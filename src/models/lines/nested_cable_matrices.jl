struct NestedCableFrequencyState
    frequency_hz::Float64
    angular_frequency_rad_s::Float64
    unreduced_series_impedance_matrix_ohm_per_m::Matrix{ComplexF64}
    series_impedance_matrix_ohm_per_m::Matrix{ComplexF64}
    unreduced_shunt_admittance_matrix_s_per_m::Matrix{ComplexF64}
    shunt_admittance_matrix_s_per_m::Matrix{ComplexF64}
    phase_zy_matrix_per_m2::Matrix{ComplexF64}
    modal_solution::LineModalSolution
    modal_series_impedance_ohm_per_m::Vector{ComplexF64}
    modal_shunt_admittance_s_per_m::Vector{ComplexF64}
    modal_characteristic_impedance_ohm::Vector{ComplexF64}
    modal_characteristic_admittance_s::Vector{ComplexF64}
    modal_attenuation_db_per_km::Vector{Float64}
    modal_velocity_m_per_s::Vector{Float64}
    grounded_conductor_reduction_executed::Bool
    internal_impedance_executed::Bool
    earth_return_impedance_executed::Bool
    electrostatic_admittance_executed::Bool
    modal_solution_executed::Bool
    complex_symmetry_max_abs_error::Float64
    capacitance_minimum_eigenvalue_f_per_m::Float64
    modal_diagonalization_max_abs_error::Float64
    grounded_admittance_balance_passed::Bool
    physical_checks_passed::Bool
    conductor_count::Int
    active_conductor_count::Int
    phase_count::Int
end

mutable struct NestedCableTransientLineState
    frequency_states::Vector{NestedCableFrequencyState}
    recursive_state::FrequencyDependentLineRecursiveConvolutionState
    response_fit::LineRecursiveConvolutionFitResult
    sample_rows::Vector{Vector{LineFrequencyPoint}}
    line_length_m::Float64
    dt_s::Float64
    modal_diagonalization_max_abs_error::Float64
    fit_relative_max_abs_error::Float64
    stable_poles::Bool
    physical_checks_passed::Bool
end

function cable_underground_homogeneous_earth_return_impedance(
    direct_distance_m::Real,
    image_distance_m::Real,
    image_angle_rad::Real,
    earth_resistivity_ohm_m::Real,
    frequency_hz::Real,
)
    direct_distance = _checked_line_positive(direct_distance_m, "earth-return direct_distance_m")
    image_distance = _checked_line_positive(image_distance_m, "earth-return image_distance_m")
    image_distance > direct_distance ||
        throw(ArgumentError("earth-return image distance must exceed direct distance"))
    earth_resistivity = _checked_line_positive(earth_resistivity_ohm_m, "earth_resistivity_ohm_m")
    frequency = _checked_line_positive(frequency_hz, "earth-return frequency_hz")
    omega = 2.0 * pi * frequency
    sqrt_jw = sqrt(ComplexF64(0.0, omega))
    diffusion_scale = sqrt(LINE_VACUUM_PERMEABILITY_H_PER_M / earth_resistivity)
    direct_argument = direct_distance * diffusion_scale * sqrt_jw
    image_argument = image_distance * diffusion_scale * sqrt_jw
    scaled = abs(image_argument) > 10.0
    _, _, direct_k0, _ = _line_skin_effect_bessel_i0_i1_k0_k1(direct_argument, scaled)
    _, _, image_k0, _ = _line_skin_effect_bessel_i0_i1_k0_k1(image_argument, scaled)
    impedance_scale = ComplexF64(0.0, omega) *
        (LINE_VACUUM_PERMEABILITY_H_PER_M / (2.0 * pi))
    correction = if abs(image_argument) > 100.0
        abs(direct_argument) > LINE_SKIN_EFFECT_EXPONENT_LIMIT ?
            0.0 + 0.0im :
            impedance_scale * direct_k0 / exp(direct_argument)
    elseif scaled
        exponent = image_argument - direct_argument
        impedance_scale * (direct_k0 - image_k0 / exp(exponent)) / exp(direct_argument)
    else
        impedance_scale * (direct_k0 - image_k0)
    end
    return ComplexF64(
        correction + cable_homogeneous_earth_return_impedance(
            image_distance,
            image_angle_rad,
            earth_resistivity,
            frequency,
        ),
    )
end

function _nested_cable_hollow_conductor_impedances(
    inner_diffusion_factor::Real,
    outer_diffusion_factor::Real,
    relative_permeability::Real,
    resistivity_ohm_m::Real,
    inner_radius_m::Real,
    outer_radius_m::Real,
    frequency_hz::Real,
)
    inner = _checked_line_positive(inner_diffusion_factor, "hollow-conductor inner diffusion factor")
    outer = _checked_line_positive(outer_diffusion_factor, "hollow-conductor outer diffusion factor")
    inner < outer || throw(ArgumentError("hollow-conductor inner diffusion factor must be less than outer"))
    permeability = _checked_line_positive(relative_permeability, "hollow-conductor relative_permeability")
    resistivity = _checked_line_positive(resistivity_ohm_m, "hollow-conductor resistivity_ohm_m")
    inner_radius = _checked_line_positive(inner_radius_m, "hollow-conductor inner_radius_m")
    outer_radius = _checked_line_positive(outer_radius_m, "hollow-conductor outer_radius_m")
    inner_radius < outer_radius || throw(ArgumentError("hollow-conductor inner radius must be less than outer radius"))
    frequency = _checked_line_positive(frequency_hz, "hollow-conductor frequency_hz")

    omega = 2.0 * pi * frequency
    jw = ComplexF64(0.0, omega)
    sqrt_jw = sqrt(jw)
    inner_argument = ComplexF64(inner) * sqrt_jw
    outer_argument = ComplexF64(outer) * sqrt_jw
    scaled = abs(inner_argument) > 10.0
    inner_i0, inner_i1, inner_k0, inner_k1 =
        _line_skin_effect_bessel_i0_i1_k0_k1(inner_argument, scaled)
    outer_i0, outer_i1, outer_k0, outer_k1 =
        _line_skin_effect_bessel_i0_i1_k0_k1(outer_argument, scaled)
    surface_scale = jw * (LINE_VACUUM_PERMEABILITY_H_PER_M / (2.0 * pi)) * permeability

    if scaled
        exponent = outer_argument - inner_argument
        if abs(exponent) > LINE_SKIN_EFFECT_EXPONENT_LIMIT
            inner_impedance = surface_scale / inner_argument * inner_k0 / inner_k1
            outer_impedance = surface_scale / outer_argument * outer_i0 / outer_i1
            return ComplexF64(inner_impedance), ComplexF64(outer_impedance), 0.0 + 0.0im
        end
        scale = exp(exponent)
        denominator = outer_i1 * inner_k1 * scale - inner_i1 * outer_k1 / scale
        inner_numerator = inner_i0 * outer_k1 / scale + inner_k0 * outer_i1 * scale
        outer_numerator = outer_i0 * inner_k1 * scale + outer_k0 * inner_i1 / scale
    else
        denominator = outer_i1 * inner_k1 - inner_i1 * outer_k1
        inner_numerator = inner_i0 * outer_k1 + inner_k0 * outer_i1
        outer_numerator = outer_i0 * inner_k1 + outer_k0 * inner_i1
    end
    denominator != 0.0 + 0.0im || throw(ArgumentError("hollow-conductor Bessel denominator must be nonzero"))
    inner_impedance = surface_scale / inner_argument * inner_numerator / denominator
    outer_impedance = surface_scale / outer_argument * outer_numerator / denominator
    mutual_impedance = resistivity / (2.0 * pi * inner_radius * outer_radius) / denominator
    return (
        ComplexF64(inner_impedance),
        ComplexF64(outer_impedance),
        ComplexF64(mutual_impedance),
    )
end

function _nested_cable_nonnegative_resistance(value::ComplexF64)
    return real(value) < 0.0 ? ComplexF64(0.0, imag(value)) : value
end

function _nested_cable_internal_impedance_matrix(
    state::CablePipeSheathDerivedState,
    radii::Matrix{Float64},
    resistivity::Matrix{Float64},
    conductor_permeability::Matrix{Float64},
    insulation_permeability::Matrix{Float64},
    frequency::Float64,
)
    phases = state.phase_count
    total = state.conductor_count
    internal = zeros(ComplexF64, total, total)
    omega = 2.0 * pi * frequency
    dielectric_scale = ComplexF64(0.0, omega * LINE_VACUUM_PERMEABILITY_H_PER_M / (2.0 * pi))
    for phase in 1:phases
        core = phase
        sheath = phases + phase
        armor = 2 * phases + phase
        core_impedance = cable_skin_effect_internal_impedance(
            state.core_inner_diffusion_factors[phase],
            state.core_outer_diffusion_factors[phase],
            conductor_permeability[phase, 1],
            frequency,
        )
        core_insulation = dielectric_scale * insulation_permeability[phase, 1] *
            state.core_insulation_log_ratios[phase]
        sheath_inner, sheath_outer, sheath_mutual =
            _nested_cable_hollow_conductor_impedances(
                state.sheath_inner_diffusion_factors[phase],
                state.sheath_outer_diffusion_factors[phase],
                conductor_permeability[phase, 2],
                resistivity[phase, 2],
                radii[phase, 3],
                radii[phase, 4],
                frequency,
            )
        sheath_insulation = dielectric_scale * insulation_permeability[phase, 2] *
            state.sheath_insulation_log_ratios[phase]
        armor_inner, armor_outer, armor_mutual =
            _nested_cable_hollow_conductor_impedances(
                state.armor_inner_diffusion_factors[phase],
                state.armor_outer_diffusion_factors[phase],
                conductor_permeability[phase, 3],
                resistivity[phase, 3],
                radii[phase, 5],
                radii[phase, 6],
                frequency,
            )
        inner_path = core_impedance + core_insulation + sheath_inner
        middle_path = sheath_outer + sheath_insulation + armor_inner
        # The external field starts at the outer cable boundary used by the
        # underground earth-return term, so the armor-to-boundary insulation
        # must not be counted again in the concentric internal impedance.
        outer_path = armor_outer
        armor_loop = outer_path - 2.0 * armor_mutual
        core_armor = _nested_cable_nonnegative_resistance(outer_path - armor_mutual)
        sheath_loop = _nested_cable_nonnegative_resistance(middle_path + armor_loop)
        core_sheath = _nested_cable_nonnegative_resistance(sheath_loop - sheath_mutual)

        internal[core, core] = inner_path + middle_path - 2.0 * sheath_mutual + armor_loop
        internal[sheath, sheath] = sheath_loop
        internal[core, sheath] = core_sheath
        internal[sheath, core] = core_sheath
        internal[armor, armor] = outer_path
        internal[core, armor] = core_armor
        internal[armor, core] = core_armor
        internal[sheath, armor] = core_armor
        internal[armor, sheath] = core_armor
    end
    return internal
end

function _nested_cable_earth_impedance_matrix(
    geometry::CableGeometryConstants,
    total_conductor_count::Int,
    frequency::Float64,
    earth_resistivity::Float64,
)
    phase_earth = Matrix{ComplexF64}(undef, geometry.conductor_count, geometry.conductor_count)
    for row in axes(phase_earth, 1), column in axes(phase_earth, 2)
        phase_earth[row, column] = cable_underground_homogeneous_earth_return_impedance(
            geometry.direct_distance_m[row, column],
            geometry.image_distance_m[row, column],
            geometry.angle_rad[row, column],
            earth_resistivity,
            frequency,
        )
    end
    phases = geometry.phase_count
    earth = Matrix{ComplexF64}(undef, total_conductor_count, total_conductor_count)
    for row in 1:total_conductor_count, column in 1:total_conductor_count
        earth[row, column] = phase_earth[mod1(row, phases), mod1(column, phases)]
    end
    return earth
end

function _nested_cable_admittance_matrix(
    radii::Matrix{Float64},
    insulation_permittivity::Matrix{Float64},
    frequency::Float64,
)
    phases = size(radii, 1)
    total = 3 * phases
    capacitance = zeros(Float64, total, total)
    for phase in 1:phases
        core = phase
        sheath = phases + phase
        armor = 2 * phases + phase
        layer_capacitances = (
            2.0 * pi * LINE_VACUUM_PERMITTIVITY_F_PER_M * insulation_permittivity[phase, 1] /
                log(radii[phase, 3] / radii[phase, 2]),
            2.0 * pi * LINE_VACUUM_PERMITTIVITY_F_PER_M * insulation_permittivity[phase, 2] /
                log(radii[phase, 5] / radii[phase, 4]),
            2.0 * pi * LINE_VACUUM_PERMITTIVITY_F_PER_M * insulation_permittivity[phase, 3] /
                log(radii[phase, 7] / radii[phase, 6]),
        )
        for (first, second, value) in (
            (core, sheath, layer_capacitances[1]),
            (sheath, armor, layer_capacitances[2]),
        )
            capacitance[first, first] += value
            capacitance[second, second] += value
            capacitance[first, second] -= value
            capacitance[second, first] -= value
        end
        capacitance[armor, armor] += layer_capacitances[3]
    end
    return ComplexF64.(0.0, 2.0 * pi * frequency .* capacitance)
end

function _nested_cable_modal_quantities(
    series::Matrix{ComplexF64},
    shunt::Matrix{ComplexF64},
    frequency::Float64,
)
    zy = series * shunt
    raw_solution = line_modal_solution(zy, frequency)
    modal_to_phase = copy(raw_solution.transform.modal_to_phase)
    for mode in axes(modal_to_phase, 2)
        maximum_component = maximum(abs.(modal_to_phase[:, mode]))
        maximum_component > 0.0 ||
            throw(ArgumentError("nested cable modal vector must be nonzero"))
        modal_to_phase[:, mode] ./= maximum_component
    end
    transform = LineModalTransform(inv(modal_to_phase), modal_to_phase)
    eigen_residual = zy * modal_to_phase -
        modal_to_phase * Diagonal(raw_solution.eigenvalues)
    inverse_product = transform.modal_to_phase * transform.phase_to_modal
    solution = LineModalSolution(
        transform,
        raw_solution.eigenvalues,
        raw_solution.propagation_roots,
        raw_solution.modal_eigen_order,
        raw_solution.mode_sequence,
        raw_solution.mode_metrics,
        raw_solution.ordered_mode_metrics,
        maximum(abs.(eigen_residual)),
        maximum(abs.(inverse_product - I)),
        raw_solution.update_count,
    )
    modal_shunt_matrix = transpose(modal_to_phase) * shunt * modal_to_phase
    modal_shunt = ComplexF64[modal_shunt_matrix[index, index] for index in axes(modal_shunt_matrix, 1)]
    modal_series = ComplexF64[
        solution.eigenvalues[index] / modal_shunt[index] for index in eachindex(modal_shunt)
    ]
    characteristic_admittance = sqrt.(modal_shunt ./ modal_series)
    characteristic_impedance = inv.(characteristic_admittance)
    attenuation = 20.0 * 1000.0 / log(10.0) .* real.(solution.propagation_roots)
    omega = 2.0 * pi * frequency
    velocity = omega ./ imag.(solution.propagation_roots)
    diagonalization_error = _line_offdiagonal_max_abs(modal_shunt_matrix)
    return (
        zy,
        solution,
        modal_series,
        modal_shunt,
        ComplexF64.(characteristic_impedance),
        ComplexF64.(characteristic_admittance),
        Float64.(attenuation),
        Float64.(velocity),
        diagonalization_error,
    )
end

function nested_cable_frequency_state(
    state::CablePipeSheathDerivedState,
    geometry::CableGeometryConstants,
    boundary_radii_m::AbstractMatrix,
    resistivity_ohm_m::AbstractMatrix,
    conductor_relative_permeability::AbstractMatrix,
    insulation_relative_permeability::AbstractMatrix,
    insulation_relative_permittivity::AbstractMatrix,
    frequency_hz::Real,
    earth_resistivity_ohm_m::Real;
    layer_counts::AbstractVector{<:Integer} = fill(3, state.phase_count),
)
    state.cable_kind == :non_pipe_cable ||
        throw(ArgumentError("nested cable frequency state currently requires single-core concentric cables"))
    phases = state.phase_count
    total = state.conductor_count
    total == 3 * phases || throw(ArgumentError("nested cable frequency state requires core, sheath, and armor per phase"))
    state.selected_grounded_conductor_count == 2 * phases ||
        throw(ArgumentError("nested cable frequency state requires grounded armor conductors"))
    length(layer_counts) == phases && all(==(3), layer_counts) ||
        throw(ArgumentError("nested cable frequency state requires three layers per phase"))
    geometry.phase_count == phases && geometry.conductor_count == phases ||
        throw(ArgumentError("nested cable outer geometry must contain one conductor per phase"))
    frequency = _checked_line_positive(frequency_hz, "nested cable frequency_hz")
    earth_resistivity = _checked_line_positive(earth_resistivity_ohm_m, "nested cable earth_resistivity_ohm_m")
    radii = _checked_line_real_matrix(boundary_radii_m, phases, 7, "boundary_radii_m")
    resistivity = _checked_line_real_matrix(resistivity_ohm_m, phases, 3, "resistivity_ohm_m")
    conductor_permeability = _checked_line_real_matrix(
        conductor_relative_permeability,
        phases,
        3,
        "conductor_relative_permeability",
    )
    insulation_permeability = _checked_line_real_matrix(
        insulation_relative_permeability,
        phases,
        3,
        "insulation_relative_permeability",
    )
    insulation_permittivity = _checked_line_real_matrix(
        insulation_relative_permittivity,
        phases,
        3,
        "insulation_relative_permittivity",
    )
    all(>(0.0), radii) && all(>(0.0), resistivity) &&
        all(>(0.0), conductor_permeability) && all(>(0.0), insulation_permeability) &&
        all(>(0.0), insulation_permittivity) ||
        throw(ArgumentError("nested cable radii and material values must be positive"))
    all(phase -> all(diff(radii[phase, :]) .> 0.0), 1:phases) ||
        throw(ArgumentError("nested cable boundary radii must be strictly increasing"))

    internal = _nested_cable_internal_impedance_matrix(
        state,
        radii,
        resistivity,
        conductor_permeability,
        insulation_permeability,
        frequency,
    )
    earth = _nested_cable_earth_impedance_matrix(
        geometry,
        total,
        frequency,
        earth_resistivity,
    )
    unreduced_series = internal + earth
    active = state.selected_grounded_conductor_count
    active_indices = 1:active
    grounded_indices = (active + 1):total
    reduced_series = Matrix{ComplexF64}(
        unreduced_series[active_indices, active_indices] -
        unreduced_series[active_indices, grounded_indices] *
        (unreduced_series[grounded_indices, grounded_indices] \
            unreduced_series[grounded_indices, active_indices]),
    )
    unreduced_shunt = _nested_cable_admittance_matrix(radii, insulation_permittivity, frequency)
    reduced_shunt = Matrix{ComplexF64}(unreduced_shunt[active_indices, active_indices])
    zy, solution, modal_series, modal_shunt, characteristic_impedance,
        characteristic_admittance, attenuation, velocity, diagonalization_error =
        _nested_cable_modal_quantities(reduced_series, reduced_shunt, frequency)

    symmetry_error = maximum((
        maximum(abs.(unreduced_series - transpose(unreduced_series))),
        maximum(abs.(reduced_series - transpose(reduced_series))),
        maximum(abs.(reduced_shunt - transpose(reduced_shunt))),
    ))
    capacitance = Matrix{Float64}(imag.(reduced_shunt) ./ (2.0 * pi * frequency))
    capacitance_minimum_eigenvalue = minimum(eigvals(Symmetric(capacitance)))
    core_row_sum_error = maximum(abs.(vec(sum(capacitance[1:phases, :]; dims = 2))))
    sheath_row_sums = vec(sum(capacitance[(phases + 1):active, :]; dims = 2))
    admittance_balance_passed = core_row_sum_error <= 1.0e-18 && all(>(0.0), sheath_row_sums)
    finite_modal = all(isfinite, attenuation) && all(isfinite, velocity) && all(>(0.0), velocity)
    physical_checks = symmetry_error <= 1.0e-12 &&
        capacitance_minimum_eigenvalue >= -1.0e-18 &&
        solution.eigenvector_residual_max_abs_error <= 1.0e-12 &&
        diagonalization_error <= 1.0e-12 && admittance_balance_passed && finite_modal
    return NestedCableFrequencyState(
        frequency,
        2.0 * pi * frequency,
        unreduced_series,
        reduced_series,
        unreduced_shunt,
        reduced_shunt,
        zy,
        solution,
        modal_series,
        modal_shunt,
        characteristic_impedance,
        characteristic_admittance,
        attenuation,
        velocity,
        true,
        true,
        true,
        true,
        true,
        symmetry_error,
        capacitance_minimum_eigenvalue,
        diagonalization_error,
        admittance_balance_passed,
        physical_checks,
        total,
        active,
        phases,
    )
end

function _nested_cable_transient_frequency_states(
    frequency_states::AbstractVector{NestedCableFrequencyState},
)
    length(frequency_states) >= 2 ||
        throw(ArgumentError("nested cable transient fitting requires at least two frequency states"))
    frequencies = Float64[state.frequency_hz for state in frequency_states]
    order = sortperm(frequencies)
    sorted = NestedCableFrequencyState[frequency_states[index] for index in order]
    sorted_frequencies = frequencies[order]
    all(diff(sorted_frequencies) .> 0.0) ||
        throw(ArgumentError("nested cable transient frequencies must be unique"))
    active_count = first(sorted).active_conductor_count
    phase_count = first(sorted).phase_count
    all(
        state -> state.active_conductor_count == active_count &&
            state.phase_count == phase_count &&
            size(state.phase_zy_matrix_per_m2) == (active_count, active_count) &&
            state.physical_checks_passed,
        sorted,
    ) || throw(ArgumentError(
        "nested cable transient frequency states must be physically accepted and dimensionally consistent",
    ))
    return sorted, sorted_frequencies, active_count
end

function _nested_cable_transient_modal_samples(
    states::Vector{NestedCableFrequencyState},
    frequencies::Vector{Float64},
    line_length_m::Float64,
)
    matrices = Matrix{ComplexF64}[copy(state.phase_zy_matrix_per_m2) for state in states]
    scan = line_modal_solution_scan(matrices, frequencies)
    mode_count = scan.mode_count
    sample_rows = Vector{Vector{LineFrequencyPoint}}(undef, length(states))
    response_samples = Matrix{ComplexF64}(undef, mode_count, length(states))
    maximum_diagonalization_error = 0.0
    for frequency_index in eachindex(states)
        state = states[frequency_index]
        solution = scan.solutions[frequency_index]
        sequence = solution.mode_sequence
        maximum_diagonalization_error = max(
            maximum_diagonalization_error,
            state.modal_diagonalization_max_abs_error,
        )
        rows = LineFrequencyPoint[]
        sizehint!(rows, mode_count)
        for mode in 1:mode_count
            source_mode = sequence[mode]
            characteristic_impedance =
                state.modal_characteristic_impedance_ohm[source_mode]
            propagation_constant =
                state.modal_solution.propagation_roots[source_mode]
            propagation_factor = exp(-propagation_constant * line_length_m)
            push!(
                rows,
                LineFrequencyPoint(
                    frequencies[frequency_index],
                    characteristic_impedance,
                    propagation_constant,
                    propagation_factor,
                ),
            )
            response_samples[mode, frequency_index] =
                propagation_factor / characteristic_impedance
        end
        sample_rows[frequency_index] = rows
    end
    return matrices, sample_rows, response_samples, maximum_diagonalization_error, scan
end

function nested_cable_transient_line_state(
    frequency_states::AbstractVector{NestedCableFrequencyState},
    line_length_m::Real,
    pole_decay::AbstractMatrix,
    dt_s::Real;
    initial_frequency_hz::Real = minimum(state.frequency_hz for state in frequency_states),
)
    line_length = _checked_line_positive_finite(line_length_m, "nested cable line_length_m")
    dt = _checked_line_positive_finite(dt_s, "nested cable dt_s")
    states, frequencies, mode_count =
        _nested_cable_transient_frequency_states(frequency_states)
    matrices, sample_rows, response_samples, diagonalization_error, _ =
        _nested_cable_transient_modal_samples(states, frequencies, line_length)
    decay = _checked_line_complex_rectangular_matrix(
        pole_decay,
        mode_count,
        "nested cable pole decay",
    )
    fit = frequency_dependent_line_recursive_convolution_fit(
        frequencies,
        response_samples,
        decay,
        dt,
    )
    response_scale = maximum(abs.(response_samples); init = 0.0)
    relative_error = response_scale == 0.0 ? fit.max_abs_error :
        fit.max_abs_error / response_scale
    stable_poles = all(value -> abs(value) <= 1.0 + 64.0 * eps(Float64), decay)
    recursive_state = FrequencyDependentLineRecursiveConvolutionState(
        matrices,
        sample_rows,
        initial_frequency_hz,
        line_length,
        fit;
        skin_effect_internal_impedance_executed = true,
        earth_return_impedance_executed = true,
        frequency_dependent_fitting_executed = true,
        frequency_loop_executed = true,
        pipe_sheath_side_effects_executed = true,
    )
    physical_checks = all(state -> state.physical_checks_passed, states) &&
        stable_poles && isfinite(relative_error) &&
        diagonalization_error <= 1.0e-10
    return NestedCableTransientLineState(
        states,
        recursive_state,
        fit,
        sample_rows,
        line_length,
        dt,
        diagonalization_error,
        relative_error,
        stable_poles,
        physical_checks,
    )
end

function nested_cable_semlyen_frequency_dependent_line_from_fit(
    frequency_states::AbstractVector{NestedCableFrequencyState},
    line_length_m::Real,
    from_nodes::AbstractVector{<:Integer},
    to_nodes::AbstractVector{<:Integer},
    timestep_s::Real,
    propagation_fits::AbstractVector{LineStepResponseExponentialFitResult},
    admittance_fits::AbstractVector{LineStepResponseExponentialFitResult},
    characteristic_admittance_s::AbstractVector{<:Real};
    phasor_frequency_hz::Real,
    propagation_relative_tolerance::Real = 0.10,
    admittance_relative_tolerance::Real = 0.10,
)
    line_length = _checked_line_positive_finite(line_length_m, "nested cable line_length_m")
    states, frequencies, _ = _nested_cable_transient_frequency_states(frequency_states)
    _, sample_rows, _, diagonalization_error, modal_scan =
        _nested_cable_transient_modal_samples(states, frequencies, line_length)
    phasor_frequency = _checked_line_positive_finite(
        phasor_frequency_hz,
        "nested cable phasor_frequency_hz",
    )
    phasor_index = findfirst(
        frequency -> isapprox(
            frequency,
            phasor_frequency;
            atol = _line_frequency_row_tolerance(frequency),
            rtol = 0.0,
        ),
        frequencies,
    )
    phasor_index === nothing && throw(ArgumentError(
        "nested cable Semlyen phasor frequency must be an exact generated scan row",
    ))
    diagonalization_error <= 1.0e-10 || throw(ArgumentError(
        "nested cable modal scan exceeds the accepted diagonalization tolerance",
    ))
    voltage_modal_to_phase =
        modal_scan.solutions[phasor_index].transform.modal_to_phase
    current_modal_to_phase = transpose(inv(voltage_modal_to_phase))
    return semlyen_frequency_dependent_line_from_scan_fit(
        from_nodes,
        to_nodes,
        sample_rows,
        line_length,
        voltage_modal_to_phase,
        current_modal_to_phase,
        timestep_s,
        propagation_fits,
        admittance_fits,
        characteristic_admittance_s;
        phasor_frequency_hz = phasor_frequency,
        propagation_relative_tolerance = propagation_relative_tolerance,
        admittance_relative_tolerance = admittance_relative_tolerance,
    )
end

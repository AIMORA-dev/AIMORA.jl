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

function _nested_cable_layer_layout(
    layer_counts::AbstractVector{<:Integer};
    include_pipe::Bool=false,
)
    counts = Int.(layer_counts)
    all(count -> 1 <= count <= 3, counts) ||
        throw(ArgumentError("nested cable layer counts must be between one and three"))
    indices = [zeros(Int, count) for count in counts]
    phase_by_conductor = Int[]
    next_index = 1
    for layer in 1:3, phase in eachindex(counts)
        counts[phase] >= layer || continue
        indices[phase][layer] = next_index
        push!(phase_by_conductor, phase)
        next_index += 1
    end
    pipe_index = include_pipe ? next_index : 0
    include_pipe && push!(phase_by_conductor, 0)
    return (
        layer_indices = indices,
        phase_by_conductor = phase_by_conductor,
        cable_conductor_count = next_index - 1,
        total_conductor_count = next_index - 1 + (include_pipe ? 1 : 0),
        pipe_index = pipe_index,
    )
end

function _nested_cable_internal_impedance_matrix(
    state::CablePipeSheathDerivedState,
    radii::Matrix{Float64},
    resistivity::Matrix{Float64},
    conductor_permeability::Matrix{Float64},
    insulation_permeability::Matrix{Float64},
    frequency::Float64,
    layer_counts::Vector{Int},
    layout,
)
    internal = zeros(
        ComplexF64,
        layout.total_conductor_count,
        layout.total_conductor_count,
    )
    omega = 2.0 * pi * frequency
    magnetic_scale =
        ComplexF64(0.0, omega * LINE_VACUUM_PERMEABILITY_H_PER_M / (2.0 * pi))
    for phase in 1:state.phase_count
        count = layer_counts[phase]
        core = layout.layer_indices[phase][1]
        core_impedance = cable_skin_effect_internal_impedance(
            state.core_inner_diffusion_factors[phase],
            state.core_outer_diffusion_factors[phase],
            conductor_permeability[phase, 1],
            frequency,
        )
        core_insulation = magnetic_scale * insulation_permeability[phase, 1] *
            state.core_insulation_log_ratios[phase]
        if count == 1
            internal[core, core] = core_impedance + core_insulation
            continue
        end
        sheath = layout.layer_indices[phase][2]
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
        sheath_insulation = magnetic_scale * insulation_permeability[phase, 2] *
            state.sheath_insulation_log_ratios[phase]
        if count == 2
            outer_path = sheath_outer + sheath_insulation
            core_sheath =
                _nested_cable_nonnegative_resistance(outer_path - sheath_mutual)
            internal[core, core] =
                core_impedance + core_insulation + sheath_inner +
                outer_path - 2.0 * sheath_mutual
            internal[sheath, sheath] = outer_path
            internal[core, sheath] = core_sheath
            internal[sheath, core] = core_sheath
            continue
        end
        armor = layout.layer_indices[phase][3]
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
    phase_by_conductor::Vector{Int},
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
    total_conductor_count = length(phase_by_conductor)
    earth = Matrix{ComplexF64}(undef, total_conductor_count, total_conductor_count)
    for row in 1:total_conductor_count, column in 1:total_conductor_count
        row_phase = phase_by_conductor[row]
        column_phase = phase_by_conductor[column]
        row_geometry =
            geometry.conductor_count == 1 ? 1 :
            row_phase == 0 ? geometry.conductor_count : row_phase
        column_geometry =
            geometry.conductor_count == 1 ? 1 :
            column_phase == 0 ? geometry.conductor_count : column_phase
        1 <= row_geometry <= geometry.conductor_count &&
            1 <= column_geometry <= geometry.conductor_count ||
            throw(ArgumentError("nested cable earth-return geometry does not cover every conductor phase"))
        earth[row, column] = phase_earth[row_geometry, column_geometry]
    end
    return earth
end

function _nested_pipe_cavity_log_matrix(
    state::CablePipeSheathDerivedState,
    radii::Matrix{Float64},
    layer_counts::Vector{Int},
)
    phase_count = state.phase_count
    pipe_radius = state.pipe_radii_m[1]
    coordinates = ComplexF64[
        state.conductor_pipe_center_distances_m[phase] *
        cis(state.conductor_angles_rad[phase])
        for phase in 1:phase_count
    ]
    outer_boundaries = Float64[
        radii[phase, 2 * layer_counts[phase] + 1] for phase in 1:phase_count
    ]
    all(
        abs(coordinates[phase]) + outer_boundaries[phase] < pipe_radius
        for phase in 1:phase_count
    ) || throw(ArgumentError("cable outer boundaries must lie inside the pipe"))
    values = zeros(Float64, phase_count, phase_count)
    for row in 1:phase_count, column in 1:phase_count
        if row == column
            numerator = pipe_radius^2 - abs2(coordinates[row])
            denominator = pipe_radius * outer_boundaries[row]
        else
            numerator =
                abs(pipe_radius^2 - coordinates[row] * conj(coordinates[column]))
            denominator =
                pipe_radius * abs(coordinates[row] - coordinates[column])
        end
        numerator > denominator > 0.0 ||
            throw(ArgumentError("pipe-cavity potential coefficient must be positive"))
        values[row, column] = log(numerator / denominator)
    end
    return 0.5 .* (values .+ transpose(values))
end

function _nested_pipe_wall_proximity_impedance_matrix(
    state::CablePipeSheathDerivedState,
    frequency::Float64;
    term_count::Int=19,
)
    term_count >= 1 ||
        throw(ArgumentError("pipe-wall proximity term count must be positive"))
    pipe_radius = state.pipe_radii_m[1]
    normalized_distances =
        state.conductor_pipe_center_distances_m ./ pipe_radius
    all(distance -> 0.0 <= distance < 1.0, normalized_distances) ||
        throw(ArgumentError("pipe-wall proximity distances must lie inside the pipe"))
    omega = 2.0 * pi * frequency
    sqrt_jw = sqrt(ComplexF64(0.0, omega))
    inner_diffusion_factor =
        pipe_radius *
        sqrt(
            LINE_VACUUM_PERMEABILITY_H_PER_M /
            state.pipe_resistivity_ohm_m *
            state.pipe_relative_permeability,
        )
    inner_argument = ComplexF64(inner_diffusion_factor) * sqrt_jw
    scaled = abs(inner_argument) > 10.0
    _, _, k0, k1 =
        _line_skin_effect_bessel_i0_i1_k0_k1(inner_argument, scaled)
    k1 != 0.0 + 0.0im ||
        throw(ArgumentError("pipe-wall proximity Bessel K1 must be nonzero"))
    bessel_ratios = Vector{ComplexF64}(undef, term_count)
    ratio = k0 / k1
    for order in 1:term_count
        bessel_ratios[order] = ratio
        ratio = inv(ratio + 2.0 * order / inner_argument)
    end

    phase_count = state.phase_count
    correction = zeros(ComplexF64, phase_count, phase_count)
    for row in 1:phase_count, column in row:phase_count
        radial_product =
            normalized_distances[row] * normalized_distances[column]
        angle =
            state.conductor_angles_rad[column] -
            state.conductor_angles_rad[row]
        value = 0.0 + 0.0im
        radial_power = radial_product
        for order in 1:term_count
            denominator =
                order * (state.pipe_relative_permeability + 1.0) +
                inner_argument * bessel_ratios[order]
            denominator != 0.0 + 0.0im ||
                throw(ArgumentError("pipe-wall proximity denominator must be nonzero"))
            value += radial_power * cos(order * angle) / denominator
            radial_power *= radial_product
        end
        correction[row, column] = value
        correction[column, row] = value
    end
    magnetic_scale =
        ComplexF64(
            0.0,
            omega * LINE_VACUUM_PERMEABILITY_H_PER_M / (2.0 * pi),
        )
    return (
        2.0 *
        state.pipe_relative_permeability *
        magnetic_scale .*
        correction
    )
end

function _nested_pipe_series_coupling!(
    series::Matrix{ComplexF64},
    state::CablePipeSheathDerivedState,
    radii::Matrix{Float64},
    resistivity::Matrix{Float64},
    conductor_permeability::Matrix{Float64},
    insulation_permeability::Matrix{Float64},
    frequency::Float64,
    layer_counts::Vector{Int},
    layout,
)
    cavity = _nested_pipe_cavity_log_matrix(state, radii, layer_counts)
    wall_proximity =
        _nested_pipe_wall_proximity_impedance_matrix(state, frequency)
    omega = 2.0 * pi * frequency
    magnetic_scale =
        ComplexF64(0.0, omega * LINE_VACUUM_PERMEABILITY_H_PER_M / (2.0 * pi))
    for row in 1:layout.cable_conductor_count,
        column in 1:layout.cable_conductor_count
        row_phase = layout.phase_by_conductor[row]
        column_phase = layout.phase_by_conductor[column]
        value = cavity[row_phase, column_phase]
        if row_phase == column_phase && layer_counts[row_phase] == 3
            layer = layer_counts[row_phase]
            value += insulation_permeability[row_phase, layer] *
                log(
                    radii[row_phase, 2 * layer + 1] /
                    radii[row_phase, 2 * layer],
                )
        end
        series[row, column] +=
            magnetic_scale * value +
            wall_proximity[row_phase, column_phase]
    end
    pipe_inner_factor =
        state.pipe_radii_m[1] *
        sqrt(
            LINE_VACUUM_PERMEABILITY_H_PER_M /
            state.pipe_resistivity_ohm_m *
            state.pipe_relative_permeability,
        )
    pipe_outer_factor =
        state.pipe_radii_m[2] *
        sqrt(
            LINE_VACUUM_PERMEABILITY_H_PER_M /
            state.pipe_resistivity_ohm_m *
            state.pipe_relative_permeability,
        )
    pipe_inner, pipe_outer, pipe_mutual =
        _nested_cable_hollow_conductor_impedances(
            pipe_inner_factor,
            pipe_outer_factor,
            state.pipe_relative_permeability,
            state.pipe_resistivity_ohm_m,
            state.pipe_radii_m[1],
            state.pipe_radii_m[2],
            frequency,
        )
    pipe_outer_insulation =
        magnetic_scale *
        log(state.pipe_radii_m[3] / state.pipe_radii_m[2])
    cable_common =
        pipe_inner + pipe_outer - 2.0 * pipe_mutual +
        pipe_outer_insulation
    cable_pipe =
        pipe_outer - pipe_mutual + pipe_outer_insulation
    for row in 1:layout.cable_conductor_count,
        column in 1:layout.cable_conductor_count
        series[row, column] += cable_common
    end
    pipe_index = layout.pipe_index
    for conductor in 1:layout.cable_conductor_count
        series[conductor, pipe_index] += cable_pipe
        series[pipe_index, conductor] += cable_pipe
    end
    series[pipe_index, pipe_index] += pipe_outer + pipe_outer_insulation
    return series
end

function _nested_cable_admittance_matrix(
    state::CablePipeSheathDerivedState,
    radii::Matrix{Float64},
    insulation_permittivity::Matrix{Float64},
    frequency::Float64,
    layer_counts::Vector{Int},
    layout,
)
    capacitance =
        zeros(Float64, layout.total_conductor_count, layout.total_conductor_count)
    outer_indices = Int[]
    for phase in 1:state.phase_count
        count = layer_counts[phase]
        phase_indices = layout.layer_indices[phase]
        for layer in 1:(count - 1)
            first = phase_indices[layer]
            second = phase_indices[layer + 1]
            value =
                2.0 * pi * LINE_VACUUM_PERMITTIVITY_F_PER_M *
                insulation_permittivity[phase, layer] /
                log(radii[phase, 2 * layer + 1] / radii[phase, 2 * layer])
            capacitance[first, first] += value
            capacitance[second, second] += value
            capacitance[first, second] -= value
            capacitance[second, first] -= value
        end
        push!(outer_indices, phase_indices[end])
    end
    if state.cable_kind == :pipe_type_cable
        potential =
            _nested_pipe_cavity_log_matrix(state, radii, layer_counts) ./
            state.pipe_inner_insulator_relative_permittivity
        for phase in 1:state.phase_count
            layer = layer_counts[phase]
            potential[phase, phase] +=
                log(
                    radii[phase, 2 * layer + 1] /
                    radii[phase, 2 * layer],
                ) / insulation_permittivity[phase, layer]
        end
        potential ./= 2.0 * pi * LINE_VACUUM_PERMITTIVITY_F_PER_M
        outer_capacitance = inv(Symmetric(potential))
        pipe_index = layout.pipe_index
        for row in 1:state.phase_count, column in 1:state.phase_count
            value = outer_capacitance[row, column]
            capacitance[outer_indices[row], outer_indices[column]] += value
            capacitance[outer_indices[row], pipe_index] -= value
            capacitance[pipe_index, outer_indices[column]] -= value
            capacitance[pipe_index, pipe_index] += value
        end
        capacitance[pipe_index, pipe_index] +=
            2.0 * pi * LINE_VACUUM_PERMITTIVITY_F_PER_M *
            state.pipe_outer_insulator_relative_permittivity /
            log(state.pipe_radii_m[3] / state.pipe_radii_m[2])
    else
        for phase in 1:state.phase_count
            layer = layer_counts[phase]
            value =
                2.0 * pi * LINE_VACUUM_PERMITTIVITY_F_PER_M *
                insulation_permittivity[phase, layer] /
                log(
                    radii[phase, 2 * layer + 1] /
                    radii[phase, 2 * layer],
                )
            capacitance[outer_indices[phase], outer_indices[phase]] += value
        end
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
    state.cable_kind in (:non_pipe_cable, :pipe_type_cable) ||
        throw(ArgumentError("nested cable frequency state requires a non-pipe or pipe-type cable"))
    phases = state.phase_count
    counts = Int.(layer_counts)
    length(counts) == phases && all(count -> 1 <= count <= 3, counts) ||
        throw(ArgumentError("nested cable layer counts must cover every phase with one to three layers"))
    include_pipe = state.cable_kind == :pipe_type_cable
    layout = _nested_cable_layer_layout(counts; include_pipe = include_pipe)
    total = state.conductor_count
    layout.total_conductor_count == total ||
        throw(ArgumentError("nested cable conductor count does not match its layer and pipe layout"))
    if include_pipe
        state.pipe_return_included && state.pipe_count == 1 ||
            throw(ArgumentError("pipe-type cable frequency state requires one owned return pipe"))
        geometry.phase_count == 1 && geometry.conductor_count == 1 ||
            throw(ArgumentError("pipe-type cable outer geometry must contain the single enclosing pipe"))
    else
        geometry.phase_count == phases && geometry.conductor_count == phases ||
            throw(ArgumentError("non-pipe cable outer geometry must contain one conductor per phase"))
    end
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
    for phase in 1:phases
        layer_count = counts[phase]
        active_radii = @view radii[phase, 1:(2 * layer_count + 1)]
        first(active_radii) >= 0.0 && all(>(0.0), @view(active_radii[2:end])) &&
            all(diff(active_radii) .> 0.0) ||
            throw(ArgumentError("active nested cable boundary radii must be nonnegative and strictly increasing"))
        all(>(0.0), @view(resistivity[phase, 1:layer_count])) &&
            all(>(0.0), @view(conductor_permeability[phase, 1:layer_count])) &&
            all(>(0.0), @view(insulation_permeability[phase, 1:layer_count])) &&
            all(>(0.0), @view(insulation_permittivity[phase, 1:layer_count])) ||
            throw(ArgumentError("active nested cable material values must be positive"))
    end

    internal = _nested_cable_internal_impedance_matrix(
        state,
        radii,
        resistivity,
        conductor_permeability,
        insulation_permeability,
        frequency,
        counts,
        layout,
    )
    include_pipe && _nested_pipe_series_coupling!(
        internal,
        state,
        radii,
        resistivity,
        conductor_permeability,
        insulation_permeability,
        frequency,
        counts,
        layout,
    )
    earth = _nested_cable_earth_impedance_matrix(
        geometry,
        layout.phase_by_conductor,
        frequency,
        earth_resistivity,
    )
    unreduced_series = internal + earth
    active = state.selected_grounded_conductor_count
    1 <= active <= total ||
        throw(ArgumentError("nested cable frequency state requires at least one active conductor"))
    active_indices = 1:active
    grounded_indices = (active + 1):total
    reduced_series = if isempty(grounded_indices)
        Matrix{ComplexF64}(unreduced_series[active_indices, active_indices])
    else
        Matrix{ComplexF64}(
            unreduced_series[active_indices, active_indices] -
            unreduced_series[active_indices, grounded_indices] *
            (unreduced_series[grounded_indices, grounded_indices] \
                unreduced_series[grounded_indices, active_indices]),
        )
    end
    unreduced_shunt = _nested_cable_admittance_matrix(
        state,
        radii,
        insulation_permittivity,
        frequency,
        counts,
        layout,
    )
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
    capacitance_row_sums = vec(sum(capacitance; dims = 2))
    balance_tolerance = max(
        1.0e-18,
        100.0 * eps(Float64) * maximum(abs, capacitance; init = 0.0),
    )
    admittance_balance_passed =
        all(>=(-balance_tolerance), capacitance_row_sums) &&
        sum(capacitance_row_sums) > balance_tolerance
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
        !isempty(grounded_indices),
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

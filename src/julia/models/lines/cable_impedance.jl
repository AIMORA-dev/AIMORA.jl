
function cable_homogeneous_earth_return_impedance_matrix(
    constants::CableGeometryConstants,
    frequency_hz::Real,
    earth_resistivity_ohm_m::Real,
)
    frequency = _checked_line_positive(frequency_hz, "earth-return frequency_hz")
    earth_resistivity = _checked_line_positive(earth_resistivity_ohm_m, "earth_resistivity_ohm_m")
    matrix = Matrix{ComplexF64}(undef, constants.conductor_count, constants.conductor_count)
    for row in 1:constants.conductor_count, col in 1:constants.conductor_count
        matrix[row, col] = cable_homogeneous_earth_return_impedance(
            constants.image_distance_m[row, col],
            constants.angle_rad[row, col],
            earth_resistivity,
            frequency,
        )
    end
    return matrix
end
function cable_bounded_earth_return_impedance_matrix(
    constants::CableGeometryConstants,
    frequency_hz::Real,
    earth_resistivity_ohm_m::Real,
)
    frequency = _checked_line_positive(frequency_hz, "earth-return frequency_hz")
    earth_resistivity = _checked_line_positive(earth_resistivity_ohm_m, "earth_resistivity_ohm_m")
    omega = 2.0 * pi * frequency
    equivalent_depth = sqrt(2.0 * earth_resistivity / (omega * LINE_VACUUM_PERMEABILITY_H_PER_M))
    matrix = Matrix{ComplexF64}(undef, constants.conductor_count, constants.conductor_count)
    for row in 1:constants.conductor_count, col in 1:constants.conductor_count
        dx = constants.horizontal_position_m[row] - constants.horizontal_position_m[col]
        image_distance = constants.image_distance_m[row, col]
        effective_image_distance = hypot(
            dx,
            constants.depth_m[row] + constants.depth_m[col] + 2.0 * equivalent_depth,
        )
        loss = earth_resistivity / (2.0 * pi * effective_image_distance^2)
        reactance = omega * (LINE_VACUUM_PERMEABILITY_H_PER_M / (2.0 * pi)) *
            log(effective_image_distance / image_distance)
        matrix[row, col] = ComplexF64(loss, reactance)
    end
    return matrix
end

function _checked_optional_cable_earth_resistivity(value)
    value === nothing && return nothing
    return _checked_line_positive(value, "earth_resistivity_ohm_m")
end

function _checked_cable_earth_return_model(model::Symbol)
    model in (:bounded, :homogeneous) ||
        throw(ArgumentError("earth_return_model must be :bounded or :homogeneous"))
    return model
end

function _cable_skin_effect_diffusion_factor_value(value)
    inner = nothing
    outer = nothing
    if value isa NamedTuple
        if haskey(value, :inner_diffusion_factor)
            inner = getproperty(value, :inner_diffusion_factor)
        elseif haskey(value, :inner)
            inner = getproperty(value, :inner)
        end
        if haskey(value, :outer_diffusion_factor)
            outer = getproperty(value, :outer_diffusion_factor)
        elseif haskey(value, :outer)
            outer = getproperty(value, :outer)
        end
    elseif value isa Pair
        inner = first(value)
        outer = last(value)
    elseif value isa Tuple || value isa AbstractVector
        length(value) == 2 ||
            throw(ArgumentError("skin-effect diffusion factor entries must have two values"))
        inner = value[1]
        outer = value[2]
    end
    (inner === nothing || outer === nothing) &&
        throw(ArgumentError("skin-effect diffusion factor requires inner and outer values"))
    inner_value = _checked_line_nonnegative(inner, "skin-effect inner diffusion factor")
    outer_value = _checked_line_positive(outer, "skin-effect outer diffusion factor")
    inner_value <= outer_value ||
        throw(ArgumentError("skin-effect inner diffusion factor must not exceed the outer factor"))
    return (
        inner_diffusion_factor = inner_value,
        outer_diffusion_factor = outer_value,
    )
end

function _cable_skin_effect_diffusion_factor_pair(
    skin_effect_diffusion_factors,
    conductor_index::Int,
    conductor_count::Int,
)
    skin_effect_diffusion_factors === nothing && return nothing
    factors = skin_effect_diffusion_factors
    if factors isa AbstractMatrix
        if size(factors, 1) == conductor_count && size(factors, 2) >= 2
            return _cable_skin_effect_diffusion_factor_value((
                factors[conductor_index, 1],
                factors[conductor_index, 2],
            ))
        elseif size(factors, 2) == conductor_count && size(factors, 1) >= 2
            return _cable_skin_effect_diffusion_factor_value((
                factors[1, conductor_index],
                factors[2, conductor_index],
            ))
        end
        throw(ArgumentError("skin-effect diffusion factor matrix must have two values per conductor"))
    elseif factors isa NamedTuple &&
           haskey(factors, :inner_diffusion_factors) &&
           haskey(factors, :outer_diffusion_factors)
        length(factors.inner_diffusion_factors) == conductor_count &&
            length(factors.outer_diffusion_factors) == conductor_count ||
            throw(ArgumentError("skin-effect diffusion factor vectors must match conductor count"))
        return _cable_skin_effect_diffusion_factor_value((
            factors.inner_diffusion_factors[conductor_index],
            factors.outer_diffusion_factors[conductor_index],
        ))
    elseif factors isa AbstractVector || factors isa Tuple
        length(factors) == conductor_count ||
            throw(ArgumentError("skin-effect diffusion factor count must match conductor count"))
        return _cable_skin_effect_diffusion_factor_value(factors[conductor_index])
    end
    throw(ArgumentError("skin-effect diffusion factors must be a vector, tuple, matrix, or named tuple"))
end

function _cable_active_conductor_count(constants::CableGeometryConstants)
    return constants.conductor_count - constants.grounded_conductor_count
end

function _cable_grounded_reduced_series_impedance(
    conductor_series::AbstractMatrix{ComplexF64},
    constants::CableGeometryConstants,
    reduce_grounded_conductors::Bool,
)
    active_count = _cable_active_conductor_count(constants)
    active_count == sum(constants.phase_conductor_counts) ||
        throw(ArgumentError("active cable conductor count must match phase conductor counts"))
    active_series = Matrix{ComplexF64}(conductor_series[1:active_count, 1:active_count])
    grounded = constants.grounded_conductor_count
    executed = reduce_grounded_conductors && grounded > 0
    if executed
        grounded_indices = (active_count + 1):constants.conductor_count
        zaa = conductor_series[1:active_count, 1:active_count]
        zag = conductor_series[1:active_count, grounded_indices]
        zgg = conductor_series[grounded_indices, grounded_indices]
        zga = conductor_series[grounded_indices, 1:active_count]
        active_series .= zaa .- zag * (zgg \ zga)
    end
    return active_series, executed
end

function cable_phase_electrostatic_admittance(
    constants::CableGeometryConstants,
    frequency_hz::Real;
    relative_permittivity = nothing,
)
    conductor_admittance = cable_geometry_electrostatic_admittance(
        constants,
        frequency_hz;
        relative_permittivity = relative_permittivity,
    )
    incidence = _cable_phase_conductor_matrix(constants)
    phase_capacitance = incidence *
        conductor_admittance.capacitance_matrix_f_per_m *
        transpose(incidence)
    phase_shunt = incidence *
        conductor_admittance.shunt_admittance_matrix_s_per_m *
        transpose(incidence)
    all(isfinite, phase_capacitance) ||
        throw(ArgumentError("cable phase capacitance matrix entries must be finite"))
    all(value -> isfinite(real(value)) && isfinite(imag(value)), phase_shunt) ||
        throw(ArgumentError("cable phase shunt-admittance matrix entries must be finite"))
    return CablePhaseElectrostaticAdmittance(
        conductor_admittance.frequency_hz,
        conductor_admittance.angular_frequency_rad_s,
        conductor_admittance.relative_permittivity,
        incidence,
        Matrix{Float64}(phase_capacitance),
        Matrix{ComplexF64}(phase_shunt),
        conductor_admittance,
        constants.conductor_count,
        constants.phase_count,
    )
end

function cable_phase_series_impedance(
    constants::CableGeometryConstants,
    frequency_hz::Real;
    include_external_inductance::Bool = true,
    include_internal_inductance::Bool = true,
    include_bounded_skin_effect::Bool = false,
    skin_effect_diffusion_factors = nothing,
    earth_resistivity_ohm_m = nothing,
    earth_return_model::Symbol = :bounded,
    reduce_grounded_conductors::Bool = true,
)
    frequency = Float64(frequency_hz)
    isfinite(frequency) && frequency > 0.0 ||
        throw(ArgumentError("cable series impedance frequency_hz must be finite and positive"))
    earth_resistivity = _checked_optional_cable_earth_resistivity(earth_resistivity_ohm_m)
    earth_model = _checked_cable_earth_return_model(earth_return_model)
    omega = 2.0 * pi * frequency
    conductor_count = constants.conductor_count
    external = zeros(ComplexF64, conductor_count, conductor_count)
    if include_external_inductance
        external .= ComplexF64.(
            0.0,
            omega * (LINE_VACUUM_PERMEABILITY_H_PER_M / (2.0 * pi)) .*
                constants.potential_log_matrix,
        )
    end
    earth_return = earth_resistivity === nothing ?
        zeros(ComplexF64, conductor_count, conductor_count) :
        earth_model == :homogeneous ?
        cable_homogeneous_earth_return_impedance_matrix(constants, frequency, earth_resistivity) :
        cable_bounded_earth_return_impedance_matrix(constants, frequency, earth_resistivity)
    earth_return_impedance_executed =
        earth_resistivity !== nothing && earth_model == :homogeneous
    internal = ComplexF64[]
    skin_delta = ComplexF64[]
    sizehint!(internal, conductor_count)
    sizehint!(skin_delta, conductor_count)
    skin_effect_internal_impedance_executed = false
    for idx in 1:conductor_count
        area_m2 = pi * constants.radius_m[idx]^2
        dc_impedance = ComplexF64(
            constants.resistivity_ohm_m[idx] / area_m2,
            include_internal_inductance ?
            omega * constants.permeability_h_per_m[idx] / (8.0 * pi) :
            0.0,
        )
        skin_factor = _cable_skin_effect_diffusion_factor_pair(
            skin_effect_diffusion_factors,
            idx,
            conductor_count,
        )
        skin_impedance =
            include_bounded_skin_effect && skin_factor !== nothing ?
            cable_skin_effect_internal_impedance(
                skin_factor.inner_diffusion_factor,
                skin_factor.outer_diffusion_factor,
                constants.relative_permeability[idx],
                frequency,
            ) :
            include_bounded_skin_effect ?
            cable_bounded_skin_effect_internal_impedance(
                constants.radius_m[idx],
                constants.resistivity_ohm_m[idx],
                constants.permeability_h_per_m[idx],
                frequency,
            ) :
            dc_impedance
        skin_effect_internal_impedance_executed |=
            include_bounded_skin_effect && skin_factor !== nothing
        if include_bounded_skin_effect && !include_internal_inductance
            skin_impedance = ComplexF64(real(skin_impedance), 0.0)
        end
        push!(internal, skin_impedance)
        push!(skin_delta, skin_impedance - dc_impedance)
    end
    conductor_series = external .+ earth_return
    for idx in 1:conductor_count
        conductor_series[idx, idx] += internal[idx]
    end
    average = _cable_phase_average_matrix(constants)
    active_series, grounded_reduction_executed =
        _cable_grounded_reduced_series_impedance(conductor_series, constants, reduce_grounded_conductors)
    active_count = _cable_active_conductor_count(constants)
    active_average = average[:, 1:active_count]
    phase_unreduced = average * conductor_series * transpose(average)
    phase_series = active_average * active_series * transpose(active_average)
    all(value -> isfinite(real(value)) && isfinite(imag(value)), conductor_series) ||
        throw(ArgumentError("cable conductor series impedance entries must be finite"))
    all(value -> isfinite(real(value)) && isfinite(imag(value)), active_series) ||
        throw(ArgumentError("cable active conductor series impedance entries must be finite"))
    all(value -> isfinite(real(value)) && isfinite(imag(value)), phase_series) ||
        throw(ArgumentError("cable phase series impedance entries must be finite"))
    return CableConductorSeriesImpedance(
        frequency,
        omega,
        external,
        internal,
        skin_delta,
        Matrix{ComplexF64}(earth_return),
        Matrix{ComplexF64}(conductor_series),
        Matrix{ComplexF64}(active_series),
        average,
        Matrix{ComplexF64}(phase_unreduced),
        Matrix{ComplexF64}(phase_series),
        include_bounded_skin_effect,
        skin_effect_internal_impedance_executed,
        earth_resistivity !== nothing,
        earth_return_impedance_executed,
        grounded_reduction_executed,
        earth_resistivity === nothing ? NaN : earth_resistivity,
        constants.grounded_conductor_count,
        conductor_count,
        constants.phase_count,
    )
end

function cable_frequency_scan_series_impedance(
    constants::CableGeometryConstants,
    frequencies_hz::AbstractVector;
    include_external_inductance::Bool = true,
    include_internal_inductance::Bool = true,
    include_bounded_skin_effect::Bool = false,
    skin_effect_diffusion_factors = nothing,
    earth_resistivity_ohm_m = nothing,
    earth_return_model::Symbol = :bounded,
    reduce_grounded_conductors::Bool = true,
)
    isempty(frequencies_hz) &&
        throw(ArgumentError("cable frequency scan series impedance must contain at least one frequency"))
    frequencies = Float64.(frequencies_hz)
    all(frequency -> isfinite(frequency) && frequency > 0.0, frequencies) ||
        throw(ArgumentError("cable frequency scan series impedance frequencies must be finite and positive"))
    sorted_order = sortperm(frequencies)
    sorted_frequencies = frequencies[sorted_order]
    for idx in 2:length(sorted_frequencies)
        sorted_frequencies[idx] != sorted_frequencies[idx - 1] ||
            throw(ArgumentError("cable frequency scan series impedance frequencies must be unique"))
    end
    rows = CableConductorSeriesImpedance[]
    matrices = Matrix{ComplexF64}[]
    sizehint!(rows, length(sorted_frequencies))
    sizehint!(matrices, length(sorted_frequencies))
    for source_index in sorted_order
        row = cable_phase_series_impedance(
            constants,
            frequencies[source_index];
            include_external_inductance = include_external_inductance,
            include_internal_inductance = include_internal_inductance,
            include_bounded_skin_effect = include_bounded_skin_effect,
            skin_effect_diffusion_factors = skin_effect_diffusion_factors,
            earth_resistivity_ohm_m = earth_resistivity_ohm_m,
            earth_return_model = earth_return_model,
            reduce_grounded_conductors = reduce_grounded_conductors,
        )
        push!(rows, row)
        push!(matrices, copy(row.phase_series_impedance_matrix_ohm_per_m))
    end
    return CableFrequencyScanSeriesImpedance(
        sorted_frequencies,
        2.0 .* pi .* sorted_frequencies,
        rows,
        matrices,
        collect(sorted_order),
        length(sorted_frequencies),
        constants.phase_count,
    )
end

function _line_offdiagonal_max_abs(matrix::AbstractMatrix{<:Complex})
    rows, cols = size(matrix)
    rows == cols || throw(ArgumentError("line modal matrix must be square"))
    rows <= 1 && return 0.0
    max_value = 0.0
    for col in 1:cols, row in 1:rows
        row == col && continue
        max_value = max(max_value, abs(matrix[row, col]))
    end
    return max_value
end

function cable_phase_line_constants(
    constants::CableGeometryConstants,
    frequency_hz::Real,
    phase_series_impedance_matrix::AbstractMatrix;
    line_length_m::Real = 0.0,
    relative_permittivity = nothing,
    unwind_state::Union{Nothing,LineModeUnwindState} = nothing,
    ntol::Integer = 1,
    nrp::Integer = 0,
)
    frequency = Float64(frequency_hz)
    isfinite(frequency) && frequency > 0.0 ||
        throw(ArgumentError("frequency_hz must be finite and positive"))
    line_length = _checked_line_length(line_length_m)
    phase_admittance = cable_phase_electrostatic_admittance(
        constants,
        frequency;
        relative_permittivity = relative_permittivity,
    )
    series = _checked_line_complex_square_matrix(
        phase_series_impedance_matrix,
        "cable phase series impedance matrix",
    )
    size(series, 1) == phase_admittance.phase_count ||
        throw(ArgumentError("cable phase series impedance matrix dimension must match phase count"))
    shunt = phase_admittance.phase_shunt_admittance_matrix_s_per_m
    phase_zy = series * shunt
    solution = line_modal_solution(
        phase_zy,
        frequency;
        unwind_state = unwind_state,
        ntol = ntol,
        nrp = nrp,
    )
    modal_series_matrix = solution.transform.phase_to_modal *
        series *
        solution.transform.modal_to_phase
    modal_shunt_matrix = solution.transform.phase_to_modal *
        shunt *
        solution.transform.modal_to_phase
    mode_count = phase_admittance.phase_count
    modal_series = ComplexF64[modal_series_matrix[idx, idx] for idx in 1:mode_count]
    modal_shunt = ComplexF64[modal_shunt_matrix[idx, idx] for idx in 1:mode_count]
    points = LineFrequencyPoint[]
    sizehint!(points, mode_count)
    root_errors = Float64[]
    for idx in 1:mode_count
        modal_shunt[idx] != 0.0 + 0.0im ||
            throw(ArgumentError("cable modal shunt admittance must be nonzero"))
        characteristic_impedance = sqrt(modal_series[idx] / modal_shunt[idx])
        propagation_constant = solution.propagation_roots[idx]
        push!(
            points,
            LineFrequencyPoint(
                frequency,
                characteristic_impedance,
                propagation_constant,
                exp(-propagation_constant * line_length),
            ),
        )
        push!(
            root_errors,
            abs(propagation_constant - sqrt(modal_series[idx] * modal_shunt[idx])),
        )
    end
    return CablePhaseLineConstants(
        frequency,
        2.0 * pi * frequency,
        line_length,
        Matrix{ComplexF64}(series),
        Matrix{ComplexF64}(shunt),
        Matrix{ComplexF64}(phase_zy),
        solution,
        modal_series,
        modal_shunt,
        points,
        max(
            _line_offdiagonal_max_abs(modal_series_matrix),
            _line_offdiagonal_max_abs(modal_shunt_matrix),
        ),
        maximum(root_errors),
        phase_admittance,
        mode_count,
    )
end

function cable_frequency_scan_line_constants(
    constants::CableGeometryConstants,
    frequencies_hz::AbstractVector,
    phase_series_impedance_matrices::AbstractVector;
    line_length_m::Real = 0.0,
    relative_permittivity = nothing,
    skin_effect_internal_impedance_executed::Bool = false,
    earth_return_impedance_executed::Bool = false,
    grounded_conductor_reduction_executed::Bool = false,
    ntol::Integer = 2,
    nrp::Integer = 0,
)
    isempty(frequencies_hz) &&
        throw(ArgumentError("cable frequency scan must contain at least one frequency"))
    length(phase_series_impedance_matrices) == length(frequencies_hz) ||
        throw(ArgumentError("cable frequency scan matrix count must match frequency count"))
    frequencies = Float64.(frequencies_hz)
    all(frequency -> isfinite(frequency) && frequency > 0.0, frequencies) ||
        throw(ArgumentError("cable frequency scan frequencies must be finite and positive"))
    sorted_order = sortperm(frequencies)
    sorted_frequencies = frequencies[sorted_order]
    for idx in 2:length(sorted_frequencies)
        sorted_frequencies[idx] != sorted_frequencies[idx - 1] ||
            throw(ArgumentError("cable frequency scan frequencies must be unique"))
    end

    line_length = _checked_line_length(line_length_m)
    mode_state = LineModeUnwindState(constants.phase_count)
    frequency_constants = CablePhaseLineConstants[]
    sizehint!(frequency_constants, length(sorted_frequencies))
    for source_index in sorted_order
        push!(
            frequency_constants,
            cable_phase_line_constants(
                constants,
                frequencies[source_index],
                phase_series_impedance_matrices[source_index];
                line_length_m = line_length,
                relative_permittivity = relative_permittivity,
                unwind_state = mode_state,
                ntol = ntol,
                nrp = nrp,
            ),
        )
    end
    sample_rows = [copy(row.frequency_points) for row in frequency_constants]
    return CableFrequencyScanLineConstants(
        sorted_frequencies,
        2.0 .* pi .* sorted_frequencies,
        line_length,
        frequency_constants,
        sample_rows,
        collect(sorted_order),
        length(sorted_frequencies),
        constants.phase_count,
        skin_effect_internal_impedance_executed,
        earth_return_impedance_executed,
        grounded_conductor_reduction_executed,
    )
end

function cable_frequency_scan_line_constants_from_geometry(
    constants::CableGeometryConstants,
    frequencies_hz::AbstractVector;
    line_length_m::Real = 0.0,
    relative_permittivity = nothing,
    include_external_inductance::Bool = true,
    include_internal_inductance::Bool = true,
    include_bounded_skin_effect::Bool = false,
    skin_effect_diffusion_factors = nothing,
    earth_resistivity_ohm_m = nothing,
    earth_return_model::Symbol = :bounded,
    reduce_grounded_conductors::Bool = true,
    ntol::Integer = 2,
    nrp::Integer = 0,
)
    series = cable_frequency_scan_series_impedance(
        constants,
        frequencies_hz;
        include_external_inductance = include_external_inductance,
        include_internal_inductance = include_internal_inductance,
        include_bounded_skin_effect = include_bounded_skin_effect,
        skin_effect_diffusion_factors = skin_effect_diffusion_factors,
        earth_resistivity_ohm_m = earth_resistivity_ohm_m,
        earth_return_model = earth_return_model,
        reduce_grounded_conductors = reduce_grounded_conductors,
    )
    skin_effect_internal_impedance_executed =
        any(row -> row.skin_effect_internal_impedance_executed, series.frequency_rows)
    earth_return_impedance_executed =
        any(row -> row.earth_return_impedance_executed, series.frequency_rows)
    grounded_conductor_reduction_executed =
        any(row -> row.grounded_conductor_reduction_executed, series.frequency_rows)
    line_constants = cable_frequency_scan_line_constants(
        constants,
        series.frequencies_hz,
        series.phase_series_impedance_matrices_ohm_per_m;
        line_length_m = line_length_m,
        relative_permittivity = relative_permittivity,
        skin_effect_internal_impedance_executed = skin_effect_internal_impedance_executed,
        earth_return_impedance_executed = earth_return_impedance_executed,
        grounded_conductor_reduction_executed = grounded_conductor_reduction_executed,
        ntol = ntol,
        nrp = nrp,
    )
    return CableFrequencyScanGeneratedLineConstants(series, line_constants)
end

function _cable_frequency_scan_schedule_earth_resistivity(
    schedule::CableFrequencyScanLoopSchedule,
    supplied_value,
)
    supplied_value !== nothing && return supplied_value
    return schedule.final_earth_resistivity_ohm_m > 0.0 ?
        schedule.final_earth_resistivity_ohm_m :
        nothing
end

function _cable_pipe_sheath_side_effects_executed(state)
    return state isa CablePipeSheathDerivedState && state.derived_state_executed
end

function cable_frequency_scan_line_constants_from_geometry(
    constants::CableGeometryConstants,
    schedule::CableFrequencyScanLoopSchedule;
    line_length_m::Real = 0.0,
    relative_permittivity = nothing,
    include_external_inductance::Bool = true,
    include_internal_inductance::Bool = true,
    include_bounded_skin_effect::Bool = false,
    skin_effect_diffusion_factors = nothing,
    earth_resistivity_ohm_m = nothing,
    earth_return_model::Symbol = :bounded,
    reduce_grounded_conductors::Bool = true,
    ntol::Integer = 2,
    nrp::Integer = 0,
)
    schedule.loop_executed ||
        throw(ArgumentError("cable frequency scan schedule must be executed"))
    return cable_frequency_scan_line_constants_from_geometry(
        constants,
        schedule.frequencies_hz;
        line_length_m = line_length_m,
        relative_permittivity = relative_permittivity,
        include_external_inductance = include_external_inductance,
        include_internal_inductance = include_internal_inductance,
        include_bounded_skin_effect = include_bounded_skin_effect,
        skin_effect_diffusion_factors = skin_effect_diffusion_factors,
        earth_resistivity_ohm_m =
            _cable_frequency_scan_schedule_earth_resistivity(schedule, earth_resistivity_ohm_m),
        earth_return_model = earth_return_model,
        reduce_grounded_conductors = reduce_grounded_conductors,
        ntol = ntol,
        nrp = nrp,
    )
end

function cable_frequency_scan_runtime_update_from_geometry(
    constants::CableGeometryConstants,
    frequencies_hz::AbstractVector,
    target_frequency_hz::Real,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector;
    line_length_m::Real = 0.0,
    relative_permittivity = nothing,
    include_external_inductance::Bool = true,
    include_internal_inductance::Bool = true,
    include_bounded_skin_effect::Bool = false,
    skin_effect_diffusion_factors = nothing,
    earth_resistivity_ohm_m = nothing,
    earth_return_model::Symbol = :bounded,
    reduce_grounded_conductors::Bool = true,
    pipe_sheath_state = nothing,
    initial_frequency_hz::Real = target_frequency_hz,
    initial_ntol::Integer = 1,
    initial_nrp::Integer = 0,
    ntol::Integer = 2,
    nrp::Integer = 0,
)
    generated = cable_frequency_scan_line_constants_from_geometry(
        constants,
        frequencies_hz;
        line_length_m = line_length_m,
        relative_permittivity = relative_permittivity,
        include_external_inductance = include_external_inductance,
        include_internal_inductance = include_internal_inductance,
        include_bounded_skin_effect = include_bounded_skin_effect,
        skin_effect_diffusion_factors = skin_effect_diffusion_factors,
        earth_resistivity_ohm_m = earth_resistivity_ohm_m,
        earth_return_model = earth_return_model,
        reduce_grounded_conductors = reduce_grounded_conductors,
        ntol = ntol,
        nrp = nrp,
    )
    result = cable_frequency_scan_runtime_update(
        generated.line_constants,
        target_frequency_hz,
        from_phase_voltage,
        to_phase_voltage;
        initial_frequency_hz = initial_frequency_hz,
        initial_ntol = initial_ntol,
        initial_nrp = initial_nrp,
        ntol = ntol,
        nrp = nrp,
    )
    skin_effect_internal_impedance_executed =
        generated.line_constants.skin_effect_internal_impedance_executed
    earth_return_impedance_executed =
        generated.line_constants.earth_return_impedance_executed
    pipe_sheath_side_effects_executed =
        _cable_pipe_sheath_side_effects_executed(pipe_sheath_state)
    runtime_deferred_calls = Symbol[]
    skin_effect_internal_impedance_executed ||
        push!(runtime_deferred_calls, :full_over47_skin_effect)
    earth_return_impedance_executed ||
        push!(runtime_deferred_calls, :full_over47_earth_return)
    append!(
        runtime_deferred_calls,
        (
            (pipe_sheath_side_effects_executed ? () : (:full_over47_pipe_or_sheath_coupling,))...,
            :full_over47_frequency_loop_side_effects,
            :full_dceign_lr_modal_eigenvector_bulk_package,
            :full_bpa_frequency_dependent_fitting,
            :recursive_convolution_line_runtime_oracle,
            :line_timestep_bulk_oracle,
        ),
    )
    return merge(
        result,
        (
            source = :cable_frequency_scan_runtime_update_from_geometry,
            cable_series_impedance = generated.series_impedance,
            bounded_cable_phase_series_impedance_generation_executed = true,
            bounded_cable_frequency_scan_series_impedance_generation_executed = true,
            bounded_over47_skin_effect_executed =
                any(row -> row.bounded_skin_effect_executed, generated.series_impedance.frequency_rows),
            skin_effect_internal_impedance_executed =
                skin_effect_internal_impedance_executed,
            bounded_over47_earth_return_executed =
                any(row -> row.bounded_earth_return_executed, generated.series_impedance.frequency_rows),
            earth_return_impedance_executed =
                earth_return_impedance_executed,
            bounded_over47_grounded_conductor_reduction_executed =
                any(row -> row.grounded_conductor_reduction_executed, generated.series_impedance.frequency_rows),
            full_over47_skin_effect_executed =
                skin_effect_internal_impedance_executed,
            full_over47_earth_return_executed =
                earth_return_impedance_executed,
            cable_pipe_sheath_state = pipe_sheath_state,
            pipe_sheath_side_effects_executed =
                pipe_sheath_side_effects_executed,
            full_over47_pipe_or_sheath_coupling_executed =
                pipe_sheath_side_effects_executed,
            deferred_calls = Tuple(runtime_deferred_calls),
        ),
    )
end

function cable_frequency_scan_runtime_update_from_geometry(
    constants::CableGeometryConstants,
    schedule::CableFrequencyScanLoopSchedule,
    target_frequency_hz::Real,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector;
    line_length_m::Real = 0.0,
    relative_permittivity = nothing,
    include_external_inductance::Bool = true,
    include_internal_inductance::Bool = true,
    include_bounded_skin_effect::Bool = false,
    skin_effect_diffusion_factors = nothing,
    earth_resistivity_ohm_m = nothing,
    earth_return_model::Symbol = :bounded,
    reduce_grounded_conductors::Bool = true,
    pipe_sheath_state = nothing,
    initial_frequency_hz::Real = target_frequency_hz,
    initial_ntol::Integer = 1,
    initial_nrp::Integer = 0,
    ntol::Integer = 2,
    nrp::Integer = 0,
)
    generated = cable_frequency_scan_line_constants_from_geometry(
        constants,
        schedule;
        line_length_m = line_length_m,
        relative_permittivity = relative_permittivity,
        include_external_inductance = include_external_inductance,
        include_internal_inductance = include_internal_inductance,
        include_bounded_skin_effect = include_bounded_skin_effect,
        skin_effect_diffusion_factors = skin_effect_diffusion_factors,
        earth_resistivity_ohm_m = earth_resistivity_ohm_m,
        earth_return_model = earth_return_model,
        reduce_grounded_conductors = reduce_grounded_conductors,
        ntol = ntol,
        nrp = nrp,
    )
    result = cable_frequency_scan_runtime_update(
        generated.line_constants,
        target_frequency_hz,
        from_phase_voltage,
        to_phase_voltage;
        initial_frequency_hz = initial_frequency_hz,
        initial_ntol = initial_ntol,
        initial_nrp = initial_nrp,
        ntol = ntol,
        nrp = nrp,
    )
    skin_effect_internal_impedance_executed =
        generated.line_constants.skin_effect_internal_impedance_executed
    earth_return_impedance_executed =
        generated.line_constants.earth_return_impedance_executed
    pipe_sheath_side_effects_executed =
        _cable_pipe_sheath_side_effects_executed(pipe_sheath_state)
    runtime_deferred_calls = Symbol[]
    skin_effect_internal_impedance_executed ||
        push!(runtime_deferred_calls, :full_over47_skin_effect)
    earth_return_impedance_executed ||
        push!(runtime_deferred_calls, :full_over47_earth_return)
    append!(
        runtime_deferred_calls,
        (
            (pipe_sheath_side_effects_executed ? () : (:full_over47_pipe_or_sheath_coupling,))...,
            :full_dceign_lr_modal_eigenvector_bulk_package,
            :full_bpa_frequency_dependent_fitting,
            :recursive_convolution_line_runtime_oracle,
            :line_timestep_bulk_oracle,
        ),
    )
    return merge(
        result,
        (
            source = :cable_frequency_scan_runtime_update_from_geometry,
            cable_frequency_scan_loop_schedule = schedule,
            cable_frequency_scan_loop_frequencies_hz = copy(schedule.frequencies_hz),
            cable_frequency_scan_loop_output_row_count = schedule.output_row_count,
            cable_frequency_scan_loop_retry_requested = schedule.retry_requested,
            cable_frequency_scan_loop_unit_distance_record_count =
                schedule.unit_distance_record_count,
            cable_series_impedance = generated.series_impedance,
            bounded_cable_phase_series_impedance_generation_executed = true,
            bounded_cable_frequency_scan_series_impedance_generation_executed = true,
            bounded_over47_skin_effect_executed =
                any(row -> row.bounded_skin_effect_executed, generated.series_impedance.frequency_rows),
            skin_effect_internal_impedance_executed =
                skin_effect_internal_impedance_executed,
            bounded_over47_earth_return_executed =
                any(row -> row.bounded_earth_return_executed, generated.series_impedance.frequency_rows),
            earth_return_impedance_executed =
                earth_return_impedance_executed,
            bounded_over47_grounded_conductor_reduction_executed =
                any(row -> row.grounded_conductor_reduction_executed, generated.series_impedance.frequency_rows),
            full_over47_skin_effect_executed =
                skin_effect_internal_impedance_executed,
            full_over47_earth_return_executed =
                earth_return_impedance_executed,
            cable_pipe_sheath_state = pipe_sheath_state,
            pipe_sheath_side_effects_executed =
                pipe_sheath_side_effects_executed,
            full_over47_pipe_or_sheath_coupling_executed =
                pipe_sheath_side_effects_executed,
            frequency_loop_executed = true,
            full_over47_frequency_loop_executed = true,
            full_over47_frequency_loop_side_effects_executed = true,
            deferred_calls = Tuple(runtime_deferred_calls),
        ),
    )
end

function cable_frequency_scan_recursive_convolution_update_from_geometry(
    constants::CableGeometryConstants,
    frequencies_hz::AbstractVector,
    pole_decay::AbstractMatrix,
    dt_s::Real,
    target_frequency_hz::Real,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector;
    line_length_m::Real = 0.0,
    relative_permittivity = nothing,
    include_external_inductance::Bool = true,
    include_internal_inductance::Bool = true,
    include_bounded_skin_effect::Bool = false,
    skin_effect_diffusion_factors = nothing,
    earth_resistivity_ohm_m = nothing,
    earth_return_model::Symbol = :bounded,
    reduce_grounded_conductors::Bool = true,
    pipe_sheath_state = nothing,
    response_kind::Symbol = :propagated_modal_admittance,
    initial_frequency_hz::Real = target_frequency_hz,
    initial_ntol::Integer = 1,
    initial_nrp::Integer = 0,
    ntol::Integer = 2,
    nrp::Integer = 0,
)
    generated = cable_frequency_scan_line_constants_from_geometry(
        constants,
        frequencies_hz;
        line_length_m = line_length_m,
        relative_permittivity = relative_permittivity,
        include_external_inductance = include_external_inductance,
        include_internal_inductance = include_internal_inductance,
        include_bounded_skin_effect = include_bounded_skin_effect,
        skin_effect_diffusion_factors = skin_effect_diffusion_factors,
        earth_resistivity_ohm_m = earth_resistivity_ohm_m,
        earth_return_model = earth_return_model,
        reduce_grounded_conductors = reduce_grounded_conductors,
        ntol = ntol,
        nrp = nrp,
    )
    pipe_sheath_side_effects_executed =
        _cable_pipe_sheath_side_effects_executed(pipe_sheath_state)
    result = cable_frequency_scan_recursive_convolution_update_from_scan_fit(
        generated.line_constants,
        pole_decay,
        dt_s,
        target_frequency_hz,
        from_phase_voltage,
        to_phase_voltage;
        response_kind = response_kind,
        initial_frequency_hz = initial_frequency_hz,
        initial_ntol = initial_ntol,
        initial_nrp = initial_nrp,
        ntol = ntol,
        nrp = nrp,
        pipe_sheath_side_effects_executed =
            pipe_sheath_side_effects_executed,
    )
    skin_effect_internal_impedance_executed =
        generated.line_constants.skin_effect_internal_impedance_executed
    earth_return_impedance_executed =
        generated.line_constants.earth_return_impedance_executed
    runtime_deferred_calls = Symbol[]
    skin_effect_internal_impedance_executed ||
        push!(runtime_deferred_calls, :full_over47_skin_effect)
    earth_return_impedance_executed ||
        push!(runtime_deferred_calls, :full_over47_earth_return)
    append!(
        runtime_deferred_calls,
        (
            (pipe_sheath_side_effects_executed ? () : (:full_over47_pipe_or_sheath_coupling,))...,
            :full_over47_frequency_loop_side_effects,
            :full_bpa_frequency_dependent_fitting,
            :pole_zero_fit_oracle,
            :recursive_convolution_line_runtime_oracle,
            :line_timestep_bulk_oracle,
        ),
    )
    return merge(
        result,
        (
            source = :cable_frequency_scan_recursive_convolution_update_from_geometry,
            cable_series_impedance = generated.series_impedance,
            bounded_cable_phase_series_impedance_generation_executed = true,
            bounded_cable_frequency_scan_series_impedance_generation_executed = true,
            bounded_cable_frequency_scan_response_fit_executed = true,
            bounded_cable_frequency_scan_recursive_convolution_runtime_executed = true,
            bounded_over47_skin_effect_executed =
                any(row -> row.bounded_skin_effect_executed, generated.series_impedance.frequency_rows),
            skin_effect_internal_impedance_executed =
                skin_effect_internal_impedance_executed,
            bounded_over47_earth_return_executed =
                any(row -> row.bounded_earth_return_executed, generated.series_impedance.frequency_rows),
            earth_return_impedance_executed =
                earth_return_impedance_executed,
            bounded_over47_grounded_conductor_reduction_executed =
                any(row -> row.grounded_conductor_reduction_executed, generated.series_impedance.frequency_rows),
            full_over47_skin_effect_executed =
                skin_effect_internal_impedance_executed,
            full_over47_earth_return_executed =
                earth_return_impedance_executed,
            cable_pipe_sheath_state = pipe_sheath_state,
            pipe_sheath_side_effects_executed =
                pipe_sheath_side_effects_executed,
            full_over47_pipe_or_sheath_coupling_executed =
                pipe_sheath_side_effects_executed,
            deferred_calls = Tuple(runtime_deferred_calls),
        ),
    )
end

function cable_frequency_scan_recursive_convolution_update_from_geometry(
    constants::CableGeometryConstants,
    schedule::CableFrequencyScanLoopSchedule,
    pole_decay::AbstractMatrix,
    dt_s::Real,
    target_frequency_hz::Real,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector;
    line_length_m::Real = 0.0,
    relative_permittivity = nothing,
    include_external_inductance::Bool = true,
    include_internal_inductance::Bool = true,
    include_bounded_skin_effect::Bool = false,
    skin_effect_diffusion_factors = nothing,
    earth_resistivity_ohm_m = nothing,
    earth_return_model::Symbol = :bounded,
    reduce_grounded_conductors::Bool = true,
    pipe_sheath_state = nothing,
    response_kind::Symbol = :propagated_modal_admittance,
    initial_frequency_hz::Real = target_frequency_hz,
    initial_ntol::Integer = 1,
    initial_nrp::Integer = 0,
    ntol::Integer = 2,
    nrp::Integer = 0,
)
    generated = cable_frequency_scan_line_constants_from_geometry(
        constants,
        schedule;
        line_length_m = line_length_m,
        relative_permittivity = relative_permittivity,
        include_external_inductance = include_external_inductance,
        include_internal_inductance = include_internal_inductance,
        include_bounded_skin_effect = include_bounded_skin_effect,
        skin_effect_diffusion_factors = skin_effect_diffusion_factors,
        earth_resistivity_ohm_m = earth_resistivity_ohm_m,
        earth_return_model = earth_return_model,
        reduce_grounded_conductors = reduce_grounded_conductors,
        ntol = ntol,
        nrp = nrp,
    )
    pipe_sheath_side_effects_executed =
        _cable_pipe_sheath_side_effects_executed(pipe_sheath_state)
    result = cable_frequency_scan_recursive_convolution_update_from_scan_fit(
        generated.line_constants,
        pole_decay,
        dt_s,
        target_frequency_hz,
        from_phase_voltage,
        to_phase_voltage;
        response_kind = response_kind,
        initial_frequency_hz = initial_frequency_hz,
        initial_ntol = initial_ntol,
        initial_nrp = initial_nrp,
        ntol = ntol,
        nrp = nrp,
        frequency_loop_executed = true,
        pipe_sheath_side_effects_executed =
            pipe_sheath_side_effects_executed,
    )
    skin_effect_internal_impedance_executed =
        generated.line_constants.skin_effect_internal_impedance_executed
    earth_return_impedance_executed =
        generated.line_constants.earth_return_impedance_executed
    runtime_deferred_calls = Symbol[]
    skin_effect_internal_impedance_executed ||
        push!(runtime_deferred_calls, :full_over47_skin_effect)
    earth_return_impedance_executed ||
        push!(runtime_deferred_calls, :full_over47_earth_return)
    append!(
        runtime_deferred_calls,
        (
            (pipe_sheath_side_effects_executed ? () : (:full_over47_pipe_or_sheath_coupling,))...,
            :full_bpa_frequency_dependent_fitting,
            :pole_zero_fit_oracle,
            :recursive_convolution_line_runtime_oracle,
            :line_timestep_bulk_oracle,
        ),
    )
    return merge(
        result,
        (
            source = :cable_frequency_scan_recursive_convolution_update_from_geometry,
            cable_frequency_scan_loop_schedule = schedule,
            cable_frequency_scan_loop_frequencies_hz = copy(schedule.frequencies_hz),
            cable_frequency_scan_loop_output_row_count = schedule.output_row_count,
            cable_frequency_scan_loop_retry_requested = schedule.retry_requested,
            cable_frequency_scan_loop_unit_distance_record_count =
                schedule.unit_distance_record_count,
            cable_series_impedance = generated.series_impedance,
            bounded_cable_phase_series_impedance_generation_executed = true,
            bounded_cable_frequency_scan_series_impedance_generation_executed = true,
            bounded_cable_frequency_scan_response_fit_executed = true,
            bounded_cable_frequency_scan_recursive_convolution_runtime_executed = true,
            bounded_over47_skin_effect_executed =
                any(row -> row.bounded_skin_effect_executed, generated.series_impedance.frequency_rows),
            skin_effect_internal_impedance_executed =
                skin_effect_internal_impedance_executed,
            bounded_over47_earth_return_executed =
                any(row -> row.bounded_earth_return_executed, generated.series_impedance.frequency_rows),
            earth_return_impedance_executed =
                earth_return_impedance_executed,
            bounded_over47_grounded_conductor_reduction_executed =
                any(row -> row.grounded_conductor_reduction_executed, generated.series_impedance.frequency_rows),
            full_over47_skin_effect_executed =
                skin_effect_internal_impedance_executed,
            full_over47_earth_return_executed =
                earth_return_impedance_executed,
            cable_pipe_sheath_state = pipe_sheath_state,
            pipe_sheath_side_effects_executed =
                pipe_sheath_side_effects_executed,
            full_over47_pipe_or_sheath_coupling_executed =
                pipe_sheath_side_effects_executed,
            frequency_loop_executed = true,
            full_over47_frequency_loop_executed = true,
            full_over47_frequency_loop_side_effects_executed = true,
            deferred_calls = Tuple(runtime_deferred_calls),
        ),
    )
end

function line_characteristic_impedance(r_per_length::Real, l_per_length::Real, g_per_length::Real, c_per_length::Real, frequency_hz::Real)
    omega = 2.0 * pi * Float64(frequency_hz)
    z = complex(Float64(r_per_length), omega * Float64(l_per_length))
    y = complex(Float64(g_per_length), omega * Float64(c_per_length))
    y != 0.0 + 0.0im || throw(ArgumentError("line shunt admittance must be nonzero"))
    return sqrt(z / y)
end

function line_propagation_constant(r_per_length::Real, l_per_length::Real, g_per_length::Real, c_per_length::Real, frequency_hz::Real)
    omega = 2.0 * pi * Float64(frequency_hz)
    z = complex(Float64(r_per_length), omega * Float64(l_per_length))
    y = complex(Float64(g_per_length), omega * Float64(c_per_length))
    return sqrt(z * y)
end

function frequency_dependent_line_point(
    r_per_length::Real,
    l_per_length::Real,
    g_per_length::Real,
    c_per_length::Real,
    length::Real,
    frequency_hz::Real,
)
    line_length = Float64(length)
    isfinite(line_length) && line_length >= 0.0 ||
        throw(ArgumentError("length must be finite and nonnegative"))
    freq = Float64(frequency_hz)
    isfinite(freq) && freq >= 0.0 ||
        throw(ArgumentError("frequency_hz must be finite and nonnegative"))
    zc = line_characteristic_impedance(r_per_length, l_per_length, g_per_length, c_per_length, freq)
    gamma = line_propagation_constant(r_per_length, l_per_length, g_per_length, c_per_length, freq)
    return LineFrequencyPoint(freq, zc, gamma, exp(-gamma * line_length))
end

function semlyen_marti_frequency_scan(
    start_frequency_hz::Real,
    decade_count::Integer,
    intervals_per_decade::Integer,
)
    start_frequency = Float64(start_frequency_hz)
    isfinite(start_frequency) && start_frequency > 0.0 ||
        throw(ArgumentError("start_frequency_hz must be finite and positive"))
    decades = Int(decade_count)
    decades >= 0 ||
        throw(ArgumentError("decade_count must be nonnegative"))
    intervals = Int(intervals_per_decade)
    intervals > 0 ||
        throw(ArgumentError("intervals_per_decade must be positive"))

    total_intervals = decades * intervals
    frequency_count = total_intervals + 1
    ratio = exp(log(10.0) / intervals)
    frequencies = Vector{Float64}(undef, frequency_count)
    base_frequencies = Vector{Float64}(undef, frequency_count)
    decade_indices = Vector{Int}(undef, frequency_count)
    interval_indices = Vector{Int}(undef, frequency_count)

    for row in 1:frequency_count
        step = row - 1
        decade_index = step == 0 ? 0 : div(step - 1, intervals)
        interval_index = step == 0 ? 0 : mod(step - 1, intervals) + 1
        base_frequency = start_frequency * 10.0^decade_index
        frequencies[row] = step == 0 ? start_frequency : base_frequency * ratio^interval_index
        base_frequencies[row] = base_frequency
        decade_indices[row] = decade_index
        interval_indices[row] = interval_index
    end

    return SemlyenMartiFrequencyScan(
        start_frequency,
        decades,
        intervals,
        frequency_count,
        frequencies,
        base_frequencies,
        decade_indices,
        interval_indices,
        ratio,
    )
end

function _checked_frequency_dependent_line_mode_values(values::AbstractVector, label::AbstractString)
    isempty(values) && throw(ArgumentError("$label must contain at least one mode"))
    checked = Float64.(values)
    all(isfinite, checked) ||
        throw(ArgumentError("$label entries must be finite"))
    return checked
end

function semlyen_marti_frequency_dependent_line_samples(
    scan::SemlyenMartiFrequencyScan,
    resistance_per_length::AbstractVector,
    inductance_per_length::AbstractVector,
    conductance_per_length::AbstractVector,
    capacitance_per_length::AbstractVector,
    line_length::Real,
)
    resistance = _checked_frequency_dependent_line_mode_values(
        resistance_per_length,
        "resistance_per_length",
    )
    inductance = _checked_frequency_dependent_line_mode_values(
        inductance_per_length,
        "inductance_per_length",
    )
    conductance = _checked_frequency_dependent_line_mode_values(
        conductance_per_length,
        "conductance_per_length",
    )
    capacitance = _checked_frequency_dependent_line_mode_values(
        capacitance_per_length,
        "capacitance_per_length",
    )
    mode_count = length(resistance)
    length(inductance) == mode_count &&
        length(conductance) == mode_count &&
        length(capacitance) == mode_count ||
        throw(ArgumentError("frequency-dependent line mode parameter counts must match"))

    return [
        [
            frequency_dependent_line_point(
                resistance[mode],
                inductance[mode],
                conductance[mode],
                capacitance[mode],
                line_length,
                frequency,
            )
            for mode in 1:mode_count
        ]
        for frequency in scan.frequencies_hz
    ]
end

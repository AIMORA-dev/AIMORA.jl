export hysteresis_runtime_table

function _deck_hysteretic_inductor_point_groups(
    rows::Vector{DeckParser.DeckHystereticInductorRow},
    points::Vector{DeckParser.DeckHystereticInductorPointRow},
)
    groups = [DeckParser.DeckHystereticInductorPointRow[] for _ in rows]
    for point in points
        1 <= point.hysteretic_inductor_row_index <= length(rows) ||
            throw(ArgumentError("type-96 point references an unknown hysteretic inductor"))
        group = groups[point.hysteretic_inductor_row_index]
        point.point_index == length(group) + 1 ||
            throw(ArgumentError("type-96 current/flux points must be consecutive"))
        push!(group, point)
    end
    for (row_index, row) in enumerate(rows)
        if row.source_kind == :copy_reference
            1 <= row.reference_index < row_index ||
                throw(ArgumentError("type-96 COPY rows must reference a prior owner"))
            isempty(groups[row_index]) ||
                throw(ArgumentError("type-96 COPY rows must not duplicate a major-loop table"))
            groups[row_index] = groups[row.reference_index]
        end
        length(groups[row_index]) >= 3 ||
            throw(ArgumentError("every type-96 hysteretic inductor requires at least three points"))
    end
    return groups
end

function _deck_hysteretic_inductor_table(
    points::Vector{DeckParser.DeckHystereticInductorPointRow},
    delta2::Float64,
)
    return _hysteretic_inductor_table(
        Float64[point.current_A for point in points],
        Float64[point.flux_Wb for point in points],
        delta2,
    )
end

function _hysteretic_inductor_table(
    currents::Vector{Float64},
    fluxes::Vector{Float64},
    delta2::Float64,
)
    point_count = length(currents)
    point_count == length(fluxes) >= 3 ||
        throw(ArgumentError("type-96 runtime table requires at least three aligned points"))
    all(isfinite, currents) && all(isfinite, fluxes) ||
        throw(ArgumentError("type-96 runtime points must be finite"))
    state_start = 1
    major_loop_start = state_start + 6
    table_length = 2 * point_count + 8
    cchar = zeros(Float64, table_length)
    vchar = zeros(Float64, table_length)
    gslope = zeros(Float64, table_length)

    major_currents = vcat(-currents[end - 1], currents)
    major_fluxes = vcat(-fluxes[end - 1], fluxes)
    major_stop = major_loop_start + point_count
    cchar[major_loop_start:major_stop] .= major_currents
    vchar[major_loop_start:major_stop] .= major_fluxes

    extension_start = major_stop + 1
    for segment in 1:point_count
        left = major_loop_start + segment - 1
        right = left + 1
        extension = extension_start + segment - 1
        current_delta = cchar[right] - cchar[left]
        flux_delta = vchar[right] - vchar[left]
        current_delta > 0.0 && flux_delta > 0.0 ||
            throw(ArgumentError("type-96 major-loop current and flux must increase"))
        flux_per_current = flux_delta / current_delta
        flux_intercept = vchar[right] - flux_per_current * cchar[right]
        vchar[extension] = flux_per_current
        cchar[extension] = flux_intercept
        gslope[left] = 1.0 / flux_per_current
        gslope[extension] = -flux_intercept / flux_per_current
    end
    final_extension = extension_start + point_count
    cchar[final_extension] = -cchar[extension_start]
    vchar[final_extension] = vchar[extension_start]
    gslope[major_stop] = gslope[major_stop - 1]
    gslope[final_extension] = -gslope[extension_start]

    zero_point = findfirst(==(0.0), currents)
    zero_point === nothing &&
        throw(ArgumentError("type-96 major loop must contain zero current"))
    trace_index = major_loop_start + zero_point
    cchar[state_start] = point_count
    cchar[state_start + 1] = 1.0
    cchar[state_start + 2] = trace_index
    cchar[state_start + 3] = 0.0
    cchar[state_start + 4] = -1.0
    cchar[state_start + 5] = 0.0
    vchar[state_start + 5] = fluxes[end - 1]
    gslope[state_start] = 0.0
    gslope[state_start + 1] = delta2 * gslope[trace_index]
    gslope[state_start + 2] = 1.0
    gslope[state_start + 3] = trace_index
    gslope[state_start + 4] = 0.0
    gslope[state_start + 5] = currents[end - 1]
    return cchar, vchar, gslope, state_start, major_loop_start
end

function _hysteretic_inductor_major_loop_trace_index(
    cchar::Vector{Float64},
    major_loop_start::Int,
    active_point_count::Int,
    current_a::Float64,
)
    major_loop_stop = major_loop_start + active_point_count - 1
    for index in major_loop_start:major_loop_stop
        current_a <= cchar[index] && return index
    end
    return major_loop_stop + 1
end

function _hysteretic_inductor_major_loop_flux_bounds(
    cchar::Vector{Float64},
    vchar::Vector{Float64},
    state_start::Int,
    major_loop_start::Int,
    current_a::Float64,
)
    active_point_count = round(Int, cchar[state_start])
    upper_trace = _hysteretic_inductor_major_loop_trace_index(
        cchar,
        major_loop_start,
        active_point_count,
        current_a,
    )
    lower_trace = _hysteretic_inductor_major_loop_trace_index(
        cchar,
        major_loop_start,
        active_point_count,
        -current_a,
    )
    upper_extension = upper_trace + active_point_count + 1
    lower_extension = lower_trace + active_point_count + 1
    upper_flux_wb =
        vchar[upper_extension] * current_a + cchar[upper_extension]
    lower_flux_wb =
        vchar[lower_extension] * current_a - cchar[lower_extension]
    return (
        lower_wb=min(lower_flux_wb, upper_flux_wb),
        upper_wb=max(lower_flux_wb, upper_flux_wb),
    )
end

function _hysteretic_inductor_steady_state_table(
    cchar::AbstractVector{<:Real},
    vchar::AbstractVector{<:Real},
    gslope::AbstractVector{<:Real};
    state_start_index::Int,
    major_loop_start_index::Int,
    branch_current_a::Real,
    branch_flux_wb::Real,
    branch_voltage_v::Real,
    residual_flux_wb::Real,
    half_timestep_s::Real,
    flux_tolerance_wb::Real,
)
    cchar_values = Float64.(cchar)
    vchar_values = Float64.(vchar)
    gslope_values = Float64.(gslope)
    state_start = state_start_index
    major_loop_start = major_loop_start_index
    current_a = Float64(branch_current_a)
    harmonic_flux_wb = Float64(branch_flux_wb)
    voltage_v = Float64(branch_voltage_v)
    residual_flux = Float64(residual_flux_wb)
    delta2 = Float64(half_timestep_s)
    tolerance = Float64(flux_tolerance_wb)
    all(isfinite, (current_a, harmonic_flux_wb, voltage_v, residual_flux, delta2, tolerance)) ||
        throw(ArgumentError("hysteretic-inductor initialization values must be finite"))
    delta2 > 0.0 || throw(ArgumentError(
        "hysteretic-inductor initialization half timestep must be positive",
    ))
    tolerance >= 0.0 || throw(ArgumentError(
        "hysteretic-inductor initialization flux tolerance must be nonnegative",
    ))
    state_start >= 1 && state_start + 5 <= length(cchar_values) ||
        throw(ArgumentError("hysteretic-inductor state table is incomplete"))
    active_point_count = round(Int, cchar_values[state_start])
    active_point_count > 0 || throw(ArgumentError(
        "hysteretic-inductor active major-loop point count must be positive",
    ))
    major_loop_stop = major_loop_start + active_point_count - 1
    final_extension = major_loop_stop + active_point_count + 2
    final_extension <= length(cchar_values) &&
        final_extension <= length(vchar_values) &&
        final_extension <= length(gslope_values) ||
        throw(ArgumentError("hysteretic-inductor major-loop table is incomplete"))

    flux_wb = abs(harmonic_flux_wb) <= tolerance ? residual_flux : harmonic_flux_wb
    bounds = _hysteretic_inductor_major_loop_flux_bounds(
        cchar_values,
        vchar_values,
        state_start,
        major_loop_start,
        current_a,
    )
    bounds.lower_wb - tolerance <= flux_wb <= bounds.upper_wb + tolerance ||
        throw(ArgumentError(
            "hysteretic-inductor harmonic flux-current point lies outside the major loop",
        ))
    flux_wb = clamp(flux_wb, bounds.lower_wb, bounds.upper_wb)

    trace_slot = state_start + 2
    cchar_values[state_start + 3] = current_a
    cchar_values[state_start + 4] = -1.0
    vchar_values[state_start + 4] = flux_wb
    gslope_values[state_start + 4] = current_a
    direction = voltage_v < 0.0 ? -1.0 : 1.0
    cchar_values[state_start + 1] = direction
    boundary_index = direction > 0.0 ? major_loop_stop : major_loop_start
    boundary_flux_wb = vchar_values[boundary_index]
    boundary_current_a = cchar_values[boundary_index]
    vchar_values[state_start + 5] = boundary_flux_wb
    gslope_values[state_start + 5] = boundary_current_a

    trace_index = if abs(boundary_flux_wb - flux_wb) <= tolerance
        vchar_values[state_start] = 0.0
        vchar_values[state_start + 1] = 0.0
        major_loop_stop + 1
    else
        candidate = _hysteretic_inductor_major_loop_trace_index(
            cchar_values,
            major_loop_start,
            active_point_count,
            direction > 0.0 ? current_a : -current_a,
        )
        extension = candidate + active_point_count + 1
        major_flux_wb = if direction > 0.0
            vchar_values[extension] * current_a + cchar_values[extension]
        else
            vchar_values[extension] * current_a - cchar_values[extension]
        end
        trajectory_offset = direction > 0.0 ?
            flux_wb - major_flux_wb :
            major_flux_wb - flux_wb
        trajectory_denominator = flux_wb - boundary_flux_wb
        abs(trajectory_denominator) > tolerance || throw(ArgumentError(
            "hysteretic-inductor initial trajectory denominator is zero",
        ))
        vchar_values[state_start] =
            trajectory_offset / trajectory_denominator
        vchar_values[state_start + 1] =
            trajectory_offset - vchar_values[state_start] * flux_wb
        candidate
    end
    cchar_values[trace_slot] = Float64(trace_index)
    extension_index = trace_index + active_point_count + 1
    inverse_incremental_slope = if direction > 0.0
        denominator = gslope_values[trace_index] *
            (1.0 - vchar_values[state_start])
        denominator != 0.0 || throw(ArgumentError(
            "hysteretic-inductor upper initial slope denominator is zero",
        ))
        inv(denominator)
    else
        denominator = gslope_values[trace_index] *
            (1.0 + vchar_values[state_start])
        denominator != 0.0 || throw(ArgumentError(
            "hysteretic-inductor lower initial slope denominator is zero",
        ))
        inv(denominator)
    end
    isfinite(inverse_incremental_slope) && extension_index <= length(gslope_values) ||
        throw(ArgumentError("hysteretic-inductor initial incremental slope is invalid"))
    companion_admittance_s = delta2 / inverse_incremental_slope
    companion_current_a = current_a - voltage_v * companion_admittance_s
    gslope_values[state_start] = companion_current_a
    gslope_values[state_start + 1] = companion_admittance_s
    vchar_values[state_start + 2] = flux_wb
    vchar_values[state_start + 3] = current_a
    gslope_values[state_start + 2] = direction
    gslope_values[state_start + 3] = Float64(trace_index)
    runtime_flux_wb = flux_wb - voltage_v * delta2
    return (
        cchar=cchar_values,
        vchar=vchar_values,
        gslope=gslope_values,
        current_a=current_a,
        flux_wb=flux_wb,
        runtime_flux_wb=runtime_flux_wb,
        companion_current_a=companion_current_a,
        companion_admittance_s=companion_admittance_s,
        direction=Int(direction),
        trace_index=trace_index,
        lower_major_loop_flux_wb=bounds.lower_wb,
        upper_major_loop_flux_wb=bounds.upper_wb,
    )
end

function hysteresis_runtime_table(
    characteristic::HysteresisLoopPreprocessResult;
    half_timestep_s::Real,
)
    characteristic.physical_checks_passed ||
        throw(ArgumentError("hysteresis preprocessing must pass physical checks"))
    delta2 = Float64(half_timestep_s)
    isfinite(delta2) && delta2 > 0.0 ||
        throw(ArgumentError("hysteresis half_timestep_s must be finite and positive"))
    cchar, vchar, gslope, state_start, major_loop_start =
        _hysteretic_inductor_table(
            characteristic.runtime_current_a,
            characteristic.runtime_flux_wb,
            delta2,
        )
    return (
        cchar,
        vchar,
        gslope,
        state_start_index = state_start,
        major_loop_start_index = major_loop_start,
        table_end_index = length(cchar),
        active_runtime_owner = :hysteretic_inductor_current_update,
        physical_checks_passed =
            gslope[state_start + 1] > 0.0 &&
            characteristic.energy_loss_j > 0.0,
    )
end

function deck_hysteretic_inductor_nonlinear_current_config(
    parsed::DeckParser.DeckParseResult;
    delta2::Real = 1.0,
    ncomp::Union{Nothing,Int} = nothing,
    t::Real = 0.0,
    deltat::Real = 2.0 * Float64(delta2),
)
    DeckParser.assert_deck_valid!(parsed)
    rows = DeckParser.deck_hysteretic_inductor_rows(parsed)
    count = length(rows)
    count > 0 || throw(ArgumentError("deck has no type-96 hysteretic inductors"))
    delta = Float64(delta2)
    isfinite(delta) && delta > 0.0 ||
        throw(ArgumentError("delta2 must be finite and positive"))
    groups = _deck_hysteretic_inductor_point_groups(
        rows,
        DeckParser.deck_hysteretic_inductor_point_rows(parsed),
    )

    cchar = Float64[]
    vchar = Float64[]
    gslope = Float64[]
    state_start_indices = zeros(Int, count)
    major_loop_start_indices = zeros(Int, count)
    table_end_indices = zeros(Int, count)
    for index in eachindex(rows)
        local_cchar, local_vchar, local_gslope, local_state, local_major =
            _deck_hysteretic_inductor_table(groups[index], delta)
        offset = length(cchar)
        append!(cchar, local_cchar)
        append!(vchar, local_vchar)
        append!(gslope, local_gslope)
        state_start_indices[index] = offset + local_state
        major_loop_start_indices[index] = offset + local_major
        table_end_indices[index] = length(cchar)
        cchar[state_start_indices[index] + 2] += offset
        gslope[state_start_indices[index] + 3] += offset
    end

    deck_node_count = length(parsed.node_map)
    shifted_from_nodes = [
        _deck_reference_shifted_runtime_node(row.from_node_index, deck_node_count)
        for row in rows
    ]
    shifted_to_nodes = [
        _deck_reference_shifted_runtime_node(row.to_node_index, deck_node_count)
        for row in rows
    ]
    deck_to_runtime_node_indices = [
        _deck_reference_shifted_runtime_node(index, deck_node_count)
        for index in 1:deck_node_count
    ]
    required_node_count = max(
        deck_node_count + 1,
        maximum(vcat(shifted_from_nodes, shifted_to_nodes); init = 1),
    )
    source_activity_flags = zeros(Int, required_node_count)
    for node in vcat(shifted_from_nodes, shifted_to_nodes)
        node > 1 && (source_activity_flags[node] = 1)
    end
    last_slot = 1 + 5 * (count - 1)
    source_next_indices = zeros(Int, last_slot)
    source_from_nodes = zeros(Int, last_slot)
    source_to_nodes = zeros(Int, last_slot)
    source_begin_indices = zeros(Int, count)
    subsystem_owner_rows = zeros(Int, last_slot)
    subsystem_simultaneous_flags = zeros(Int, last_slot)
    for index in 1:count
        slot = 1 + 5 * (index - 1)
        source_next_indices[slot] = slot
        source_from_nodes[slot] = shifted_from_nodes[index]
        source_to_nodes[slot] = shifted_to_nodes[index]
        source_begin_indices[index] = slot
        subsystem_owner_rows[slot] = index
    end
    component_count = ncomp === nothing ? count : ncomp
    component_count > 0 || throw(ArgumentError("ncomp must be positive"))
    steady_flux = Float64[row.steady_state_flux_Wb for row in rows]
    steady_current = Float64[row.steady_state_current_A for row in rows]

    return (
        source = :deck_hysteretic_inductor_nonlinear_current_config,
        outcome = :nonlinear_current_config,
        nonlinear_types = fill(-96, count),
        nonlinear_from_nodes = shifted_from_nodes,
        nonlinear_to_nodes = shifted_to_nodes,
        nonlinear_admittance_nodes = state_start_indices,
        nonlinear_table_end_indices = table_end_indices,
        nonlinear_subsystem_indices = collect(1:count),
        subsystem_begin_indices = source_begin_indices,
        subsystem_owner_rows = subsystem_owner_rows,
        subsystem_simultaneous_flags = subsystem_simultaneous_flags,
        cchar = cchar,
        vchar = vchar,
        gslope = gslope,
        delta2 = delta,
        flzero = 1.0e-12,
        ncomp = component_count,
        fltinf = Inf,
        t = Float64(t),
        deltat = Float64(deltat),
        deck_owned_voltage_context = true,
        reference_shifted_voltage_context = true,
        deck_to_runtime_node_indices = deck_to_runtime_node_indices,
        reference_node_index = 1,
        nonlinear_required_node_count = required_node_count,
        initialize_nonlinear_state = true,
        initial_companion_current_values = steady_current,
        initial_characteristic_current_values = steady_current,
        initial_stored_voltage_values = steady_flux,
        initial_runtime_voltage_values = steady_flux,
        initial_current_segment_values = steady_current,
        initial_table_index_values = major_loop_start_indices,
        initial_cursub_values = zeros(Float64, required_node_count),
        nonlinear_current_segments = steady_current,
        nonlinear_output_codes = [row.output_code for row in rows],
        nonlinear_owner_names = [row.name for row in rows],
        nonlinear_owner_line_numbers = [row.line_no for row in rows],
        nonlinear_reference_indices = [row.reference_index for row in rows],
        nonlinear_source_kinds = [row.source_kind for row in rows],
        nonlinear_characteristic_owner_indices = [
            row.source_kind == :copy_reference ? row.reference_index : index
            for (index, row) in enumerate(rows)
        ],
        dense_primary_nonlinear_compensation = true,
        use_state_nonlinear_inverse_columns = false,
        nonlinear_inverse_config = (
            require_fortran_sparse_factor_result = true,
            ntot = required_node_count,
            ncomp = component_count,
            source_begin_indices = source_begin_indices,
            source_next_indices = source_next_indices,
            source_from_nodes = source_from_nodes,
            source_to_nodes = source_to_nodes,
            source_activity_flags = source_activity_flags,
            nonlinear_types = fill(-96, count),
            nonlinear_admittance_nodes = state_start_indices,
            nonlinear_from_nodes = shifted_from_nodes,
            nonlinear_to_nodes = shifted_to_nodes,
            nonlinear_source_flags = zeros(Int, count),
            partition_boundary = max(required_node_count - 1, 1),
            delta2 = delta,
            fltinf = Inf,
        ),
        nonlinear_source_begin_indices = source_begin_indices,
        nonlinear_source_next_indices = source_next_indices,
        nonlinear_source_from_nodes = source_from_nodes,
        nonlinear_source_to_nodes = source_to_nodes,
        nonlinear_source_activity_flags = source_activity_flags,
        nonlinear_source_flags = zeros(Int, count),
        nonlinear_inverse_partition_boundary = max(required_node_count - 1, 1),
        nonlinear_deck_from_nodes = [row.from_node_index for row in rows],
        nonlinear_deck_to_nodes = [row.to_node_index for row in rows],
        nonlinear_deck_from_node_names = [row.from_node for row in rows],
        nonlinear_deck_to_node_names = [row.to_node for row in rows],
        fortran_nonlinear_admittance_nodes = state_start_indices,
        fortran_nonlinear_state_start_indices = state_start_indices,
        fortran_gap_status_values = [row.residual_flux_Wb for row in rows],
        table_entry_count = length(cchar),
        nonlinear_row_count = count,
        fortran_files = (:OVER2_FOR, :OVER8_FOR, :OVER16_FOR),
        fortran_labels = (4221, 73418, 187, 188, 189, 5300, 5380, 4317, 1116, 1118, 1110, 1123, 1125, 1140, 1141, 1149, 1150, 1212, 1312, 1315, 1322, 4372, 3950),
        mutation_order = (
            :fixed_card_intake,
            :major_loop_table,
            :steady_state_seed,
            :branch_voltage,
            :flux_current_trajectory,
            :admittance_source_restamp,
            :accepted_nodal_solve,
        ),
        deferred_calls = Symbol[],
        complete_nonlinear_source_loop = true,
        replacement_ready = true,
    )
end

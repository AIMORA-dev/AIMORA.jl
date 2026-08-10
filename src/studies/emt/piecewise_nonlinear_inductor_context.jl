export deck_piecewise_nonlinear_inductor_current_config

function _piecewise_nonlinear_inductor_point_groups(
    rows::Vector{DeckParser.DeckPiecewiseNonlinearInductorRow},
    points::Vector{DeckParser.DeckPiecewiseNonlinearInductorPointRow},
)
    groups = [DeckParser.DeckPiecewiseNonlinearInductorPointRow[] for _ in rows]
    for point in points
        1 <= point.inductor_row_index <= length(rows) ||
            throw(ArgumentError("nonlinear-inductor point references an unknown owner"))
        push!(groups[point.inductor_row_index], point)
    end
    for index in eachindex(rows)
        if rows[index].source_kind == :copy_reference
            reference = rows[index].reference_index
            1 <= reference < index ||
                throw(ArgumentError("nonlinear-inductor copy must reference an earlier owner"))
            groups[index] = groups[reference]
        else
            sort!(groups[index]; by = point -> point.point_index)
            length(groups[index]) >= 2 ||
                throw(ArgumentError("nonlinear inductor requires at least two points"))
            all(
                point.point_index == position
                for (position, point) in enumerate(groups[index])
            ) || throw(ArgumentError("nonlinear-inductor point indices must be consecutive"))
        end
    end
    return groups
end

function _piecewise_nonlinear_inductor_initial_segment(
    current_a::Float64,
    points::Vector{DeckParser.DeckPiecewiseNonlinearInductorPointRow},
)
    current_a <= first(points).current_a && return 1
    current_a >= last(points).current_a && return length(points) - 1
    return clamp(
        searchsortedlast([point.current_a for point in points], current_a),
        1,
        length(points) - 1,
    )
end

function _append_pseudo_nonlinear_inductor_table!(
    cchar::Vector{Float64},
    vchar::Vector{Float64},
    gslope::Vector{Float64},
    points::Vector{DeckParser.DeckPiecewiseNonlinearInductorPointRow},
    delta2::Float64,
)
    previous_current = 0.0
    previous_flux = 0.0
    for point in points
        point.current_a > previous_current ||
            throw(ArgumentError("type-98 nonlinear-inductor current points must increase from zero"))
        point.flux_wb > previous_flux ||
            throw(ArgumentError("type-98 nonlinear-inductor flux points must increase from zero"))
        current_delta = point.current_a - previous_current
        flux_delta = point.flux_wb - previous_flux
        flux_per_current = flux_delta / current_delta
        push!(cchar, previous_flux - flux_per_current * previous_current)
        push!(vchar, point.flux_wb)
        push!(gslope, delta2 * current_delta / flux_delta)
        previous_current = point.current_a
        previous_flux = point.flux_wb
    end
    return cchar
end

function deck_piecewise_nonlinear_inductor_current_config(
    parsed::DeckParser.DeckParseResult;
    delta2::Real=1.0,
)
    DeckParser.assert_deck_valid!(parsed)
    rows = DeckParser.deck_piecewise_nonlinear_inductor_rows(parsed)
    isempty(rows) && throw(ArgumentError("deck has no piecewise nonlinear inductors"))
    groups = _piecewise_nonlinear_inductor_point_groups(
        rows,
        DeckParser.deck_piecewise_nonlinear_inductor_point_rows(parsed),
    )
    delta = Float64(delta2)
    isfinite(delta) && delta > 0.0 ||
        throw(ArgumentError("delta2 must be finite and positive"))

    cchar = Float64[]
    vchar = Float64[]
    gslope = Float64[]
    table_starts = zeros(Int, length(rows))
    table_ends = zeros(Int, length(rows))
    for index in eachindex(rows)
        if rows[index].source_kind == :copy_reference
            reference = rows[index].reference_index
            table_starts[index] = table_starts[reference]
            table_ends[index] = table_ends[reference]
            continue
        end
        table_starts[index] = length(cchar) + 1
        if rows[index].nonlinear_type == PIECEWISE_NONLINEAR_INDUCTOR_TYPE
            append!(cchar, (point.current_a for point in groups[index]))
            append!(vchar, (point.flux_wb for point in groups[index]))
            append!(gslope, zeros(Float64, length(groups[index])))
        elseif rows[index].nonlinear_type == PSEUDO_NONLINEAR_INDUCTOR_TYPE
            _append_pseudo_nonlinear_inductor_table!(
                cchar,
                vchar,
                gslope,
                groups[index],
                delta,
            )
        else
            throw(ArgumentError("unsupported nonlinear-inductor input type"))
        end
        table_ends[index] = length(cchar)
    end
    initial_segments = [
        rows[index].nonlinear_type == PIECEWISE_NONLINEAR_INDUCTOR_TYPE ?
            table_starts[index] - 1 + _piecewise_nonlinear_inductor_initial_segment(
                rows[index].steady_state_current_a,
                groups[index],
            ) :
            table_starts[index]
        for index in eachindex(rows)
    ]

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
    count = length(rows)
    last_slot = 1 + 5 * (count - 1)
    source_next_indices = zeros(Int, last_slot)
    source_from_nodes = zeros(Int, last_slot)
    source_to_nodes = zeros(Int, last_slot)
    subsystem_begin_indices = zeros(Int, count)
    subsystem_owner_rows = zeros(Int, last_slot)
    subsystem_simultaneous_flags = zeros(Int, last_slot)
    for index in 1:count
        slot = 1 + 5 * (index - 1)
        source_next_indices[slot] = index == count ? 1 : slot + 5
        source_from_nodes[slot] = shifted_from_nodes[index]
        source_to_nodes[slot] = shifted_to_nodes[index]
        subsystem_begin_indices[index] = slot
        subsystem_owner_rows[slot] = index
    end
    initial_currents = [row.steady_state_current_a for row in rows]
    initial_fluxes = [row.steady_state_flux_wb for row in rows]
    input_types = [row.nonlinear_type for row in rows]
    runtime_segments = [
        input_type == PSEUDO_NONLINEAR_INDUCTOR_TYPE ? 1.0 : initial_currents[index]
        for (index, input_type) in enumerate(input_types)
    ]
    runtime_stored_fluxes = [
        input_type == PSEUDO_NONLINEAR_INDUCTOR_TYPE ? 0.0 : initial_fluxes[index]
        for (index, input_type) in enumerate(input_types)
    ]
    initial_characteristic_values = [
        input_type == PSEUDO_NONLINEAR_INDUCTOR_TYPE ?
            initial_currents[index] :
            initial_fluxes[index]
        for (index, input_type) in enumerate(input_types)
    ]
    return (
        source = :deck_piecewise_nonlinear_inductor_current_config,
        outcome = :nonlinear_current_config,
        nonlinear_types = input_types,
        nonlinear_from_nodes = shifted_from_nodes,
        nonlinear_to_nodes = shifted_to_nodes,
        nonlinear_admittance_nodes = table_starts,
        nonlinear_table_end_indices = table_ends,
        nonlinear_subsystem_indices = collect(1:count),
        subsystem_begin_indices = subsystem_begin_indices,
        subsystem_owner_rows = subsystem_owner_rows,
        subsystem_simultaneous_flags = subsystem_simultaneous_flags,
        cchar = cchar,
        vchar = vchar,
        gslope = gslope,
        flzero = 1.0e-12,
        delta2 = delta,
        deltat = 2.0 * delta,
        ncomp = count,
        deck_owned_voltage_context = true,
        reference_shifted_voltage_context = true,
        deck_to_runtime_node_indices = deck_to_runtime_node_indices,
        reference_node_index = 1,
        nonlinear_required_node_count = required_node_count,
        initialize_nonlinear_state = true,
        initial_companion_current_values = zeros(Float64, count),
        nonlinear_steady_state_current_values = initial_currents,
        nonlinear_steady_state_flux_values = initial_fluxes,
        initial_characteristic_current_values = initial_characteristic_values,
        initial_stored_voltage_values = runtime_stored_fluxes,
        initial_runtime_voltage_values = copy(runtime_stored_fluxes),
        initial_current_segment_values = runtime_segments,
        initial_table_index_values = initial_segments,
        initial_cursub_values = zeros(Float64, max(required_node_count, count + 1)),
        nonlinear_current_segments = runtime_segments,
        nonlinear_output_codes = [row.output_code for row in rows],
        nonlinear_owner_names = [row.name for row in rows],
        nonlinear_owner_line_numbers = [row.line_no for row in rows],
        nonlinear_source_begin_indices = [1],
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
        subnetwork_next_indices = Int[],
        subnetwork_from_nodes = Int[],
        subnetwork_to_nodes = Int[],
        subnetwork_nonlinear_indices = Int[],
        subnetwork_element_types = Int[],
        table_entry_count = length(cchar),
        nonlinear_row_count = count,
        true_nonlinear_inductor_count =
            Base.count(==(PIECEWISE_NONLINEAR_INDUCTOR_TYPE), input_types),
        pseudo_nonlinear_inductor_count =
            Base.count(==(PSEUDO_NONLINEAR_INDUCTOR_TYPE), input_types),
        seed_initial_nonlinear_state = any(==(PSEUDO_NONLINEAR_INDUCTOR_TYPE), input_types),
        dense_primary_nonlinear_compensation = true,
        mutation_order = (
            :flux_current_input,
            :network_response_columns,
            :coupled_piecewise_solution,
            :current_injection,
            :corrected_network_solve,
            :current_flux_output,
        ),
        deferred_calls = (),
        complete_nonlinear_source_loop = true,
    )
end

function _piecewise_nonlinear_inductor_report(run, config::NamedTuple)
    types = Int.(config.nonlinear_types)
    indices = findall(
        type_code ->
            type_code in (
                PIECEWISE_NONLINEAR_INDUCTOR_TYPE,
                PSEUDO_NONLINEAR_INDUCTOR_TYPE,
            ),
        types,
    )
    names = Symbol.(config.nonlinear_owner_names)
    line_numbers = Int.(config.nonlinear_owner_line_numbers)
    output_codes = Int.(config.nonlinear_output_codes)
    from_nodes = Int.(config.nonlinear_from_nodes)
    to_nodes = abs.(Int.(config.nonlinear_to_nodes))
    node_mapping = Int.(config.deck_to_runtime_node_indices)
    times = Float64[]
    currents = zeros(Float64, length(indices), length(run.over16_updates))
    fluxes = similar(currents)
    voltages = similar(currents)
    segments = zeros(Int, size(currents))
    for (step, update) in enumerate(run.over16_updates)
        push!(times, Float64(update.t_s))
        result = update.over16_update.nonlinear_current_result
        result === nothing && throw(ArgumentError(
            "nonlinear-inductor report requires a nonlinear timestep result",
        ))
        reference_voltage = zeros(Float64, Int(config.nonlinear_required_node_count))
        for deck_index in eachindex(node_mapping)
            reference_voltage[node_mapping[deck_index]] = update.voltage_pu[deck_index]
        end
        for (position, nonlinear_index) in enumerate(indices)
            if types[nonlinear_index] == PIECEWISE_NONLINEAR_INDUCTOR_TYPE
                currents[position, step] =
                    result.piecewise_nonlinear_inductor_accepted_currents[nonlinear_index]
                fluxes[position, step] =
                    result.piecewise_nonlinear_inductor_accepted_fluxes[nonlinear_index]
                voltages[position, step] =
                    result.piecewise_nonlinear_inductor_accepted_voltages[nonlinear_index]
                segments[position, step] = result.ilast[nonlinear_index]
            else
                voltage = reference_voltage[from_nodes[nonlinear_index]] -
                    reference_voltage[to_nodes[nonlinear_index]]
                segment = round(Int, result.curr[nonlinear_index])
                table_index =
                    config.nonlinear_admittance_nodes[nonlinear_index] +
                    abs(segment) - 1
                config.nonlinear_admittance_nodes[nonlinear_index] <= table_index <=
                    config.nonlinear_table_end_indices[nonlinear_index] ||
                    throw(ArgumentError("type-98 report segment is outside its characteristic"))
                voltages[position, step] = voltage
                currents[position, step] =
                    result.pseudo_nonlinear_inductor_accepted_currents[
                        nonlinear_index
                    ]
                fluxes[position, step] = result.vnonl[nonlinear_index]
                segments[position, step] = segment
            end
        end
    end
    return (
        source = :piecewise_nonlinear_inductor_report,
        outcome = :report_output,
        names = names[indices],
        line_numbers = line_numbers[indices],
        nonlinear_types = types[indices],
        output_codes = output_codes[indices],
        output_requested_flags = output_codes[indices] .> 0,
        time_s = times,
        voltage_v = voltages,
        current_a = currents,
        flux_wb = fluxes,
        active_segments = segments,
        deferred_effects = (),
    )
end

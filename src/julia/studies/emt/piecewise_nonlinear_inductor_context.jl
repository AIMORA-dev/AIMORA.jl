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

    currents_a = Float64[]
    fluxes_wb = Float64[]
    table_starts = zeros(Int, length(rows))
    table_ends = zeros(Int, length(rows))
    for index in eachindex(rows)
        if rows[index].source_kind == :copy_reference
            reference = rows[index].reference_index
            table_starts[index] = table_starts[reference]
            table_ends[index] = table_ends[reference]
            continue
        end
        table_starts[index] = length(currents_a) + 1
        append!(currents_a, (point.current_a for point in groups[index]))
        append!(fluxes_wb, (point.flux_wb for point in groups[index]))
        table_ends[index] = length(currents_a)
    end
    initial_segments = [
        table_starts[index] - 1 + _piecewise_nonlinear_inductor_initial_segment(
            rows[index].steady_state_current_a,
            groups[index],
        )
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
    return (
        source = :deck_piecewise_nonlinear_inductor_current_config,
        outcome = :nonlinear_current_config,
        nonlinear_types = fill(PIECEWISE_NONLINEAR_INDUCTOR_TYPE, count),
        nonlinear_from_nodes = shifted_from_nodes,
        nonlinear_to_nodes = shifted_to_nodes,
        nonlinear_admittance_nodes = table_starts,
        nonlinear_table_end_indices = table_ends,
        nonlinear_subsystem_indices = collect(1:count),
        subsystem_begin_indices = subsystem_begin_indices,
        subsystem_owner_rows = subsystem_owner_rows,
        subsystem_simultaneous_flags = subsystem_simultaneous_flags,
        cchar = currents_a,
        vchar = fluxes_wb,
        gslope = zeros(Float64, length(currents_a)),
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
        initial_characteristic_current_values = copy(initial_fluxes),
        initial_stored_voltage_values = copy(initial_fluxes),
        initial_runtime_voltage_values = copy(initial_fluxes),
        initial_current_segment_values = copy(initial_currents),
        initial_table_index_values = initial_segments,
        initial_cursub_values = zeros(Float64, max(required_node_count, count + 1)),
        nonlinear_current_segments = copy(initial_currents),
        nonlinear_output_codes = [row.output_code for row in rows],
        nonlinear_owner_names = [row.name for row in rows],
        nonlinear_owner_line_numbers = [row.line_no for row in rows],
        dense_primary_nonlinear_compensation = true,
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
        table_entry_count = length(currents_a),
        nonlinear_row_count = count,
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
    indices = findall(==(PIECEWISE_NONLINEAR_INDUCTOR_TYPE), types)
    names = Symbol.(config.nonlinear_owner_names)
    line_numbers = Int.(config.nonlinear_owner_line_numbers)
    output_codes = Int.(config.nonlinear_output_codes)
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
        for (position, nonlinear_index) in enumerate(indices)
            currents[position, step] =
                result.piecewise_nonlinear_inductor_accepted_currents[nonlinear_index]
            fluxes[position, step] =
                result.piecewise_nonlinear_inductor_accepted_fluxes[nonlinear_index]
            voltages[position, step] =
                result.piecewise_nonlinear_inductor_accepted_voltages[nonlinear_index]
            segments[position, step] = result.ilast[nonlinear_index]
        end
    end
    return (
        source = :piecewise_nonlinear_inductor_report,
        outcome = :report_output,
        names = names[indices],
        line_numbers = line_numbers[indices],
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

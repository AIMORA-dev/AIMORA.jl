export deck_triggered_timed_resistance_current_config

function _triggered_timed_resistance_point_groups(
    rows::Vector{DeckParser.DeckTriggeredTimedResistanceRow},
    points::Vector{DeckParser.DeckTriggeredTimedResistancePointRow},
)
    groups = [DeckParser.DeckTriggeredTimedResistancePointRow[] for _ in rows]
    for point in points
        1 <= point.resistance_row_index <= length(rows) ||
            throw(ArgumentError("timed-resistance point references an unknown row"))
        push!(groups[point.resistance_row_index], point)
    end
    for (row_index, row) in enumerate(rows)
        if row.source_kind == :copy_reference
            1 <= row.reference_index < row_index ||
                throw(ArgumentError("timed-resistance COPY rows must reference a prior owner"))
            isempty(groups[row_index]) ||
                throw(ArgumentError("timed-resistance COPY rows must not duplicate a schedule"))
            groups[row_index] = groups[row.reference_index]
        end
        group = groups[row_index]
        isempty(group) && throw(ArgumentError("timed resistance requires a schedule"))
        sort!(group; by = point -> point.point_index)
        for (index, point) in enumerate(group)
            point.point_index == index ||
                throw(ArgumentError("timed-resistance point indices must be consecutive"))
        end
    end
    return groups
end

function deck_triggered_timed_resistance_current_config(
    parsed::DeckParser.DeckParseResult;
    delta2::Real=1.0,
    fltinf::Real=1.0e99,
)
    DeckParser.assert_deck_valid!(parsed)
    rows = DeckParser.deck_triggered_timed_resistance_rows(parsed)
    isempty(rows) && throw(ArgumentError("deck has no triggered timed-resistance rows"))
    point_groups = _triggered_timed_resistance_point_groups(
        rows,
        DeckParser.deck_triggered_timed_resistance_point_rows(parsed),
    )
    delta = Float64(delta2)
    delta > 0.0 && isfinite(delta) ||
        throw(ArgumentError("delta2 must be finite and positive"))
    infinity = Float64(fltinf)
    infinity > 0.0 || throw(ArgumentError("fltinf must be positive"))

    elapsed_times_s = Float64[]
    resistances_ohm = Float64[]
    conductances_s = Float64[]
    table_start_indices = Int[]
    table_end_indices = Int[]
    for (row_index, points) in enumerate(point_groups)
        row = rows[row_index]
        if row.source_kind == :copy_reference
            push!(table_start_indices, table_start_indices[row.reference_index])
            push!(table_end_indices, table_end_indices[row.reference_index])
            continue
        end
        push!(table_start_indices, length(elapsed_times_s) + 1)
        previous_elapsed = -Inf
        for point in points
            point.elapsed_time_s >= 0.0 && point.elapsed_time_s > previous_elapsed ||
                throw(ArgumentError("timed-resistance elapsed times must increase"))
            point.resistance_ohm > 0.0 && isfinite(point.resistance_ohm) ||
                throw(ArgumentError("timed resistance must be finite and positive"))
            push!(elapsed_times_s, point.elapsed_time_s)
            push!(resistances_ohm, point.resistance_ohm)
            push!(conductances_s, inv(point.resistance_ohm))
            previous_elapsed = point.elapsed_time_s
        end
        push!(table_end_indices, length(elapsed_times_s))
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
    count = length(rows)
    last_slot = 1 + 5 * (count - 1)
    source_next_indices = zeros(Int, last_slot)
    source_from_nodes = zeros(Int, last_slot)
    source_to_nodes = zeros(Int, last_slot)
    subsystem_begin_indices = Int[]
    subsystem_owner_rows = zeros(Int, last_slot)
    subsystem_simultaneous_flags = zeros(Int, last_slot)
    for index in 1:count
        slot = 1 + 5 * (index - 1)
        source_next_indices[slot] = index == count ? 1 : slot + 5
        source_from_nodes[slot] = shifted_from_nodes[index]
        source_to_nodes[slot] = shifted_to_nodes[index]
        push!(subsystem_begin_indices, slot)
        subsystem_owner_rows[slot] = index
    end
    initial_segments = [row.arm_time_s < 0.0 ? 1.0 : 0.0 for row in rows]
    return (
        source = :deck_triggered_timed_resistance_current_config,
        outcome = :nonlinear_current_config,
        nonlinear_types = fill(TRIGGERED_TIMED_RESISTANCE_TYPE, count),
        nonlinear_from_nodes = shifted_from_nodes,
        nonlinear_to_nodes = shifted_to_nodes,
        nonlinear_admittance_nodes = table_start_indices,
        nonlinear_table_end_indices = table_end_indices,
        nonlinear_subsystem_indices = collect(1:count),
        subsystem_begin_indices = subsystem_begin_indices,
        subsystem_owner_rows = subsystem_owner_rows,
        subsystem_simultaneous_flags = subsystem_simultaneous_flags,
        cchar = elapsed_times_s,
        vchar = resistances_ohm,
        gslope = conductances_s,
        delta2 = delta,
        deltat = 2.0 * delta,
        fltinf = infinity,
        ncomp = count,
        deck_owned_voltage_context = true,
        reference_shifted_voltage_context = true,
        deck_to_runtime_node_indices = deck_to_runtime_node_indices,
        reference_node_index = 1,
        nonlinear_required_node_count = required_node_count,
        initialize_nonlinear_state = true,
        initial_companion_current_values = zeros(Float64, count),
        initial_characteristic_current_values = zeros(Float64, count),
        initial_stored_voltage_values = [row.trigger_voltage_v for row in rows],
        initial_runtime_voltage_values = [row.trigger_voltage_v for row in rows],
        initial_current_segment_values = initial_segments,
        initial_table_index_values = copy(table_start_indices),
        initial_cursub_values = zeros(Float64, max(required_node_count, count + 1)),
        nonlinear_current_segments = copy(initial_segments),
        timed_resistance_arm_time_values = [row.arm_time_s for row in rows],
        nonlinear_output_codes = [row.output_code for row in rows],
        nonlinear_owner_names = [row.name for row in rows],
        nonlinear_owner_line_numbers = [row.line_no for row in rows],
        nonlinear_reference_indices = [row.reference_index for row in rows],
        nonlinear_source_kinds = [row.source_kind for row in rows],
        dense_primary_nonlinear_compensation = true,
        use_state_nonlinear_inverse_columns = true,
        nonlinear_inverse_config = (
            ntot = required_node_count,
            ncomp = count,
            source_begin_indices = [1],
            source_next_indices = source_next_indices,
            source_from_nodes = source_from_nodes,
            source_to_nodes = source_to_nodes,
            source_activity_flags = source_activity_flags,
            nonlinear_types = fill(TRIGGERED_TIMED_RESISTANCE_TYPE, count),
            nonlinear_admittance_nodes = table_start_indices,
            nonlinear_from_nodes = shifted_from_nodes,
            nonlinear_to_nodes = shifted_to_nodes,
            nonlinear_source_flags = zeros(Int, count),
            partition_boundary = max(required_node_count - 1, 1),
            delta2 = delta,
            fltinf = infinity,
        ),
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
        table_entry_count = length(elapsed_times_s),
        nonlinear_row_count = count,
        mutation_order = (
            :timed_resistance_header,
            :resistance_schedule,
            :voltage_trigger,
            :conductance_restamp,
            :corrected_solve,
            :current_output,
        ),
        deferred_calls = (),
        complete_nonlinear_source_loop = true,
        replacement_ready = true,
    )
end

function _triggered_timed_resistance_report(run, config::NamedTuple)
    types = Int.(get(config, :nonlinear_types, Int[]))
    indices = findall(==(TRIGGERED_TIMED_RESISTANCE_TYPE), types)
    names = Symbol.(get(config, :nonlinear_owner_names, fill(Symbol(""), length(types))))
    line_numbers = Int.(get(config, :nonlinear_owner_line_numbers, zeros(Int, length(types))))
    output_codes = Int.(get(config, :nonlinear_output_codes, zeros(Int, length(types))))
    from_nodes = Int.(get(config, :nonlinear_from_nodes, Int[]))
    to_nodes = abs.(Int.(get(config, :nonlinear_to_nodes, Int[])))
    table_starts = Int.(get(config, :nonlinear_admittance_nodes, Int[]))
    nodal_to_reference = Int.(get(config, :deck_to_runtime_node_indices, Int[]))
    times = Float64[]
    currents = zeros(Float64, length(indices), length(run.over16_updates))
    voltages = similar(currents)
    segments = zeros(Int, size(currents))
    conductances = similar(currents)
    for (step, update) in enumerate(run.over16_updates)
        push!(times, Float64(update.t_s))
        result = update.over16_update.nonlinear_current_result
        result === nothing && throw(ArgumentError(
            "timed-resistance report requires a nonlinear timestep result",
        ))
        reference_voltage = zeros(Float64, Int(config.nonlinear_required_node_count))
        for nodal_index in eachindex(nodal_to_reference)
            reference_voltage[nodal_to_reference[nodal_index]] = update.voltage_pu[nodal_index]
        end
        for (row_position, nonlinear_index) in enumerate(indices)
            voltage = reference_voltage[from_nodes[nonlinear_index]] -
                      reference_voltage[to_nodes[nonlinear_index]]
            segment = round(Int, result.curr[nonlinear_index])
            conductance = segment == 0 ? 0.0 :
                result.gslope[table_starts[nonlinear_index] + segment - 1]
            voltages[row_position, step] = voltage
            conductances[row_position, step] = conductance
            currents[row_position, step] = conductance * voltage
            segments[row_position, step] = segment
        end
    end
    return (
        source = :triggered_timed_resistance_report,
        outcome = :report_output,
        names = names[indices],
        line_numbers = line_numbers[indices],
        output_codes = output_codes[indices],
        output_requested_flags = output_codes[indices] .> 0,
        time_s = times,
        voltage_v = voltages,
        current_a = currents,
        conductance_s = conductances,
        active_segments = segments,
        deferred_effects = (),
    )
end

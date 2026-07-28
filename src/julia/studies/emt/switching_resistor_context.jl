export SwitchingNonlinearResistorSafetyShunt,
       deck_switching_nonlinear_resistor_current_config,
       deck_switching_nonlinear_resistor_safety_shunts

struct SaturatedTransformerSafetyShunt
    owner_name::Symbol
    from_node_index::Int
    to_node_index::Int
    resistance_ohm::Float64
    conductance_s::Float64
end

struct SwitchingNonlinearResistorSafetyShunt
    resistor_name::Symbol
    from_node::Symbol
    to_node::Symbol
    from_node_index::Int
    to_node_index::Int
    resistance_ohm::Float64
    conductance_s::Float64
end

function _switching_resistor_parallel_branch_exists(
    parsed::DeckParser.DeckParseResult,
    from_node::Int,
    to_node::Int,
)
    target = minmax(from_node, to_node)
    for element in parsed.elements
        hasproperty(element, :a) && hasproperty(element, :b) || continue
        a = getproperty(element, :a)
        b = getproperty(element, :b)
        a isa Integer && b isa Integer || continue
        minmax(Int(a), Int(b)) == target && return true
    end
    return false
end

function _switching_resistor_isolation_resistance_ohm(
    parsed::DeckParser.DeckParseResult,
    numerical_zero::Float64,
)
    numerical_zero > 0.0 && isfinite(numerical_zero) ||
        throw(ArgumentError("switching resistor numerical zero must be finite and positive"))
    resistance = inv(100.0 * numerical_zero)
    for request in DeckParser.deck_study_option_request_rows(parsed)
        request.request_kind == :high_resistance_exponent || continue
        isempty(request.numeric_values) || (resistance = Float64(first(request.numeric_values)))
    end
    resistance > 0.0 && isfinite(resistance) ||
        throw(ArgumentError("switching resistor isolation resistance must be finite and positive"))
    return resistance
end

function deck_switching_nonlinear_resistor_safety_shunts(
    parsed::DeckParser.DeckParseResult;
    numerical_zero::Real=1.0e-12,
)
    DeckParser.assert_deck_valid!(parsed)
    isempty(DeckParser.deck_synchronous_machine_terminal_voltage_rows(parsed)) &&
        return SwitchingNonlinearResistorSafetyShunt[]
    resistance = _switching_resistor_isolation_resistance_ohm(
        parsed,
        Float64(numerical_zero),
    )
    shunts = SwitchingNonlinearResistorSafetyShunt[]
    for row in DeckParser.deck_switching_nonlinear_resistor_rows(parsed)
        _switching_resistor_parallel_branch_exists(
            parsed,
            row.from_node_index,
            row.to_node_index,
        ) && continue
        push!(
            shunts,
            SwitchingNonlinearResistorSafetyShunt(
                row.name,
                row.from_node,
                row.to_node,
                row.from_node_index,
                row.to_node_index,
                resistance,
                inv(resistance),
            ),
        )
    end
    return shunts
end

function _append_switching_nonlinear_resistor_safety_shunts!(
    elements::Vector,
    element_names::Vector{Symbol},
    parsed::DeckParser.DeckParseResult,
)
    shunts = deck_switching_nonlinear_resistor_safety_shunts(parsed)
    for shunt in shunts
        push!(
            elements,
            ConductanceBranch(
                shunt.from_node_index,
                shunt.to_node_index,
                shunt.conductance_s,
            ),
        )
        push!(element_names, Symbol(shunt.resistor_name, "_isolation_shunt"))
    end
    return shunts
end

function _saturated_transformer_parallel_branch_exists(
    branch_assembly,
    from_node::Int,
    to_node::Int,
)
    target = minmax(from_node, to_node)
    for (from, to) in zip(
        branch_assembly.magnetizing_branch_from_node_indices,
        branch_assembly.magnetizing_branch_to_node_indices,
    )
        minmax(Int(from), Int(to)) == target && return true
    end
    return false
end

function _saturated_transformer_safety_shunts(
    parsed::DeckParser.DeckParseResult,
    current_config::NamedTuple,
)
    isempty(DeckParser.deck_synchronous_machine_terminal_voltage_rows(parsed)) &&
        return SaturatedTransformerSafetyShunt[]
    resistance = _switching_resistor_isolation_resistance_ohm(parsed, 1.0e-12)
    types = Int.(get(current_config, :nonlinear_types, Int[]))
    from_nodes = Int.(get(current_config, :nonlinear_from_nodes, Int[]))
    to_nodes = abs.(Int.(get(current_config, :nonlinear_to_nodes, Int[])))
    length(types) == length(from_nodes) == length(to_nodes) ||
        throw(ArgumentError("saturated transformer safety-shunt arrays must align"))
    branch_assembly = current_config.saturated_transformer_branch_assembly
    shunts = SaturatedTransformerSafetyShunt[]
    for index in eachindex(types)
        types[index] == SATURATED_TRANSFORMER_NONLINEAR_TYPE || continue
        _saturated_transformer_parallel_branch_exists(
            branch_assembly,
            from_nodes[index],
            to_nodes[index],
        ) && continue
        push!(
            shunts,
            SaturatedTransformerSafetyShunt(
                Symbol("saturated_transformer_", index),
                from_nodes[index],
                to_nodes[index],
                resistance,
                inv(resistance),
            ),
        )
    end
    return shunts
end

function _append_saturated_transformer_safety_shunts!(
    elements::Vector,
    element_names::Vector{Symbol},
    parsed::DeckParser.DeckParseResult,
    current_config::NamedTuple,
)
    shunts = _saturated_transformer_safety_shunts(parsed, current_config)
    for shunt in shunts
        push!(
            elements,
            ConductanceBranch(
                shunt.from_node_index,
                shunt.to_node_index,
                shunt.conductance_s,
            ),
        )
        push!(element_names, Symbol(shunt.owner_name, "_isolation_shunt"))
    end
    return shunts
end

function _switching_nonlinear_resistor_point_groups(
    rows::Vector{DeckParser.DeckSwitchingNonlinearResistorRow},
    points::Vector{DeckParser.DeckSwitchingNonlinearResistorPointRow},
)
    groups = [DeckParser.DeckSwitchingNonlinearResistorPointRow[] for _ in rows]
    for point in points
        1 <= point.resistor_row_index <= length(rows) ||
            throw(ArgumentError("switching resistor point references an unknown row"))
        push!(groups[point.resistor_row_index], point)
    end
    for (index, row) in enumerate(rows)
        if row.source_kind == :copy_reference
            1 <= row.reference_index < index ||
                throw(ArgumentError("switching resistor COPY row must reference a prior owner"))
            isempty(groups[index]) ||
                throw(ArgumentError("switching resistor COPY row must not duplicate a V-I table"))
            groups[index] = groups[row.reference_index]
        end
        group = groups[index]
        isempty(group) && throw(ArgumentError("switching resistor requires a V-I table"))
        sort!(group; by = point -> point.point_index)
        for (index, point) in enumerate(group)
            point.point_index == index ||
                throw(ArgumentError("switching resistor point indices must be consecutive"))
        end
    end
    return groups
end

function deck_switching_nonlinear_resistor_current_config(
    parsed::DeckParser.DeckParseResult;
    delta2::Real=1.0,
    fltinf::Real=1.0e99,
    flzero::Real=1.0e-12,
)
    DeckParser.assert_deck_valid!(parsed)
    rows = DeckParser.deck_switching_nonlinear_resistor_rows(parsed)
    isempty(rows) && throw(ArgumentError("deck has no switching nonlinear resistor rows"))
    point_groups = _switching_nonlinear_resistor_point_groups(
        rows,
        DeckParser.deck_switching_nonlinear_resistor_point_rows(parsed),
    )
    delta = Float64(delta2)
    delta > 0.0 && isfinite(delta) ||
        throw(ArgumentError("delta2 must be finite and positive"))
    infinity = Float64(fltinf)
    infinity > 0.0 || throw(ArgumentError("fltinf must be positive"))
    zero_tolerance = Float64(flzero)
    zero_tolerance >= 0.0 && isfinite(zero_tolerance) ||
        throw(ArgumentError("flzero must be finite and nonnegative"))

    cchar = Float64[]
    gslope = Float64[]
    vchar = Float64[]
    table_start_indices = Int[]
    table_end_indices = Int[]
    for (index, points) in enumerate(point_groups)
        row = rows[index]
        if row.source_kind == :copy_reference
            push!(table_start_indices, table_start_indices[row.reference_index])
            push!(table_end_indices, table_end_indices[row.reference_index])
            continue
        end
        push!(table_start_indices, length(cchar) + 1)
        previous_current = 0.0
        previous_voltage = 0.0
        for point in points
            voltage_delta = point.voltage_v - previous_voltage
            voltage_delta > 0.0 ||
                throw(ArgumentError("switching resistor voltages must increase"))
            slope = (point.current_a - previous_current) / voltage_delta
            slope >= 0.0 && isfinite(slope) ||
                throw(ArgumentError("switching resistor segment conductance must be finite and nonnegative"))
            push!(gslope, slope)
            push!(cchar, previous_current - slope * previous_voltage)
            push!(vchar, point.voltage_v)
            previous_current = point.current_a
            previous_voltage = point.voltage_v
        end
        push!(table_end_indices, length(cchar))
    end
    rearm_time_state_indices = Int[]
    for _ in rows
        push!(cchar, 0.0)
        push!(gslope, 0.0)
        push!(vchar, 0.0)
        push!(rearm_time_state_indices, length(vchar))
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
    safety_shunts = deck_switching_nonlinear_resistor_safety_shunts(
        parsed;
        numerical_zero = max(zero_tolerance, 1.0e-12),
    )
    characteristic_owner_indices = collect(eachindex(rows))
    for (index, row) in enumerate(rows)
        row.source_kind == :copy_reference || continue
        characteristic_owner_indices[index] =
            characteristic_owner_indices[row.reference_index]
    end
    return (
        source = :deck_switching_nonlinear_resistor_current_config,
        outcome = :nonlinear_current_config,
        nonlinear_types = fill(SWITCHING_NONLINEAR_RESISTOR_TYPE, count),
        nonlinear_from_nodes = shifted_from_nodes,
        nonlinear_to_nodes = shifted_to_nodes,
        nonlinear_admittance_nodes = table_start_indices,
        nonlinear_table_end_indices = table_end_indices,
        nonlinear_subsystem_indices = collect(1:count),
        subsystem_begin_indices = subsystem_begin_indices,
        subsystem_owner_rows = subsystem_owner_rows,
        subsystem_simultaneous_flags = subsystem_simultaneous_flags,
        cchar = cchar,
        vchar = vchar,
        gslope = gslope,
        delta2 = delta,
        deltat = 2.0 * delta,
        fltinf = infinity,
        flzero = zero_tolerance,
        ncomp = count,
        deck_owned_voltage_context = true,
        reference_shifted_voltage_context = true,
        deck_to_runtime_node_indices = deck_to_runtime_node_indices,
        reference_node_index = 1,
        nonlinear_required_node_count = required_node_count,
        initialize_nonlinear_state = true,
        initial_companion_current_values = [
            Float64(row.activation_segment_count) for row in rows
        ],
        initial_characteristic_current_values = [row.turn_off_voltage for row in rows],
        initial_stored_voltage_values = [row.turn_on_voltage for row in rows],
        initial_runtime_voltage_values = [row.turn_on_voltage for row in rows],
        initial_current_segment_values = zeros(Float64, count),
        initial_table_index_values = copy(table_start_indices),
        initial_cursub_values = zeros(Float64, max(required_node_count, count + 1)),
        nonlinear_current_segments = zeros(Float64, count),
        minimum_on_time_values = [row.minimum_on_time_s for row in rows],
        rearm_time_state_indices = rearm_time_state_indices,
        single_flash_flags = [row.single_flash for row in rows],
        nonlinear_output_codes = [row.output_code for row in rows],
        nonlinear_owner_names = [row.name for row in rows],
        nonlinear_owner_line_numbers = [row.line_no for row in rows],
        nonlinear_characteristic_owner_indices = characteristic_owner_indices,
        safety_shunt_resistor_names = [shunt.resistor_name for shunt in safety_shunts],
        safety_shunt_resistance_ohm_values = [shunt.resistance_ohm for shunt in safety_shunts],
        safety_shunt_conductance_s_values = [shunt.conductance_s for shunt in safety_shunts],
        dense_primary_nonlinear_compensation = true,
        use_state_nonlinear_inverse_columns = true,
        nonlinear_inverse_config = (
            require_fortran_sparse_factor_result = true,
            ntot = required_node_count,
            ncomp = count,
            source_begin_indices = [1],
            source_next_indices = source_next_indices,
            source_from_nodes = source_from_nodes,
            source_to_nodes = source_to_nodes,
            source_activity_flags = source_activity_flags,
            nonlinear_types = fill(SWITCHING_NONLINEAR_RESISTOR_TYPE, count),
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
        table_entry_count = length(cchar),
        nonlinear_row_count = count,
        mutation_order = (
            :switching_resistor_header,
            :piecewise_vi_table,
            :machine_triggered_isolation_shunt,
            :segment_and_admittance_update,
            :companion_current_rhs,
            :current_output,
        ),
        deferred_calls = (),
    )
end

function _switching_nonlinear_resistor_report(run, config::NamedTuple)
    types = Int.(get(config, :nonlinear_types, Int[]))
    indices = findall(==(SWITCHING_NONLINEAR_RESISTOR_TYPE), types)
    names = Symbol.(get(config, :nonlinear_owner_names, fill(Symbol(""), length(types))))
    line_numbers = Int.(get(config, :nonlinear_owner_line_numbers, zeros(Int, length(types))))
    output_codes = Int.(get(config, :nonlinear_output_codes, zeros(Int, length(types))))
    from_nodes = Int.(get(config, :nonlinear_from_nodes, Int[]))
    to_nodes = abs.(Int.(get(config, :nonlinear_to_nodes, Int[])))
    table_starts = Int.(get(config, :nonlinear_admittance_nodes, Int[]))
    subsystem_indices = Int.(get(config, :nonlinear_subsystem_indices, Int[]))
    subsystem_heads = Int.(get(config, :subsystem_begin_indices, Int[]))
    nodal_to_reference = Int.(get(config, :deck_to_runtime_node_indices, Int[]))
    times = Float64[]
    currents = zeros(Float64, length(indices), length(run.over16_updates))
    voltages = similar(currents)
    segments = zeros(Int, size(currents))
    for (step, update) in enumerate(run.over16_updates)
        push!(times, Float64(update.t_s))
        result = update.over16_update.nonlinear_current_result
        result === nothing && throw(ArgumentError(
            "switching nonlinear resistor report requires a nonlinear timestep result",
        ))
        reference_voltage = zeros(Float64, Int(config.nonlinear_required_node_count))
        for nodal_index in eachindex(nodal_to_reference)
            reference_voltage[nodal_to_reference[nodal_index]] = update.voltage_pu[nodal_index]
        end
        for (row_position, nonlinear_index) in enumerate(indices)
            voltage = reference_voltage[from_nodes[nonlinear_index]] -
                      reference_voltage[to_nodes[nonlinear_index]]
            segment = round(Int, result.curr[nonlinear_index])
            subsystem = subsystem_indices[nonlinear_index]
            head = subsystem_heads[subsystem]
            companion_index = div(head, 5) + 1
            companion = result.cursub[companion_index]
            current = if segment == 0
                0.0
            else
                table_index = table_starts[nonlinear_index] + abs(segment) - 1
                result.gslope[table_index] * voltage + companion
            end
            voltages[row_position, step] = voltage
            currents[row_position, step] = current
            segments[row_position, step] = segment
        end
    end
    return (
        source = :switching_nonlinear_resistor_report,
        outcome = :report_output,
        names = names[indices],
        line_numbers = line_numbers[indices],
        output_codes = output_codes[indices],
        output_requested_flags = output_codes[indices] .> 0,
        time_s = times,
        voltage_v = voltages,
        current_a = currents,
        active_segments = segments,
        deferred_effects = (),
    )
end

function _nonlinear_isolation_shunt_report(config::NamedTuple)
    return (
        source = :nonlinear_isolation_shunt_report,
        outcome = :report_output,
        owner_names = Symbol.(get(config, :safety_shunt_resistor_names, Symbol[])),
        resistance_ohm =
            Float64.(get(config, :safety_shunt_resistance_ohm_values, Float64[])),
        conductance_s =
            Float64.(get(config, :safety_shunt_conductance_s_values, Float64[])),
        deferred_effects = (),
    )
end

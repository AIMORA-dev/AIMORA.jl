function _terminal_node_name(node_names::AbstractVector{Symbol}, index::Int)
    index == 0 && return :ground
    return 1 <= index <= length(node_names) ? node_names[index] : Symbol("node_", index)
end

function _terminal_branch_states(context::EMTStepContext)
    rows = EMTTerminalBranchState[]
    voltage = context.system.v
    for (index, element) in enumerate(context.system.elements)
        snapshot = branch_companion_snapshot(element, voltage, context.dt_s)
        snapshot === nothing && continue
        accepted_current = branch_current_value(element, voltage, context.dt_s)
        accepted_history =
            accepted_current - snapshot.conductance * snapshot.branch_voltage
        name = index <= length(context.element_names) ?
            context.element_names[index] : Symbol("branch_", index)
        energy = index in context.branch_power_output_branch_indices ?
            context.branch_energy_values[index] : NaN
        resistance_ohm = element isa SeriesRLCBranch ? element.r : NaN
        inductance_h = element isa SeriesRLCBranch ? element.l : NaN
        capacitance_f = element isa SeriesRLCBranch ? element.c : NaN
        push!(
            rows,
            EMTTerminalBranchState(
                name,
                snapshot.kind,
                _terminal_node_name(context.node_names, snapshot.a),
                _terminal_node_name(context.node_names, snapshot.b),
                snapshot.conductance,
                accepted_history,
                snapshot.branch_voltage,
                accepted_current,
                snapshot.previous_current,
                snapshot.previous_voltage,
                energy,
                resistance_ohm,
                inductance_h,
                capacitance_f,
            ),
        )
    end
    return rows
end

function _terminal_switch_states(context::EMTStepContext)
    count = context.deck_time_switch_count
    rows = Vector{EMTTerminalSwitchState}(undef, count)
    for index in 1:count
        rows[index] = EMTTerminalSwitchState(
            context.deck_time_switch_names[index],
            _terminal_node_name(
                context.node_names,
                context.deck_time_switch_from_node_indices[index],
            ),
            _terminal_node_name(
                context.node_names,
                context.deck_time_switch_to_node_indices[index],
            ),
            context.switch_closed_step_flags[index] != 0,
            context.switch_conductance_step_values[index],
            context.switch_voltage_step_values[index],
            context.switch_current_step_values[index],
            context.switch_power_step_values[index],
            index <= length(context.deck_over5_switch_output_codes) &&
            context.deck_over5_switch_output_codes[index] > 3 ?
                context.switch_energy_values[index] : NaN,
            context.deck_time_switch_close_time_s_values[index],
            context.deck_time_switch_open_time_s_values[index],
        )
    end
    return rows
end

function _nonlinear_kind(type_code::Int)
    type_code == SWITCHING_NONLINEAR_RESISTOR_TYPE &&
        return :switching_nonlinear_resistor
    type_code == TRIGGERED_TIMED_RESISTANCE_TYPE &&
        return :triggered_timed_resistance
    type_code == PIECEWISE_NONLINEAR_INDUCTOR_TYPE &&
        return :piecewise_nonlinear_inductor
    type_code == HYSTERETIC_INDUCTOR_NONLINEAR_TYPE &&
        return :hysteretic_inductor
    type_code == SATURATED_TRANSFORMER_NONLINEAR_TYPE &&
        return :saturated_transformer
    type_code in (94, -97) && return :surge_arrester
    return :nonlinear_element
end

function _terminal_nonlinear_owner_rows(parsed::DeckParser.DeckParseResult)
    rows = Dict{Symbol,NamedTuple{(:from_node, :to_node),Tuple{Symbol,Symbol}}}()
    families = (
        DeckParser.deck_nonlinear_resistance_rows(parsed),
        DeckParser.deck_triggered_timed_resistance_rows(parsed),
        DeckParser.deck_switching_nonlinear_resistor_rows(parsed),
        DeckParser.deck_piecewise_nonlinear_inductor_rows(parsed),
        DeckParser.deck_hysteretic_inductor_rows(parsed),
        DeckParser.deck_arrester_nonlinear_rows(parsed),
    )
    for family in families, row in family
        rows[row.name] = (from_node = row.from_node, to_node = row.to_node)
    end
    return rows
end

function _terminal_state_complete(
    parsed::DeckParser.DeckParseResult,
    branches,
    nonlinear_elements,
    switches,
)
    expected_branch_names = Set{Symbol}()
    for (index, element) in enumerate(parsed.elements)
        element isa Union{
            ConductanceBranch,
            SeriesRLBranch,
            SeriesRLCBranch,
            CoupledInductiveBranch,
            CoupledSeriesRLBranch,
            CapacitorBranch,
        } || continue
        push!(expected_branch_names, parsed.element_names[index])
    end
    actual_branch_names = Set(row.name for row in branches)
    issubset(expected_branch_names, actual_branch_names) || return false

    expected_switch_names = Set(DeckParser.deck_time_switch_names(parsed))
    actual_switch_names = Set(row.name for row in switches)
    expected_switch_names == actual_switch_names || return false

    expected_nonlinear_names = Set(keys(_terminal_nonlinear_owner_rows(parsed)))
    actual_nonlinear_names = Set(row.name for row in nonlinear_elements)
    union!(actual_nonlinear_names, actual_branch_names)
    issubset(expected_nonlinear_names, actual_nonlinear_names) || return false
    return true
end

function _terminal_energy_j(
    times_s::AbstractVector{<:Real},
    voltages_v::AbstractVector{<:Real},
    currents_a::AbstractVector{<:Real},
)
    length(times_s) == length(voltages_v) == length(currents_a) ||
        throw(ArgumentError("nonlinear terminal-state histories must have equal lengths"))
    energy = 0.0
    @inbounds for index in 2:length(times_s)
        dt = Float64(times_s[index] - times_s[index - 1])
        previous_power = Float64(voltages_v[index - 1] * currents_a[index - 1])
        current_power = Float64(voltages_v[index] * currents_a[index])
        energy += 0.5 * dt * (previous_power + current_power)
    end
    return energy
end

function _append_terminal_nonlinear_report!(
    rows::Vector{EMTTerminalNonlinearState},
    seen::Set{Symbol},
    report,
    type_codes::AbstractVector{<:Integer},
    owners,
)
    isempty(report.names) && return rows
    size(report.voltage_v) == size(report.current_a) ||
        throw(ArgumentError("nonlinear terminal-state voltage/current shapes must match"))
    size(report.voltage_v, 1) == length(report.names) ||
        throw(ArgumentError("nonlinear terminal-state rows must match their names"))
    length(type_codes) == length(report.names) ||
        throw(ArgumentError("nonlinear terminal-state type codes must match their names"))
    for index in eachindex(report.names)
        name = report.names[index]
        name in seen && continue
        owner = get(owners, name, (from_node = :unknown, to_node = :unknown))
        voltages = view(report.voltage_v, index, :)
        currents = view(report.current_a, index, :)
        flux = hasproperty(report, :flux_wb) ?
            Float64(report.flux_wb[index, end]) : NaN
        segment = hasproperty(report, :active_segments) ?
            Int(report.active_segments[index, end]) : 0
        push!(
            rows,
            EMTTerminalNonlinearState(
                name,
                _nonlinear_kind(Int(type_codes[index])),
                owner.from_node,
                owner.to_node,
                Float64(last(voltages)),
                Float64(last(currents)),
                flux,
                segment,
                _terminal_energy_j(report.time_s, voltages, currents),
            ),
        )
        push!(seen, name)
    end
    return rows
end

function _terminal_nonlinear_states(parsed::DeckParser.DeckParseResult, run)
    owners = _terminal_nonlinear_owner_rows(parsed)
    rows = EMTTerminalNonlinearState[]
    seen = Set{Symbol}()
    reports = (
        (:switching_nonlinear_resistor_report, SWITCHING_NONLINEAR_RESISTOR_TYPE),
        (:triggered_timed_resistance_report, TRIGGERED_TIMED_RESISTANCE_TYPE),
        (:piecewise_nonlinear_inductor_report, PIECEWISE_NONLINEAR_INDUCTOR_TYPE),
    )
    for (property, type_code) in reports
        hasproperty(run, property) || continue
        report = getproperty(run, property)
        _append_terminal_nonlinear_report!(
            rows,
            seen,
            report,
            fill(type_code, length(report.names)),
            owners,
        )
    end
    if hasproperty(run, :nonlinear_branch_report)
        report = run.nonlinear_branch_report
        _append_terminal_nonlinear_report!(
            rows,
            seen,
            report,
            report.nonlinear_types,
            owners,
        )
    end
    return rows
end

function _terminal_state_checks(state::EMTTerminalState)
    isfinite(state.time_s) || return false
    all(row -> isfinite(row.voltage_v), state.nodes) || return false
    all(
        row ->
            all(isfinite, (
                row.conductance_s,
                row.history_current_a,
                row.voltage_v,
                row.current_a,
                row.previous_current_a,
                row.previous_voltage_v,
            )) &&
            (isfinite(row.energy_j) || isnan(row.energy_j)) &&
            (
                row.kind == :series_rlc ?
                (
                    isfinite(row.resistance_ohm) &&
                    row.resistance_ohm >= 0.0 &&
                    isfinite(row.inductance_h) &&
                    row.inductance_h >= 0.0 &&
                    isfinite(row.capacitance_f) &&
                    row.capacitance_f > 0.0
                ) :
                (
                    isnan(row.resistance_ohm) &&
                    isnan(row.inductance_h) &&
                    isnan(row.capacitance_f)
                )
            ) &&
            isapprox(
                row.current_a,
                row.conductance_s * row.voltage_v + row.history_current_a;
                atol = 1.0e-10,
                rtol = 1.0e-10,
            ),
        state.branches,
    ) || return false
    all(
        row ->
            all(isfinite, (
                row.conductance_s,
                row.voltage_v,
                row.current_a,
                row.power_w,
            )) &&
            (isfinite(row.energy_j) || isnan(row.energy_j)) &&
            isapprox(row.power_w, row.voltage_v * row.current_a; atol = 1.0e-8, rtol = 1.0e-10),
        state.switches,
    ) || return false
    all(
        row ->
            all(isfinite, (row.voltage_v, row.current_a, row.energy_j)) &&
            (isfinite(row.flux_wb) || isnan(row.flux_wb)),
        state.nonlinear_elements,
    ) || return false
    return true
end

function electromagnetic_terminal_state(
    parsed::DeckParser.DeckParseResult,
    trace::DeckEMTTrace;
    context::Union{Nothing,EMTStepContext}=nothing,
    nonlinear_run=nothing,
)
    node_names = trace.node_names
    node_values = view(trace.voltage_pu, :, size(trace.voltage_pu, 2))
    nodes = EMTTerminalNodeState[
        EMTTerminalNodeState(node_names[index], Float64(node_values[index]))
        for index in eachindex(node_names)
    ]
    branches = context === nothing ? EMTTerminalBranchState[] :
        _terminal_branch_states(context)
    switches = context === nothing ? EMTTerminalSwitchState[] :
        _terminal_switch_states(context)
    nonlinear_elements = nonlinear_run === nothing ?
        EMTTerminalNonlinearState[] :
        _terminal_nonlinear_states(parsed, nonlinear_run)
    provisional = EMTTerminalState(
        parsed.source,
        Float64(last(trace.time_s)),
        nodes,
        branches,
        nonlinear_elements,
        switches,
        false,
    )
    return EMTTerminalState(
        provisional.source,
        provisional.time_s,
        provisional.nodes,
        provisional.branches,
        provisional.nonlinear_elements,
        provisional.switches,
        _terminal_state_checks(provisional) &&
        _terminal_state_complete(
            parsed,
            provisional.branches,
            provisional.nonlinear_elements,
            provisional.switches,
        ),
    )
end

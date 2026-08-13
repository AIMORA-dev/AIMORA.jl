@enum PowerSemiconductorTopologyFault begin
    BRIDGE_POSITION_HEALTHY
    BRIDGE_POSITION_STUCK_OPEN
    BRIDGE_POSITION_STUCK_CLOSED
end

"""Accepted aggregate terminal, valve, passive, fault, KCL, power, and energy state."""
struct PowerSemiconductorBridgeTopologyState
    schema::Symbol
    family::Symbol
    topology_signature::String
    terminal_names::Vector{Symbol}
    terminal_voltage_v::Vector{Float64}
    terminal_current_a::Vector{Float64}
    valve_names::Vector{Symbol}
    requested_gate_state::BitVector
    applied_gate_state::BitVector
    conducting_state::BitVector
    passive_names::Vector{Symbol}
    passive_voltage_v::Vector{Float64}
    passive_current_a::Vector{Float64}
    passive_charge_c::Vector{Float64}
    passive_flux_wb::Vector{Float64}
    terminal_kcl_residual_a::Float64
    terminal_power_w::Float64
    semiconductor_loss_w::Float64
    passive_dissipated_power_w::Float64
    stored_energy_j::Float64
    dissipated_energy_j::Float64
    blocked::Bool
    position_faults::Vector{PowerSemiconductorTopologyFault}
    transition_count::Int
    refusal_count::Int
    deterministic_signature::String
end

"""One executable aggregate that composes canonical switch, bridge-leg, and passive owners."""
mutable struct PowerSemiconductorBridgeTopology{T,V,P,L} <: EMTElement
    topology::T
    valves::V
    passives::P
    bridge_legs::L
    blocked::Bool
    position_faults::Vector{PowerSemiconductorTopologyFault}
    transition_count::Int
    refusal_count::Int
    last_terminal_voltage_v::Vector{Float64}
    last_terminal_current_a::Vector{Float64}
    last_step_s::Float64
    dissipated_energy_j::Float64
end

function _bridge_topology_position_index(
    topology::BridgeTopologyDescriptor,
    name::Symbol,
)
    index = findfirst(position -> position.name === name, topology.valve_positions)
    index === nothing && throw(ArgumentError("bridge topology has no valve position $name"))
    return index
end

function _validate_bridge_topology_valve(position::BridgeValvePosition, valve::PowerSemiconductorSwitch)
    valve.a == position.from_node && valve.b == position.to_node || throw(ArgumentError(
        "bridge valve $(position.name) terminal orientation does not match its topology contract",
    ))
    kind = power_semiconductor_kind(valve)
    compatible = if position.valve_class === :diode
        kind === :diode
    elseif position.valve_class === :thyristor
        kind in (:thyristor, :gto)
    else
        kind in (:igbt, :mosfet)
    end
    compatible || throw(ArgumentError(
        "bridge valve $(position.name) kind $kind is incompatible with $(position.valve_class)",
    ))
    return valve
end

function _validate_bridge_topology_passive(position::BridgePassivePosition, passive::EMTElement)
    hasproperty(passive, :a) && hasproperty(passive, :b) || throw(ArgumentError(
        "bridge passive $(position.name) must expose oriented scalar terminals",
    ))
    passive.a == position.from_node && passive.b == position.to_node || throw(ArgumentError(
        "bridge passive $(position.name) terminal orientation does not match its topology contract",
    ))
    kind = passive isa ConductanceBranch ? :conductance :
        passive isa SeriesRLBranch ? :series_rl :
        passive isa SeriesRLCBranch ? :series_rlc :
        passive isa CapacitorBranch ? :capacitor : :unsupported
    kind === position.kind || throw(ArgumentError(
        "bridge passive $(position.name) owner does not match declared kind $(position.kind)",
    ))
    return passive
end

function PowerSemiconductorBridgeTopology(
    topology::BridgeTopologyDescriptor,
    valves::AbstractVector{<:PowerSemiconductorSwitch};
    passives::AbstractVector{<:EMTElement}=EMTElement[],
    bridge_legs::AbstractVector{<:PowerSemiconductorBridgeLeg}=PowerSemiconductorBridgeLeg[],
    initially_blocked::Bool=false,
)
    bridge_topology_signature(topology)
    length(valves) == length(topology.valve_positions) || throw(DimensionMismatch(
        "bridge topology requires one canonical switch owner per valve position",
    ))
    length(passives) == length(topology.passive_positions) || throw(DimensionMismatch(
        "bridge topology requires one canonical branch owner per passive position",
    ))
    valve_tuple = Tuple(valves)
    passive_tuple = Tuple(passives)
    leg_tuple = Tuple(bridge_legs)
    for (position, valve) in zip(topology.valve_positions, valve_tuple)
        _validate_bridge_topology_valve(position, valve)
    end
    for (position, passive) in zip(topology.passive_positions, passive_tuple)
        _validate_bridge_topology_passive(position, passive)
    end
    length(unique(objectid.(valve_tuple))) == length(valve_tuple) || throw(ArgumentError(
        "bridge topology valve positions must own distinct switch objects",
    ))
    for leg in leg_tuple
        any(valve -> valve === leg.upper_switch, valve_tuple) &&
            any(valve -> valve === leg.lower_switch, valve_tuple) || throw(ArgumentError(
            "bridge leg must alias two exact switch objects owned by the topology aggregate",
        ))
    end
    terminal_count = count(node -> node.role === :external, topology.nodes)
    runtime = PowerSemiconductorBridgeTopology(
        topology,
        valve_tuple,
        passive_tuple,
        leg_tuple,
        initially_blocked,
        fill(BRIDGE_POSITION_HEALTHY, length(valve_tuple)),
        0,
        0,
        zeros(Float64, terminal_count),
        zeros(Float64, terminal_count),
        0.0,
        0.0,
    )
    initially_blocked && block_power_semiconductor_topology!(runtime, 0.0)
    return runtime
end

power_semiconductor_bridge_topology_valves(bridge::PowerSemiconductorBridgeTopology) = bridge.valves
power_semiconductor_bridge_topology_linear_valves(bridge::PowerSemiconductorBridgeTopology) =
    tuple((valve for valve in bridge.valves if !power_semiconductor_has_extended_fidelity(valve))...)
power_semiconductor_bridge_topology_nonlinear_valves(bridge::PowerSemiconductorBridgeTopology) =
    tuple((valve for valve in bridge.valves if power_semiconductor_has_extended_fidelity(valve))...)

function _bridge_topology_requested_state(bridge::PowerSemiconductorBridgeTopology)
    return BitVector(map(bridge.valves) do valve
        valve.gate_driver === nothing ? false : valve.gate_driver.commanded_on
    end)
end

function _bridge_topology_applied_state(bridge::PowerSemiconductorBridgeTopology)
    return BitVector(map(bridge.valves) do valve
        valve.gate_driver === nothing ? false : valve.gate_driver.applied_on
    end)
end

function _bridge_topology_conducting_state(bridge::PowerSemiconductorBridgeTopology)
    return BitVector(map(bridge.valves) do valve
        valve.closed || valve.reverse_diode_conducting
    end)
end

function _validate_bridge_topology_requested_state!(bridge::PowerSemiconductorBridgeTopology)
    for group in bridge.topology.state_groups
        group.kind === :natural_conduction && continue
        state_is_allowed = false
        for column in axes(group.admitted_states, 2)
            state_matches = true
            for (row, position_index) in enumerate(group.position_indices)
                driver = bridge.valves[position_index].gate_driver
                requested = driver === nothing ? false : driver.commanded_on
                if requested != group.admitted_states[row, column]
                    state_matches = false
                    break
                end
            end
            if state_matches
                state_is_allowed = true
                break
            end
        end
        state_is_allowed || begin
            bridge.refusal_count += 1
            throw(ArgumentError("bridge gate request violates an admitted topology state"))
        end
    end
    return bridge
end

"""Request an exact ordered gate vector, validating all topology interlocks before mutation."""
function request_power_semiconductor_topology_gates!(
    bridge::PowerSemiconductorBridgeTopology,
    requested_state::AbstractVector{Bool},
    time_s::Real,
)
    length(requested_state) == length(bridge.valves) || throw(DimensionMismatch(
        "bridge gate request length must match valve-position count",
    ))
    bridge.blocked && any(requested_state) && begin
        bridge.refusal_count += 1
        throw(ArgumentError("a blocked bridge topology refuses turn-on requests"))
    end
    bridge_topology_state_is_allowed(bridge.topology, requested_state) || begin
        bridge.refusal_count += 1
        throw(ArgumentError("bridge gate request violates an admitted topology state"))
    end
    previous = _bridge_topology_requested_state(bridge)
    for index in eachindex(bridge.valves)
        valve = bridge.valves[index]
        valve.gate_driver === nothing && continue
        fault = bridge.position_faults[index]
        command = fault === BRIDGE_POSITION_STUCK_OPEN ? false :
            fault === BRIDGE_POSITION_STUCK_CLOSED ? true : requested_state[index]
        request_power_semiconductor_gate!(valve, command, time_s)
    end
    bridge.transition_count += count(previous .!= _bridge_topology_requested_state(bridge))
    return bridge
end

function block_power_semiconductor_topology!(bridge::PowerSemiconductorBridgeTopology, time_s::Real)
    changed = !bridge.blocked
    for (index, valve) in enumerate(bridge.valves)
        valve.gate_driver === nothing && continue
        bridge.position_faults[index] === BRIDGE_POSITION_STUCK_CLOSED && begin
            bridge.refusal_count += 1
            throw(ArgumentError("bridge block cannot override a declared stuck-closed position"))
        end
    end
    bridge.blocked = true
    for valve in bridge.valves
        valve.gate_driver === nothing && continue
        request_power_semiconductor_gate!(valve, false, time_s)
    end
    changed && (bridge.transition_count += 1)
    return changed
end

function restart_power_semiconductor_topology!(bridge::PowerSemiconductorBridgeTopology, time_s::Real)
    bridge.blocked || return false
    all(bridge.valves) do valve
        valve.gate_driver === nothing || (!valve.gate_driver.applied_on && valve.gate_driver.pending_state !== true)
    end || throw(ArgumentError("bridge restart requires every controlled position safely off"))
    bridge.blocked = false
    bridge.transition_count += 1
    return true
end

function apply_power_semiconductor_topology_fault!(
    bridge::PowerSemiconductorBridgeTopology,
    position_name::Symbol,
    fault::PowerSemiconductorTopologyFault,
    time_s::Real,
)
    index = _bridge_topology_position_index(bridge.topology, position_name)
    previous = bridge.position_faults[index]
    previous === fault && return false
    valve = bridge.valves[index]
    if fault === BRIDGE_POSITION_STUCK_OPEN
        valve.gate_driver === nothing || request_power_semiconductor_gate!(valve, false, time_s)
        valve.closed = false
        valve.reverse_diode_conducting = false
    elseif fault === BRIDGE_POSITION_STUCK_CLOSED
        valve.gate_driver === nothing && throw(ArgumentError(
            "stuck-closed forcing requires a controlled valve with an explicit gate owner",
        ))
        requested = _bridge_topology_requested_state(bridge)
        requested[index] = true
        bridge_topology_state_is_allowed(bridge.topology, requested) || begin
            bridge.refusal_count += 1
            throw(ArgumentError("stuck-closed fault would violate an admitted topology state"))
        end
        request_power_semiconductor_gate!(valve, true, time_s)
    end
    bridge.position_faults[index] = fault
    bridge.transition_count += 1
    return true
end

clear_power_semiconductor_topology_fault!(bridge, position_name, time_s) =
    apply_power_semiconductor_topology_fault!(
        bridge, position_name, BRIDGE_POSITION_HEALTHY, time_s,
    )

function power_semiconductor_event_localization!(bridge::PowerSemiconductorBridgeTopology)
    foreach(power_semiconductor_event_localization!, bridge.valves)
    return bridge
end

function _bridge_topology_apply_due_gate_transitions!(bridge::PowerSemiconductorBridgeTopology, time_s)
    count = 0
    for valve in bridge.valves
        transition_time = power_semiconductor_gate_transition_time(valve)
        transition_time === nothing && continue
        transition_time <= time_s + 16 * eps(max(1.0, abs(time_s))) &&
            (count += apply_power_semiconductor_gate_transition!(valve, time_s))
    end
    return count
end

function power_semiconductor_bridge_gate_transition_time(bridge::PowerSemiconductorBridgeTopology)
    next_time = Inf
    for valve in bridge.valves
        transition_time = power_semiconductor_gate_transition_time(valve)
        transition_time === nothing || (next_time = min(next_time, transition_time))
    end
    return isfinite(next_time) ? next_time : nothing
end

apply_power_semiconductor_bridge_gate_transitions!(
    bridge::PowerSemiconductorBridgeTopology,
    time_s::Real,
) = _bridge_topology_apply_due_gate_transitions!(bridge, Float64(time_s))

function _stamp_power_semiconductor_bridge_topology!(admittance, rhs,
    bridge::PowerSemiconductorBridgeTopology, time_s::Float64, dt_s::Float64, method)
    _bridge_topology_apply_due_gate_transitions!(bridge, time_s)
    _validate_bridge_topology_requested_state!(bridge)
    for valve in bridge.valves
        power_semiconductor_has_extended_fidelity(valve) && continue
        stamp!(admittance, rhs, valve, time_s, dt_s, method)
    end
    for passive in bridge.passives
        stamp!(admittance, rhs, passive, time_s, dt_s, method)
    end
    bridge.last_step_s = dt_s
    return nothing
end

function stamp!(admittance, rhs, bridge::PowerSemiconductorBridgeTopology,
    time_s::Float64, dt_s::Float64)
    return _stamp_power_semiconductor_bridge_topology!(
        admittance, rhs, bridge, time_s, dt_s, Val(TrapezoidalCompanion),
    )
end


function stamp!(admittance, rhs, bridge::PowerSemiconductorBridgeTopology,
    time_s::Float64, dt_s::Float64, method::Val{TrapezoidalCompanion})
    return _stamp_power_semiconductor_bridge_topology!(
        admittance, rhs, bridge, time_s, dt_s, method,
    )
end


function stamp!(admittance, rhs, bridge::PowerSemiconductorBridgeTopology,
    time_s::Float64, dt_s::Float64, method::Val{BackwardEulerCompanion})
    return _stamp_power_semiconductor_bridge_topology!(
        admittance, rhs, bridge, time_s, dt_s, method,
    )
end

Branches.backward_euler_companion_supported(bridge::PowerSemiconductorBridgeTopology) =
    all(Branches.backward_euler_companion_supported, bridge.passives) &&
    all(valve -> power_semiconductor_has_extended_fidelity(valve) ||
        Branches.backward_euler_companion_supported(valve), bridge.valves)

function _update_power_semiconductor_bridge_terminal_state!(
    bridge::PowerSemiconductorBridgeTopology,
    voltages::AbstractVector{Float64},
    dt_s::Float64,
)
    external_index = 0
    for node in bridge.topology.nodes
        node.role === :external || continue
        external_index += 1
        bridge.last_terminal_voltage_v[external_index] =
            node.node == 0 ? 0.0 : voltages[node.node]
        bridge.last_terminal_current_a[external_index] = 0.0
    end
    for valve in bridge.valves
        terminal_state = power_semiconductor_terminal_state(valve)
        external_index = 0
        for node in bridge.topology.nodes
            node.role === :external || continue
            external_index += 1
            node.node == valve.a &&
                (bridge.last_terminal_current_a[external_index] +=
                    terminal_state.terminal_current_a)
            node.node == valve.b &&
                (bridge.last_terminal_current_a[external_index] -=
                    terminal_state.terminal_current_a)
        end
    end
    for passive in bridge.passives
        passive_state = _bridge_topology_passive_state(passive, voltages, dt_s)
        external_index = 0
        for node in bridge.topology.nodes
            node.role === :external || continue
            external_index += 1
            node.node == passive.a &&
                (bridge.last_terminal_current_a[external_index] += passive_state.current)
            node.node == passive.b &&
                (bridge.last_terminal_current_a[external_index] -= passive_state.current)
        end
    end
    return bridge
end

function update!(bridge::PowerSemiconductorBridgeTopology, voltages, dt_s::Float64)
    for valve in bridge.valves
        power_semiconductor_has_extended_fidelity(valve) && continue
        update!(valve, voltages, dt_s)
    end
    for passive in bridge.passives
        update!(passive, voltages, dt_s)
    end
    passive_loss_w = sum(bridge.passives; init=0.0) do passive
        _bridge_topology_passive_state(passive, voltages, dt_s).loss
    end
    bridge.dissipated_energy_j += dt_s * passive_loss_w
    _update_power_semiconductor_bridge_terminal_state!(bridge, voltages, dt_s)
    bridge.last_step_s = dt_s
    return nothing
end

function _bridge_topology_passive_state(passive, voltages, dt_s)
    voltage = Branches.branch_voltage(voltages, passive.a, passive.b)
    current = branch_current_value(passive, voltages, dt_s)
    charge = passive isa CapacitorBranch ? passive.c * voltage :
        passive isa SeriesRLCBranch ? passive.c * passive.capacitor_voltage_prev : 0.0
    flux = passive isa SeriesRLBranch ? passive.l * current :
        passive isa SeriesRLCBranch ? passive.l * current : 0.0
    resistance = passive isa ConductanceBranch ? inv(passive.g) :
        hasproperty(passive, :r) ? passive.r : 0.0
    loss = resistance * current^2
    stored = 0.5 * (charge * (passive isa SeriesRLCBranch ? passive.capacitor_voltage_prev : voltage) +
        flux * current)
    return (; voltage, current, charge, flux, loss, stored)
end

function _bridge_topology_state_signature(topology_signature, requested, applied, conducting,
    faults, terminal_voltage, terminal_current, passive_voltage, passive_current, transitions, refusals)
    io = IOBuffer()
    print(io, topology_signature, '|', join(Int.(requested)), '|', join(Int.(applied)), '|',
        join(Int.(conducting)), '|', join(Int.(faults)), '|')
    for values in (terminal_voltage, terminal_current, passive_voltage, passive_current)
        for value in values
            print(io, bitstring(value), ';')
        end
        print(io, '|')
    end
    print(io, transitions, '|', refusals)
    return bytes2hex(SHA.sha256(take!(io)))
end

function power_semiconductor_bridge_topology_state(
    bridge::PowerSemiconductorBridgeTopology,
    voltages::AbstractVector{Float64},
    dt_s::Real=bridge.last_step_s,
)
    step = Float64(dt_s)
    isfinite(step) && step > 0.0 || throw(ArgumentError(
        "bridge topology state requires a finite positive accepted timestep",
    ))
    topology = bridge.topology
    topology_signature = bridge_topology_signature(topology)
    external_nodes = [node for node in topology.nodes if node.role === :external]
    terminal_voltage = Float64[node.node == 0 ? 0.0 : voltages[node.node] for node in external_nodes]
    terminal_current = zeros(Float64, length(external_nodes))
    terminal_lookup = Dict(node.node => index for (index, node) in enumerate(external_nodes))
    semiconductor_loss = 0.0
    semiconductor_dissipated_energy = 0.0
    for valve in bridge.valves
        state = power_semiconductor_terminal_state(valve)
        haskey(terminal_lookup, valve.a) && (terminal_current[terminal_lookup[valve.a]] += state.terminal_current_a)
        haskey(terminal_lookup, valve.b) && (terminal_current[terminal_lookup[valve.b]] -= state.terminal_current_a)
        semiconductor_loss += state.semiconductor_loss_w + state.snubber_resistor_loss_w
        semiconductor_dissipated_energy += state.semiconductor_dissipated_energy_j +
            state.snubber_dissipated_energy_j
    end
    passive_state = [_bridge_topology_passive_state(passive, voltages, step)
        for passive in bridge.passives]
    for (passive, state) in zip(bridge.passives, passive_state)
        haskey(terminal_lookup, passive.a) && (terminal_current[terminal_lookup[passive.a]] += state.current)
        haskey(terminal_lookup, passive.b) && (terminal_current[terminal_lookup[passive.b]] -= state.current)
    end
    requested = _bridge_topology_requested_state(bridge)
    applied = _bridge_topology_applied_state(bridge)
    conducting = _bridge_topology_conducting_state(bridge)
    passive_voltage = Float64[state.voltage for state in passive_state]
    passive_current = Float64[state.current for state in passive_state]
    passive_charge = Float64[state.charge for state in passive_state]
    passive_flux = Float64[state.flux for state in passive_state]
    passive_loss = sum(getproperty.(passive_state, :loss); init=0.0)
    stored_energy = sum(getproperty.(passive_state, :stored); init=0.0) +
        sum(valve -> begin
            state = power_semiconductor_terminal_state(valve)
            state.snubber_capacitor_energy_j +
                (power_semiconductor_has_extended_fidelity(valve) ?
                    power_semiconductor_extended_state(valve).junction_stored_energy_j : 0.0)
        end, bridge.valves; init=0.0)
    terminal_power = sum(terminal_voltage .* terminal_current)
    kcl_residual = sum(terminal_current)
    signature = _bridge_topology_state_signature(topology_signature, requested, applied,
        conducting, bridge.position_faults, terminal_voltage, terminal_current,
        passive_voltage, passive_current, bridge.transition_count, bridge.refusal_count)
    return PowerSemiconductorBridgeTopologyState(
        :aimora_bridge_result_v1,
        topology.family,
        topology_signature,
        getfield.(external_nodes, :name),
        terminal_voltage,
        terminal_current,
        getfield.(topology.valve_positions, :name),
        requested,
        applied,
        conducting,
        getfield.(topology.passive_positions, :name),
        passive_voltage,
        passive_current,
        passive_charge,
        passive_flux,
        kcl_residual,
        terminal_power,
        semiconductor_loss,
        passive_loss,
        stored_energy,
        semiconductor_dissipated_energy + bridge.dissipated_energy_j,
        bridge.blocked,
        copy(bridge.position_faults),
        bridge.transition_count + sum(power_semiconductor_bridge_topology_transition_count,
            bridge.bridge_legs; init=0),
        bridge.refusal_count,
        signature,
    )
end

trace_output_channel_count(bridge::PowerSemiconductorBridgeTopology) =
    6 + 3 * length(bridge.valves) + 2 * length(bridge.passives) +
    2 * count(node -> node.role === :external, bridge.topology.nodes)
trace_output_is_public(::PowerSemiconductorBridgeTopology) = true

function trace_output_channel_names!(names, element_name, bridge::PowerSemiconductorBridgeTopology)
    for node in bridge.topology.nodes
        node.role === :external || continue
        push!(names, Symbol(element_name, :_, node.name, :_voltage_v))
        push!(names, Symbol(element_name, :_, node.name, :_current_a))
    end
    for position in bridge.topology.valve_positions
        push!(names, Symbol(element_name, :_, position.name, :_requested))
        push!(names, Symbol(element_name, :_, position.name, :_applied))
        push!(names, Symbol(element_name, :_, position.name, :_conducting))
    end
    for position in bridge.topology.passive_positions
        push!(names, Symbol(element_name, :_, position.name, :_voltage_v))
        push!(names, Symbol(element_name, :_, position.name, :_current_a))
    end
    append!(names, Symbol[
        Symbol(element_name, :_terminal_kcl_residual_a),
        Symbol(element_name, :_terminal_power_w),
        Symbol(element_name, :_semiconductor_loss_w),
        Symbol(element_name, :_passive_loss_w),
        Symbol(element_name, :_stored_energy_j),
        Symbol(element_name, :_blocked),
    ])
    return names
end

function trace_output_values!(output, first_channel, sample, bridge::PowerSemiconductorBridgeTopology, voltage)
    state = power_semiconductor_bridge_topology_state(bridge, voltage, bridge.last_step_s)
    channel = first_channel
    for index in eachindex(state.terminal_names)
        output[channel, sample] = state.terminal_voltage_v[index]; channel += 1
        output[channel, sample] = state.terminal_current_a[index]; channel += 1
    end
    for index in eachindex(state.valve_names)
        output[channel, sample] = state.requested_gate_state[index]; channel += 1
        output[channel, sample] = state.applied_gate_state[index]; channel += 1
        output[channel, sample] = state.conducting_state[index]; channel += 1
    end
    for index in eachindex(state.passive_names)
        output[channel, sample] = state.passive_voltage_v[index]; channel += 1
        output[channel, sample] = state.passive_current_a[index]; channel += 1
    end
    for value in (
        state.terminal_kcl_residual_a,
        state.terminal_power_w,
        state.semiconductor_loss_w,
        state.passive_dissipated_power_w,
        state.stored_energy_j,
        state.blocked ? 1.0 : 0.0,
    )
        output[channel, sample] = value
        channel += 1
    end
    return channel
end

@enum PowerSemiconductorBridgeGateDisposition begin
    BRIDGE_GATE_ACCEPTED
    BRIDGE_GATE_REDUNDANT
    BRIDGE_GATE_BLOCKED
    BRIDGE_GATE_SHOOT_THROUGH_REJECTED
end

"""Typed three-terminal electrical, gating, topology, loss, and energy state of one self-commutated bridge leg."""
struct PowerSemiconductorBridgeTerminalState
    dc_positive_voltage_v::Float64
    ac_terminal_voltage_v::Float64
    dc_negative_voltage_v::Float64
    dc_positive_current_a::Float64
    ac_terminal_current_a::Float64
    dc_negative_current_a::Float64
    upper_switch::PowerSemiconductorTerminalState
    lower_switch::PowerSemiconductorTerminalState
    terminal_kcl_residual_a::Float64
    terminal_power_w::Float64
    semiconductor_loss_w::Float64
    snubber_resistor_loss_w::Float64
    dissipated_energy_j::Float64
    stored_energy_j::Float64
    blocked::Bool
    shoot_through_rejection_count::Int
    topology_transition_count::Int
end

"""Reusable two-position self-commutated bridge leg with natural antiparallel paths, complementary interlock, commutation dead time, ideal zero-recovery diode commutation, block/restart state, and typed terminal accounting."""
mutable struct PowerSemiconductorBridgeLeg{
    U<:PowerSemiconductorSwitch,
    L<:PowerSemiconductorSwitch,
} <: EMTElement
    upper_switch::U
    lower_switch::L
    commutation_dead_time_s::Float64
    blocked::Bool
    requested_upper_on::Bool
    requested_lower_on::Bool
    last_command_time_s::Float64
    command_count::Int
    shoot_through_rejection_count::Int
    block_count::Int
    restart_count::Int
    last_dc_positive_voltage_v::Float64
    last_ac_terminal_voltage_v::Float64
    last_dc_negative_voltage_v::Float64
end

function _self_commutated_bridge_switch(device::PowerSemiconductorSwitch)
    kind = power_semiconductor_kind(device)
    kind in (:igbt, :mosfet) || throw(ArgumentError(
        "a complementary bridge leg requires self-commutated IGBT or MOSFET switches",
    ))
    device.gate_driver !== nothing || throw(ArgumentError(
        "a complementary bridge switch requires a gate driver",
    ))
    return device
end

function PowerSemiconductorBridgeLeg(
    upper_switch::U,
    lower_switch::L;
    commutation_dead_time_s::Real=0.0,
    initially_blocked::Bool=false,
) where {
    U<:PowerSemiconductorSwitch,
    L<:PowerSemiconductorSwitch,
}
    _self_commutated_bridge_switch(upper_switch)
    _self_commutated_bridge_switch(lower_switch)
    upper_switch.b == lower_switch.a || throw(ArgumentError(
        "bridge upper cathode/source-side node must equal the lower anode/drain-side AC node",
    ))
    upper_switch.a != lower_switch.b || throw(ArgumentError(
        "bridge positive and negative DC terminals must be distinct",
    ))
    upper_switch.a != upper_switch.b && lower_switch.a != lower_switch.b ||
        throw(ArgumentError("bridge switches require distinct terminals"))
    dead_time = Float64(commutation_dead_time_s)
    isfinite(dead_time) && dead_time >= 0.0 || throw(ArgumentError(
        "bridge commutation dead time must be finite and nonnegative",
    ))
    upper_driver = something(upper_switch.gate_driver)
    lower_driver = something(lower_switch.gate_driver)
    !(upper_driver.applied_on && lower_driver.applied_on) || throw(ArgumentError(
        "bridge switches cannot both be initially applied",
    ))
    !(upper_driver.commanded_on && lower_driver.commanded_on) || throw(ArgumentError(
        "bridge switches cannot both be initially commanded",
    ))
    if initially_blocked
        !upper_driver.applied_on && !lower_driver.applied_on &&
            !upper_driver.commanded_on && !lower_driver.commanded_on &&
            upper_driver.pending_state !== true &&
            lower_driver.pending_state !== true || throw(ArgumentError(
            "an initially blocked bridge requires both switch gates to be commanded and applied off",
        ))
    end
    return PowerSemiconductorBridgeLeg{U,L}(
        upper_switch,
        lower_switch,
        dead_time,
        initially_blocked,
        upper_driver.commanded_on,
        lower_driver.commanded_on,
        0.0,
        0,
        0,
        0,
        0,
        0.0,
        0.0,
        0.0,
    )
end

function power_semiconductor_bridge_switch(
    bridge::PowerSemiconductorBridgeLeg,
    position::Symbol,
)
    position === :upper && return bridge.upper_switch
    position === :lower && return bridge.lower_switch
    throw(ArgumentError("bridge switch position must be :upper or :lower"))
end

power_semiconductor_bridge_topology_transition_count(
    bridge::PowerSemiconductorBridgeLeg,
) = bridge.upper_switch.topology_transition_count +
    bridge.lower_switch.topology_transition_count

function power_semiconductor_event_localization!(bridge::PowerSemiconductorBridgeLeg)
    power_semiconductor_event_localization!(bridge.upper_switch)
    power_semiconductor_event_localization!(bridge.lower_switch)
    return bridge
end

function _bridge_event_localization_enabled(
    bridge::PowerSemiconductorBridgeLeg,
)
    upper_enabled = bridge.upper_switch.event_localization_enabled
    lower_enabled = bridge.lower_switch.event_localization_enabled
    upper_enabled == lower_enabled || throw(ArgumentError(
        "bridge event localization must be enabled through the bridge owner so both complementary switches share one mutation boundary",
    ))
    return upper_enabled
end

function power_semiconductor_bridge_gate_transition_time(
    bridge::PowerSemiconductorBridgeLeg,
)
    upper_time = power_semiconductor_gate_transition_time(bridge.upper_switch)
    lower_time = power_semiconductor_gate_transition_time(bridge.lower_switch)
    upper_time === nothing && return lower_time
    lower_time === nothing && return upper_time
    return min(upper_time, lower_time)
end

function _bridge_gate_transition_due(
    device::PowerSemiconductorSwitch,
    transition_time_s::Float64,
)
    pending_time = power_semiconductor_gate_transition_time(device)
    pending_time === nothing && return false
    tolerance = 16 * eps(max(1.0, abs(transition_time_s)))
    return pending_time <= transition_time_s + tolerance
end

function _bridge_gate_state_after_due_transition(
    device::PowerSemiconductorSwitch,
    transition_due::Bool,
)
    driver = something(device.gate_driver)
    return transition_due ? something(driver.pending_state) : driver.applied_on
end

function _bridge_extinguish_opposing_freewheel_path!(
    outgoing::PowerSemiconductorSwitch,
    transition_time_s::Float64,
)
    outgoing.reverse_diode_conducting || return outgoing
    apply_power_semiconductor_reverse_extinction!(outgoing, transition_time_s)
    return outgoing
end

"""Apply every bridge gate edge due at `time_s`, always committing turn-off edges before turn-on edges and rejecting a simultaneous-on result before mutation."""
function apply_power_semiconductor_bridge_gate_transitions!(
    bridge::PowerSemiconductorBridgeLeg,
    time_s::Real,
)
    transition_time = _bridge_command_time(time_s)
    upper_due = _bridge_gate_transition_due(bridge.upper_switch, transition_time)
    lower_due = _bridge_gate_transition_due(bridge.lower_switch, transition_time)
    upper_due || lower_due || return 0
    upper_final = _bridge_gate_state_after_due_transition(
        bridge.upper_switch,
        upper_due,
    )
    lower_final = _bridge_gate_state_after_due_transition(
        bridge.lower_switch,
        lower_due,
    )
    upper_final && lower_final && throw(ArgumentError(
        "bridge gate transitions would apply both complementary positions",
    ))
    transition_count = 0
    for device in (bridge.upper_switch, bridge.lower_switch)
        driver = something(device.gate_driver)
        if _bridge_gate_transition_due(device, transition_time) &&
           driver.pending_state === false
            transition_count += apply_power_semiconductor_gate_transition!(
                device,
                transition_time,
            )
        end
    end
    for device in (bridge.upper_switch, bridge.lower_switch)
        driver = something(device.gate_driver)
        if _bridge_gate_transition_due(device, transition_time) &&
           driver.pending_state === true
            outgoing = device === bridge.upper_switch ?
                bridge.lower_switch : bridge.upper_switch
            _bridge_extinguish_opposing_freewheel_path!(
                outgoing,
                transition_time,
            )
            transition_count += apply_power_semiconductor_gate_transition!(
                device,
                transition_time,
            )
        end
    end
    _bridge_assert_interlock(bridge)
    return transition_count
end

function _bridge_command_time(time_s::Real)
    command_time = Float64(time_s)
    isfinite(command_time) && command_time >= 0.0 || throw(ArgumentError(
        "bridge command time must be finite and nonnegative",
    ))
    return command_time
end

function _bridge_gate_off_boundary!(
    device::PowerSemiconductorSwitch,
    command_time_s::Float64,
)
    driver = something(device.gate_driver)
    request_power_semiconductor_gate!(device, false, command_time_s)
    transition_time = power_semiconductor_gate_transition_time(device)
    if transition_time !== nothing && driver.pending_state === false
        return max(command_time_s, transition_time)
    end
    return max(command_time_s, driver.last_turn_off_time_s)
end

function _request_bridge_all_off!(
    bridge::PowerSemiconductorBridgeLeg,
    command_time_s::Float64,
)
    request_power_semiconductor_gate!(bridge.upper_switch, false, command_time_s)
    request_power_semiconductor_gate!(bridge.lower_switch, false, command_time_s)
    bridge.requested_upper_on = false
    bridge.requested_lower_on = false
    return bridge
end

function _request_bridge_position!(
    bridge::PowerSemiconductorBridgeLeg,
    incoming::PowerSemiconductorSwitch,
    outgoing::PowerSemiconductorSwitch,
    command_time_s::Float64,
)
    outgoing_off_time = _bridge_gate_off_boundary!(outgoing, command_time_s)
    incoming_transition_time =
        outgoing_off_time + bridge.commutation_dead_time_s
    tolerance = 16 * eps(max(1.0, abs(command_time_s)))
    incoming_transition_time <= command_time_s + tolerance &&
        _bridge_extinguish_opposing_freewheel_path!(
            outgoing,
            command_time_s,
        )
    request_power_semiconductor_gate!(
        incoming,
        true,
        command_time_s;
        earliest_transition_time_s = incoming_transition_time,
    )
    return bridge
end

"""Request the upper and lower bridge gates through the complementary interlock and return the typed disposition."""
function request_power_semiconductor_bridge_gates!(
    bridge::PowerSemiconductorBridgeLeg,
    upper_on::Bool,
    lower_on::Bool,
    time_s::Real,
)
    command_time = _bridge_command_time(time_s)
    bridge.command_count += 1
    bridge.last_command_time_s = command_time
    if bridge.blocked
        _request_bridge_all_off!(bridge, command_time)
        return BRIDGE_GATE_BLOCKED
    end
    if upper_on && lower_on
        bridge.shoot_through_rejection_count += 1
        _request_bridge_all_off!(bridge, command_time)
        return BRIDGE_GATE_SHOOT_THROUGH_REJECTED
    end
    if upper_on == bridge.requested_upper_on &&
       lower_on == bridge.requested_lower_on
        return BRIDGE_GATE_REDUNDANT
    end
    if upper_on
        _request_bridge_position!(
            bridge,
            bridge.upper_switch,
            bridge.lower_switch,
            command_time,
        )
    elseif lower_on
        _request_bridge_position!(
            bridge,
            bridge.lower_switch,
            bridge.upper_switch,
            command_time,
        )
    else
        _request_bridge_all_off!(bridge, command_time)
    end
    bridge.requested_upper_on = upper_on
    bridge.requested_lower_on = lower_on
    return BRIDGE_GATE_ACCEPTED
end

"""Request one bridge pole state, where `upper_on=true` selects the upper switch and `false` selects the lower switch."""
request_power_semiconductor_bridge_pole!(
    bridge::PowerSemiconductorBridgeLeg,
    upper_on::Bool,
    time_s::Real,
) = request_power_semiconductor_bridge_gates!(
    bridge,
    upper_on,
    !upper_on,
    time_s,
)

"""Block a bridge leg by cancelling every pending turn-on and commanding both switches off."""
function block_power_semiconductor_bridge!(
    bridge::PowerSemiconductorBridgeLeg,
    time_s::Real,
)
    command_time = _bridge_command_time(time_s)
    changed = !bridge.blocked
    bridge.blocked = true
    changed && (bridge.block_count += 1)
    bridge.last_command_time_s = command_time
    _request_bridge_all_off!(bridge, command_time)
    return changed
end

function _bridge_gate_is_safely_off(device::PowerSemiconductorSwitch)
    driver = something(device.gate_driver)
    return !driver.applied_on && driver.pending_state !== true
end

"""Release a blocked bridge only after both applied and pending gate states are safely off."""
function restart_power_semiconductor_bridge!(
    bridge::PowerSemiconductorBridgeLeg,
    time_s::Real,
)
    command_time = _bridge_command_time(time_s)
    bridge.blocked || return false
    _bridge_gate_is_safely_off(bridge.upper_switch) &&
        _bridge_gate_is_safely_off(bridge.lower_switch) || throw(ArgumentError(
            "bridge restart requires both applied and pending switch gates to be off",
        ))
    bridge.blocked = false
    bridge.restart_count += 1
    bridge.last_command_time_s = command_time
    bridge.requested_upper_on = false
    bridge.requested_lower_on = false
    return true
end

function _bridge_assert_interlock(bridge::PowerSemiconductorBridgeLeg)
    upper_driver = something(bridge.upper_switch.gate_driver)
    lower_driver = something(bridge.lower_switch.gate_driver)
    !(upper_driver.applied_on && lower_driver.applied_on) || throw(ArgumentError(
        "bridge complementary interlock detected simultaneous applied gates",
    ))
    return bridge
end

function stamp!(
    admittance::AbstractMatrix{Float64},
    rhs::AbstractVector{Float64},
    bridge::PowerSemiconductorBridgeLeg,
    time_s::Float64,
    dt_s::Float64,
)
    if !_bridge_event_localization_enabled(bridge)
        apply_power_semiconductor_bridge_gate_transitions!(bridge, time_s)
    end
    _bridge_assert_interlock(bridge)
    stamp!(admittance, rhs, bridge.upper_switch, time_s, dt_s)
    stamp!(admittance, rhs, bridge.lower_switch, time_s, dt_s)
    return nothing
end

function _bridge_node_voltage(
    voltages::AbstractVector{Float64},
    node::Int,
)
    return node == 0 ? 0.0 : voltages[node]
end

function update!(
    bridge::PowerSemiconductorBridgeLeg,
    voltages::AbstractVector{Float64},
    dt_s::Float64,
)
    update!(bridge.upper_switch, voltages, dt_s)
    update!(bridge.lower_switch, voltages, dt_s)
    bridge.last_dc_positive_voltage_v = _bridge_node_voltage(
        voltages,
        bridge.upper_switch.a,
    )
    bridge.last_ac_terminal_voltage_v = _bridge_node_voltage(
        voltages,
        bridge.upper_switch.b,
    )
    bridge.last_dc_negative_voltage_v = _bridge_node_voltage(
        voltages,
        bridge.lower_switch.b,
    )
    _bridge_assert_interlock(bridge)
    return nothing
end

"""Return the three-terminal KCL, switch state, loss, stored-energy, and topology snapshot of a bridge leg."""
function power_semiconductor_bridge_terminal_state(
    bridge::PowerSemiconductorBridgeLeg,
)
    upper = power_semiconductor_terminal_state(bridge.upper_switch)
    lower = power_semiconductor_terminal_state(bridge.lower_switch)
    dc_positive_current = upper.terminal_current_a
    ac_terminal_current = lower.terminal_current_a - upper.terminal_current_a
    dc_negative_current = -lower.terminal_current_a
    kcl_residual = dc_positive_current + ac_terminal_current + dc_negative_current
    terminal_power =
        bridge.last_dc_positive_voltage_v * dc_positive_current +
        bridge.last_ac_terminal_voltage_v * ac_terminal_current +
        bridge.last_dc_negative_voltage_v * dc_negative_current
    semiconductor_loss = upper.semiconductor_loss_w + lower.semiconductor_loss_w
    snubber_loss = upper.snubber_resistor_loss_w + lower.snubber_resistor_loss_w
    dissipated_energy =
        upper.semiconductor_dissipated_energy_j +
        lower.semiconductor_dissipated_energy_j +
        upper.snubber_dissipated_energy_j +
        lower.snubber_dissipated_energy_j
    stored_energy =
        upper.snubber_capacitor_energy_j + lower.snubber_capacitor_energy_j
    return PowerSemiconductorBridgeTerminalState(
        bridge.last_dc_positive_voltage_v,
        bridge.last_ac_terminal_voltage_v,
        bridge.last_dc_negative_voltage_v,
        dc_positive_current,
        ac_terminal_current,
        dc_negative_current,
        upper,
        lower,
        kcl_residual,
        terminal_power,
        semiconductor_loss,
        snubber_loss,
        dissipated_energy,
        stored_energy,
        bridge.blocked,
        bridge.shoot_through_rejection_count,
        power_semiconductor_bridge_topology_transition_count(bridge),
    )
end

trace_output_channel_count(::PowerSemiconductorBridgeLeg) = 24
trace_output_is_public(::PowerSemiconductorBridgeLeg) = true

function trace_output_channel_names!(
    names::Vector{Symbol},
    element_name::Symbol,
    ::PowerSemiconductorBridgeLeg,
)
    append!(
        names,
        Symbol[
            Symbol(element_name, :_dc_positive_current_a),
            Symbol(element_name, :_ac_terminal_current_a),
            Symbol(element_name, :_dc_negative_current_a),
            Symbol(element_name, :_upper_voltage_v),
            Symbol(element_name, :_upper_current_a),
            Symbol(element_name, :_lower_voltage_v),
            Symbol(element_name, :_lower_current_a),
            Symbol(element_name, :_upper_forward_conducting),
            Symbol(element_name, :_upper_reverse_diode_conducting),
            Symbol(element_name, :_upper_gate_commanded),
            Symbol(element_name, :_upper_gate_applied),
            Symbol(element_name, :_lower_forward_conducting),
            Symbol(element_name, :_lower_reverse_diode_conducting),
            Symbol(element_name, :_lower_gate_commanded),
            Symbol(element_name, :_lower_gate_applied),
            Symbol(element_name, :_semiconductor_loss_w),
            Symbol(element_name, :_snubber_resistor_loss_w),
            Symbol(element_name, :_dissipated_energy_j),
            Symbol(element_name, :_stored_energy_j),
            Symbol(element_name, :_blocked),
            Symbol(element_name, :_shoot_through_rejection_count),
            Symbol(element_name, :_topology_transition_count),
            Symbol(element_name, :_terminal_kcl_residual_a),
            Symbol(element_name, :_terminal_power_w),
        ],
    )
    return names
end

function trace_output_values!(
    output::AbstractMatrix{Float64},
    first_channel::Int,
    sample::Int,
    bridge::PowerSemiconductorBridgeLeg,
    _voltage::AbstractVector{Float64},
)
    terminal = power_semiconductor_bridge_terminal_state(bridge)
    upper = terminal.upper_switch
    lower = terminal.lower_switch
    values = (
        terminal.dc_positive_current_a,
        terminal.ac_terminal_current_a,
        terminal.dc_negative_current_a,
        upper.terminal_voltage_v,
        upper.terminal_current_a,
        lower.terminal_voltage_v,
        lower.terminal_current_a,
        upper.forward_conducting ? 1.0 : 0.0,
        upper.reverse_diode_conducting ? 1.0 : 0.0,
        upper.gate_commanded_on ? 1.0 : 0.0,
        upper.gate_applied_on ? 1.0 : 0.0,
        lower.forward_conducting ? 1.0 : 0.0,
        lower.reverse_diode_conducting ? 1.0 : 0.0,
        lower.gate_commanded_on ? 1.0 : 0.0,
        lower.gate_applied_on ? 1.0 : 0.0,
        terminal.semiconductor_loss_w,
        terminal.snubber_resistor_loss_w,
        terminal.dissipated_energy_j,
        terminal.stored_energy_j,
        terminal.blocked ? 1.0 : 0.0,
        Float64(terminal.shoot_through_rejection_count),
        Float64(terminal.topology_transition_count),
        terminal.terminal_kcl_residual_a,
        terminal.terminal_power_w,
    )
    for offset in eachindex(values)
        output[first_channel + offset - 1, sample] = values[offset]
    end
    return first_channel + length(values)
end

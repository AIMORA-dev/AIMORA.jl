abstract type PowerSemiconductorJunction end

import ..NonlinearNetwork: _validate_nonlinear_current_correction_buffers

struct DiodeJunction <: PowerSemiconductorJunction end
struct ThyristorJunction <: PowerSemiconductorJunction end
struct TriacJunction <: PowerSemiconductorJunction end
struct GateTurnOffThyristorJunction <: PowerSemiconductorJunction end
struct InsulatedGateBipolarTransistorJunction <: PowerSemiconductorJunction end
struct MetalOxideSemiconductorFieldEffectTransistorJunction <: PowerSemiconductorJunction end

"""Inertial gate-command state with separate turn-on/turn-off delay, minimum off-time dead time, and minimum applied pulse width."""
mutable struct PowerSemiconductorGateDriver
    turn_on_delay_s::Float64
    turn_off_delay_s::Float64
    dead_time_s::Float64
    minimum_pulse_width_s::Float64
    commanded_on::Bool
    applied_on::Bool
    pending_state::Union{Nothing,Bool}
    pending_transition_time_s::Float64
    last_command_time_s::Float64
    last_turn_on_time_s::Float64
    last_turn_off_time_s::Float64
    command_count::Int
    transition_count::Int
    filtered_pulse_count::Int
end

function PowerSemiconductorGateDriver(;
    turn_on_delay_s::Real=0.0,
    turn_off_delay_s::Real=0.0,
    dead_time_s::Real=0.0,
    minimum_pulse_width_s::Real=0.0,
    initially_on::Bool=false,
)
    turn_on_delay = Float64(turn_on_delay_s)
    turn_off_delay = Float64(turn_off_delay_s)
    dead_time = Float64(dead_time_s)
    minimum_pulse_width = Float64(minimum_pulse_width_s)
    all(isfinite, (turn_on_delay, turn_off_delay, dead_time, minimum_pulse_width)) ||
        throw(ArgumentError("gate-driver times must be finite"))
    all(>=(0.0), (turn_on_delay, turn_off_delay, dead_time, minimum_pulse_width)) ||
        throw(ArgumentError("gate-driver times must be nonnegative"))
    return PowerSemiconductorGateDriver(
        turn_on_delay,
        turn_off_delay,
        dead_time,
        minimum_pulse_width,
        initially_on,
        initially_on,
        nothing,
        Inf,
        0.0,
        initially_on ? 0.0 : -Inf,
        initially_on ? -Inf : 0.0,
        0,
        0,
        0,
    )
end

"""A gate-transition time and the deterministic physical constraint that owns it."""
struct PowerSemiconductorGateTimingSelection
    transition_time_s::Float64
    active_constraint::Symbol
end

function _later_power_semiconductor_gate_timing(
    retained::PowerSemiconductorGateTimingSelection,
    candidate_time_s::Float64,
    candidate_constraint::Symbol,
)
    isnan(candidate_time_s) && throw(ArgumentError(
        "gate-timing constraint must not be NaN",
    ))
    candidate_time_s > retained.transition_time_s || return retained
    return PowerSemiconductorGateTimingSelection(
        candidate_time_s,
        candidate_constraint,
    )
end

"""
Select the active gate-transition timing constraint in a fixed physical order.

An exact tie retains the first constraint: commanded delay, prior gate-state
constraint, then an external interlock. This direction-independent selection
defines one valid semismooth chart for batched tangents without changing the
primal transition time.
"""
function power_semiconductor_gate_transition_timing(
    driver::PowerSemiconductorGateDriver,
    commanded_on::Bool,
    time_s::Real;
    earliest_transition_time_s::Union{Nothing,Real}=nothing,
)
    command_time = Float64(time_s)
    isfinite(command_time) && command_time >= 0.0 || throw(ArgumentError(
        "gate command time must be finite and nonnegative",
    ))
    earliest_transition_time = if earliest_transition_time_s === nothing
        nothing
    else
        boundary = Float64(earliest_transition_time_s)
        isfinite(boundary) && boundary >= command_time || throw(ArgumentError(
            "earliest gate-transition time must be finite and no earlier than the command",
        ))
        boundary
    end
    selection = if commanded_on
        PowerSemiconductorGateTimingSelection(
            command_time + driver.turn_on_delay_s,
            :turn_on_delay,
        )
    else
        PowerSemiconductorGateTimingSelection(
            command_time + driver.turn_off_delay_s,
            :turn_off_delay,
        )
    end
    if commanded_on
        selection = _later_power_semiconductor_gate_timing(
            selection,
            driver.last_turn_off_time_s + driver.dead_time_s,
            :minimum_off_time,
        )
    else
        selection = _later_power_semiconductor_gate_timing(
            selection,
            driver.last_turn_on_time_s + driver.minimum_pulse_width_s,
            :minimum_pulse_width,
        )
    end
    if earliest_transition_time !== nothing
        selection = _later_power_semiconductor_gate_timing(
            selection,
            earliest_transition_time,
            :external_interlock,
        )
    end
    isfinite(selection.transition_time_s) || throw(ArgumentError(
        "active gate-transition timing constraint must be finite",
    ))
    return selection
end

"""Select the deterministic command or completed/pending turn-off boundary used by a complementary bridge interlock."""
function power_semiconductor_gate_off_boundary_timing(
    driver::PowerSemiconductorGateDriver,
    time_s::Real,
)
    command_time = Float64(time_s)
    isfinite(command_time) && command_time >= 0.0 || throw(ArgumentError(
        "gate command time must be finite and nonnegative",
    ))
    selection = PowerSemiconductorGateTimingSelection(
        command_time,
        :gate_command_time,
    )
    if driver.pending_state === false
        return _later_power_semiconductor_gate_timing(
            selection,
            driver.pending_transition_time_s,
            :pending_turn_off_transition,
        )
    end
    return _later_power_semiconductor_gate_timing(
        selection,
        driver.last_turn_off_time_s,
        :last_turn_off_transition,
    )
end

"""Generic piecewise-linear antiparallel-diode parameters; current is positive from the semiconductor cathode-side terminal back toward its anode-side terminal."""
struct AntiparallelDiodeParameters
    forward_voltage_v::Float64
    holding_current_a::Float64
    on_conductance_s::Float64
end

function AntiparallelDiodeParameters(;
    forward_voltage_v::Real=0.0,
    holding_current_a::Real=0.0,
    on_conductance_s::Real=1.0e3,
)
    forward_voltage = Float64(forward_voltage_v)
    holding_current = Float64(holding_current_a)
    on_conductance = Float64(on_conductance_s)
    isfinite(forward_voltage) && forward_voltage >= 0.0 || throw(ArgumentError(
        "antiparallel-diode forward voltage must be finite and nonnegative",
    ))
    isfinite(holding_current) && holding_current >= 0.0 || throw(ArgumentError(
        "antiparallel-diode holding current must be finite and nonnegative",
    ))
    isfinite(on_conductance) && on_conductance > 0.0 || throw(ArgumentError(
        "antiparallel-diode on conductance must be finite and positive",
    ))
    return AntiparallelDiodeParameters(forward_voltage, holding_current, on_conductance)
end

"""Causal RC accepted companion at original fixed kernel/storage precision."""
struct AcceptedRCCompanionHistory
    incoming_capacitor_voltage_v::Float64
    incoming_current_a::Float64
    physical_voltage_v::BigFloat
    step_s::Float64
    method::Symbol
    current_a::Float64
    outgoing_capacitor_voltage_v::Float64
    arithmetic_policy::Symbol
end

struct RCBoundaryCurrentAnchor
    physical_voltage_v::BigFloat
    outgoing_capacitor_voltage_v::Float64
    current_a::Float64
    arithmetic_policy::Symbol
end

"""Trapezoidal/backward-Euler companion state for a series RC semiconductor snubber."""
mutable struct SeriesRCSnubber
    resistance_ohm::Float64
    capacitance_f::Float64
    previous_current_a::Float64
    capacitor_voltage_v::Float64
    last_branch_voltage_v::Float64
    last_current_a::Float64
    last_resistor_loss_w::Float64
    dissipated_energy_j::Float64
    accepted_companion_history::Union{Nothing,AcceptedRCCompanionHistory}
    boundary_voltage_v::Union{Nothing,BigFloat}
    boundary_anchor::Union{Nothing,RCBoundaryCurrentAnchor}
end

function SeriesRCSnubber(
    resistance_ohm::Real,
    capacitance_f::Real;
    initial_current_a::Real=0.0,
    initial_capacitor_voltage_v::Real=0.0,
)
    resistance = Float64(resistance_ohm)
    capacitance = Float64(capacitance_f)
    current = Float64(initial_current_a)
    capacitor_voltage = Float64(initial_capacitor_voltage_v)
    isfinite(resistance) && resistance > 0.0 ||
        throw(ArgumentError("snubber resistance must be finite and positive"))
    isfinite(capacitance) && capacitance > 0.0 ||
        throw(ArgumentError("snubber capacitance must be finite and positive"))
    isfinite(current) && isfinite(capacitor_voltage) || throw(ArgumentError(
        "snubber initial current and capacitor voltage must be finite",
    ))
    initial_loss = resistance * current^2
    return SeriesRCSnubber(
        resistance,
        capacitance,
        current,
        capacitor_voltage,
        capacitor_voltage + resistance * current,
        current,
        initial_loss,
        0.0,
        nothing,
        nothing,
        nothing,
    )
end

"""Typed terminal, path-state, loss, and stored-energy snapshot for one accepted semiconductor state."""
struct PowerSemiconductorTerminalState
    device_kind::Symbol
    forward_conducting::Bool
    reverse_diode_conducting::Bool
    gate_commanded_on::Bool
    gate_applied_on::Bool
    terminal_voltage_v::Float64
    terminal_current_a::Float64
    forward_current_a::Float64
    reverse_diode_current_a::Float64
    snubber_current_a::Float64
    semiconductor_loss_w::Float64
    snubber_resistor_loss_w::Float64
    semiconductor_dissipated_energy_j::Float64
    snubber_dissipated_energy_j::Float64
    snubber_capacitor_energy_j::Float64
end

"""
Evidence-bounded semiconductor switch whose unchanged baseline owns piecewise-linear diode, thyristor, IGBT, and MOSFET behavior. Explicit optional fidelity components add TRIAC/GTO semantics, recovered charge, nonlinear junction charge, tail current, event-energy maps, and passive lumped thermal state without changing the baseline constructor or trace.

Validity requires an EMT timestep fine enough to resolve the external circuit and every admitted dynamic component. Generic extended parameters remain AIMORA-authored physical contracts and never establish manufacturer prediction, arbitrary compact-model compatibility, standard conformance, ATP/PSCAD equivalence, or certification.
"""
mutable struct PowerSemiconductorSwitch{K<:PowerSemiconductorJunction} <: EMTElement
    a::Int
    b::Int
    threshold_v::Float64
    holding_current::Float64
    on_conductance::Float64
    off_conductance::Float64
    closed::Bool
    last_voltage::Float64
    last_current::Float64
    last_conductance::Float64
    forward_voltage_drop_v::Float64
    gate_driver::Union{Nothing,PowerSemiconductorGateDriver}
    antiparallel_diode::Union{Nothing,AntiparallelDiodeParameters}
    snubber::Union{Nothing,SeriesRCSnubber}
    reverse_diode_conducting::Bool
    event_localization_enabled::Bool
    last_evaluation_time_s::Float64
    last_history_current_a::Float64
    last_forward_current_a::Float64
    last_reverse_diode_current_a::Float64
    last_snubber_current_a::Float64
    last_semiconductor_loss_w::Float64
    previous_semiconductor_loss_w::Float64
    semiconductor_dissipated_energy_j::Float64
    topology_transition_count::Int
    last_transition_time_s::Float64
    extended_fidelity::Union{Nothing,PowerSemiconductorExtendedFidelity}
    conduction_direction::Int8
    gate_turn_off_current_limit_a::Float64
    gate_turn_off_policy::Symbol
    gate_turn_off_disposition::Symbol
end

"""Piecewise-linear naturally commutated diode specialization."""
const DiodeValveSwitch = PowerSemiconductorSwitch{DiodeJunction}
"""Piecewise-linear gate-triggered, holding-current-commutated thyristor specialization."""
const ThyristorValveSwitch = PowerSemiconductorSwitch{ThyristorJunction}
"""Bidirectional latching TRIAC composed by the accepted semiconductor terminal owner."""
const TriacSwitch = PowerSemiconductorSwitch{TriacJunction}
"""Forward latching gate-turn-off thyristor with an explicit interruptible-current policy."""
const GateTurnOffThyristorSwitch = PowerSemiconductorSwitch{GateTurnOffThyristorJunction}
"""Piecewise-linear forward-channel IGBT specialization with an optional antiparallel diode."""
const IGBTSwitch = PowerSemiconductorSwitch{InsulatedGateBipolarTransistorJunction}
"""Piecewise-linear bidirectional-on-channel MOSFET specialization with an optional body-diode path."""
const MOSFETSwitch = PowerSemiconductorSwitch{MetalOxideSemiconductorFieldEffectTransistorJunction}

"""Minimal accepted chart state needed to audit one semiconductor topology transition."""
struct PowerSemiconductorTopologyChartState
    forward_conducting::Bool
    reverse_diode_conducting::Bool
    topology_transition_count::Int
    last_transition_time_s::Float64
    gate_applied_on::Union{Nothing,Bool}
    gate_transition_count::Int
    tail_active::Bool
    last_forward_current_a::Float64
    last_voltage_v::Float64
end

"""One ordered physical topology transition observed between numerical callbacks."""
struct PowerSemiconductorTopologyTransitionObservation{D<:PowerSemiconductorSwitch}
    device::D
    transition_owner::Symbol
    before::PowerSemiconductorTopologyChartState
    after::PowerSemiconductorTopologyChartState
end

function power_semiconductor_topology_chart_state(
    device::PowerSemiconductorSwitch,
)
    driver = device.gate_driver
    fidelity = device.extended_fidelity
    tail = fidelity === nothing ? nothing : fidelity.turn_off_tail
    return PowerSemiconductorTopologyChartState(
        device.closed,
        device.reverse_diode_conducting,
        device.topology_transition_count,
        device.last_transition_time_s,
        driver === nothing ? nothing : driver.applied_on,
        driver === nothing ? 0 : driver.transition_count,
        tail !== nothing && tail.active,
        device.last_forward_current_a,
        device.last_voltage,
    )
end

"""Append exactly one observed topology transition, rejecting a phase that hides multiple transitions."""
function record_power_semiconductor_topology_transition_observation!(
    journal::AbstractVector,
    device::PowerSemiconductorSwitch,
    before::PowerSemiconductorTopologyChartState,
    transition_owner::Symbol,
    ; tail_history=nothing, after=power_semiconductor_topology_chart_state(device),
)
    transition_owner === :none && throw(ArgumentError(
        "a semiconductor topology transition requires a physical owner",
    ))
    topology_increment = after.topology_transition_count -
        before.topology_transition_count
    topology_increment >= 0 || throw(ArgumentError(
        "semiconductor topology transition count moved backward",
    ))
    gate_increment = after.gate_transition_count - before.gate_transition_count
    gate_increment in (0, 1) || throw(ArgumentError(
        "one observed semiconductor phase crossed more than one gate transition",
    ))
    if topology_increment == 0
        before.forward_conducting == after.forward_conducting &&
            before.reverse_diode_conducting ==
                after.reverse_diode_conducting || throw(ArgumentError(
            "semiconductor conduction changed without a topology-transition count",
        ))
        return false
    end
    topology_increment == 1 || throw(ArgumentError(
        "one observed semiconductor phase crossed more than one topology transition",
    ))
    (before.forward_conducting != after.forward_conducting ||
     before.reverse_diode_conducting != after.reverse_diode_conducting) || throw(
        ArgumentError(
            "a semiconductor topology-transition count changed without a conduction-state change",
        ),
    )
    isfinite(after.last_transition_time_s) || throw(ArgumentError(
        "an observed semiconductor topology transition requires a finite time",
    ))
    push!(journal, (;
        device, transition_owner, before, after, tail_history,
        time_s=after.last_transition_time_s,
    ))
    return true
end

"""Retain the causal tail history before a topology phase for its derivative.

This observation does not change state. A new turn-off rebases and adds to a
surviving exponential, so the previous amplitude and clock cannot be recovered
from the replaced history alone.
"""
function power_semiconductor_topology_tail_history(device::PowerSemiconductorSwitch)
    fidelity = device.extended_fidelity
    tail = fidelity === nothing ? nothing : fidelity.turn_off_tail
    return tail === nothing ? nothing : (;
        active=tail.active, initial_current_a=tail.initial_current_a,
        turn_off_time_s=tail.turn_off_time_s,
    )
end

NonlinearNetwork.nonlinear_device_event_transition_snapshot(device::PowerSemiconductorSwitch) = (;
    chart=power_semiconductor_topology_chart_state(device),
    tail_history=power_semiconductor_topology_tail_history(device),
)

function NonlinearNetwork.record_nonlinear_device_event_transition_observation!(
    journal::AbstractVector, device::PowerSemiconductorSwitch, before, after,
    surface::NonlinearNetwork.NonlinearDeviceEventSurface,
    time_s::Union{Float64,ContinuousEvaluationTime}, endpoint_owner::Symbol,
)
    endpoint_owner in (:left, :right) || throw(ArgumentError(
        "semiconductor event observation has no substep endpoint owner",
    ))
    first_index = length(journal) + 1
    record_power_semiconductor_topology_transition_observation!(
        journal, device, before.chart, surface.name;
        tail_history=before.tail_history, after=after.chart,
    )
    for index in first_index:length(journal)
        observation = journal[index]
        journal[index] = (;
            device, transition_owner=observation.transition_owner,
            before=observation.before, after=observation.after,
            tail_history=before.tail_history, time_s, localized_event=true, endpoint_owner,
        )
    end
    # Reconstruction occurs after transition application and before the next
    # algebraic solve. Capture that exact stage, including non-topology events.
    push!(journal, (;
        device, boundary_reconstruction=true, time_s,
        voltage_v=device.last_voltage,
        channel_current_a=_power_semiconductor_accepted_boundary_channel_current_jacobian(
            device, device.last_voltage)[1],
        tail_current_a=device.extended_fidelity === nothing ||
            device.extended_fidelity.turn_off_tail === nothing ? 0.0 :
            device.extended_fidelity.turn_off_tail.current_a,
        prepare_pending_event=false, localized_event=true, endpoint_owner,
    ))
    return nothing
end

power_semiconductor_kind(::PowerSemiconductorSwitch{DiodeJunction}) = :diode
power_semiconductor_kind(::PowerSemiconductorSwitch{ThyristorJunction}) = :thyristor
power_semiconductor_kind(::PowerSemiconductorSwitch{TriacJunction}) = :triac
power_semiconductor_kind(::PowerSemiconductorSwitch{GateTurnOffThyristorJunction}) = :gto
power_semiconductor_kind(::PowerSemiconductorSwitch{InsulatedGateBipolarTransistorJunction}) = :igbt
power_semiconductor_kind(::PowerSemiconductorSwitch{MetalOxideSemiconductorFieldEffectTransistorJunction}) = :mosfet

function _power_semiconductor_switch(
    ::Type{K},
    a::Int,
    b::Int;
    threshold_v::Real=0.0,
    forward_voltage_drop_v::Real=0.0,
    holding_current::Real=0.0,
    on_conductance::Real=1.0e3,
    off_conductance::Real=0.0,
    initially_closed::Bool=false,
    gate_driver::Union{Nothing,PowerSemiconductorGateDriver}=nothing,
    antiparallel_diode::Union{Nothing,AntiparallelDiodeParameters}=nothing,
    snubber::Union{Nothing,SeriesRCSnubber}=nothing,
    extended_fidelity::Union{Nothing,PowerSemiconductorExtendedFidelity}=nothing,
    gate_turn_off_current_limit_a::Real=Inf,
    gate_turn_off_policy::Symbol=:not_applicable,
) where {K<:PowerSemiconductorJunction}
    a >= 0 && b >= 0 ||
        throw(ArgumentError("power-semiconductor nodes must be nonnegative"))
    a != b || throw(ArgumentError("power-semiconductor terminals must be distinct"))
    threshold = Float64(threshold_v)
    forward_drop = Float64(forward_voltage_drop_v)
    holding = Float64(holding_current)
    on = Float64(on_conductance)
    off = Float64(off_conductance)
    isfinite(threshold) && threshold >= 0.0 || throw(ArgumentError(
        "power-semiconductor turn-on threshold must be finite and nonnegative",
    ))
    isfinite(forward_drop) && forward_drop >= 0.0 || throw(ArgumentError(
        "power-semiconductor forward voltage drop must be finite and nonnegative",
    ))
    isfinite(holding) && holding >= 0.0 || throw(ArgumentError(
        "power-semiconductor holding current must be finite and nonnegative",
    ))
    isfinite(on) && on > 0.0 || throw(ArgumentError(
        "power-semiconductor on conductance must be finite and positive",
    ))
    isfinite(off) && off >= 0.0 && off <= on || throw(ArgumentError(
        "power-semiconductor off conductance must be finite, nonnegative, and no greater than on conductance",
    ))
    K === DiodeJunction && gate_driver !== nothing && throw(ArgumentError(
        "a diode junction cannot own a gate driver",
    ))
    K !== DiodeJunction && gate_driver === nothing && throw(ArgumentError(
        "a controlled power semiconductor requires a gate driver",
    ))
    if K === TriacJunction
        antiparallel_diode === nothing || throw(ArgumentError(
            "a TRIAC owns its two oriented latching channels and cannot add a baseline antiparallel diode",
        ))
    end
    turn_off_current_limit = Float64(gate_turn_off_current_limit_a)
    if K === GateTurnOffThyristorJunction
        isfinite(turn_off_current_limit) && turn_off_current_limit >= 0.0 ||
            throw(ArgumentError(
                "GTO interruptible-current limit must be finite and nonnegative",
            ))
        gate_turn_off_policy in (:refuse, :remain_latched) || throw(ArgumentError(
            "GTO turn-off policy must be :refuse or :remain_latched",
        ))
    else
        turn_off_current_limit == Inf || throw(ArgumentError(
            "only a GTO may declare an interruptible-current limit",
        ))
        gate_turn_off_policy === :not_applicable || throw(ArgumentError(
            "only a GTO may declare a gate-turn-off policy",
        ))
    end
    if extended_fidelity !== nothing
        recovery = extended_fidelity.recovered_charge
        recovery === nothing || K === DiodeJunction || antiparallel_diode !== nothing ||
            throw(ArgumentError(
                "recovered charge requires a diode junction or an explicit antiparallel diode path",
            ))
        tail = extended_fidelity.turn_off_tail
        tail === nothing || K === InsulatedGateBipolarTransistorJunction ||
            K === GateTurnOffThyristorJunction || throw(ArgumentError(
            "turn-off tail fidelity is admitted only for IGBT or GTO junctions",
        ))
    end
    if K === InsulatedGateBipolarTransistorJunction ||
       K === MetalOxideSemiconductorFieldEffectTransistorJunction
        gate_driver.applied_on == initially_closed || throw(ArgumentError(
            "initial gated-channel state must match the applied gate-driver state",
        ))
    end
    conductance = initially_closed ? on : off
    return PowerSemiconductorSwitch{K}(
        a,
        b,
        threshold,
        holding,
        on,
        off,
        initially_closed,
        0.0,
        0.0,
        conductance,
        forward_drop,
        gate_driver,
        antiparallel_diode,
        snubber,
        false,
        false,
        0.0,
        initially_closed ? -on * forward_drop : 0.0,
        0.0,
        0.0,
        snubber === nothing ? 0.0 : snubber.last_current_a,
        0.0,
        0.0,
        0.0,
        0,
        0.0,
        extended_fidelity,
        initially_closed ? Int8(1) : Int8(0),
        turn_off_current_limit,
        gate_turn_off_policy,
        :not_requested,
    )
end

function PowerSemiconductorSwitch{DiodeJunction}(
    a::Int,
    b::Int;
    threshold_v::Real=0.0,
    forward_voltage_drop_v::Real=0.0,
    holding_current::Real=0.0,
    on_conductance::Real=1.0e9,
    off_conductance::Real=0.0,
    initially_closed::Bool=false,
    snubber::Union{Nothing,SeriesRCSnubber}=nothing,
    extended_fidelity::Union{Nothing,PowerSemiconductorExtendedFidelity}=nothing,
)
    return _power_semiconductor_switch(
        DiodeJunction,
        a,
        b;
        threshold_v,
        forward_voltage_drop_v,
        holding_current,
        on_conductance,
        off_conductance,
        initially_closed,
        snubber,
        extended_fidelity,
    )
end

function PowerSemiconductorSwitch{ThyristorJunction}(
    a::Int,
    b::Int;
    gate_driver::PowerSemiconductorGateDriver=PowerSemiconductorGateDriver(),
    threshold_v::Real=0.0,
    forward_voltage_drop_v::Real=0.0,
    holding_current::Real=0.0,
    on_conductance::Real=1.0e3,
    off_conductance::Real=0.0,
    initially_closed::Bool=false,
    antiparallel_diode::Union{Nothing,AntiparallelDiodeParameters}=nothing,
    snubber::Union{Nothing,SeriesRCSnubber}=nothing,
    extended_fidelity::Union{Nothing,PowerSemiconductorExtendedFidelity}=nothing,
)
    return _power_semiconductor_switch(
        ThyristorJunction,
        a,
        b;
        threshold_v,
        forward_voltage_drop_v,
        holding_current,
        on_conductance,
        off_conductance,
        initially_closed,
        gate_driver,
        antiparallel_diode,
        snubber,
        extended_fidelity,
    )
end

function PowerSemiconductorSwitch{TriacJunction}(
    a::Int,
    b::Int;
    gate_driver::PowerSemiconductorGateDriver=PowerSemiconductorGateDriver(),
    threshold_v::Real=0.0,
    forward_voltage_drop_v::Real=0.0,
    holding_current::Real=0.0,
    on_conductance::Real=1.0e3,
    off_conductance::Real=0.0,
    initially_closed::Bool=false,
    initial_direction::Integer=initially_closed ? 1 : 0,
    snubber::Union{Nothing,SeriesRCSnubber}=nothing,
    extended_fidelity::Union{Nothing,PowerSemiconductorExtendedFidelity}=nothing,
)
    direction = Int(initial_direction)
    direction in (-1, 0, 1) || throw(ArgumentError(
        "TRIAC initial direction must be -1, 0, or 1",
    ))
    initially_closed == (direction != 0) || throw(ArgumentError(
        "TRIAC initial direction must be nonzero exactly when initially closed",
    ))
    device = _power_semiconductor_switch(
        TriacJunction,
        a,
        b;
        threshold_v,
        forward_voltage_drop_v,
        holding_current,
        on_conductance,
        off_conductance,
        initially_closed,
        gate_driver,
        snubber,
        extended_fidelity,
    )
    device.conduction_direction = Int8(direction)
    device.last_history_current_a = initially_closed ?
        -device.on_conductance * device.forward_voltage_drop_v * direction : 0.0
    return device
end

function PowerSemiconductorSwitch{GateTurnOffThyristorJunction}(
    a::Int,
    b::Int;
    gate_driver::PowerSemiconductorGateDriver=PowerSemiconductorGateDriver(),
    threshold_v::Real=0.0,
    forward_voltage_drop_v::Real=0.0,
    holding_current::Real=0.0,
    on_conductance::Real=1.0e3,
    off_conductance::Real=0.0,
    initially_closed::Bool=false,
    interruptible_current_a::Real,
    turn_off_policy::Symbol=:refuse,
    antiparallel_diode::Union{Nothing,AntiparallelDiodeParameters}=nothing,
    snubber::Union{Nothing,SeriesRCSnubber}=nothing,
    extended_fidelity::Union{Nothing,PowerSemiconductorExtendedFidelity}=nothing,
)
    return _power_semiconductor_switch(
        GateTurnOffThyristorJunction,
        a,
        b;
        threshold_v,
        forward_voltage_drop_v,
        holding_current,
        on_conductance,
        off_conductance,
        initially_closed,
        gate_driver,
        antiparallel_diode,
        snubber,
        extended_fidelity,
        gate_turn_off_current_limit_a=interruptible_current_a,
        gate_turn_off_policy=turn_off_policy,
    )
end

function PowerSemiconductorSwitch{InsulatedGateBipolarTransistorJunction}(
    a::Int,
    b::Int;
    gate_driver::PowerSemiconductorGateDriver=PowerSemiconductorGateDriver(),
    threshold_v::Real=0.0,
    forward_voltage_drop_v::Real=0.0,
    holding_current::Real=0.0,
    on_conductance::Real=1.0e3,
    off_conductance::Real=0.0,
    initially_closed::Bool=false,
    antiparallel_diode::Union{Nothing,AntiparallelDiodeParameters}=nothing,
    snubber::Union{Nothing,SeriesRCSnubber}=nothing,
    extended_fidelity::Union{Nothing,PowerSemiconductorExtendedFidelity}=nothing,
)
    return _power_semiconductor_switch(
        InsulatedGateBipolarTransistorJunction,
        a,
        b;
        threshold_v,
        forward_voltage_drop_v,
        holding_current,
        on_conductance,
        off_conductance,
        initially_closed,
        gate_driver,
        antiparallel_diode,
        snubber,
        extended_fidelity,
    )
end

function PowerSemiconductorSwitch{MetalOxideSemiconductorFieldEffectTransistorJunction}(
    a::Int,
    b::Int;
    gate_driver::PowerSemiconductorGateDriver=PowerSemiconductorGateDriver(),
    on_conductance::Real=1.0e3,
    off_conductance::Real=0.0,
    initially_closed::Bool=false,
    antiparallel_diode::Union{Nothing,AntiparallelDiodeParameters}=nothing,
    snubber::Union{Nothing,SeriesRCSnubber}=nothing,
    extended_fidelity::Union{Nothing,PowerSemiconductorExtendedFidelity}=nothing,
)
    return _power_semiconductor_switch(
        MetalOxideSemiconductorFieldEffectTransistorJunction,
        a,
        b;
        threshold_v=0.0,
        forward_voltage_drop_v=0.0,
        holding_current=0.0,
        on_conductance,
        off_conductance,
        initially_closed,
        gate_driver,
        antiparallel_diode,
        snubber,
        extended_fidelity,
    )
end

function diode_next_closed(
    device::DiodeValveSwitch,
    voltage::Real,
    current::Real,
)::Bool
    return device.closed ? Float64(current) >= device.holding_current :
        Float64(voltage) >= power_semiconductor_forward_turn_on_voltage(device)
end

diode_conductance(device::DiodeValveSwitch)::Float64 =
    device.closed ? device.on_conductance : device.off_conductance

power_semiconductor_forward_turn_on_voltage(device::PowerSemiconductorSwitch) =
    max(
        device.threshold_v,
        device.forward_voltage_drop_v +
        device.holding_current / device.on_conductance,
    )

function power_semiconductor_reverse_turn_on_voltage(
    device::PowerSemiconductorSwitch,
)
    diode = something(device.antiparallel_diode)
    return diode.forward_voltage_v +
        diode.holding_current_a / diode.on_conductance_s
end

function power_semiconductor_event_localization!(device::PowerSemiconductorSwitch)
    device.event_localization_enabled = true
    return device
end

function _record_power_semiconductor_topology_transition!(
    device::PowerSemiconductorSwitch,
    previous_forward_state::Bool,
    previous_reverse_state::Bool,
    time_s::Real,
)
    if previous_forward_state != device.closed ||
       previous_reverse_state != device.reverse_diode_conducting
        device.topology_transition_count += 1
        device.last_transition_time_s = Float64(time_s)
    end
    return device
end

function _apply_power_semiconductor_gate_state!(
    device::PowerSemiconductorSwitch{DiodeJunction},
    _time_s::Float64,
)
    return device
end

function _apply_power_semiconductor_gate_state!(
    device::PowerSemiconductorSwitch{ThyristorJunction},
    time_s::Float64,
)
    previous_forward = device.closed
    previous_reverse = device.reverse_diode_conducting
    if device.gate_driver.applied_on &&
       !device.reverse_diode_conducting &&
       device.last_voltage >= power_semiconductor_forward_turn_on_voltage(device)
        device.closed = true
    end
    return _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
    )
end

function _apply_power_semiconductor_gate_state!(
    device::PowerSemiconductorSwitch{TriacJunction},
    time_s::Float64,
)
    previous_forward = device.closed
    previous_reverse = false
    if device.gate_driver.applied_on && !device.closed &&
       abs(device.last_voltage) >= power_semiconductor_forward_turn_on_voltage(device)
        device.closed = true
        device.conduction_direction = device.last_voltage < 0.0 ? Int8(-1) : Int8(1)
    end
    return _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
    )
end

function _apply_power_semiconductor_gate_state!(
    device::PowerSemiconductorSwitch{GateTurnOffThyristorJunction},
    time_s::Float64,
)
    previous_forward = device.closed
    previous_reverse = device.reverse_diode_conducting
    if device.gate_driver.applied_on && !device.closed &&
       !device.reverse_diode_conducting &&
       device.last_voltage >= power_semiconductor_forward_turn_on_voltage(device)
        device.closed = true
        device.conduction_direction = Int8(1)
        device.gate_turn_off_disposition = :not_requested
    elseif !device.gate_driver.applied_on && device.closed
        if abs(device.last_forward_current_a) <= device.gate_turn_off_current_limit_a
            _activate_power_semiconductor_tail!(device, time_s)
            device.closed = false
            device.conduction_direction = Int8(0)
            device.gate_turn_off_disposition = :accepted
        elseif device.gate_turn_off_policy === :remain_latched
            device.gate_turn_off_disposition = :remained_latched
        else
            throw(DomainError(
                device.last_forward_current_a,
                "GTO gate turn-off current exceeds its declared interruptible domain",
            ))
        end
    end
    return _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
    )
end

function _apply_power_semiconductor_gate_state!(
    device::PowerSemiconductorSwitch{InsulatedGateBipolarTransistorJunction},
    time_s::Float64,
)
    previous_forward = device.closed
    previous_reverse = device.reverse_diode_conducting
    if device.gate_driver.applied_on &&
       !device.reverse_diode_conducting &&
       device.last_voltage >= power_semiconductor_forward_turn_on_voltage(device)
        device.closed = true
        device.conduction_direction = Int8(1)
    elseif !device.gate_driver.applied_on
        device.closed && _activate_power_semiconductor_tail!(device, time_s)
        device.closed = false
        device.conduction_direction = Int8(0)
    end
    return _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
    )
end

function _apply_power_semiconductor_gate_state!(
    device::PowerSemiconductorSwitch{MetalOxideSemiconductorFieldEffectTransistorJunction},
    time_s::Float64,
)
    previous_forward = device.closed
    previous_reverse = device.reverse_diode_conducting
    device.closed = device.gate_driver.applied_on
    device.closed && (device.reverse_diode_conducting = false)
    return _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
    )
end

function request_power_semiconductor_gate!(
    device::PowerSemiconductorSwitch,
    commanded_on::Bool,
    time_s::Real,
    ;
    earliest_transition_time_s::Union{Nothing,Real}=nothing,
)
    driver = device.gate_driver
    driver === nothing && throw(ArgumentError("a diode junction has no gate command"))
    command_time = Float64(time_s)
    timing = power_semiconductor_gate_transition_timing(
        driver,
        commanded_on,
        command_time;
        earliest_transition_time_s,
    )
    if device isa PowerSemiconductorSwitch{GateTurnOffThyristorJunction} &&
       !commanded_on && device.closed &&
       abs(device.last_forward_current_a) > device.gate_turn_off_current_limit_a &&
       device.gate_turn_off_policy === :refuse
        device.gate_turn_off_disposition = :refused
        throw(DomainError(
            device.last_forward_current_a,
            "GTO gate turn-off current exceeds its declared interruptible domain",
        ))
    end
    commanded_on == driver.commanded_on && return false
    driver.commanded_on = commanded_on
    driver.last_command_time_s = command_time
    driver.command_count += 1
    if commanded_on == driver.applied_on
        driver.pending_state !== nothing && (driver.filtered_pulse_count += 1)
        driver.pending_state = nothing
        driver.pending_transition_time_s = Inf
        return false
    end
    transition_time = timing.transition_time_s
    driver.pending_state = commanded_on
    driver.pending_transition_time_s = transition_time
    if transition_time <= command_time + 16 * eps(max(1.0, command_time))
        apply_power_semiconductor_gate_transition!(device, command_time)
        return true
    end
    return false
end

function power_semiconductor_gate_transition_time(device::PowerSemiconductorSwitch)
    driver = device.gate_driver
    driver === nothing && return nothing
    driver.pending_state === nothing && return nothing
    return driver.pending_transition_time_s
end

function power_semiconductor_gate_transition_residual(
    device::PowerSemiconductorSwitch,
    time_s::Real,
)
    transition_time_s = power_semiconductor_gate_transition_time(device)
    transition_time_s === nothing && return nothing
    return transition_time_s - Float64(time_s)
end

"""Compare the unchanged gate deadline with the complete physical endpoint."""
function power_semiconductor_gate_transition_residual(
    device::PowerSemiconductorSwitch,
    time::ContinuousEvaluationTime,
)
    transition_time_s = power_semiconductor_gate_transition_time(device)
    transition_time_s === nothing && return nothing
    return setprecision(BigFloat, 106) do
        -continuous_time_elapsed(time, transition_time_s)
    end
end

function _apply_power_semiconductor_gate_event!(
    device::PowerSemiconductorSwitch,
    time_s::Real,
)
    fidelity = something(device.extended_fidelity)
    fidelity.pending_event_current_a = abs(device.last_forward_current_a)
    fidelity.pending_event_blocking_voltage_v = abs(device.last_voltage)
    transitioned = apply_power_semiconductor_gate_transition!(device, time_s)
    transitioned || return device
    event_energy_j = _deposit_power_semiconductor_switching_energy!(
        device,
        device.last_voltage,
        device.last_forward_current_a,
        Float64(time_s),
    )
    device.semiconductor_dissipated_energy_j += event_energy_j
    return device
end

function apply_power_semiconductor_gate_transition!(
    device::PowerSemiconductorSwitch,
    time_s::Real,
    ; defer_conduction::Bool=false,
)
    driver = device.gate_driver
    driver === nothing && throw(ArgumentError("a diode junction has no gate transition"))
    driver.pending_state === nothing && return false
    defer_conduction && !(power_semiconductor_kind(device) in (:igbt, :mosfet)) &&
        throw(ArgumentError(
            "deferred gate conduction requires a self-commutated bridge switch",
        ))
    transition_time = Float64(time_s)
    tolerance = 16 * eps(max(1.0, abs(transition_time)))
    transition_time + tolerance >= driver.pending_transition_time_s || throw(ArgumentError(
        "gate transition was requested before its delay/dead-time/minimum-pulse boundary",
    ))
    applied_on = something(driver.pending_state)
    driver.applied_on = applied_on
    driver.pending_state = nothing
    driver.pending_transition_time_s = Inf
    driver.transition_count += 1
    if applied_on
        driver.last_turn_on_time_s = transition_time
    else
        driver.last_turn_off_time_s = transition_time
    end
    if !defer_conduction
        _apply_power_semiconductor_gate_state!(device, transition_time)
    end
    return true
end

function power_semiconductor_forward_turn_on_residual(
    device::PowerSemiconductorSwitch{DiodeJunction},
)
    device.closed || device.reverse_diode_conducting ? nothing :
        device.last_voltage - power_semiconductor_forward_turn_on_voltage(device)
end

function power_semiconductor_forward_turn_on_residual(
    device::PowerSemiconductorSwitch{ThyristorJunction},
)
    driver = device.gate_driver
    !device.closed && !device.reverse_diode_conducting && driver.applied_on ?
        device.last_voltage - power_semiconductor_forward_turn_on_voltage(device) : nothing
end

function power_semiconductor_forward_turn_on_residual(
    device::PowerSemiconductorSwitch{TriacJunction},
)
    driver = device.gate_driver
    !device.closed && driver.applied_on ?
        abs(device.last_voltage) - power_semiconductor_forward_turn_on_voltage(device) : nothing
end

function power_semiconductor_forward_turn_on_residual(
    device::PowerSemiconductorSwitch{GateTurnOffThyristorJunction},
)
    driver = device.gate_driver
    !device.closed && !device.reverse_diode_conducting && driver.applied_on ?
        device.last_voltage - power_semiconductor_forward_turn_on_voltage(device) : nothing
end

function power_semiconductor_forward_turn_on_residual(
    device::PowerSemiconductorSwitch{InsulatedGateBipolarTransistorJunction},
)
    driver = device.gate_driver
    !device.closed && !device.reverse_diode_conducting && driver.applied_on ?
        device.last_voltage - power_semiconductor_forward_turn_on_voltage(device) : nothing
end

power_semiconductor_forward_turn_on_residual(
    ::PowerSemiconductorSwitch{MetalOxideSemiconductorFieldEffectTransistorJunction},
) = nothing

function power_semiconductor_forward_extinction_residual(
    device::PowerSemiconductorSwitch,
)
    device.closed ? device.last_forward_current_a - device.holding_current : nothing
end

function power_semiconductor_forward_extinction_residual(
    device::PowerSemiconductorSwitch{TriacJunction},
)
    device.closed ?
        device.conduction_direction * device.last_forward_current_a - device.holding_current :
        nothing
end

"""Linearize a localized forward-conduction event on its pre-event chart.

The accepted channel law owns holding-current extinction; the open-device
voltage law owns forward turn-on. These slopes are with respect to terminal
voltage only. A zero-voltage triac corner has no unique first derivative and
must be rejected by a tangent caller rather than assigned one arbitrarily.
"""
function power_semiconductor_forward_event_voltage_linearization(
    device::PowerSemiconductorSwitch,
    surface_name::Symbol,
)
    if surface_name === :forward_holding_extinction
        residual = power_semiconductor_forward_extinction_residual(device)
        residual === nothing && throw(ArgumentError(
            "forward-extinction linearization requires a closed device",
        ))
        _, conductance = _power_semiconductor_accepted_boundary_channel_current_jacobian(
            device, device.last_voltage)
        direction = device isa PowerSemiconductorSwitch{TriacJunction} ?
            Float64(device.conduction_direction) : 1.0
        slope = direction * conductance
    elseif surface_name === :forward_turn_on
        residual = power_semiconductor_forward_turn_on_residual(device)
        residual === nothing && throw(ArgumentError(
            "forward-turn-on linearization requires an eligible open device",
        ))
        if device isa PowerSemiconductorSwitch{TriacJunction}
            iszero(device.last_voltage) && throw(DomainError(
                device.last_voltage,
                "zero-voltage triac turn-on has no unique tangent",
            ))
            slope = sign(device.last_voltage)
        else
            slope = 1.0
        end
    else
        throw(ArgumentError("unsupported forward event surface $surface_name"))
    end
    isfinite(residual) && isfinite(slope) && !iszero(slope) || throw(
        DomainError((residual, slope),
            "forward event surface lost a finite transversal voltage slope"),
    )
    return (; residual, voltage_slope=slope)
end

power_semiconductor_forward_extinction_residual(
    ::PowerSemiconductorSwitch{MetalOxideSemiconductorFieldEffectTransistorJunction},
) = nothing

function power_semiconductor_reverse_turn_on_residual(
    device::PowerSemiconductorSwitch,
)
    diode = device.antiparallel_diode
    diode === nothing && return nothing
    device.closed || device.reverse_diode_conducting ? nothing :
        -device.last_voltage - power_semiconductor_reverse_turn_on_voltage(device)
end

function power_semiconductor_reverse_extinction_residual(
    device::PowerSemiconductorSwitch,
)
    diode = device.antiparallel_diode
    diode === nothing && return nothing
    device.reverse_diode_conducting ?
        device.last_reverse_diode_current_a - diode.holding_current_a : nothing
end

"""Return the accepted stored-charge exhaustion surface without mutating device state."""
function power_semiconductor_recovery_zero_residual(
    device::PowerSemiconductorSwitch,
)
    fidelity = device.extended_fidelity
    fidelity === nothing && return nothing
    recovery = fidelity.recovered_charge
    recovery === nothing && return nothing
    if recovery.stored_charge_c > 0.0 ||
       recovery.previous_stored_charge_c > 0.0 &&
       recovery.last_recovery_zero_time_s < fidelity.candidate_time_s
        return fidelity.candidate_recovery_charge_c
    end
    return nothing
end

"""Return the within-companion exhaustion boundary after an endpoint recovery probe."""
function power_semiconductor_recovery_zero_time(device::PowerSemiconductorSwitch)
    fidelity = device.extended_fidelity
    fidelity === nothing && return nothing
    recovery = fidelity.recovered_charge
    recovery === nothing && return nothing
    return nothing
end

"""Accept one localized stored-charge exhaustion transition exactly once."""
function apply_power_semiconductor_recovery_zero!(
    device::PowerSemiconductorSwitch,
    time_s::Union{Real,ContinuousEvaluationTime},
)
    fidelity = device.extended_fidelity
    fidelity === nothing && throw(ArgumentError(
        "baseline power semiconductor has no recovered-charge transition",
    ))
    recovery = fidelity.recovered_charge
    recovery === nothing && throw(ArgumentError(
        "power semiconductor has no recovered-charge component",
    ))
    time = time_s isa ContinuousEvaluationTime ? time_s.high_s : Float64(time_s)
    isfinite(time_s) || throw(ArgumentError("recovery-zero event time must be finite"))
    recovery.last_recovery_zero_time_s < time || throw(ArgumentError(
        "recovery-zero event must advance beyond the last accepted occurrence",
    ))
    localized_charge_c = fidelity.candidate_recovery_charge_c
    tolerance_c = _power_semiconductor_recovery_zero_tolerance(
        recovery.previous_stored_charge_c,
        recovery.last_recovery_current_a,
        fidelity.candidate_step_s,
    )
    abs(localized_charge_c) <= tolerance_c || throw(DomainError(
        localized_charge_c,
        "recovery-zero transition requires localized exhausted charge",
    ))
    # Validate full chronology and compute physical duration before any reset.
    duration = if isfinite(recovery.recovery_start_time_s)
        exact_age = _semiconductor_time_exact(time_s)-_semiconductor_time_exact(recovery.recovery_start_time_s)
        exact_age >= 0 || throw(ArgumentError("recovery-zero physical origin exceeds event time"))
        Float64(_semiconductor_history_elapsed(time_s,recovery.recovery_start_time_s))
    else
        0.0
    end
    recovery.stored_charge_c = 0.0
    recovery.previous_stored_charge_c = 0.0
    recovery.recovery_active = false
    recovery.last_recovery_current_a = 0.0
    recovery.last_recovery_duration_s = duration
    recovery.recovery_zero_event_count += 1
    recovery.last_recovery_zero_time_s = time
    return device
end

"""Return the analytic tail-current cutoff surface at one candidate time."""
function power_semiconductor_tail_cutoff_residual(
    device::PowerSemiconductorSwitch,
    time_s::Real,
)
    fidelity = device.extended_fidelity
    fidelity === nothing && return nothing
    tail = fidelity.turn_off_tail
    tail === nothing && return nothing
    time = Float64(time_s)
    isfinite(time) || throw(ArgumentError("tail-cutoff evaluation time must be finite"))
    tail.initial_current_a > tail.cutoff_current_a || return nothing
    isfinite(tail.turn_off_time_s) || return nothing
    tail.last_cutoff_time_s < _semiconductor_history_event_time(tail.turn_off_time_s) || return nothing
    elapsed_time_s = max(0.0, _semiconductor_tail_elapsed(time, tail.turn_off_time_s))
    return tail.initial_current_a * exp(-elapsed_time_s / tail.decay_time_s) -
        tail.cutoff_current_a
end

function power_semiconductor_tail_cutoff_residual(device::PowerSemiconductorSwitch,
    time::ContinuousEvaluationTime)
    fidelity = device.extended_fidelity
    fidelity===nothing && return nothing
    tail = fidelity.turn_off_tail
    tail===nothing && return nothing
    tail.initial_current_a>tail.cutoff_current_a && isfinite(tail.turn_off_time_s) &&
        tail.last_cutoff_time_s<_semiconductor_history_event_time(tail.turn_off_time_s) || return nothing
    setprecision(BigFloat,106) do
        elapsed = max(BigFloat(0),_semiconductor_tail_elapsed(time,tail.turn_off_time_s))
        BigFloat(tail.initial_current_a)*exp(-elapsed/BigFloat(tail.decay_time_s))-BigFloat(tail.cutoff_current_a)
    end
end

"""Return the exact positive-cutoff time for one active exponential tail."""
function power_semiconductor_tail_cutoff_time(device::PowerSemiconductorSwitch)
    fidelity = device.extended_fidelity
    fidelity === nothing && return nothing
    tail = fidelity.turn_off_tail
    tail === nothing && return nothing
    tail.active || return nothing
    tail.cutoff_current_a > 0.0 || return nothing
    tail.initial_current_a > tail.cutoff_current_a || return nothing
    tail.last_cutoff_time_s < _semiconductor_history_event_time(tail.turn_off_time_s) || return nothing
    duration = tail.decay_time_s * log(tail.initial_current_a / tail.cutoff_current_a)
    return Float64(_semiconductor_time_storage(
        _semiconductor_time_exact(tail.turn_off_time_s) + _semiconductor_time_exact(duration)))
end

"""Accept one localized exponential-tail cutoff exactly once."""
function apply_power_semiconductor_tail_cutoff!(
    device::PowerSemiconductorSwitch,
    time_s::Union{Real,ContinuousEvaluationTime},
)
    fidelity = device.extended_fidelity
    fidelity === nothing && throw(ArgumentError(
        "baseline power semiconductor has no tail-cutoff transition",
    ))
    tail = fidelity.turn_off_tail
    tail === nothing && throw(ArgumentError(
        "power semiconductor has no turn-off-tail component",
    ))
    time = time_s isa ContinuousEvaluationTime ? time_s.high_s : Float64(time_s)
    isfinite(time) || throw(ArgumentError("tail-cutoff event time must be finite"))
    tail.last_cutoff_time_s < time || throw(ArgumentError(
        "tail-cutoff event must advance beyond the last accepted occurrence",
    ))
    residual = power_semiconductor_tail_cutoff_residual(device, time_s)
    residual === nothing && throw(ArgumentError(
        "turn-off tail has no pending localized cutoff",
    ))
    root_time_tolerance_s = 1.0e-12
    tolerance = max(
        8192.0 * eps(Float64) * max(1.0, abs(tail.cutoff_current_a)),
        abs(tail.cutoff_current_a) *
            expm1(root_time_tolerance_s / tail.decay_time_s),
    )
    abs(residual) <= tolerance || throw(DomainError(
        residual,
        "tail-cutoff transition requires its localized current surface",
    ))
    tail.active = false
    tail.current_a = 0.0
    tail.initial_current_a = 0.0
    tail.last_duration_s = Float64(max(0.0, _semiconductor_tail_elapsed(time_s, tail.turn_off_time_s)))
    tail.cutoff_event_count += 1
    tail.last_cutoff_time_s = time
    return device
end

const _POWER_SEMICONDUCTOR_RECOVERY_ZERO_EVENT = NonlinearDeviceEventSurface(
    :stored_charge_exhaustion,
    (device, _time_s) -> power_semiconductor_recovery_zero_residual(device),
    apply_power_semiconductor_recovery_zero!;
    direction=:falling,
    priority=10,
    candidate_time=power_semiconductor_recovery_zero_time,
)

"""Localize recovered charge at least as tightly as its physical reset guard.

The reset allowance always includes `8192eps(Float64)` coulombs; using this
conservative floor prevents a generic event policy from accepting a root that
the unchanged physical transition must reject.
"""
function NonlinearNetwork.nonlinear_device_event_value_tolerance(
    surface::NonlinearDeviceEventSurface,
    device::PowerSemiconductorSwitch,
    policy_value_tolerance::Real,
)
    tolerance = invoke(
        NonlinearNetwork.nonlinear_device_event_value_tolerance,
        Tuple{NonlinearDeviceEventSurface,Any,Real},
        surface,
        device,
        policy_value_tolerance,
    )
    return surface === _POWER_SEMICONDUCTOR_RECOVERY_ZERO_EVENT ?
        min(tolerance, 8192eps(Float64)) : tolerance
end

const _POWER_SEMICONDUCTOR_GATE_TRANSITION_EVENT = NonlinearDeviceEventSurface(
    :gate_transition,
    power_semiconductor_gate_transition_residual,
    _apply_power_semiconductor_gate_event!;
    direction=:falling,
    priority=0,
    topology_invalidating=true,
    candidate_time=power_semiconductor_gate_transition_time,
)

const _POWER_SEMICONDUCTOR_FORWARD_EXTINCTION_EVENT = NonlinearDeviceEventSurface(
    :forward_holding_extinction,
    (device, _time_s) -> power_semiconductor_forward_extinction_residual(device),
    (device, time_s) -> apply_power_semiconductor_forward_extinction!(device, time_s);
    direction=:falling,
    priority=5,
    topology_invalidating=true,
)

const _POWER_SEMICONDUCTOR_FORWARD_TURN_ON_EVENT = NonlinearDeviceEventSurface(
    :forward_turn_on,
    (device, _time_s) -> power_semiconductor_forward_turn_on_residual(device),
    (device, time_s) -> apply_power_semiconductor_forward_turn_on!(device, time_s);
    direction=:rising,
    priority=6,
    topology_invalidating=true,
)

const _POWER_SEMICONDUCTOR_TAIL_CUTOFF_EVENT = NonlinearDeviceEventSurface(
    :turn_off_tail_cutoff,
    power_semiconductor_tail_cutoff_residual,
    apply_power_semiconductor_tail_cutoff!;
    direction=:falling,
    priority=20,
    candidate_time=power_semiconductor_tail_cutoff_time,
)

function nonlinear_device_event_surfaces(device::PowerSemiconductorSwitch)
    fidelity = device.extended_fidelity
    fidelity === nothing && return ()
    gate_event = power_semiconductor_gate_transition_time(device) === nothing ? () :
        (_POWER_SEMICONDUCTOR_GATE_TRANSITION_EVENT,)
    extinction_event = device.event_localization_enabled && device.closed &&
        !(device isa PowerSemiconductorSwitch{MetalOxideSemiconductorFieldEffectTransistorJunction}) ?
        (_POWER_SEMICONDUCTOR_FORWARD_EXTINCTION_EVENT,) : ()
    turn_on_event = device.event_localization_enabled && !device.closed &&
        power_semiconductor_forward_turn_on_residual(device) !== nothing ?
        (_POWER_SEMICONDUCTOR_FORWARD_TURN_ON_EVENT,) : ()
    recovery_event = fidelity.recovered_charge === nothing ? () :
        (_POWER_SEMICONDUCTOR_RECOVERY_ZERO_EVENT,)
    tail_event = fidelity.turn_off_tail === nothing ? () :
        (_POWER_SEMICONDUCTOR_TAIL_CUTOFF_EVENT,)
    return (gate_event..., extinction_event..., turn_on_event...,
        recovery_event..., tail_event...)
end

function apply_power_semiconductor_forward_turn_on!(
    device::PowerSemiconductorSwitch,
    time_s::Real,
)
    previous_forward = device.closed
    previous_reverse = device.reverse_diode_conducting
    device.reverse_diode_conducting && throw(ArgumentError(
        "forward conduction cannot begin while the antiparallel diode conducts",
    ))
    device.closed = true
    device.conduction_direction = device isa PowerSemiconductorSwitch{TriacJunction} ?
        (device.last_voltage < 0.0 ? Int8(-1) : Int8(1)) : Int8(1)
    _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
    )
    return device
end

function apply_power_semiconductor_forward_turn_on!(device::PowerSemiconductorSwitch,
    time::ContinuousEvaluationTime)
    # Turn-on reads the accepted voltage/conduction owner, not elapsed time.
    # Its stored transition timestamp is the existing discrete bookkeeping ID.
    apply_power_semiconductor_forward_turn_on!(device, time.high_s)
end

function apply_power_semiconductor_forward_extinction!(
    device::PowerSemiconductorSwitch,
    time_s::Union{Real,ContinuousEvaluationTime},
)
    previous_forward = device.closed
    previous_reverse = device.reverse_diode_conducting
    _activate_power_semiconductor_tail!(device, time_s)
    device.closed = false
    device.conduction_direction = Int8(0)
    _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s isa ContinuousEvaluationTime ? time_s.high_s : Float64(time_s),
    )
    return device
end

function apply_power_semiconductor_reverse_turn_on!(
    device::PowerSemiconductorSwitch,
    time_s::Real,
)
    device.antiparallel_diode === nothing && throw(ArgumentError(
        "power semiconductor has no antiparallel diode",
    ))
    previous_forward = device.closed
    previous_reverse = device.reverse_diode_conducting
    device.closed && throw(ArgumentError(
        "antiparallel conduction cannot begin while the forward path conducts",
    ))
    device.reverse_diode_conducting = true
    _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
    )
    return device
end

function apply_power_semiconductor_reverse_extinction!(
    device::PowerSemiconductorSwitch,
    time_s::Real,
)
    previous_forward = device.closed
    previous_reverse = device.reverse_diode_conducting
    device.reverse_diode_conducting = false
    _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
    )
    return device
end

# These transition owners store Float64 occurrence timestamps, as does the
# existing forward-turn-on overload. Continuous network evaluation retains
# the endpoint low part; it must not make these discrete callbacks unavailable
# or change their physical reset guards and accepted transition bookkeeping.
apply_power_semiconductor_reverse_turn_on!(device::PowerSemiconductorSwitch,
    time::ContinuousEvaluationTime) =
    apply_power_semiconductor_reverse_turn_on!(device, time.high_s)
apply_power_semiconductor_reverse_extinction!(device::PowerSemiconductorSwitch,
    time::ContinuousEvaluationTime) =
    apply_power_semiconductor_reverse_extinction!(device, time.high_s)

function _series_rc_snubber_companion(snubber::SeriesRCSnubber, dt_s::Float64)
    dt_s > 0.0 && isfinite(dt_s) ||
        throw(ArgumentError("snubber timestep must be finite and positive"))
    capacitive_resistance = dt_s / (2.0 * snubber.capacitance_f)
    conductance = inv(snubber.resistance_ohm + capacitive_resistance)
    history_current = conductance * (
        -capacitive_resistance * snubber.previous_current_a -
        snubber.capacitor_voltage_v
    )
    return conductance, history_current
end

function _series_rc_snubber_backward_euler_companion(
    snubber::SeriesRCSnubber,
    dt_s::Float64,
)
    dt_s > 0.0 && isfinite(dt_s) || throw(ArgumentError(
        "snubber timestep must be finite and positive",
    ))
    capacitive_resistance = dt_s / snubber.capacitance_f
    conductance = inv(snubber.resistance_ohm + capacitive_resistance)
    history_current = -conductance * snubber.capacitor_voltage_v
    return conductance, history_current
end

function _stamp_power_semiconductor_gate_due!(
    device::PowerSemiconductorSwitch,
    time_s::Float64,
)
    device.event_localization_enabled && return device
    transition_time = power_semiconductor_gate_transition_time(device)
    transition_time === nothing && return device
    transition_time <= time_s + 16 * eps(max(1.0, abs(time_s))) &&
        apply_power_semiconductor_gate_transition!(device, time_s)
    return device
end

"""Return the exact accepted piecewise-linear channel companion for a semiconductor's discrete conduction state."""
function _power_semiconductor_boundary_channel_companion(
    device::PowerSemiconductorSwitch,
)
    if device.closed
        direction = device isa PowerSemiconductorSwitch{TriacJunction} ?
            Float64(device.conduction_direction) : 1.0
        return (
            conductance_s=device.on_conductance,
            history_current_a=-device.on_conductance *
                device.forward_voltage_drop_v * direction,
        )
    elseif device.reverse_diode_conducting
        diode = something(device.antiparallel_diode)
        return (
            conductance_s=diode.on_conductance_s,
            history_current_a=diode.on_conductance_s * diode.forward_voltage_v,
        )
    end
    return (
        conductance_s=device.off_conductance,
        history_current_a=0.0,
    )
end

function stamp!(
    admittance::AbstractMatrix{<:Real},
    rhs::AbstractVector{<:Real},
    device::PowerSemiconductorSwitch,
    time_s::Float64,
    dt_s::Float64,
)
    _stamp_power_semiconductor_gate_due!(device, time_s)
    device.last_evaluation_time_s = time_s
    companion = _power_semiconductor_boundary_channel_companion(device)
    device.last_conductance = companion.conductance_s
    device.last_history_current_a = companion.history_current_a
    stamp_conductance!(admittance, device.a, device.b, companion.conductance_s)
    stamp_history_current!(rhs, device.a, device.b, companion.history_current_a)
    if device.snubber !== nothing
        snubber_conductance, snubber_history_current =
            _series_rc_snubber_companion(device.snubber, dt_s)
        stamp_conductance!(admittance, device.a, device.b, snubber_conductance)
        stamp_history_current!(rhs, device.a, device.b, snubber_history_current)
    end
    return nothing
end

Branches.backward_euler_companion_supported(::PowerSemiconductorSwitch) = true

function stamp!(
    admittance::AbstractMatrix{<:Real},
    rhs::AbstractVector{<:Real},
    device::PowerSemiconductorSwitch,
    time_s::Float64,
    dt_s::Float64,
    ::Val{TrapezoidalCompanion},
)
    return stamp!(admittance, rhs, device, time_s, dt_s)
end

function stamp!(
    admittance::AbstractMatrix{<:Real},
    rhs::AbstractVector{<:Real},
    device::PowerSemiconductorSwitch,
    time_s::Float64,
    dt_s::Float64,
    ::Val{BackwardEulerCompanion},
)
    _stamp_power_semiconductor_gate_due!(device, time_s)
    device.last_evaluation_time_s = time_s
    companion = _power_semiconductor_boundary_channel_companion(device)
    device.last_conductance = companion.conductance_s
    device.last_history_current_a = companion.history_current_a
    stamp_conductance!(admittance, device.a, device.b, companion.conductance_s)
    stamp_history_current!(rhs, device.a, device.b, companion.history_current_a)
    if device.snubber !== nothing
        snubber_conductance, snubber_history_current =
            _series_rc_snubber_backward_euler_companion(device.snubber, dt_s)
        stamp_conductance!(admittance, device.a, device.b, snubber_conductance)
        stamp_history_current!(rhs, device.a, device.b, snubber_history_current)
    end
    return nothing
end

function _update_series_rc_snubber!(
    snubber::SeriesRCSnubber,
    terminal_voltage_v::Float64,
    dt_s::Float64,
    ;
    accepted_current_a::Union{Nothing,Float64}=nothing,
)
    previous_current = snubber.previous_current_a
    previous_loss = snubber.last_resistor_loss_w
    conductance, history_current = _series_rc_snubber_companion(snubber, dt_s)
    current = accepted_current_a === nothing ?
        conductance * terminal_voltage_v + history_current : accepted_current_a
    capacitor_voltage = snubber.capacitor_voltage_v +
        dt_s / (2.0 * snubber.capacitance_f) * (current + previous_current)
    resistor_loss = snubber.resistance_ohm * current^2
    snubber.dissipated_energy_j += 0.5 * dt_s * (previous_loss + resistor_loss)
    snubber.previous_current_a = current
    snubber.capacitor_voltage_v = capacitor_voltage
    snubber.last_branch_voltage_v = terminal_voltage_v
    snubber.last_current_a = current
    snubber.last_resistor_loss_w = resistor_loss
    return snubber
end

function _update_series_rc_snubber_backward_euler!(
    snubber::SeriesRCSnubber,
    terminal_voltage_v::Float64,
    dt_s::Float64,
    ;
    accepted_current_a::Union{Nothing,Float64}=nothing,
)
    conductance, history_current =
        _series_rc_snubber_backward_euler_companion(snubber, dt_s)
    current = accepted_current_a === nothing ?
        conductance * terminal_voltage_v + history_current : accepted_current_a
    capacitor_voltage = snubber.capacitor_voltage_v +
        dt_s * current / snubber.capacitance_f
    resistor_loss = snubber.resistance_ohm * current^2
    snubber.dissipated_energy_j += dt_s * resistor_loss
    snubber.previous_current_a = current
    snubber.capacitor_voltage_v = capacitor_voltage
    snubber.last_branch_voltage_v = terminal_voltage_v
    snubber.last_current_a = current
    snubber.last_resistor_loss_w = resistor_loss
    return snubber
end

function _update_power_semiconductor_sampled_state!(device::DiodeValveSwitch)
    device.closed = diode_next_closed(
        device,
        device.last_voltage,
        device.last_forward_current_a,
    )
    return device
end

function _update_power_semiconductor_sampled_state!(
    device::PowerSemiconductorSwitch{ThyristorJunction},
)
    if device.closed
        device.last_forward_current_a < device.holding_current &&
            (device.closed = false)
    elseif device.gate_driver.applied_on &&
           !device.reverse_diode_conducting &&
           device.last_voltage >= power_semiconductor_forward_turn_on_voltage(device)
        device.closed = true
    end
    return device
end

function _update_power_semiconductor_sampled_state!(
    device::PowerSemiconductorSwitch{TriacJunction},
)
    if device.closed
        device.conduction_direction * device.last_forward_current_a < device.holding_current &&
            (device.closed = false; device.conduction_direction = Int8(0))
    elseif device.gate_driver.applied_on &&
           abs(device.last_voltage) >= power_semiconductor_forward_turn_on_voltage(device)
        device.closed = true
        device.conduction_direction = device.last_voltage < 0.0 ? Int8(-1) : Int8(1)
    end
    return device
end

function _update_power_semiconductor_sampled_state!(
    device::PowerSemiconductorSwitch{GateTurnOffThyristorJunction},
    evaluation_time::Union{Real,ContinuousEvaluationTime}=device.last_evaluation_time_s,
)
    if !device.gate_driver.applied_on && device.closed
        if abs(device.last_forward_current_a) <= device.gate_turn_off_current_limit_a
            _activate_power_semiconductor_tail!(device, evaluation_time)
            device.closed = false
            device.conduction_direction = Int8(0)
            device.gate_turn_off_disposition = :accepted
        elseif device.gate_turn_off_policy === :remain_latched
            device.gate_turn_off_disposition = :remained_latched
        else
            throw(DomainError(
                device.last_forward_current_a,
                "GTO gate turn-off current exceeds its declared interruptible domain",
            ))
        end
    elseif device.closed
        device.last_forward_current_a < device.holding_current &&
            (device.closed = false; device.conduction_direction = Int8(0))
    elseif device.gate_driver.applied_on && !device.reverse_diode_conducting &&
           device.last_voltage >= power_semiconductor_forward_turn_on_voltage(device)
        device.closed = true
        device.conduction_direction = Int8(1)
    end
    return device
end

function _update_power_semiconductor_sampled_state!(
    device::PowerSemiconductorSwitch{InsulatedGateBipolarTransistorJunction},
)
    if !device.gate_driver.applied_on
        device.closed = false
    elseif device.closed
        device.last_forward_current_a < device.holding_current &&
            (device.closed = false)
    elseif !device.reverse_diode_conducting &&
           device.last_voltage >= power_semiconductor_forward_turn_on_voltage(device)
        device.closed = true
    end
    return device
end

function _update_power_semiconductor_sampled_state!(
    device::PowerSemiconductorSwitch{MetalOxideSemiconductorFieldEffectTransistorJunction},
)
    device.closed = device.gate_driver.applied_on
    device.closed && (device.reverse_diode_conducting = false)
    return device
end

function _update_power_semiconductor_reverse_state!(device::PowerSemiconductorSwitch)
    diode = device.antiparallel_diode
    diode === nothing && return device
    if device.reverse_diode_conducting
        device.last_reverse_diode_current_a < diode.holding_current_a &&
            (device.reverse_diode_conducting = false)
    elseif !device.closed &&
           -device.last_voltage >= power_semiconductor_reverse_turn_on_voltage(device)
        device.reverse_diode_conducting = true
    end
    return device
end

function update!(
    device::PowerSemiconductorSwitch,
    voltages::AbstractVector{Float64},
    dt_s::Float64,
)
    terminal_voltage = Branches.branch_voltage(voltages, device.a, device.b)
    semiconductor_current = device.last_conductance * terminal_voltage +
        device.last_history_current_a
    forward_current = device.closed ? semiconductor_current : 0.0
    reverse_current = device.reverse_diode_conducting ? -semiconductor_current : 0.0
    snubber_current = 0.0
    if device.snubber !== nothing
        _update_series_rc_snubber!(device.snubber, terminal_voltage, dt_s)
        snubber_current = device.snubber.last_current_a
    end
    semiconductor_loss = max(0.0, terminal_voltage * semiconductor_current)
    device.semiconductor_dissipated_energy_j += 0.5 * dt_s * (
        device.previous_semiconductor_loss_w + semiconductor_loss
    )
    device.last_voltage = terminal_voltage
    device.last_forward_current_a = forward_current
    device.last_reverse_diode_current_a = reverse_current
    device.last_snubber_current_a = snubber_current
    device.last_current = semiconductor_current + snubber_current
    device.last_semiconductor_loss_w = semiconductor_loss
    device.previous_semiconductor_loss_w = semiconductor_loss
    if !device.event_localization_enabled
        previous_forward = device.closed
        previous_reverse = device.reverse_diode_conducting
        _update_power_semiconductor_sampled_state!(device)
        _update_power_semiconductor_reverse_state!(device)
        _record_power_semiconductor_topology_transition!(
            device,
            previous_forward,
            previous_reverse,
            device.last_evaluation_time_s,
        )
    end
    return nothing
end


nonlinear_device_formulation(::PowerSemiconductorSwitch) = PhysicalConstitutiveCurrent

function nonlinear_device_provenance(device::PowerSemiconductorSwitch)
    fidelity = device.extended_fidelity
    fidelity === nothing && throw(ArgumentError(
        "only a power semiconductor with explicit extended fidelity may enter the nonlinear-device owner list",
    ))
    return fidelity.provenance
end

nonlinear_terminal_nodes(device::PowerSemiconductorSwitch) = (device.a, device.b)

function _power_semiconductor_candidate_method(method::Symbol)
    if method === :TrapezoidalCompanion || method === :trapezoidal
        return :trapezoidal
    elseif method === :BackwardEulerCompanion || method === :backward_euler
        return :backward_euler
    end
    throw(ArgumentError(
        "extended semiconductor companion method must be trapezoidal or backward Euler",
    ))
end

function prepare_nonlinear_device_step!(
    device::PowerSemiconductorSwitch,
    time_s::Float64,
    step_s::Float64,
    companion_method::Symbol,
)
    fidelity = device.extended_fidelity
    fidelity === nothing && throw(ArgumentError(
        "baseline power semiconductor is linear and cannot be prepared as a nonlinear device",
    ))
    isfinite(time_s) || throw(ArgumentError(
        "extended semiconductor candidate time must be finite",
    ))
    isfinite(step_s) && step_s > 0.0 || throw(ArgumentError(
        "extended semiconductor candidate step must be finite and positive",
    ))
    method = _power_semiconductor_candidate_method(companion_method)
    fidelity.candidate_time_s = time_s
    fidelity.candidate_step_s = step_s
    fidelity.candidate_method = method
    fidelity.candidate_prepared = true
    if device.topology_transition_count > fidelity.accepted_topology_transition_count
        fidelity.pending_event_current_a = abs(device.last_forward_current_a)
        fidelity.pending_event_blocking_voltage_v = abs(device.last_voltage)
    end
    return nothing
end

function _power_semiconductor_terminal_voltage(
    terminal_voltage_v::AbstractVector{Float64},
)
    length(terminal_voltage_v) == 2 || throw(DimensionMismatch(
        "power semiconductor requires two terminal voltages",
    ))
    voltage_a_v = terminal_voltage_v[1]
    voltage_b_v = terminal_voltage_v[2]
    isfinite(voltage_a_v) && isfinite(voltage_b_v) || throw(ArgumentError(
        "power-semiconductor terminal voltages must be finite",
    ))
    return voltage_a_v - voltage_b_v
end

function _power_semiconductor_accepted_boundary_channel_current_jacobian(
    device::PowerSemiconductorSwitch,
    voltage_v::Float64,
)
    if device.closed
        direction = device isa PowerSemiconductorSwitch{TriacJunction} ?
            Float64(device.conduction_direction) : 1.0
        return device.on_conductance * (
            voltage_v - direction * device.forward_voltage_drop_v
        ), device.on_conductance
    elseif device.reverse_diode_conducting
        diode = something(device.antiparallel_diode)
        return diode.on_conductance_s * (
            voltage_v + diode.forward_voltage_v
        ), diode.on_conductance_s
    end
    fidelity = device.extended_fidelity
    recovery = fidelity === nothing ? nothing : fidelity.recovered_charge
    if recovery !== nothing && recovery.stored_charge_c > 0.0 && voltage_v < 0.0
        return device.on_conductance * voltage_v, device.on_conductance
    end
    return device.off_conductance * voltage_v, device.off_conductance
end

function _power_semiconductor_channel_current_jacobian(
    device::PowerSemiconductorSwitch,
    voltage_v::Float64,
)
    fidelity = device.extended_fidelity
    recovery = fidelity === nothing ? nothing : fidelity.recovered_charge
    recovery_conducting = !device.closed &&
        !device.reverse_diode_conducting && recovery !== nothing &&
        recovery.stored_charge_c > 0.0 && voltage_v < 0.0
    recovery_conducting && !fidelity.candidate_prepared && throw(ArgumentError(
        "recovered-charge evaluation requires one prepared candidate step",
    ))
    return _power_semiconductor_accepted_boundary_channel_current_jacobian(
        device,
        voltage_v,
    )
end

function _power_semiconductor_dynamic_current_jacobian(
    fidelity::PowerSemiconductorExtendedFidelity,
    voltage_v::Float64,
    time_s::Float64,
    localized_tail::Bool=false,
)
    fidelity.candidate_prepared || throw(ArgumentError(
        "extended semiconductor must be prepared for an exact candidate step before evaluation",
    ))
    tolerance = 16 * eps(max(1.0, abs(time_s), abs(fidelity.candidate_time_s)))
    abs(time_s - fidelity.candidate_time_s) <= tolerance || throw(ArgumentError(
        "extended semiconductor candidate time does not match the nonlinear evaluation time",
    ))
    step_s = fidelity.candidate_step_s
    current_a = 0.0
    jacobian_s = 0.0
    junction = fidelity.junction_charge
    if junction !== nothing
        charge, capacitance = _power_semiconductor_junction_charge_capacitance(
            junction,
            voltage_v,
        )
        if fidelity.candidate_method === :trapezoidal
            current_a += 2.0 * (charge - junction.previous_charge_c) / step_s -
                junction.last_displacement_current_a
            jacobian_s += 2.0 * capacitance / step_s
        else
            current_a += (charge - junction.previous_charge_c) / step_s
            jacobian_s += capacitance / step_s
        end
    end
    tail = fidelity.turn_off_tail
    if tail !== nothing && tail.active
        duration_s = max(0.0, _semiconductor_tail_elapsed(time_s, tail.turn_off_time_s))
        tail_current = tail.initial_current_a * exp(-duration_s / tail.decay_time_s)
        # Localized integration accepts the incoming chart before the event
        # owner resets it. Endpoint roundoff must not select a second chart.
        (localized_tail || tail_current > tail.cutoff_current_a) &&
            (current_a += tail_current)
    end
    return current_a, jacobian_s
end

function _power_semiconductor_snubber_current_jacobian(
    device::PowerSemiconductorSwitch,
    voltage_v::Float64,
)
    snubber = device.snubber
    snubber === nothing && return 0.0, 0.0
    fidelity = something(device.extended_fidelity)
    conductance, history_current = if fidelity.candidate_method === :trapezoidal
        _series_rc_snubber_companion(snubber, fidelity.candidate_step_s)
    else
        _series_rc_snubber_backward_euler_companion(snubber, fidelity.candidate_step_s)
    end
    return conductance * voltage_v + history_current, conductance
end

function nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    device::PowerSemiconductorSwitch,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::Float64,
)
    length(terminal_current_a) == 2 || throw(DimensionMismatch(
        "power semiconductor requires two terminal currents",
    ))
    size(terminal_jacobian_s) == (2, 2) || throw(DimensionMismatch(
        "power semiconductor requires a two-by-two terminal Jacobian",
    ))
    voltage_v = _power_semiconductor_terminal_voltage(terminal_voltage_v)
    channel_current_a, channel_jacobian_s =
        _power_semiconductor_channel_current_jacobian(device, voltage_v)
    fidelity = something(device.extended_fidelity)
    dynamic_current_a, dynamic_jacobian_s =
        _power_semiconductor_dynamic_current_jacobian(
            fidelity, voltage_v, time_s, device.event_localization_enabled)
    snubber_current_a, snubber_jacobian_s =
        _power_semiconductor_snubber_current_jacobian(device, voltage_v)
    current_a = channel_current_a + dynamic_current_a + snubber_current_a
    jacobian_s = channel_jacobian_s + dynamic_jacobian_s + snubber_jacobian_s
    terminal_current_a[1] = current_a
    terminal_current_a[2] = -current_a
    terminal_jacobian_s[1, 1] = jacobian_s
    terminal_jacobian_s[1, 2] = -jacobian_s
    terminal_jacobian_s[2, 1] = -jacobian_s
    terminal_jacobian_s[2, 2] = jacobian_s
    return nothing
end

function _power_semiconductor_compensated_channel_current_jacobian(
    device::PowerSemiconductorSwitch,
    voltage_v::BigFloat,
)
    if device.closed
        direction = device isa PowerSemiconductorSwitch{TriacJunction} ?
            BigFloat(device.conduction_direction) : BigFloat(1)
        conductance = BigFloat(device.on_conductance)
        return conductance * (
            voltage_v - direction * BigFloat(device.forward_voltage_drop_v)
        ), conductance
    elseif device.reverse_diode_conducting
        diode = something(device.antiparallel_diode)
        conductance = BigFloat(diode.on_conductance_s)
        return conductance * (
            voltage_v + BigFloat(diode.forward_voltage_v)
        ), conductance
    end
    fidelity = device.extended_fidelity
    recovery = fidelity === nothing ? nothing : fidelity.recovered_charge
    conductance = if recovery !== nothing &&
                     recovery.stored_charge_c > 0.0 && voltage_v < 0
        BigFloat(device.on_conductance)
    else
        BigFloat(device.off_conductance)
    end
    return conductance * voltage_v, conductance
end

function _power_semiconductor_compensated_junction_charge_capacitance(
    junction::NonlinearJunctionChargeFidelity,
    voltage_v::BigFloat,
)
    BigFloat(junction.minimum_voltage_v) <= voltage_v <=
        BigFloat(junction.maximum_voltage_v) || throw(DomainError(
            Float64(voltage_v),
            "junction voltage lies outside its declared domain",
        ))
    # Q(V) and dQ/dV are constitutive quantities, not accepted-history
    # lookups. A Float64 history is not an exact106-bit evaluation: returning
    # it at one matching voltage creates a jump and a different Jacobian at
    # that point. The dynamic companion consumes causal charge separately.
    capacitance = BigFloat(junction.zero_bias_capacitance_f)
    voltage_v >= 0 && return capacitance * voltage_v, capacitance
    junction_voltage = BigFloat(junction.junction_voltage_v)
    exponent = BigFloat(junction.grading_exponent)
    normalized = 1 - voltage_v / junction_voltage
    normalized_negative_power = normalized^(-exponent)
    charge = -capacitance * junction_voltage * (
        normalized * normalized_negative_power - 1
    ) / (1 - exponent)
    return charge, capacitance * normalized_negative_power
end

function _power_semiconductor_compensated_dynamic_current_jacobian(
    fidelity::PowerSemiconductorExtendedFidelity,
    voltage_v::BigFloat,
    time_s::Union{Float64,ContinuousEvaluationTime},
    localized_tail::Bool=false,
)
    fidelity.candidate_prepared || throw(ArgumentError(
        "extended semiconductor must be prepared for an exact candidate step before evaluation",
    ))
    bookkeeping_time_s = time_s isa ContinuousEvaluationTime ? time_s.high_s : time_s
    tolerance = 16 * eps(max(1.0, abs(bookkeeping_time_s), abs(fidelity.candidate_time_s)))
    abs(bookkeeping_time_s - fidelity.candidate_time_s) <= tolerance || throw(ArgumentError(
        "extended semiconductor candidate time does not match the nonlinear evaluation time",
    ))
    current_a = BigFloat(0)
    jacobian_s = BigFloat(0)
    junction = fidelity.junction_charge
    if junction !== nothing
        charge, capacitance =
            _power_semiconductor_compensated_junction_charge_capacitance(
                junction,
                voltage_v,
            )
        step_s = BigFloat(fidelity.candidate_step_s)
        if fidelity.candidate_method === :trapezoidal
            current_a += 2 * (charge - BigFloat(junction.previous_charge_c)) /
                step_s - BigFloat(junction.last_displacement_current_a)
            jacobian_s += 2 * capacitance / step_s
        else
            current_a += (charge - BigFloat(junction.previous_charge_c)) / step_s
            jacobian_s += capacitance / step_s
        end
    end
    tail = fidelity.turn_off_tail
    if tail !== nothing && tail.active
        tail_current = if time_s isa ContinuousEvaluationTime
            duration_s = max(BigFloat(0), _semiconductor_tail_elapsed(time_s, tail.turn_off_time_s))
            BigFloat(tail.initial_current_a) * exp(-duration_s / BigFloat(tail.decay_time_s))
        else
            duration_s = max(0.0, _semiconductor_tail_elapsed(time_s, tail.turn_off_time_s))
            tail.initial_current_a * exp(-duration_s / tail.decay_time_s)
        end
        (localized_tail || tail_current > tail.cutoff_current_a) &&
            (current_a += BigFloat(tail_current))
    end
    return current_a, jacobian_s
end

function _power_semiconductor_compensated_snubber_current_jacobian(
    device::PowerSemiconductorSwitch,
    voltage_v::BigFloat,
)
    snubber = device.snubber
    snubber === nothing && return BigFloat(0), BigFloat(0)
    fidelity = something(device.extended_fidelity)
    step_s = BigFloat(fidelity.candidate_step_s)
    capacitance = BigFloat(snubber.capacitance_f)
    resistance = BigFloat(snubber.resistance_ohm)
    if fidelity.candidate_method === :trapezoidal
        capacitive_resistance = step_s / (2 * capacitance)
        conductance = inv(resistance + capacitive_resistance)
        current = conductance * (
            voltage_v -
            capacitive_resistance * BigFloat(snubber.previous_current_a) -
            BigFloat(snubber.capacitor_voltage_v)
        )
        return current, conductance
    end
    capacitive_resistance = step_s / capacitance
    conductance = inv(resistance + capacitive_resistance)
    return conductance * (
        voltage_v - BigFloat(snubber.capacitor_voltage_v)
    ), conductance
end

function nonlinear_current_jacobian_compensated!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    device::PowerSemiconductorSwitch,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_voltage_correction_v::AbstractVector{Float64},
    time_s::Union{Float64,ContinuousEvaluationTime},
    ; terminal_current_correction_a=nothing,
)
    terminal_current_correction_a === nothing ||
        _validate_nonlinear_current_correction_buffers(terminal_current_a,
            terminal_current_correction_a, terminal_voltage_v, terminal_voltage_correction_v)
    length(terminal_current_a) == 2 || throw(DimensionMismatch(
        "power semiconductor requires two terminal currents",
    ))
    size(terminal_jacobian_s) == (2, 2) || throw(DimensionMismatch(
        "power semiconductor requires a two-by-two terminal Jacobian",
    ))
    length(terminal_voltage_v) == 2 &&
        length(terminal_voltage_correction_v) == 2 || throw(DimensionMismatch(
            "power semiconductor requires two compensated terminal voltages",
        ))
    all(isfinite, terminal_voltage_v) &&
        all(isfinite, terminal_voltage_correction_v) || throw(ArgumentError(
            "power-semiconductor compensated terminal voltages must be finite",
        ))
    current_a, jacobian_s = setprecision(BigFloat, 106) do
        voltage_v =
            BigFloat(terminal_voltage_v[1]) +
            BigFloat(terminal_voltage_correction_v[1]) -
            BigFloat(terminal_voltage_v[2]) -
            BigFloat(terminal_voltage_correction_v[2])
        channel_current_a, channel_jacobian_s =
            _power_semiconductor_compensated_channel_current_jacobian(
                device,
                voltage_v,
            )
        fidelity = something(device.extended_fidelity)
        dynamic_current_a, dynamic_jacobian_s =
            _power_semiconductor_compensated_dynamic_current_jacobian(
                fidelity,
                voltage_v,
                time_s,
                device.event_localization_enabled,
            )
        snubber_current_a, snubber_jacobian_s =
            _power_semiconductor_compensated_snubber_current_jacobian(
                device,
                voltage_v,
            )
        return (
            channel_current_a + dynamic_current_a + snubber_current_a,
            channel_jacobian_s + dynamic_jacobian_s + snubber_jacobian_s,
        )
    end
    all(isfinite, (current_a, jacobian_s)) || throw(ArgumentError(
        "power-semiconductor compensated current and Jacobian must be finite",
    ))
    high_current, high_jacobian = Float64(current_a), Float64(jacobian_s)
    all(isfinite, (high_current, high_jacobian)) || throw(ArgumentError(
        "power-semiconductor compensated current and Jacobian exceed their Float64 output range",
    ))
    terminal_current_a[1] = high_current
    terminal_current_a[2] = -high_current
    terminal_jacobian_s[1, 1] = high_jacobian
    terminal_jacobian_s[1, 2] = -high_jacobian
    terminal_jacobian_s[2, 1] = -high_jacobian
    terminal_jacobian_s[2, 2] = high_jacobian
    if terminal_current_correction_a !== nothing
        correction = setprecision(BigFloat, 106) do
            Float64(current_a - BigFloat(high_current))
        end
        terminal_current_correction_a[1] = correction
        terminal_current_correction_a[2] = -correction
    end
    return nothing
end

function nonlinear_current_jacobian_compensated!(
    current::AbstractVector{Float64}, jacobian::AbstractMatrix{Float64},
    current_correction::AbstractVector{Float64}, device::PowerSemiconductorSwitch,
    voltage::AbstractVector{Float64}, voltage_correction::AbstractVector{Float64},
    time::Union{Float64,ContinuousEvaluationTime},
)
    nonlinear_current_jacobian_compensated!(current, jacobian, device,
        voltage, voltage_correction, time; terminal_current_correction_a=current_correction)
end

"""Accept compensated semiconductor histories through their constitutive owner.

A zero voltage low part does not discard a nonzero physical current low part.
Charge/capacitance and companion current are reconstructed by the shared device
owner, never inferred by cancellation in the rounded terminal balance.
"""
function accept_nonlinear_device_compensated_state!(
    device::PowerSemiconductorSwitch, voltage::AbstractVector{Float64},
    voltage_correction::AbstractVector{Float64}, current::AbstractVector{Float64},
    current_correction::AbstractVector{Float64}, jacobian::AbstractMatrix{Float64},
    time::Union{Float64,ContinuousEvaluationTime},
)
    _validate_nonlinear_current_correction_buffers(current, current_correction,
        voltage, voltage_correction)
    all(isfinite, current_correction) &&
        all(i -> current[i] + current_correction[i] == current[i], eachindex(current)) ||
        throw(ArgumentError("semiconductor current correction must be finite and nonoverlapping"))
    if device.extended_fidelity !== nothing
        bookkeeping = time isa ContinuousEvaluationTime ? time.high_s : time
        return _accept_power_semiconductor_state!(device, voltage, current,
            bookkeeping, jacobian, voltage_correction;
            evaluation_time=time, terminal_current_correction_a=current_correction)
    end
    accept_nonlinear_device_compensated_state!(device, voltage, voltage_correction,
        current, jacobian, time)
end

function nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    device::PowerSemiconductorSwitch,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::ContinuousEvaluationTime,
)
    nonlinear_current_jacobian_compensated!(terminal_current_a, terminal_jacobian_s,
        device, terminal_voltage_v, zeros(length(terminal_voltage_v)), time_s)
end

function accept_nonlinear_device_compensated_state!(
    device::PowerSemiconductorSwitch,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_voltage_correction_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    time_s::Float64,
)
    length(terminal_voltage_correction_v) == 2 || throw(DimensionMismatch(
        "power semiconductor requires two accepted voltage corrections",
    ))
    @inbounds for index in eachindex(
        terminal_voltage_v,
        terminal_voltage_correction_v,
    )
        correction = terminal_voltage_correction_v[index]
        high = terminal_voltage_v[index]
        high + correction == high || throw(ArgumentError(
            "power-semiconductor accepted voltage correction must remain a nonoverlapping low part",
        ))
    end
    return _accept_power_semiconductor_state!(
        device,
        terminal_voltage_v,
        terminal_current_a,
        time_s,
        terminal_jacobian_s,
        terminal_voltage_correction_v,
    )
end

function _power_semiconductor_switching_temperature_k(
    fidelity::PowerSemiconductorExtendedFidelity,
)
    thermal = fidelity.thermal
    thermal === nothing && return 293.15
    return thermal.node_temperature_k[1]
end

"""Validate an accepted electrical boundary before changing counter coordinates.

The common offset preserves every equality and pending-transition predicate.
An unpaid transition must be accepted by its owner, never erased by a rebase.
"""
function _power_semiconductor_periodic_counter_offset(device::PowerSemiconductorSwitch)
    fidelity = something(device.extended_fidelity)
    offset = device.topology_transition_count
    offset >= 0 && fidelity.accepted_topology_transition_count == offset ||
        throw(ArgumentError("periodic semiconductor boundary has an unaccepted transition"))
    for frame in (fidelity.accepted_electrical_history,
        fidelity.boundary_electrical_anchor, fidelity.boundary_electrical_state)
        frame !== nothing && frame.topology_transition_count == offset ||
            throw(ArgumentError("periodic semiconductor boundary has an unmatched electrical witness"))
    end
    energy = fidelity.switching_energy
    energy === nothing || energy.last_event_transition_count == offset ||
        throw(ArgumentError("periodic semiconductor boundary has an unpaid switching transition"))
    return offset
end

function _rebase_power_semiconductor_boundary_counts!(device::PowerSemiconductorSwitch)
    offset = _power_semiconductor_periodic_counter_offset(device)
    iszero(offset) && return device
    fidelity = something(device.extended_fidelity)
    rebase(frame) = AcceptedSemiconductorElectricalHistory(
        frame.voltage_v, frame.terminal_current_a, frame.channel_current_a,
        frame.displacement_current_a, frame.tail_current_a, frame.snubber_current_a,
        frame.topology_transition_count - offset, frame.arithmetic_policy)
    # Immutable past witnesses remain available to their original owners.
    fidelity.accepted_electrical_history = rebase(fidelity.accepted_electrical_history)
    fidelity.boundary_electrical_anchor = rebase(fidelity.boundary_electrical_anchor)
    fidelity.boundary_electrical_state = rebase(fidelity.boundary_electrical_state)
    device.topology_transition_count -= offset
    fidelity.accepted_topology_transition_count -= offset
    fidelity.switching_energy === nothing ||
        (fidelity.switching_energy.last_event_transition_count -= offset)
    return device
end

function _deposit_power_semiconductor_switching_energy!(
    device::PowerSemiconductorSwitch,
    voltage_v::Float64,
    current_a::Float64,
    time_s::Float64,
)
    fidelity = something(device.extended_fidelity)
    energy_table = fidelity.switching_energy
    energy_table === nothing && return 0.0
    event_kind = :none
    recovery = fidelity.recovered_charge
    if recovery !== nothing && _semiconductor_history_event_time(recovery.recovery_start_time_s) == time_s &&
       energy_table.last_reverse_recovery_start_time_s != time_s
        event_kind = :reverse_recovery
    elseif device.topology_transition_count > energy_table.last_event_transition_count
        event_kind = device.closed ? :turn_on : :turn_off
    end
    event_kind === :none && return 0.0
    interpolation_current_a = event_kind in (:turn_on, :turn_off) ?
        fidelity.pending_event_current_a : abs(current_a)
    interpolation_voltage_v = event_kind in (:turn_on, :turn_off) ?
        fidelity.pending_event_blocking_voltage_v : abs(voltage_v)
    event_energy_j = power_semiconductor_switching_energy(
        energy_table,
        event_kind,
        interpolation_current_a,
        interpolation_voltage_v,
        _power_semiconductor_switching_temperature_k(fidelity),
    )
    if event_kind === :turn_on
        energy_table.cumulative_turn_on_energy_j += event_energy_j
        energy_table.last_event_transition_count = device.topology_transition_count
    elseif event_kind === :turn_off
        energy_table.cumulative_turn_off_energy_j += event_energy_j
        energy_table.last_event_transition_count = device.topology_transition_count
    else
        energy_table.cumulative_reverse_recovery_energy_j += event_energy_j
        energy_table.last_reverse_recovery_start_time_s = time_s
    end
    energy_table.last_event_kind = event_kind
    energy_table.last_event_energy_j = event_energy_j
    fidelity.thermal === nothing ||
        _deposit_cauer_event_energy!(fidelity.thermal, event_energy_j)
    return event_energy_j
end

function _accept_power_semiconductor_state!(
    device::PowerSemiconductorSwitch,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    time_s::Float64,
    terminal_jacobian_s::Union{Nothing,AbstractMatrix{Float64}},
    terminal_voltage_correction_v::Union{Nothing,AbstractVector{Float64}}=nothing,
    ; evaluation_time::Union{Float64,ContinuousEvaluationTime}=time_s,
      terminal_current_correction_a::Union{Nothing,AbstractVector{Float64}}=nothing,
)
    fidelity = something(device.extended_fidelity)
    fidelity.candidate_prepared || throw(ArgumentError(
        "extended semiconductor acceptance requires one prepared candidate step",
    ))
    voltage_v = _power_semiconductor_terminal_voltage(terminal_voltage_v)
    # A supplied compensated contract remains compensated at zero low voltage:
    # its current can still contain physical contributions below the high part.
    compensated_voltage = terminal_voltage_correction_v !== nothing
    exact_voltage_v = if compensated_voltage
        length(terminal_voltage_correction_v) == 2 || throw(DimensionMismatch(
            "power semiconductor requires two accepted voltage corrections",
        ))
        setprecision(BigFloat, 106) do
            BigFloat(terminal_voltage_v[1]) +
                BigFloat(terminal_voltage_correction_v[1]) -
                BigFloat(terminal_voltage_v[2]) -
                BigFloat(terminal_voltage_correction_v[2])
        end
    else
        BigFloat(voltage_v)
    end
    compensated_voltage && (voltage_v = Float64(exact_voltage_v))
    length(terminal_current_a) == 2 || throw(DimensionMismatch(
        "power semiconductor requires two accepted terminal currents",
    ))
    current_a = terminal_current_a[1]
    isfinite(current_a) && isfinite(time_s) || throw(ArgumentError(
        "accepted power-semiconductor current and time must be finite",
    ))
    terminal_jacobian_s === nothing || size(terminal_jacobian_s) == (2, 2) ||
        throw(DimensionMismatch(
            "power semiconductor requires a two-by-two accepted terminal Jacobian",
        ))
    step_s = fidelity.candidate_step_s
    previous_terminal_voltage_v = fidelity.previous_terminal_voltage_v
    previous_terminal_current_a = fidelity.previous_terminal_current_a
    incoming_snubber_capacitor_v = device.snubber === nothing ? 0.0 : device.snubber.capacitor_voltage_v
    incoming_snubber_history_current_a = device.snubber === nothing ? 0.0 : device.snubber.previous_current_a
    previous_snubber_current_a = device.snubber === nothing ? 0.0 :
        device.snubber.last_current_a
    previous_snubber_stored_energy_j = device.snubber === nothing ? 0.0 :
        0.5 * device.snubber.capacitance_f * device.snubber.capacitor_voltage_v^2
    previous_snubber_dissipated_energy_j = device.snubber === nothing ? 0.0 :
        device.snubber.dissipated_energy_j
    previous_semiconductor_loss_w = device.previous_semiconductor_loss_w
    channel_current_a, channel_jacobian_s = if compensated_voltage
        current, jacobian = setprecision(BigFloat, 106) do
            _power_semiconductor_compensated_channel_current_jacobian(
                device,
                exact_voltage_v,
            )
        end
        Float64(current), Float64(jacobian)
    else
        _power_semiconductor_channel_current_jacobian(device, voltage_v)
    end
    tail = fidelity.turn_off_tail
    evaluated_tail_current_a = if tail !== nothing && tail.active
        if evaluation_time isa ContinuousEvaluationTime
            setprecision(BigFloat,106) do
                duration_s = max(BigFloat(0),_semiconductor_tail_elapsed(evaluation_time,tail.turn_off_time_s))
                Float64(BigFloat(tail.initial_current_a)*exp(-duration_s/BigFloat(tail.decay_time_s)))
            end
        else
            duration_s = max(0.0,_semiconductor_tail_elapsed(time_s,tail.turn_off_time_s))
            tail.initial_current_a*exp(-duration_s/tail.decay_time_s)
        end
    else
        0.0
    end
    accepted_tail_current_a = tail !== nothing && tail.active &&
        (device.event_localization_enabled ||
         evaluated_tail_current_a > tail.cutoff_current_a) ? evaluated_tail_current_a : 0.0
    snubber_current_a, snubber_jacobian_s = if compensated_voltage
        current, jacobian = setprecision(BigFloat, 106) do
            _power_semiconductor_compensated_snubber_current_jacobian(
                device,
                exact_voltage_v,
            )
        end
        Float64(current), Float64(jacobian)
    else
        _power_semiconductor_snubber_current_jacobian(device, voltage_v)
    end
    junction = fidelity.junction_charge
    if junction !== nothing
        previous_charge_c = junction.previous_charge_c
        previous_voltage_v = junction.previous_voltage_v
        previous_energy_j = power_semiconductor_junction_stored_energy(
            junction,
            previous_voltage_v,
        )
        previous_displacement_current_a = junction.last_displacement_current_a
        displacement_current_a = current_a - channel_current_a -
            accepted_tail_current_a - snubber_current_a
        if compensated_voltage
            # The junction owns Q(V), C(V), and its companion current. Recovering
            # them from differences of rounded terminal/channel currents or
            # Jacobians loses small dynamics beside a large conducting channel.
            charge, capacitance, displacement_current_a = setprecision(BigFloat, 106) do
                exact_charge, exact_capacitance =
                    _power_semiconductor_compensated_junction_charge_capacitance(
                        junction, exact_voltage_v,
                    )
                factor = fidelity.candidate_method === :trapezoidal ? 2 : 1
                displacement = factor * (exact_charge - BigFloat(previous_charge_c)) /
                    BigFloat(step_s)
                fidelity.candidate_method === :trapezoidal &&
                    (displacement -= BigFloat(previous_displacement_current_a))
                Float64(exact_charge), Float64(exact_capacitance), Float64(displacement)
            end
        elseif terminal_jacobian_s === nothing
            charge, capacitance = _power_semiconductor_junction_charge_capacitance(
                junction,
                voltage_v,
            )
        else
            dynamic_jacobian_s = (
                terminal_jacobian_s[1, 1] -
                channel_jacobian_s -
                snubber_jacobian_s
            )
            if fidelity.candidate_method === :trapezoidal
                capacitance = 0.5 * dynamic_jacobian_s * step_s
                charge = previous_charge_c + 0.5 * step_s * (
                    junction.last_displacement_current_a + displacement_current_a
                )
            else
                capacitance = dynamic_jacobian_s * step_s
                charge = previous_charge_c + step_s * displacement_current_a
            end
        end
        junction.previous_voltage_v = voltage_v
        junction.previous_charge_c = charge
        junction.last_capacitance_f = capacitance
        junction.last_charge_c = charge
        junction.last_displacement_current_a = displacement_current_a
        current_energy_j = power_semiconductor_junction_stored_energy(
            junction,
            voltage_v,
        )
        junction_work_j = if fidelity.candidate_method === :trapezoidal
            0.5 * step_s * (
                previous_voltage_v * previous_displacement_current_a +
                voltage_v * displacement_current_a
            )
        else
            step_s * voltage_v * displacement_current_a
        end
        fidelity.companion_energy_residual_j +=
            junction_work_j - (current_energy_j - previous_energy_j)
    end
    recovery = fidelity.recovered_charge
    if recovery !== nothing
        previous_charge_c = recovery.stored_charge_c
        candidate_charge_c = _power_semiconductor_recovery_charge(
            recovery, channel_current_a, step_s,
        )
        fidelity.candidate_recovery_charge_c = candidate_charge_c
        accepted_charge_c = max(0.0, candidate_charge_c)
        charge_change_c = accepted_charge_c - previous_charge_c
        if channel_current_a >= 0.0
            recovery.last_recovery_current_a = 0.0
            recovery.recovery_active = false
        else
            recovered_charge_c = max(0.0, -charge_change_c)
            recovery.last_recovery_current_a = min(0.0, channel_current_a)
            if recovered_charge_c > 0.0 && !recovery.recovery_active
                recovery.recovery_start_time_s = SemiconductorHistoryTime(evaluation_time)
            end
            recovery.recovery_active = accepted_charge_c > 0.0 && recovered_charge_c > 0.0
            recovery.peak_reverse_current_a = min(
                recovery.peak_reverse_current_a,
                recovery.last_recovery_current_a,
            )
            recovery.cumulative_recovered_charge_c += recovered_charge_c
            recovery.last_recovery_duration_s = isfinite(recovery.recovery_start_time_s) ?
                max(0.0, _semiconductor_history_elapsed(evaluation_time, recovery.recovery_start_time_s)) : 0.0
        end
        recovery.previous_stored_charge_c = previous_charge_c
        recovery.stored_charge_c = max(0.0, accepted_charge_c)
    end
    if tail !== nothing && tail.active
        duration_s = evaluation_time isa ContinuousEvaluationTime ?
            Float64(max(BigFloat(0),_semiconductor_tail_elapsed(evaluation_time,tail.turn_off_time_s))) :
            max(0.0,_semiconductor_tail_elapsed(time_s,tail.turn_off_time_s))
        if !device.event_localization_enabled &&
           evaluated_tail_current_a <= tail.cutoff_current_a
            tail.active = false
            tail.current_a = 0.0
        else
            tail.current_a = evaluated_tail_current_a
        end
        tail.last_duration_s = duration_s
    end
    tail_current_a = tail === nothing ? 0.0 : tail.current_a
    instantaneous_loss_w = max(0.0, setprecision(BigFloat, 106) do
        Float64(exact_voltage_v * BigFloat(channel_current_a + tail_current_a))
    end)
    thermal = fidelity.thermal
    device.last_voltage = voltage_v
    device.last_current = current_a
    device.last_forward_current_a = device.closed ? channel_current_a : 0.0
    device.last_reverse_diode_current_a =
        device.reverse_diode_conducting ? -channel_current_a : 0.0
    if device.snubber !== nothing
        if fidelity.candidate_method === :trapezoidal
            _update_series_rc_snubber!(
                device.snubber,
                voltage_v,
                step_s;
                accepted_current_a=compensated_voltage ? snubber_current_a : nothing,
            )
        else
            _update_series_rc_snubber_backward_euler!(
                device.snubber,
                voltage_v,
                step_s;
                accepted_current_a=compensated_voltage ? snubber_current_a : nothing,
            )
        end
        device.last_snubber_current_a = device.snubber.last_current_a
    end
    if !device.event_localization_enabled
        previous_forward = device.closed
        previous_reverse = device.reverse_diode_conducting
        if device isa PowerSemiconductorSwitch{GateTurnOffThyristorJunction}
            _update_power_semiconductor_sampled_state!(device, evaluation_time)
        else
            _update_power_semiconductor_sampled_state!(device)
        end
        if previous_forward && !device.closed
            _activate_power_semiconductor_tail!(device, evaluation_time)
        end
        _update_power_semiconductor_reverse_state!(device)
        _record_power_semiconductor_topology_transition!(
            device,
            previous_forward,
            previous_reverse,
            time_s,
        )
    end
    event_energy_j = _deposit_power_semiconductor_switching_energy!(
        device,
        voltage_v,
        channel_current_a,
        time_s,
    )
    thermal === nothing || _accept_cauer_thermal_step!(
        thermal,
        instantaneous_loss_w,
        step_s,
    )
    conduction_energy_increment_j = if fidelity.candidate_method === :trapezoidal
        0.5 * step_s * (previous_semiconductor_loss_w + instantaneous_loss_w)
    else
        step_s * instantaneous_loss_w
    end
    device.semiconductor_dissipated_energy_j += conduction_energy_increment_j + event_energy_j
    previous_junction_current_a = junction === nothing ? 0.0 :
        previous_displacement_current_a
    channel_work_j = if fidelity.candidate_method === :trapezoidal
        0.5 * step_s * (
            previous_terminal_voltage_v * (
                previous_terminal_current_a - previous_junction_current_a -
                previous_snubber_current_a
            ) +
            voltage_v * (channel_current_a + tail_current_a)
        )
    else
        step_s * voltage_v * (channel_current_a + tail_current_a)
    end
    semiconductor_energy_increment_j = conduction_energy_increment_j + event_energy_j
    fidelity.companion_energy_residual_j +=
        channel_work_j - semiconductor_energy_increment_j
    if device.snubber !== nothing
        snubber = device.snubber
        snubber_work_j = if fidelity.candidate_method === :trapezoidal
            0.5 * step_s * (previous_terminal_voltage_v * previous_snubber_current_a +
                voltage_v * snubber.last_current_a)
        else
            step_s * voltage_v * snubber.last_current_a
        end
        snubber_stored_energy_j =
            0.5 * snubber.capacitance_f * snubber.capacitor_voltage_v^2
        fidelity.companion_energy_residual_j += snubber_work_j -
            (snubber.dissipated_energy_j - previous_snubber_dissipated_energy_j) -
            (snubber_stored_energy_j - previous_snubber_stored_energy_j)
    end
    device.last_semiconductor_loss_w = instantaneous_loss_w
    device.previous_semiconductor_loss_w = instantaneous_loss_w
    fidelity.previous_terminal_voltage_v = voltage_v
    fidelity.previous_terminal_current_a = current_a
    fidelity.accepted_topology_transition_count = device.topology_transition_count
    physical_current = setprecision(BigFloat,106) do
        BigFloat(current_a) + (terminal_current_correction_a === nothing ? BigFloat(0) : BigFloat(terminal_current_correction_a[1]))
    end
    fidelity.accepted_electrical_history = AcceptedSemiconductorElectricalHistory(
        setprecision(BigFloat,106) do; BigFloat(exact_voltage_v); end,
        physical_current, channel_current_a,
        junction === nothing ? 0.0 : junction.last_displacement_current_a,
        tail_current_a, device.last_snubber_current_a, device.topology_transition_count,
        compensated_voltage ? :compensated_binary : :binary64)
    fidelity.boundary_electrical_anchor=fidelity.accepted_electrical_history
    fidelity.boundary_electrical_state=fidelity.accepted_electrical_history
    if device.snubber !== nothing
        snubber = device.snubber
        snubber.boundary_voltage_v=setprecision(BigFloat,106) do;BigFloat(exact_voltage_v);end
        snubber.accepted_companion_history = AcceptedRCCompanionHistory(
            incoming_snubber_capacitor_v, incoming_snubber_history_current_a,
            setprecision(BigFloat,106) do; BigFloat(exact_voltage_v); end,
            step_s, fidelity.candidate_method, snubber.last_current_a, snubber.capacitor_voltage_v,
            compensated_voltage ? :compensated_binary : :binary64)
        snubber.boundary_anchor=RCBoundaryCurrentAnchor(snubber.boundary_voltage_v,snubber.capacitor_voltage_v,
            snubber.last_current_a,compensated_voltage ? :compensated_binary : :binary64)
    end
    return nothing
end


function accept_nonlinear_device_state!(
    device::PowerSemiconductorSwitch,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    time_s::Float64,
)
    return _accept_power_semiconductor_state!(
        device,
        terminal_voltage_v,
        terminal_current_a,
        time_s,
        nothing,
    )
end

function accept_nonlinear_device_state!(
    device::PowerSemiconductorSwitch,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    time_s::Float64,
)
    return _accept_power_semiconductor_state!(
        device,
        terminal_voltage_v,
        terminal_current_a,
        time_s,
        terminal_jacobian_s,
    )
end

function finish_nonlinear_device_step!(device::PowerSemiconductorSwitch)
    fidelity = something(device.extended_fidelity)
    fidelity.candidate_prepared || throw(ArgumentError(
        "extended semiconductor candidate step was already finished or never prepared",
    ))
    companion = _power_semiconductor_boundary_channel_companion(device)
    device.last_conductance = companion.conductance_s
    device.last_history_current_a = companion.history_current_a
    device.last_evaluation_time_s = fidelity.candidate_time_s
    fidelity.candidate_prepared = false
    return nothing
end

function accept_nonlinear_device_compensated_state!(device::PowerSemiconductorSwitch,
    voltage::AbstractVector{Float64},correction::AbstractVector{Float64},
    current::AbstractVector{Float64},jacobian::AbstractMatrix{Float64},time::ContinuousEvaluationTime)
    length(voltage)==length(correction)==2 || throw(DimensionMismatch("semiconductor compensated terminal count"))
    all(i->voltage[i]+correction[i]==voltage[i],eachindex(voltage)) ||
        throw(ArgumentError("semiconductor voltage correction must be nonoverlapping"))
    _accept_power_semiconductor_state!(device,voltage,current,time.high_s,jacobian,correction;
        evaluation_time=time)
end

function accept_nonlinear_device_state!(device::PowerSemiconductorSwitch,
    voltage::AbstractVector{Float64},current::AbstractVector{Float64},
    jacobian::AbstractMatrix{Float64},time::ContinuousEvaluationTime)
    _accept_power_semiconductor_state!(device,voltage,current,time.high_s,jacobian;
        evaluation_time=time)
end

# Localized gate ownership is fixed by the event owner, not by timed stamping.
Branches.continuous_time_invariant(device::PowerSemiconductorSwitch) = device.event_localization_enabled

function update!(
    device::PowerSemiconductorSwitch,
    voltages::AbstractVector{Float64},
    dt_s::Float64,
    ::Val{TrapezoidalCompanion},
)
    return update!(device, voltages, dt_s)
end

function update!(
    device::PowerSemiconductorSwitch,
    voltages::AbstractVector{Float64},
    dt_s::Float64,
    ::Val{BackwardEulerCompanion},
)
    terminal_voltage = Branches.branch_voltage(voltages, device.a, device.b)
    semiconductor_current = device.last_conductance * terminal_voltage +
        device.last_history_current_a
    forward_current = device.closed ? semiconductor_current : 0.0
    reverse_current = device.reverse_diode_conducting ? -semiconductor_current : 0.0
    snubber_current = 0.0
    if device.snubber !== nothing
        _update_series_rc_snubber_backward_euler!(device.snubber, terminal_voltage, dt_s)
        snubber_current = device.snubber.last_current_a
    end
    semiconductor_loss = max(0.0, terminal_voltage * semiconductor_current)
    device.semiconductor_dissipated_energy_j += dt_s * semiconductor_loss
    device.last_voltage = terminal_voltage
    device.last_forward_current_a = forward_current
    device.last_reverse_diode_current_a = reverse_current
    device.last_snubber_current_a = snubber_current
    device.last_current = semiconductor_current + snubber_current
    device.last_semiconductor_loss_w = semiconductor_loss
    device.previous_semiconductor_loss_w = semiconductor_loss
    if !device.event_localization_enabled
        previous_forward = device.closed
        previous_reverse = device.reverse_diode_conducting
        _update_power_semiconductor_sampled_state!(device)
        _update_power_semiconductor_reverse_state!(device)
        _record_power_semiconductor_topology_transition!(
            device,
            previous_forward,
            previous_reverse,
            device.last_evaluation_time_s,
        )
    end
    return nothing
end

function power_semiconductor_terminal_state(device::PowerSemiconductorSwitch)
    driver = device.gate_driver
    snubber = device.snubber
    return PowerSemiconductorTerminalState(
        power_semiconductor_kind(device),
        device.closed,
        device.reverse_diode_conducting,
        driver === nothing ? false : driver.commanded_on,
        driver === nothing ? false : driver.applied_on,
        device.last_voltage,
        device.last_current,
        device.last_forward_current_a,
        device.last_reverse_diode_current_a,
        device.last_snubber_current_a,
        device.last_semiconductor_loss_w,
        snubber === nothing ? 0.0 : snubber.last_resistor_loss_w,
        device.semiconductor_dissipated_energy_j,
        snubber === nothing ? 0.0 : snubber.dissipated_energy_j,
        snubber === nothing ? 0.0 :
            0.5 * snubber.capacitance_f * snubber.capacitor_voltage_v^2,
    )
end

power_semiconductor_has_extended_fidelity(device::PowerSemiconductorSwitch) =
    device.extended_fidelity !== nothing

function initialize_power_semiconductor_junction_state!(
    device::PowerSemiconductorSwitch,
    voltage_v::Real,
)
    fidelity = something(device.extended_fidelity)
    junction = something(fidelity.junction_charge)
    voltage = Float64(voltage_v)
    charge, capacitance = _power_semiconductor_junction_charge_capacitance(
        junction,
        voltage,
    )
    junction.previous_voltage_v = voltage
    junction.previous_charge_c = charge
    junction.last_capacitance_f = capacitance
    junction.last_charge_c = charge
    junction.last_displacement_current_a = 0.0
    fidelity.companion_energy_residual_j = 0.0
    fidelity.previous_terminal_voltage_v = voltage
    return device
end

function _activate_power_semiconductor_tail!(
    device::PowerSemiconductorSwitch,
    time_s::Union{Real,ContinuousEvaluationTime},
)
    fidelity = device.extended_fidelity
    fidelity === nothing && return device
    tail = fidelity.turn_off_tail
    tail === nothing && return device
    # Consecutive turn-off events do not erase carriers from an earlier tail.
    # Both contributions have this device's same decay constant, so their sum
    # is represented exactly by one exponential rebased at the new event time.
    # Activation/rebase has always been binary64, including localized events.
    elapsed = tail.active ? _semiconductor_tail_elapsed_binary64(time_s, tail.turn_off_time_s) : 0.0
    surviving_current = tail.active ? tail.initial_current_a * exp(
        -max(0.0, elapsed) / tail.decay_time_s,
    ) : 0.0
    initial_current = surviving_current + max(0.0, device.last_forward_current_a)
    if initial_current > tail.cutoff_current_a
        tail.active = true
        tail.current_a = initial_current
        tail.initial_current_a = initial_current
        tail.turn_off_time_s = SemiconductorHistoryTime(time_s)
        tail.last_duration_s = 0.0
    else
        tail.active = false
        tail.current_a = 0.0
        tail.initial_current_a = 0.0
        tail.turn_off_time_s = SemiconductorHistoryTime(time_s)
        tail.last_duration_s = 0.0
    end
    return device
end

"""Reconstruct analytically dependent semiconductor state at an EMT boundary.

The causal turn-off-tail history is its activation state, initial current, and
turn-off time. Instantaneous tail current and elapsed duration are derived from
that history and must never be optimized or restored independently. Active
recovery duration is likewise derived from its accepted start time.
"""
function reconstruct_power_semiconductor_boundary_state!(
    device::PowerSemiconductorSwitch,
    boundary_time_s::Union{Real,ContinuousEvaluationTime},
)
    boundary_time = boundary_time_s isa ContinuousEvaluationTime ?
        boundary_time_s.high_s : Float64(boundary_time_s)
    isfinite(boundary_time) || throw(ArgumentError(
        "a semiconductor boundary time must be finite",
    ))
    fidelity = device.extended_fidelity
    fidelity === nothing && return device
    recovery = fidelity.recovered_charge
    tail = fidelity.turn_off_tail
    function active_elapsed_age(start_time_s, history_name)
        isfinite(start_time_s) || throw(ArgumentError(
            "active $history_name requires a finite start time",
        ))
        exact_elapsed = _semiconductor_time_exact(boundary_time_s) -
            _semiconductor_time_exact(start_time_s)
        exact_elapsed >= 0 || throw(ArgumentError(
            "active $history_name start time exceeds its boundary time: " *
            "kind=$(power_semiconductor_kind(device)), nodes=($(device.a),$(device.b)), " *
            "boundary_high=$(boundary_time), boundary_low=$(boundary_time_s isa ContinuousEvaluationTime ? boundary_time_s.low_s : 0.0), " *
            "origin_type=$(nameof(typeof(start_time_s))), origin=$(first(repr(start_time_s),192)), " *
            "origin_high=$(start_time_s isa SemiconductorHistoryTime ? start_time_s.high_s : start_time_s), " *
            "origin_low=$(start_time_s isa SemiconductorHistoryTime ? first(string(start_time_s.low_s),192) : "0"), " *
            "event_time=$(_semiconductor_history_event_time(start_time_s)), " *
            "topology=$(device.topology_transition_count), accepted_topology=$(fidelity.accepted_topology_transition_count), " *
            "frame_topology=$(fidelity.boundary_electrical_state === nothing ? nothing : fidelity.boundary_electrical_state.topology_transition_count), " *
            "exact_elapsed=$(first(string(exact_elapsed),192))",
        ))
        return history_name == "turn-off tail" ?
            _semiconductor_tail_elapsed(boundary_time_s, start_time_s) :
            start_time_s isa SemiconductorHistoryTime ?
                _semiconductor_history_elapsed(boundary_time_s, start_time_s) :
                boundary_time_s isa ContinuousEvaluationTime ?
                    continuous_time_elapsed(boundary_time_s, start_time_s) :
                    boundary_time - start_time_s
    end
    # Validate both histories before rebuilding either dependent duration.
    recovery_age = recovery !== nothing && recovery.recovery_active ?
        active_elapsed_age(recovery.recovery_start_time_s, "recovered charge") : nothing
    tail_age = if tail !== nothing && tail.active
        if boundary_time_s isa ContinuousEvaluationTime
            setprecision(BigFloat, 106) do
                active_elapsed_age(tail.turn_off_time_s, "turn-off tail")
            end
        else
            active_elapsed_age(tail.turn_off_time_s, "turn-off tail")
        end
    else
        nothing
    end
    if recovery_age !== nothing
        recovery.last_recovery_duration_s = Float64(recovery_age)
    end
    if tail !== nothing
        if tail_age !== nothing
            if boundary_time_s isa ContinuousEvaluationTime
                setprecision(BigFloat, 106) do
                    tail.last_duration_s = Float64(tail_age)
                    tail.current_a = Float64(BigFloat(tail.initial_current_a) *
                        exp(-tail_age / BigFloat(tail.decay_time_s)))
                end
            else
                tail.last_duration_s = tail_age
                tail.current_a = tail.initial_current_a * exp(-tail_age / tail.decay_time_s)
            end
        else
            tail.current_a = 0.0
        end
    end
    return device
end


function _accepted_rc_companion_current(snubber::SeriesRCSnubber, frame::AcceptedRCCompanionHistory)
    if frame.arithmetic_policy===:binary64
        voltage=Float64(frame.physical_voltage_v)
        BigFloat(voltage)==frame.physical_voltage_v || throw(ArgumentError("binary64 RC frame contains an unconsumed voltage low"))
        alpha=frame.method===:trapezoidal ? frame.step_s/(2.0*snubber.capacitance_f) : frame.step_s/snubber.capacitance_f
        conductance=inv(snubber.resistance_ohm+alpha)
        history=frame.method===:trapezoidal ? conductance*(-alpha*frame.incoming_current_a-frame.incoming_capacitor_voltage_v) : -conductance*frame.incoming_capacitor_voltage_v
        return conductance*voltage+history
    end
    frame.arithmetic_policy===:compensated_binary || throw(ArgumentError("unknown accepted RC arithmetic policy"))
    setprecision(BigFloat,106) do
        step=BigFloat(frame.step_s);capacitance=BigFloat(snubber.capacitance_f)
        alpha=frame.method===:trapezoidal ? step/(2capacitance) : step/capacitance
        conductance=inv(BigFloat(snubber.resistance_ohm)+alpha)
        drive=frame.method===:trapezoidal ? frame.physical_voltage_v-alpha*BigFloat(frame.incoming_current_a)-BigFloat(frame.incoming_capacitor_voltage_v) :
            frame.physical_voltage_v-BigFloat(frame.incoming_capacitor_voltage_v)
        Float64(conductance*drive)
    end
end

# The fixed anchor is a represented numerical frame, not another physical DOF.
# Trial/backtrack composition always references this same immutable endpoint.
function _semiconductor_boundary_terminal_voltage(device,voltage,correction)
    fidelity=device.extended_fidelity
    anchor=fidelity===nothing ? nothing : fidelity.boundary_electrical_anchor===nothing ? fidelity.accepted_electrical_history : fidelity.boundary_electrical_anchor
    high(node)=node==0 ? 0.0 : voltage[node]
    low(node)=node==0 ? 0.0 : correction[node]
    if anchor!==nothing && anchor.arithmetic_policy===:binary64
        iszero(low(device.a))&&iszero(low(device.b)) || throw(ArgumentError("binary64 producer cannot consume supplied terminal voltage lows"))
        value=high(device.a)-high(device.b)
        return setprecision(BigFloat,106) do;BigFloat(value);end
    end
    setprecision(BigFloat,106) do
        BigFloat(high(device.a))+BigFloat(low(device.a))-BigFloat(high(device.b))-BigFloat(low(device.b))
    end
end

function _guard_rc_boundary_markov_current(snubber,voltage,capacitor,current)
    all(isfinite,(voltage,capacitor,current)) || throw(ArgumentError("nonfinite RC boundary frame"))
    exact=Rational{BigInt}
    v,q,i,r=exact(voltage),exact(capacitor),exact(current),exact(snubber.resistance_ohm)
    defect=r*i-(v-q)
    abs(defect/r)<=exact(1e-8) || throw(DomainError(defect,"RC boundary exceeds original current/KCL bound"))
    scale=max(abs(v*i),abs(r*i*i),abs(q*i),exact(eps(Float64)))
    abs(defect*i)<=exact(2e-6)*scale || throw(DomainError(defect,"RC boundary exceeds original power/energy bound"))
    return current
end

function _rc_boundary_markov_current(snubber::SeriesRCSnubber,anchor::Union{AcceptedRCCompanionHistory,RCBoundaryCurrentAnchor},voltage::Real,capacitor::Float64)
    if anchor.arithmetic_policy===:binary64
        high=Float64(voltage)
        BigFloat(high)==BigFloat(voltage) || throw(ArgumentError("binary64 boundary contains an unconsumed voltage low"))
        current=anchor.current_a+((high-Float64(anchor.physical_voltage_v))-(capacitor-anchor.outgoing_capacitor_voltage_v))/snubber.resistance_ohm
        return _guard_rc_boundary_markov_current(snubber,voltage,capacitor,current)
    end
    anchor.arithmetic_policy===:compensated_binary || throw(ArgumentError("unknown RC boundary arithmetic"))
    current=setprecision(BigFloat,106) do
        Float64(BigFloat(anchor.current_a)+((BigFloat(voltage)-anchor.physical_voltage_v)-
            (BigFloat(capacitor)-BigFloat(anchor.outgoing_capacitor_voltage_v)))/BigFloat(snubber.resistance_ohm))
    end
    return _guard_rc_boundary_markov_current(snubber,voltage,capacitor,current)
end

_rc_boundary_markov_tangent(snubber::SeriesRCSnubber,voltage_direction,capacitor_direction) =
    (voltage_direction .- capacitor_direction) ./ snubber.resistance_ohm

function _accepted_rc_companion_capacitor(snubber::SeriesRCSnubber, frame::AcceptedRCCompanionHistory, current::Float64)
    frame.method===:trapezoidal ? frame.incoming_capacitor_voltage_v+
        frame.step_s/(2.0*snubber.capacitance_f)*(current+frame.incoming_current_a) :
        frame.incoming_capacitor_voltage_v+frame.step_s*current/snubber.capacitance_f
end

"""Reconstruct every instantaneous semiconductor boundary dependent from its accepted owners.

The discrete conduction state, recovered charge, junction displacement-current
history, snubber capacitor voltage, tail activation history, thermal nodes, and
event observations remain canonical owners. Terminal, companion, junction,
snubber, tail, and thermal endpoint fields are rebuilt deterministically from
those owners and the accepted terminal voltage.
"""
function reconstruct_power_semiconductor_boundary_state!(
    device::PowerSemiconductorSwitch,
    boundary_time_s::Union{Real,ContinuousEvaluationTime},
    terminal_voltage_v::Real,
)
    # Validate all temporal, endpoint and causal companion guards on an isolated
    # copy first. No destination/alias changes occur if reconstruction rejects.
    _reconstruct_power_semiconductor_boundary_state!(deepcopy(device),boundary_time_s,terminal_voltage_v)
    return _reconstruct_power_semiconductor_boundary_state!(device,boundary_time_s,terminal_voltage_v)
end

function _reconstruct_power_semiconductor_boundary_state!(
    device::PowerSemiconductorSwitch,
    boundary_time_s::Union{Real,ContinuousEvaluationTime},
    terminal_voltage_v::Real,
)
    boundary_time = boundary_time_s isa ContinuousEvaluationTime ?
        boundary_time_s.high_s : Float64(boundary_time_s)
    voltage_v = Float64(terminal_voltage_v)
    all(isfinite, (boundary_time, voltage_v)) || throw(ArgumentError(
        "semiconductor boundary time and terminal voltage must be finite",
    ))
    reconstruct_power_semiconductor_boundary_state!(device, boundary_time_s)
    companion = _power_semiconductor_boundary_channel_companion(device)
    device.last_evaluation_time_s = boundary_time
    device.last_conductance = companion.conductance_s
    device.last_history_current_a = companion.history_current_a
    fidelity = device.extended_fidelity
    endpoint = fidelity===nothing ? nothing : fidelity.boundary_electrical_anchor===nothing ? fidelity.accepted_electrical_history : fidelity.boundary_electrical_anchor
    physical_voltage = setprecision(BigFloat,106) do;BigFloat(terminal_voltage_v);end
    if endpoint !== nothing
        endpoint.topology_transition_count==device.topology_transition_count ||
            throw(ArgumentError("semiconductor boundary changed its fixed topology anchor"))
    end
    channel_current_a = endpoint === nothing || endpoint.arithmetic_policy===:binary64 ? first(_power_semiconductor_accepted_boundary_channel_current_jacobian(device,voltage_v)) :
        setprecision(BigFloat,106) do
            Float64(first(_power_semiconductor_compensated_channel_current_jacobian(device,physical_voltage)))
        end
    displacement_current_a = 0.0
    tail_current_a = 0.0
    if fidelity !== nothing
        junction = fidelity.junction_charge
        if junction !== nothing
            displacement_current_a = junction.last_displacement_current_a
            charge_c, capacitance_f = if endpoint===nothing || endpoint.arithmetic_policy===:binary64
                _power_semiconductor_junction_charge_capacitance(junction,voltage_v)
            else
                setprecision(BigFloat,106) do
                    charge,capacitance=_power_semiconductor_compensated_junction_charge_capacitance(junction,physical_voltage)
                    Float64(charge),Float64(capacitance)
                end
            end
            junction.previous_voltage_v = voltage_v
            junction.previous_charge_c = charge_c
            junction.last_capacitance_f = capacitance_f
            junction.last_charge_c = charge_c
            junction.last_displacement_current_a = displacement_current_a
        end
        tail = fidelity.turn_off_tail
        tail_current_a = tail === nothing || !tail.active ? 0.0 : tail.current_a
    end
    snubber = device.snubber
    companion_frame=snubber === nothing ? nothing : snubber.accepted_companion_history
    snubber_current_a = if snubber===nothing
        0.0
    elseif companion_frame===nothing
        if snubber.boundary_anchor!==nothing
            _rc_boundary_markov_current(snubber,snubber.boundary_anchor,physical_voltage,snubber.capacitor_voltage_v)
        else
            endpoint===nothing || throw(ArgumentError("RC current boundary anchor is missing"))
            (voltage_v-snubber.capacitor_voltage_v)/snubber.resistance_ohm
        end
    else
        witnessed=_accepted_rc_companion_current(snubber,companion_frame)
        witnessed==companion_frame.current_a && _accepted_rc_companion_capacitor(snubber,companion_frame,witnessed)==companion_frame.outgoing_capacitor_voltage_v ||
            throw(ArgumentError("immutable RC past witness fails its original rounded law"))
        anchor=snubber.boundary_anchor===nothing ? companion_frame : snubber.boundary_anchor
        _rc_boundary_markov_current(snubber,anchor,physical_voltage,snubber.capacitor_voltage_v)
    end
    physical_current=if endpoint===nothing
        setprecision(BigFloat,106) do;BigFloat(channel_current_a+displacement_current_a+tail_current_a+snubber_current_a);end
    elseif endpoint.arithmetic_policy===:binary64
        setprecision(BigFloat,106) do
            BigFloat(Float64(endpoint.terminal_current_a)+(channel_current_a-endpoint.channel_current_a)+
                (displacement_current_a-endpoint.displacement_current_a)+(tail_current_a-endpoint.tail_current_a)+
                (snubber_current_a-endpoint.snubber_current_a))
        end
    else
        setprecision(BigFloat,106) do
            endpoint.terminal_current_a+(BigFloat(channel_current_a)-BigFloat(endpoint.channel_current_a))+
                (BigFloat(displacement_current_a)-BigFloat(endpoint.displacement_current_a))+
                (BigFloat(tail_current_a)-BigFloat(endpoint.tail_current_a))+
                (BigFloat(snubber_current_a)-BigFloat(endpoint.snubber_current_a))
        end
    end
    terminal_current_a=Float64(physical_current)
    if snubber !== nothing
        snubber.boundary_voltage_v=physical_voltage
        snubber.previous_current_a = snubber_current_a
        snubber.last_branch_voltage_v = voltage_v
        snubber.last_current_a = snubber_current_a
        snubber.last_resistor_loss_w =
            snubber.resistance_ohm * snubber_current_a^2
    end
    device.last_voltage = voltage_v
    device.last_current = terminal_current_a
    device.last_forward_current_a = device.closed ? channel_current_a : 0.0
    device.last_reverse_diode_current_a = device.reverse_diode_conducting ?
        -channel_current_a : 0.0
    device.last_snubber_current_a = snubber_current_a
    device.last_semiconductor_loss_w = setprecision(BigFloat,106) do
        max(0.0,Float64(physical_voltage*BigFloat(channel_current_a+tail_current_a)))
    end
    device.previous_semiconductor_loss_w = device.last_semiconductor_loss_w
    if fidelity !== nothing
        fidelity.previous_terminal_voltage_v = voltage_v
        fidelity.previous_terminal_current_a = terminal_current_a
        if endpoint!==nothing
            fidelity.boundary_electrical_anchor=endpoint
            fidelity.boundary_electrical_state=AcceptedSemiconductorElectricalHistory(physical_voltage,physical_current,
                channel_current_a,displacement_current_a,tail_current_a,snubber_current_a,device.topology_transition_count,endpoint.arithmetic_policy)
        end
        thermal = fidelity.thermal
        if thermal !== nothing
            copyto!(thermal.trial_temperature_k, thermal.node_temperature_k)
            thermal.last_stored_energy_j =
                power_semiconductor_thermal_stored_energy(thermal)
            thermal.last_ambient_heat_flow_w =
                thermal.accepted_temperature_rise_k[end] /
                thermal.resistance_k_per_w[end]
        end
    end
    return device
end

struct SemiconductorEventBoundaryContext{D,T}
    device::D
    time_s::T
    chart::PowerSemiconductorTopologyChartState
    frame::AcceptedSemiconductorElectricalHistory
    accepted_topology_transition_count::Int
end

function NonlinearNetwork.nonlinear_device_event_boundary_context(
    device::PowerSemiconductorSwitch, time_s,
)
    fidelity = device.extended_fidelity
    fidelity === nothing && return nothing
    frame = fidelity.boundary_electrical_state === nothing ?
        fidelity.accepted_electrical_history : fidelity.boundary_electrical_state
    frame === nothing && return nothing
    topology_matches = frame.topology_transition_count == device.topology_transition_count
    voltage_matches = Float64(frame.voltage_v) == device.last_voltage
    current_matches = Float64(frame.terminal_current_a) == device.last_current
    topology_matches && voltage_matches && current_matches || throw(ArgumentError(
        "device event requires a matching pre-transition current frame: " *
        "kind=$(power_semiconductor_kind(device)), nodes=($(device.a),$(device.b)), time=$(time_s), " *
        "topology=$(topology_matches) (frame=$(frame.topology_transition_count), device=$(device.topology_transition_count)), " *
        "rounded_voltage=$(voltage_matches) (frame=$(Float64(frame.voltage_v)), report=$(device.last_voltage)), " *
        "rounded_current=$(current_matches) (frame=$(Float64(frame.terminal_current_a)), report=$(device.last_current))",
    ))
    return SemiconductorEventBoundaryContext(device, time_s,
        power_semiconductor_topology_chart_state(device), frame,
        fidelity.accepted_topology_transition_count)
end

function NonlinearNetwork.refresh_nonlinear_device_event_dependents!(
    device::PowerSemiconductorSwitch,
    time_s::Union{Float64,ContinuousEvaluationTime},
    context::SemiconductorEventBoundaryContext,
)
    fidelity = something(device.extended_fidelity)
    frame = context.frame
    current = fidelity.boundary_electrical_state === nothing ?
        fidelity.accepted_electrical_history : fidelity.boundary_electrical_state
    changed = context.chart.forward_conducting != device.closed ||
        context.chart.reverse_diode_conducting != device.reverse_diode_conducting
    context.device === device && isequal(context.time_s, time_s) && current === frame &&
        device.topology_transition_count == context.chart.topology_transition_count + Int(changed) &&
        fidelity.accepted_topology_transition_count == context.accepted_topology_transition_count ||
        throw(ArgumentError("device event refresh lost its pre-transition boundary context"))
    event_time = time_s isa ContinuousEvaluationTime ? time_s.high_s : time_s
    !changed || device.last_transition_time_s == event_time || throw(ArgumentError(
        "device event topology transition has a different physical clock",
    ))
    channel = frame.arithmetic_policy === :binary64 ?
        first(_power_semiconductor_accepted_boundary_channel_current_jacobian(device, Float64(frame.voltage_v))) :
        setprecision(BigFloat, 106) do
            Float64(first(_power_semiconductor_compensated_channel_current_jacobian(device, frame.voltage_v)))
        end
    displacement = fidelity.junction_charge === nothing ? 0.0 : fidelity.junction_charge.last_displacement_current_a
    tail = fidelity.turn_off_tail === nothing || !fidelity.turn_off_tail.active ? 0.0 : fidelity.turn_off_tail.current_a
    snubber = device.snubber === nothing ? 0.0 : device.snubber.last_current_a
    physical_current = setprecision(BigFloat, 106) do
        if frame.arithmetic_policy === :binary64
            BigFloat(Float64(frame.terminal_current_a) + (channel - frame.channel_current_a) +
                (displacement - frame.displacement_current_a) + (tail - frame.tail_current_a) +
                (snubber - frame.snubber_current_a))
        else
            frame.arithmetic_policy === :compensated_binary || throw(ArgumentError("unknown event boundary arithmetic"))
            frame.terminal_current_a + (BigFloat(channel) - BigFloat(frame.channel_current_a)) +
                (BigFloat(displacement) - BigFloat(frame.displacement_current_a)) +
                (BigFloat(tail) - BigFloat(frame.tail_current_a)) +
                (BigFloat(snubber) - BigFloat(frame.snubber_current_a))
        end
    end
    post_event = AcceptedSemiconductorElectricalHistory(frame.voltage_v, physical_current,
        channel, displacement, tail, snubber, device.topology_transition_count, frame.arithmetic_policy)
    # Check the rebuilt current chart before publishing it. The enclosing event
    # transaction also owns the preceding physical transition and any projection.
    candidate = deepcopy(device)
    candidate.extended_fidelity.boundary_electrical_anchor = post_event
    candidate.extended_fidelity.boundary_electrical_state = post_event
    _reconstruct_power_semiconductor_boundary_state!(candidate, time_s, frame.voltage_v)
    fidelity.boundary_electrical_anchor = post_event
    fidelity.boundary_electrical_state = post_event
    _reconstruct_power_semiconductor_boundary_state!(device, time_s, frame.voltage_v)
    return nothing
end

function refresh_nonlinear_device_event_dependents!(
    device::PowerSemiconductorSwitch,
    time_s::Union{Float64,ContinuousEvaluationTime},
)
    device.extended_fidelity === nothing && return nothing
    frame=device.extended_fidelity.boundary_electrical_state===nothing ? device.extended_fidelity.accepted_electrical_history : device.extended_fidelity.boundary_electrical_state
    physical_voltage=frame===nothing ? device.last_voltage : frame.voltage_v
    reconstruct_power_semiconductor_boundary_state!(device,time_s,physical_voltage)
    return nothing
end

function NonlinearNetwork.nonlinear_device_event_capacitance(
    device::PowerSemiconductorSwitch,
)
    fidelity = device.extended_fidelity
    fidelity === nothing && return nothing
    junction = fidelity.junction_charge
    junction === nothing && return nothing
    return (device.a, device.b, junction.last_capacitance_f)
end

NonlinearNetwork.nonlinear_device_event_current(device::PowerSemiconductorSwitch) =
    device.last_current

function NonlinearNetwork.rebase_nonlinear_device_event_derivative!(
    device::PowerSemiconductorSwitch,
    current_increment_a::Float64,
)
    fidelity = device.extended_fidelity
    fidelity === nothing && throw(ArgumentError(
        "semiconductor event derivative rebase requires extended fidelity",
    ))
    junction = fidelity.junction_charge
    junction === nothing && throw(ArgumentError(
        "semiconductor event derivative rebase requires junction charge",
    ))
    isfinite(current_increment_a) || throw(ArgumentError(
        "semiconductor event derivative increment must be finite",
    ))
    displacement=junction.last_displacement_current_a+current_increment_a
    frame=fidelity.boundary_electrical_state===nothing ? fidelity.accepted_electrical_history : fidelity.boundary_electrical_state
    if frame!==nothing
        anchor=fidelity.boundary_electrical_anchor===nothing ? fidelity.accepted_electrical_history : fidelity.boundary_electrical_anchor
        channel = frame.arithmetic_policy === :binary64 ?
            first(_power_semiconductor_accepted_boundary_channel_current_jacobian(device, Float64(frame.voltage_v))) :
            setprecision(BigFloat, 106) do
                Float64(first(_power_semiconductor_compensated_channel_current_jacobian(device, frame.voltage_v)))
            end
        frame.topology_transition_count==device.topology_transition_count &&
            anchor!==nothing && anchor.topology_transition_count==frame.topology_transition_count &&
            fidelity.accepted_topology_transition_count<=frame.topology_transition_count &&
            frame.arithmetic_policy in (:binary64, :compensated_binary) &&
            frame.channel_current_a==channel &&
            frame.displacement_current_a==junction.last_displacement_current_a &&
            Float64(frame.voltage_v)==device.last_voltage &&
            Float64(frame.terminal_current_a)==device.last_current ||
            throw(ArgumentError("event derivative rebase requires matching topology and accepted physical frame"))
    end
    updated_frame=if frame===nothing
        nothing
    else
        physical_current=setprecision(BigFloat,106) do
            frame.arithmetic_policy===:binary64 ? BigFloat(Float64(frame.terminal_current_a)+current_increment_a) :
                frame.terminal_current_a+BigFloat(current_increment_a)
        end
        AcceptedSemiconductorElectricalHistory(frame.voltage_v,physical_current,
            frame.channel_current_a,displacement,frame.tail_current_a,frame.snubber_current_a,frame.topology_transition_count,frame.arithmetic_policy)
    end
    terminal=frame===nothing ? device.last_current+current_increment_a : Float64(updated_frame.terminal_current_a)
    previous=frame===nothing ? fidelity.previous_terminal_current_a+current_increment_a : terminal
    all(isfinite,(displacement,terminal,previous)) || throw(ArgumentError("event derivative rebase overflowed physical history"))
    junction.last_displacement_current_a=displacement
    device.last_current=terminal
    fidelity.previous_terminal_current_a=previous
    fidelity.boundary_electrical_anchor=updated_frame
    fidelity.boundary_electrical_state=updated_frame
    return nothing
end

function _power_semiconductor_boundary_maximum_error(current, expected)
    if current isa AbstractArray || current isa Tuple
        length(current) == length(expected) || return Inf
        return maximum((
            _power_semiconductor_boundary_maximum_error(left, right)
            for (left, right) in zip(current, expected)
        ); init=0.0)
    end
    return abs(Float64(current) - Float64(expected))
end

"""Return exact errors for every reconstructed semiconductor boundary family."""
function power_semiconductor_complete_boundary_state_errors(
    device::PowerSemiconductorSwitch,
    boundary_time_s::Real,
    terminal_voltage_v::Real,
)
    expected = deepcopy(device)
    reconstruct_power_semiconductor_boundary_state!(
        expected,
        boundary_time_s,
        terminal_voltage_v,
    )
    fidelity = device.extended_fidelity
    expected_fidelity = expected.extended_fidelity
    companion = maximum((
        abs(device.last_conductance - expected.last_conductance),
        abs(device.last_history_current_a - expected.last_history_current_a),
        abs(device.last_evaluation_time_s - expected.last_evaluation_time_s),
    ))
    terminal = maximum((
        abs(device.last_voltage - expected.last_voltage),
        abs(device.last_current - expected.last_current),
        abs(device.last_forward_current_a - expected.last_forward_current_a),
        abs(device.last_reverse_diode_current_a -
            expected.last_reverse_diode_current_a),
        abs(device.last_snubber_current_a - expected.last_snubber_current_a),
        fidelity === nothing ? 0.0 : abs(
            fidelity.previous_terminal_voltage_v -
            something(expected_fidelity).previous_terminal_voltage_v,
        ),
        fidelity === nothing ? 0.0 : abs(
            fidelity.previous_terminal_current_a -
            something(expected_fidelity).previous_terminal_current_a,
        ),
    ))
    junction = if fidelity === nothing || fidelity.junction_charge === nothing
        0.0
    else
        current = something(fidelity.junction_charge)
        target = something(something(expected_fidelity).junction_charge)
        maximum((
            abs(current.previous_voltage_v - target.previous_voltage_v),
            abs(current.previous_charge_c - target.previous_charge_c),
            abs(current.last_capacitance_f - target.last_capacitance_f),
            abs(current.last_charge_c - target.last_charge_c),
            abs(current.last_displacement_current_a -
                target.last_displacement_current_a),
        ))
    end
    snubber = if device.snubber === nothing
        0.0
    else
        current = something(device.snubber)
        target = something(expected.snubber)
        maximum((
            abs(current.previous_current_a - target.previous_current_a),
            abs(current.capacitor_voltage_v - target.capacitor_voltage_v),
            abs(current.last_branch_voltage_v - target.last_branch_voltage_v),
            abs(current.last_current_a - target.last_current_a),
            abs(current.last_resistor_loss_w - target.last_resistor_loss_w),
        ))
    end
    thermal = if fidelity === nothing || fidelity.thermal === nothing
        0.0
    else
        current = something(fidelity.thermal)
        target = something(something(expected_fidelity).thermal)
        maximum((
            _power_semiconductor_boundary_maximum_error(
                current.trial_temperature_k,
                target.trial_temperature_k,
            ),
            abs(current.last_stored_energy_j - target.last_stored_energy_j),
            abs(current.last_ambient_heat_flow_w -
                target.last_ambient_heat_flow_w),
        ))
    end
    return (; companion, terminal, junction, snubber, thermal)
end

"""Return exact invariant errors for dependent semiconductor boundary state."""
function power_semiconductor_boundary_state_errors(
    device::PowerSemiconductorSwitch,
    boundary_time_s::Real,
)
    boundary_time = Float64(boundary_time_s)
    isfinite(boundary_time) || throw(ArgumentError(
        "a semiconductor boundary time must be finite",
    ))
    fidelity = device.extended_fidelity
    fidelity === nothing && return (
        recovery_start_presence=0.0,
        recovery_duration=0.0,
        tail_turn_off_presence=0.0,
        tail_current=0.0,
        tail_duration=0.0,
    )
    recovery = fidelity.recovered_charge
    recovery_start_presence = recovery !== nothing &&
        recovery.recovery_active && !isfinite(recovery.recovery_start_time_s) ?
        1.0 : 0.0
    recovery_duration = if recovery !== nothing && recovery.recovery_active &&
                           isfinite(recovery.recovery_start_time_s)
        abs(recovery.last_recovery_duration_s - Float64(max(
            0.0,
            _semiconductor_history_elapsed(boundary_time, recovery.recovery_start_time_s),
        )))
    else
        0.0
    end
    tail = fidelity.turn_off_tail
    tail_turn_off_presence = tail !== nothing && tail.active &&
        !isfinite(tail.turn_off_time_s) ? 1.0 : 0.0
    expected_tail_current = if tail !== nothing && tail.active &&
                               isfinite(tail.turn_off_time_s)
        tail.initial_current_a * exp(
            -max(0.0, _semiconductor_tail_elapsed(boundary_time, tail.turn_off_time_s)) /
                tail.decay_time_s,
        )
    else
        0.0
    end
    tail_current = tail === nothing ? 0.0 :
        abs(tail.current_a - expected_tail_current)
    tail_duration = if tail !== nothing && tail.active &&
                       isfinite(tail.turn_off_time_s)
        abs(tail.last_duration_s - max(
            0.0,
            _semiconductor_tail_elapsed(boundary_time, tail.turn_off_time_s),
        ))
    else
        0.0
    end
    return (;
        recovery_start_presence,
        recovery_duration,
        tail_turn_off_presence,
        tail_current,
        tail_duration,
    )
end

function _power_semiconductor_fidelity_components(
    fidelity::PowerSemiconductorExtendedFidelity,
)
    components = Symbol[]
    fidelity.recovered_charge === nothing || push!(components, :recovered_charge)
    fidelity.junction_charge === nothing || push!(components, :nonlinear_junction_charge)
    fidelity.turn_off_tail === nothing || push!(components, :turn_off_tail)
    fidelity.switching_energy === nothing || push!(components, :event_energy_map)
    fidelity.thermal === nothing || push!(components, :lumped_electrothermal)
    fidelity.declared_model === nothing || push!(components, :declared_compact_model)
    return Tuple(components)
end

function power_semiconductor_extended_state(device::PowerSemiconductorSwitch)
    fidelity = device.extended_fidelity
    fidelity === nothing && throw(ArgumentError(
        "baseline power semiconductor has no extended state",
    ))
    recovery = fidelity.recovered_charge
    junction = fidelity.junction_charge
    tail = fidelity.turn_off_tail
    energy = fidelity.switching_energy
    thermal = fidelity.thermal
    return PowerSemiconductorExtendedState(
        1,
        _power_semiconductor_fidelity_components(fidelity),
        device.conduction_direction,
        device.gate_turn_off_disposition,
        recovery === nothing ? false : recovery.recovery_active,
        recovery === nothing ? 0.0 : recovery.stored_charge_c,
        recovery === nothing ? 0.0 : recovery.last_recovery_current_a,
        recovery === nothing ? 0.0 : recovery.peak_reverse_current_a,
        recovery === nothing ? 0.0 : recovery.cumulative_recovered_charge_c,
        recovery === nothing ? 0 : recovery.recovery_zero_event_count,
        recovery === nothing ? -Inf : recovery.last_recovery_zero_time_s,
        junction === nothing ? 0.0 : junction.last_capacitance_f,
        junction === nothing ? 0.0 : junction.last_charge_c,
        junction === nothing ? 0.0 : junction.last_displacement_current_a,
        junction === nothing ? 0.0 :
            power_semiconductor_junction_stored_energy(junction, junction.last_charge_c == junction.previous_charge_c ? junction.previous_voltage_v : device.last_voltage),
        fidelity.companion_energy_residual_j,
        tail === nothing ? false : tail.active,
        tail === nothing ? 0.0 : tail.current_a,
        tail === nothing ? 0 : tail.cutoff_event_count,
        tail === nothing ? -Inf : tail.last_cutoff_time_s,
        energy === nothing ? :none : energy.last_event_kind,
        energy === nothing ? 0.0 : energy.last_event_energy_j,
        energy === nothing ? 0.0 : energy.cumulative_turn_on_energy_j,
        energy === nothing ? 0.0 : energy.cumulative_turn_off_energy_j,
        energy === nothing ? 0.0 : energy.cumulative_reverse_recovery_energy_j,
        thermal === nothing ? 0.0 : thermal.node_temperature_k[1],
        thermal === nothing ? Float64[] : copy(thermal.node_temperature_k),
        thermal === nothing ? 0.0 : thermal.last_ambient_heat_flow_w,
        thermal === nothing ? 0.0 : power_semiconductor_thermal_stored_energy(thermal),
    )
end

trace_output_channel_count(::PowerSemiconductorSwitch) = 10
trace_output_is_public(::PowerSemiconductorSwitch) = true

function trace_output_channel_names!(
    names::Vector{Symbol},
    element_name::Symbol,
    ::PowerSemiconductorSwitch,
)
    append!(
        names,
        Symbol[
            Symbol(element_name, :_terminal_voltage_v),
            Symbol(element_name, :_terminal_current_a),
            Symbol(element_name, :_forward_current_a),
            Symbol(element_name, :_reverse_diode_current_a),
            Symbol(element_name, :_gate_commanded),
            Symbol(element_name, :_gate_applied),
            Symbol(element_name, :_semiconductor_loss_w),
            Symbol(element_name, :_snubber_resistor_loss_w),
            Symbol(element_name, :_semiconductor_dissipated_energy_j),
            Symbol(element_name, :_snubber_dissipated_energy_j),
        ],
    )
    return names
end

function trace_output_values!(
    output::AbstractMatrix{Float64},
    first_channel::Int,
    sample::Int,
    device::PowerSemiconductorSwitch,
    _voltage::AbstractVector{Float64},
)
    terminal = power_semiconductor_terminal_state(device)
    values = (
        terminal.terminal_voltage_v,
        terminal.terminal_current_a,
        terminal.forward_current_a,
        terminal.reverse_diode_current_a,
        terminal.gate_commanded_on ? 1.0 : 0.0,
        terminal.gate_applied_on ? 1.0 : 0.0,
        terminal.semiconductor_loss_w,
        terminal.snubber_resistor_loss_w,
        terminal.semiconductor_dissipated_energy_j,
        terminal.snubber_dissipated_energy_j,
    )
    for offset in eachindex(values)
        output[first_channel + offset - 1, sample] = values[offset]
    end
    return first_channel + length(values)
end

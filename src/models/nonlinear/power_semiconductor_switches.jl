abstract type PowerSemiconductorJunction end

struct DiodeJunction <: PowerSemiconductorJunction end
struct ThyristorJunction <: PowerSemiconductorJunction end
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

"""Trapezoidal companion state for a series resistance-capacitance snubber connected across a semiconductor device."""
mutable struct SeriesRCSnubber
    resistance_ohm::Float64
    capacitance_f::Float64
    previous_current_a::Float64
    capacitor_voltage_v::Float64
    last_branch_voltage_v::Float64
    last_current_a::Float64
    last_resistor_loss_w::Float64
    dissipated_energy_j::Float64
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
Evidence-bounded piecewise-linear diode, thyristor, IGBT, or MOSFET electrical switch. The baseline owns gating, latching/extinction, optional antiparallel conduction, a series-RC snubber, and conduction/snubber energy; it deliberately excludes reverse recovery, nonlinear capacitance, tail current, package parasitics, switching-energy maps, and thermal state.

Validity requires an EMT timestep fine enough to resolve the external circuit, localized switching surfaces, commanded gate timing, and series-RC dynamics. The model claims no separate frequency-domain limit; reported dissipated energy integrates piecewise-linear terminal conduction and snubber resistance only, so it is not a manufacturer switching-loss or thermal estimate.
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
end

"""Piecewise-linear naturally commutated diode specialization."""
const DiodeValveSwitch = PowerSemiconductorSwitch{DiodeJunction}
"""Piecewise-linear gate-triggered, holding-current-commutated thyristor specialization."""
const ThyristorValveSwitch = PowerSemiconductorSwitch{ThyristorJunction}
"""Piecewise-linear forward-channel IGBT specialization with an optional antiparallel diode."""
const IGBTSwitch = PowerSemiconductorSwitch{InsulatedGateBipolarTransistorJunction}
"""Piecewise-linear bidirectional-on-channel MOSFET specialization with an optional body-diode path."""
const MOSFETSwitch = PowerSemiconductorSwitch{MetalOxideSemiconductorFieldEffectTransistorJunction}

power_semiconductor_kind(::PowerSemiconductorSwitch{DiodeJunction}) = :diode
power_semiconductor_kind(::PowerSemiconductorSwitch{ThyristorJunction}) = :thyristor
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
    device::PowerSemiconductorSwitch{InsulatedGateBipolarTransistorJunction},
    time_s::Float64,
)
    previous_forward = device.closed
    previous_reverse = device.reverse_diode_conducting
    if device.gate_driver.applied_on &&
       !device.reverse_diode_conducting &&
       device.last_voltage >= power_semiconductor_forward_turn_on_voltage(device)
        device.closed = true
    elseif !device.gate_driver.applied_on
        device.closed = false
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
    transition_time = if commanded_on
        max(
            command_time + driver.turn_on_delay_s,
            driver.last_turn_off_time_s + driver.dead_time_s,
        )
    else
        max(
            command_time + driver.turn_off_delay_s,
            driver.last_turn_on_time_s + driver.minimum_pulse_width_s,
        )
    end
    if earliest_transition_time !== nothing
        transition_time = max(transition_time, earliest_transition_time)
    end
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

function apply_power_semiconductor_gate_transition!(
    device::PowerSemiconductorSwitch,
    time_s::Real,
)
    driver = device.gate_driver
    driver === nothing && throw(ArgumentError("a diode junction has no gate transition"))
    driver.pending_state === nothing && return false
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
    _apply_power_semiconductor_gate_state!(device, transition_time)
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
    _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
    )
    return device
end

function apply_power_semiconductor_forward_extinction!(
    device::PowerSemiconductorSwitch,
    time_s::Real,
)
    previous_forward = device.closed
    previous_reverse = device.reverse_diode_conducting
    device.closed = false
    _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
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

function stamp!(
    admittance::AbstractMatrix{Float64},
    rhs::AbstractVector{Float64},
    device::PowerSemiconductorSwitch,
    time_s::Float64,
    dt_s::Float64,
)
    _stamp_power_semiconductor_gate_due!(device, time_s)
    device.last_evaluation_time_s = time_s
    conductance = device.off_conductance
    history_current = 0.0
    if device.closed
        conductance = device.on_conductance
        history_current = -conductance * device.forward_voltage_drop_v
    elseif device.reverse_diode_conducting
        diode = something(device.antiparallel_diode)
        conductance = diode.on_conductance_s
        history_current = conductance * diode.forward_voltage_v
    end
    device.last_conductance = conductance
    device.last_history_current_a = history_current
    stamp_conductance!(admittance, device.a, device.b, conductance)
    stamp_history_current!(rhs, device.a, device.b, history_current)
    if device.snubber !== nothing
        snubber_conductance, snubber_history_current =
            _series_rc_snubber_companion(device.snubber, dt_s)
        stamp_conductance!(admittance, device.a, device.b, snubber_conductance)
        stamp_history_current!(rhs, device.a, device.b, snubber_history_current)
    end
    return nothing
end

function _update_series_rc_snubber!(
    snubber::SeriesRCSnubber,
    terminal_voltage_v::Float64,
    dt_s::Float64,
)
    previous_current = snubber.previous_current_a
    previous_loss = snubber.last_resistor_loss_w
    conductance, history_current = _series_rc_snubber_companion(snubber, dt_s)
    current = conductance * terminal_voltage_v + history_current
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

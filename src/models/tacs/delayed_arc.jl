"""
Mutable state for the controlled-switch delayed-arc continuation.

The four immutable parameters correspond to the accepted `A8SW` continuation:
current coefficient, current exponent, time scale, and cutoff current. Mutable
fields own zero-crossing prediction and the decaying post-opening current.
"""
mutable struct ControlledSwitchDelayedArcState
    current_coefficient::Float64
    current_exponent::Float64
    time_scale_s::Float64
    cutoff_current_a::Float64
    previous_current_a::Float64
    tail_amplitude_a::Float64
    tail_time_constant_s::Float64
    scheduled_open_time_s::Float64
    opening_requested::Bool
    tail_active::Bool
    tail_current_a::Float64
    decay_factor::Float64
    transition_count::Int
    absolute_tail_energy_j::Float64
end

function ControlledSwitchDelayedArcState(
    current_coefficient::Real,
    current_exponent::Real,
    time_scale_s::Real,
    cutoff_current_a::Real,
)
    coefficient = Float64(current_coefficient)
    exponent = Float64(current_exponent)
    time_scale = Float64(time_scale_s)
    cutoff = Float64(cutoff_current_a)
    coefficient >= 0.0 && isfinite(coefficient) ||
        throw(ArgumentError("delayed-arc current coefficient must be finite and nonnegative"))
    isfinite(exponent) ||
        throw(ArgumentError("delayed-arc current exponent must be finite"))
    time_scale > 0.0 && isfinite(time_scale) ||
        throw(ArgumentError("delayed-arc time scale must be finite and positive"))
    cutoff >= 0.0 && isfinite(cutoff) ||
        throw(ArgumentError("delayed-arc cutoff current must be finite and nonnegative"))
    return ControlledSwitchDelayedArcState(
        coefficient,
        exponent,
        time_scale,
        cutoff,
        0.0,
        0.0,
        1.0,
        Inf,
        false,
        false,
        0.0,
        0.0,
        0,
        0.0,
    )
end

function _defer_controlled_switch_arc_opening!(
    switch,
    previously_closed::Bool,
)
    arc = switch.delayed_arc
    if arc !== nothing && previously_closed && !switch.closed
        switch.closed = true
        arc.opening_requested = true
    end
    return switch
end

function _advance_controlled_switch_delayed_arc!(
    switch,
    time_s::Float64,
    dt_s::Float64,
)
    arc = switch.delayed_arc
    (arc === nothing || !arc.opening_requested || !switch.closed) &&
        return nothing
    update = over16_a8sw_delayed_open_step(
        switch.last_current,
        arc.previous_current_a,
        arc.tail_amplitude_a,
        arc.tail_time_constant_s,
        arc.scheduled_open_time_s,
        time_s,
        dt_s,
        arc.current_coefficient,
        arc.current_exponent,
        arc.time_scale_s,
    )
    arc.previous_current_a = update.previous_current
    arc.tail_amplitude_a = update.shape_current
    arc.tail_time_constant_s = update.shape_delay
    arc.scheduled_open_time_s = update.scheduled_open_time
    if update.opens
        switch.closed = false
        arc.opening_requested = false
        arc.tail_active = true
        arc.transition_count += 1
    end
    return update
end

function _stamp_controlled_switch_arc_tail!(
    rhs::AbstractVector{Float64},
    switch,
    time_s::Float64,
    dt_s::Float64,
)
    arc = switch.delayed_arc
    arc === nothing && return 0.0
    if !arc.tail_active
        arc.tail_current_a = 0.0
        arc.decay_factor = 0.0
        return 0.0
    end
    current_from, current_to, decay_factor, clear_tail =
        over16_switch_tail_current_injection(
            arc.scheduled_open_time_s,
            time_s,
            dt_s,
            arc.tail_amplitude_a,
            arc.tail_time_constant_s,
            arc.cutoff_current_a,
        )
    switch.a == 0 || (rhs[switch.a] += current_from)
    switch.b == 0 || (rhs[switch.b] += current_to)
    arc.tail_current_a = current_from
    arc.decay_factor = decay_factor
    if clear_tail
        arc.tail_active = false
    end
    return current_from
end

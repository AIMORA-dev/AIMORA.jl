export SWITCHING_NONLINEAR_RESISTOR_TYPE,
       SwitchingNonlinearResistorStep,
       switching_nonlinear_resistor_step

const SWITCHING_NONLINEAR_RESISTOR_TYPE = -99

struct SwitchingNonlinearResistorStep
    current_segment::Float64
    companion_current_a::Float64
    turn_on_voltage_v::Float64
    rearm_time_s::Float64
    reported_current_a::Float64
    accepted_current_a::Float64
    admittance_delta_s::Float64
    companion_current_delta_a::Float64
    segment_changed::Bool
    activated::Bool
    deactivated::Bool
    polarity_reversed::Bool
end

"""
Advance one symmetric piecewise-linear switching resistor.

`cchar` contains segment current intercepts in ampere, `gslope` segment
conductances in siemens, and `vchar` the positive-voltage upper boundary of
each segment. The returned admittance and companion-current deltas are applied
to the nodal matrix and RHS in that order.
"""
function switching_nonlinear_resistor_step(
    current_segment::Real,
    companion_current_a::Real,
    branch_voltage_v::Real,
    table_start::Int,
    table_end::Int,
    cchar::AbstractVector{<:Real},
    gslope::AbstractVector{<:Real},
    vchar::AbstractVector{<:Real};
    turn_on_voltage_v::Real,
    turn_off_voltage_v::Real,
    activation_segment_count::Int,
    minimum_on_time_s::Real,
    rearm_time_s::Real,
    time_s::Real,
    single_flash::Bool=false,
    infinity::Real=Inf,
    voltage_tolerance::Real=0.0,
)
    1 <= table_start <= table_end <= length(cchar) ||
        throw(ArgumentError("switching resistor table range must address cchar"))
    table_end <= length(gslope) && table_end <= length(vchar) ||
        throw(ArgumentError("switching resistor table range must address gslope and vchar"))
    segment_count = table_end - table_start + 1
    1 <= activation_segment_count <= segment_count ||
        throw(ArgumentError("activation segment must lie inside the switching resistor table"))
    segment_value = Float64(current_segment)
    isinteger(segment_value) ||
        throw(ArgumentError("switching resistor current segment must be integral"))
    segment = Int(segment_value)
    abs(segment) <= segment_count ||
        throw(ArgumentError("switching resistor current segment is outside its table"))
    voltage = Float64(branch_voltage_v)
    companion = Float64(companion_current_a)
    turn_on = Float64(turn_on_voltage_v)
    turn_off = Float64(turn_off_voltage_v)
    minimum_on_time = Float64(minimum_on_time_s)
    deadline = Float64(rearm_time_s)
    time = Float64(time_s)
    tolerance = Float64(voltage_tolerance)
    all(isfinite, (voltage, companion, turn_off, minimum_on_time, deadline, time, tolerance)) ||
        throw(ArgumentError("switching resistor finite inputs contain a non-finite value"))
    (isfinite(turn_on) || turn_on == Inf) && turn_on >= 0.0 ||
        throw(ArgumentError("switching resistor turn-on voltage must be nonnegative"))
    turn_off >= 0.0 && minimum_on_time >= 0.0 && tolerance >= 0.0 ||
        throw(ArgumentError("switching resistor limits must be nonnegative"))
    infinity_value = Float64(infinity)
    infinity_value > 0.0 || throw(ArgumentError("switching resistor infinity must be positive"))

    voltage_sign = voltage < 0.0 ? -1 : 1
    voltage_magnitude = abs(voltage)
    function selected_segment_magnitude(seed::Int)
        magnitude = seed
        while magnitude < segment_count &&
              voltage_magnitude > Float64(vchar[table_start + magnitude - 1])
            magnitude += 1
        end
        while magnitude > 1 &&
              voltage_magnitude < Float64(vchar[table_start + magnitude - 2])
            magnitude -= 1
        end
        return magnitude
    end
    if segment == 0
        if voltage_magnitude <= turn_on
            return SwitchingNonlinearResistorStep(
                0.0, 0.0, turn_on, deadline, 0.0, 0.0, 0.0, -companion,
                false, false, false, false,
            )
        end
        activated_segment = selected_segment_magnitude(activation_segment_count)
        index = table_start + activated_segment - 1
        new_companion = voltage_sign * Float64(cchar[index])
        accepted_current = Float64(gslope[index]) * voltage + new_companion
        return SwitchingNonlinearResistorStep(
            Float64(voltage_sign * activated_segment),
            new_companion,
            turn_on,
            time + minimum_on_time,
            0.0,
            accepted_current,
            Float64(gslope[index]),
            new_companion - companion,
            true,
            true,
            false,
            false,
        )
    end

    magnitude = abs(segment)
    old_index = table_start + magnitude - 1
    old_slope = Float64(gslope[old_index])
    reported_current = old_slope * voltage + companion
    old_sign = sign(segment)
    polarity_reversed = old_sign != voltage_sign
    if time >= deadline && voltage_magnitude < turn_off
        disable_forever =
            single_flash && magnitude == 1 &&
            abs(Float64(vchar[table_start]) - Float64(vchar[old_index])) <= tolerance
        next_turn_on = disable_forever ? infinity_value : turn_on
        return SwitchingNonlinearResistorStep(
            0.0,
            0.0,
            next_turn_on,
            deadline,
            reported_current,
            0.0,
            -old_slope,
            -companion,
            true,
            false,
            true,
            polarity_reversed,
        )
    end

    if polarity_reversed && time >= deadline && voltage_magnitude <= turn_on
        next_turn_on = single_flash ? infinity_value : turn_on
        return SwitchingNonlinearResistorStep(
            0.0, 0.0, next_turn_on, deadline, reported_current, 0.0,
            -old_slope, -companion, true, false, true, true,
        )
    end

    new_magnitude = selected_segment_magnitude(magnitude)
    new_index = table_start + new_magnitude - 1
    new_slope = Float64(gslope[new_index])
    new_companion = voltage_sign * Float64(cchar[new_index])
    accepted_current = new_slope * voltage + new_companion
    new_segment = voltage_sign * new_magnitude
    return SwitchingNonlinearResistorStep(
        Float64(new_segment),
        new_companion,
        turn_on,
        deadline,
        reported_current,
        accepted_current,
        new_slope - old_slope,
        new_companion - companion,
        new_segment != segment,
        false,
        false,
        polarity_reversed,
    )
end

export TRIGGERED_TIMED_RESISTANCE_TYPE,
       TriggeredTimedResistanceStep,
       triggered_timed_resistance_step

const TRIGGERED_TIMED_RESISTANCE_TYPE = -97

struct TriggeredTimedResistanceStep
    active_segment::Int
    activation_time_s::Float64
    conductance_s::Float64
    admittance_delta_s::Float64
    accepted_current_a::Float64
    activated::Bool
    segment_change_count::Int
end

"""
Advance a voltage-triggered staircase resistance schedule.

The schedule contains elapsed times in seconds and strictly positive
resistances represented by their conductances. Once armed and triggered, the
first resistance is connected. All schedule points due at the accepted time
are applied in order; their combined nodal update is the final-minus-initial
conductance, which is algebraically identical to repeated restamping.
"""
function triggered_timed_resistance_step(
    active_segment::Integer,
    activation_time_s::Real,
    branch_voltage_v::Real,
    table_start::Int,
    table_end::Int,
    elapsed_times_s::AbstractVector{<:Real},
    conductances_s::AbstractVector{<:Real};
    trigger_voltage_v::Real,
    arm_time_s::Real,
    time_s::Real,
)
    1 <= table_start <= table_end <= length(elapsed_times_s) ||
        throw(ArgumentError("timed-resistance table range must address elapsed times"))
    table_end <= length(conductances_s) ||
        throw(ArgumentError("timed-resistance table range must address conductances"))
    segment_count = table_end - table_start + 1
    segment = Int(active_segment)
    0 <= segment <= segment_count ||
        throw(ArgumentError("timed-resistance active segment must address its schedule"))
    voltage = Float64(branch_voltage_v)
    trigger = Float64(trigger_voltage_v)
    arm = Float64(arm_time_s)
    time = Float64(time_s)
    activation = Float64(activation_time_s)
    all(isfinite, (voltage, trigger, arm, time, activation)) ||
        throw(ArgumentError("timed-resistance state and inputs must be finite"))
    trigger >= 0.0 || throw(ArgumentError("trigger voltage must be nonnegative"))

    previous_elapsed = -Inf
    for table_index in table_start:table_end
        elapsed = Float64(elapsed_times_s[table_index])
        conductance = Float64(conductances_s[table_index])
        elapsed >= 0.0 && elapsed > previous_elapsed ||
            throw(ArgumentError("timed-resistance elapsed times must be nonnegative and increasing"))
        conductance > 0.0 && isfinite(conductance) ||
            throw(ArgumentError("timed-resistance conductances must be finite and positive"))
        previous_elapsed = elapsed
    end

    old_conductance = segment == 0 ? 0.0 :
        Float64(conductances_s[table_start + segment - 1])
    activated = false
    changes = 0
    if segment == 0
        if time >= arm && abs(voltage) > trigger
            segment = 1
            activation = time
            activated = true
            changes = 1
        end
    end
    if segment > 0
        elapsed_since_activation = time - activation
        while segment < segment_count &&
              elapsed_since_activation >=
              Float64(elapsed_times_s[table_start + segment])
            segment += 1
            changes += 1
        end
    end
    conductance = segment == 0 ? 0.0 :
        Float64(conductances_s[table_start + segment - 1])
    return TriggeredTimedResistanceStep(
        segment,
        activation,
        conductance,
        conductance - old_conductance,
        conductance * voltage,
        activated,
        changes,
    )
end

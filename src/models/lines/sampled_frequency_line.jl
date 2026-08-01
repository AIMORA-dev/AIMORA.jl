export LineWeightingSamples,
       SampledLineWeightingCoefficients,
       SampledFrequencyDependentLine,
       SampledFrequencyDependentLineGroup,
       line_weighting_samples,
       sampled_line_weighting_coefficients,
       sampled_frequency_dependent_line,
       sampled_frequency_dependent_line_group,
       sampled_line_steady_state_terminal_admittance,
       sampled_line_history_convolution,
       sampled_line_history_convolution!

import ..Branches: trace_output_channel_count,
                   trace_output_channel_names!,
                   trace_output_values!

struct LineWeightingSamples
    time_s::Vector{Float64}
    amplitude::Vector{Float64}
end

function line_weighting_samples(time_s::AbstractVector, amplitude::AbstractVector)
    times = Float64.(time_s)
    values = Float64.(amplitude)
    length(times) == length(values) ||
        throw(ArgumentError("line weighting times and amplitudes must have equal length"))
    length(times) >= 2 ||
        throw(ArgumentError("a line weighting function requires at least two samples"))
    all(isfinite, times) && all(isfinite, values) ||
        throw(ArgumentError("line weighting samples must be finite"))
    all(values .>= 0.0) ||
        throw(ArgumentError("passive line weighting amplitudes must be nonnegative"))
    all(diff(times) .> 0.0) ||
        throw(ArgumentError("line weighting sample times must be strictly increasing"))
    times[1] >= 0.0 ||
        throw(ArgumentError("line weighting sample times must be nonnegative"))
    return LineWeightingSamples(times, values)
end

struct SampledLineWeightingCoefficients
    timestep_s::Float64
    characteristic_impedance_ohm::Float64
    eta::Float64
    propagation_delay_steps::Int
    history_sample_count::Int
    admittance_weights::Vector{Float64}
    propagation_weights::Vector{Float64}
    admittance_tail::NTuple{3,Float64}
    propagation_tail::NTuple{3,Float64}
    input_integral::Float64
    normalized::Bool
    total_resistance_ohm::Float64
    tail_relative_tolerance::Float64
    coefficient_sum::Float64
    dc_gain::Float64
end

function _trapezoid_integral(samples::LineWeightingSamples)
    total = 0.0
    for index in 2:length(samples.time_s)
        total += 0.5 *
                 (samples.time_s[index] - samples.time_s[index - 1]) *
                 (samples.amplitude[index] + samples.amplitude[index - 1])
    end
    return total
end

function _discrete_tail_coefficients(raw_tail::NTuple{3,Float64}, timestep_s::Float64)
    tail_amplitude, tail_rate, tail_start = raw_tail
    tail_start == 0.0 && return (0.0, 0.0, 0.0)
    tail_rate > 0.0 || throw(ArgumentError("line weighting exponential-tail rate must be positive"))
    scaled_rate = tail_rate * timestep_s
    decay = exp(-scaled_rate)
    integral_ratio = abs(scaled_rate) <= sqrt(eps(Float64)) ?
        -1.0 + scaled_rate / 2.0 - scaled_rate^2 / 6.0 :
        expm1(-scaled_rate) / scaled_rate
    amplitude_scale = tail_amplitude / tail_rate
    return (
        amplitude_scale * (1.0 + integral_ratio),
        -amplitude_scale * (integral_ratio + decay),
        decay,
    )
end

function _sampled_weighting_panels(
    samples::LineWeightingSamples,
    timestep_s::Float64,
    cutoff_index::Int,
    cutoff_value::Float64,
    cutoff_direction::Symbol,
    scale::Union{Nothing,Float64},
    eta::Union{Nothing,Float64},
)
    times = samples.time_s
    amplitudes = samples.amplitude
    point_count = length(times)
    1 <= cutoff_index <= point_count ||
        throw(ArgumentError("line weighting cutoff index is outside the sample table"))

    weights = [0.0]
    index = 1
    previous_time = times[1]
    previous_amplitude = amplitudes[1]
    initial_step = trunc(Int, times[1] / timestep_s)
    panel_end = initial_step * timestep_s
    panel_index = initial_step - 1
    panel_area = 0.0
    panel_moment = 0.0
    additional_integral = 0.0
    tail_mode = false
    tail_start = 0.0
    tail_amplitude = 0.0
    raw_tail = (0.0, 0.0, 0.0)
    local_scale = scale
    local_eta = eta

    # The first pass through the original state machine allocates the first
    # coefficient and advances to the first timestep boundary.
    panel_end += timestep_s
    panel_index += 1
    index += 1

    while true
        current_time = times[index]
        current_amplitude = amplitudes[index]
        boundary_hit = false
        endpoint_index = index

        if !tail_mode && current_time > panel_end
            slope = (current_amplitude - previous_amplitude) /
                    (current_time - previous_time)
            current_amplitude = previous_amplitude + slope * (panel_end - previous_time)
            current_time = panel_end
            endpoint_index = index - 1
            boundary_hit = true
        elseif index == point_count
            boundary_hit = true
            endpoint_index = -1
        end

        half_width = 0.5 * (current_time - previous_time)
        segment_area = half_width * (previous_amplitude + current_amplitude)
        panel_area += segment_area
        panel_moment += half_width *
                        (previous_time * previous_amplitude +
                         current_time * current_amplitude)
        previous_time = current_time
        previous_amplitude = current_amplitude

        if tail_mode
            endpoint_index < 0 && begin
                panel_area != 0.0 ||
                    throw(ArgumentError("line weighting exponential tail has zero integral"))
                mean_time = panel_moment / panel_area
                mean_time > tail_start ||
                    throw(ArgumentError("line weighting exponential tail has nonpositive decay span"))
                tail_rate = inv(mean_time - tail_start)
                raw_tail = (tail_rate * panel_area * local_scale, tail_rate, tail_start)
                break
            end
            index += 1
            continue
        end

        additional_integral += segment_area
        if !boundary_hit
            index += 1
            continue
        end

        panel_moment /= timestep_s
        if local_scale === nothing
            zero_area = panel_area
            correction = zero_area - panel_moment
            upper = 1.0 + correction
            lower = 1.0 - correction
            upper != 0.0 && lower != 0.0 ||
                throw(ArgumentError("line weighting first-panel normalization is singular"))
            local_eta = upper / lower
            local_scale = inv(upper)
            weights[end] = panel_moment * local_scale
        else
            leading = panel_moment - panel_index * panel_area
            weights[end] += (panel_area - leading) * local_scale
            push!(weights, leading * local_scale)
        end

        cutoff_reached = if cutoff_direction == :rising
            endpoint_index > cutoff_index && current_amplitude >= cutoff_value
        elseif cutoff_direction == :falling
            endpoint_index > cutoff_index && current_amplitude <= cutoff_value
        else
            throw(ArgumentError("line weighting cutoff direction must be :rising or :falling"))
        end
        if cutoff_reached
            tail_mode = true
            tail_start = current_time
            tail_amplitude = current_amplitude
        end

        endpoint_index < 0 && break
        panel_area = 0.0
        panel_moment = 0.0
        panel_end += timestep_s
        panel_index += 1
        if boundary_hit
            # Revisit the raw endpoint after the interpolated boundary.
            continue
        end
        index += 1
    end

    if raw_tail[3] != 0.0 && tail_amplitude == 0.0
        throw(ArgumentError("line weighting tail amplitude must be nonzero"))
    end
    return (
        weights = weights,
        raw_tail = raw_tail,
        scale = something(local_scale),
        eta = something(local_eta),
        additional_integral = additional_integral,
        tail_amplitude = tail_amplitude,
    )
end

function _resistance_adjusted_line_tail(
    raw_tail::NTuple{3,Float64},
    additional_integral::Float64,
    cutoff_amplitude::Float64,
    target_integral::Float64,
    scale::Float64,
    maximum_iterations::Int,
    relative_tolerance::Float64,
)
    _, tail_rate, tail_start = raw_tail
    tail_start > 0.0 && tail_rate > 0.0 ||
        throw(ArgumentError("resistance-adjusted line weighting requires an exponential tail"))
    cutoff_amplitude > 0.0 ||
        throw(ArgumentError("resistance-adjusted line weighting cutoff amplitude must be positive"))
    correction_integral = target_integral - additional_integral
    correction_integral != 0.0 || return (0.0, tail_rate, tail_start)
    log_ratio = log(abs(correction_integral / cutoff_amplitude))
    correction = 0.0
    for iteration in 1:maximum_iterations
        candidate_rate = tail_rate + correction
        if iteration > 1 && abs(correction) / tail_rate <= relative_tolerance
            return (
                correction_integral * tail_rate * scale,
                tail_rate,
                tail_start,
            )
        end
        tail_rate = candidate_rate
        tail_rate > 0.0 ||
            throw(ArgumentError("resistance-adjusted line weighting produced a nonpositive tail rate"))
        residual = -log(tail_rate) - tail_rate * tail_start - log_ratio
        derivative = inv(tail_rate) + tail_start
        correction = residual / derivative
        for _ in 1:6
            abs(correction) < tail_rate && break
            correction /= 5.0
        end
        abs(correction) < tail_rate ||
            throw(ArgumentError("resistance-adjusted line weighting Newton step left its positive domain"))
    end
    throw(ArgumentError(
        "resistance-adjusted line weighting did not converge in $maximum_iterations iterations",
    ))
end

function sampled_line_weighting_coefficients(
    propagation::LineWeightingSamples,
    admittance::LineWeightingSamples,
    timestep_s::Real,
    characteristic_impedance_ohm::Real;
    propagation_peak_index::Integer,
    admittance_rise_index::Integer,
    propagation_cutoff_fraction::Real = 0.01,
    admittance_cutoff_fraction::Real = 0.1,
    total_resistance_ohm::Real = 0.0,
    maximum_tail_iterations::Integer = 100,
    tail_relative_tolerance::Real = 1.0e-8,
)
    dt = Float64(timestep_s)
    zinf = Float64(characteristic_impedance_ohm)
    resistance = Float64(total_resistance_ohm)
    propagation_cutoff = Float64(propagation_cutoff_fraction)
    admittance_cutoff = Float64(admittance_cutoff_fraction)
    relative_tolerance = Float64(tail_relative_tolerance)
    isfinite(dt) && dt > 0.0 || throw(ArgumentError("line weighting timestep must be positive"))
    isfinite(zinf) && zinf > 0.0 ||
        throw(ArgumentError("line weighting characteristic impedance must be positive"))
    isfinite(resistance) && resistance >= 0.0 ||
        throw(ArgumentError("line weighting total resistance must be finite and nonnegative"))
    propagation_cutoff > 0.0 && admittance_cutoff > 0.0 ||
        throw(ArgumentError("line weighting cutoff fractions must be positive"))
    Int(maximum_tail_iterations) > 0 ||
        throw(ArgumentError("maximum_tail_iterations must be positive"))
    isfinite(relative_tolerance) && relative_tolerance > 0.0 ||
        throw(ArgumentError("tail_relative_tolerance must be positive and finite"))

    combined_integral = _trapezoid_integral(propagation) +
                        _trapezoid_integral(admittance)
    combined_integral > 0.0 ||
        throw(ArgumentError("combined line weighting integral must be positive"))
    normalized = resistance <= 0.0
    normalization = normalized ? inv(combined_integral) : 1.0
    normalized_propagation = line_weighting_samples(
        propagation.time_s,
        propagation.amplitude .* normalization,
    )
    normalized_admittance = line_weighting_samples(
        admittance.time_s,
        admittance.amplitude .* normalization,
    )

    rise_index = Int(admittance_rise_index)
    peak_index = Int(propagation_peak_index)
    admittance_result = _sampled_weighting_panels(
        normalized_admittance,
        dt,
        rise_index,
        admittance_cutoff * normalized_admittance.amplitude[rise_index],
        :rising,
        nothing,
        nothing,
    )
    propagation_result = _sampled_weighting_panels(
        normalized_propagation,
        dt,
        peak_index,
        propagation_cutoff * normalized_propagation.amplitude[peak_index],
        :falling,
        admittance_result.scale,
        admittance_result.eta,
    )

    admittance_raw_tail = admittance_result.raw_tail
    propagation_raw_tail = propagation_result.raw_tail
    if resistance > 0.0
        denominator = 2.0 * zinf + resistance
        admittance_raw_tail = _resistance_adjusted_line_tail(
            admittance_raw_tail,
            admittance_result.additional_integral,
            admittance_result.tail_amplitude,
            resistance / denominator,
            admittance_result.scale,
            Int(maximum_tail_iterations),
            relative_tolerance,
        )
        propagation_raw_tail = _resistance_adjusted_line_tail(
            propagation_raw_tail,
            propagation_result.additional_integral,
            propagation_result.tail_amplitude,
            2.0 * zinf / denominator,
            propagation_result.scale,
            Int(maximum_tail_iterations),
            relative_tolerance,
        )
    end
    skip_steps = trunc(Int, propagation.time_s[1] / dt)
    history_count = max(
        length(admittance_result.weights),
        length(propagation_result.weights) + max(skip_steps - 1, 0),
    )
    discrete_admittance_tail =
        _discrete_tail_coefficients(admittance_raw_tail, dt)
    discrete_propagation_tail =
        _discrete_tail_coefficients(propagation_raw_tail, dt)
    passive_tolerance = 1.0e-12
    all(weight -> weight >= -passive_tolerance, admittance_result.weights) &&
        all(weight -> weight >= -passive_tolerance, propagation_result.weights) &&
        all(value -> value >= -passive_tolerance, discrete_admittance_tail[1:2]) &&
        all(value -> value >= -passive_tolerance, discrete_propagation_tail[1:2]) ||
        throw(ArgumentError("sampled line weighting response contains an active coefficient"))
    coefficient_sum = sum(admittance_result.weights) +
                      sum(propagation_result.weights)
    tail_dc_gain(tail) = tail[3] == 1.0 ? Inf :
        (tail[1] + tail[2]) / (1.0 - tail[3])
    dc_gain = coefficient_sum +
              tail_dc_gain(discrete_admittance_tail) +
              tail_dc_gain(discrete_propagation_tail)
    isfinite(dc_gain) && dc_gain >= -1.0e-12 && dc_gain <= 1.0 + 1.0e-10 ||
        throw(ArgumentError("sampled line weighting response is not passive at zero frequency"))
    return SampledLineWeightingCoefficients(
        dt,
        zinf,
        admittance_result.eta,
        skip_steps,
        history_count,
        admittance_result.weights,
        propagation_result.weights,
        discrete_admittance_tail,
        discrete_propagation_tail,
        combined_integral,
        normalized,
        resistance,
        relative_tolerance,
        coefficient_sum,
        dc_gain,
    )
end

mutable struct SampledFrequencyDependentLine <: EMTElement
    a::Int
    b::Int
    coefficients::SampledLineWeightingCoefficients
    loss_factor::Float64
    from_wave_history::Vector{Float64}
    to_wave_history::Vector{Float64}
    write_index::Int
    admittance_tail_from::Float64
    admittance_tail_to::Float64
    propagation_tail_from::Float64
    propagation_tail_to::Float64
    history_current_from::Float64
    history_current_to::Float64
    terminal_voltage_from::Float64
    terminal_voltage_to::Float64
    terminal_current_from::Float64
    terminal_current_to::Float64
    convolution_from::Float64
    convolution_to::Float64
    update_count::Int
end

function sampled_frequency_dependent_line(
    a::Integer,
    b::Integer,
    coefficients::SampledLineWeightingCoefficients;
    loss_factor::Real = 1.0e-5,
)
    from_node = Int(a)
    to_node = Int(b)
    from_node >= 0 && to_node >= 0 ||
        throw(ArgumentError("sampled frequency-dependent line nodes must be nonnegative"))
    loss = Float64(loss_factor)
    isfinite(loss) && 0.0 <= loss < 1.0 ||
        throw(ArgumentError("sampled line loss factor must be finite and in [0, 1)"))
    history_count = max(coefficients.history_sample_count, 2)
    return SampledFrequencyDependentLine(
        from_node,
        to_node,
        coefficients,
        loss,
        zeros(history_count),
        zeros(history_count),
        1,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0,
    )
end

function _sampled_line_history_value(history::AbstractVector{Float64}, newest::Int, lag::Int)
    index = mod1(newest - lag, length(history))
    return history[index]
end

function _sampled_line_weighted_history(
    weights::AbstractVector{Float64},
    history::AbstractVector{Float64},
    newest::Int,
    initial_lag::Int,
)
    total = 0.0
    for (offset, weight) in enumerate(weights)
        total += weight *
                 _sampled_line_history_value(history, newest, initial_lag + offset - 1)
    end
    return total
end

function _sampled_line_tail_update(
    previous_tail::Float64,
    coefficients::NTuple{3,Float64},
    history::AbstractVector{Float64},
    newest::Int,
    initial_lag::Int,
)
    current_gain, previous_gain, decay = coefficients
    return current_gain * _sampled_line_history_value(history, newest, initial_lag) +
           previous_gain * _sampled_line_history_value(history, newest, initial_lag + 1) +
           decay * previous_tail
end

function sampled_line_history_convolution!(line::SampledFrequencyDependentLine)
    coefficients = line.coefficients
    newest = line.write_index
    propagation_lag = max(coefficients.propagation_delay_steps - 1, 0)
    admittance_from = _sampled_line_weighted_history(
        coefficients.admittance_weights,
        line.from_wave_history,
        newest,
        0,
    )
    admittance_to = _sampled_line_weighted_history(
        coefficients.admittance_weights,
        line.to_wave_history,
        newest,
        0,
    )
    propagation_from = _sampled_line_weighted_history(
        coefficients.propagation_weights,
        line.from_wave_history,
        newest,
        propagation_lag,
    )
    propagation_to = _sampled_line_weighted_history(
        coefficients.propagation_weights,
        line.to_wave_history,
        newest,
        propagation_lag,
    )
    line.admittance_tail_from = _sampled_line_tail_update(
        line.admittance_tail_from,
        coefficients.admittance_tail,
        line.from_wave_history,
        newest,
        0,
    )
    line.admittance_tail_to = _sampled_line_tail_update(
        line.admittance_tail_to,
        coefficients.admittance_tail,
        line.to_wave_history,
        newest,
        0,
    )
    line.propagation_tail_from = _sampled_line_tail_update(
        line.propagation_tail_from,
        coefficients.propagation_tail,
        line.from_wave_history,
        newest,
        propagation_lag,
    )
    line.propagation_tail_to = _sampled_line_tail_update(
        line.propagation_tail_to,
        coefficients.propagation_tail,
        line.to_wave_history,
        newest,
        propagation_lag,
    )
    # The first table cross-couples terminal histories; the delayed second
    # table propagates the same-terminal wave. This is the single-mode form
    # of the accepted line-history mutation order.
    line.convolution_from = admittance_to + line.admittance_tail_to +
                            propagation_from + line.propagation_tail_from
    line.convolution_to = admittance_from + line.admittance_tail_from +
                          propagation_to + line.propagation_tail_to
    return (from = line.convolution_from, to = line.convolution_to)
end

function sampled_line_history_convolution(
    coefficients::SampledLineWeightingCoefficients,
    from_wave_history::AbstractVector,
    to_wave_history::AbstractVector,
    newest_index::Integer;
    admittance_tail_from::Real = 0.0,
    admittance_tail_to::Real = 0.0,
    propagation_tail_from::Real = 0.0,
    propagation_tail_to::Real = 0.0,
)
    line = sampled_frequency_dependent_line(1, 2, coefficients)
    length(from_wave_history) == length(line.from_wave_history) ||
        throw(ArgumentError("from-wave history length must match coefficient history count"))
    length(to_wave_history) == length(line.to_wave_history) ||
        throw(ArgumentError("to-wave history length must match coefficient history count"))
    line.from_wave_history .= from_wave_history
    line.to_wave_history .= to_wave_history
    line.write_index = Int(newest_index)
    line.admittance_tail_from = Float64(admittance_tail_from)
    line.admittance_tail_to = Float64(admittance_tail_to)
    line.propagation_tail_from = Float64(propagation_tail_from)
    line.propagation_tail_to = Float64(propagation_tail_to)
    result = sampled_line_history_convolution!(line)
    return merge(result, (
        admittance_tail_from = line.admittance_tail_from,
        admittance_tail_to = line.admittance_tail_to,
        propagation_tail_from = line.propagation_tail_from,
        propagation_tail_to = line.propagation_tail_to,
    ))
end

function line_surge_admittance(line::SampledFrequencyDependentLine)::Float64
    return inv(line.coefficients.characteristic_impedance_ohm * line.coefficients.eta)
end

function _sampled_line_harmonic_weighting(
    weights::AbstractVector{<:Real},
    tail::NTuple{3,Float64},
    initial_lag::Int,
    step_rotation::ComplexF64,
)
    inverse_rotation = inv(step_rotation)
    finite_response = zero(ComplexF64)
    lag_rotation = inverse_rotation^initial_lag
    for weight in weights
        finite_response += Float64(weight) * lag_rotation
        lag_rotation *= inverse_rotation
    end
    current_gain, previous_gain, decay = tail
    tail_denominator = 1.0 - decay * inverse_rotation
    abs(tail_denominator) > eps(Float64) || throw(ArgumentError(
        "sampled-line steady-state tail recurrence is singular",
    ))
    tail_response = inverse_rotation^initial_lag *
        (current_gain + previous_gain * inverse_rotation) /
        tail_denominator
    return finite_response + tail_response
end

"""
    sampled_line_steady_state_terminal_admittance(line, frequency_hz)

Return the exact complex terminal admittance of the sampled line's discrete
history recurrence at `frequency_hz`. The transform follows the production
mutation order: the current solve uses the preceding history source, then the
current wave slot and both weighting recurrences are updated for the next
timestep. Phasors use peak-value cosine convention.
"""
function sampled_line_steady_state_terminal_admittance(
    line::SampledFrequencyDependentLine,
    frequency_hz::Real,
)
    frequency = Float64(frequency_hz)
    isfinite(frequency) && frequency >= 0.0 || throw(ArgumentError(
        "sampled-line steady-state frequency must be finite and nonnegative",
    ))
    coefficients = line.coefficients
    step_rotation = cis(2.0 * pi * frequency * coefficients.timestep_s)
    admittance_response = _sampled_line_harmonic_weighting(
        coefficients.admittance_weights,
        coefficients.admittance_tail,
        0,
        step_rotation,
    )
    propagation_lag = max(coefficients.propagation_delay_steps - 1, 0)
    propagation_response = _sampled_line_harmonic_weighting(
        coefficients.propagation_weights,
        coefficients.propagation_tail,
        propagation_lag,
        step_rotation,
    )
    ring_rotation = step_rotation^(-length(line.from_wave_history))
    delayed_loss = line.loss_factor * ring_rotation
    wave_denominator = 1.0 - delayed_loss^2
    abs(wave_denominator) > eps(Float64) || throw(ArgumentError(
        "sampled-line steady-state wave recurrence is singular",
    ))
    history_rotation = inv(step_rotation)
    conductance = line_surge_admittance(line)
    response_scale = conductance *
        coefficients.characteristic_impedance_ohm *
        history_rotation / wave_denominator
    self_admittance = conductance - response_scale *
        (admittance_response - delayed_loss * propagation_response)
    mutual_admittance = -response_scale *
        (propagation_response - delayed_loss * admittance_response)
    return ComplexF64[
        self_admittance mutual_admittance
        mutual_admittance self_admittance
    ]
end

line_terminal_voltages(line::SampledFrequencyDependentLine) =
    (from = line.terminal_voltage_from, to = line.terminal_voltage_to)
line_terminal_currents(line::SampledFrequencyDependentLine) =
    (from = line.terminal_current_from, to = line.terminal_current_to)
line_history_currents(line::SampledFrequencyDependentLine) =
    (from = line.history_current_from, to = line.history_current_to)

function stamp!(
    y::AbstractMatrix{Float64},
    rhs::AbstractVector{Float64},
    line::SampledFrequencyDependentLine,
    _time_s::Float64,
    _timestep_s::Float64,
)
    conductance = line_surge_admittance(line)
    stamp_conductance!(y, line.a, 0, conductance)
    stamp_conductance!(y, line.b, 0, conductance)
    stamp_history_current!(rhs, line.a, 0, line.history_current_from)
    stamp_history_current!(rhs, line.b, 0, line.history_current_to)
    return nothing
end

function update!(
    line::SampledFrequencyDependentLine,
    voltages::Union{AbstractVector{Float64},NTuple{2,Float64}},
    timestep_s::Float64,
)
    abs(timestep_s - line.coefficients.timestep_s) <=
        64.0 * eps(Float64) * max(timestep_s, line.coefficients.timestep_s) ||
        throw(ArgumentError("sampled line timestep does not match its weighting coefficients"))
    conductance = line_surge_admittance(line)
    from_voltage = line.a == 0 ? 0.0 : voltages[line.a]
    to_voltage = line.b == 0 ? 0.0 : voltages[line.b]
    from_current = conductance * from_voltage + line.history_current_from
    to_current = conductance * to_voltage + line.history_current_to
    previous_from_wave = line.from_wave_history[line.write_index]
    previous_to_wave = line.to_wave_history[line.write_index]
    impedance = line.coefficients.characteristic_impedance_ohm
    # Each endpoint history stores the wave arriving from the opposite
    # terminal.  OVER16 overwrites the current ring slot before applying the
    # two weighting tables, so voltage and attenuation both cross terminals.
    line.from_wave_history[line.write_index] =
        impedance * to_voltage - line.loss_factor * previous_to_wave
    line.to_wave_history[line.write_index] =
        impedance * from_voltage - line.loss_factor * previous_from_wave
    convolution = sampled_line_history_convolution!(line)
    line.history_current_from = -conductance * convolution.from
    line.history_current_to = -conductance * convolution.to
    line.terminal_voltage_from = from_voltage
    line.terminal_voltage_to = to_voltage
    line.terminal_current_from = from_current
    line.terminal_current_to = to_current
    line.write_index = line.write_index == length(line.from_wave_history) ?
        1 : line.write_index + 1
    line.update_count += 1
    return nothing
end

trace_output_channel_count(::SampledFrequencyDependentLine) = 4

function trace_output_channel_names!(
    names::Vector{Symbol},
    element_name::Symbol,
    ::SampledFrequencyDependentLine,
)
    append!(names, (
        Symbol(element_name, :_from_current_a),
        Symbol(element_name, :_to_current_a),
        Symbol(element_name, :_from_history_a),
        Symbol(element_name, :_to_history_a),
    ))
    return names
end

function trace_output_values!(
    output::AbstractMatrix{Float64},
    first_channel::Int,
    sample::Int,
    line::SampledFrequencyDependentLine,
    _voltage::AbstractVector{Float64},
)
    output[first_channel, sample] = line.terminal_current_from
    output[first_channel + 1, sample] = line.terminal_current_to
    output[first_channel + 2, sample] = line.history_current_from
    output[first_channel + 3, sample] = line.history_current_to
    return first_channel + 4
end

mutable struct SampledFrequencyDependentLineGroup <: EMTElement
    from_nodes::Vector{Int}
    to_nodes::Vector{Int}
    modal_lines::Vector{SampledFrequencyDependentLine}
    phase_to_modal::Matrix{Float64}
    modal_to_phase::Matrix{Float64}
    phase_admittance::Matrix{Float64}
    history_current_from::Vector{Float64}
    history_current_to::Vector{Float64}
    terminal_voltage_from::Vector{Float64}
    terminal_voltage_to::Vector{Float64}
    terminal_current_from::Vector{Float64}
    terminal_current_to::Vector{Float64}
    modal_voltage_from::Vector{Float64}
    modal_voltage_to::Vector{Float64}
    modal_current_from::Vector{Float64}
    modal_current_to::Vector{Float64}
    update_count::Int
end

function sampled_frequency_dependent_line_group(
    from_nodes::AbstractVector{<:Integer},
    to_nodes::AbstractVector{<:Integer},
    coefficients::AbstractVector{SampledLineWeightingCoefficients};
    loss_factors::AbstractVector{<:Real}=fill(1.0e-5, length(coefficients)),
    modal_to_phase::AbstractMatrix{<:Real}=_transposed_line_modal_transform(),
    phase_to_modal::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
)
    phase_count = length(coefficients)
    phase_count >= 2 ||
        throw(ArgumentError("sampled frequency-dependent line group requires at least two modes"))
    length(from_nodes) == phase_count && length(to_nodes) == phase_count ||
        throw(ArgumentError("sampled line terminal count must match its modal coefficient count"))
    length(loss_factors) == phase_count ||
        throw(ArgumentError("sampled line loss-factor count must match its modal coefficient count"))
    from_indices = Int.(from_nodes)
    to_indices = Int.(to_nodes)
    all(>=(0), from_indices) && all(>=(0), to_indices) ||
        throw(ArgumentError("sampled line terminal indices must be nonnegative"))
    length(unique(from_indices)) == phase_count && length(unique(to_indices)) == phase_count ||
        throw(ArgumentError("sampled line phase terminals must be distinct"))
    transform = Matrix{Float64}(modal_to_phase)
    size(transform) == (phase_count, phase_count) ||
        throw(ArgumentError("modal-to-phase transform size must match the sampled line group"))
    all(isfinite, transform) ||
        throw(ArgumentError("modal-to-phase transform entries must be finite"))
    inverse_transform = phase_to_modal === nothing ?
        Matrix{Float64}(inv(transform)) : Matrix{Float64}(phase_to_modal)
    size(inverse_transform) == (phase_count, phase_count) ||
        throw(ArgumentError("phase-to-modal transform size must match the sampled line group"))
    all(isfinite, inverse_transform) ||
        throw(ArgumentError("phase-to-modal transform entries must be finite"))
    maximum(abs.(inverse_transform * transform - I(phase_count))) <= 1.0e-10 ||
        throw(ArgumentError("sampled line modal transforms must be inverses"))
    maximum(abs.(transpose(transform) * transform - I(phase_count))) <= 1.0e-10 ||
        throw(ArgumentError("sampled line modal transform must preserve phase power"))

    modal_lines = SampledFrequencyDependentLine[
        sampled_frequency_dependent_line(
            1,
            2,
            coefficient;
            loss_factor = loss_factor,
        )
        for (coefficient, loss_factor) in zip(coefficients, loss_factors)
    ]
    modal_admittance = Float64[line_surge_admittance(line) for line in modal_lines]
    phase_admittance = transform * Diagonal(modal_admittance) * inverse_transform
    symmetry_error = maximum(abs.(phase_admittance - transpose(phase_admittance)))
    symmetry_error <= 1.0e-10 * max(maximum(abs.(phase_admittance)), 1.0) ||
        throw(ArgumentError("sampled line phase admittance must be reciprocal"))
    minimum(eigvals(Symmetric(phase_admittance))) >= -1.0e-12 ||
        throw(ArgumentError("sampled line phase admittance must be passive"))
    zeros_phase = zeros(phase_count)
    return SampledFrequencyDependentLineGroup(
        from_indices,
        to_indices,
        modal_lines,
        inverse_transform,
        transform,
        Matrix{Float64}(phase_admittance),
        copy(zeros_phase),
        copy(zeros_phase),
        copy(zeros_phase),
        copy(zeros_phase),
        copy(zeros_phase),
        copy(zeros_phase),
        copy(zeros_phase),
        copy(zeros_phase),
        copy(zeros_phase),
        copy(zeros_phase),
        0,
    )
end

function sampled_line_steady_state_terminal_admittance(
    line::SampledFrequencyDependentLineGroup,
    frequency_hz::Real,
)
    modal_terminal_admittances = [
        sampled_line_steady_state_terminal_admittance(modal_line, frequency_hz)
        for modal_line in line.modal_lines
    ]
    modal_self = Diagonal([
        admittance[1, 1] for admittance in modal_terminal_admittances
    ])
    modal_mutual = Diagonal([
        admittance[1, 2] for admittance in modal_terminal_admittances
    ])
    phase_self = line.modal_to_phase * modal_self * line.phase_to_modal
    phase_mutual = line.modal_to_phase * modal_mutual * line.phase_to_modal
    return ComplexF64[
        phase_self phase_mutual
        phase_mutual phase_self
    ]
end

line_terminal_voltages(line::SampledFrequencyDependentLineGroup) =
    (from = copy(line.terminal_voltage_from), to = copy(line.terminal_voltage_to))
line_terminal_currents(line::SampledFrequencyDependentLineGroup) =
    (from = copy(line.terminal_current_from), to = copy(line.terminal_current_to))
line_history_currents(line::SampledFrequencyDependentLineGroup) =
    (from = copy(line.history_current_from), to = copy(line.history_current_to))

function stamp!(
    y::AbstractMatrix{Float64},
    rhs::AbstractVector{Float64},
    line::SampledFrequencyDependentLineGroup,
    _time_s::Float64,
    _timestep_s::Float64,
)
    for column in eachindex(line.from_nodes)
        for row in eachindex(line.from_nodes)
            value = line.phase_admittance[row, column]
            stamp_admittance_entry!(
                y,
                line.from_nodes[row],
                0,
                line.from_nodes[column],
                0,
                value,
            )
            stamp_admittance_entry!(
                y,
                line.to_nodes[row],
                0,
                line.to_nodes[column],
                0,
                value,
            )
        end
    end
    for phase in eachindex(line.from_nodes)
        stamp_history_current!(rhs, line.from_nodes[phase], 0, line.history_current_from[phase])
        stamp_history_current!(rhs, line.to_nodes[phase], 0, line.history_current_to[phase])
    end
    return nothing
end

function update!(
    line::SampledFrequencyDependentLineGroup,
    voltages::AbstractVector{Float64},
    timestep_s::Float64,
)
    for phase in eachindex(line.from_nodes)
        from_node = line.from_nodes[phase]
        to_node = line.to_nodes[phase]
        line.terminal_voltage_from[phase] = from_node == 0 ? 0.0 : voltages[from_node]
        line.terminal_voltage_to[phase] = to_node == 0 ? 0.0 : voltages[to_node]
    end
    mul!(line.modal_voltage_from, line.phase_to_modal, line.terminal_voltage_from)
    mul!(line.modal_voltage_to, line.phase_to_modal, line.terminal_voltage_to)
    for mode in eachindex(line.modal_lines)
        modal_line = line.modal_lines[mode]
        update!(
            modal_line,
            (line.modal_voltage_from[mode], line.modal_voltage_to[mode]),
            timestep_s,
        )
        line.modal_current_from[mode] = modal_line.terminal_current_from
        line.modal_current_to[mode] = modal_line.terminal_current_to
    end
    mul!(line.terminal_current_from, line.modal_to_phase, line.modal_current_from)
    mul!(line.terminal_current_to, line.modal_to_phase, line.modal_current_to)
    for mode in eachindex(line.modal_lines)
        line.modal_current_from[mode] = line.modal_lines[mode].history_current_from
        line.modal_current_to[mode] = line.modal_lines[mode].history_current_to
    end
    mul!(line.history_current_from, line.modal_to_phase, line.modal_current_from)
    mul!(line.history_current_to, line.modal_to_phase, line.modal_current_to)
    line.update_count += 1
    return nothing
end

trace_output_channel_count(line::SampledFrequencyDependentLineGroup) =
    4 * length(line.from_nodes)

function trace_output_channel_names!(
    names::Vector{Symbol},
    element_name::Symbol,
    line::SampledFrequencyDependentLineGroup,
)
    for phase in eachindex(line.from_nodes)
        append!(names, (
            Symbol(element_name, :_phase_, phase, :_from_current_a),
            Symbol(element_name, :_phase_, phase, :_to_current_a),
            Symbol(element_name, :_phase_, phase, :_from_history_a),
            Symbol(element_name, :_phase_, phase, :_to_history_a),
        ))
    end
    return names
end

function trace_output_values!(
    output::AbstractMatrix{Float64},
    first_channel::Int,
    sample::Int,
    line::SampledFrequencyDependentLineGroup,
    _voltage::AbstractVector{Float64},
)
    channel = first_channel
    for phase in eachindex(line.from_nodes)
        output[channel, sample] = line.terminal_current_from[phase]
        output[channel + 1, sample] = line.terminal_current_to[phase]
        output[channel + 2, sample] = line.history_current_from[phase]
        output[channel + 3, sample] = line.history_current_to[phase]
        channel += 4
    end
    return channel
end

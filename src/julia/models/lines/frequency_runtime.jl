
function frequency_dependent_line_sample_fit(
    sample_rows::AbstractVector,
    target_frequency_hz::Real,
    line_length::Real,
)
    target_frequency = _checked_line_target_frequency(target_frequency_hz)
    line_length_value = _checked_line_length(line_length)
    frequencies, rows, mode_count = _checked_line_frequency_sample_rows(sample_rows)
    return _line_frequency_sample_fit_from_sorted(
        frequencies,
        rows,
        mode_count,
        target_frequency,
        line_length_value,
    )
end

function FrequencyDependentLineModalState(
    transform::LineModalTransform,
    frequency_points::AbstractVector{LineFrequencyPoint},
)
    mode_count = _line_modal_dimension(transform)
    points = _checked_line_frequency_points(frequency_points, mode_count)
    return FrequencyDependentLineModalState(
        transform,
        points,
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        0,
    )
end

FrequencyDependentLineModalState(
    transform::LineModalTransform,
    fit::LineFrequencySampleFitResult,
) = FrequencyDependentLineModalState(transform, fit.frequency_points)

function FrequencyDependentLineRuntimeState(
    transform::LineModalTransform,
    frequency_points::AbstractVector{LineFrequencyPoint},
)
    modal_state = FrequencyDependentLineModalState(transform, frequency_points)
    mode_count = length(modal_state.frequency_points)
    return FrequencyDependentLineRuntimeState(
        modal_state,
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        0,
    )
end

FrequencyDependentLineRuntimeState(
    transform::LineModalTransform,
    fit::LineFrequencySampleFitResult,
) = FrequencyDependentLineRuntimeState(transform, fit.frequency_points)

function _checked_line_modal_yz_sample_matrices(
    yz_sample_matrices::AbstractVector,
    order::AbstractVector{Int},
    mode_count::Int,
)
    length(yz_sample_matrices) == length(order) ||
        throw(ArgumentError("line modal YZ sample matrix count must match the frequency sample row count"))
    matrices = Vector{Matrix{ComplexF64}}(undef, length(order))
    for (target_index, source_index) in pairs(order)
        matrix = _checked_line_complex_square_matrix(
            yz_sample_matrices[source_index],
            "line modal YZ sample matrix",
        )
        size(matrix, 1) == mode_count ||
            throw(ArgumentError("line modal YZ sample matrix dimension must match the frequency sample mode count"))
        matrices[target_index] = matrix
    end
    return matrices
end

function FrequencyDependentLineSampleRuntimeState(
    transform::LineModalTransform,
    sample_rows::AbstractVector,
    initial_frequency_hz::Real,
    line_length::Real,
)
    initial_frequency = _checked_line_target_frequency(initial_frequency_hz)
    line_length_value = _checked_line_length(line_length)
    frequencies, rows, mode_count = _checked_line_frequency_sample_rows(sample_rows)
    fit = _line_frequency_sample_fit_from_sorted(
        frequencies,
        rows,
        mode_count,
        initial_frequency,
        line_length_value,
    )
    runtime_state = FrequencyDependentLineRuntimeState(transform, fit)
    return FrequencyDependentLineSampleRuntimeState(
        collect(frequencies),
        [copy(row) for row in rows],
        line_length_value,
        fit,
        runtime_state,
        fit.frequency_hz,
        fit.frequency_hz,
        fit.interpolation_weight,
        fit.interpolation_weight,
        0,
    )
end

function _line_modal_sample_matrix_from_sorted(
    frequencies::AbstractVector{Float64},
    matrices::AbstractVector{<:AbstractMatrix{ComplexF64}},
    target_frequency::Float64,
)
    exact_index = findfirst(
        frequency -> abs(frequency - target_frequency) <= _line_frequency_row_tolerance(frequency),
        frequencies,
    )
    exact_index !== nothing && return copy(matrices[exact_index])
    first(frequencies) < target_frequency < last(frequencies) ||
        throw(ArgumentError("target_frequency_hz must be within the sampled frequency range"))
    upper_index = findfirst(frequency -> frequency > target_frequency, frequencies)
    lower_index = upper_index - 1
    weight = _line_frequency_interpolation_weight(
        frequencies[lower_index],
        frequencies[upper_index],
        target_frequency,
    )
    return matrices[lower_index] .+
        weight .* (matrices[upper_index] .- matrices[lower_index])
end

function FrequencyDependentLineModalSampleRuntimeState(
    yz_sample_matrices::AbstractVector,
    sample_rows::AbstractVector,
    initial_frequency_hz::Real,
    line_length::Real;
    ntol::Integer = 1,
    nrp::Integer = 0,
)
    initial_frequency = _checked_line_target_frequency(initial_frequency_hz)
    line_length_value = _checked_line_length(line_length)
    frequencies, rows, mode_count, order = _checked_line_frequency_sample_rows_with_order(sample_rows)
    matrices = _checked_line_modal_yz_sample_matrices(yz_sample_matrices, order, mode_count)
    mode_order_state = LineModeUnwindState(mode_count)
    initial_matrix = _line_modal_sample_matrix_from_sorted(frequencies, matrices, initial_frequency)
    solution = line_modal_solution(
        initial_matrix,
        initial_frequency;
        unwind_state = mode_order_state,
        ntol = ntol,
        nrp = nrp,
    )
    sample_runtime_state = FrequencyDependentLineSampleRuntimeState(
        solution.transform,
        rows,
        initial_frequency,
        line_length_value,
    )
    return FrequencyDependentLineModalSampleRuntimeState(
        collect(frequencies),
        [copy(row) for row in rows],
        matrices,
        line_length_value,
        mode_order_state,
        solution,
        sample_runtime_state,
        initial_frequency,
        initial_frequency,
        0,
    )
end

function _checked_line_recursive_convolution_matrices(
    pole_decay::AbstractMatrix,
    residue::AbstractMatrix,
    mode_count::Int,
)
    decay = _checked_line_complex_rectangular_matrix(
        pole_decay,
        mode_count,
        "line recursive convolution pole decay",
    )
    terms = _checked_line_complex_rectangular_matrix(
        residue,
        mode_count,
        "line recursive convolution residue",
    )
    size(decay) == size(terms) ||
        throw(ArgumentError("line recursive convolution decay/residue matrices must have identical dimensions"))
    all(value -> abs(value) <= 1.0 + 64.0 * eps(Float64), decay) ||
        throw(ArgumentError("line recursive convolution pole decays must be stable"))
    return decay, terms
end

function _checked_line_recursive_fit_samples(
    sample_frequencies_hz::AbstractVector,
    modal_response_samples::AbstractMatrix,
)
    isempty(sample_frequencies_hz) &&
        throw(ArgumentError("line recursive convolution fit samples must be nonempty"))
    frequencies = Vector{Float64}(undef, length(sample_frequencies_hz))
    for idx in eachindex(sample_frequencies_hz)
        frequency = Float64(sample_frequencies_hz[idx])
        isfinite(frequency) && frequency >= 0.0 ||
            throw(ArgumentError("line recursive convolution fit frequencies must be finite and nonnegative"))
        frequencies[idx] = frequency
    end
    rows, cols = size(modal_response_samples)
    rows > 0 && cols == length(frequencies) ||
        throw(ArgumentError("line recursive convolution modal response samples must be mode x frequency"))
    samples = Matrix{ComplexF64}(undef, rows, cols)
    for mode in 1:rows, sample in 1:cols
        value = ComplexF64(modal_response_samples[mode, sample])
        isfinite(real(value)) && isfinite(imag(value)) ||
            throw(ArgumentError("line recursive convolution modal response samples must be finite"))
        samples[mode, sample] = value
    end
    order = sortperm(frequencies)
    sorted_frequencies = frequencies[order]
    for idx in 2:length(sorted_frequencies)
        abs(sorted_frequencies[idx] - sorted_frequencies[idx - 1]) >
            _line_frequency_row_tolerance(sorted_frequencies[idx]) ||
            throw(ArgumentError("line recursive convolution fit frequencies must not contain duplicates"))
    end
    return sorted_frequencies, samples[:, order]
end

function _line_recursive_convolution_basis(
    frequency_hz::Float64,
    dt_s::Float64,
    pole_decay::ComplexF64,
)
    phase = cis(-2.0 * pi * frequency_hz * dt_s)
    denominator = 1.0 + 0.0im - pole_decay * phase
    abs(denominator) > 0.0 ||
        throw(ArgumentError("line recursive convolution fit basis is singular"))
    return inv(denominator)
end

function _checked_line_complex(value, label::AbstractString)
    complex_value = ComplexF64(value)
    isfinite(real(complex_value)) && isfinite(imag(complex_value)) ||
        throw(ArgumentError("$label must be finite"))
    return complex_value
end

function _checked_line_positive_finite(value::Real, label::AbstractString)
    finite_value = Float64(value)
    isfinite(finite_value) && finite_value > 0.0 ||
        throw(ArgumentError("$label must be finite and positive"))
    return finite_value
end

function _semlyen_line_convolution_general_coefficients(
    pole::ComplexF64,
    residue::ComplexF64,
    dt::Float64,
)
    scaled_pole = pole * dt
    pole_norm2 = abs2(scaled_pole)
    pole_norm2 > 0.0 ||
        throw(ArgumentError("Semlyen exponential convolution pole must be nonzero"))
    decay = exp(-real(scaled_pole)) * ComplexF64(cos(-imag(scaled_pole)), sin(-imag(scaled_pole)))
    real_ratio = ((1.0 - real(decay)) * real(scaled_pole) - imag(decay) * imag(scaled_pole)) / pole_norm2
    imag_ratio = (-(1.0 - real(decay)) * imag(scaled_pole) - imag(decay) * real(scaled_pole)) / pole_norm2
    ratio = ComplexF64(real_ratio, imag_ratio)
    current_gain = residue * (1.0 - ratio)
    delayed_gain = residue * (ratio - decay)
    return decay, current_gain, delayed_gain
end

function _semlyen_line_convolution_taylor_coefficients(
    pole::ComplexF64,
    residue::ComplexF64,
    dt::Float64,
    convergence_tolerance::Float64,
)
    scaled_real = real(pole) * dt
    scaled_imag = imag(pole) * dt
    tolerance2 = convergence_tolerance * convergence_tolerance
    iteration = 0
    denominator = 1
    numerator_real = -1.0
    numerator_imag = 0.0
    first_sum_real = 0.0
    first_sum_imag = 0.0
    second_sum_real = 0.0
    second_sum_imag = 0.0
    while true
        iteration += 1
        iteration < 15 || return nothing
        previous_denominator = denominator
        denominator *= iteration + 1
        next_numerator_real = -numerator_real * scaled_real + numerator_imag * scaled_imag
        numerator_imag = -numerator_real * scaled_imag - numerator_imag * scaled_real
        numerator_real = next_numerator_real
        first_real = numerator_real / denominator
        first_imag = numerator_imag / denominator
        first_sum_real += first_real
        first_sum_imag += first_imag
        first_error = (first_real * first_real + first_imag * first_imag) /
            (first_sum_real * first_sum_real + first_sum_imag * first_sum_imag)
        factor = (denominator - previous_denominator) /
            (denominator * previous_denominator)
        second_real = numerator_real * factor
        second_imag = numerator_imag * factor
        second_sum_real += second_real
        second_sum_imag += second_imag
        second_error = (second_real * second_real + second_imag * second_imag) /
            (second_sum_real * second_sum_real + second_sum_imag * second_sum_imag)
        (second_error <= tolerance2 && first_error <= tolerance2) && break
    end
    decay = exp(-scaled_real) * ComplexF64(cos(scaled_imag), -sin(scaled_imag))
    current_gain = residue * ComplexF64(first_sum_real, first_sum_imag)
    delayed_gain = residue * ComplexF64(second_sum_real, second_sum_imag)
    return decay, current_gain, delayed_gain, iteration
end

function semlyen_line_exponential_convolution_coefficients(
    pole,
    residue,
    dt_s::Real;
    convergence_tolerance::Real = 1.0e-12,
)
    pole_value = _checked_line_complex(pole, "Semlyen exponential convolution pole")
    residue_value = _checked_line_complex(residue, "Semlyen exponential convolution residue")
    dt = _checked_line_positive_finite(dt_s, "Semlyen exponential convolution dt_s")
    tolerance = _checked_line_positive_finite(
        convergence_tolerance,
        "Semlyen exponential convolution convergence tolerance",
    )
    scaled_norm2 = abs2(pole_value * dt)
    if scaled_norm2 <= 1.0e-4
        taylor = _semlyen_line_convolution_taylor_coefficients(
            pole_value,
            residue_value,
            dt,
            tolerance,
        )
        if taylor !== nothing
            decay, current_gain, delayed_gain, iteration = taylor
            return SemlyenLineExponentialConvolutionCoefficients(
                pole_value,
                residue_value,
                dt,
                decay,
                current_gain,
                delayed_gain,
                true,
                iteration,
            )
        end
    end
    decay, current_gain, delayed_gain = _semlyen_line_convolution_general_coefficients(
        pole_value,
        residue_value,
        dt,
    )
    return SemlyenLineExponentialConvolutionCoefficients(
        pole_value,
        residue_value,
        dt,
        decay,
        current_gain,
        delayed_gain,
        false,
        0,
    )
end

function semlyen_line_harmonic_history_update(
    pole,
    residue,
    angular_frequency_rad_s::Real,
    from_voltage_phasor,
    to_voltage_phasor;
    previous_from_history = 0.0 + 0.0im,
    previous_to_history = 0.0 + 0.0im,
)
    pole_value = _checked_line_complex(pole, "Semlyen harmonic history pole")
    residue_value = _checked_line_complex(residue, "Semlyen harmonic history residue")
    omega = Float64(angular_frequency_rad_s)
    isfinite(omega) && omega >= 0.0 ||
        throw(ArgumentError("Semlyen harmonic history angular frequency must be finite and nonnegative"))
    from_voltage = _checked_line_complex(from_voltage_phasor, "from voltage phasor")
    to_voltage = _checked_line_complex(to_voltage_phasor, "to voltage phasor")
    previous_from = _checked_line_complex(previous_from_history, "previous from history")
    previous_to = _checked_line_complex(previous_to_history, "previous to history")
    scale = 1.0e-6
    angular_frequency = omega * scale
    pole_real = real(pole_value) * scale
    pole_imag = -imag(pole_value) * scale
    residue_real = real(residue_value) * scale
    residue_imag = imag(residue_value) * scale
    numerator_real = residue_real * pole_real - residue_imag * pole_imag
    numerator_imag = residue_real * angular_frequency
    denominator_real = pole_real * pole_real + pole_imag * pole_imag -
        angular_frequency * angular_frequency
    denominator_imag = 2.0 * pole_real * angular_frequency
    denominator_norm2 = denominator_real * denominator_real +
        denominator_imag * denominator_imag
    denominator_norm2 > 0.0 ||
        throw(ArgumentError("Semlyen harmonic history transfer denominator is singular"))
    inverse_denominator_real = denominator_real / denominator_norm2
    inverse_denominator_imag = denominator_imag / denominator_norm2
    transfer_real = inverse_denominator_real * numerator_real +
        inverse_denominator_imag * numerator_imag
    transfer_imag = inverse_denominator_real * numerator_imag -
        inverse_denominator_imag * numerator_real
    quadrature_numerator_real = residue_imag * pole_real + residue_real * pole_imag
    quadrature_numerator_imag = residue_imag * angular_frequency
    quadrature_real = inverse_denominator_real * quadrature_numerator_real +
        inverse_denominator_imag * quadrature_numerator_imag
    quadrature_imag = inverse_denominator_real * quadrature_numerator_imag -
        inverse_denominator_imag * quadrature_numerator_real
    transfer = ComplexF64(transfer_real, transfer_imag)
    quadrature_transfer = ComplexF64(quadrature_real, quadrature_imag)
    from_increment = ComplexF64(
        transfer_real * real(from_voltage) - transfer_imag * imag(from_voltage),
        quadrature_real * real(from_voltage) + quadrature_imag * imag(from_voltage),
    )
    to_increment = ComplexF64(
        transfer_real * real(to_voltage) - transfer_imag * imag(to_voltage),
        quadrature_real * real(to_voltage) + quadrature_imag * imag(to_voltage),
    )
    return SemlyenLineHarmonicHistoryUpdate(
        pole_value,
        residue_value,
        omega,
        transfer,
        quadrature_transfer,
        from_voltage,
        to_voltage,
        previous_from,
        previous_to,
        from_increment,
        to_increment,
        from_increment - previous_from,
        to_increment - previous_to,
        true,
    )
end

function _checked_line_step_response_frequency_samples(
    angular_frequencies_rad_s::AbstractVector,
    frequency_response_values::AbstractVector,
)
    angular = Float64.(collect(angular_frequencies_rad_s))
    response = Float64.(collect(frequency_response_values))
    length(angular) == length(response) ||
        throw(ArgumentError("step-response angular frequency and response counts must match"))
    length(angular) >= 4 ||
        throw(ArgumentError("step-response fitting needs at least four frequency samples"))
    angular[1] == 0.0 ||
        throw(ArgumentError("step-response first angular frequency must be zero"))
    all(isfinite, angular) && all(isfinite, response) ||
        throw(ArgumentError("step-response samples must be finite"))
    all(diff(angular) .> 0.0) ||
        throw(ArgumentError("step-response angular frequencies must be strictly increasing"))
    ratios = [angular[idx + 1] / angular[idx] for idx in 2:(length(angular) - 1)]
    ratio = first(ratios)
    ratio_tolerance = max(
        64.0 * eps(Float64) * max(1.0, abs(ratio)),
        1.0e-10 * max(1.0, abs(ratio)),
    )
    all(value -> abs(value - ratio) <= ratio_tolerance, ratios) ||
        throw(ArgumentError("step-response angular frequencies must follow a geometric ratio"))
    return angular, response, ratio
end

function _line_step_response_inverse_fourier_geometric(
    time_s::Float64,
    angular_frequencies::Vector{Float64},
    frequency_response::Vector{Float64},
    frequency_ratio::Float64,
)
    time_s > 0.0 ||
        throw(ArgumentError("step-response inverse transform time must be positive"))
    nfr = length(angular_frequencies) - 1
    twopi = 2.0 * pi
    pi2 = 4.0 / twopi
    dmin = frequency_ratio - 1.0
    dplu = frequency_ratio + 1.0
    abs(dmin) > 0.0 ||
        throw(ArgumentError("step-response frequency ratio must not be one"))
    c = Vector{Float64}(undef, nfr)
    time2 = time_s * time_s
    for idx in 1:nfr
        phase = time_s * angular_frequencies[idx + 1]
        cycles = trunc(Int, phase / twopi)
        c[idx] = sin(phase - cycles * twopi)
    end
    first_phase = angular_frequencies[2] * time_s
    first_phase != 0.0 ||
        throw(ArgumentError("step-response first nonzero frequency produces zero phase"))
    transformed = first_phase * 0.5 * (
        frequency_response[1] * time2 +
        frequency_response[2] * time_s * (2.0 * cos(first_phase) / first_phase + c[1])
    )
    transformed -= frequency_ratio / dmin * frequency_response[2] *
        (c[2] - c[1]) / angular_frequencies[3]
    for idx in 3:nfr
        transformed += frequency_response[idx] *
            (dplu * c[idx - 1] - frequency_ratio * c[idx - 2] - c[idx]) /
            (dmin * angular_frequencies[idx])
    end
    return pi2 * transformed / time2
end

function _line_step_response_initial_values(
    values::Vector{Float64},
    final_value::Float64,
    fit2z::Float64,
    control_mode::Int,
    npoint::Int,
    no::Int,
    pivot_threshold::Float64,
)
    threshold = final_value * fit2z
    control_mode < 0 && (threshold *= 0.1)
    start_index = findfirst(value -> value >= threshold, values)
    start_index === nothing &&
        throw(ArgumentError("step-response fit did not reach the initial threshold"))
    ncount = 0
    for idx in start_index:npoint
        if idx == 1
            if values[idx + 1] - values[idx] <= values[idx]
                ncount = idx
                break
            end
        elseif idx < npoint && values[idx + 1] - values[idx] <= values[idx] - values[idx - 1]
            ncount = idx
            break
        end
    end
    ncount > 0 ||
        throw(ArgumentError("step-response fit did not find a valid curvature anchor"))
    ntrd = div(npoint, 3)
    ntrd > 0 || throw(ArgumentError("step-response fit needs at least three points"))
    dy = (values[npoint] - values[npoint - ntrd]) * no / ntrd
    dy <= 0.0 && (dy = pivot_threshold)
    second_rate = abs(dy / (final_value - values[npoint]))
    ncount == 1 && (ncount = 2)
    dy = (values[ncount] - values[ncount - 1]) * no
    first_rate = dy / values[ncount] / 2.0
    delay_fraction = (ncount - 0.5) / no -
        0.5 * (values[ncount] + values[ncount - 1]) / dy
    if control_mode < 0
        first_rate = no * values[2] * 0.5 / maximum(values)
        delay_fraction = 0.0
        return first_rate, second_rate / 1.2, 0.0, 2
    end
    return first_rate, second_rate, delay_fraction, 3
end

function _line_step_response_time_domain_sequence(
    angular::Vector{Float64},
    response::Vector{Float64},
    ratio::Float64,
    final_value::Float64,
    fit_span_s::Float64,
    time_step_s::Float64,
    time_start_s::Float64,
    control_mode::Int,
)
    span = fit_span_s
    start_time = time_start_s
    point_count = trunc(Int, span / time_step_s + 0.5)
    point_count >= 4 ||
        throw(ArgumentError("step-response fitting needs at least four time samples"))
    critical_peak_index = trunc(Int, point_count * 0.8)
    time_zero_shift_count = 0
    span_division_count = 0
    values = Vector{Float64}(undef, point_count)
    while true
        time_step = span / point_count
        zero_index = 0
        peak_index = 0
        fill_start = 1
        time = control_mode < 0 ? 5.0 * time_step : time_step
        peak_value = -Inf
        while true
            for idx in fill_start:point_count
                value = _line_step_response_inverse_fourier_geometric(
                    time,
                    angular,
                    response,
                    ratio,
                )
                values[idx] = value
                if value > peak_value
                    peak_value = value
                    peak_index = idx
                end
                if value < final_value / 1000.0
                    zero_index = idx
                    peak_value = -Inf
                end
                time += time_step
            end
            if control_mode < 0 || zero_index <= 0
                break
            end
            time_zero_shift_count += 1
            time_zero_shift_count <= 5 ||
                throw(ArgumentError("step-response fit discarded invalid inverse-transform samples more than five times"))
            kept_count = point_count - zero_index
            peak_index -= zero_index
            if kept_count > 0
                values[1:kept_count] .= values[(zero_index + 1):point_count]
            end
            fill_start = kept_count + 1
            start_time += zero_index * time_step
            zero_index = 0
        end
        if control_mode >= 0 || peak_index >= critical_peak_index
            return values, span, start_time, time_step, point_count
        end
        span = peak_index * time_step / 2.0
        span_division_count += 1
        span_division_count <= 5 ||
            throw(ArgumentError("step-response fit reduced the fitting span more than five times"))
    end
end

function _line_step_response_fit_error_and_amplitude!(
    hhm::Vector{Float64},
    hhn::Vector{Float64},
    values::Vector{Float64},
    rates::Vector{Float64},
    final_value::Float64,
    no::Int,
    npoint::Int,
    variable_count::Int,
    with_jacobian::Bool,
    hac::Matrix{Float64},
    rhs::Vector{Float64},
)
    delay_fraction = max(rates[3], 0.0)
    first_rate = rates[1]
    second_rate = rates[2]
    first_rate > 0.0 && second_rate > 0.0 ||
        return NaN, NaN, 1
    first_index = max(1, trunc(Int, delay_fraction * no + 1.0))
    first_index <= npoint || return NaN, NaN, first_index
    numerator = 0.0
    denominator = 0.0
    for idx in first_index:npoint
        t = idx / no - delay_fraction
        hhm[idx] = exp(-min(first_rate * t, 80.0))
        hhn[idx] = exp(-min(second_rate * t, 80.0))
        delta = hhm[idx] - hhn[idx]
        numerator += delta * (final_value * (1.0 - hhn[idx]) - values[idx])
        denominator += delta * delta
    end
    denominator > 0.0 || return NaN, NaN, first_index
    first_amplitude = numerator / denominator
    second_amplitude = final_value - first_amplitude
    error = 0.0
    if with_jacobian
        fill!(hac, 0.0)
        fill!(rhs, 0.0)
    end
    for idx in first_index:npoint
        t = idx / no - delay_fraction
        fitted = final_value - first_amplitude * hhm[idx] - second_amplitude * hhn[idx]
        residual = fitted - values[idx]
        error += residual * residual
        if with_jacobian
            b1 = t * hhm[idx] * first_amplitude
            b2 = t * hhn[idx] * second_amplitude
            b3 = -hhm[idx] * first_amplitude * first_rate -
                hhn[idx] * second_amplitude * second_rate
            db11 = -t * t * hhm[idx] * first_amplitude
            db13 = -hhm[idx] * (1.0 - first_rate * t) * first_amplitude
            db22 = -t * t * hhn[idx] * second_amplitude
            db23 = -hhn[idx] * (1.0 - second_rate * t) * second_amplitude
            db33 = -first_rate * first_rate * first_amplitude * hhm[idx] -
                second_rate * second_rate * second_amplitude * hhn[idx]
            b = (b1, b2, b3)
            db = (
                (db11, 0.0, db13),
                (0.0, db22, db23),
                (db13, db23, db33),
            )
            for row in 1:variable_count
                rhs[row] -= residual * b[row]
                for col in 1:variable_count
                    hac[row, col] += b[row] * b[col] + residual * db[row][col]
                end
            end
        end
    end
    return error / (npoint * final_value * final_value), first_amplitude, first_index
end

function _line_step_response_solve_increment!(
    increment::Vector{Float64},
    matrix::Matrix{Float64},
    rhs::Vector{Float64},
    variable_count::Int,
    pivot_threshold::Float64,
)
    a = copy(@view matrix[1:variable_count, 1:variable_count])
    b = copy(@view rhs[1:variable_count])
    for pivot in 1:variable_count
        if abs(a[pivot, pivot]) < pivot_threshold
            swap = findfirst(row -> abs(a[row, pivot]) >= pivot_threshold, (pivot + 1):variable_count)
            swap === nothing && return false
            swap_row = pivot + swap
            a[pivot, :], a[swap_row, :] = copy(a[swap_row, :]), copy(a[pivot, :])
        end
        value = a[pivot, pivot]
        b[pivot] /= value
        a[pivot, pivot:variable_count] ./= value
        for row in 1:variable_count
            row == pivot && continue
            factor = a[row, pivot]
            b[row] -= factor * b[pivot]
            a[row, pivot:variable_count] .-= factor .* a[pivot, pivot:variable_count]
        end
    end
    increment[1:variable_count] .= b
    return true
end

function line_step_response_exponential_fit(
    angular_frequencies_rad_s::AbstractVector,
    frequency_response_values::AbstractVector;
    final_value::Real = 1.0,
    fit_span_s::Real,
    time_step_s::Real,
    time_start_s::Real = 0.0,
    fit_control_mode::Integer = 0,
    steady_state_angular_frequency_rad_s::Real = 2.0 * pi * 60.0,
    steady_state_shift = 1.0 + 0.0im,
    fitting_error_tolerance::Real = 0.5e-4,
    relative_increment_tolerance::Real = 0.5e-2,
    initial_rate_scan_step::Real = 0.1,
    pivot_threshold::Real = 1.0e-5,
    steady_state_error_tolerance::Real = 1.0e-4,
    pivot_tolerance::Real = 1.0e-16,
    maximum_iterations::Integer = 10,
    steady_state_maximum_iterations::Integer = 10,
)
    angular, response, ratio =
        _checked_line_step_response_frequency_samples(angular_frequencies_rad_s, frequency_response_values)
    final = _checked_line_positive_finite(final_value, "step-response final_value")
    fit_span = _checked_line_positive_finite(fit_span_s, "step-response fit_span_s")
    time_step = _checked_line_positive_finite(time_step_s, "step-response time_step_s")
    start_time = Float64(time_start_s)
    isfinite(start_time) && start_time >= 0.0 ||
        throw(ArgumentError("step-response time_start_s must be finite and nonnegative"))
    control_mode = Int(fit_control_mode)
    omega = Float64(steady_state_angular_frequency_rad_s)
    isfinite(omega) && omega >= 0.0 ||
        throw(ArgumentError("steady_state_angular_frequency_rad_s must be finite and nonnegative"))
    shift = _checked_line_complex(steady_state_shift, "steady-state shift")
    eps_fit = _checked_line_positive_finite(fitting_error_tolerance, "fitting_error_tolerance")
    eps_increment = _checked_line_positive_finite(relative_increment_tolerance, "relative_increment_tolerance")
    fit2z = _checked_line_positive_finite(initial_rate_scan_step, "initial_rate_scan_step")
    pivthr = _checked_line_positive_finite(pivot_threshold, "pivot_threshold")
    ft2emx = _checked_line_positive_finite(steady_state_error_tolerance, "steady_state_error_tolerance")
    epspv2 = _checked_line_positive_finite(pivot_tolerance, "pivot_tolerance")
    niter = Int(maximum_iterations)
    niter1 = Int(steady_state_maximum_iterations)
    niter > 0 && niter1 > 0 ||
        throw(ArgumentError("step-response iteration limits must be positive"))
    values, fit_span, start_time, time_step, npoint =
        _line_step_response_time_domain_sequence(
            angular,
            response,
            ratio,
            final,
            fit_span,
            time_step,
            start_time,
            control_mode,
        )
    first_rate, second_rate, delay_fraction, variable_count =
        _line_step_response_initial_values(values, final, fit2z, control_mode, npoint, npoint, pivthr)
    rates = [first_rate, second_rate, delay_fraction]
    first_rate_anchor = control_mode < 0 ? second_rate : first_rate
    best_anchor = first_rate_anchor
    previous_error = 100.0
    error = previous_error
    first_amplitude = 0.0
    iteration = 0
    gradient_scan_count = 0
    hhm = zeros(Float64, npoint)
    hhn = zeros(Float64, npoint)
    hac = zeros(Float64, 3, 3)
    rhs = zeros(Float64, 3)
    increment = zeros(Float64, 3)
    old_rates = copy(rates)
    old_increment = ones(Float64, 3)
    while true
        iteration += 1
        local_error = previous_error
        recovery_needed = false
        converged_by_increment = false
        advanced_by_increment = false
        while true
            gradient_scan_count += 1
            if gradient_scan_count <= 20
                if control_mode >= 0
                    rates[1] = first_rate_anchor * (1.0 + fit2z * (gradient_scan_count - 1))
                else
                    rates[2] = first_rate_anchor * (1.0 + fit2z * (gradient_scan_count - 1))
                end
            elseif gradient_scan_count == 21
                if control_mode >= 0
                    rates[1] = best_anchor
                else
                    rates[2] = best_anchor
                end
            end
            rates[3] = max(rates[3], 0.0)
            error, first_amplitude, _ = _line_step_response_fit_error_and_amplitude!(
                hhm,
                hhn,
                values,
                rates,
                final,
                npoint,
                npoint,
                variable_count,
                gradient_scan_count > 20,
                hac,
                rhs,
            )
            if !isfinite(error)
                throw(ArgumentError("step-response exponential fit encountered a singular trial"))
            end
            if gradient_scan_count > 20 && error > previous_error
                recovery_needed = true
                break
            end
            if error < previous_error
                previous_error = error
                best_anchor = control_mode < 0 ? rates[2] : rates[1]
            end
            gradient_scan_count <= 20 && continue
            iteration >= niter && break
            error < eps_fit && break
            _line_step_response_solve_increment!(increment, hac, rhs, variable_count, epspv2) ||
                throw(ArgumentError("step-response exponential fit Jacobian is singular"))
            old_rates .= rates
            old_increment .= increment
            rates[1:variable_count] .+= increment[1:variable_count]
            if all(idx -> abs(increment[idx] / rates[idx]) <= eps_increment, 1:variable_count)
                converged_by_increment = true
                break
            end
            local_error = error
            advanced_by_increment = true
            break
        end
        if iteration >= niter || error < eps_fit ||
           converged_by_increment
            break
        end
        advanced_by_increment && continue
        recovery_needed || continue
        scale = 0.1
        recovered = false
        for _ in 1:niter
            trial = copy(old_rates)
            for idx in 1:variable_count
                trial[idx] = old_rates[idx] + scale * old_rates[idx] *
                    old_increment[idx] / abs(old_increment[idx])
            end
            trial_error, trial_amplitude, _ = _line_step_response_fit_error_and_amplitude!(
                hhm,
                hhn,
                values,
                trial,
                final,
                npoint,
                npoint,
                variable_count,
                true,
                hac,
                rhs,
            )
            if isfinite(trial_error) && trial_error <= local_error
                rates .= trial
                error = trial_error
                first_amplitude = trial_amplitude
                if error < previous_error
                    previous_error = error
                    best_anchor = control_mode < 0 ? rates[2] : rates[1]
                end
                _line_step_response_solve_increment!(increment, hac, rhs, variable_count, epspv2) ||
                    throw(ArgumentError("step-response exponential fit Jacobian is singular"))
                old_rates .= rates
                old_increment .= increment
                rates[1:variable_count] .+= increment[1:variable_count]
                converged_by_increment =
                    all(idx -> abs(increment[idx] / rates[idx]) <= eps_increment, 1:variable_count)
                recovered = true
                break
            end
            scale /= 10.0
        end
        recovered || throw(ArgumentError("step-response exponential fit failed to recover decreasing error"))
        converged_by_increment && break
    end
    delay = rates[3] * fit_span + start_time
    first_time_constant = fit_span / rates[1]
    control_mode < 0 && (control_mode = -control_mode - 1)
    delay_phase = omega * delay
    rotated_shift = ComplexF64(
        real(shift) * cos(delay_phase) + imag(shift) * sin(delay_phase),
        real(shift) * sin(delay_phase) - imag(shift) * cos(delay_phase),
    )
    second_amplitude = final - first_amplitude
    second_time_constant = fit_span / rates[2]
    steady_state_adjusted = false
    if control_mode != 0
        trial_error = fit_span / rates[2]
        if trial_error / first_time_constant >= 10.0
            first_trial = first_amplitude
            second_trial = second_amplitude
            for _ in 1:niter1
                ddt = omega * first_time_constant
                zp = 1.0 + ddt * ddt
                ddt /= zp
                zr = real(rotated_shift) - first_trial / zp
                zi = imag(rotated_shift) + first_trial * ddt
                d1 = omega * trial_error
                d2 = 1.0 + d1 * d1
                d3 = zr - second_trial / d2
                d4 = zi + second_trial * d1 / d2
                steady_error = hypot(d3, d4)
                steady_error < ft2emx * abs(rotated_shift) && break
                new_error = -zi / zr / omega
                new_error <= 0.0 && break
                second_trial = (zr * zr + zi * zi) / zr
                first_trial = final - second_trial
                trial_error = new_error
                steady_state_adjusted = true
            end
            if steady_state_adjusted
                first_amplitude = first_trial
                second_amplitude = second_trial
                second_time_constant = trial_error
            end
        end
    end
    final_error = 0.0
    for idx in 1:npoint
        shifted_time = idx * time_step - (delay - start_time)
        fitted = 0.0
        if shifted_time >= 0.0
            first_decay = exp(-min(shifted_time / first_time_constant, 30.0))
            second_decay = exp(-min(shifted_time / second_time_constant, 30.0))
            fitted = first_amplitude * (1.0 - first_decay) +
                second_amplitude * (1.0 - second_decay)
        end
        residual = fitted - values[idx]
        final_error += residual * residual
    end
    final_error /= npoint * final * final
    return LineStepResponseExponentialFitResult(
        angular,
        response,
        values,
        final,
        fit_span,
        start_time,
        time_step,
        first_amplitude,
        first_time_constant,
        second_amplitude,
        second_time_constant,
        delay,
        rotated_shift,
        final_error,
        iteration,
        variable_count,
        true,
        steady_state_adjusted,
    )
end

function frequency_dependent_line_recursive_convolution_fit_from_step_response(
    sample_frequencies_hz::AbstractVector,
    modal_response_samples::AbstractMatrix,
    step_response_angular_frequencies_rad_s::AbstractVector,
    step_response_frequency_values::AbstractMatrix,
    dt_s::Real;
    final_values = ones(size(step_response_frequency_values, 1)),
    fit_span_s::Real,
    time_step_s::Real,
    time_start_s::Real = 0.0,
    fit_control_mode::Integer = 0,
    steady_state_angular_frequency_rad_s::Real = 2.0 * pi * 60.0,
    steady_state_shift = 1.0 + 0.0im,
)
    frequencies, samples = _checked_line_recursive_fit_samples(
        sample_frequencies_hz,
        modal_response_samples,
    )
    mode_count, _ = size(samples)
    response_matrix = Matrix{Float64}(step_response_frequency_values)
    size(response_matrix, 1) == mode_count ||
        throw(ArgumentError("step-response frequency value rows must match mode count"))
    finals = Float64.(collect(final_values))
    length(finals) == mode_count ||
        throw(ArgumentError("step-response final value count must match mode count"))
    fits = Vector{LineStepResponseExponentialFitResult}(undef, mode_count)
    pole_decay = Matrix{ComplexF64}(undef, mode_count, 2)
    dt = _checked_line_positive_finite(dt_s, "dt_s")
    for mode in 1:mode_count
        fit = line_step_response_exponential_fit(
            step_response_angular_frequencies_rad_s,
            view(response_matrix, mode, :);
            final_value = finals[mode],
            fit_span_s = fit_span_s,
            time_step_s = time_step_s,
            time_start_s = time_start_s,
            fit_control_mode = fit_control_mode,
            steady_state_angular_frequency_rad_s = steady_state_angular_frequency_rad_s,
            steady_state_shift = steady_state_shift,
        )
        fits[mode] = fit
        pole_decay[mode, 1] = exp(-dt / fit.first_time_constant_s)
        pole_decay[mode, 2] = exp(-dt / fit.second_time_constant_s)
    end
    recursive_fit = frequency_dependent_line_recursive_convolution_fit(
        frequencies,
        samples,
        pole_decay,
        dt,
    )
    return (
        recursive_fit = recursive_fit,
        step_response_fits = fits,
        pole_decay = pole_decay,
        full_bpa_frequency_dependent_fitting_executed = true,
        tdfit_step_response_fit_executed = true,
    )
end

function frequency_dependent_line_recursive_convolution_fit(
    sample_frequencies_hz::AbstractVector,
    modal_response_samples::AbstractMatrix,
    pole_decay::AbstractMatrix,
    dt_s::Real,
)
    frequencies, samples = _checked_line_recursive_fit_samples(
        sample_frequencies_hz,
        modal_response_samples,
    )
    mode_count, sample_count = size(samples)
    decay = _checked_line_complex_rectangular_matrix(
        pole_decay,
        mode_count,
        "line recursive convolution fit pole decay",
    )
    all(value -> abs(value) <= 1.0 + 64.0 * eps(Float64), decay) ||
        throw(ArgumentError("line recursive convolution fit pole decays must be stable"))
    dt = Float64(dt_s)
    isfinite(dt) && dt > 0.0 ||
        throw(ArgumentError("dt_s must be finite and positive"))
    term_count = size(decay, 2)
    sample_count >= term_count ||
        throw(ArgumentError("line recursive convolution fit needs at least one sample per pole term"))
    residue = Matrix{ComplexF64}(undef, mode_count, term_count)
    fitted = zeros(ComplexF64, mode_count, sample_count)
    for mode in 1:mode_count
        basis = Matrix{ComplexF64}(undef, sample_count, term_count)
        for sample in 1:sample_count, term in 1:term_count
            basis[sample, term] = _line_recursive_convolution_basis(
                frequencies[sample],
                dt,
                decay[mode, term],
            )
        end
        coefficients = basis \ vec(samples[mode, :])
        all(value -> isfinite(real(value)) && isfinite(imag(value)), coefficients) ||
            throw(ArgumentError("line recursive convolution fit residue solution must be finite"))
        residue[mode, :] .= coefficients
        fitted[mode, :] .= basis * coefficients
    end
    residual = samples .- fitted
    max_abs_error = isempty(residual) ? 0.0 : maximum(abs.(residual))
    return LineRecursiveConvolutionFitResult(
        collect(frequencies),
        decay,
        residue,
        samples,
        fitted,
        residual,
        max_abs_error,
        dt,
        mode_count,
        term_count,
    )
end

function FrequencyDependentLineRecursiveConvolutionState(
    yz_sample_matrices::AbstractVector,
    sample_rows::AbstractVector,
    initial_frequency_hz::Real,
    line_length::Real,
    pole_decay::AbstractMatrix,
    residue::AbstractMatrix;
    initial_ntol::Integer = 1,
    initial_nrp::Integer = 0,
    skin_effect_internal_impedance_executed::Bool = false,
    earth_return_impedance_executed::Bool = false,
    frequency_dependent_fitting_executed::Bool = false,
    frequency_loop_executed::Bool = false,
    pipe_sheath_side_effects_executed::Bool = false,
)
    modal_sample_state = FrequencyDependentLineModalSampleRuntimeState(
        yz_sample_matrices,
        sample_rows,
        initial_frequency_hz,
        line_length;
        ntol = initial_ntol,
        nrp = initial_nrp,
    )
    mode_count = length(modal_sample_state.sample_rows[1])
    decay, terms = _checked_line_recursive_convolution_matrices(
        pole_decay,
        residue,
        mode_count,
    )
    term_count = size(decay, 2)
    return FrequencyDependentLineRecursiveConvolutionState(
        modal_sample_state,
        decay,
        terms,
        zeros(ComplexF64, mode_count, term_count),
        zeros(ComplexF64, mode_count, term_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        skin_effect_internal_impedance_executed,
        earth_return_impedance_executed,
        frequency_dependent_fitting_executed,
        frequency_loop_executed,
        pipe_sheath_side_effects_executed,
        0,
    )
end

function FrequencyDependentLineRecursiveConvolutionState(
    yz_sample_matrices::AbstractVector,
    sample_rows::AbstractVector,
    initial_frequency_hz::Real,
    line_length::Real,
    fit::LineRecursiveConvolutionFitResult;
    initial_ntol::Integer = 1,
    initial_nrp::Integer = 0,
    skin_effect_internal_impedance_executed::Bool = false,
    earth_return_impedance_executed::Bool = false,
    frequency_dependent_fitting_executed::Bool = false,
    frequency_loop_executed::Bool = false,
    pipe_sheath_side_effects_executed::Bool = false,
)
    return FrequencyDependentLineRecursiveConvolutionState(
        yz_sample_matrices,
        sample_rows,
        initial_frequency_hz,
        line_length,
        fit.pole_decay,
        fit.residue;
        initial_ntol = initial_ntol,
        initial_nrp = initial_nrp,
        skin_effect_internal_impedance_executed = skin_effect_internal_impedance_executed,
        earth_return_impedance_executed = earth_return_impedance_executed,
        frequency_dependent_fitting_executed = frequency_dependent_fitting_executed,
        frequency_loop_executed = frequency_loop_executed,
        pipe_sheath_side_effects_executed = pipe_sheath_side_effects_executed,
    )
end

function _checked_line_complex_vector(
    values::AbstractVector,
    expected_count::Int,
    label::AbstractString,
)
    length(values) == expected_count ||
        throw(ArgumentError("$label count must match the line phase/mode count"))
    checked = Vector{ComplexF64}(undef, expected_count)
    for idx in 1:expected_count
        value = ComplexF64(values[idx])
        all(isfinite, (real(value), imag(value))) ||
            throw(ArgumentError("$label entries must be finite"))
        checked[idx] = value
    end
    return checked
end

function _apply_line_transform!(
    output::AbstractVector{ComplexF64},
    matrix::AbstractMatrix{ComplexF64},
    input::AbstractVector,
    label::AbstractString,
)
    n = size(matrix, 1)
    length(input) == n ||
        throw(ArgumentError("$label input vector length must match the transform dimension"))
    length(output) >= n ||
        throw(ArgumentError("$label output vector is shorter than the transform dimension"))
    for row in 1:n
        total = 0.0 + 0.0im
        for col in 1:n
            value = ComplexF64(input[col])
            isfinite(real(value)) && isfinite(imag(value)) ||
                throw(ArgumentError("$label input vector entries must be finite"))
            total += matrix[row, col] * value
        end
        output[row] = total
    end
    return output
end

function line_modal_transform!(
    modal_values::AbstractVector{ComplexF64},
    transform::LineModalTransform,
    phase_values::AbstractVector,
)
    return _apply_line_transform!(modal_values, transform.phase_to_modal, phase_values, "line_modal_transform")
end

function line_modal_transform(transform::LineModalTransform, phase_values::AbstractVector)
    modal_values = Vector{ComplexF64}(undef, size(transform.phase_to_modal, 1))
    return line_modal_transform!(modal_values, transform, phase_values)
end

function line_phase_transform!(
    phase_values::AbstractVector{ComplexF64},
    transform::LineModalTransform,
    modal_values::AbstractVector,
)
    return _apply_line_transform!(phase_values, transform.modal_to_phase, modal_values, "line_phase_transform")
end

function line_phase_transform(transform::LineModalTransform, modal_values::AbstractVector)
    phase_values = Vector{ComplexF64}(undef, size(transform.modal_to_phase, 1))
    return line_phase_transform!(phase_values, transform, modal_values)
end

function _checked_line_complex_square_matrix(matrix::AbstractMatrix, label::AbstractString)
    rows, cols = size(matrix)
    rows == cols && rows > 0 ||
        throw(ArgumentError("$label must be a nonempty square matrix"))
    checked = Matrix{ComplexF64}(undef, rows, cols)
    for row in 1:rows, col in 1:cols
        value = ComplexF64(matrix[row, col])
        isfinite(real(value)) && isfinite(imag(value)) ||
            throw(ArgumentError("$label entries must be finite"))
        checked[row, col] = value
    end
    return checked
end

function _checked_line_complex_rectangular_matrix(
    matrix::AbstractMatrix,
    expected_rows::Int,
    label::AbstractString,
)
    rows, cols = size(matrix)
    rows == expected_rows && cols > 0 ||
        throw(ArgumentError("$label must have one row per line mode and at least one column"))
    checked = Matrix{ComplexF64}(undef, rows, cols)
    for row in 1:rows, col in 1:cols
        value = ComplexF64(matrix[row, col])
        isfinite(real(value)) && isfinite(imag(value)) ||
            throw(ArgumentError("$label entries must be finite"))
        checked[row, col] = value
    end
    return checked
end

function _normalize_line_modal_vectors!(vectors::Matrix{ComplexF64})
    for col in axes(vectors, 2)
        norm_value = sqrt(sum(abs2, vectors[:, col]))
        norm_value > 0.0 || throw(ArgumentError("line modal eigenvectors must be nonzero"))
        max_row = argmax(abs.(vectors[:, col]))
        phase = angle(vectors[max_row, col])
        vectors[:, col] .*= cis(-phase) / norm_value
        real(vectors[max_row, col]) < 0.0 && (vectors[:, col] .*= -1.0)
    end
    return vectors
end

function line_modal_eigen_order(eigenvalues::AbstractVector)
    isempty(eigenvalues) && throw(ArgumentError("line modal eigenvalue vector must be nonempty"))
    values = ComplexF64.(eigenvalues)
    all(value -> isfinite(real(value)) && isfinite(imag(value)), values) ||
        throw(ArgumentError("line modal eigenvalues must be finite"))
    order = collect(eachindex(values))
    magnitudes = abs2.(values)
    for pass in 1:length(values)
        limit = length(values) - pass
        for idx in 1:limit
            if magnitudes[idx] < magnitudes[idx + 1]
                magnitudes[idx], magnitudes[idx + 1] = magnitudes[idx + 1], magnitudes[idx]
                order[idx], order[idx + 1] = order[idx + 1], order[idx]
            end
        end
    end
    return order
end

function _line_modal_solution_metrics(eigenvalues::AbstractVector{ComplexF64}, frequency_hz::Float64)
    metrics = Matrix{Float64}(undef, 4, length(eigenvalues))
    omega = 2.0 * pi * frequency_hz
    for (idx, value) in pairs(eigenvalues)
        root = sqrt(value)
        magnitude = abs(value)
        angle_degrees = angle(value) * 180.0 / pi
        velocity = imag(root) == 0.0 ? copysign(Inf, omega) : omega / imag(root)
        attenuation = real(root) / sqrt(frequency_hz)
        isfinite(velocity) && isfinite(attenuation) ||
            throw(ArgumentError("line modal solution metrics require finite propagation roots"))
        metrics[:, idx] .= (magnitude / frequency_hz, angle_degrees, velocity, attenuation)
    end
    return metrics
end

function line_modal_solution(
    yz_matrix::AbstractMatrix,
    frequency_hz::Real;
    unwind_state::Union{Nothing,LineModeUnwindState} = nothing,
    ntol::Integer = 1,
    nrp::Integer = 0,
)
    frequency = Float64(frequency_hz)
    isfinite(frequency) && frequency > 0.0 ||
        throw(ArgumentError("frequency_hz must be finite and positive"))
    matrix = _checked_line_complex_square_matrix(yz_matrix, "line modal YZ matrix")
    mode_count = size(matrix, 1)
    factorization = eigen(matrix)
    raw_values = ComplexF64.(factorization.values)
    raw_vectors = Matrix{ComplexF64}(factorization.vectors)
    eigen_order = line_modal_eigen_order(raw_values)
    values = raw_values[eigen_order]
    vectors = raw_vectors[:, eigen_order]
    _normalize_line_modal_vectors!(vectors)
    metrics = _line_modal_solution_metrics(values, frequency)
    if unwind_state === nothing
        sequence = collect(1:mode_count)
        ordered_metrics = copy(metrics)
        update_count = 0
    else
        order_result = line_mode_unwind!(unwind_state, metrics; ntol = ntol, nrp = nrp)
        sequence = order_result.sequence
        ordered_metrics = order_result.ordered_metrics
        update_count = order_result.update_count
    end
    ordered_values = values[sequence]
    ordered_vectors = vectors[:, sequence]
    modal_eigen_order = eigen_order[sequence]
    transform = LineModalTransform(inv(ordered_vectors), ordered_vectors)
    eigen_residual = matrix * ordered_vectors - ordered_vectors * Diagonal(ordered_values)
    inverse_product = transform.modal_to_phase * transform.phase_to_modal
    return LineModalSolution(
        transform,
        ordered_values,
        sqrt.(ordered_values),
        copy(modal_eigen_order),
        copy(sequence),
        copy(metrics),
        copy(ordered_metrics),
        maximum(abs.(eigen_residual)),
        maximum(abs.(inverse_product - I)),
        update_count,
    )
end

function _checked_line_modal_solution_scan_inputs(
    yz_matrices::AbstractVector,
    frequencies_hz::AbstractVector,
)
    isempty(yz_matrices) && throw(ArgumentError("line modal solution scan must contain at least one matrix"))
    length(yz_matrices) == length(frequencies_hz) ||
        throw(ArgumentError("line modal solution scan matrix count must match frequency count"))
    frequencies = Float64.(frequencies_hz)
    all(frequency -> isfinite(frequency) && frequency > 0.0, frequencies) ||
        throw(ArgumentError("line modal solution scan frequencies must be finite and positive"))
    order = sortperm(frequencies)
    sorted_frequencies = frequencies[order]
    for idx in 2:length(sorted_frequencies)
        abs(sorted_frequencies[idx] - sorted_frequencies[idx - 1]) >
            _line_frequency_row_tolerance(sorted_frequencies[idx]) ||
            throw(ArgumentError("line modal solution scan frequencies must be unique"))
    end
    first_matrix = _checked_line_complex_square_matrix(
        yz_matrices[order[1]],
        "line modal solution scan YZ matrix",
    )
    mode_count = size(first_matrix, 1)
    matrices = Vector{Matrix{ComplexF64}}(undef, length(order))
    matrices[1] = first_matrix
    for idx in 2:length(order)
        matrix = _checked_line_complex_square_matrix(
            yz_matrices[order[idx]],
            "line modal solution scan YZ matrix",
        )
        size(matrix, 1) == mode_count ||
            throw(ArgumentError("line modal solution scan matrices must share one mode count"))
        matrices[idx] = matrix
    end
    return sorted_frequencies, matrices, collect(order), mode_count
end

function line_modal_solution_scan!(
    mode_order_state::LineModeUnwindState,
    yz_matrices::AbstractVector,
    frequencies_hz::AbstractVector;
    ntol::Integer = 2,
    nrp::Integer = 0,
)
    frequencies, matrices, order, mode_count =
        _checked_line_modal_solution_scan_inputs(yz_matrices, frequencies_hz)
    _check_line_mode_state!(mode_order_state, mode_count)
    solutions = LineModalSolution[]
    sizehint!(solutions, length(frequencies))
    for idx in eachindex(frequencies)
        push!(
            solutions,
            line_modal_solution(
                matrices[idx],
                frequencies[idx];
                unwind_state = mode_order_state,
                ntol = ntol,
                nrp = nrp,
            ),
        )
    end
    residuals = [solution.eigenvector_residual_max_abs_error for solution in solutions]
    inverse_errors = [solution.inverse_product_max_abs_error for solution in solutions]
    return LineModalSolutionScan(
        collect(frequencies),
        2.0 .* pi .* frequencies,
        [copy(matrix) for matrix in matrices],
        solutions,
        [copy(solution.eigenvalues) for solution in solutions],
        [copy(solution.modal_eigen_order) for solution in solutions],
        [copy(solution.mode_sequence) for solution in solutions],
        order,
        mode_order_state.update_count,
        length(frequencies),
        mode_count,
        maximum(residuals),
        maximum(inverse_errors),
    )
end

function line_modal_solution_scan(
    yz_matrices::AbstractVector,
    frequencies_hz::AbstractVector;
    ntol::Integer = 2,
    nrp::Integer = 0,
)
    _, _, _, mode_count =
        _checked_line_modal_solution_scan_inputs(yz_matrices, frequencies_hz)
    mode_order_state = LineModeUnwindState(mode_count)
    return line_modal_solution_scan!(
        mode_order_state,
        yz_matrices,
        frequencies_hz;
        ntol = ntol,
        nrp = nrp,
    )
end

function frequency_dependent_line_modal_response!(
    state::FrequencyDependentLineModalState,
    phase_voltage::AbstractVector,
)
    mode_count = _line_modal_dimension(state.transform)
    state.frequency_points = _checked_line_frequency_points(state.frequency_points, mode_count)
    length(phase_voltage) == mode_count ||
        throw(ArgumentError("phase voltage count must match the modal transform dimension"))
    for idx in 1:mode_count
        value = ComplexF64(phase_voltage[idx])
        all(isfinite, (real(value), imag(value))) ||
            throw(ArgumentError("phase voltage values must be finite"))
        state.phase_voltage[idx] = value
    end
    line_modal_transform!(state.modal_voltage, state.transform, state.phase_voltage)
    for idx in 1:mode_count
        point = state.frequency_points[idx]
        admittance = inv(point.characteristic_impedance)
        state.modal_admittance[idx] = admittance
        state.modal_current[idx] = admittance * state.modal_voltage[idx]
        state.propagated_modal_current[idx] = point.propagation_factor * state.modal_current[idx]
    end
    line_phase_transform!(state.phase_current, state.transform, state.modal_current)
    line_phase_transform!(state.propagated_phase_current, state.transform, state.propagated_modal_current)
    state.update_count += 1
    return FrequencyDependentLineModalResponse(
        copy(state.phase_voltage),
        copy(state.modal_voltage),
        copy(state.modal_admittance),
        copy(state.modal_current),
        copy(state.propagated_modal_current),
        copy(state.phase_current),
        copy(state.propagated_phase_current),
    )
end

function frequency_dependent_line_modal_response(
    transform::LineModalTransform,
    frequency_points::AbstractVector{LineFrequencyPoint},
    phase_voltage::AbstractVector,
)
    state = FrequencyDependentLineModalState(transform, frequency_points)
    return frequency_dependent_line_modal_response!(state, phase_voltage)
end

function frequency_dependent_line_runtime_update!(
    state::FrequencyDependentLineRuntimeState,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector,
)
    mode_count = _line_modal_dimension(state.modal_state.transform)
    length(state.from_phase_voltage) == mode_count ||
        throw(ArgumentError("runtime state phase-voltage storage must match the transform dimension"))
    from_values = _checked_line_complex_vector(from_phase_voltage, mode_count, "from_phase_voltage")
    to_values = _checked_line_complex_vector(to_phase_voltage, mode_count, "to_phase_voltage")
    sending_before = copy(state.sending_phase_current)
    receiving_before = copy(state.receiving_phase_current)
    update_count_before = state.update_count

    state.previous_sending_phase_current .= state.sending_phase_current
    state.previous_receiving_phase_current .= state.receiving_phase_current
    state.from_phase_voltage .= from_values
    state.to_phase_voltage .= to_values
    state.phase_voltage_difference .= state.from_phase_voltage .- state.to_phase_voltage
    response = frequency_dependent_line_modal_response!(
        state.modal_state,
        state.phase_voltage_difference,
    )
    state.sending_phase_current .= response.phase_current
    state.receiving_phase_current .= .-response.propagated_phase_current
    state.update_count += 1

    sending_mutated = state.sending_phase_current != sending_before
    receiving_mutated = state.receiving_phase_current != receiving_before
    update_count_mutated = state.update_count != update_count_before
    return (
        source = :frequency_dependent_line_runtime_update,
        outcome = :state_mutation,
        fortran_files = (:OVER10_FOR, :OVER11_FOR, :OVER12_FOR, :OVER13_FOR),
        fortran_routines = (:OVER10, :SSEQIV, :OVER11, :OVER12, :OVER13),
        phase_voltage_difference = copy(state.phase_voltage_difference),
        modal_response = response,
        sending_phase_current = copy(state.sending_phase_current),
        receiving_phase_current = copy(state.receiving_phase_current),
        previous_sending_phase_current = copy(state.previous_sending_phase_current),
        previous_receiving_phase_current = copy(state.previous_receiving_phase_current),
        update_count = state.update_count,
        sending_phase_current_mutated = sending_mutated,
        receiving_phase_current_mutated = receiving_mutated,
        update_count_mutated = update_count_mutated,
        runtime_state_mutated = sending_mutated || receiving_mutated || update_count_mutated,
        legacy_fortran_in_loop = false,
        full_frequency_dependent_line_runtime_executed = false,
        deferred_calls = (:frequency_dependent_fitting, :line_timestep_bulk_oracle),
    )
end

function frequency_dependent_line_runtime_update(
    transform::LineModalTransform,
    frequency_points::AbstractVector{LineFrequencyPoint},
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector,
)
    state = FrequencyDependentLineRuntimeState(transform, frequency_points)
    return frequency_dependent_line_runtime_update!(state, from_phase_voltage, to_phase_voltage)
end

function frequency_dependent_line_runtime_update_from_fit!(
    state::FrequencyDependentLineRuntimeState,
    fit::LineFrequencySampleFitResult,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector,
)
    mode_count = _line_modal_dimension(state.modal_state.transform)
    fit.mode_count == mode_count ||
        throw(ArgumentError("frequency fit mode count must match the modal transform dimension"))
    state.modal_state.frequency_points = _checked_line_frequency_points(fit.frequency_points, mode_count)
    result = frequency_dependent_line_runtime_update!(state, from_phase_voltage, to_phase_voltage)
    return merge(
        result,
        (
            source = :frequency_dependent_line_runtime_update_from_fit,
            frequency_hz = fit.frequency_hz,
            frequency_fit_lower_hz = fit.lower_frequency_hz,
            frequency_fit_upper_hz = fit.upper_frequency_hz,
            frequency_fit_interpolation_weight = fit.interpolation_weight,
            frequency_fit_exact_match = fit.exact_frequency_match,
            frequency_fit_line_length = fit.line_length,
            frequency_fit_mode_count = fit.mode_count,
            frequency_fit = fit,
            bounded_frequency_sample_fit_executed = true,
            full_frequency_dependent_line_fitting_executed = false,
            deferred_calls = (:full_bpa_frequency_dependent_fitting, :line_timestep_bulk_oracle),
        ),
    )
end

function frequency_dependent_line_runtime_update_from_fit(
    transform::LineModalTransform,
    fit::LineFrequencySampleFitResult,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector,
)
    state = FrequencyDependentLineRuntimeState(transform, fit)
    return frequency_dependent_line_runtime_update_from_fit!(state, fit, from_phase_voltage, to_phase_voltage)
end

function frequency_dependent_line_sample_runtime_update!(
    state::FrequencyDependentLineSampleRuntimeState,
    target_frequency_hz::Real,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector,
)
    target_frequency = _checked_line_target_frequency(target_frequency_hz)
    previous_frequency = state.current_frequency_hz
    previous_weight = state.current_interpolation_weight
    fit = _line_frequency_sample_fit_from_sorted(
        state.sample_frequencies_hz,
        state.sample_rows,
        length(state.sample_rows[1]),
        target_frequency,
        state.line_length,
    )
    state.previous_frequency_hz = previous_frequency
    state.previous_interpolation_weight = previous_weight
    state.current_fit = fit
    state.current_frequency_hz = fit.frequency_hz
    state.current_interpolation_weight = fit.interpolation_weight
    state.frequency_update_count += 1
    result = frequency_dependent_line_runtime_update_from_fit!(
        state.runtime_state,
        fit,
        from_phase_voltage,
        to_phase_voltage,
    )
    return merge(
        result,
        (
            source = :frequency_dependent_line_sample_runtime_update,
            fortran_files = (:OVER10_FOR, :OVER11_FOR, :OVER12_FOR, :OVER13_FOR, :OVER44_FOR, :OVER45_FOR),
            fortran_routines = (:OVER10, :SSEQIV, :OVER11, :OVER12, :OVER13, :GUTS45),
            sample_frequency_hz = copy(state.sample_frequencies_hz),
            sample_row_count = length(state.sample_rows),
            previous_frequency_hz = state.previous_frequency_hz,
            current_frequency_hz = state.current_frequency_hz,
            previous_frequency_fit_interpolation_weight = state.previous_interpolation_weight,
            current_frequency_fit_interpolation_weight = state.current_interpolation_weight,
            frequency_changed = abs(state.current_frequency_hz - previous_frequency) >
                _line_frequency_row_tolerance(state.current_frequency_hz),
            frequency_update_count = state.frequency_update_count,
            runtime_update_count = state.runtime_state.update_count,
            bounded_frequency_scan_runtime_executed = true,
            full_frequency_dependent_line_fitting_executed = false,
            full_frequency_dependent_line_runtime_executed = false,
            deferred_calls = (
                :full_bpa_frequency_dependent_fitting,
                :recursive_convolution_line_runtime_oracle,
                :line_timestep_bulk_oracle,
            ),
        ),
    )
end

function frequency_dependent_line_sample_runtime_update(
    transform::LineModalTransform,
    sample_rows::AbstractVector,
    target_frequency_hz::Real,
    line_length::Real,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector;
    initial_frequency_hz::Real = target_frequency_hz,
)
    state = FrequencyDependentLineSampleRuntimeState(
        transform,
        sample_rows,
        initial_frequency_hz,
        line_length,
    )
    return frequency_dependent_line_sample_runtime_update!(
        state,
        target_frequency_hz,
        from_phase_voltage,
        to_phase_voltage,
    )
end

function frequency_dependent_line_modal_sample_runtime_update!(
    state::FrequencyDependentLineModalSampleRuntimeState,
    target_frequency_hz::Real,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector;
    ntol::Integer = 2,
    nrp::Integer = 0,
)
    target_frequency = _checked_line_target_frequency(target_frequency_hz)
    previous_frequency = state.current_frequency_hz
    matrix = _line_modal_sample_matrix_from_sorted(
        state.sample_frequencies_hz,
        state.modal_yz_matrices,
        target_frequency,
    )
    solution = line_modal_solution(
        matrix,
        target_frequency;
        unwind_state = state.mode_order_state,
        ntol = ntol,
        nrp = nrp,
    )
    state.previous_frequency_hz = previous_frequency
    state.current_frequency_hz = target_frequency
    state.current_modal_solution = solution
    state.modal_solution_update_count += 1
    state.sample_runtime_state.runtime_state.modal_state.transform = solution.transform
    result = frequency_dependent_line_sample_runtime_update!(
        state.sample_runtime_state,
        target_frequency,
        from_phase_voltage,
        to_phase_voltage,
    )
    return merge(
        result,
        (
            source = :frequency_dependent_line_modal_sample_runtime_update,
            fortran_files = (:OVER10_FOR, :OVER11_FOR, :OVER12_FOR, :OVER13_FOR, :OVER44_FOR, :OVER45_FOR),
            fortran_routines = (:OVER10, :SSEQIV, :OVER11, :OVER12, :OVER13, :MODAL, :DCEIGN, :COMLR2, :UNWIND, :GUTS45),
            modal_solution = solution,
            modal_eigenvalues = copy(solution.eigenvalues),
            modal_propagation_roots = copy(solution.propagation_roots),
            modal_eigen_order = copy(solution.modal_eigen_order),
            mode_sequence = copy(solution.mode_sequence),
            modal_ordered_metrics = copy(solution.ordered_mode_metrics),
            modal_eigenvector_residual_max_abs_error =
                solution.eigenvector_residual_max_abs_error,
            modal_inverse_product_max_abs_error =
                solution.inverse_product_max_abs_error,
            previous_modal_frequency_hz = state.previous_frequency_hz,
            current_modal_frequency_hz = state.current_frequency_hz,
            modal_solution_update_count = state.modal_solution_update_count,
            mode_order_update_count = state.mode_order_state.update_count,
            bounded_modal_sample_runtime_executed = true,
            full_dceign_lr_eigenvector_bulk_executed = false,
            full_frequency_dependent_line_fitting_executed = false,
            full_frequency_dependent_line_runtime_executed = false,
            deferred_calls = (
                :full_dceign_lr_modal_eigenvector_bulk_package,
                :full_bpa_frequency_dependent_fitting,
                :recursive_convolution_line_runtime_oracle,
                :line_timestep_bulk_oracle,
            ),
        ),
    )
end

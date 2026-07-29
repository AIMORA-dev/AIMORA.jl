export SaturationPreprocessResult,
       HysteresisLoopPreprocessResult,
       ZincOxidePowerSegment,
       ZincOxideFitResult,
       saturation_from_incremental_inductance,
       saturation_from_rms_excitation,
       saturation_runtime_table,
       normalized_hysteresis_loop,
       zinc_oxide_piecewise_fit,
       zinc_oxide_runtime_table,
       zinc_oxide_fitted_current,
       zinc_oxide_fitted_current_and_derivative

struct SaturationPreprocessResult
    input_kind::Symbol
    current_a::Vector{Float64}
    flux_wb::Vector{Float64}
    reconstructed_rms_current_a::Vector{Float64}
    maximum_relative_reconstruction_error::Float64
    monotone::Bool
    physical_checks_passed::Bool
end

struct HysteresisLoopPreprocessResult
    material_profile::Symbol
    level::Int
    saturation_current_a::Float64
    saturation_flux_wb::Float64
    current_a::Vector{Float64}
    flux_wb::Vector{Float64}
    runtime_current_a::Vector{Float64}
    runtime_flux_wb::Vector{Float64}
    closed_loop_current_a::Vector{Float64}
    closed_loop_flux_wb::Vector{Float64}
    energy_loss_j::Float64
    closure_error::Float64
    monotone_ascending_branch::Bool
    physical_checks_passed::Bool
end

struct ZincOxidePowerSegment
    first_sample_index::Int
    last_sample_index::Int
    current_coefficient_a::Float64
    voltage_exponent::Float64
    minimum_voltage_pu::Float64
    maximum_voltage_pu::Float64
    maximum_relative_current_error::Float64
end

struct ZincOxideFitResult
    fit_mode::Symbol
    reference_voltage_v::Float64
    sample_current_a::Vector{Float64}
    sample_voltage_v::Vector{Float64}
    fitted_current_a::Vector{Float64}
    segments::Vector{ZincOxidePowerSegment}
    maximum_relative_current_error::Float64
    normalized_rms_current_error::Float64
    continuity_error_a::Float64
    positive::Bool
    continuous::Bool
    fit_checks_passed::Bool
end

function _checked_strictly_increasing_positive(
    values::AbstractVector{<:Real},
    label::AbstractString,
)
    checked = Float64.(values)
    !isempty(checked) ||
        throw(ArgumentError("$label must not be empty"))
    all(value -> isfinite(value) && value > 0.0, checked) ||
        throw(ArgumentError("$label entries must be finite and positive"))
    all(index -> checked[index] > checked[index - 1], 2:length(checked)) ||
        throw(ArgumentError("$label entries must be strictly increasing"))
    return checked
end

function _checked_saturation_result(
    input_kind::Symbol,
    current_a::Vector{Float64},
    flux_wb::Vector{Float64},
    reconstructed_rms_current_a::Vector{Float64},
    maximum_relative_error::Float64,
)
    length(current_a) == length(flux_wb) >= 2 ||
        throw(ArgumentError("saturation characteristic requires at least two aligned points"))
    all(isfinite, current_a) && all(isfinite, flux_wb) ||
        throw(ArgumentError("saturation characteristic must be finite"))
    current_a[1] == 0.0 && flux_wb[1] == 0.0 ||
        throw(ArgumentError("saturation characteristic must begin at the origin"))
    monotone =
        all(index -> current_a[index] > current_a[index - 1], 2:length(current_a)) &&
        all(index -> flux_wb[index] > flux_wb[index - 1], 2:length(flux_wb))
    monotone ||
        throw(ArgumentError("saturation current and flux must increase strictly"))
    physical = isfinite(maximum_relative_error) && maximum_relative_error >= 0.0
    return SaturationPreprocessResult(
        input_kind,
        current_a,
        flux_wb,
        reconstructed_rms_current_a,
        maximum_relative_error,
        monotone,
        physical,
    )
end

"""
    saturation_from_incremental_inductance(current_a, incremental_inductance_h)

Integrate an incremental-inductance curve with the trapezoidal rule into the
current-versus-flux characteristic consumed by the existing piecewise
nonlinear-inductor runtime. A missing zero-current point is supplied using the
first measured incremental inductance.
"""
function saturation_from_incremental_inductance(
    current_a::AbstractVector{<:Real},
    incremental_inductance_h::AbstractVector{<:Real},
)
    length(current_a) == length(incremental_inductance_h) ||
        throw(ArgumentError("incremental saturation current and inductance counts must match"))
    !isempty(current_a) ||
        throw(ArgumentError("incremental saturation samples must not be empty"))
    currents = Float64.(current_a)
    inductances = Float64.(incremental_inductance_h)
    all(isfinite, currents) &&
        all(value -> isfinite(value) && value > 0.0, inductances) ||
        throw(ArgumentError("incremental saturation samples must be finite with positive inductance"))
    all(value -> value >= 0.0, currents) ||
        throw(ArgumentError("incremental saturation currents must be nonnegative"))
    all(index -> currents[index] > currents[index - 1], 2:length(currents)) ||
        throw(ArgumentError("incremental saturation currents must increase strictly"))
    if first(currents) > 0.0
        pushfirst!(currents, 0.0)
        pushfirst!(inductances, first(inductances))
    end
    length(currents) >= 2 ||
        throw(ArgumentError("incremental saturation curve requires a nonzero-current sample"))
    fluxes = zeros(Float64, length(currents))
    for index in 2:length(currents)
        fluxes[index] =
            fluxes[index - 1] +
            0.5 * (inductances[index - 1] + inductances[index]) *
            (currents[index] - currents[index - 1])
    end
    return _checked_saturation_result(
        :incremental_inductance,
        currents,
        fluxes,
        Float64[],
        0.0,
    )
end

function _piecewise_saturation_current(
    flux_wb::Float64,
    flux_knots::AbstractVector{Float64},
    current_knots::AbstractVector{Float64},
)
    flux_wb <= 0.0 && return 0.0
    flux_wb >= last(flux_knots) && return last(current_knots)
    right = searchsortedfirst(flux_knots, flux_wb)
    left = right - 1
    fraction =
        (flux_wb - flux_knots[left]) /
        (flux_knots[right] - flux_knots[left])
    return current_knots[left] +
           fraction * (current_knots[right] - current_knots[left])
end

function _saturation_rms_current(
    peak_flux_wb::Float64,
    flux_knots::AbstractVector{Float64},
    current_knots::AbstractVector{Float64},
    quadrature_steps::Int,
)
    square_sum = 0.0
    final_square = 0.0
    for step in 1:quadrature_steps
        angle = (pi / 2.0) * step / quadrature_steps
        current = _piecewise_saturation_current(
            peak_flux_wb * sin(angle),
            flux_knots,
            current_knots,
        )
        final_square = current * current
        square_sum += final_square
    end
    return sqrt(max(square_sum - 0.5 * final_square, 0.0) / quadrature_steps)
end

function _saturation_endpoint_current(
    target_rms_current_a::Float64,
    peak_flux_wb::Float64,
    flux_knots::Vector{Float64},
    current_knots::Vector{Float64},
    quadrature_steps::Int,
)
    lower = last(current_knots)
    trial_flux = [flux_knots; peak_flux_wb]
    function residual(endpoint_current::Float64)
        trial_current = [current_knots; endpoint_current]
        return _saturation_rms_current(
            peak_flux_wb,
            trial_flux,
            trial_current,
            quadrature_steps,
        ) - target_rms_current_a
    end
    lower_residual = residual(lower)
    lower_residual <= 1.0e-12 * max(target_rms_current_a, 1.0) ||
        throw(ArgumentError(
            "RMS saturation samples are inconsistent with a monotone peak-current curve",
        ))
    upper = max(sqrt(2.0) * target_rms_current_a, 2.0 * lower, 1.0)
    upper_residual = residual(upper)
    expansions = 0
    while upper_residual < 0.0 && expansions < 80
        upper *= 2.0
        upper_residual = residual(upper)
        expansions += 1
    end
    upper_residual >= 0.0 ||
        throw(ArgumentError("RMS saturation endpoint current could not be bracketed"))
    for _ in 1:100
        midpoint = 0.5 * (lower + upper)
        midpoint_residual = residual(midpoint)
        if midpoint_residual > 0.0
            upper = midpoint
        else
            lower = midpoint
        end
    end
    return 0.5 * (lower + upper)
end

"""
    saturation_from_rms_excitation(rms_current_a, rms_voltage_v; frequency_hz)

Convert a sinusoidal RMS excitation curve to the peak current-versus-flux
curve used by the transient nonlinear-inductor owner. The inverse conversion
is recomputed numerically and retained as an independent acceptance check.
"""
function saturation_from_rms_excitation(
    rms_current_a::AbstractVector{<:Real},
    rms_voltage_v::AbstractVector{<:Real};
    frequency_hz::Real,
    quadrature_steps::Integer = 90,
)
    length(rms_current_a) == length(rms_voltage_v) ||
        throw(ArgumentError("RMS saturation current and voltage counts must match"))
    currents_rms = _checked_strictly_increasing_positive(
        rms_current_a,
        "RMS saturation currents",
    )
    voltages_rms = _checked_strictly_increasing_positive(
        rms_voltage_v,
        "RMS saturation voltages",
    )
    length(currents_rms) >= 2 ||
        throw(ArgumentError("RMS saturation conversion requires at least two samples"))
    frequency = Float64(frequency_hz)
    isfinite(frequency) && frequency > 0.0 ||
        throw(ArgumentError("RMS saturation frequency_hz must be finite and positive"))
    steps = Int(quadrature_steps)
    steps >= 16 ||
        throw(ArgumentError("RMS saturation quadrature_steps must be at least 16"))
    peak_fluxes = [0.0; sqrt(2.0) .* voltages_rms ./ (2.0 * pi * frequency)]
    peak_currents = Float64[0.0, sqrt(2.0) * first(currents_rms)]
    for sample_index in 2:length(currents_rms)
        endpoint = _saturation_endpoint_current(
            currents_rms[sample_index],
            peak_fluxes[sample_index + 1],
            peak_fluxes[1:sample_index],
            peak_currents,
            steps,
        )
        endpoint > last(peak_currents) ||
            throw(ArgumentError("derived saturation peak current must increase strictly"))
        push!(peak_currents, endpoint)
    end
    reconstructed = Float64[
        _saturation_rms_current(
            peak_fluxes[index + 1],
            peak_fluxes,
            peak_currents,
            steps,
        )
        for index in eachindex(currents_rms)
    ]
    relative_errors = abs.(reconstructed .- currents_rms) ./ currents_rms
    return _checked_saturation_result(
        :sinusoidal_rms_excitation,
        peak_currents,
        peak_fluxes,
        reconstructed,
        maximum(relative_errors; init = 0.0),
    )
end

function saturation_runtime_table(result::SaturationPreprocessResult)
    result.physical_checks_passed ||
        throw(ArgumentError("saturation preprocessing must pass physical checks"))
    return (
        currents_a = copy(result.current_a),
        fluxes_wb = copy(result.flux_wb),
        table_start_index = 1,
        table_end_index = length(result.current_a),
        initial_segment = 1,
        physical_checks_passed = result.monotone,
    )
end

const _NORMALIZED_HYSTERESIS_FLUX = Float64[
    -15.6, 14.6, 17.0, 17.1,
    -16.5, -15.63, -14.4, -11.2, 9.7, 12.8, 14.7, 16.0, 17.0, 17.1,
    -16.6, -16.4, -15.9, -15.4, -14.2, -12.0, 8.6, 12.3, 13.8, 15.0,
    15.8, 16.4, 17.0, 17.1,
    -16.6, -16.5, -16.15, -15.8, -15.5, -14.9, -14.2, -13.0, -11.0,
    -8.0, 5.35, 7.4, 10.0, 12.0, 13.0, 14.0, 14.9, 15.6, 16.1, 16.6,
    17.0, 17.1,
]

const _NORMALIZED_HYSTERESIS_CURRENT = Float64[
    0.04, 0.17, 1.60, 2.20,
    -0.40, -0.05, 0.03, 0.07, 0.135, 0.21, 0.36, 0.665, 1.60, 2.20,
    -0.60, -0.30, -0.10, -0.02, 0.035, 0.066, 0.12, 0.19, 0.27, 0.40,
    0.59, 0.92, 1.60, 2.20,
    -0.60, -0.40, -0.18, -0.08, -0.03, 0.01, 0.035, 0.058, 0.07, 0.08,
    0.10, 0.11, 0.14, 0.18, 0.218, 0.285, 0.39, 0.535, 0.70, 1.00,
    1.60, 2.20,
]

const _NORMALIZED_HYSTERESIS_RANGES = (1:4, 5:14, 15:28, 29:50)

function _hysteresis_runtime_origin(
    current_a::Vector{Float64},
    flux_wb::Vector{Float64},
)
    zero_index = findfirst(==(0.0), current_a)
    zero_index !== nothing && return copy(current_a), copy(flux_wb)
    right = findfirst(>(0.0), current_a)
    if right === nothing
        throw(ArgumentError("normalized hysteresis branch must reach positive current"))
    elseif right == 1
        current_slope =
            (flux_wb[2] - flux_wb[1]) /
            (current_a[2] - current_a[1])
        zero_flux = flux_wb[1] - current_slope * current_a[1]
        return [0.0; current_a], [zero_flux; flux_wb]
    end
    left = right - 1
    fraction = -current_a[left] / (current_a[right] - current_a[left])
    zero_flux = flux_wb[left] +
                fraction * (flux_wb[right] - flux_wb[left])
    return (
        [current_a[1:left]; 0.0; current_a[right:end]],
        [flux_wb[1:left]; zero_flux; flux_wb[right:end]],
    )
end

function _closed_hysteresis_energy(
    current_a::Vector{Float64},
    flux_wb::Vector{Float64},
)
    closed_current = [current_a; -current_a; first(current_a)]
    closed_flux = [flux_wb; -flux_wb; first(flux_wb)]
    energy = 0.0
    for index in 2:length(closed_current)
        energy += 0.5 * (closed_current[index - 1] + closed_current[index]) *
                  (closed_flux[index] - closed_flux[index - 1])
    end
    closure_error = hypot(
        last(closed_current) - first(closed_current),
        last(closed_flux) - first(closed_flux),
    )
    return closed_current, closed_flux, abs(energy), closure_error
end

"""
    normalized_hysteresis_loop(level, saturation_current_a, saturation_flux_wb)

Scale one of the four established normalized material-loop resolutions to a
physical saturation point. The returned ascending points feed the existing
type-96 runtime; an exact zero-current interpolation and a closed symmetric
loop are retained for runtime construction and energy-loss checks.
"""
function normalized_hysteresis_loop(
    level::Integer,
    saturation_current_a::Real,
    saturation_flux_wb::Real,
)
    selected_level = Int(level)
    1 <= selected_level <= length(_NORMALIZED_HYSTERESIS_RANGES) ||
        throw(ArgumentError("normalized hysteresis level must be 1 through 4"))
    saturation_current = Float64(saturation_current_a)
    saturation_flux = Float64(saturation_flux_wb)
    all(value -> isfinite(value) && value > 0.0, (saturation_current, saturation_flux)) ||
        throw(ArgumentError("hysteresis saturation current and flux must be finite and positive"))
    indices = _NORMALIZED_HYSTERESIS_RANGES[selected_level]
    normalized_current = _NORMALIZED_HYSTERESIS_CURRENT[indices]
    normalized_flux = _NORMALIZED_HYSTERESIS_FLUX[indices]
    current_scale = saturation_current / normalized_current[end - 1]
    flux_scale = saturation_flux / normalized_flux[end - 1]
    current = current_scale .* normalized_current
    flux = flux_scale .* normalized_flux
    monotone =
        all(index -> current[index] > current[index - 1], 2:length(current)) &&
        all(index -> flux[index] > flux[index - 1], 2:length(flux))
    monotone ||
        throw(ArgumentError("normalized hysteresis ascending branch is not monotone"))
    runtime_current, runtime_flux =
        _hysteresis_runtime_origin(current, flux)
    closed_current, closed_flux, energy_loss, closure_error =
        _closed_hysteresis_energy(current, flux)
    physical =
        closure_error <= 64.0 * eps(Float64) &&
        energy_loss > 0.0 &&
        monotone
    return HysteresisLoopPreprocessResult(
        :normalized_electrical_steel,
        selected_level,
        saturation_current,
        saturation_flux,
        current,
        flux,
        runtime_current,
        runtime_flux,
        closed_current,
        closed_flux,
        energy_loss,
        closure_error,
        monotone,
        physical,
    )
end

function _zno_segment_cost(
    log_voltage::Vector{Float64},
    log_current::Vector{Float64},
    first_index::Int,
    last_index::Int,
)
    first_index < last_index ||
        throw(ArgumentError("a ZnO fit segment requires two distinct samples"))
    exponent =
        (log_current[last_index] - log_current[first_index]) /
        (log_voltage[last_index] - log_voltage[first_index])
    exponent > 0.0 ||
        return Inf
    intercept =
        log_current[first_index] - exponent * log_voltage[first_index]
    cost = 0.0
    for index in first_index:last_index
        residual =
            intercept + exponent * log_voltage[index] - log_current[index]
        cost += residual * residual
    end
    return cost
end

function _zno_optimal_breakpoints(
    log_voltage::Vector{Float64},
    log_current::Vector{Float64},
    segment_count::Int,
)
    sample_count = length(log_voltage)
    1 <= segment_count <= sample_count - 1 ||
        throw(ArgumentError("ZnO segment count must be between one and sample_count - 1"))
    costs = fill(Inf, segment_count, sample_count)
    predecessor = zeros(Int, segment_count, sample_count)
    for last_index in 2:sample_count
        costs[1, last_index] =
            _zno_segment_cost(log_voltage, log_current, 1, last_index)
    end
    for segment in 2:segment_count
        for last_index in (segment + 1):sample_count
            for split_index in segment:last_index - 1
                prior = costs[segment - 1, split_index]
                isfinite(prior) || continue
                candidate =
                    prior +
                    _zno_segment_cost(
                        log_voltage,
                        log_current,
                        split_index,
                        last_index,
                    )
                if candidate < costs[segment, last_index]
                    costs[segment, last_index] = candidate
                    predecessor[segment, last_index] = split_index
                end
            end
        end
    end
    isfinite(costs[segment_count, sample_count]) ||
        throw(ArgumentError("ZnO samples do not admit positive power-law segments"))
    breakpoints = zeros(Int, segment_count + 1)
    breakpoints[end] = sample_count
    last_index = sample_count
    for segment in segment_count:-1:2
        split_index = predecessor[segment, last_index]
        split_index > 0 || error("unreachable ZnO breakpoint reconstruction failure")
        breakpoints[segment] = split_index
        last_index = split_index
    end
    breakpoints[1] = 1
    return breakpoints
end

function _zno_fit_for_breakpoints(
    currents::Vector{Float64},
    voltages::Vector{Float64},
    reference_voltage::Float64,
    breakpoints::Vector{Int},
    fit_mode::Symbol,
)
    normalized_voltage = voltages ./ reference_voltage
    log_voltage = log.(normalized_voltage)
    log_current = log.(currents)
    segments = ZincOxidePowerSegment[]
    fitted = zeros(Float64, length(currents))
    continuity_error = 0.0
    for segment_index in 1:(length(breakpoints) - 1)
        first_index = breakpoints[segment_index]
        last_index = breakpoints[segment_index + 1]
        exponent =
            (log_current[last_index] - log_current[first_index]) /
            (log_voltage[last_index] - log_voltage[first_index])
        exponent > 0.0 ||
            throw(ArgumentError("ZnO voltage exponent must be positive"))
        coefficient =
            exp(log_current[first_index] - exponent * log_voltage[first_index])
        local_errors = Float64[]
        for index in first_index:last_index
            prediction = coefficient * normalized_voltage[index]^exponent
            if fitted[index] != 0.0
                continuity_error =
                    max(continuity_error, abs(prediction - fitted[index]))
            end
            fitted[index] = prediction
            push!(local_errors, abs(prediction - currents[index]) / currents[index])
        end
        push!(
            segments,
            ZincOxidePowerSegment(
                first_index,
                last_index,
                coefficient,
                exponent,
                normalized_voltage[first_index],
                normalized_voltage[last_index],
                maximum(local_errors; init = 0.0),
            ),
        )
    end
    relative_errors = abs.(fitted .- currents) ./ currents
    maximum_error = maximum(relative_errors; init = 0.0)
    normalized_rms_error =
        sqrt(sum(abs2, relative_errors) / length(relative_errors))
    continuity_scale = maximum(currents; init = 1.0)
    continuous = continuity_error <= 64.0 * eps(Float64) * continuity_scale
    positive = all(
        segment -> segment.current_coefficient_a > 0.0 &&
                   segment.voltage_exponent > 0.0,
        segments,
    )
    return ZincOxideFitResult(
        fit_mode,
        reference_voltage,
        currents,
        voltages,
        fitted,
        segments,
        maximum_error,
        normalized_rms_error,
        continuity_error,
        positive,
        continuous,
        positive && continuous,
    )
end

"""
    zinc_oxide_piecewise_fit(current_a, voltage_v; ...)

Fit a continuous piecewise power law `I = c(V/Vref)^p`. Supplying
`segment_count` selects a fixed optimal segmentation. Omitting it chooses the
smallest segment count satisfying `relative_error_tolerance`.
"""
function zinc_oxide_piecewise_fit(
    current_a::AbstractVector{<:Real},
    voltage_v::AbstractVector{<:Real};
    reference_voltage_v::Union{Nothing,Real} = nothing,
    segment_count::Union{Nothing,Integer} = nothing,
    relative_error_tolerance::Real = 0.05,
    maximum_segments::Integer = 18,
)
    length(current_a) == length(voltage_v) ||
        throw(ArgumentError("ZnO current and voltage sample counts must match"))
    currents = _checked_strictly_increasing_positive(current_a, "ZnO currents")
    voltages = _checked_strictly_increasing_positive(voltage_v, "ZnO voltages")
    length(currents) >= 2 ||
        throw(ArgumentError("ZnO fitting requires at least two samples"))
    reference_voltage = reference_voltage_v === nothing ?
        last(voltages) : Float64(reference_voltage_v)
    isfinite(reference_voltage) && reference_voltage > 0.0 ||
        throw(ArgumentError("ZnO reference_voltage_v must be finite and positive"))
    tolerance = Float64(relative_error_tolerance)
    isfinite(tolerance) && tolerance > 0.0 ||
        throw(ArgumentError("ZnO relative_error_tolerance must be finite and positive"))
    maximum_count = min(Int(maximum_segments), length(currents) - 1)
    maximum_count >= 1 ||
        throw(ArgumentError("ZnO maximum_segments must be positive"))
    log_voltage = log.(voltages ./ reference_voltage)
    log_current = log.(currents)
    if segment_count !== nothing
        count = Int(segment_count)
        count <= maximum_count ||
            throw(ArgumentError("ZnO fixed segment_count exceeds maximum_segments or sample support"))
        breakpoints =
            _zno_optimal_breakpoints(log_voltage, log_current, count)
        result = _zno_fit_for_breakpoints(
            currents,
            voltages,
            reference_voltage,
            breakpoints,
            :fixed_segment_count,
        )
        return ZincOxideFitResult(
            result.fit_mode,
            result.reference_voltage_v,
            result.sample_current_a,
            result.sample_voltage_v,
            result.fitted_current_a,
            result.segments,
            result.maximum_relative_current_error,
            result.normalized_rms_current_error,
            result.continuity_error_a,
            result.positive,
            result.continuous,
            result.fit_checks_passed &&
                result.maximum_relative_current_error <= tolerance,
        )
    end
    retained = nothing
    for count in 1:maximum_count
        breakpoints =
            _zno_optimal_breakpoints(log_voltage, log_current, count)
        candidate = _zno_fit_for_breakpoints(
            currents,
            voltages,
            reference_voltage,
            breakpoints,
            :automatic_error_bounded,
        )
        retained = candidate
        candidate.maximum_relative_current_error <= tolerance && break
    end
    result = something(retained)
    return ZincOxideFitResult(
        result.fit_mode,
        result.reference_voltage_v,
        result.sample_current_a,
        result.sample_voltage_v,
        result.fitted_current_a,
        result.segments,
        result.maximum_relative_current_error,
        result.normalized_rms_current_error,
        result.continuity_error_a,
        result.positive,
        result.continuous,
        result.fit_checks_passed &&
            result.maximum_relative_current_error <= tolerance,
    )
end

function zinc_oxide_runtime_table(result::ZincOxideFitResult)
    result.fit_checks_passed ||
        throw(ArgumentError("ZnO fit must pass its requested error and continuity checks"))
    first_segment = first(result.segments)
    first_voltage = first_segment.minimum_voltage_pu
    initial_conductance =
        first_segment.current_coefficient_a *
        first_voltage^first_segment.voltage_exponent /
        (first_voltage * result.reference_voltage_v)
    return (
        reference_voltage_v = result.reference_voltage_v,
        cchar = [initial_conductance; getfield.(result.segments, :current_coefficient_a)],
        gslope = [0.0; getfield.(result.segments, :voltage_exponent)],
        vchar = [0.0; getfield.(result.segments, :minimum_voltage_pu)],
        table_start_index = 1,
        table_end_index = length(result.segments) + 1,
        initial_table_index = length(result.segments) + 1,
        physical_checks_passed =
            result.positive && result.continuous && initial_conductance > 0.0,
    )
end

function zinc_oxide_fitted_current_and_derivative(
    result::ZincOxideFitResult,
    voltage_v::Real,
)
    voltage = Float64(voltage_v)
    isfinite(voltage) ||
        throw(ArgumentError("ZnO evaluation voltage must be finite"))
    table = zinc_oxide_runtime_table(result)
    current, derivative, search_steps = _over16_zno_current_and_derivative(
        voltage,
        table.reference_voltage_v,
        table.initial_table_index,
        table.table_start_index,
        table.table_end_index,
        table.cchar,
        table.gslope,
        table.vchar,
    )
    return (
        voltage_v = voltage,
        current_a = current,
        derivative_s = derivative,
        search_steps,
        active_runtime_owner = :over16_simultaneous_zno_current_and_derivative,
    )
end

function zinc_oxide_fitted_current(
    result::ZincOxideFitResult,
    voltage_v::Real,
)
    return zinc_oxide_fitted_current_and_derivative(result, voltage_v).current_a
end

export SemlyenRationalTerm,
       SemlyenModeParameters,
       SemlyenFrequencyDependentLine,
       SemlyenFrequencyScanFitResult,
       PoleResidueTransfer,
       PoleResidueFitResult,
       PoleResidueReductionResult,
       RationalLineModeConversion,
       RationalLineFrequencyFitResult,
       semlyen_rational_terms,
       pole_residue_transfer_value,
       pole_residue_transfer_fit,
       pole_residue_transfer_for_timestep,
       rational_frequency_dependent_mode_parameters,
       semlyen_frequency_dependent_line,
       semlyen_frequency_dependent_line_from_scan_fit,
       semlyen_frequency_dependent_line_from_rational_fit,
       semlyen_line_physical_checks,
       semlyen_line_steady_state_terminal_admittance,
       initialize_semlyen_line_steady_state!

struct SemlyenRationalTerm
    pole::ComplexF64
    residue::ComplexF64
    conjugate_pair::Bool
end

function SemlyenRationalTerm(pole, residue; conjugate_pair::Bool=false)
    pole_value = ComplexF64(pole)
    residue_value = ComplexF64(residue)
    all(isfinite, (real(pole_value), imag(pole_value), real(residue_value), imag(residue_value))) ||
        throw(ArgumentError("Semlyen pole and residue must be finite"))
    real(pole_value) > 0.0 ||
        throw(ArgumentError("Semlyen poles must have positive real parts"))
    !conjugate_pair && (imag(pole_value) != 0.0 || imag(residue_value) != 0.0) &&
        throw(ArgumentError("a complex Semlyen term must represent a conjugate pair"))
    return SemlyenRationalTerm(pole_value, residue_value, conjugate_pair)
end

"""Convert OVER2 `(kind, pole, residue)` triples to physical rational terms."""
function semlyen_rational_terms(triples::AbstractVector)
    terms = SemlyenRationalTerm[]
    index = 1
    while index <= length(triples)
        triple = triples[index]
        kind = Float64(triple[1])
        pole = Float64(triple[2])
        residue = Float64(triple[3])
        if kind == 0.0
            push!(terms, SemlyenRationalTerm(pole, residue))
            index += 1
        elseif kind > 0.0
            index < length(triples) ||
                throw(ArgumentError("Semlyen complex term is missing its quadrature row"))
            quadrature = triples[index + 1]
            Float64(quadrature[1]) < 0.0 ||
                throw(ArgumentError("Semlyen complex term quadrature marker must be negative"))
            push!(
                terms,
                SemlyenRationalTerm(
                    complex(pole, Float64(quadrature[2])),
                    complex(residue, Float64(quadrature[3]));
                    conjugate_pair = true,
                ),
            )
            index += 2
        else
            throw(ArgumentError("Semlyen quadrature term must follow a positive marker"))
        end
    end
    return terms
end

"""Convert a two-exponential step-response fit to continuous Semlyen terms."""
function semlyen_rational_terms(
    fit::LineStepResponseExponentialFitResult;
    response_kind::Symbol = :propagation,
)
    response_kind in (:propagation, :characteristic_admittance) ||
        throw(ArgumentError("unsupported Semlyen fitted response kind $response_kind"))
    fit.fit_executed || throw(ArgumentError("Semlyen conversion requires an executed fit"))
    time_constants = (fit.first_time_constant_s, fit.second_time_constant_s)
    amplitudes = (fit.first_amplitude, fit.second_amplitude)
    sign = response_kind == :propagation ? 1.0 : -1.0
    return SemlyenRationalTerm[
        SemlyenRationalTerm(inv(time_constant), sign * amplitude)
        for (time_constant, amplitude) in zip(time_constants, amplitudes)
    ]
end

struct SemlyenModeParameters
    characteristic_admittance_s::Float64
    travel_time_s::Float64
    phasor_series_impedance::ComplexF64
    phasor_characteristic_admittance::ComplexF64
    phasor_frequency_hz::Float64
    propagation_terms::Vector{SemlyenRationalTerm}
    admittance_terms::Vector{SemlyenRationalTerm}
end

"""A stable real-coefficient transfer `d + sum(r[k] / (s + p[k]))`."""
struct PoleResidueTransfer
    direct_term::Float64
    residues::Vector{Float64}
    poles::Vector{Float64}
end

function PoleResidueTransfer(
    direct_term::Real,
    residues::AbstractVector{<:Real},
    poles::AbstractVector{<:Real},
)
    direct = Float64(direct_term)
    residue_values = Float64.(residues)
    pole_values = Float64.(poles)
    length(residue_values) == length(pole_values) ||
        throw(ArgumentError("pole-residue transfer residue and pole counts must match"))
    isfinite(direct) && all(isfinite, residue_values) && all(isfinite, pole_values) ||
        throw(ArgumentError("pole-residue transfer values must be finite"))
    all(>(0.0), pole_values) ||
        throw(ArgumentError("pole-residue transfer poles must be positive in the s+p convention"))
    length(unique(pole_values)) == length(pole_values) ||
        throw(ArgumentError("pole-residue transfer poles must be distinct"))
    return PoleResidueTransfer(direct, residue_values, pole_values)
end

function pole_residue_transfer_value(response::PoleResidueTransfer, s::Number)
    frequency = ComplexF64(s)
    value = complex(response.direct_term)
    for index in eachindex(response.residues, response.poles)
        value += response.residues[index] / (frequency + response.poles[index])
    end
    return value
end

struct PoleResidueFitResult
    response::PoleResidueTransfer
    frequencies_hz::Vector{Float64}
    sample_values::Vector{ComplexF64}
    fitted_values::Vector{ComplexF64}
    maximum_absolute_error::Float64
    relative_maximum_absolute_error::Float64
    normalized_root_mean_square_error::Float64
    pole_refinement_sweeps::Int
    maximum_frequency_response_gain::Float64
    passivity_projection_scale::Float64
    stable_poles::Bool
    passivity_checks_passed::Bool
    fit_checks_passed::Bool
end

function _pole_residue_weighted_system(
    frequencies_hz::Vector{Float64},
    sample_values::Vector{ComplexF64},
    poles::Vector{Float64},
    fixed_direct_term::Union{Nothing,Float64},
)
    sample_scale = maximum(abs, sample_values; init = 0.0)
    weight_floor = max(sample_scale * 1.0e-6, eps(Float64))
    direct_unknown = fixed_direct_term === nothing
    column_count = length(poles) + (direct_unknown ? 1 : 0)
    matrix = zeros(Float64, 2 * length(frequencies_hz), column_count)
    target = zeros(Float64, 2 * length(frequencies_hz))
    for sample_index in eachindex(frequencies_hz)
        frequency = 2.0im * pi * frequencies_hz[sample_index]
        sample = sample_values[sample_index] -
            (direct_unknown ? 0.0 : something(fixed_direct_term))
        weight = inv(max(abs(sample_values[sample_index]), weight_floor))
        real_row = 2 * sample_index - 1
        imaginary_row = real_row + 1
        column = 1
        if direct_unknown
            matrix[real_row, column] = weight
            column += 1
        end
        for pole in poles
            basis = pole / (frequency + pole)
            matrix[real_row, column] = weight * real(basis)
            matrix[imaginary_row, column] = weight * imag(basis)
            column += 1
        end
        target[real_row] = weight * real(sample)
        target[imaginary_row] = weight * imag(sample)
    end
    return matrix, target
end

function _pole_residue_nonnegative_least_squares(
    matrix::Matrix{Float64},
    target::Vector{Float64};
    maximum_sweeps::Int = 20_000,
)
    column_count = size(matrix, 2)
    column_norms = vec(sum(abs2, matrix; dims = 1))
    all(>(0.0), column_norms) ||
        throw(ArgumentError("pole-residue fit basis contains a zero column"))
    values = max.(matrix \ target, 0.0)
    residual = target - matrix * values
    for sweep in 1:maximum_sweeps
        maximum_change = 0.0
        for column in 1:column_count
            old_value = values[column]
            numerator =
                dot(@view(matrix[:, column]), residual) +
                column_norms[column] * old_value
            new_value = max(0.0, numerator / column_norms[column])
            change = new_value - old_value
            change == 0.0 && continue
            residual .-= change .* @view(matrix[:, column])
            values[column] = new_value
            maximum_change = max(maximum_change, abs(change))
        end
        scale = max(maximum(abs, values; init = 0.0), 1.0)
        maximum_change <= 1.0e-12 * scale && return values, sweep
    end
    return values, maximum_sweeps
end

function _pole_residue_linear_fit(
    frequencies_hz::Vector{Float64},
    sample_values::Vector{ComplexF64},
    poles::Vector{Float64},
    fixed_direct_term::Union{Nothing,Float64},
    enforce_nonnegative_coefficients::Bool,
)
    matrix, target = _pole_residue_weighted_system(
        frequencies_hz,
        sample_values,
        poles,
        fixed_direct_term,
    )
    coefficients, sweeps = if enforce_nonnegative_coefficients
        _pole_residue_nonnegative_least_squares(matrix, target)
    else
        matrix \ target, 1
    end
    if fixed_direct_term === nothing
        direct_term = coefficients[1]
        normalized_residues = coefficients[2:end]
    else
        direct_term = fixed_direct_term
        normalized_residues = coefficients
    end
    residues = normalized_residues .* poles
    return PoleResidueTransfer(direct_term, residues, poles), sweeps
end

function _pole_residue_fit_objective(
    response::PoleResidueTransfer,
    frequencies_hz::Vector{Float64},
    sample_values::Vector{ComplexF64},
)
    fitted = ComplexF64[
        pole_residue_transfer_value(response, 2.0im * pi * frequency)
        for frequency in frequencies_hz
    ]
    sample_scale = maximum(abs, sample_values; init = 0.0)
    floor_scale = max(sample_scale * 1.0e-6, eps(Float64))
    normalized_errors = Float64[
        abs(fitted[index] - sample_values[index]) /
        max(abs(sample_values[index]), floor_scale)
        for index in eachindex(sample_values)
    ]
    return maximum(normalized_errors; init = 0.0)
end

function _pole_residue_passivity_grid(frequencies_hz::Vector{Float64})
    positive = filter(>(0.0), frequencies_hz)
    grid = copy(frequencies_hz)
    for index in 1:(length(positive) - 1)
        lower = positive[index]
        upper = positive[index + 1]
        for fraction in 1:7
            push!(
                grid,
                exp(
                    log(lower) +
                    (fraction / 8.0) * (log(upper) - log(lower)),
                ),
            )
        end
    end
    return sort!(unique(grid))
end

function _pole_residue_global_rate_grid(
    response::PoleResidueTransfer,
    frequencies_hz::Vector{Float64},
)
    sample_rates = 2.0 .* pi .* filter(>(0.0), frequencies_hz)
    anchors = sort!(unique(vcat(sample_rates, response.poles)))
    isempty(anchors) &&
        throw(ArgumentError("pole-residue global response grid requires a positive rate"))
    lower = first(anchors) / 1024.0
    upper = last(anchors) * 1024.0
    knots = sort!(unique(vcat(lower, anchors, upper)))
    rates = Float64[0.0]
    for fraction in 1:16
        push!(rates, lower * fraction / 16.0)
    end
    for index in 1:(length(knots) - 1)
        first_rate = knots[index]
        last_rate = knots[index + 1]
        for fraction in 0:16
            push!(
                rates,
                exp(
                    log(first_rate) +
                    (fraction / 16.0) * (log(last_rate) - log(first_rate)),
                ),
            )
        end
    end
    return sort!(unique(rates))
end

function _pole_residue_maximum_gain(
    response::PoleResidueTransfer,
    frequencies_hz::Vector{Float64},
)
    return maximum(
        rate -> abs(pole_residue_transfer_value(response, complex(0.0, rate))),
        _pole_residue_global_rate_grid(response, frequencies_hz);
        init = abs(response.direct_term),
    )
end

function _pole_residue_passivity_projection(
    response::PoleResidueTransfer,
    frequencies_hz::Vector{Float64},
)
    maximum_gain = _pole_residue_maximum_gain(response, frequencies_hz)
    maximum_gain > 1.0 || return response, maximum_gain, 1.0
    scale = prevfloat(1.0) / maximum_gain
    projected = PoleResidueTransfer(
        response.direct_term * scale,
        response.residues .* scale,
        response.poles,
    )
    return projected, maximum_gain, scale
end

"""
    pole_residue_transfer_fit(frequencies_hz, sample_values; order, ...)

Fit a stable real-coefficient rational response directly to complex-frequency
samples. Stable logarithmic poles are refined in the physical frequency band;
the dimensionless `r/p` coefficients and optional direct term are recomputed
at every sweep in a conditioned `p/(s+p)` basis.
Characteristic-impedance fits use a nonnegative Stieltjes basis, which makes
the fitted impedance positive real and gives its reciprocal stable,
real-interlacing poles.
"""
function pole_residue_transfer_fit(
    frequencies_hz::AbstractVector,
    sample_values::AbstractVector;
    order::Integer,
    fixed_direct_term::Union{Nothing,Real} = nothing,
    response_kind::Symbol = :generic,
    relative_tolerance::Real = 0.10,
    passivity_tolerance::Real = 1.0e-8,
    maximum_refinement_sweeps::Integer = 6,
)
    frequencies = Float64.(frequencies_hz)
    samples = ComplexF64.(sample_values)
    length(frequencies) == length(samples) && !isempty(frequencies) ||
        throw(ArgumentError("pole-residue fit frequencies and samples must be nonempty and aligned"))
    all(value -> isfinite(value) && value >= 0.0, frequencies) &&
        all(isfinite, real.(samples)) && all(isfinite, imag.(samples)) ||
        throw(ArgumentError("pole-residue fit samples must be finite at nonnegative frequencies"))
    issorted(frequencies) && all(diff(frequencies) .> 0.0) ||
        throw(ArgumentError("pole-residue fit frequencies must be strictly increasing"))
    response_kind in (:generic, :characteristic_impedance, :propagation) ||
        throw(ArgumentError("unsupported pole-residue response kind $response_kind"))
    term_count = Int(order)
    term_count > 0 || throw(ArgumentError("pole-residue fit order must be positive"))
    fixed_direct = fixed_direct_term === nothing ?
        nothing : Float64(fixed_direct_term)
    fixed_direct === nothing || isfinite(fixed_direct) ||
        throw(ArgumentError("pole-residue fixed direct term must be finite"))
    unknown_count = term_count + (fixed_direct === nothing ? 1 : 0)
    2 * length(frequencies) >= unknown_count ||
        throw(ArgumentError("pole-residue fit has fewer real equations than coefficients"))
    tolerance = Float64(relative_tolerance)
    passivity_limit = Float64(passivity_tolerance)
    isfinite(tolerance) && tolerance > 0.0 &&
        isfinite(passivity_limit) && passivity_limit > 0.0 ||
        throw(ArgumentError("pole-residue fit tolerances must be finite and positive"))
    refinement_sweeps = Int(maximum_refinement_sweeps)
    refinement_sweeps >= 0 ||
        throw(ArgumentError("pole-residue refinement sweep count must be nonnegative"))
    positive_frequencies = filter(>(0.0), frequencies)
    isempty(positive_frequencies) &&
        throw(ArgumentError("pole-residue fit requires at least one positive frequency"))
    minimum_rate = 2.0 * pi * first(positive_frequencies)
    maximum_rate = 2.0 * pi * last(positive_frequencies)
    lower_pole = max(0.1 * minimum_rate, sqrt(eps(Float64)))
    upper_pole = max(10.0 * maximum_rate, 10.0 * lower_pole)
    initial_grid = range(
        log(lower_pole),
        log(upper_pole);
        length = term_count + 2,
    )
    poles = exp.(initial_grid[2:(end - 1)])
    enforce_nonnegative = response_kind == :characteristic_impedance
    response, _ = _pole_residue_linear_fit(
        frequencies,
        samples,
        poles,
        fixed_direct,
        enforce_nonnegative,
    )
    unprojected_maximum_gain = _pole_residue_maximum_gain(response, frequencies)
    passivity_projection_scale = 1.0
    if response_kind == :propagation
        response, unprojected_maximum_gain, passivity_projection_scale =
            _pole_residue_passivity_projection(response, frequencies)
    end
    objective = _pole_residue_fit_objective(response, frequencies, samples)
    completed_sweeps = 0
    for sweep in 1:refinement_sweeps
        changed = false
        for index in eachindex(poles)
            retained_pole = poles[index]
            retained_response = response
            retained_objective = objective
            lower_bound =
                index == firstindex(poles) ? lower_pole / 100.0 :
                poles[index - 1] * (1.0 + 1.0e-6)
            upper_bound =
                index == lastindex(poles) ? upper_pole * 100.0 :
                poles[index + 1] * (1.0 - 1.0e-6)
            for multiplier in (0.5, 0.75, 0.9, 0.95, 1.05, 1.1, 1.25, 2.0)
                candidate_pole = clamp(retained_pole * multiplier, lower_bound, upper_bound)
                candidate_pole == retained_pole && continue
                candidate_poles = copy(poles)
                candidate_poles[index] = candidate_pole
                candidate_response, _ = _pole_residue_linear_fit(
                    frequencies,
                    samples,
                    candidate_poles,
                    fixed_direct,
                    enforce_nonnegative,
                )
                candidate_maximum_gain =
                    _pole_residue_maximum_gain(candidate_response, frequencies)
                candidate_projection_scale = 1.0
                if response_kind == :propagation
                    candidate_response,
                    candidate_maximum_gain,
                    candidate_projection_scale =
                        _pole_residue_passivity_projection(
                            candidate_response,
                            frequencies,
                        )
                end
                candidate_objective = _pole_residue_fit_objective(
                    candidate_response,
                    frequencies,
                    samples,
                )
                if candidate_objective < retained_objective
                    retained_pole = candidate_pole
                    retained_response = candidate_response
                    retained_objective = candidate_objective
                    unprojected_maximum_gain = candidate_maximum_gain
                    passivity_projection_scale = candidate_projection_scale
                end
            end
            if retained_pole != poles[index]
                poles[index] = retained_pole
                response = retained_response
                objective = retained_objective
                changed = true
            end
        end
        completed_sweeps = sweep
        changed || break
    end
    fitted = ComplexF64[
        pole_residue_transfer_value(response, 2.0im * pi * frequency)
        for frequency in frequencies
    ]
    errors = abs.(fitted .- samples)
    maximum_error = maximum(errors; init = 0.0)
    sample_scale = maximum(abs, samples; init = 0.0)
    floor_scale = max(sample_scale * 1.0e-6, eps(Float64))
    relative_errors = Float64[
        errors[index] / max(abs(samples[index]), floor_scale)
        for index in eachindex(samples)
    ]
    relative_maximum_error = maximum(relative_errors; init = 0.0)
    normalized_rms_error =
        sqrt(sum(abs2, relative_errors) / length(relative_errors))
    stable_poles = all(>(0.0), response.poles)
    passivity_values = ComplexF64[
        pole_residue_transfer_value(response, complex(0.0, rate))
        for rate in _pole_residue_global_rate_grid(response, frequencies)
    ]
    maximum_gain = maximum(abs, passivity_values; init = abs(response.direct_term))
    passivity_checks = if response_kind == :characteristic_impedance
        response.direct_term > 0.0 &&
            all(value -> real(value) >= -passivity_limit, passivity_values)
    elseif response_kind == :propagation
        abs(response.direct_term) <= passivity_limit &&
            all(value -> abs(value) <= 1.0 + passivity_limit, passivity_values)
    else
        true
    end
    fit_checks = stable_poles && passivity_checks &&
        relative_maximum_error <= tolerance
    return PoleResidueFitResult(
        response,
        frequencies,
        samples,
        fitted,
        maximum_error,
        relative_maximum_error,
        normalized_rms_error,
        completed_sweeps,
        maximum_gain,
        passivity_projection_scale,
        stable_poles,
        passivity_checks,
        fit_checks,
    )
end

struct PoleResidueReductionResult
    response::PoleResidueTransfer
    original_term_count::Int
    effective_threshold::Float64
    reduction_applied::Bool
    inverse_first_moment_error::Float64
    inverse_second_moment_error::Float64
end

"""
    pole_residue_transfer_for_timestep(response, timestep_s, minimum_pole_step)

Apply the DISTR2 high-pole tail equivalence. The merged term preserves
`sum(r/p)` and `sum(r/p^2)`, so the zero-frequency value and first derivative
of the reduced tail are unchanged.
"""
function pole_residue_transfer_for_timestep(
    response::PoleResidueTransfer,
    timestep_s::Real,
    minimum_pole_step::Real,
)
    timestep = Float64(timestep_s)
    threshold = Float64(minimum_pole_step)
    isfinite(threshold) && threshold >= 0.0 ||
        throw(ArgumentError("pole-residue reduction threshold must be nonnegative"))
    term_count = length(response.poles)
    if threshold == 0.0
        return PoleResidueReductionResult(response, term_count, threshold, false, 0.0, 0.0)
    end
    isfinite(timestep) && timestep > 0.0 ||
        throw(ArgumentError("pole-residue reduction timestep must be positive"))
    if term_count < 3
        return PoleResidueReductionResult(response, term_count, threshold, false, 0.0, 0.0)
    end
    issorted(response.poles) || throw(ArgumentError(
        "pole-residue reduction requires poles in ascending order",
    ))
    effective_threshold = max(threshold, 1.0)
    while true
        first_tail = findfirst(pole -> pole * timestep >= effective_threshold, response.poles)
        if first_tail === nothing
            return PoleResidueReductionResult(
                response,
                term_count,
                effective_threshold,
                false,
                0.0,
                0.0,
            )
        end
        first_tail == 1 && (first_tail = 2)
        if first_tail == term_count
            return PoleResidueReductionResult(
                response,
                term_count,
                effective_threshold,
                false,
                0.0,
                0.0,
            )
        end
        first_moment = 0.0
        second_moment = 0.0
        for index in first_tail:term_count
            residue = response.residues[index]
            pole = response.poles[index]
            first_moment += residue / pole
            second_moment += residue / (pole * pole)
        end
        abs(first_moment) > eps(Float64) && abs(second_moment) > eps(Float64) ||
            throw(ArgumentError("pole-residue tail moments are singular"))
        equivalent_pole = first_moment / second_moment
        equivalent_residue = first_moment * equivalent_pole
        if equivalent_pole > 0.0 && equivalent_pole * timestep >= effective_threshold
            reduced = PoleResidueTransfer(
                response.direct_term,
                vcat(response.residues[1:(first_tail - 1)], equivalent_residue),
                vcat(response.poles[1:(first_tail - 1)], equivalent_pole),
            )
            reduced_first = equivalent_residue / equivalent_pole
            reduced_second = equivalent_residue / (equivalent_pole * equivalent_pole)
            return PoleResidueReductionResult(
                reduced,
                term_count,
                effective_threshold,
                true,
                abs(first_moment - reduced_first),
                abs(second_moment - reduced_second),
            )
        end
        effective_threshold += 1.0
    end
end

struct RationalLineModeConversion
    parameters::SemlyenModeParameters
    characteristic_impedance::PoleResidueTransfer
    propagation_without_delay::PoleResidueTransfer
    characteristic_admittance_terms::Vector{SemlyenRationalTerm}
    reciprocal_reconstruction_max_abs_error::Float64
    physical_checks_passed::Bool
end

function _ascending_polynomial_product(left::Vector{Float64}, right::Vector{Float64})
    result = zeros(length(left) + length(right) - 1)
    for left_index in eachindex(left), right_index in eachindex(right)
        result[left_index + right_index - 1] += left[left_index] * right[right_index]
    end
    return result
end

function _pole_denominator(poles::Vector{Float64}; omitted_index::Int=0)
    coefficients = [1.0]
    for (index, pole) in pairs(poles)
        index == omitted_index && continue
        coefficients = _ascending_polynomial_product(coefficients, [pole, 1.0])
    end
    return coefficients
end

function _ascending_polynomial_value(coefficients::Vector{Float64}, value::ComplexF64)
    result = 0.0 + 0.0im
    for coefficient in Iterators.reverse(coefficients)
        result = muladd(result, value, coefficient)
    end
    return result
end

function _ascending_polynomial_roots(coefficients::Vector{Float64})
    degree = length(coefficients) - 1
    degree == 0 && return ComplexF64[]
    leading = coefficients[end]
    abs(leading) > 0.0 || throw(ArgumentError("rational transfer numerator is singular"))
    companion = zeros(ComplexF64, degree, degree)
    for index in 2:degree
        companion[index, index - 1] = 1.0
    end
    companion[:, end] .= -coefficients[1:degree] ./ leading
    return eigvals(companion)
end

function _reciprocal_rational_terms(
    impedance::PoleResidueTransfer;
    root_tolerance::Float64=1.0e-8,
)
    impedance.direct_term > 0.0 || throw(ArgumentError(
        "characteristic-impedance infinity value must be positive",
    ))
    isempty(impedance.poles) && return SemlyenRationalTerm[]
    frequency_scale = exp(sum(log, impedance.poles) / length(impedance.poles))
    normalized_poles = impedance.poles ./ frequency_scale
    normalized_residues = impedance.residues ./ frequency_scale
    denominator = _pole_denominator(normalized_poles)
    numerator = impedance.direct_term .* denominator
    for index in eachindex(normalized_poles)
        quotient = _pole_denominator(normalized_poles; omitted_index = index)
        numerator[1:length(quotient)] .+= normalized_residues[index] .* quotient
    end
    numerator_scale = maximum(abs, numerator; init = 1.0)
    while length(numerator) > 1 && abs(numerator[end]) <= eps(Float64) * numerator_scale
        pop!(numerator)
    end
    roots = _ascending_polynomial_roots(numerator)
    derivative = [index * numerator[index + 1] for index in 1:(length(numerator) - 1)]
    terms = SemlyenRationalTerm[]
    ordered_roots = sort(roots; by = value -> (real(value), imag(value)))
    consumed = falses(length(ordered_roots))
    for root_index in eachindex(ordered_roots)
        consumed[root_index] && continue
        root = ordered_roots[root_index]
        normalized_pole = -root
        scale = max(1.0, abs(root))
        real(normalized_pole) > root_tolerance * scale || throw(ArgumentError(
            "characteristic-admittance reciprocal has a non-stable pole",
        ))
        conventional_residue =
            _ascending_polynomial_value(denominator, root) /
            _ascending_polynomial_value(derivative, root)
        if abs(imag(root)) <= root_tolerance * scale
            residue_scale = max(1.0, abs(conventional_residue))
            abs(imag(conventional_residue)) <= root_tolerance * residue_scale ||
                throw(ArgumentError("characteristic-admittance residue is not real"))
            real_pole = real(normalized_pole) * frequency_scale
            push!(
                terms,
                SemlyenRationalTerm(
                    real_pole,
                    real(conventional_residue) / real(normalized_pole),
                ),
            )
            consumed[root_index] = true
            continue
        end
        conjugate_index = findfirst(
            candidate_index ->
                !consumed[candidate_index] &&
                candidate_index != root_index &&
                abs(
                    ordered_roots[candidate_index] - conj(root),
                ) <= root_tolerance * scale,
            eachindex(ordered_roots),
        )
        conjugate_index === nothing && throw(ArgumentError(
            "characteristic-admittance reciprocal has an unpaired complex pole",
        ))
        conjugate_root = ordered_roots[conjugate_index]
        conjugate_residue =
            _ascending_polynomial_value(denominator, conjugate_root) /
            _ascending_polynomial_value(derivative, conjugate_root)
        residue_scale = max(1.0, abs(conventional_residue))
        abs(conjugate_residue - conj(conventional_residue)) <=
            root_tolerance * residue_scale || throw(ArgumentError(
            "characteristic-admittance reciprocal has inconsistent conjugate residues",
        ))
        retained_root, retained_residue =
            imag(normalized_pole) >= 0.0 ?
            (root, conventional_residue) :
            (conjugate_root, conjugate_residue)
        retained_pole = -retained_root
        push!(
            terms,
            SemlyenRationalTerm(
                retained_pole * frequency_scale,
                retained_residue / retained_pole;
                conjugate_pair = true,
            ),
        )
        consumed[root_index] = true
        consumed[conjugate_index] = true
    end
    return terms
end

function _rational_terms_value(terms::AbstractVector{SemlyenRationalTerm}, s::ComplexF64)
    value = 0.0 + 0.0im
    for term in terms
        value += term.residue * term.pole / (s + term.pole)
        term.conjugate_pair || continue
        conjugate_pole = conj(term.pole)
        conjugate_residue = conj(term.residue)
        value += conjugate_residue * conjugate_pole / (s + conjugate_pole)
    end
    return value
end

"""
    rational_frequency_dependent_mode_parameters(...)

Convert a characteristic-impedance and propagation pole-residue basis into
the characteristic-admittance realization used by the executable fitted-line
solver. Inputs use SI ohms, seconds, hertz, and poles in radians per second.
"""
function rational_frequency_dependent_mode_parameters(
    characteristic_impedance::PoleResidueTransfer,
    propagation_without_delay::PoleResidueTransfer,
    travel_time_s::Real,
    phasor_frequency_hz::Real;
    reciprocal_tolerance::Real=5.0e-10,
)
    propagation_without_delay.direct_term == 0.0 || throw(ArgumentError(
        "propagation pole-residue response must have zero direct term",
    ))
    delay = Float64(travel_time_s)
    frequency_hz = Float64(phasor_frequency_hz)
    isfinite(delay) && delay > 0.0 ||
        throw(ArgumentError("rational line travel time must be finite and positive"))
    isfinite(frequency_hz) && frequency_hz > 0.0 ||
        throw(ArgumentError("rational line phasor frequency must be finite and positive"))
    tolerance = Float64(reciprocal_tolerance)
    isfinite(tolerance) && tolerance > 0.0 ||
        throw(ArgumentError("rational line reciprocal tolerance must be positive"))

    admittance_terms = _reciprocal_rational_terms(characteristic_impedance)
    propagation_terms = SemlyenRationalTerm[
        SemlyenRationalTerm(pole, residue / pole)
        for (residue, pole) in zip(
            propagation_without_delay.residues,
            propagation_without_delay.poles,
        )
    ]
    characteristic_admittance_s = inv(characteristic_impedance.direct_term)
    angular_frequency = 2.0 * pi * frequency_hz
    phasor_s = ComplexF64(0.0, angular_frequency)
    impedance = pole_residue_transfer_value(characteristic_impedance, phasor_s)
    propagation = pole_residue_transfer_value(propagation_without_delay, phasor_s)
    abs(propagation) > eps(Float64) ||
        throw(ArgumentError("rational line propagation response is zero at the phasor frequency"))
    abs(propagation) <= 1.0 + 64.0 * eps(Float64) ||
        throw(ArgumentError("rational line propagation response is active rather than passive"))
    propagation_constant = -log(propagation) + complex(0.0, angular_frequency * delay)
    real(propagation_constant) >= -64.0 * eps(Float64) ||
        throw(ArgumentError("rational line propagation constant has negative attenuation"))
    series_impedance = impedance * propagation_constant
    shunt_admittance = propagation_constant / impedance
    real(impedance) > 0.0 ||
        throw(ArgumentError("rational line characteristic impedance is not passive"))
    real(inv(impedance)) > 0.0 ||
        throw(ArgumentError("rational line characteristic admittance is not passive"))
    real(series_impedance) >= -tolerance * max(abs(series_impedance), 1.0) ||
        throw(ArgumentError("rational line modal series impedance is not passive"))
    real(shunt_admittance) >= -tolerance * max(abs(shunt_admittance), 1.0) ||
        throw(ArgumentError("rational line modal shunt admittance is not passive"))

    check_rates = unique(sort(vcat(
        0.0,
        angular_frequency,
        characteristic_impedance.poles,
        real.(getfield.(admittance_terms, :pole)),
    )))
    reconstruction_error = 0.0
    for rate in check_rates
        point = complex(0.0, rate)
        impedance_point = pole_residue_transfer_value(characteristic_impedance, point)
        exact = inv(impedance_point)
        real(impedance_point) >= -tolerance * characteristic_impedance.direct_term ||
            throw(ArgumentError("rational line characteristic impedance fails passivity"))
        real(exact) >= -tolerance * characteristic_admittance_s ||
            throw(ArgumentError("rational line characteristic admittance fails passivity"))
        reconstructed = characteristic_admittance_s +
            _rational_terms_value(admittance_terms, point)
        reconstruction_error = max(reconstruction_error, abs(exact - reconstructed))
    end
    propagation_check_rates = unique(sort(vcat(check_rates, propagation_without_delay.poles)))
    all(
        rate -> abs(pole_residue_transfer_value(
            propagation_without_delay,
            complex(0.0, rate),
        )) <= 1.0 + tolerance,
        propagation_check_rates,
    ) || throw(ArgumentError("rational line propagation response fails passive gain limits"))
    scale = max(characteristic_admittance_s, eps(Float64))
    reconstruction_error <= tolerance * scale || throw(ArgumentError(
        "characteristic-admittance reciprocal exceeds its reconstruction tolerance",
    ))
    parameters = SemlyenModeParameters(
        characteristic_admittance_s,
        delay,
        series_impedance,
        shunt_admittance,
        frequency_hz,
        propagation_terms,
        admittance_terms,
    )
    return RationalLineModeConversion(
        parameters,
        characteristic_impedance,
        propagation_without_delay,
        copy(admittance_terms),
        reconstruction_error,
        true,
    )
end

function SemlyenModeParameters(
    characteristic_admittance_s::Real,
    travel_time_s::Real,
    phasor_series_impedance,
    phasor_characteristic_admittance,
    phasor_frequency_hz::Real,
    propagation_terms::AbstractVector{SemlyenRationalTerm},
    admittance_terms::AbstractVector{SemlyenRationalTerm},
)
    admittance = Float64(characteristic_admittance_s)
    delay = Float64(travel_time_s)
    frequency = Float64(phasor_frequency_hz)
    isfinite(admittance) && admittance > 0.0 ||
        throw(ArgumentError("Semlyen characteristic admittance must be finite and positive"))
    isfinite(delay) && delay > 0.0 ||
        throw(ArgumentError("Semlyen travel time must be finite and positive"))
    isfinite(frequency) && frequency > 0.0 ||
        throw(ArgumentError("Semlyen phasor frequency must be finite and positive"))
    series = ComplexF64(phasor_series_impedance)
    phasor_admittance = ComplexF64(phasor_characteristic_admittance)
    all(isfinite, (real(series), imag(series), real(phasor_admittance), imag(phasor_admittance))) ||
        throw(ArgumentError("Semlyen phasor parameters must be finite"))
    return SemlyenModeParameters(
        admittance,
        delay,
        series,
        phasor_admittance,
        frequency,
        collect(propagation_terms),
        collect(admittance_terms),
    )
end

mutable struct SemlyenModeRuntimeState
    parameters::SemlyenModeParameters
    propagation_coefficients::Vector{SemlyenLineExponentialConvolutionCoefficients}
    admittance_coefficients::Vector{SemlyenLineExponentialConvolutionCoefficients}
    propagation_from_state::Vector{ComplexF64}
    propagation_to_state::Vector{ComplexF64}
    admittance_from_state::Vector{ComplexF64}
    admittance_to_state::Vector{ComplexF64}
    outgoing_from_history::Vector{Float64}
    outgoing_to_history::Vector{Float64}
    write_index::Int
    delay_steps::Int
    delay_fraction::Float64
    history_current_from::Float64
    history_current_to::Float64
end

mutable struct SemlyenFrequencyDependentLine <: EMTElement
    from_nodes::Vector{Int}
    to_nodes::Vector{Int}
    modes::Vector{SemlyenModeRuntimeState}
    voltage_modal_to_phase::Matrix{ComplexF64}
    current_modal_to_phase::Matrix{ComplexF64}
    runtime_current_modal_to_phase::Matrix{Float64}
    phase_to_modal_voltage::Matrix{Float64}
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
    modal_history_current_from::Vector{Float64}
    modal_history_current_to::Vector{Float64}
    timestep_s::Float64
    update_count::Int
end

struct SemlyenFrequencyScanFitResult
    line::SemlyenFrequencyDependentLine
    sample_frequencies_hz::Vector{Float64}
    phasor_frequency_hz::Float64
    line_length::Float64
    propagation_fits::Vector{LineStepResponseExponentialFitResult}
    admittance_fits::Vector{LineStepResponseExponentialFitResult}
    propagation_max_abs_error::Vector{Float64}
    propagation_relative_max_abs_error::Vector{Float64}
    admittance_max_abs_error_s::Vector{Float64}
    admittance_relative_max_abs_error::Vector{Float64}
    maximum_tdfit_normalized_square_error::Float64
    stable_poles::Bool
    fit_checks_passed::Bool
    physical_checks_passed::Bool
end

struct RationalLineFrequencyFitResult
    line::SemlyenFrequencyDependentLine
    sample_frequencies_hz::Vector{Float64}
    phasor_frequency_hz::Float64
    requested_travel_times_s::Vector{Float64}
    travel_times_s::Vector{Float64}
    travel_time_refinement_evaluations::Vector{Int}
    maximum_travel_time_adjustment_s::Float64
    characteristic_impedance_orders::Vector{Int}
    characteristic_impedance_fits::Vector{PoleResidueFitResult}
    propagation_without_delay_fits::Vector{PoleResidueFitResult}
    mode_conversions::Vector{RationalLineModeConversion}
    maximum_relative_fit_error::Float64
    stable_poles::Bool
    passivity_checks_passed::Bool
    fit_checks_passed::Bool
    physical_checks_passed::Bool
end

function _pole_residue_propagation_fit(
    frequencies::Vector{Float64},
    propagation_samples::Vector{ComplexF64},
    requested_travel_time_s::Float64,
    timestep_s::Float64;
    order::Int,
    relative_tolerance::Float64,
    passivity_tolerance::Float64,
    maximum_refinement_sweeps::Int,
    refine_travel_time::Bool,
)
    maximum_frequency = last(frequencies)
    maximum_frequency > 0.0 ||
        throw(ArgumentError("propagation fitting requires a positive maximum frequency"))
    evaluation_count = 0
    function fit_at_delay(delay_s::Float64)
        evaluation_count += 1
        without_delay = ComplexF64[
            propagation_samples[index] *
            cis(2.0 * pi * frequencies[index] * delay_s)
            for index in eachindex(frequencies)
        ]
        return pole_residue_transfer_fit(
            frequencies,
            without_delay;
            order,
            fixed_direct_term = 0.0,
            response_kind = :propagation,
            relative_tolerance,
            passivity_tolerance,
            maximum_refinement_sweeps,
        )
    end
    best_delay = requested_travel_time_s
    best_fit = fit_at_delay(best_delay)
    refine_travel_time || return best_fit, best_delay, evaluation_count
    half_span = min(
        1.0e-3 / maximum_frequency,
        max(1.0e-3 * requested_travel_time_s, 32.0 * eps(requested_travel_time_s)),
    )
    absolute_lower = max(timestep_s, requested_travel_time_s - half_span)
    absolute_upper = requested_travel_time_s + half_span
    for _ in 1:3
        lower = max(absolute_lower, best_delay - half_span)
        upper = min(absolute_upper, best_delay + half_span)
        candidates = range(lower, upper; length = 17)
        for candidate in candidates
            delay = Float64(candidate)
            delay == best_delay && continue
            fit = fit_at_delay(delay)
            fit_key = (
                !fit.passivity_checks_passed,
                fit.relative_maximum_absolute_error,
                fit.normalized_root_mean_square_error,
            )
            best_key = (
                !best_fit.passivity_checks_passed,
                best_fit.relative_maximum_absolute_error,
                best_fit.normalized_root_mean_square_error,
            )
            if fit_key < best_key
                best_fit = fit
                best_delay = delay
            end
        end
        half_span /= 8.0
    end
    return best_fit, best_delay, evaluation_count
end

_semlyen_term_multiplier(term::SemlyenRationalTerm) = term.conjugate_pair ? 2.0 : 1.0

function _semlyen_mode_state(parameters::SemlyenModeParameters, timestep_s::Float64)
    propagation_coefficients = [
        semlyen_line_exponential_convolution_coefficients(term.pole, term.residue, timestep_s)
        for term in parameters.propagation_terms
    ]
    admittance_coefficients = [
        semlyen_line_exponential_convolution_coefficients(term.pole, term.residue, timestep_s)
        for term in parameters.admittance_terms
    ]
    delay_ratio = parameters.travel_time_s / timestep_s
    delay_steps = floor(Int, delay_ratio)
    delay_steps >= 1 ||
        throw(ArgumentError("Semlyen travel time must be at least one timestep"))
    delay_fraction = delay_ratio - delay_steps
    history_length = delay_steps + 3
    return SemlyenModeRuntimeState(
        parameters,
        propagation_coefficients,
        admittance_coefficients,
        zeros(ComplexF64, length(propagation_coefficients)),
        zeros(ComplexF64, length(propagation_coefficients)),
        zeros(ComplexF64, length(admittance_coefficients)),
        zeros(ComplexF64, length(admittance_coefficients)),
        zeros(history_length),
        zeros(history_length),
        1,
        delay_steps,
        delay_fraction,
        0.0,
        0.0,
    )
end

function semlyen_line_physical_checks(
    modes::AbstractVector{SemlyenModeParameters},
    voltage_modal_to_phase::AbstractMatrix,
    current_modal_to_phase::AbstractMatrix;
    transform_tolerance::Real=2.5e-2,
    inactive_phase_indices::AbstractVector{<:Integer}=Int[],
)
    mode_count = length(modes)
    voltage_transform = Matrix{ComplexF64}(voltage_modal_to_phase)
    current_transform = Matrix{ComplexF64}(current_modal_to_phase)
    size(voltage_transform) == (mode_count, mode_count) ||
        throw(ArgumentError("Semlyen voltage transform size must match its mode count"))
    size(current_transform) == (mode_count, mode_count) ||
        throw(ArgumentError("Semlyen current transform size must match its mode count"))
    all(isfinite, voltage_transform) && all(isfinite, current_transform) ||
        throw(ArgumentError("Semlyen transforms must be finite"))
    inverse_error = maximum(abs.(transpose(current_transform) * voltage_transform - I))
    tolerance = Float64(transform_tolerance)
    inactive = Int.(inactive_phase_indices)
    all(index -> 1 <= index <= mode_count, inactive) ||
        throw(ArgumentError("Semlyen inactive phase index is outside the transform"))
    (inverse_error <= tolerance || !isempty(inactive)) ||
        throw(ArgumentError("Semlyen voltage/current transforms are not power duals"))
    runtime_current_transform = real.(current_transform)
    phase_admittance = runtime_current_transform *
        Diagonal(getfield.(modes, :characteristic_admittance_s)) *
        transpose(runtime_current_transform)
    symmetry_error = maximum(abs.(phase_admittance - transpose(phase_admittance)))
    symmetry_scale = max(maximum(abs.(phase_admittance)), eps(Float64))
    symmetry_error / symmetry_scale <= 1.0e-3 ||
        throw(ArgumentError("Semlyen phase companion admittance must be reciprocal"))
    minimum_eigenvalue = minimum(eigvals(Symmetric(phase_admittance)))
    minimum_eigenvalue >= -64.0 * eps(Float64) ||
        throw(ArgumentError("Semlyen phase companion admittance must be passive"))
    return (
        transform_inverse_max_abs_error = inverse_error,
        phase_admittance_symmetry_max_abs_error = symmetry_error,
        phase_admittance_minimum_eigenvalue = minimum_eigenvalue,
        physical_checks_passed = true,
    )
end

function semlyen_frequency_dependent_line(
    from_nodes::AbstractVector{<:Integer},
    to_nodes::AbstractVector{<:Integer},
    modes::AbstractVector{SemlyenModeParameters},
    voltage_modal_to_phase::AbstractMatrix,
    current_modal_to_phase::AbstractMatrix,
    timestep_s::Real,
)
    mode_count = length(modes)
    mode_count > 0 || throw(ArgumentError("Semlyen line requires at least one mode"))
    length(from_nodes) == mode_count && length(to_nodes) == mode_count ||
        throw(ArgumentError("Semlyen terminal count must match its mode count"))
    from_indices = Int.(from_nodes)
    to_indices = Int.(to_nodes)
    all(>=(0), from_indices) && all(>=(0), to_indices) ||
        throw(ArgumentError("Semlyen terminal indices must be nonnegative"))
    timestep = Float64(timestep_s)
    isfinite(timestep) && timestep > 0.0 ||
        throw(ArgumentError("Semlyen timestep must be finite and positive"))
    inactive_phases = [
        index for index in eachindex(from_indices)
        if from_indices[index] == 0 && to_indices[index] == 0
    ]
    checks = semlyen_line_physical_checks(
        modes,
        voltage_modal_to_phase,
        current_modal_to_phase;
        inactive_phase_indices = inactive_phases,
    )
    voltage_transform = Matrix{ComplexF64}(voltage_modal_to_phase)
    current_transform = Matrix{ComplexF64}(current_modal_to_phase)
    runtime_current_transform = real.(current_transform)
    phase_to_modal = transpose(runtime_current_transform)
    phase_admittance = runtime_current_transform *
        Diagonal(getfield.(modes, :characteristic_admittance_s)) *
        phase_to_modal
    checks.physical_checks_passed || error("unreachable Semlyen physical check failure")
    zeros_phase = zeros(mode_count)
    return SemlyenFrequencyDependentLine(
        from_indices,
        to_indices,
        [_semlyen_mode_state(mode, timestep) for mode in modes],
        voltage_transform,
        current_transform,
        Matrix(runtime_current_transform),
        Matrix(phase_to_modal),
        Matrix(phase_admittance),
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
        copy(zeros_phase),
        copy(zeros_phase),
        timestep,
        0,
    )
end

function _semlyen_fitted_transfer(
    terms::AbstractVector{SemlyenRationalTerm},
    angular_frequency_rad_s::Float64,
)
    frequency = ComplexF64(0.0, angular_frequency_rad_s)
    return sum(
        term.residue * term.pole / (term.pole + frequency)
        for term in terms;
        init = 0.0 + 0.0im,
    )
end

function _semlyen_fit_frequency_alignment(
    fit::LineStepResponseExponentialFitResult,
    frequencies_hz::AbstractVector{Float64},
    label::AbstractString,
)
    angular = fit.angular_frequencies_rad_s
    length(angular) == length(frequencies_hz) + 1 || throw(ArgumentError(
        "$label TDFIT frequency count must equal the scan count plus the zero-frequency row",
    ))
    angular[1] == 0.0 ||
        throw(ArgumentError("$label TDFIT frequency sequence must start at zero"))
    for index in eachindex(frequencies_hz)
        expected = 2.0 * pi * frequencies_hz[index]
        isapprox(angular[index + 1], expected; atol = 1.0e-10, rtol = 1.0e-10) ||
            throw(ArgumentError("$label TDFIT frequencies must match the generated scan"))
    end
    return nothing
end

function _semlyen_relative_fit_error(error::Float64, samples)
    scale = maximum(abs.(samples); init = 0.0)
    return scale == 0.0 ? error : error / scale
end

"""
    semlyen_frequency_dependent_line_from_scan_fit(...)

Build the executable fitted-line owner from frequency-domain modal samples and
the two-exponential fits produced by `line_step_response_exponential_fit`.
Propagation fits retain their delay and amplitudes. Characteristic-admittance
fits use the sign convention `Yc(s) = Y0 - sum(Ak*pk/(s+pk))`, matching the
coefficients consumed by the timestep convolution.
"""
function semlyen_frequency_dependent_line_from_scan_fit(
    from_nodes::AbstractVector{<:Integer},
    to_nodes::AbstractVector{<:Integer},
    sample_rows::AbstractVector,
    line_length::Real,
    voltage_modal_to_phase::AbstractMatrix,
    current_modal_to_phase::AbstractMatrix,
    timestep_s::Real,
    propagation_fits::AbstractVector{LineStepResponseExponentialFitResult},
    admittance_fits::AbstractVector{LineStepResponseExponentialFitResult},
    characteristic_admittance_s::AbstractVector{<:Real};
    phasor_frequency_hz::Real,
    propagation_relative_tolerance::Real = 0.10,
    admittance_relative_tolerance::Real = 0.10,
)
    line_length_value = _checked_line_positive_finite(line_length, "Semlyen line_length")
    timestep = _checked_line_positive_finite(timestep_s, "Semlyen timestep_s")
    frequencies, rows, mode_count = _checked_line_frequency_sample_rows(sample_rows)
    length(from_nodes) == mode_count && length(to_nodes) == mode_count ||
        throw(ArgumentError("Semlyen scan terminal count must match its modal samples"))
    length(propagation_fits) == mode_count && length(admittance_fits) == mode_count ||
        throw(ArgumentError("Semlyen scan fit count must match its modal samples"))
    length(characteristic_admittance_s) == mode_count ||
        throw(ArgumentError("Semlyen characteristic-admittance count must match its modal samples"))
    propagation_tolerance = _checked_line_positive_finite(
        propagation_relative_tolerance,
        "Semlyen propagation_relative_tolerance",
    )
    admittance_tolerance = _checked_line_positive_finite(
        admittance_relative_tolerance,
        "Semlyen admittance_relative_tolerance",
    )
    phasor_frequency = _checked_line_positive_finite(
        phasor_frequency_hz,
        "Semlyen phasor_frequency_hz",
    )
    phasor_fit = frequency_dependent_line_sample_fit(
        rows,
        phasor_frequency,
        line_length_value,
    )

    modes = SemlyenModeParameters[]
    propagation_abs_errors = Float64[]
    propagation_relative_errors = Float64[]
    admittance_abs_errors = Float64[]
    admittance_relative_errors = Float64[]
    stable_poles = true
    scan_physical_checks_passed = true
    maximum_tdfit_error = 0.0
    for mode in 1:mode_count
        propagation_fit = propagation_fits[mode]
        admittance_fit = admittance_fits[mode]
        _semlyen_fit_frequency_alignment(propagation_fit, frequencies, "propagation mode $mode")
        _semlyen_fit_frequency_alignment(admittance_fit, frequencies, "admittance mode $mode")
        propagation_fit.delay_s >= timestep || throw(ArgumentError(
            "Semlyen fitted propagation delay must be at least one timestep",
        ))
        abs(admittance_fit.delay_s) <= 1.0e-14 || throw(ArgumentError(
            "Semlyen characteristic-admittance fit must have zero delay",
        ))
        propagation_terms = semlyen_rational_terms(
            propagation_fit;
            response_kind = :propagation,
        )
        admittance_terms = semlyen_rational_terms(
            admittance_fit;
            response_kind = :characteristic_admittance,
        )
        stable_poles &= all(
            term -> real(term.pole) > 0.0,
            Iterators.flatten((propagation_terms, admittance_terms)),
        )
        initial_admittance = Float64(characteristic_admittance_s[mode])
        isfinite(initial_admittance) && initial_admittance > 0.0 ||
            throw(ArgumentError("Semlyen characteristic admittance must be finite and positive"))

        propagation_samples = ComplexF64[row[mode].propagation_factor for row in rows]
        propagation_consistency = ComplexF64[
            exp(-row[mode].propagation_constant * line_length_value)
            for row in rows
        ]
        consistency_scale = max(maximum(abs.(propagation_samples); init = 0.0), 1.0)
        maximum(abs.(propagation_samples .- propagation_consistency); init = 0.0) <=
            1.0e-10 * consistency_scale || throw(ArgumentError(
            "Semlyen propagation samples are inconsistent with propagation constants and line length",
        ))
        scan_physical_checks_passed &=
            all(value -> abs(value) <= 1.0 + 1.0e-10, propagation_samples)
        fitted_propagation = ComplexF64[
            cis(-2.0 * pi * frequency * propagation_fit.delay_s) *
            _semlyen_fitted_transfer(propagation_terms, 2.0 * pi * frequency)
            for frequency in frequencies
        ]
        admittance_samples = ComplexF64[inv(row[mode].characteristic_impedance) for row in rows]
        admittance_scale = max(maximum(abs.(admittance_samples); init = 0.0), 1.0)
        scan_physical_checks_passed &= all(
            value -> real(value) >= -1.0e-10 * admittance_scale,
            admittance_samples,
        )
        fitted_admittance = ComplexF64[
            initial_admittance +
            _semlyen_fitted_transfer(admittance_terms, 2.0 * pi * frequency)
            for frequency in frequencies
        ]
        propagation_error = maximum(abs.(propagation_samples .- fitted_propagation); init = 0.0)
        admittance_error = maximum(abs.(admittance_samples .- fitted_admittance); init = 0.0)
        push!(propagation_abs_errors, propagation_error)
        push!(
            propagation_relative_errors,
            _semlyen_relative_fit_error(propagation_error, propagation_samples),
        )
        push!(admittance_abs_errors, admittance_error)
        push!(
            admittance_relative_errors,
            _semlyen_relative_fit_error(admittance_error, admittance_samples),
        )
        maximum_tdfit_error = max(
            maximum_tdfit_error,
            propagation_fit.normalized_square_error,
            admittance_fit.normalized_square_error,
        )

        point = phasor_fit.frequency_points[mode]
        propagation = point.propagation_constant * line_length_value
        series_impedance = point.characteristic_impedance * propagation
        shunt_admittance = propagation / point.characteristic_impedance
        push!(
            modes,
            SemlyenModeParameters(
                initial_admittance,
                propagation_fit.delay_s,
                series_impedance,
                shunt_admittance,
                phasor_frequency,
                propagation_terms,
                admittance_terms,
            ),
        )
    end
    fit_checks_passed = stable_poles && scan_physical_checks_passed &&
        all(<=(propagation_tolerance), propagation_relative_errors) &&
        all(<=(admittance_tolerance), admittance_relative_errors)
    fit_checks_passed || throw(ArgumentError(
        "Semlyen frequency-scan fit exceeds its accepted relative error tolerance",
    ))
    physical_checks = semlyen_line_physical_checks(
        modes,
        voltage_modal_to_phase,
        current_modal_to_phase,
    )
    line = semlyen_frequency_dependent_line(
        from_nodes,
        to_nodes,
        modes,
        voltage_modal_to_phase,
        current_modal_to_phase,
        timestep,
    )
    return SemlyenFrequencyScanFitResult(
        line,
        collect(frequencies),
        phasor_frequency,
        line_length_value,
        collect(propagation_fits),
        collect(admittance_fits),
        propagation_abs_errors,
        propagation_relative_errors,
        admittance_abs_errors,
        admittance_relative_errors,
        maximum_tdfit_error,
        stable_poles,
        fit_checks_passed,
        physical_checks.physical_checks_passed && scan_physical_checks_passed,
    )
end

"""
    semlyen_frequency_dependent_line_from_rational_fit(...)

Fit general stable Marti pole-residue responses to modal complex-frequency
samples and construct the existing executable rational-line owner. The
propagation delay remains explicit; only the attenuation/dispersion response
is fitted, preserving the runtime history-ring mutation order.
"""
function semlyen_frequency_dependent_line_from_rational_fit(
    from_nodes::AbstractVector{<:Integer},
    to_nodes::AbstractVector{<:Integer},
    sample_rows::AbstractVector,
    travel_times_s::AbstractVector,
    voltage_modal_to_phase::AbstractMatrix,
    current_modal_to_phase::AbstractMatrix,
    timestep_s::Real;
    phasor_frequency_hz::Real,
    characteristic_impedance_order::Integer = 10,
    minimum_characteristic_impedance_order::Integer = 2,
    propagation_order::Integer = 6,
    relative_fit_tolerance::Real = 0.10,
    passivity_tolerance::Real = 1.0e-8,
    reciprocal_tolerance::Real = 5.0e-8,
    maximum_refinement_sweeps::Integer = 6,
    refine_travel_times::Bool = true,
)
    frequencies, rows, mode_count =
        _checked_line_frequency_sample_rows(sample_rows)
    length(from_nodes) == mode_count && length(to_nodes) == mode_count ||
        throw(ArgumentError("rational scan terminal count must match its modal samples"))
    travel_times = Float64.(travel_times_s)
    length(travel_times) == mode_count ||
        throw(ArgumentError("rational scan travel-time count must match its modes"))
    timestep = _checked_line_positive_finite(timestep_s, "rational scan timestep_s")
    all(value -> isfinite(value) && value >= timestep, travel_times) ||
        throw(ArgumentError("rational scan travel times must be finite and at least one timestep"))
    phasor_frequency = _checked_line_positive_finite(
        phasor_frequency_hz,
        "rational scan phasor_frequency_hz",
    )
    fit_tolerance = _checked_line_positive_finite(
        relative_fit_tolerance,
        "rational scan relative_fit_tolerance",
    )
    passivity_limit = _checked_line_positive_finite(
        passivity_tolerance,
        "rational scan passivity_tolerance",
    )
    refinement_sweeps = Int(maximum_refinement_sweeps)
    refinement_sweeps >= 0 ||
        throw(ArgumentError("rational scan refinement sweep count must be nonnegative"))
    requested_characteristic_order = Int(characteristic_impedance_order)
    minimum_characteristic_order = Int(minimum_characteristic_impedance_order)
    requested_characteristic_order >= minimum_characteristic_order >= 1 ||
        throw(ArgumentError(
            "rational scan characteristic-impedance order range must be positive and ordered",
        ))
    characteristic_fits = PoleResidueFitResult[]
    propagation_fits = PoleResidueFitResult[]
    conversions = RationalLineModeConversion[]
    characteristic_orders = Int[]
    fitted_travel_times = Float64[]
    travel_time_evaluations = Int[]
    maximum_relative_error = 0.0
    for mode in 1:mode_count
        characteristic_samples = ComplexF64[
            row[mode].characteristic_impedance for row in rows
        ]
        propagation_samples =
            ComplexF64[row[mode].propagation_factor for row in rows]
        propagation_fit, fitted_travel_time, evaluation_count =
            _pole_residue_propagation_fit(
            frequencies,
            propagation_samples,
            travel_times[mode],
            timestep;
            order = Int(propagation_order),
            relative_tolerance = fit_tolerance,
            passivity_tolerance = passivity_limit,
                maximum_refinement_sweeps = refinement_sweeps,
                refine_travel_time = refine_travel_times,
            )
        propagation_fit.fit_checks_passed || throw(ArgumentError(
            "mode $mode propagation rational fit exceeds its accepted tolerance or passive-gain boundary",
        ))
        characteristic_fit = nothing
        conversion = nothing
        selected_characteristic_order = 0
        conversion_error = nothing
        for candidate_order in
            requested_characteristic_order:-1:minimum_characteristic_order
            candidate_fit = pole_residue_transfer_fit(
                frequencies,
                characteristic_samples;
                order = candidate_order,
                response_kind = :characteristic_impedance,
                relative_tolerance = fit_tolerance,
                passivity_tolerance = passivity_limit,
                maximum_refinement_sweeps = refinement_sweeps,
            )
            candidate_fit.fit_checks_passed || continue
            candidate_conversion = try
                rational_frequency_dependent_mode_parameters(
                    candidate_fit.response,
                    propagation_fit.response,
                    fitted_travel_time,
                    phasor_frequency;
                    reciprocal_tolerance,
                )
            catch error
                error isa ArgumentError || rethrow()
                conversion_error = error
                nothing
            end
            candidate_conversion === nothing && continue
            characteristic_fit = candidate_fit
            conversion = candidate_conversion
            selected_characteristic_order = candidate_order
            break
        end
        characteristic_fit === nothing && throw(ArgumentError(
            "mode $mode has no accepted characteristic-impedance order in " *
            "$minimum_characteristic_order:$requested_characteristic_order" *
            (conversion_error === nothing ? "" :
             "; last conversion error: $(sprint(showerror, conversion_error))"),
        ))
        push!(characteristic_fits, characteristic_fit)
        push!(propagation_fits, propagation_fit)
        push!(conversions, conversion)
        push!(characteristic_orders, selected_characteristic_order)
        push!(fitted_travel_times, fitted_travel_time)
        push!(travel_time_evaluations, evaluation_count)
        maximum_relative_error = max(
            maximum_relative_error,
            characteristic_fit.relative_maximum_absolute_error,
            propagation_fit.relative_maximum_absolute_error,
        )
    end
    stable_poles = all(
        fit -> fit.stable_poles,
        Iterators.flatten((characteristic_fits, propagation_fits)),
    )
    passivity_checks = all(
        fit -> fit.passivity_checks_passed,
        Iterators.flatten((characteristic_fits, propagation_fits)),
    )
    fit_checks = stable_poles && passivity_checks &&
        maximum_relative_error <= fit_tolerance
    physical = semlyen_line_physical_checks(
        getfield.(conversions, :parameters),
        voltage_modal_to_phase,
        current_modal_to_phase,
    )
    line = semlyen_frequency_dependent_line(
        from_nodes,
        to_nodes,
        getfield.(conversions, :parameters),
        voltage_modal_to_phase,
        current_modal_to_phase,
        timestep,
    )
    return RationalLineFrequencyFitResult(
        line,
        collect(frequencies),
        phasor_frequency,
        travel_times,
        fitted_travel_times,
        travel_time_evaluations,
        maximum(abs.(fitted_travel_times .- travel_times); init = 0.0),
        characteristic_orders,
        characteristic_fits,
        propagation_fits,
        conversions,
        maximum_relative_error,
        stable_poles,
        passivity_checks,
        fit_checks,
        physical.physical_checks_passed,
    )
end

function _semlyen_modal_terminal_admittance(parameters::SemlyenModeParameters)
    series_impedance = parameters.phasor_series_impedance
    shunt_admittance = parameters.phasor_characteristic_admittance
    abs(series_impedance) > 0.0 && abs(shunt_admittance) > 0.0 ||
        throw(ArgumentError("Semlyen steady-state series impedance and shunt admittance must be nonzero"))
    propagation = sqrt(series_impedance * shunt_admittance)
    characteristic_impedance = sqrt(series_impedance / shunt_admittance)
    abs(sinh(propagation)) > eps(Float64) ||
        throw(ArgumentError("Semlyen steady-state propagation factor is singular"))
    self_admittance = coth(propagation) / characteristic_impedance
    mutual_admittance = -csch(propagation) / characteristic_impedance
    return self_admittance, mutual_admittance
end

function _semlyen_modal_pi_parameters(
    parameters::SemlyenModeParameters,
    frequency_hz::Float64,
)
    series_impedance, shunt_admittance = if isapprox(
        parameters.phasor_frequency_hz,
        frequency_hz;
        atol = 1.0e-9,
        rtol = 1.0e-9,
    )
        (
            parameters.phasor_series_impedance,
            parameters.phasor_characteristic_admittance,
        )
    else
        angular_frequency = 2.0 * pi * frequency_hz
        point = ComplexF64(0.0, angular_frequency)
        propagation_without_delay =
            _rational_terms_value(parameters.propagation_terms, point)
        abs(propagation_without_delay) > eps(Float64) ||
            throw(ArgumentError("Semlyen frequency response has zero propagation gain"))
        abs(propagation_without_delay) <= 1.0 + 1.0e-8 ||
            throw(ArgumentError("Semlyen frequency response has active propagation gain"))
        characteristic_admittance =
            parameters.characteristic_admittance_s +
            _rational_terms_value(parameters.admittance_terms, point)
        real(characteristic_admittance) > 0.0 ||
            throw(ArgumentError("Semlyen frequency response has non-passive characteristic admittance"))
        characteristic_impedance = inv(characteristic_admittance)
        propagation =
            -log(propagation_without_delay) +
            ComplexF64(0.0, angular_frequency * parameters.travel_time_s)
        real(propagation) >= -1.0e-8 ||
            throw(ArgumentError("Semlyen frequency response has negative attenuation"))
        (
            characteristic_impedance * propagation,
            propagation / characteristic_impedance,
        )
    end
    propagation = sqrt(series_impedance * shunt_admittance)
    characteristic_impedance = sqrt(series_impedance / shunt_admittance)
    series_branch_impedance = characteristic_impedance * sinh(propagation)
    abs(series_branch_impedance) > eps(Float64) ||
        throw(ArgumentError("Semlyen steady-state series branch is singular"))
    terminal_shunt_admittance = (cosh(propagation) - 1.0) /
        series_branch_impedance
    return series_branch_impedance, terminal_shunt_admittance
end

"""
    semlyen_line_steady_state_terminal_admittance(line, frequency_hz)

Return the complex phase-domain terminal admittance at the supplied phasor
frequency. The modal `Z` and `Y` data are converted to the same equivalent-pi
series impedance and terminal shunt admittance assembled by `OVER8`; voltage
and current transforms therefore retain their distinct physical roles.
"""
function semlyen_line_steady_state_terminal_admittance(
    line::SemlyenFrequencyDependentLine,
    frequency_hz::Real,
)
    frequency = Float64(frequency_hz)
    isfinite(frequency) && frequency > 0.0 ||
        throw(ArgumentError("Semlyen steady-state frequency must be finite and positive"))
    series_modal = ComplexF64[]
    shunt_modal = ComplexF64[]
    for state in line.modes
        series_impedance, shunt_admittance =
            _semlyen_modal_pi_parameters(state.parameters, frequency)
        push!(series_modal, series_impedance)
        push!(shunt_modal, shunt_admittance)
    end
    phase_count = length(line.modes)
    phase_series_impedance = line.voltage_modal_to_phase *
        Diagonal(series_modal) * transpose(line.voltage_modal_to_phase)
    phase_series_admittance = phase_series_impedance \
        Matrix{ComplexF64}(I, phase_count, phase_count)
    phase_shunt_admittance = line.current_modal_to_phase *
        Diagonal(shunt_modal) * transpose(line.current_modal_to_phase)
    self_phase = phase_series_admittance + phase_shunt_admittance
    mutual_phase = -phase_series_admittance
    return [self_phase mutual_phase; mutual_phase self_phase]
end

function _semlyen_seed_ring!(
    history::Vector{Float64},
    write_index::Int,
    phasor::ComplexF64,
    angular_frequency::Float64,
    timestep_s::Float64,
)
    fill!(history, 0.0)
    for lag in 1:(length(history) - 1)
        time_s = -lag * timestep_s
        history[mod1(write_index - lag, length(history))] =
            real(phasor * cis(angular_frequency * time_s))
    end
    return history
end

function _semlyen_correct_state_sum!(
    values::Vector{ComplexF64},
    terms::Vector{SemlyenRationalTerm},
    target::Float64,
    label::AbstractString,
)
    if isempty(values)
        abs(target) <= 1.0e-10 ||
            throw(ArgumentError("Semlyen $label cannot be represented without rational terms"))
        return values
    end
    # OVER13 labels 14920-14924 correct only a final real exponential.
    # A final conjugate-pair row takes the STYPE > 0 branch and retains both
    # quadrature states exactly as evaluated from the harmonic transfer.
    last(terms).conjugate_pair && return values
    correction = (target - _semlyen_state_sum(terms, values)) /
        _semlyen_term_multiplier(last(terms))
    values[end] += correction
    return values
end

"""
    initialize_semlyen_line_steady_state!(line, from_phasors, to_phasors, frequency_hz)

Initialize a fitted line from terminal voltage phasors. This owns the `OVER13`
traveling-wave ring, propagation and characteristic-admittance convolution
states, and stages the Norton history source for the first solve at `t = dt`.
Phasors use peak-value cosine convention and SI seconds/hertz.
"""
function initialize_semlyen_line_steady_state!(
    line::SemlyenFrequencyDependentLine,
    from_voltage_phasors::AbstractVector{<:Complex},
    to_voltage_phasors::AbstractVector{<:Complex},
    frequency_hz::Real,
)
    phase_count = length(line.from_nodes)
    length(from_voltage_phasors) == phase_count &&
        length(to_voltage_phasors) == phase_count ||
        throw(ArgumentError("Semlyen steady-state terminal phasor counts must match the line"))
    frequency = Float64(frequency_hz)
    terminal_admittance = semlyen_line_steady_state_terminal_admittance(line, frequency)
    phase_voltage = vcat(
        ComplexF64.(from_voltage_phasors),
        ComplexF64.(to_voltage_phasors),
    )
    phase_current = terminal_admittance * phase_voltage
    phase_to_modal = transpose(line.current_modal_to_phase)
    modal_voltage_from = phase_to_modal * phase_voltage[1:phase_count]
    modal_voltage_to = phase_to_modal * phase_voltage[(phase_count + 1):end]
    phase_to_modal_current = transpose(line.voltage_modal_to_phase)
    modal_current_from =
        phase_to_modal_current * phase_current[1:phase_count]
    modal_current_to =
        phase_to_modal_current * phase_current[(phase_count + 1):end]
    angular_frequency = 2.0 * pi * frequency

    for mode in eachindex(line.modes)
        state = line.modes[mode]
        state.write_index = 1
        parameters = state.parameters
        characteristic_impedance = sqrt(
            parameters.phasor_series_impedance /
            parameters.phasor_characteristic_admittance,
        )
        outgoing_from = 0.5 * (
            modal_voltage_from[mode] + characteristic_impedance * modal_current_from[mode]
        )
        outgoing_to = 0.5 * (
            modal_voltage_to[mode] + characteristic_impedance * modal_current_to[mode]
        )
        incoming_from = modal_voltage_from[mode] - outgoing_from
        incoming_to = modal_voltage_to[mode] - outgoing_to
        _semlyen_seed_ring!(
            state.outgoing_from_history,
            state.write_index,
            outgoing_from,
            angular_frequency,
            line.timestep_s,
        )
        _semlyen_seed_ring!(
            state.outgoing_to_history,
            state.write_index,
            outgoing_to,
            angular_frequency,
            line.timestep_s,
        )
        delayed_to = outgoing_to * cis(-angular_frequency * parameters.travel_time_s)
        delayed_from = outgoing_from * cis(-angular_frequency * parameters.travel_time_s)
        for index in eachindex(parameters.propagation_terms)
            term = parameters.propagation_terms[index]
            seeded = semlyen_line_harmonic_history_update(
                term.pole,
                term.residue,
                angular_frequency,
                delayed_to,
                delayed_from,
            )
            state.propagation_from_state[index] = seeded.from_history
            state.propagation_to_state[index] = seeded.to_history
        end
        _semlyen_correct_state_sum!(
            state.propagation_from_state,
            parameters.propagation_terms,
            real(incoming_from),
            "from-terminal propagation",
        )
        _semlyen_correct_state_sum!(
            state.propagation_to_state,
            parameters.propagation_terms,
            real(incoming_to),
            "to-terminal propagation",
        )

        effective_from = modal_voltage_from[mode] - 2.0 * incoming_from
        effective_to = modal_voltage_to[mode] - 2.0 * incoming_to
        for index in eachindex(parameters.admittance_terms)
            term = parameters.admittance_terms[index]
            seeded = semlyen_line_harmonic_history_update(
                term.pole,
                term.residue,
                angular_frequency,
                effective_from,
                effective_to,
            )
            state.admittance_from_state[index] = seeded.from_history
            state.admittance_to_state[index] = seeded.to_history
        end
        _semlyen_correct_state_sum!(
            state.admittance_from_state,
            parameters.admittance_terms,
            real(modal_current_from[mode]) -
            parameters.characteristic_admittance_s * real(effective_from),
            "from-terminal admittance",
        )
        _semlyen_correct_state_sum!(
            state.admittance_to_state,
            parameters.admittance_terms,
            real(modal_current_to[mode]) -
            parameters.characteristic_admittance_s * real(effective_to),
            "to-terminal admittance",
        )
        next_step_phase = cis(angular_frequency * line.timestep_s)
        state.history_current_from = real(
            (modal_current_from[mode] -
             parameters.characteristic_admittance_s * modal_voltage_from[mode]) *
            next_step_phase,
        )
        state.history_current_to = real(
            (modal_current_to[mode] -
             parameters.characteristic_admittance_s * modal_voltage_to[mode]) *
            next_step_phase,
        )
        line.modal_history_current_from[mode] = state.history_current_from
        line.modal_history_current_to[mode] = state.history_current_to
    end

    line.terminal_voltage_from .= real.(phase_voltage[1:phase_count])
    line.terminal_voltage_to .= real.(phase_voltage[(phase_count + 1):end])
    line.terminal_current_from .= real.(phase_current[1:phase_count])
    line.terminal_current_to .= real.(phase_current[(phase_count + 1):end])
    line.modal_voltage_from .= real.(modal_voltage_from)
    line.modal_voltage_to .= real.(modal_voltage_to)
    line.modal_current_from .= real.(modal_current_from)
    line.modal_current_to .= real.(modal_current_to)
    mul!(
        line.history_current_from,
        line.runtime_current_modal_to_phase,
        line.modal_history_current_from,
    )
    mul!(
        line.history_current_to,
        line.runtime_current_modal_to_phase,
        line.modal_history_current_to,
    )
    line.update_count = 0
    return line
end

function _semlyen_state_sum(terms, values)
    total = 0.0
    for (term, value) in zip(terms, values)
        total += _semlyen_term_multiplier(term) * real(value)
    end
    return total
end

function _semlyen_history_value(history::Vector{Float64}, newest::Int, lag::Int)
    return history[mod1(newest - lag, length(history))]
end

function _semlyen_propagation_inputs(state::SemlyenModeRuntimeState, history::Vector{Float64})
    lag = max(state.delay_steps - 1, 0)
    x0 = _semlyen_history_value(history, state.write_index, lag)
    x1 = _semlyen_history_value(history, state.write_index, lag + 1)
    x2 = _semlyen_history_value(history, state.write_index, lag + 2)
    return x0, x1, x2
end

function _semlyen_update_propagation_state!(
    values::Vector{ComplexF64},
    coefficients,
    fraction::Float64,
    inputs,
)
    x0, x1, x2 = inputs
    complement = 1.0 - fraction
    for index in eachindex(values, coefficients)
        coefficient = coefficients[index]
        values[index] = coefficient.decay * values[index] +
            (complement * coefficient.current_gain) * x0 +
            (fraction * coefficient.current_gain + complement * coefficient.delayed_gain) * x1 +
            (fraction * coefficient.delayed_gain) * x2
    end
    return values
end

function _semlyen_update_admittance_state!(
    values::Vector{ComplexF64},
    coefficients,
    voltage::Float64,
    propagated_next::Float64,
    effective_previous::Float64,
)
    for index in eachindex(values, coefficients)
        coefficient = coefficients[index]
        partial = values[index] + coefficient.current_gain * voltage
        values[index] = coefficient.decay * partial -
            2.0 * coefficient.current_gain * propagated_next +
            coefficient.delayed_gain * effective_previous
    end
    return values
end

function _semlyen_mode_update!(
    state::SemlyenModeRuntimeState,
    from_voltage::Float64,
    to_voltage::Float64,
)
    parameters = state.parameters
    propagated_from_previous =
        _semlyen_state_sum(parameters.propagation_terms, state.propagation_from_state)
    propagated_to_previous =
        _semlyen_state_sum(parameters.propagation_terms, state.propagation_to_state)
    state.outgoing_from_history[state.write_index] = from_voltage - propagated_from_previous
    state.outgoing_to_history[state.write_index] = to_voltage - propagated_to_previous
    _semlyen_update_propagation_state!(
        state.propagation_from_state,
        state.propagation_coefficients,
        state.delay_fraction,
        _semlyen_propagation_inputs(state, state.outgoing_to_history),
    )
    _semlyen_update_propagation_state!(
        state.propagation_to_state,
        state.propagation_coefficients,
        state.delay_fraction,
        _semlyen_propagation_inputs(state, state.outgoing_from_history),
    )
    propagated_from_next =
        _semlyen_state_sum(parameters.propagation_terms, state.propagation_from_state)
    propagated_to_next =
        _semlyen_state_sum(parameters.propagation_terms, state.propagation_to_state)
    effective_from_previous = from_voltage - 2.0 * propagated_from_previous
    effective_to_previous = to_voltage - 2.0 * propagated_to_previous
    _semlyen_update_admittance_state!(
        state.admittance_from_state,
        state.admittance_coefficients,
        from_voltage,
        propagated_from_next,
        effective_from_previous,
    )
    _semlyen_update_admittance_state!(
        state.admittance_to_state,
        state.admittance_coefficients,
        to_voltage,
        propagated_to_next,
        effective_to_previous,
    )
    state.history_current_from =
        -2.0 * parameters.characteristic_admittance_s * propagated_from_next +
        _semlyen_state_sum(parameters.admittance_terms, state.admittance_from_state)
    state.history_current_to =
        -2.0 * parameters.characteristic_admittance_s * propagated_to_next +
        _semlyen_state_sum(parameters.admittance_terms, state.admittance_to_state)
    state.write_index = state.write_index == length(state.outgoing_from_history) ?
        1 : state.write_index + 1
    return state
end

line_terminal_voltages(line::SemlyenFrequencyDependentLine) =
    (from = copy(line.terminal_voltage_from), to = copy(line.terminal_voltage_to))
line_terminal_currents(line::SemlyenFrequencyDependentLine) =
    (from = copy(line.terminal_current_from), to = copy(line.terminal_current_to))
line_history_currents(line::SemlyenFrequencyDependentLine) =
    (from = copy(line.history_current_from), to = copy(line.history_current_to))

function stamp!(
    y::AbstractMatrix{Float64},
    rhs::AbstractVector{Float64},
    line::SemlyenFrequencyDependentLine,
    _time_s::Float64,
    timestep_s::Float64,
)
    abs(timestep_s - line.timestep_s) <= 64.0 * eps(Float64) * max(timestep_s, line.timestep_s) ||
        throw(ArgumentError("Semlyen line timestep changed after initialization"))
    for column in eachindex(line.from_nodes), row in eachindex(line.from_nodes)
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
    for phase in eachindex(line.from_nodes)
        stamp_history_current!(rhs, line.from_nodes[phase], 0, line.history_current_from[phase])
        stamp_history_current!(rhs, line.to_nodes[phase], 0, line.history_current_to[phase])
    end
    return nothing
end

function update!(
    line::SemlyenFrequencyDependentLine,
    voltages::AbstractVector{Float64},
    timestep_s::Float64,
)
    abs(timestep_s - line.timestep_s) <= 64.0 * eps(Float64) * max(timestep_s, line.timestep_s) ||
        throw(ArgumentError("Semlyen line timestep changed after initialization"))
    for phase in eachindex(line.from_nodes)
        from_node = line.from_nodes[phase]
        to_node = line.to_nodes[phase]
        line.terminal_voltage_from[phase] = from_node == 0 ? 0.0 : voltages[from_node]
        line.terminal_voltage_to[phase] = to_node == 0 ? 0.0 : voltages[to_node]
    end
    mul!(line.modal_voltage_from, line.phase_to_modal_voltage, line.terminal_voltage_from)
    mul!(line.modal_voltage_to, line.phase_to_modal_voltage, line.terminal_voltage_to)
    for mode in eachindex(line.modes)
        state = line.modes[mode]
        line.modal_current_from[mode] =
            state.parameters.characteristic_admittance_s * line.modal_voltage_from[mode] +
            state.history_current_from
        line.modal_current_to[mode] =
            state.parameters.characteristic_admittance_s * line.modal_voltage_to[mode] +
            state.history_current_to
        _semlyen_mode_update!(state, line.modal_voltage_from[mode], line.modal_voltage_to[mode])
        line.modal_history_current_from[mode] = state.history_current_from
        line.modal_history_current_to[mode] = state.history_current_to
    end
    mul!(
        line.terminal_current_from,
        line.runtime_current_modal_to_phase,
        line.modal_current_from,
    )
    mul!(
        line.terminal_current_to,
        line.runtime_current_modal_to_phase,
        line.modal_current_to,
    )
    mul!(
        line.history_current_from,
        line.runtime_current_modal_to_phase,
        line.modal_history_current_from,
    )
    mul!(
        line.history_current_to,
        line.runtime_current_modal_to_phase,
        line.modal_history_current_to,
    )
    line.update_count += 1
    return nothing
end

trace_output_channel_count(line::SemlyenFrequencyDependentLine) = 4 * length(line.from_nodes)

function trace_output_channel_names!(
    names::Vector{Symbol},
    element_name::Symbol,
    line::SemlyenFrequencyDependentLine,
)
    for phase in eachindex(line.from_nodes)
        append!(names, (
            Symbol(element_name, :_from_current_, phase, :_a),
            Symbol(element_name, :_to_current_, phase, :_a),
            Symbol(element_name, :_from_history_, phase, :_a),
            Symbol(element_name, :_to_history_, phase, :_a),
        ))
    end
    return names
end

function trace_output_values!(
    output::AbstractMatrix{Float64},
    first_channel::Int,
    sample::Int,
    line::SemlyenFrequencyDependentLine,
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

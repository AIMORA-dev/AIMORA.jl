export SemlyenRationalTerm,
       SemlyenModeParameters,
       SemlyenFrequencyDependentLine,
       SemlyenFrequencyScanFitResult,
       PoleResidueTransfer,
       PoleResidueReductionResult,
       RationalLineModeConversion,
       semlyen_rational_terms,
       pole_residue_transfer_value,
       pole_residue_transfer_for_timestep,
       rational_frequency_dependent_mode_parameters,
       semlyen_frequency_dependent_line,
       semlyen_frequency_dependent_line_from_scan_fit,
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
    for root in sort(roots; by = value -> (real(value), imag(value)))
        scale = max(1.0, abs(real(root)))
        abs(imag(root)) <= root_tolerance * scale || throw(ArgumentError(
            "characteristic-admittance poles must be real for the real-time line owner",
        ))
        normalized_pole = -real(root)
        pole = normalized_pole * frequency_scale
        pole > 0.0 || throw(ArgumentError(
            "characteristic-admittance reciprocal has a non-stable pole",
        ))
        residue = _ascending_polynomial_value(denominator, root) /
            _ascending_polynomial_value(derivative, root)
        abs(imag(residue)) <= root_tolerance * max(1.0, abs(real(residue))) ||
            throw(ArgumentError("characteristic-admittance residue is not real"))
        push!(terms, SemlyenRationalTerm(pole, real(residue) / normalized_pole))
    end
    return terms
end

function _rational_terms_value(terms::AbstractVector{SemlyenRationalTerm}, s::ComplexF64)
    value = 0.0 + 0.0im
    for term in terms
        value += term.residue * term.pole / (s + term.pole)
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

function _semlyen_modal_pi_parameters(parameters::SemlyenModeParameters)
    series_impedance = parameters.phasor_series_impedance
    shunt_admittance = parameters.phasor_characteristic_admittance
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
    all(
        mode -> isapprox(
            mode.parameters.phasor_frequency_hz,
            frequency;
            atol = 1.0e-9,
            rtol = 1.0e-9,
        ),
        line.modes,
    ) || throw(ArgumentError(
        "Semlyen supplied phasor parameters do not match the steady-state frequency",
    ))
    series_modal = ComplexF64[]
    shunt_modal = ComplexF64[]
    for state in line.modes
        series_impedance, shunt_admittance =
            _semlyen_modal_pi_parameters(state.parameters)
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

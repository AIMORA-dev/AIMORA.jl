module CoupledLineFitting

using LinearAlgebra
using SHA
using TOML

export CoupledLineFrequencyResponse,
       CoupledLineFitSettings,
       CoupledLineFitAlternative,
       CoupledLineFitRequest,
       CoupledLineModalPreparation,
       CoupledLineRationalModel,
       CoupledLinePassivityCertificate,
       CoupledLineEnforcementReport,
       CoupledLineFitAttempt,
       CoupledLineFitErrorDiagnostics,
       CoupledLineFitUncertainty,
       CoupledLineFitResult,
       coupled_line_terminal_response,
       coupled_line_response_from_scattering,
       coupled_line_modal_preparation,
       coupled_line_admittance_to_scattering,
       coupled_line_scattering_to_admittance,
       coupled_line_rational_model,
       coupled_line_rational_model_from_poles,
       coupled_line_rational_model_from_state_space,
       coupled_line_model_value,
       coupled_line_passivity_certificate,
       coupled_line_fit_signature,
       coupled_line_fit_error_diagnostics,
       coupled_line_fit_uncertainty,
       write_coupled_line_fit,
       read_coupled_line_fit,
       coupled_line_fit_report_text,
       write_coupled_line_fit_report

const COUPLED_LINE_FIT_SCHEMA_VERSION = 2
const COUPLED_LINE_MINIMUM_FREQUENCY_HZ = 0.1
const COUPLED_LINE_MAXIMUM_FREQUENCY_HZ = 1.0e6
const COUPLED_LINE_MAXIMUM_FREQUENCY_COUNT = 401
const COUPLED_LINE_MAXIMUM_PHASE_COUNT = 12
const COUPLED_LINE_MAXIMUM_LENGTH_M = 2.0e6
const COUPLED_LINE_FIT_REQUIRED_TOP_LEVEL_KEYS = Set((
    "schema",
    "schema_version",
    "source_signature_sha256",
    "response_signature_sha256",
    "settings_signature_sha256",
    "settings",
    "deterministic_signature_sha256",
    "maximum_absolute_fit_error",
    "maximum_relative_fit_error",
    "source_response",
    "model",
    "certificate_before_enforcement",
    "certificate_after_enforcement",
    "enforcement",
    "attempts",
    "fit_error_before_enforcement",
    "fit_error_after_enforcement",
    "fitted_scattering_matrices",
    "fitted_terminal_admittance_matrices_s",
))
const COUPLED_LINE_FIT_OPTIONAL_TOP_LEVEL_KEYS = Set((
    "modal_preparation",
    "uncertainty",
))

function _checked_signature(value::AbstractString, label::AbstractString)
    signature = lowercase(String(value))
    occursin(r"^[0-9a-f]{64}$", signature) || throw(ArgumentError(
        "$label must be lowercase 64-hex",
    ))
    return signature
end

function _checked_frequencies(values)
    frequencies = Float64.(values)
    2 <= length(frequencies) <= COUPLED_LINE_MAXIMUM_FREQUENCY_COUNT ||
        throw(ArgumentError("coupled line fitting requires 2 through 401 frequencies"))
    all(value -> isfinite(value) &&
        COUPLED_LINE_MINIMUM_FREQUENCY_HZ <= value <= COUPLED_LINE_MAXIMUM_FREQUENCY_HZ,
        frequencies) || throw(ArgumentError(
        "coupled line fitting frequencies must be finite inside the released band",
    ))
    issorted(frequencies) && all(diff(frequencies) .> 0.0) || throw(ArgumentError(
        "coupled line fitting frequencies must be strictly increasing and unique",
    ))
    return frequencies
end

function _checked_complex_matrix(matrix, dimension::Int, label::AbstractString)
    size(matrix) == (dimension, dimension) || throw(ArgumentError(
        "$label dimensions must be $dimension by $dimension",
    ))
    result = Matrix{ComplexF64}(matrix)
    all(value -> isfinite(real(value)) && isfinite(imag(value)), result) ||
        throw(ArgumentError("$label entries must be finite"))
    return result
end

function _checked_reference_impedance(values, port_count::Int)
    reference = Float64.(values)
    length(reference) == port_count || throw(ArgumentError(
        "reference impedance count must match coupled line ports",
    ))
    all(value -> isfinite(value) && 0.01 <= value <= 100_000.0, reference) ||
        throw(ArgumentError(
            "reference impedances must be finite from 0.01 through 100000 ohms",
        ))
    return reference
end

function _line_response_signature(
    source_signature,
    segment_id,
    segment_kind,
    phases,
    frequencies,
    length_m,
    reference,
    terminal,
)
    io = IOBuffer()
    println(io, COUPLED_LINE_FIT_SCHEMA_VERSION, '|', source_signature)
    println(io, segment_id, '|', segment_kind, "|uniform_segment")
    println(io, join(phases, ','), '|', bitstring(length_m))
    for value in reference
        println(io, bitstring(value))
    end
    for index in eachindex(frequencies)
        println(io, bitstring(frequencies[index]))
        for value in terminal[index]
            println(io, bitstring(real(value)), '|', bitstring(imag(value)))
        end
    end
    return bytes2hex(sha256(take!(io)))
end

"""Immutable frequency samples of one physical coupled line two-terminal response."""
struct CoupledLineFrequencyResponse
    schema_version::Int
    source_signature_sha256::String
    response_signature_sha256::String
    segment_id::Symbol
    segment_kind::Symbol
    source_scope::Symbol
    phase_order::Vector{Symbol}
    port_order::Vector{Symbol}
    frequencies_hz::Vector{Float64}
    length_m::Float64
    reference_impedance_ohm::Vector{Float64}
    terminal_admittance_matrices_s::Vector{Matrix{ComplexF64}}
    scattering_matrices::Vector{Matrix{ComplexF64}}
    construction_residuals::Vector{Float64}
    reciprocity_errors::Vector{Float64}
    minimum_physical_loss_eigenvalues_s::Vector{Float64}
    deterministic_signature_sha256::String
end

function _checked_uniform_segment_identity(segment_id::Symbol, segment_kind::Symbol)
    isempty(String(segment_id)) && throw(ArgumentError(
        "coupled line uniform segment identity must not be empty",
    ))
    segment_kind in (:overhead, :cable, :imported_uniform, :manufactured_uniform) ||
        throw(ArgumentError(
            "coupled line fitting accepts one typed uniform segment, not a mixed-route average",
        ))
    return segment_id, segment_kind
end

struct CoupledLineModalPreparation
    frequencies_hz::Vector{Float64}
    phase_order::Vector{Symbol}
    length_m::Float64
    propagation_constants_per_m::Vector{Vector{ComplexF64}}
    modal_to_phase_matrices::Vector{Matrix{ComplexF64}}
    phase_to_modal_matrices::Vector{Matrix{ComplexF64}}
    mode_assignments::Vector{Vector{Int}}
    minimum_mode_overlaps::Vector{Float64}
    maximum_principal_angles_rad::Vector{Float64}
    extracted_delays_s::Vector{Float64}
    delay_phase_residuals_rad::Vector{Float64}
    constant_transform_maximum_relative_variation::Float64
    deterministic_signature_sha256::String
end

function _causal_propagation_root(value::ComplexF64)
    root = sqrt(value)
    if real(root) < 0.0 || (iszero(real(root)) && imag(root) < 0.0)
        root = -root
    end
    return root
end

function _maximum_overlap_assignment(overlaps::Matrix{Float64})
    count = size(overlaps, 1)
    size(overlaps) == (count, count) || throw(ArgumentError(
        "mode-overlap matrix must be square",
    ))
    state_count = 1 << count
    scores = fill(-Inf, state_count)
    predecessors = fill(-1, state_count)
    selected = fill(0, state_count)
    scores[1] = 0.0
    for mask in 0:(state_count - 1)
        score = scores[mask + 1]
        isfinite(score) || continue
        row = count_ones(mask) + 1
        row > count && continue
        for column in 1:count
            bit = 1 << (column - 1)
            iszero(mask & bit) || continue
            next_mask = mask | bit
            candidate = score + overlaps[row, column]
            if candidate > scores[next_mask + 1] + 16.0 * eps(Float64) ||
                (abs(candidate - scores[next_mask + 1]) <= 16.0 * eps(Float64) &&
                    (selected[next_mask + 1] == 0 || column < selected[next_mask + 1]))
                scores[next_mask + 1] = candidate
                predecessors[next_mask + 1] = mask
                selected[next_mask + 1] = column
            end
        end
    end
    assignment = zeros(Int, count)
    mask = state_count - 1
    for row in count:-1:1
        assignment[row] = selected[mask + 1]
        mask = predecessors[mask + 1]
        mask >= 0 || error("mode-overlap assignment failed")
    end
    return assignment
end

function _normalized_modal_vectors(vectors::AbstractMatrix{<:Number})
    result = Matrix{ComplexF64}(vectors)
    for column in axes(result, 2)
        column_norm = norm(@view result[:, column])
        isfinite(column_norm) && column_norm > 0.0 || throw(ArgumentError(
            "modal eigenvector is nonfinite or singular",
        ))
        result[:, column] ./= column_norm
        pivot = argmax(abs.(@view result[:, column]))
        pivot_value = result[pivot, column]
        iszero(pivot_value) || (result[:, column] .*= exp(-1.0im * angle(pivot_value)))
    end
    return result
end

function _modal_degenerate_components(
    previous_roots::Vector{ComplexF64},
    current_roots::Vector{ComplexF64},
    relative_tolerance::Float64,
)
    count = length(current_roots)
    visited = falses(count)
    components = Vector{Vector{Int}}()
    for first_index in 1:count
        visited[first_index] && continue
        component = Int[first_index]
        visited[first_index] = true
        cursor = 1
        while cursor <= length(component)
            left = component[cursor]
            for right in 1:count
                visited[right] && continue
                scale = max(
                    abs(previous_roots[left]),
                    abs(previous_roots[right]),
                    abs(current_roots[left]),
                    abs(current_roots[right]),
                    1.0,
                )
                close_previous = abs(previous_roots[left] - previous_roots[right]) <=
                    relative_tolerance * scale
                close_current = abs(current_roots[left] - current_roots[right]) <=
                    relative_tolerance * scale
                if close_previous || close_current
                    push!(component, right)
                    visited[right] = true
                end
            end
            cursor += 1
        end
        length(component) > 1 && push!(components, sort!(component))
    end
    return components
end

function _maximum_subspace_angle(previous, current, indices)
    previous_basis = Matrix(qr(previous[:, indices]).Q)[:, 1:length(indices)]
    current_basis = Matrix(qr(current[:, indices]).Q)[:, 1:length(indices)]
    singular_values = svdvals(adjoint(previous_basis) * current_basis)
    minimum_cosine = clamp(minimum(singular_values; init=1.0), 0.0, 1.0)
    return acos(minimum_cosine)
end

function _unwrapped_phase(values::Vector{ComplexF64})
    phases = angle.(values)
    for index in 2:length(phases)
        delta = phases[index] - phases[index - 1]
        while delta > pi
            phases[index] -= 2.0 * pi
            delta -= 2.0 * pi
        end
        while delta < -pi
            phases[index] += 2.0 * pi
            delta += 2.0 * pi
        end
    end
    return phases
end

function _modal_delay(frequencies_hz, propagation_constants, length_m)
    angular_frequencies = 2.0 .* pi .* frequencies_hz
    phases = -imag.(propagation_constants) .* length_m
    frequency_mean = sum(angular_frequencies) / length(angular_frequencies)
    phase_mean = sum(phases) / length(phases)
    denominator = sum(abs2, angular_frequencies .- frequency_mean)
    denominator > 0.0 || throw(ArgumentError("modal delay extraction requires frequency spread"))
    slope = sum(
        (angular_frequencies .- frequency_mean) .*
        (phases .- phase_mean),
    ) / denominator
    delay = -slope
    delay >= -128.0 * eps(Float64) * max(abs(delay), 1.0) || throw(ArgumentError(
        "modal delay extraction produced a negative delay",
    ))
    delay = max(delay, 0.0)
    intercept = phase_mean + delay * frequency_mean
    residual = maximum(abs, phases .- (intercept .- delay .* angular_frequencies))
    return delay, residual
end

function _modal_preparation_signature(frequencies, phases, roots, transforms, delays)
    io = IOBuffer()
    println(io, join(phases, ','))
    for index in eachindex(frequencies)
        println(io, bitstring(frequencies[index]))
        for value in roots[index]
            println(io, bitstring(real(value)), '|', bitstring(imag(value)))
        end
        for value in transforms[index]
            println(io, bitstring(real(value)), '|', bitstring(imag(value)))
        end
    end
    for delay in delays
        println(io, bitstring(delay))
    end
    return bytes2hex(sha256(take!(io)))
end

"""Track coupled modes and extract immutable nonnegative propagation delays."""
function coupled_line_modal_preparation(
    frequencies_hz,
    series_impedance_matrices_ohm_per_m,
    shunt_admittance_matrices_s_per_m,
    length_m::Real;
    phase_order,
    minimum_mode_overlap::Real=0.05,
    degeneracy_relative_tolerance::Real=1.0e-5,
)
    frequencies = _checked_frequencies(frequencies_hz)
    phases = Symbol.(phase_order)
    1 <= length(phases) <= COUPLED_LINE_MAXIMUM_PHASE_COUNT || throw(ArgumentError(
        "coupled modal preparation supports one through twelve phases",
    ))
    length(unique(phases)) == length(phases) || throw(ArgumentError(
        "coupled modal phase identities must be unique",
    ))
    length(series_impedance_matrices_ohm_per_m) == length(frequencies) ==
        length(shunt_admittance_matrices_s_per_m) || throw(ArgumentError(
        "coupled modal Z and Y rows must match the frequency grid",
    ))
    length_value = Float64(length_m)
    isfinite(length_value) && 0.0 < length_value <= COUPLED_LINE_MAXIMUM_LENGTH_M ||
        throw(ArgumentError("coupled modal length must be positive inside the released domain"))
    overlap_limit = Float64(minimum_mode_overlap)
    degeneracy_tolerance = Float64(degeneracy_relative_tolerance)
    isfinite(overlap_limit) && 0.0 <= overlap_limit <= 1.0 &&
        isfinite(degeneracy_tolerance) && degeneracy_tolerance > 0.0 ||
        throw(ArgumentError("coupled modal tracking tolerances are invalid"))
    phase_count = length(phases)
    roots_by_frequency = Vector{ComplexF64}[]
    modal_to_phase = Matrix{ComplexF64}[]
    phase_to_modal = Matrix{ComplexF64}[]
    assignments = Vector{Int}[]
    minimum_overlaps = Float64[]
    maximum_angles = Float64[]
    previous_roots = ComplexF64[]
    previous_vectors = zeros(ComplexF64, phase_count, phase_count)
    for frequency_index in eachindex(frequencies)
        series = _checked_complex_matrix(
            series_impedance_matrices_ohm_per_m[frequency_index],
            phase_count,
            "modal series impedance matrix",
        )
        shunt = _checked_complex_matrix(
            shunt_admittance_matrices_s_per_m[frequency_index],
            phase_count,
            "modal shunt admittance matrix",
        )
        decomposition = eigen(series * shunt)
        roots = _causal_propagation_root.(ComplexF64.(decomposition.values))
        vectors = _normalized_modal_vectors(decomposition.vectors)
        if frequency_index == firstindex(frequencies)
            order = sortperm(eachindex(roots); by=index -> (
                real(roots[index]),
                imag(roots[index]),
                index,
            ))
            roots = roots[order]
            vectors = vectors[:, order]
            assignment = collect(1:phase_count)
            minimum_overlap_value = 1.0
            maximum_angle = 0.0
        else
            overlaps = abs.(adjoint(previous_vectors) * vectors)
            assignment = _maximum_overlap_assignment(overlaps)
            roots = roots[assignment]
            vectors = vectors[:, assignment]
            components = _modal_degenerate_components(
                previous_roots,
                roots,
                degeneracy_tolerance,
            )
            maximum_angle = maximum(
                component -> _maximum_subspace_angle(
                    previous_vectors,
                    vectors,
                    component,
                ),
                components;
                init=0.0,
            )
            for component in components
                previous_basis = Matrix(qr(previous_vectors[:, component]).Q)[:, 1:length(component)]
                current_basis = Matrix(qr(vectors[:, component]).Q)[:, 1:length(component)]
                alignment = svd(adjoint(current_basis) * previous_basis)
                vectors[:, component] .= current_basis * alignment.U * alignment.Vt
            end
            overlap_values = Float64[]
            for mode in 1:phase_count
                alignment = dot(previous_vectors[:, mode], vectors[:, mode])
                push!(overlap_values, abs(alignment))
                iszero(alignment) || (vectors[:, mode] .*= exp(-1.0im * angle(alignment)))
            end
            minimum_overlap_value = minimum(overlap_values; init=1.0)
            minimum_overlap_value >= overlap_limit || throw(ArgumentError(
                "coupled modal tracking overlap is below the declared bound",
            ))
        end
        reciprocal_condition = inv(cond(vectors))
        reciprocal_condition > 64.0 * eps(Float64) || throw(ArgumentError(
            "coupled modal transform is singular or unresolved",
        ))
        push!(roots_by_frequency, roots)
        push!(modal_to_phase, vectors)
        push!(phase_to_modal, inv(vectors))
        push!(assignments, assignment)
        push!(minimum_overlaps, minimum_overlap_value)
        push!(maximum_angles, maximum_angle)
        previous_roots = roots
        previous_vectors = vectors
    end
    delays = Float64[]
    delay_residuals = Float64[]
    for mode in 1:phase_count
        delay, residual = _modal_delay(
            frequencies,
            ComplexF64[roots[mode] for roots in roots_by_frequency],
            length_value,
        )
        push!(delays, delay)
        push!(delay_residuals, residual)
    end
    reference_transform = first(modal_to_phase)
    transform_scale = max(opnorm(reference_transform), 1.0)
    maximum_transform_variation = maximum(
        matrix -> opnorm(matrix - reference_transform) / transform_scale,
        modal_to_phase;
        init=0.0,
    )
    signature = _modal_preparation_signature(
        frequencies,
        phases,
        roots_by_frequency,
        modal_to_phase,
        delays,
    )
    return CoupledLineModalPreparation(
        frequencies,
        phases,
        length_value,
        roots_by_frequency,
        modal_to_phase,
        phase_to_modal,
        assignments,
        minimum_overlaps,
        maximum_angles,
        delays,
        delay_residuals,
        maximum_transform_variation,
        signature,
    )
end

function coupled_line_admittance_to_scattering(
    terminal_admittance,
    reference_impedance_ohm,
)
    dimension = size(terminal_admittance, 1)
    size(terminal_admittance) == (dimension, dimension) || throw(ArgumentError(
        "terminal admittance must be square",
    ))
    reference = _checked_reference_impedance(reference_impedance_ohm, dimension)
    admittance = _checked_complex_matrix(
        terminal_admittance,
        dimension,
        "terminal admittance",
    )
    square_root = Diagonal(sqrt.(reference))
    normalized = square_root * admittance * square_root
    identity_matrix = Matrix{ComplexF64}(I, dimension, dimension)
    denominator = identity_matrix + normalized
    reciprocal_condition = inv(cond(denominator))
    reciprocal_condition > 32.0 * eps(Float64) || throw(ArgumentError(
        "terminal admittance to scattering transformation is singular",
    ))
    return denominator \ (identity_matrix - normalized)
end

function coupled_line_scattering_to_admittance(
    scattering,
    reference_impedance_ohm,
)
    dimension = size(scattering, 1)
    size(scattering) == (dimension, dimension) || throw(ArgumentError(
        "scattering response must be square",
    ))
    reference = _checked_reference_impedance(reference_impedance_ohm, dimension)
    response = _checked_complex_matrix(scattering, dimension, "scattering response")
    identity_matrix = Matrix{ComplexF64}(I, dimension, dimension)
    denominator = identity_matrix + response
    reciprocal_condition = inv(cond(denominator))
    reciprocal_condition > 32.0 * eps(Float64) || throw(ArgumentError(
        "scattering to terminal admittance transformation is singular",
    ))
    normalized = denominator \ (identity_matrix - response)
    inverse_square_root = Diagonal(inv.(sqrt.(reference)))
    return inverse_square_root * normalized * inverse_square_root
end

function _coupled_line_port_order(phases)
    return vcat(
        [Symbol("sending_", phase) for phase in phases],
        [Symbol("receiving_", phase) for phase in phases],
    )
end

function _coupled_line_terminal_matrix(
    series_impedance,
    shunt_admittance,
    length_m::Float64,
)
    phase_count = size(series_impedance, 1)
    zero_block = zeros(ComplexF64, phase_count, phase_count)
    propagation_matrix = [
        zero_block -series_impedance
        -shunt_admittance zero_block
    ]
    transition = exp(propagation_matrix * length_m)
    half_transition = exp(propagation_matrix * (length_m / 2.0))
    first_voltage = @view transition[1:phase_count, 1:phase_count]
    first_current = @view transition[1:phase_count, (phase_count + 1):(2 * phase_count)]
    second_voltage = @view transition[(phase_count + 1):(2 * phase_count), 1:phase_count]
    second_current = @view transition[(phase_count + 1):(2 * phase_count), (phase_count + 1):(2 * phase_count)]
    reciprocal_condition = inv(cond(first_current))
    reciprocal_condition > 32.0 * eps(Float64) || throw(ArgumentError(
        "coupled line terminal response block is singular or unresolved",
    ))
    identity_matrix = Matrix{ComplexF64}(I, phase_count, phase_count)
    sending_current_map = first_current \ hcat(-first_voltage, identity_matrix)
    receiving_current_map = -hcat(second_voltage, second_current) * vcat(
        hcat(identity_matrix, zero_block),
        sending_current_map,
    )
    terminal = vcat(sending_current_map, receiving_current_map)
    transition_residual = opnorm(transition - half_transition * half_transition, Inf)
    block_residual = opnorm(
        first_current * sending_current_map - hcat(-first_voltage, identity_matrix),
        Inf,
    )
    scale = max(opnorm(transition, Inf), opnorm(terminal, Inf), 1.0)
    return Matrix{ComplexF64}(terminal), max(transition_residual, block_residual) / scale
end

"""
    coupled_line_terminal_response(frequencies, Z, Y, length; ...)

Construct the complete inward-current two-terminal response of one uniform
multiconductor segment from SI per-metre matrices. No rational fit or runtime
history is created.
"""
function coupled_line_terminal_response(
    frequencies_hz,
    series_impedance_matrices_ohm_per_m,
    shunt_admittance_matrices_s_per_m,
    length_m::Real;
    phase_order,
    reference_impedance_ohm,
    source_signature_sha256::AbstractString,
    segment_id::Symbol,
    segment_kind::Symbol,
)
    frequencies = _checked_frequencies(frequencies_hz)
    phases = Symbol.(phase_order)
    1 <= length(phases) <= COUPLED_LINE_MAXIMUM_PHASE_COUNT || throw(ArgumentError(
        "coupled line fitting supports one through twelve phases",
    ))
    length(unique(phases)) == length(phases) || throw(ArgumentError(
        "coupled line phase identities must be unique",
    ))
    length(series_impedance_matrices_ohm_per_m) == length(frequencies) ==
        length(shunt_admittance_matrices_s_per_m) || throw(ArgumentError(
        "coupled line Z and Y matrix rows must match the frequency grid",
    ))
    length_value = Float64(length_m)
    isfinite(length_value) && 0.0 < length_value <= COUPLED_LINE_MAXIMUM_LENGTH_M ||
        throw(ArgumentError("coupled line length must be positive inside the released domain"))
    source_signature = _checked_signature(source_signature_sha256, "line source signature")
    uniform_segment_id, uniform_segment_kind = _checked_uniform_segment_identity(
        segment_id,
        segment_kind,
    )
    phase_count = length(phases)
    port_count = 2 * phase_count
    reference = _checked_reference_impedance(reference_impedance_ohm, port_count)
    terminal = Matrix{ComplexF64}[]
    scattering = Matrix{ComplexF64}[]
    residuals = Float64[]
    reciprocity = Float64[]
    minimum_loss = Float64[]
    for index in eachindex(frequencies)
        series = _checked_complex_matrix(
            series_impedance_matrices_ohm_per_m[index],
            phase_count,
            "series impedance matrix",
        )
        shunt = _checked_complex_matrix(
            shunt_admittance_matrices_s_per_m[index],
            phase_count,
            "shunt admittance matrix",
        )
        terminal_matrix, residual = _coupled_line_terminal_matrix(
            series,
            shunt,
            length_value,
        )
        response_scale = max(opnorm(terminal_matrix, Inf), 1.0e-15)
        reciprocity_error = opnorm(terminal_matrix - transpose(terminal_matrix), Inf) /
            response_scale
        reciprocity_error <= 1.0e-8 || throw(ArgumentError(
            "coupled line terminal response violates reciprocity",
        ))
        loss = eigmin(Hermitian((terminal_matrix + adjoint(terminal_matrix)) / 2.0))
        loss >= -1.0e-9 * response_scale || throw(ArgumentError(
            "coupled line terminal response violates passive physical loss",
        ))
        push!(terminal, terminal_matrix)
        push!(scattering, coupled_line_admittance_to_scattering(terminal_matrix, reference))
        push!(residuals, residual)
        push!(reciprocity, reciprocity_error)
        push!(minimum_loss, loss)
    end
    response_signature = _line_response_signature(
        source_signature,
        uniform_segment_id,
        uniform_segment_kind,
        phases,
        frequencies,
        length_value,
        reference,
        terminal,
    )
    return CoupledLineFrequencyResponse(
        COUPLED_LINE_FIT_SCHEMA_VERSION,
        source_signature,
        response_signature,
        uniform_segment_id,
        uniform_segment_kind,
        :uniform_segment,
        phases,
        _coupled_line_port_order(phases),
        frequencies,
        length_value,
        reference,
        terminal,
        scattering,
        residuals,
        reciprocity,
        minimum_loss,
        response_signature,
    )
end

function coupled_line_response_from_scattering(
    frequencies_hz,
    scattering_matrices;
    phase_order,
    reference_impedance_ohm,
    source_signature_sha256::AbstractString,
    segment_id::Symbol,
    segment_kind::Symbol,
    length_m::Real=1.0,
)
    frequencies = _checked_frequencies(frequencies_hz)
    phases = Symbol.(phase_order)
    1 <= length(phases) <= COUPLED_LINE_MAXIMUM_PHASE_COUNT || throw(ArgumentError(
        "coupled line fitting supports one through twelve phases",
    ))
    length(unique(phases)) == length(phases) || throw(ArgumentError(
        "coupled line phase identities must be unique",
    ))
    port_count = 2 * length(phases)
    length(scattering_matrices) == length(frequencies) || throw(ArgumentError(
        "scattering matrix rows must match the frequency grid",
    ))
    reference = _checked_reference_impedance(reference_impedance_ohm, port_count)
    scattering = [
        _checked_complex_matrix(matrix, port_count, "scattering response")
        for matrix in scattering_matrices
    ]
    terminal = [
        coupled_line_scattering_to_admittance(matrix, reference)
        for matrix in scattering
    ]
    length_value = Float64(length_m)
    isfinite(length_value) && 0.0 < length_value <= COUPLED_LINE_MAXIMUM_LENGTH_M ||
        throw(ArgumentError("coupled line length must be positive inside the released domain"))
    source_signature = _checked_signature(source_signature_sha256, "line source signature")
    uniform_segment_id, uniform_segment_kind = _checked_uniform_segment_identity(
        segment_id,
        segment_kind,
    )
    reciprocity = [
        opnorm(matrix - transpose(matrix), Inf) / max(opnorm(matrix, Inf), 1.0e-15)
        for matrix in terminal
    ]
    minimum_loss = [
        eigmin(Hermitian((matrix + adjoint(matrix)) / 2.0))
        for matrix in terminal
    ]
    signature = _line_response_signature(
        source_signature,
        uniform_segment_id,
        uniform_segment_kind,
        phases,
        frequencies,
        length_value,
        reference,
        terminal,
    )
    return CoupledLineFrequencyResponse(
        COUPLED_LINE_FIT_SCHEMA_VERSION,
        source_signature,
        signature,
        uniform_segment_id,
        uniform_segment_kind,
        :uniform_segment,
        phases,
        _coupled_line_port_order(phases),
        frequencies,
        length_value,
        reference,
        terminal,
        scattering,
        zeros(length(frequencies)),
        reciprocity,
        minimum_loss,
        signature,
    )
end

struct CoupledLineFitSettings
    candidate_orders::Vector{Int}
    relative_fit_tolerance::Float64
    maximum_weighted_condition_number::Float64
    minimum_effective_rank_fraction::Float64
    maximum_direct_term_singular_value::Float64
    passivity_tolerance::Float64
    stable_pole_margin_per_s::Float64
    maximum_enforcement_relative_perturbation::Float64
    maximum_relocation_sweeps::Int
    enforce_passivity::Bool
    initial_poles_per_s::Vector{ComplexF64}
    deterministic_signature_sha256::String
end

function _checked_conjugate_poles_and_residues(
    poles_per_s,
    residue_matrices_per_s,
    dimension::Int,
)
    poles = ComplexF64.(poles_per_s)
    residues = Matrix{ComplexF64}.(residue_matrices_per_s)
    !isempty(poles) && length(residues) == length(poles) || throw(ArgumentError(
        "coupled line pole and residue counts must be nonzero and equal",
    ))
    all(pole -> isfinite(real(pole)) && isfinite(imag(pole)) && real(pole) < 0.0, poles) ||
        throw(ArgumentError("coupled line poles must be finite in the open left half-plane"))
    all(matrix -> size(matrix) == (dimension, dimension) &&
        all(value -> isfinite(real(value)) && isfinite(imag(value)), matrix), residues) ||
        throw(ArgumentError("coupled line residue matrices must be finite and square"))
    tolerance = 256.0 * eps(Float64)
    used = falses(length(poles))
    real_indices = Int[]
    pair_indices = Tuple{Int,Int}[]
    for index in eachindex(poles)
        used[index] && continue
        pole = poles[index]
        scale = max(abs(pole), 1.0)
        if abs(imag(pole)) <= tolerance * scale
            residue_scale = max(opnorm(residues[index]), 1.0)
            opnorm(imag.(residues[index])) <= tolerance * residue_scale || throw(ArgumentError(
                "a real coupled line pole must have a real residue matrix",
            ))
            poles[index] = complex(real(pole), 0.0)
            residues[index] = complex.(real.(residues[index]), 0.0)
            push!(real_indices, index)
            used[index] = true
            continue
        end
        imag(pole) > 0.0 || throw(ArgumentError(
            "coupled line conjugate pole pairs must list the positive-imaginary member first",
        ))
        partner = findfirst(eachindex(poles)) do candidate
            !used[candidate] && candidate != index &&
                abs(poles[candidate] - conj(pole)) <= tolerance * scale
        end
        partner === nothing && throw(ArgumentError(
            "every complex coupled line pole must have an explicit conjugate partner",
        ))
        residue_scale = max(opnorm(residues[index]), opnorm(residues[partner]), 1.0)
        opnorm(residues[partner] - conj.(residues[index])) <= tolerance * residue_scale ||
            throw(ArgumentError(
                "complex coupled line pole residues must be conjugate pairs",
            ))
        push!(pair_indices, (index, partner))
        used[index] = true
        used[partner] = true
    end
    canonical_indices = vcat(
        sort(real_indices; by=index -> (real(poles[index]), index)),
        reduce(vcat, (
            collect(pair) for pair in sort(pair_indices; by=pair -> (
                real(poles[first(pair)]),
                imag(poles[first(pair)]),
                first(pair),
            ))
        ); init=Int[]),
    )
    return poles[canonical_indices], residues[canonical_indices]
end

function _checked_initial_poles(values)
    isempty(values) && return ComplexF64[]
    poles = ComplexF64.(values)
    dimension = 1
    dummy_residues = [zeros(ComplexF64, dimension, dimension) for _ in poles]
    checked, _ = _checked_conjugate_poles_and_residues(poles, dummy_residues, dimension)
    return checked
end

function CoupledLineFitSettings(;
    candidate_orders=collect(2:2:12),
    relative_fit_tolerance::Real=0.05,
    maximum_weighted_condition_number::Real=1.0e8,
    minimum_effective_rank_fraction::Real=0.85,
    maximum_direct_term_singular_value::Real=0.92,
    passivity_tolerance::Real=1.0e-8,
    stable_pole_margin_per_s::Real=1.0e-9,
    maximum_enforcement_relative_perturbation::Real=0.02,
    maximum_relocation_sweeps::Integer=64,
    enforce_passivity::Bool=true,
    initial_poles_per_s=ComplexF64[],
)
    orders = sort!(unique(Int.(candidate_orders)))
    !isempty(orders) && all(order -> 1 <= order <= 80, orders) || throw(ArgumentError(
        "coupled line fit candidate orders must be unique from one through eighty",
    ))
    fit_tolerance = Float64(relative_fit_tolerance)
    condition_limit = Float64(maximum_weighted_condition_number)
    rank_fraction_limit = Float64(minimum_effective_rank_fraction)
    direct_term_limit = Float64(maximum_direct_term_singular_value)
    passivity_limit = Float64(passivity_tolerance)
    pole_margin = Float64(stable_pole_margin_per_s)
    perturbation_limit = Float64(maximum_enforcement_relative_perturbation)
    all(value -> isfinite(value) && value > 0.0, (
        fit_tolerance,
        condition_limit,
        passivity_limit,
        pole_margin,
        perturbation_limit,
    )) || throw(ArgumentError("coupled line fit tolerances must be finite and positive"))
    condition_limit >= 1.0 || throw(ArgumentError(
        "coupled line fit maximum condition number must be at least one",
    ))
    0.0 < rank_fraction_limit <= 1.0 || throw(ArgumentError(
        "coupled line fit minimum effective-rank fraction must be in (0, 1]",
    ))
    isfinite(direct_term_limit) && 0.0 < direct_term_limit < 1.0 ||
        throw(ArgumentError(
            "coupled line fit direct-term singular-value limit must be in (0, 1)",
        ))
    sweeps = Int(maximum_relocation_sweeps)
    sweeps >= 0 || throw(ArgumentError("pole relocation sweep count must be nonnegative"))
    initial_poles = _checked_initial_poles(initial_poles_per_s)
    isempty(initial_poles) || length(initial_poles) in orders || throw(ArgumentError(
        "the explicit initial pole count must be one of the candidate orders",
    ))
    io = IOBuffer()
    println(io, join(orders, ','))
    for value in (
        fit_tolerance,
        condition_limit,
        rank_fraction_limit,
        direct_term_limit,
        passivity_limit,
        pole_margin,
        perturbation_limit,
    )
        println(io, bitstring(value))
    end
    println(io, sweeps, '|', enforce_passivity)
    for pole in initial_poles
        println(io, bitstring(real(pole)), '|', bitstring(imag(pole)))
    end
    return CoupledLineFitSettings(
        orders,
        fit_tolerance,
        condition_limit,
        rank_fraction_limit,
        direct_term_limit,
        passivity_limit,
        pole_margin,
        perturbation_limit,
        sweeps,
        enforce_passivity,
        initial_poles,
        bytes2hex(sha256(take!(io))),
    )
end

struct CoupledLineFitAlternative
    response::CoupledLineFrequencyResponse
    modal_preparation::Union{Nothing,CoupledLineModalPreparation}

    function CoupledLineFitAlternative(
        response::CoupledLineFrequencyResponse,
        modal_preparation::Union{Nothing,CoupledLineModalPreparation}=nothing,
    )
        _check_modal_response_alignment(response, modal_preparation)
        return new(response, modal_preparation)
    end
end

function _check_modal_response_alignment(response, modal_preparation)
    modal_preparation === nothing && return nothing
    modal_preparation.frequencies_hz == response.frequencies_hz || throw(ArgumentError(
        "coupled line modal preparation frequency grid must match its response",
    ))
    modal_preparation.phase_order == response.phase_order || throw(ArgumentError(
        "coupled line modal preparation phase order must match its response",
    ))
    modal_preparation.length_m == response.length_m || throw(ArgumentError(
        "coupled line modal preparation length must match its response",
    ))
    return nothing
end

struct CoupledLineFitRequest
    response::CoupledLineFrequencyResponse
    settings::CoupledLineFitSettings
    modal_preparation::Union{Nothing,CoupledLineModalPreparation}
    uncertainty_alternatives::Vector{CoupledLineFitAlternative}
    uncertainty_set_complete::Bool

    function CoupledLineFitRequest(
        response::CoupledLineFrequencyResponse,
        settings::CoupledLineFitSettings,
        modal_preparation::Union{Nothing,CoupledLineModalPreparation},
        uncertainty_alternatives::AbstractVector{CoupledLineFitAlternative},
        uncertainty_set_complete::Bool,
    )
        _check_modal_response_alignment(response, modal_preparation)
        alternatives = collect(uncertainty_alternatives)
        signatures = String[]
        for alternative in alternatives
            candidate = alternative.response
            candidate.frequencies_hz == response.frequencies_hz || throw(ArgumentError(
                "coupled line uncertainty frequency grid must match the nominal response",
            ))
            candidate.phase_order == response.phase_order &&
                candidate.port_order == response.port_order || throw(ArgumentError(
                "coupled line uncertainty phase and port order must match the nominal response",
            ))
            candidate.length_m == response.length_m || throw(ArgumentError(
                "coupled line uncertainty length must match the nominal response",
            ))
            candidate.reference_impedance_ohm == response.reference_impedance_ohm ||
                throw(ArgumentError(
                    "coupled line uncertainty reference impedance must match the nominal response",
                ))
            candidate.source_signature_sha256 != response.source_signature_sha256 ||
                throw(ArgumentError(
                    "coupled line uncertainty alternatives require distinct source signatures",
                ))
            push!(signatures, candidate.source_signature_sha256)
        end
        length(unique(signatures)) == length(signatures) || throw(ArgumentError(
            "coupled line uncertainty source signatures must be unique",
        ))
        uncertainty_set_complete && isempty(alternatives) && throw(ArgumentError(
            "a complete coupled line uncertainty set must contain an alternative",
        ))
        return new(
            response,
            settings,
            modal_preparation,
            alternatives,
            uncertainty_set_complete,
        )
    end
end

CoupledLineFitRequest(
    response::CoupledLineFrequencyResponse,
    settings::CoupledLineFitSettings,
) = CoupledLineFitRequest(response, settings, nothing, CoupledLineFitAlternative[], false)

CoupledLineFitRequest(
    response::CoupledLineFrequencyResponse,
    settings::CoupledLineFitSettings,
    modal_preparation::Union{Nothing,CoupledLineModalPreparation},
) = CoupledLineFitRequest(
    response,
    settings,
    modal_preparation,
    CoupledLineFitAlternative[],
    false,
)

struct CoupledLineRationalModel
    poles_per_s::Vector{ComplexF64}
    decay_rates_per_s::Vector{Float64}
    direct_term::Matrix{Float64}
    residue_matrices_per_s::Vector{Matrix{ComplexF64}}
    state_matrix_per_s::Matrix{Float64}
    input_matrix::Matrix{Float64}
    output_matrix_per_s::Matrix{Float64}
    port_order::Vector{Symbol}
    reference_impedance_ohm::Vector{Float64}
    deterministic_signature_sha256::String
end

function _rational_model_signature(state, input, output, direct, ports, reference)
    io = IOBuffer()
    println(io, join(ports, ','))
    for value in reference
        println(io, bitstring(value))
    end
    for matrix in (state, input, output, direct)
        for value in matrix
            println(io, bitstring(Float64(value)))
        end
    end
    return bytes2hex(sha256(take!(io)))
end

function coupled_line_rational_model(
    decay_rates_per_s,
    direct_term,
    residue_matrices_per_s;
    port_order,
    reference_impedance_ohm,
)
    rates = Float64.(decay_rates_per_s)
    !isempty(rates) && all(value -> isfinite(value) && value > 0.0, rates) ||
        throw(ArgumentError("coupled line fit decay rates must be finite and positive"))
    issorted(rates) && length(unique(rates)) == length(rates) || throw(ArgumentError(
        "coupled line fit decay rates must be strictly increasing and unique",
    ))
    return coupled_line_rational_model_from_poles(
        complex.(-rates, 0.0),
        direct_term,
        complex.(Matrix{Float64}.(residue_matrices_per_s));
        port_order,
        reference_impedance_ohm,
    )
end

"""Build a real state-space realization from stable real or conjugate-pair poles."""
function coupled_line_rational_model_from_poles(
    poles_per_s,
    direct_term,
    residue_matrices_per_s;
    port_order,
    reference_impedance_ohm,
)
    dimension = size(direct_term, 1)
    direct = Matrix{Float64}(direct_term)
    size(direct) == (dimension, dimension) && dimension > 0 || throw(ArgumentError(
        "coupled line fit direct term must be nonempty and square",
    ))
    all(isfinite, direct) || throw(ArgumentError("coupled line fit direct term must be finite"))
    poles, residues = _checked_conjugate_poles_and_residues(
        poles_per_s,
        residue_matrices_per_s,
        dimension,
    )
    ports = Symbol.(port_order)
    length(ports) == dimension && length(unique(ports)) == dimension || throw(ArgumentError(
        "coupled line fit port order must contain one unique identity per port",
    ))
    reference = _checked_reference_impedance(reference_impedance_ohm, dimension)
    state_count = dimension * length(poles)
    state_matrix = zeros(Float64, state_count, state_count)
    input_matrix = zeros(Float64, state_count, dimension)
    output_matrix = zeros(Float64, dimension, state_count)
    identity_block = Matrix{Float64}(I, dimension, dimension)
    state_cursor = 1
    pole_index = 1
    while pole_index <= length(poles)
        pole = poles[pole_index]
        if iszero(imag(pole))
            state_range = state_cursor:(state_cursor + dimension - 1)
            state_matrix[state_range, state_range] .= real(pole) .* identity_block
            input_matrix[state_range, :] .= identity_block
            output_matrix[:, state_range] .= real.(residues[pole_index])
            state_cursor += dimension
            pole_index += 1
        else
            conjugate_index = pole_index + 1
            conjugate_index <= length(poles) && poles[conjugate_index] == conj(pole) ||
                error("canonical coupled line pole pairing failed")
            real_range = state_cursor:(state_cursor + dimension - 1)
            imaginary_range = (state_cursor + dimension):(state_cursor + 2 * dimension - 1)
            damping = real(pole)
            oscillation = imag(pole)
            state_matrix[real_range, real_range] .= damping .* identity_block
            state_matrix[real_range, imaginary_range] .= -oscillation .* identity_block
            state_matrix[imaginary_range, real_range] .= oscillation .* identity_block
            state_matrix[imaginary_range, imaginary_range] .= damping .* identity_block
            input_matrix[real_range, :] .= identity_block
            residue = residues[pole_index]
            output_matrix[:, real_range] .= 2.0 .* real.(residue)
            output_matrix[:, imaginary_range] .= -2.0 .* imag.(residue)
            state_cursor += 2 * dimension
            pole_index += 2
        end
    end
    signature = _rational_model_signature(
        state_matrix,
        input_matrix,
        output_matrix,
        direct,
        ports,
        reference,
    )
    return CoupledLineRationalModel(
        poles,
        -real.(poles),
        direct,
        residues,
        state_matrix,
        input_matrix,
        output_matrix,
        ports,
        reference,
        signature,
    )
end

function _state_space_pole_residues(state, input, output)
    decomposition = eigen(state)
    vectors = decomposition.vectors
    reciprocal_vectors = inv(vectors)
    state_poles = ComplexF64.(decomposition.values)
    state_residues = Matrix{ComplexF64}[]
    for index in eachindex(state_poles)
        right_output = output * @view(vectors[:, index:index])
        left_input = @view(reciprocal_vectors[index:index, :]) * input
        push!(state_residues, Matrix{ComplexF64}(right_output * left_input))
    end
    order = sortperm(eachindex(state_poles); by=index -> (
        real(state_poles[index]),
        -abs(imag(state_poles[index])),
        -imag(state_poles[index]),
        index,
    ))
    poles = ComplexF64[]
    residues = Matrix{ComplexF64}[]
    for index in order
        pole = state_poles[index]
        match = findfirst(existing ->
            abs(existing - pole) <= 1.0e-8 * max(abs(existing), abs(pole), 1.0),
            poles,
        )
        if match === nothing
            push!(poles, pole)
            push!(residues, state_residues[index])
        else
            residues[match] .+= state_residues[index]
        end
    end
    canonical_order = sortperm(eachindex(poles); by=index -> (
        real(poles[index]),
        -imag(poles[index]),
        index,
    ))
    return poles[canonical_order], residues[canonical_order]
end

"""Build an immutable coupled rational model from an explicit real state-space realization."""
function coupled_line_rational_model_from_state_space(
    state_matrix_per_s,
    input_matrix,
    output_matrix_per_s,
    direct_term;
    port_order,
    reference_impedance_ohm,
)
    state = Matrix{Float64}(state_matrix_per_s)
    input = Matrix{Float64}(input_matrix)
    output = Matrix{Float64}(output_matrix_per_s)
    direct = Matrix{Float64}(direct_term)
    state_count = size(state, 1)
    port_count = size(direct, 1)
    state_count > 0 && size(state) == (state_count, state_count) || throw(ArgumentError(
        "coupled line state matrix must be nonempty and square",
    ))
    port_count > 0 && size(direct) == (port_count, port_count) || throw(ArgumentError(
        "coupled line direct term must be nonempty and square",
    ))
    size(input) == (state_count, port_count) && size(output) == (port_count, state_count) ||
        throw(DimensionMismatch("coupled line state-space input and output dimensions disagree"))
    all(matrix -> all(isfinite, matrix), (state, input, output, direct)) || throw(ArgumentError(
        "coupled line state-space matrices must be finite",
    ))
    maximum(real, eigvals(state); init=-Inf) < 0.0 || throw(ArgumentError(
        "coupled line state-space realization must be strictly stable",
    ))
    ports = Symbol.(port_order)
    length(ports) == port_count && length(unique(ports)) == port_count || throw(ArgumentError(
        "coupled line state-space port order must contain unique identities",
    ))
    reference = _checked_reference_impedance(reference_impedance_ohm, port_count)
    poles, residues = _state_space_pole_residues(state, input, output)
    signature = _rational_model_signature(state, input, output, direct, ports, reference)
    return CoupledLineRationalModel(
        poles,
        -real.(poles),
        direct,
        residues,
        state,
        input,
        output,
        ports,
        reference,
        signature,
    )
end

function coupled_line_model_value(model::CoupledLineRationalModel, s::Number)
    frequency = ComplexF64(s)
    value = ComplexF64.(model.direct_term)
    for index in eachindex(model.poles_per_s)
        value .+= model.residue_matrices_per_s[index] ./
            (frequency - model.poles_per_s[index])
    end
    return value
end

struct CoupledLinePassivityCertificate
    method::Symbol
    gamma::Float64
    stable_realization::Bool
    direct_term_bounded::Bool
    zero_frequency_bounded::Bool
    hamiltonian_crossing_free::Bool
    continuous_passivity_passed::Bool
    maximum_real_pole_per_s::Float64
    direct_maximum_singular_value::Float64
    zero_frequency_maximum_singular_value::Float64
    hamiltonian_minimum_real_separation_per_s::Float64
    hamiltonian_crossing_frequencies_hz::Vector{Float64}
    worst_diagnostic_frequency_hz::Float64
    worst_diagnostic_singular_value::Float64
    worst_incident_direction::Vector{ComplexF64}
    minimum_physical_loss_eigenvalue_s::Float64
    diagnostic_frequencies_hz::Vector{Float64}
end

function _passivity_diagnostic_grid(model, frequencies)
    positive = sort!(unique(Float64.(filter(>(0.0), frequencies))))
    anchors = sort!(unique(filter(>(0.0), vcat(
        2.0 .* pi .* positive,
        -real.(eigvals(model.state_matrix_per_s)),
        abs.(imag.(eigvals(model.state_matrix_per_s))),
        abs.(eigvals(model.state_matrix_per_s)),
    ))))
    isempty(anchors) && throw(ArgumentError("passivity certificate requires a positive rate"))
    rates = Float64[0.0]
    lower = first(anchors) / 1024.0
    upper = last(anchors) * 1024.0
    knots = sort!(unique(vcat(lower, anchors, upper)))
    for index in 1:(length(knots) - 1), fraction in 0:8
        push!(rates, exp(
            log(knots[index]) +
            (fraction / 8.0) * (log(knots[index + 1]) - log(knots[index])),
        ))
    end
    push!(rates, upper)
    return sort!(unique(rates ./ (2.0 * pi)))
end

function _bounded_real_hamiltonian(model::CoupledLineRationalModel, gamma::Float64)
    state = model.state_matrix_per_s
    input = model.input_matrix
    output = model.output_matrix_per_s
    direct = model.direct_term
    port_count = size(direct, 1)
    residual = gamma^2 .* Matrix{Float64}(I, port_count, port_count) -
        transpose(direct) * direct
    eigmin(Symmetric(residual)) > 64.0 * eps(Float64) * max(opnorm(residual), 1.0) ||
        return nothing
    state_shift = state + input * (residual \ (transpose(direct) * output))
    positive_block = input * (residual \ transpose(input))
    identity_ports = Matrix{Float64}(I, port_count, port_count)
    negative_block = transpose(output) *
        (identity_ports + direct * (residual \ transpose(direct))) * output
    return [state_shift positive_block; -negative_block -transpose(state_shift)]
end

function coupled_line_passivity_certificate(
    model::CoupledLineRationalModel;
    frequencies_hz,
    passivity_tolerance::Real=1.0e-8,
    stable_pole_margin_per_s::Real=1.0e-9,
)
    tolerance = Float64(passivity_tolerance)
    pole_margin = Float64(stable_pole_margin_per_s)
    isfinite(tolerance) && tolerance > 0.0 && isfinite(pole_margin) && pole_margin > 0.0 ||
        throw(ArgumentError("passivity certificate tolerances must be finite and positive"))
    gamma = 1.0 + tolerance
    maximum_real_pole = maximum(real, eigvals(model.state_matrix_per_s); init=-Inf)
    stable = maximum_real_pole <= -pole_margin
    direct_gain = opnorm(model.direct_term)
    zero_response = coupled_line_model_value(model, 0.0)
    zero_gain = opnorm(zero_response)
    direct_bounded = direct_gain < gamma
    zero_bounded = zero_gain < gamma
    hamiltonian = direct_bounded ? _bounded_real_hamiltonian(model, gamma) : nothing
    balanced_hamiltonian = hamiltonian === nothing ? nothing : copy(hamiltonian)
    balanced_hamiltonian === nothing ||
        LinearAlgebra.LAPACK.gebal!('B', balanced_hamiltonian)
    hamiltonian_eigenvalues = balanced_hamiltonian === nothing ? ComplexF64[] :
        ComplexF64.(eigvals(balanced_hamiltonian))
    minimum_separation = if isempty(hamiltonian_eigenvalues)
        0.0
    else
        minimum(abs ∘ real, hamiltonian_eigenvalues; init=Inf)
    end
    separation_floor = balanced_hamiltonian === nothing ? Inf :
        256.0 * eps(Float64) * max(opnorm(balanced_hamiltonian), 1.0)
    crossing_free = balanced_hamiltonian !== nothing &&
        minimum_separation > separation_floor
    crossing_frequencies = crossing_free ? Float64[] : sort!(unique((
        abs(imag(value)) / (2.0 * pi)
        for value in hamiltonian_eigenvalues
        if abs(real(value)) <= 8.0 * separation_floor
    )))
    diagnostic_grid = _passivity_diagnostic_grid(model, frequencies_hz)
    worst_gain = -Inf
    worst_frequency = 0.0
    worst_direction = ComplexF64[]
    minimum_loss = Inf
    for frequency in diagnostic_grid
        response = coupled_line_model_value(model, 2.0im * pi * frequency)
        decomposition = svd(response)
        gain = isempty(decomposition.S) ? 0.0 : first(decomposition.S)
        if gain > worst_gain
            worst_gain = gain
            worst_frequency = frequency
            worst_direction = Vector{ComplexF64}(decomposition.V[:, 1])
        end
        admittance = try
            coupled_line_scattering_to_admittance(
                response,
                model.reference_impedance_ohm,
            )
        catch
            nothing
        end
        if admittance === nothing
            minimum_loss = -Inf
        elseif minimum_loss != -Inf
            minimum_loss = min(
                minimum_loss,
                eigmin(Hermitian((admittance + adjoint(admittance)) / 2.0)),
            )
        end
    end
    passed = stable && direct_bounded && zero_bounded && crossing_free
    return CoupledLinePassivityCertificate(
        :continuous_bounded_real_hamiltonian,
        gamma,
        stable,
        direct_bounded,
        zero_bounded,
        crossing_free,
        passed,
        maximum_real_pole,
        direct_gain,
        zero_gain,
        minimum_separation,
        crossing_frequencies,
        worst_frequency,
        worst_gain,
        worst_direction,
        minimum_loss,
        diagnostic_grid,
    )
end

struct CoupledLineEnforcementReport
    applied::Bool
    method::Symbol
    enforcement_parameter::Float64
    scattering_maximum_absolute_perturbation::Float64
    scattering_maximum_relative_perturbation::Float64
    admittance_maximum_absolute_perturbation_s::Float64
    admittance_maximum_relative_perturbation::Float64
    within_budget::Bool
end

struct CoupledLineFitAttempt
    order::Int
    maximum_relative_fit_error::Float64
    maximum_weighted_condition_number::Float64
    effective_rank::Int
    weighted_basis_column_count::Int
    effective_rank_fraction::Float64
    regularization::Float64
    pole_policy::Symbol
    pole_relocation_sweeps::Int
    pole_relocation_converged::Bool
    pole_relocation_outcome::Symbol
    pole_relocation_maximum_relative_change::Float64
    pole_relocation_condition_number::Float64
    pole_relocation_effective_rank::Int
    low_frequency_anchor_applied::Bool
    passivity_passed::Bool
    enforcement_relative_perturbation::Float64
    accepted::Bool
    outcome::Symbol
end

struct CoupledLineFitUncertainty
    alternative_source_signatures_sha256::Vector{String}
    alternative_fit_signatures_sha256::Vector{String}
    maximum_source_scattering_relative_deviation::Float64
    maximum_fit_scattering_relative_deviation::Float64
    maximum_fit_error::Float64
    minimum_passivity_margin::Float64
    maximum_enforcement_relative_perturbation::Float64
    maximum_delay_deviation_s::Float64
    pole_association_status::Symbol
    mode_association_status::Symbol
    complete_set::Bool
end

struct CoupledLineFitErrorDiagnostics
    frequencies_hz::Vector{Float64}
    matrix_absolute_errors::Vector{Float64}
    matrix_relative_errors::Vector{Float64}
    singular_value_absolute_errors::Matrix{Float64}
    singular_value_relative_errors::Matrix{Float64}
    port_absolute_errors::Matrix{Float64}
    port_relative_errors::Matrix{Float64}
    maximum_absolute_error::Float64
    maximum_relative_error::Float64
end

function coupled_line_fit_error_diagnostics(
    frequencies_hz,
    source_matrices,
    fitted_matrices,
)
    frequencies = Float64.(frequencies_hz)
    length(frequencies) == length(source_matrices) == length(fitted_matrices) ||
        throw(DimensionMismatch("coupled line fit error rows do not align"))
    isempty(frequencies) && throw(ArgumentError(
        "coupled line fit error diagnostics require frequency rows",
    ))
    port_count = size(first(source_matrices), 1)
    matrix_absolute = zeros(Float64, length(frequencies))
    matrix_relative = zeros(Float64, length(frequencies))
    singular_absolute = zeros(Float64, length(frequencies), port_count)
    singular_relative = zeros(Float64, length(frequencies), port_count)
    port_absolute = zeros(Float64, length(frequencies), port_count)
    port_relative = zeros(Float64, length(frequencies), port_count)
    for frequency_index in eachindex(frequencies)
        source = _checked_complex_matrix(
            source_matrices[frequency_index],
            port_count,
            "source fit-error matrix",
        )
        fitted = _checked_complex_matrix(
            fitted_matrices[frequency_index],
            port_count,
            "fitted fit-error matrix",
        )
        difference = fitted - source
        matrix_absolute[frequency_index] = opnorm(difference)
        matrix_relative[frequency_index] = matrix_absolute[frequency_index] /
            max(opnorm(source), 1.0e-12)
        source_singular_values = svdvals(source)
        fitted_singular_values = svdvals(fitted)
        singular_absolute[frequency_index, :] .=
            abs.(fitted_singular_values - source_singular_values)
        singular_relative[frequency_index, :] .=
            singular_absolute[frequency_index, :] ./
            max.(source_singular_values, 1.0e-12)
        for port_index in 1:port_count
            source_port = @view source[port_index, :]
            difference_port = @view difference[port_index, :]
            port_absolute[frequency_index, port_index] = norm(difference_port)
            port_relative[frequency_index, port_index] =
                port_absolute[frequency_index, port_index] /
                max(norm(source_port), 1.0e-12)
        end
    end
    all(isfinite, vcat(
        matrix_absolute,
        matrix_relative,
        vec(singular_absolute),
        vec(singular_relative),
        vec(port_absolute),
        vec(port_relative),
    )) || throw(ArgumentError("coupled line fit error diagnostics are nonfinite"))
    return CoupledLineFitErrorDiagnostics(
        frequencies,
        matrix_absolute,
        matrix_relative,
        singular_absolute,
        singular_relative,
        port_absolute,
        port_relative,
        maximum(matrix_absolute),
        maximum(matrix_relative),
    )
end

function coupled_line_fit_uncertainty(nominal, alternatives; complete_set::Bool=false)
    results = CoupledLineFitResult[alternative for alternative in alternatives]
    isempty(results) && return CoupledLineFitUncertainty(
        String[], String[], 0.0, 0.0, nominal.maximum_relative_fit_error,
        nominal.certificate_after_enforcement.gamma -
            nominal.certificate_after_enforcement.worst_diagnostic_singular_value,
        nominal.enforcement.admittance_maximum_relative_perturbation,
        0.0,
        :not_applicable_without_alternatives,
        :not_applicable_without_alternatives,
        false,
    )
    source_deviation = 0.0
    fit_deviation = 0.0
    maximum_fit_error = nominal.maximum_relative_fit_error
    minimum_margin = nominal.certificate_after_enforcement.gamma -
        nominal.certificate_after_enforcement.worst_diagnostic_singular_value
    maximum_perturbation = nominal.enforcement.admittance_maximum_relative_perturbation
    maximum_delay_deviation = 0.0
    for alternative in results
        length(alternative.source_response.scattering_matrices) ==
            length(nominal.source_response.scattering_matrices) || throw(ArgumentError(
                "coupled line uncertainty results must share the nominal frequency grid",
            ))
        for index in eachindex(nominal.source_response.scattering_matrices)
            source_scale = max(
                opnorm(nominal.source_response.scattering_matrices[index]),
                1.0e-12,
            )
            fit_scale = max(opnorm(nominal.fitted_scattering_matrices[index]), 1.0e-12)
            source_deviation = max(
                source_deviation,
                opnorm(alternative.source_response.scattering_matrices[index] -
                    nominal.source_response.scattering_matrices[index]) / source_scale,
            )
            fit_deviation = max(
                fit_deviation,
                opnorm(alternative.fitted_scattering_matrices[index] -
                    nominal.fitted_scattering_matrices[index]) / fit_scale,
            )
        end
        maximum_fit_error = max(maximum_fit_error, alternative.maximum_relative_fit_error)
        minimum_margin = min(
            minimum_margin,
            alternative.certificate_after_enforcement.gamma -
                alternative.certificate_after_enforcement.worst_diagnostic_singular_value,
        )
        maximum_perturbation = max(
            maximum_perturbation,
            alternative.enforcement.admittance_maximum_relative_perturbation,
        )
        if nominal.modal_preparation !== nothing && alternative.modal_preparation !== nothing
            maximum_delay_deviation = max(
                maximum_delay_deviation,
                maximum(abs.(alternative.modal_preparation.extracted_delays_s .-
                    nominal.modal_preparation.extracted_delays_s); init=0.0),
            )
        end
    end
    return CoupledLineFitUncertainty(
        getfield.(results, :source_signature_sha256),
        getfield.(results, :deterministic_signature_sha256),
        source_deviation,
        fit_deviation,
        maximum_fit_error,
        minimum_margin,
        maximum_perturbation,
        maximum_delay_deviation,
        :not_claimed_across_alternatives,
        :not_claimed_across_alternatives,
        complete_set,
    )
end

struct CoupledLineFitResult
    schema_version::Int
    source_signature_sha256::String
    response_signature_sha256::String
    settings_signature_sha256::String
    settings::CoupledLineFitSettings
    modal_preparation::Union{Nothing,CoupledLineModalPreparation}
    model::CoupledLineRationalModel
    certificate_before_enforcement::CoupledLinePassivityCertificate
    certificate_after_enforcement::CoupledLinePassivityCertificate
    enforcement::CoupledLineEnforcementReport
    attempts::Vector{CoupledLineFitAttempt}
    fit_error_before_enforcement::CoupledLineFitErrorDiagnostics
    fit_error_after_enforcement::CoupledLineFitErrorDiagnostics
    source_response::CoupledLineFrequencyResponse
    fitted_scattering_matrices::Vector{Matrix{ComplexF64}}
    fitted_terminal_admittance_matrices_s::Vector{Matrix{ComplexF64}}
    maximum_absolute_fit_error::Float64
    maximum_relative_fit_error::Float64
    uncertainty::Union{Nothing,CoupledLineFitUncertainty}
    deterministic_signature_sha256::String
end

function coupled_line_fit_signature(
    response::CoupledLineFrequencyResponse,
    settings::CoupledLineFitSettings,
    modal_preparation::Union{Nothing,CoupledLineModalPreparation},
    model::CoupledLineRationalModel,
    certificate_before::CoupledLinePassivityCertificate,
    certificate_after::CoupledLinePassivityCertificate,
    enforcement::CoupledLineEnforcementReport,
    fit_error_before::CoupledLineFitErrorDiagnostics,
    fit_error_after::CoupledLineFitErrorDiagnostics,
    attempts::AbstractVector{CoupledLineFitAttempt}=CoupledLineFitAttempt[],
    uncertainty::Union{Nothing,CoupledLineFitUncertainty}=nothing,
)
    return _coupled_line_fit_signature(
        response.deterministic_signature_sha256,
        settings.deterministic_signature_sha256,
        modal_preparation,
        model,
        certificate_before,
        certificate_after,
        enforcement,
        fit_error_before,
        fit_error_after,
        attempts,
        uncertainty,
    )
end

function _coupled_line_fit_signature(
    response_signature::AbstractString,
    settings_signature::AbstractString,
    modal_preparation::Union{Nothing,CoupledLineModalPreparation},
    model::CoupledLineRationalModel,
    certificate_before::CoupledLinePassivityCertificate,
    certificate_after::CoupledLinePassivityCertificate,
    enforcement::CoupledLineEnforcementReport,
    fit_error_before::CoupledLineFitErrorDiagnostics,
    fit_error_after::CoupledLineFitErrorDiagnostics,
    attempts::AbstractVector{CoupledLineFitAttempt}=CoupledLineFitAttempt[],
    uncertainty::Union{Nothing,CoupledLineFitUncertainty}=nothing,
)
    io = IOBuffer()
    println(io, response_signature)
    println(io, settings_signature)
    println(io, modal_preparation === nothing ? "modal=none" :
        "modal=$(modal_preparation.deterministic_signature_sha256)")
    println(io, model.deterministic_signature_sha256)
    for (label, certificate) in (
        (:before, certificate_before),
        (:after, certificate_after),
    )
        println(io, label, '|', certificate.method)
        for value in (
            certificate.gamma,
            certificate.maximum_real_pole_per_s,
            certificate.direct_maximum_singular_value,
            certificate.zero_frequency_maximum_singular_value,
            certificate.hamiltonian_minimum_real_separation_per_s,
            certificate.worst_diagnostic_frequency_hz,
            certificate.worst_diagnostic_singular_value,
            certificate.minimum_physical_loss_eigenvalue_s,
        )
            println(io, bitstring(value))
        end
        println(io, certificate.stable_realization, '|', certificate.direct_term_bounded,
            '|', certificate.zero_frequency_bounded, '|',
            certificate.hamiltonian_crossing_free, '|',
            certificate.continuous_passivity_passed)
        for value in certificate.worst_incident_direction
            println(io, bitstring(real(value)), '|', bitstring(imag(value)))
        end
        for value in certificate.hamiltonian_crossing_frequencies_hz
            println(io, "crossing|", bitstring(value))
        end
        for value in certificate.diagnostic_frequencies_hz
            println(io, bitstring(value))
        end
    end
    println(io, enforcement.applied, '|', enforcement.method, '|', enforcement.within_budget)
    for value in (
        enforcement.enforcement_parameter,
        enforcement.scattering_maximum_absolute_perturbation,
        enforcement.scattering_maximum_relative_perturbation,
        enforcement.admittance_maximum_absolute_perturbation_s,
        enforcement.admittance_maximum_relative_perturbation,
    )
        println(io, bitstring(value))
    end
    for (label, diagnostics) in (
        (:fit_error_before, fit_error_before),
        (:fit_error_after, fit_error_after),
    )
        println(io, label)
        for field in (
            :frequencies_hz,
            :matrix_absolute_errors,
            :matrix_relative_errors,
            :singular_value_absolute_errors,
            :singular_value_relative_errors,
            :port_absolute_errors,
            :port_relative_errors,
            :maximum_absolute_error,
            :maximum_relative_error,
        )
            values = getfield(diagnostics, field)
            if values isa Number
                println(io, bitstring(Float64(values)))
            else
                for value in values
                    println(io, bitstring(value))
                end
            end
        end
    end
    for attempt in attempts
        println(io, attempt.order, '|', bitstring(attempt.maximum_relative_fit_error), '|',
            bitstring(attempt.maximum_weighted_condition_number), '|', attempt.effective_rank,
            '|', attempt.weighted_basis_column_count, '|',
            bitstring(attempt.effective_rank_fraction), '|',
            bitstring(attempt.regularization), '|', attempt.pole_policy, '|',
            attempt.pole_relocation_sweeps, '|', attempt.pole_relocation_converged, '|',
            attempt.pole_relocation_outcome, '|',
            bitstring(attempt.pole_relocation_maximum_relative_change), '|',
            bitstring(attempt.pole_relocation_condition_number), '|',
            attempt.pole_relocation_effective_rank, '|',
            attempt.low_frequency_anchor_applied, '|', attempt.passivity_passed, '|',
            bitstring(attempt.enforcement_relative_perturbation), '|', attempt.accepted,
            '|', attempt.outcome)
    end
    if uncertainty === nothing
        println(io, "uncertainty=unknown")
    else
        println(io, "uncertainty_complete=", uncertainty.complete_set)
        for index in eachindex(uncertainty.alternative_source_signatures_sha256)
            println(io, uncertainty.alternative_source_signatures_sha256[index], '|',
                uncertainty.alternative_fit_signatures_sha256[index])
        end
        for value in (
            uncertainty.maximum_source_scattering_relative_deviation,
            uncertainty.maximum_fit_scattering_relative_deviation,
            uncertainty.maximum_fit_error,
            uncertainty.minimum_passivity_margin,
            uncertainty.maximum_enforcement_relative_perturbation,
            uncertainty.maximum_delay_deviation_s,
        )
            println(io, bitstring(value))
        end
        println(io, uncertainty.pole_association_status)
        println(io, uncertainty.mode_association_status)
    end
    return bytes2hex(sha256(take!(io)))
end

function _complex_matrix_data(matrix)
    return Dict{String,Any}(
        "rows" => size(matrix, 1),
        "columns" => size(matrix, 2),
        "real" => vec(real.(matrix)),
        "imaginary" => vec(imag.(matrix)),
    )
end

function _real_matrix_data(matrix)
    return Dict{String,Any}(
        "rows" => size(matrix, 1),
        "columns" => size(matrix, 2),
        "values" => vec(Float64.(matrix)),
    )
end

function _complex_matrix_from_data(data, label)
    rows = Int(data["rows"])
    columns = Int(data["columns"])
    real_values = Float64.(data["real"])
    imaginary_values = Float64.(data["imaginary"])
    length(real_values) == length(imaginary_values) == rows * columns ||
        throw(ArgumentError("$label matrix storage is malformed"))
    return reshape(complex.(real_values, imaginary_values), rows, columns)
end

function _real_matrix_from_data(data, label)
    rows = Int(data["rows"])
    columns = Int(data["columns"])
    values = Float64.(data["values"])
    length(values) == rows * columns || throw(ArgumentError(
        "$label matrix storage is malformed",
    ))
    return reshape(values, rows, columns)
end

function _fit_error_data(diagnostics::CoupledLineFitErrorDiagnostics)
    return Dict{String,Any}(
        "frequencies_hz" => diagnostics.frequencies_hz,
        "matrix_absolute_errors" => diagnostics.matrix_absolute_errors,
        "matrix_relative_errors" => diagnostics.matrix_relative_errors,
        "singular_value_absolute_errors" =>
            _real_matrix_data(diagnostics.singular_value_absolute_errors),
        "singular_value_relative_errors" =>
            _real_matrix_data(diagnostics.singular_value_relative_errors),
        "port_absolute_errors" => _real_matrix_data(diagnostics.port_absolute_errors),
        "port_relative_errors" => _real_matrix_data(diagnostics.port_relative_errors),
        "maximum_absolute_error" => diagnostics.maximum_absolute_error,
        "maximum_relative_error" => diagnostics.maximum_relative_error,
    )
end

function _fit_error_from_data(data)
    diagnostics = CoupledLineFitErrorDiagnostics(
        Float64.(data["frequencies_hz"]),
        Float64.(data["matrix_absolute_errors"]),
        Float64.(data["matrix_relative_errors"]),
        _real_matrix_from_data(
            data["singular_value_absolute_errors"],
            "fit singular-value absolute error",
        ),
        _real_matrix_from_data(
            data["singular_value_relative_errors"],
            "fit singular-value relative error",
        ),
        _real_matrix_from_data(data["port_absolute_errors"], "fit port absolute error"),
        _real_matrix_from_data(data["port_relative_errors"], "fit port relative error"),
        Float64(data["maximum_absolute_error"]),
        Float64(data["maximum_relative_error"]),
    )
    row_count = length(diagnostics.frequencies_hz)
    length(diagnostics.matrix_absolute_errors) == row_count ==
        length(diagnostics.matrix_relative_errors) || throw(ArgumentError(
        "fit error frequency and matrix rows do not align",
    ))
    all(size(matrix, 1) == row_count for matrix in (
        diagnostics.singular_value_absolute_errors,
        diagnostics.singular_value_relative_errors,
        diagnostics.port_absolute_errors,
        diagnostics.port_relative_errors,
    )) || throw(ArgumentError("fit error frequency and port rows do not align"))
    all(isfinite, vcat(
        diagnostics.frequencies_hz,
        diagnostics.matrix_absolute_errors,
        diagnostics.matrix_relative_errors,
        vec(diagnostics.singular_value_absolute_errors),
        vec(diagnostics.singular_value_relative_errors),
        vec(diagnostics.port_absolute_errors),
        vec(diagnostics.port_relative_errors),
        diagnostics.maximum_absolute_error,
        diagnostics.maximum_relative_error,
    )) || throw(ArgumentError("fit error diagnostics are nonfinite"))
    return diagnostics
end

function _fit_error_diagnostics_match(left, right)
    left.frequencies_hz == right.frequencies_hz || return false
    for field in (
        :matrix_absolute_errors,
        :matrix_relative_errors,
        :singular_value_absolute_errors,
        :singular_value_relative_errors,
        :port_absolute_errors,
        :port_relative_errors,
        :maximum_absolute_error,
        :maximum_relative_error,
    )
        isapprox(getfield(left, field), getfield(right, field);
            rtol=2.0e-12, atol=2.0e-13) || return false
    end
    return true
end

function _frequency_response_data(response::CoupledLineFrequencyResponse)
    return Dict{String,Any}(
        "schema_version" => response.schema_version,
        "source_signature_sha256" => response.source_signature_sha256,
        "response_signature_sha256" => response.response_signature_sha256,
        "segment_id" => String(response.segment_id),
        "segment_kind" => String(response.segment_kind),
        "source_scope" => String(response.source_scope),
        "phase_order" => String.(response.phase_order),
        "port_order" => String.(response.port_order),
        "frequencies_hz" => response.frequencies_hz,
        "length_m" => response.length_m,
        "reference_impedance_ohm" => response.reference_impedance_ohm,
        "terminal_admittance_matrices_s" =>
            _complex_matrix_data.(response.terminal_admittance_matrices_s),
        "scattering_matrices" => _complex_matrix_data.(response.scattering_matrices),
        "construction_residuals" => response.construction_residuals,
        "reciprocity_errors" => response.reciprocity_errors,
        "minimum_physical_loss_eigenvalues_s" =>
            response.minimum_physical_loss_eigenvalues_s,
        "deterministic_signature_sha256" => response.deterministic_signature_sha256,
    )
end

function _fit_settings_data(settings::CoupledLineFitSettings)
    return Dict{String,Any}(
        "candidate_orders" => settings.candidate_orders,
        "relative_fit_tolerance" => settings.relative_fit_tolerance,
        "maximum_weighted_condition_number" =>
            settings.maximum_weighted_condition_number,
        "minimum_effective_rank_fraction" => settings.minimum_effective_rank_fraction,
        "maximum_direct_term_singular_value" =>
            settings.maximum_direct_term_singular_value,
        "passivity_tolerance" => settings.passivity_tolerance,
        "stable_pole_margin_per_s" => settings.stable_pole_margin_per_s,
        "maximum_enforcement_relative_perturbation" =>
            settings.maximum_enforcement_relative_perturbation,
        "maximum_relocation_sweeps" => settings.maximum_relocation_sweeps,
        "enforce_passivity" => settings.enforce_passivity,
        "initial_poles_per_s" => Dict(
            "real" => real.(settings.initial_poles_per_s),
            "imaginary" => imag.(settings.initial_poles_per_s),
        ),
        "deterministic_signature_sha256" => settings.deterministic_signature_sha256,
    )
end

function _fit_settings_from_data(data)
    poles_data = data["initial_poles_per_s"]
    real_values = Float64.(poles_data["real"])
    imaginary_values = Float64.(poles_data["imaginary"])
    length(real_values) == length(imaginary_values) || throw(ArgumentError(
        "fit settings initial pole storage is malformed",
    ))
    settings = CoupledLineFitSettings(
        candidate_orders=Int.(data["candidate_orders"]),
        relative_fit_tolerance=Float64(data["relative_fit_tolerance"]),
        maximum_weighted_condition_number=
            Float64(data["maximum_weighted_condition_number"]),
        minimum_effective_rank_fraction=
            Float64(data["minimum_effective_rank_fraction"]),
        maximum_direct_term_singular_value=
            Float64(data["maximum_direct_term_singular_value"]),
        passivity_tolerance=Float64(data["passivity_tolerance"]),
        stable_pole_margin_per_s=Float64(data["stable_pole_margin_per_s"]),
        maximum_enforcement_relative_perturbation=
            Float64(data["maximum_enforcement_relative_perturbation"]),
        maximum_relocation_sweeps=Int(data["maximum_relocation_sweeps"]),
        enforce_passivity=Bool(data["enforce_passivity"]),
        initial_poles_per_s=complex.(real_values, imaginary_values),
    )
    settings.deterministic_signature_sha256 == data["deterministic_signature_sha256"] ||
        throw(ArgumentError("fit settings signature mismatch"))
    return settings
end

function _frequency_response_from_data(data)
    response = CoupledLineFrequencyResponse(
        Int(data["schema_version"]),
        _checked_signature(data["source_signature_sha256"], "source response signature"),
        _checked_signature(data["response_signature_sha256"], "response identity"),
        Symbol(data["segment_id"]),
        Symbol(data["segment_kind"]),
        Symbol(data["source_scope"]),
        Symbol.(data["phase_order"]),
        Symbol.(data["port_order"]),
        Float64.(data["frequencies_hz"]),
        Float64(data["length_m"]),
        Float64.(data["reference_impedance_ohm"]),
        [_complex_matrix_from_data(row, "source terminal admittance")
            for row in data["terminal_admittance_matrices_s"]],
        [_complex_matrix_from_data(row, "source scattering")
            for row in data["scattering_matrices"]],
        Float64.(data["construction_residuals"]),
        Float64.(data["reciprocity_errors"]),
        Float64.(data["minimum_physical_loss_eigenvalues_s"]),
        _checked_signature(data["deterministic_signature_sha256"],
            "source deterministic signature"),
    )
    response.response_signature_sha256 == response.deterministic_signature_sha256 ||
        throw(ArgumentError("source response signatures disagree"))
    _checked_uniform_segment_identity(response.segment_id, response.segment_kind)
    response.source_scope == :uniform_segment || throw(ArgumentError(
        "source response is not one typed uniform segment",
    ))
    _checked_frequencies(response.frequencies_hz) == response.frequencies_hz ||
        throw(ArgumentError("source response frequency grid is malformed"))
    1 <= length(response.phase_order) <= COUPLED_LINE_MAXIMUM_PHASE_COUNT &&
        length(unique(response.phase_order)) == length(response.phase_order) ||
        throw(ArgumentError("source response phase order is malformed"))
    response.port_order == _coupled_line_port_order(response.phase_order) ||
        throw(ArgumentError("source response port order is malformed"))
    port_count = length(response.port_order)
    _checked_reference_impedance(response.reference_impedance_ohm, port_count) ==
        response.reference_impedance_ohm || throw(ArgumentError(
            "source response reference impedance is malformed",
        ))
    row_count = length(response.frequencies_hz)
    all(length(values) == row_count for values in (
        response.terminal_admittance_matrices_s,
        response.scattering_matrices,
        response.construction_residuals,
        response.reciprocity_errors,
        response.minimum_physical_loss_eigenvalues_s,
    )) || throw(ArgumentError("source response rows do not align"))
    all(isfinite, vcat(
        response.construction_residuals,
        response.reciprocity_errors,
        response.minimum_physical_loss_eigenvalues_s,
    )) || throw(ArgumentError("source response diagnostics are nonfinite"))
    for matrix in response.terminal_admittance_matrices_s
        _checked_complex_matrix(matrix, port_count, "source terminal admittance")
    end
    for matrix in response.scattering_matrices
        _checked_complex_matrix(matrix, port_count, "source scattering")
    end
    expected_signature = _line_response_signature(
        response.source_signature_sha256,
        response.segment_id,
        response.segment_kind,
        response.phase_order,
        response.frequencies_hz,
        response.length_m,
        response.reference_impedance_ohm,
        response.terminal_admittance_matrices_s,
    )
    expected_signature == response.response_signature_sha256 || throw(ArgumentError(
        "source response content signature mismatch",
    ))
    all(eachindex(response.scattering_matrices)) do index
        isapprox(
            coupled_line_admittance_to_scattering(
            response.terminal_admittance_matrices_s[index],
            response.reference_impedance_ohm,
            ),
            response.scattering_matrices[index];
            rtol=1.0e-12,
            atol=1.0e-13,
        )
    end || throw(ArgumentError("source terminal and scattering responses disagree"))
    return response
end

function _certificate_data(certificate::CoupledLinePassivityCertificate)
    return Dict{String,Any}(
        "method" => String(certificate.method),
        "gamma" => certificate.gamma,
        "stable_realization" => certificate.stable_realization,
        "direct_term_bounded" => certificate.direct_term_bounded,
        "zero_frequency_bounded" => certificate.zero_frequency_bounded,
        "hamiltonian_crossing_free" => certificate.hamiltonian_crossing_free,
        "continuous_passivity_passed" => certificate.continuous_passivity_passed,
        "maximum_real_pole_per_s" => certificate.maximum_real_pole_per_s,
        "direct_maximum_singular_value" => certificate.direct_maximum_singular_value,
        "zero_frequency_maximum_singular_value" => certificate.zero_frequency_maximum_singular_value,
        "hamiltonian_minimum_real_separation_per_s" => certificate.hamiltonian_minimum_real_separation_per_s,
        "hamiltonian_crossing_frequencies_hz" =>
            certificate.hamiltonian_crossing_frequencies_hz,
        "worst_diagnostic_frequency_hz" => certificate.worst_diagnostic_frequency_hz,
        "worst_diagnostic_singular_value" => certificate.worst_diagnostic_singular_value,
        "worst_incident_direction_real" => real.(certificate.worst_incident_direction),
        "worst_incident_direction_imaginary" => imag.(certificate.worst_incident_direction),
        "minimum_physical_loss_eigenvalue_s" => certificate.minimum_physical_loss_eigenvalue_s,
        "diagnostic_frequencies_hz" => certificate.diagnostic_frequencies_hz,
    )
end

function _certificate_from_data(data)
    direction_real = Float64.(data["worst_incident_direction_real"])
    direction_imaginary = Float64.(data["worst_incident_direction_imaginary"])
    length(direction_real) == length(direction_imaginary) || throw(ArgumentError(
        "coupled line certificate incident direction is malformed",
    ))
    return CoupledLinePassivityCertificate(
        Symbol(data["method"]),
        Float64(data["gamma"]),
        Bool(data["stable_realization"]),
        Bool(data["direct_term_bounded"]),
        Bool(data["zero_frequency_bounded"]),
        Bool(data["hamiltonian_crossing_free"]),
        Bool(data["continuous_passivity_passed"]),
        Float64(data["maximum_real_pole_per_s"]),
        Float64(data["direct_maximum_singular_value"]),
        Float64(data["zero_frequency_maximum_singular_value"]),
        Float64(data["hamiltonian_minimum_real_separation_per_s"]),
        Float64.(data["hamiltonian_crossing_frequencies_hz"]),
        Float64(data["worst_diagnostic_frequency_hz"]),
        Float64(data["worst_diagnostic_singular_value"]),
        complex.(direction_real, direction_imaginary),
        Float64(data["minimum_physical_loss_eigenvalue_s"]),
        Float64.(data["diagnostic_frequencies_hz"]),
    )
end

function _modal_data(modal::CoupledLineModalPreparation)
    return Dict{String,Any}(
        "frequencies_hz" => modal.frequencies_hz,
        "phase_order" => String.(modal.phase_order),
        "length_m" => modal.length_m,
        "propagation_constants" => [
            Dict(
                "real" => real.(values),
                "imaginary" => imag.(values),
            ) for values in modal.propagation_constants_per_m
        ],
        "modal_to_phase" => _complex_matrix_data.(modal.modal_to_phase_matrices),
        "phase_to_modal" => _complex_matrix_data.(modal.phase_to_modal_matrices),
        "mode_assignments" => modal.mode_assignments,
        "minimum_mode_overlaps" => modal.minimum_mode_overlaps,
        "maximum_principal_angles_rad" => modal.maximum_principal_angles_rad,
        "extracted_delays_s" => modal.extracted_delays_s,
        "delay_phase_residuals_rad" => modal.delay_phase_residuals_rad,
        "constant_transform_maximum_relative_variation" =>
            modal.constant_transform_maximum_relative_variation,
        "deterministic_signature_sha256" => modal.deterministic_signature_sha256,
    )
end

function _modal_from_data(data)
    roots = Vector{ComplexF64}[
        complex.(Float64.(row["real"]), Float64.(row["imaginary"]))
        for row in data["propagation_constants"]
    ]
    modal = CoupledLineModalPreparation(
        Float64.(data["frequencies_hz"]),
        Symbol.(data["phase_order"]),
        Float64(data["length_m"]),
        roots,
        [_complex_matrix_from_data(row, "modal-to-phase") for row in data["modal_to_phase"]],
        [_complex_matrix_from_data(row, "phase-to-modal") for row in data["phase_to_modal"]],
        [Int.(row) for row in data["mode_assignments"]],
        Float64.(data["minimum_mode_overlaps"]),
        Float64.(data["maximum_principal_angles_rad"]),
        Float64.(data["extracted_delays_s"]),
        Float64.(data["delay_phase_residuals_rad"]),
        Float64(data["constant_transform_maximum_relative_variation"]),
        _checked_signature(
            data["deterministic_signature_sha256"],
            "modal preparation signature",
        ),
    )
    signature = _modal_preparation_signature(
        modal.frequencies_hz,
        modal.phase_order,
        modal.propagation_constants_per_m,
        modal.modal_to_phase_matrices,
        modal.extracted_delays_s,
    )
    signature == modal.deterministic_signature_sha256 || throw(ArgumentError(
        "coupled line modal preparation signature mismatch",
    ))
    return modal
end

"""Write a complete deterministic coupled-line fit interchange artifact."""
function write_coupled_line_fit(path::AbstractString, result::CoupledLineFitResult)
    model = result.model
    data = Dict{String,Any}(
        "schema" => "aimora.coupled_line_fit.v1",
        "schema_version" => result.schema_version,
        "source_signature_sha256" => result.source_signature_sha256,
        "response_signature_sha256" => result.response_signature_sha256,
        "settings_signature_sha256" => result.settings_signature_sha256,
        "settings" => _fit_settings_data(result.settings),
        "deterministic_signature_sha256" => result.deterministic_signature_sha256,
        "maximum_absolute_fit_error" => result.maximum_absolute_fit_error,
        "maximum_relative_fit_error" => result.maximum_relative_fit_error,
        "source_response" => _frequency_response_data(result.source_response),
        "model" => Dict(
            "poles_per_s" => Dict(
                "real" => real.(model.poles_per_s),
                "imaginary" => imag.(model.poles_per_s),
            ),
            "decay_rates_per_s" => model.decay_rates_per_s,
            "direct_term" => _real_matrix_data(model.direct_term),
            "residue_matrices_per_s" => _complex_matrix_data.(model.residue_matrices_per_s),
            "state_matrix_per_s" => _real_matrix_data(model.state_matrix_per_s),
            "input_matrix" => _real_matrix_data(model.input_matrix),
            "output_matrix_per_s" => _real_matrix_data(model.output_matrix_per_s),
            "port_order" => String.(model.port_order),
            "reference_impedance_ohm" => model.reference_impedance_ohm,
            "deterministic_signature_sha256" => model.deterministic_signature_sha256,
        ),
        "certificate_before_enforcement" =>
            _certificate_data(result.certificate_before_enforcement),
        "certificate_after_enforcement" =>
            _certificate_data(result.certificate_after_enforcement),
        "enforcement" => Dict(
            "applied" => result.enforcement.applied,
            "method" => String(result.enforcement.method),
            "enforcement_parameter" => result.enforcement.enforcement_parameter,
            "scattering_maximum_absolute_perturbation" =>
                result.enforcement.scattering_maximum_absolute_perturbation,
            "scattering_maximum_relative_perturbation" =>
                result.enforcement.scattering_maximum_relative_perturbation,
            "admittance_maximum_absolute_perturbation_s" =>
                result.enforcement.admittance_maximum_absolute_perturbation_s,
            "admittance_maximum_relative_perturbation" =>
                result.enforcement.admittance_maximum_relative_perturbation,
            "within_budget" => result.enforcement.within_budget,
        ),
        "attempts" => [
            Dict(
                "order" => attempt.order,
                "maximum_relative_fit_error" => attempt.maximum_relative_fit_error,
                "maximum_weighted_condition_number" =>
                    attempt.maximum_weighted_condition_number,
                "effective_rank" => attempt.effective_rank,
                "weighted_basis_column_count" => attempt.weighted_basis_column_count,
                "effective_rank_fraction" => attempt.effective_rank_fraction,
                "regularization" => attempt.regularization,
                "pole_policy" => String(attempt.pole_policy),
                "pole_relocation_sweeps" => attempt.pole_relocation_sweeps,
                "pole_relocation_converged" => attempt.pole_relocation_converged,
                "pole_relocation_outcome" => String(attempt.pole_relocation_outcome),
                "pole_relocation_maximum_relative_change" =>
                    attempt.pole_relocation_maximum_relative_change,
                "pole_relocation_condition_number" =>
                    attempt.pole_relocation_condition_number,
                "pole_relocation_effective_rank" =>
                    attempt.pole_relocation_effective_rank,
                "low_frequency_anchor_applied" =>
                    attempt.low_frequency_anchor_applied,
                "passivity_passed" => attempt.passivity_passed,
                "enforcement_relative_perturbation" =>
                    attempt.enforcement_relative_perturbation,
                "accepted" => attempt.accepted,
                "outcome" => String(attempt.outcome),
            ) for attempt in result.attempts
        ],
        "fit_error_before_enforcement" =>
            _fit_error_data(result.fit_error_before_enforcement),
        "fit_error_after_enforcement" =>
            _fit_error_data(result.fit_error_after_enforcement),
        "fitted_scattering_matrices" =>
            _complex_matrix_data.(result.fitted_scattering_matrices),
        "fitted_terminal_admittance_matrices_s" =>
            _complex_matrix_data.(result.fitted_terminal_admittance_matrices_s),
    )
    result.modal_preparation === nothing ||
        (data["modal_preparation"] = _modal_data(result.modal_preparation))
    result.uncertainty === nothing ||
        (data["uncertainty"] = _uncertainty_data(result.uncertainty))
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return path
end

function _enforcement_from_data(data)
    return CoupledLineEnforcementReport(
        Bool(data["applied"]),
        Symbol(data["method"]),
        Float64(data["enforcement_parameter"]),
        Float64(data["scattering_maximum_absolute_perturbation"]),
        Float64(data["scattering_maximum_relative_perturbation"]),
        Float64(data["admittance_maximum_absolute_perturbation_s"]),
        Float64(data["admittance_maximum_relative_perturbation"]),
        Bool(data["within_budget"]),
    )
end

function _uncertainty_data(uncertainty::CoupledLineFitUncertainty)
    return Dict{String,Any}(
        "alternative_source_signatures_sha256" =>
            uncertainty.alternative_source_signatures_sha256,
        "alternative_fit_signatures_sha256" => uncertainty.alternative_fit_signatures_sha256,
        "maximum_source_scattering_relative_deviation" =>
            uncertainty.maximum_source_scattering_relative_deviation,
        "maximum_fit_scattering_relative_deviation" =>
            uncertainty.maximum_fit_scattering_relative_deviation,
        "maximum_fit_error" => uncertainty.maximum_fit_error,
        "minimum_passivity_margin" => uncertainty.minimum_passivity_margin,
        "maximum_enforcement_relative_perturbation" =>
            uncertainty.maximum_enforcement_relative_perturbation,
        "maximum_delay_deviation_s" => uncertainty.maximum_delay_deviation_s,
        "pole_association_status" => String(uncertainty.pole_association_status),
        "mode_association_status" => String(uncertainty.mode_association_status),
        "complete_set" => uncertainty.complete_set,
    )
end

function _uncertainty_from_data(data)
    source_signatures = _checked_signature.(
        String.(data["alternative_source_signatures_sha256"]),
        Ref("uncertainty source signature"),
    )
    fit_signatures = _checked_signature.(
        String.(data["alternative_fit_signatures_sha256"]),
        Ref("uncertainty fit signature"),
    )
    length(source_signatures) == length(fit_signatures) || throw(ArgumentError(
        "coupled line uncertainty signature counts disagree",
    ))
    return CoupledLineFitUncertainty(
        source_signatures,
        fit_signatures,
        Float64(data["maximum_source_scattering_relative_deviation"]),
        Float64(data["maximum_fit_scattering_relative_deviation"]),
        Float64(data["maximum_fit_error"]),
        Float64(data["minimum_passivity_margin"]),
        Float64(data["maximum_enforcement_relative_perturbation"]),
        Float64(data["maximum_delay_deviation_s"]),
        Symbol(data["pole_association_status"]),
        Symbol(data["mode_association_status"]),
        Bool(data["complete_set"]),
    )
end

"""Read, hash-check, and continuously recertify a coupled-line fit artifact."""
function read_coupled_line_fit(
    path::AbstractString;
    expected_source_signature_sha256::Union{Nothing,AbstractString}=nothing,
    expected_response_signature_sha256::Union{Nothing,AbstractString}=nothing,
)
    data = TOML.parsefile(path)
    supplied_keys = Set(String.(keys(data)))
    missing_keys = setdiff(COUPLED_LINE_FIT_REQUIRED_TOP_LEVEL_KEYS, supplied_keys)
    unknown_keys = setdiff(
        supplied_keys,
        union(
            COUPLED_LINE_FIT_REQUIRED_TOP_LEVEL_KEYS,
            COUPLED_LINE_FIT_OPTIONAL_TOP_LEVEL_KEYS,
        ),
    )
    isempty(missing_keys) || throw(ArgumentError(
        "coupled line fit is missing required top-level fields: " *
        join(sort!(collect(missing_keys)), ','),
    ))
    isempty(unknown_keys) || throw(ArgumentError(
        "coupled line fit contains unsupported top-level fields: " *
        join(sort!(collect(unknown_keys)), ','),
    ))
    get(data, "schema", "") == "aimora.coupled_line_fit.v1" || throw(ArgumentError(
        "unsupported coupled line fit schema",
    ))
    Int(data["schema_version"]) == COUPLED_LINE_FIT_SCHEMA_VERSION ||
        throw(ArgumentError("unsupported coupled line fit schema version"))
    source_signature = _checked_signature(
        data["source_signature_sha256"],
        "fit source signature",
    )
    response_signature = _checked_signature(
        data["response_signature_sha256"],
        "fit response signature",
    )
    if expected_source_signature_sha256 !== nothing
        source_signature == _checked_signature(
            expected_source_signature_sha256,
            "expected fit source signature",
        ) || throw(ArgumentError("coupled line fit source signature is stale"))
    end
    if expected_response_signature_sha256 !== nothing
        response_signature == _checked_signature(
            expected_response_signature_sha256,
            "expected fit response signature",
        ) || throw(ArgumentError("coupled line fit response signature is stale"))
    end
    model_data = data["model"]
    model = coupled_line_rational_model_from_state_space(
        _real_matrix_from_data(model_data["state_matrix_per_s"], "fit state matrix"),
        _real_matrix_from_data(model_data["input_matrix"], "fit input matrix"),
        _real_matrix_from_data(model_data["output_matrix_per_s"], "fit output matrix"),
        _real_matrix_from_data(model_data["direct_term"], "fit direct term"),
        ;
        port_order=Symbol.(model_data["port_order"]),
        reference_impedance_ohm=Float64.(model_data["reference_impedance_ohm"]),
    )
    model.deterministic_signature_sha256 == model_data["deterministic_signature_sha256"] ||
        throw(ArgumentError("coupled line fit model signature mismatch"))
    certificate_before = _certificate_from_data(data["certificate_before_enforcement"])
    certificate_after = _certificate_from_data(data["certificate_after_enforcement"])
    recertification_source_frequencies = Float64.(
        data["source_response"]["frequencies_hz"],
    )
    settings = _fit_settings_from_data(data["settings"])
    settings.deterministic_signature_sha256 == data["settings_signature_sha256"] ||
        throw(ArgumentError("fit settings identity disagrees with stored settings"))
    recertified = coupled_line_passivity_certificate(
        model;
        frequencies_hz=recertification_source_frequencies,
        passivity_tolerance=settings.passivity_tolerance,
        stable_pole_margin_per_s=settings.stable_pole_margin_per_s,
    )
    for field in (
        :method,
        :stable_realization,
        :direct_term_bounded,
        :zero_frequency_bounded,
        :hamiltonian_crossing_free,
        :continuous_passivity_passed,
    )
        getfield(recertified, field) == getfield(certificate_after, field) ||
            throw(ArgumentError(
                "coupled line fit passivity recertification disagrees in $(field)",
            ))
    end
    for field in (
        :gamma,
        :maximum_real_pole_per_s,
        :direct_maximum_singular_value,
        :zero_frequency_maximum_singular_value,
        :hamiltonian_minimum_real_separation_per_s,
        :hamiltonian_crossing_frequencies_hz,
        :worst_diagnostic_frequency_hz,
        :worst_diagnostic_singular_value,
        :minimum_physical_loss_eigenvalue_s,
        :diagnostic_frequencies_hz,
    )
        isapprox(getfield(recertified, field), getfield(certificate_after, field);
            rtol=2.0e-10, atol=2.0e-12) || throw(ArgumentError(
                "coupled line fit passivity recertification disagrees in $(field)",
            ))
    end
    length(recertified.worst_incident_direction) ==
        length(certificate_after.worst_incident_direction) || throw(ArgumentError(
            "coupled line fit passivity incident direction dimension changed",
        ))
    direction_overlap = abs(dot(
        recertified.worst_incident_direction,
        certificate_after.worst_incident_direction,
    ))
    isapprox(direction_overlap, 1.0; rtol=2.0e-10, atol=2.0e-12) ||
        throw(ArgumentError(
            "coupled line fit passivity incident direction changed",
        ))
    recertified.stable_realization && recertified.continuous_passivity_passed ||
        throw(ArgumentError("imported coupled line fit is not stable and globally passive"))
    enforcement = _enforcement_from_data(data["enforcement"])
    attempts = CoupledLineFitAttempt[
        CoupledLineFitAttempt(
            Int(row["order"]),
            Float64(row["maximum_relative_fit_error"]),
            Float64(row["maximum_weighted_condition_number"]),
            Int(row["effective_rank"]),
            Int(row["weighted_basis_column_count"]),
            Float64(row["effective_rank_fraction"]),
            Float64(row["regularization"]),
            Symbol(row["pole_policy"]),
            Int(row["pole_relocation_sweeps"]),
            Bool(row["pole_relocation_converged"]),
            Symbol(row["pole_relocation_outcome"]),
            Float64(row["pole_relocation_maximum_relative_change"]),
            Float64(row["pole_relocation_condition_number"]),
            Int(row["pole_relocation_effective_rank"]),
            Bool(row["low_frequency_anchor_applied"]),
            Bool(row["passivity_passed"]),
            Float64(row["enforcement_relative_perturbation"]),
            Bool(row["accepted"]),
            Symbol(row["outcome"]),
        ) for row in data["attempts"]
    ]
    modal = haskey(data, "modal_preparation") ?
        _modal_from_data(data["modal_preparation"]) : nothing
    uncertainty = haskey(data, "uncertainty") ?
        _uncertainty_from_data(data["uncertainty"]) : nothing
    source_response = _frequency_response_from_data(data["source_response"])
    source_response.source_signature_sha256 == source_signature || throw(ArgumentError(
        "fit source response signature is stale",
    ))
    source_response.response_signature_sha256 == response_signature ||
        throw(ArgumentError("fit response identity disagrees with source response"))
    result = CoupledLineFitResult(
        Int(data["schema_version"]),
        source_signature,
        response_signature,
        _checked_signature(data["settings_signature_sha256"], "fit settings signature"),
        settings,
        modal,
        model,
        certificate_before,
        certificate_after,
        enforcement,
        attempts,
        _fit_error_from_data(data["fit_error_before_enforcement"]),
        _fit_error_from_data(data["fit_error_after_enforcement"]),
        source_response,
        [
            _complex_matrix_from_data(row, "fitted scattering")
            for row in data["fitted_scattering_matrices"]
        ],
        [
            _complex_matrix_from_data(row, "fitted terminal admittance")
            for row in data["fitted_terminal_admittance_matrices_s"]
        ],
        Float64(data["maximum_absolute_fit_error"]),
        Float64(data["maximum_relative_fit_error"]),
        uncertainty,
        _checked_signature(data["deterministic_signature_sha256"], "fit result signature"),
    )
    length(result.fitted_scattering_matrices) ==
        length(result.source_response.frequencies_hz) ==
        length(result.fitted_terminal_admittance_matrices_s) || throw(ArgumentError(
            "fit source and fitted response rows do not align",
        ))
    recomputed_scattering = [
        coupled_line_model_value(result.model, 2.0im * pi * frequency)
        for frequency in result.source_response.frequencies_hz
    ]
    all(eachindex(recomputed_scattering)) do index
        isapprox(recomputed_scattering[index], result.fitted_scattering_matrices[index];
            rtol=2.0e-12, atol=2.0e-13)
    end || throw(ArgumentError("stored fitted scattering disagrees with the realization"))
    recomputed_admittance = [
        coupled_line_scattering_to_admittance(
            matrix,
            result.model.reference_impedance_ohm,
        ) for matrix in recomputed_scattering
    ]
    all(eachindex(recomputed_admittance)) do index
        isapprox(recomputed_admittance[index], result.fitted_terminal_admittance_matrices_s[index];
            rtol=2.0e-12, atol=2.0e-13)
    end || throw(ArgumentError("stored fitted admittance disagrees with the realization"))
    for diagnostics in (
        result.fit_error_before_enforcement,
        result.fit_error_after_enforcement,
    )
        diagnostics.frequencies_hz == result.source_response.frequencies_hz ||
            throw(ArgumentError("stored fit-error frequencies disagree with the source"))
    end
    recomputed_after_error = coupled_line_fit_error_diagnostics(
        result.source_response.frequencies_hz,
        result.source_response.scattering_matrices,
        recomputed_scattering,
    )
    _fit_error_diagnostics_match(
        recomputed_after_error,
        result.fit_error_after_enforcement,
    ) || throw(ArgumentError(
        "stored post-enforcement fit errors disagree with the realization",
    ))
    recomputed_absolute = maximum((
        opnorm(recomputed_scattering[index] -
            result.source_response.scattering_matrices[index])
        for index in eachindex(recomputed_scattering)
    ); init=0.0)
    recomputed_relative = maximum((
        opnorm(recomputed_scattering[index] -
            result.source_response.scattering_matrices[index]) /
            max(opnorm(result.source_response.scattering_matrices[index]), 1.0e-12)
        for index in eachindex(recomputed_scattering)
    ); init=0.0)
    isapprox(recomputed_absolute, result.maximum_absolute_fit_error;
        rtol=2.0e-12, atol=2.0e-13) || throw(ArgumentError(
            "stored absolute fit error disagrees with the realization",
        ))
    isapprox(recomputed_relative, result.maximum_relative_fit_error;
        rtol=2.0e-12, atol=2.0e-13) || throw(ArgumentError(
            "stored relative fit error disagrees with the realization",
        ))
    expected_signature = _coupled_line_fit_signature(
        result.response_signature_sha256,
        result.settings_signature_sha256,
        result.modal_preparation,
        result.model,
        result.certificate_before_enforcement,
        result.certificate_after_enforcement,
        result.enforcement,
        result.fit_error_before_enforcement,
        result.fit_error_after_enforcement,
        result.attempts,
        result.uncertainty,
    )
    expected_signature == result.deterministic_signature_sha256 ||
        throw(ArgumentError("coupled line fit result signature mismatch"))
    return result
end

function coupled_line_fit_report_text(result::CoupledLineFitResult)
    certificate = result.certificate_after_enforcement
    lines = String[
        "AIMORA coupled line fitting and global passivity report",
        "schema_version=$(result.schema_version)",
        "source_signature_sha256=$(result.source_signature_sha256)",
        "response_signature_sha256=$(result.response_signature_sha256)",
        "settings_signature_sha256=$(result.settings_signature_sha256)",
        "candidate_orders=$(join(result.settings.candidate_orders, ','))",
        "relative_fit_tolerance=$(result.settings.relative_fit_tolerance)",
        "maximum_weighted_condition_number=$(result.settings.maximum_weighted_condition_number)",
        "minimum_effective_rank_fraction=$(result.settings.minimum_effective_rank_fraction)",
        "maximum_direct_term_singular_value=$(result.settings.maximum_direct_term_singular_value)",
        "passivity_tolerance=$(result.settings.passivity_tolerance)",
        "maximum_enforcement_relative_perturbation=$(result.settings.maximum_enforcement_relative_perturbation)",
        "result_signature_sha256=$(result.deterministic_signature_sha256)",
        "segment_id=$(result.source_response.segment_id)",
        "segment_kind=$(result.source_response.segment_kind)",
        "source_scope=$(result.source_response.source_scope)",
        "port_order=$(join(result.model.port_order, ','))",
        "source_frequency_band_hz=$(first(result.source_response.frequencies_hz)),$(last(result.source_response.frequencies_hz))",
        "source_frequency_count=$(length(result.source_response.frequencies_hz))",
        "realization_pole_count=$(length(result.model.poles_per_s))",
        "selected_order=$(last(result.attempts).order)",
        "maximum_absolute_fit_error=$(result.maximum_absolute_fit_error)",
        "maximum_relative_fit_error=$(result.maximum_relative_fit_error)",
        "pre_enforcement_matrix_relative_error=$(result.fit_error_before_enforcement.maximum_relative_error)",
        "post_enforcement_matrix_relative_error=$(result.fit_error_after_enforcement.maximum_relative_error)",
        "post_enforcement_singular_value_relative_error=$(maximum(result.fit_error_after_enforcement.singular_value_relative_errors))",
        "post_enforcement_port_relative_error=$(maximum(result.fit_error_after_enforcement.port_relative_errors))",
        "stable_realization=$(certificate.stable_realization)",
        "continuous_global_passivity=$(certificate.continuous_passivity_passed)",
        "certificate_method=$(certificate.method)",
        "hamiltonian_crossing_frequencies_hz=$(join(certificate.hamiltonian_crossing_frequencies_hz, ','))",
        "worst_diagnostic_frequency_hz=$(certificate.worst_diagnostic_frequency_hz)",
        "worst_diagnostic_singular_value=$(certificate.worst_diagnostic_singular_value)",
        "minimum_physical_loss_eigenvalue_s=$(certificate.minimum_physical_loss_eigenvalue_s)",
        "passivity_enforcement_applied=$(result.enforcement.applied)",
        "passivity_enforcement_method=$(result.enforcement.method)",
        "passivity_enforcement_parameter=$(result.enforcement.enforcement_parameter)",
        "scattering_relative_perturbation=$(result.enforcement.scattering_maximum_relative_perturbation)",
        "admittance_relative_perturbation=$(result.enforcement.admittance_maximum_relative_perturbation)",
        "enforcement_within_budget=$(result.enforcement.within_budget)",
        "selected_weighted_condition_number=$(last(result.attempts).maximum_weighted_condition_number)",
        "selected_effective_rank=$(last(result.attempts).effective_rank)",
        "selected_weighted_basis_column_count=$(last(result.attempts).weighted_basis_column_count)",
        "selected_effective_rank_fraction=$(last(result.attempts).effective_rank_fraction)",
        "selected_regularization=$(last(result.attempts).regularization)",
        "selected_pole_policy=$(last(result.attempts).pole_policy)",
        "selected_pole_relocation_sweeps=$(last(result.attempts).pole_relocation_sweeps)",
        "selected_pole_relocation_converged=$(last(result.attempts).pole_relocation_converged)",
        "selected_pole_relocation_outcome=$(last(result.attempts).pole_relocation_outcome)",
        "selected_pole_relocation_maximum_relative_change=$(last(result.attempts).pole_relocation_maximum_relative_change)",
        "selected_pole_relocation_condition_number=$(last(result.attempts).pole_relocation_condition_number)",
        "selected_pole_relocation_effective_rank=$(last(result.attempts).pole_relocation_effective_rank)",
        "low_frequency_anchor_applied=$(last(result.attempts).low_frequency_anchor_applied)",
        "runtime_executed=false",
        "ulm_file_compatibility_claimed=false",
        "atp_pscad_equivalence_claimed=false",
    ]
    if result.modal_preparation !== nothing
        append!(lines, [
            "modal_delays_s=$(join(result.modal_preparation.extracted_delays_s, ','))",
            "minimum_mode_overlap=$(minimum(result.modal_preparation.minimum_mode_overlaps))",
            "maximum_principal_angle_rad=$(maximum(result.modal_preparation.maximum_principal_angles_rad))",
        ])
    end
    if result.uncertainty === nothing ||
        isempty(result.uncertainty.alternative_source_signatures_sha256)
        push!(lines, "uncertainty_status=unknown")
    else
        append!(lines, [
            "uncertainty_status=$(result.uncertainty.complete_set ? "complete_declared_set" : "bounded_alternatives_not_complete")",
            "uncertainty_alternative_count=$(length(result.uncertainty.alternative_source_signatures_sha256))",
            "uncertainty_source_scattering_relative_deviation=$(result.uncertainty.maximum_source_scattering_relative_deviation)",
            "uncertainty_fit_scattering_relative_deviation=$(result.uncertainty.maximum_fit_scattering_relative_deviation)",
            "uncertainty_minimum_passivity_margin=$(result.uncertainty.minimum_passivity_margin)",
            "uncertainty_maximum_delay_deviation_s=$(result.uncertainty.maximum_delay_deviation_s)",
            "uncertainty_pole_association_status=$(result.uncertainty.pole_association_status)",
            "uncertainty_mode_association_status=$(result.uncertainty.mode_association_status)",
        ])
    end
    return join(lines, '\n') * "\n"
end

function write_coupled_line_fit_report(path::AbstractString, result::CoupledLineFitResult)
    write(path, coupled_line_fit_report_text(result))
    return path
end

end

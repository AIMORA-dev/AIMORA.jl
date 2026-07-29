export FaultSequenceBranch,
       FaultGeneratorEquivalentFitResult,
       fault_sequence_admittance_matrix,
       fit_fault_generator_equivalent,
       fitted_generator_equivalent_history_injection

struct FaultSequenceBranch
    from_node::Int
    to_node::Int
    reactance_ohm::Float64

    function FaultSequenceBranch(
        from_node::Integer,
        to_node::Integer,
        reactance_ohm::Real,
    )
        from = Int(from_node)
        to = Int(to_node)
        from >= 0 && to >= 0 && from != to ||
            throw(ArgumentError("fault-sequence branch terminals must be distinct nonnegative indices"))
        reactance = Float64(reactance_ohm)
        isfinite(reactance) && reactance > 0.0 ||
            throw(ArgumentError("fault-sequence branch reactance must be finite and positive"))
        return new(from, to, reactance)
    end
end

struct FaultGeneratorEquivalentFitResult
    sequence_kind::Symbol
    frequency_hz::Float64
    network_admittance_s::Matrix{Float64}
    target_thevenin_reactance_ohm::Vector{Float64}
    target_x_over_r::Vector{Float64}
    generator_admittance_s::Vector{Float64}
    generator_reactance_ohm::Vector{Float64}
    generator_resistance_ohm::Vector{Float64}
    fitted_modal_branches::Vector{GeneratorEquivalentModalBranch}
    reconstructed_thevenin_reactance_ohm::Vector{Float64}
    residual_ohm::Vector{Float64}
    maximum_residual_ohm::Float64
    iteration_count::Int
    continuation_levels::Vector{Float64}
    correction_history::Vector{Float64}
    converged::Bool
    reciprocal::Bool
    passive::Bool
    physical_checks_passed::Bool
end

function fault_sequence_admittance_matrix(
    node_count::Integer,
    branches::AbstractVector{FaultSequenceBranch},
)
    count = Int(node_count)
    count > 0 ||
        throw(ArgumentError("fault-sequence node_count must be positive"))
    !isempty(branches) ||
        throw(ArgumentError("fault-sequence network requires at least one branch"))
    admittance = zeros(Float64, count, count)
    for branch in branches
        branch.from_node <= count && branch.to_node <= count ||
            throw(ArgumentError("fault-sequence branch terminal exceeds node_count"))
        susceptance = inv(branch.reactance_ohm)
        if branch.from_node != 0
            admittance[branch.from_node, branch.from_node] += susceptance
        end
        if branch.to_node != 0
            admittance[branch.to_node, branch.to_node] += susceptance
        end
        if branch.from_node != 0 && branch.to_node != 0
            admittance[branch.from_node, branch.to_node] -= susceptance
            admittance[branch.to_node, branch.from_node] -= susceptance
        end
    end
    return admittance
end

function _fault_equivalent_checked_admittance(
    network_admittance_s::AbstractMatrix{<:Real},
)
    rows, columns = size(network_admittance_s)
    rows == columns > 0 ||
        throw(ArgumentError("fault-equivalent network admittance must be nonempty and square"))
    admittance = Matrix{Float64}(network_admittance_s)
    all(isfinite, admittance) ||
        throw(ArgumentError("fault-equivalent network admittance must be finite"))
    symmetry_error = maximum(abs, admittance - transpose(admittance); init = 0.0)
    scale = max(maximum(abs, admittance; init = 0.0), 1.0)
    symmetry_error <= 1.0e-11 * scale ||
        throw(ArgumentError("fault-equivalent network admittance must be reciprocal"))
    minimum(eigvals(Symmetric(admittance))) >= -1.0e-10 * scale ||
        throw(ArgumentError("fault-equivalent network admittance must be positive semidefinite"))
    return admittance
end

function _fault_equivalent_driving_reactance(
    network_admittance_s::Matrix{Float64},
    generator_admittance_s::Vector{Float64},
)
    total = network_admittance_s + Diagonal(generator_admittance_s)
    factor = cholesky(Symmetric(total); check = true)
    impedance = Matrix(factor \ I)
    driving = diag(impedance)
    all(value -> isfinite(value) && value > 0.0, driving) ||
        throw(ArgumentError("fault-equivalent driving reactances must be finite and positive"))
    return driving, impedance
end

function _fault_equivalent_level_solve(
    network_admittance_s::Matrix{Float64},
    target_admittance_s::Vector{Float64},
    initial_generator_admittance_s::Vector{Float64};
    absolute_tolerance_ohm::Float64,
    relative_tolerance::Float64,
    maximum_iterations::Int,
)
    generator_admittance = copy(initial_generator_admittance_s)
    correction_history = Float64[]
    completed_iterations = 0
    converged = false
    for iteration in 1:maximum_iterations
        completed_iterations = iteration
        driving_reactance, impedance =
            _fault_equivalent_driving_reactance(
                network_admittance_s,
                generator_admittance,
            )
        driving_admittance = inv.(driving_reactance)
        residual = driving_admittance - target_admittance_s
        target_reactance = inv.(target_admittance_s)
        reconstructed_reactance = inv.(driving_admittance)
        reactance_residual = reconstructed_reactance - target_reactance
        allowance =
            absolute_tolerance_ohm +
            relative_tolerance *
            max(maximum(abs, target_reactance; init = 0.0), 1.0)
        if maximum(abs, reactance_residual; init = 0.0) <= allowance
            converged = true
            break
        end
        diagonal_impedance = diag(impedance)
        jacobian = Matrix{Float64}(
            undef,
            length(generator_admittance),
            length(generator_admittance),
        )
        for row in axes(jacobian, 1), column in axes(jacobian, 2)
            jacobian[row, column] =
                impedance[row, column]^2 / diagonal_impedance[row]^2
        end
        correction = try
            -(jacobian \ residual)
        catch error
            error isa SingularException || rethrow()
            -(pinv(jacobian) * residual)
        end
        all(isfinite, correction) ||
            throw(ArgumentError("fault-equivalent Newton correction is nonfinite"))
        correction_norm = maximum(abs, correction; init = 0.0)
        push!(correction_history, correction_norm)
        residual_norm = maximum(abs, residual; init = 0.0)
        accepted = false
        for damping in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625)
            candidate = generator_admittance + damping * correction
            all(>(0.0), candidate) || continue
            candidate_reactance, _ =
                _fault_equivalent_driving_reactance(
                    network_admittance_s,
                    candidate,
                )
            candidate_residual =
                inv.(candidate_reactance) - target_admittance_s
            if maximum(abs, candidate_residual; init = 0.0) < residual_norm
                generator_admittance = candidate
                accepted = true
                break
            end
        end
        accepted ||
            throw(ArgumentError("fault-equivalent Newton step did not reduce its residual"))
    end
    converged ||
        throw(ArgumentError(
            "fault-equivalent fit did not converge after $maximum_iterations iterations",
        ))
    return generator_admittance, completed_iterations, correction_history
end

"""
    fit_fault_generator_equivalent(network_admittance_s, target_x_ohm; ...)

Fit one passive generator shunt per sequence-network node so the complete
network reproduces the requested driving-point fault reactances. Continuation
from an uncoupled network and an analytic Newton Jacobian replace the legacy
`NTOT^2` overlay workspace while preserving the same physical equations.
"""
function fit_fault_generator_equivalent(
    network_admittance_s::AbstractMatrix{<:Real},
    target_thevenin_reactance_ohm::AbstractVector{<:Real};
    sequence_kind::Symbol = :positive,
    frequency_hz::Real = 60.0,
    target_x_over_r::AbstractVector{<:Real} =
        fill(Inf, length(target_thevenin_reactance_ohm)),
    continuation_fraction::Real = 0.1,
    absolute_tolerance_ohm::Real = 1.0e-10,
    relative_tolerance::Real = 1.0e-10,
    maximum_iterations::Integer = 40,
)
    sequence_kind in (:zero, :positive, :negative) ||
        throw(ArgumentError("fault-equivalent sequence_kind must be zero, positive, or negative"))
    admittance = _fault_equivalent_checked_admittance(network_admittance_s)
    target_reactance = Float64.(target_thevenin_reactance_ohm)
    length(target_reactance) == size(admittance, 1) ||
        throw(ArgumentError("fault-equivalent target count must match the network"))
    all(value -> isfinite(value) && value > 0.0, target_reactance) ||
        throw(ArgumentError("fault-equivalent target reactances must be finite and positive"))
    ratios = Float64.(target_x_over_r)
    length(ratios) == length(target_reactance) ||
        throw(ArgumentError("fault-equivalent X/R count must match the network"))
    all(value -> value > 0.0, ratios) ||
        throw(ArgumentError("fault-equivalent X/R values must be positive"))
    frequency = Float64(frequency_hz)
    isfinite(frequency) && frequency > 0.0 ||
        throw(ArgumentError("fault-equivalent frequency_hz must be finite and positive"))
    continuation = Float64(continuation_fraction)
    isfinite(continuation) && 0.0 < continuation <= 1.0 ||
        throw(ArgumentError("fault-equivalent continuation_fraction must be in (0, 1]"))
    absolute_tolerance = Float64(absolute_tolerance_ohm)
    relative_limit = Float64(relative_tolerance)
    all(value -> isfinite(value) && value > 0.0, (absolute_tolerance, relative_limit)) ||
        throw(ArgumentError("fault-equivalent tolerances must be finite and positive"))
    iteration_limit = Int(maximum_iterations)
    iteration_limit > 0 ||
        throw(ArgumentError("fault-equivalent maximum_iterations must be positive"))

    levels = Float64[]
    level = continuation
    while level < 1.0
        push!(levels, level)
        level = min(1.0, level + continuation)
    end
    if isempty(levels) || last(levels) < 1.0
        push!(levels, 1.0)
    end
    generator_admittance = inv.(target_reactance)
    total_iterations = 0
    correction_history = Float64[]
    target_admittance = inv.(target_reactance)
    for level in levels
        generator_admittance, iterations, corrections =
            _fault_equivalent_level_solve(
                level .* admittance,
                target_admittance,
                generator_admittance;
                absolute_tolerance_ohm = absolute_tolerance,
                relative_tolerance = relative_limit,
                maximum_iterations = iteration_limit,
            )
        total_iterations += iterations
        append!(correction_history, corrections)
    end

    reconstructed, _ =
        _fault_equivalent_driving_reactance(
            admittance,
            generator_admittance,
        )
    residual = reconstructed - target_reactance
    maximum_residual = maximum(abs, residual; init = 0.0)
    generator_reactance = inv.(generator_admittance)
    generator_resistance = Float64[
        isfinite(ratios[index]) ?
        generator_reactance[index] / ratios[index] :
        0.0
        for index in eachindex(ratios)
    ]
    angular_frequency = 2.0 * pi * frequency
    modal_branches = GeneratorEquivalentModalBranch[
        GeneratorEquivalentModalBranch(
            generator_resistance[index],
            generator_reactance[index] / angular_frequency,
            0.0,
            0.0,
        )
        for index in eachindex(generator_reactance)
    ]
    symmetry_error = maximum(abs, admittance - transpose(admittance); init = 0.0)
    reciprocal =
        symmetry_error <=
        1.0e-11 * max(maximum(abs, admittance; init = 0.0), 1.0)
    passive =
        all(>(0.0), generator_admittance) &&
        all(>=(0.0), generator_resistance) &&
        minimum(eigvals(Symmetric(admittance))) >=
        -1.0e-10 * max(maximum(abs, admittance; init = 0.0), 1.0)
    allowance =
        absolute_tolerance +
        relative_limit * max(maximum(abs, target_reactance; init = 0.0), 1.0)
    converged = maximum_residual <= allowance
    return FaultGeneratorEquivalentFitResult(
        sequence_kind,
        frequency,
        admittance,
        target_reactance,
        ratios,
        generator_admittance,
        generator_reactance,
        generator_resistance,
        modal_branches,
        reconstructed,
        residual,
        maximum_residual,
        total_iterations,
        levels,
        correction_history,
        converged,
        reciprocal,
        passive,
        converged && reciprocal && passive,
    )
end

function fit_fault_generator_equivalent(
    node_count::Integer,
    branches::AbstractVector{FaultSequenceBranch},
    target_thevenin_reactance_ohm::AbstractVector{<:Real};
    kwargs...,
)
    return fit_fault_generator_equivalent(
        fault_sequence_admittance_matrix(node_count, branches),
        target_thevenin_reactance_ohm;
        kwargs...,
    )
end

function fitted_generator_equivalent_history_injection(
    from_nodes::AbstractVector{<:Integer},
    to_nodes::AbstractVector{<:Integer},
    zero_sequence_fit::FaultGeneratorEquivalentFitResult,
    positive_sequence_fit::FaultGeneratorEquivalentFitResult;
    generator_index::Integer = 1,
    initial_phase_voltage::AbstractVector{<:Real} =
        zeros(Float64, length(from_nodes)),
)
    zero_sequence_fit.sequence_kind == :zero ||
        throw(ArgumentError("zero_sequence_fit must own the zero-sequence fit"))
    positive_sequence_fit.sequence_kind in (:positive, :negative) ||
        throw(ArgumentError("positive_sequence_fit must own a nonzero-sequence fit"))
    zero_sequence_fit.physical_checks_passed &&
        positive_sequence_fit.physical_checks_passed ||
        throw(ArgumentError("generator-equivalent fits must pass physical checks"))
    index = Int(generator_index)
    1 <= index <= length(zero_sequence_fit.fitted_modal_branches) &&
        index <= length(positive_sequence_fit.fitted_modal_branches) ||
        throw(ArgumentError("generator_index is outside the fitted modal branches"))
    return generator_equivalent_history_injection(
        from_nodes,
        to_nodes,
        [zero_sequence_fit.fitted_modal_branches[index]],
        [positive_sequence_fit.fitted_modal_branches[index]];
        initial_phase_voltage,
    )
end

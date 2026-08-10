module NonlinearNetwork

using ..StudyCore: ParameterNature,
                   PhysicalModelParameter,
                   ScalingBasisParameter,
                   NumericalPolicyParameter,
                   ParameterProvenance

export AbstractNonlinearCurrentDevice,
       NonlinearDeviceFormulation,
       NonlinearParameterNature,
       PhysicalModelParameter,
       ScalingBasisParameter,
       NumericalPolicyParameter,
       NonlinearParameterProvenance,
       PhysicalConstitutiveCurrent,
       ComplementarityCurrent,
       SemismoothCurrent,
       BoundedRegularizedCurrent,
       ExponentialCurrentBranch,
       CubicCurrentBranch,
       IdealVoltageConstraint,
       NonlinearNetworkScales,
       NonlinearSolveOptions,
       NonlinearChatterOptions,
       NonlinearChatterObservation,
       NonlinearChatterDecision,
       NonlinearResidualDiagnostic,
       NonlinearPowerDiagnostic,
       NonlinearSolveDiagnostics,
       NonlinearSolveFailure,
       NonlinearSolveResult,
       nonlinear_device_formulation,
       nonlinear_device_provenance,
       nonlinear_terminal_nodes,
       nonlinear_current_jacobian!,
       accept_nonlinear_device_state!,
       classify_numerical_chatter,
       validate_nonlinear_solve_options,
       nonlinear_result_accepted

"""Physical device contract for instantaneous terminal-current and analytic-Jacobian evaluation.

Implementations must evaluate trial states without mutating accepted device state. Terminal currents are positive leaving each declared terminal, terminal voltages are relative to ground, and the Jacobian entry `(row, column)` is the derivative of terminal current `row` with respect to terminal voltage `column` in siemens.
"""
abstract type AbstractNonlinearCurrentDevice end

const NonlinearParameterNature = ParameterNature
const NonlinearParameterProvenance = ParameterProvenance

function _caller_physical_parameter_provenance(
    units::AbstractString,
    validity_domain::AbstractString,
)
    return NonlinearParameterProvenance(
        "caller-supplied constructor values",
        units,
        "converted to Float64 without an implicit physical transformation",
        "not supplied; caller retains uncertainty ownership",
        validity_domain,
        PhysicalModelParameter,
    )
end

function _caller_scaling_provenance()
    return NonlinearParameterProvenance(
        "caller-declared physical scaling bases",
        "volt and ampere",
        "converted to positive Float64 SI magnitudes",
        "not supplied; scales are numerical bases rather than measured estimates",
        "declared nodes and ideal constraints in this nonlinear network",
        ScalingBasisParameter,
    )
end

function _aimora_numerical_policy_provenance(policy::AbstractString)
    return NonlinearParameterProvenance(
        "AIMORA documented nonlinear-network numerical policy",
        "field-specific SI or dimensionless units",
        "direct typed policy values; no conversion into device physics",
        "deterministic numerical limits, not physical uncertainty",
        policy,
        NumericalPolicyParameter,
    )
end

"""Declared mathematical formulation used by a nonlinear terminal-current device."""
@enum NonlinearDeviceFormulation begin
    PhysicalConstitutiveCurrent
    ComplementarityCurrent
    SemismoothCurrent
    BoundedRegularizedCurrent
end

"""Return the device's declared physical or numerical nonlinear formulation."""
function nonlinear_device_formulation(device::AbstractNonlinearCurrentDevice)
    throw(MethodError(nonlinear_device_formulation, (device,)))
end

"""Return the scientific provenance of a nonlinear device's physical or numerical parameters."""
function nonlinear_device_provenance(device::AbstractNonlinearCurrentDevice)
    throw(MethodError(nonlinear_device_provenance, (device,)))
end

"""Return the ordered physical terminal nodes of a nonlinear current device; node zero is ground."""
function nonlinear_terminal_nodes(device::AbstractNonlinearCurrentDevice)
    throw(MethodError(nonlinear_terminal_nodes, (device,)))
end

"""Fill terminal currents in amperes and their analytic voltage Jacobian in siemens without accepting state."""
function nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    device::AbstractNonlinearCurrentDevice,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::Float64,
)
    throw(MethodError(
        nonlinear_current_jacobian!,
        (terminal_current_a, terminal_jacobian_s, device, terminal_voltage_v, time_s),
    ))
end

"""Accept device state only after the coupled network solution has converged; stateless devices use the default no-op."""
function accept_nonlinear_device_state!(
    device::AbstractNonlinearCurrentDevice,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    time_s::Float64,
)
    return nothing
end

function _validate_two_terminal_nodes(positive_node::Integer, negative_node::Integer)
    positive = Int(positive_node)
    negative = Int(negative_node)
    positive >= 0 || throw(ArgumentError("positive terminal node must be nonnegative"))
    negative >= 0 || throw(ArgumentError("negative terminal node must be nonnegative"))
    positive != negative || throw(ArgumentError("nonlinear branch terminals must be distinct"))
    return positive, negative
end

"""Passive two-terminal exponential current law `i = g*v + Is*expm1(v/Vs)` from positive to negative terminal."""
struct ExponentialCurrentBranch <: AbstractNonlinearCurrentDevice
    positive_node::Int
    negative_node::Int
    saturation_current_a::Float64
    voltage_scale_v::Float64
    parallel_conductance_s::Float64
    provenance::NonlinearParameterProvenance

    function ExponentialCurrentBranch(
        positive_node::Integer,
        negative_node::Integer;
        saturation_current_a::Real,
        voltage_scale_v::Real,
        parallel_conductance_s::Real=0.0,
        provenance::NonlinearParameterProvenance=
            _caller_physical_parameter_provenance(
                "ampere, volt, and siemens",
                "finite positive saturation current/voltage scale and nonnegative parallel conductance",
            ),
    )
        positive, negative = _validate_two_terminal_nodes(positive_node, negative_node)
        saturation_current = Float64(saturation_current_a)
        voltage_scale = Float64(voltage_scale_v)
        parallel_conductance = Float64(parallel_conductance_s)
        isfinite(saturation_current) && saturation_current > 0.0 || throw(ArgumentError(
            "exponential saturation current must be finite and positive",
        ))
        isfinite(voltage_scale) && voltage_scale > 0.0 || throw(ArgumentError(
            "exponential voltage scale must be finite and positive",
        ))
        isfinite(parallel_conductance) && parallel_conductance >= 0.0 || throw(ArgumentError(
            "exponential parallel conductance must be finite and nonnegative",
        ))
        provenance.nature === PhysicalModelParameter || throw(ArgumentError(
            "exponential branch provenance must describe physical model parameters",
        ))
        return new(
            positive,
            negative,
            saturation_current,
            voltage_scale,
            parallel_conductance,
            provenance,
        )
    end
end

nonlinear_terminal_nodes(device::ExponentialCurrentBranch) =
    (device.positive_node, device.negative_node)

nonlinear_device_formulation(::ExponentialCurrentBranch) =
    PhysicalConstitutiveCurrent

nonlinear_device_provenance(device::ExponentialCurrentBranch) = device.provenance

function nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    device::ExponentialCurrentBranch,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::Float64,
)
    length(terminal_current_a) >= 2 || throw(DimensionMismatch(
        "exponential branch current workspace must contain two terminals",
    ))
    size(terminal_jacobian_s, 1) >= 2 && size(terminal_jacobian_s, 2) >= 2 ||
        throw(DimensionMismatch("exponential branch Jacobian workspace must be at least 2x2"))
    length(terminal_voltage_v) >= 2 || throw(DimensionMismatch(
        "exponential branch voltage workspace must contain two terminals",
    ))
    isfinite(time_s) || throw(ArgumentError("nonlinear evaluation time must be finite"))
    branch_voltage = terminal_voltage_v[1] - terminal_voltage_v[2]
    normalized_voltage = branch_voltage / device.voltage_scale_v
    exponential_value = exp(normalized_voltage)
    branch_current = device.parallel_conductance_s * branch_voltage +
        device.saturation_current_a * (exponential_value - 1.0)
    differential_conductance = device.parallel_conductance_s +
        device.saturation_current_a * exponential_value / device.voltage_scale_v
    terminal_current_a[1] = branch_current
    terminal_current_a[2] = -branch_current
    terminal_jacobian_s[1, 1] = differential_conductance
    terminal_jacobian_s[1, 2] = -differential_conductance
    terminal_jacobian_s[2, 1] = -differential_conductance
    terminal_jacobian_s[2, 2] = differential_conductance
    return nothing
end

"""Passive two-terminal polynomial current law `i = g*v + k*v^3` from positive to negative terminal."""
struct CubicCurrentBranch <: AbstractNonlinearCurrentDevice
    positive_node::Int
    negative_node::Int
    linear_conductance_s::Float64
    cubic_coefficient_a_per_v3::Float64
    provenance::NonlinearParameterProvenance

    function CubicCurrentBranch(
        positive_node::Integer,
        negative_node::Integer;
        linear_conductance_s::Real,
        cubic_coefficient_a_per_v3::Real,
        provenance::NonlinearParameterProvenance=
            _caller_physical_parameter_provenance(
                "siemens and ampere per volt cubed",
                "finite nonnegative coefficients with at least one positive constitutive term",
            ),
    )
        positive, negative = _validate_two_terminal_nodes(positive_node, negative_node)
        linear_conductance = Float64(linear_conductance_s)
        cubic_coefficient = Float64(cubic_coefficient_a_per_v3)
        isfinite(linear_conductance) && linear_conductance >= 0.0 || throw(ArgumentError(
            "cubic branch linear conductance must be finite and nonnegative",
        ))
        isfinite(cubic_coefficient) && cubic_coefficient >= 0.0 || throw(ArgumentError(
            "cubic branch coefficient must be finite and nonnegative",
        ))
        linear_conductance > 0.0 || cubic_coefficient > 0.0 || throw(ArgumentError(
            "cubic branch must have a positive linear or cubic coefficient",
        ))
        provenance.nature === PhysicalModelParameter || throw(ArgumentError(
            "cubic branch provenance must describe physical model parameters",
        ))
        return new(
            positive,
            negative,
            linear_conductance,
            cubic_coefficient,
            provenance,
        )
    end
end

nonlinear_terminal_nodes(device::CubicCurrentBranch) =
    (device.positive_node, device.negative_node)

nonlinear_device_formulation(::CubicCurrentBranch) =
    PhysicalConstitutiveCurrent

nonlinear_device_provenance(device::CubicCurrentBranch) = device.provenance

function nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    device::CubicCurrentBranch,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::Float64,
)
    length(terminal_current_a) >= 2 || throw(DimensionMismatch(
        "cubic branch current workspace must contain two terminals",
    ))
    size(terminal_jacobian_s, 1) >= 2 && size(terminal_jacobian_s, 2) >= 2 ||
        throw(DimensionMismatch("cubic branch Jacobian workspace must be at least 2x2"))
    length(terminal_voltage_v) >= 2 || throw(DimensionMismatch(
        "cubic branch voltage workspace must contain two terminals",
    ))
    isfinite(time_s) || throw(ArgumentError("nonlinear evaluation time must be finite"))
    branch_voltage = terminal_voltage_v[1] - terminal_voltage_v[2]
    branch_voltage_squared = branch_voltage * branch_voltage
    branch_current = device.linear_conductance_s * branch_voltage +
        device.cubic_coefficient_a_per_v3 * branch_voltage_squared * branch_voltage
    differential_conductance = device.linear_conductance_s +
        3.0 * device.cubic_coefficient_a_per_v3 * branch_voltage_squared
    terminal_current_a[1] = branch_current
    terminal_current_a[2] = -branch_current
    terminal_jacobian_s[1, 1] = differential_conductance
    terminal_jacobian_s[1, 2] = -differential_conductance
    terminal_jacobian_s[2, 1] = -differential_conductance
    terminal_jacobian_s[2, 2] = differential_conductance
    return nothing
end

"""Linear ideal-voltage equation `sum(coefficients .* terminal_voltage) = value_v` with an MNA current unknown."""
struct IdealVoltageConstraint{N}
    terminal_nodes::NTuple{N,Int}
    terminal_coefficients::NTuple{N,Float64}
    value_v::Float64
    provenance::NonlinearParameterProvenance

    function IdealVoltageConstraint(
        terminal_nodes::NTuple{N,<:Integer},
        terminal_coefficients::NTuple{N,<:Real},
        value_v::Real,
        ;
        provenance::NonlinearParameterProvenance=
            _caller_physical_parameter_provenance(
                "volt and dimensionless terminal coefficients",
                "declared linear ideal-voltage source or closed ideal-switch constraint",
            ),
    ) where {N}
        N > 0 || throw(ArgumentError("ideal voltage constraint must have a terminal"))
        nodes = ntuple(index -> Int(terminal_nodes[index]), Val(N))
        coefficients = ntuple(index -> Float64(terminal_coefficients[index]), Val(N))
        all(node -> node >= 0, nodes) || throw(ArgumentError(
            "ideal voltage constraint nodes must be nonnegative",
        ))
        length(unique(nodes)) == N || throw(ArgumentError(
            "ideal voltage constraint terminal nodes must be unique",
        ))
        all(isfinite, coefficients) || throw(ArgumentError(
            "ideal voltage constraint coefficients must be finite",
        ))
        any(index -> nodes[index] != 0 && coefficients[index] != 0.0, 1:N) ||
            throw(ArgumentError("ideal voltage constraint must address a non-ground voltage"))
        value = Float64(value_v)
        isfinite(value) || throw(ArgumentError("ideal voltage constraint value must be finite"))
        provenance.nature === PhysicalModelParameter || throw(ArgumentError(
            "ideal voltage-constraint provenance must describe physical model parameters",
        ))
        return new{N}(nodes, coefficients, value, provenance)
    end
end

"""Declared physical magnitudes used to scale nodal voltage/current and ideal-constraint equations."""
struct NonlinearNetworkScales
    node_voltage_v::Vector{Float64}
    node_current_a::Vector{Float64}
    constraint_voltage_v::Vector{Float64}
    constraint_current_a::Vector{Float64}
    provenance::NonlinearParameterProvenance

    function NonlinearNetworkScales(
        node_voltage_v,
        node_current_a,
        constraint_voltage_v,
        constraint_current_a,
        ;
        provenance::NonlinearParameterProvenance=_caller_scaling_provenance(),
    )
        node_voltage = Float64.(node_voltage_v)
        node_current = Float64.(node_current_a)
        constraint_voltage = Float64.(constraint_voltage_v)
        constraint_current = Float64.(constraint_current_a)
        length(node_voltage) == length(node_current) || throw(DimensionMismatch(
            "node voltage/current scale lengths must match",
        ))
        length(constraint_voltage) == length(constraint_current) ||
            throw(DimensionMismatch("constraint voltage/current scale lengths must match"))
        for (name, values) in (
            ("node voltage", node_voltage),
            ("node current", node_current),
            ("constraint voltage", constraint_voltage),
            ("constraint current", constraint_current),
        )
            all(value -> isfinite(value) && value > 0.0, values) || throw(ArgumentError(
                "$name scales must be finite and positive",
            ))
        end
        provenance.nature === ScalingBasisParameter || throw(ArgumentError(
            "nonlinear network scale provenance must describe scaling bases",
        ))
        return new(
            node_voltage,
            node_current,
            constraint_voltage,
            constraint_current,
            provenance,
        )
    end
end

function NonlinearNetworkScales(
    node_count::Integer,
    constraint_count::Integer;
    nominal_voltage_v::Real=1.0,
    nominal_current_a::Real=1.0,
)
    nodes = Int(node_count)
    constraints = Int(constraint_count)
    nodes > 0 || throw(ArgumentError("nonlinear network must contain at least one node"))
    constraints >= 0 || throw(ArgumentError("constraint count must be nonnegative"))
    voltage = Float64(nominal_voltage_v)
    current = Float64(nominal_current_a)
    return NonlinearNetworkScales(
        fill(voltage, nodes),
        fill(current, nodes),
        fill(voltage, constraints),
        fill(current, constraints),
    )
end

"""Numerical policy for scaled safeguarded Newton solution of physical nonlinear KCL/MNA equations."""
Base.@kwdef struct NonlinearSolveOptions
    maximum_iterations::Int = 30
    current_absolute_tolerance_a::Float64 = 1.0e-12
    current_relative_tolerance::Float64 = 1.0e-10
    voltage_absolute_tolerance_v::Float64 = 1.0e-12
    voltage_relative_tolerance::Float64 = 1.0e-10
    scaled_step_tolerance::Float64 = 1.0e-10
    armijo_sufficient_decrease::Float64 = 1.0e-4
    minimum_step_fraction::Float64 = 9.5367431640625e-7
    maximum_condition_estimate::Float64 = 1.0e12
    rank_threshold_multiplier::Float64 = 1.0
    maximum_refinement_steps::Int = 3
    refinement_residual_multiplier::Float64 = 32.0
    maximum_regularized_trials::Int = 8
    initial_regularization::Float64 = 1.0e-10
    maximum_scaled_step::Float64 = 10.0
    sparse_dimension_threshold::Int = 64
    maximum_condition_estimation_steps::Int = 5
    provenance::NonlinearParameterProvenance = _aimora_numerical_policy_provenance(
        "scaled safeguarded Newton, factorization, condition estimation, and refinement",
    )
end

function validate_nonlinear_solve_options(options::NonlinearSolveOptions)
    options.maximum_iterations > 0 || throw(ArgumentError("maximum iterations must be positive"))
    isfinite(options.current_absolute_tolerance_a) &&
        options.current_absolute_tolerance_a > 0.0 || throw(ArgumentError(
        "current absolute tolerance must be finite and positive",
    ))
    isfinite(options.current_relative_tolerance) &&
        options.current_relative_tolerance >= 0.0 || throw(ArgumentError(
        "current relative tolerance must be finite and nonnegative",
    ))
    isfinite(options.voltage_absolute_tolerance_v) &&
        options.voltage_absolute_tolerance_v > 0.0 || throw(ArgumentError(
        "voltage absolute tolerance must be finite and positive",
    ))
    isfinite(options.voltage_relative_tolerance) &&
        options.voltage_relative_tolerance >= 0.0 || throw(ArgumentError(
        "voltage relative tolerance must be finite and nonnegative",
    ))
    isfinite(options.scaled_step_tolerance) &&
        options.scaled_step_tolerance > 0.0 || throw(ArgumentError(
        "scaled step tolerance must be finite and positive",
    ))
    0.0 < options.armijo_sufficient_decrease < 1.0 || throw(ArgumentError(
        "Armijo sufficient-decrease coefficient must lie strictly between zero and one",
    ))
    0.0 < options.minimum_step_fraction <= 1.0 || throw(ArgumentError(
        "minimum step fraction must lie in (0, 1]",
    ))
    isfinite(options.maximum_condition_estimate) &&
        options.maximum_condition_estimate >= 1.0 || throw(ArgumentError(
        "maximum condition estimate must be finite and at least one",
    ))
    isfinite(options.rank_threshold_multiplier) &&
        options.rank_threshold_multiplier > 0.0 || throw(ArgumentError(
        "rank threshold multiplier must be finite and positive",
    ))
    options.maximum_refinement_steps >= 0 || throw(ArgumentError(
        "maximum refinement steps must be nonnegative",
    ))
    options.refinement_residual_multiplier > 0.0 &&
        isfinite(options.refinement_residual_multiplier) || throw(ArgumentError(
            "refinement residual multiplier must be finite and positive",
        ))
    options.maximum_regularized_trials >= 0 || throw(ArgumentError(
        "maximum regularized trials must be nonnegative",
    ))
    isfinite(options.initial_regularization) &&
        options.initial_regularization > 0.0 || throw(ArgumentError(
        "initial regularization must be finite and positive",
    ))
    isfinite(options.maximum_scaled_step) &&
        options.maximum_scaled_step > 0.0 || throw(ArgumentError(
        "maximum scaled step must be finite and positive",
    ))
    options.sparse_dimension_threshold > 0 || throw(ArgumentError(
        "sparse factorization dimension threshold must be positive",
    ))
    options.maximum_condition_estimation_steps > 0 || throw(ArgumentError(
        "maximum condition-estimation steps must be positive",
    ))
    options.provenance.nature === NumericalPolicyParameter || throw(ArgumentError(
        "nonlinear solve-option provenance must describe numerical policy",
    ))
    return options
end

"""Numerical-chatter thresholds expressed against a declared physical channel scale."""
Base.@kwdef struct NonlinearChatterOptions
    minimum_increment_pu::Float64 = 1.0e-8
    adjacent_ratio_minimum::Float64 = 0.8
    adjacent_ratio_maximum::Float64 = 1.25
    provenance::NonlinearParameterProvenance = _aimora_numerical_policy_provenance(
        "trapezoidal-chatter discrimination and critical-damping veto policy",
    )
end

"""Three accepted signed increments and the physical/event facts needed to distinguish numerical chatter from resolved behavior."""
struct NonlinearChatterObservation
    accepted_increments::NTuple{3,Float64}
    physical_scale::Float64
    topology_unchanged::Bool
    task_calendar_unchanged::Bool
    control_mode_unchanged::Bool
    localized_event_present::Bool
    pwm_transition_present::Bool
    controller_transition_present::Bool
    resolved_physical_oscillation::Bool
    residual_or_energy_supports_physical_mode::Bool

    function NonlinearChatterObservation(
        accepted_increments::NTuple{3,<:Real},
        physical_scale::Real;
        topology_unchanged::Bool,
        task_calendar_unchanged::Bool,
        control_mode_unchanged::Bool,
        localized_event_present::Bool=false,
        pwm_transition_present::Bool=false,
        controller_transition_present::Bool=false,
        resolved_physical_oscillation::Bool=false,
        residual_or_energy_supports_physical_mode::Bool=false,
    )
        increments = ntuple(index -> Float64(accepted_increments[index]), Val(3))
        all(isfinite, increments) || throw(ArgumentError(
            "accepted chatter-observation increments must be finite",
        ))
        scale = Float64(physical_scale)
        isfinite(scale) && scale > 0.0 || throw(ArgumentError(
            "chatter-observation physical scale must be finite and positive",
        ))
        return new(
            increments,
            scale,
            topology_unchanged,
            task_calendar_unchanged,
            control_mode_unchanged,
            localized_event_present,
            pwm_transition_present,
            controller_transition_present,
            resolved_physical_oscillation,
            residual_or_energy_supports_physical_mode,
        )
    end
end

"""Auditable decision describing whether critical-damping treatment is permitted for an observed channel."""
struct NonlinearChatterDecision
    numerical_chatter::Bool
    critical_damping_allowed::Bool
    classification::Symbol
    maximum_increment_pu::Float64
    minimum_adjacent_ratio::Float64
    maximum_adjacent_ratio::Float64
end

function _validate_nonlinear_chatter_options(options::NonlinearChatterOptions)
    isfinite(options.minimum_increment_pu) && options.minimum_increment_pu > 0.0 ||
        throw(ArgumentError("minimum chatter increment must be finite and positive"))
    isfinite(options.adjacent_ratio_minimum) &&
        options.adjacent_ratio_minimum > 0.0 || throw(ArgumentError(
            "minimum adjacent chatter ratio must be finite and positive",
        ))
    isfinite(options.adjacent_ratio_maximum) &&
        options.adjacent_ratio_maximum >= options.adjacent_ratio_minimum ||
        throw(ArgumentError(
            "maximum adjacent chatter ratio must be finite and no smaller than the minimum",
        ))
    options.provenance.nature === NumericalPolicyParameter || throw(ArgumentError(
        "nonlinear chatter-option provenance must describe numerical policy",
    ))
    return options
end

function _chatter_decision(
    classification::Symbol,
    numerical_chatter::Bool,
    increments_pu::NTuple{3,Float64},
    adjacent_ratios::NTuple{2,Float64},
)
    return NonlinearChatterDecision(
        numerical_chatter,
        numerical_chatter,
        classification,
        maximum(abs, increments_pu),
        minimum(adjacent_ratios),
        maximum(adjacent_ratios),
    )
end

"""Classify a three-increment alternating pattern without damping physical oscillation, PWM, controller, event, or topology behavior."""
function classify_numerical_chatter(
    observation::NonlinearChatterObservation;
    options::NonlinearChatterOptions=NonlinearChatterOptions(),
)
    _validate_nonlinear_chatter_options(options)
    increments_pu = ntuple(
        index -> observation.accepted_increments[index] / observation.physical_scale,
        Val(3),
    )
    adjacent_ratios = (
        abs(increments_pu[2]) / max(abs(increments_pu[1]), eps(Float64)),
        abs(increments_pu[3]) / max(abs(increments_pu[2]), eps(Float64)),
    )
    alternating = increments_pu[1] * increments_pu[2] < 0.0 &&
        increments_pu[2] * increments_pu[3] < 0.0
    above_floor = all(
        increment -> abs(increment) >= options.minimum_increment_pu,
        increments_pu,
    )
    ratios_in_band = all(
        ratio -> options.adjacent_ratio_minimum <= ratio <= options.adjacent_ratio_maximum,
        adjacent_ratios,
    )
    if observation.localized_event_present
        return _chatter_decision(:localized_event_veto, false, increments_pu, adjacent_ratios)
    elseif observation.pwm_transition_present
        return _chatter_decision(:pwm_transition_veto, false, increments_pu, adjacent_ratios)
    elseif observation.controller_transition_present
        return _chatter_decision(
            :controller_transition_veto,
            false,
            increments_pu,
            adjacent_ratios,
        )
    elseif !observation.topology_unchanged
        return _chatter_decision(:topology_change_veto, false, increments_pu, adjacent_ratios)
    elseif !observation.task_calendar_unchanged
        return _chatter_decision(:task_transition_veto, false, increments_pu, adjacent_ratios)
    elseif !observation.control_mode_unchanged
        return _chatter_decision(:control_mode_veto, false, increments_pu, adjacent_ratios)
    end
    if !(alternating && above_floor && ratios_in_band)
        return _chatter_decision(
            :insufficient_chatter_evidence,
            false,
            increments_pu,
            adjacent_ratios,
        )
    end
    if observation.resolved_physical_oscillation ||
       observation.residual_or_energy_supports_physical_mode
        return _chatter_decision(
            :physical_oscillation_veto,
            false,
            increments_pu,
            adjacent_ratios,
        )
    end
    return _chatter_decision(:numerical_chatter, true, increments_pu, adjacent_ratios)
end

"""Typed numerical failure that leaves the last accepted network and device state unchanged."""
struct NonlinearSolveFailure <: Exception
    code::Symbol
    message::String
    iteration::Int
    condition_estimate::Float64
    numerical_rank::Int
    system_dimension::Int
    topology_signature::UInt64
end

function Base.showerror(io::IO, failure::NonlinearSolveFailure)
    print(io, failure.message, " [", failure.code, ", iteration=", failure.iteration,
        ", rank=", failure.numerical_rank, "/", failure.system_dimension,
        ", condition=", failure.condition_estimate, "]")
end

"""One physical and tolerance-normalized residual sample from a declared nonlinear substep and iteration."""
struct NonlinearResidualDiagnostic
    substep_index::Int
    iteration_index::Int
    maximum_kcl_residual_a::Float64
    maximum_constraint_residual_v::Float64
    maximum_tolerance_normalized_residual::Float64
end

"""Instantaneous power accounting at the final evaluated nonlinear candidate; a rejected result does not imply this candidate was accepted."""
struct NonlinearPowerDiagnostic
    nonlinear_device_absorbed_power_w::Float64
    ideal_constraint_absorbed_power_w::Float64
    algebraic_power_balance_residual_w::Float64
end

"""Physical, numerical, factorization, and discontinuity diagnostics for one nonlinear network solve."""
struct NonlinearSolveDiagnostics
    iteration_count::Int
    maximum_kcl_residual_a::Float64
    maximum_constraint_residual_v::Float64
    maximum_tolerance_normalized_residual::Float64
    residual_history::Vector{NonlinearResidualDiagnostic}
    power::NonlinearPowerDiagnostic
    maximum_scaled_step::Float64
    condition_estimate::Float64
    numerical_rank::Int
    system_dimension::Int
    line_search_backtrack_count::Int
    regularized_fallback_count::Int
    symbolic_factorization_count::Int
    numeric_factorization_count::Int
    factor_reuse_count::Int
    iterative_refinement_count::Int
    jacobian_nonzero_count::Int
    linear_solver::Symbol
    topology_signature::UInt64
    companion_method::Symbol
    accepted_substep_count::Int
    chatter_classification::Symbol
    discontinuity_reason::Symbol
end

"""Accepted nonlinear voltage/constraint-current state or a typed failure with unchanged accepted state."""
struct NonlinearSolveResult
    accepted::Bool
    voltage_v::Vector{Float64}
    constraint_current_a::Vector{Float64}
    diagnostics::NonlinearSolveDiagnostics
    failure::Union{Nothing,NonlinearSolveFailure}
end

nonlinear_result_accepted(result::NonlinearSolveResult) = result.accepted

end

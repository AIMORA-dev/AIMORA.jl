using AIMORA.NonlinearNetwork

function synthetic_physical_parameter_provenance()
    return NonlinearParameterProvenance(
        "AIMORA synthetic nonlinear-network test fixture",
        "siemens",
        "direct deterministic fixture parameter",
        "exact synthetic value without measurement uncertainty",
        "bounded algebraic test networks",
        PhysicalModelParameter,
    )
end

@testset "nonlinear terminal-current and ideal-constraint public contracts" begin
    exponential = ExponentialCurrentBranch(
        1,
        0;
        saturation_current_a=2.0e-3,
        voltage_scale_v=0.25,
        parallel_conductance_s=0.1,
    )
    @test nonlinear_device_formulation(exponential) === PhysicalConstitutiveCurrent
    @test nonlinear_device_provenance(exponential).nature === PhysicalModelParameter
    @test occursin("not supplied", nonlinear_device_provenance(exponential).uncertainty)
    terminal_current = zeros(Float64, 2)
    terminal_jacobian = zeros(Float64, 2, 2)
    terminal_voltage = [0.4, 0.0]
    nonlinear_current_jacobian!(
        terminal_current,
        terminal_jacobian,
        exponential,
        terminal_voltage,
        0.0,
    )
    perturbation = 1.0e-7
    positive_current = zeros(Float64, 2)
    negative_current = zeros(Float64, 2)
    scratch_jacobian = zeros(Float64, 2, 2)
    nonlinear_current_jacobian!(
        positive_current,
        scratch_jacobian,
        exponential,
        [terminal_voltage[1] + perturbation, 0.0],
        0.0,
    )
    nonlinear_current_jacobian!(
        negative_current,
        scratch_jacobian,
        exponential,
        [terminal_voltage[1] - perturbation, 0.0],
        0.0,
    )
    central_derivative = (positive_current[1] - negative_current[1]) /
        (2.0 * perturbation)
    @test terminal_current[1] == -terminal_current[2]
    @test terminal_jacobian[1, 1] ≈ central_derivative rtol=5.0e-9
    @test terminal_jacobian[1, 1] == terminal_jacobian[2, 2]
    @test terminal_jacobian[1, 2] == -terminal_jacobian[1, 1]

    cubic = CubicCurrentBranch(
        2,
        1;
        linear_conductance_s=0.2,
        cubic_coefficient_a_per_v3=0.1,
    )
    @test nonlinear_device_formulation(cubic) === PhysicalConstitutiveCurrent
    nonlinear_current_jacobian!(
        terminal_current,
        terminal_jacobian,
        cubic,
        [2.0, 1.0],
        0.0,
    )
    @test terminal_current ≈ [0.3, -0.3] atol=1.0e-16
    @test terminal_jacobian == [0.5 -0.5; -0.5 0.5]

    constraint = IdealVoltageConstraint((1, 2), (1.0, -1.0), 1.0)
    @test constraint.terminal_nodes == (1, 2)
    @test constraint.terminal_coefficients == (1.0, -1.0)
    @test constraint.provenance.nature === PhysicalModelParameter
    @test_throws ArgumentError IdealVoltageConstraint((1, 1), (1.0, -1.0), 0.0)
    @test_throws ArgumentError ExponentialCurrentBranch(
        1,
        0;
        saturation_current_a=0.0,
        voltage_scale_v=0.25,
    )
    @test_throws ArgumentError NonlinearParameterProvenance(
        " ",
        "ampere",
        "direct",
        "unknown",
        "positive current",
        PhysicalModelParameter,
    )
    scale_parameter_provenance = NonlinearParameterProvenance(
        "synthetic scale",
        "volt",
        "direct",
        "not applicable",
        "bounded test network",
        ScalingBasisParameter,
    )
    @test_throws ArgumentError IdealVoltageConstraint(
        (1, 2),
        (1.0, -1.0),
        0.0;
        provenance=scale_parameter_provenance,
    )
    @test_throws ArgumentError ExponentialCurrentBranch(
        1,
        0;
        saturation_current_a=1.0e-3,
        voltage_scale_v=0.25,
        provenance=scale_parameter_provenance,
    )
    @test validate_nonlinear_solve_options(NonlinearSolveOptions()) isa
        NonlinearSolveOptions
    for options in (
        NonlinearSolveOptions(current_absolute_tolerance_a=Inf),
        NonlinearSolveOptions(current_relative_tolerance=Inf),
        NonlinearSolveOptions(voltage_absolute_tolerance_v=Inf),
        NonlinearSolveOptions(voltage_relative_tolerance=Inf),
        NonlinearSolveOptions(scaled_step_tolerance=Inf),
        NonlinearSolveOptions(maximum_condition_estimate=Inf),
        NonlinearSolveOptions(rank_threshold_multiplier=Inf),
        NonlinearSolveOptions(initial_regularization=Inf),
        NonlinearSolveOptions(maximum_scaled_step=Inf),
        NonlinearSolveOptions(maximum_condition_estimation_steps=0),
    )
        @test_throws ArgumentError validate_nonlinear_solve_options(options)
    end
    @test_throws ArgumentError validate_nonlinear_solve_options(NonlinearSolveOptions(
        provenance=scale_parameter_provenance,
    ))
    @test_throws ArgumentError classify_numerical_chatter(
        NonlinearChatterObservation(
            (0.10, -0.09, 0.095),
            1.0;
            topology_unchanged=true,
            task_calendar_unchanged=true,
            control_mode_unchanged=true,
        );
        options=NonlinearChatterOptions(provenance=scale_parameter_provenance),
    )

    chatter_observation = NonlinearChatterObservation(
        (0.10, -0.09, 0.095),
        1.0;
        topology_unchanged=true,
        task_calendar_unchanged=true,
        control_mode_unchanged=true,
    )
    chatter_decision = classify_numerical_chatter(chatter_observation)
    @test chatter_decision.numerical_chatter
    @test chatter_decision.critical_damping_allowed
    @test chatter_decision.classification == :numerical_chatter
    pwm_decision = classify_numerical_chatter(NonlinearChatterObservation(
        chatter_observation.accepted_increments,
        chatter_observation.physical_scale;
        topology_unchanged=true,
        task_calendar_unchanged=true,
        control_mode_unchanged=true,
        pwm_transition_present=true,
    ))
    @test !pwm_decision.critical_damping_allowed
    @test pwm_decision.classification == :pwm_transition_veto
    physical_decision = classify_numerical_chatter(NonlinearChatterObservation(
        chatter_observation.accepted_increments,
        chatter_observation.physical_scale;
        topology_unchanged=true,
        task_calendar_unchanged=true,
        control_mode_unchanged=true,
        resolved_physical_oscillation=true,
    ))
    @test physical_decision.classification == :physical_oscillation_veto
end

if AIMORA.solver_available()
    using AIMORA.Branches
    using AIMORA.EMTStudy
    using AIMORA.Nodal
    using AIMORA.Nonlinear
    using AIMORA.NonlinearNodal
    using AIMORA.Switches

    struct DiagonalConditioningNetwork <: AIMORA.Branches.EMTElement
        conductance_s::Vector{Float64}
        source_current_a::Vector{Float64}
    end

    function AIMORA.Branches.stamp!(
        admittance,
        rhs,
        network::DiagonalConditioningNetwork,
        time_s::Float64,
        step_s::Float64,
    )
        isfinite(time_s) || throw(ArgumentError(
            "diagonal-conditioning time must be finite",
        ))
        isfinite(step_s) && step_s > 0.0 || throw(ArgumentError(
            "diagonal-conditioning step must be finite and positive",
        ))
        @inbounds for node in eachindex(network.conductance_s, network.source_current_a)
            admittance[node, node] += network.conductance_s[node]
            rhs[node] += network.source_current_a[node]
        end
        return nothing
    end

    AIMORA.Branches.update!(
        network::DiagonalConditioningNetwork,
        voltage,
        step_s::Float64,
    ) = nothing

    mutable struct TimeLimitedLinearCurrentBranch <: AbstractNonlinearCurrentDevice
        positive_node::Int
        negative_node::Int
        failure_time_s::Float64
        accepted_state_count::Int
    end

    AIMORA.NonlinearNetwork.nonlinear_terminal_nodes(
        device::TimeLimitedLinearCurrentBranch,
    ) = (device.positive_node, device.negative_node)

    AIMORA.NonlinearNetwork.nonlinear_device_formulation(
        ::TimeLimitedLinearCurrentBranch,
    ) = PhysicalConstitutiveCurrent

    AIMORA.NonlinearNetwork.nonlinear_device_provenance(
        ::TimeLimitedLinearCurrentBranch,
    ) = synthetic_physical_parameter_provenance()

    function AIMORA.NonlinearNetwork.nonlinear_current_jacobian!(
        terminal_current_a::AbstractVector{Float64},
        terminal_jacobian_s::AbstractMatrix{Float64},
        device::TimeLimitedLinearCurrentBranch,
        terminal_voltage_v::AbstractVector{Float64},
        time_s::Float64,
    )
        conductance = time_s >= device.failure_time_s ? NaN : 0.1
        branch_current = conductance * (terminal_voltage_v[1] - terminal_voltage_v[2])
        terminal_current_a[1] = branch_current
        terminal_current_a[2] = -branch_current
        terminal_jacobian_s[1, 1] = conductance
        terminal_jacobian_s[1, 2] = -conductance
        terminal_jacobian_s[2, 1] = -conductance
        terminal_jacobian_s[2, 2] = conductance
        return nothing
    end

    function AIMORA.NonlinearNetwork.accept_nonlinear_device_state!(
        device::TimeLimitedLinearCurrentBranch,
        terminal_voltage_v::AbstractVector{Float64},
        terminal_current_a::AbstractVector{Float64},
        time_s::Float64,
    )
        device.accepted_state_count += 1
        return nothing
    end

    mutable struct FailingAcceptanceCurrentBranch <: AbstractNonlinearCurrentDevice
        positive_node::Int
        negative_node::Int
        accepted_state_count::Int
    end

    AIMORA.NonlinearNetwork.nonlinear_terminal_nodes(
        device::FailingAcceptanceCurrentBranch,
    ) = (device.positive_node, device.negative_node)

    AIMORA.NonlinearNetwork.nonlinear_device_formulation(
        ::FailingAcceptanceCurrentBranch,
    ) = PhysicalConstitutiveCurrent

    AIMORA.NonlinearNetwork.nonlinear_device_provenance(
        ::FailingAcceptanceCurrentBranch,
    ) = synthetic_physical_parameter_provenance()

    function AIMORA.NonlinearNetwork.nonlinear_current_jacobian!(
        terminal_current_a::AbstractVector{Float64},
        terminal_jacobian_s::AbstractMatrix{Float64},
        device::FailingAcceptanceCurrentBranch,
        terminal_voltage_v::AbstractVector{Float64},
        time_s::Float64,
    )
        conductance_s = 0.1
        branch_current_a = conductance_s *
            (terminal_voltage_v[1] - terminal_voltage_v[2])
        terminal_current_a .= (branch_current_a, -branch_current_a)
        terminal_jacobian_s .= [
            conductance_s -conductance_s
            -conductance_s conductance_s
        ]
        return nothing
    end

    function AIMORA.NonlinearNetwork.accept_nonlinear_device_state!(
        device::FailingAcceptanceCurrentBranch,
        terminal_voltage_v::AbstractVector{Float64},
        terminal_current_a::AbstractVector{Float64},
        time_s::Float64,
    )
        device.accepted_state_count += 1
        error("deliberate nonlinear device acceptance failure")
    end

    struct DuplicateTerminalCurrentBranch <: AbstractNonlinearCurrentDevice end

    struct UndeclaredFormulationCurrentBranch <: AbstractNonlinearCurrentDevice end

    struct UnsupportedFormulationCurrentBranch <: AbstractNonlinearCurrentDevice end

    struct MissingProvenanceCurrentBranch <: AbstractNonlinearCurrentDevice end

    struct NumericalPolicyTaggedCurrentBranch <: AbstractNonlinearCurrentDevice end

    struct UndeclaredCompanionHistoryBranch <: AIMORA.Branches.EMTElement end

    AIMORA.NonlinearNetwork.nonlinear_terminal_nodes(
        ::DuplicateTerminalCurrentBranch,
    ) = (1, 1)

    AIMORA.NonlinearNetwork.nonlinear_device_formulation(
        ::DuplicateTerminalCurrentBranch,
    ) = PhysicalConstitutiveCurrent

    AIMORA.NonlinearNetwork.nonlinear_device_provenance(
        ::DuplicateTerminalCurrentBranch,
    ) = synthetic_physical_parameter_provenance()

    AIMORA.NonlinearNetwork.nonlinear_terminal_nodes(
        ::UndeclaredFormulationCurrentBranch,
    ) = (1, 0)

    AIMORA.NonlinearNetwork.nonlinear_terminal_nodes(
        ::UnsupportedFormulationCurrentBranch,
    ) = (1, 0)

    AIMORA.NonlinearNetwork.nonlinear_device_formulation(
        ::UnsupportedFormulationCurrentBranch,
    ) = SemismoothCurrent

    AIMORA.NonlinearNetwork.nonlinear_device_provenance(
        ::UnsupportedFormulationCurrentBranch,
    ) = synthetic_physical_parameter_provenance()

    AIMORA.NonlinearNetwork.nonlinear_terminal_nodes(
        ::MissingProvenanceCurrentBranch,
    ) = (1, 0)

    AIMORA.NonlinearNetwork.nonlinear_device_formulation(
        ::MissingProvenanceCurrentBranch,
    ) = PhysicalConstitutiveCurrent

    AIMORA.NonlinearNetwork.nonlinear_terminal_nodes(
        ::NumericalPolicyTaggedCurrentBranch,
    ) = (1, 0)

    AIMORA.NonlinearNetwork.nonlinear_device_formulation(
        ::NumericalPolicyTaggedCurrentBranch,
    ) = PhysicalConstitutiveCurrent

    AIMORA.NonlinearNetwork.nonlinear_device_provenance(
        ::NumericalPolicyTaggedCurrentBranch,
    ) = NonlinearParameterProvenance(
        "synthetic numerical policy",
        "dimensionless",
        "direct",
        "not applicable",
        "negative-control constructor",
        NumericalPolicyParameter,
    )

    function AIMORA.NonlinearNetwork.nonlinear_current_jacobian!(
        terminal_current_a::AbstractVector{Float64},
        terminal_jacobian_s::AbstractMatrix{Float64},
        ::DuplicateTerminalCurrentBranch,
        terminal_voltage_v::AbstractVector{Float64},
        time_s::Float64,
    )
        terminal_current_a .= 0.0
        terminal_jacobian_s .= 0.0
        return nothing
    end

    @testset "scaled nonlinear KCL and ideal-constraint solution" begin
        exponential_linear = NodalSystem(
            1,
            [
                ConductanceBranch(1, 0, 0.5),
                CurrentInjection(1, _time_s -> 1.0),
            ],
        )
        exponential = ExponentialCurrentBranch(
            1,
            0;
            saturation_current_a=0.01,
            voltage_scale_v=0.2,
            parallel_conductance_s=0.05,
        )
        exponential_system = NonlinearNodalSystem(
            exponential_linear,
            [exponential];
            scales=NonlinearNetworkScales(
                [1.0],
                [1.0],
                Float64[],
                Float64[],
            ),
        )
        exponential_result = solve_nonlinear_algebraic_state!(
            exponential_system,
            0.0,
            1.0e-5,
        )
        voltage = only(exponential_result.voltage_v)
        @test nonlinear_result_accepted(exponential_result)
        @test abs(0.55 * voltage + 0.01 * expm1(voltage / 0.2) - 1.0) <= 2.0e-12
        @test exponential_result.diagnostics.maximum_kcl_residual_a <= 2.0e-12
        @test exponential_result.diagnostics.numerical_rank == 1
        @test exponential_result.diagnostics.condition_estimate == 1.0
        exponential_device_current_a =
            0.05 * voltage + 0.01 * expm1(voltage / 0.2)
        exponential_power = exponential_result.diagnostics.power
        @test exponential_power.nonlinear_device_absorbed_power_w ≈
            voltage * exponential_device_current_a atol=2.0e-15
        @test exponential_power.ideal_constraint_absorbed_power_w == 0.0
        @test abs(exponential_power.algebraic_power_balance_residual_w) <= 4.0e-12
        @test !isempty(exponential_result.diagnostics.residual_history)
        @test getproperty(
            last(exponential_result.diagnostics.residual_history),
            :maximum_tolerance_normalized_residual,
        ) ==
            exponential_result.diagnostics.maximum_tolerance_normalized_residual

        constrained_linear = NodalSystem(
            2,
            [
                ConductanceBranch(1, 0, 1.0),
                ConductanceBranch(2, 0, 1.0),
                CurrentInjection(1, _time_s -> 2.7),
                CurrentInjection(2, _time_s -> 0.3),
            ],
        )
        constrained_system = NonlinearNodalSystem(
            constrained_linear,
            [CubicCurrentBranch(
                1,
                2;
                linear_conductance_s=0.2,
                cubic_coefficient_a_per_v3=0.1,
            )];
            ideal_constraints=[IdealVoltageConstraint((1, 2), (1.0, -1.0), 1.0)],
            scales=NonlinearNetworkScales([2.0, 1.0], [3.0, 1.0], [1.0], [1.0]),
            options=NonlinearSolveOptions(sparse_dimension_threshold=1),
        )
        constrained_result = solve_nonlinear_algebraic_state!(
            constrained_system,
            0.0,
            1.0e-5,
        )
        @test constrained_result.accepted
        @test constrained_result.voltage_v ≈ [2.0, 1.0] atol=1.0e-12
        @test constrained_result.constraint_current_a ≈ [0.4] atol=1.0e-12
        @test constrained_result.diagnostics.maximum_constraint_residual_v <= 1.0e-12
        constrained_power = constrained_result.diagnostics.power
        @test constrained_power.nonlinear_device_absorbed_power_w ≈ 0.3 atol=2.0e-15
        @test constrained_power.ideal_constraint_absorbed_power_w ≈ 0.4 atol=2.0e-15
        @test abs(constrained_power.algebraic_power_balance_residual_w) <= 4.0e-12
        @test constrained_result.diagnostics.linear_solver == :sparse_umfpack
        @test constrained_result.diagnostics.symbolic_factorization_count == 1
        @test constrained_result.diagnostics.factor_reuse_count >= 1
        @test constrained_result.diagnostics.jacobian_nonzero_count == 8
        fill!(constrained_linear.v, 0.0)
        repeated_constrained_result = solve_nonlinear_algebraic_state!(
            constrained_system,
            0.0,
            1.0e-5,
        )
        @test repeated_constrained_result.accepted
        @test repeated_constrained_result.diagnostics.symbolic_factorization_count == 0
        @test repeated_constrained_result.diagnostics.factor_reuse_count >= 1

        zno_fit = zinc_oxide_piecewise_fit(
            [1.0e-3, 8.0e-3, 64.0e-3],
            [100.0, 200.0, 400.0];
            reference_voltage_v=100.0,
            segment_count=1,
            relative_error_tolerance=1.0e-10,
        )
        zno_device = FittedZincOxideCurrentBranch(1, 0, zno_fit)
        zno_current = zeros(Float64, 2)
        zno_jacobian = zeros(Float64, 2, 2)
        nonlinear_current_jacobian!(
            zno_current,
            zno_jacobian,
            zno_device,
            [200.0, 0.0],
            0.0,
        )
        @test zno_current[1] ≈ 8.0e-3 rtol=1.0e-12
        @test zno_current[1] == -zno_current[2]
        @test zno_jacobian[1, 1] ≈ 3.0 * zno_current[1] / 200.0 rtol=1.0e-12

        accepted_device = TimeLimitedLinearCurrentBranch(1, 0, Inf, 0)
        accepted_system = NonlinearNodalSystem(
            NodalSystem(
                1,
                [
                    ConductanceBranch(1, 0, 1.0),
                    CurrentInjection(1, _time_s -> 1.0),
                ],
            ),
            [accepted_device];
            scales=NonlinearNetworkScales([1.0], [1.0], Float64[], Float64[]),
        )
        accepted_result = solve_nonlinear_algebraic_state!(
            accepted_system,
            0.0,
            1.0e-5,
        )
        @test accepted_result.accepted
        @test accepted_device.accepted_state_count == 1

        failing_device = FailingAcceptanceCurrentBranch(1, 0, 0)
        failing_linear = NodalSystem(
            1,
            [
                ConductanceBranch(1, 0, 1.0),
                CurrentInjection(1, _time_s -> 1.0),
            ],
        )
        failing_voltage_identity = objectid(failing_linear.v)
        failing_system = NonlinearNodalSystem(
            failing_linear,
            [failing_device];
            scales=NonlinearNetworkScales([1.0], [1.0], Float64[], Float64[]),
        )
        failing_result = solve_nonlinear_algebraic_state!(
            failing_system,
            0.0,
            1.0e-5,
        )
        @test !failing_result.accepted
        @test failing_result.failure.code == :state_acceptance_failure
        @test failing_linear.v == [0.0]
        @test objectid(failing_linear.v) == failing_voltage_identity
        @test failing_device.accepted_state_count == 0

        @test_throws ArgumentError NonlinearNodalSystem(
            NodalSystem(1, [ConductanceBranch(1, 0, 1.0)]),
            [DuplicateTerminalCurrentBranch()];
            scales=NonlinearNetworkScales([1.0], [1.0], Float64[], Float64[]),
        )
        @test_throws ArgumentError NonlinearNodalSystem(
            NodalSystem(1, [ConductanceBranch(1, 0, 1.0)]),
            [MissingProvenanceCurrentBranch()];
            scales=NonlinearNetworkScales([1.0], [1.0], Float64[], Float64[]),
        )
        @test_throws ArgumentError NonlinearNodalSystem(
            NodalSystem(1, [ConductanceBranch(1, 0, 1.0)]),
            [NumericalPolicyTaggedCurrentBranch()];
            scales=NonlinearNetworkScales([1.0], [1.0], Float64[], Float64[]),
        )
        @test_throws ArgumentError NonlinearNodalSystem(
            NodalSystem(1, [ConductanceBranch(1, 0, 1.0)]),
            [UndeclaredFormulationCurrentBranch()];
            scales=NonlinearNetworkScales([1.0], [1.0], Float64[], Float64[]),
        )
        @test_throws ArgumentError NonlinearNodalSystem(
            NodalSystem(1, [ConductanceBranch(1, 0, 1.0)]),
            [UnsupportedFormulationCurrentBranch()];
            scales=NonlinearNetworkScales([1.0], [1.0], Float64[], Float64[]),
        )
    end

    @testset "rank diagnosis and rollback-safe discontinuity treatment" begin
        sparse_dimension = 64
        sparse_conductance_s = 10.0 .^
            range(0.0, -8.0; length=sparse_dimension)
        sparse_exact_voltage_v = collect(range(0.5, 1.5; length=sparse_dimension))
        sparse_source_current_a = sparse_conductance_s .* sparse_exact_voltage_v
        sparse_conditioning_system = NonlinearNodalSystem(
            NodalSystem(
                sparse_dimension,
                [DiagonalConditioningNetwork(
                    sparse_conductance_s,
                    sparse_source_current_a,
                )],
            ),
            ();
            scales=NonlinearNetworkScales(
                ones(sparse_dimension),
                ones(sparse_dimension),
                Float64[],
                Float64[],
            ),
        )
        sparse_conditioning_result = solve_nonlinear_algebraic_state!(
            sparse_conditioning_system,
            0.0,
            1.0e-5,
        )
        @test sparse_conditioning_result.accepted
        @test sparse_conditioning_result.voltage_v ≈ sparse_exact_voltage_v atol=2.0e-15
        @test sparse_conditioning_result.diagnostics.condition_estimate ≈ 1.0e8 rtol=2.0e-15
        @test sparse_conditioning_result.diagnostics.numerical_rank == sparse_dimension
        @test sparse_conditioning_result.diagnostics.linear_solver == :sparse_umfpack

        outside_sparse_conductance_s = 10.0 .^
            range(0.0, -13.0; length=sparse_dimension)
        outside_sparse_system = NonlinearNodalSystem(
            NodalSystem(
                sparse_dimension,
                [DiagonalConditioningNetwork(
                    outside_sparse_conductance_s,
                    outside_sparse_conductance_s,
                )],
            ),
            ();
            scales=NonlinearNetworkScales(
                ones(sparse_dimension),
                ones(sparse_dimension),
                Float64[],
                Float64[],
            ),
        )
        outside_sparse_result = solve_nonlinear_algebraic_state!(
            outside_sparse_system,
            0.0,
            1.0e-5,
        )
        @test !outside_sparse_result.accepted
        @test outside_sparse_result.failure.code == :ill_conditioned_network
        @test outside_sparse_result.diagnostics.condition_estimate ≈ 1.0e13 rtol=2.0e-15
        @test outside_sparse_result.diagnostics.numerical_rank == sparse_dimension

        singular_sparse_conductance_s = copy(sparse_conductance_s)
        singular_sparse_conductance_s[end] = 0.0
        singular_sparse_system = NonlinearNodalSystem(
            NodalSystem(
                sparse_dimension,
                [DiagonalConditioningNetwork(
                    singular_sparse_conductance_s,
                    zeros(sparse_dimension),
                )],
            ),
            ();
            scales=NonlinearNetworkScales(
                ones(sparse_dimension),
                ones(sparse_dimension),
                Float64[],
                Float64[],
            ),
        )
        singular_sparse_result = solve_nonlinear_algebraic_state!(
            singular_sparse_system,
            0.0,
            1.0e-5,
        )
        @test !singular_sparse_result.accepted
        @test singular_sparse_result.failure.code == :structural_rank_deficiency
        @test singular_sparse_result.diagnostics.numerical_rank == sparse_dimension - 1

        undeclared_companion_system = NonlinearNodalSystem(
            NodalSystem(1, [UndeclaredCompanionHistoryBranch()]),
            ();
            scales=NonlinearNetworkScales([1.0], [1.0], Float64[], Float64[]),
        )
        @test_throws ArgumentError advance_nonlinear_step!(
            undeclared_companion_system,
            1.0e-3,
            1.0e-3;
            discontinuity_treatment=:two_backward_euler_half_steps,
            discontinuity_reason=:localized_event,
        )

        singular_linear = NodalSystem(1, [ConductanceBranch(1, 0, 1.0)])
        duplicate_constraints = [
            IdealVoltageConstraint((1,), (1.0,), 1.0),
            IdealVoltageConstraint((1,), (1.0,), 1.0),
        ]
        singular_system = NonlinearNodalSystem(
            singular_linear,
            ();
            ideal_constraints=duplicate_constraints,
            scales=NonlinearNetworkScales(
                [1.0],
                [1.0],
                [1.0, 1.0],
                [1.0, 1.0],
            ),
        )
        singular_result = solve_nonlinear_algebraic_state!(singular_system, 0.0, 1.0e-5)
        @test !singular_result.accepted
        @test singular_result.failure.code == :structural_rank_deficiency
        @test singular_linear.v == [0.0]

        branch = SeriesRLBranch(1, 0, 1.0, 1.0e-3)
        rollback_linear = NodalSystem(
            1,
            [branch, CurrentInjection(1, _time_s -> 1.0)],
        )
        time_limited_device = TimeLimitedLinearCurrentBranch(1, 0, 0.75, 0)
        rollback_system = NonlinearNodalSystem(
            rollback_linear,
            [time_limited_device];
            scales=NonlinearNetworkScales(
                [1.0],
                [1.0],
                Float64[],
                Float64[],
            ),
        )
        voltage_identity = objectid(rollback_linear.v)
        accepted_voltage = copy(rollback_linear.v)
        accepted_branch_state = (branch.i_prev, branch.v_prev, branch.i_last)
        rollback_result = advance_nonlinear_step!(
            rollback_system,
            1.0,
            1.0;
            discontinuity_treatment=:two_backward_euler_half_steps,
            discontinuity_reason=:localized_event,
        )
        @test !rollback_result.accepted
        @test rollback_result.failure.code == :nonfinite_device_current
        @test objectid(rollback_linear.v) == voltage_identity
        @test rollback_linear.v == accepted_voltage
        @test rollback_result.voltage_v == accepted_voltage
        @test (branch.i_prev, branch.v_prev, branch.i_last) == accepted_branch_state
        @test time_limited_device.accepted_state_count == 0
        @test unique(getproperty.(
            rollback_result.diagnostics.residual_history,
            :substep_index,
        )) == [1, 2]

        rlc_branch = SeriesRLCBranch(1, 0, 0.2, 2.0e-3, 200.0e-6)
        discontinuity_linear = NodalSystem(
            1,
            [rlc_branch, CurrentInjection(1, _time_s -> 2.0)],
        )
        discontinuity_system = NonlinearNodalSystem(
            discontinuity_linear,
            [CubicCurrentBranch(
                1,
                0;
                linear_conductance_s=0.05,
                cubic_coefficient_a_per_v3=0.002,
            )];
            scales=NonlinearNetworkScales(
                [10.0],
                [2.0],
                Float64[],
                Float64[],
            ),
        )
        discontinuity_result = advance_nonlinear_step!(
            discontinuity_system,
            1.0e-3,
            1.0e-3;
            discontinuity_treatment=:two_backward_euler_half_steps,
            discontinuity_reason=:localized_event,
        )
        @test discontinuity_result.accepted
        @test discontinuity_result.diagnostics.accepted_substep_count == 2
        @test discontinuity_result.diagnostics.companion_method ==
            :two_backward_euler_half_steps
        @test unique(getproperty.(
            discontinuity_result.diagnostics.residual_history,
            :substep_index,
        )) == [1, 2]
        @test all(isfinite, (
            rlc_branch.i_prev,
            rlc_branch.inductor_voltage_prev,
            rlc_branch.capacitor_voltage_prev,
        ))
        checkpoint = nonlinear_nodal_checkpoint(discontinuity_system)
        checkpoint_voltage = copy(discontinuity_linear.v)
        checkpoint_branch_state = (
            rlc_branch.i_prev,
            rlc_branch.inductor_voltage_prev,
            rlc_branch.capacitor_voltage_prev,
        )
        unsupported_schema_checkpoint = NonlinearNodalCheckpoint(
            checkpoint.schema_version + 1,
            checkpoint.accepted_state,
            checkpoint.topology_signature,
            checkpoint.solve_count,
        )
        @test_throws ArgumentError restore_nonlinear_nodal_checkpoint!(
            discontinuity_system,
            unsupported_schema_checkpoint,
        )
        @test discontinuity_linear.v == checkpoint_voltage
        uninterrupted_result = advance_nonlinear_step!(
            discontinuity_system,
            2.0e-3,
            1.0e-3,
        )
        @test uninterrupted_result.accepted
        uninterrupted_voltage = copy(uninterrupted_result.voltage_v)
        restore_nonlinear_nodal_checkpoint!(discontinuity_system, checkpoint)
        @test discontinuity_linear.v == checkpoint_voltage
        @test (
            rlc_branch.i_prev,
            rlc_branch.inductor_voltage_prev,
            rlc_branch.capacitor_voltage_prev,
        ) == checkpoint_branch_state
        restarted_result = advance_nonlinear_step!(
            discontinuity_system,
            2.0e-3,
            1.0e-3,
        )
        @test restarted_result.accepted
        @test reinterpret(UInt64, restarted_result.voltage_v) ==
            reinterpret(UInt64, uninterrupted_voltage)
        mismatched_voltage = copy(rollback_linear.v)
        topology_failure = try
            restore_nonlinear_nodal_checkpoint!(rollback_system, checkpoint)
            nothing
        catch error
            error
        end
        @test topology_failure isa NonlinearSolveFailure
        @test topology_failure.code == :invalid_topology_signature
        @test rollback_linear.v == mismatched_voltage
        @test_throws ArgumentError advance_nonlinear_step!(
            discontinuity_system,
            2.0e-3,
            1.0e-3;
            discontinuity_treatment=:two_backward_euler_half_steps,
        )
    end

    @testset "typed nonlinear EMT schedule and discontinuity order" begin
        accepted_breaker_event = EMTHybridEventOccurrence(
            :breaker_open,
            3.0e-3,
            0.0,
            -10,
            true,
            2,
            1.0e-12,
        )
        schedule = NonlinearEMTStudySchedule(
            1.0e-3,
            3.0e-3;
            discontinuities=[
                NonlinearEMTDiscontinuity(1.0e-3, :localized_event),
                NonlinearEMTDiscontinuity(accepted_breaker_event),
                NonlinearEMTDiscontinuity(
                    3.0e-3,
                    :localized_event;
                    event_name=:controller_limit,
                    priority=-5,
                ),
            ],
            accepted_tasks=[EMTSampledTaskOccurrence(
                :sampled_control,
                2.0e-3,
                2,
                -1,
                1,
            )],
        )
        scheduled_branch = SeriesRLCBranch(1, 0, 0.2, 2.0e-3, 200.0e-6)
        scheduled_linear = NodalSystem(
            1,
            [scheduled_branch, CurrentInjection(1, _time_s -> 2.0)],
        )
        scheduled_system = NonlinearNodalSystem(
            scheduled_linear,
            [CubicCurrentBranch(
                1,
                0;
                linear_conductance_s=0.05,
                cubic_coefficient_a_per_v3=0.002,
            )];
            scales=NonlinearNetworkScales(
                [10.0],
                [2.0],
                Float64[],
                Float64[],
            ),
        )
        trace = evaluate_nonlinear_emt_network!(scheduled_system, schedule)
        @test trace.time_s == [0.0, 1.0e-3, 2.0e-3, 3.0e-3]
        @test size(trace.voltage_v) == (1, 4)
        @test size(trace.constraint_current_a) == (0, 4)
        @test length(trace.diagnostics) == 3
        @test getproperty.(trace.diagnostics, :discontinuity_reason) ==
            [:localized_event, :none, :topology_change]
        @test getproperty.(trace.diagnostics, :companion_method) == [
            :two_backward_euler_half_steps,
            :trapezoidal,
            :two_backward_euler_half_steps,
        ]
        @test getproperty.(trace.diagnostics, :accepted_substep_count) == [2, 1, 2]
        @test all(isfinite, trace.voltage_v)
        @test getproperty.(trace.accepted_discontinuities, :event_name) ==
            [:localized_event, :breaker_open, :controller_limit]
        @test getproperty.(trace.accepted_discontinuities, :priority) == [0, -10, -5]
        @test getproperty.(trace.accepted_tasks, :name) == [:sampled_control]

        @test_throws ArgumentError NonlinearEMTStudySchedule(
            1.0e-3,
            3.0e-3;
            discontinuities=[
                NonlinearEMTDiscontinuity(
                    1.0e-3,
                    :localized_event;
                    event_name=:duplicate,
                    priority=-2,
                ),
                NonlinearEMTDiscontinuity(
                    1.0e-3,
                    :topology_change;
                    event_name=:duplicate,
                    priority=-1,
                ),
            ],
        )
        @test_throws ArgumentError NonlinearEMTStudySchedule(
            1.0e-3,
            3.0e-3;
            accepted_tasks=[
                EMTSampledTaskOccurrence(:late_task, 1.0e-3, 1, 1, 1),
                EMTSampledTaskOccurrence(:early_task, 1.0e-3, 1, -1, 1),
            ],
        )
        @test_throws ArgumentError NonlinearEMTStudySchedule(
            1.0e-3,
            3.0e-3;
            discontinuities=[
                NonlinearEMTDiscontinuity(
                    1.0e-3,
                    :localized_event;
                    event_name=:late_priority,
                    priority=1,
                ),
                NonlinearEMTDiscontinuity(
                    1.0e-3,
                    :topology_change;
                    event_name=:early_priority,
                    priority=-1,
                ),
            ],
        )
        @test_throws ArgumentError NonlinearEMTStudySchedule(
            1.0e-3,
            3.0e-3;
            discontinuities=[
                NonlinearEMTDiscontinuity(1.5e-3, :localized_event),
            ],
        )
        @test_throws ArgumentError NonlinearEMTStudySchedule(1.0e-3, 3.5e-3)

        task_boundary_linear = NodalSystem(
            1,
            [
                ConductanceBranch(1, 0, 1.0),
                CurrentInjection(1, _time_s -> 1.0),
            ],
        )
        task_boundary_system = NonlinearNodalSystem(
            task_boundary_linear,
            ();
            scales=NonlinearNetworkScales([1.0], [1.0], Float64[], Float64[]),
        )
        task_boundary_schedule = NonlinearEMTStudySchedule(
            1.0e-3,
            1.0e-3;
            accepted_tasks=[EMTSampledTaskOccurrence(
                :control_tick,
                1.0e-3,
                1,
                0,
                1,
            )],
        )
        critical_damping_decision = classify_numerical_chatter(
            NonlinearChatterObservation(
                (0.10, -0.09, 0.095),
                1.0;
                topology_unchanged=true,
                task_calendar_unchanged=true,
                control_mode_unchanged=true,
            ),
        )
        @test_throws ArgumentError evaluate_nonlinear_emt_network!(
            task_boundary_system,
            task_boundary_schedule;
            chatter_decisions=Dict(1 => critical_damping_decision),
        )
        @test task_boundary_linear.v == [0.0]
    end

    @testset "sparse topology invalidation and numeric reuse" begin
        switching_time_s = 2.0e-3
        switching_conductance_s = 0.25
        open_voltage_v = [1.0, 2.0]
        switched_voltage_v = [1.2, 1.8]
        retained_voltage_v = [1.25, 1.75]
        switch = TimeSwitch(
            1,
            2;
            close_time_s=switching_time_s,
            on_conductance=switching_conductance_s,
            off_conductance=0.0,
        )
        source_one = CurrentInjection(
            1,
            time_s -> begin
                voltage = time_s < switching_time_s ? open_voltage_v :
                    time_s < 3.0e-3 ? switched_voltage_v : retained_voltage_v
                voltage[1] + (time_s >= switching_time_s ?
                    switching_conductance_s * (voltage[1] - voltage[2]) : 0.0)
            end,
        )
        source_two = CurrentInjection(
            2,
            time_s -> begin
                voltage = time_s < switching_time_s ? open_voltage_v :
                    time_s < 3.0e-3 ? switched_voltage_v : retained_voltage_v
                node_voltage_v = voltage[2]
                node_voltage_v + 0.1 * node_voltage_v + 0.01 * node_voltage_v^3 +
                    (time_s >= switching_time_s ?
                        switching_conductance_s * (voltage[2] - voltage[1]) : 0.0)
            end,
        )
        topology_linear = NodalSystem(
            2,
            [
                ConductanceBranch(1, 0, 1.0),
                ConductanceBranch(2, 0, 1.0),
                source_one,
                source_two,
                switch,
            ],
        )
        topology_system = NonlinearNodalSystem(
            topology_linear,
            [CubicCurrentBranch(
                2,
                0;
                linear_conductance_s=0.1,
                cubic_coefficient_a_per_v3=0.01,
            )];
            scales=NonlinearNetworkScales(
                retained_voltage_v,
                [1.0, 2.5],
                Float64[],
                Float64[],
            ),
            options=NonlinearSolveOptions(sparse_dimension_threshold=1),
        )
        open_result = advance_nonlinear_step!(topology_system, 1.0e-3, 1.0e-3)
        @test open_result.accepted
        @test open_result.voltage_v ≈ open_voltage_v atol=2.0e-12
        open_signature = open_result.diagnostics.topology_signature

        rejected_topology_linear = NodalSystem(
            2,
            [
                ConductanceBranch(1, 0, 1.0),
                ConductanceBranch(2, 0, 1.0),
                source_one,
                source_two,
                TimeSwitch(
                    1,
                    2;
                    close_time_s=switching_time_s,
                    on_conductance=switching_conductance_s,
                    off_conductance=0.0,
                ),
            ],
        )
        rejected_topology_device = TimeLimitedLinearCurrentBranch(
            2,
            0,
            switching_time_s,
            0,
        )
        rejected_topology_system = NonlinearNodalSystem(
            rejected_topology_linear,
            [rejected_topology_device];
            scales=NonlinearNetworkScales(
                retained_voltage_v,
                [1.0, 2.5],
                Float64[],
                Float64[],
            ),
            options=NonlinearSolveOptions(sparse_dimension_threshold=1),
        )
        accepted_before_topology_failure = advance_nonlinear_step!(
            rejected_topology_system,
            1.0e-3,
            1.0e-3,
        )
        @test accepted_before_topology_failure.accepted
        accepted_topology_voltage = copy(rejected_topology_linear.v)
        accepted_topology_signature = rejected_topology_system.topology_signature
        accepted_previous_topology_signature =
            rejected_topology_system.previous_topology_signature
        rejected_topology_result = advance_nonlinear_step!(
            rejected_topology_system,
            switching_time_s,
            1.0e-3;
            discontinuity_treatment=:two_backward_euler_half_steps,
            discontinuity_reason=:topology_change,
        )
        @test !rejected_topology_result.accepted
        @test rejected_topology_result.failure.code == :nonfinite_device_current
        @test rejected_topology_linear.v == accepted_topology_voltage
        @test rejected_topology_system.topology_signature ==
            accepted_topology_signature
        @test rejected_topology_system.previous_topology_signature ==
            accepted_previous_topology_signature
        @test rejected_topology_system.workspace.sparse_factor === nothing
        @test rejected_topology_system.workspace.sparse_pattern_signature == UInt64(0)

        topology_result = advance_nonlinear_step!(
            topology_system,
            switching_time_s,
            1.0e-3;
            discontinuity_treatment=:two_backward_euler_half_steps,
            discontinuity_reason=:topology_change,
        )
        @test topology_result.accepted
        @test topology_result.voltage_v ≈ switched_voltage_v atol=2.0e-12
        @test topology_result.diagnostics.topology_signature != open_signature
        @test topology_result.diagnostics.symbolic_factorization_count == 1
        retained_result = advance_nonlinear_step!(
            topology_system,
            3.0e-3,
            1.0e-3;
            discontinuity_treatment=:two_backward_euler_half_steps,
            discontinuity_reason=:localized_event,
        )
        @test retained_result.accepted
        @test retained_result.voltage_v ≈ retained_voltage_v atol=2.0e-12
        @test retained_result.diagnostics.symbolic_factorization_count == 0
        @test retained_result.diagnostics.factor_reuse_count >= 1
    end
end

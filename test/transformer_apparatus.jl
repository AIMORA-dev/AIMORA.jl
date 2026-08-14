using Test
using AIMORA
using LinearAlgebra
using AIMORA.NonlinearNetwork

const TransformerPublic = AIMORA.TransformerApparatus

function synthetic_transformer_source()
    provenance = AIMORA.StudyCore.ParameterProvenance(
        "AIMORA synthetic transformer apparatus fixture",
        "SI",
        "direct deterministic fixture values",
        "exact synthetic values without measurement uncertainty",
        "public transformer apparatus contract tests",
        AIMORA.StudyCore.PhysicalModelParameter,
    )
    return TransformerPublic.TransformerSourceRecord(
        :synthetic_transformer_apparatus,
        repeat("a", 64),
        provenance,
    )
end

function test_model_consistent_transformer_initialization()
@testset "model-consistent transformer initialization without hidden settling" begin
    source = synthetic_transformer_source()
    connection = grounded_transformer_connection()
    graph = linear_transformer_magnetic_graph(source)
    settings = TransformerPublic.TransformerRuntimeSettings(
        timestep_s=10.0e-6,
        initialization_frequency_hz=60.0,
    )
    physical_frequency = 2.0 * pi * settings.initialization_frequency_hz
    discrete_frequency = (2.0 / settings.timestep_s) *
        tan(0.5 * physical_frequency * settings.timestep_s)
    initial_time = 0.0
    next_time = settings.timestep_s
    terminal_voltage_phasor = ComplexF64[10.0 + 2.0im, -4.0 + 1.0im]
    sample(phasor, time_s) = real.(phasor .* cis(physical_frequency * time_s))
    derivative(phasor, time_s) =
        real.(im * physical_frequency .* phasor .* cis(physical_frequency * time_s))

    terminal_matrices = TransformerPublic.TransformerTerminalMatrices(
        [0.4 0.0; 0.0 0.5],
        [0.12 0.02; 0.02 0.10];
        capacitance_f=[2.0e-9 0.0; 0.0 1.0e-9],
        conductance_s=[1.0e-7 0.0; 0.0 2.0e-7],
    )
    bctran = TransformerPublic.BCTRANTransformerModel(
        terminal_matrices;
        positive_pair_reconstruction_residual=0.0,
        zero_pair_reconstruction_residual=0.0,
        inverse_reconstruction_residual=0.0,
    )
    magnetic_inductance =
        TransformerPublic.transformer_magnetic_linear_inductance(graph)
    branch_flux_per_current = hcat([
        TransformerPublic.transformer_magnetic_linear_response(
            graph,
            [coil == basis ? 1.0 : 0.0 for coil in 1:2],
        ).branch_flux_wb for basis in 1:2
    ]...)
    winding_loss_state_matrix = -1_000.0 .* Matrix{Float64}(I, 2, 2)
    winding_loss_input_matrix = Matrix{Float64}(I, 2, 2)
    winding_loss_output_matrix = Matrix{Float64}(I, 2, 2)
    winding_loss_direct = 0.05 .* Matrix{Float64}(I, 2, 2)
    hybrid = TransformerPublic.HybridTransformerModel(
        terminal_matrices,
        graph;
        winding_loss_state_matrix_per_s=winding_loss_state_matrix,
        winding_loss_input_matrix=winding_loss_input_matrix,
        winding_loss_output_matrix_ohm_per_s=winding_loss_output_matrix,
        winding_loss_direct_ohm=winding_loss_direct,
        winding_loss_storage_matrix_j=Matrix{Float64}(I, 2, 2),
        continuous_passivity_margin_ohm=0.05,
    )
    mec = TransformerPublic.MagneticEquivalentCircuitModel(
        terminal_matrices.resistance_ohm,
        terminal_matrices.inductance_h,
        graph;
        terminal_capacitance_f=terminal_matrices.capacitance_f,
        terminal_conductance_s=terminal_matrices.conductance_s,
    )
    winding_loss_impedance = winding_loss_output_matrix * (
        (
            im * discrete_frequency .* Matrix{ComplexF64}(I, 2, 2) .-
            ComplexF64.(winding_loss_state_matrix)
        ) \ ComplexF64.(winding_loss_input_matrix)
    ) .+ ComplexF64.(winding_loss_direct)

    terminal_models = (
        (
            tier=TransformerPublic.LowFrequencyTerminalTier,
            model=TransformerPublic.LowFrequencyTransformerModel(terminal_matrices),
            impedance=ComplexF64.(terminal_matrices.resistance_ohm) .+
                im * discrete_frequency .* ComplexF64.(terminal_matrices.inductance_h),
        ),
        (
            tier=TransformerPublic.BCTRANTerminalTier,
            model=bctran,
            impedance=ComplexF64.(terminal_matrices.resistance_ohm) .+
                im * discrete_frequency .* ComplexF64.(terminal_matrices.inductance_h),
        ),
    )
    for row in terminal_models
        current_phasor = row.impedance \ terminal_voltage_phasor
        specification = transformer_specification(
            row.tier,
            connection,
            row.model,
            source,
        )
        preparations = [
            TransformerPublic.prepare_transformer_apparatus(
                specification;
                initialization_mode=
                    TransformerPublic.SinusoidalTransformerOperatingPoint,
                initial_time_s=time_s,
                initial_node_voltage_v=sample(terminal_voltage_phasor, time_s),
                initial_node_voltage_derivative_v_per_s=
                    derivative(terminal_voltage_phasor, time_s),
                initial_coil_current_a=sample(current_phasor, time_s),
                initial_coil_current_derivative_a_per_s=
                    derivative(current_phasor, time_s),
            ) for time_s in (initial_time, next_time)
        ]
        @test preparations[1].operating_point_electrical_residual_v <= 1.0e-7
        runtime = TransformerPublic.transformer_apparatus_runtime(
            preparations[1],
            [1, 2],
        )
        expected = TransformerPublic.transformer_apparatus_runtime(
            preparations[2],
            [1, 2],
        )
        accept_transformer_test_step!(
            runtime,
            sample(terminal_voltage_phasor, next_time),
            next_time,
        )
        @test runtime.accepted_state.coil_current_a ≈
            expected.accepted_state.coil_current_a rtol=1.0e-10 atol=1.0e-11
        @test runtime.accepted_state.capacitor_current_a ≈
            expected.accepted_state.capacitor_current_a rtol=1.0e-10 atol=1.0e-11
        @test runtime.accepted_state.terminal_current_a ≈
            expected.accepted_state.terminal_current_a rtol=1.0e-10 atol=1.0e-11
    end

    magnetic_models = (
        (
            tier=TransformerPublic.HybridTransformerTier,
            model=hybrid,
            impedance=ComplexF64.(terminal_matrices.resistance_ohm) .+
                winding_loss_impedance .+
                im * discrete_frequency .* ComplexF64.(
                    terminal_matrices.inductance_h .+ magnetic_inductance,
                ),
        ),
        (
            tier=TransformerPublic.MagneticEquivalentCircuitTier,
            model=mec,
            impedance=ComplexF64.(terminal_matrices.resistance_ohm) .+
                im * discrete_frequency .* ComplexF64.(
                    terminal_matrices.inductance_h .+ magnetic_inductance,
                ),
        ),
    )
    for row in magnetic_models
        current_phasor = row.impedance \ terminal_voltage_phasor
        branch_flux_phasor = branch_flux_per_current * current_phasor
        specification = transformer_specification(
            row.tier,
            connection,
            row.model,
            source,
        )
        preparations = [
            TransformerPublic.prepare_transformer_apparatus(
                specification;
                initialization_mode=
                    TransformerPublic.SinusoidalTransformerOperatingPoint,
                initial_time_s=time_s,
                initial_node_voltage_v=sample(terminal_voltage_phasor, time_s),
                initial_node_voltage_derivative_v_per_s=
                    derivative(terminal_voltage_phasor, time_s),
                initial_coil_current_a=sample(current_phasor, time_s),
                initial_coil_current_derivative_a_per_s=
                    derivative(current_phasor, time_s),
                initial_branch_flux_wb=sample(branch_flux_phasor, time_s),
                initial_branch_flux_derivative_wb_per_s=
                    derivative(branch_flux_phasor, time_s),
            ) for time_s in (initial_time, next_time)
        ]
        @test preparations[1].operating_point_electrical_residual_v <= 1.0e-7
        runtime = TransformerPublic.transformer_apparatus_runtime(
            preparations[1],
            [1, 2],
        )
        expected = TransformerPublic.transformer_apparatus_runtime(
            preparations[2],
            [1, 2],
        )
        accept_transformer_test_step!(
            runtime,
            sample(terminal_voltage_phasor, next_time),
            next_time,
        )
        @test runtime.accepted_state.coil_current_a ≈
            expected.accepted_state.coil_current_a rtol=2.0e-9 atol=1.0e-10
        @test runtime.accepted_state.branch_flux_wb ≈
            expected.accepted_state.branch_flux_wb rtol=2.0e-9 atol=1.0e-12
        @test runtime.accepted_state.winding_loss_state ≈
            expected.accepted_state.winding_loss_state rtol=2.0e-9 atol=1.0e-11
        @test runtime.accepted_state.terminal_current_a ≈
            expected.accepted_state.terminal_current_a rtol=2.0e-9 atol=1.0e-10
    end

    wideband = TransformerPublic.WidebandTransformerModel(
        -1_000.0 .* Matrix{Float64}(I, 2, 2),
        Matrix{Float64}(I, 2, 2),
        Matrix{Float64}(I, 2, 2),
        0.1 .* Matrix{Float64}(I, 2, 2);
        storage_matrix_j=Matrix{Float64}(I, 2, 2),
        port_order=connection.node_order,
        frequency_band_hz=(10.0, 10_000.0),
        continuous_passivity_margin_s=0.0,
        enforcement_perturbation_relative=0.0,
        source_response_sha256=repeat("e", 64),
    )
    grey_box = TransformerPublic.GreyBoxTransformerModel(
        node_order=[:primary_terminal, :secondary_terminal, :internal_node],
        terminal_node_indices=[1, 2],
        branches=[
            TransformerPublic.TransformerLadderBranch(
                :primary_section,
                1,
                3;
                resistance_ohm=0.2,
                inductance_h=2.0e-3,
            ),
            TransformerPublic.TransformerLadderBranch(
                :secondary_section,
                3,
                2;
                resistance_ohm=0.3,
                inductance_h=3.0e-3,
            ),
        ],
        capacitance_f=Diagonal([1.0e-9, 1.0e-9, 2.0e-9]),
        conductance_s=Diagonal([1.0e-7, 1.0e-7, 2.0e-7]),
        source_response_sha256=repeat("f", 64),
        identification_residual_relative=1.0e-4,
        parameter_nonuniqueness="synthetic initialization fixture",
    )
    white_box = TransformerPublic.WhiteBoxTransformerModel(
        winding_order=[:primary_winding, :secondary_winding],
        section_count_per_winding=[2, 3],
        section_length_m=[0.5, 0.5, 0.5],
        series_resistance_ohm_per_m=fill([0.2 0.02; 0.02 0.3], 3),
        series_inductance_h_per_m=fill([2.0e-3 0.2e-3; 0.2e-3 3.0e-3], 3),
        shunt_conductance_s_per_m=
            fill(1.0e-7 .* Matrix{Float64}(I, 2, 2), 3),
        shunt_capacitance_f_per_m=
            fill(1.0e-9 .* Matrix{Float64}(I, 2, 2), 3),
        geometry_sha256=repeat("1", 64),
        frequency_band_hz=(10.0, 10_000.0),
        section_refinement_residual_relative=1.0e-3,
    )
    represented_rows = (
        (
            tier=TransformerPublic.WidebandBlackBoxTier,
            model=wideband,
            connection=connection,
            voltage_phasor=terminal_voltage_phasor,
            nodes=[1, 2],
        ),
        (
            tier=TransformerPublic.GreyBoxLadderTier,
            model=grey_box,
            connection=connection,
            voltage_phasor=terminal_voltage_phasor,
            nodes=[1, 2],
        ),
        (
            tier=TransformerPublic.WhiteBoxWindingTier,
            model=white_box,
            connection=white_box_transformer_connection(),
            voltage_phasor=ComplexF64[
                10.0 + 2.0im,
                0.0 + 0.0im,
                -4.0 + 1.0im,
                0.0 + 0.0im,
            ],
            nodes=[1, 2, 3, 4],
        ),
    )
    for row in represented_rows
        specification = transformer_specification(
            row.tier,
            row.connection,
            row.model,
            source,
        )
        preparations = [
            TransformerPublic.prepare_transformer_apparatus(
                specification;
                initialization_mode=
                    TransformerPublic.SinusoidalTransformerOperatingPoint,
                initial_time_s=time_s,
                initial_node_voltage_v=sample(row.voltage_phasor, time_s),
                initial_node_voltage_derivative_v_per_s=
                    derivative(row.voltage_phasor, time_s),
            ) for time_s in (initial_time, next_time)
        ]
        runtime = TransformerPublic.transformer_apparatus_runtime(
            preparations[1],
            row.nodes,
        )
        expected = TransformerPublic.transformer_apparatus_runtime(
            preparations[2],
            row.nodes,
        )
        accept_transformer_test_step!(
            runtime,
            sample(row.voltage_phasor, next_time),
            next_time,
        )
        state = runtime.accepted_state
        expected_state = expected.accepted_state
        if hasproperty(state, :rational_state)
            @test state.rational_state ≈
                expected_state.rational_state rtol=1.0e-10 atol=1.0e-11
        else
            @test state.represented_node_voltage_v ≈
                expected_state.represented_node_voltage_v rtol=2.0e-9 atol=1.0e-10
            @test state.branch_current_a ≈
                expected_state.branch_current_a rtol=2.0e-9 atol=1.0e-10
            @test state.capacitor_current_a ≈
                expected_state.capacitor_current_a rtol=2.0e-9 atol=1.0e-11
            @test state.maximum_internal_kcl_residual_a <= 1.0e-9
        end
        @test state.terminal_current_a ≈
            expected_state.terminal_current_a rtol=2.0e-9 atol=1.0e-10
        @test state.accepted_time_s == next_time
    end

    hysteretic_model = TransformerPublic.HybridTransformerModel(
        TransformerPublic.TransformerTerminalMatrices(
            [0.2 0.0; 0.0 0.2],
            [1.0e-3 0.0; 0.0 1.0e-3],
        ),
        hysteretic_transformer_magnetic_graph(source),
    )
    hysteretic_specification = transformer_specification(
        TransformerPublic.HybridTransformerTier,
        connection,
        hysteretic_model,
        source,
    )
    residual_flux = [1.0e-3, 1.0e-3]
    residual_preparation = TransformerPublic.prepare_transformer_apparatus(
        hysteretic_specification;
        initialization_mode=TransformerPublic.DeenergizedTransformerInitialization,
        initial_branch_flux_wb=residual_flux,
    )
    residual_runtime = TransformerPublic.transformer_apparatus_runtime(
        residual_preparation,
        [1, 2],
    )
    @test residual_runtime.accepted_state.branch_flux_wb == residual_flux
    @test all(
        state -> state.field_strength_a_per_m == 0.0,
        residual_runtime.accepted_state.tellinen_state,
    )
    @test residual_runtime.accepted_state.maximum_magnetic_continuity_residual_wb <=
        settings.magnetic_continuity_absolute_tolerance_wb
    @test_throws TransformerPublic.TransformerApparatusRefusal begin
        TransformerPublic.transformer_apparatus_runtime(
            TransformerPublic.prepare_transformer_apparatus(
                transformer_specification(
                    TransformerPublic.HybridTransformerTier,
                    connection,
                    TransformerPublic.HybridTransformerModel(
                        terminal_matrices,
                        nonlinear_transformer_magnetic_graph(source),
                    ),
                    source,
                );
                initialization_mode=
                    TransformerPublic.DeenergizedTransformerInitialization,
                initial_branch_flux_wb=residual_flux,
            ),
            [1, 2],
        )
    end
    @test_throws TransformerPublic.TransformerApparatusRefusal begin
        TransformerPublic.prepare_transformer_apparatus(
            hysteretic_specification;
            initialization_mode=TransformerPublic.DeenergizedTransformerInitialization,
            initial_branch_flux_wb=residual_flux,
            initial_branch_flux_derivative_wb_per_s=[1.0, 1.0],
        )
    end
end
end

function grounded_transformer_connection()
    return TransformerPublic.TransformerConnectionTopology(
        node_order=[:primary_terminal, :secondary_terminal],
        coil_order=[:primary_coil, :secondary_coil],
        winding_order=[:primary_winding, :secondary_winding],
        phase_order=[:phase_a],
        coil_winding=[:primary_winding, :secondary_winding],
        coil_phase=[:phase_a, :phase_a],
        incidence=Matrix{Float64}(I, 2, 2),
        vector_group="Ii0",
        clock_number=0,
        phase_shift_rad=0.0,
    )
end

function single_coil_reactor_connection()
    return TransformerPublic.TransformerConnectionTopology(
        node_order=[:reactor_line_terminal, :reactor_neutral_terminal],
        coil_order=[:reactor_coil],
        winding_order=[:reactor_winding],
        phase_order=[:phase_a],
        coil_winding=[:reactor_winding],
        coil_phase=[:phase_a],
        incidence=reshape([1.0, -1.0], 2, 1),
        vector_group="I0",
        clock_number=0,
        phase_shift_rad=0.0,
    )
end

function split_winding_reactor_connection()
    return TransformerPublic.TransformerConnectionTopology(
        node_order=[
            :reactor_line_terminal,
            :reactor_section_junction,
            :reactor_neutral_terminal,
        ],
        coil_order=[:reactor_line_section, :reactor_neutral_section],
        winding_order=[:reactor_line_winding, :reactor_neutral_winding],
        phase_order=[:phase_a],
        coil_winding=[:reactor_line_winding, :reactor_neutral_winding],
        coil_phase=[:phase_a, :phase_a],
        incidence=[
            1.0 0.0
            -1.0 1.0
            0.0 -1.0
        ],
        vector_group="Ii0",
        clock_number=0,
        phase_shift_rad=0.0,
    )
end

function white_box_transformer_connection()
    return TransformerPublic.TransformerConnectionTopology(
        node_order=[
            :primary_start,
            :primary_end,
            :secondary_start,
            :secondary_end,
        ],
        coil_order=[:primary_coil, :secondary_coil],
        winding_order=[:primary_winding, :secondary_winding],
        phase_order=[:phase_a],
        coil_winding=[:primary_winding, :secondary_winding],
        coil_phase=[:phase_a, :phase_a],
        incidence=[
            1.0 0.0
            -1.0 0.0
            0.0 1.0
            0.0 -1.0
        ],
        vector_group="Ii0",
        clock_number=0,
        phase_shift_rad=0.0,
    )
end

function linear_transformer_magnetic_graph(source)
    material = TransformerPublic.LinearTransformerMagneticMaterial(
        2_000.0,
        2.0,
        source,
    )
    return TransformerPublic.TransformerMagneticGraph(
        node_order=[:magnetic_node],
        branches=[
            TransformerPublic.MagneticBranchGeometry(
                :core_limb;
                length_m=0.8,
                cross_section_m2=0.01,
            ),
            TransformerPublic.MagneticBranchGeometry(
                :return_limb;
                length_m=0.8,
                cross_section_m2=0.01,
            ),
        ],
        incidence=reshape([1.0, -1.0], 1, 2),
        winding_turns=[100.0 0.0; 0.0 100.0],
        materials=[material],
    )
end

function nonlinear_transformer_magnetic_graph(source)
    material = TransformerPublic.PiecewiseLinearTransformerMagneticMaterial(
        [0.0, 0.5, 1.5, 2.0],
        [0.0, 100.0, 2_000.0, 10_000.0],
        source,
    )
    return TransformerPublic.TransformerMagneticGraph(
        node_order=[:magnetic_node],
        branches=[
            TransformerPublic.MagneticBranchGeometry(
                :nonlinear_core_limb;
                length_m=0.8,
                cross_section_m2=0.01,
            ),
            TransformerPublic.MagneticBranchGeometry(
                :nonlinear_return_limb;
                length_m=0.8,
                cross_section_m2=0.01,
            ),
        ],
        incidence=reshape([1.0, -1.0], 1, 2),
        winding_turns=[100.0 0.0; 0.0 100.0],
        materials=[material],
    )
end

function hysteretic_transformer_magnetic_graph(source)
    lower_curve = TransformerPublic.TellinenLimitingCurve(
        [-1_000.0, -500.0, 0.0, 500.0, 1_000.0],
        [-1.5, -1.0, -0.2, 0.5, 1.0],
    )
    upper_curve = TransformerPublic.TellinenLimitingCurve(
        [-1_000.0, -500.0, 0.0, 500.0, 1_000.0],
        [-1.0, -0.5, 0.2, 1.0, 1.5],
    )
    material = TransformerPublic.TellinenTransformerMagneticMaterial(
        lower_curve,
        upper_curve,
        source;
        integration_field_increment_a_per_m=2.0,
    )
    return TransformerPublic.TransformerMagneticGraph(
        node_order=[:magnetic_node],
        branches=[
            TransformerPublic.MagneticBranchGeometry(
                :hysteretic_core_limb;
                length_m=0.8,
                cross_section_m2=0.01,
            ),
            TransformerPublic.MagneticBranchGeometry(
                :hysteretic_return_limb;
                length_m=0.8,
                cross_section_m2=0.01,
            ),
        ],
        incidence=reshape([1.0, -1.0], 1, 2),
        winding_turns=[100.0 0.0; 0.0 100.0],
        materials=[material],
    )
end

function transformer_specification(
    tier,
    connection,
    model,
    source;
    reactor_definition=nothing,
)
    return TransformerPublic.TransformerApparatusSpecification(
        Symbol(TransformerPublic.transformer_apparatus_contract(tier).id),
        tier,
        connection,
        model,
        TransformerPublic.TransformerRuntimeSettings(
            timestep_s=10.0e-6,
            initialization_frequency_hz=60.0,
        );
        phase_count=1,
        rated_power_va=1.0e6,
        rated_voltage_v=13.8e3,
        rated_frequency_hz=60.0,
        reactor_definition,
        sources=[source],
        uncertainty="exact synthetic contract fixture",
        validity_domain="one-phase public apparatus contract test",
    )
end

function accept_transformer_test_step!(
    runtime,
    voltage,
    time_s=10.0e-6;
    step_s=10.0e-6,
    companion_method=:trapezoidal,
)
    terminal_count = length(voltage)
    current = zeros(Float64, terminal_count)
    jacobian = zeros(Float64, terminal_count, terminal_count)
    prepare_nonlinear_device_step!(
        runtime,
        time_s,
        step_s,
        companion_method,
    )
    nonlinear_current_jacobian!(
        current,
        jacobian,
        runtime,
        Float64.(voltage),
        time_s,
    )
    accept_nonlinear_device_state!(
        runtime,
        Float64.(voltage),
        current,
        jacobian,
        time_s,
    )
    finish_nonlinear_device_step!(runtime)
    return current, jacobian
end


@testset "explicit air-core and iron-core reactor families" begin
    source = synthetic_transformer_source()
    connection = single_coil_reactor_connection()
    terminal_matrices = TransformerPublic.TransformerTerminalMatrices(
        reshape([0.25], 1, 1),
        reshape([8.0e-3], 1, 1);
        capacitance_f=reshape([2.0e-8], 1, 1),
        conductance_s=reshape([1.0e-6], 1, 1),
    )
    air_core_definition = TransformerPublic.ReactorApparatusDefinition(
        TransformerPublic.ShuntReactorApplication,
        TransformerPublic.AirCoreReactorConstruction;
        control_mode=TransformerPublic.BreakerSwitchedReactorControl,
    )
    air_core_specification = transformer_specification(
        TransformerPublic.LowFrequencyTerminalTier,
        connection,
        TransformerPublic.LowFrequencyTransformerModel(terminal_matrices),
        source;
        reactor_definition=air_core_definition,
    )
    air_core_runtime = TransformerPublic.transformer_apparatus_runtime(
        TransformerPublic.prepare_transformer_apparatus(air_core_specification),
        [1, 2],
    )
    TransformerPublic.queue_transformer_apparatus_event!(
        air_core_runtime,
        TransformerPublic.TransformerApparatusEventCommand(
            :shunt_reactor_breaker_opens,
            TransformerPublic.TransformerBreakerOpenEvent,
            10.0e-6;
            target_indices=[1, 2],
        ),
    )
    air_core_surface = only(nonlinear_device_event_surfaces(air_core_runtime))
    air_core_current, _ = accept_transformer_test_step!(
        air_core_runtime,
        [100.0, 0.0],
    )
    @test air_core_current[1] > 0.0
    @test air_core_current[2] < 0.0
    apply_nonlinear_device_event!(air_core_surface, air_core_runtime, 10.0e-6)
    opened_current, opened_jacobian = accept_transformer_test_step!(
        air_core_runtime,
        [100.0, 0.0],
        20.0e-6,
        companion_method=:backward_euler,
    )
    @test opened_current == zeros(2)
    @test opened_jacobian == zeros(2, 2)
    air_core_result = TransformerPublic.transformer_apparatus_result(air_core_runtime)
    @test air_core_result.reactor_definition == air_core_definition
    @test !TransformerPublic.transformer_result_quantity_available(
        air_core_result.magnetic_branch_flux_wb,
    )
    application_signatures = String[]
    for application in instances(TransformerPublic.ReactorApplication)
        definition = TransformerPublic.ReactorApparatusDefinition(
            application,
            TransformerPublic.AirCoreReactorConstruction,
        )
        specification = transformer_specification(
            TransformerPublic.LowFrequencyTerminalTier,
            connection,
            TransformerPublic.LowFrequencyTransformerModel(terminal_matrices),
            source;
            reactor_definition=definition,
        )
        runtime = TransformerPublic.transformer_apparatus_runtime(
            TransformerPublic.prepare_transformer_apparatus(specification),
            [1, 2],
        )
        accept_transformer_test_step!(runtime, [10.0, 0.0])
        result = TransformerPublic.transformer_apparatus_result(runtime)
        @test result.reactor_definition.application == application
        push!(application_signatures, result.deterministic_signature_sha256)
    end
    @test length(unique(application_signatures)) == 4

    split_connection = split_winding_reactor_connection()
    split_terminal_matrices = TransformerPublic.TransformerTerminalMatrices(
        [0.12 0.0; 0.0 0.13],
        [4.0e-3 0.0; 0.0 4.5e-3];
        capacitance_f=[2.0e-9 0.0; 0.0 2.0e-9],
        conductance_s=[1.0e-7 0.0; 0.0 1.0e-7],
    )
    split_definition = TransformerPublic.ReactorApparatusDefinition(
        TransformerPublic.SeriesReactorApplication,
        TransformerPublic.AirCoreReactorConstruction;
        winding_configuration=TransformerPublic.SplitWindingReactorWinding,
    )
    split_specification = transformer_specification(
        TransformerPublic.LowFrequencyTerminalTier,
        split_connection,
        TransformerPublic.LowFrequencyTransformerModel(split_terminal_matrices),
        source;
        reactor_definition=split_definition,
    )
    split_runtime = TransformerPublic.transformer_apparatus_runtime(
        TransformerPublic.prepare_transformer_apparatus(split_specification),
        [1, 2, 3],
    )
    split_current, split_jacobian = accept_transformer_test_step!(
        split_runtime,
        [100.0, 40.0, 0.0],
    )
    @test maximum(abs, split_current) > 0.0
    @test sum(split_current) ≈ 0.0 atol=2.0e-12
    @test split_jacobian ≈ transpose(split_jacobian) atol=2.0e-12
    @test TransformerPublic.transformer_apparatus_result(split_runtime).
        reactor_definition == split_definition

    mutually_coupled_terminal_matrices =
        TransformerPublic.TransformerTerminalMatrices(
            [0.12 0.0; 0.0 0.13],
            [4.0e-3 1.0e-3; 1.0e-3 4.5e-3];
            capacitance_f=[2.0e-9 0.0; 0.0 2.0e-9],
            conductance_s=[1.0e-7 0.0; 0.0 1.0e-7],
        )
    mutually_coupled_definition = TransformerPublic.ReactorApparatusDefinition(
        TransformerPublic.NeutralReactorApplication,
        TransformerPublic.AirCoreReactorConstruction;
        winding_configuration=TransformerPublic.MutuallyCoupledReactorWinding,
    )
    mutually_coupled_specification = transformer_specification(
        TransformerPublic.LowFrequencyTerminalTier,
        split_connection,
        TransformerPublic.LowFrequencyTransformerModel(
            mutually_coupled_terminal_matrices,
        ),
        source;
        reactor_definition=mutually_coupled_definition,
    )
    mutually_coupled_runtime = TransformerPublic.transformer_apparatus_runtime(
        TransformerPublic.prepare_transformer_apparatus(
            mutually_coupled_specification,
        ),
        [1, 2, 3],
    )
    mutually_coupled_current, mutually_coupled_jacobian =
        accept_transformer_test_step!(
            mutually_coupled_runtime,
            [100.0, 40.0, 0.0],
        )
    @test mutually_coupled_current != split_current
    @test sum(mutually_coupled_current) ≈ 0.0 atol=2.0e-12
    @test mutually_coupled_jacobian ≈
        transpose(mutually_coupled_jacobian) atol=2.0e-12
    @test mutually_coupled_specification.deterministic_signature_sha256 !=
        split_specification.deterministic_signature_sha256
    @test TransformerPublic.transformer_apparatus_result(
        mutually_coupled_runtime,
    ).reactor_definition == mutually_coupled_definition
    @test_throws ArgumentError transformer_specification(
        TransformerPublic.LowFrequencyTerminalTier,
        split_connection,
        TransformerPublic.LowFrequencyTransformerModel(split_terminal_matrices),
        source;
        reactor_definition=mutually_coupled_definition,
    )

    iron_material = TransformerPublic.PiecewiseLinearTransformerMagneticMaterial(
        [0.0, 0.5, 1.5, 2.0],
        [0.0, 100.0, 2_000.0, 10_000.0],
        source,
    )
    iron_graph = TransformerPublic.TransformerMagneticGraph(
        node_order=[:reactor_magnetic_node],
        branches=[
            TransformerPublic.MagneticBranchGeometry(
                :gapped_reactor_limb;
                length_m=0.8,
                cross_section_m2=0.01,
                air_gap_length_m=2.0e-3,
                air_gap_effective_area_factor=1.2,
            ),
            TransformerPublic.MagneticBranchGeometry(
                :reactor_return_limb;
                length_m=0.8,
                cross_section_m2=0.01,
            ),
        ],
        incidence=reshape([1.0, -1.0], 1, 2),
        winding_turns=reshape([100.0, 0.0], 2, 1),
        materials=[iron_material],
    )
    iron_definition = TransformerPublic.ReactorApparatusDefinition(
        TransformerPublic.SmoothingReactorApplication,
        TransformerPublic.IronCoreReactorConstruction;
        gap_model=TransformerPublic.EffectiveAreaReactorAirGap,
        total_air_gap_length_m=2.0e-3,
        air_gap_effective_area_factor=1.2,
    )
    iron_model = TransformerPublic.MagneticEquivalentCircuitModel(
        reshape([0.25], 1, 1),
        reshape([1.0e-3], 1, 1),
        iron_graph;
        terminal_capacitance_f=reshape([2.0e-9], 1, 1),
        terminal_conductance_s=reshape([1.0e-7], 1, 1),
    )
    iron_specification = transformer_specification(
        TransformerPublic.MagneticEquivalentCircuitTier,
        connection,
        iron_model,
        source;
        reactor_definition=iron_definition,
    )
    iron_runtime = TransformerPublic.transformer_apparatus_runtime(
        TransformerPublic.prepare_transformer_apparatus(iron_specification),
        [1, 2],
    )
    iron_current, iron_jacobian = accept_transformer_test_step!(
        iron_runtime,
        [25.0, 0.0],
    )
    @test all(isfinite, iron_current)
    @test iron_jacobian[1, 1] > 0.0
    @test iron_jacobian[1, 2] < 0.0
    iron_result = TransformerPublic.transformer_apparatus_result(iron_runtime)
    @test iron_result.reactor_definition == iron_definition
    @test TransformerPublic.transformer_result_quantity_available(
        iron_result.magnetic_branch_flux_wb,
    )
    @test iron_result.energy.stored_energy_j > 0.0

    @test_throws ArgumentError TransformerPublic.ReactorApparatusDefinition(
        TransformerPublic.SeriesReactorApplication,
        TransformerPublic.AirCoreReactorConstruction;
        gap_model=TransformerPublic.UniformReactorAirGap,
        total_air_gap_length_m=1.0e-3,
    )
    @test_throws ArgumentError transformer_specification(
        TransformerPublic.LowFrequencyTerminalTier,
        grounded_transformer_connection(),
        TransformerPublic.LowFrequencyTransformerModel(TransformerPublic.TransformerTerminalMatrices(
            [0.4 0.0; 0.0 0.5],
            [2.0e-3 0.0; 0.0 3.0e-3],
        )),
        source;
        reactor_definition=air_core_definition,
    )
    @test_throws ArgumentError transformer_specification(
        TransformerPublic.LowFrequencyTerminalTier,
        connection,
        TransformerPublic.LowFrequencyTransformerModel(terminal_matrices),
        source;
        reactor_definition=iron_definition,
    )
end

@testset "nonlinear transformer magnetic companion and rollback boundary" begin
    source = synthetic_transformer_source()
    connection = grounded_transformer_connection()
    graph = nonlinear_transformer_magnetic_graph(source)
    leakage = TransformerPublic.TransformerTerminalMatrices(
        [0.4 0.0; 0.0 0.5],
        [2.0e-3 0.1e-3; 0.1e-3 3.0e-3];
        capacitance_f=[2.0e-9 0.0; 0.0 1.0e-9],
        conductance_s=[1.0e-7 0.0; 0.0 2.0e-7],
    )
    model = TransformerPublic.HybridTransformerModel(leakage, graph)
    specification = transformer_specification(
        TransformerPublic.HybridTransformerTier,
        connection,
        model,
        source,
    )
    preparation = TransformerPublic.prepare_transformer_apparatus(specification)
    @test preparation.effective_linear_inductance_h === nothing
    runtime = TransformerPublic.transformer_apparatus_runtime(preparation, [1, 2])
    accepted_flux = copy(runtime.accepted_state.branch_flux_wb)
    current = zeros(2)
    jacobian = zeros(2, 2)
    prepare_nonlinear_device_step!(runtime, 10.0e-6, 10.0e-6, :trapezoidal)
    nonlinear_current_jacobian!(
        current,
        jacobian,
        runtime,
        [10.0, -4.0],
        10.0e-6,
    )
    @test runtime.accepted_state.branch_flux_wb == accepted_flux
    @test all(isfinite, current)
    @test jacobian ≈ transpose(jacobian) atol=1.0e-10
    @test minimum(eigvals(Symmetric(jacobian))) >= -1.0e-12
    perturbation_v = 1.0e-5
    positive_current = zeros(2)
    negative_current = zeros(2)
    scratch_jacobian = zeros(2, 2)
    for terminal in 1:2
        positive_voltage = [10.0, -4.0]
        negative_voltage = [10.0, -4.0]
        positive_voltage[terminal] += perturbation_v
        negative_voltage[terminal] -= perturbation_v
        nonlinear_current_jacobian!(
            positive_current,
            scratch_jacobian,
            runtime,
            positive_voltage,
            10.0e-6,
        )
        nonlinear_current_jacobian!(
            negative_current,
            scratch_jacobian,
            runtime,
            negative_voltage,
            10.0e-6,
        )
        finite_difference = (positive_current .- negative_current) ./
            (2.0 * perturbation_v)
        @test finite_difference ≈ jacobian[:, terminal] rtol=2.0e-6 atol=1.0e-9
    end
    nonlinear_current_jacobian!(
        current,
        jacobian,
        runtime,
        [10.0, -4.0],
        10.0e-6,
    )
    accept_nonlinear_device_state!(
        runtime,
        [10.0, -4.0],
        current,
        jacobian,
        10.0e-6,
    )
    finish_nonlinear_device_step!(runtime)
    state = runtime.accepted_state
    @test maximum(abs, state.branch_flux_wb) > 0.0
    @test state.maximum_magnetic_continuity_residual_wb <= 1.0e-12
    @test state.maximum_magnetic_constitutive_residual_at <= 1.0e-6
    @test state.maximum_local_nonlinear_iterations <= 10
    @test state.stored_leakage_energy_j > 0.0
    @test state.stored_magnetic_energy_j > 0.0
    snapshot = TransformerPublic.transformer_apparatus_runtime_snapshot(runtime)
    restored = TransformerPublic.transformer_apparatus_runtime(preparation, [1, 2])
    TransformerPublic.restore_transformer_apparatus_runtime_snapshot!(restored, snapshot)
    @test TransformerPublic.transformer_apparatus_runtime_snapshot(restored).
        deterministic_signature_sha256 == snapshot.deterministic_signature_sha256

    mec_model = TransformerPublic.MagneticEquivalentCircuitModel(
        leakage.resistance_ohm,
        leakage.inductance_h,
        graph;
        terminal_capacitance_f=leakage.capacitance_f,
        terminal_conductance_s=leakage.conductance_s,
    )
    mec_specification = transformer_specification(
        TransformerPublic.MagneticEquivalentCircuitTier,
        connection,
        mec_model,
        source,
    )
    mec_runtime = TransformerPublic.transformer_apparatus_runtime(
        TransformerPublic.prepare_transformer_apparatus(mec_specification),
        [1, 2],
    )
    mec_current, mec_jacobian = accept_transformer_test_step!(
        mec_runtime,
        [10.0, -4.0],
    )
    @test mec_current ≈ current rtol=1.0e-10 atol=1.0e-12
    @test mec_jacobian ≈ jacobian rtol=1.0e-10 atol=1.0e-12
end

@testset "Tellinen transformer hysteresis accepted-state trajectory" begin
    source = synthetic_transformer_source()
    connection = grounded_transformer_connection()
    graph = hysteretic_transformer_magnetic_graph(source)
    leakage = TransformerPublic.TransformerTerminalMatrices(
        [0.2 0.0; 0.0 0.2],
        [1.0e-3 0.0; 0.0 1.0e-3],
    )
    model = TransformerPublic.HybridTransformerModel(leakage, graph)
    specification = transformer_specification(
        TransformerPublic.HybridTransformerTier,
        connection,
        model,
        source,
    )
    preparation = TransformerPublic.prepare_transformer_apparatus(specification)
    runtime = TransformerPublic.transformer_apparatus_runtime(preparation, [1, 2])
    initial_states = copy(runtime.accepted_state.tellinen_state)
    accept_transformer_test_step!(runtime, [5.0, 5.0], 10.0e-6)
    rising_states = copy(runtime.accepted_state.tellinen_state)
    @test all(state -> state !== nothing, rising_states)
    @test all(
        index -> rising_states[index].flux_density_t >
            initial_states[index].flux_density_t,
        eachindex(rising_states),
    )
    accept_transformer_test_step!(runtime, [-10.0, -10.0], 20.0e-6)
    falling_states = runtime.accepted_state.tellinen_state
    @test all(
        index -> falling_states[index].direction == -1,
        eachindex(falling_states),
    )
    @test all(
        index -> falling_states[index].reversal_count >= 1,
        eachindex(falling_states),
    )
    @test all(isfinite, (
        runtime.accepted_state.stored_magnetic_energy_j,
        runtime.accepted_state.hysteresis_loss_energy_j,
        runtime.accepted_state.supplied_energy_j,
    ))
    snapshot = TransformerPublic.transformer_apparatus_runtime_snapshot(runtime)
    restored = TransformerPublic.transformer_apparatus_runtime(preparation, [1, 2])
    TransformerPublic.restore_transformer_apparatus_runtime_snapshot!(restored, snapshot)
    @test TransformerPublic.transformer_apparatus_runtime_snapshot(restored).
        deterministic_signature_sha256 == snapshot.deterministic_signature_sha256
end

@testset "separate classical eddy and excess transformer core loss" begin
    source = synthetic_transformer_source()
    connection = grounded_transformer_connection()
    @test_throws ArgumentError TransformerPublic.TransformerDynamicCoreLossModel(
        -1.0,
        0.0,
        source,
    )
    dynamic_loss = TransformerPublic.TransformerDynamicCoreLossModel(
        0.02,
        0.05,
        source;
        excess_rate_regularization_t_per_s=1.0e-3,
    )
    branches = [
        TransformerPublic.MagneticBranchGeometry(
            :dynamic_core_limb;
            length_m=0.8,
            cross_section_m2=0.01,
        ),
        TransformerPublic.MagneticBranchGeometry(
            :dynamic_return_limb;
            length_m=0.8,
            cross_section_m2=0.01,
        ),
    ]
    material = TransformerPublic.LinearTransformerMagneticMaterial(
        2_000.0,
        2.0,
        source,
    )
    graph = TransformerPublic.TransformerMagneticGraph(
        node_order=[:magnetic_node],
        branches=branches,
        incidence=reshape([1.0, -1.0], 1, 2),
        winding_turns=[100.0 0.0; 0.0 100.0],
        materials=[material],
        dynamic_core_loss=[dynamic_loss, dynamic_loss],
    )
    @test graph.dynamic_core_loss == [dynamic_loss, dynamic_loss]
    @test TransformerPublic.transformer_magnetic_graph_signature(graph) !=
        TransformerPublic.transformer_magnetic_graph_signature(
            linear_transformer_magnetic_graph(source),
        )
    terminal_matrices = TransformerPublic.TransformerTerminalMatrices(
        [0.2 0.0; 0.0 0.2],
        [1.0e-3 0.0; 0.0 1.0e-3],
    )
    model = TransformerPublic.HybridTransformerModel(terminal_matrices, graph)
    specification = transformer_specification(
        TransformerPublic.HybridTransformerTier,
        connection,
        model,
        source,
    )
    runtime = TransformerPublic.transformer_apparatus_runtime(
        TransformerPublic.prepare_transformer_apparatus(specification),
        [1, 2],
    )
    accept_transformer_test_step!(runtime, [5.0, 5.0], 10.0e-6)
    state = runtime.accepted_state
    @test state.hysteresis_loss_energy_j == 0.0
    @test state.classical_eddy_core_loss_energy_j > 0.0
    @test state.excess_core_loss_energy_j > 0.0
    @test maximum(abs, state.last_classical_dynamic_mmf_at) > 0.0
    @test maximum(abs, state.last_excess_dynamic_mmf_at) > 0.0
    @test abs(
        state.supplied_energy_j - state.stored_leakage_energy_j -
        state.stored_magnetic_energy_j - state.stored_electric_energy_j -
        state.winding_loss_energy_j - state.dielectric_loss_energy_j -
        state.classical_eddy_core_loss_energy_j - state.excess_core_loss_energy_j,
    ) <= 1.0e-9
    diagnostics = TransformerPublic.transformer_apparatus_runtime_diagnostics(runtime)
    @test diagnostics.classical_eddy_core_loss_energy_j ==
        state.classical_eddy_core_loss_energy_j
    @test diagnostics.excess_core_loss_energy_j == state.excess_core_loss_energy_j
    snapshot = TransformerPublic.transformer_apparatus_runtime_snapshot(runtime)
    restored = TransformerPublic.transformer_apparatus_runtime(
        TransformerPublic.prepare_transformer_apparatus(specification),
        [1, 2],
    )
    TransformerPublic.restore_transformer_apparatus_runtime_snapshot!(restored, snapshot)
    @test restored.accepted_state.classical_eddy_core_loss_energy_j ==
        state.classical_eddy_core_loss_energy_j
    @test restored.accepted_state.excess_core_loss_energy_j ==
        state.excess_core_loss_energy_j
end

@testset "transformer apparatus tiers, topology, and magnetic owners" begin
    source = synthetic_transformer_source()
    @test length(TransformerPublic.TRANSFORMER_APPARATUS_TIERS) == 7
    @test length(unique(TransformerPublic.TRANSFORMER_APPARATUS_TIERS)) == 7
    for tier in TransformerPublic.TRANSFORMER_APPARATUS_TIERS
        contract = TransformerPublic.transformer_apparatus_contract(tier)
        @test contract.maturity === :implemented
        @test :atp_or_pscad_equivalence in
            contract.validity_domain.unsupported_phenomena
        @test contract.mutation_order[end] ===
            :accept_state_history_energy_and_output_once_or_restore
    end

    connection = grounded_transformer_connection()
    @test TransformerPublic.transformer_connection_rank(connection) == 2
    @test length(TransformerPublic.transformer_connection_signature(connection)) == 64
    node_voltage = [11.0, -7.0]
    coil_voltage = zeros(2)
    TransformerPublic.transformer_coil_voltages!(
        coil_voltage,
        connection,
        node_voltage,
    )
    @test coil_voltage == node_voltage
    coil_current = [3.0, -2.0]
    terminal_current = zeros(2)
    TransformerPublic.transformer_terminal_currents!(
        terminal_current,
        connection,
        coil_current,
    )
    @test terminal_current == coil_current
    @test dot(node_voltage, terminal_current) == dot(coil_voltage, coil_current)
    delta_star_incidence = [
        1.0 0.0 -1.0 0.0 0.0 0.0
        -1.0 1.0 0.0 0.0 0.0 0.0
        0.0 -1.0 1.0 0.0 0.0 0.0
        0.0 0.0 0.0 1.0 0.0 0.0
        0.0 0.0 0.0 0.0 1.0 0.0
        0.0 0.0 0.0 0.0 0.0 1.0
    ]
    delta_star = TransformerPublic.TransformerConnectionTopology(
        node_order=[:primary_a, :primary_b, :primary_c, :secondary_a, :secondary_b, :secondary_c],
        node_phase=[:phase_a, :phase_b, :phase_c, :phase_a, :phase_b, :phase_c],
        coil_order=[:primary_a, :primary_b, :primary_c, :secondary_a, :secondary_b, :secondary_c],
        winding_order=[:primary, :secondary],
        phase_order=[:phase_a, :phase_b, :phase_c],
        coil_winding=[:primary, :primary, :primary, :secondary, :secondary, :secondary],
        coil_phase=[:phase_a, :phase_b, :phase_c, :phase_a, :phase_b, :phase_c],
        incidence=delta_star_incidence,
        vector_group="Dyn11",
        clock_number=11,
        phase_shift_rad=-pi / 6.0,
    )
    @test TransformerPublic.transformer_connection_phase_shift(delta_star) == -pi / 6.0
    @test_throws ArgumentError TransformerPublic.TransformerConnectionTopology(
        node_order=[:primary_a, :primary_b, :primary_c, :secondary_a, :secondary_b, :secondary_c],
        node_phase=[:phase_a, :phase_b, :phase_c, :phase_a, :phase_b, :phase_c],
        coil_order=[:primary_a, :primary_b, :primary_c, :secondary_a, :secondary_b, :secondary_c],
        winding_order=[:primary, :secondary],
        phase_order=[:phase_a, :phase_b, :phase_c],
        coil_winding=[:primary, :primary, :primary, :secondary, :secondary, :secondary],
        coil_phase=[:phase_a, :phase_b, :phase_c, :phase_a, :phase_b, :phase_c],
        incidence=delta_star_incidence,
        vector_group="Dyn0",
        clock_number=0,
        phase_shift_rad=0.0,
    )
    @test_throws ArgumentError TransformerPublic.TransformerConnectionTopology(
        node_order=[:terminal],
        coil_order=[:coil],
        winding_order=[:winding],
        phase_order=[:phase_a],
        coil_winding=[:winding],
        coil_phase=[:phase_a],
        incidence=reshape([1.0], 1, 1),
        vector_group="Ii1",
        clock_number=1,
        phase_shift_rad=0.0,
    )

    graph = linear_transformer_magnetic_graph(source)
    response = TransformerPublic.transformer_magnetic_linear_response(
        graph,
        [1.0, 1.0],
    )
    @test response.continuity_residual_wb <= 1.0e-14
    @test response.constitutive_residual_at <= 1.0e-10
    @test response.stored_energy_j > 0.0
    inductance = TransformerPublic.transformer_magnetic_linear_inductance(graph)
    @test inductance ≈ transpose(inductance) atol=1.0e-14
    @test minimum(eigvals(Symmetric(inductance))) >= -1.0e-12

    piecewise = TransformerPublic.PiecewiseLinearTransformerMagneticMaterial(
        [0.0, 1.0, 1.8],
        [0.0, 100.0, 1_000.0],
        source,
    )
    @test TransformerPublic.magnetic_material_field(piecewise, -0.5) == -50.0
    @test TransformerPublic.magnetic_material_differential_reluctivity(
        piecewise,
        1.4,
    ) ≈ 1_125.0

    lower_curve = TransformerPublic.TellinenLimitingCurve(
        [-1_000.0, -500.0, 0.0, 500.0, 1_000.0],
        [-1.5, -1.0, -0.2, 0.5, 1.0],
    )
    upper_curve = TransformerPublic.TellinenLimitingCurve(
        [-1_000.0, -500.0, 0.0, 500.0, 1_000.0],
        [-1.0, -0.5, 0.2, 1.0, 1.5],
    )
    tellinen = TransformerPublic.TellinenTransformerMagneticMaterial(
        lower_curve,
        upper_curve,
        source;
        integration_field_increment_a_per_m=2.0,
    )
    initial_tellinen = TransformerPublic.tellinen_state(tellinen)
    tellinen_trial = TransformerPublic.tellinen_trial_from_flux_density(
        tellinen,
        initial_tellinen,
        0.05,
    )
    @test tellinen_trial.state.flux_density_t == 0.05
    @test tellinen_trial.state.field_strength_a_per_m > 0.0
    @test tellinen_trial.differential_reluctivity_m_per_h > 0.0
    @test initial_tellinen.reversal_count == 0
end

@testset "seven transformer apparatus runtime tiers" begin
    source = synthetic_transformer_source()
    connection = grounded_transformer_connection()
    graph = linear_transformer_magnetic_graph(source)
    terminal_matrices = TransformerPublic.TransformerTerminalMatrices(
        [0.4 0.0; 0.0 0.5],
        [0.12 0.02; 0.02 0.10];
        capacitance_f=[2.0e-9 0.0; 0.0 1.0e-9],
        conductance_s=[1.0e-7 0.0; 0.0 2.0e-7],
    )
    models = (
        TransformerPublic.LowFrequencyTransformerModel(terminal_matrices),
        TransformerPublic.BCTRANTransformerModel(
            terminal_matrices;
            positive_pair_reconstruction_residual=0.0,
            zero_pair_reconstruction_residual=0.0,
            inverse_reconstruction_residual=0.0,
        ),
        TransformerPublic.HybridTransformerModel(terminal_matrices, graph),
        TransformerPublic.MagneticEquivalentCircuitModel(
            terminal_matrices.resistance_ohm,
            terminal_matrices.inductance_h,
            graph;
            terminal_capacitance_f=terminal_matrices.capacitance_f,
            terminal_conductance_s=terminal_matrices.conductance_s,
        ),
    )
    for (tier, model) in zip(
        TransformerPublic.TRANSFORMER_APPARATUS_TIERS[1:4],
        models,
    )
        specification = transformer_specification(tier, connection, model, source)
        preparation = TransformerPublic.prepare_transformer_apparatus(specification)
        readiness = TransformerPublic.transformer_apparatus_readiness(preparation)
        @test readiness.ready
        @test !readiness.production_backend_available
        runtime = TransformerPublic.transformer_apparatus_runtime(
            preparation,
            [1, 2],
        )
        current, jacobian = accept_transformer_test_step!(runtime, [10.0, -4.0])
        @test all(isfinite, current)
        @test jacobian ≈ transpose(jacobian) atol=1.0e-12
        diagnostics = TransformerPublic.transformer_apparatus_runtime_diagnostics(runtime)
        @test diagnostics.accepted_step_count == 1
        @test diagnostics.accepted_time_s == 10.0e-6
        result = TransformerPublic.transformer_apparatus_result(runtime)
        @test result.schema_version == 1
        @test result.tier == tier
        @test result.accepted_step_count == 1
        @test result.terminal_voltage_v == (10.0, -4.0)
        @test length(result.deterministic_signature_sha256) == 64
        @test TransformerPublic.transformer_result_quantity_available(
            result.coil_current_a,
        )
        magnetic_output_available = tier in (
            TransformerPublic.HybridTransformerTier,
            TransformerPublic.MagneticEquivalentCircuitTier,
        )
        @test TransformerPublic.transformer_result_quantity_available(
            result.magnetic_branch_flux_wb,
        ) == magnetic_output_available
        @test abs(result.energy.unexplained_balance_residual_j) <=
            specification.settings.energy_absolute_tolerance_j
        @test TransformerPublic.transformer_apparatus_result(runtime).
            deterministic_signature_sha256 == result.deterministic_signature_sha256
        state = runtime.accepted_state
        @test state.winding_loss_energy_j >= 0.0
        @test state.dielectric_loss_energy_j >= 0.0
        stored_inductive_energy_j = hasproperty(state, :stored_leakage_energy_j) ?
            state.stored_leakage_energy_j + state.stored_magnetic_energy_j :
            state.stored_magnetic_energy_j
        @test abs(
            state.supplied_energy_j - stored_inductive_energy_j -
            state.stored_electric_energy_j - state.winding_loss_energy_j -
            state.dielectric_loss_energy_j,
        ) <= 1.0e-10
        snapshot = TransformerPublic.transformer_apparatus_runtime_snapshot(runtime)
        restored = TransformerPublic.transformer_apparatus_runtime(preparation, [1, 2])
        @test TransformerPublic.restore_transformer_apparatus_runtime_snapshot!(
            restored,
            snapshot,
        ) === restored
        @test TransformerPublic.transformer_apparatus_runtime_snapshot(restored).
            deterministic_signature_sha256 == snapshot.deterministic_signature_sha256
    end

    @test_throws ArgumentError TransformerPublic.HybridTransformerModel(
        terminal_matrices,
        graph;
        winding_loss_state_matrix_per_s=-1_000.0 .* Matrix{Float64}(I, 2, 2),
        winding_loss_input_matrix=Matrix{Float64}(I, 2, 2),
        winding_loss_output_matrix_ohm_per_s=Matrix{Float64}(I, 2, 2),
        winding_loss_direct_ohm=0.05 .* Matrix{Float64}(I, 2, 2),
        winding_loss_storage_matrix_j=2.0 .* Matrix{Float64}(I, 2, 2),
        continuous_passivity_margin_ohm=0.05,
    )
    dynamic_loss_hybrid = TransformerPublic.HybridTransformerModel(
        terminal_matrices,
        graph;
        winding_loss_state_matrix_per_s=-1_000.0 .* Matrix{Float64}(I, 2, 2),
        winding_loss_input_matrix=Matrix{Float64}(I, 2, 2),
        winding_loss_output_matrix_ohm_per_s=Matrix{Float64}(I, 2, 2),
        winding_loss_direct_ohm=0.05 .* Matrix{Float64}(I, 2, 2),
        winding_loss_storage_matrix_j=Matrix{Float64}(I, 2, 2),
        continuous_passivity_margin_ohm=0.05,
    )
    dynamic_loss_specification = transformer_specification(
        TransformerPublic.HybridTransformerTier,
        connection,
        dynamic_loss_hybrid,
        source,
    )
    dynamic_loss_runtime = TransformerPublic.transformer_apparatus_runtime(
        TransformerPublic.prepare_transformer_apparatus(dynamic_loss_specification),
        [1, 2],
    )
    dynamic_loss_current, dynamic_loss_jacobian = accept_transformer_test_step!(
        dynamic_loss_runtime,
        [10.0, -4.0],
    )
    @test all(isfinite, dynamic_loss_current)
    @test dynamic_loss_jacobian ≈ transpose(dynamic_loss_jacobian) atol=1.0e-12
    dynamic_loss_state = dynamic_loss_runtime.accepted_state
    @test maximum(abs, dynamic_loss_state.winding_loss_state) > 0.0
    @test dynamic_loss_state.stored_frequency_dependent_winding_energy_j > 0.0
    @test dynamic_loss_state.frequency_dependent_winding_loss_energy_j > 0.0
    @test abs(
        dynamic_loss_state.supplied_energy_j -
        dynamic_loss_state.stored_leakage_energy_j -
        dynamic_loss_state.stored_magnetic_energy_j -
        dynamic_loss_state.stored_frequency_dependent_winding_energy_j -
        dynamic_loss_state.stored_electric_energy_j -
        dynamic_loss_state.winding_loss_energy_j -
        dynamic_loss_state.frequency_dependent_winding_loss_energy_j -
        dynamic_loss_state.dielectric_loss_energy_j,
    ) <= 1.0e-10

    @test_throws ArgumentError TransformerPublic.WidebandTransformerModel(
        -1_000.0 .* Matrix{Float64}(I, 2, 2),
        Matrix{Float64}(I, 2, 2),
        Matrix{Float64}(I, 2, 2),
        0.1 .* Matrix{Float64}(I, 2, 2);
        storage_matrix_j=2.0 .* Matrix{Float64}(I, 2, 2),
        port_order=connection.node_order,
        frequency_band_hz=(10.0, 10_000.0),
        continuous_passivity_margin_s=0.0,
        enforcement_perturbation_relative=0.0,
        source_response_sha256=repeat("b", 64),
    )
    wideband = TransformerPublic.WidebandTransformerModel(
        -1_000.0 .* Matrix{Float64}(I, 2, 2),
        Matrix{Float64}(I, 2, 2),
        Matrix{Float64}(I, 2, 2),
        0.1 .* Matrix{Float64}(I, 2, 2);
        storage_matrix_j=Matrix{Float64}(I, 2, 2),
        port_order=connection.node_order,
        frequency_band_hz=(10.0, 10_000.0),
        continuous_passivity_margin_s=0.0,
        enforcement_perturbation_relative=0.0,
        source_response_sha256=repeat("b", 64),
    )
    wideband_specification = transformer_specification(
        TransformerPublic.WidebandBlackBoxTier,
        connection,
        wideband,
        source,
    )
    wideband_runtime = TransformerPublic.transformer_apparatus_runtime(
        TransformerPublic.prepare_transformer_apparatus(wideband_specification),
        [1, 2],
    )
    wideband_current, wideband_jacobian =
        accept_transformer_test_step!(wideband_runtime, [1.0, -1.0])
    @test dot([1.0, -1.0], wideband_current) > 0.0
    @test wideband_jacobian ≈ transpose(wideband_jacobian) atol=1.0e-12
    @test wideband_runtime.accepted_state.stored_energy_j > 0.0
    @test wideband_runtime.accepted_state.dissipated_energy_j > 0.0
    @test wideband_runtime.accepted_state.maximum_energy_residual_j <= 1.0e-10
    @test abs(
        wideband_runtime.accepted_state.supplied_energy_j -
        wideband_runtime.accepted_state.stored_energy_j -
        wideband_runtime.accepted_state.dissipated_energy_j,
    ) <= 1.0e-10
    wideband_result = TransformerPublic.transformer_apparatus_result(wideband_runtime)
    @test TransformerPublic.transformer_result_quantity_available(
        wideband_result.passive_rational_state,
    )
    @test !TransformerPublic.transformer_result_quantity_available(
        wideband_result.coil_current_a,
    )
    @test_throws ArgumentError TransformerPublic.transformer_result_quantity_value(
        wideband_result.coil_current_a,
    )

    grey_box = TransformerPublic.GreyBoxTransformerModel(
        node_order=[:primary_terminal, :secondary_terminal, :internal_node],
        terminal_node_indices=[1, 2],
        branches=[
            TransformerPublic.TransformerLadderBranch(
                :primary_section,
                1,
                3;
                resistance_ohm=0.2,
                inductance_h=2.0e-3,
            ),
            TransformerPublic.TransformerLadderBranch(
                :secondary_section,
                3,
                2;
                resistance_ohm=0.3,
                inductance_h=3.0e-3,
            ),
        ],
        capacitance_f=Diagonal([1.0e-9, 1.0e-9, 2.0e-9]),
        conductance_s=Diagonal([1.0e-7, 1.0e-7, 2.0e-7]),
        source_response_sha256=repeat("c", 64),
        identification_residual_relative=1.0e-4,
        parameter_nonuniqueness="synthetic topology selected from equivalent fits",
    )
    grey_specification = transformer_specification(
        TransformerPublic.GreyBoxLadderTier,
        connection,
        grey_box,
        source,
    )
    grey_runtime = TransformerPublic.transformer_apparatus_runtime(
        TransformerPublic.prepare_transformer_apparatus(grey_specification),
        [1, 2],
    )
    grey_current, grey_jacobian = accept_transformer_test_step!(
        grey_runtime,
        [10.0, 0.0],
    )
    @test all(isfinite, grey_current)
    @test grey_jacobian ≈ transpose(grey_jacobian) atol=1.0e-12
    @test grey_runtime.accepted_state.maximum_internal_kcl_residual_a <= 1.0e-12
    @test grey_runtime.accepted_state.dielectric_loss_energy_j >= 0.0
    grey_result = TransformerPublic.transformer_apparatus_result(grey_runtime)
    @test TransformerPublic.transformer_result_quantity_available(
        grey_result.ladder_node_voltage_v,
    )
    @test !TransformerPublic.transformer_result_quantity_available(
        grey_result.winding_section_voltage_v,
    )

    @test_throws ArgumentError TransformerPublic.WhiteBoxTransformerModel(
        winding_order=[:primary_winding, :secondary_winding],
        section_count_per_winding=[1, 3],
        section_length_m=[0.5, 0.5, 0.5],
        series_resistance_ohm_per_m=fill(Matrix{Float64}(I, 2, 2), 3),
        series_inductance_h_per_m=fill(Matrix{Float64}(I, 2, 2), 3),
        shunt_conductance_s_per_m=fill(zeros(2, 2), 3),
        shunt_capacitance_f_per_m=fill(zeros(2, 2), 3),
        geometry_sha256=repeat("d", 64),
        frequency_band_hz=(10.0, 10_000.0),
        section_refinement_residual_relative=0.0,
    )
    series_resistance = fill([0.2 0.02; 0.02 0.3], 3)
    series_inductance = fill([2.0e-3 0.2e-3; 0.2e-3 3.0e-3], 3)
    white_box = TransformerPublic.WhiteBoxTransformerModel(
        winding_order=[:primary_winding, :secondary_winding],
        section_count_per_winding=[2, 3],
        section_length_m=[0.5, 0.5, 0.5],
        series_resistance_ohm_per_m=series_resistance,
        series_inductance_h_per_m=series_inductance,
        shunt_conductance_s_per_m=fill(1.0e-7 .* Matrix{Float64}(I, 2, 2), 3),
        shunt_capacitance_f_per_m=fill(1.0e-9 .* Matrix{Float64}(I, 2, 2), 3),
        geometry_sha256=repeat("d", 64),
        frequency_band_hz=(10.0, 10_000.0),
        section_refinement_residual_relative=1.0e-3,
    )
    white_connection = white_box_transformer_connection()
    white_specification = transformer_specification(
        TransformerPublic.WhiteBoxWindingTier,
        white_connection,
        white_box,
        source,
    )
    white_runtime = TransformerPublic.transformer_apparatus_runtime(
        TransformerPublic.prepare_transformer_apparatus(white_specification),
        [1, 2, 3, 4],
    )
    white_current, white_jacobian = accept_transformer_test_step!(
        white_runtime,
        [10.0, 0.0, 5.0, 0.0],
    )
    @test all(isfinite, white_current)
    @test white_jacobian ≈ transpose(white_jacobian) atol=1.0e-12
    @test minimum(eigvals(Symmetric(white_jacobian))) >= -1.0e-12
    @test white_runtime.accepted_state.maximum_internal_kcl_residual_a <= 1.0e-12
    @test white_runtime.accepted_state.dielectric_loss_energy_j >= 0.0
    white_result = TransformerPublic.transformer_apparatus_result(white_runtime)
    @test TransformerPublic.transformer_result_quantity_available(
        white_result.winding_section_voltage_v,
    )
    @test TransformerPublic.transformer_result_quantity_available(
        white_result.winding_section_current_a,
    )
    @test !TransformerPublic.transformer_result_quantity_available(
        white_result.ladder_node_voltage_v,
    )
end

test_model_consistent_transformer_initialization()

@testset "transformer apparatus localized topology and fault events" begin
    source = synthetic_transformer_source()
    connection = grounded_transformer_connection()
    matrices = TransformerPublic.TransformerTerminalMatrices(
        [0.4 0.0; 0.0 0.5],
        [2.0e-3 0.0; 0.0 3.0e-3];
        capacitance_f=[2.0e-7 0.0; 0.0 1.0e-7],
        conductance_s=[1.0e-5 0.0; 0.0 2.0e-5],
    )
    specification = transformer_specification(
        TransformerPublic.LowFrequencyTerminalTier,
        connection,
        TransformerPublic.LowFrequencyTransformerModel(matrices),
        source,
    )
    runtime = TransformerPublic.transformer_apparatus_runtime(
        TransformerPublic.prepare_transformer_apparatus(specification),
        [1, 2],
    )
    open_command = TransformerPublic.TransformerApparatusEventCommand(
        :primary_breaker_opens,
        TransformerPublic.TransformerBreakerOpenEvent,
        10.0e-6;
        target_indices=[1, 2],
    )
    TransformerPublic.queue_transformer_apparatus_event!(runtime, open_command)
    surfaces = nonlinear_device_event_surfaces(runtime)
    @test length(surfaces) == 1
    @test nonlinear_device_event_candidate_time(only(surfaces), runtime) == 10.0e-6
    @test nonlinear_device_event_value(only(surfaces), runtime, 0.0) === nothing
    accepted_current, _ = accept_transformer_test_step!(runtime, [10.0, -4.0])
    @test maximum(abs, accepted_current) > 0.0
    apply_nonlinear_device_event!(only(surfaces), runtime, 10.0e-6)
    before_open_state = copy(runtime.accepted_state.coil_current_a)
    open_current, open_jacobian = accept_transformer_test_step!(
        runtime,
        [50.0, -20.0],
        20.0e-6,
        companion_method=:backward_euler,
    )
    @test open_current == zeros(2)
    @test open_jacobian == zeros(2, 2)
    @test runtime.accepted_state.coil_current_a != before_open_state
    @test runtime.accepted_state.accepted_time_s == 20.0e-6

    close_command = TransformerPublic.TransformerApparatusEventCommand(
        :primary_breaker_recloses,
        TransformerPublic.TransformerBreakerCloseEvent,
        30.0e-6;
        target_indices=[1, 2],
    )
    TransformerPublic.queue_transformer_apparatus_event!(runtime, close_command)
    close_surface = only(nonlinear_device_event_surfaces(runtime))
    accept_transformer_test_step!(runtime, [0.0, 0.0], 30.0e-6)
    apply_nonlinear_device_event!(close_surface, runtime, 30.0e-6)
    reclosed_current, _ = accept_transformer_test_step!(
        runtime,
        [10.0, -4.0],
        40.0e-6,
        companion_method=:backward_euler,
    )
    @test maximum(abs, reclosed_current) > 0.0

    fault_command = TransformerPublic.TransformerApparatusEventCommand(
        :primary_terminal_ground_fault,
        TransformerPublic.TransformerTerminalFaultApplyEvent,
        50.0e-6;
        target_indices=[1],
        conductance_s=2.0,
    )
    clear_command = TransformerPublic.TransformerApparatusEventCommand(
        :primary_terminal_ground_fault_clears,
        TransformerPublic.TransformerTerminalFaultClearEvent,
        60.0e-6;
        reference_id=:primary_terminal_ground_fault,
    )
    TransformerPublic.queue_transformer_apparatus_event!(runtime, fault_command)
    TransformerPublic.queue_transformer_apparatus_event!(runtime, clear_command)
    event_surfaces = nonlinear_device_event_surfaces(runtime)
    fault_surface = only(filter(row -> row.name == :primary_terminal_ground_fault, event_surfaces))
    accept_transformer_test_step!(runtime, [10.0, -4.0], 50.0e-6)
    apply_nonlinear_device_event!(fault_surface, runtime, 50.0e-6)
    fault_current, fault_jacobian = accept_transformer_test_step!(
        runtime,
        [10.0, -4.0],
        60.0e-6,
        companion_method=:backward_euler,
    )
    @test fault_jacobian[1, 1] >= 2.0
    @test fault_current[1] >= 20.0
    event_state = TransformerPublic.transformer_apparatus_event_state(runtime)
    @test event_state.event_energy_j > 0.0
    @test event_state.external_fault_energy_j == event_state.event_energy_j
    @test event_state.internal_fault_energy_j == 0.0
    @test event_state.numerical_dissipation_energy_j >= 0.0
    @test event_state.maximum_energy_balance_residual_j <=
        specification.settings.energy_absolute_tolerance_j + 1024.0 * eps(Float64)
    @test event_state.topology_transition_count == 3
    @test getfield.(
        TransformerPublic.transformer_apparatus_event_occurrences(runtime),
        :id,
    ) == [
        :primary_breaker_opens,
        :primary_breaker_recloses,
        :primary_terminal_ground_fault,
    ]

    snapshot = TransformerPublic.transformer_apparatus_runtime_snapshot(runtime)
    restored = TransformerPublic.transformer_apparatus_runtime(
        TransformerPublic.prepare_transformer_apparatus(specification),
        [1, 2],
    )
    TransformerPublic.restore_transformer_apparatus_runtime_snapshot!(restored, snapshot)
    @test TransformerPublic.transformer_apparatus_runtime_snapshot(restored).
        deterministic_signature_sha256 == snapshot.deterministic_signature_sha256
    @test TransformerPublic.transformer_apparatus_event_state(restored) == event_state

    clear_surface = only(filter(
        row -> row.name == :primary_terminal_ground_fault_clears,
        nonlinear_device_event_surfaces(runtime),
    ))
    apply_nonlinear_device_event!(clear_surface, runtime, 60.0e-6)
    cleared_current, cleared_jacobian = accept_transformer_test_step!(
        runtime,
        [10.0, -4.0],
        70.0e-6,
        companion_method=:backward_euler,
    )
    @test cleared_jacobian[1, 1] < fault_jacobian[1, 1] - 1.0
    @test cleared_current[1] < fault_current[1]

    @test_throws TransformerPublic.TransformerApparatusRefusal begin
        TransformerPublic.queue_transformer_apparatus_event!(
            runtime,
            TransformerPublic.TransformerApparatusEventCommand(
                :unrepresented_internal_fault,
                TransformerPublic.TransformerInternalFaultApplyEvent,
                80.0e-6;
                target_indices=[1, 2],
                conductance_s=1.0,
            ),
        )
    end
    @test_throws ArgumentError TransformerPublic.TransformerApparatusEventCommand(
        :invalid_clear,
        TransformerPublic.TransformerTerminalFaultClearEvent,
        80.0e-6,
    )

    backward_euler_runtime = TransformerPublic.transformer_apparatus_runtime(
        TransformerPublic.prepare_transformer_apparatus(specification),
        [1, 2],
    )
    backward_euler_current, backward_euler_jacobian =
        accept_transformer_test_step!(
            backward_euler_runtime,
            [10.0, -4.0],
            5.0e-6;
            step_s=5.0e-6,
            companion_method=:backward_euler,
        )
    @test all(isfinite, backward_euler_current)
    @test backward_euler_jacobian ≈ transpose(backward_euler_jacobian) atol=1.0e-12
    @test backward_euler_runtime.accepted_state.accepted_time_s == 5.0e-6

    represented_model = TransformerPublic.GreyBoxTransformerModel(
        node_order=[:primary_terminal, :secondary_terminal, :internal_node],
        terminal_node_indices=[1, 2],
        branches=[
            TransformerPublic.TransformerLadderBranch(
                :primary_section,
                1,
                3;
                resistance_ohm=0.2,
                inductance_h=2.0e-3,
            ),
            TransformerPublic.TransformerLadderBranch(
                :secondary_section,
                3,
                2;
                resistance_ohm=0.3,
                inductance_h=3.0e-3,
            ),
        ],
        capacitance_f=Diagonal([1.0e-9, 1.0e-9, 2.0e-9]),
        conductance_s=Diagonal([1.0e-7, 1.0e-7, 2.0e-7]),
        source_response_sha256=repeat("8", 64),
        identification_residual_relative=1.0e-4,
        parameter_nonuniqueness="synthetic represented-event regression topology",
    )
    represented_specification = transformer_specification(
        TransformerPublic.GreyBoxLadderTier,
        connection,
        represented_model,
        source,
    )
    represented_runtime = TransformerPublic.transformer_apparatus_runtime(
        TransformerPublic.prepare_transformer_apparatus(represented_specification),
        [1, 2],
    )
    represented_fault = TransformerPublic.TransformerApparatusEventCommand(
        :represented_terminal_fault,
        TransformerPublic.TransformerTerminalFaultApplyEvent,
        10.0e-6;
        target_indices=[1],
        conductance_s=0.01,
    )
    TransformerPublic.queue_transformer_apparatus_event!(
        represented_runtime,
        represented_fault,
    )
    represented_surface = only(nonlinear_device_event_surfaces(represented_runtime))
    accept_transformer_test_step!(
        represented_runtime,
        [10.0, -4.0],
        10.0e-6;
        companion_method=:backward_euler,
    )
    apply_nonlinear_device_event!(
        represented_surface,
        represented_runtime,
        10.0e-6,
    )
    post_event_current, post_event_jacobian = accept_transformer_test_step!(
        represented_runtime,
        [8.0, -3.0],
        20.0e-6;
        companion_method=:trapezoidal,
    )
    @test all(isfinite, post_event_current)
    @test post_event_jacobian ≈ transpose(post_event_jacobian) atol=1.0e-12
    represented_event_state =
        TransformerPublic.transformer_apparatus_event_state(represented_runtime)
    @test represented_event_state.external_fault_energy_j > 0.0
    @test represented_event_state.maximum_energy_balance_residual_j <=
        represented_specification.settings.energy_absolute_tolerance_j +
        1024.0 * eps(Float64)
    @test abs(
        TransformerPublic.transformer_apparatus_result(represented_runtime).energy.
            unexplained_balance_residual_j,
    ) <= represented_specification.settings.energy_absolute_tolerance_j +
        1024.0 * eps(Float64)
end

@testset "represented winding fault and three-phase shift events" begin
    source = synthetic_transformer_source()
    white_connection = white_box_transformer_connection()
    white_model = TransformerPublic.WhiteBoxTransformerModel(
        winding_order=[:primary_winding, :secondary_winding],
        section_count_per_winding=[2, 2],
        section_length_m=[0.5, 0.5],
        series_resistance_ohm_per_m=fill([0.2 0.02; 0.02 0.3], 2),
        series_inductance_h_per_m=fill([2.0e-3 0.2e-3; 0.2e-3 3.0e-3], 2),
        shunt_conductance_s_per_m=
            fill(1.0e-7 .* Matrix{Float64}(I, 2, 2), 2),
        shunt_capacitance_f_per_m=
            fill(1.0e-9 .* Matrix{Float64}(I, 2, 2), 2),
        geometry_sha256=repeat("5", 64),
        frequency_band_hz=(10.0, 10_000.0),
        section_refinement_residual_relative=1.0e-3,
    )
    white_specification = transformer_specification(
        TransformerPublic.WhiteBoxWindingTier,
        white_connection,
        white_model,
        source,
    )
    white_runtime = TransformerPublic.transformer_apparatus_runtime(
        TransformerPublic.prepare_transformer_apparatus(white_specification),
        [1, 2, 3, 4],
    )
    baseline_runtime = deepcopy(white_runtime)
    internal_fault = TransformerPublic.TransformerApparatusEventCommand(
        :represented_section_fault,
        TransformerPublic.TransformerInternalFaultApplyEvent,
        10.0e-6;
        target_indices=[5, 6],
        conductance_s=5.0,
    )
    TransformerPublic.queue_transformer_apparatus_event!(white_runtime, internal_fault)
    internal_surface = only(nonlinear_device_event_surfaces(white_runtime))
    first_voltage = [10.0, 0.0, 5.0, 0.0]
    accept_transformer_test_step!(white_runtime, first_voltage)
    accept_transformer_test_step!(baseline_runtime, first_voltage)
    apply_nonlinear_device_event!(internal_surface, white_runtime, 10.0e-6)
    second_voltage = [8.0, 0.0, 4.0, 0.0]
    fault_current, fault_jacobian = accept_transformer_test_step!(
        white_runtime,
        second_voltage,
        20.0e-6,
        companion_method=:backward_euler,
    )
    baseline_current, baseline_jacobian = accept_transformer_test_step!(
        baseline_runtime,
        second_voltage,
        20.0e-6,
        companion_method=:backward_euler,
    )
    @test fault_current != baseline_current
    @test fault_jacobian != baseline_jacobian
    @test fault_jacobian ≈ transpose(fault_jacobian) atol=2.0e-12
    @test minimum(eigvals(Symmetric(fault_jacobian))) >= -2.0e-12
    @test white_runtime.accepted_state.maximum_internal_kcl_residual_a <= 2.0e-9
    @test haskey(
        TransformerPublic.transformer_apparatus_event_state(white_runtime).
            active_internal_faults,
        :represented_section_fault,
    )
    internal_event_state =
        TransformerPublic.transformer_apparatus_event_state(white_runtime)
    @test internal_event_state.internal_fault_energy_j > 0.0
    @test internal_event_state.external_fault_energy_j == 0.0
    @test internal_event_state.event_energy_j ==
        internal_event_state.internal_fault_energy_j
    @test internal_event_state.numerical_dissipation_energy_j >= 0.0
    @test internal_event_state.maximum_energy_balance_residual_j <=
        white_specification.settings.energy_absolute_tolerance_j +
        1024.0 * eps(Float64)

    phases = [:phase_a, :phase_b, :phase_c]
    three_phase_connection = TransformerPublic.TransformerConnectionTopology(
        node_order=[:terminal_a, :terminal_b, :terminal_c],
        node_phase=phases,
        coil_order=[:coil_a, :coil_b, :coil_c],
        winding_order=[:winding],
        phase_order=phases,
        coil_winding=fill(:winding, 3),
        coil_phase=phases,
        incidence=Matrix{Float64}(I, 3, 3),
        vector_group="Y",
    )
    three_phase_model = TransformerPublic.LowFrequencyTransformerModel(
        TransformerPublic.TransformerTerminalMatrices(
            Diagonal([0.4, 0.4, 0.4]),
            [0.12 0.01 0.01; 0.01 0.12 0.01; 0.01 0.01 0.12];
            capacitance_f=1.0e-9 .* Matrix{Float64}(I, 3, 3),
            conductance_s=1.0e-7 .* Matrix{Float64}(I, 3, 3),
        ),
    )
    three_phase_specification = TransformerPublic.TransformerApparatusSpecification(
        :three_phase_shift_transformer,
        TransformerPublic.LowFrequencyTerminalTier,
        three_phase_connection,
        three_phase_model,
        TransformerPublic.TransformerRuntimeSettings(
            timestep_s=10.0e-6,
            initialization_frequency_hz=60.0,
        );
        phase_count=3,
        rated_power_va=1.0e6,
        rated_voltage_v=13.8e3,
        rated_frequency_hz=60.0,
        sources=[source],
        uncertainty="exact synthetic phase-shift event values",
        validity_domain="balanced and unbalanced synthetic three-phase transform",
    )
    phase_runtime = TransformerPublic.transformer_apparatus_runtime(
        TransformerPublic.prepare_transformer_apparatus(three_phase_specification),
        [1, 2, 3],
    )
    phase_angle_rad = pi / 12.0
    phase_transform = [
        1.0 / 3.0 + 2.0 / 3.0 * cos(
            phase_angle_rad + 2.0 * pi * (row - column) / 3.0,
        ) for row in 1:3, column in 1:3
    ]
    @test phase_transform * transpose(phase_transform) ≈
        Matrix{Float64}(I, 3, 3) atol=2.0e-15
    phase_command = TransformerPublic.TransformerApparatusEventCommand(
        :continuous_phase_shift_changes,
        TransformerPublic.TransformerPhaseShiftChangeEvent,
        10.0e-6;
        terminal_transform=phase_transform,
    )
    TransformerPublic.queue_transformer_apparatus_event!(phase_runtime, phase_command)
    phase_surface = only(nonlinear_device_event_surfaces(phase_runtime))
    accept_transformer_test_step!(phase_runtime, [10.0, -5.0, -5.0])
    apply_nonlinear_device_event!(phase_surface, phase_runtime, 10.0e-6)
    network_voltage = [8.0, -2.0, -6.0]
    network_current, network_jacobian = accept_transformer_test_step!(
        phase_runtime,
        network_voltage,
        20.0e-6,
        companion_method=:backward_euler,
    )
    @test network_jacobian ≈ transpose(network_jacobian) atol=2.0e-12
    @test phase_runtime.accepted_state.terminal_voltage_v ≈
        phase_transform * network_voltage atol=2.0e-13
    @test dot(network_voltage, network_current) ≈ dot(
        phase_runtime.accepted_state.terminal_voltage_v,
        phase_runtime.accepted_state.terminal_current_a,
    ) atol=2.0e-12
    @test TransformerPublic.transformer_apparatus_event_state(phase_runtime).
        phase_shift_change_count == 1
end

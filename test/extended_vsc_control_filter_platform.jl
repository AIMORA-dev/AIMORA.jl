using Test, AIMORA

const ExtendedVSC = AIMORA.SwitchDetailedVSC

@testset "VSC balance uses the actual accepted companion quadrature" begin
    frame(voltage, current) = (
        power=(;external_power_w=current, dissipated_power_w=current^2,
            bridge=(;companion_energy_residual_j=0.0)),
        stored_j=0.5voltage^2,
        linear=(;capacitor_power_w=voltage * current,
            capacitor_stored_energy_j=0.5voltage^2, inductor_power_w=0.0,
            inductor_dissipated_power_w=0.0, inductor_stored_energy_j=0.0),
    )
    increment(left, right, step, method) = AIMORA.EMTStudy._extended_vsc_balance_increment(
        left.power, right.power, left.stored_j, right.stored_j,
        left.linear, right.linear, step, method)
    # Analytic unit R-C apparatus driven by a constant unit source: the
    # resistor's physical loss and the method's storage defect stay distinct.
    for step in (0.1, 0.01), method in
        (AIMORA.Branches.TrapezoidalCompanion, AIMORA.Branches.BackwardEulerCompanion)
        implicit = method === AIMORA.Branches.BackwardEulerCompanion
        left = frame(0.0, 1.0)
        voltage = implicit ? step / (1 + step) : step / (1 + step / 2)
        right = frame(voltage, 1 - voltage)
        balance = increment(left, right, step, method)
        @test abs(balance.residual_j) <= 4eps(Float64)
        expected_external = step * (implicit ? right.power.external_power_w :
            0.5 * (left.power.external_power_w + right.power.external_power_w))
        expected_loss = step * (implicit ? right.power.dissipated_power_w :
            0.5 * (left.power.dissipated_power_w + right.power.dissipated_power_w))
        @test balance.external_work_j == expected_external
        @test balance.dissipated_work_j == expected_loss
        expected_defect = implicit ? 0.5voltage^2 :
            step * voltage * (right.power.external_power_w - 1.0) / 4
        @test balance.linear_defect_j ≈ expected_defect atol=4eps(Float64) rtol=0
        if implicit
            @test balance.linear_defect_j > 0
            @test balance.external_work_j - balance.dissipated_work_j - right.stored_j > 0
        end
    end
    left = frame(0.0, 1.0)
    middle = frame(0.05 / 1.05, 1 / 1.05)
    final_voltage = (middle.stored_j == 0 ? 0.0 : sqrt(2middle.stored_j))
    final_voltage = (final_voltage + 0.05) / 1.05
    right = frame(final_voltage, 1 - final_voltage)
    first_half = increment(left, middle, 0.05, AIMORA.Branches.BackwardEulerCompanion)
    second_half = increment(middle, right, 0.05, AIMORA.Branches.BackwardEulerCompanion)
    @test abs(first_half.residual_j + second_half.residual_j) <= 4eps(Float64)
    @test_throws ArgumentError increment(left, right, 0.1, :undeclared_companion)
end

function sequence_phase_sample(
    angle_rad;
    positive_magnitude_v=0.0,
    positive_angle_rad=0.0,
    negative_magnitude_v=0.0,
    negative_angle_rad=0.0,
    zero_magnitude_v=0.0,
    zero_angle_rad=0.0,
)
    positive_angle = angle_rad + positive_angle_rad
    negative_angle = angle_rad + negative_angle_rad
    zero = zero_magnitude_v * cos(angle_rad + zero_angle_rad)
    return (
        positive_magnitude_v * cos(positive_angle) +
            negative_magnitude_v * cos(negative_angle) + zero,
        positive_magnitude_v * cos(positive_angle - 2.0 * pi / 3.0) +
            negative_magnitude_v * cos(negative_angle + 2.0 * pi / 3.0) + zero,
        positive_magnitude_v * cos(positive_angle + 2.0 * pi / 3.0) +
            negative_magnitude_v * cos(negative_angle - 2.0 * pi / 3.0) + zero,
    )
end

@testset "extended VSC causal rotating sequence extraction" begin
    parameters = ExtendedVSC.ExtendedVSCParameters()
    sequence_state = ExtendedVSC.VSCRotatingSequenceState(parameters)
    sample_count = length(sequence_state.samples)
    omega = 2.0 * pi * parameters.frequency_hz
    extracted = nothing
    for sample_index in 0:(sample_count - 1)
        angle = omega * sample_index * parameters.controller.control_period_s
        phase_voltage = sequence_phase_sample(
            angle;
            positive_magnitude_v=325.0,
            positive_angle_rad=0.23,
            negative_magnitude_v=32.5,
            negative_angle_rad=-0.41,
            zero_magnitude_v=16.25,
            zero_angle_rad=0.67,
        )
        extracted = ExtendedVSC.advance_vsc_sequence_extractor!(
            sequence_state,
            phase_voltage,
            angle,
        )
    end
    @test !isnothing(extracted)
    @test extracted.settled
    @test extracted.positive_magnitude_v ≈ 325.0 atol=2.0e-12
    @test extracted.negative_magnitude_v ≈ 32.5 atol=2.0e-12
    @test extracted.zero_magnitude_v ≈ 16.25 atol=2.0e-12
    @test atan(extracted.positive.quadrature, extracted.positive.direct) ≈
        0.23 atol=2.0e-14
    @test atan(extracted.negative.quadrature, extracted.negative.direct) ≈
        0.41 atol=2.0e-14

    balanced_state = ExtendedVSC.VSCRotatingSequenceState(parameters)
    for sample_index in 0:(sample_count - 1)
        angle = omega * sample_index * parameters.controller.control_period_s
        extracted = ExtendedVSC.advance_vsc_sequence_extractor!(
            balanced_state,
            sequence_phase_sample(angle; positive_magnitude_v=325.0),
            angle,
        )
    end
    @test extracted.positive_magnitude_v ≈ 325.0 atol=2.0e-12
    @test extracted.negative_magnitude_v <= 2.0e-12
    @test extracted.zero_magnitude_v <= 2.0e-12

    reversed_state = ExtendedVSC.VSCRotatingSequenceState(parameters)
    for sample_index in 0:(sample_count - 1)
        angle = omega * sample_index * parameters.controller.control_period_s
        phase_voltage = sequence_phase_sample(angle; positive_magnitude_v=325.0)
        extracted = ExtendedVSC.advance_vsc_sequence_extractor!(
            reversed_state,
            (phase_voltage[1], phase_voltage[3], phase_voltage[2]),
            angle,
        )
    end
    @test extracted.positive_magnitude_v <= 2.0e-12
    @test extracted.negative_magnitude_v ≈ 325.0 atol=2.0e-12

    checkpoint_state = ExtendedVSC.VSCRotatingSequenceState(parameters)
    for sample_index in 0:(sample_count ÷ 2 - 1)
        angle = omega * sample_index * parameters.controller.control_period_s
        ExtendedVSC.advance_vsc_sequence_extractor!(
            checkpoint_state,
            sequence_phase_sample(
                angle;
                positive_magnitude_v=325.0,
                negative_magnitude_v=16.25,
            ),
            angle,
        )
    end
    restored_state = deepcopy(checkpoint_state)
    for sample_index in (sample_count ÷ 2):(sample_count + 10)
        angle = omega * sample_index * parameters.controller.control_period_s
        phase_voltage = sequence_phase_sample(
            angle;
            positive_magnitude_v=325.0,
            negative_magnitude_v=16.25,
        )
        original_result = ExtendedVSC.advance_vsc_sequence_extractor!(
            checkpoint_state,
            phase_voltage,
            angle,
        )
        restored_result = ExtendedVSC.advance_vsc_sequence_extractor!(
            restored_state,
            phase_voltage,
            angle,
        )
        @test original_result == restored_result
    end
    @test checkpoint_state.samples == restored_state.samples
    @test checkpoint_state.next_sample_index == restored_state.next_sample_index
end

@testset "extended VSC public control and refusal contracts" begin
    @test ExtendedVSC.extended_vsc_contract().id ==
        :extended_switch_detailed_vsc_platform

    active_priority = ExtendedVSC.project_vsc_current_reference(
        200.0,
        100.0,
        180.0,
        ExtendedVSC.ActiveCurrentPriority,
    )
    @test active_priority == (direct=180.0, quadrature=0.0, limited=true)
    reactive_priority = ExtendedVSC.project_vsc_current_reference(
        200.0,
        100.0,
        180.0,
        ExtendedVSC.ReactiveCurrentPriority,
    )
    @test reactive_priority.direct ≈ sqrt(180.0^2 - 100.0^2)
    @test reactive_priority.quadrature == 100.0
    vector_priority = ExtendedVSC.project_vsc_current_reference(
        200.0,
        100.0,
        180.0,
        ExtendedVSC.VectorMagnitudePriority,
    )
    @test hypot(vector_priority.direct, vector_priority.quadrature) ≈ 180.0
    @test vector_priority.direct / vector_priority.quadrature ≈ 2.0

    three_wire_duties = ExtendedVSC.extended_vsc_modulation_duties(
        (240.0, -120.0, -120.0),
        70.0,
        800.0,
        ExtendedVSC.ThreeWireForm,
    )
    @test collect(three_wire_duties) ≈ [0.725, 0.275, 0.275, 0.5]
    four_wire_duties = ExtendedVSC.extended_vsc_modulation_duties(
        (240.0, -120.0, -120.0),
        70.0,
        800.0,
        ExtendedVSC.FourWireForm,
    )
    @test all(0.02 .<= collect(four_wire_duties) .<= 0.98)
    @test four_wire_duties[1] - four_wire_duties[4] ≈ (240.0 - 70.0) / 800.0

    @test_throws ArgumentError ExtendedVSC.validate_extended_vsc_parameters(
        ExtendedVSC.ExtendedVSCParameters(
            scenario=ExtendedVSC.ExtendedVSCScenarioParameters(
                zero_sequence_voltage_ratio=0.05,
            ),
        ),
    )
    @test_throws ArgumentError ExtendedVSC.validate_extended_vsc_parameters(
        ExtendedVSC.ExtendedVSCParameters(
            controller=ExtendedVSC.ExtendedVSCControlParameters(
                selected_harmonic_orders=(1, 5, 5),
            ),
        ),
    )
    @test_throws ArgumentError ExtendedVSC.validate_extended_vsc_parameters(
        ExtendedVSC.ExtendedVSCParameters(
            controller=ExtendedVSC.ExtendedVSCControlParameters(
                virtual_inertia_w_s2_per_rad=0.0,
            ),
        ),
    )

    stale_parameters = ExtendedVSC.ExtendedVSCParameters()
    stale_state = ExtendedVSC.ExtendedVSCControlState(stale_parameters)
    stale_measurement = ExtendedVSC.ExtendedVSCMeasurement(
        (325.0, -162.5, -162.5),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (325.0, -162.5, -162.5),
        0.0,
        800.0,
        2.0e-3,
    )
    stale_request = ExtendedVSC.ExtendedVSCPlantRequest(
        timestamp_s=0.0,
        valid_until_s=1.0e-3,
    )
    stale_command = ExtendedVSC.compute_extended_vsc_control!(
        stale_state,
        stale_measurement,
        stale_request,
        stale_parameters,
    )
    @test stale_command.mode === ExtendedVSC.VSCBlockedOperation
    @test stale_command.request_disposition === ExtendedVSC.VSCPlantRequestStale
    @test stale_command.duties == (0.5, 0.5, 0.5, 0.5)
    @test stale_state.refusal_count == 1
    @test stale_state.sequence.retained_sample_count == 0
end

@testset "extended VSC controller-family equations" begin
    request = ExtendedVSC.ExtendedVSCPlantRequest(
        active_power_reference_w=0.0,
        reactive_power_reference_var=0.0,
    )
    for family in (
        ExtendedVSC.SynchronousPLLGridFollowing,
        ExtendedVSC.StationaryResonantGridFollowing,
        ExtendedVSC.PowerDroopGridForming,
        ExtendedVSC.VirtualSynchronousGridForming,
    )
        parameters = ExtendedVSC.ExtendedVSCParameters(
            controller=ExtendedVSC.ExtendedVSCControlParameters(family=family),
        )
        state = ExtendedVSC.ExtendedVSCControlState(parameters)
        sample_count = length(state.sequence.samples)
        omega = 2.0 * pi * parameters.frequency_hz
        command = nothing
        for sample_index in 0:(sample_count - 1)
            angle = state.angle_rad
            phase_voltage = sequence_phase_sample(
                angle;
                positive_magnitude_v=325.0,
                negative_magnitude_v=16.25,
            )
            measurement = ExtendedVSC.ExtendedVSCMeasurement(
                phase_voltage,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                phase_voltage,
                0.0,
                800.0,
                sample_index * parameters.controller.control_period_s,
            )
            command = ExtendedVSC.compute_extended_vsc_control!(
                state,
                measurement,
                request,
                parameters,
            )
        end
        @test !isnothing(command)
        @test command.controller_family === family
        @test command.sequence_extractor_settled
        @test command.positive_sequence_voltage_v ≈ 325.0 atol=2.0e-10
        @test command.negative_sequence_voltage_v ≈ 16.25 atol=2.0e-10
        @test command.zero_sequence_voltage_v <= 2.0e-12
        @test all(isfinite, command.duties)
        @test all(parameters.minimum_duty .<= collect(command.duties) .<=
            parameters.maximum_duty)
        @test state.sample_count == sample_count
        if family in (
            ExtendedVSC.SynchronousPLLGridFollowing,
            ExtendedVSC.StationaryResonantGridFollowing,
        )
            @test command.frequency_hz ≈ parameters.frequency_hz atol=1.0e-11
            @test state.pll_locked
        else
            @test !state.pll_locked
        end
    end

    droop_control = ExtendedVSC.ExtendedVSCControlParameters(
        family=ExtendedVSC.PowerDroopGridForming,
    )
    droop_parameters = ExtendedVSC.ExtendedVSCParameters(controller=droop_control)
    droop_state = ExtendedVSC.ExtendedVSCControlState(droop_parameters)
    droop_request = ExtendedVSC.ExtendedVSCPlantRequest(
        active_power_reference_w=10.0e3,
        reactive_power_reference_var=0.0,
    )
    droop_measurement = ExtendedVSC.ExtendedVSCMeasurement(
        (325.0, -162.5, -162.5),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (325.0, -162.5, -162.5),
        0.0,
        800.0,
        0.0,
    )
    droop_command = ExtendedVSC.compute_extended_vsc_control!(
        droop_state,
        droop_measurement,
        droop_request,
        droop_parameters,
    )
    expected_droop_frequency = droop_request.frequency_reference_hz +
        droop_control.active_power_frequency_droop_rad_per_ws *
        droop_request.active_power_reference_w / (2.0 * pi)
    @test droop_command.frequency_hz ≈ expected_droop_frequency atol=1.0e-12

    swing_control = ExtendedVSC.ExtendedVSCControlParameters(
        family=ExtendedVSC.VirtualSynchronousGridForming,
    )
    swing_parameters = ExtendedVSC.ExtendedVSCParameters(controller=swing_control)
    swing_state = ExtendedVSC.ExtendedVSCControlState(swing_parameters)
    swing_command = ExtendedVSC.compute_extended_vsc_control!(
        swing_state,
        droop_measurement,
        droop_request,
        swing_parameters,
    )
    expected_swing_frequency = droop_request.frequency_reference_hz +
        swing_control.control_period_s * droop_request.active_power_reference_w /
        swing_control.virtual_inertia_w_s2_per_rad / (2.0 * pi)
    @test swing_command.frequency_hz ≈ expected_swing_frequency atol=1.0e-12

    three_wire_state = ExtendedVSC.ExtendedVSCControlState(
        ExtendedVSC.ExtendedVSCParameters(),
    )
    three_wire_state.zero_current_integral_v = 3.25
    three_wire_measurement = ExtendedVSC.ExtendedVSCMeasurement(
        droop_measurement.phase_voltage_v,
        droop_measurement.converter_current_a,
        (1.0, 1.0, 1.0),
        droop_measurement.filter_voltage_v,
        droop_measurement.neutral_current_a,
        droop_measurement.dc_link_voltage_v,
        droop_measurement.time_s,
    )
    ExtendedVSC.compute_extended_vsc_control!(
        three_wire_state,
        three_wire_measurement,
        request,
        ExtendedVSC.ExtendedVSCParameters(),
    )
    @test three_wire_state.zero_current_integral_v == 3.25

    pr_parameters = ExtendedVSC.ExtendedVSCParameters(
        controller=ExtendedVSC.ExtendedVSCControlParameters(
            family=ExtendedVSC.StationaryResonantGridFollowing,
            selected_harmonic_orders=(1,),
        ),
    )
    pr_state = ExtendedVSC.ExtendedVSCControlState(pr_parameters)
    pr_request = ExtendedVSC.ExtendedVSCPlantRequest(active_power_reference_w=20.0e3)
    ExtendedVSC.compute_extended_vsc_control!(
        pr_state,
        droop_measurement,
        pr_request,
        pr_parameters,
    )
    @test pr_state.resonant_alpha[1].derivative_state != 0.0
    @test all(
        resonator -> resonator.derivative_state == 0.0 && resonator.integral_state == 0.0,
        pr_state.resonant_alpha[2:end],
    )

    four_wire_parameters = ExtendedVSC.ExtendedVSCParameters(
        wire_form=ExtendedVSC.FourWireForm,
        scenario=ExtendedVSC.ExtendedVSCScenarioParameters(
            zero_sequence_voltage_ratio=0.05,
        ),
    )
    four_wire_state = ExtendedVSC.ExtendedVSCControlState(four_wire_parameters)
    four_wire_state.zero_current_integral_v = 3.25
    four_wire_voltage = sequence_phase_sample(
        0.0;
        positive_magnitude_v=325.0,
        zero_magnitude_v=16.25,
    )
    four_wire_command = ExtendedVSC.compute_extended_vsc_control!(
        four_wire_state,
        ExtendedVSC.ExtendedVSCMeasurement(
            four_wire_voltage,
            (1.0, 1.0, 1.0),
            (0.0, 0.0, 0.0),
            four_wire_voltage,
            0.0,
            800.0,
            0.0,
        ),
        request,
        four_wire_parameters,
    )
    @test abs(sum(four_wire_command.phase_current_reference_a)) > 0.0
    @test four_wire_command.duties[4] != 0.5
    @test four_wire_state.zero_current_integral_v != 3.25
end
@testset "declared VSC source endpoints preserve boundary ownership" begin
    emt = AIMORA.EMTStudy
    endpoint = 0.7e-3
    plan = [(:commanded_block, endpoint, true),
        (:island_begins, endpoint, true), (:grid_reconnects, 0.8e-3, true)]
    runtime = (; boundary_plan=plan, next_boundary_index=1)
    @test !emt._extended_vsc_source_boundary_at(runtime, prevfloat(endpoint))
    @test emt._extended_vsc_source_boundary_at(runtime, endpoint)
    @test !emt._extended_vsc_source_boundary_at(runtime, nextfloat(endpoint))
    @test emt._extended_vsc_source_boundary_at(runtime, 0.8e-3)
    @test runtime.next_boundary_index == 1
    @test runtime.boundary_plan == plan
    consumed = (; boundary_plan=plan, next_boundary_index=3)
    @test !emt._extended_vsc_source_boundary_at(consumed, endpoint)
    @test emt._extended_vsc_source_boundary_at(consumed, 0.8e-3)
    @test !emt._extended_vsc_is_source_boundary(:commanded_block)
    @test !emt._extended_vsc_is_source_boundary(:commanded_restart)
    for name in (:fault_begins, :fault_clears, :dc_sag_begins,
                 :dc_sag_clears, :island_begins, :grid_reconnects)
        @test emt._extended_vsc_is_source_boundary(name)
    end
end

@testset "VSC accounting preserves compensated cancellation and numeric transitions" begin
    E=AIMORA.EMTStudy;B=AIMORA.Branches
    power(x)=(;external_power_w=x,dissipated_power_w=zero(x),bridge=(;companion_energy_residual_j=0.0))
    linear(x)=(;capacitor_power_w=zero(x),capacitor_stored_energy_j=zero(x),
        inductor_power_w=x,inductor_dissipated_power_w=x,inductor_stored_energy_j=zero(x))
    for bits in (106,256), method in (B.TrapezoidalCompanion,B.BackwardEulerCompanion)
        setprecision(BigFloat,bits) do
            delta=BigFloat(2)^(-70);initial=BigFloat(1);final=initial+delta
            @test Float64(final)-Float64(initial)==0.0
            account=E.ExtendedVSCBalanceAccounting(power(delta),initial,linear(BigFloat(0)),0.0,0.0)
            setprecision(BigFloat,53) do
                E._accumulate_extended_vsc_balance!(account,power(delta),final,linear(BigFloat(0)),1.0,method)
                @test precision(BigFloat)==53
            end
            @test account.energy_residual_j==0
            @test account.previous_stored_energy_j==final
            @test precision(account.previous_stored_energy_j)==bits
            @test isconcretetype(typeof(account.previous_power))
            wrong=E.ExtendedVSCBalanceAccounting(power(2delta),initial,linear(BigFloat(0)),0.0,0.0)
            E._accumulate_extended_vsc_balance!(wrong,power(2delta),final,linear(BigFloat(0)),1.0,method)
            @test wrong.energy_residual_j==delta
            report=E._extended_vsc_reported_energy_ratio(wrong.energy_residual_j,eps(Float64))
            @test report.exact_ratio==delta/BigFloat(eps(Float64))
            @test report.exact_ratio>1e-6 && report.reported>1e-6
            @test isfinite(report.reported)
            @test BigFloat(prevfloat(report.reported))<=report.exact_ratio<=BigFloat(nextfloat(report.reported))
            # A Float64 baseline must admit a compensated accepted substep
            # without freezing either power/linear record to its initial type.
            mixed=E.ExtendedVSCBalanceAccounting(power(0.0),0.0,linear(0.0),0.0,0.0)
            factor=method===B.TrapezoidalCompanion ? 2 : 1
            compensated_power=power(factor*final)
            compensated_linear=linear(final)
            setprecision(BigFloat,53) do
                E._accumulate_extended_vsc_balance!(mixed,compensated_power,final,compensated_linear,1.0,method)
            end
            @test mixed.energy_residual_j==0
            @test mixed.previous_power.external_power_w==factor*final
            @test mixed.previous_power.external_power_w isa BigFloat
            @test mixed.previous_linear.inductor_power_w==final
            @test mixed.previous_linear.inductor_power_w isa BigFloat
            E._accumulate_extended_vsc_balance!(mixed,power(0.0),0.0,linear(0.0),1.0,method)
            @test mixed.energy_residual_j==(method===B.TrapezoidalCompanion ? 2final : final)
            @test mixed.previous_stored_energy_j isa Float64
            @test mixed.previous_power.external_power_w isa Float64
            @test mixed.previous_linear.inductor_power_w isa Float64
        end
    end
end

@testset "VSC live conservation survives window cancellation rollback and restart" begin
    E=AIMORA.EMTStudy
    integrator=E.configure_extended_vsc(AIMORA.SwitchDetailedVSC.ExtendedVSCParameters())
    runtime=integrator.runtime
    delta=setprecision(BigFloat,256) do;BigFloat(2)^(-70);end
    transaction=E.ExtendedVSCRuntimeTransaction(runtime)
    E._begin_extended_vsc_transaction!(transaction)
    setprecision(BigFloat,53) do
        for increment in (1.0,delta,-1.0)
            E._accumulate_extended_vsc_conservation!(runtime,increment,increment)
        end
    end
    @test runtime.cumulative_energy_residual_j==delta
    @test runtime.cumulative_linear_companion_energy_residual_j==delta
    @test precision(runtime.cumulative_energy_residual_j)>=256
    E._restore_extended_vsc_transaction!(transaction)
    @test runtime.cumulative_energy_residual_j===0.0
    @test runtime.cumulative_linear_companion_energy_residual_j===0.0
    for increment in (1.0,delta,-1.0)
        E._accumulate_extended_vsc_conservation!(runtime,increment,increment)
    end
    E._begin_extended_vsc_transaction!(transaction)
    runtime.cumulative_energy_residual_j=0.0
    runtime.cumulative_linear_companion_energy_residual_j=0.0
    E._restore_extended_vsc_transaction!(transaction)
    @test runtime.cumulative_energy_residual_j==delta
    @test runtime.cumulative_linear_companion_energy_residual_j==delta
    mktempdir() do directory
        path=joinpath(directory,"retained-conservation.checkpoint")
        E.write_extended_vsc_checkpoint(path,integrator)
        restored=E.read_extended_vsc_checkpoint(path)
        @test restored.runtime.cumulative_energy_residual_j==delta
        @test restored.runtime.cumulative_linear_companion_energy_residual_j==delta
        @test precision(restored.runtime.cumulative_energy_residual_j)==precision(runtime.cumulative_energy_residual_j)
        @test precision(restored.runtime.cumulative_linear_companion_energy_residual_j)==precision(runtime.cumulative_linear_companion_energy_residual_j)
    end
end

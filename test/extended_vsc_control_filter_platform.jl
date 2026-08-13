const ExtendedVSC = AIMORA.SwitchDetailedVSC

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
    four_wire_voltage = sequence_phase_sample(
        0.0;
        positive_magnitude_v=325.0,
        zero_magnitude_v=16.25,
    )
    four_wire_command = ExtendedVSC.compute_extended_vsc_control!(
        four_wire_state,
        ExtendedVSC.ExtendedVSCMeasurement(
            four_wire_voltage,
            (0.0, 0.0, 0.0),
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
end

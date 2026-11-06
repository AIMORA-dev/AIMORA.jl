using Test
using AIMORA
using AIMORA.NonlinearNetwork

@testset "localized tail cutoff owns the incoming current and thermal chart" begin
    nonlinear = AIMORA.Nonlinear
    for localized in (false, true), compensated in (false, true),
        method in (:BackwardEulerCompanion, :TrapezoidalCompanion),
        endpoint in (prevfloat(1e-3), 1e-3, nextfloat(1e-3))
        tail = nonlinear.TurnOffTailFidelity(1e-3; cutoff_current_a=exp(-1.0))
        tail.active = true
        tail.initial_current_a = tail.current_a = 1.0
        tail.turn_off_time_s = 0.0
        thermal = nonlinear.CauerThermalFidelity([2.0], [3.0])
        fidelity = nonlinear.PowerSemiconductorExtendedFidelity(
            turn_off_tail=tail, thermal=thermal)
        device = nonlinear.IGBTSwitch(1, 0; off_conductance=1e-6,
            extended_fidelity=fidelity)
        localized && nonlinear.power_semiconductor_event_localization!(device)
        step = 1e-5
        prepare_nonlinear_device_step!(device, endpoint, step, method)
        voltage, correction = [10.0, 0.0], [1e-17, 0.0]
        current, jacobian = zeros(2), zeros(2, 2)
        if compensated
            nonlinear.nonlinear_current_jacobian_compensated!(current, jacobian,
                device, voltage, correction, endpoint)
        else
            nonlinear_current_jacobian!(current, jacobian, device, voltage, endpoint)
        end
        exponential = exp(-endpoint / tail.decay_time_s)
        expected_tail = localized || exponential > tail.cutoff_current_a ? exponential : 0.0
        @test current[1] ≈ 1e-5 + expected_tail rtol=2eps(Float64)
        initial_temperature = thermal.node_temperature_k[1]
        nonlinear._accept_power_semiconductor_state!(device, voltage, current,
            endpoint, jacobian, compensated ? correction : nothing)
        expected_loss = 10.0 * (1e-5 + expected_tail)
        @test device.last_semiconductor_loss_w ≈ expected_loss rtol=2eps(Float64)
        # Independent one-stage backward-Euler thermal balance (both electrical methods).
        rise = initial_temperature - thermal.ambient_temperature_k
        expected_temperature = thermal.ambient_temperature_k +
            (2.0rise + step * expected_loss) / (2.0 + step / 3.0)
        @test thermal.node_temperature_k[1] ≈ expected_temperature rtol=2eps(Float64)
        @test thermal.cumulative_input_energy_j ≈ step * expected_loss rtol=2eps(Float64)
        @test tail.active == (localized || exponential > tail.cutoff_current_a)
        if localized
            @test tail.cutoff_event_count == 0
            @test tail.current_a == exponential
            nonlinear.apply_power_semiconductor_tail_cutoff!(device, endpoint)
            @test !tail.active && tail.current_a == 0.0 && tail.initial_current_a == 0.0
            @test tail.cutoff_event_count == 1
            @test_throws ArgumentError nonlinear.apply_power_semiconductor_tail_cutoff!(device, endpoint)
        end
        finish_nonlinear_device_step!(device)
    end
end

@testset "semiconductor work follows the admitted integration method" begin
    for (method,expected_energy) in ((:BackwardEulerCompanion,0.4),(:TrapezoidalCompanion,0.65))
        fidelity = AIMORA.Nonlinear.PowerSemiconductorExtendedFidelity(
            thermal=AIMORA.Nonlinear.CauerThermalFidelity([1.0],[1.0]))
        device = AIMORA.Nonlinear.DiodeValveSwitch(1,0;on_conductance=1.0,
            off_conductance=1.0,extended_fidelity=fidelity)
        # Unit conductance: v=i. A backward-Euler substep uses the right-limit
        # endpoint power 2*2 W; trapezoidal retains both 3*3 and 2*2 W samples.
        fidelity.previous_terminal_voltage_v = 3.0
        fidelity.previous_terminal_current_a = 3.0
        device.previous_semiconductor_loss_w = 9.0
        device.last_semiconductor_loss_w = 9.0
        step = 0.1
        prepare_nonlinear_device_step!(device,step,step,method)
        current,jacobian = zeros(2),zeros(2,2)
        nonlinear_current_jacobian!(current,jacobian,device,[2.0,0.0],step)
        @test current == [2.0,-2.0]
        accept_nonlinear_device_state!(device,[2.0,0.0],current,jacobian,step)
        @test device.semiconductor_dissipated_energy_j ≈ expected_energy atol=1e-15 rtol=0
        @test fidelity.companion_energy_residual_j ≈ 0.0 atol=1e-15
        @test device.last_semiconductor_loss_w == 4.0
        @test fidelity.thermal.cumulative_input_energy_j ≈ 0.4 atol=1e-15
        finish_nonlinear_device_step!(device)
        simple = AIMORA.Nonlinear.DiodeValveSwitch(1,0;on_conductance=1.0,off_conductance=1.0)
        simple.last_conductance = 1.0
        simple.previous_semiconductor_loss_w = 9.0
        integration = method === :BackwardEulerCompanion ?
            Val(AIMORA.Branches.BackwardEulerCompanion) : Val(AIMORA.Branches.TrapezoidalCompanion)
        AIMORA.Branches.update!(simple,[2.0],step,integration)
        @test simple.semiconductor_dissipated_energy_j ≈ expected_energy atol=1e-15 rtol=0
    end
end

@testset "series RC snubber work follows its companion method" begin
    for implicit in (false, true), initial_current in (0.0, 2.0), voltage in (-100.0, 30.0, 100.0)
        resistance, capacitance, step = 10.0, 1e-6, 5e-6
        initial_voltage = 20.0
        snubber = AIMORA.Nonlinear.SeriesRCSnubber(resistance, capacitance;
            initial_current_a=initial_current, initial_capacitor_voltage_v=initial_voltage)
        old_power = resistance * initial_current^2
        old_branch_voltage = initial_voltage + resistance * initial_current
        if implicit
            AIMORA.Nonlinear._update_series_rc_snubber_backward_euler!(snubber, voltage, step)
        else
            AIMORA.Nonlinear._update_series_rc_snubber!(snubber, voltage, step)
        end
        current = snubber.last_current_a
        expected_loss = step * resistance * current^2
        implicit || (expected_loss = 0.5step * (old_power + resistance * current^2))
        @test snubber.dissipated_energy_j ≈ expected_loss rtol=2e-15
        expected_voltage = initial_voltage + step / capacitance *
            (implicit ? current : 0.5 * (initial_current + current))
        @test snubber.capacitor_voltage_v ≈ expected_voltage rtol=2e-15
        # This identity distinguishes physical resistor loss from BE numerical
        # dissipation; it does not subtract numerical loss from system acceptance.
        work = implicit ? step * voltage * current :
            0.5step * (old_branch_voltage * initial_current + voltage * current)
        stored_change = 0.5capacitance * (snubber.capacitor_voltage_v^2 - initial_voltage^2)
        discretization_defect = implicit ?
            0.5capacitance * (snubber.capacitor_voltage_v - initial_voltage)^2 :
            step * (snubber.capacitor_voltage_v - initial_voltage) *
                (current - initial_current) / 4
        @test work - stored_change - expected_loss ≈ discretization_defect atol=2e-17 rtol=2e-13
        @test snubber.dissipated_energy_j >= 0.0
    end
end

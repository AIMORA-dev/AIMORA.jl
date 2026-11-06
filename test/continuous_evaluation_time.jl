module ContinuousEvaluationTimeTests
using Test, AIMORA
using AIMORA.StudyCore: ContinuousEvaluationTime, continuous_time_seconds, continuous_time_elapsed

@testset "continuous time keeps elapsed coordinates below the absolute clock lattice" begin
    for origin in (0.0, 0.0150125, 1.0, 4096.0, -1.0), elapsed in (1.6382184114200127e-6, 5e-6)
        time = ContinuousEvaluationTime(origin, elapsed)
        setprecision(BigFloat, 192) do
            @test continuous_time_seconds(time) == BigFloat(origin) + BigFloat(elapsed)
            @test continuous_time_elapsed(time, origin) == BigFloat(elapsed)
        end
        @test time.high_s == origin + elapsed
    end
    @test_throws ArgumentError ContinuousEvaluationTime(Inf, 1.0)
    @test_throws ArgumentError ContinuousEvaluationTime(1.0, NaN)
    @test_throws ArgumentError ContinuousEvaluationTime(floatmax(Float64), floatmax(Float64))
    @test_throws MethodError Float64(ContinuousEvaluationTime(1.0, 1e-20))
end

@testset "pending gate surfaces retain the continuous left endpoint" begin
    nonlinear=AIMORA.Nonlinear
    device=nonlinear.IGBTSwitch(1,0)
    driver=device.gate_driver
    driver.pending_state=true
    driver.pending_transition_time_s=1.0
    before=deepcopy(driver)
    surface=nonlinear._POWER_SEMICONDUCTOR_GATE_TRANSITION_EVENT
    for offset in (1e-20,-1e-20,0.0)
        time=ContinuousEvaluationTime(1.0,offset)
        expected=setprecision(BigFloat,192) do
            BigFloat(driver.pending_transition_time_s)-continuous_time_seconds(time)
        end
        @test nonlinear.power_semiconductor_gate_transition_residual(device,time)==expected
        @test AIMORA.NonlinearNetwork.nonlinear_device_event_value(surface,device,time)==Float64(expected)
        @test all(f->getfield(driver,f)==getfield(before,f),fieldnames(typeof(driver)))
    end
    @test nonlinear.power_semiconductor_gate_transition_residual(device,1.0)==0.0
    driver.pending_state=nothing
    @test nonlinear.power_semiconductor_gate_transition_residual(device,ContinuousEvaluationTime(1.0,1e-20))===nothing
end

@testset "scheduled switch stamps keep both sides of exact deadlines" begin
    switches=AIMORA.Switches;branches=AIMORA.Branches
    switch=switches.TimeSwitch(1,0;close_time_s=1.0,open_time_s=2.0,
        on_conductance=20.0,off_conductance=1e-6)
    before=deepcopy(switch)
    for (deadline,low,expected) in ((1.0,-1e-20,false),(1.0,0.0,true),(1.0,1e-20,true),
        (2.0,-1e-20,true),(2.0,0.0,false),(2.0,1e-20,false))
        time=ContinuousEvaluationTime(deadline,low)
        @test time.high_s==deadline
        @test switches.switch_closed(switch,time)==expected
        for method in (branches.TrapezoidalCompanion,branches.BackwardEulerCompanion)
            matrix=zeros(1,1);rhs=zeros(1)
            branches.stamp!(matrix,rhs,switch,time,1e-5,Val(method))
            @test matrix[1,1]==(expected ? 20.0 : 1e-6)
            @test iszero(rhs[1])
        end
    end
    @test all(k->isequal(getfield(switch,k),getfield(before,k)),fieldnames(typeof(switch)))
    simultaneous=switches.TimeSwitch(1,0;close_time_s=1.0,open_time_s=1.0)
    @test !switches.switch_closed(simultaneous,ContinuousEvaluationTime(1.0,1e-20))
    initially=switches.TimeSwitch(1,0;initially_closed=true,open_time_s=1.0)
    @test switches.switch_closed(initially,ContinuousEvaluationTime(1.0,-1e-20))
    @test !switches.switch_closed(initially,ContinuousEvaluationTime(1.0,1e-20))
    for owner in (switches.IdealSwitch(1,0,true),switches.CurrentZeroSwitch(initially))
        @test switches.switch_closed(owner,ContinuousEvaluationTime(1.0,1e-20))
        matrix=zeros(1,1);branches.stamp!(matrix,zeros(1),owner,ContinuousEvaluationTime(1.0,-1e-20),1e-5)
        @test matrix[1,1]==owner.on_conductance
    end
end

if AIMORA.solver_available()
    @testset "continuous source stamps retain precision until RHS storage" begin
        branches = AIMORA.Branches
        time = ContinuousEvaluationTime(0.0150125, 1.6382184114200127e-6)
        signal(endpoint::ContinuousEvaluationTime) = continuous_time_seconds(endpoint)
        method = Val(branches.BackwardEulerCompanion)
        setprecision(BigFloat, 192) do
            value = signal(time)
            for source in (branches.TheveninSource(1, 2.0, signal),
                           branches.TwoTerminalTheveninSource(1, 2, 2.0, signal),
                           branches.CurrentInjection(1, signal))
                y = zeros(BigFloat, 2, 2)
                rhs = zeros(BigFloat, 2)
                branches.stamp!(y, rhs, source, time, 5e-6, method)
                expected = source isa branches.CurrentInjection ? value : 2value
                @test rhs[1] == expected
                @test rhs[2] == (source isa branches.TwoTerminalTheveninSource ? -expected : 0)
                @test rhs[1] != BigFloat(Float64(expected))
            end
        end
        unsupported = branches.CurrentInjection(1, (_::Float64) -> 1.0)
        @test_throws ArgumentError branches.stamp!(zeros(1, 1), zeros(1),
            unsupported, time, 5e-6, method)
        rhs = zeros(1)
        branches.stamp!(zeros(1, 1), rhs, unsupported,
            ContinuousEvaluationTime(1.0, 0.0), 5e-6, method)
        @test rhs == [1.0]
    end

    @testset "source and semiconductor consume the same continuous evaluation endpoint" begin
        origin = 0.0150125
        elapsed = 1.6382184114200127e-6
        time = ContinuousEvaluationTime(origin, elapsed)
        nonlinear = AIMORA.Nonlinear
        tail = nonlinear.TurnOffTailFidelity(1e-3)
        tail.active = true
        tail.initial_current_a = 100.0
        tail.turn_off_time_s = origin - 1e-5
        device = nonlinear.IGBTSwitch(1, 0;
            extended_fidelity=nonlinear.PowerSemiconductorExtendedFidelity(turn_off_tail=tail))
        AIMORA.NonlinearNetwork.prepare_nonlinear_device_step!(device, time.high_s, elapsed, :BackwardEulerCompanion)
        fidelity = device.extended_fidelity
        setprecision(BigFloat, 106) do
            current, jacobian = nonlinear._power_semiconductor_compensated_dynamic_current_jacobian(
                fidelity, BigFloat(0), time, true)
            expected = setprecision(BigFloat, 256) do
                exact_time = BigFloat(origin) + BigFloat(elapsed)
                BigFloat(tail.initial_current_a) * exp(-(exact_time - BigFloat(tail.turn_off_time_s))/BigFloat(tail.decay_time_s))
            end
            @test abs(current - expected) <= BigFloat(2)^(-70)
            @test iszero(jacobian)
            legacy, _ = nonlinear._power_semiconductor_compensated_dynamic_current_jacobian(
                fidelity, BigFloat(0), time.high_s, true)
            @test legacy == BigFloat(tail.initial_current_a * exp(-max(0.0, time.high_s-tail.turn_off_time_s)/tail.decay_time_s))
            @test current != legacy
            terminal_current = zeros(2);terminal_jacobian = zeros(2,2)
            AIMORA.NonlinearNetwork.nonlinear_current_jacobian!(terminal_current,terminal_jacobian,
                device,zeros(2),time)
            @test terminal_current[1] == -terminal_current[2]
            @test terminal_current[1] == Float64(current)
        end
        AIMORA.NonlinearNetwork.refresh_nonlinear_device_event_dependents!(device, time)
        expected_tail = setprecision(BigFloat, 256) do
            Float64(BigFloat(tail.initial_current_a) * exp(-
                (BigFloat(origin) + BigFloat(elapsed) - BigFloat(tail.turn_off_time_s)) /
                BigFloat(tail.decay_time_s)))
        end
        @test tail.current_a == expected_tail
        @test device.last_current == expected_tail
        @test device.extended_fidelity.previous_terminal_current_a == expected_tail
        diode = nonlinear.DiodeValveSwitch(1, 0;
            extended_fidelity=nonlinear.PowerSemiconductorExtendedFidelity(
                recovered_charge=nonlinear.RecoveredChargeFidelity(1e-3)))
        diode.last_voltage = 1e-11
        surface = nonlinear._POWER_SEMICONDUCTOR_FORWARD_TURN_ON_EVENT
        @test 0 <= AIMORA.NonlinearNetwork.nonlinear_device_event_value(surface, diode, time) <= 1e-10
        AIMORA.NonlinearNetwork.apply_nonlinear_device_event!(surface, diode, time)
        AIMORA.NonlinearNetwork.refresh_nonlinear_device_event_dependents!(diode, time)
        @test diode.closed
        @test diode.last_transition_time_s == time.high_s
        @test diode.last_evaluation_time_s == time.high_s
        journal = Any[]
        before = AIMORA.NonlinearNetwork.nonlinear_device_event_transition_snapshot(diode)
        AIMORA.NonlinearNetwork.record_nonlinear_device_event_transition_observation!(
            journal, diode, before, before, surface, time, :left)
        @test only(journal).time_s === time
        parameters = AIMORA.SwitchDetailedVSC.ExtendedVSCParameters()
        for phase in 1:3
            source = AIMORA.EMTStudy.ExtendedVSCGridSource(phase,0,phase,0.2,parameters)
            voltage = setprecision(BigFloat, 106) do
                AIMORA.EMTStudy._extended_vsc_grid_source_voltage(source,time)
            end
            expected = setprecision(BigFloat, 256) do
                exact_time = BigFloat(origin) + BigFloat(elapsed)
                angle = BigFloat(2pi*parameters.frequency_hz)*exact_time
                offset = BigFloat((phase-1)*2pi/3)
                scenario = parameters.scenario
                # Independent complex-phasor superposition includes every
                # configured sequence/harmonic, not a fundamental-only oracle.
                waveform = real(exp(im*angle)*exp(-im*offset)) +
                    BigFloat(scenario.negative_sequence_voltage_ratio)*real(exp(im*angle)*exp(im*offset)) +
                    BigFloat(scenario.zero_sequence_voltage_ratio)*real(exp(im*angle))
                for (order,ratio) in zip(scenario.harmonic_orders,scenario.harmonic_voltage_ratios)
                    offset_sign = mod(order,6)==5 ? -1 : 1
                    waveform += BigFloat(ratio)*real(exp(im*order*angle)*exp(-im*offset_sign*offset))
                end
                voltage = BigFloat(sqrt(2/3)*parameters.grid_line_line_rms_v)*waveform
                faulted = scenario.fault_kind===:three_phase || phase in scenario.faulted_phases
                scenario.fault_start_s<=exact_time<scenario.fault_end_s && faulted &&
                    (voltage *= BigFloat(scenario.fault_voltage_factor))
                voltage
            end
            @test abs(voltage-expected) <= BigFloat(2)^(-60)
        end
    end
end
end

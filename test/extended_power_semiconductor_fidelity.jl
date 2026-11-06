using Test
using AIMORA.NonlinearNetwork

@testset "continuous semiconductor event callbacks preserve discrete ownership" begin
    nonlinear = AIMORA.Nonlinear
    time = AIMORA.StudyCore.ContinuousEvaluationTime(0.0133875, 5.16352205790695e-10)
    @test !iszero(time.low_s)
    tail = nonlinear.TurnOffTailFidelity(1.0e-3; cutoff_current_a=0.1)
    device = nonlinear.IGBTSwitch(1, 0;
        antiparallel_diode=nonlinear.AntiparallelDiodeParameters(),
        extended_fidelity=nonlinear.PowerSemiconductorExtendedFidelity(turn_off_tail=tail))
    device.closed=true; device.last_forward_current_a=0.5
    reference=deepcopy(device)
    nonlinear.apply_power_semiconductor_forward_extinction!(reference, time.high_s)
    nonlinear.apply_power_semiconductor_forward_extinction!(device, time)
    @test !device.closed && device.topology_transition_count == reference.topology_transition_count
    @test device.last_transition_time_s == reference.last_transition_time_s == time.high_s
    @test device.extended_fidelity.turn_off_tail.current_a == reference.extended_fidelity.turn_off_tail.current_a
    physical_time=device.extended_fidelity.turn_off_tail.turn_off_time_s
    rounded_time=reference.extended_fidelity.turn_off_tail.turn_off_time_s
    @test Rational{BigInt}(physical_time.high_s)+Rational{BigInt}(physical_time.low_s) ==
        Rational{BigInt}(time.high_s)+Rational{BigInt}(time.low_s)
    @test Rational{BigInt}(rounded_time.high_s)+Rational{BigInt}(rounded_time.low_s) == Rational{BigInt}(time.high_s)
    @test physical_time.event_time_s == rounded_time.event_time_s == time.high_s
    @test physical_time != rounded_time
    nonlinear.apply_power_semiconductor_reverse_turn_on!(reference, time.high_s)
    nonlinear.apply_power_semiconductor_reverse_turn_on!(device, time)
    @test device.reverse_diode_conducting && device.topology_transition_count == reference.topology_transition_count
    nonlinear.apply_power_semiconductor_reverse_extinction!(reference, time.high_s)
    nonlinear.apply_power_semiconductor_reverse_extinction!(device, time)
    @test !device.reverse_diode_conducting && device.topology_transition_count == reference.topology_transition_count
    for callback in (nonlinear.apply_power_semiconductor_recovery_zero!, nonlinear.apply_power_semiconductor_tail_cutoff!)
        @test applicable(callback, device, time)
    end
    @test_throws ArgumentError nonlinear.apply_power_semiconductor_recovery_zero!(device, time)
    cutoff = nonlinear.power_semiconductor_tail_cutoff_time(device)
    @test cutoff !== nothing
    nonlinear.apply_power_semiconductor_tail_cutoff!(device, AIMORA.StudyCore.ContinuousEvaluationTime(cutoff, 0.0))
    @test !device.extended_fidelity.turn_off_tail.active
    @test device.extended_fidelity.turn_off_tail.cutoff_event_count == 1
    @test_throws ArgumentError nonlinear.apply_power_semiconductor_tail_cutoff!(device, AIMORA.StudyCore.ContinuousEvaluationTime(cutoff, 0.0))
end

@testset "two-level bridge selection uses mapped channel current" begin
    nonlinear = AIMORA.Nonlinear
    controlled = [nonlinear.IGBTSwitch(1, 0) for _ in 1:6]
    diodes = [nonlinear.DiodeValveSwitch(1, 0) for _ in 1:6]
    diodes[6].last_voltage = nonlinear.power_semiconductor_forward_turn_on_voltage(diodes[6])
    mapped = deepcopy(vcat(controlled, diodes))
    nonlinear.apply_power_semiconductor_forward_turn_on!(mapped[12], 0.0)
    mapped[12].last_forward_current_a = -0.15129709758684884
    @test !AIMORA.EMTStudy._synchronize_two_level_bridge_conduction!(
        controlled, diodes, (0.0, 0.0, 50.183771907642864), 0.0;
        mapped_devices=mapped,
    )
    @test !diodes[6].closed
    @test diodes[6].topology_transition_count == 0
    mapped[12].last_forward_current_a = 0.25
    @test AIMORA.EMTStudy._synchronize_two_level_bridge_conduction!(
        controlled, diodes, (0.0, 0.0, 50.183771907642864), 0.0;
        mapped_devices=mapped,
    )
    @test diodes[6].closed
    neutral_controlled = [nonlinear.IGBTSwitch(1, 0) for _ in 1:6]
    neutral_diodes = [nonlinear.DiodeValveSwitch(1, 0) for _ in 1:6]
    neutral_diodes[1].last_voltage = nonlinear.power_semiconductor_forward_turn_on_voltage(neutral_diodes[1])
    @test !AIMORA.EMTStudy._synchronize_two_level_bridge_conduction!(
        neutral_controlled, neutral_diodes, (-1.0, 0.0, 0.0), 0.0;
        excluded_positions=(7,),
    )
    @test all(!device.closed for device in neutral_diodes)
    @test all(iszero(device.topology_transition_count) for device in neutral_diodes)
    @test AIMORA.EMTStudy._synchronize_two_level_bridge_conduction!(
        neutral_controlled, neutral_diodes, (-1.0, 0.0, 0.0), 0.0,
    )
    @test neutral_diodes[1].closed
    @test_throws BoundsError AIMORA.EMTStudy._synchronize_two_level_bridge_conduction!(
        neutral_controlled, neutral_diodes, (0.0, 0.0, 0.0), 0.0;
        excluded_positions=(13,),
    )
end

@testset "forward completion preserves unaffected two-level bridge paths" begin
    nonlinear = AIMORA.Nonlinear
    emt = AIMORA.EMTStudy
    controlled = [nonlinear.IGBTSwitch(1, 0) for _ in 1:6]
    diodes = [nonlinear.DiodeValveSwitch(1, 0) for _ in 1:6]
    runtime = (; controlled_valves=controlled, freewheel_diodes=diodes)
    chart = ntuple(position -> position == 10, 12)
    @test emt._complete_two_level_bridge_conduction_chart!(controlled, diodes, chart, 0.0)
    @test !emt._complete_two_level_bridge_conduction_chart!(controlled, diodes, chart, 0.0)
    transitions = getfield.(vcat(controlled, diodes), :topology_transition_count)
    extended_chart = ntuple(position -> position in (8, 10), 12)
    @test emt._complete_two_level_bridge_conduction_chart!(
        controlled, diodes, extended_chart, 0.0)
    @test Tuple(device.closed for device in vcat(controlled, diodes)) ==
        ntuple(position -> position in (8, 10), 12)
    @test diodes[4].topology_transition_count == transitions[10]
    @test diodes[5].topology_transition_count == 0
    @test all(!device.closed for device in controlled)
    before = deepcopy(vcat(controlled, diodes))
    conflicting = ntuple(position -> position in (3, 10), 12)
    @test_throws ArgumentError emt._complete_two_level_bridge_conduction_chart!(
        controlled, diodes, conflicting, 0.0)
    @test_throws DimensionMismatch emt._complete_two_level_bridge_conduction_chart!(
        controlled, diodes, chart[1:11], 0.0)
    @test_throws ArgumentError emt._complete_two_level_bridge_conduction_chart!(
        controlled, diodes, Int.(chart), 0.0)
    @test getfield.(vcat(controlled, diodes), :closed) == getfield.(before, :closed)
    @test getfield.(vcat(controlled, diodes), :topology_transition_count) ==
        getfield.(before, :topology_transition_count)
    @test_throws BoundsError emt._complete_two_level_bridge_conduction_chart!(
        controlled, diodes, chart, 0.0; extinguished_positions=(13,))
    @test_throws ArgumentError emt._complete_two_level_bridge_conduction_chart!(
        controlled, diodes, chart, 0.0; extinguished_positions=(7,))
    # Two failed holding paths do not authorize changing the valid third leg.
    mixed = ntuple(position -> position in (7, 10, 12), 12)
    held_transition = diodes[4].topology_transition_count
    @test emt._complete_two_level_bridge_conduction_chart!(
        controlled, diodes, mixed, 0.0; extinguished_positions=(7, 12))
    @test Tuple(device.closed for device in vcat(controlled, diodes)) == chart
    @test diodes[4].topology_transition_count == held_transition
    @test !emt._complete_two_level_bridge_conduction_chart!(
        controlled, diodes, mixed, 0.0; extinguished_positions=(7, 12))
end

@testset "semiconductor event derivative rebase preserves physical storage" begin
    nonlinear = AIMORA.Nonlinear
    junction = nonlinear.NonlinearJunctionChargeFidelity(
        1.0e-6, 10.0, 0.5; voltage_domain_v=(-100.0, 100.0))
    device = nonlinear.IGBTSwitch(1, 0;
        extended_fidelity=nonlinear.PowerSemiconductorExtendedFidelity(
            junction_charge=junction))
    nonlinear.reconstruct_power_semiconductor_boundary_state!(device, 0.0, 2.0)
    before_charge = junction.last_charge_c
    before_energy = nonlinear.power_semiconductor_extended_state(device).
        junction_stored_energy_j
    before_terminal_current = device.last_current
    before_displacement_current = junction.last_displacement_current_a
    before_history_current = device.extended_fidelity.previous_terminal_current_a
    @test nonlinear_device_event_capacitance(device) ==
        (1, 0, junction.last_capacitance_f)
    @test nonlinear_device_event_current(device) == before_terminal_current
    rebase_nonlinear_device_event_derivative!(device, 0.75)
    @test junction.last_displacement_current_a ==
        before_displacement_current + 0.75
    @test device.last_current == before_terminal_current + 0.75
    @test device.extended_fidelity.previous_terminal_current_a ==
        before_history_current + 0.75
    @test junction.last_charge_c == before_charge
    @test nonlinear.power_semiconductor_extended_state(device).
        junction_stored_energy_j == before_energy
end

@testset "self-commutated gate and conduction ownership" begin
    nonlinear = AIMORA.Nonlinear
    switch = nonlinear.IGBTSwitch(
        1, 0; gate_driver=nonlinear.PowerSemiconductorGateDriver(
            turn_on_delay_s=1.0e-6,
        ),
    )
    nonlinear.request_power_semiconductor_gate!(switch, true, 0.0)
    ordinary = deepcopy(switch)
    @test nonlinear.apply_power_semiconductor_gate_transition!(
        switch, 1.0e-6; defer_conduction=true,
    )
    @test switch.gate_driver.applied_on
    @test !switch.closed
    @test switch.topology_transition_count == 0
    @test nonlinear.apply_power_semiconductor_gate_transition!(ordinary, 1.0e-6)
    @test ordinary.closed
    @test ordinary.topology_transition_count == 1
    nonlinear.apply_power_semiconductor_forward_turn_on!(switch, 1.0e-6)
    @test switch.closed
    @test switch.topology_transition_count == 1
end

@testset "bridge stamp preserves localized gate-event ownership" begin
    nonlinear = AIMORA.Nonlinear
    descriptor = AIMORA.BridgeTopologies.two_level_bridge_topology(
        [3, 4, 5], 1, 2,
    )
    valves = [nonlinear.IGBTSwitch(
        position.from_node,
        position.to_node;
        gate_driver=nonlinear.PowerSemiconductorGateDriver(
            turn_on_delay_s=1.0e-6,
        ),
    ) for position in descriptor.valve_positions]
    bridge = nonlinear.PowerSemiconductorBridgeTopology(descriptor, valves)
    nonlinear.power_semiconductor_event_localization!(valves[1])
    nonlinear.request_power_semiconductor_gate!(valves[1], true, 0.0)
    nonlinear.request_power_semiconductor_gate!(valves[3], true, 0.0)
    AIMORA.Branches.stamp!(
        zeros(5, 5), zeros(5), bridge, 1.0e-6, 1.0e-6,
        Val(AIMORA.Branches.BackwardEulerCompanion),
    )
    @test valves[1].gate_driver.pending_state === true
    @test !valves[1].gate_driver.applied_on
    @test !valves[1].closed
    @test valves[3].gate_driver.pending_state === nothing
    @test valves[3].gate_driver.applied_on
    @test valves[3].closed
    @test nonlinear.apply_power_semiconductor_bridge_gate_transitions!(
        bridge, 1.0e-6,
    ) == 1
    @test valves[1].gate_driver.applied_on
end

@testset "recovery-root allowance uses matching history stages" begin
    tolerance = AIMORA.Nonlinear._power_semiconductor_recovery_zero_tolerance
    for charge in (0.0,4.514383086639529e-7,2.0),
        current in (0.0,-11922.081314353127,1.0), step in (3.7867575883865364e-11,1e-6)
        old_formula = max(8192eps(Float64)*max(1.0,abs(charge),abs(current)*step,eps(Float64)),
            abs(current)*1e-12)
        @test tolerance(charge,current,step) == old_formula
    end
    charge = 4.514383086639529e-7
    current = -11922.081314353127
    step = 3.7867575883865364e-11
    residual = -2.2010192592720176e-11
    @test abs(residual) > tolerance(0.0,0.0,step)
    @test abs(residual) <= tolerance(charge,current,step)
    surface = AIMORA.Nonlinear._POWER_SEMICONDUCTOR_RECOVERY_ZERO_EVENT
    device = AIMORA.Nonlinear.IGBTSwitch(1, 0)
    @test nonlinear_device_event_value_tolerance(surface, device, 1.0e-10) ==
        8192eps(Float64)
    @test nonlinear_device_event_value_tolerance(surface, device, 1.0e-13) ==
        1.0e-13
    @test nonlinear_device_event_value_tolerance(
        AIMORA.Nonlinear._POWER_SEMICONDUCTOR_FORWARD_TURN_ON_EVENT,
        device,
        1.0e-10,
    ) == 1.0e-10
end

@testset "recovered-charge discrete surface linearization" begin
    nonlinear = AIMORA.Nonlinear
    lifetime_s = 1.0e-6
    charge_surface(q, current, step) =
        (q + step * current) / (1 + step / BigFloat(lifetime_s))
    for stored_charge_c in (0.0, 5.0e-7, 5.0e-7 + 1.5e-12),
        current_a in (-1.0, 0.0, 2.0), step_s in (5.0e-7, 2.0e-6)
        recovery = nonlinear.RecoveredChargeFidelity(
            lifetime_s; initial_charge_c=stored_charge_c,
        )
        linearization = nonlinear._power_semiconductor_recovery_charge_linearization(
            recovery, current_a, step_s,
        )
        @test isequal(linearization.charge_c,
            (stored_charge_c + step_s * current_a) / (1 + step_s / lifetime_s))
        values = BigFloat.([stored_charge_c, current_a, step_s])
        epsilon = big"1e-24"
        partials = map(eachindex(values)) do index
            positive, negative = copy(values), copy(values)
            positive[index] += epsilon
            negative[index] -= epsilon
            Float64((charge_surface(positive...) - charge_surface(negative...)) /
                (2epsilon))
        end
        @test linearization.stored_charge_derivative ≈ partials[1] rtol=1.0e-14
        @test linearization.current_derivative_s ≈ partials[2] rtol=1.0e-14
        @test linearization.step_derivative_a ≈ partials[3] rtol=1.0e-14
        @test recovery.stored_charge_c == stored_charge_c
        @test recovery.recovery_zero_event_count == 0
    end
end

function synthetic_semiconductor_energy_table()
    current_axis_a = [0.0, 10.0]
    voltage_axis_v = [0.0, 100.0]
    temperature_axis_k = [250.0, 350.0]
    values = Array{Float64}(undef, 2, 2, 2)
    for current_index in 1:2, voltage_index in 1:2, temperature_index in 1:2
        values[current_index, voltage_index, temperature_index] =
            1.0e-6 * (
                current_axis_a[current_index] +
                voltage_axis_v[voltage_index] / 10.0 +
                (temperature_axis_k[temperature_index] - 250.0) / 10.0
            )
    end
    return AIMORA.Nonlinear.SwitchingEnergyTable(
        current_axis_a,
        voltage_axis_v,
        temperature_axis_k;
        turn_on_energy_j=values,
        turn_off_energy_j=2.0 .* values,
        reverse_recovery_energy_j=3.0 .* values,
    )
end

@testset "extended power-semiconductor constitutive fidelity" begin
    nonlinear_charge = AIMORA.Nonlinear.NonlinearJunctionChargeFidelity(
        2.0e-6,
        10.0,
        0.4;
        voltage_domain_v=(-100.0, 100.0),
    )
    for voltage_v in (-50.0, -1.0, 0.0, 20.0)
        perturbation_v = 1.0e-6
        central_charge_derivative = (
            AIMORA.Nonlinear.power_semiconductor_junction_charge(
                nonlinear_charge,
                voltage_v + perturbation_v,
            ) -
            AIMORA.Nonlinear.power_semiconductor_junction_charge(
                nonlinear_charge,
                voltage_v - perturbation_v,
            )
        ) / (2.0 * perturbation_v)
        @test central_charge_derivative ≈
            AIMORA.Nonlinear.power_semiconductor_junction_capacitance(
                nonlinear_charge,
                voltage_v,
            ) rtol=1.0e-7
    end

    energy_table = synthetic_semiconductor_energy_table()
    @test AIMORA.Nonlinear.power_semiconductor_switching_energy(
        energy_table,
        :turn_on,
        10.0,
        100.0,
        350.0,
    ) == energy_table.turn_on_energy_j[2, 2, 2]
    @test_throws DomainError AIMORA.Nonlinear.power_semiconductor_switching_energy(
        energy_table,
        :turn_on,
        11.0,
        100.0,
        350.0,
    )

    fidelity = AIMORA.Nonlinear.PowerSemiconductorExtendedFidelity(
        recovered_charge=AIMORA.Nonlinear.RecoveredChargeFidelity(
            2.0e-6;
            initial_charge_c=2.0e-6,
        ),
        junction_charge=nonlinear_charge,
        thermal=AIMORA.Nonlinear.CauerThermalFidelity([1.0, 2.0], [0.5, 1.0]),
    )
    recovered_diode = AIMORA.Nonlinear.DiodeValveSwitch(
        1,
        0;
        on_conductance=1.0,
        off_conductance=0.0,
        extended_fidelity=fidelity,
    )
    prepare_nonlinear_device_step!(
        recovered_diode,
        1.0e-6,
        1.0e-6,
        :TrapezoidalCompanion,
    )
    terminal_current_a = zeros(Float64, 2)
    terminal_jacobian_s = zeros(Float64, 2, 2)
    nonlinear_current_jacobian!(
        terminal_current_a,
        terminal_jacobian_s,
        recovered_diode,
        [-1.0, 0.0],
        1.0e-6,
    )
    perturbation_v = 1.0e-7
    positive_current_a = zeros(Float64, 2)
    negative_current_a = zeros(Float64, 2)
    nonlinear_current_jacobian!(
        positive_current_a,
        zeros(Float64, 2, 2),
        recovered_diode,
        [-1.0 + perturbation_v, 0.0],
        1.0e-6,
    )
    nonlinear_current_jacobian!(
        negative_current_a,
        zeros(Float64, 2, 2),
        recovered_diode,
        [-1.0 - perturbation_v, 0.0],
        1.0e-6,
    )
    @test terminal_jacobian_s[1, 1] ≈
        (positive_current_a[1] - negative_current_a[1]) / (2.0 * perturbation_v) rtol=1.0e-7
    fidelity_snapshot = sprint(show, fidelity)
    nonlinear_current_jacobian!(
        terminal_current_a,
        terminal_jacobian_s,
        recovered_diode,
        [-1.0, 0.0],
        1.0e-6,
    )
    @test sprint(show, fidelity) == fidelity_snapshot
    fallback_acceptance_diode = deepcopy(recovered_diode)
    jacobian_acceptance_diode = deepcopy(recovered_diode)
    accept_nonlinear_device_state!(
        fallback_acceptance_diode,
        [-1.0, 0.0],
        terminal_current_a,
        1.0e-6,
    )
    accept_nonlinear_device_state!(
        jacobian_acceptance_diode,
        [-1.0, 0.0],
        terminal_current_a,
        terminal_jacobian_s,
        1.0e-6,
    )
    finish_nonlinear_device_step!(fallback_acceptance_diode)
    finish_nonlinear_device_step!(jacobian_acceptance_diode)
    fallback_acceptance_state = AIMORA.Nonlinear.power_semiconductor_extended_state(
        fallback_acceptance_diode,
    )
    jacobian_acceptance_state = AIMORA.Nonlinear.power_semiconductor_extended_state(
        jacobian_acceptance_diode,
    )
    for field in fieldnames(typeof(fallback_acceptance_state))
        fallback_value = getfield(fallback_acceptance_state, field)
        jacobian_value = getfield(jacobian_acceptance_state, field)
        if fallback_value isa AbstractArray
            @test fallback_value ≈ jacobian_value rtol=16 * eps(Float64) atol=16 * eps(Float64)
        elseif fallback_value isa AbstractFloat
            @test isequal(fallback_value, jacobian_value) ||
                isapprox(
                    fallback_value,
                    jacobian_value;
                    rtol=16 * eps(Float64),
                    atol=16 * eps(Float64),
                )
        else
            @test fallback_value == jacobian_value
        end
    end
    @test fallback_acceptance_diode.last_voltage ≈ jacobian_acceptance_diode.last_voltage
    @test fallback_acceptance_diode.last_current ≈ jacobian_acceptance_diode.last_current
    @test fallback_acceptance_diode.semiconductor_dissipated_energy_j ≈
        jacobian_acceptance_diode.semiconductor_dissipated_energy_j
    fallback_fidelity = something(fallback_acceptance_diode.extended_fidelity)
    @test fallback_fidelity.recovered_charge.stored_charge_c ≈ 2.0e-6 / 3.0
    @test fallback_fidelity.recovered_charge.recovery_active
    @test fallback_fidelity.thermal.node_temperature_k[1] >=
        fallback_fidelity.thermal.ambient_temperature_k
    @test !fallback_fidelity.candidate_prepared
    @test_throws ArgumentError AIMORA.Nonlinear._power_semiconductor_channel_current_jacobian(
        fallback_acceptance_diode,
        -1.0,
    )
    accepted_boundary_snapshot = sprint(show, fallback_acceptance_diode)
    accepted_boundary_current_a, accepted_boundary_jacobian_s =
        AIMORA.Nonlinear._power_semiconductor_accepted_boundary_channel_current_jacobian(
            fallback_acceptance_diode,
            -1.0,
        )
    @test accepted_boundary_current_a ==
        -fallback_acceptance_diode.on_conductance
    @test accepted_boundary_jacobian_s ==
        fallback_acceptance_diode.on_conductance
    @test sprint(show, fallback_acceptance_diode) == accepted_boundary_snapshot

    triac = AIMORA.Nonlinear.TriacSwitch(
        1,
        0;
        gate_driver=AIMORA.Nonlinear.PowerSemiconductorGateDriver(initially_on=true),
        threshold_v=0.5,
        holding_current=0.1,
    )
    triac.last_voltage = -1.0
    AIMORA.Nonlinear.apply_power_semiconductor_forward_turn_on!(triac, 0.0)
    @test triac.closed
    @test triac.conduction_direction == -1
    triac.last_forward_current_a = 0.05
    AIMORA.Nonlinear.apply_power_semiconductor_forward_extinction!(triac, 1.0e-6)
    @test !triac.closed
    @test triac.conduction_direction == 0

    refusing_gto = AIMORA.Nonlinear.GateTurnOffThyristorSwitch(
        1,
        0;
        gate_driver=AIMORA.Nonlinear.PowerSemiconductorGateDriver(initially_on=true),
        initially_closed=true,
        interruptible_current_a=1.0,
        turn_off_policy=:refuse,
    )
    refusing_gto.last_forward_current_a = 2.0
    refusing_driver_state = sprint(show, refusing_gto.gate_driver)
    @test_throws DomainError AIMORA.Nonlinear.request_power_semiconductor_gate!(
        refusing_gto,
        false,
        0.0,
    )
    @test refusing_gto.closed
    @test sprint(show, refusing_gto.gate_driver) == refusing_driver_state

    remaining_gto = AIMORA.Nonlinear.GateTurnOffThyristorSwitch(
        1,
        0;
        gate_driver=AIMORA.Nonlinear.PowerSemiconductorGateDriver(initially_on=true),
        initially_closed=true,
        interruptible_current_a=1.0,
        turn_off_policy=:remain_latched,
    )
    remaining_gto.last_forward_current_a = 2.0
    @test AIMORA.Nonlinear.request_power_semiconductor_gate!(
        remaining_gto,
        false,
        0.0,
    )
    @test remaining_gto.closed
    @test remaining_gto.gate_turn_off_disposition == :remained_latched
end

@testset "compensated nonlinear-device voltage ownership" begin
    terminal_voltage_v = [1.0000000000000002, 0.0]
    terminal_voltage_correction_v = [-9.997429242503131e-17, 0.0]
    @test terminal_voltage_v[1] + terminal_voltage_correction_v[1] ==
        terminal_voltage_v[1]
    time_s = 1.0e-6
    device = AIMORA.Nonlinear.IGBTSwitch(
        1,
        0;
        initially_closed=true,
        gate_driver=AIMORA.Nonlinear.PowerSemiconductorGateDriver(
            initially_on=true,
        ),
        on_conductance=1.0e12,
        forward_voltage_drop_v=1.0,
        extended_fidelity=AIMORA.Nonlinear.PowerSemiconductorExtendedFidelity(
            turn_off_tail=AIMORA.Nonlinear.TurnOffTailFidelity(1.0e-3),
        ),
    )
    prepare_nonlinear_device_step!(
        device,
        time_s,
        time_s,
        :TrapezoidalCompanion,
    )
    ordinary_current_a = zeros(2)
    ordinary_jacobian_s = zeros(2, 2)
    nonlinear_current_jacobian!(
        ordinary_current_a,
        ordinary_jacobian_s,
        device,
        terminal_voltage_v,
        time_s,
    )
    compensated_current_a = zeros(2)
    compensated_jacobian_s = zeros(2, 2)
    nonlinear_current_jacobian_compensated!(
        compensated_current_a,
        compensated_jacobian_s,
        device,
        terminal_voltage_v,
        terminal_voltage_correction_v,
        time_s,
    )
    expected_current_a = setprecision(BigFloat, 106) do
        Float64(
            BigFloat(device.on_conductance) * (
                BigFloat(terminal_voltage_v[1]) +
                BigFloat(terminal_voltage_correction_v[1]) -
                BigFloat(device.forward_voltage_drop_v)
            ),
        )
    end
    @test compensated_current_a == [expected_current_a, -expected_current_a]
    @test compensated_current_a != ordinary_current_a
    @test compensated_jacobian_s == ordinary_jacobian_s

    zero_correction_current_a = zeros(2)
    zero_correction_jacobian_s = zeros(2, 2)
    nonlinear_current_jacobian_compensated!(
        zero_correction_current_a,
        zero_correction_jacobian_s,
        device,
        terminal_voltage_v,
        zeros(2),
        time_s,
    )
    @test zero_correction_current_a == ordinary_current_a
    @test zero_correction_jacobian_s == ordinary_jacobian_s

    accept_nonlinear_device_compensated_state!(
        device,
        terminal_voltage_v,
        terminal_voltage_correction_v,
        compensated_current_a,
        compensated_jacobian_s,
        time_s,
    )
    finish_nonlinear_device_step!(device)
    @test device.last_current == expected_current_a
    @test device.last_forward_current_a == expected_current_a
    @test device.last_voltage == terminal_voltage_v[1]
    @test device.last_conductance == device.on_conductance
    @test device.last_history_current_a ==
        -device.on_conductance * device.forward_voltage_drop_v
    @test device.last_evaluation_time_s == time_s
    @test something(device.extended_fidelity).previous_terminal_current_a ==
        expected_current_a
    @test !something(device.extended_fidelity).candidate_prepared

    unsupported = AIMORA.NonlinearNetwork.ExponentialCurrentBranch(
        1,
        0;
        saturation_current_a=0.1,
        voltage_scale_v=1.0,
        parallel_conductance_s=1.0,
    )
    @test_throws ArgumentError nonlinear_current_jacobian_compensated!(
        zeros(2),
        zeros(2, 2),
        unsupported,
        [1.0, 0.0],
        [eps(1.0) / 4, 0.0],
        0.0,
    )
    @test_throws ArgumentError accept_nonlinear_device_compensated_state!(
        unsupported,
        [1.0, 0.0],
        [eps(1.0) / 4, 0.0],
        [1.1, -1.1],
        [1.3 -1.3; -1.3 1.3],
        0.0,
    )
end

@testset "private nonlinear ownership and accepted integration" begin
    fidelity = AIMORA.Nonlinear.PowerSemiconductorExtendedFidelity(
        junction_charge=AIMORA.Nonlinear.NonlinearJunctionChargeFidelity(
            1.0e-6,
            10.0,
            0.5;
            voltage_domain_v=(-100.0, 100.0),
        ),
    )
    diode = AIMORA.Nonlinear.DiodeValveSwitch(
        1,
        0;
        on_conductance=10.0,
        off_conductance=0.01,
        extended_fidelity=fidelity,
    )
    linear_system = AIMORA.Nodal.NodalSystem(
        1,
        [AIMORA.Branches.CurrentInjection(1, _time_s -> 1.0)],
    )
    nonlinear_system = AIMORA.NonlinearNodal.NonlinearNodalSystem(
        linear_system,
        [diode];
        scales=NonlinearNetworkScales([1.0], [1.0], Float64[], Float64[]),
    )
    result = AIMORA.NonlinearNodal.advance_nonlinear_step!(
        nonlinear_system,
        1.0e-6,
        1.0e-6,
    )
    @test result.accepted
    @test diode.last_current ≈ 1.0 atol=1.0e-10
    @test !fidelity.candidate_prepared
    accepted_voltage_v = copy(linear_system.v)
    accepted_charge_c = fidelity.junction_charge.previous_charge_c
    # Force the next 0.249 V candidate outside the admitted domain.  A 0.5 V
    # ceiling admits that candidate and therefore cannot exercise rollback.
    fidelity.junction_charge.maximum_voltage_v = 0.2
    rejected_result = AIMORA.NonlinearNodal.advance_nonlinear_step!(
        nonlinear_system,
        2.0e-6,
        1.0e-6,
    )
    @test !rejected_result.accepted
    @test linear_system.v == accepted_voltage_v
    @test fidelity.junction_charge.previous_charge_c == accepted_charge_c
    @test_throws ArgumentError AIMORA.NonlinearNodal.NonlinearNodalSystem(
        AIMORA.Nodal.NodalSystem(1, [diode]),
        [diode];
        scales=NonlinearNetworkScales([1.0], [1.0], Float64[], Float64[]),
    )
end

@testset "localized recovered-charge exhaustion and transactional retry" begin
    recovery = AIMORA.Nonlinear.RecoveredChargeFidelity(
        1.0e-3;
        initial_charge_c=0.5e-3,
    )
    fidelity = AIMORA.Nonlinear.PowerSemiconductorExtendedFidelity(
        recovered_charge=recovery,
    )
    diode = AIMORA.Nonlinear.DiodeValveSwitch(
        1,
        0;
        on_conductance=1.0,
        off_conductance=0.1,
        extended_fidelity=fidelity,
    )
    linear_system = AIMORA.Nodal.NodalSystem(
        1,
        [AIMORA.Branches.CurrentInjection(1, _time_s -> -1.0)],
    )
    nonlinear_system = AIMORA.NonlinearNodal.NonlinearNodalSystem(
        linear_system,
        [diode];
        scales=NonlinearNetworkScales([1.0], [1.0], Float64[], Float64[]),
    )
    result = AIMORA.NonlinearNodal.advance_nonlinear_step!(
        nonlinear_system,
        1.0e-3,
        1.0e-3;
        event_policy=AIMORA.OVER16TimestepIntegration.HybridEventPolicy(
            root_time_tolerance_s=1.0e-10,
            root_value_tolerance=1.0e-12,
            simultaneity_tolerance_s=1.0e-10,
        ),
    )
    @test result.accepted
    @test result.diagnostics.discontinuity_reason == :localized_event
    @test result.diagnostics.accepted_substep_count == 2
    @test recovery.recovery_zero_event_count == 1
    @test recovery.last_recovery_zero_time_s ≈ 0.5e-3 atol=1.0e-8
    @test recovery.stored_charge_c == 0.0
    @test fidelity.candidate_prepared == false

    checkpoint = AIMORA.NonlinearNodal.nonlinear_nodal_checkpoint(nonlinear_system)
    accepted_energy_j = diode.semiconductor_dissipated_energy_j
    accepted_event_count = recovery.recovery_zero_event_count
    retry = AIMORA.NonlinearNodal.advance_nonlinear_step!(
        nonlinear_system,
        2.0e-3,
        1.0e-3,
    )
    @test retry.accepted
    @test recovery.recovery_zero_event_count == accepted_event_count
    AIMORA.NonlinearNodal.restore_nonlinear_nodal_checkpoint!(
        nonlinear_system,
        checkpoint,
    )
    replay = AIMORA.NonlinearNodal.advance_nonlinear_step!(
        nonlinear_system,
        2.0e-3,
        1.0e-3,
    )
    @test replay.accepted
    @test recovery.recovery_zero_event_count == accepted_event_count
    @test diode.semiconductor_dissipated_energy_j >= accepted_energy_j
end

@testset "semiconductor boundary-dependent state reconstruction" begin
    tail = AIMORA.Nonlinear.TurnOffTailFidelity(
        2.0e-3;
        cutoff_current_a=1.0e-6,
    )
    tail.active = true
    tail.current_a = 99.0
    tail.initial_current_a = 4.0
    tail.turn_off_time_s = 0.25e-3
    tail.last_duration_s = 99.0
    recovery = AIMORA.Nonlinear.RecoveredChargeFidelity(
        1.0e-3;
        initial_charge_c=1.0e-3,
    )
    recovery.recovery_active = true
    recovery.recovery_start_time_s = 0.5e-3
    recovery.last_recovery_duration_s = 99.0
    device = AIMORA.Nonlinear.IGBTSwitch(
        1,
        0;
        antiparallel_diode=AIMORA.Nonlinear.AntiparallelDiodeParameters(),
        extended_fidelity=AIMORA.Nonlinear.PowerSemiconductorExtendedFidelity(
            recovered_charge=recovery,
            turn_off_tail=tail,
        ),
    )
    before = AIMORA.Nonlinear.power_semiconductor_boundary_state_errors(
        device,
        1.25e-3,
    )
    @test before.tail_current > 0.0
    @test before.tail_duration > 0.0
    @test before.recovery_duration > 0.0
    AIMORA.Nonlinear.reconstruct_power_semiconductor_boundary_state!(
        device,
        1.25e-3,
    )
    @test tail.last_duration_s == 1.0e-3
    @test tail.current_a == 4.0 * exp(-0.5)
    @test recovery.last_recovery_duration_s == 0.75e-3
    @test all(iszero, values(
        AIMORA.Nonlinear.power_semiconductor_boundary_state_errors(
            device,
            1.25e-3,
        ),
    ))
    tail.active = false
    tail.current_a = 3.0
    retained_duration_s = tail.last_duration_s
    AIMORA.Nonlinear.reconstruct_power_semiconductor_boundary_state!(
        device,
        1.25e-3,
    )
    @test tail.current_a == 0.0
    @test tail.last_duration_s == retained_duration_s

    junction = AIMORA.Nonlinear.NonlinearJunctionChargeFidelity(
        1.0e-6,
        10.0,
        0.5;
        voltage_domain_v=(-100.0, 100.0),
    )
    junction.last_displacement_current_a = 0.25
    thermal = AIMORA.Nonlinear.CauerThermalFidelity(
        [1.0, 2.0],
        [0.5, 1.0];
        initial_temperature_k=[310.0, 305.0],
    )
    complete_fidelity = AIMORA.Nonlinear.PowerSemiconductorExtendedFidelity(
        junction_charge=junction,
        thermal=thermal,
    )
    complete_fidelity.pending_event_current_a = 17.0
    complete_fidelity.pending_event_blocking_voltage_v = 700.0
    complete_device = AIMORA.Nonlinear.IGBTSwitch(
        1,
        0;
        initially_closed=true,
        gate_driver=AIMORA.Nonlinear.PowerSemiconductorGateDriver(
            initially_on=true,
        ),
        on_conductance=2.0,
        forward_voltage_drop_v=1.5,
        snubber=AIMORA.Nonlinear.SeriesRCSnubber(10.0, 1.0e-6),
        extended_fidelity=complete_fidelity,
    )
    complete_device.snubber.capacitor_voltage_v = 4.0
    complete_device.last_conductance = 99.0
    complete_device.last_history_current_a = 99.0
    complete_device.last_evaluation_time_s = -1.0
    complete_device.last_voltage = -99.0
    junction.previous_voltage_v = -99.0
    junction.previous_charge_c = -99.0
    junction.last_capacitance_f = 99.0
    junction.last_charge_c = 99.0
    complete_device.snubber.previous_current_a = 99.0
    complete_device.snubber.last_branch_voltage_v = -99.0
    complete_device.snubber.last_current_a = 99.0
    complete_device.snubber.last_resistor_loss_w = 99.0
    thermal.trial_temperature_k .= 99.0
    thermal.last_stored_energy_j = 99.0
    thermal.last_ambient_heat_flow_w = 99.0
    complete_before =
        AIMORA.Nonlinear.power_semiconductor_complete_boundary_state_errors(
            complete_device,
            1.25e-3,
            10.0,
        )
    @test all(>(0.0), values(complete_before))
    AIMORA.Nonlinear.reconstruct_power_semiconductor_boundary_state!(
        complete_device,
        1.25e-3,
        10.0,
    )
    @test all(iszero, values(
        AIMORA.Nonlinear.power_semiconductor_complete_boundary_state_errors(
            complete_device,
            1.25e-3,
            10.0,
        ),
    ))
    @test complete_fidelity.pending_event_current_a == 17.0
    @test complete_fidelity.pending_event_blocking_voltage_v == 700.0
end

@testset "gate timing has deterministic physical constraint ownership" begin
    nonlinear = AIMORA.Nonlinear
    driver = nonlinear.PowerSemiconductorGateDriver(
        turn_on_delay_s=2.0,
        turn_off_delay_s=2.0,
        dead_time_s=1.0,
        minimum_pulse_width_s=1.0,
    )
    driver.last_turn_off_time_s = 2.0
    tied_turn_on = nonlinear.power_semiconductor_gate_transition_timing(
        driver,
        true,
        1.0;
        earliest_transition_time_s=3.0,
    )
    @test tied_turn_on.transition_time_s == 3.0
    @test tied_turn_on.active_constraint === :turn_on_delay
    replayed_tie = nonlinear.power_semiconductor_gate_transition_timing(
        deepcopy(driver),
        true,
        1.0;
        earliest_transition_time_s=3.0,
    )
    @test replayed_tie.transition_time_s == tied_turn_on.transition_time_s
    @test replayed_tie.active_constraint === tied_turn_on.active_constraint

    driver.last_turn_off_time_s = 4.0
    prior_off_time = nonlinear.power_semiconductor_gate_transition_timing(
        driver,
        true,
        1.0,
    )
    @test prior_off_time.transition_time_s == 5.0
    @test prior_off_time.active_constraint === :minimum_off_time
    driver.last_turn_off_time_s = 2.0
    interlocked = nonlinear.power_semiconductor_gate_transition_timing(
        driver,
        true,
        1.0;
        earliest_transition_time_s=4.0,
    )
    @test interlocked.transition_time_s == 4.0
    @test interlocked.active_constraint === :external_interlock

    driver.last_turn_on_time_s = 2.0
    tied_turn_off = nonlinear.power_semiconductor_gate_transition_timing(
        driver,
        false,
        1.0,
    )
    @test tied_turn_off.transition_time_s == 3.0
    @test tied_turn_off.active_constraint === :turn_off_delay

    requested_driver = deepcopy(driver)
    requested_driver.last_turn_off_time_s = 2.0
    requested_device = nonlinear.IGBTSwitch(
        1,
        0;
        gate_driver=requested_driver,
    )
    @test !nonlinear.request_power_semiconductor_gate!(
        requested_device,
        true,
        1.0;
        earliest_transition_time_s=3.0,
    )
    @test requested_driver.pending_transition_time_s ==
        tied_turn_on.transition_time_s

    requested_driver.pending_state = false
    requested_driver.pending_transition_time_s = 5.0
    tied_off_boundary = nonlinear.power_semiconductor_gate_off_boundary_timing(
        requested_driver,
        5.0,
    )
    @test tied_off_boundary.transition_time_s == 5.0
    @test tied_off_boundary.active_constraint === :gate_command_time
    requested_driver.pending_transition_time_s = 6.0
    pending_off_boundary =
        nonlinear.power_semiconductor_gate_off_boundary_timing(
            requested_driver,
            5.0,
        )
    @test pending_off_boundary.transition_time_s == 6.0
    @test pending_off_boundary.active_constraint ===
        :pending_turn_off_transition
    requested_driver.pending_state = nothing
    requested_driver.last_turn_off_time_s = 7.0
    completed_off_boundary =
        nonlinear.power_semiconductor_gate_off_boundary_timing(
            requested_driver,
            5.0,
        )
    @test completed_off_boundary.transition_time_s == 7.0
    @test completed_off_boundary.active_constraint ===
        :last_turn_off_transition
end

@testset "analytic tail cutoff localization" begin
    tail = AIMORA.Nonlinear.TurnOffTailFidelity(
        1.0e-3;
        cutoff_current_a=exp(-1.0),
    )
    tail.active = true
    tail.current_a = 1.0
    tail.initial_current_a = 1.0
    tail.turn_off_time_s = 0.0
    fidelity = AIMORA.Nonlinear.PowerSemiconductorExtendedFidelity(
        turn_off_tail=tail,
    )
    igbt = AIMORA.Nonlinear.IGBTSwitch(
        1,
        0;
        gate_driver=AIMORA.Nonlinear.PowerSemiconductorGateDriver(),
        off_conductance=0.2,
        extended_fidelity=fidelity,
    )
    system = AIMORA.NonlinearNodal.NonlinearNodalSystem(
        AIMORA.Nodal.NodalSystem(
            1,
            [AIMORA.Branches.CurrentInjection(1, _time_s -> 1.0)],
        ),
        [igbt];
        scales=NonlinearNetworkScales([10.0], [1.0], Float64[], Float64[]),
    )
    result = AIMORA.NonlinearNodal.advance_nonlinear_step!(system, 2.0e-3, 2.0e-3)
    @test result.accepted
    @test result.diagnostics.discontinuity_reason == :localized_event
    @test tail.cutoff_event_count == 1
    @test tail.last_cutoff_time_s ≈ 1.0e-3 atol=1.0e-12
    @test !tail.active
    @test tail.current_a == 0.0

    rounded_tail = AIMORA.Nonlinear.TurnOffTailFidelity(
        1.0e-3;
        cutoff_current_a=exp(-1.0),
    )
    rounded_tail.active = true
    rounded_tail.current_a = 1.0
    rounded_tail.initial_current_a = 1.0
    rounded_tail.turn_off_time_s = 0.0
    rounded_igbt = AIMORA.Nonlinear.IGBTSwitch(
        1,
        0;
        extended_fidelity=AIMORA.Nonlinear.PowerSemiconductorExtendedFidelity(
            turn_off_tail=rounded_tail,
        ),
    )
    AIMORA.Nonlinear.apply_power_semiconductor_tail_cutoff!(
        rounded_igbt,
        1.0e-3 + 0.5e-12,
    )
    @test rounded_tail.cutoff_event_count == 1
    @test !rounded_tail.active
end

@testset "overlapping turn-off tail current is conserved at commutation" begin
    tail = AIMORA.Nonlinear.TurnOffTailFidelity(
        2.0e-3;
        cutoff_current_a=1.0e-6,
    )
    tail.active = true
    tail.initial_current_a = 4.0
    tail.current_a = 4.0
    tail.turn_off_time_s = 0.0
    device = AIMORA.Nonlinear.IGBTSwitch(
        1,
        0;
        extended_fidelity=AIMORA.Nonlinear.PowerSemiconductorExtendedFidelity(
            turn_off_tail=tail,
        ),
        initially_closed=true,
        gate_driver=AIMORA.Nonlinear.PowerSemiconductorGateDriver(
            initially_on=true,
        ),
    )
    device.last_forward_current_a = 0.5
    first_surviving_a = 4.0 * exp(-0.5)
    AIMORA.Nonlinear.apply_power_semiconductor_forward_extinction!(device, 1.0e-3)
    @test tail.active
    @test tail.initial_current_a ≈ first_surviving_a + 0.5 atol=1.0e-12
    @test tail.current_a == tail.initial_current_a
    @test tail.turn_off_time_s == 1.0e-3
    AIMORA.Nonlinear.apply_power_semiconductor_forward_turn_on!(device, 1.0e-3)
    device.last_forward_current_a = 0.25
    second_surviving_a = tail.initial_current_a * exp(-0.5)
    AIMORA.Nonlinear.apply_power_semiconductor_forward_extinction!(device, 2.0e-3)
    @test tail.initial_current_a ≈ second_surviving_a + 0.25 atol=1.0e-12
    @test tail.current_a == tail.initial_current_a
    @test tail.turn_off_time_s == 2.0e-3
end

@testset "physical event energy rollback and single retry deposition" begin
    energy_table = synthetic_semiconductor_energy_table()
    fidelity = AIMORA.Nonlinear.PowerSemiconductorExtendedFidelity(
        switching_energy=energy_table,
        thermal=AIMORA.Nonlinear.CauerThermalFidelity(
            [1.0],
            [1.0];
            initial_temperature_k=[300.0],
        ),
    )
    igbt = AIMORA.Nonlinear.IGBTSwitch(
        1,
        0;
        gate_driver=AIMORA.Nonlinear.PowerSemiconductorGateDriver(
            initially_on=true,
        ),
        initially_closed=true,
        on_conductance=1.0,
        off_conductance=0.1,
        extended_fidelity=fidelity,
    )
    AIMORA.Nonlinear.power_semiconductor_event_localization!(igbt)
    igbt.last_voltage = 1.0
    igbt.last_forward_current_a = 1.0
    AIMORA.Nonlinear.request_power_semiconductor_gate!(
        igbt,
        false,
        0.0;
        earliest_transition_time_s=0.4,
    )
    failure_enabled = Ref(true)
    failure = ToggleableFailureCurrentBranch(1, 0, failure_enabled, 0.75, 0)
    linear_system = AIMORA.Nodal.NodalSystem(
        1,
        [AIMORA.Branches.CurrentInjection(1, _time_s -> 1.0)],
    )
    system = AIMORA.NonlinearNodal.NonlinearNodalSystem(
        linear_system,
        [igbt, failure];
        scales=NonlinearNetworkScales([10.0], [1.0], Float64[], Float64[]),
    )
    baseline_voltage_v = copy(linear_system.v)
    baseline_temperature_k = copy(fidelity.thermal.node_temperature_k)
    baseline_energy_j = igbt.semiconductor_dissipated_energy_j
    baseline_transition_count = igbt.topology_transition_count
    baseline_driver_transition_count = igbt.gate_driver.transition_count

    rejected = AIMORA.NonlinearNodal.advance_nonlinear_step!(system, 1.0, 1.0)
    @test !rejected.accepted
    @test rejected.failure.code == :nonfinite_device_current
    @test linear_system.v == baseline_voltage_v
    @test igbt.closed
    @test igbt.gate_driver.pending_state === false
    @test igbt.gate_driver.transition_count == baseline_driver_transition_count
    @test igbt.topology_transition_count == baseline_transition_count
    @test energy_table.cumulative_turn_off_energy_j == 0.0
    @test energy_table.last_event_kind == :none
    @test fidelity.thermal.node_temperature_k == baseline_temperature_k
    @test igbt.semiconductor_dissipated_energy_j == baseline_energy_j
    @test failure.accepted_state_count == 0

    failure_enabled[] = false
    accepted = AIMORA.NonlinearNodal.advance_nonlinear_step!(system, 1.0, 1.0)
    @test accepted.accepted
    @test accepted.diagnostics.discontinuity_reason == :localized_event
    @test !igbt.closed
    @test igbt.gate_driver.pending_state === nothing
    @test igbt.gate_driver.transition_count == baseline_driver_transition_count + 1
    @test igbt.topology_transition_count == baseline_transition_count + 1
    expected_event_energy_j = energy_table.last_event_energy_j
    @test expected_event_energy_j > 0.0
    @test energy_table.cumulative_turn_off_energy_j == expected_event_energy_j
    @test energy_table.last_event_kind == :turn_off
    @test energy_table.last_event_energy_j == expected_event_energy_j
    @test igbt.semiconductor_dissipated_energy_j >= expected_event_energy_j
    @test isfinite(fidelity.thermal.node_temperature_k[1])
    @test fidelity.thermal.node_temperature_k[1] != baseline_temperature_k[1]
end

@testset "ordered power-semiconductor topology-transition observations" begin
    nonlinear = AIMORA.Nonlinear
    device = nonlinear.IGBTSwitch(
        1,
        0;
        gate_driver=nonlinear.PowerSemiconductorGateDriver(
            initially_on=true,
        ),
        antiparallel_diode=nonlinear.AntiparallelDiodeParameters(),
        initially_closed=true,
    )
    device.last_voltage = -1.0
    device.last_forward_current_a = 1.0
    initial = nonlinear.power_semiconductor_topology_chart_state(device)
    journal = Any[]

    @test nonlinear.request_power_semiconductor_gate!(device, false, 1.0)
    @test nonlinear.record_power_semiconductor_topology_transition_observation!(
        journal,
        device,
        initial,
        :gate_release,
    )
    after_gate = nonlinear.power_semiconductor_topology_chart_state(device)
    nonlinear.apply_power_semiconductor_reverse_turn_on!(device, 1.0)
    @test nonlinear.record_power_semiconductor_topology_transition_observation!(
        journal,
        device,
        after_gate,
        :conduction_selection,
    )

    @test length(journal) == 2
    @test getfield.(journal, :transition_owner) ==
        [:gate_release, :conduction_selection]
    @test journal[1].before == initial
    @test journal[1].after == journal[2].before == after_gate
    @test journal[2].after ==
        nonlinear.power_semiconductor_topology_chart_state(device)
    @test journal[1].after.topology_transition_count ==
        initial.topology_transition_count + 1
    @test journal[2].after.topology_transition_count ==
        initial.topology_transition_count + 2
    @test !journal[2].after.forward_conducting
    @test journal[2].after.reverse_diode_conducting

    hidden = deepcopy(device)
    hidden_initial = nonlinear.power_semiconductor_topology_chart_state(hidden)
    nonlinear.apply_power_semiconductor_reverse_extinction!(hidden, 2.0)
    nonlinear.apply_power_semiconductor_forward_turn_on!(hidden, 2.0)
    @test_throws ArgumentError begin
        nonlinear.record_power_semiconductor_topology_transition_observation!(
            Any[],
            hidden,
            hidden_initial,
            :hidden_multiple_transitions,
        )
    end
    @test_throws ArgumentError begin
        nonlinear.record_power_semiconductor_topology_transition_observation!(
            Any[],
            device,
            nonlinear.power_semiconductor_topology_chart_state(device),
            :none,
        )
    end
end

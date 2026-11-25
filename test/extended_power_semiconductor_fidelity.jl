using AIMORA.NonlinearNetwork

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


@testset "Declared semiconductor schema version" begin
    device_models = AIMORA.Nonlinear
    provenance = device_models.NonlinearJunctionChargeFidelity(200e-9, 50.0, 0.45).provenance
    source_hash = repeat("a", 64)
    function identity(version)
        return device_models.DeclaredSemiconductorModelIdentity(
            :aimora_generic_semiconductor_v1, version, "synthetic schema fixture",
            source_hash, "CC0-1.0", :public, provenance,
        )
    end
    fidelity = device_models.PowerSemiconductorExtendedFidelity(declared_model=identity(v"1.0.0"))
    diode = device_models.DiodeValveSwitch(1, 0; extended_fidelity=fidelity)
    @test :declared_compact_model in device_models.power_semiconductor_extended_state(diode).fidelity_components
    @test fidelity.declared_model.version == v"1.0.0"
    @test_throws ArgumentError identity(v"0.9.0")
    @test_throws ArgumentError identity(v"1.1.0")
    @test_throws ArgumentError identity(v"9.0.0")
end

@testset "Inclusive signed terminal current at accepted endpoints" begin
    models = AIMORA.Nonlinear
    domain = models.SemiconductorTerminalCurrentDomain
    @test_throws ArgumentError domain(-Inf, 3)
    @test_throws ArgumentError domain(-2, NaN)
    @test_throws ArgumentError domain(3, -2)

    function prepared_device(minimum_a, maximum_a; snubber=nothing, nodes=(1, 0))
        fidelity = models.PowerSemiconductorExtendedFidelity(
            recovered_charge=models.RecoveredChargeFidelity(1e-3),
            terminal_current_domain=domain(minimum_a, maximum_a),
        )
        device = models.DiodeValveSwitch(
            nodes...; on_conductance=1.0, off_conductance=1.0, snubber,
            extended_fidelity=fidelity,
        )
        prepare_nonlinear_device_step!(device, 1.0, 1.0, :BackwardEulerCompanion)
        return device
    end

    for i in (-2.0, 0.0, 3.0, prevfloat(-2.0), nextfloat(3.0))
        device = prepared_device(-2, 3)
        voltages = [i, 0.0]
        currents = zeros(2)
        jacobian = zeros(2, 2)
        nonlinear_current_jacobian!(currents, jacobian, device, voltages, 1.0)
        @test currents == [i, -i]
        if -2 <= i <= 3
            accept_nonlinear_device_state!(device, voltages, currents, jacobian, 1.0)
            @test device.last_current == i
        else
            @test_throws DomainError accept_nonlinear_device_state!(device, voltages, currents, 1.0)
            @test device.last_current == 0.0
            @test device.extended_fidelity.previous_terminal_current_a == 0.0
            @test device.extended_fidelity.candidate_prepared
        end
    end

    snubber = models.SeriesRCSnubber(1, 1; initial_capacitor_voltage_v=-3)
    cancelling = prepared_device(0, 0; snubber, nodes=(0, 1))
    currents, jacobian = zeros(2), zeros(2, 2)
    nonlinear_current_jacobian!(currents, jacobian, cancelling, [0.0, 1.0], 1.0)
    @test currents == [0.0, 0.0]
    accept_nonlinear_device_state!(cancelling, [0.0, 1.0], currents, jacobian, 1.0)
    @test cancelling.last_snubber_current_a ≈ 1.0
    @test cancelling.last_current == 0.0

    adding = prepared_device(-1.25, 1.25; snubber=models.SeriesRCSnubber(1, 1))
    nonlinear_current_jacobian!(currents, jacobian, adding, [-1.0, 0.0], 1.0)
    @test currents[1] ≈ -1.5
    @test_throws DomainError accept_nonlinear_device_state!(adding, [-1.0, 0.0], currents, 1.0)
    @test adding.snubber.last_current_a == 0.0
    @test_throws DomainError accept_nonlinear_device_state!(adding, [-1.0, 0.0], zeros(2), 1.0)
    @test adding.last_current == 0.0
end

@testset "Recovery event probes and transactional terminal-current refusal" begin
    models, branches, solver = AIMORA.Nonlinear, AIMORA.Branches, AIMORA.NonlinearNodal

    function recovery_network(minimum_a, off_conductance_s)
        stored_charge = models.RecoveredChargeFidelity(1e-3; initial_charge_c=0.5e-3)
        fidelity = models.PowerSemiconductorExtendedFidelity(
            recovered_charge=stored_charge,
            terminal_current_domain=models.SemiconductorTerminalCurrentDomain(minimum_a, 0.1),
        )
        valve = models.DiodeValveSwitch(
            1, 0; on_conductance=1.0, off_conductance=off_conductance_s,
            extended_fidelity=fidelity,
        )
        upstream_valve = models.DiodeValveSwitch(
            2, 0; off_conductance=0.25,
            extended_fidelity=models.PowerSemiconductorExtendedFidelity(
                recovered_charge=models.RecoveredChargeFidelity(1e-3),
            ),
        )
        series_branch = branches.SeriesRLCBranch(2, 0, 1.0, 1e-3, 1e-3)
        circuit = AIMORA.Nodal.NodalSystem(2, [
            branches.CurrentInjection(1, t -> -4000t),
            branches.ConductanceBranch(1, 0, 1.0),
            branches.CurrentInjection(2, _ -> 1.0),
            series_branch,
        ])
        circuit_owner = solver.NonlinearNodalSystem(
            circuit, [upstream_valve, valve];
            scales=NonlinearNetworkScales([1.0, 1.0], [1.0, 1.0], Float64[], Float64[]),
        )
        return circuit_owner, valve, stored_charge, upstream_valve, series_branch
    end

    network, diode, recovery, _, _ = recovery_network(-1.2, 0.1)
    observations = Tuple{Float64,Float64}[]
    observe = (_, time, _step, _method, _index) -> push!(observations, (time, diode.last_current))
    pending = solver._pending_nonlinear_device_events(network, 0.0)
    _, probe = solver._probe_right_nonlinear_device_event_values!(
        network, pending, 0.0, 1e-3, network.topology_signature, network.previous_topology_signature,
    )
    @test probe.accepted
    @test probe.voltage_v[1] ≈ -2.0 atol=1e-9
    @test probe.voltage_v[1] * diode.on_conductance < -1.2
    @test diode.last_current == 0.0
    @test recovery.stored_charge_c == 0.5e-3
    result = solver.advance_nonlinear_step!(
        network, 1e-3, 1e-3; accepted_substep_observer=observe,
    )
    # The uncommitted right-hand probe in recovery sees -4 A divided equally
    # by the 1 S diode and external 1 S shunt. The event side is -1 A.
    @test result.accepted
    @test length(observations) == 2
    @test recovery.recovery_zero_event_count == 1
    @test observations[1][1] ≈ 0.5e-3 atol=1e-8
    @test observations[1][2] ≈ -1.0 atol=1e-5
    @test diode.last_current ≈ -4 / 11 atol=1e-9
    @test observations[end][2] == diode.last_current

    rejected_network, rejected_diode, rejected_recovery, earlier, storage = recovery_network(-1.2, 0.9)
    before = (
        voltage=copy(solver.nonlinear_linear_system(rejected_network).v),
        charge=rejected_recovery.stored_charge_c,
        current=rejected_diode.last_current,
        work=rejected_diode.semiconductor_dissipated_energy_j,
        events=rejected_recovery.recovery_zero_event_count,
        earlier=(earlier.last_current, earlier.semiconductor_dissipated_energy_j),
        storage=(storage.i_prev, storage.v_prev, storage.i_last),
    )
    observed = Ref(0)
    rejected = solver.advance_nonlinear_step!(
        rejected_network, 1e-3, 1e-3;
        accepted_substep_observer=(_system, _time, _step, _method, _index) -> (observed[] += 1),
    )
    @test !rejected.accepted
    @test occursin("terminal current from node 1 to 0", sprint(showerror, rejected.failure))
    @test observed[] == 0
    @test solver.nonlinear_linear_system(rejected_network).v == before.voltage
    @test rejected_recovery.stored_charge_c == before.charge
    @test rejected_diode.last_current == before.current
    @test rejected_diode.semiconductor_dissipated_energy_j == before.work
    @test rejected_recovery.recovery_zero_event_count == before.events
    @test !rejected_diode.extended_fidelity.candidate_prepared
    @test (earlier.last_current, earlier.semiconductor_dissipated_energy_j) == before.earlier
    @test (storage.i_prev, storage.v_prev, storage.i_last) == before.storage
    @test solver.advance_nonlinear_step!(rejected_network, 0.4e-3, 0.4e-3).accepted
    @test rejected_recovery.recovery_zero_event_count == 0
    @test earlier.last_current != before.earlier[1]
    @test storage.i_last != before.storage[3]

    checkpoint = solver.nonlinear_nodal_checkpoint(network)
    saved_current = diode.last_current
    @test saved_current ≈ -4 / 11 atol=1e-9
    @test solver.advance_nonlinear_step!(network, 1.05e-3, 0.05e-3).accepted
    continued_current = diode.last_current
    solver.restore_nonlinear_nodal_checkpoint!(network, checkpoint)
    @test diode.last_current == saved_current
    @test solver.advance_nonlinear_step!(network, 1.05e-3, 0.05e-3).accepted
    @test diode.last_current == continued_current
    diode.extended_fidelity.terminal_current_domain = nothing
    @test_throws ArgumentError solver.restore_nonlinear_nodal_checkpoint!(network, checkpoint)
    @test diode.extended_fidelity.terminal_current_domain === nothing
    @test diode.last_current == continued_current
    @info "Terminal-current domain recovery probe" held_probe_a=probe.voltage_v[1] * diode.on_conductance event_time_s=observations[1][1] event_a=observations[1][2] accepted_endpoint_a=saved_current refused=sprint(showerror, rejected.failure)
end

@testset "Signed depletion current limits preserve admitted charge history" begin
    models = AIMORA.Nonlinear

    function depletion_endpoint(path, domain)
        junction = models.NonlinearJunctionChargeFidelity(
            40e-6, 50, 0.45; initial_voltage_v=first(path), voltage_domain_v=(-50, -10),
        )
        fidelity = models.PowerSemiconductorExtendedFidelity(
            junction_charge=junction, terminal_current_domain=domain,
        )
        device = models.DiodeValveSwitch(1, 0; off_conductance=0, extended_fidelity=fidelity)
        currents, jacobian = zeros(2), zeros(2, 2)
        old_charge = junction.previous_charge_c
        prepare_nonlinear_device_step!(device, 1e-9, 1e-9, :BackwardEulerCompanion)
        voltage = [last(path), 0.0]
        nonlinear_current_jacobian!(currents, jacobian, device, voltage, 1e-9)
        return device, junction, old_charge, voltage, currents, jacobian
    end

    bounds = models.SemiconductorTerminalCurrentDomain(-10_000, 10_000)
    for sign in (-1, 1)
        path = sign > 0 ? (-10.2, -10.0) : (-10.0, -10.2)
        switch, charge, initial_q, voltage, currents, jacobian = depletion_endpoint(path, bounds)
        @test sign * currents[1] ≈ 7364.327728395467 atol=1e-6
        @test abs(currents[1]) < 10_000
        accept_nonlinear_device_state!(switch, voltage, currents, jacobian, 1e-9)
        finish_nonlinear_device_step!(switch)
        @test sign * (charge.previous_charge_c - initial_q) > 0
    end
    for sign in (-1, 1), domain in (nothing, bounds)
        path = sign > 0 ? (-10.4, -10.0) : (-10.0, -10.4)
        switch, charge, initial_q, voltage, currents, jacobian = depletion_endpoint(path, domain)
        @test sign * currents[1] ≈ 14717.65386554496 atol=1e-6
        @test abs(currents[1]) > 10_000
        if domain === nothing
            accept_nonlinear_device_state!(switch, voltage, currents, jacobian, 1e-9)
            @test sign * (charge.previous_charge_c - initial_q) > 0
        else
            @test_throws DomainError accept_nonlinear_device_state!(switch, voltage, currents, jacobian, 1e-9)
            @test charge.previous_charge_c == initial_q
            @test switch.last_current == 0.0
        end
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
        # The baseline retains its established endpoint-average loss in both methods.
        simple = AIMORA.Nonlinear.DiodeValveSwitch(1,0;on_conductance=1.0,off_conductance=1.0)
        simple.last_conductance = 1.0
        simple.previous_semiconductor_loss_w = 9.0
        integration = method === :BackwardEulerCompanion ?
            Val(AIMORA.Branches.BackwardEulerCompanion) : Val(AIMORA.Branches.TrapezoidalCompanion)
        AIMORA.Branches.update!(simple,[2.0],step,integration)
        @test simple.semiconductor_dissipated_energy_j ≈ 0.65 atol=1e-15 rtol=0
    end
end

@testset "Backward-Euler signed channel work is separate from heat" begin
    models = AIMORA.Nonlinear
    fidelity = models.PowerSemiconductorExtendedFidelity(
        thermal=models.CauerThermalFidelity([1.0], [1.0]),
    )
    device = models.DiodeValveSwitch(
        1, 0; on_conductance=1.0, initially_closed=true,
        forward_voltage_drop_v=3.0, extended_fidelity=fidelity,
    )
    models.power_semiconductor_event_localization!(device)
    voltage = [2.0, 0.0]
    current, jacobian = zeros(2), zeros(2, 2)
    prepare_nonlinear_device_step!(device, 0.1, 0.1, :BackwardEulerCompanion)
    nonlinear_current_jacobian!(current, jacobian, device, voltage, 0.1)
    @test current == [-1.0, 1.0]
    accept_nonlinear_device_state!(device, voltage, current, jacobian, 0.1)
    @test fidelity.companion_energy_residual_j ≈ -0.2 atol=1e-15
    @test device.semiconductor_dissipated_energy_j == 0.0
    @test fidelity.thermal.cumulative_input_energy_j == 0.0
    finish_nonlinear_device_step!(device)
end

@testset "Backward-Euler channel work excludes the owned RC snubber" begin
    models = AIMORA.Nonlinear
    snubber = models.SeriesRCSnubber(1.0, 1.0)
    fidelity = models.PowerSemiconductorExtendedFidelity(
        thermal=models.CauerThermalFidelity([1.0], [1.0]),
    )
    device = models.DiodeValveSwitch(
        1, 0; off_conductance=1.0, snubber, extended_fidelity=fidelity,
    )
    voltage = [2.0, 0.0]
    current, jacobian = zeros(2), zeros(2, 2)
    prepare_nonlinear_device_step!(device, 0.1, 0.1, :BackwardEulerCompanion)
    nonlinear_current_jacobian!(current, jacobian, device, voltage, 0.1)
    snubber_current = 2.0 / 1.1
    @test current[1] ≈ 2.0 + snubber_current atol=1e-15
    accept_nonlinear_device_state!(device, voltage, current, jacobian, 0.1)
    @test device.semiconductor_dissipated_energy_j ≈ 0.4 atol=1e-15
    @test fidelity.thermal.cumulative_input_energy_j ≈ 0.4 atol=1e-15
    @test snubber.capacitor_voltage_v ≈ 0.1 * snubber_current atol=1e-15
    @test snubber.dissipated_energy_j ≈ 0.05 * snubber_current^2 atol=1e-15
    expected_snubber_residual = 0.05 * 2.0 * snubber_current -
        0.05 * snubber_current^2 - 0.5 * (0.1 * snubber_current)^2
    @test fidelity.companion_energy_residual_j ≈ expected_snubber_residual atol=1e-15
    finish_nonlinear_device_step!(device)
end

const ModernMachineTest = AIMORA.ModernMachines

function modern_machine_test_electrical(family)
    synchronous = family in (
        ModernMachineTest.WoundFieldSynchronousMachine,
        ModernMachineTest.SynchronousCondenserMachine,
    )
    permanent = family === ModernMachineTest.PermanentMagnetSynchronousMachine
    return ModernMachineTest.MachineElectricalParameters(
        stator_resistance_ohm=0.18,
        zero_sequence_inductance_h=0.012,
        stator_d_leakage_inductance_h=0.006,
        stator_q_leakage_inductance_h=0.007,
        d_axis_magnetizing_inductance_h=0.085,
        q_axis_magnetizing_inductance_h=0.072,
        field_resistance_ohm=synchronous ? 1.2 : 0.0,
        field_leakage_inductance_h=synchronous ? 0.018 : 0.0,
        d_damper_resistance_ohm=synchronous ? 0.45 : 0.0,
        d_damper_leakage_inductance_h=synchronous ? 0.011 : 0.0,
        q_damper_resistance_ohm=synchronous ? 0.52 : 0.0,
        q_damper_leakage_inductance_h=synchronous ? 0.013 : 0.0,
        permanent_magnet_flux_wb=permanent ? 0.48 : 0.0,
    )
end

function modern_machine_test_rotor_branches(family; branch_count=1)
    induction = family in (
        ModernMachineTest.CageInductionMachine,
        ModernMachineTest.WoundRotorInductionMachine,
        ModernMachineTest.DoublyFedInductionMachine,
    )
    induction || return ModernMachineTest.MachineRotorBranch[]
    exposed_family = family in (
        ModernMachineTest.WoundRotorInductionMachine,
        ModernMachineTest.DoublyFedInductionMachine,
    )
    return ModernMachineTest.MachineRotorBranch[
        ModernMachineTest.MachineRotorBranch(
            Symbol("rotor_branch_$(index)");
            resistance_ohm=0.12 + 0.05 * index,
            leakage_inductance_h=0.004 + 0.002 * index,
            terminal_exposed=exposed_family && index == 1,
        ) for index in 1:branch_count
    ]
end

function modern_machine_test_specification(
    family;
    id=Symbol(lowercase(string(family))),
    timestep_s=1.0e-5,
    branch_count=1,
    mass_count=1,
    controls=false,
)
    synchronous_speed = 2.0 * pi * 60.0 / 2.0
    masses = ModernMachineTest.MachineShaftMass[
        ModernMachineTest.MachineShaftMass(
            Symbol("shaft_mass_$(index)");
            inertia_kg_m2=1.5 + 0.25 * index,
            damping_nm_s_per_rad=0.002 * index,
            initial_angle_rad=0.0,
            initial_speed_rad_s=synchronous_speed,
        ) for index in 1:mass_count
    ]
    couplings = ModernMachineTest.MachineShaftCoupling[
        ModernMachineTest.MachineShaftCoupling(
            Symbol("shaft_section_$(index)"),
            masses[index].id,
            masses[index + 1].id;
            stiffness_nm_per_rad=2.0e4 + 1.0e3 * index,
            damping_nm_s_per_rad=2.0 + index,
        ) for index in 1:(mass_count - 1)
    ]
    operating_mode = family === ModernMachineTest.SynchronousCondenserMachine ?
        ModernMachineTest.MachineCondenserMode :
        family in (
            ModernMachineTest.WoundFieldSynchronousMachine,
            ModernMachineTest.PermanentMagnetSynchronousMachine,
            ModernMachineTest.DoublyFedInductionMachine,
        ) ? ModernMachineTest.MachineGeneratorMode : ModernMachineTest.MachineMotorMode
    control_parameters = ModernMachineTest.MachineControlParameters(
        enabled=controls,
        task_period_s=5.0 * timestep_s,
        voltage_reference_v=230.0,
        excitation_gain=0.02,
        excitation_time_constant_s=0.02,
        field_voltage_min_v=0.0,
        field_voltage_max_v=20.0,
        speed_reference_rad_s=synchronous_speed,
        governor_droop_rad_s_per_nm=0.1,
        governor_time_constant_s=0.04,
        torque_min_nm=-200.0,
        torque_max_nm=200.0,
        stabilizer_gain=0.01,
        stabilizer_washout_s=0.1,
        stabilizer_lead_s=0.02,
        stabilizer_lag_s=0.05,
    )
    field_voltage = family in (
        ModernMachineTest.WoundFieldSynchronousMachine,
        ModernMachineTest.SynchronousCondenserMachine,
    ) ? 5.0 : 0.0
    return ModernMachineTest.ModernMachineSpecification(
        id,
        family;
        operating_mode=operating_mode,
        pole_pairs=2,
        electrical=modern_machine_test_electrical(family),
        rotor_branches=modern_machine_test_rotor_branches(
            family;
            branch_count=branch_count,
        ),
        saturation=ModernMachineTest.MachineMagneticCoenergyLaw(
            radial_coefficient_per_wb2_h=0.02,
            cross_coefficient_per_wb2_h=0.01,
            maximum_flux_wb=20.0,
        ),
        shaft_masses=masses,
        shaft_couplings=couplings,
        electromagnetic_mass=first(masses).id,
        controls=control_parameters,
        settings=ModernMachineTest.MachineRuntimeSettings(
            timestep_s=timestep_s,
            nonlinear_tolerance=1.0e-13,
            maximum_nonlinear_iterations=12,
            energy_tolerance_j=1.0e-7,
        ),
        initialization_mode=ModernMachineTest.SpecifiedMachineInitialization,
        initial_field_voltage_v=field_voltage,
        initial_rotor_voltage_dq_v=(0.0, 0.0),
        initial_mechanical_torque_nm=0.0,
    )
end

function modern_machine_balanced_voltage(time_s; amplitude_v=325.0, frequency_hz=60.0)
    angle = 2.0 * pi * frequency_hz * time_s
    return [
        amplitude_v * sin(angle),
        amplitude_v * sin(angle - 2.0 * pi / 3.0),
        amplitude_v * sin(angle + 2.0 * pi / 3.0),
        0.0,
    ]
end

@testset "modern machine phase transform and coenergy" begin
    phase_voltage = [120.0, -35.0, 18.0]
    phase_current = [2.0, -0.75, 0.4]
    angle = 0.37
    transform = ModernMachineTest.machine_phase_transform(angle)
    @test transform * transpose(transform) ≈ Matrix{Float64}(I, 3, 3) atol=2.0e-15
    rotor_voltage = ModernMachineTest.machine_phase_to_rotor(phase_voltage, angle)
    rotor_current = ModernMachineTest.machine_phase_to_rotor(phase_current, angle)
    @test dot(phase_voltage, phase_current) ≈ dot(rotor_voltage, rotor_current) atol=2.0e-13
    @test ModernMachineTest.machine_rotor_to_phase(rotor_voltage, angle) ≈
        phase_voltage atol=2.0e-13

    specification = modern_machine_test_specification(
        ModernMachineTest.PermanentMagnetSynchronousMachine,
    )
    preparation = ModernMachineTest.prepare_modern_machine(specification)
    flux = preparation.initial_flux_wb .+ collect(range(-0.02, 0.03; length=3))
    current, hessian, coenergy, margin =
        ModernMachineTest.machine_coenergy_current_hessian(preparation, flux)
    @test hessian ≈ transpose(hessian) atol=1.0e-14
    @test minimum(eigvals(Symmetric(hessian))) > 0.0
    @test coenergy >= 0.0
    @test margin > 0.0
    @test margin <= inv(maximum(eigvals(Symmetric(hessian))))
    epsilon = 1.0e-7
    for column in eachindex(flux)
        plus = copy(flux)
        minus = copy(flux)
        plus[column] += epsilon
        minus[column] -= epsilon
        current_plus = ModernMachineTest.machine_coenergy_current_hessian(
            preparation,
            plus,
        )[1]
        current_minus = ModernMachineTest.machine_coenergy_current_hessian(
            preparation,
            minus,
        )[1]
        @test (current_plus - current_minus) / (2.0 * epsilon) ≈
            hessian[:, column] rtol=2.0e-7 atol=2.0e-8
    end
    @test length(current) == 3
end

@testset "modern machine family preparation and refusal" begin
    families = (
        ModernMachineTest.WoundFieldSynchronousMachine,
        ModernMachineTest.CageInductionMachine,
        ModernMachineTest.WoundRotorInductionMachine,
        ModernMachineTest.PermanentMagnetSynchronousMachine,
        ModernMachineTest.DoublyFedInductionMachine,
        ModernMachineTest.SynchronousCondenserMachine,
    )
    for family in families
        branch_count = family === ModernMachineTest.CageInductionMachine ? 4 : 1
        specification = modern_machine_test_specification(
            family;
            branch_count=branch_count,
        )
        readiness = ModernMachineTest.modern_machine_readiness(specification)
        @test readiness.ready
        @test readiness.state_count >= 3
        @test length(readiness.deterministic_signature_sha256) == 64
        preparation = ModernMachineTest.prepare_modern_machine(specification)
        @test preparation.deterministic_signature_sha256 ==
            readiness.deterministic_signature_sha256
        @test minimum(eigvals(Symmetric(preparation.inductance_h))) > 0.0
    end

    invalid_cage = ModernMachineTest.ModernMachineSpecification(
        :invalid_cage,
        ModernMachineTest.CageInductionMachine;
        pole_pairs=2,
        electrical=modern_machine_test_electrical(
            ModernMachineTest.CageInductionMachine,
        ),
        rotor_branches=ModernMachineTest.MachineRotorBranch[],
        shaft_masses=[ModernMachineTest.MachineShaftMass(
            :rotor;
            inertia_kg_m2=1.0,
        )],
        settings=ModernMachineTest.MachineRuntimeSettings(timestep_s=1.0e-5),
    )
    @test !ModernMachineTest.modern_machine_readiness(invalid_cage).ready
    @test_throws ModernMachineTest.ModernMachineRefusal ModernMachineTest.prepare_modern_machine(
        invalid_cage,
    )
    @test_throws ArgumentError ModernMachineTest.MachineMagneticCoenergyLaw(
        radial_coefficient_per_wb2_h=0.01,
        cross_coefficient_per_wb2_h=0.021,
        maximum_flux_wb=2.0,
    )
end

@testset "modern machine analytic terminal companion" begin
    preparation = ModernMachineTest.prepare_modern_machine(
        modern_machine_test_specification(
            ModernMachineTest.CageInductionMachine;
            branch_count=4,
        ),
    )
    runtime = ModernMachineTest.modern_machine_runtime(preparation)
    step = preparation.specification.settings.timestep_s
    ModernMachineTest.prepare_nonlinear_device_step!(
        runtime,
        step,
        step,
        :trapezoidal,
    )
    voltage = modern_machine_balanced_voltage(step)
    current = zeros(4)
    jacobian = zeros(4, 4)
    ModernMachineTest.nonlinear_current_jacobian!(
        current,
        jacobian,
        runtime,
        voltage,
        step,
    )
    candidate_flux = copy(runtime.candidate.flux_wb)
    candidate_hessian = copy(runtime.candidate.current_flux_jacobian_per_h)
    candidate_flux_sensitivity = copy(
        runtime.candidate.terminal_flux_sensitivity_wb_per_v,
    )
    flux_epsilon = 1.0e-7
    for column in eachindex(candidate_flux)
        plus_flux = copy(candidate_flux)
        minus_flux = copy(candidate_flux)
        plus_flux[column] += flux_epsilon
        minus_flux[column] -= flux_epsilon
        plus_winding_current = ModernMachineTest.machine_coenergy_current_hessian(
            preparation,
            plus_flux,
        )[1]
        minus_winding_current = ModernMachineTest.machine_coenergy_current_hessian(
            preparation,
            minus_flux,
        )[1]
        @test (plus_winding_current - minus_winding_current) / (2.0 * flux_epsilon) ≈
            candidate_hessian[:, column] rtol=2.0e-7 atol=2.0e-7
    end
    epsilon = 1.0e-1
    mass_index = ModernMachineTest._machine_electromagnetic_mass_index(
        preparation.specification,
    )
    electrical_angle = preparation.specification.pole_pairs *
        runtime.shaft_state.angle_rad[mass_index]
    speed_matrix = ModernMachineTest._machine_electrical_speed_matrix(
        preparation,
        runtime.shaft_state.speed_rad_s[mass_index],
    )
    tangent = Matrix{Float64}(I, length(candidate_flux), length(candidate_flux)) .-
        0.5 * step .* (
            -Diagonal(preparation.resistance_ohm) * candidate_hessian + speed_matrix
        )
    input_map = ModernMachineTest._machine_terminal_input_matrix(
        preparation,
        electrical_angle,
    )
    analytic_flux_sensitivity = tangent \ (0.5 * step .* input_map)
    @test candidate_flux_sensitivity ≈ analytic_flux_sensitivity atol=2.0e-15
    for column in 1:4
        plus_voltage = copy(voltage)
        minus_voltage = copy(voltage)
        plus_voltage[column] += epsilon
        minus_voltage[column] -= epsilon
        plus_current = zeros(4)
        minus_current = zeros(4)
        workspace = zeros(4, 4)
        ModernMachineTest.nonlinear_current_jacobian!(
            plus_current,
            workspace,
            runtime,
            plus_voltage,
            step,
        )
        plus_flux = copy(runtime.candidate.flux_wb)
        ModernMachineTest.nonlinear_current_jacobian!(
            minus_current,
            workspace,
            runtime,
            minus_voltage,
            step,
        )
        minus_flux = copy(runtime.candidate.flux_wb)
        @test (plus_flux - minus_flux) / (2.0 * epsilon) ≈
            analytic_flux_sensitivity[:, column] rtol=2.0e-5 atol=2.0e-10
        @test (plus_current - minus_current) / (2.0 * epsilon) ≈
            jacobian[:, column] rtol=2.0e-5 atol=2.0e-8
    end
    @test maximum(abs, sum(jacobian; dims=1)) <= 1.0e-12
    @test maximum(abs, sum(jacobian; dims=2)) <= 1.0e-12
    @test abs(sum(current)) <= 1.0e-12
    ModernMachineTest.accept_nonlinear_device_state!(
        runtime,
        voltage,
        current,
        jacobian,
        step,
    )
    ModernMachineTest.finish_nonlinear_device_step!(runtime)
    @test runtime.accepted_state.accepted_step_count == 1
    @test runtime.accepted_state.maximum_flux_residual_wb <= 1.0e-10
end

@testset "modern machine families execute state, controls, events, and restart" begin
    families = (
        ModernMachineTest.WoundFieldSynchronousMachine,
        ModernMachineTest.CageInductionMachine,
        ModernMachineTest.WoundRotorInductionMachine,
        ModernMachineTest.PermanentMagnetSynchronousMachine,
        ModernMachineTest.DoublyFedInductionMachine,
        ModernMachineTest.SynchronousCondenserMachine,
    )
    for (family_index, family) in enumerate(families)
        branch_count = family === ModernMachineTest.CageInductionMachine ? 4 : 1
        specification = modern_machine_test_specification(
            family;
            id=Symbol("executable_machine_$(family_index)"),
            branch_count=branch_count,
            mass_count=family_index == 1 ? 8 : 1,
            controls=family_index in (1, 6),
        )
        preparation = ModernMachineTest.prepare_modern_machine(specification)
        step = specification.settings.timestep_s
        events = ModernMachineTest.ModernMachineEvent[
            ModernMachineTest.ModernMachineEvent(
                :torque_step,
                10 * step,
                ModernMachineTest.MachineMechanicalTorqueEvent;
                value=12.0,
            ),
        ]
        trace = ModernMachineTest.simulate_modern_machine(
            preparation,
            modern_machine_balanced_voltage;
            duration_s=20 * step,
            events=events,
        )
        @test length(trace.time_s) == 21
        @test all(isfinite, trace.phase_current_a)
        @test all(isfinite, trace.electromagnetic_torque_nm)
        @test trace.event_count[end] == 1
        @test maximum(abs, vec(sum(trace.phase_current_a; dims=1))) <= 2.0e-12
        @test trace.result.status == :accepted
        @test trace.result.diagnostics.accepted_step_count == 20
        @test trace.result.diagnostics.maximum_flux_residual_wb <= 1.0e-9
        @test trace.result.diagnostics.maximum_kcl_residual_a <= 2.0e-12
        @test trace.dissipated_energy_j[end] >= 0.0
        if family_index in (1, 6)
            @test trace.control_sample_count[end] >= 4
        end
        if family === ModernMachineTest.DoublyFedInductionMachine
            @test ModernMachineTest.modern_machine_result_quantity(
                trace.result,
                :rotor_current_dq_a,
            ) isa Matrix{Float64}
        end
        if family !== ModernMachineTest.PermanentMagnetSynchronousMachine
            @test !ModernMachineTest.modern_machine_result_quantity(
                trace.result,
                :permanent_magnet_flux_wb,
            ).available
        end
    end

    preparation = ModernMachineTest.prepare_modern_machine(
        modern_machine_test_specification(
            ModernMachineTest.PermanentMagnetSynchronousMachine,
        ),
    )
    step = preparation.specification.settings.timestep_s
    uninterrupted = ModernMachineTest.modern_machine_runtime(preparation)
    for index in 1:20
        ModernMachineTest.advance_modern_machine!(
            uninterrupted,
            modern_machine_balanced_voltage(index * step);
            time_s=index * step,
        )
    end
    split = ModernMachineTest.modern_machine_runtime(preparation)
    for index in 1:9
        ModernMachineTest.advance_modern_machine!(
            split,
            modern_machine_balanced_voltage(index * step);
            time_s=index * step,
        )
    end
    snapshot = ModernMachineTest.modern_machine_runtime_snapshot(split)
    restored = ModernMachineTest.modern_machine_runtime(preparation)
    ModernMachineTest.restore_modern_machine_runtime_snapshot!(restored, snapshot)
    for index in 10:20
        ModernMachineTest.advance_modern_machine!(
            restored,
            modern_machine_balanced_voltage(index * step);
            time_s=index * step,
        )
    end
    @test restored.accepted_state.flux_wb == uninterrupted.accepted_state.flux_wb
    @test restored.accepted_state.terminal_current_a ==
        uninterrupted.accepted_state.terminal_current_a
    @test restored.shaft_state.angle_rad == uninterrupted.shaft_state.angle_rad
    @test restored.shaft_state.speed_rad_s == uninterrupted.shaft_state.speed_rad_s
    @test restored.accepted_state.supplied_electrical_energy_j ==
        uninterrupted.accepted_state.supplied_electrical_energy_j
end

@testset "modern machine result and report" begin
    preparation = ModernMachineTest.prepare_modern_machine(
        modern_machine_test_specification(
            ModernMachineTest.WoundFieldSynchronousMachine;
            controls=true,
        ),
    )
    trace = ModernMachineTest.simulate_modern_machine(
        preparation,
        modern_machine_balanced_voltage;
        duration_s=5 * preparation.specification.settings.timestep_s,
    )
    report = ModernMachineTest.modern_machine_report_text(trace)
    @test occursin("wound_field_synchronous", report)
    @test occursin(trace.result.deterministic_signature_sha256, report)
    mktempdir() do directory
        path = ModernMachineTest.write_modern_machine_trace_csv(
            joinpath(directory, "machine.csv"),
            trace,
        )
        text = read(path, String)
        @test occursin("electromagnetic_torque_nm", text)
        @test count(==('\n'), text) == length(trace.time_s) + 1
    end
end

@testset "modern machine simultaneous event ordering" begin
    preparation = ModernMachineTest.prepare_modern_machine(
        modern_machine_test_specification(
            ModernMachineTest.PermanentMagnetSynchronousMachine,
        ),
    )
    step = preparation.specification.settings.timestep_s
    runtime = ModernMachineTest.modern_machine_runtime(
        preparation;
        events=ModernMachineTest.ModernMachineEvent[
            ModernMachineTest.ModernMachineEvent(
                :second_torque,
                step,
                ModernMachineTest.MachineMechanicalTorqueEvent;
                value=12.0,
                priority=2,
            ),
            ModernMachineTest.ModernMachineEvent(
                :first_torque,
                step,
                ModernMachineTest.MachineMechanicalTorqueEvent;
                value=4.0,
                priority=1,
            ),
        ],
    )
    ModernMachineTest.advance_modern_machine!(
        runtime,
        modern_machine_balanced_voltage(step);
        time_s=step,
    )
    @test runtime.accepted_state.event_count == 2
    @test runtime.inputs.mechanical_torque_nm == 12.0
end

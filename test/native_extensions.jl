using UUIDs

const NativeExtensionAPI = AIMORA.NativeExtensions

@testset "caller-owned native extension registry" begin
    registry = NativeExtensionAPI.ExtensionRegistry()
    types = (
        NativeExtensionAPI.SampledSaturatingLag,
        NativeExtensionAPI.NativeCubicCurrentBranch,
        NativeExtensionAPI.NativeSeriesRLCompanion,
    )
    registrations = [NativeExtensionAPI.register_extension!(registry, type) for type in types]
    @test length(registrations) == 3
    @test NativeExtensionAPI.registered_extension_identities(registry) == sort(
        getfield.(registrations, :identity);
        by = identity -> (
            identity.namespace,
            identity.semantic_type,
            identity.semantic_version,
            string(identity.package_uuid),
            identity.content_sha256,
        ),
    )
    for registration in registrations
        @test NativeExtensionAPI.resolve_extension(registry, registration.identity) === registration
        @test registration.identity.package_uuid ==
            UUID("bb4e3f50-0ac0-4ed3-9c38-7355ab8ea495")
    end
    unknown = NativeExtensionAPI.ExtensionIdentity(
        UUID("4cf129db-b3fb-41c2-8667-cf22679bd5ba"),
        "example.extensions",
        "unknown_component",
        v"1.0.0",
        :UnknownComponent,
        v"1.0.0",
        repeat("a", 64),
    )
    @test_throws NativeExtensionAPI.ExtensionFailure NativeExtensionAPI.resolve_extension(
        registry,
        unknown,
    )
end

@testset "sampled saturating lag task checkpoint and event" begin
    control = NativeExtensionAPI.SampledSaturatingLag(
        2.0,
        2.0e-3,
        -1.0,
        1.0,
        10;
        delay_ticks = 2,
    )
    registry = NativeExtensionAPI.ExtensionRegistry()
    NativeExtensionAPI.register_extension!(registry, typeof(control))
    checkpoint = NativeExtensionAPI.extension_checkpoint(control)
    candidate = NativeExtensionAPI.sample_extension_task!(control, 0.75, 0, 1.0e-4)
    expected_state = -expm1(-1.0e-3 / 2.0e-3) * 0.75
    @test control.state ≈ expected_state atol = 1.0e-15
    @test candidate ≈ 2.0 * expected_state atol = 1.0e-15
    @test control.held_output == 0.0
    @test NativeExtensionAPI.release_extension_task_output!(control, 2) == candidate
    @test control.sample_count == 1
    @test control.write_count == 1
    outputs = NativeExtensionAPI.extension_outputs(control, 2.0e-4)
    @test getfield.(outputs, :name) == (:state, :output)
    @test NativeExtensionAPI.extension_source_value(control, 2.0e-4) == candidate
    execution = NativeExtensionAPI.ExtensionExecutionResult(
        true,
        NativeExtensionAPI.extension_contract(control),
        (1,),
        repeat("1", 64),
        outputs,
        1,
        1,
        checkpoint.state_sha256,
        (control_error = 0.0,),
        ("declared_control_unit",),
    )
    @test execution.accepted
    @test execution.output_cursor == 2
    @test length(execution.deterministic_sha256) == 64
    @test NativeExtensionAPI.validate_extension_result(
        execution,
        NativeExtensionAPI.extension_contract(control),
        repeat("1", 64),
        checkpoint.state_sha256,
        ("declared_control_unit", "declared_control_unit"),
    ) === execution
    stale_project_failure = try
        NativeExtensionAPI.validate_extension_result(
            execution,
            NativeExtensionAPI.extension_contract(control),
            repeat("2", 64),
            checkpoint.state_sha256,
            ("declared_control_unit", "declared_control_unit"),
        )
        nothing
    catch exception
        exception
    end
    @test stale_project_failure isa NativeExtensionAPI.ExtensionFailure
    @test stale_project_failure.code == :stale_extension_project_signature
    wrong_unit_failure = try
        NativeExtensionAPI.validate_extension_result(
            execution,
            NativeExtensionAPI.extension_contract(control),
            repeat("1", 64),
            checkpoint.state_sha256,
            ("A", "A"),
        )
        nothing
    catch exception
        exception
    end
    @test wrong_unit_failure isa NativeExtensionAPI.ExtensionFailure
    @test wrong_unit_failure.code == :extension_output_unit_mismatch
    NativeExtensionAPI.restore_extension_checkpoint!(control, checkpoint)
    @test control.state == 0.0
    @test control.sample_count == 0
    @test NativeExtensionAPI.extension_state_signature(checkpoint.state) ==
        checkpoint.state_sha256

    event = NativeExtensionAPI.DirectedExtensionEvent(
        :output,
        0.25;
        direction = :rising,
        maximum_occurrences = 1,
    )
    @test NativeExtensionAPI.extension_event_crossed(event, 0.0, candidate)
    @test NativeExtensionAPI.accept_extension_event!(event) == 1
    @test_throws NativeExtensionAPI.ExtensionFailure NativeExtensionAPI.accept_extension_event!(event)
end

@testset "passive cubic extension residual and analytic Jacobian" begin
    device = NativeExtensionAPI.NativeCubicCurrentBranch(1, 0, 0.2, 0.03)
    currents = zeros(2)
    jacobian = zeros(2, 2)
    voltages = [1.5, 0.0]
    AIMORA.NonlinearNetwork.nonlinear_current_jacobian!(
        currents,
        jacobian,
        device,
        voltages,
        0.0,
    )
    @test currents == [0.2 * 1.5 + 0.03 * 1.5^3, -(0.2 * 1.5 + 0.03 * 1.5^3)]
    derivative = 0.2 + 3.0 * 0.03 * 1.5^2
    @test jacobian == [derivative -derivative; -derivative derivative]
    perturbation = 1.0e-6
    shifted = zeros(2)
    AIMORA.NonlinearNetwork.nonlinear_current_jacobian!(
        shifted,
        zeros(2, 2),
        device,
        [1.5 + perturbation, 0.0],
        0.0,
    )
    @test (shifted[1] - currents[1]) / perturbation ≈ derivative rtol = 1.0e-6
    outputs = NativeExtensionAPI.extension_outputs(device, voltages, 0.0)
    @test outputs[2].value >= 0.0
    checkpoint = NativeExtensionAPI.extension_checkpoint(device)
    @test NativeExtensionAPI.restore_extension_checkpoint!(device, checkpoint) === device
end

@testset "stateful series R-L extension trial acceptance and restore" begin
    component = NativeExtensionAPI.NativeSeriesRLCompanion(
        1,
        0,
        2.0,
        10.0e-3;
        initial_current_a = 0.5,
        initial_voltage_v = 1.0,
    )
    step_s = 100.0e-6
    checkpoint = NativeExtensionAPI.extension_checkpoint(component)
    companion = NativeExtensionAPI.extension_companion(component, step_s)
    expected_g = inv(2.0 + 2.0 * 10.0e-3 / step_s)
    expected_history = expected_g * (1.0 + (2.0 * 10.0e-3 / step_s - 2.0) * 0.5)
    @test companion.conductance_s == expected_g
    @test companion.history_current_a == expected_history
    @test component.previous_current_a == 0.5
    accepted_current = NativeExtensionAPI.accept_extension_state!(component, [4.0, 0.0], step_s)
    @test accepted_current == expected_g * 4.0 + expected_history
    @test component.previous_voltage_v == 4.0
    NativeExtensionAPI.restore_extension_checkpoint!(component, checkpoint)
    @test component.previous_current_a == 0.5
    @test component.previous_voltage_v == 1.0
    @test getfield.(NativeExtensionAPI.extension_outputs(component, 0.0), :name) ==
        (:terminal_current, :stored_energy, :dissipated_power)
end

@testset "deterministic one-step extension migrations" begin
    registry = NativeExtensionAPI.ExtensionRegistry()
    from = NativeExtensionAPI.extension_identity(NativeExtensionAPI.SampledSaturatingLag)
    to = NativeExtensionAPI.ExtensionIdentity(
        from.package_uuid,
        from.namespace,
        from.semantic_type,
        v"1.1.0",
        from.implementation_symbol,
        from.api_version,
        repeat("d", 64),
    )
    migrate_state(state) = (state = state.state, added_limit = 1.0)
    NativeExtensionAPI.register_extension_migration!(registry, from, to, migrate_state)
    source = (state = 0.25,)
    @test NativeExtensionAPI.migrate_extension_state(registry, from, to, source) ==
        (state = 0.25, added_limit = 1.0)
    @test source == (state = 0.25,)
end

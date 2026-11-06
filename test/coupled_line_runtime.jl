using AIMORA.CoupledLineFitting
using AIMORA.CoupledLineRuntime


function manufactured_coupled_line_runtime_preparation(;
    timestep_s=10.0e-6,
    maximum_frequency_hz=1_000.0,
    uncertainty_alternative_fit_signatures_sha256=[repeat("d", 64)],
    uncertainty_set_complete=false,
)
    rate_per_s = 2.0 * pi * 500.0
    direct_term = [
        0.050 0.010 0.020 0.005
        0.010 0.055 0.005 0.018
        0.020 0.005 0.050 0.010
        0.005 0.018 0.010 0.055
    ]
    normalized_residue = [
        0.180 0.020 0.040 0.010
        0.020 0.160 0.010 0.035
        0.040 0.010 0.180 0.020
        0.010 0.035 0.020 0.160
    ]
    model = coupled_line_rational_model(
        [rate_per_s],
        direct_term,
        [rate_per_s .* normalized_residue];
        port_order=[:sending_a, :sending_b, :receiving_a, :receiving_b],
        reference_impedance_ohm=[40.0, 55.0, 40.0, 55.0],
    )
    frequencies_hz = [1.0, 100.0, maximum_frequency_hz]
    certificate = coupled_line_passivity_certificate(model; frequencies_hz)
    @assert certificate.continuous_passivity_passed
    settings = CoupledLineRuntimeSettings(; timestep_s)
    preparation = CoupledLineRuntime._coupled_line_runtime_preparation(
        model,
        settings;
        source_signature_sha256=repeat("a", 64),
        response_signature_sha256=repeat("b", 64),
        fit_signature_sha256=repeat("c", 64),
        phase_order=[:a, :b],
        frequencies_hz,
        continuous_passivity_passed=true,
        uncertainty_alternative_fit_signatures_sha256,
        uncertainty_set_complete,
    )
    return (; preparation, model, rate_per_s, frequencies_hz)
end

@testset "coupled phase-domain line runtime preparation and recurrence" begin
    fixture = manufactured_coupled_line_runtime_preparation()
    preparation = fixture.preparation
    settings = preparation.settings
    expected_pole = (preparation.bilinear_alpha_per_s - fixture.rate_per_s) /
        (preparation.bilinear_alpha_per_s + fixture.rate_per_s)
    @test eigvals(preparation.state_transition) ≈ fill(expected_pole, 4) atol=1.0e-14
    @test preparation.maximum_discrete_pole_magnitude < 1.0
    @test preparation.sampled_maximum_scattering_singular_value < 1.0
    @test preparation.minimum_companion_conductance_eigenvalue_s > 0.0
    @test preparation.uncertainty_alternative_fit_signatures_sha256 ==
        [repeat("d", 64)]
    @test !preparation.uncertainty_set_complete
    @test preparation.companion_admittance_s ≈
        transpose(preparation.companion_admittance_s) atol=1.0e-14
    @test CoupledLineRuntimeSettings(timestep_s=settings.timestep_s).
        deterministic_signature_sha256 == settings.deterministic_signature_sha256
    @test_throws ArgumentError CoupledLineRuntimeSettings(timestep_s=0.0)
    @test_throws ArgumentError CoupledLineRuntimeSettings(
        timestep_s=settings.timestep_s,
        prewarp_frequency_hz=1.0 / (2.0 * settings.timestep_s),
    )
    @test_throws ArgumentError manufactured_coupled_line_runtime_preparation(
        uncertainty_alternative_fit_signatures_sha256=String[],
        uncertainty_set_complete=true,
    )
    @test_throws ArgumentError manufactured_coupled_line_runtime_preparation(
        uncertainty_alternative_fit_signatures_sha256=[repeat("c", 64)],
    )
    @test_throws ArgumentError manufactured_coupled_line_runtime_preparation(
        uncertainty_alternative_fit_signatures_sha256=[
            repeat("d", 64),
            repeat("d", 64),
        ],
    )

    frequency_hz = 500.0
    mapped_s = 1.0im * preparation.bilinear_alpha_per_s *
        tan(pi * frequency_hz * settings.timestep_s)
    @test coupled_line_runtime_discrete_response(preparation, frequency_hz) ≈
        coupled_line_model_value(fixture.model, mapped_s) atol=2.0e-14
    terminal_admittance = coupled_line_runtime_terminal_admittance(
        preparation,
        frequency_hz,
    )
    @test terminal_admittance ≈ transpose(terminal_admittance) atol=2.0e-14
    @test_throws ArgumentError coupled_line_runtime_discrete_response(
        preparation,
        2_000.0,
    )

    state = coupled_line_runtime_state(preparation)
    terminal_voltage_v = [100.0, -50.0, 20.0, -10.0]
    expected_incident = preparation.incident_from_voltage * terminal_voltage_v
    expected_state = preparation.endpoint_input * expected_incident
    expected_outgoing = preparation.output_matrix * expected_state +
        preparation.continuous_direct_term * expected_incident
    expected_current = preparation.reference_impedance_inverse_sqrt_per_ohm_sqrt .*
        (expected_incident - expected_outgoing)
    initial_snapshot = coupled_line_runtime_snapshot(state)
    @test accept_coupled_line_runtime_step!(state, terminal_voltage_v) === state
    @test state.rational_state ≈ expected_state atol=1.0e-14
    @test state.incident_wave ≈ expected_incident atol=1.0e-14
    @test state.outgoing_wave ≈ expected_outgoing atol=1.0e-14
    @test state.terminal_current_a ≈ expected_current atol=1.0e-14
    @test state.accepted_step_count == 1
    @test state.accepted_time_s == settings.timestep_s
    @test state.maximum_kcl_residual_a <= settings.kcl_absolute_tolerance_a
    @test state.cumulative_supplied_energy_j >= -settings.energy_absolute_tolerance_j
    @test coupled_line_runtime_diagnostics(state).passive_energy_balance_passed

    accepted_snapshot = coupled_line_runtime_snapshot(state)
    mktempdir() do directory
        snapshot_path = write_coupled_line_runtime_snapshot(
            joinpath(directory, "runtime_snapshot.toml"),
            accepted_snapshot,
        )
        restored_snapshot = read_coupled_line_runtime_snapshot(snapshot_path)
        @test restored_snapshot.deterministic_signature_sha256 ==
            accepted_snapshot.deterministic_signature_sha256
        replay = coupled_line_runtime_state(preparation)
        restore_coupled_line_runtime_snapshot!(replay, restored_snapshot)
        @test replay.rational_state == state.rational_state
        malformed_path = joinpath(directory, "malformed_snapshot.toml")
        write(malformed_path, "schema = \"unknown\"\n")
        @test_throws ArgumentError read_coupled_line_runtime_snapshot(malformed_path)
    end
    @test_throws ArgumentError accept_coupled_line_runtime_step!(
        state,
        [NaN, 0.0, 0.0, 0.0],
    )
    @test coupled_line_runtime_snapshot(state).deterministic_signature_sha256 ==
        accepted_snapshot.deterministic_signature_sha256
    rejected_state = coupled_line_runtime_state(preparation)
    rejected_state.previous_terminal_power_w = -1.0e9
    rejected_snapshot = coupled_line_runtime_snapshot(rejected_state)
    @test_throws ArgumentError accept_coupled_line_runtime_step!(
        rejected_state,
        zeros(length(preparation.port_order)),
    )
    @test coupled_line_runtime_snapshot(rejected_state).
        deterministic_signature_sha256 ==
        rejected_snapshot.deterministic_signature_sha256
    restore_coupled_line_runtime_snapshot!(state, initial_snapshot)
    @test state.accepted_step_count == 0
    @test all(iszero, state.rational_state)

    too_coarse = CoupledLineRuntimeSettings(timestep_s=1.0e-3)
    @test_throws ArgumentError CoupledLineRuntime._coupled_line_runtime_preparation(
        fixture.model,
        too_coarse;
        source_signature_sha256=repeat("a", 64),
        response_signature_sha256=repeat("b", 64),
        fit_signature_sha256=repeat("c", 64),
        phase_order=[:a, :b],
        frequencies_hz=fixture.frequencies_hz,
        continuous_passivity_passed=true,
    )
    @test_throws ArgumentError CoupledLineRuntime._coupled_line_runtime_preparation(
        fixture.model,
        settings;
        source_signature_sha256=repeat("a", 64),
        response_signature_sha256=repeat("b", 64),
        fit_signature_sha256=repeat("c", 64),
        phase_order=[:a, :b],
        frequencies_hz=fixture.frequencies_hz,
        continuous_passivity_passed=false,
    )
end

@testset "coupled phase-domain line sinusoidal initialization and exact restart" begin
    fixture = manufactured_coupled_line_runtime_preparation()
    preparation = fixture.preparation
    frequency_hz = 60.0
    voltage_phasor_v = ComplexF64[
        120.0 + 5.0im,
        -55.0 + 20.0im,
        35.0 - 10.0im,
        -15.0 - 8.0im,
    ]
    current_phasor_a = coupled_line_runtime_terminal_admittance(
        preparation,
        frequency_hz,
    ) * voltage_phasor_v
    state = coupled_line_runtime_state(preparation)
    initialize_coupled_line_runtime_sinusoidal!(
        state,
        voltage_phasor_v,
        frequency_hz,
    )
    @test state.initialization_kind == :sinusoidal_discrete_operating_point
    @test state.terminal_current_a ≈ real.(current_phasor_a) atol=2.0e-14
    angle = 2.0 * pi * frequency_hz * preparation.settings.timestep_s
    restart_snapshot = nothing
    for step in 1:20
        voltage_v = real.(voltage_phasor_v .* cis(step * angle))
        accept_coupled_line_runtime_step!(state, voltage_v)
        @test state.terminal_current_a ≈
            real.(current_phasor_a .* cis(step * angle)) atol=2.0e-11
        step == 10 && (restart_snapshot = coupled_line_runtime_snapshot(state))
    end
    uninterrupted_snapshot = coupled_line_runtime_snapshot(state)

    replay = coupled_line_runtime_state(preparation)
    restore_coupled_line_runtime_snapshot!(replay, something(restart_snapshot))
    for step in 11:20
        voltage_v = real.(voltage_phasor_v .* cis(step * angle))
        accept_coupled_line_runtime_step!(replay, voltage_v)
    end
    replay_snapshot = coupled_line_runtime_snapshot(replay)
    @test replay_snapshot.deterministic_signature_sha256 ==
        uninterrupted_snapshot.deterministic_signature_sha256
    @test replay.rational_state == state.rational_state
    @test replay.history_current_a == state.history_current_a

    stale_fixture = manufactured_coupled_line_runtime_preparation(timestep_s=5.0e-6)
    @test_throws ArgumentError restore_coupled_line_runtime_snapshot!(
        coupled_line_runtime_state(stale_fixture.preparation),
        something(restart_snapshot),
    )
    snapshot = something(restart_snapshot)
    values = Any[getfield(snapshot, name) for name in fieldnames(typeof(snapshot))]
    state_index = findfirst(==(:rational_state), fieldnames(typeof(snapshot)))
    values[state_index] = copy(snapshot.rational_state)
    values[state_index][1] += 1.0
    tampered = CoupledLineRuntimeSnapshot(values...)
    @test_throws ArgumentError restore_coupled_line_runtime_snapshot!(replay, tampered)

    diagnostics = coupled_line_runtime_diagnostics(state)
    @test diagnostics.source_signature_sha256 == repeat("a", 64)
    @test diagnostics.response_signature_sha256 == repeat("b", 64)
    @test diagnostics.fit_signature_sha256 == repeat("c", 64)
    @test diagnostics.uncertainty_alternative_fit_signatures_sha256 ==
        [repeat("d", 64)]
    @test diagnostics.unknown_uncertainty_explicit
    @test diagnostics.coupled_phase_domain_runtime_executed
    @test !diagnostics.ulm_file_compatibility_claimed
    @test !diagnostics.atp_or_pscad_equivalence_claimed
    report = coupled_line_runtime_report_text(state)
    @test occursin("runtime_signature_sha256=", report)
    @test occursin("accepted_step_count=20", report)
end

@testset "coupled line event-step realization and clock" begin
    fixture = manufactured_coupled_line_runtime_preparation()
    preparation, model = fixture.preparation, fixture.model
    state = coupled_line_runtime_state(preparation)
    ns, np = length(state.rational_state), length(state.incident_wave)
    for (index, (duration, implicit)) in enumerate(
        [(2.5e-6, true), (7.5e-6, true), (10.0e-6, false),
            (1.0e-10, true), (4.0e-6, false), (6.0e-6, false)])
        old_state = copy(state.rational_state)
        old_incident = copy(state.previous_incident_wave)
        voltage = [100sin(index), -30cos(index), 20sin(index), -10cos(index)]
        theta = implicit ? 1.0 : 0.5
        # Solve independent continuous state+wave equations simultaneously,
        # not the production companion or its history reconstruction.
        block = vcat(
            hcat(Matrix{Float64}(I, ns, ns) - theta * duration * model.state_matrix_per_s,
                -theta * duration * model.input_matrix),
            hcat(model.output_matrix_per_s, Matrix{Float64}(I, np, np) + model.direct_term),
        )
        rhs_state = implicit ? old_state :
            (Matrix{Float64}(I, ns, ns) + theta * duration * model.state_matrix_per_s) *
                old_state + theta * duration * model.input_matrix * old_incident
        reference = block \ vcat(rhs_state,
            preparation.reference_impedance_inverse_sqrt_per_ohm_sqrt .* voltage)
        expected_state = reference[1:ns]
        incident = reference[ns+1:end]
        outgoing = model.output_matrix_per_s * expected_state + model.direct_term * incident
        current = preparation.reference_impedance_inverse_sqrt_per_ohm_sqrt .* (incident - outgoing)
        step = coupled_line_runtime_step(preparation, duration; backward_euler=implicit)
        before = coupled_line_runtime_snapshot(state).deterministic_signature_sha256
        @test step.preparation_signature_sha256 == preparation.deterministic_signature_sha256
        @test coupled_line_runtime_snapshot(state).deterministic_signature_sha256 == before
        @test minimum(eigvals(Symmetric(step.companion_admittance_s))) > 0.0
        @test accept_coupled_line_runtime_step!(state, voltage, step) === state
        @test state.rational_state ≈ expected_state rtol=2e-13 atol=2e-14
        @test state.incident_wave ≈ incident rtol=2e-13 atol=2e-14
        @test state.terminal_current_a ≈ current rtol=2e-13 atol=2e-14
        @test state.maximum_kcl_residual_a <= preparation.settings.kcl_absolute_tolerance_a
        @test state.cumulative_supplied_energy_j >= -preparation.settings.energy_absolute_tolerance_j
    end
    @test state.accepted_step_count == 1
    @test state.accepted_substep_count == 5
    @test state.accepted_time_s ≈ 30e-6 + 1e-10 atol=1e-20
    snapshot = coupled_line_runtime_snapshot(state)
    @test snapshot.schema_version == 2
    mktempdir() do directory
        path = write_coupled_line_runtime_snapshot(joinpath(directory, "event_steps.toml"), snapshot)
        replay = coupled_line_runtime_state(preparation)
        restore_coupled_line_runtime_snapshot!(replay, read_coupled_line_runtime_snapshot(path))
        @test coupled_line_runtime_snapshot(replay).deterministic_signature_sha256 == snapshot.deterministic_signature_sha256
        for _ in 1:3
            accept_coupled_line_runtime_step!(state, fill(10.0, np))
            accept_coupled_line_runtime_step!(replay, fill(10.0, np))
        end
        @test coupled_line_runtime_snapshot(replay).deterministic_signature_sha256 ==
            coupled_line_runtime_snapshot(state).deterministic_signature_sha256
    end
    for invalid in (0.0, -1e-6, Inf, NaN, 11e-6)
        @test_throws ArgumentError coupled_line_runtime_step(preparation, invalid)
    end
    a, b = coupled_line_runtime_state(preparation), coupled_line_runtime_state(preparation)
    nominal = coupled_line_runtime_step(preparation, preparation.settings.timestep_s)
    for index in 1:20
        voltage = fill(sin(index), np)
        accept_coupled_line_runtime_step!(a, voltage)
        accept_coupled_line_runtime_step!(b, voltage, nominal)
    end
    @test coupled_line_runtime_snapshot(a).deterministic_signature_sha256 ==
        coupled_line_runtime_snapshot(b).deterministic_signature_sha256
    @test coupled_line_runtime_snapshot(a).schema_version == 1
    old_dictionary = CoupledLineRuntime._coupled_line_runtime_snapshot_dictionary(coupled_line_runtime_snapshot(a))
    delete!(old_dictionary, "accepted_substep_count")
    delete!(old_dictionary, "accepted_substep_duration_s")
    old_snapshot = CoupledLineRuntime._coupled_line_runtime_snapshot_from_dictionary(old_dictionary)
    @test restore_coupled_line_runtime_snapshot!(b, old_snapshot) === b
end

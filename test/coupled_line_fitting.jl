using AIMORA.CoupledLineFitting

@testset "coupled line response and global passivity contracts" begin
    frequencies_hz = [10.0, 100.0, 1_000.0]
    series = [
        ComplexF64[
            0.08 + 2.0im * pi * frequency * 1.0e-3 0.01 + 2.0im * pi * frequency * 0.2e-3
            0.01 + 2.0im * pi * frequency * 0.2e-3 0.09 + 2.0im * pi * frequency * 1.1e-3
        ] for frequency in frequencies_hz
    ]
    shunt = [
        2.0im * pi * frequency .* ComplexF64[
            1.2e-6 -0.2e-6
            -0.2e-6 1.3e-6
        ] for frequency in frequencies_hz
    ]
    response = coupled_line_terminal_response(
        frequencies_hz,
        series,
        shunt,
        10.0;
        phase_order=[:a, :b],
        reference_impedance_ohm=fill(50.0, 4),
        source_signature_sha256=repeat("a", 64),
        segment_id=:engine_contract_uniform_segment,
        segment_kind=:manufactured_uniform,
    )
    @test response.port_order == [
        :sending_a,
        :sending_b,
        :receiving_a,
        :receiving_b,
    ]
    @test maximum(response.reciprocity_errors) <= 1.0e-12
    @test minimum(response.minimum_physical_loss_eigenvalues_s) >= -1.0e-12
    for index in eachindex(frequencies_hz)
        reconstructed = coupled_line_scattering_to_admittance(
            response.scattering_matrices[index],
            response.reference_impedance_ohm,
        )
        @test reconstructed ≈ response.terminal_admittance_matrices_s[index] rtol=1.0e-11
    end

    modal = coupled_line_modal_preparation(
        frequencies_hz,
        series,
        shunt,
        10.0;
        phase_order=[:a, :b],
    )
    @test length(modal.extracted_delays_s) == 2
    @test all(>=(0.0), modal.extracted_delays_s)
    @test minimum(modal.minimum_mode_overlaps) >= 0.05
    @test all(isfinite, modal.delay_phase_residuals_rad)
    default_fit_settings = CoupledLineFitSettings(candidate_orders=[1])
    @test default_fit_settings.maximum_direct_term_singular_value == 0.92
    @test default_fit_settings.maximum_relocation_sweeps == 64
    @test CoupledLineFitSettings(
        candidate_orders=[1],
        maximum_direct_term_singular_value=0.91,
    ).deterministic_signature_sha256 != default_fit_settings.deterministic_signature_sha256
    @test_throws ArgumentError CoupledLineFitSettings(
        candidate_orders=[1],
        maximum_direct_term_singular_value=0.0,
    )
    @test_throws ArgumentError CoupledLineFitSettings(
        candidate_orders=[1],
        maximum_direct_term_singular_value=1.0,
    )
    @test CoupledLineFitRequest(
        response,
        default_fit_settings,
        modal,
    ).modal_preparation === modal
    long_modal = coupled_line_modal_preparation(
        frequencies_hz,
        series,
        shunt,
        1_000.0;
        phase_order=[:a, :b],
    )
    @test all(>=(0.0), long_modal.extracted_delays_s)

    passive_model = coupled_line_rational_model(
        [10.0],
        zeros(2, 2),
        [2.0 .* Matrix{Float64}(I, 2, 2)];
        port_order=[:sending_a, :receiving_a],
        reference_impedance_ohm=fill(50.0, 2),
    )
    passive_certificate = coupled_line_passivity_certificate(
        passive_model;
        frequencies_hz,
    )
    @test passive_certificate.stable_realization
    @test passive_certificate.hamiltonian_crossing_free
    @test passive_certificate.continuous_passivity_passed
    @test passive_certificate.zero_frequency_maximum_singular_value ≈ 0.2
    @test isfinite(passive_certificate.minimum_physical_loss_eigenvalue_s)
    @test passive_certificate.minimum_physical_loss_eigenvalue_s > 0.0

    state_scaling = Diagonal([1.0e8, 1.0e-8])
    scaled_passive_model = coupled_line_rational_model_from_state_space(
        state_scaling \ passive_model.state_matrix_per_s * state_scaling,
        state_scaling \ passive_model.input_matrix,
        passive_model.output_matrix_per_s * state_scaling,
        passive_model.direct_term;
        port_order=passive_model.port_order,
        reference_impedance_ohm=passive_model.reference_impedance_ohm,
    )
    scaled_passive_certificate = coupled_line_passivity_certificate(
        scaled_passive_model;
        frequencies_hz,
    )
    @test scaled_passive_certificate.continuous_passivity_passed
    @test scaled_passive_certificate.zero_frequency_maximum_singular_value ≈
        passive_certificate.zero_frequency_maximum_singular_value

    active_model = coupled_line_rational_model(
        [10.0],
        zeros(2, 2),
        [12.0 .* Matrix{Float64}(I, 2, 2)];
        port_order=[:sending_a, :receiving_a],
        reference_impedance_ohm=fill(50.0, 2),
    )
    active_certificate = coupled_line_passivity_certificate(
        active_model;
        frequencies_hz,
    )
    @test active_certificate.stable_realization
    @test !active_certificate.zero_frequency_bounded
    @test !active_certificate.continuous_passivity_passed

    oscillatory_poles = ComplexF64[-20.0 + 100.0im, -20.0 - 100.0im]
    oscillatory_residue = ComplexF64[
        2.0 + 1.0im 0.2 + 0.1im
        0.2 + 0.1im 1.5 + 0.5im
    ]
    oscillatory_model = coupled_line_rational_model_from_poles(
        oscillatory_poles,
        [0.05 0.01; 0.01 0.04],
        [oscillatory_residue, conj.(oscillatory_residue)];
        port_order=[:sending_a, :receiving_a],
        reference_impedance_ohm=fill(50.0, 2),
    )
    @test oscillatory_model.poles_per_s == oscillatory_poles
    @test eltype(oscillatory_model.state_matrix_per_s) == Float64
    for frequency_hz in frequencies_hz
        s = 2.0im * pi * frequency_hz
        state_space_value = ComplexF64.(oscillatory_model.direct_term) +
            oscillatory_model.output_matrix_per_s *
            ((s * I - oscillatory_model.state_matrix_per_s) \
                oscillatory_model.input_matrix)
        @test state_space_value ≈ coupled_line_model_value(oscillatory_model, s) atol=1.0e-13
    end
    @test_throws ArgumentError coupled_line_rational_model_from_poles(
        oscillatory_poles[1:1],
        [0.05 0.01; 0.01 0.04],
        [oscillatory_residue];
        port_order=[:sending_a, :receiving_a],
        reference_impedance_ohm=fill(50.0, 2),
    )
end

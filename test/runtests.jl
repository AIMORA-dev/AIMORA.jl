using Test
using AIMORA
using LinearAlgebra
using Printf

@testset "public package isolation" begin
    for module_name in (
        :LegacyEMTP,
        :ExecutableComparison,
        :TranslationProgress,
        :ValidationBenchmarks,
    )
        @test !isdefined(AIMORA, module_name)
    end
end

@testset "open engineering core" begin
    issue = AIMORA.ValidationCore.missing_data("line", "length is required")
    result = AIMORA.ValidationCore.validation_result(source = "public package test")
    AIMORA.ValidationCore.add_issue!(result, issue)
    @test !AIMORA.ValidationCore.is_valid(result)

    scenario = AIMORA.ProjectData.Scenario(id = :base, name = "Base")
    row = AIMORA.InverterAssets.inverter_row(
        id = :inv1,
        bus = :bus1,
        rated_kva = 1000.0,
        v_ll_rms_v = 4160.0,
    )
    AIMORA.InverterAssets.set_inverter_table!(scenario, [row])
    @test only(AIMORA.InverterAssets.inverter_table(scenario))[:id] == :inv1

    profile = AIMORA.StudyInputProfiles.input_profile(:power_flow)
    @test :buses in AIMORA.StudyInputs.required_keys(profile)
    @test :emt in (study.id for study in AIMORA.StudyCatalog.available_studies())
end

if AIMORA.solver_available()
    @testset "grounded scalar branch references" begin
        fixed_control = [
            "BEGIN NEW DATA CASE",
            "POWER FREQUENCY                     60.0",
            "ABSOLUTE U.M. DIMENSIONS              20       2      50      60",
            "C PRINTED NUMBER WIDTH  13  2",
            " .000200    .100",
            "       1       1       1       1       1      -1                               1",
            "       5       5      20      20     100     100",
            "C BRANCHES",
        ]
        fixed_terminators = [
            "BLANK card ending branch cards",
            "BLANK card ending nonexistent switch cards",
            "BLANK card ending all electric-network sources",
            "BLANK card terminating output variable requests",
            "BLANK card terminating plot cards",
            "BEGIN NEW DATA CASE",
        ]
        parsed = AIMORA.DeckParser.parse_deck_lines(
            vcat(
                fixed_control,
                [
                    "  BUSAS2                  1.0E+6",
                    "  BUSBS2      BUSAS2",
                    "  BUSCS2      BUSAS2",
                ],
                fixed_terminators,
            );
            source = "grounded-scalar-reference-contract",
        )
        @test AIMORA.ValidationCore.is_valid(parsed.validation)
        @test length(parsed.elements) == 3
        @test all(element -> element isa AIMORA.Branches.ConductanceBranch, parsed.elements)
        @test getfield.(parsed.elements, :g) == fill(1.0e-6, 3)
        @test get(
            parsed.card_counts,
            :fixed_grounded_scalar_branch_reference,
            0,
        ) == 2
        @test get(
            parsed.card_counts,
            :deferred_single_terminal_capacitance_group,
            0,
        ) == 0
        @test getfield.(parsed.over2_branch_rows, :layout_kind) == [
            :fixed_sparse_numeric,
            :fixed_grounded_scalar_branch_reference,
            :fixed_grounded_scalar_branch_reference,
        ]
        @test getfield.(parsed.over2_branch_rows, :reference_name) == [
            :none,
            :branch_fixed_1,
            :branch_fixed_1,
        ]

        missing_reference = AIMORA.DeckParser.parse_deck_lines(
            vcat(
                fixed_control,
                ["  BUSBS2      ABSENT"],
                fixed_terminators,
            );
            source = "missing-grounded-scalar-reference-contract",
        )
        @test !AIMORA.ValidationCore.is_valid(missing_reference.validation)
        @test any(
            issue -> occursin(
                "does not match a prior accepted scalar branch owner",
                issue.message,
            ),
            missing_reference.validation.issues,
        )
    end
end

@testset "public inverter model" begin
    rows = AIMORA.Inverter.simulate_inverter(t_end = 2.0e-4, dt = 20.0e-6)
    summary = AIMORA.Inverter.inverter_summary(rows)
    @test summary.samples == 11
    @test all(isfinite, rows[end])
end

@testset "transformer parameter studies" begin
    short_circuit_lines = [
        "XFORMER",
        "22.        700.",
        "139.4      13.6     2100.     12.       700.",
        "\$PUNCH",
        "BRANCH  NAME1 NAME2 NAME3 NAME4 NAME5 NAME6",
        "3.3        83.3",
        "132.8      250.     6.7       83.3",
        "66.4       56.8     5.1       18.96",
        "13.2       56.8     3.2       18.96",
        "BLANK card ending XFORMER cases",
    ]
    short_deck =
        AIMORA.TransformerParameterInput.parse_over41_transformer_parameter_lines(
            short_circuit_lines;
            source = "short-circuit-contract",
        )
    @test length(short_deck.cases) == 2
    short_study =
        AIMORA.TransformerParameterStudy.run_transformer_parameter_study(
            short_deck,
        )
    @test short_study.generated_branch_count == 5
    @test short_study.physical_checks_passed
    two_winding = short_study.case_results[1]
    three_winding = short_study.case_results[2]
    @test two_winding.admittance_matrix_s *
          two_winding.impedance_matrix_ohm ≈
          -I atol = 2.0e-12
    @test three_winding.admittance_matrix_s *
          three_winding.impedance_matrix_ohm ≈
          -I atol = 1.0e-10
    @test two_winding.impedance_matrix_ohm ≈
          transpose(two_winding.impedance_matrix_ohm) atol = 1.0e-12
    @test eigmin(Symmetric(real.(two_winding.impedance_matrix_ohm))) >=
          -1.0e-12
    @test two_winding.generated_branches[2].resistance_values_ohm[1] ≈
          -0.002028822586562 atol = 5.0e-13
    @test two_winding.generated_branches[2].
          inductance_or_reactance_values[1] ≈
          135.3359666262 atol = 5.0e-10
    @test getfield.(two_winding.generated_branches, :from_node) ==
          [Symbol(""), Symbol("")]
    @test getfield.(two_winding.generated_branches, :to_node) ==
          [Symbol(""), Symbol("")]
    @test getfield.(
        three_winding.generated_branches,
        :from_node,
    ) == [:NAME1, :NAME3, :NAME5]

    saturable_header =
        rpad("XFORMER", 32) * @sprintf("%8.1f%8.1f", 60.0, 0.0)
    saturable_case =
        @sprintf(
            "  %-6s%18s%6.2f%6.2f%6s%6.2f",
            "TRANSF",
            "",
            0.5,
            1.2,
            "",
            1000.0,
        )
    saturable_winding_1 =
        @sprintf(
            "%2d%-6s%-6s%12s%6.2f%6.2f%6.2f",
            1,
            "H1",
            "H0",
            "",
            0.2,
            5.0,
            132.0,
        )
    saturable_winding_2 =
        @sprintf(
            "%2d%-6s%-6s%12s%6.2f%6.2f%6.2f",
            2,
            "L1",
            "L0",
            "",
            0.1,
            2.0,
            33.0,
        )
    saturable_study =
        AIMORA.TransformerParameterStudy.run_transformer_parameter_study_lines(
            [
                saturable_header,
                saturable_case,
                saturable_winding_1,
                saturable_winding_2,
                "BLANK",
            ];
            source = "saturable-contract",
        )
    saturable = only(saturable_study.case_results)
    @test saturable.magnetizing_parallel_inductance_h ≈ 2.4
    @test saturable.branch_resistance_matrix_ohm ≈
          transpose(saturable.branch_resistance_matrix_ohm)
    @test saturable.branch_inductance_or_reactance_matrix ≈
          1000.0 .* saturable.physical_inductance_matrix_h
    @test eigmin(Symmetric(saturable.physical_inductance_matrix_h)) >=
          -1.0e-12
    @test saturable.physical_checks_passed

    bctran_lines = [
        "ACCESS MODULE BCTRAN",
        " 360.       .428      300.      135.73    .428      300.      135.73       1 3 1",
        "  1132.79056 .2054666   H-1         H-2         H-3",
        "  263.393059 .0742333   L-1         L-2         L-3",
        "  350.       .0822      T-1   T-2   T-2               T-1",
        " 1 20.        8.74      300.      7.3431941 300.       3 1",
        " 1 30.        8.68      76.       26.258183 300.",
        " 2 30.        5.31      76.       18.552824 300.",
        "BLANK card that terminates short-circuit test data of BCTRAN",
        "BEGIN NEW DATA CASE",
    ]
    bctran_study =
        AIMORA.TransformerParameterStudy.run_transformer_parameter_study_lines(
            bctran_lines;
            source = "bctran-contract",
        )
    bctran = only(bctran_study.case_results)
    @test bctran.input.output_representation == :reactance
    @test size(bctran.reactance_matrix_ohm) == (9, 9)
    @test bctran.winding_resistances_ohm ≈
          [0.2054666, 0.0742333, 0.0822]
    @test bctran.inverse_inductance_matrix_per_h[1, 1] ≈
          26.512692374898 atol = 5.0e-11
    @test bctran.inverse_inductance_matrix_per_h[2, 1] ≈
          -59.57848438329 atol = 5.0e-11
    @test bctran.reactance_matrix_ohm[1, 1] ≈
          41432.097487193 atol = 2.0e-8
    @test bctran.reactance_matrix_ohm[4, 1] ≈
          -0.0533106395901 atol = 5.0e-8
    @test only(bctran.magnetizing_shunts).self_resistance_ohm ≈
          55098.277352343 atol = 5.0e-8
    @test bctran.positive_pair_reconstruction_residual <= 1.0e-14
    @test bctran.zero_pair_reconstruction_residual <= 1.0e-14
    @test bctran.inverse_inductance_matrix_per_h ≈
          transpose(bctran.inverse_inductance_matrix_per_h) atol = 1.0e-12
    @test eigmin(Symmetric(bctran.inverse_inductance_matrix_per_h)) >
          0.0
    @test bctran.physical_checks_passed

    legacy_bctran_lines = copy(bctran_lines)
    legacy_bctran_lines[1] = "XFORMER, 44."
    legacy_bctran =
        only(
            AIMORA.TransformerParameterStudy.run_transformer_parameter_study_lines(
                legacy_bctran_lines;
                source = "legacy-bctran-contract",
            ).case_results,
        )
    @test legacy_bctran.reactance_matrix_ohm ≈
          bctran.reactance_matrix_ohm atol = 0.0

    report =
        AIMORA.TransformerParameterReport.transformer_parameter_report_text(
            bctran_study,
        )
    @test report ==
          AIMORA.TransformerParameterReport.transformer_parameter_report_text(
        bctran_study,
    )
    @test occursin("CASE 1 MULTIPHASE_TRANSFORMER", report)
    @test occursin("GENERATED_BRANCH_COUNT 9", report)
    mktempdir() do directory
        path = joinpath(directory, "transformer_parameters.txt")
        @test AIMORA.TransformerParameterReport.write_transformer_parameter_report(
            path,
            bctran_study,
        ) == abspath(path)
        @test read(path, String) == report
    end

    @test_throws ArgumentError AIMORA.TransformerParameterInput.
                               parse_over41_transformer_parameter_lines(
        bctran_lines[1:6],
    )
    @test_throws ArgumentError AIMORA.TransformerParameters.
                               TransformerShortCircuitCase(
        1,
        [132.0, 33.0],
        [10000.0],
        [1.0],
        [100.0],
        1.0,
        100.0,
        Tuple[],
    ) |>
                               AIMORA.TransformerParameters.
                               transformer_short_circuit_parameters
end

if AIMORA.solver_available()
    @testset "sampled line steady-state terminal admittance" begin
        propagation = AIMORA.Lines.line_weighting_samples(
            (2:7) .* 1.0e-6,
            [0.01, 0.15, 0.50, 0.20, 0.04, 0.01],
        )
        admittance = AIMORA.Lines.line_weighting_samples(
            (0:5) .* 1.0e-6,
            [0.01, 0.10, 0.40, 0.20, 0.05, 0.01],
        )
        coefficients = AIMORA.Lines.sampled_line_weighting_coefficients(
            propagation,
            admittance,
            1.0e-6,
            300.0;
            propagation_peak_index = 3,
            admittance_rise_index = 2,
        )
        line = AIMORA.Lines.sampled_frequency_dependent_line(
            1,
            2,
            coefficients,
        )
        terminal_admittance =
            AIMORA.Lines.sampled_line_steady_state_terminal_admittance(
                line,
                60.0,
            )
        @test size(terminal_admittance) == (2, 2)
        @test all(isfinite, terminal_admittance)
        @test terminal_admittance ≈ transpose(terminal_admittance)
        @test_throws ArgumentError AIMORA.Lines.
                                   sampled_line_steady_state_terminal_admittance(
            line,
            -60.0,
        )
    end

    @testset "retired input-converter dispositions" begin
        parsed = AIMORA.DeckParser.parse_deck_lines(
            [
                "BEGIN NEW DATA CASE",
                "CS",
                "CZ",
                "BLANK CARD TERMINATING THE CASE",
            ];
            source = "retired-input-converter-contract",
        )
        @test AIMORA.ValidationCore.is_valid(parsed.validation)
        auxiliary = AIMORA.EMTStudy.run_deck_auxiliary_studies(parsed)
        @test auxiliary.compatibility_exclusions == [
            :pre_m37_zinc_oxide_card_converter,
            :pre_m37_switch_pseudononlinear_card_converter,
        ]
        @test isempty(auxiliary.deferred_requests)
    end

    @testset "private solver integration" begin
        system = AIMORA.Nodal.NodalSystem(2, [
            AIMORA.Branches.TheveninSource(1, 1.0e9, _ -> 1.0),
            AIMORA.Branches.ConductanceBranch(1, 2, 1.0),
            AIMORA.Branches.ConductanceBranch(2, 0, 1.0),
        ])
        voltage = AIMORA.Nodal.solve_step!(system, 0.0, 20.0e-6)
        @test voltage[2] ≈ 0.5 atol = 1.0e-6
        @test AIMORA.solver_status().mode == :full_engine
    end
else
    @testset "public checkout has no solver source" begin
        @test AIMORA.solver_status().mode == :open_core
        @test !isdefined(AIMORA, :Nodal)
        @test_throws ErrorException AIMORA.require_solver()
    end
end

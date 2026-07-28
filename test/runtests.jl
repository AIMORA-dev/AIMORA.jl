using Test
using AIMORA

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

@testset "public inverter model" begin
    rows = AIMORA.Inverter.simulate_inverter(t_end = 2.0e-4, dt = 20.0e-6)
    summary = AIMORA.Inverter.inverter_summary(rows)
    @test summary.samples == 11
    @test all(isfinite, rows[end])
end

if AIMORA.solver_available()
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

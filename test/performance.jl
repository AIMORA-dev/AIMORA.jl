using Test
using AIMORA

function measure_nodal_hot_path(system, step_count::Int)
    dt_s = 20.0e-6
    AIMORA.Nodal.solve_step!(system, 0.0, dt_s)
    allocations = @allocated for step in 1:step_count
        AIMORA.Nodal.solve_step!(system, step * dt_s, dt_s)
    end
    elapsed_s = @elapsed for step in 1:step_count
        AIMORA.Nodal.solve_step!(system, step * dt_s, dt_s)
    end
    return (; allocations, elapsed_s, seconds_per_step = elapsed_s / step_count)
end

function nodal_hot_path_metrics(step_count::Int = 100_000)
    system = AIMORA.Nodal.NodalSystem(2, [
        AIMORA.Branches.TheveninSource(1, 1.0e9, _ -> 1.0),
        AIMORA.Branches.ConductanceBranch(1, 2, 1.0),
        AIMORA.Branches.ConductanceBranch(2, 0, 1.0),
    ])
    return measure_nodal_hot_path(system, step_count)
end

function coupled_inductive_step!(branch, admittance, rhs, voltage, dt_s)
    fill!(admittance, 0.0)
    fill!(rhs, 0.0)
    AIMORA.Branches.stamp!(admittance, rhs, branch, 0.0, dt_s)
    snapshot = AIMORA.Branches.branch_companion_snapshot(branch, voltage, dt_s)
    AIMORA.Branches.update!(branch, voltage, dt_s)
    return snapshot
end

function coupled_inductive_hot_path_metrics()
    branch = AIMORA.Branches.CoupledInductiveBranch(
        [1, 2],
        [0, 0],
        [-2.0 0.5; 0.5 -1.5],
        2.0 * pi * 60.0;
        series_resistance = 0.25,
    )
    branch.previous_current .= [0.4, -0.2]
    branch.previous_voltage .= [1.5, -0.5]
    admittance = zeros(2, 2)
    rhs = zeros(2)
    voltage = [2.0, -1.0]
    dt_s = 20.0e-6
    snapshot = coupled_inductive_step!(branch, admittance, rhs, voltage, dt_s)
    allocations = @allocated coupled_inductive_step!(
        branch,
        admittance,
        rhs,
        voltage,
        dt_s,
    )
    return (; allocations, snapshot)
end

metrics = nodal_hot_path_metrics()
@testset "measured nodal hot path" begin
    @info "Nodal hot-path performance" metrics.allocations metrics.elapsed_s metrics.seconds_per_step
    @test metrics.allocations == 0
    @test metrics.seconds_per_step <= 2.0e-6
end

coupled_metrics = coupled_inductive_hot_path_metrics()
@testset "measured coupled-inductive hot path" begin
    @info "Coupled-inductive hot-path performance" coupled_metrics.allocations
    @test coupled_metrics.allocations == 0
    @test isfinite(coupled_metrics.snapshot.branch_current)
end

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

metrics = nodal_hot_path_metrics()
@testset "measured nodal hot path" begin
    @info "Nodal hot-path performance" metrics.allocations metrics.elapsed_s metrics.seconds_per_step
    @test metrics.allocations == 0
    @test metrics.seconds_per_step <= 2.0e-6
end

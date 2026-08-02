using Test
using AIMORA
using LinearAlgebra
using Statistics

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

function coupled_inductive_hot_path_metrics(step_count::Int = 100_000)
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
    elapsed_s = @elapsed for step in 1:step_count
        voltage[1] = 2.0 + (step & 1) * eps(Float64)
        coupled_inductive_step!(branch, admittance, rhs, voltage, dt_s)
    end
    return (
        ;
        allocations,
        elapsed_s,
        seconds_per_step = elapsed_s / step_count,
        snapshot,
        branch,
        admittance,
        rhs,
    )
end

function _measure_matrix_vector_kernel!(
    kernel!,
    destination,
    matrix,
    source,
    iteration_count::Int,
)
    kernel!(destination, matrix, source)
    allocations = @allocated kernel!(destination, matrix, source)
    elapsed_samples_s = Float64[]
    for _ in 1:7
        elapsed_s = @elapsed for iteration in 1:iteration_count
            source[1] = 1.0 + (iteration & 1) * eps(Float64)
            kernel!(destination, matrix, source)
        end
        push!(elapsed_samples_s, elapsed_s)
    end
    return (
        ;
        allocations,
        median_elapsed_s = median(elapsed_samples_s),
        seconds_per_call = median(elapsed_samples_s) / iteration_count,
    )
end

function coupled_matrix_vector_hot_path_metrics(iteration_count::Int = 1_000_000)
    matrix = [-2.0 0.5; 0.5 -1.5]
    source = [1.0, -0.5]
    actual = zeros(2)
    reference = zeros(2)
    inferred_result = @inferred AIMORA.Branches._coupled_matrix_vector_mul!(
        actual,
        matrix,
        source,
    )
    mul!(reference, matrix, source)
    reference_error = maximum(abs.(actual .- reference))
    manual = _measure_matrix_vector_kernel!(
        AIMORA.Branches._coupled_matrix_vector_mul!,
        actual,
        matrix,
        source,
        iteration_count,
    )
    blas = _measure_matrix_vector_kernel!(
        mul!,
        reference,
        matrix,
        source,
        iteration_count,
    )
    return (
        ;
        inferred_result,
        actual,
        reference,
        reference_error,
        manual,
        blas,
        speedup = blas.median_elapsed_s / manual.median_elapsed_s,
    )
end

metrics = nodal_hot_path_metrics()
@testset "measured nodal hot path" begin
    @info "Nodal hot-path performance" metrics.allocations metrics.elapsed_s metrics.seconds_per_step
    @test metrics.allocations == 0
    @test metrics.seconds_per_step <= 2.0e-6
end

coupled_metrics = coupled_inductive_hot_path_metrics()
@testset "measured coupled-inductive hot path" begin
    @info "Coupled-inductive hot-path performance" coupled_metrics.allocations coupled_metrics.elapsed_s coupled_metrics.seconds_per_step
    @test coupled_metrics.allocations == 0
    @test coupled_metrics.seconds_per_step <= 2.0e-6
    @test coupled_metrics.snapshot.branch_current == 0.42765671552407436
    @test coupled_metrics.snapshot.history_current == 0.41072403255015133
    @test coupled_metrics.branch.last_current == [9.000000000000208, -4.999999999999971]
    @test coupled_metrics.admittance == [
        0.007525636877299119 -0.0018814092193247798
        -0.0018814092193247798 0.005644227657974339
    ]
    @test all(isfinite, coupled_metrics.rhs)
end

matrix_vector_metrics = coupled_matrix_vector_hot_path_metrics()
@testset "small coupled matrix-vector kernel" begin
    @info "Small coupled matrix-vector performance" matrix_vector_metrics.manual.seconds_per_call matrix_vector_metrics.blas.seconds_per_call matrix_vector_metrics.speedup
    @test matrix_vector_metrics.inferred_result === matrix_vector_metrics.actual
    @test matrix_vector_metrics.reference_error == 0.0
    @test matrix_vector_metrics.actual == matrix_vector_metrics.reference
    @test matrix_vector_metrics.manual.allocations == 0
    @test matrix_vector_metrics.manual.median_elapsed_s <=
          matrix_vector_metrics.blas.median_elapsed_s
    @test_throws DimensionMismatch AIMORA.Branches._coupled_matrix_vector_mul!(
        zeros(1),
        zeros(2, 2),
        zeros(2),
    )
    @test_throws DimensionMismatch AIMORA.Branches._coupled_matrix_vector_mul!(
        zeros(2),
        zeros(2, 1),
        zeros(2),
    )
end

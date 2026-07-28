#!/usr/bin/env julia

using Printf

ROOT = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "AIMORA.jl"))

using .AIMORA.Inverter
using .AIMORA.Figures

function main()
    output_dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "outputs")
    mkpath(output_dir)

    t0 = time()
    rows = simulate_inverter(t_end = 0.150, dt = 20e-6)
    elapsed = time() - t0
    summary = inverter_summary(rows; elapsed_s = elapsed)

    write_inverter_csv(joinpath(output_dir, "inverter_timeseries.csv"), rows)
    write_inverter_summary(joinpath(output_dir, "summary.json"), summary)
    write_inverter_power_svg(joinpath(output_dir, "inverter_power.svg"), rows)
    write_inverter_current_svg(joinpath(output_dir, "inverter_current.svg"), rows)

    @printf("Output: %s\n", abspath(output_dir))
    @printf("Engine: Julia single-model EMT timestep loop\n")
    @printf("Legacy Fortran in loop: no\n")
    @printf("Final P/Q: %.4f / %.4f pu\n", summary.final_p_pu, summary.final_q_pu)
    @printf("Simulation loop time: %.6f s\n", elapsed)
end

main()

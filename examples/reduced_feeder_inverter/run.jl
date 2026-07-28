#!/usr/bin/env julia

using Printf

ROOT = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "AIMORA.jl"))

using .AIMORA.EMTStudy
using .AIMORA.Figures

function main()
    output_dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "outputs")
    mkpath(output_dir)

    cfg = UnifiedEMTConfig(
        t_end_s = 0.150,
        dt_s = 20e-6,
        initial_bus_pu = 0.972,
        inverter_p0_pu = 0.35,
        inverter_p1_pu = 0.80,
        inverter_q_pu = 0.05,
    )

    t0 = time()
    rows = run_reduced_feeder_inverter(cfg)
    elapsed = time() - t0

    write_unified_csv(joinpath(output_dir, "reduced_feeder_timeseries.csv"), rows)
    write_unified_summary(joinpath(output_dir, "summary.json"), rows, cfg; elapsed_s = elapsed)
    write_unified_voltage_svg(joinpath(output_dir, "reduced_feeder_voltage.svg"), rows)
    write_unified_power_svg(joinpath(output_dir, "reduced_feeder_inverter_power.svg"), rows)

    last = rows[end]
    @printf("Output: %s\n", abspath(output_dir))
    @printf("Engine: Julia unified EMT timestep loop\n")
    @printf("Legacy Fortran in loop: no\n")
    @printf("Samples: %d, dt: %.1f us\n", length(rows), cfg.dt_s * 1e6)
    @printf("Final bus voltage: %.4f pu\n", last[2])
    @printf("Final inverter P/Q: %.4f / %.4f pu\n", last[6], last[7])
    @printf("Simulation loop time: %.6f s\n", elapsed)
end

main()

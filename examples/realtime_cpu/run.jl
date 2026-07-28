#!/usr/bin/env julia

using Printf

ROOT = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "AIMORA.jl"))

using .AIMORA.RealtimeLoop
using .AIMORA.Inverter

mutable struct BenchState
    x::InverterState
    params::InverterParams
    sink::Float64
end

function main()
    output_dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "outputs")
    mkpath(output_dir)

    # This benchmark measures whether the CPU can execute the model step under
    # a 20 us budget. It is a real-time feasibility check, not RTDS hardware.
    params = InverterParams(p_ref_0_pu = 0.35, p_ref_1_pu = 0.80, q_ref_pu = 0.05)
    state = BenchState(initial_inverter_state(params), params, 0.0)
    report = run_fixed_step(
        (s, t, dt) -> begin
            row = inverter_row(s.x, t, s.params)
            s.sink += row[6]
            s.x = step_inverter(s.x, t, dt, s.params)
            return nothing
        end;
        state = state,
        dt_s = 20e-6,
        duration_s = 0.050,
        realtime = false,
    )
    write_realtime_summary(joinpath(output_dir, "realtime_summary.json"), report)

    @printf("Output: %s\n", abspath(output_dir))
    @printf("Steps: %d, dt: %.1f us\n", report.steps, report.dt_s * 1e6)
    @printf("Mean step: %.3f us, max step: %.3f us, overruns: %d\n",
        report.mean_step_s * 1e6,
        report.max_step_s * 1e6,
        report.overruns,
    )
end

main()

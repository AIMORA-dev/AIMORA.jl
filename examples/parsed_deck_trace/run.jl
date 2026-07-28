#!/usr/bin/env julia

using Printf

ROOT = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "AIMORA.jl"))

using .AIMORA.EMTStudy
using .AIMORA.ReportArtifacts

const DECK_LINES = [
    # This is the accepted small Julia deck subset, not full BPA EMTP deck execution.
    "source src bus1 1.0e9 0.0 60.0 0.0 1.0",
    "resistor line bus1 bus2 1.0",
    "resistor load bus2 0 1.0",
    "time_switch tie bus2 bus3 0.00004 1.0 false 1.0e9 0.0",
    "resistor load3 bus3 0 1.0",
]

function main()
    output_dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "outputs")
    mkpath(output_dir)

    dt_s = 20.0e-6
    t_end_s = 100.0e-6
    trace = run_deck_emt(DECK_LINES; dt_s = dt_s, t_end_s = t_end_s, source = "parsed_deck_trace example")
    artifacts = write_deck_trace_report_artifacts(
        output_dir,
        trace;
        basename = "parsed_deck_trace",
        title = "Parsed Deck Trace",
    )
    summary_path = write_deck_trace_summary(joinpath(output_dir, "parsed_deck_trace_summary.json"), trace)
    result = deck_trace_result(trace)

    @printf("Output CSV: %s\n", abspath(artifacts.csv_path))
    @printf("Output summary: %s\n", abspath(summary_path))
    @printf("Output EMTP-style report: %s\n", abspath(artifacts.text_path))
    @printf("Output SVG: %s\n", abspath(artifacts.svg_path))
    @printf("Output report manifest: %s\n", abspath(artifacts.manifest_path))
    @printf("Output OVER20 finalization CSV: %s\n", abspath(artifacts.finalization_csv_path))
    @printf("Output OVER20 finalization report: %s\n", abspath(artifacts.finalization_text_path))
    @printf("Output OVER20 finalization manifest: %s\n", abspath(artifacts.finalization_manifest_path))
    @printf("Engine: Julia fixed-step parsed-deck subset runner\n")
    @printf("Legacy Fortran in loop: no\n")
    @printf("Full BPA EMTP deck execution: no\n")
    @printf("Result status: %s\n", String(result.status))
    @printf("Samples: %d, dt: %.1f us\n", result.quantities[:samples].value, trace.dt_s * 1.0e6)
    @printf("Final bus2 voltage: %.6f pu\n", final_voltage_pu(trace, :bus2))
    @printf("Final bus3 voltage: %.6f pu\n", final_voltage_pu(trace, :bus3))
end

main()

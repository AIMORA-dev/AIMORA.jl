#!/usr/bin/env julia

using Printf

ROOT = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "AIMORA.jl"))

using .AIMORA.EMTStudy
using .AIMORA.Lines

const DT = 20.0e-6

const DECK_LINES = [
    "source src source 1.0e9 0.0 60.0 0.0 1.0",
    "bergeron_line line source load 1.0 $(2.0 * DT) $(DT)",
    "resistor matched load 0 1.0",
]

function write_line_trace_csv(path::AbstractString, trace::DeckEMTTrace)
    mkpath(dirname(path))
    source_index = trace.node_map[:source]
    load_index = trace.node_map[:load]
    open(path, "w") do io
        println(io, "time_s,source_v_pu,load_v_pu")
        for sample in eachindex(trace.time_s)
            @printf(
                io,
                "%.9f,%.9f,%.9f\n",
                trace.time_s[sample],
                trace.voltage_pu[source_index, sample],
                trace.voltage_pu[load_index, sample],
            )
        end
    end
    return path
end

function main()
    output_dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "outputs")
    trace = run_deck_emt(DECK_LINES; dt_s = DT, t_end_s = 5.0 * DT, source = "line_models example")
    csv_path = write_line_trace_csv(joinpath(output_dir, "bergeron_line_trace.csv"), trace)
    point = frequency_dependent_line_point(0.0, 2.0e-3, 0.0, 0.5e-6, 10.0, 60.0)

    @printf("Output CSV: %s\n", abspath(csv_path))
    @printf("Engine: Julia fixed-step Bergeron line slice\n")
    @printf("Legacy Fortran in loop: no\n")
    @printf("Samples: %d, dt: %.1f us\n", length(trace.time_s), trace.dt_s * 1.0e6)
    @printf("Final load voltage: %.6f pu\n", final_voltage_pu(trace, :load))
    @printf("Frequency point Zc: %.6f%+.6fim\n", real(point.characteristic_impedance), imag(point.characteristic_impedance))
    @printf("Frequency point gamma: %.9f%+.9fim\n", real(point.propagation_constant), imag(point.propagation_constant))
end

main()

#!/usr/bin/env julia

using Printf

ROOT = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "AIMORA.jl"))

using .AIMORA.Companion
using .AIMORA.Machines

function write_transfer_csv(path::AbstractString, preview)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "request,kind,etac_index,value")
        for row in eachindex(preview.ismtac_requests)
            @printf(
                io,
                "%d,%s,%d,%.9f\n",
                preview.ismtac_requests[row],
                String(preview.request_kinds[row]),
                preview.etac_storage_indices[row],
                preview.request_values[row],
            )
        end
    end
    return path
end

function run_case()
    past_state = past_machine_history_state(
        histq = [11.0, 22.0, 33.0, 44.0],
        shp = [0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 2.0],
    )
    state = SynchronousMachineTACSInterfaceState(fill(-99.0, 8); lmset = 1)
    preview = synchronous_machine_tacs_transfer_from_past_update!(
        state,
        past_state,
        [-1, -8, -12, -13, 1, 5];
        cu_values = [1.25],
        cv1 = -0.75,
        a3 = 2.0,
        a4 = -1.0,
        c1 = 0.5,
        c2 = 1.5,
        c3 = -0.5,
        c4 = 2.5,
        elp_i26_21 = 0.25,
        elp_i75_4 = 2.0,
        numask = 2,
        n22 = 0,
    )
    return preview, state
end

function main()
    output_dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "outputs")
    preview, state = run_case()
    csv_path = write_transfer_csv(joinpath(output_dir, "machine_tacs_transfer.csv"), preview)

    @printf("Output CSV: %s\n", abspath(csv_path))
    @printf("Engine: Julia OVER16 synchronous-machine TACS interface slice\n")
    @printf("Legacy Fortran in loop: no\n")
    @printf("Transferred ETAC values: %d\n", preview.transfer_count)
    @printf("Final LMSET: %d\n", state.lmset)
    @printf("First transferred value: %.6f\n", preview.request_values[1])
    @printf("Last transferred value: %.6f\n", preview.request_values[end])
end

main()

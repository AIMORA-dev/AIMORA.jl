#!/usr/bin/env julia

using Printf

ROOT = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "AIMORA.jl"))

using .AIMORA.Nodal
using .AIMORA.Sources
using .AIMORA.Nonlinear
using .AIMORA.TACS

function write_trace_csv(path::AbstractString, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "step,time_s,control,load_v_pu,valve_closed,tacs_switch_closed,sat_inductor_current_pu")
        for row in rows
            @printf(
                io,
                "%d,%.9f,%.9f,%.9f,%s,%s,%.9f\n",
                row.step,
                row.time_s,
                row.control,
                row.load_v_pu,
                string(row.valve_closed),
                string(row.tacs_switch_closed),
                row.sat_inductor_current_pu,
            )
        end
    end
    return path
end

function run_case()
    dt = 20.0e-6
    block = FirstOrderTACSBlock(1.0, 5.0 * dt)
    control = Ref(0.0)
    valve = DiodeValveSwitch(1, 2; threshold_v = 0.2, holding_current = 0.0)
    tacs_switch = TACSControlledSwitch(2, 3, control; threshold = 0.5)
    saturable = SaturableInductorBranch(3, 0, 0.2, 1.0e-3, 0.2e-3, 1.0; i_prev = 0.25)
    system = NodalSystem(3, [
        constant_thevenin_source(1, 1.0e9, 1.0),
        valve,
        tacs_switch,
        ConductanceBranch(2, 0, 1.0e-9),
        ConductanceBranch(3, 0, 1.0),
        saturable,
    ])

    rows = NamedTuple[]
    for step in 1:10
        control[] = tacs_step!(block, 1.0, dt)
        voltages = solve_step!(system, (step - 1) * dt, dt)
        push!(rows, (;
            step,
            time_s = (step - 1) * dt,
            control = control[],
            load_v_pu = voltages[3],
            valve_closed = valve.closed,
            tacs_switch_closed = tacs_switch.closed,
            sat_inductor_current_pu = saturable.i_last,
        ))
    end
    return rows
end

function main()
    output_dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "outputs")
    rows = run_case()
    csv_path = write_trace_csv(joinpath(output_dir, "nonlinear_controls_trace.csv"), rows)
    last = rows[end]

    @printf("Output CSV: %s\n", abspath(csv_path))
    @printf("Engine: Julia fixed-step nonlinear/control EMT slice\n")
    @printf("Legacy Fortran in loop: no\n")
    @printf("Samples: %d, dt: %.1f us\n", length(rows), 20.0)
    @printf("Final control: %.6f\n", last.control)
    @printf("Final load voltage: %.6f pu\n", last.load_v_pu)
    @printf("Final saturable-inductor current: %.6f pu\n", last.sat_inductor_current_pu)
end

main()

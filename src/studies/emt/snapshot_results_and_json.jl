
function inverter_params(cfg::UnifiedEMTConfig, bus_pu::Float64)
    return InverterParams(
        v_ll_rms_v = cfg.v_ll_base_v * max(bus_pu, 0.2),
        p_ref_0_pu = cfg.inverter_p0_pu,
        p_ref_1_pu = cfg.inverter_p1_pu,
        q_ref_pu = cfg.inverter_q_pu,
    )
end

function make_network(cfg::UnifiedEMTConfig, inverter_current_pu::Base.RefValue{Float64})
    # Node 1: stiff feeder source behind a conductance.
    # Node 2: inverter/load bus.
    elements = [
        TheveninSource(1, cfg.source_stiffness_pu, t -> cfg.source_pu),
        SeriesRLBranch(1, 2, cfg.feeder_r_pu, cfg.feeder_l_pu_s),
        ConductanceBranch(2, 0, cfg.load_conductance_pu),
        CurrentInjection(2, t -> inverter_current_pu[]),
    ]
    return NodalSystem(2, elements)
end

function reduced_feeder_step_context(cfg::UnifiedEMTConfig, inverter_current_pu::Base.RefValue{Float64})
    network = make_network(cfg, inverter_current_pu)
    network.v[1] = cfg.source_pu
    network.v[2] = cfg.initial_bus_pu
    return initialize_step_context(
        network;
        node_map = Dict(:source => 1, :bus => 2),
        element_names = [:source, :feeder, :load, :inverter_current],
        dt_s = cfg.dt_s,
        t_end_s = cfg.t_end_s,
        source = "reduced feeder inverter",
    )
end

function run_reduced_feeder_inverter_with_trace(cfg::UnifiedEMTConfig = UnifiedEMTConfig())
    inverter_current_pu = Ref(0.0)
    context = reduced_feeder_step_context(cfg, inverter_current_pu)
    steps = context.step_count
    rows = Vector{NTuple{13, Float64}}(undef, steps + 1)

    bus_pu = cfg.initial_bus_pu
    p0 = inverter_params(cfg, bus_pu)
    inv_state = initial_inverter_state(p0)
    bus_index = context.node_map[:bus]
    source_index = context.node_map[:source]

    for n in 0:steps
        t = context.t_s

        inv = inverter_row(inv_state, t, inverter_params(cfg, bus_pu))
        for _ in 1:2
            inverter_current_pu[] = inv[6] / max(abs(bus_pu), 0.20)
            solved_v = solve_step!(context.system, t, context.dt_s)
            bus_pu = solved_v[bus_index]
            inv = inverter_row(inv_state, t, inverter_params(cfg, bus_pu))
        end

        p_inv, q_inv = inv[6], inv[7]
        source_node_pu = context.system.v[source_index]
        record_step!(context, context.system.v)

        rows[n + 1] = (
            t,
            bus_pu,
            source_node_pu,
            cfg.source_pu,
            cfg.load_conductance_pu,
            p_inv,
            q_inv,
            inv[2],
            inv[3],
            inv[4],
            inv[5],
            inv[8],
            inv[9],
        )

        if n < steps
            inv_state = step_inverter(inv_state, t, context.dt_s, inverter_params(cfg, bus_pu))
        end
    end
    return rows, deck_trace(context)
end

function run_reduced_feeder_inverter(cfg::UnifiedEMTConfig = UnifiedEMTConfig())
    rows, _ = run_reduced_feeder_inverter_with_trace(cfg)
    return rows
end

function reduced_feeder_result(rows, cfg::UnifiedEMTConfig; elapsed_s::Float64 = 0.0)
    last = rows[end]
    return study_result(
        :emt;
        status = :warning,
        quantities = [
            result_quantity(:samples, length(rows); unit = "count", description = "Recorded EMT samples."),
            result_quantity(:dt_s, cfg.dt_s; unit = "s", description = "Fixed EMT timestep."),
            result_quantity(:t_end_s, cfg.t_end_s; unit = "s", description = "Simulation duration."),
            result_quantity(:elapsed_s, elapsed_s; unit = "s", description = "Wall-clock runtime."),
            result_quantity(:final_bus_v_pu, last[2]; unit = "pu", base = "$(cfg.v_ll_base_v) V LL RMS", description = "Final reduced feeder load/inverter bus voltage."),
            result_quantity(:final_inverter_p_pu, last[6]; unit = "pu", description = "Final inverter active power."),
            result_quantity(:final_inverter_q_pu, last[7]; unit = "pu", description = "Final inverter reactive power."),
        ],
        assumptions = [
            study_assumption(:model, "reduced_feeder_with_inverter"; description = "Reduced two-node feeder prototype, not full IEEE 13 bus."),
            study_assumption(:external_reference_in_loop, false; description = "No external reference executable participates in the timestep loop."),
            study_assumption(:network_coupling, "scalar active-power current injection only"; description = "The inverter affects the network through scalar active current."),
            study_assumption(:timestep_context, "EMTStepContext"; description = "The reduced feeder network advances through the explicit fixed-step EMT context."),
        ],
        warnings = [
            study_warning(:reactive_power_not_coupled, "Reactive inverter power is reported but does not affect network voltage in this reduced prototype."),
        ],
        metadata = Dict{Symbol,Any}(
            :engine => "Julia unified EMT timestep engine",
            :source => "src/studies/emt.jl",
        ),
    )
end

function ensure_dir(path::AbstractString)
    isdir(path) || mkpath(path)
end

function write_unified_csv(path::AbstractString, rows)
    ensure_dir(dirname(path))
    open(path, "w") do io
        println(io, "time_s,bus_v_pu,source_node_v_pu,source_v_pu,load_conductance_pu,inverter_p_pu,inverter_q_pu,id_pu,iq_pu,id_ref_pu,iq_ref_pu,p_ref_pu,q_ref_pu")
        for r in rows
            @printf(io, "%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n", r...)
        end
    end
end

function write_deck_trace_summary(path::AbstractString, trace::DeckEMTTrace; elapsed_s::Float64 = 0.0)
    ensure_dir(dirname(path))
    result = deck_trace_result(trace; elapsed_s = elapsed_s)
    open(path, "w") do io
        println(io, "{")
        println(io, "  \"schema\": \"$(DECK_TRACE_SUMMARY_SCHEMA)\",")
        println(io, "  \"engine\": \"$(result.metadata[:engine])\",")
        println(io, "  \"source\": \"$(result.metadata[:source])\",")
        println(io, "  \"status\": \"$(result.status)\",")
        println(io, "  \"external_reference_in_loop\": $(result.assumptions[3].value),")
        println(io, "  \"admitted_deck_execution\": $(result.assumptions[4].value),")
        @printf(io, "  \"samples\": %d,\n", result.quantities[:samples].value)
        @printf(io, "  \"dt_s\": %.9f,\n", result.quantities[:dt_s].value)
        @printf(io, "  \"t_end_s\": %.9f,\n", result.quantities[:t_end_s].value)
        @printf(io, "  \"elapsed_s\": %.9f,\n", result.quantities[:elapsed_s].value)
        @printf(io, "  \"node_count\": %d,\n", result.quantities[:node_count].value)
        @printf(io, "  \"element_count\": %d,\n", result.quantities[:element_count].value)
        @printf(io, "  \"output_channel_count\": %d,\n", result.quantities[:output_channel_count].value)
        println(io, "  \"final_voltages_pu\": {")
        for (idx, node) in enumerate(trace.node_names)
            suffix = idx == length(trace.node_names) ? "" : ","
            @printf(io, "    \"%s\": %.9f%s\n", String(node), final_voltage_pu(trace, node), suffix)
        end
        println(io, "  },")
        println(io, "  \"final_outputs_pu\": {")
        for (idx, channel) in enumerate(trace.output_channel_names)
            suffix = idx == length(trace.output_channel_names) ? "" : ","
            @printf(io, "    \"%s\": %.9f%s\n", String(channel), final_output_pu(trace, channel), suffix)
        end
        println(io, "  }")
        println(io, "}")
    end
    return path
end

function write_unified_summary(
    path::AbstractString,
    rows,
    cfg::UnifiedEMTConfig;
    elapsed_s::Float64 = 0.0,
)
    ensure_dir(dirname(path))
    result = reduced_feeder_result(rows, cfg; elapsed_s = elapsed_s)
    open(path, "w") do io
        println(io, "{")
        println(io, "  \"engine\": \"$(result.metadata[:engine])\",")
        println(io, "  \"model\": \"$(result.assumptions[1].value)\",")
        println(io, "  \"status\": \"$(result.status)\",")
        println(io, "  \"external_reference_in_loop\": $(result.assumptions[2].value),")
        println(io, "  \"network_coupling\": \"$(result.assumptions[3].value)\",")
        println(io, "  \"timestep_context\": \"$(result.assumptions[4].value)\",")
        println(io, "  \"reactive_power_affects_network_voltage\": false,")
        @printf(io, "  \"samples\": %d,\n", result.quantities[:samples].value)
        @printf(io, "  \"dt_s\": %.9f,\n", result.quantities[:dt_s].value)
        @printf(io, "  \"t_end_s\": %.9f,\n", result.quantities[:t_end_s].value)
        @printf(io, "  \"elapsed_s\": %.9f,\n", result.quantities[:elapsed_s].value)
        @printf(io, "  \"final_bus_v_pu\": %.9f,\n", result.quantities[:final_bus_v_pu].value)
        @printf(io, "  \"final_inverter_p_pu\": %.9f,\n", result.quantities[:final_inverter_p_pu].value)
        @printf(io, "  \"final_inverter_q_pu\": %.9f\n", result.quantities[:final_inverter_q_pu].value)
        println(io, "}")
    end
    return path
end

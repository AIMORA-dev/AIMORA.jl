function deck_trace_result(trace::DeckEMTTrace; elapsed_s::Float64 = 0.0)
    quantities = ResultQuantity[
        result_quantity(:samples, length(trace.time_s); unit = "count", description = "Recorded fixed-step EMT samples."),
        result_quantity(:dt_s, trace.dt_s; unit = "s", description = "Fixed EMT timestep."),
        result_quantity(:t_end_s, trace.t_end_s; unit = "s", description = "Simulation duration."),
        result_quantity(:elapsed_s, elapsed_s; unit = "s", description = "Wall-clock runtime."),
        result_quantity(:node_count, length(trace.node_names); unit = "count", description = "Parsed deck node count."),
        result_quantity(:element_count, length(trace.element_names); unit = "count", description = "Parsed deck element count."),
        result_quantity(:output_channel_count, length(trace.output_channel_names); unit = "count", description = "Accepted non-voltage trace output channel count."),
    ]
    for node in trace.node_names
        push!(
            quantities,
            result_quantity(
                Symbol("final_", String(node), "_v_pu"),
                final_voltage_pu(trace, node);
                unit = "pu",
                description = "Final parsed-deck node voltage for $(String(node)).",
            ),
        )
    end
    for channel in trace.output_channel_names
        push!(
            quantities,
            result_quantity(
                Symbol("final_", String(channel)),
                final_output_pu(trace, channel);
                unit = "pu",
                description = "Final accepted parsed-deck output channel $(String(channel)).",
            ),
        )
    end

    return study_result(
        :emt;
        status = :ok,
        quantities = quantities,
        assumptions = [
            study_assumption(:model, "parsed_deck_emt"; description = "Typed Julia electromagnetic-transient deck execution."),
            study_assumption(:timestep_context, "EMTStepContext"; description = "The fixed-step Julia context owns timestep state and trace buffers."),
            study_assumption(:external_reference_in_loop, false; description = "No external reference executable participates in the study."),
            study_assumption(:admitted_deck_execution, true; description = "The accepted deck ran through its admitted Julia owners."),
        ],
        metadata = Dict{Symbol,Any}(
            :engine => "AIMORA Julia EMT",
            :source => trace.source,
            :source_path => "src/studies/emt.jl",
            :node_names => copy(trace.node_names),
            :element_names => copy(trace.element_names),
            :output_channel_names => copy(trace.output_channel_names),
        ),
    )
end

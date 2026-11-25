function _deck_control_system_expression_runtimes(
    parsed::DeckParser.DeckParseResult,
    state::ControlSystemExecutionState,
)
    rows = DeckParser.deck_control_system_expression_rows(parsed)
    runtimes = map(rows) do row
        program = compile_control_expression(
            row.name, row.expression; bpa_factor_arithmetic = true,
        )
        return ControlExpressionRuntime(program)
    end
    if any(runtime -> runtime.random_draw_count > 0, runtimes)
        DeckParser.deck_output_schedule_options(parsed).energization_count == 0 ||
            throw(ArgumentError("TACS RAN with statistical energization requires a shared switching stream"))
        # An output's reserved initial slot does not make its producing expression ready.
        pending_outputs = Set(runtime.program.output_name for runtime in runtimes)
        for runtime in runtimes
            program = runtime.program
            for name in program.input_names
                name in pending_outputs || continue
                throw(ArgumentError(
                    "TACS RAN requires producers before consumers: $(program.output_name) depends on $name",
                ))
            end
            delete!(pending_outputs, program.output_name)
        end
    end
    for (row, runtime) in zip(rows, runtimes)
        for input_name in runtime.program.input_names
            haskey(state.values, input_name) || throw(ArgumentError(
                "control expression $(row.name) on line $(row.line_no) references unknown signal $input_name",
            ))
        end
        state.values[row.name] = 0.0
    end
    return runtimes
end

function _advance_control_system_expressions!(
    runtimes::AbstractVector{ControlExpressionRuntime},
    state::ControlSystemExecutionState,
    step::Int,
)
    draws_per_step = sum(runtime -> runtime.random_draw_count, runtimes; init = 0)
    # Every expression runs once per accepted tick in deck order. Reconstructing
    # the stream from the absolute tick makes repeated trials and restart agree.
    stream = step > 0 && draws_per_step > 0 ?
        _bpa_control_random_stream(step, draws_per_step) : nothing
    for runtime in runtimes
        step == 0 && runtime.random_draw_count > 0 && continue
        state.values[runtime.program.output_name] =
            evaluate_control_expression!(runtime, state.values; bpa_random_stream = stream)
    end
    return length(runtimes)
end

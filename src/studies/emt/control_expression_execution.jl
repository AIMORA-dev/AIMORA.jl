function _deck_control_system_expression_runtimes(
    parsed::DeckParser.DeckParseResult,
    state::ControlSystemExecutionState,
)
    runtimes = ControlExpressionRuntime[]
    for row in DeckParser.deck_control_system_expression_rows(parsed)
        program = compile_control_expression(
            row.name, row.expression; bpa_factor_arithmetic = true,
        )
        for input_name in program.input_names
            haskey(state.values, input_name) ||
                throw(ArgumentError(
                    "control expression $(row.name) on line $(row.line_no) references unknown signal $input_name",
                ))
        end
        runtime = ControlExpressionRuntime(program)
        runtime.random_draw_count <= 1 || throw(ArgumentError(
            "multiple RAN calls in one BPA expression require legacy operand evaluation order",
        ))
        if runtime.random_draw_count > 0 &&
           DeckParser.deck_output_schedule_options(parsed).energization_count != 0
            throw(ArgumentError("TACS RAN with statistical energization requires a shared switching stream"))
        end
        state.values[row.name] = 0.0
        push!(runtimes, runtime)
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

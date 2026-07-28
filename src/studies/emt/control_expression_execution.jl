function _deck_control_system_expression_runtimes(
    parsed::DeckParser.DeckParseResult,
    state::ControlSystemExecutionState,
)
    runtimes = ControlExpressionRuntime[]
    for row in DeckParser.deck_control_system_expression_rows(parsed)
        program = compile_control_expression(row.name, row.expression)
        for input_name in program.input_names
            haskey(state.values, input_name) ||
                throw(ArgumentError(
                    "control expression $(row.name) on line $(row.line_no) references unknown signal $input_name",
                ))
        end
        runtime = ControlExpressionRuntime(program)
        state.values[row.name] = 0.0
        push!(runtimes, runtime)
    end
    return runtimes
end

function _advance_control_system_expressions!(
    runtimes::AbstractVector{ControlExpressionRuntime},
    state::ControlSystemExecutionState,
)
    for runtime in runtimes
        state.values[runtime.program.output_name] =
            evaluate_control_expression!(runtime, state.values)
    end
    return length(runtimes)
end

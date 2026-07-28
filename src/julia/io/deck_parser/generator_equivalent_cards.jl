mutable struct GeneratorEquivalentParseState
    header_line_numbers::Vector{Int}
    from_node_names::Vector{Symbol}
    to_node_names::Vector{Symbol}
    from_node_indices::Vector{Int}
    to_node_indices::Vector{Int}
    output_codes::Vector{Int}
    marker_line_no::Int
    stage::Symbol
    zero_mode_branches::Vector{DeckGeneratorEquivalentModalBranchRow}
    positive_mode_branches::Vector{DeckGeneratorEquivalentModalBranchRow}
end

function generator_equivalent_header_card(line::AbstractString)::Bool
    image = fixed_image(line)
    fixed_int_value(image, 1, 2) == 1 || return false
    fixed_float_value(image, 27, 32) == -6666.0 || return false
    return fixed_float_value(image, 33, 38) == -6666.0
end

function generator_equivalent_marker_card(line::AbstractString)::Bool
    image = fixed_image(line)
    isempty(fixed_field(image, 1, 2)) || return false
    uppercase(fixed_field(image, 3, 8)) == "BRANCH" || return false
    return uppercase(fixed_field(image, 9, 14)) == "ES"
end

function _generator_equivalent_phase_terminals!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
)
    from_node = fixed_field(image, 3, 8)
    to_node = fixed_field(image, 9, 14)
    if isempty(from_node)
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "expected frequency-dependent generator-equivalent phase node in columns 3-8",
            ),
        )
        return nothing
    end
    isempty(to_node) && (to_node = "0")
    return (
        from_name = Symbol(from_node),
        to_name = Symbol(to_node),
        from_index = node_id!(result, from_node),
        to_index = node_id!(result, to_node),
    )
end

function start_generator_equivalent_parse!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
)
    image = fixed_image(line)
    terminals = _generator_equivalent_phase_terminals!(result, image, line_no)
    terminals === nothing && return nothing
    record_card!(result, :fixed_field)
    record_card!(result, :generator_equivalent_header)
    return GeneratorEquivalentParseState(
        [line_no],
        [terminals.from_name],
        [terminals.to_name],
        [terminals.from_index],
        [terminals.to_index],
        [bpa_fixed_branch_output_code_from_image(image)],
        0,
        :phase_rows,
        DeckGeneratorEquivalentModalBranchRow[],
        DeckGeneratorEquivalentModalBranchRow[],
    )
end

function _parse_generator_equivalent_phase_row!(
    result::DeckParseResult,
    state::GeneratorEquivalentParseState,
    line::AbstractString,
    line_no::Int,
)
    image = fixed_image(line)
    expected_phase = length(state.header_line_numbers) + 1
    phase = fixed_int_value(image, 1, 2)
    if phase != expected_phase
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "frequency-dependent generator-equivalent phase rows must be consecutive; expected type $expected_phase",
            ),
        )
        return false
    end
    terminals = _generator_equivalent_phase_terminals!(result, image, line_no)
    terminals === nothing && return false
    push!(state.header_line_numbers, line_no)
    push!(state.from_node_names, terminals.from_name)
    push!(state.to_node_names, terminals.to_name)
    push!(state.from_node_indices, terminals.from_index)
    push!(state.to_node_indices, terminals.to_index)
    push!(state.output_codes, bpa_fixed_branch_output_code_from_image(image))
    record_card!(result, :fixed_field)
    record_card!(result, :generator_equivalent_phase_row)
    return true
end

function _generator_equivalent_modal_values!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
)
    image = fixed_image(line)
    values = Float64[]
    for (first_col, last_col, field) in (
        (1, 16, "resistance"),
        (17, 32, "inductance"),
        (33, 48, "capacitance"),
        (49, 64, "damping resistance"),
    )
        raw = fixed_field(image, first_col, last_col)
        if isempty(raw)
            push!(values, 0.0)
            continue
        end
        value = tryparse_deck_float(raw)
        if value === nothing || !isfinite(value)
            add_issue!(
                result.validation,
                invalid_value(
                    "line $line_no",
                    "generator-equivalent $field must be a finite 16-column numeric value",
                ),
            )
            return nothing
        end
        push!(values, Float64(value))
    end
    return values
end

function _push_generator_equivalent_modal_branch!(
    result::DeckParseResult,
    state::GeneratorEquivalentParseState,
    line::AbstractString,
    line_no::Int,
)
    values = _generator_equivalent_modal_values!(result, line, line_no)
    values === nothing && return false
    raw_resistance, raw_inductance, raw_capacitance, raw_damping = values
    sequence_kind = state.stage == :zero_mode ? :zero : :positive
    rows = state.stage == :zero_mode ?
        state.zero_mode_branches : state.positive_mode_branches
    if raw_resistance == 9999.0
        if isempty(rows)
            add_issue!(
                result.validation,
                missing_data(
                    "line $line_no",
                    "frequency-dependent generator-equivalent $sequence_kind mode requires at least one branch",
                ),
            )
            return false
        end
        state.stage = state.stage == :zero_mode ? :positive_mode : :complete
        record_card!(result, :generator_equivalent_modal_sentinel)
        return true
    end
    all(>=(0.0), values) || begin
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "generator-equivalent modal R/L/C/damping values must be nonnegative",
            ),
        )
        return false
    end
    any(>(0.0), values) || begin
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "generator-equivalent modal branch must contain a positive component",
            ),
        )
        return false
    end

    resistance = raw_resistance
    inductance = raw_inductance
    damping = raw_damping
    if inductance <= 0.0
        resistance += damping
        inductance = 0.0
        damping = 0.0
    end
    branch = GeneratorEquivalentModalBranch(
        resistance,
        fixed_card_branch_timestep_inductance(result, inductance),
        fixed_card_branch_timestep_capacitance(result, raw_capacitance),
        damping,
    )
    push!(
        rows,
        DeckGeneratorEquivalentModalBranchRow(
            line_no,
            sequence_kind,
            raw_resistance,
            raw_inductance,
            raw_capacitance,
            raw_damping,
            branch,
        ),
    )
    record_card!(result, :fixed_field)
    record_card!(result, :generator_equivalent_modal_branch)
    return true
end

function _finish_generator_equivalent_parse!(
    result::DeckParseResult,
    state::GeneratorEquivalentParseState,
)
    zero_branches = GeneratorEquivalentModalBranch[
        row.branch for row in state.zero_mode_branches
    ]
    positive_branches = GeneratorEquivalentModalBranch[
        row.branch for row in state.positive_mode_branches
    ]
    element = generator_equivalent_history_injection(
        state.from_node_indices,
        state.to_node_indices,
        zero_branches,
        positive_branches,
    )
    name = Symbol("generator_equivalent_", length(result.generator_equivalent_rows) + 1)
    push!(result.elements, element)
    push!(result.element_names, name)
    push!(result.element_line_numbers, first(state.header_line_numbers))
    push!(
        result.generator_equivalent_rows,
        DeckGeneratorEquivalentRow(
            name,
            copy(state.header_line_numbers),
            state.marker_line_no,
            copy(state.from_node_names),
            copy(state.to_node_names),
            copy(state.from_node_indices),
            copy(state.to_node_indices),
            copy(state.zero_mode_branches),
            copy(state.positive_mode_branches),
            copy(state.output_codes),
        ),
    )
    record_card!(result, :generator_equivalent)
    return nothing
end

function parse_generator_equivalent_card!(
    result::DeckParseResult,
    state::GeneratorEquivalentParseState,
    line::AbstractString,
    line_no::Int,
)
    if state.stage == :phase_rows
        if generator_equivalent_marker_card(line)
            if length(state.header_line_numbers) < 2
                add_issue!(
                    result.validation,
                    missing_data(
                        "line $line_no",
                        "frequency-dependent generator equivalent requires at least two phase rows",
                    ),
                )
                return true
            end
            state.marker_line_no = line_no
            state.stage = :zero_mode
            record_card!(result, :generator_equivalent_modal_marker)
            return false
        end
        _parse_generator_equivalent_phase_row!(result, state, line, line_no)
        return false
    end
    _push_generator_equivalent_modal_branch!(result, state, line, line_no)
    if state.stage == :complete
        _finish_generator_equivalent_parse!(result, state)
        return true
    end
    return false
end

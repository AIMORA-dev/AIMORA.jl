
function parse_line_constants_frequency_card!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
)::Bool
    image = fixed_image(line)
    initial_issues = length(result.validation.issues)
    earth_resistivity = line_constants_float_or_default!(
        result,
        image,
        line_no,
        1,
        8,
        2,
        "line_constants_earth_resistivity",
        0.0,
    )
    frequency = line_constants_float_or_default!(
        result,
        image,
        line_no,
        9,
        18,
        2,
        "line_constants_frequency",
        0.0,
    )
    correction = line_constants_float_or_default!(
        result,
        image,
        line_no,
        19,
        28,
        6,
        "line_constants_carson_correction_factor",
        0.0,
    )
    if line_constants_blank(image, 19, 28) && correction == 0.0
        correction = 1.0e-6
    end
    capacitance_flags = (
        line_constants_int_or_default!(result, image, line_no, 29, 30, "line_constants_capacitance_print_flag_1", 0),
        line_constants_int_or_default!(result, image, line_no, 31, 31, "line_constants_capacitance_print_flag_2", 0),
        line_constants_int_or_default!(result, image, line_no, 32, 32, "line_constants_capacitance_print_flag_3", 0),
        line_constants_int_or_default!(result, image, line_no, 33, 33, "line_constants_capacitance_print_flag_4", 0),
        line_constants_int_or_default!(result, image, line_no, 34, 34, "line_constants_capacitance_print_flag_5", 0),
        line_constants_int_or_default!(result, image, line_no, 35, 35, "line_constants_capacitance_print_flag_6", 0),
    )
    impedance_flags = (
        line_constants_int_or_default!(result, image, line_no, 36, 37, "line_constants_impedance_print_flag_1", 0),
        line_constants_int_or_default!(result, image, line_no, 38, 38, "line_constants_impedance_print_flag_2", 0),
        line_constants_int_or_default!(result, image, line_no, 39, 39, "line_constants_impedance_print_flag_3", 0),
        line_constants_int_or_default!(result, image, line_no, 40, 40, "line_constants_impedance_print_flag_4", 0),
        line_constants_int_or_default!(result, image, line_no, 41, 41, "line_constants_impedance_print_flag_5", 0),
        line_constants_int_or_default!(result, image, line_no, 42, 42, "line_constants_impedance_print_flag_6", 0),
    )
    matrix_selector = line_constants_int_or_default!(
        result,
        image,
        line_no,
        43,
        44,
        "line_constants_matrix_output_selector",
        0,
    )
    distance = line_constants_float_or_default!(
        result,
        image,
        line_no,
        45,
        52,
        3,
        "line_constants_distance",
        0.0,
    )
    punch_request = line_constants_int_or_default!(
        result,
        image,
        line_no,
        53,
        54,
        "line_constants_punch_request",
        0,
    )
    alternate_flags = (
        line_constants_int_or_default!(result, image, line_no, 55, 55, "line_constants_alternate_punch_flag_1", 0),
        line_constants_int_or_default!(result, image, line_no, 56, 56, "line_constants_alternate_punch_flag_2", 0),
        line_constants_int_or_default!(result, image, line_no, 57, 57, "line_constants_alternate_punch_flag_3", 0),
        line_constants_int_or_default!(result, image, line_no, 58, 58, "line_constants_segmentation_flag", 0),
        line_constants_int_or_default!(result, image, line_no, 59, 59, "line_constants_mutual_flag", 0),
    )
    frequency_decade_count = line_constants_int_or_default!(
        result,
        image,
        line_no,
        60,
        62,
        "line_constants_frequency_decade_count",
        0,
    )
    points_per_decade = line_constants_int_or_default!(
        result,
        image,
        line_no,
        63,
        65,
        "line_constants_points_per_decade",
        0,
    )
    line_model_punch_request = line_constants_int_or_default!(
        result,
        image,
        line_no,
        66,
        68,
        "line_constants_line_model_punch_request",
        0,
    )
    modal_output_flag = line_constants_int_or_default!(
        result,
        image,
        line_no,
        69,
        70,
        "line_constants_modal_output_flag",
        0,
    )
    transform_output_flag = line_constants_int_or_default!(
        result,
        image,
        line_no,
        71,
        72,
        "line_constants_transform_output_flag",
        0,
    )
    conductance = line_constants_float_or_default!(
        result,
        image,
        line_no,
        73,
        80,
        2,
        "line_constants_conductance",
        0.0,
    )
    if conductance == 0.0
        conductance = 3.22e-9
    end
    if length(result.validation.issues) != initial_issues
        record_fixed_blocker!(
            result,
            :line_constants_frequency_card_blocked,
            :line_constants_frequency_card_parse_error,
        )
        return true
    end
    if earth_resistivity == 0.0
        return false
    end

    push!(
        result.line_constants_frequency_cards,
        DeckLineConstantsFrequencyCard(
            line_no,
            earth_resistivity,
            frequency,
            correction,
            capacitance_flags,
            impedance_flags,
            matrix_selector,
            distance,
            punch_request,
            alternate_flags,
            frequency_decade_count,
            points_per_decade,
            line_model_punch_request,
            modal_output_flag,
            transform_output_flag,
            conductance,
        ),
    )
    record_card!(result, :fixed_field)
    record_card!(result, :line_constants_frequency_card)
    return true
end

const EMPTY_SATURATED_TRANSFORMER_REFERENCE = Symbol("")

function record_saturated_transformer_card!(counts::Dict{Symbol,Int}, card::Symbol)
    counts[card] = get(counts, card, 0) + 1
    return counts
end

function saturated_transformer_float_or_issue!(
    validation::ValidationResult,
    raw::AbstractString,
    line_no::Int,
    field::AbstractString,
)
    value = tryparse_deck_float(raw)
    if value === nothing
        add_issue!(
            validation,
            invalid_value("line $line_no", "$field=$(String(raw)): expected Float64"),
        )
    end
    return value
end

function saturated_transformer_symbol(raw::AbstractString)
    return Symbol(strip(String(raw)))
end

function parse_saturated_transformer_header_tokens!(
    transformers::Vector{DeckSaturatedTransformerHeaderRow},
    card_counts::Dict{Symbol,Int},
    validation::ValidationResult,
    tokens,
    line_no::Int,
)
    length(tokens) >= 4 || begin
        add_issue!(
            validation,
            missing_data("line $line_no", "expected TRANSFORMER current, flux, and name fields"),
        )
        return nothing
    end
    current = saturated_transformer_float_or_issue!(validation, tokens[2], line_no, "initial_current")
    flux = saturated_transformer_float_or_issue!(validation, tokens[3], line_no, "initial_flux")
    current === nothing && return nothing
    flux === nothing && return nothing
    name = saturated_transformer_symbol(tokens[4])
    isempty(String(name)) && begin
        add_issue!(validation, missing_data("line $line_no", "expected transformer name"))
        return nothing
    end
    magnetizing_resistance = missing
    if length(tokens) >= 5
        parsed = saturated_transformer_float_or_issue!(
            validation,
            tokens[5],
            line_no,
            "magnetizing_resistance",
        )
        parsed === nothing && return nothing
        magnetizing_resistance = parsed
    end
    push!(
        transformers,
        DeckSaturatedTransformerHeaderRow(
            name,
            EMPTY_SATURATED_TRANSFORMER_REFERENCE,
            line_no,
            current,
            flux,
            magnetizing_resistance,
        ),
    )
    record_saturated_transformer_card!(card_counts, :saturated_transformer_header)
    record_saturated_transformer_card!(card_counts, :saturated_transformer_intake)
    return name
end

function parse_saturated_transformer_copy_tokens!(
    transformers::Vector{DeckSaturatedTransformerHeaderRow},
    copies::Vector{DeckSaturatedTransformerCopyRow},
    card_counts::Dict{Symbol,Int},
    validation::ValidationResult,
    line::AbstractString,
    tokens,
    line_no::Int,
)
    image = fixed_image(line)
    reference_text = fixed_field(image, 15, 20)
    target_text = fixed_field(image, 39, 44)
    if isempty(reference_text) && length(tokens) >= 2
        reference_text = String(tokens[2])
    end
    if isempty(target_text) && length(tokens) >= 3
        target_text = String(tokens[3])
    end
    if isempty(reference_text) || isempty(target_text)
        add_issue!(
            validation,
            missing_data("line $line_no", "expected TRANSFORMER copy reference and target names"),
        )
        return nothing
    end
    reference = saturated_transformer_symbol(reference_text)
    target = saturated_transformer_symbol(target_text)
    push!(
        copies,
        DeckSaturatedTransformerCopyRow(reference, target, line_no),
    )
    push!(
        transformers,
        DeckSaturatedTransformerHeaderRow(
            target,
            reference,
            line_no,
            missing,
            missing,
            missing,
        ),
    )
    record_saturated_transformer_card!(card_counts, :saturated_transformer_copy)
    record_saturated_transformer_card!(card_counts, :saturated_transformer_intake)
    return target
end

function parse_saturated_transformer_breakpoint_tokens!(
    breakpoints::Vector{DeckSaturatedTransformerBreakpointRow},
    card_counts::Dict{Symbol,Int},
    validation::ValidationResult,
    active_name::Symbol,
    tokens,
    line_no::Int,
)
    length(tokens) >= 2 || begin
        add_issue!(
            validation,
            missing_data("line $line_no", "expected saturated transformer current and flux"),
        )
        return nothing
    end
    current = saturated_transformer_float_or_issue!(validation, tokens[1], line_no, "current")
    flux = saturated_transformer_float_or_issue!(validation, tokens[2], line_no, "flux")
    current === nothing && return nothing
    flux === nothing && return nothing
    if !isempty(breakpoints)
        previous = breakpoints[end]
        if previous.transformer_name == active_name
            if current <= previous.current
                add_issue!(
                    validation,
                    invalid_value(
                        "line $line_no",
                        "saturated transformer current must increase monotonically",
                    ),
                )
            end
            if flux <= previous.flux
                add_issue!(
                    validation,
                    invalid_value(
                        "line $line_no",
                        "saturated transformer flux must increase monotonically",
                    ),
                )
            end
        end
    end
    push!(breakpoints, DeckSaturatedTransformerBreakpointRow(active_name, line_no, current, flux))
    record_saturated_transformer_card!(card_counts, :saturated_transformer_breakpoint)
    record_saturated_transformer_card!(card_counts, :saturated_transformer_intake)
    return active_name
end

function parse_saturated_transformer_winding_line!(
    windings::Vector{DeckSaturatedTransformerWindingRow},
    card_counts::Dict{Symbol,Int},
    active_name::Symbol,
    line::AbstractString,
    line_no::Int,
)
    image = fixed_image(line)
    winding_number = fixed_int_value(image, 1, 2)
    winding_number === nothing && return false
    1 <= winding_number <= 9 || return false
    from_node_text = fixed_field(image, 3, 8)
    isempty(from_node_text) && return false
    to_node_text = fixed_field(image, 9, 14)
    resistance = fixed_float_value(image, 27, 32)
    inductance = fixed_float_value(image, 33, 38)
    turns = fixed_float_value(image, 39, 44)
    inherited = resistance === nothing && inductance === nothing && turns === nothing
    stored_resistance = inherited ? missing : something(resistance, 0.0)
    stored_inductance = inherited ? missing : something(inductance, 0.0)
    stored_turns = inherited ? missing : something(turns, 0.0)
    push!(
        windings,
        DeckSaturatedTransformerWindingRow(
            active_name,
            line_no,
            winding_number,
            saturated_transformer_symbol(from_node_text),
            saturated_transformer_symbol(to_node_text),
            stored_resistance,
            stored_inductance,
            stored_turns,
            inherited,
        ),
    )
    record_saturated_transformer_card!(card_counts, :saturated_transformer_winding)
    record_saturated_transformer_card!(card_counts, :saturated_transformer_intake)
    return true
end

function parse_saturated_transformer_intake_lines(lines; source::AbstractString="deck")
    source_text = String(source)
    validation = validation_result(source=source_text)
    transformers = DeckSaturatedTransformerHeaderRow[]
    breakpoints = DeckSaturatedTransformerBreakpointRow[]
    windings = DeckSaturatedTransformerWindingRow[]
    copies = DeckSaturatedTransformerCopyRow[]
    card_counts = Dict{Symbol,Int}()
    active_name = nothing
    reading_breakpoints = false

    for (line_no, raw_line) in enumerate(lines)
        line = strip_deck_line(raw_line)
        isempty(line) && continue
        tokens = deck_tokens(line)
        isempty(tokens) && continue
        first_token = normalized_deck_token(tokens[1])

        if first_token == "transformer"
            header_current = length(tokens) >= 2 ? tryparse_deck_float(tokens[2]) : nothing
            if header_current !== nothing
                parsed_name = parse_saturated_transformer_header_tokens!(
                    transformers,
                    card_counts,
                    validation,
                    tokens,
                    line_no,
                )
                active_name = parsed_name
                reading_breakpoints = parsed_name !== nothing
            else
                parsed_name = parse_saturated_transformer_copy_tokens!(
                    transformers,
                    copies,
                    card_counts,
                    validation,
                    line,
                    tokens,
                    line_no,
                )
                active_name = parsed_name
                reading_breakpoints = false
            end
            continue
        end

        if reading_breakpoints && active_name !== nothing
            current = tryparse_deck_float(tokens[1])
            current === nothing && continue
            if current == 9999.0
                record_saturated_transformer_card!(
                    card_counts,
                    :saturated_transformer_breakpoint_termination,
                )
                record_saturated_transformer_card!(card_counts, :saturated_transformer_intake)
                reading_breakpoints = false
                continue
            end
            parse_saturated_transformer_breakpoint_tokens!(
                breakpoints,
                card_counts,
                validation,
                active_name,
                tokens,
                line_no,
            )
            continue
        end

        active_name === nothing && continue
        parse_saturated_transformer_winding_line!(
            windings,
            card_counts,
            active_name,
            line,
            line_no,
        )
    end

    return DeckSaturatedTransformerIntake(
        source_text,
        transformers,
        breakpoints,
        windings,
        copies,
        card_counts,
        validation,
    )
end

function parse_saturated_transformer_intake_file(path::AbstractString)
    return parse_saturated_transformer_intake_lines(readlines(path); source=String(path))
end

function saturated_transformer_branch_section_lines(lines)
    selected = String[]
    in_transformer_branch = false
    for (line_no, raw_line) in enumerate(lines)
        line = strip_deck_line(raw_line)
        tokens = deck_tokens(line)
        section_marker = fixed_card_section_marker(tokens)
        named_section_marker = fixed_card_named_section_marker(tokens)
        if in_transformer_branch &&
           (section_marker !== nothing || named_section_marker !== nothing)
            break
        end
        if !in_transformer_branch &&
           fixed_card_saturated_transformer_header_card(tokens)
            append!(selected, fill("", line_no - 1))
            in_transformer_branch = true
        end
        in_transformer_branch || continue
        push!(selected, String(raw_line))
    end
    return selected
end

function parse_saturated_transformer_branch_section_intake_lines(
    lines;
    source::AbstractString="deck",
)
    return parse_saturated_transformer_intake_lines(
        saturated_transformer_branch_section_lines(lines);
        source = source,
    )
end

function parse_saturated_transformer_branch_section_intake_file(path::AbstractString)
    return parse_saturated_transformer_branch_section_intake_lines(
        readlines(path);
        source = String(path),
    )
end

function parse_saturated_transformer_branch_section_shunt_capacitance_rows(
    lines;
    source::AbstractString="deck",
)
    rows = DeckTransformerBranchShuntCapacitanceRow[]
    in_transformer_branch = false
    for (line_no, raw_line) in enumerate(lines)
        line = strip_deck_line(raw_line)
        image = fixed_image(line)
        tokens = deck_tokens(line)
        section_marker = fixed_card_section_marker(tokens)
        named_section_marker = fixed_card_named_section_marker(tokens)
        if in_transformer_branch &&
           (section_marker !== nothing || named_section_marker !== nothing)
            break
        end
        if !in_transformer_branch &&
           fixed_card_saturated_transformer_header_card(tokens)
            in_transformer_branch = true
        end
        in_transformer_branch || continue
        isempty(line) && continue
        from_node = fixed_field(image, 3, 8)
        to_node = fixed_field(image, 9, 14)
        !isempty(from_node) && isempty(to_node) && (to_node = "0")
        aux_1 = fixed_field(image, 15, 20)
        aux_2 = fixed_field(image, 21, 26)
        single_terminal_capacitance_row(image, from_node, to_node, aux_1, aux_2) ||
            continue
        capacitance = fixed_float_value(image, 39, 44)
        capacitance === nothing && continue
        push!(
            rows,
            DeckTransformerBranchShuntCapacitanceRow(
                Symbol(from_node),
                line_no,
                Float64(capacitance),
            ),
        )
    end
    return rows
end

function bpa_fixed_branch_triplet!(result::DeckParseResult, image::AbstractString,
                                   line_no::Int)
    compact = (
        fixed_float_value(image, 27, 32),
        fixed_float_value(image, 33, 38),
        fixed_float_value(image, 39, 44),
    )
    if all(value -> value !== nothing, compact)
        return compact[1], compact[2], compact[3], :bpa_fixed_branch_compact
    end

    wide = (
        fixed_float_value(image, 27, 42),
        fixed_float_value(image, 43, 58),
        fixed_float_value(image, 59, 74),
    )
    if all(value -> value !== nothing, wide)
        return wide[1], wide[2], wide[3], :bpa_fixed_branch_wide
    end

    implied_decimal = fixed_card_implied_decimal_triplet(image)
    if implied_decimal !== nothing
        return (
            implied_decimal[1],
            implied_decimal[2],
            implied_decimal[3],
            :fixed_card_branch_implied_decimal_layout,
        )
    end

    compact_present = [value !== nothing for value in compact]
    if any(compact_present)
        return (
            compact_present[1] ? Float64(compact[1]) : 0.0,
            compact_present[2] ? Float64(compact[2]) : 0.0,
            compact_present[3] ? Float64(compact[3]) : 0.0,
            :fixed_card_branch_sparse_numeric_layout,
        )
    end

    whitespace_values = branch_numeric_tail_values(image)
    if length(whitespace_values) >= 2
        return (
            whitespace_values[1],
            whitespace_values[2],
            length(whitespace_values) >= 3 ? whitespace_values[3] : 0.0,
            :fixed_card_branch_sparse_numeric_layout,
        )
    end

    add_issue!(result.validation,
               invalid_value("line $line_no",
                             "expected OVER2 fixed-field R/L/C in compact columns 27-44 or wide columns 27-74"))
    return nothing, nothing, nothing, :bpa_fixed_branch_unknown_layout
end

function fixed_card_implied_decimal_value(raw::AbstractString)
    text = replace(String(raw), ' ' => "")
    isempty(text) && return 0.0
    text in (".", "+", "-", "+.", "-.") && return 0.0
    value = tryparse_deck_float(text)
    value === nothing && return nothing
    uppercase_text = uppercase(text)
    return occursin('.', text) || occursin('E', uppercase_text) ?
        Float64(value) : Float64(value) / 100.0
end

function fixed_card_implied_decimal_triplet(image::AbstractString)
    raw_fields = (
        image[27:32],
        image[33:38],
        image[39:44],
    )
    all(field -> isempty(strip(field)), raw_fields) && return nothing
    requires_implied_decimal = any(raw_fields) do field
        text = strip(field)
        isempty(text) && return false
        compact = replace(text, ' ' => "")
        isolated_sign_or_decimal = compact in (".", "+", "-", "+.", "-.")
        split_token = occursin(' ', text)
        lacks_explicit_decimal =
            !occursin('.', compact) && !occursin('E', uppercase(compact))
        return isolated_sign_or_decimal || split_token || lacks_explicit_decimal
    end
    requires_implied_decimal || return nothing
    values = map(fixed_card_implied_decimal_value, raw_fields)
    any(value -> value === nothing, values) && return nothing
    return Float64(values[1]), Float64(values[2]), Float64(values[3])
end

function fixed_card_branch_triplet_slices(image::AbstractString)
    values = Float64[]
    for first_col in 27:6:75
        value = fixed_float_value(image, first_col, min(first_col + 5, 80))
        value === nothing && continue
        push!(values, Float64(value))
    end
    triplet_count = length(values) ÷ 3
    return (
        resistance = Float64[values[3 * index - 2] for index in 1:triplet_count],
        inductance = Float64[values[3 * index - 1] for index in 1:triplet_count],
        capacitance = Float64[values[3 * index] for index in 1:triplet_count],
    )
end

function fixed_card_compact_matrix_triplet_slices(image::AbstractString)
    tail = fixed_field(image, 27, 80)
    if occursin('E', uppercase(tail)) || occursin('D', uppercase(tail))
        return (resistance = Float64[], inductance = Float64[], capacitance = Float64[])
    end
    return fixed_card_branch_triplet_slices(image)
end

function coupled_phase_pi_numeric_continuation_row(image::AbstractString)::Bool
    isempty(fixed_field(image, 1, 26)) || return false
    slices = fixed_card_compact_matrix_triplet_slices(image)
    return !isempty(slices.resistance)
end

function append_coupled_phase_pi_continuation_row!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    initial_issues::Int,
)::Bool
    row_index = findlast(
        row -> length(row.raw_resistance_values) < row.phase_index,
        result.coupled_phase_pi_section_rows,
    )
    if row_index === nothing
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "OVER2 coupled phase PI continuation row has no incomplete owner row",
            ),
        )
        return true
    end
    row = result.coupled_phase_pi_section_rows[row_index]
    slices = fixed_card_branch_triplet_slices(image)
    append!(row.raw_resistance_values, slices.resistance)
    append!(row.raw_inductance_values, slices.inductance)
    append!(row.raw_capacitance_values, slices.capacitance)
    push!(row.continuation_line_numbers, line_no)
    record_fixed_card!(
        result,
        :bpa_fixed_branch,
        :bpa_fixed_coupled_phase_pi_continuation,
        initial_issues,
    )
    return true
end

function branch_numeric_tail_values(image::AbstractString)
    tail = fixed_field(image, 27, 74)
    fields = split(tail)
    length(fields) >= 2 || return Float64[]
    values = Float64[]
    for field in fields
        value = tryparse_deck_float(field)
        value === nothing && continue
        push!(values, Float64(value))
    end
    if length(values) >= 3 && isinteger(values[end])
        pop!(values)
    end
    return values
end

ground_terminal_text(node::AbstractString)::Bool =
    lowercase(String(node)) in ("", "0", "gnd", "ground", "ref")

function single_terminal_capacitance_row(image::AbstractString,
                                         from_node::AbstractString,
                                         to_node::AbstractString,
                                         aux_1::AbstractString,
                                         aux_2::AbstractString)::Bool
    !isempty(from_node) || return false
    ground_terminal_text(to_node) || return false
    isempty(aux_1) || return false
    isempty(aux_2) || return false
    compact = (
        fixed_float_value(image, 27, 32),
        fixed_float_value(image, 33, 38),
        fixed_float_value(image, 39, 44),
    )
    compact[1] === nothing && compact[2] === nothing && compact[3] !== nothing ||
        return false
    return true
end

function grounded_scalar_branch_reference_row(image::AbstractString,
                                              from_node::AbstractString,
                                              to_node::AbstractString,
                                              aux_1::AbstractString,
                                              aux_2::AbstractString)::Bool
    !isempty(from_node) || return false
    ground_terminal_text(to_node) || return false
    isempty(aux_1) && return false
    isempty(aux_2) || return false
    return all(
        value -> value === nothing,
        (
            fixed_float_value(image, 27, 32),
            fixed_float_value(image, 33, 38),
            fixed_float_value(image, 39, 44),
            fixed_float_value(image, 45, 50),
        ),
    )
end

function parse_fixed_grounded_scalar_branch_reference_row!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    from_node::AbstractString,
    to_node::AbstractString,
    reference_node::AbstractString,
    initial_issues::Int,
)::Bool
    reference_index =
        fixed_branch_reference_by_existing_node_pair(result, reference_node, "0")
    if reference_index === nothing
        record_fixed_blocker!(
            result,
            :bpa_fixed_branch_blocked,
            :bpa_fixed_branch_missing_copy_reference,
        )
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "grounded scalar branch reference $(String(reference_node))-0 does not match a prior accepted scalar branch owner",
            ),
        )
        return true
    end
    branch_type = fixed_card_branch_type!(result, image, line_no, from_node, to_node)
    branch_type === nothing && return true
    reference_kind =
        result.elements[reference_index] isa CapacitorBranch ?
        :fixed_single_terminal_capacitance_reference :
        :fixed_grounded_scalar_branch_reference
    return parse_bpa_fixed_branch_copy_reference!(
        result,
        image,
        line_no,
        from_node,
        to_node,
        branch_type,
        reference_index,
        reference_kind,
        initial_issues;
        branch_layout_kind = reference_kind,
    )
end

function parse_bpa_fixed_single_terminal_capacitance_row!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    from_node::AbstractString,
    to_node::AbstractString,
    initial_issues::Int,
)::Bool
    capacitance = fixed_float_field!(result, image, line_no, 39, 44, "branch_capacitance")
    capacitance === nothing && return true
    branch_type = fixed_card_branch_type!(result, image, line_no, from_node, to_node)
    branch_type === nothing && return true
    return push_bpa_fixed_branch_scalar_row!(
        result,
        from_node,
        to_node,
        branch_type,
        0.0,
        0.0,
        Float64(capacitance),
        line_no,
        initial_issues;
        branch_layout_count = :bpa_fixed_branch_single_terminal_capacitance,
        output_image = image,
    )
end

function bpa_fixed_existing_node_index(result::DeckParseResult, node::AbstractString)
    normalized = lowercase(String(node))
    if normalized in ("0", "gnd", "ground", "ref")
        return 0
    end
    name = Symbol(String(node))
    if haskey(result.node_map, name)
        return result.node_map[name]
    end
    return nothing
end

function bpa_fixed_existing_node_index!(result::DeckParseResult,
                                        node::AbstractString,
                                        line_no::Int,
                                        field::AbstractString)
    index = bpa_fixed_existing_node_index(result, node)
    index === nothing || return index
    add_issue!(
        result.validation,
        invalid_value(
            "line $line_no",
            "$field=$(String(node)): expected an existing OVER2 branch reference node",
        ),
    )
    return nothing
end

function bpa_fixed_branch_reference_by_name!(result::DeckParseResult,
                                             reference_name::AbstractString,
                                             line_no::Int)
    name = Symbol(String(reference_name))
    element_index = findfirst(==(name), result.element_names)
    if element_index === nothing
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "OVER2 fixed-field branch COPY reference $(String(reference_name)) does not match a prior branch owner",
            ),
        )
        return nothing
    elseif !accepted_branch_output_element(result.elements[element_index])
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "OVER2 fixed-field branch COPY reference $(String(reference_name)) does not target an accepted scalar branch owner",
            ),
        )
        return nothing
    end
    return element_index
end

function fixed_branch_reference_by_existing_node_pair(
    result::DeckParseResult,
    reference_from::AbstractString,
    reference_to::AbstractString,
)
    from_index = bpa_fixed_existing_node_index(result, reference_from)
    to_index = bpa_fixed_existing_node_index(result, reference_to)
    if from_index === nothing || to_index === nothing
        return nothing
    end
    for (element_index, element) in enumerate(result.elements)
        accepted_branch_output_element(element) || continue
        if abs(element.a) == abs(from_index) && abs(element.b) == abs(to_index)
            return element_index
        end
    end
    return nothing
end

function bpa_fixed_branch_reference_by_node_pair!(result::DeckParseResult,
                                                  reference_from::AbstractString,
                                                  reference_to::AbstractString,
                                                  line_no::Int)
    from_index = bpa_fixed_existing_node_index!(result, reference_from, line_no, "reference_from_node")
    to_index = bpa_fixed_existing_node_index!(result, reference_to, line_no, "reference_to_node")
    if from_index === nothing || to_index === nothing
        return nothing
    end
    element_index =
        fixed_branch_reference_by_existing_node_pair(result, reference_from, reference_to)
    element_index === nothing || return element_index
    add_issue!(
        result.validation,
        invalid_value(
            "line $line_no",
            "OVER2 fixed-field branch reference node pair $(String(reference_from))-$(String(reference_to)) does not match a prior accepted scalar branch owner",
        ),
    )
    return nothing
end

function coupled_phase_pi_reference_by_node_pair(
    result::DeckParseResult,
    reference_from::AbstractString,
    reference_to::AbstractString,
)
    from_index = bpa_fixed_existing_node_index(result, reference_from)
    to_index = bpa_fixed_existing_node_index(result, reference_to)
    if from_index === nothing || to_index === nothing
        return nothing
    end
    for row in Iterators.reverse(result.coupled_phase_pi_section_rows)
        row.phase_index == 1 || continue
        if row.from_node_value == from_index && row.to_node_value == to_index
            return row
        end
    end
    return nothing
end

function push_bpa_fixed_coupled_phase_pi_row!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    from_node::AbstractString,
    to_node::AbstractString,
    branch_type::Int,
    initial_issues::Int;
    reference_kind::Symbol = :none,
    reference_from_node::Union{Nothing,String} = nothing,
    reference_to_node::Union{Nothing,String} = nothing,
)::Bool
    from_index = node_id!(result, from_node)
    to_index = node_id!(result, to_node)
    ref_from_name = reference_from_node === nothing ? missing : Symbol(reference_from_node)
    ref_to_name = reference_to_node === nothing ? missing : Symbol(reference_to_node)
    ref_from_index =
        reference_from_node === nothing ? missing :
        bpa_fixed_existing_node_index(result, reference_from_node)
    ref_to_index =
        reference_to_node === nothing ? missing :
        bpa_fixed_existing_node_index(result, reference_to_node)
    ref_from_index === nothing && (ref_from_index = missing)
    ref_to_index === nothing && (ref_to_index = missing)
    slices = fixed_card_branch_triplet_slices(image)
    output_code = bpa_fixed_branch_output_code_from_image(image)
    push!(
        result.coupled_phase_pi_section_rows,
        DeckCoupledPhasePiSectionRow(
            Symbol("coupled_phase_pi_row_", length(result.coupled_phase_pi_section_rows) + 1),
            Symbol(from_node),
            Symbol(to_node),
            from_index,
            to_index,
            line_no,
            branch_type,
            branch_type,
            reference_kind,
            ref_from_name,
            ref_to_name,
            ref_from_index,
            ref_to_index,
            slices.resistance,
            slices.inductance,
            slices.capacitance,
            Int[],
            output_code,
        ),
    )
    record_fixed_card!(result, :bpa_fixed_branch, :bpa_fixed_coupled_phase_pi_row, initial_issues)
    reference_kind == :none ||
        record_card!(result, :bpa_fixed_coupled_phase_pi_copy_reference)
    length(slices.resistance) < branch_type &&
        record_card!(result, :bpa_fixed_coupled_phase_pi_continuation_expected)
    return true
end

function coupled_phase_pi_copy_continuation_expected(
    result::DeckParseResult,
    branch_type::Int,
)::Bool
    branch_type > 1 || return false
    !isempty(result.coupled_phase_pi_section_rows) || return false
    result.coupled_phase_pi_section_rows[end].phase_index == branch_type - 1 ||
        return false
    for row in Iterators.reverse(result.coupled_phase_pi_section_rows)
        if row.phase_index == 1
            return row.reference_kind != :none
        end
    end
    return false
end

function bpa_fixed_nonlinear_reference_by_name!(
    result::DeckParseResult,
    rows::AbstractVector,
    family::AbstractString,
    reference_name::AbstractString,
    line_no::Int,
)
    name = Symbol(String(reference_name))
    row_index = findfirst(row -> row.name == name, rows)
    if row_index === nothing
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "OVER2 fixed-field nonlinear $family COPY reference $(String(reference_name)) does not match a prior nonlinear owner",
            ),
        )
        return nothing
    end
    return row_index
end

function bpa_fixed_nonlinear_reference_by_node_pair!(
    result::DeckParseResult,
    rows::AbstractVector,
    family::AbstractString,
    reference_from::AbstractString,
    reference_to::AbstractString,
    line_no::Int,
)
    from_index = bpa_fixed_existing_node_index!(
        result,
        reference_from,
        line_no,
        "nonlinear_reference_from_node",
    )
    to_index = bpa_fixed_existing_node_index!(
        result,
        reference_to,
        line_no,
        "nonlinear_reference_to_node",
    )
    if from_index === nothing || to_index === nothing
        return nothing
    end
    for (row_index, row) in enumerate(rows)
        if abs(row.from_node_index) == abs(from_index) &&
           abs(row.to_node_index) == abs(to_index)
            return row_index
        end
    end
    add_issue!(
        result.validation,
        invalid_value(
            "line $line_no",
            "OVER2 fixed-field nonlinear $family reference node pair $(String(reference_from))-$(String(reference_to)) does not match a prior nonlinear owner",
        ),
    )
    return nothing
end

function bpa_fixed_nonlinear_reference_index!(
    result::DeckParseResult,
    rows::AbstractVector,
    family::AbstractString,
    reference_kind::Symbol,
    reference_name::Union{Nothing,String},
    reference_from::Union{Nothing,String},
    reference_to::Union{Nothing,String},
    line_no::Int,
)
    if reference_kind == :bpa_fixed_branch_copy_name_reference
        reference_name === nothing && return nothing
        return bpa_fixed_nonlinear_reference_by_name!(
            result,
            rows,
            family,
            reference_name,
            line_no,
        )
    elseif reference_kind == :bpa_fixed_branch_copy_node_pair_reference
        (reference_from === nothing || reference_to === nothing) && return nothing
        return bpa_fixed_nonlinear_reference_by_node_pair!(
            result,
            rows,
            family,
            reference_from,
            reference_to,
            line_no,
        )
    end
    return nothing
end

function bpa_fixed_branch_copy_element(reference, from_index::Int, to_index::Int)
    if reference isa ConductanceBranch
        return ConductanceBranch(from_index, to_index, reference.g),
               :bpa_fixed_branch_copy_conductance
    elseif reference isa SeriesRLBranch
        return SeriesRLBranch(
                   from_index,
                   to_index,
                   reference.r,
                   reference.l,
                   reference.i_prev,
                   reference.v_prev,
                   reference.i_last,
               ),
               :bpa_fixed_branch_copy_series_rl
    elseif reference isa SeriesRLCBranch
        return SeriesRLCBranch(
                   from_index,
                   to_index,
                   reference.r,
                   reference.l,
                   reference.c,
                   reference.i_prev,
                   reference.inductor_voltage_prev,
                   reference.capacitor_voltage_prev,
                   reference.v_prev,
                   reference.i_last,
               ),
               :bpa_fixed_branch_copy_series_rlc
    elseif reference isa CapacitorBranch
        return CapacitorBranch(
                   from_index,
                   to_index,
                   reference.c,
                   reference.i_prev,
                   reference.v_prev,
                   reference.i_last,
               ),
               :bpa_fixed_branch_copy_capacitor
    end
    return nothing, :bpa_fixed_branch_copy_unsupported
end

function record_bpa_branch_output_requests_by_code!(result::DeckParseResult,
                                                    line_no::Int,
                                                    name::AbstractString,
                                                    output_code::Int)
    signal_output_code = output_code > 3 ? 3 : output_code
    if signal_output_code >= 2
        output_name = string("branch_voltage_", name)
        output_initial_issues = length(result.validation.issues)
        parse_over16_branch_voltage_output!(
            result,
            ["over16_branch_voltage_output", output_name, name],
            line_no,
        )
        if length(result.validation.issues) == output_initial_issues
            record_card!(result, :bpa_fixed_branch_voltage_output)
            record_card!(result, :output_branch_voltage)
        end
    end
    if signal_output_code == 1 || signal_output_code == 3
        output_name = string("branch_current_", name)
        output_initial_issues = length(result.validation.issues)
        parse_over16_branch_current_output!(
            result,
            ["over16_branch_current_output", output_name, name],
            line_no,
        )
        if length(result.validation.issues) == output_initial_issues
            record_card!(result, :bpa_fixed_branch_current_output)
            record_card!(result, :output_branch_current)
        end
    end
    if output_code > 3
        output_name = string("branch_power_", name)
        output_initial_issues = length(result.validation.issues)
        parse_over16_branch_power_output!(
            result,
            ["over16_branch_power_output", output_name, name],
            line_no,
        )
        if length(result.validation.issues) == output_initial_issues
            record_card!(result, :bpa_fixed_branch_power_output)
            record_card!(result, :output_branch_power_energy)
        end
    end
    return result
end

function record_bpa_fixed_branch_output_requests!(result::DeckParseResult,
                                                  image::AbstractString,
                                                  line_no::Int,
                                                  name::AbstractString)
    output_code = fixed_int_or_default!(result, image, line_no, 80, 80, "branch_output_code", 0)
    output_code === nothing && return result
    return record_bpa_branch_output_requests_by_code!(result, line_no, name, output_code)
end

function bpa_fixed_branch_output_code_from_image(image::AbstractString)::Int
    raw = fixed_field(image, 80, 80)
    isempty(raw) && return 0
    value = tryparse(Int, raw)
    return value === nothing ? 0 : value
end

function bpa_fixed_branch_layout_kind(branch_layout_count::Union{Nothing,Symbol})::Symbol
    branch_layout_count === :bpa_fixed_branch_compact && return :fixed_compact
    branch_layout_count === :bpa_fixed_branch_wide && return :fixed_wide
    branch_layout_count === :bpa_fixed_branch_free_field && return :free_field
    branch_layout_count === :bpa_fixed_branch_single_terminal_capacitance &&
        return :fixed_single_terminal_capacitance
    branch_layout_count === :fixed_card_branch_sparse_numeric_layout &&
        return :fixed_sparse_numeric
    branch_layout_count === :fixed_card_branch_implied_decimal_layout &&
        return :fixed_compact_implied_decimal
    return :unknown
end

function fixed_card_branch_type!(result::DeckParseResult, image::AbstractString,
                                 line_no::Int, from_node::AbstractString,
                                 to_node::AbstractString)
    raw = fixed_field(image, 1, 2)
    if isempty(raw) && !isempty(from_node) && !isempty(to_node)
        record_card!(result, :fixed_card_blank_branch_type_default)
        return 0
    end
    return fixed_int_field!(result, image, line_no, 1, 2, "branch_type")
end

function fixed_card_branch_terminal_pair!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
)
    from_node = fixed_field(image, 3, 8)
    to_node = fixed_field(image, 9, 14)
    if !isempty(from_node) && isempty(to_node)
        record_card!(result, :fixed_card_branch_reference_terminal_default)
        to_node = "0"
    end
    return from_node, to_node
end

function zinc_oxide_table_values(line::AbstractString)
    values = Float64[]
    for field in split(String(line))
        value = tryparse_deck_float(field)
        value === nothing && return nothing
        push!(values, Float64(value))
    end
    return values
end

function zinc_oxide_table_sentinel(line::AbstractString)::Bool
    values = zinc_oxide_table_values(line)
    values === nothing && return false
    return length(values) == 1 && values[1] == 9999.0
end

function zinc_oxide_three_value_row!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    row_kind::AbstractString,
)
    values = zinc_oxide_table_values(line)
    if values === nothing || length(values) != 3
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "expected zinc-oxide $row_kind row with three numeric fields",
            ),
        )
        return nothing
    end
    return values
end

function push_zinc_oxide_nonlinear_row!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    from_node::AbstractString,
    to_node::AbstractString,
    nonlinear_type::Int,
    initial_issues::Int,
    ;
    inline_name::Union{Nothing,AbstractString}=nothing,
    reference_kind::Symbol=:none,
    reference_index::Int=0,
)::Bool
    wide_resistance = fixed_float_value(image, 27, 42)
    wide_inductance = fixed_float_value(image, 43, 58)
    wide_capacitance_marker = fixed_float_value(image, 59, 74)
    if wide_resistance !== nothing && wide_inductance !== nothing &&
       wide_capacitance_marker !== nothing
        raw_resistance = wide_resistance
        raw_inductance = wide_inductance
        raw_capacitance_marker = wide_capacitance_marker
    else
        raw_resistance = fixed_float_or_default!(
            result,
            image,
            line_no,
            27,
            32,
            "zinc_oxide_reference_resistance",
            0.0,
        )
        raw_inductance = fixed_float_or_default!(
            result,
            image,
            line_no,
            33,
            38,
            "zinc_oxide_reference_inductance",
            0.0,
        )
        raw_capacitance_marker = fixed_float_or_default!(
            result,
            image,
            line_no,
            39,
            44,
            "zinc_oxide_capacitance_marker",
            0.0,
        )
    end
    output_code = fixed_int_or_default!(
        result,
        image,
        line_no,
        80,
        80,
        "zinc_oxide_output_code",
        0,
    )
    if raw_resistance === nothing || raw_inductance === nothing ||
       raw_capacitance_marker === nothing || output_code === nothing
        return true
    end
    if raw_resistance == 5555.0
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "zinc-oxide type-92 TR field must not be 5555"),
        )
        return true
    end
    if raw_capacitance_marker != 5555.0
        record_fixed_blocker!(result, :bpa_fixed_branch_blocked,
                              :bpa_fixed_branch_blocked_type)
        add_issue!(
            result.validation,
            unknown_field(
                "line $line_no",
                "Unsupported OVER2 branch type 92 with nonlinear table marker $(raw_capacitance_marker); only zinc-oxide marker 5555 is parser-owned",
            ),
        )
        return true
    end

    copied = reference_index > 0
    if copied
        1 <= reference_index <= length(result.zinc_oxide_nonlinear_rows) ||
            throw(ArgumentError("zinc-oxide nonlinear reference index is outside prior rows"))
    end
    name = bpa_fixed_nonlinear_owner_name(
        result,
        string("zinc_oxide_nonlinear_", length(result.zinc_oxide_nonlinear_rows) + 1);
        explicit_name = inline_name,
    )
    from_index = node_id!(result, from_node)
    to_index = node_id!(result, to_node)
    first_characteristic_index =
        copied ?
        result.zinc_oxide_nonlinear_rows[reference_index].first_characteristic_index :
        length(result.zinc_oxide_initialization_rows) +
        length(result.zinc_oxide_breakpoint_rows) + 1
    reference_name =
        copied ? result.zinc_oxide_nonlinear_rows[reference_index].name : :none
    reference_line_no =
        copied ? result.zinc_oxide_nonlinear_rows[reference_index].line_no : 0
    push!(
        result.zinc_oxide_nonlinear_rows,
        DeckZincOxideNonlinearRow(
            name,
            Symbol(String(from_node)),
            Symbol(String(to_node)),
            from_index,
            to_index,
            line_no,
            nonlinear_type,
            Float64(raw_resistance),
            Float64(raw_inductance),
            Float64(raw_capacitance_marker),
            Int(output_code),
            first_characteristic_index,
            copied ? :copy_reference : :characteristic_table,
            reference_kind,
            reference_index,
            reference_name,
            reference_line_no,
            String(image),
        ),
    )
    record_fixed_card!(result, :zinc_oxide_nonlinear, :zinc_oxide_nonlinear_row, initial_issues)
    if copied && length(result.validation.issues) == initial_issues
        record_card!(result, :zinc_oxide_copy_reference)
        record_card!(result, reference_kind)
    end
    inline_name === nothing || record_card!(result, :zinc_oxide_inline_name)
    record_card!(result, :zinc_oxide_type_92_row)
    return true
end

function parse_zinc_oxide_initialization_row!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    nonlinear_row_index::Int,
)::Bool
    initial_issues = length(result.validation.issues)
    values = zinc_oxide_three_value_row!(result, line, line_no, "initialization")
    values === nothing && return true
    push!(
        result.zinc_oxide_initialization_rows,
        DeckZincOxideInitializationRow(
            nonlinear_row_index,
            line_no,
            values[1],
            values[2],
            values[3],
            String(line),
        ),
    )
    record_fixed_card!(result, :zinc_oxide_nonlinear, :zinc_oxide_initialization, initial_issues)
    return true
end

function parse_zinc_oxide_breakpoint_row!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    nonlinear_row_index::Int,
    breakpoint_index::Int,
)::Bool
    initial_issues = length(result.validation.issues)
    values = zinc_oxide_three_value_row!(result, line, line_no, "breakpoint")
    values === nothing && return true
    push!(
        result.zinc_oxide_breakpoint_rows,
        DeckZincOxideBreakpointRow(
            nonlinear_row_index,
            breakpoint_index,
            line_no,
            values[1],
            values[2],
            values[3],
            String(line),
        ),
    )
    record_fixed_card!(result, :zinc_oxide_nonlinear, :zinc_oxide_breakpoint, initial_issues)
    return true
end

function nonlinear_resistance_table_values(line::AbstractString)
    values = Float64[]
    for field in split(String(line))
        value = tryparse_deck_float(field)
        value === nothing && return nothing
        push!(values, Float64(value))
    end
    return values
end

function nonlinear_resistance_table_sentinel(line::AbstractString)::Bool
    values = nonlinear_resistance_table_values(line)
    values === nothing && return false
    return length(values) == 1 && values[1] == 9999.0
end

function switching_nonlinear_resistor_table_sentinel(line::AbstractString)
    fields = split(strip(String(line)))
    isempty(fields) && return nothing
    marker = tryparse_deck_float(first(fields))
    marker == 9999.0 || return nothing
    words = uppercase.(String.(fields[2:end]))
    single_flash = "SINGLE" in words && "FLASH" in words
    isempty(words) || single_flash || return nothing
    return (single_flash = single_flash,)
end

triggered_timed_resistance_table_sentinel(line::AbstractString)::Bool =
    nonlinear_resistance_table_sentinel(line)

function push_triggered_timed_resistance_row!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    from_node::AbstractString,
    to_node::AbstractString,
    initial_issues::Int;
    inline_name::Union{Nothing,AbstractString}=nothing,
    reference_kind::Symbol=:none,
    reference_index::Int=0,
)::Bool
    trigger_voltage_v = fixed_float_or_default!(
        result, image, line_no, 27, 38,
        "triggered_timed_resistance_trigger_voltage_v", 0.0,
    )
    arm_time_s = fixed_float_or_default!(
        result, image, line_no, 39, 50,
        "triggered_timed_resistance_arm_time_s", 0.0,
    )
    output_code = fixed_int_or_default!(
        result, image, line_no, 80, 80,
        "triggered_timed_resistance_output_code", 0,
    )
    any(isnothing, (trigger_voltage_v, arm_time_s, output_code)) && return true
    Float64(trigger_voltage_v) >= 0.0 || begin
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "triggered timed-resistance voltage threshold must be nonnegative",
        ))
        return true
    end
    isfinite(Float64(arm_time_s)) || begin
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "triggered timed-resistance arm time must be finite",
        ))
        return true
    end
    copied = reference_index > 0
    reference_row = if copied
        1 <= reference_index <= length(result.triggered_timed_resistance_rows) ||
            throw(ArgumentError("timed-resistance reference index is outside prior rows"))
        result.triggered_timed_resistance_rows[reference_index]
    else
        nothing
    end
    name = bpa_fixed_nonlinear_owner_name(
        result,
        string(
            "triggered_timed_resistance_",
            length(result.triggered_timed_resistance_rows) + 1,
        );
        explicit_name = inline_name,
    )
    push!(
        result.triggered_timed_resistance_rows,
        DeckTriggeredTimedResistanceRow(
            name,
            Symbol(String(from_node)),
            Symbol(String(to_node)),
            node_id!(result, from_node),
            node_id!(result, to_node),
            line_no,
            Float64(trigger_voltage_v),
            reference_row === nothing ? Float64(arm_time_s) : reference_row.arm_time_s,
            Int(output_code),
            copied ? :copy_reference : :schedule_table,
            reference_kind,
            reference_index,
            reference_row === nothing ? :none : reference_row.name,
            reference_row === nothing ? 0 : reference_row.line_no,
            String(image),
        ),
    )
    if copied && length(result.validation.issues) == initial_issues
        record_card!(result, :triggered_timed_resistance_copy_reference)
        record_card!(result, reference_kind)
    end
    record_fixed_card!(
        result,
        :triggered_timed_resistance,
        :triggered_timed_resistance_row,
        initial_issues,
    )
    inline_name === nothing || record_card!(result, :triggered_timed_resistance_inline_name)
    return true
end

function parse_triggered_timed_resistance_point_row!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    resistance_row_index::Int,
    point_index::Int,
)::Bool
    initial_issues = length(result.validation.issues)
    values = nonlinear_resistance_table_values(line)
    if values === nothing || length(values) != 2
        image = fixed_image(line)
        elapsed = fixed_float_value(image, 1, 16)
        resistance = fixed_float_value(image, 17, 32)
        values = elapsed === nothing || resistance === nothing ?
            nothing : [Float64(elapsed), Float64(resistance)]
    end
    if values === nothing || length(values) != 2
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "expected timed-resistance schedule point with elapsed time and resistance",
        ))
        return true
    end
    elapsed_time_s = Float64(values[1])
    resistance_ohm = Float64(values[2])
    elapsed_time_s >= 0.0 || begin
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "timed-resistance elapsed time must be nonnegative",
        ))
        return true
    end
    resistance_ohm > 0.0 && isfinite(resistance_ohm) || begin
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "timed resistance must be finite and positive",
        ))
        return true
    end
    if point_index > 1
        prior = findlast(
            point -> point.resistance_row_index == resistance_row_index,
            result.triggered_timed_resistance_point_rows,
        )
        prior === nothing && throw(ArgumentError("timed-resistance point sequence is not contiguous"))
        elapsed_time_s > result.triggered_timed_resistance_point_rows[prior].elapsed_time_s || begin
            add_issue!(result.validation, invalid_value(
                "line $line_no",
                "timed-resistance elapsed times must strictly increase",
            ))
            return true
        end
    end
    push!(
        result.triggered_timed_resistance_point_rows,
        DeckTriggeredTimedResistancePointRow(
            resistance_row_index,
            point_index,
            line_no,
            elapsed_time_s,
            resistance_ohm,
            String(line),
        ),
    )
    record_fixed_card!(
        result,
        :triggered_timed_resistance,
        :triggered_timed_resistance_point,
        initial_issues,
    )
    return true
end

function push_switching_nonlinear_resistor_row!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    from_node::AbstractString,
    to_node::AbstractString,
    initial_issues::Int;
    inline_name::Union{Nothing,AbstractString}=nothing,
    reference_kind::Symbol=:none,
    reference_index::Int=0,
)::Bool
    turn_on_voltage = fixed_float_or_default!(
        result, image, line_no, 27, 38,
        "switching_nonlinear_resistor_turn_on_voltage", 0.0,
    )
    minimum_on_time_s = fixed_float_or_default!(
        result, image, line_no, 39, 50,
        "switching_nonlinear_resistor_minimum_on_time", 0.0,
    )
    raw_segment_count = fixed_float_or_default!(
        result, image, line_no, 51, 62,
        "switching_nonlinear_resistor_activation_segment_count", 1.0,
    )
    turn_off_voltage = fixed_float_or_default!(
        result, image, line_no, 63, 74,
        "switching_nonlinear_resistor_turn_off_voltage", 0.0,
    )
    output_code = fixed_int_or_default!(
        result, image, line_no, 80, 80,
        "switching_nonlinear_resistor_output_code", 0,
    )
    any(isnothing, (
        turn_on_voltage,
        minimum_on_time_s,
        raw_segment_count,
        turn_off_voltage,
        output_code,
    )) && return true
    Float64(turn_on_voltage) >= 0.0 || begin
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "switching nonlinear resistor turn-on voltage must be nonnegative",
        ))
        return true
    end
    Float64(minimum_on_time_s) >= 0.0 || begin
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "switching nonlinear resistor minimum on-time must be nonnegative",
        ))
        return true
    end
    Float64(turn_off_voltage) >= 0.0 || begin
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "switching nonlinear resistor turn-off voltage must be nonnegative",
        ))
        return true
    end
    segment_count = Float64(raw_segment_count) == 0.0 ? 1 : trunc(Int, raw_segment_count)
    segment_count >= 1 && Float64(raw_segment_count) == segment_count || begin
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "switching nonlinear resistor activation segment count must be a positive integer",
        ))
        return true
    end
    name = bpa_fixed_nonlinear_owner_name(
        result,
        string(
            "switching_nonlinear_resistor_",
            length(result.switching_nonlinear_resistor_rows) + 1,
        );
        explicit_name = inline_name,
    )
    copied = reference_index > 0
    if copied
        1 <= reference_index <= length(result.switching_nonlinear_resistor_rows) ||
            throw(ArgumentError("switching nonlinear resistor reference index is outside prior rows"))
    end
    reference_name =
        copied ? result.switching_nonlinear_resistor_rows[reference_index].name : :none
    reference_line_no =
        copied ? result.switching_nonlinear_resistor_rows[reference_index].line_no : 0
    push!(
        result.switching_nonlinear_resistor_rows,
        DeckSwitchingNonlinearResistorRow(
            name,
            Symbol(String(from_node)),
            Symbol(String(to_node)),
            node_id!(result, from_node),
            node_id!(result, to_node),
            line_no,
            Float64(turn_on_voltage),
            Float64(minimum_on_time_s),
            segment_count,
            Float64(turn_off_voltage) > Float64(turn_on_voltage) ?
                0.0 : Float64(turn_off_voltage),
            Int(output_code),
            false,
            copied ? :copy_reference : :vi_table,
            reference_kind,
            reference_index,
            reference_name,
            reference_line_no,
            String(image),
        ),
    )
    record_fixed_card!(
        result,
        :switching_nonlinear_resistor,
        :switching_nonlinear_resistor_row,
        initial_issues,
    )
    inline_name === nothing || record_card!(result, :switching_nonlinear_resistor_inline_name)
    if copied && length(result.validation.issues) == initial_issues
        record_card!(result, :switching_nonlinear_resistor_copy_reference)
        record_card!(result, reference_kind)
    end
    return true
end

function parse_switching_nonlinear_resistor_point_row!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    resistor_row_index::Int,
    point_index::Int,
)::Bool
    initial_issues = length(result.validation.issues)
    values = nonlinear_resistance_table_values(line)
    if values === nothing || length(values) != 2
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "expected switching nonlinear resistor point with current and voltage",
        ))
        return true
    end
    current_a, voltage_v = values
    previous_current = 0.0
    previous_voltage = 0.0
    if point_index > 1
        prior = findlast(
            row -> row.resistor_row_index == resistor_row_index,
            result.switching_nonlinear_resistor_point_rows,
        )
        prior === nothing && throw(ArgumentError("switching resistor point sequence is not contiguous"))
        previous = result.switching_nonlinear_resistor_point_rows[prior]
        previous_current = previous.current_a
        previous_voltage = previous.voltage_v
    end
    if !(current_a > previous_current && voltage_v > previous_voltage)
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "switching nonlinear resistor current and voltage points must both increase from the origin",
        ))
        return true
    end
    push!(
        result.switching_nonlinear_resistor_point_rows,
        DeckSwitchingNonlinearResistorPointRow(
            resistor_row_index,
            point_index,
            line_no,
            current_a,
            voltage_v,
            String(line),
        ),
    )
    record_fixed_card!(
        result,
        :switching_nonlinear_resistor,
        :switching_nonlinear_resistor_point,
        initial_issues,
    )
    return true
end

function finish_switching_nonlinear_resistor_table!(
    result::DeckParseResult,
    row_index::Int,
    single_flash::Bool,
)
    row = result.switching_nonlinear_resistor_rows[row_index]
    point_count = count(
        point -> point.resistor_row_index == row_index,
        result.switching_nonlinear_resistor_point_rows,
    )
    row.activation_segment_count <= point_count || begin
        add_issue!(result.validation, invalid_value(
            "line $(row.line_no)",
            "switching nonlinear resistor activation segment exceeds its V-I table",
        ))
        return result
    end
    result.switching_nonlinear_resistor_rows[row_index] =
        DeckSwitchingNonlinearResistorRow(
            row.name,
            row.from_node,
            row.to_node,
            row.from_node_index,
            row.to_node_index,
            row.line_no,
            row.turn_on_voltage,
            row.minimum_on_time_s,
            row.activation_segment_count,
            row.turn_off_voltage,
            row.output_code,
            single_flash,
            row.source_kind,
            row.reference_kind,
            row.reference_index,
            row.reference_name,
            row.reference_line_no,
            row.raw_text,
        )
    single_flash && record_card!(result, :switching_nonlinear_resistor_single_flash)
    return result
end

function nonlinear_resistance_initialization_row!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
)
    values = nonlinear_resistance_table_values(line)
    if values === nothing || length(values) != 3
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "expected nonlinear resistance initialization row with three numeric fields",
            ),
        )
        return nothing
    end
    return values
end

function nonlinear_resistance_point_row!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
)
    values = nonlinear_resistance_table_values(line)
    if values === nothing || length(values) != 2
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "expected nonlinear resistance characteristic row with two numeric fields",
            ),
        )
        return nothing
    end
    return values
end

function push_nonlinear_resistance_row!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    from_node::AbstractString,
    to_node::AbstractString,
    nonlinear_type::Int,
    initial_issues::Int,
    ;
    inline_name::Union{Nothing,AbstractString}=nothing,
    reference_kind::Symbol=:none,
    reference_index::Int=0,
)::Bool
    wide_reference = fixed_float_value(image, 27, 42)
    wide_secondary = fixed_float_value(image, 43, 58)
    wide_marker = fixed_float_value(image, 59, 74)
    if wide_reference !== nothing && wide_secondary !== nothing && wide_marker !== nothing
        steady_state_reference = wide_reference
        secondary_reference = wide_secondary
        table_marker = wide_marker
    else
        steady_state_reference = fixed_float_or_default!(
            result,
            image,
            line_no,
            27,
            32,
            "nonlinear_resistance_steady_state_reference",
            0.0,
        )
        secondary_reference = fixed_float_or_default!(
            result,
            image,
            line_no,
            33,
            38,
            "nonlinear_resistance_secondary_reference",
            0.0,
        )
        table_marker = fixed_float_or_default!(
            result,
            image,
            line_no,
            39,
            44,
            "nonlinear_resistance_table_marker",
            0.0,
        )
    end
    output_code = fixed_int_or_default!(
        result,
        image,
        line_no,
        80,
        80,
        "nonlinear_resistance_output_code",
        0,
    )
    if steady_state_reference === nothing || secondary_reference === nothing ||
       table_marker === nothing || output_code === nothing
        return true
    end
    expected_marker = nonlinear_type == 91 ? 3333.0 : 4444.0
    if Float64(table_marker) != expected_marker
        record_fixed_blocker!(result, :bpa_fixed_branch_blocked,
                              :bpa_fixed_branch_blocked_type)
        add_issue!(
            result.validation,
            unknown_field(
                "line $line_no",
                "Unsupported OVER2 branch type $nonlinear_type nonlinear resistance marker $(table_marker)",
            ),
        )
        return true
    end
    if nonlinear_type == 92 && Float64(steady_state_reference) == 5555.0
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "piecewise resistance type-92 TR field must not be 5555",
            ),
        )
        return true
    end

    element_kind =
        nonlinear_type == 91 ? :time_varying_resistance : :piecewise_resistance
    copied = reference_index > 0
    if copied
        1 <= reference_index <= length(result.nonlinear_resistance_rows) ||
            throw(ArgumentError("nonlinear resistance reference index is outside prior rows"))
        result.nonlinear_resistance_rows[reference_index].nonlinear_type == nonlinear_type ||
            throw(ArgumentError("nonlinear resistance COPY reference type must match the current row"))
    end
    name = bpa_fixed_nonlinear_owner_name(
        result,
        string(
            element_kind,
            "_nonlinear_",
            length(result.nonlinear_resistance_rows) + 1,
        );
        explicit_name = inline_name,
    )
    from_index = node_id!(result, from_node)
    to_index = node_id!(result, to_node)
    first_characteristic_index =
        copied ?
        result.nonlinear_resistance_rows[reference_index].first_characteristic_index :
        length(result.nonlinear_resistance_initialization_rows) +
        length(result.nonlinear_resistance_point_rows) + 1
    reference_name =
        copied ? result.nonlinear_resistance_rows[reference_index].name : :none
    reference_line_no =
        copied ? result.nonlinear_resistance_rows[reference_index].line_no : 0
    push!(
        result.nonlinear_resistance_rows,
        DeckNonlinearResistanceRow(
            name,
            Symbol(String(from_node)),
            Symbol(String(to_node)),
            from_index,
            to_index,
            line_no,
            nonlinear_type,
            element_kind,
            Float64(steady_state_reference),
            Float64(secondary_reference),
            Float64(table_marker),
            Int(output_code),
            first_characteristic_index,
            copied ? :copy_reference : :characteristic_table,
            reference_kind,
            reference_index,
            reference_name,
            reference_line_no,
            String(image),
        ),
    )
    record_fixed_card!(result, :nonlinear_resistance, :nonlinear_resistance_row, initial_issues)
    if copied && length(result.validation.issues) == initial_issues
        record_card!(result, :nonlinear_resistance_copy_reference)
        record_card!(result, reference_kind)
    end
    inline_name === nothing || record_card!(result, :nonlinear_resistance_inline_name)
    record_card!(result, nonlinear_type == 91 ?
                 :time_varying_resistance_type_91_row :
                 :piecewise_resistance_type_92_row)
    return true
end

function parse_nonlinear_resistance_initialization_row!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    nonlinear_row_index::Int,
)::Bool
    initial_issues = length(result.validation.issues)
    values = nonlinear_resistance_initialization_row!(result, line, line_no)
    values === nothing && return true
    push!(
        result.nonlinear_resistance_initialization_rows,
        DeckNonlinearResistanceInitializationRow(
            nonlinear_row_index,
            line_no,
            values[1],
            values[2],
            values[3],
            String(line),
        ),
    )
    record_fixed_card!(result, :nonlinear_resistance, :nonlinear_resistance_initialization, initial_issues)
    return true
end

function parse_nonlinear_resistance_point_row!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    nonlinear_row_index::Int,
    point_index::Int,
)::Bool
    initial_issues = length(result.validation.issues)
    values = nonlinear_resistance_point_row!(result, line, line_no)
    values === nothing && return true
    push!(
        result.nonlinear_resistance_point_rows,
        DeckNonlinearResistancePointRow(
            nonlinear_row_index,
            point_index,
            line_no,
            values[1],
            values[2],
            String(line),
        ),
    )
    record_fixed_card!(result, :nonlinear_resistance, :nonlinear_resistance_point, initial_issues)
    return true
end

function push_piecewise_nonlinear_inductor_row!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    from_node::AbstractString,
    to_node::AbstractString,
    nonlinear_type::Int,
    initial_issues::Int;
    inline_name::Union{Nothing,AbstractString}=nothing,
    reference_kind::Symbol=:none,
    reference_index::Int=0,
)::Bool
    nonlinear_type in (93, 98) ||
        throw(ArgumentError("nonlinear-inductor type must be 93 or 98"))
    steady_state_current = fixed_float_or_default!(
        result,
        image,
        line_no,
        27,
        38,
        "piecewise_nonlinear_inductor_steady_state_current_a",
        0.0,
    )
    steady_state_flux = fixed_float_or_default!(
        result,
        image,
        line_no,
        39,
        50,
        "piecewise_nonlinear_inductor_steady_state_flux_wb",
        0.0,
    )
    output_code = fixed_int_or_default!(
        result,
        image,
        line_no,
        80,
        80,
        "piecewise_nonlinear_inductor_output_code",
        0,
    )
    if steady_state_current === nothing || steady_state_flux === nothing ||
       output_code === nothing
        return true
    end
    isfinite(steady_state_current) || add_issue!(
        result.validation,
        invalid_value("line $line_no", "nonlinear-inductor steady-state current must be finite"),
    )
    isfinite(steady_state_flux) || add_issue!(
        result.validation,
        invalid_value(
            "line $line_no",
            "nonlinear-inductor steady-state flux must be finite",
        ),
    )
    fallback_name = nonlinear_type == 93 ?
        string(
            "piecewise_nonlinear_inductor_",
            length(result.piecewise_nonlinear_inductor_rows) + 1,
        ) :
        string(
            "pseudo_nonlinear_inductor_",
            length(result.piecewise_nonlinear_inductor_rows) + 1,
        )
    name = bpa_fixed_nonlinear_owner_name(
        result,
        fallback_name;
        explicit_name = inline_name,
    )
    source_kind = reference_index == 0 ? :direct : :copy_reference
    reference_row = reference_index == 0 ? nothing :
        result.piecewise_nonlinear_inductor_rows[reference_index]
    reference_row === nothing ||
        reference_row.nonlinear_type == nonlinear_type ||
        throw(ArgumentError("nonlinear-inductor COPY reference type must match the current row"))
    owned_steady_state_current = reference_row === nothing ?
        Float64(steady_state_current) : reference_row.steady_state_current_a
    owned_steady_state_flux = reference_row === nothing ?
        Float64(steady_state_flux) : reference_row.steady_state_flux_wb
    push!(
        result.piecewise_nonlinear_inductor_rows,
        DeckPiecewiseNonlinearInductorRow(
            name,
            Symbol(String(from_node)),
            Symbol(String(to_node)),
            node_id!(result, from_node),
            node_id!(result, to_node),
            line_no,
            nonlinear_type,
            owned_steady_state_current,
            owned_steady_state_flux,
            Int(output_code),
            source_kind,
            reference_kind,
            reference_index,
            reference_row === nothing ? Symbol("") : reference_row.name,
            reference_row === nothing ? 0 : reference_row.line_no,
            String(image),
        ),
    )
    record_fixed_card!(
        result,
        :piecewise_nonlinear_inductor,
        :piecewise_nonlinear_inductor_row,
        initial_issues,
    )
    record_card!(
        result,
        nonlinear_type == 93 ?
            :piecewise_nonlinear_inductor_type_93_row :
            :pseudo_nonlinear_inductor_type_98_row,
    )
    inline_name === nothing || record_card!(result, :piecewise_nonlinear_inductor_inline_name)
    return true
end

function parse_piecewise_nonlinear_inductor_point_row!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    row_index::Int,
    point_index::Int,
)::Bool
    initial_issues = length(result.validation.issues)
    values = nonlinear_resistance_table_values(line)
    if values === nothing || length(values) < 2
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "expected nonlinear-inductor current/flux point with at least two numeric fields",
            ),
        )
        return true
    end
    current = Float64(values[1])
    flux = Float64(values[2])
    if !isfinite(current) || !isfinite(flux)
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "nonlinear-inductor points must be finite"),
        )
        return true
    end
    if point_index > 1
        prior = last(filter(
            row -> row.inductor_row_index == row_index,
            result.piecewise_nonlinear_inductor_point_rows,
        ))
        if current <= prior.current_a || flux <= prior.flux_wb
            add_issue!(
                result.validation,
                invalid_value(
                    "line $line_no",
                    "nonlinear-inductor current and flux must both increase strictly",
                ),
            )
            return true
        end
    end
    push!(
        result.piecewise_nonlinear_inductor_point_rows,
        DeckPiecewiseNonlinearInductorPointRow(
            row_index,
            point_index,
            line_no,
            current,
            flux,
            String(line),
        ),
    )
    record_fixed_card!(
        result,
        :piecewise_nonlinear_inductor,
        :piecewise_nonlinear_inductor_point,
        initial_issues,
    )
    return true
end

function push_hysteretic_inductor_row!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    from_node::AbstractString,
    to_node::AbstractString,
    initial_issues::Int;
    inline_name::Union{Nothing,AbstractString}=nothing,
    reference_kind::Symbol=:none,
    reference_index::Int=0,
)::Bool
    steady_state_current = fixed_float_or_default!(
        result,
        image,
        line_no,
        27,
        38,
        "hysteretic_inductor_steady_state_current_A",
        0.0,
    )
    steady_state_flux = fixed_float_or_default!(
        result,
        image,
        line_no,
        39,
        50,
        "hysteretic_inductor_steady_state_flux_Wb",
        0.0,
    )
    residual_flux = fixed_float_or_default!(
        result,
        image,
        line_no,
        51,
        62,
        "hysteretic_inductor_residual_flux_Wb",
        0.0,
    )
    output_code = fixed_int_or_default!(
        result,
        image,
        line_no,
        80,
        80,
        "hysteretic_inductor_output_code",
        0,
    )
    if steady_state_flux === nothing || steady_state_current === nothing ||
       residual_flux === nothing || output_code === nothing
        return true
    end
    copied = reference_index > 0
    reference_row = if copied
        1 <= reference_index <= length(result.hysteretic_inductor_rows) ||
            throw(ArgumentError("hysteretic-inductor reference index is outside prior rows"))
        result.hysteretic_inductor_rows[reference_index]
    else
        nothing
    end
    name = bpa_fixed_nonlinear_owner_name(
        result,
        string("hysteretic_inductor_", length(result.hysteretic_inductor_rows) + 1);
        explicit_name = inline_name,
    )
    push!(
        result.hysteretic_inductor_rows,
        DeckHystereticInductorRow(
            name,
            Symbol(String(from_node)),
            Symbol(String(to_node)),
            node_id!(result, from_node),
            node_id!(result, to_node),
            line_no,
            Float64(steady_state_flux),
            Float64(steady_state_current),
            Float64(residual_flux),
            Int(output_code),
            reference_row === nothing ?
                length(result.hysteretic_inductor_point_rows) + 1 :
                reference_row.first_point_index,
            copied ? :copy_reference : :major_loop_table,
            reference_kind,
            reference_index,
            reference_row === nothing ? :none : reference_row.name,
            reference_row === nothing ? 0 : reference_row.line_no,
            String(image),
        ),
    )
    record_fixed_card!(result, :hysteretic_inductor, :hysteretic_inductor_row, initial_issues)
    record_card!(result, :hysteretic_inductor_type_96_row)
    if copied && length(result.validation.issues) == initial_issues
        record_card!(result, :hysteretic_inductor_copy_reference)
        record_card!(result, reference_kind)
    end
    inline_name === nothing || record_card!(result, :hysteretic_inductor_inline_name)
    return true
end

function parse_hysteretic_inductor_point_row!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    row_index::Int,
    point_index::Int,
)::Bool
    initial_issues = length(result.validation.issues)
    values = nonlinear_resistance_table_values(line)
    if values === nothing || length(values) != 2
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "expected type-96 current/flux point with two numeric fields",
            ),
        )
        return true
    end
    current = Float64(values[1])
    flux = Float64(values[2])
    if point_index == 1
        if current != 0.0 || flux != 0.0
            add_issue!(
                result.validation,
                invalid_value(
                    "line $line_no",
                    "type-96 major-loop data must start at zero current and zero flux",
                ),
            )
            return true
        end
    else
        prior = last(filter(
            row -> row.hysteretic_inductor_row_index == row_index,
            result.hysteretic_inductor_point_rows,
        ))
        if current <= prior.current_A || flux <= prior.flux_Wb
            add_issue!(
                result.validation,
                invalid_value(
                    "line $line_no",
                    "type-96 major-loop current and flux must both increase strictly",
                ),
            )
            return true
        end
    end
    push!(
        result.hysteretic_inductor_point_rows,
        DeckHystereticInductorPointRow(
            row_index,
            point_index,
            line_no,
            current,
            flux,
            String(line),
        ),
    )
    record_fixed_card!(result, :hysteretic_inductor, :hysteretic_inductor_point, initial_issues)
    return true
end

function arrester_constant_values(line::AbstractString)
    image = fixed_image(line)
    values = Float64[]
    fixed_field_count = 0
    for first_col in 1:16:65
        value = fixed_float_value(image, first_col, min(first_col + 15, 80))
        value === nothing && continue
        push!(values, value)
        fixed_field_count += 1
    end
    fixed_field_count > 0 && return values
    for field in split(String(line))
        value = tryparse_deck_float(field)
        value === nothing && return nothing
        push!(values, Float64(value))
    end
    return values
end

function push_arrester_nonlinear_row!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    from_node::AbstractString,
    to_node::AbstractString,
    nonlinear_type::Int,
    initial_issues::Int,
    ;
    inline_name::Union{Nothing,AbstractString}=nothing,
    reference_kind::Symbol=:none,
    reference_index::Int=0,
)::Bool
    wide_flashover_voltage = fixed_float_value(image, 27, 42)
    wide_voltage_division = fixed_float_value(image, 43, 58)
    wide_current_division = fixed_float_value(image, 59, 74)
    if wide_flashover_voltage !== nothing && wide_voltage_division !== nothing &&
       wide_current_division !== nothing
        flashover_voltage = wide_flashover_voltage
        voltage_division_factor = wide_voltage_division
        current_division_factor = wide_current_division
    else
        flashover_voltage = fixed_float_or_default!(
            result,
            image,
            line_no,
            27,
            32,
            "arrester_flashover_voltage",
            0.0,
        )
        voltage_division_factor = fixed_float_or_default!(
            result,
            image,
            line_no,
            33,
            38,
            "arrester_voltage_division_factor",
            1.0,
        )
        current_division_factor = fixed_float_or_default!(
            result,
            image,
            line_no,
            39,
            44,
            "arrester_current_division_factor",
            1.0,
        )
    end
    output_code = fixed_int_or_default!(
        result,
        image,
        line_no,
        80,
        80,
        "arrester_output_code",
        0,
    )
    if flashover_voltage === nothing || voltage_division_factor === nothing ||
       current_division_factor === nothing || output_code === nothing
        return true
    end
    voltage_division_factor == 0.0 && (voltage_division_factor = 1.0)
    current_division_factor == 0.0 && (current_division_factor = 1.0)

    copied = reference_index > 0
    if copied
        1 <= reference_index <= length(result.arrester_nonlinear_rows) ||
            throw(ArgumentError("arrester nonlinear reference index is outside prior rows"))
    end
    name = bpa_fixed_nonlinear_owner_name(
        result,
        string("arrester_nonlinear_", length(result.arrester_nonlinear_rows) + 1);
        explicit_name = inline_name,
    )
    from_index = node_id!(result, from_node)
    to_index = node_id!(result, to_node)
    first_constant_index =
        copied ?
        result.arrester_nonlinear_rows[reference_index].first_constant_index :
        sum((length(row.values) for row in result.arrester_constant_rows); init = 0) + 1
    reference_name =
        copied ? result.arrester_nonlinear_rows[reference_index].name : :none
    reference_line_no =
        copied ? result.arrester_nonlinear_rows[reference_index].line_no : 0
    push!(
        result.arrester_nonlinear_rows,
        DeckArresterNonlinearRow(
            name,
            Symbol(String(from_node)),
            Symbol(String(to_node)),
            from_index,
            to_index,
            line_no,
            nonlinear_type,
            Float64(flashover_voltage),
            Float64(voltage_division_factor),
            Float64(current_division_factor),
            Int(output_code),
            first_constant_index,
            copied ? :copy_reference : :constant_table,
            reference_kind,
            reference_index,
            reference_name,
            reference_line_no,
            String(image),
        ),
    )
    record_fixed_card!(result, :arrester_nonlinear, :arrester_nonlinear_row, initial_issues)
    if copied && length(result.validation.issues) == initial_issues
        record_card!(result, :arrester_copy_reference)
        record_card!(result, reference_kind)
    end
    inline_name === nothing || record_card!(result, :arrester_inline_name)
    record_card!(result, :arrester_type_94_row)
    return true
end

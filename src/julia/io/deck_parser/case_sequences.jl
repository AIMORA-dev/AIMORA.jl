"""
One independently parsed, executable data case from a physical deck stream.

`start_line_no` and `end_line_no` retain positions in the original stream, and
the contained parse result retains those same line numbers.
"""
struct DeckExecutableCase
    case_index::Int
    start_line_no::Int
    end_line_no::Int
    parsed::DeckParseResult
end

"""
Ordered data cases and lifecycle mutations accepted from one deck stream.

Aborted cases are deliberately absent from `cases`; their count and discarded
card count remain explicit so orchestration cannot accidentally execute them.
"""
struct DeckCaseSequence
    source::String
    boundaries::Vector{DeckCaseBoundaryRow}
    cases::Vector{DeckExecutableCase}
    aborted_case_count::Int
    discarded_card_count::Int
    run_terminated::Bool
end

function _case_sequence_raw_text(tokens)
    return join(token_strings(tokens), " ")
end

function _case_sequence_boundary_kind(tokens)
    begin_new_data_case_card(tokens) && return :new_data_case
    abort_data_case_card(tokens) && return :abort_data_case
    end_last_data_case_card(tokens) && return :end_data_case
    if !isempty(tokens) &&
       normalized_deck_token(tokens[1]) == "end" &&
       length(tokens) >= 4 &&
       normalized_deck_token(tokens[2]) == "new" &&
       normalized_deck_token(tokens[3]) == "data" &&
       normalized_deck_token(tokens[4]) == "case"
        return :end_data_case
    end
    marker = fixed_card_section_marker(tokens)
    marker !== nothing && marker.kind == :fixed_card_case_terminator &&
        return :run_termination
    return nothing
end

function _case_sequence_payload_line(tokens)
    isempty(tokens) && return false
    first_token = normalized_deck_token(tokens[1])
    first_token in ("c", "\$", "comment") && return false
    _case_sequence_boundary_kind(tokens) === nothing || return false
    return true
end

function _parse_sequence_case(
    lines::Vector{String},
    source::String,
    case_index::Int,
    start_line_no::Int,
    end_line_no::Int,
)
    padded = vcat(fill("", start_line_no - 1), lines[start_line_no:end_line_no])
    parsed = parse_deck_lines(
        padded;
        source="$(source)#case$(case_index):$(start_line_no)-$(end_line_no)",
    )
    return DeckExecutableCase(case_index, start_line_no, end_line_no, parsed)
end

"""
    parse_deck_case_sequence(lines; source="deck")

Split a physical deck stream at DATAIN lifecycle boundaries, discard aborted
cases, and parse every surviving case with isolated model/state ownership.
"""
function parse_deck_case_sequence(lines; source::AbstractString="deck")
    physical_lines = String.(collect(lines))
    source_text = String(source)
    boundaries = DeckCaseBoundaryRow[]
    cases = DeckExecutableCase[]
    case_start = 0
    case_has_payload = false
    aborted = false
    aborted_case_count = 0
    discarded_card_count = 0
    run_terminated = false

    function finish_case!(end_line_no::Int)
        if case_start > 0 && !aborted && case_has_payload && end_line_no >= case_start
            push!(
                cases,
                _parse_sequence_case(
                    physical_lines,
                    source_text,
                    length(cases) + 1,
                    case_start,
                    end_line_no,
                ),
            )
        end
        return nothing
    end

    for (line_no, line) in enumerate(physical_lines)
        tokens = deck_tokens(line)
        isempty(tokens) && continue
        boundary_kind = _case_sequence_boundary_kind(tokens)

        if aborted
            if boundary_kind == :new_data_case
                aborted = false
                case_start = line_no
                case_has_payload = false
                push!(
                    boundaries,
                    DeckCaseBoundaryRow(
                        line_no,
                        :new_data_case,
                        _case_sequence_raw_text(tokens),
                    ),
                )
            else
                discarded_card_count += 1
                push!(
                    boundaries,
                    DeckCaseBoundaryRow(
                        line_no,
                        :aborted_case_discarded_card,
                        _case_sequence_raw_text(tokens),
                    ),
                )
            end
            continue
        end

        if boundary_kind == :new_data_case
            finish_case!(line_no - 1)
            case_start = line_no
            case_has_payload = false
        elseif case_start == 0
            case_start = line_no
        end

        boundary_kind === nothing ||
            push!(
                boundaries,
                DeckCaseBoundaryRow(
                    line_no,
                    boundary_kind,
                    _case_sequence_raw_text(tokens),
                ),
            )

        if boundary_kind == :abort_data_case
            aborted_case_count += 1
            aborted = true
            case_has_payload = false
            continue
        end

        case_has_payload |= _case_sequence_payload_line(tokens)

        if boundary_kind == :end_data_case
            finish_case!(line_no)
            case_start = 0
            case_has_payload = false
            if end_last_data_case_card(tokens)
                run_terminated = true
                break
            end
        elseif boundary_kind == :run_termination
            finish_case!(line_no)
            case_start = 0
            case_has_payload = false
            run_terminated = true
            break
        end
    end

    !aborted && !run_terminated && finish_case!(length(physical_lines))
    return DeckCaseSequence(
        source_text,
        boundaries,
        cases,
        aborted_case_count,
        discarded_card_count,
        run_terminated,
    )
end

function parse_deck_file_sequence(path::AbstractString)
    return parse_deck_case_sequence(readlines(path); source=String(path))
end

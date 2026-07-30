
function switch_family_card(tokens)::Bool
    length(tokens) >= 2 || return false
    normalized_deck_token(tokens[2]) in ("ideal", "time", "timed", "time_switch")
end
function record_control_card!(result::DeckParseResult, kind::Symbol, tokens, line_no::Int)
    push!(result.control_cards, DeckControlCard(kind, deck_control_label(tokens), line_no,
                                                token_strings(tokens)))
    record_card!(result, kind)
    return record_card!(result, :section)
end

function record_case_boundary!(
    result::DeckParseResult,
    boundary_kind::Symbol,
    tokens,
    line_no::Int,
)
    push!(
        result.case_boundary_rows,
        DeckCaseBoundaryRow(line_no, boundary_kind, join(token_strings(tokens), " ")),
    )
    return result
end

function parse_control_card!(result::DeckParseResult, tokens, line_no::Int)
    card = lowercase(String(tokens[1]))
    tail = [lowercase(String(token)) for token in tokens[2:end]]
    kind = if begin_new_data_case_card(tokens)
        :begin_new_data_case
    elseif card == "end" && length(tail) >= 3 && tail[1:3] == ["new", "data", "case"]
        :end_new_data_case
    elseif card == "eldc" || card == "end" && length(tail) >= 3 && tail[1:3] == ["last", "data", "case"]
        :end_last_data_case
    elseif abort_data_case_card(tokens)
        :abort_data_case
    else
        Symbol(card)
    end
    if kind == :begin_new_data_case
        record_case_boundary!(result, :new_data_case, tokens, line_no)
    elseif kind in (:end_new_data_case, :end_last_data_case)
        record_case_boundary!(result, :end_data_case, tokens, line_no)
    elseif kind == :abort_data_case
        record_case_boundary!(result, :abort_data_case, tokens, line_no)
    end
    return record_control_card!(result, kind, tokens, line_no)
end

function parse_blank_control_card!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 2) || return result
    category = normalized_deck_token(tokens[2])
    kind = get(BPA_BLANK_SECTION_KINDS, category, nothing)
    if kind === nothing
        add_issue!(result.validation,
                   unknown_field("line $line_no",
                                 "Unsupported BLANK deck section $(String(tokens[2]))"))
        return result
    end
    return record_control_card!(result, kind, tokens, line_no)
end

function parse_branch_family_card!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 6) || return result
    kind = normalized_deck_token(tokens[2])
    name = require_family_value!(result, tokens, line_no, ("name", "id"), 3, "name")
    from_node = require_family_value!(result, tokens, line_no,
                                      ("from", "from_node", "bus1", "a"), 4, "from")
    to_node = require_family_value!(result, tokens, line_no,
                                    ("to", "to_node", "bus2", "b"), 5, "to")
    if name === nothing || from_node === nothing || to_node === nothing
        return result
    end
    initial_issues = length(result.validation.issues)
    if kind in ("resistor", "resistance", "r")
        resistance = require_family_value!(result, tokens, line_no,
                                           ("resistance", "r", "value"), 6, "resistance")
        resistance === nothing && return result
        parse_resistor!(result, ["resistor", name, from_node, to_node, resistance], line_no)
        return record_family_card!(result, :branch_family, :branch_resistor, initial_issues)
    elseif kind in ("conductance", "g")
        conductance = require_family_value!(result, tokens, line_no,
                                            ("conductance", "g", "value"), 6, "conductance")
        conductance === nothing && return result
        parse_conductance!(result, ["conductance", name, from_node, to_node, conductance], line_no)
        return record_family_card!(result, :branch_family, :branch_conductance, initial_issues)
    elseif kind in ("rl", "series_rl")
        resistance = require_family_value!(result, tokens, line_no,
                                           ("resistance", "r"), 6, "resistance")
        inductance = require_family_value!(result, tokens, line_no,
                                           ("inductance", "l"), 7, "inductance")
        if resistance === nothing || inductance === nothing
            return result
        end
        previous_current = deck_family_value(tokens, ("previous_current", "i_prev", "iprev"), 8)
        previous_voltage = deck_family_value(tokens, ("previous_voltage", "v_prev", "vprev"), 9)
        mapped = ["rl", name, from_node, to_node, resistance, inductance]
        previous_current === nothing || push!(mapped, previous_current)
        previous_voltage === nothing || push!(mapped, previous_voltage)
        parse_series_rl!(result, mapped, line_no)
        return record_family_card!(result, :branch_family, :branch_rl, initial_issues)
    elseif kind in ("inductor", "inductance", "l")
        inductance = require_family_value!(result, tokens, line_no,
                                           ("inductance", "l", "value"), 6, "inductance")
        inductance === nothing && return result
        previous_current = deck_family_value(tokens, ("previous_current", "i_prev", "iprev"), 7)
        previous_voltage = deck_family_value(tokens, ("previous_voltage", "v_prev", "vprev"), 8)
        mapped = ["inductor", name, from_node, to_node, inductance]
        previous_current === nothing || push!(mapped, previous_current)
        previous_voltage === nothing || push!(mapped, previous_voltage)
        parse_inductor!(result, mapped, line_no)
        return record_family_card!(result, :branch_family, :branch_inductor, initial_issues)
    elseif kind in ("capacitor", "capacitance", "c")
        capacitance = require_family_value!(result, tokens, line_no,
                                            ("capacitance", "c", "value"), 6, "capacitance")
        capacitance === nothing && return result
        previous_current = deck_family_value(tokens, ("previous_current", "i_prev", "iprev"), 7)
        previous_voltage = deck_family_value(tokens, ("previous_voltage", "v_prev", "vprev"), 8)
        mapped = ["capacitor", name, from_node, to_node, capacitance]
        previous_current === nothing || push!(mapped, previous_current)
        previous_voltage === nothing || push!(mapped, previous_voltage)
        parse_capacitor!(result, mapped, line_no)
        return record_family_card!(result, :branch_family, :branch_capacitor, initial_issues)
    end
    add_issue!(result.validation,
               unknown_field("line $line_no",
                             "Unsupported BRANCH family kind $(String(tokens[2]))"))
    return result
end

function parse_source_family_card!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 6) || return result
    kind = normalized_deck_token(tokens[2])
    if kind in ("current", "isource", "i")
        return parse_current_source_family_card!(result, tokens, line_no)
    end
    kind in ("voltage", "thevenin", "vsource", "v") || begin
        add_issue!(result.validation,
                   unknown_field("line $line_no",
                                 "Unsupported SOURCE family kind $(String(tokens[2]))"))
        return result
    end
    name = require_family_value!(result, tokens, line_no, ("name", "id"), 3, "name")
    node = require_family_value!(result, tokens, line_no, ("node", "bus", "bus1"), 4, "node")
    conductance = require_family_value!(result, tokens, line_no,
                                        ("conductance", "g"), 5, "conductance")
    amplitude = require_family_value!(result, tokens, line_no,
                                      ("amplitude", "amplitude_pu", "crest"), 6, "amplitude")
    frequency = require_family_value!(result, tokens, line_no,
                                      ("frequency", "frequency_hz", "sfreq"), 7, "frequency")
    if name === nothing || node === nothing || conductance === nothing ||
       amplitude === nothing || frequency === nothing
        return result
    end
    phase = deck_family_value(tokens, ("phase", "phase_rad"), 8)
    offset = deck_family_value(tokens, ("offset", "offset_pu"), 9)
    mapped = ["source", name, node, conductance, amplitude, frequency]
    phase === nothing || push!(mapped, phase)
    offset === nothing || push!(mapped, offset)
    initial_issues = length(result.validation.issues)
    parse_source!(result, mapped, line_no)
    return record_family_card!(result, :source_family, :source_voltage, initial_issues)
end

function parse_current_source_family_card!(result::DeckParseResult, tokens, line_no::Int)
    start_index = normalized_deck_token(tokens[1]) in ("current_source", "source_current",
                                                       "bpa_current_source") ? 1 : 2
    require_fields!(result, tokens, line_no, start_index + 4) || return result
    name = require_family_value!(result, tokens, line_no, ("name", "id"), start_index + 1, "name")
    node = require_family_value!(result, tokens, line_no, ("node", "bus", "bus1"), start_index + 2, "node")
    amplitude = require_family_value!(result, tokens, line_no,
                                      ("amplitude", "amplitude_pu", "crest"), start_index + 3, "amplitude")
    frequency = require_family_value!(result, tokens, line_no,
                                      ("frequency", "frequency_hz", "sfreq"), start_index + 4, "frequency")
    if name === nothing || node === nothing || amplitude === nothing || frequency === nothing
        return result
    end
    phase = deck_family_value(tokens, ("phase", "phase_rad"), start_index + 5)
    offset = deck_family_value(tokens, ("offset", "offset_pu"), start_index + 6)
    mapped = ["current", name, node, amplitude, frequency]
    phase === nothing || push!(mapped, phase)
    offset === nothing || push!(mapped, offset)
    initial_issues = length(result.validation.issues)
    parse_current!(result, mapped, line_no)
    return record_family_card!(result, :source_family, :source_current, initial_issues)
end

function parse_switch_family_card!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 6) || return result
    kind = normalized_deck_token(tokens[2])
    name = require_family_value!(result, tokens, line_no, ("name", "id"), 3, "name")
    from_node = require_family_value!(result, tokens, line_no,
                                      ("from", "from_node", "bus1", "a"), 4, "from")
    to_node = require_family_value!(result, tokens, line_no,
                                    ("to", "to_node", "bus2", "b"), 5, "to")
    if name === nothing || from_node === nothing || to_node === nothing
        return result
    end
    initial_issues = length(result.validation.issues)
    if kind == "ideal"
        state = require_family_value!(result, tokens, line_no,
                                      ("state", "closed", "initially_closed"), 6, "state")
        state === nothing && return result
        on_conductance = deck_family_value(tokens, ("on_conductance", "gon", "closed_conductance"), 7)
        off_conductance = deck_family_value(tokens, ("off_conductance", "goff", "open_conductance"), 8)
        mapped = ["switch", name, from_node, to_node, state]
        on_conductance === nothing || push!(mapped, on_conductance)
        off_conductance === nothing || push!(mapped, off_conductance)
        parse_switch!(result, mapped, line_no)
        return record_family_card!(result, :switch_family, :switch_ideal, initial_issues)
    elseif kind in ("time", "timed", "time_switch")
        close_time = require_family_value!(result, tokens, line_no,
                                           ("close_time", "close_time_s", "tclose"), 6, "close_time")
        open_time = require_family_value!(result, tokens, line_no,
                                          ("open_time", "open_time_s", "topen"), 7, "open_time")
        initial_state = require_family_value!(result, tokens, line_no,
                                              ("initially_closed", "state", "initial_state"), 8,
                                              "initially_closed")
        if close_time === nothing || open_time === nothing || initial_state === nothing
            return result
        end
        on_conductance = deck_family_value(tokens, ("on_conductance", "gon", "closed_conductance"), 9)
        off_conductance = deck_family_value(tokens, ("off_conductance", "goff", "open_conductance"), 10)
        mapped = ["time_switch", name, from_node, to_node, close_time, open_time, initial_state]
        on_conductance === nothing || push!(mapped, on_conductance)
        off_conductance === nothing || push!(mapped, off_conductance)
        parse_time_switch!(result, mapped, line_no)
        return record_family_card!(result, :switch_family, :switch_time, initial_issues)
    end
    add_issue!(result.validation,
               unknown_field("line $line_no",
                             "Unsupported SWITCH family kind $(String(tokens[2]))"))
    return result
end

function parse_output_family_card!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 4) || return result
    kind = normalized_deck_token(tokens[2])
    if kind in ("voltage", "over16", "over16_voltage", "node_voltage")
        name = require_family_value!(result, tokens, line_no, ("name", "id", "channel"), 3, "name")
        node = require_family_value!(result, tokens, line_no, ("node", "bus", "bus1"), 4, "node")
        if name === nothing || node === nothing
            return result
        end
        initial_issues = length(result.validation.issues)
        parse_over16_output!(result, ["over16_output", name, node], line_no)
        return record_family_card!(result, :output_family, :output_voltage, initial_issues)
    elseif kind in (
        "branch_voltage",
        "branch_volt",
        "over16_branch_voltage",
        "over16_branch_volt",
    )
        name = require_family_value!(result, tokens, line_no, ("name", "id", "channel"), 3, "name")
        branch = require_family_value!(result, tokens, line_no, ("branch", "element", "line"), 4, "branch")
        if name === nothing || branch === nothing
            return result
        end
        initial_issues = length(result.validation.issues)
        parse_over16_branch_voltage_output!(result, ["over16_branch_voltage_output", name, branch], line_no)
        return record_family_card!(result, :output_family, :output_branch_voltage, initial_issues)
    elseif kind in (
        "branch_current",
        "branch_curr",
        "over16_branch_current",
        "over16_branch_curr",
    )
        name = require_family_value!(result, tokens, line_no, ("name", "id", "channel"), 3, "name")
        branch = require_family_value!(result, tokens, line_no, ("branch", "element", "line"), 4, "branch")
        if name === nothing || branch === nothing
            return result
        end
        initial_issues = length(result.validation.issues)
        parse_over16_branch_current_output!(result, ["over16_branch_current_output", name, branch], line_no)
        return record_family_card!(result, :output_family, :output_branch_current, initial_issues)
    elseif kind in (
        "branch_power",
        "branch_energy",
        "branch_power_energy",
        "over16_branch_power",
        "over16_branch_energy",
        "over16_branch_power_energy",
    )
        name = require_family_value!(result, tokens, line_no, ("name", "id", "channel"), 3, "name")
        branch = require_family_value!(result, tokens, line_no, ("branch", "element", "line"), 4, "branch")
        if name === nothing || branch === nothing
            return result
        end
        initial_issues = length(result.validation.issues)
        parse_over16_branch_power_output!(result, ["over16_branch_power_output", name, branch], line_no)
        return record_family_card!(result, :output_family, :output_branch_power_energy, initial_issues)
    else
        add_issue!(result.validation,
                   unknown_field("line $line_no",
                                 "Unsupported OUTPUT family kind $(String(tokens[2]))"))
        return result
    end
end

function record_card!(result::DeckParseResult, card::Symbol)
    result.card_counts[card] = get(result.card_counts, card, 0) + 1
    return result
end

function require_fields!(result::DeckParseResult, tokens, line_no::Int, count::Int)::Bool
    if length(tokens) < count
        add_issue!(result.validation,
                   missing_data("line $line_no", "expected at least $count fields"))
        return false
    end
    return true
end

function parse_float!(result::DeckParseResult, token, line_no::Int, field::AbstractString)
    value = tryparse_deck_float(deck_token_value(token))
    if value === nothing
        add_issue!(result.validation,
                   invalid_value("line $line_no",
                                 "$field=$(String(token)): expected Float64"))
    end
    return value
end

function parse_int!(result::DeckParseResult, token, line_no::Int, field::AbstractString)
    value = tryparse(Int, deck_token_value(token))
    if value === nothing
        add_issue!(result.validation,
                   invalid_value("line $line_no",
                                 "$field=$(String(token)): expected Int"))
    end
    return value
end

function parse_positive_float!(result::DeckParseResult, token, line_no::Int, field::AbstractString)
    value = parse_float!(result, token, line_no, field)
    value === nothing && return nothing
    if !isfinite(value) || value <= 0.0
        add_issue!(result.validation,
                   invalid_value("line $line_no",
                                 "$field=$(String(token)): expected a positive finite Float64"))
        return nothing
    end
    return value
end

function parse_power_frequency_request!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 3) || return result
    frequency = parse_positive_float!(result, tokens[3], line_no, "frequency_hz")
    frequency === nothing && return result
    push!(
        result.power_frequency_request_rows,
        DeckPowerFrequencyRequestRow(line_no, frequency, join(token_strings(tokens), " ")),
    )
    record_card!(result, :power_frequency_request)
    return record_card!(result, :special_request)
end

function parse_universal_machine_dimension_request!(
    result::DeckParseResult,
    tokens,
    line_no::Int,
)
    first_value_index = compact_deck_keyword(tokens[1]) == "aumd" ? 2 : 4
    require_fields!(result, tokens, line_no, first_value_index + 3) || return result
    values = Vector{Int}(undef, 4)
    for index in 1:4
        value = parse_int!(
            result,
            tokens[first_value_index + index - 1],
            line_no,
            "universal_machine_dimension_$index",
        )
        value === nothing && return result
        if value <= 0
            add_issue!(
                result.validation,
                invalid_value(
                    "line $line_no",
                    "universal_machine_dimension_$index must be positive; got $value",
                ),
            )
            return result
        end
        values[index] = value
    end
    push!(
        result.universal_machine_dimension_request_rows,
        DeckUniversalMachineDimensionRequestRow(
            line_no,
            values[1],
            values[2],
            values[3],
            values[4],
            join(token_strings(tokens), " "),
        ),
    )
    record_card!(result, :universal_machine_dimension_request)
    return record_card!(result, :special_request)
end

function parse_output_width_request!(result::DeckParseResult, tokens, line_no::Int)
    column_width = output_width_request_columns(tokens)
    if !(column_width in (80, 132))
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "output width must be 80 or 132"),
        )
        return result
    end
    push!(
        result.output_width_request_rows,
        DeckOutputWidthRequestRow(line_no, column_width, join(token_strings(tokens), " ")),
    )
    record_card!(result, :output_width_request)
    record_card!(result, Symbol("output_width_$(column_width)_request"))
    return record_card!(result, :special_request)
end

function request_numeric_values!(result::DeckParseResult, tokens, line_no::Int,
                                 field::AbstractString, required_count::Int)
    values = Float64[]
    for token in tokens
        value = tryparse_deck_float(deck_token_value(token))
        value === nothing && continue
        push!(values, value)
    end
    if length(values) < required_count
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "expected at least $required_count numeric values for $field",
            ),
        )
        return nothing
    end
    return values
end

function parse_peak_voltage_monitor_request!(result::DeckParseResult, tokens, line_no::Int)
    push!(
        result.peak_voltage_monitor_request_rows,
        DeckPeakVoltageMonitorRequestRow(line_no, join(token_strings(tokens), " ")),
    )
    record_card!(result, :peak_voltage_monitor_request)
    return record_card!(result, :special_request)
end

function parse_diagnostic_print_request!(result::DeckParseResult, tokens, line_no::Int)
    values = request_numeric_values!(result, tokens, line_no, "diagnostic_print_control", 4)
    values === nothing && return result
    controls = ntuple(index -> trunc(Int, values[index]), 4)
    push!(
        result.diagnostic_print_request_rows,
        DeckDiagnosticPrintRequestRow(line_no, controls, join(token_strings(tokens), " ")),
    )
    record_card!(result, :diagnostic_print_request)
    return record_card!(result, :special_request)
end

function parse_tacs_warning_limit_request!(result::DeckParseResult, tokens, line_no::Int)
    values = request_numeric_values!(result, tokens, line_no, "tacs_warning_limit", 2)
    values === nothing && return result
    push!(
        result.tacs_warning_limit_request_rows,
        DeckTACSWarningLimitRequestRow(
            line_no,
            trunc(Int, values[1]),
            values[2],
            join(token_strings(tokens), " "),
        ),
    )
    record_card!(result, :tacs_warning_limit_request)
    return record_card!(result, :special_request)
end

function next_plot_file_request_mode(result::DeckParseResult)
    current = isempty(result.plot_file_request_rows) ? 2 :
              last(result.plot_file_request_rows).plot_file_mode
    current == 0 && return 2
    current == 2 && return 0
    return current
end

function parse_plot_file_request!(result::DeckParseResult, tokens, line_no::Int)
    plot_file_mode = next_plot_file_request_mode(result)
    push!(
        result.plot_file_request_rows,
        DeckPlotFileRequestRow(line_no, plot_file_mode, join(token_strings(tokens), " ")),
    )
    record_card!(result, :plot_file_request)
    record_card!(result, Symbol("plot_file_mode_$(plot_file_mode)_request"))
    return record_card!(result, :special_request)
end

function next_switch_logic_request_value(result::DeckParseResult)
    current = isempty(result.switch_logic_request_rows) ? 0 :
              last(result.switch_logic_request_rows).control_value
    value = current + 1
    return value >= 2 ? 0 : value
end

function parse_switch_logic_request!(result::DeckParseResult, tokens, line_no::Int)
    control_value = next_switch_logic_request_value(result)
    push!(
        result.switch_logic_request_rows,
        DeckSwitchLogicRequestRow(line_no, control_value, join(token_strings(tokens), " ")),
    )
    record_card!(result, :switch_logic_request)
    record_card!(result, Symbol("switch_logic_value_$(control_value)_request"))
    return record_card!(result, :special_request)
end

function simulation_control_request_kind(tokens)
    first_token = compact_deck_keyword(tokens[1])
    first_token == "obc" && return :omit_base_case
    first_token == "cs" && return :switch_pseudononlinear_conversion
    first_token == "mdc" && return :miscellaneous_data_cards
    first_token == "rte" && return :tolerance_constant
    first_token == "bpvs" && return :peak_value_search_start
    first_token == "todr" && return :statistics_table_save_time
    first_token == "usst" && return :random_switch_time_file
    deck_phrase_match(tokens, ("omit", "base", "case")) && return :omit_base_case
    deck_phrase_match(tokens, ("change", "switch")) &&
        return :switch_pseudononlinear_conversion
    deck_phrase_match(tokens, ("miscellaneous", "data", "cards")) &&
        return :miscellaneous_data_cards
    deck_phrase_match(tokens, ("redefine", "tolerance", "epsiln")) &&
        return :tolerance_constant
    deck_phrase_match(tokens, ("time", "of", "dice", "roll")) &&
        return :statistics_table_save_time
    deck_phrase_match(tokens, ("user", "supplied", "switch", "times")) &&
        return :random_switch_time_file
    return nothing
end

function parse_simulation_control_request!(result::DeckParseResult, tokens, line_no::Int)
    request_kind = simulation_control_request_kind(tokens)
    if request_kind === nothing
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "unsupported simulation control request"),
        )
        return result
    end
    numeric_value = missing
    integer_value = missing
    if request_kind in (
        :tolerance_constant,
        :peak_value_search_start,
        :statistics_table_save_time,
    )
        values = request_numeric_values!(result, tokens, line_no, String(request_kind), 1)
        values === nothing && return result
        numeric_value = values[1]
    elseif request_kind == :random_switch_time_file
        values = request_numeric_values!(result, tokens, line_no, "random_switch_time_file", 1)
        values === nothing && return result
        file_unit = trunc(Int, values[1])
        integer_value = file_unit == 0 ? 24 : file_unit
    end
    push!(
        result.simulation_control_request_rows,
        DeckSimulationControlRequestRow(
            line_no,
            request_kind,
            numeric_value,
            integer_value,
            join(token_strings(tokens), " "),
        ),
    )
    record_card!(result, :simulation_control_request)
    record_card!(result, Symbol(String(request_kind), "_request"))
    return record_card!(result, :special_request)
end

function request_integer_values(values::AbstractVector{<:Real}, count::Int)
    return Int[trunc(Int, values[index]) for index in 1:count]
end

function diagnostic_code_values!(result::DeckParseResult, tokens, line_no::Int)
    payload = join(String.(tokens[2:end]), "")
    if isempty(payload)
        add_issue!(
            result.validation,
            missing_data("line $line_no", "expected diagnostic code payload"),
        )
        return nothing
    end
    codes = Int[]
    index = firstindex(payload)
    while index <= lastindex(payload)
        next = nextind(payload, index)
        if next > lastindex(payload)
            add_issue!(
                result.validation,
                invalid_value("line $line_no", "diagnostic code payload has odd width"),
            )
            return nothing
        end
        value = tryparse(Int, payload[index:next])
        if value === nothing
            add_issue!(
                result.validation,
                invalid_value("line $line_no", "diagnostic code payload must be I2 values"),
            )
            return nothing
        end
        push!(codes, value)
        index = nextind(payload, next)
    end
    return codes
end

function automatic_branch_naming_value(result::DeckParseResult)
    prior_count = count(
        row -> row.request_kind == :automatic_branch_naming,
        result.study_option_request_rows,
    )
    return isodd(prior_count) ? 0 : 1
end

function parse_study_option_request!(result::DeckParseResult, tokens, line_no::Int)
    request_kind = study_option_request_kind(tokens)
    if request_kind === nothing
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "unsupported study option request"),
        )
        return result
    end
    numeric_values = Float64[]
    integer_values = Int[]
    text_values = String[]
    if request_kind in (
        :type99_message_limit,
        :printer_lines_per_inch,
        :mode_voltage_output,
        :plot_oscillation_limit,
        :statistics_output_salvage,
    )
        values = request_numeric_values!(result, tokens, line_no, String(request_kind), 1)
        values === nothing && return result
        integer_values = request_integer_values(values, 1)
    elseif request_kind == :plotter_paper_height
        values = request_numeric_values!(result, tokens, line_no, "plotter_paper_height", 1)
        values === nothing && return result
        numeric_values = [values[1]]
    elseif request_kind in (:frequency_scan, :line_model_frequency_scan)
        values = request_numeric_values!(result, tokens, line_no, String(request_kind), 4)
        values === nothing && return result
        numeric_values = values[1:3]
        integer_values = request_integer_values(values[4:4], 1)
    elseif request_kind == :diagnostic_codes
        codes = diagnostic_code_values!(result, tokens, line_no)
        codes === nothing && return result
        integer_values = codes
    elseif request_kind == :user_identification
        text_values = token_strings(tokens[2:end])
        if isempty(text_values)
            add_issue!(
                result.validation,
                missing_data("line $line_no", "expected user identification text"),
            )
            return result
        end
    elseif request_kind == :free_format_characters
        text_values = token_strings(tokens[2:end])
        if length(text_values) < 2
            add_issue!(
                result.validation,
                missing_data("line $line_no", "expected free-format separator and continuation characters"),
            )
            return result
        end
    elseif request_kind == :high_resistance_exponent
        values = request_numeric_values!(result, tokens, line_no, "high_resistance_exponent", 1)
        values === nothing && return result
        exponent = trunc(Int, values[1])
        integer_values = [exponent]
        numeric_values = [10.0 ^ exponent]
    elseif request_kind == :automatic_branch_naming
        integer_values = [automatic_branch_naming_value(result)]
    elseif request_kind == :zinc_oxide_constants
        values = request_numeric_values!(result, tokens, line_no, "zinc_oxide_constants", 6)
        values === nothing && return result
        integer_values = request_integer_values(values, 1)
        numeric_values = values[2:6]
    end
    push!(
        result.study_option_request_rows,
        DeckStudyOptionRequestRow(
            line_no,
            request_kind,
            numeric_values,
            integer_values,
            text_values,
            join(token_strings(tokens), " "),
        ),
    )
    record_card!(result, :study_option_request)
    record_card!(result, Symbol(String(request_kind), "_request"))
    return record_card!(result, :special_request)
end

function printout_frequency_active_pair_count(change_steps::NTuple{5,Int})
    for index in 1:5
        change_steps[index] == 0 && return index - 1
    end
    return 5
end

function parse_printout_frequency_change_payload!(
    result::DeckParseResult,
    request_line_no::Int,
    request_text::AbstractString,
    tokens,
    payload_line_no::Int,
)
    values = request_numeric_values!(
        result,
        tokens,
        payload_line_no,
        "printout_frequency_change",
        10,
    )
    values === nothing && return result
    integers = ntuple(index -> trunc(Int, values[index]), 10)
    change_steps = ntuple(index -> integers[2 * index - 1], 5)
    multipliers = ntuple(index -> integers[2 * index], 5)
    push!(
        result.printout_frequency_change_rows,
        DeckPrintoutFrequencyChangeRow(
            request_line_no,
            payload_line_no,
            change_steps,
            multipliers,
            printout_frequency_active_pair_count(change_steps),
            length(values),
            String(request_text),
            join(token_strings(tokens), " "),
        ),
    )
    record_card!(result, :printout_frequency_change_payload)
    return record_card!(result, :special_request_payload)
end

function parse_tacs_dimension_payload_values!(
    result::DeckParseResult,
    tokens,
    line_no::Int,
)
    values = Float64[]
    for (index, token) in enumerate(tokens)
        value = parse_float!(result, token, line_no, "tacs_dimension_value_$index")
        value === nothing && continue
        push!(values, value)
    end
    if length(values) < TACS_DIMENSION_VALUE_COUNT
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "expected at least $TACS_DIMENSION_VALUE_COUNT TACS dimension values",
            ),
        )
        return nothing
    end
    return values
end

function absolute_tacs_dimension_list_sizes(values::NTuple{8,Float64})
    sizes = ntuple(index -> trunc(Int, values[index]), TACS_DIMENSION_VALUE_COUNT)
    return (max(sizes[1], 5), sizes[2], sizes[3], sizes[4],
            sizes[5], sizes[6], sizes[7], sizes[8])
end

function relative_tacs_dimension_list_sizes(
    values::NTuple{8,Float64};
    total_storage_cells::Int = DEFAULT_TACS_TOTAL_STORAGE_CELLS,
    numeric_byte_counts::NTuple{4,Int} = DEFAULT_TACS_NUMERIC_BYTE_COUNTS,
)
    total_weight = sum(values)
    if total_weight <= 0.0
        return nothing
    end
    _byte1, _byte2, byte3, byte4 = numeric_byte_counts
    scale = total_storage_cells * byte3 / total_weight
    denominators = (
        4 * byte3 + 8 * byte4,
        2 * byte3 + byte4,
        2 * byte4,
        5 * byte3 + byte4,
        3 * byte4,
        byte4,
        byte3,
        6 * byte3 + 2 * byte4,
    )
    return ntuple(
        index -> trunc(Int, values[index] * scale / denominators[index]),
        TACS_DIMENSION_VALUE_COUNT,
    )
end

function parse_tacs_dimension_request_payload!(
    result::DeckParseResult,
    allocation_kind::Symbol,
    request_line_no::Int,
    request_text::AbstractString,
    tokens,
    payload_line_no::Int,
)
    values = parse_tacs_dimension_payload_values!(result, tokens, payload_line_no)
    values === nothing && return result
    request_values = ntuple(index -> values[index], TACS_DIMENSION_VALUE_COUNT)
    list_sizes = if allocation_kind == :absolute_tacs_dimension_request
        absolute_tacs_dimension_list_sizes(request_values)
    elseif allocation_kind == :relative_tacs_dimension_request
        sizes = relative_tacs_dimension_list_sizes(request_values)
        if sizes === nothing
            add_issue!(
                result.validation,
                invalid_value(
                    "line $payload_line_no",
                    "relative TACS dimension values must have positive total weight",
                ),
            )
            return result
        end
        sizes
    else
        add_issue!(
            result.validation,
            invalid_value("line $request_line_no", "unsupported TACS dimension request kind"),
        )
        return result
    end
    push!(
        result.tacs_dimension_request_rows,
        DeckTACSDimensionRequestRow(
            request_line_no,
            payload_line_no,
            allocation_kind,
            request_values,
            list_sizes,
            length(values),
            String(request_text),
            join(token_strings(tokens), " "),
        ),
    )
    record_card!(result, Symbol(String(allocation_kind), "_payload"))
    return record_card!(result, :tacs_dimension_request_payload)
end

function parse_bool!(result::DeckParseResult, token, line_no::Int, field::AbstractString)
    value = lowercase(deck_token_value(token))
    if value in ("true", "closed", "on", "1")
        return true
    elseif value in ("false", "open", "off", "0")
        return false
    end
    add_issue!(result.validation,
               invalid_value("line $line_no",
                             "$field=$(String(token)): expected true/false, closed/open, on/off, or 1/0"))
    return nothing
end

function optional_float!(result::DeckParseResult, tokens, line_no::Int, index::Int,
                         field::AbstractString, default::Float64)
    index > length(tokens) && return default
    value = parse_float!(result, tokens[index], line_no, field)
    return value === nothing ? default : value
end

function node_id!(result::DeckParseResult, token)::Int
    value = lowercase(String(token))
    if value in ("0", "gnd", "ground", "ref")
        return 0
    end
    name = Symbol(String(token))
    return get!(result.node_map, name, length(result.node_map) + 1)
end

function push_element!(result::DeckParseResult, tokens, element, line_no::Int)
    push!(result.elements, element)
    push!(result.element_line_numbers, line_no)
    push!(result.element_names, Symbol(String(tokens[2])))
    return result
end

function parse_bus!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 2) || return result
    node_id!(result, tokens[2])
    return record_card!(result, :bus)
end

function parse_source!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 6) || return result
    node = node_id!(result, tokens[3])
    conductance = parse_float!(result, tokens[4], line_no, "conductance")
    amplitude = parse_float!(result, tokens[5], line_no, "amplitude_pu")
    frequency = parse_float!(result, tokens[6], line_no, "frequency_hz")
    phase = optional_float!(result, tokens, line_no, 7, "phase_rad", 0.0)
    offset = optional_float!(result, tokens, line_no, 8, "offset_pu", 0.0)
    if conductance === nothing || amplitude === nothing || frequency === nothing
        return result
    end
    push_element!(result, tokens,
                  sinusoidal_thevenin_source(node, conductance, amplitude,
                                             frequency; phase_rad=phase,
                                             offset_pu=offset), line_no)
    return record_card!(result, :source)
end

function parse_current!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 5) || return result
    node = node_id!(result, tokens[3])
    amplitude = parse_float!(result, tokens[4], line_no, "amplitude_pu")
    frequency = parse_float!(result, tokens[5], line_no, "frequency_hz")
    phase = optional_float!(result, tokens, line_no, 6, "phase_rad", 0.0)
    offset = optional_float!(result, tokens, line_no, 7, "offset_pu", 0.0)
    if amplitude === nothing || frequency === nothing
        return result
    end
    amplitude_value = Float64(amplitude)
    frequency_value = Float64(frequency)
    phase_value = Float64(phase)
    offset_value = Float64(offset)
    push_element!(
        result,
        tokens,
        CurrentInjection(
            node,
            t -> sinusoidal_value(
                t,
                amplitude_value,
                frequency_value;
                phase_rad=phase_value,
                offset_pu=offset_value,
            ),
        ),
        line_no,
    )
    return record_card!(result, :current)
end

function parse_conductance!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 5) || return result
    a = node_id!(result, tokens[3])
    b = node_id!(result, tokens[4])
    conductance = parse_float!(result, tokens[5], line_no, "conductance")
    conductance === nothing && return result
    push_element!(result, tokens, ConductanceBranch(a, b, conductance), line_no)
    return record_card!(result, :conductance)
end

function parse_resistor!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 5) || return result
    a = node_id!(result, tokens[3])
    b = node_id!(result, tokens[4])
    resistance = parse_positive_float!(result, tokens[5], line_no, "resistance")
    resistance === nothing && return result
    push_element!(result, tokens, ConductanceBranch(a, b, inv(resistance)), line_no)
    return record_card!(result, :resistor)
end

function parse_series_rl!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 6) || return result
    a = node_id!(result, tokens[3])
    b = node_id!(result, tokens[4])
    resistance = parse_float!(result, tokens[5], line_no, "resistance")
    inductance = parse_float!(result, tokens[6], line_no, "inductance")
    previous_current = optional_float!(result, tokens, line_no, 7, "previous_current", 0.0)
    previous_voltage = optional_float!(result, tokens, line_no, 8, "previous_voltage", 0.0)
    if resistance === nothing || inductance === nothing
        return result
    end
    if !isfinite(resistance)
        add_issue!(result.validation,
                   invalid_value("line $line_no",
                                 "resistance=$(String(tokens[5])): expected a finite Float64"))
        return result
    end
    if !isfinite(inductance) || inductance <= 0.0
        add_issue!(result.validation,
                   invalid_value("line $line_no",
                                 "inductance=$(String(tokens[6])): expected a positive finite Float64"))
        return result
    end
    if !isfinite(previous_current) || !isfinite(previous_voltage)
        add_issue!(result.validation,
                   invalid_value("line $line_no",
                                 "series R-L initial current and voltage must be finite"))
        return result
    end
    push_element!(
        result,
        tokens,
        SeriesRLBranch(
            a,
            b,
            resistance,
            inductance,
            previous_current,
            previous_voltage,
            previous_current,
        ),
        line_no,
    )
    return record_card!(result, :rl)
end

function parse_inductor!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 5) || return result
    a = node_id!(result, tokens[3])
    b = node_id!(result, tokens[4])
    inductance = parse_positive_float!(result, tokens[5], line_no, "inductance")
    previous_current = optional_float!(result, tokens, line_no, 6, "previous_current", 0.0)
    previous_voltage = optional_float!(result, tokens, line_no, 7, "previous_voltage", 0.0)
    inductance === nothing && return result
    push_element!(
        result,
        tokens,
        SeriesRLBranch(a, b, 0.0, inductance, previous_current, previous_voltage, previous_current),
        line_no,
    )
    return record_card!(result, :inductor)
end

function parse_capacitor!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 5) || return result
    a = node_id!(result, tokens[3])
    b = node_id!(result, tokens[4])
    capacitance = parse_float!(result, tokens[5], line_no, "capacitance")
    previous_current = optional_float!(result, tokens, line_no, 6, "previous_current", 0.0)
    previous_voltage = optional_float!(result, tokens, line_no, 7, "previous_voltage", 0.0)
    capacitance === nothing && return result
    push_element!(
        result,
        tokens,
        CapacitorBranch(a, b, capacitance, previous_current, previous_voltage, previous_current),
        line_no,
    )
    return record_card!(result, :capacitor)
end

function parse_breqiv!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 11) || return result
    a = node_id!(result, tokens[3])
    b = node_id!(result, tokens[4])
    resistance = parse_float!(result, tokens[5], line_no, "resistance")
    inductance = parse_float!(result, tokens[6], line_no, "inductance")
    capacitance = parse_float!(result, tokens[7], line_no, "capacitance")
    damping_resistance = parse_float!(result, tokens[8], line_no, "damping_resistance")
    previous_current = parse_float!(result, tokens[9], line_no, "previous_current")
    capacitor_voltage = parse_float!(result, tokens[10], line_no, "capacitor_voltage")
    initial_voltage = parse_float!(result, tokens[11], line_no, "initial_voltage")
    if resistance === nothing || inductance === nothing || capacitance === nothing ||
       damping_resistance === nothing || previous_current === nothing ||
       capacitor_voltage === nothing || initial_voltage === nothing
        return result
    end
    push_element!(
        result,
        tokens,
        single_phase_breqiv_history_injection(
            a,
            b,
            resistance,
            inductance,
            capacitance,
            damping_resistance,
            previous_current,
            capacitor_voltage,
            initial_voltage,
        ),
        line_no,
    )
    return record_card!(result, :breqiv)
end

function parse_breqiv3!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 19) || return result
    a1 = node_id!(result, tokens[3])
    b1 = node_id!(result, tokens[4])
    a2 = node_id!(result, tokens[5])
    b2 = node_id!(result, tokens[6])
    a3 = node_id!(result, tokens[7])
    b3 = node_id!(result, tokens[8])
    zero_r = parse_float!(result, tokens[9], line_no, "zero_resistance")
    zero_l = parse_float!(result, tokens[10], line_no, "zero_inductance")
    zero_c = parse_float!(result, tokens[11], line_no, "zero_capacitance")
    zero_rl = parse_float!(result, tokens[12], line_no, "zero_damping_resistance")
    positive_r = parse_float!(result, tokens[13], line_no, "positive_resistance")
    positive_l = parse_float!(result, tokens[14], line_no, "positive_inductance")
    positive_c = parse_float!(result, tokens[15], line_no, "positive_capacitance")
    positive_rl = parse_float!(result, tokens[16], line_no, "positive_damping_resistance")
    initial_v1 = parse_float!(result, tokens[17], line_no, "initial_v1")
    initial_v2 = parse_float!(result, tokens[18], line_no, "initial_v2")
    initial_v3 = parse_float!(result, tokens[19], line_no, "initial_v3")
    zero_previous_current = optional_float!(result, tokens, line_no, 20, "zero_previous_current", 0.0)
    zero_capacitor_voltage = optional_float!(result, tokens, line_no, 21, "zero_capacitor_voltage", 0.0)
    positive_previous_current_1 = optional_float!(result, tokens, line_no, 22, "positive_previous_current_1", 0.0)
    positive_capacitor_voltage_1 = optional_float!(result, tokens, line_no, 23, "positive_capacitor_voltage_1", 0.0)
    positive_previous_current_2 = optional_float!(result, tokens, line_no, 24, "positive_previous_current_2", 0.0)
    positive_capacitor_voltage_2 = optional_float!(result, tokens, line_no, 25, "positive_capacitor_voltage_2", 0.0)
    if zero_r === nothing || zero_l === nothing || zero_c === nothing ||
       zero_rl === nothing || positive_r === nothing || positive_l === nothing ||
       positive_c === nothing || positive_rl === nothing || initial_v1 === nothing ||
       initial_v2 === nothing || initial_v3 === nothing
        return result
    end
    push_element!(
        result,
        tokens,
        three_phase_breqiv_history_injection(
            a1,
            b1,
            a2,
            b2,
            a3,
            b3,
            zero_r,
            zero_l,
            zero_c,
            zero_rl,
            positive_r,
            positive_l,
            positive_c,
            positive_rl,
            initial_v1,
            initial_v2,
            initial_v3,
            zero_previous_current=zero_previous_current,
            zero_capacitor_voltage=zero_capacitor_voltage,
            positive_previous_current_1=positive_previous_current_1,
            positive_capacitor_voltage_1=positive_capacitor_voltage_1,
            positive_previous_current_2=positive_previous_current_2,
            positive_capacitor_voltage_2=positive_capacitor_voltage_2,
        ),
        line_no,
    )
    return record_card!(result, :breqiv3)
end

function parse_bergeron_line!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 7) || return result
    a = node_id!(result, tokens[3])
    b = node_id!(result, tokens[4])
    zc = parse_positive_float!(result, tokens[5], line_no, "surge_impedance")
    travel_time = parse_positive_float!(result, tokens[6], line_no, "travel_time_s")
    dt = parse_positive_float!(result, tokens[7], line_no, "dt_s")
    attenuation = optional_float!(result, tokens, line_no, 8, "attenuation", 1.0)
    if zc === nothing || travel_time === nothing || dt === nothing
        return result
    end
    try
        push_element!(result, tokens, BergeronLine(a, b, zc, travel_time, dt; attenuation = attenuation), line_no)
        push_bergeron_line_row!(result, length(result.elements), line_no)
    catch err
        add_issue!(result.validation, invalid_value("line $line_no", sprint(showerror, err)))
        return result
    end
    return record_card!(result, :bergeron_line)
end

function parse_switch!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 5) || return result
    a = node_id!(result, tokens[3])
    b = node_id!(result, tokens[4])
    closed = parse_bool!(result, tokens[5], line_no, "state")
    closed === nothing && return result
    on_conductance = optional_float!(result, tokens, line_no, 6, "on_conductance", 1.0e9)
    off_conductance = optional_float!(result, tokens, line_no, 7, "off_conductance", 0.0)
    push_element!(result, tokens,
                  IdealSwitch(a, b, closed; on_conductance=on_conductance,
                              off_conductance=off_conductance), line_no)
    return record_card!(result, :switch)
end

function parse_time_switch!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 7) || return result
    a = node_id!(result, tokens[3])
    b = node_id!(result, tokens[4])
    close_time = parse_float!(result, tokens[5], line_no, "close_time_s")
    open_time = parse_float!(result, tokens[6], line_no, "open_time_s")
    initially_closed = parse_bool!(result, tokens[7], line_no, "initially_closed")
    if close_time === nothing || open_time === nothing || initially_closed === nothing
        return result
    end
    on_conductance = optional_float!(result, tokens, line_no, 8, "on_conductance", 1.0e9)
    off_conductance = optional_float!(result, tokens, line_no, 9, "off_conductance", 0.0)
    push_element!(result, tokens,
                  TimeSwitch(a, b; close_time_s=close_time, open_time_s=open_time,
                             initially_closed=initially_closed,
                             on_conductance=on_conductance,
                             off_conductance=off_conductance), line_no)
    return record_card!(result, :time_switch)
end

function parse_over16_output!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 3) || return result
    push!(
        result.over16_output_channels,
        DeckOVER16OutputChannel(Symbol(String(tokens[2])), Symbol(String(tokens[3])), line_no),
    )
    return record_card!(result, :over16_output)
end

function parse_over16_branch_voltage_output!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 3) || return result
    push!(
        result.over16_branch_voltage_outputs,
        DeckOVER16BranchVoltageOutput(Symbol(String(tokens[2])), Symbol(String(tokens[3])), line_no),
    )
    return record_card!(result, :over16_branch_voltage_output)
end

function parse_over16_branch_current_output!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 3) || return result
    push!(
        result.over16_branch_current_outputs,
        DeckOVER16BranchCurrentOutput(Symbol(String(tokens[2])), Symbol(String(tokens[3])), line_no),
    )
    return record_card!(result, :over16_branch_current_output)
end

function parse_over16_branch_power_output!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 3) || return result
    push!(
        result.over16_branch_power_outputs,
        DeckOVER16BranchPowerOutput(Symbol(String(tokens[2])), Symbol(String(tokens[3])), line_no),
    )
    return record_card!(result, :over16_branch_power_output)
end

function parse_over16_source_value_row!(result::DeckParseResult, tokens, line_no::Int,
                                        card_name::AbstractString)
    require_fields!(result, tokens, line_no, 2) || return nothing
    if length(tokens) > 11
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "$card_name rows accept at most ten numeric values"),
        )
        return nothing
    end
    values = Float64[]
    provided_value_count = length(tokens) - 1
    for index in 2:length(tokens)
        value = parse_float!(result, tokens[index], line_no, "$(card_name)_value_$(index - 1)")
        value === nothing && return nothing
        push!(values, Float64(value))
    end
    while length(values) < 10
        push!(values, 0.0)
    end
    return values, provided_value_count
end

function parse_over16_source_card!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 3) || return result
    if length(tokens) > 12
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "OVER16 source-card rows accept at most ten numeric values"),
        )
        return result
    end
    card_kind = normalized_deck_token(tokens[2])
    kind = if card_kind in ("fixed", "fixed_field")
        :fixed_field
    elseif card_kind in ("free", "free_field")
        :free_field
    else
        add_issue!(
            result.validation,
            unknown_field(
                "line $line_no",
                "Unsupported OVER16 source-card kind $(String(tokens[2])); expected fixed_field or free_field",
            ),
        )
        return result
    end
    value_row = parse_over16_source_value_row!(
        result,
        tokens[2:end],
        line_no,
        "source_card",
    )
    value_row === nothing && return result
    values, provided_value_count = value_row
    push!(
        result.over16_source_card_rows,
        DeckOVER16SourceCardRow(kind, values, provided_value_count, line_no),
    )
    record_card!(result, :over16_source_card)
    provided_value_count == 10 || record_card!(result, :over16_source_card_zero_filled)
    return record_card!(result, Symbol("over16_source_card_", String(kind)))
end

function parse_over16_source_interpolation!(result::DeckParseResult, tokens, line_no::Int)
    value_row = parse_over16_source_value_row!(
        result,
        tokens,
        line_no,
        "source_interpolation",
    )
    value_row === nothing && return result
    values, provided_value_count = value_row
    push!(
        result.over16_source_interpolation_rows,
        DeckOVER16SourceInterpolationRow(values, provided_value_count, line_no),
    )
    record_card!(result, :over16_source_interpolation)
    provided_value_count == 10 || record_card!(result, :over16_source_interpolation_zero_filled)
    return result
end

function parse_over16_source_tacs_override!(result::DeckParseResult, tokens, line_no::Int)
    require_fields!(result, tokens, line_no, 3) || return result
    if length(tokens) > 3
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "OVER16 source TACS override rows accept position and XTCS index only"),
        )
        return result
    end
    position = parse_int!(result, tokens[2], line_no, "source_tacs_override_position")
    xtcs_index = parse_int!(result, tokens[3], line_no, "source_tacs_override_xtcs_index")
    position === nothing && return result
    xtcs_index === nothing && return result
    if position < 1 || position > 10
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "source_tacs_override_position=$(String(tokens[2])): expected 1 through 10"),
        )
        return result
    end
    if xtcs_index < 0
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "source_tacs_override_xtcs_index=$(String(tokens[3])): expected nonnegative Int"),
        )
        return result
    end
    push!(
        result.over16_source_tacs_override_rows,
        DeckOVER16SourceTACSOverrideRow(position, xtcs_index, line_no),
    )
    return record_card!(result, :over16_source_tacs_override)
end

function parse_over16_source_analytic!(result::DeckParseResult, tokens, line_no::Int)
    value_row = parse_over16_source_value_row!(
        result,
        tokens,
        line_no,
        "source_analytic",
    )
    value_row === nothing && return result
    values, provided_value_count = value_row
    push!(
        result.over16_source_analytic_rows,
        DeckOVER16SourceAnalyticRow(values, provided_value_count, line_no),
    )
    record_card!(result, :over16_source_analytic)
    provided_value_count == 10 || record_card!(result, :over16_source_analytic_zero_filled)
    return result
end

function validate_over16_output_channels!(result::DeckParseResult)
    for channel in result.over16_output_channels
        if !haskey(result.node_map, channel.node)
            add_issue!(
                result.validation,
                invalid_value(
                    "line $(channel.line_no)",
                    "over16_output node $(String(channel.node)) is not declared by any supported bus or element card",
                ),
            )
        end
    end
    return result
end

function accepted_branch_output_element(element)::Bool
    return element isa ConductanceBranch ||
           element isa SeriesRLBranch ||
           element isa SeriesRLCBranch ||
           element isa CapacitorBranch
end

function validate_over16_branch_voltage_outputs!(result::DeckParseResult)
    for request in result.over16_branch_voltage_outputs
        element_index = findfirst(==(request.branch), result.element_names)
        if element_index === nothing
            add_issue!(
                result.validation,
                invalid_value(
                    "line $(request.line_no)",
                    "over16 branch voltage output branch $(String(request.branch)) is not declared by any supported branch card",
                ),
            )
        elseif !accepted_branch_output_element(result.elements[element_index])
            add_issue!(
                result.validation,
                invalid_value(
                    "line $(request.line_no)",
                    "over16 branch voltage output branch $(String(request.branch)) does not have an accepted scalar branch companion owner",
                ),
            )
        end
    end
    return result
end

function validate_over16_branch_current_outputs!(result::DeckParseResult)
    for request in result.over16_branch_current_outputs
        element_index = findfirst(==(request.branch), result.element_names)
        if element_index === nothing
            add_issue!(
                result.validation,
                invalid_value(
                    "line $(request.line_no)",
                    "over16 branch current output branch $(String(request.branch)) is not declared by any supported branch card",
                ),
            )
        elseif !accepted_branch_output_element(result.elements[element_index])
            add_issue!(
                result.validation,
                invalid_value(
                    "line $(request.line_no)",
                    "over16 branch current output branch $(String(request.branch)) does not have an accepted scalar branch companion owner",
                ),
            )
        end
    end
    return result
end

function validate_over16_branch_power_outputs!(result::DeckParseResult)
    for request in result.over16_branch_power_outputs
        element_index = findfirst(==(request.branch), result.element_names)
        if element_index === nothing
            add_issue!(
                result.validation,
                invalid_value(
                    "line $(request.line_no)",
                    "over16 branch power output branch $(String(request.branch)) is not declared by any supported branch card",
                ),
            )
        elseif !accepted_branch_output_element(result.elements[element_index])
            add_issue!(
                result.validation,
                invalid_value(
                    "line $(request.line_no)",
                    "over16 branch power output branch $(String(request.branch)) does not have an accepted scalar branch companion owner",
                ),
            )
        end
    end
    return result
end

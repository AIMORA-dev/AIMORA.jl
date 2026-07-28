
function validate_coupled_phase_pi_section_rows!(result::DeckParseResult)
    groups = _coupled_phase_pi_row_groups(result.coupled_phase_pi_section_rows)
    accepted_sections = CoupledLumpedPhasePiSection[]
    for (index, group) in enumerate(groups)
        first_row = first(group)
        try
            section =
                _coupled_phase_pi_section_from_rows(result, group, index, accepted_sections)
            push!(accepted_sections, section)
            record_card!(result, :fixed_card_coupled_phase_pi_section)
            first_row.reference_kind == :none ||
                record_card!(result, :fixed_card_coupled_phase_pi_copy_section)
        catch err
            record_fixed_blocker!(
                result,
                :bpa_fixed_branch_blocked,
                :fixed_card_coupled_phase_pi_section_deferred,
            )
            add_issue!(
                result.validation,
                unknown_field(
                    "line $(first_row.line_no)",
                    "Deferred fixed-card coupled phase PI section: $(sprint(showerror, err))",
                ),
            )
        end
    end
    return result
end

function validate_cascaded_pi_requests!(result::DeckParseResult)
    for (index, row) in enumerate(result.cascaded_pi_request_rows)
        try
            _deck_cascaded_phase_pi_equivalent(result, row, index)
        catch err
            add_issue!(result.validation, unknown_field(
                "line $(row.header_line_no)",
                "Invalid cascaded PI request: $(sprint(showerror, err))",
            ))
        end
    end
    return result
end

function push_bpa_fixed_switch_time_row!(
    result::DeckParseResult,
    switch_type::Int,
    from_node::AbstractString,
    to_node::AbstractString,
    close_time::Float64,
    open_time::Float64,
    closed_text::AbstractString,
    line_no::Int,
    initial_issues::Int;
    switch_layout_count::Union{Nothing,Symbol}=nothing,
    switch_layout_kind::Symbol=:fixed_field,
    marker_text::AbstractString="",
    output_code::Int=0,
)::Bool
    if switch_type != 0
        blocker = switch_type in (11, 12, 13) ?
                  :bpa_fixed_switch_blocked_tacs_controlled :
                  :bpa_fixed_switch_blocked_other
        record_fixed_blocker!(result, :bpa_fixed_switch_blocked, blocker)
        add_issue!(result.validation,
                   unknown_field("line $line_no",
                                 switch_type in (11, 12, 13) ?
                                 "Unsupported OVER5 switch type $switch_type: switch types 11-13 require TACS/controlled-switch state ownership" :
                                 "Unsupported OVER5 switch type $switch_type"))
        return true
    end
    closed_marker = uppercase(strip(String(closed_text)))
    normalized_marker_text = uppercase(strip(String(marker_text)))
    measuring_marker =
        closed_marker == "MEASURING" ||
        string(closed_marker, normalized_marker_text) == "MEASURING"
    initially_closed = measuring_marker || closed_marker == "CLOSED" || close_time < 0.0
    accepted_close_time = (measuring_marker || close_time < 0.0) ? 0.0 : close_time
    accepted_open_time = measuring_marker ? Inf : open_time
    name = bpa_fixed_owner_name(result, "switch")
    element_count_before = length(result.elements)
    parse_time_switch!(
        result,
        [
            "time_switch",
            name,
            from_node,
            to_node,
            string(accepted_close_time),
            string(accepted_open_time),
            string(initially_closed),
        ],
        line_no,
    )
    record_fixed_card!(result, :bpa_fixed_switch, :bpa_fixed_switch_time, initial_issues)
    if length(result.validation.issues) == initial_issues
        switch_layout_count === nothing || record_card!(result, switch_layout_count)
        initially_closed && record_card!(result, :bpa_fixed_switch_initially_closed)
        measuring_marker && record_card!(result, :bpa_fixed_switch_measuring)
        if length(result.elements) == element_count_before + 1 && result.elements[end] isa TimeSwitch
            switch = result.elements[end]
            node_names = _deck_node_names_for_indices(result, [switch.a, switch.b])
            push!(
                result.over5_switch_rows,
                DeckOVER5SwitchRow(
                    result.element_names[end],
                    node_names[1],
                    node_names[2],
                    switch.a,
                    switch.b,
                    line_no,
                    switch_type,
                    output_code,
                    switch_layout_kind,
                    close_time,
                    open_time,
                    switch.close_time_s,
                    switch.open_time_s,
                    switch.initially_closed,
                    measuring_marker,
                    closed_marker,
                    normalized_marker_text,
                    switch.on_conductance,
                    switch.off_conductance,
                ),
            )
        end
    end
    return true
end

function bpa_fixed_switch_free_field_row_candidate(line::AbstractString)::Bool
    occursin(',', String(line)) || return false
    fields = bpa_fixed_switch_free_field_fields(line)
    isempty(fields) && return false
    uppercase(fields[1]) == "NAME" && return true
    return tryparse(Int, fields[1]) !== nothing
end

function parse_bpa_fixed_switch_free_field_card!(result::DeckParseResult,
                                                 line::AbstractString,
                                                 line_no::Int)::Bool
    fields = bpa_fixed_switch_free_field_fields(line)
    initial_issues = length(result.validation.issues)
    if isempty(fields)
        return true
    elseif uppercase(fields[1]) == "NAME"
        if length(fields) < 2 || isempty(strip(fields[2]))
            add_issue!(result.validation,
                       missing_data("line $line_no",
                                    "expected OVER5 free-field switch moniker in field 2"))
            return true
        end
        enqueue_fixed_owner_name!(result, :switch, fields[2])
        record_card!(result, :fixed_field)
        record_card!(result, :bpa_fixed_switch_name)
        record_card!(result, :bpa_fixed_switch_free_field)
        return true
    end
    switch_type = free_field_int_or_default!(result, fields, 1, line_no, "switch_type", 0)
    switch_type === nothing && return true
    maximum_field_count = switch_type in (11, 12, 13) ? 13 : 7
    if length(fields) > maximum_field_count
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "OVER5 free-field switch type $switch_type accepts at most $maximum_field_count fields",
            ),
        )
        return true
    end
    from_node = length(fields) >= 2 ? strip(fields[2]) : ""
    to_node = length(fields) >= 3 ? strip(fields[3]) : ""
    if isempty(from_node) || isempty(to_node)
        add_issue!(result.validation,
                   missing_data("line $line_no", "expected OVER5 free-field switch bus names in fields 2 and 3"))
        return true
    end
    if switch_type in (11, 12, 13)
        optional_float(index, field) = begin
            index > length(fields) && return missing
            raw = strip(fields[index])
            isempty(raw) && return missing
            value = tryparse_deck_float(raw)
            if value === nothing
                add_issue!(
                    result.validation,
                    invalid_value("line $line_no", "$field=$raw: expected free-field Float64"),
                )
                return missing
            end
            return Float64(value)
        end
        ignition_voltage = optional_float(4, "switch_ignition_voltage")
        holding_current = optional_float(5, "switch_holding_current")
        deionization_time_s = optional_float(6, "switch_deionization_time_s")
        switching_delay_or_model = optional_float(7, "switch_switching_delay_or_model")
        length(result.validation.issues) == initial_issues || return true
        initial_state = length(fields) >= 8 ? fields[8] : ""
        parameter_marker = length(fields) >= 9 ? fields[9] : ""
        gate_signal = length(fields) >= 10 ? fields[10] : ""
        clamp_signal = length(fields) >= 11 ? fields[11] : ""
        event_output_code = free_field_int_or_default!(
            result,
            fields,
            12,
            line_no,
            "switch_event_output_code",
            0,
        )
        output_code = free_field_int_or_default!(
            result,
            fields,
            13,
            line_no,
            "switch_output_code",
            0,
        )
        event_output_code === nothing && return true
        output_code === nothing && return true
        return push_control_system_switch_row!(
            result,
            line,
            line_no,
            switch_type,
            from_node,
            to_node,
            ignition_voltage,
            holding_current,
            deionization_time_s,
            switching_delay_or_model,
            initial_state,
            parameter_marker,
            gate_signal,
            clamp_signal,
            event_output_code,
            output_code,
            :free_field,
        )
    end
    close_time = free_field_float_or_default!(result, fields, 4, line_no, "switch_close_time", 0.0)
    open_time = free_field_float_or_default!(result, fields, 5, line_no, "switch_open_time", 0.0)
    if close_time === nothing || open_time === nothing
        return true
    end
    closed_text = length(fields) >= 6 ? fields[6] : ""
    marker_text = length(fields) >= 7 ? fields[7] : ""
    return push_bpa_fixed_switch_time_row!(
        result,
        switch_type,
        from_node,
        to_node,
        close_time,
        open_time,
        closed_text,
        line_no,
        initial_issues;
        switch_layout_count = :bpa_fixed_switch_free_field,
        switch_layout_kind = :free_field,
        marker_text = marker_text,
    )
end

function fixed_card_switch_type!(result::DeckParseResult, image::AbstractString,
                                 line_no::Int, from_node::AbstractString,
                                 to_node::AbstractString)
    raw_type = fixed_field(image, 1, 2)
    if !isempty(raw_type)
        return fixed_int_field!(result, image, line_no, 1, 2, "switch_type")
    end
    if isempty(from_node) || isempty(to_node)
        return fixed_int_field!(result, image, line_no, 1, 2, "switch_type")
    end
    record_card!(result, :fixed_card_blank_switch_type_default)
    return 0
end

function fixed_card_switch_output_code!(result::DeckParseResult, image::AbstractString,
                                        line_no::Int)
    output_code = fixed_int_or_default!(result, image, line_no, 80, 80,
                                        "switch_output_code", 0)
    output_code === nothing && return nothing
    output_code != 0 && record_card!(result, :fixed_card_switch_output_code)
    return output_code
end

function fixed_control_system_switch_optional_float(image::AbstractString,
                                                    first_col::Int,
                                                    last_col::Int)
    raw = fixed_field(image, first_col, last_col)
    isempty(raw) && return missing
    value = tryparse_deck_float(raw)
    return value === nothing ? missing : value
end

function push_control_system_switch_row!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    switch_type::Int,
    from_node::AbstractString,
    to_node::AbstractString,
    ignition_voltage::Union{Missing,Float64},
    holding_current::Union{Missing,Float64},
    deionization_time_s::Union{Missing,Float64},
    switching_delay_or_model::Union{Missing,Float64},
    initial_state::AbstractString,
    parameter_marker::AbstractString,
    gate_signal_text::AbstractString,
    clamp_signal_text::AbstractString,
    event_output_code::Int,
    output_code::Int,
    layout_kind::Symbol,
)::Bool
    if isempty(from_node) || isempty(to_node)
        record_fixed_blocker!(
            result,
            :bpa_fixed_switch_blocked,
            :bpa_fixed_switch_blocked_tacs_controlled,
        )
        add_issue!(
            result.validation,
            missing_data("line $line_no", "expected control-system switch endpoints"),
        )
        return true
    end
    parameter_source_kind = :card_values
    parameter_reference_index = 0
    parameter_reference_line_no = 0
    if uppercase(strip(String(parameter_marker))) == "SAME"
        if switch_type != 11 || isempty(result.control_system_switch_rows)
            record_fixed_blocker!(
                result,
                :bpa_fixed_switch_blocked,
                :bpa_fixed_switch_missing_same_reference,
            )
            add_issue!(
                result.validation,
                invalid_value(
                    "line $line_no",
                    "SAME controlled-valve parameters require a preceding type-11/12/13 switch",
                ),
            )
            return true
        end
        reference_index = length(result.control_system_switch_rows)
        reference = result.control_system_switch_rows[reference_index]
        ignition_voltage = reference.ignition_voltage
        holding_current = reference.holding_current
        deionization_time_s = reference.deionization_time_s
        parameter_source_kind = :prior_switch
        parameter_reference_index = reference_index
        parameter_reference_line_no = reference.line_no
        record_card!(result, :control_system_switch_same_parameters)
    end

    gate_signal = if switch_type == 13 || isempty(strip(gate_signal_text))
        missing
    else
        Symbol(strip(gate_signal_text))
    end
    clamp_signal = isempty(strip(clamp_signal_text)) ?
        missing : Symbol(strip(clamp_signal_text))
    control_signal = if switch_type == 13
        clamp_signal === missing ? :ALWAYS_DISABLED : clamp_signal
    elseif gate_signal !== missing
        gate_signal
    elseif switch_type == 11
        :ALWAYS_ENABLED
    else
        :ALWAYS_DISABLED
    end
    push!(
        result.control_system_switch_rows,
        DeckControlSystemSwitchRow(
            line_no,
            switch_type,
            Symbol(from_node),
            Symbol(to_node),
            ignition_voltage,
            holding_current,
            deionization_time_s,
            String(initial_state),
            control_signal,
            gate_signal,
            clamp_signal,
            switching_delay_or_model,
            nothing,
            parameter_source_kind,
            parameter_reference_index,
            parameter_reference_line_no,
            layout_kind,
            event_output_code,
            output_code,
            String(line),
        ),
    )
    record_card!(result, :control_system_switch_card)
    record_card!(result, layout_kind == :free_field ?
                         :bpa_fixed_switch_free_field :
                         :bpa_fixed_switch_fixed_field)
    event_output_code != 0 && record_card!(result, :control_system_switch_event_output)
    return true
end

function _delayed_arc_continuation_values(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
)
    image = fixed_image(line)
    raw_fields = occursin(',', line) ?
        strip.(split(String(line), ','; keepempty = true)) :
        [
            strip(String(SubString(image, first_col, last_col)))
            for (first_col, last_col) in ((1, 16), (17, 32), (33, 48), (49, 64))
        ]
    length(raw_fields) == 4 || begin
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "delayed-arc continuation requires four numeric fields",
            ),
        )
        return nothing
    end
    values = Float64[]
    for (field_index, text) in enumerate(raw_fields)
        value = tryparse_deck_float(text)
        if value === nothing
            add_issue!(
                result.validation,
                invalid_value(
                    "line $line_no",
                    "delayed-arc field $field_index='$text': expected Float64",
                ),
            )
            return nothing
        end
        push!(values, Float64(value))
    end
    coefficient, exponent, time_scale, cutoff = values
    coefficient >= 0.0 && isfinite(coefficient) || begin
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "delayed-arc current coefficient must be finite and nonnegative"),
        )
        return nothing
    end
    isfinite(exponent) || begin
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "delayed-arc current exponent must be finite"),
        )
        return nothing
    end
    time_scale > 0.0 && isfinite(time_scale) || begin
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "delayed-arc time scale must be finite and positive"),
        )
        return nothing
    end
    cutoff >= 0.0 && isfinite(cutoff) || begin
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "delayed-arc cutoff current must be finite and nonnegative"),
        )
        return nothing
    end
    return DeckDelayedArcSwitchParameters(
        line_no,
        coefficient,
        exponent,
        time_scale,
        cutoff,
        String(line),
    )
end

function parse_delayed_arc_switch_continuation!(
    result::DeckParseResult,
    row_index::Int,
    line::AbstractString,
    line_no::Int,
)
    1 <= row_index <= length(result.control_system_switch_rows) ||
        throw(ArgumentError("delayed-arc switch row index is invalid"))
    row = result.control_system_switch_rows[row_index]
    row.switching_delay_or_model !== missing &&
        row.switching_delay_or_model == 7777.0 ||
        throw(ArgumentError("delayed-arc continuation does not follow a 7777 switch row"))
    parameters = _delayed_arc_continuation_values(result, line, line_no)
    parameters === nothing && return true
    result.control_system_switch_rows[row_index] =
        DeckControlSystemSwitchRow(
            row.line_no,
            row.switch_type,
            row.from_node,
            row.to_node,
            row.ignition_voltage,
            row.holding_current,
            row.deionization_time_s,
            row.initial_state,
            row.control_signal,
            row.gate_signal,
            row.clamp_signal,
            row.switching_delay_or_model,
            parameters,
            row.parameter_source_kind,
            row.parameter_reference_index,
            row.parameter_reference_line_no,
            row.layout_kind,
            row.event_output_code,
            row.output_code,
            row.raw_text,
        )
    record_card!(result, :control_system_switch_delayed_arc_continuation)
    return true
end

function parse_control_system_switch_card!(
    result::DeckParseResult,
    image::AbstractString,
    line::AbstractString,
    line_no::Int,
    switch_type::Int,
    from_node::AbstractString,
    to_node::AbstractString,
    output_code::Int,
)::Bool
    event_output_code = fixed_int_or_default!(
        result,
        image,
        line_no,
        79,
        79,
        "switch_event_output_code",
        0,
    )
    event_output_code === nothing && return true
    return push_control_system_switch_row!(
        result,
        line,
        line_no,
        switch_type,
        from_node,
        to_node,
        fixed_control_system_switch_optional_float(image, 15, 24),
        fixed_control_system_switch_optional_float(image, 25, 34),
        fixed_control_system_switch_optional_float(image, 35, 44),
        fixed_control_system_switch_optional_float(image, 45, 54),
        fixed_field(image, 55, 60),
        fixed_field(image, 61, 64),
        fixed_field(image, 65, 70),
        fixed_field(image, 71, 76),
        event_output_code,
        output_code,
        :fixed_field,
    )
end

function parse_bpa_fixed_switch_card!(result::DeckParseResult, line::AbstractString,
                                      line_no::Int)::Bool
    if bpa_fixed_switch_free_field_row_candidate(line)
        return parse_bpa_fixed_switch_free_field_card!(result, line, line_no)
    end
    image = fixed_image(line)
    initial_issues = length(result.validation.issues)
    from_node = fixed_field(image, 3, 8)
    to_node = fixed_field(image, 9, 14)
    if !isempty(from_node) && isempty(to_node)
        to_node = "0"
        record_card!(result, :fixed_card_single_terminal_switch_reference)
    end
    if uppercase(from_node) == "NAME"
        if isempty(to_node)
            add_issue!(result.validation,
                       missing_data("line $line_no",
                                    "expected OVER5 fixed-field switch moniker in columns 9-14"))
            return true
        end
        enqueue_fixed_owner_name!(result, :switch, to_node)
        record_card!(result, :fixed_field)
        record_card!(result, :bpa_fixed_switch_name)
        return true
    end

    switch_type = fixed_card_switch_type!(result, image, line_no, from_node, to_node)
    output_code = fixed_card_switch_output_code!(result, image, line_no)
    if switch_type in (11, 12, 13) && output_code !== nothing
        return parse_control_system_switch_card!(
            result,
            image,
            line,
            line_no,
            switch_type,
            from_node,
            to_node,
            output_code,
        )
    elseif switch_type in (11, 12, 13)
        return true
    end
    close_time = fixed_float_field!(result, image, line_no, 15, 24, "switch_close_time")
    open_time =
        fixed_float_or_default!(result, image, line_no, 25, 34, "switch_open_time", Inf)
    closed_text = fixed_field(image, 55, 60)
    marker_text = fixed_field(image, 61, 64)
    if switch_type === nothing || output_code === nothing ||
       close_time === nothing || open_time === nothing
        return true
    end
    if isempty(from_node) || isempty(to_node)
        add_issue!(result.validation,
                   missing_data("line $line_no", "expected OVER5 fixed-field switch bus names in columns 3-14"))
        return true
    end
    return push_bpa_fixed_switch_time_row!(
        result,
        switch_type,
        from_node,
        to_node,
        close_time,
        open_time,
        closed_text,
        line_no,
        initial_issues,
        marker_text = marker_text,
        output_code = output_code,
    )
end

function fixed_a6_fields(image::AbstractString)
    fields = String[]
    for first_col in 1:6:73
        value = fixed_field(image, first_col, first_col + 5)
        isempty(value) || push!(fields, value)
    end
    return fields
end

function bpa_fixed_output_dummy_node_name(index::Int)
    return Symbol("CHAN", lpad(string(index), 2, '0'))
end

function parse_bpa_fixed_output_dummy_nodes!(result::DeckParseResult,
                                             image::AbstractString,
                                             line_no::Int)::Bool
    channel_count = fixed_int_field!(result, image, line_no, 9, 16,
                                     "dummy_node_voltage_output_count")
    channel_count === nothing && return true
    if channel_count <= 0 || channel_count > 99
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "OVER15 CHAN01 dummy node-voltage output count must be 1 through 99",
            ),
        )
        return true
    end
    record_card!(result, :fixed_field)
    for index in 1:channel_count
        node = bpa_fixed_output_dummy_node_name(index)
        get!(result.node_map, node, length(result.node_map) + 1)
        initial_issues = length(result.validation.issues)
        parse_over16_output!(result, ["over16_output", String(node), String(node)], line_no)
        if length(result.validation.issues) == initial_issues
            record_card!(result, :bpa_fixed_output)
            record_card!(result, :bpa_fixed_output_voltage)
            record_card!(result, :bpa_fixed_output_dummy_node_voltage)
            push_over15_node_output_request_row!(
                result,
                node,
                node,
                line_no,
                :node_voltage,
                :dummy_node_voltage,
                :fixed_field,
                0,
            )
        end
    end
    return true
end

function bpa_fixed_output_node_pair_fields(image::AbstractString)
    pairs = Tuple{String,String}[]
    for first_col in 3:12:63
        from_node = fixed_field(image, first_col, first_col + 5)
        to_node = fixed_field(image, first_col + 6, first_col + 11)
        isempty(from_node) && isempty(to_node) && continue
        push!(pairs, (from_node, to_node))
    end
    return pairs
end

function bpa_fixed_output_branch_name_fields(image::AbstractString)
    names = String[]
    for first_col in 3:6:75
        name = fixed_field(image, first_col, first_col + 5)
        isempty(name) || push!(names, name)
    end
    return names
end

bpa_fixed_output_i2_node_fields(image::AbstractString) =
    bpa_fixed_output_branch_name_fields(image)

function fixed_card_declared_node_count(result::DeckParseResult,
                                        names::AbstractVector{<:AbstractString})
    return count(name -> haskey(result.node_map, Symbol(String(name))), names)
end

function fixed_card_output_voltage_node_fields(result::DeckParseResult,
                                               image::AbstractString)
    shifted_fields = bpa_fixed_output_branch_name_fields(image)
    aligned_fields = fixed_a6_fields(image)
    isempty(shifted_fields) && return aligned_fields
    isempty(aligned_fields) && return shifted_fields

    shifted_declared = fixed_card_declared_node_count(result, shifted_fields)
    aligned_declared = fixed_card_declared_node_count(result, aligned_fields)
    shifted_missing = length(shifted_fields) - shifted_declared
    aligned_missing = length(aligned_fields) - aligned_declared

    if shifted_missing < aligned_missing
        return shifted_fields
    elseif aligned_missing < shifted_missing
        return aligned_fields
    elseif shifted_declared > aligned_declared
        return shifted_fields
    elseif aligned_declared > shifted_declared
        return aligned_fields
    end
    return isempty(fixed_field(image, 1, 2)) ? shifted_fields : aligned_fields
end

function bpa_fixed_all_node_voltage_channel_name(node::Symbol)
    return Symbol("all_node_voltage_", String(node))
end

function push_over15_node_output_request_row!(
    result::DeckParseResult,
    name::Union{Symbol,AbstractString},
    node::Union{Symbol,AbstractString},
    line_no::Int,
    output_kind::Symbol,
    request_kind::Symbol,
    layout_kind::Symbol,
    output_code::Int,
)
    node_name = Symbol(String(node))
    node_index = get(result.node_map, node_name, 0)
    push!(
        result.over15_output_request_rows,
        DeckOVER15OutputRequestRow(
            Symbol(String(name)),
            output_kind,
            request_kind,
            layout_kind,
            line_no,
            output_code,
            node_name,
            node_index,
            :none,
            0,
        ),
    )
    return result
end

function push_over15_branch_output_request_row!(
    result::DeckParseResult,
    name::Union{Symbol,AbstractString},
    branch::Union{Symbol,AbstractString},
    line_no::Int,
    output_kind::Symbol,
    request_kind::Symbol,
    layout_kind::Symbol,
    output_code::Int,
)
    branch_name = Symbol(String(branch))
    branch_index = findfirst(==(branch_name), result.element_names)
    push!(
        result.over15_output_request_rows,
        DeckOVER15OutputRequestRow(
            Symbol(String(name)),
            output_kind,
            request_kind,
            layout_kind,
            line_no,
            output_code,
            :none,
            0,
            branch_name,
            branch_index === nothing ? 0 : branch_index,
        ),
    )
    return result
end

function record_bpa_fixed_output_free_field_row!(result::DeckParseResult)
    record_card!(result, :fixed_field)
    return record_card!(result, :bpa_fixed_output_free_field)
end

function parse_bpa_free_field_output_dummy_nodes!(result::DeckParseResult,
                                                  fields::Vector{String},
                                                  line_no::Int)::Bool
    if length(fields) < 2 || isempty(strip(fields[2]))
        add_issue!(
            result.validation,
            missing_data("line $line_no", "expected OVER15 CHAN01 free-field dummy-node count"),
        )
        return true
    end
    channel_count = free_field_int_or_default!(
        result,
        fields,
        2,
        line_no,
        "dummy_node_voltage_output_count",
        0,
    )
    channel_count === nothing && return true
    if channel_count <= 0 || channel_count > 99
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "OVER15 CHAN01 dummy node-voltage output count must be 1 through 99",
            ),
        )
        return true
    end
    record_bpa_fixed_output_free_field_row!(result)
    for index in 1:channel_count
        node = bpa_fixed_output_dummy_node_name(index)
        get!(result.node_map, node, length(result.node_map) + 1)
        initial_issues = length(result.validation.issues)
        parse_over16_output!(result, ["over16_output", String(node), String(node)], line_no)
        if length(result.validation.issues) == initial_issues
            record_card!(result, :bpa_fixed_output)
            record_card!(result, :bpa_fixed_output_voltage)
            record_card!(result, :bpa_fixed_output_dummy_node_voltage)
            push_over15_node_output_request_row!(
                result,
                node,
                node,
                line_no,
                :node_voltage,
                :dummy_node_voltage,
                :free_field,
                0,
            )
        end
    end
    return true
end

function parse_bpa_free_field_output_end_sentinel!(result::DeckParseResult,
                                                   line_no::Int)::Bool
    record_bpa_fixed_output_free_field_row!(result)
    record_card!(result, :bpa_fixed_output_end_sentinel)
    record_card!(result, :bpa_fixed_output_interpolation_deferred)
    return true
end

function parse_bpa_free_field_output_all_node_voltage!(result::DeckParseResult,
                                                       fields::Vector{String},
                                                       line_no::Int)::Bool
    node_names = deck_node_names(result)
    record_bpa_fixed_output_free_field_row!(result)
    record_card!(result, :bpa_fixed_output_all_node_voltage_request)
    if isempty(node_names)
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "OVER15 all-node voltage output request requires declared nodes",
            ),
        )
        return true
    end
    if any(index -> index <= length(fields) && !isempty(strip(fields[index])),
           2:length(fields))
        record_card!(result, :bpa_fixed_output_all_node_extra_names_ignored)
    end
    record_card!(result, :bpa_fixed_output_interpolation_deferred)
    for node in node_names
        name = bpa_fixed_all_node_voltage_channel_name(node)
        initial_issues = length(result.validation.issues)
        parse_over16_output!(result, ["over16_output", String(name), String(node)], line_no)
        if length(result.validation.issues) == initial_issues
            record_card!(result, :bpa_fixed_output)
            record_card!(result, :bpa_fixed_output_voltage)
            record_card!(result, :bpa_fixed_output_all_node_voltage)
            push_over15_node_output_request_row!(
                result,
                name,
                node,
                line_no,
                :node_voltage,
                :all_node_voltage,
                :free_field,
                1,
            )
        end
    end
    return true
end

function parse_bpa_free_field_output_i2_node_voltages!(result::DeckParseResult,
                                                       fields::Vector{String},
                                                       line_no::Int,
                                                       output_code::Int)::Bool
    nodes = [strip(field) for field in fields[2:end] if !isempty(strip(field))]
    isempty(nodes) && return parse_bpa_free_field_output_end_sentinel!(result, line_no)
    record_bpa_fixed_output_free_field_row!(result)
    record_card!(result, :bpa_fixed_output_selected_node_voltage_request)
    if get(result.card_counts, :bpa_fixed_output_all_node_voltage_request, 0) > 0
        record_card!(result, :bpa_fixed_output_selected_node_ignored_after_all_node)
        return true
    end
    for node in nodes
        if !haskey(result.node_map, Symbol(node))
            record_card!(result, :bpa_fixed_output_missing_node_ignored)
            continue
        end
        name = bpa_fixed_owner_name(result, "output")
        initial_issues = length(result.validation.issues)
        parse_over16_output!(result, ["over16_output", name, node], line_no)
        if length(result.validation.issues) == initial_issues
            record_card!(result, :bpa_fixed_output)
            record_card!(result, :bpa_fixed_output_voltage)
            record_card!(result, :bpa_fixed_output_selected_node_voltage)
            push_over15_node_output_request_row!(
                result,
                name,
                node,
                line_no,
                :node_voltage,
                :selected_node_voltage,
                :free_field,
                output_code,
            )
        end
    end
    return true
end

function parse_bpa_free_field_output_branch_voltage_pairs!(result::DeckParseResult,
                                                           fields::Vector{String},
                                                           line_no::Int)::Bool
    nodes = [strip(field) for field in fields[2:end]]
    if isempty(nodes) || isodd(length(nodes)) || any(isempty, nodes)
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "expected OVER15 free-field branch-voltage node pairs after output code -5",
            ),
        )
        return true
    end
    record_bpa_fixed_output_free_field_row!(result)
    for index in 1:2:length(nodes)
        branch = bpa_fixed_branch_name_for_node_pair!(
            result,
            nodes[index],
            nodes[index + 1],
            line_no,
        )
        branch === nothing && continue
        output_name = string("branch_voltage_", branch)
        initial_issues = length(result.validation.issues)
        parse_over16_branch_voltage_output!(
            result,
            ["over16_branch_voltage_output", output_name, branch],
            line_no,
        )
        if length(result.validation.issues) == initial_issues
            record_card!(result, :bpa_fixed_output)
            record_card!(result, :bpa_fixed_output_branch_voltage_pair)
            record_card!(result, :output_branch_voltage)
            push_over15_branch_output_request_row!(
                result,
                output_name,
                branch,
                line_no,
                :branch_voltage,
                :branch_voltage_pair,
                :free_field,
                -5,
            )
        end
    end
    return true
end

function parse_bpa_free_field_output_named_branch_outputs!(result::DeckParseResult,
                                                           fields::Vector{String},
                                                           line_no::Int,
                                                           output_code::Int)::Bool
    names = [strip(field) for field in fields[2:end] if !isempty(strip(field))]
    if isempty(names)
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "expected OVER15 free-field branch-output names after negative output code",
            ),
        )
        return true
    end
    record_bpa_fixed_output_free_field_row!(result)
    for raw_name in names
        branch = bpa_fixed_output_branch_owner_name!(result, raw_name, line_no)
        branch === nothing && continue
        if output_code == -1
            parse_bpa_fixed_output_branch_current_name!(
                result,
                branch,
                line_no;
                layout_kind = :free_field,
                request_kind = :branch_current_name,
                output_code = output_code,
            )
        elseif output_code == -2
            parse_bpa_fixed_output_branch_voltage_name!(
                result,
                branch,
                line_no;
                layout_kind = :free_field,
                request_kind = :branch_voltage_name,
                output_code = output_code,
            )
        elseif output_code == -3
            parse_bpa_fixed_output_branch_voltage_name!(
                result,
                branch,
                line_no;
                layout_kind = :free_field,
                request_kind = :branch_voltage_current_name,
                output_code = output_code,
            )
            parse_bpa_fixed_output_branch_current_name!(
                result,
                branch,
                line_no;
                layout_kind = :free_field,
                request_kind = :branch_voltage_current_name,
                output_code = output_code,
            )
        elseif output_code <= -4
            parse_bpa_fixed_output_branch_voltage_name!(
                result,
                branch,
                line_no;
                layout_kind = :free_field,
                request_kind = :branch_power_name,
                output_code = output_code,
            )
            parse_bpa_fixed_output_branch_power_name!(
                result,
                branch,
                line_no;
                layout_kind = :free_field,
                request_kind = :branch_power_name,
                output_code = output_code,
            )
        end
    end
    return true
end

function bpa_fixed_output_free_field_row_candidate(line::AbstractString)::Bool
    return occursin(',', String(line))
end

function parse_bpa_fixed_output_free_field_card!(result::DeckParseResult,
                                                 line::AbstractString,
                                                 line_no::Int)::Bool
    fields = bpa_fixed_source_free_field_fields(line)
    isempty(fields) && return true
    if uppercase(strip(fields[1])) == "CHAN01"
        return parse_bpa_free_field_output_dummy_nodes!(result, fields, line_no)
    end
    output_code = free_field_int_or_default!(result, fields, 1, line_no, "output_code", 0)
    output_code === nothing && return true
    if output_code == 1
        return parse_bpa_free_field_output_all_node_voltage!(result, fields, line_no)
    elseif output_code >= 0
        return parse_bpa_free_field_output_i2_node_voltages!(result, fields, line_no, output_code)
    elseif output_code == -5
        return parse_bpa_free_field_output_branch_voltage_pairs!(result, fields, line_no)
    end
    return parse_bpa_free_field_output_named_branch_outputs!(
        result,
        fields,
        line_no,
        output_code,
    )
end

function parse_bpa_fixed_output_all_node_voltage!(result::DeckParseResult,
                                                  image::AbstractString,
                                                  line_no::Int)::Bool
    node_names = deck_node_names(result)
    record_card!(result, :fixed_field)
    record_card!(result, :bpa_fixed_output_all_node_voltage_request)
    if isempty(node_names)
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "OVER15 all-node voltage output request requires declared nodes",
            ),
        )
        return true
    end
    isempty(bpa_fixed_output_branch_name_fields(image)) ||
        record_card!(result, :bpa_fixed_output_all_node_extra_names_ignored)
    record_card!(result, :bpa_fixed_output_interpolation_deferred)
    for node in node_names
        name = bpa_fixed_all_node_voltage_channel_name(node)
        initial_issues = length(result.validation.issues)
        parse_over16_output!(result, ["over16_output", String(name), String(node)], line_no)
        if length(result.validation.issues) == initial_issues
            record_card!(result, :bpa_fixed_output)
            record_card!(result, :bpa_fixed_output_voltage)
            record_card!(result, :bpa_fixed_output_all_node_voltage)
            push_over15_node_output_request_row!(
                result,
                name,
                node,
                line_no,
                :node_voltage,
                :all_node_voltage,
                :fixed_field,
                1,
            )
        end
    end
    return true
end

function parse_bpa_fixed_output_i2_node_voltages!(result::DeckParseResult,
                                                  image::AbstractString,
                                                  line_no::Int,
                                                  output_code::Int)::Bool
    nodes = bpa_fixed_output_i2_node_fields(image)
    if isempty(nodes)
        record_card!(result, :fixed_field)
        record_card!(result, :bpa_fixed_output_end_sentinel)
        record_card!(result, :bpa_fixed_output_interpolation_deferred)
        return true
    end
    record_card!(result, :fixed_field)
    record_card!(result, :bpa_fixed_output_selected_node_voltage_request)
    if get(result.card_counts, :bpa_fixed_output_all_node_voltage_request, 0) > 0
        record_card!(result, :bpa_fixed_output_selected_node_ignored_after_all_node)
        return true
    end
    for node in nodes
        if !haskey(result.node_map, Symbol(node))
            record_card!(result, :bpa_fixed_output_missing_node_ignored)
            continue
        end
        name = bpa_fixed_owner_name(result, "output")
        initial_issues = length(result.validation.issues)
        parse_over16_output!(result, ["over16_output", name, node], line_no)
        if length(result.validation.issues) == initial_issues
            record_card!(result, :bpa_fixed_output)
            record_card!(result, :bpa_fixed_output_voltage)
            record_card!(result, :bpa_fixed_output_selected_node_voltage)
            push_over15_node_output_request_row!(
                result,
                name,
                node,
                line_no,
                :node_voltage,
                :selected_node_voltage,
                :fixed_field,
                output_code,
            )
        end
    end
    return true
end

function bpa_fixed_output_node_index(result::DeckParseResult,
                                     node::AbstractString,
                                     line_no::Int)
    normalized = lowercase(String(node))
    if normalized in ("0", "gnd", "ground", "ref")
        return 0
    end
    name = Symbol(String(node))
    if haskey(result.node_map, name)
        return result.node_map[name]
    end
    add_issue!(
        result.validation,
        invalid_value(
            "line $line_no",
            "OVER15 fixed-field branch-voltage output node $(String(node)) is not declared by any supported bus or branch card",
        ),
    )
    return nothing
end

function bpa_fixed_output_branch_owner_name!(result::DeckParseResult,
                                             branch_name::AbstractString,
                                             line_no::Int)
    name = Symbol(String(branch_name))
    element_index = findfirst(==(name), result.element_names)
    if element_index === nothing
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "OVER15 fixed-field branch-output request name $(String(branch_name)) does not match a declared branch owner",
            ),
        )
        return nothing
    end
    if !accepted_branch_output_element(result.elements[element_index])
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "OVER15 fixed-field branch-output request name $(String(branch_name)) does not target an accepted scalar branch owner",
            ),
        )
        return nothing
    end
    return String(branch_name)
end

function parse_bpa_fixed_output_branch_voltage_name!(
    result::DeckParseResult,
    branch::AbstractString,
    line_no::Int;
    layout_kind::Symbol = :fixed_field,
    request_kind::Symbol = :branch_voltage_name,
    output_code::Int = -2,
)
    output_name = string("branch_voltage_", branch)
    initial_issues = length(result.validation.issues)
    parse_over16_branch_voltage_output!(
        result,
        ["over16_branch_voltage_output", output_name, branch],
        line_no,
    )
    if length(result.validation.issues) == initial_issues
        record_card!(result, :bpa_fixed_output)
        record_card!(result, :bpa_fixed_output_branch_voltage_name)
        record_card!(result, :output_branch_voltage)
        push_over15_branch_output_request_row!(
            result,
            output_name,
            branch,
            line_no,
            :branch_voltage,
            request_kind,
            layout_kind,
            output_code,
        )
    end
    return result
end

function parse_bpa_fixed_output_branch_current_name!(
    result::DeckParseResult,
    branch::AbstractString,
    line_no::Int;
    layout_kind::Symbol = :fixed_field,
    request_kind::Symbol = :branch_current_name,
    output_code::Int = -1,
)
    output_name = string("branch_current_", branch)
    initial_issues = length(result.validation.issues)
    parse_over16_branch_current_output!(
        result,
        ["over16_branch_current_output", output_name, branch],
        line_no,
    )
    if length(result.validation.issues) == initial_issues
        record_card!(result, :bpa_fixed_output)
        record_card!(result, :bpa_fixed_output_branch_current_name)
        record_card!(result, :output_branch_current)
        push_over15_branch_output_request_row!(
            result,
            output_name,
            branch,
            line_no,
            :branch_current,
            request_kind,
            layout_kind,
            output_code,
        )
    end
    return result
end

function parse_bpa_fixed_output_branch_power_name!(
    result::DeckParseResult,
    branch::AbstractString,
    line_no::Int;
    layout_kind::Symbol = :fixed_field,
    request_kind::Symbol = :branch_power_name,
    output_code::Int = -4,
)
    output_name = string("branch_power_", branch)
    initial_issues = length(result.validation.issues)
    parse_over16_branch_power_output!(
        result,
        ["over16_branch_power_output", output_name, branch],
        line_no,
    )
    if length(result.validation.issues) == initial_issues
        record_card!(result, :bpa_fixed_output)
        record_card!(result, :bpa_fixed_output_branch_power_name)
        record_card!(result, :output_branch_power_energy)
        push_over15_branch_output_request_row!(
            result,
            output_name,
            branch,
            line_no,
            :branch_power_energy,
            request_kind,
            layout_kind,
            output_code,
        )
    end
    return result
end

function bpa_fixed_branch_name_for_node_pair!(result::DeckParseResult,
                                              from_node::AbstractString,
                                              to_node::AbstractString,
                                              line_no::Int)
    if isempty(from_node) || isempty(to_node)
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "expected OVER15 fixed-field branch-voltage node pairs in A6 columns",
            ),
        )
        return nothing
    end
    from_index = bpa_fixed_output_node_index(result, from_node, line_no)
    to_index = bpa_fixed_output_node_index(result, to_node, line_no)
    if from_index === nothing || to_index === nothing
        return nothing
    end
    for (element_index, element) in enumerate(result.elements)
        accepted_branch_output_element(element) || continue
        if element.a == from_index && element.b == to_index
            return String(result.element_names[element_index])
        end
    end
    add_issue!(
        result.validation,
        invalid_value(
            "line $line_no",
            "OVER15 fixed-field branch-voltage output node pair $(String(from_node))-$(String(to_node)) does not match an accepted scalar branch owner",
        ),
    )
    return nothing
end

function parse_bpa_fixed_output_branch_voltage_pairs!(result::DeckParseResult,
                                                      image::AbstractString,
                                                      line_no::Int)::Bool
    pairs = bpa_fixed_output_node_pair_fields(image)
    if isempty(pairs)
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "expected OVER15 fixed-field branch-voltage node pairs after output code -5",
            ),
        )
        return true
    end
    record_card!(result, :fixed_field)
    for (from_node, to_node) in pairs
        branch = bpa_fixed_branch_name_for_node_pair!(result, from_node, to_node, line_no)
        branch === nothing && continue
        output_name = string("branch_voltage_", branch)
        initial_issues = length(result.validation.issues)
        parse_over16_branch_voltage_output!(
            result,
            ["over16_branch_voltage_output", output_name, branch],
            line_no,
        )
        if length(result.validation.issues) == initial_issues
            record_card!(result, :bpa_fixed_output)
            record_card!(result, :bpa_fixed_output_branch_voltage_pair)
            record_card!(result, :output_branch_voltage)
            push_over15_branch_output_request_row!(
                result,
                output_name,
                branch,
                line_no,
                :branch_voltage,
                :branch_voltage_pair,
                :fixed_field,
                -5,
            )
        end
    end
    return true
end

function parse_bpa_fixed_output_named_branch_outputs!(result::DeckParseResult,
                                                      image::AbstractString,
                                                      line_no::Int,
                                                      output_code::Int)::Bool
    names = bpa_fixed_output_branch_name_fields(image)
    if isempty(names)
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "expected OVER15 fixed-field branch-output names after negative output code",
            ),
        )
        return true
    end
    record_card!(result, :fixed_field)
    for raw_name in names
        branch = bpa_fixed_output_branch_owner_name!(result, raw_name, line_no)
        branch === nothing && continue
        if output_code == -1
            parse_bpa_fixed_output_branch_current_name!(
                result,
                branch,
                line_no;
                layout_kind = :fixed_field,
                request_kind = :branch_current_name,
                output_code = output_code,
            )
        elseif output_code == -2
            parse_bpa_fixed_output_branch_voltage_name!(
                result,
                branch,
                line_no;
                layout_kind = :fixed_field,
                request_kind = :branch_voltage_name,
                output_code = output_code,
            )
        elseif output_code == -3
            parse_bpa_fixed_output_branch_voltage_name!(
                result,
                branch,
                line_no;
                layout_kind = :fixed_field,
                request_kind = :branch_voltage_current_name,
                output_code = output_code,
            )
            parse_bpa_fixed_output_branch_current_name!(
                result,
                branch,
                line_no;
                layout_kind = :fixed_field,
                request_kind = :branch_voltage_current_name,
                output_code = output_code,
            )
        elseif output_code <= -4
            parse_bpa_fixed_output_branch_voltage_name!(
                result,
                branch,
                line_no;
                layout_kind = :fixed_field,
                request_kind = :branch_power_name,
                output_code = output_code,
            )
            parse_bpa_fixed_output_branch_power_name!(
                result,
                branch,
                line_no;
                layout_kind = :fixed_field,
                request_kind = :branch_power_name,
                output_code = output_code,
            )
        end
    end
    return true
end

function node_initial_condition_row(line::AbstractString)::Bool
    image = fixed_image(line)
    code = fixed_int_value(image, 1, 2)
    code in (2, 3) || return false
    return !isempty(fixed_field(image, 3, 8))
end

function parse_node_initial_condition_row!(result::DeckParseResult,
                                           line::AbstractString,
                                           line_no::Int)::Bool
    image = fixed_image(line)
    code = fixed_int_value(image, 1, 2)
    node = fixed_field(image, 3, 8)
    reference_node = fixed_field(image, 9, 14)
    real_value = fixed_float_or_default!(
        result,
        image,
        line_no,
        15,
        29,
        "node_initial_condition_real_value",
        0.0,
    )
    imaginary_value = fixed_float_or_default!(
        result,
        image,
        line_no,
        30,
        44,
        "node_initial_condition_imaginary_value",
        0.0,
    )
    if real_value === nothing || imaginary_value === nothing
        return true
    end
    record_card!(result, :fixed_field)
    record_card!(result, :node_initial_condition)
    condition_kind = :unknown_initial_condition
    if code == 2
        record_card!(result, :node_voltage_initial_condition)
        condition_kind = :node_voltage_initial_condition
    elseif code == 3
        record_card!(result, :node_current_initial_condition)
        condition_kind = :node_current_initial_condition
    end
    node_index = isempty(node) ? 0 : node_id!(result, node)
    reference_node_symbol = isempty(reference_node) ? missing : Symbol(reference_node)
    reference_node_index =
        isempty(reference_node) ? missing : node_id!(result, reference_node)
    push!(
        result.node_initial_condition_rows,
        DeckNodeInitialConditionRow(
            line_no,
            condition_kind,
            Symbol(node),
            node_index,
            reference_node_symbol,
            reference_node_index,
            real_value,
            imaginary_value,
            String(line),
        ),
    )
    return true
end

function parse_bpa_fixed_output_card!(result::DeckParseResult, line::AbstractString,
                                      line_no::Int)::Bool
    if bpa_fixed_output_free_field_row_candidate(line)
        return parse_bpa_fixed_output_free_field_card!(result, line, line_no)
    end
    image = fixed_image(line)
    output_code = fixed_int_value(image, 1, 2)
    if uppercase(fixed_field(image, 3, 8)) == "CHAN01"
        return parse_bpa_fixed_output_dummy_nodes!(result, image, line_no)
    elseif output_code == 1
        return parse_bpa_fixed_output_all_node_voltage!(result, image, line_no)
    elseif output_code !== nothing && output_code >= 0
        return parse_bpa_fixed_output_i2_node_voltages!(result, image, line_no, output_code)
    elseif output_code == -5
        return parse_bpa_fixed_output_branch_voltage_pairs!(result, image, line_no)
    elseif output_code !== nothing && output_code < 0
        return parse_bpa_fixed_output_named_branch_outputs!(result, image, line_no, output_code)
    end
    nodes = fixed_card_output_voltage_node_fields(result, image)
    if isempty(nodes)
        add_issue!(result.validation,
                   missing_data("line $line_no",
                                "expected BPA fixed-field output node names in A6 columns"))
        return true
    end
    record_card!(result, :fixed_field)
    for node in nodes
        name = bpa_fixed_owner_name(result, "output")
        initial_issues = length(result.validation.issues)
        parse_over16_output!(result, ["over16_output", name, node], line_no)
        if length(result.validation.issues) == initial_issues
            record_card!(result, :bpa_fixed_output)
            record_card!(result, :bpa_fixed_output_voltage)
            push_over15_node_output_request_row!(
                result,
                name,
                node,
                line_no,
                :node_voltage,
                :fixed_a6_output_node,
                :fixed_field,
                0,
            )
        end
    end
    return true
end

function _continuation_payload(line::AbstractString)::Union{Nothing,String}
    stripped = lstrip(String(line))
    isempty(stripped) && return nothing
    first(stripped) == '+' || return nothing
    length(stripped) == 1 && return ""
    return strip(stripped[nextind(stripped, firstindex(stripped)):end])
end

function _strip_trailing_continuation_marker(line::AbstractString)
    stripped = rstrip(String(line))
    isempty(stripped) && return stripped, false, :none
    last_index = lastindex(stripped)
    if stripped[last_index] == '&'
        if firstindex(stripped) == last_index
            return "", true, :plus_card
        end
        return strip(stripped[firstindex(stripped):prevind(stripped, last_index)]), true, :plus_card
    elseif stripped[last_index] == '$'
        first_nonblank = findfirst(!isspace, stripped)
        if first_nonblank !== nothing && first_nonblank != last_index
            prefix = strip(stripped[firstindex(stripped):prevind(stripped, last_index)])
            if !isempty(prefix) && !endswith(prefix, ",")
                prefix = string(prefix, ",")
            end
            return prefix, true, :field_card
        end
    end
    return stripped, false, :none
end

function deck_logical_records!(result::DeckParseResult, lines)
    records = Vector{Tuple{Int,String}}()
    pending_line_no = 0
    pending_text = ""
    pending_requires_continuation = false
    pending_continuation_kind = :none
    for (line_no, raw_line) in enumerate(lines)
        line = strip_deck_line(raw_line)
        isempty(line) && continue
        continuation = _continuation_payload(line)
        if continuation !== nothing
            if pending_line_no == 0
                add_issue!(
                    result.validation,
                    missing_data("line $line_no", "continuation card without a preceding physical card"),
                )
                continue
            end
            continuation_text, continues, continuation_kind =
                _strip_trailing_continuation_marker(continuation)
            pending_text = strip(string(pending_text, " ", continuation_text))
            pending_requires_continuation = continues
            pending_continuation_kind = continuation_kind
            record_card!(result, :continuation)
            continue
        end
        if pending_line_no != 0 && pending_requires_continuation &&
           pending_continuation_kind == :field_card
            continuation_text, continues, continuation_kind =
                _strip_trailing_continuation_marker(line)
            pending_text = strip(string(pending_text, " ", continuation_text))
            pending_requires_continuation = continues
            pending_continuation_kind = continuation_kind
            record_card!(result, :continuation)
            continue
        end
        if pending_line_no != 0
            if pending_requires_continuation
                add_issue!(
                    result.validation,
                    missing_data("line $line_no", "expected '+' continuation after trailing '&'"),
                )
            end
            push!(records, (pending_line_no, pending_text))
        end
        pending_text, pending_requires_continuation, pending_continuation_kind =
            _strip_trailing_continuation_marker(line)
        pending_line_no = line_no
    end
    if pending_line_no != 0
        if pending_requires_continuation
            message = pending_continuation_kind == :field_card ?
                "unterminated free-field continuation after trailing '\$'" :
                "unterminated trailing '&' continuation"
            add_issue!(
                result.validation,
                missing_data("line $pending_line_no", message),
            )
        end
        push!(records, (pending_line_no, pending_text))
    end
    return records
end

const BPA_BLANK_SECTION_KINDS = Dict(
    "branch" => :blank_branch,
    "switch" => :blank_switch,
    "source" => :blank_source,
    "output" => :blank_output,
    "plot" => :blank_plot,
    "tacs" => :blank_tacs,
    "machine" => :blank_machine,
    "load" => :blank_load,
    "line" => :blank_line,
    "transformer" => :blank_transformer,
    "nonlinear" => :blank_nonlinear,
    "network" => :blank_network,
    "request" => :blank_request,
    "frequency" => :blank_frequency,
    "statistics" => :blank_statistics,
    "monitor" => :blank_monitor,
    "data" => :blank_data,
    "misc" => :blank_misc,
    "miscellaneous" => :blank_miscellaneous,
    "measurement" => :blank_measurement,
    "universal_machine" => :blank_universal_machine,
    "synchronous_machine" => :blank_synchronous_machine,
)

function parse_deck_card!(result::DeckParseResult, tokens::Vector{SubString{String}},
                          line_no::Int)
    card = normalized_deck_token(tokens[1])
    if card in ("begin", "end", "section") || end_last_data_case_card(tokens)
        parse_control_card!(result, tokens, line_no)
    elseif power_frequency_request_card(tokens)
        parse_power_frequency_request!(result, tokens, line_no)
    elseif universal_machine_dimension_request_card(tokens)
        parse_universal_machine_dimension_request!(result, tokens, line_no)
    elseif output_width_request_card(tokens)
        parse_output_width_request!(result, tokens, line_no)
    elseif peak_voltage_monitor_request_card(tokens)
        parse_peak_voltage_monitor_request!(result, tokens, line_no)
    elseif diagnostic_print_request_card(tokens)
        parse_diagnostic_print_request!(result, tokens, line_no)
    elseif tacs_warning_limit_request_card(tokens)
        parse_tacs_warning_limit_request!(result, tokens, line_no)
    elseif plot_file_request_card(tokens)
        parse_plot_file_request!(result, tokens, line_no)
    elseif switch_logic_request_card(tokens)
        parse_switch_logic_request!(result, tokens, line_no)
    elseif simulation_control_request_card(tokens)
        parse_simulation_control_request!(result, tokens, line_no)
    elseif study_option_request_card(tokens)
        parse_study_option_request!(result, tokens, line_no)
    elseif card == "blank"
        parse_blank_control_card!(result, tokens, line_no)
    elseif card in ("branch", "bpa_branch")
        parse_branch_family_card!(result, tokens, line_no)
    elseif card in ("bus", "node")
        parse_bus!(result, tokens, line_no)
    elseif card in ("source", "bpa_source") && source_family_card(tokens)
        parse_source_family_card!(result, tokens, line_no)
    elseif card in ("current_source", "source_current", "bpa_current_source")
        parse_current_source_family_card!(result, tokens, line_no)
    elseif card == "source"
        parse_source!(result, tokens, line_no)
    elseif card == "current"
        parse_current!(result, tokens, line_no)
    elseif card == "conductance"
        parse_conductance!(result, tokens, line_no)
    elseif card == "resistor"
        parse_resistor!(result, tokens, line_no)
    elseif card == "rl"
        parse_series_rl!(result, tokens, line_no)
    elseif card == "inductor"
        parse_inductor!(result, tokens, line_no)
    elseif card == "capacitor"
        parse_capacitor!(result, tokens, line_no)
    elseif card == "breqiv"
        parse_breqiv!(result, tokens, line_no)
    elseif card == "breqiv3"
        parse_breqiv3!(result, tokens, line_no)
    elseif card in ("bergeron_line", "line")
        parse_bergeron_line!(result, tokens, line_no)
    elseif card in ("switch", "bpa_switch") && switch_family_card(tokens)
        parse_switch_family_card!(result, tokens, line_no)
    elseif card == "switch"
        parse_switch!(result, tokens, line_no)
    elseif card == "time_switch"
        parse_time_switch!(result, tokens, line_no)
    elseif card in ("output", "bpa_output")
        parse_output_family_card!(result, tokens, line_no)
    elseif card in ("over16_output", "over16_voltage_output")
        parse_over16_output!(result, tokens, line_no)
    elseif card in ("over16_branch_voltage", "over16_branch_voltage_output")
        parse_over16_branch_voltage_output!(result, tokens, line_no)
    elseif card in ("over16_branch_current", "over16_branch_current_output")
        parse_over16_branch_current_output!(result, tokens, line_no)
    elseif card in (
        "over16_branch_power",
        "over16_branch_power_output",
        "over16_branch_energy",
        "over16_branch_energy_output",
    )
        parse_over16_branch_power_output!(result, tokens, line_no)
    elseif card == "over16_source_card"
        parse_over16_source_card!(result, tokens, line_no)
    elseif card in (
        "over16_source_interpolation",
        "over16_source_interpolation_values",
        "over16_source_interp",
        "over16_source_interp_values",
    )
        parse_over16_source_interpolation!(result, tokens, line_no)
    elseif card in ("over16_source_tacs_override", "over16_vstacs_override")
        parse_over16_source_tacs_override!(result, tokens, line_no)
    elseif card in ("over16_source_analytic", "over16_source_analytic_values")
        parse_over16_source_analytic!(result, tokens, line_no)
    else
        add_issue!(result.validation,
                   unknown_field("line $line_no",
                                 "Unsupported deck card $(String(tokens[1]))"))
    end
    return result
end

function compact_deck_keyword(token)::String
    return replace(normalized_deck_token(token), "." => "", "_" => "")
end

function power_frequency_request_card(tokens)::Bool
    length(tokens) >= 2 || return false
    return compact_deck_keyword(tokens[1]) == "power" &&
           compact_deck_keyword(tokens[2]) == "frequency"
end

function universal_machine_dimension_request_card(tokens)::Bool
    isempty(tokens) && return false
    compact_deck_keyword(tokens[1]) == "aumd" && return true
    length(tokens) >= 3 || return false
    return compact_deck_keyword(tokens[1]) == "absolute" &&
           compact_deck_keyword(tokens[2]) == "um" &&
           compact_deck_keyword(tokens[3]) == "dimensions"
end

function output_width_request_card(tokens)::Bool
    isempty(tokens) && return false
    first_token = compact_deck_keyword(tokens[1])
    first_token in ("ow1", "ow8") && return true
    length(tokens) >= 3 || return false
    return first_token == "output" &&
           compact_deck_keyword(tokens[2]) == "width" &&
           compact_deck_keyword(tokens[3]) in ("132", "80")
end

function output_width_request_columns(tokens)
    first_token = compact_deck_keyword(tokens[1])
    first_token == "ow1" && return 132
    first_token == "ow8" && return 80
    return parse(Int, compact_deck_keyword(tokens[3]))
end

function deck_phrase_match(tokens, words)::Bool
    length(tokens) >= length(words) || return false
    for (index, word) in enumerate(words)
        compact_deck_keyword(tokens[index]) == word || return false
    end
    return true
end

function begin_new_data_case_card(tokens)::Bool
    return deck_phrase_match(tokens, ("begin", "new", "data", "case"))
end

function end_last_data_case_card(tokens)::Bool
    isempty(tokens) && return false
    compact_deck_keyword(tokens[1]) == "eldc" && return true
    return deck_phrase_match(tokens, ("end", "last", "data", "case"))
end

function abort_data_case_card(tokens)::Bool
    isempty(tokens) && return false
    compact_deck_keyword(tokens[1]) == "adc" && return true
    return deck_phrase_match(tokens, ("abort", "data", "case"))
end

function peak_voltage_monitor_request_card(tokens)::Bool
    isempty(tokens) && return false
    compact_deck_keyword(tokens[1]) == "pvm" && return true
    return deck_phrase_match(tokens, ("peak", "voltage", "monitor"))
end

function diagnostic_print_request_card(tokens)::Bool
    isempty(tokens) && return false
    compact_deck_keyword(tokens[1]) == "adp" && return true
    return deck_phrase_match(tokens, ("alternate", "diagnostic", "printout"))
end

function tacs_warning_limit_request_card(tokens)::Bool
    isempty(tokens) && return false
    compact_deck_keyword(tokens[1]) == "twl" && return true
    return deck_phrase_match(tokens, ("tacs", "warn", "limit"))
end

function plot_file_request_card(tokens)::Bool
    return deck_phrase_match(tokens, ("custom", "plot", "file"))
end

function switch_logic_request_card(tokens)::Bool
    isempty(tokens) && return false
    compact_deck_keyword(tokens[1]) == "msl" && return true
    return deck_phrase_match(tokens, ("modify", "switch", "logic"))
end

function simulation_control_request_card(tokens)::Bool
    isempty(tokens) && return false
    first_token = compact_deck_keyword(tokens[1])
    first_token in ("obc", "cs", "mdc", "rte", "bpvs", "todr", "usst") &&
        return true
    deck_phrase_match(tokens, ("omit", "base", "case")) && return true
    deck_phrase_match(tokens, ("change", "switch")) && return true
    deck_phrase_match(tokens, ("miscellaneous", "data", "cards")) && return true
    deck_phrase_match(tokens, ("redefine", "tolerance", "epsiln")) && return true
    deck_phrase_match(tokens, ("time", "of", "dice", "roll")) && return true
    return deck_phrase_match(tokens, ("user", "supplied", "switch", "times"))
end

const STUDY_OPTION_REQUEST_KEYWORDS = Dict(
    "tl" => :type99_message_limit,
    "fr" => :file_request,
    "cz" => :zinc_oxide_format_conversion,
    "hs" => :hauer_impulse_response_setup,
    "lmfs" => :line_model_frequency_scan,
    "ff" => :free_format_characters,
    "pph" => :plotter_paper_height,
    "plpi" => :printer_lines_per_inch,
    "mvo" => :mode_voltage_output,
    "asu" => :analytic_source_usage,
    "lopo" => :plot_oscillation_limit,
    "lbu" => :linear_bias_usage,
    "an" => :automatic_branch_naming,
    "rb" => :renumbering_bypass,
    "fs" => :frequency_scan,
    "d" => :diagnostic_codes,
    "ui" => :user_identification,
    "hr" => :high_resistance_exponent,
    "ao" => :average_output,
    "sos" => :statistics_output_salvage,
    "zo" => :zinc_oxide_constants,
    "fxs" => :fixed_source_declaration,
)

function study_option_request_kind(tokens)
    isempty(tokens) && return nothing
    first_token = compact_deck_keyword(tokens[1])
    kind = get(STUDY_OPTION_REQUEST_KEYWORDS, first_token, nothing)
    kind !== nothing && return kind
    deck_phrase_match(tokens, ("file", "request")) && return :file_request
    deck_phrase_match(tokens, ("convert", "zno")) && return :zinc_oxide_format_conversion
    deck_phrase_match(tokens, ("hauer", "setup")) && return :hauer_impulse_response_setup
    deck_phrase_match(tokens, ("line", "model", "freq", "scan")) &&
        return :line_model_frequency_scan
    deck_phrase_match(tokens, ("type99", "limit")) && return :type99_message_limit
    deck_phrase_match(tokens, ("free", "format")) && return :free_format_characters
    deck_phrase_match(tokens, ("plotter", "paper", "height")) &&
        return :plotter_paper_height
    deck_phrase_match(tokens, ("printer", "lines", "per", "inch")) &&
        return :printer_lines_per_inch
    deck_phrase_match(tokens, ("mode", "voltage", "output")) &&
        return :mode_voltage_output
    deck_phrase_match(tokens, ("analytic", "sources", "usage")) &&
        return :analytic_source_usage
    deck_phrase_match(tokens, ("limit", "on", "plot", "oscillations")) &&
        return :plot_oscillation_limit
    deck_phrase_match(tokens, ("linear", "bias", "usage")) &&
        return :linear_bias_usage
    deck_phrase_match(tokens, ("auto", "name")) && return :automatic_branch_naming
    deck_phrase_match(tokens, ("renumber", "bypass")) && return :renumbering_bypass
    deck_phrase_match(tokens, ("frequency", "scan")) && return :frequency_scan
    deck_phrase_match(tokens, ("diagnostic",)) && return :diagnostic_codes
    deck_phrase_match(tokens, ("user", "identification")) && return :user_identification
    deck_phrase_match(tokens, ("high", "resistance")) && return :high_resistance_exponent
    deck_phrase_match(tokens, ("average", "output")) && return :average_output
    deck_phrase_match(tokens, ("statistics", "output", "salvage")) &&
        return :statistics_output_salvage
    deck_phrase_match(tokens, ("zinc", "oxide")) && return :zinc_oxide_constants
    deck_phrase_match(tokens, ("fix", "source")) && return :fixed_source_declaration
    return nothing
end

study_option_request_card(tokens)::Bool = study_option_request_kind(tokens) !== nothing

function token_strings(tokens)
    return String.(tokens)
end

function deck_control_label(tokens)
    length(tokens) <= 1 && return ""
    return join(String.(tokens[2:end]), " ")
end

function normalized_deck_token(token)::String
    return replace(lowercase(String(token)), '-' => '_', '/' => '_')
end

function deck_token_value(token)::String
    value = String(token)
    equals = findfirst(==('='), value)
    equals === nothing && return value
    return strip(value[nextind(value, equals):end])
end

function deck_field_value(tokens, aliases)::Union{Nothing,String}
    for token in tokens
        value = String(token)
        equals = findfirst(==('='), value)
        equals === nothing && continue
        field = normalized_deck_token(strip(value[begin:prevind(value, equals)]))
        field in aliases && return strip(value[nextind(value, equals):end])
    end
    return nothing
end

function positional_deck_value(tokens, index::Int)::Union{Nothing,String}
    index > length(tokens) && return nothing
    return deck_token_value(tokens[index])
end

function deck_family_value(tokens, aliases, fallback_index::Int)::Union{Nothing,String}
    value = deck_field_value(tokens, aliases)
    value === nothing || return value
    return positional_deck_value(tokens, fallback_index)
end

function require_family_value!(result::DeckParseResult, tokens, line_no::Int, aliases,
                               fallback_index::Int, field::AbstractString)
    value = deck_family_value(tokens, aliases, fallback_index)
    if value === nothing || isempty(value)
        add_issue!(result.validation,
                   missing_data("line $line_no", "expected $field field"))
        return nothing
    end
    return value
end

function record_family_card!(result::DeckParseResult, family::Symbol, kind::Symbol,
                             initial_issue_count::Int)
    length(result.validation.issues) == initial_issue_count || return result
    record_card!(result, family)
    return record_card!(result, kind)
end

function source_family_card(tokens)::Bool
    length(tokens) >= 2 || return false
    normalized_deck_token(tokens[2]) in
        ("voltage", "thevenin", "vsource", "v", "current", "isource", "i")
end


function record_tacs_dimension_request_marker!(
    result::DeckParseResult,
    allocation_kind::Symbol,
    tokens,
    line_no::Int,
)
    record_control_card!(result, allocation_kind, tokens, line_no)
    record_card!(result, :tacs_dimension_request)
    return record_card!(result, :special_request)
end

function printout_frequency_change_request_marker(tokens)::Bool
    return deck_phrase_match(tokens, ("change", "printout", "frequency"))
end

function record_printout_frequency_change_request_marker!(
    result::DeckParseResult,
    tokens,
    line_no::Int,
)
    record_control_card!(result, :printout_frequency_change_request, tokens, line_no)
    record_card!(result, :printout_frequency_change_request)
    return record_card!(result, :special_request)
end

function miscellaneous_control_float(values::AbstractVector, index::Int)
    index > length(values) && return 0.0
    value = values[index]
    return value === nothing ? 0.0 : Float64(value)
end

function miscellaneous_control_int(values::AbstractVector, index::Int)
    return _integer_control_value(index > length(values) ? nothing : values[index])
end

function miscellaneous_control_bool(values::AbstractVector, index::Int)
    return miscellaneous_control_int(values, index) != 0
end

function parse_fixed_card_miscellaneous_control!(result::DeckParseResult, tokens,
                                                 line_no::Int)
    record_control_card!(result, :fixed_card_miscellaneous_control, tokens, line_no)
    values = [tryparse_deck_float(String(token)) for token in tokens]
    if get(result.card_counts, :fixed_card_miscellaneous_control, 0) == 1
        push!(
            result.time_horizon_control_rows,
            DeckTimeHorizonControlRow(
                line_no,
                miscellaneous_control_float(values, 1),
                miscellaneous_control_float(values, 2),
                miscellaneous_control_float(values, 3),
                miscellaneous_control_float(values, 4),
                miscellaneous_control_float(values, 5),
                miscellaneous_control_float(values, 6),
                miscellaneous_control_float(values, 7),
                join(token_strings(tokens), " "),
            ),
        )
        record_card!(result, :fixed_card_time_horizon_control)
    else
        push!(
            result.output_schedule_control_rows,
            DeckOutputScheduleControlRow(
                line_no,
                miscellaneous_control_int(values, 1),
                miscellaneous_control_int(values, 2),
                miscellaneous_control_bool(values, 3),
                miscellaneous_control_bool(values, 4),
                miscellaneous_control_bool(values, 5),
                miscellaneous_control_bool(values, 6),
                miscellaneous_control_bool(values, 7),
                miscellaneous_control_int(values, 8),
                miscellaneous_control_int(values, 9),
                miscellaneous_control_int(values, 10),
                join(token_strings(tokens), " "),
            ),
        )
        record_card!(result, :fixed_card_output_network_control)
    end
    return result
end

function fixed_card_saturated_transformer_header_card(tokens)::Bool
    isempty(tokens) && return false
    return normalized_deck_token(tokens[1]) == "transformer"
end

mutable struct SynchronousMachineDataParseState
    machine_index::Int
    source_group_index::Int
    terminal_phase_count::Int
    expected_mass_count::Union{Nothing,Int}
    mass_count::Int
    dynamic_output_count::Int
    phase::Symbol
    input_family::Symbol
end

mutable struct UniversalMachineDataParseState
    expected_machine_count::Union{Nothing,Int}
    machine_index::Int
    definition_card_index::Int
    coil_index::Int
    terminal_branch_count::Int
    machine_type::Union{Nothing,Int}
    expected_definition_card_count::Int
    expected_coil_count::Int
    mechanical_node::Union{Nothing,Symbol}
    mechanical_slack_source::Union{Nothing,Symbol}
    field_slack_source::Union{Nothing,Symbol}
    torque_output_flag::Int
    speed_output_flag::Int
    angle_output_flag::Int
    total_output_count::Int
    tacs_transfer_count::Int
    initialization_mode::Symbol
    input_layout::Symbol
    phase::Symbol
    detailed_machine_state::Union{Nothing,SynchronousMachineDataParseState}
end

UniversalMachineDataParseState() =
    UniversalMachineDataParseState(
        nothing,
        0,
        0,
        0,
        0,
        nothing,
        0,
        0,
        nothing,
        nothing,
        nothing,
        0,
        0,
        0,
        0,
        0,
        :automatic,
        :undetermined,
        :machine_count,
        nothing,
    )

SynchronousMachineDataParseState() =
    SynchronousMachineDataParseState(
        1,
        1,
        0,
        nothing,
        0,
        0,
        :terminal_voltage,
        :standard,
    )

SynchronousMachineDataParseState(
    machine_index::Int,
    input_family::Symbol=:standard,
) =
    SynchronousMachineDataParseState(
        machine_index,
        machine_index,
        0,
        nothing,
        0,
        0,
        :terminal_voltage,
        input_family,
    )

const UNIVERSAL_MACHINE_LEAKAGE_REACTANCE_SCALE = 1.0e3
const UNIVERSAL_MACHINE_SPEED_CAPACITOR_RESISTANCE = 1.0e-8

function universal_machine_data_section_marker(tokens)::Bool
    isempty(tokens) && return false
    source_type = tryparse(Int, deck_token_value(tokens[1]))
    source_type == 19 || return false
    length(tokens) == 1 && return true
    marker = compact_deck_keyword(tokens[2])
    return marker == "um" || marker == "universalmachine"
end

function universal_machine_class_separator(tokens)::Bool
    length(tokens) >= 3 || return false
    compact_deck_keyword(tokens[1]) == "blank" || return false
    compact_deck_keyword(tokens[2]) == "card" || return false
    tail = Set(compact_deck_keyword(token) for token in tokens[3:end])
    return "class" in tail && "um" in tail && "data" in tail
end

function universal_machine_machine_table_separator(tokens)::Bool
    length(tokens) >= 3 || return false
    compact_deck_keyword(tokens[1]) == "blank" || return false
    compact_deck_keyword(tokens[2]) == "card" || return false
    tail = Set(compact_deck_keyword(token) for token in tokens[3:end])
    return "machine" in tail && "table" in tail
end

function universal_machine_data_terminator(tokens)::Bool
    length(tokens) >= 3 || return false
    compact_deck_keyword(tokens[1]) == "blank" || return false
    compact_deck_keyword(tokens[2]) == "card" || return false
    tail = Set(compact_deck_keyword(token) for token in tokens[3:end])
    return "terminating" in tail && "um" in tail && "data" in tail
end

function universal_machine_type_card_counts(machine_type::Int)
    machine_type == 1 && return (definition = 4, coil = 0)
    machine_type == 2 && return (definition = 4, coil = 0)
    machine_type == 3 && return (definition = 4, coil = 5)
    machine_type == 4 && return (definition = 4, coil = 6)
    machine_type == 5 && return (definition = 4, coil = 5)
    machine_type == 6 && return (definition = 4, coil = 4)
    machine_type == 7 && return (definition = 4, coil = 5)
    machine_type == 8 && return (definition = 4, coil = 4)
    machine_type == 9 && return (definition = 4, coil = 5)
    machine_type == 10 && return (definition = 4, coil = 5)
    machine_type == 11 && return (definition = 4, coil = 5)
    machine_type == 12 && return (definition = 4, coil = 5)
    return nothing
end

function universal_machine_fixed_node_symbol(
    image::AbstractString,
    first_col::Int,
    last_col::Int,
)
    value = strip(String(SubString(image, first_col, last_col)))
    return isempty(value) ? nothing : Symbol(value)
end

function universal_machine_node_symbol(tokens)
    for token in tokens
        tryparse_deck_float(String(token)) === nothing || continue
        matched = match(r"[A-Za-z][A-Za-z0-9_]*", String(token))
        matched === nothing && continue
        return Symbol(matched.match)
    end
    return missing
end

function universal_machine_numeric_tokens(tokens)
    values = Float64[]
    for token in tokens
        value = tryparse_deck_float(String(token))
        value === nothing && continue
        push!(values, Float64(value))
    end
    return values
end

function universal_machine_terminal_node_name(machine_index::Int, terminal_index::Int)
    terminal_index in 1:3 || throw(ArgumentError("terminal_index must be 1, 2, or 3"))
    suffix = Char(Int('A') + terminal_index - 1)
    return Symbol("UM$(machine_index)TL$(suffix)")
end

universal_machine_current_node_name(machine_index::Int) = Symbol("UM$(machine_index)MCC")

function universal_machine_fixed_output_flag(image::AbstractString, column::Int)
    value = fixed_int_value(image, column, column)
    return value === nothing ? 0 : value
end

function universal_machine_output_channel_count(flag::Int, multi_channel::Bool)
    flag == 0 && return 0
    if multi_channel
        return flag < 2 ? 1 : (flag == 2 ? 2 : 3)
    end
    return 1
end

function universal_machine_source_index(result::DeckParseResult, source::Union{Nothing,Symbol})
    source === nothing && return 0
    for (index, row) in enumerate(result.over5a_source_rows)
        row.node == source && return index
    end
    return 0
end

function universal_machine_machine_coils(result::DeckParseResult, machine_index::Int)
    return DeckUniversalMachineCoilRow[
        row for row in result.universal_machine_coil_rows if row.machine_index == machine_index
    ]
end

function universal_machine_generated_branch_reactance(coils::AbstractVector{DeckUniversalMachineCoilRow})
    length(coils) >= 3 || return missing
    leakage = min(coils[2].inductance, coils[3].inductance)
    return leakage * UNIVERSAL_MACHINE_LEAKAGE_REACTANCE_SCALE
end

function push_universal_machine_generated_branch_element!(
    result::DeckParseResult,
    row::DeckUniversalMachineGeneratedBranchRow,
)
    ismissing(row.reactance) && return result
    name = string("machine_terminal_reactance_", row.machine_index, "_", row.branch_index)
    inductance = fixed_card_branch_timestep_inductance(result, Float64(row.reactance))
    parse_inductor!(
        result,
        [
            "inductor",
            name,
            String(row.from_node),
            String(row.to_node),
            string(inductance),
        ],
        row.line_no,
    )
    return record_card!(result, :universal_machine_generated_branch_element)
end

function push_universal_machine_speed_capacitor_element!(
    result::DeckParseResult,
    row::DeckUniversalMachineSpeedCapacitorRow,
)
    name = string("machine_speed_capacitor_", row.machine_index)
    if row.resistance > 0.0 && row.capacitance == 0.0
        parse_resistor!(
            result,
            [
                "resistor",
                name,
                String(row.capacitor_node),
                String(row.mass_node),
                string(row.resistance),
            ],
            row.line_no,
        )
    elseif row.resistance == 0.0 && row.capacitance > 0.0
        parse_capacitor!(
            result,
            [
                "capacitor",
                name,
                String(row.capacitor_node),
                String(row.mass_node),
                string(row.capacitance),
            ],
            row.line_no,
        )
    else
        add_issue!(
            result.validation,
            invalid_value(
                "line $(row.line_no)",
                "unsupported universal-machine speed-capacitor resistance/capacitance combination",
            ),
        )
        return result
    end
    return record_card!(result, :universal_machine_speed_capacitor_element)
end

function parse_universal_machine_class1_card!(
    result::DeckParseResult,
    state::UniversalMachineDataParseState,
    line::AbstractString,
    line_no::Int,
)::Bool
    image = fixed_image(line)
    input_mode = something(fixed_int_value(image, 1, 1), 0)
    initialization_flag = fixed_int_value(image, 2, 2)
    detailed_machine_input =
        uppercase(strip(fixed_field(image, 3, 8))) == "SMDATA"
    maximum_shaft_mass_count =
        something(fixed_int_value(image, 9, 14), 10)
    terminal_coupling_flag = something(fixed_int_value(image, 15, 15), 0)
    if input_mode ∉ (0, 1, 2)
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "universal-machine parameter mode must be 0, 1, or 2; got $input_mode",
            ),
        )
        return false
    end
    if initialization_flag === nothing
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "expected universal-machine initialization flag in column 2",
            ),
        )
        return false
    end
    if initialization_flag ∉ (0, 1)
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "universal_machine_initialization_flag must be 0 or 1; got $initialization_flag",
            ),
        )
        return false
    end
    if maximum_shaft_mass_count <= 0
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "universal-machine maximum shaft-mass count must be positive; got $maximum_shaft_mass_count",
            ),
        )
        return false
    end
    if terminal_coupling_flag ∉ (0, 1)
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "universal-machine terminal-coupling flag must be 0 or 1; got $terminal_coupling_flag",
            ),
        )
        return false
    end
    parameter_basis =
        input_mode == 1 ? :power_frequency_normalized : :physical_si
    remanent_flux_enabled = input_mode == 2
    initialization_mode =
        detailed_machine_input || initialization_flag == 1 ? :automatic : :manual
    terminal_coupling =
        terminal_coupling_flag == 0 ? :fully_compensated : :predicted_current
    push!(
        result.universal_machine_section_rows,
        DeckUniversalMachineSectionRow(
            line_no,
            1,
            :undetermined,
            parameter_basis,
            remanent_flux_enabled,
            initialization_mode,
            detailed_machine_input,
            maximum_shaft_mass_count,
            terminal_coupling,
            String(line),
        ),
    )
    # Pure universal-machine cards do not declare a machine count. Successive
    # definition/coil groups continue until the blank section terminator.
    state.expected_machine_count = nothing
    state.machine_index = 1
    state.definition_card_index = 0
    state.coil_index = 0
    state.terminal_branch_count = 0
    state.mechanical_node = nothing
    state.mechanical_slack_source = nothing
    state.field_slack_source = nothing
    state.torque_output_flag = 0
    state.speed_output_flag = 0
    state.angle_output_flag = 0
    state.total_output_count = 0
    state.tacs_transfer_count = 0
    state.initialization_mode = initialization_mode
    state.input_layout = :undetermined
    state.phase = :layout
    record_card!(result, :universal_machine_class1_card)
    return true
end

function set_universal_machine_input_layout!(
    result::DeckParseResult,
    state::UniversalMachineDataParseState,
    input_layout::Symbol,
)
    input_layout in (:modular, :combined_tables) ||
        throw(ArgumentError("unsupported universal-machine input layout $input_layout"))
    section_index = lastindex(result.universal_machine_section_rows)
    section = result.universal_machine_section_rows[section_index]
    result.universal_machine_section_rows[section_index] =
        DeckUniversalMachineSectionRow(
            section.line_no,
            section.machine_count,
            input_layout,
            section.parameter_basis,
            section.remanent_flux_enabled,
            section.initialization_mode,
            section.detailed_machine_input,
            section.maximum_shaft_mass_count,
            section.terminal_coupling,
            section.raw_text,
        )
    state.input_layout = input_layout
    return state
end

function parse_universal_machine_definition_card!(
    result::DeckParseResult,
    state::UniversalMachineDataParseState,
    line::AbstractString,
    tokens,
    line_no::Int,
)::Bool
    state.definition_card_index += 1
    card_index = state.definition_card_index
    node = universal_machine_node_symbol(tokens)
    value1::Union{Missing,Float64} = missing
    value2::Union{Missing,Float64} = missing
    mechanical_damping_coefficient::Union{Missing,Float64} = missing
    speed_convergence_margin::Union{Missing,Float64} = missing
    d_axis_coil_count = 0
    q_axis_coil_count = 0
    torque_output_flag = 0
    speed_output_flag = 0
    angle_output_flag = 0
    pole_pair_count = 0
    rotor_mass::Union{Missing,Float64} = missing
    saturation_mode = 0
    saturated_inductance::Union{Missing,Float64} = missing
    saturation_flux::Union{Missing,Float64} = missing
    remanent_flux::Union{Missing,Float64} = missing
    row_machine_type = 0

    if card_index == 1
        image = fixed_image(line)
        machine_type = fixed_int_value(image, 1, 2)
        if machine_type === nothing
            machine_type = parse_int!(result, tokens[1], line_no, "universal_machine_type")
        end
        machine_type === nothing && return false
        counts = universal_machine_type_card_counts(machine_type)
        if counts === nothing
            record_fixed_blocker!(
                result,
                :universal_machine_data_blocked,
                :universal_machine_type_not_accepted,
            )
            add_issue!(
                result.validation,
                unknown_field(
                    "line $line_no",
                    "Universal-machine type $machine_type is not accepted by the typed deck intake yet",
                ),
            )
            return false
        end
        state.machine_type = machine_type
        state.expected_definition_card_count =
            state.initialization_mode == :manual ? 3 : counts.definition
        state.mechanical_node = node === missing ? nothing : node
        d_axis_coil_count = something(fixed_int_value(image, 3, 4), 0)
        q_axis_coil_count = something(fixed_int_value(image, 5, 6), 0)
        state.expected_coil_count = machine_type in (1, 2) ?
            3 + d_axis_coil_count + q_axis_coil_count : counts.coil
        torque_output_flag = universal_machine_fixed_output_flag(image, 7)
        speed_output_flag = universal_machine_fixed_output_flag(image, 8)
        angle_output_flag = universal_machine_fixed_output_flag(image, 9)
        state.torque_output_flag = torque_output_flag
        state.speed_output_flag = speed_output_flag
        state.angle_output_flag = angle_output_flag
        pole_pair_count = something(fixed_int_value(image, 22, 23), 0)
        rotor_mass = something(fixed_float_value(image, 24, 37), missing)
        mechanical_damping_coefficient = something(fixed_float_value(image, 38, 51), 0.0)
        speed_convergence_margin =
            something(fixed_float_value(image, 52, 65), 0.0)
        speed_convergence_margin >= 0.0 || add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "universal-machine speed convergence margin must be nonnegative",
            ),
        )
        row_machine_type = machine_type
    else
        if card_index in (2, 3)
            image = fixed_image(line)
            fixed_value1 = fixed_float_value(image, 1, 14)
            fixed_value2 = fixed_float_value(image, 15, 28)
            if fixed_value2 === nothing
                values = universal_machine_numeric_tokens(tokens)
                if length(values) == 1
                    value1 = 0.0
                    value2 = values[1]
                elseif length(values) >= 2
                    value1 = values[1]
                    value2 = values[2]
                else
                    value1 = 0.0
                    value2 = 0.0
                end
            else
                value1 = fixed_value1 === nothing ? 0.0 : fixed_value1
                value2 = fixed_value2
            end
            saturation_mode = something(fixed_int_value(image, 29, 29), 0)
            saturated_inductance =
                something(fixed_float_value(image, 30, 43), missing)
            saturation_flux =
                something(fixed_float_value(image, 44, 57), missing)
            remanent_flux =
                something(fixed_float_value(image, 58, 71), missing)
        else
            values = universal_machine_numeric_tokens(tokens)
            if length(values) == 1
                value1 = 0.0
                value2 = values[1]
            elseif length(values) >= 2
                value1 = values[1]
                value2 = values[2]
            else
                value1 = 0.0
                value2 = 0.0
            end
        end
        if card_index == 4
            image = fixed_image(line)
            field_source = universal_machine_fixed_node_symbol(image, 29, 34)
            mechanical_source = universal_machine_fixed_node_symbol(image, 35, 40)
            if state.machine_type in (1, 2, 6, 7, 8, 9, 10, 11, 12)
                state.field_slack_source = field_source
                state.mechanical_slack_source = mechanical_source
                node = mechanical_source === nothing ? missing : mechanical_source
            elseif node !== missing
                state.mechanical_slack_source = node
            end
        end
    end

    push!(
        result.universal_machine_definition_rows,
        DeckUniversalMachineDefinitionRow(
            line_no,
            state.machine_index,
            card_index,
            row_machine_type,
            d_axis_coil_count,
            q_axis_coil_count,
            torque_output_flag,
            speed_output_flag,
            angle_output_flag,
            pole_pair_count,
            rotor_mass,
            saturation_mode,
            saturated_inductance,
            saturation_flux,
            remanent_flux,
            value1,
            value2,
            mechanical_damping_coefficient,
            speed_convergence_margin,
            node,
            String(line),
        ),
    )
    record_card!(result, :universal_machine_definition_card)
    if card_index == state.expected_definition_card_count
        state.phase =
            state.input_layout == :combined_tables ?
            :combined_machine_boundary : :coil
        state.coil_index = 0
        state.terminal_branch_count = 0
    end
    return true
end

function universal_machine_output_flag(tokens)
    isempty(tokens) && return 0
    flag = tryparse(Int, deck_token_value(tokens[end]))
    return flag === nothing ? 0 : flag
end

function universal_machine_coil_numeric_values(tokens)
    values = universal_machine_numeric_tokens(tokens)
    if !isempty(tokens) && tryparse(Int, deck_token_value(tokens[end])) !== nothing &&
       !isempty(values)
        pop!(values)
    end
    return values
end

function push_universal_machine_terminal_rows!(
    result::DeckParseResult,
    state::UniversalMachineDataParseState,
    terminal_node,
    line::AbstractString,
    line_no::Int,
)
    generated_leakage_branch = state.initialization_mode == :automatic
    if generated_leakage_branch && terminal_node !== missing && state.terminal_branch_count < 3
        state.terminal_branch_count += 1
        generated_node =
            universal_machine_terminal_node_name(state.machine_index, state.terminal_branch_count)
        generated_node_value = node_id!(result, String(generated_node))
        connected_node_value = node_id!(result, String(terminal_node))
        push!(
            result.universal_machine_terminal_rows,
            DeckUniversalMachineTerminalRow(
                line_no,
                state.machine_index,
                state.coil_index,
                generated_node,
                :ground,
                generated_node_value,
                0,
                String(line),
            ),
        )
        push!(
            result.universal_machine_generated_branch_rows,
            DeckUniversalMachineGeneratedBranchRow(
                line_no,
                state.machine_index,
                state.terminal_branch_count,
                generated_node,
                terminal_node,
                generated_node_value,
                connected_node_value,
                missing,
                String(line),
            ),
        )
        record_card!(result, :universal_machine_terminal_node)
        return record_card!(result, :universal_machine_generated_branch)
    end

    connected_terminal = terminal_node === missing ? :ground : terminal_node
    connected_terminal_value = terminal_node === missing ? 0 :
        node_id!(result, String(terminal_node))
    push!(
        result.universal_machine_terminal_rows,
        DeckUniversalMachineTerminalRow(
            line_no,
            state.machine_index,
            state.coil_index,
            connected_terminal,
            :ground,
            connected_terminal_value,
            0,
            String(line),
        ),
    )
    return record_card!(result, :universal_machine_terminal_node)
end

function finalize_universal_machine_machine_rows!(
    result::DeckParseResult,
    state::UniversalMachineDataParseState,
    line::AbstractString,
    line_no::Int,
)
    coils = universal_machine_machine_coils(result, state.machine_index)
    if state.field_slack_source === nothing &&
       state.machine_type in (1, 2, 6, 7, 8) &&
       length(coils) >= 4
        field_terminal = coils[4].terminal_node
        ismissing(field_terminal) || (state.field_slack_source = field_terminal)
    end
    if state.mechanical_slack_source === nothing && state.mechanical_node !== nothing
        source_nodes = Set(row.node for row in result.over5a_source_rows)
        for row in result.over2_branch_rows
            if row.from_node == state.mechanical_node && row.to_node in source_nodes
                state.mechanical_slack_source = row.to_node
                break
            elseif row.to_node == state.mechanical_node && row.from_node in source_nodes
                state.mechanical_slack_source = row.from_node
                break
            end
        end
    end
    branch_reactance = universal_machine_generated_branch_reactance(coils)
    for (index, row) in pairs(result.universal_machine_generated_branch_rows)
        row.machine_index == state.machine_index || continue
        result.universal_machine_generated_branch_rows[index] =
            DeckUniversalMachineGeneratedBranchRow(
                row.line_no,
                row.machine_index,
                row.branch_index,
                row.from_node,
                row.to_node,
                row.from_node_value,
                row.to_node_value,
                branch_reactance,
                row.raw_text,
            )
        push_universal_machine_generated_branch_element!(
            result,
            result.universal_machine_generated_branch_rows[index],
        )
    end

    mass_node = state.mechanical_node
    if mass_node !== nothing
        mass_node_value = node_id!(result, String(mass_node))
        if state.initialization_mode == :automatic
            capacitor_node = universal_machine_current_node_name(state.machine_index)
            capacitor_node_value = node_id!(result, String(capacitor_node))
            push!(
                result.universal_machine_speed_capacitor_rows,
                DeckUniversalMachineSpeedCapacitorRow(
                    line_no,
                    state.machine_index,
                    capacitor_node,
                    mass_node,
                    capacitor_node_value,
                    mass_node_value,
                    UNIVERSAL_MACHINE_SPEED_CAPACITOR_RESISTANCE,
                    0.0,
                    String(line),
                ),
            )
            record_card!(result, :universal_machine_speed_capacitor)
            push_universal_machine_speed_capacitor_element!(
                result,
                result.universal_machine_speed_capacitor_rows[end],
            )
        end

        mechanical_slack_source =
            state.mechanical_slack_source === nothing ? missing : state.mechanical_slack_source
        field_slack_source =
            state.field_slack_source === nothing ? missing : state.field_slack_source
        push!(
            result.universal_machine_node_summary_rows,
            DeckUniversalMachineNodeSummaryRow(
                line_no,
                state.machine_index,
                mass_node,
                mass_node_value,
                mechanical_slack_source,
                universal_machine_source_index(result, state.mechanical_slack_source),
                field_slack_source,
                universal_machine_source_index(result, state.field_slack_source),
                String(line),
            ),
        )
        record_card!(result, :universal_machine_node_summary)
    end

    machine_output_count =
        universal_machine_output_channel_count(state.torque_output_flag, true) +
        universal_machine_output_channel_count(state.speed_output_flag, true) +
        universal_machine_output_channel_count(state.angle_output_flag, false) +
        count(row -> row.output_flag != 0, coils)
    state.total_output_count += machine_output_count
    machine_count = state.expected_machine_count
    if machine_count !== nothing && state.machine_index == machine_count
        push!(
            result.universal_machine_output_summary_rows,
            DeckUniversalMachineOutputSummaryRow(
                line_no,
                machine_count,
                state.total_output_count,
                state.tacs_transfer_count,
                String(line),
            ),
        )
        record_card!(result, :universal_machine_output_summary)
    end
    return state
end

function begin_next_universal_machine!(state::UniversalMachineDataParseState)
    state.machine_index += 1
    state.definition_card_index = 0
    state.coil_index = 0
    state.terminal_branch_count = 0
    state.machine_type = nothing
    state.expected_definition_card_count = 0
    state.expected_coil_count = 0
    state.mechanical_node = nothing
    state.mechanical_slack_source = nothing
    state.field_slack_source = nothing
    state.torque_output_flag = 0
    state.speed_output_flag = 0
    state.angle_output_flag = 0
    state.phase = :definition
    return state
end

function prepare_combined_table_coils!(
    result::DeckParseResult,
    state::UniversalMachineDataParseState,
    machine_index::Int,
)
    definitions = [
        row for row in result.universal_machine_definition_rows
        if row.machine_index == machine_index
    ]
    card1_index = findfirst(row -> row.card_index == 1, definitions)
    card1_index === nothing &&
        throw(ArgumentError("combined universal-machine table is missing machine $machine_index card 1"))
    card1 = definitions[card1_index]
    counts = universal_machine_type_card_counts(card1.machine_type)
    counts === nothing &&
        throw(ArgumentError("combined universal-machine table has unsupported type $(card1.machine_type)"))
    state.machine_index = machine_index
    state.definition_card_index =
        state.initialization_mode == :manual ? 3 : counts.definition
    state.coil_index = 0
    state.terminal_branch_count = 0
    state.machine_type = card1.machine_type
    state.expected_definition_card_count = state.definition_card_index
    state.expected_coil_count = card1.machine_type in (1, 2) ?
        3 + card1.d_axis_coil_count + card1.q_axis_coil_count : counts.coil
    state.mechanical_node = ismissing(card1.node) ? nothing : card1.node
    state.mechanical_slack_source = nothing
    state.field_slack_source = nothing
    state.torque_output_flag = card1.torque_output_flag
    state.speed_output_flag = card1.speed_output_flag
    state.angle_output_flag = card1.angle_output_flag
    card4_index = findfirst(row -> row.card_index == 4, definitions)
    if card4_index !== nothing
        card4 = definitions[card4_index]
        if card1.machine_type in (1, 2, 6, 7, 8)
            state.mechanical_slack_source =
                ismissing(card4.node) ? nothing : card4.node
            image = fixed_image(card4.raw_text)
            state.field_slack_source =
                universal_machine_fixed_node_symbol(image, 29, 34)
        elseif !ismissing(card4.node)
            state.mechanical_slack_source = card4.node
        end
    end
    state.phase = :coil
    return state
end

function advance_universal_machine_state!(result::DeckParseResult,
                                          state::UniversalMachineDataParseState,
                                          line::AbstractString,
                                          line_no::Int)
    if state.coil_index != state.expected_coil_count
        return state
    end
    finalize_universal_machine_machine_rows!(result, state, line, line_no)
    if state.expected_machine_count !== nothing &&
       state.machine_index < state.expected_machine_count
        if state.input_layout == :combined_tables
            prepare_combined_table_coils!(
                result,
                state,
                state.machine_index + 1,
            )
        else
            begin_next_universal_machine!(state)
        end
    else
        state.phase = :await_terminator
    end
    return state
end

function parse_universal_machine_coil_card!(
    result::DeckParseResult,
    state::UniversalMachineDataParseState,
    line::AbstractString,
    tokens,
    line_no::Int,
)::Bool
    state.coil_index += 1
    image = fixed_image(line)
    terminal_node = universal_machine_fixed_node_symbol(image, 29, 34)
    terminal_node = terminal_node === nothing ? missing : terminal_node
    control_signal = universal_machine_fixed_node_symbol(image, 41, 46)
    control_signal = control_signal === nothing ? missing : control_signal
    resistance = something(fixed_float_value(image, 1, 14), 0.0)
    inductance = something(fixed_float_value(image, 15, 28), 0.0)
    output_flag = something(fixed_int_value(image, 47, 47), 0)
    initial_history_current = something(fixed_float_value(image, 48, 61), 0.0)
    push!(
        result.universal_machine_coil_rows,
        DeckUniversalMachineCoilRow(
            line_no,
            state.machine_index,
            state.coil_index,
            resistance,
            inductance,
            terminal_node,
            control_signal,
            output_flag,
            initial_history_current,
            String(line),
        ),
    )
    record_card!(result, :universal_machine_coil_card)
    push_universal_machine_terminal_rows!(result, state, terminal_node, line, line_no)
    advance_universal_machine_state!(result, state, line, line_no)
    return true
end

function parse_universal_machine_data_card!(
    result::DeckParseResult,
    state::UniversalMachineDataParseState,
    line::AbstractString,
    tokens,
    line_no::Int,
)::Bool
    if universal_machine_data_terminator(tokens)
        accepted_terminator_phase =
            state.phase == :await_terminator ||
            state.phase == :await_detailed_machine_terminator
        if !accepted_terminator_phase
            add_issue!(
                result.validation,
                missing_data(
                    "line $line_no",
                    "universal-machine data ended before all accepted type-4 rows were read",
                ),
            )
        end
        if accepted_terminator_phase && state.expected_machine_count === nothing
            section_index = findlast(result.universal_machine_section_rows) do row
                row.initialization_mode == state.initialization_mode
            end
            section_index === nothing && state.phase == :await_detailed_machine_terminator &&
                (section_index = lastindex(result.universal_machine_section_rows))
            if section_index !== nothing
                section = result.universal_machine_section_rows[section_index]
                result.universal_machine_section_rows[section_index] =
                    DeckUniversalMachineSectionRow(
                        section.line_no,
                        state.machine_index,
                        section.input_layout,
                        section.parameter_basis,
                        section.remanent_flux_enabled,
                        section.initialization_mode,
                        section.detailed_machine_input,
                        section.maximum_shaft_mass_count,
                        section.terminal_coupling,
                        section.raw_text,
                    )
            end
            push!(
                result.universal_machine_output_summary_rows,
                DeckUniversalMachineOutputSummaryRow(
                    line_no,
                    state.machine_index,
                    state.total_output_count,
                    state.tacs_transfer_count,
                    String(line),
                ),
            )
            record_card!(result, :universal_machine_output_summary)
        end
        record_control_card!(result, :universal_machine_data_section_end, tokens, line_no)
        record_card!(result, :universal_machine_data_section_end)
        return true
    end
    if state.phase == :detailed_synchronous_machine
        detailed = state.detailed_machine_state
        detailed === nothing &&
            error("detailed synchronous-machine parser state is missing")
        previous_phase = detailed.phase
        finished = parse_synchronous_machine_data_card!(
            result,
            detailed,
            line,
            tokens,
            line_no,
        )
        if previous_phase == :control && detailed.phase == :model_option
            push!(
                result.synchronous_machine_model_parameter_rows,
                DeckSynchronousMachineModelParameterRow(
                    line_no,
                    detailed.machine_index,
                    :model_option,
                    zeros(7),
                    Union{Missing,Float64}[0.0 for _ in 1:7],
                    "implicit zero model-option row for detailed machine in universal section",
                ),
            )
            detailed.phase = :electrical_reactance
            record_card!(result, :synchronous_machine_implicit_model_option)
        end
        if finished
            state.tacs_transfer_count = count(
                row -> row.machine_index == detailed.machine_index,
                result.synchronous_machine_control_interface_rows,
            )
            state.total_output_count = detailed.dynamic_output_count
            state.phase = :await_detailed_machine_terminator
        end
        return false
    elseif state.phase == :await_detailed_machine_terminator
        add_issue!(
            result.validation,
            unknown_field(
                "line $line_no",
                "expected the universal-machine data terminator after detailed synchronous-machine FINISH",
            ),
        )
        return false
    end
    if state.phase == :machine_count
        parse_universal_machine_class1_card!(result, state, line, line_no)
        return false
    elseif state.phase == :layout
        if universal_machine_class_separator(tokens)
            set_universal_machine_input_layout!(result, state, :modular)
            state.phase = :definition
            record_control_card!(
                result,
                :universal_machine_class_section_end,
                tokens,
                line_no,
            )
            return false
        end
        set_universal_machine_input_layout!(result, state, :combined_tables)
        state.phase = :definition
        parse_universal_machine_definition_card!(result, state, line, tokens, line_no)
        return false
    elseif state.phase == :combined_machine_boundary
        if universal_machine_machine_table_separator(tokens)
            state.expected_machine_count = state.machine_index
            prepare_combined_table_coils!(result, state, 1)
            record_control_card!(
                result,
                :universal_machine_machine_table_end,
                tokens,
                line_no,
            )
            return false
        end
        begin_next_universal_machine!(state)
        state.phase = :definition
        parse_universal_machine_definition_card!(result, state, line, tokens, line_no)
        return false
    elseif state.phase == :definition
        source_type = fixed_int_value(fixed_image(line), 1, 2)
        if source_type !== nothing && source_type in (50, 52, 59)
            detailed = SynchronousMachineDataParseState(
                length(result.synchronous_machine_output_summary_rows) + 1,
                :universal_machine_section,
            )
            state.detailed_machine_state = detailed
            state.phase = :detailed_synchronous_machine
            state.initialization_mode = :automatic
            section_index = lastindex(result.universal_machine_section_rows)
            section = result.universal_machine_section_rows[section_index]
            result.universal_machine_section_rows[section_index] =
                DeckUniversalMachineSectionRow(
                    section.line_no,
                    1,
                    section.input_layout,
                    section.parameter_basis,
                    section.remanent_flux_enabled,
                    :automatic,
                    true,
                    section.maximum_shaft_mass_count,
                    section.terminal_coupling,
                    section.raw_text,
                )
            record_card!(result, :synchronous_machine_universal_section_input)
            parse_synchronous_machine_data_card!(
                result,
                detailed,
                line,
                tokens,
                line_no,
            )
            return false
        end
        parse_universal_machine_definition_card!(result, state, line, tokens, line_no)
        return false
    elseif state.phase == :coil
        parse_universal_machine_coil_card!(result, state, line, tokens, line_no)
        return false
    elseif state.phase == :await_terminator && state.expected_machine_count === nothing
        begin_next_universal_machine!(state)
        parse_universal_machine_definition_card!(result, state, line, tokens, line_no)
        return false
    end
    add_issue!(
        result.validation,
        unknown_field(
            "line $line_no",
            "Unexpected universal-machine data row after accepted machine data",
        ),
    )
    return false
end

function synchronous_machine_terminal_voltage_card(line::AbstractString)::Bool
    source_type = fixed_int_value(fixed_image(line), 1, 2)
    return source_type !== nothing && source_type in 50:59
end

function synchronous_machine_numeric_values!(
    result::DeckParseResult,
    tokens,
    line_no::Int,
    field_prefix::AbstractString,
)
    values = Float64[]
    for (index, token) in enumerate(tokens)
        value = parse_float!(result, token, line_no, "$(field_prefix)_$index")
        value === nothing && return nothing
        push!(values, Float64(value))
    end
    return values
end

function synchronous_machine_int_values!(
    result::DeckParseResult,
    tokens,
    line_no::Int,
    field_prefix::AbstractString,
)
    values = Int[]
    for (index, token) in enumerate(tokens)
        value = parse_int!(result, token, line_no, "$(field_prefix)_$index")
        value === nothing && return nothing
        push!(values, Int(value))
    end
    return values
end

function parse_synchronous_machine_terminal_voltage_card!(
    result::DeckParseResult,
    state::SynchronousMachineDataParseState,
    line::AbstractString,
    line_no::Int,
)::Bool
    image = fixed_image(line)
    source_type = fixed_int_field!(result, image, line_no, 1, 2, "source_type")
    source_type in 50:59 || return false
    node = fixed_field(image, 3, 8)
    if isempty(node)
        add_issue!(
            result.validation,
            missing_data("line $line_no", "expected type-50 through type-59 synchronous-machine terminal node in columns 3-8"),
        )
        return false
    end
    terminal_node = Symbol(node)
    terminal_node_value = node_id!(result, node)
    phase_index = state.terminal_phase_count + 1
    push!(
        result.synchronous_machine_terminal_voltage_rows,
        DeckSynchronousMachineTerminalVoltageRow(
            line_no,
            state.machine_index,
            state.source_group_index,
            source_type,
            phase_index,
            terminal_node,
            terminal_node_value,
            something(fixed_float_value(image, 11, 20), missing),
            something(fixed_float_value(image, 21, 30), missing),
            something(fixed_float_value(image, 31, 40), missing),
            String(line),
        ),
    )
    state.terminal_phase_count = phase_index
    state.phase = phase_index >= 3 ? :control : :terminal_voltage
    record_card!(result, :fixed_field)
    record_card!(result, :bpa_fixed_source)
    record_card!(result, Symbol("synchronous_machine_source_type_", source_type))
    source_type == 59 && record_card!(result, :bpa_fixed_source_type59)
    record_card!(result, :synchronous_machine_terminal_voltage)
    record_card!(result, :synchronous_machine_data_section)
    return true
end

function parse_synchronous_machine_tolerances!(
    result::DeckParseResult,
    state::SynchronousMachineDataParseState,
    tokens,
    line_no::Int,
)::Bool
    values = synchronous_machine_numeric_values!(
        result,
        tokens[2:end],
        line_no,
        "synchronous_machine_tolerance",
    )
    values === nothing && return true
    push!(
        result.synchronous_machine_tolerance_rows,
        DeckSynchronousMachineToleranceRow(
            line_no,
            state.machine_index,
            values,
            join(token_strings(tokens), " "),
        ),
    )
    record_card!(result, :fixed_field)
    record_card!(result, :synchronous_machine_tolerance)
    return true
end

function parse_synchronous_machine_parameter_fitting!(
    result::DeckParseResult,
    state::SynchronousMachineDataParseState,
    tokens,
    line_no::Int,
)::Bool
    value::Union{Missing,Float64} = missing
    if length(tokens) >= 3
        parsed = parse_float!(
            result,
            tokens[3],
            line_no,
            "synchronous_machine_parameter_fitting_value",
        )
        parsed === nothing && return true
        value = Float64(parsed)
    end
    push!(
        result.synchronous_machine_parameter_fitting_rows,
        DeckSynchronousMachineParameterFittingRow(
            line_no,
            state.machine_index,
            value,
            join(token_strings(tokens), " "),
        ),
    )
    record_card!(result, :fixed_field)
    record_card!(result, :synchronous_machine_parameter_fitting)
    return true
end

function parse_synchronous_machine_delta_connection!(
    result::DeckParseResult,
    state::SynchronousMachineDataParseState,
    line::AbstractString,
    tokens,
    line_no::Int,
)
    push!(
        result.synchronous_machine_model_parameter_rows,
        DeckSynchronousMachineModelParameterRow(
            line_no,
            state.machine_index,
            :delta_connection,
            Float64[],
            Union{Missing,Float64}[],
            String(line),
        ),
    )
    record_control_card!(result, :synchronous_machine_delta_connection, tokens, line_no)
    record_card!(result, :fixed_field)
    record_card!(result, :synchronous_machine_delta_connection)
    return true
end

function synchronous_machine_model_positional_values(
    line::AbstractString,
    parameter_kind::Symbol,
)
    image = fixed_image(line)
    ranges =
        parameter_kind == :definition ?
        ((1, 2), (3, 4), (5, 6), (7, 10), (11, 20), (21, 30),
         (31, 40), (41, 50), (51, 60), (61, 70), (71, 80)) :
        ntuple(index -> (10 * index - 9, 10 * index),
               parameter_kind == :electrical_reactance ? 8 : 7)
    return Union{Missing,Float64}[
        something(fixed_float_value(image, first, last), missing)
        for (first, last) in ranges
    ]
end

function parse_synchronous_machine_model_parameter!(
    result::DeckParseResult,
    state::SynchronousMachineDataParseState,
    line::AbstractString,
    tokens,
    line_no::Int,
    parameter_kind::Symbol,
)::Bool
    values = synchronous_machine_numeric_values!(
        result,
        tokens,
        line_no,
        "synchronous_machine_$(String(parameter_kind))",
    )
    values === nothing && return true
    positional_values = synchronous_machine_model_positional_values(line, parameter_kind)
    push!(
        result.synchronous_machine_model_parameter_rows,
        DeckSynchronousMachineModelParameterRow(
            line_no,
            state.machine_index,
            parameter_kind,
            values,
            positional_values,
            String(line),
        ),
    )
    if parameter_kind == :definition && !isempty(values)
        expected_mass_count = round(Int, values[1])
        if expected_mass_count <= 0
            add_issue!(
                result.validation,
                invalid_value(
                    "line $line_no",
                    "type-59 synchronous-machine mass count must be positive",
                ),
            )
            return true
        end
        state.expected_mass_count = expected_mass_count
        state.phase = :model_option
    elseif parameter_kind == :model_option
        state.phase = :electrical_reactance
    elseif parameter_kind == :electrical_reactance
        state.phase = :electrical_time_constant
    elseif parameter_kind == :electrical_time_constant
        state.phase = :mass
    end
    record_card!(result, :fixed_field)
    record_card!(result, :synchronous_machine_model_parameter)
    record_card!(result, Symbol("synchronous_machine_", String(parameter_kind)))
    return true
end

function parse_synchronous_machine_mass!(
    result::DeckParseResult,
    state::SynchronousMachineDataParseState,
    line::AbstractString,
    tokens,
    line_no::Int,
)::Bool
    values = synchronous_machine_numeric_values!(
        result,
        tokens,
        line_no,
        "synchronous_machine_mass",
    )
    values === nothing && return true
    mass_index = round(Int, values[1])
    if mass_index <= 0
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "type-59 synchronous-machine mass index must be positive"),
        )
        return true
    end
    image = fixed_image(line)
    mass_values = Float64[
        something(fixed_float_value(image, 10 * index + 1, 10 * index + 10), 0.0)
        for index in 1:6
    ]
    push!(
        result.synchronous_machine_mass_rows,
        DeckSynchronousMachineMassRow(
            line_no,
            state.machine_index,
            mass_index,
            values,
            mass_values...,
            String(line),
        ),
    )
    state.mass_count += 1
    if state.expected_mass_count !== nothing &&
       state.mass_count >= state.expected_mass_count
        state.phase = :output
    end
    record_card!(result, :fixed_field)
    record_card!(result, :synchronous_machine_mass)
    return true
end

function parse_synchronous_machine_output_request!(
    result::DeckParseResult,
    state::SynchronousMachineDataParseState,
    line::AbstractString,
    tokens,
    line_no::Int,
)::Bool
    values = synchronous_machine_int_values!(
        result,
        tokens,
        line_no,
        "synchronous_machine_output",
    )
    values === nothing && return true
    isempty(values) && return true
    group_index = values[1]
    output_codes = values[2:end]
    dynamic_output_count = group_index in 1:4 ? length(output_codes) : 0
    push!(
        result.synchronous_machine_output_request_rows,
        DeckSynchronousMachineOutputRequestRow(
            line_no,
            state.machine_index,
            group_index,
            output_codes,
            dynamic_output_count,
            String(line),
        ),
    )
    state.dynamic_output_count += dynamic_output_count
    record_card!(result, :fixed_field)
    record_card!(result, :synchronous_machine_output_request)
    dynamic_output_count > 0 ?
        record_card!(result, :synchronous_machine_dynamic_output_request) :
        record_card!(result, :synchronous_machine_steady_state_output_request)
    return true
end

function synchronous_machine_control_interface_code(line::AbstractString)
    code = fixed_int_value(fixed_image(line), 1, 2)
    return code !== nothing && code in 71:74 ? code : nothing
end

function parse_synchronous_machine_control_interface!(
    result::DeckParseResult,
    state::SynchronousMachineDataParseState,
    line::AbstractString,
    line_no::Int,
    interface_code::Int,
)
    image = fixed_image(line)
    signal_text = fixed_field(image, 3, 8)
    if isempty(signal_text)
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "synchronous-machine control interface requires a signal name in columns 3-8",
            ),
        )
        return true
    end
    variable_index = something(fixed_int_value(image, 15, 17), 0)
    alternative_input = state.input_family == :universal_machine_section
    if !alternative_input && interface_code != 71 && variable_index <= 0
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "synchronous-machine control interface code $interface_code requires a positive index in columns 15-17",
            ),
        )
        return true
    end
    direction = interface_code <= 72 ? :control_to_machine : :machine_to_control
    coupling_kind = if alternative_input && interface_code == 71
        :external_field_voltage_input
    elseif interface_code == 71
        :field_voltage_multiplier
    elseif alternative_input && interface_code == 72
        :total_applied_torque_input
    elseif alternative_input && interface_code == 73
        :exciter_voltage_output
    elseif alternative_input && interface_code == 74
        :exciter_current_output
    elseif interface_code == 72
        :mechanical_torque_input
    elseif interface_code == 73
        :machine_state_output
    else
        :rotor_mass_output
    end
    push!(
        result.synchronous_machine_control_interface_rows,
        DeckSynchronousMachineControlInterfaceRow(
            line_no,
            state.machine_index,
            interface_code,
            direction,
            coupling_kind,
            Symbol(signal_text),
            variable_index,
            String(line),
        ),
    )
    record_card!(result, :fixed_field)
    record_card!(result, :synchronous_machine_control_interface)
    record_card!(result, Symbol("synchronous_machine_control_interface_", interface_code))
    return true
end

function synchronous_machine_blank_card(tokens)::Bool
    length(tokens) >= 3 || return false
    compact_deck_keyword(tokens[1]) == "blank" || return false
    compact_deck_keyword(tokens[2]) == "card" || return false
    return true
end

function finish_synchronous_machine_data!(
    result::DeckParseResult,
    state::SynchronousMachineDataParseState,
    tokens,
    line_no::Int,
)::Bool
    parallel_part = length(tokens) >= 2 &&
                    compact_deck_keyword(tokens[2]) == "part"
    if state.expected_mass_count !== nothing && state.mass_count != state.expected_mass_count
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "type-59 synchronous-machine mass count $(state.mass_count) does not match expected $(state.expected_mass_count)",
            ),
        )
    end
    push!(
        result.synchronous_machine_output_summary_rows,
        DeckSynchronousMachineOutputSummaryRow(
            line_no,
            state.machine_index,
            state.terminal_phase_count,
            state.mass_count,
            state.dynamic_output_count,
            join(token_strings(tokens), " "),
        ),
    )
    control_kind = parallel_part ?
        :synchronous_machine_parallel_part :
        :synchronous_machine_data_section_end
    record_control_card!(result, control_kind, tokens, line_no)
    record_card!(result, :fixed_field)
    record_card!(result, :synchronous_machine_output_summary)
    parallel_part || return true

    source_machine_index = state.machine_index
    state.machine_index += 1
    source_rows = [
        row for row in result.synchronous_machine_terminal_voltage_rows
        if row.machine_index == source_machine_index
    ]
    length(source_rows) == 3 || add_issue!(
        result.validation,
        invalid_value(
            "line $line_no",
            "parallel synchronous-machine part requires three inherited terminal phases",
        ),
    )
    for row in source_rows
        push!(
            result.synchronous_machine_terminal_voltage_rows,
            DeckSynchronousMachineTerminalVoltageRow(
                row.line_no,
                state.machine_index,
                row.source_group_index,
                row.source_type,
                row.phase_index,
                row.terminal_node,
                row.terminal_node_value,
                row.peak_terminal_voltage,
                row.frequency_hz,
                row.angle_deg,
                row.raw_text,
            ),
        )
    end
    state.expected_mass_count = nothing
    state.mass_count = 0
    state.dynamic_output_count = 0
    state.phase = :control
    return false
end

function parse_synchronous_machine_data_card!(
    result::DeckParseResult,
    state::SynchronousMachineDataParseState,
    line::AbstractString,
    tokens,
    line_no::Int,
)::Bool
    isempty(tokens) && return false
    first_token = compact_deck_keyword(tokens[1])
    interface_code = synchronous_machine_control_interface_code(line)
    if synchronous_machine_terminal_voltage_card(line)
        parse_synchronous_machine_terminal_voltage_card!(result, state, line, line_no)
        return false
    elseif interface_code !== nothing
        parse_synchronous_machine_control_interface!(
            result,
            state,
            line,
            line_no,
            interface_code,
        )
        return false
    elseif first_token == "tolerances"
        parse_synchronous_machine_tolerances!(result, state, tokens, line_no)
        return false
    elseif deck_phrase_match(tokens, ("parameter", "fitting"))
        parse_synchronous_machine_parameter_fitting!(result, state, tokens, line_no)
        return false
    elseif deck_phrase_match(tokens, ("delta", "connection"))
        parse_synchronous_machine_delta_connection!(result, state, line, tokens, line_no)
        return false
    elseif first_token == "finish"
        return finish_synchronous_machine_data!(result, state, tokens, line_no)
    elseif synchronous_machine_blank_card(tokens)
        record_control_card!(result, :synchronous_machine_data_separator, tokens, line_no)
        record_card!(result, :fixed_field)
        record_card!(result, :synchronous_machine_data_separator)
        state.phase == :mass && (state.phase = :output)
        return false
    elseif state.phase == :control
        parse_synchronous_machine_model_parameter!(
            result,
            state,
            line,
            tokens,
            line_no,
            :definition,
        )
        return false
    elseif state.phase == :model_option
        parse_synchronous_machine_model_parameter!(
            result,
            state,
            line,
            tokens,
            line_no,
            :model_option,
        )
        return false
    elseif state.phase == :electrical_reactance
        parse_synchronous_machine_model_parameter!(
            result,
            state,
            line,
            tokens,
            line_no,
            :electrical_reactance,
        )
        return false
    elseif state.phase == :electrical_time_constant
        parse_synchronous_machine_model_parameter!(
            result,
            state,
            line,
            tokens,
            line_no,
            :electrical_time_constant,
        )
        return false
    elseif state.phase == :mass
        parse_synchronous_machine_mass!(result, state, line, tokens, line_no)
        return false
    elseif state.phase == :output
        parse_synchronous_machine_output_request!(result, state, line, tokens, line_no)
        return false
    end
    add_issue!(
        result.validation,
        unknown_field(
            "line $line_no",
            "Unexpected type-59 synchronous-machine data row in $(state.phase) section",
        ),
    )
    return false
end

function parse_fixed_field_section_card!(result::DeckParseResult, active_section,
                                         line::AbstractString, tokens,
                                         line_no::Int;
                                         branch_vintage_mode::Int = 0)::Bool
    active_section === nothing && return false
    isempty(tokens) && return false
    first_token = normalized_deck_token(tokens[1])
    if first_token in EXPLICIT_DECK_CARD_TOKENS
        if active_section == :blank_output
            first_token in ("bus", "node") || return false
        elseif isempty(String(line)) || !isspace(first(String(line)))
            return false
        end
    end
    if active_section == :blank_branch
        return parse_bpa_fixed_branch_card!(
            result,
            line,
            line_no;
            branch_vintage_mode = branch_vintage_mode,
        )
    elseif active_section == :blank_source
        return parse_bpa_fixed_source_card!(result, line, line_no)
    elseif active_section == :blank_switch
        return parse_bpa_fixed_switch_card!(result, line, line_no)
    elseif active_section == :blank_initial_condition_or_output
        if node_initial_condition_row(line)
            return parse_node_initial_condition_row!(result, line, line_no)
        end
        return parse_bpa_fixed_output_card!(result, line, line_no)
    elseif active_section == :blank_output
        return parse_bpa_fixed_output_card!(result, line, line_no)
    elseif active_section == :line_constants_conductor
        return parse_line_constants_conductor_card!(result, line, tokens, line_no)
    elseif active_section == :line_constants_frequency
        return parse_line_constants_frequency_card!(result, line, line_no)
    elseif active_section == :control_system_hybrid
        return parse_control_system_hybrid_card!(result, line, tokens, line_no)
    end
    if haskey(BPA_FIXED_UNSUPPORTED_SECTION_BLOCKERS, active_section)
        blocker, message = BPA_FIXED_UNSUPPORTED_SECTION_BLOCKERS[active_section]
        record_fixed_blocker!(result, :bpa_fixed_unsupported_section, blocker)
        add_issue!(result.validation, unknown_field("line $line_no", message))
        return true
    end
    return false
end

fixed_image(line::AbstractString) = rpad(String(line), 80)

function fixed_field(image::AbstractString, first_col::Int, last_col::Int)::String
    first_col <= last_col || return ""
    first_col > lastindex(image) && return ""
    stop_col = min(last_col, lastindex(image))
    return strip(image[first_col:stop_col])
end

function normalized_fortran_float_text(raw::AbstractString)
    text = replace(strip(String(raw)), 'D' => 'E', 'd' => 'e')
    isempty(text) && return text
    occursin('E', uppercase(text)) && return text
    sign_index = nothing
    current = nextind(text, firstindex(text), 1)
    while current <= lastindex(text)
        char = text[current]
        if char == '+' || char == '-'
            sign_index = current
        end
        current = nextind(text, current)
    end
    sign_index === nothing && return text
    mantissa = text[firstindex(text):prevind(text, sign_index)]
    exponent = text[sign_index:lastindex(text)]
    any(isdigit, mantissa) && any(isdigit, exponent) || return text
    return string(mantissa, "E", exponent)
end

function tryparse_deck_float(raw::AbstractString)
    text = strip(String(raw))
    isempty(text) && return nothing
    value = tryparse(Float64, text)
    value === nothing || return value
    normalized = normalized_fortran_float_text(text)
    normalized == text && return nothing
    return tryparse(Float64, normalized)
end

function fixed_float_field!(result::DeckParseResult, image::AbstractString, line_no::Int,
                            first_col::Int, last_col::Int, field::AbstractString)
    raw = fixed_field(image, first_col, last_col)
    if isempty(raw)
        add_issue!(result.validation,
                   missing_data("line $line_no", "expected fixed-field $field in columns $first_col-$last_col"))
        return nothing
    end
    value = tryparse_deck_float(raw)
    if value === nothing
        add_issue!(result.validation,
                   invalid_value("line $line_no",
                                 "$field=$(raw): expected fixed-field Float64 in columns $first_col-$last_col"))
    end
    return value
end

function fixed_int_field!(result::DeckParseResult, image::AbstractString, line_no::Int,
                          first_col::Int, last_col::Int, field::AbstractString)
    raw = fixed_field(image, first_col, last_col)
    if isempty(raw)
        add_issue!(result.validation,
                   missing_data("line $line_no", "expected fixed-field $field in columns $first_col-$last_col"))
        return nothing
    end
    value = tryparse(Int, raw)
    if value === nothing
        add_issue!(result.validation,
                   invalid_value("line $line_no",
                                 "$field=$(raw): expected fixed-field Int in columns $first_col-$last_col"))
    end
    return value
end

function fixed_float_value(image::AbstractString, first_col::Int, last_col::Int)
    raw = fixed_field(image, first_col, last_col)
    isempty(raw) && return nothing
    return tryparse_deck_float(raw)
end

function fixed_int_value(image::AbstractString, first_col::Int, last_col::Int)
    raw = fixed_field(image, first_col, last_col)
    isempty(raw) && return nothing
    return tryparse(Int, raw)
end

function fixed_source_node_name(image::AbstractString, first_col::Int, last_col::Int)
    name = strip(fixed_field(image, first_col, last_col))
    return isempty(name) ? Symbol("") : Symbol(name)
end

function parse_fixed_source_load_flow_card!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
)::Bool
    image = fixed_image(line)
    constraint_code = fixed_int_value(image, 1, 2)
    if constraint_code !== nothing
        constraint_kind = get(
            Dict(
                0 => :active_reactive_power,
                1 => :active_power_voltage,
                2 => :angle_reactive_power,
            ),
            constraint_code,
            nothing,
        )
        if constraint_kind === nothing
            add_issue!(
                result.validation,
                invalid_value(
                    "line $line_no",
                    "FIX SOURCE constraint code must be 0 (P/Q), 1 (P/V), or 2 (angle/Q)",
                ),
            )
            return false
        end
        node_names = (
            fixed_source_node_name(image, 3, 8),
            fixed_source_node_name(image, 9, 14),
            fixed_source_node_name(image, 15, 20),
        )
        if isempty(String(node_names[1]))
            add_issue!(
                result.validation,
                missing_data("line $line_no", "FIX SOURCE requires a first source node name"),
            )
            return false
        end
        first_target = fixed_float_or_default!(
            result, image, line_no, 21, 36, "FIX SOURCE active power or angle", 0.0,
        )
        second_target = fixed_float_or_default!(
            result, image, line_no, 37, 52, "FIX SOURCE reactive power or voltage", 0.0,
        )
        minimum_voltage = fixed_float_or_default!(
            result, image, line_no, 53, 60, "FIX SOURCE minimum voltage", 0.0,
        )
        maximum_voltage = fixed_float_or_default!(
            result, image, line_no, 61, 68, "FIX SOURCE maximum voltage", 0.0,
        )
        minimum_angle_deg = fixed_float_or_default!(
            result, image, line_no, 69, 74, "FIX SOURCE minimum angle", 0.0,
        )
        maximum_angle_deg = fixed_float_or_default!(
            result, image, line_no, 75, 80, "FIX SOURCE maximum angle", 0.0,
        )
        any(isnothing, (
            first_target,
            second_target,
            minimum_voltage,
            maximum_voltage,
            minimum_angle_deg,
            maximum_angle_deg,
        )) && return false
        push!(
            result.fixed_source_constraint_rows,
            DeckFixedSourceConstraintRow(
                line_no,
                constraint_kind,
                node_names,
                constraint_kind == :angle_reactive_power ? missing : first_target,
                constraint_kind == :active_power_voltage ? missing : second_target,
                constraint_kind == :active_power_voltage ? second_target : missing,
                constraint_kind == :angle_reactive_power ? first_target : missing,
                minimum_voltage,
                maximum_voltage == 0.0 ? Inf : maximum_voltage,
                minimum_angle_deg == 0.0 ? -Inf : minimum_angle_deg,
                maximum_angle_deg == 0.0 ? Inf : maximum_angle_deg,
                String(line),
            ),
        )
        record_card!(result, :fixed_source_constraint)
        return false
    end

    print_voltage_changes = fixed_int_or_default!(
        result, image, line_no, 9, 16, "FIX SOURCE voltage-change print flag", 0,
    )
    maximum_iterations = fixed_int_or_default!(
        result, image, line_no, 17, 24, "FIX SOURCE maximum iterations", 0,
    )
    report_interval = fixed_int_or_default!(
        result, image, line_no, 25, 32, "FIX SOURCE report interval", 0,
    )
    print_final_sources = fixed_int_or_default!(
        result, image, line_no, 33, 40, "FIX SOURCE final-source print flag", 0,
    )
    relative_tolerance = fixed_float_or_default!(
        result, image, line_no, 41, 48, "FIX SOURCE relative power tolerance", 0.0,
    )
    voltage_factor = fixed_float_or_default!(
        result, image, line_no, 49, 56, "FIX SOURCE voltage correction factor", 0.0,
    )
    angle_factor = fixed_float_or_default!(
        result, image, line_no, 57, 64, "FIX SOURCE angle correction factor", 0.0,
    )
    any(isnothing, (
        print_voltage_changes,
        maximum_iterations,
        report_interval,
        print_final_sources,
        relative_tolerance,
        voltage_factor,
        angle_factor,
    )) && return true
    isempty(result.fixed_source_constraint_rows) && add_issue!(
        result.validation,
        missing_data("line $line_no", "FIX SOURCE requires at least one constraint card"),
    )
    push!(
        result.fixed_source_control_rows,
        DeckFixedSourceControlRow(
            line_no,
            print_voltage_changes != 0,
            maximum_iterations == 0 ? 500 : maximum_iterations,
            report_interval,
            print_final_sources != 0,
            relative_tolerance == 0.0 ? 0.01 : relative_tolerance,
            voltage_factor == 0.0 ? 0.2 : voltage_factor,
            angle_factor == 0.0 ? 2.5 : angle_factor,
            String(line),
        ),
    )
    record_card!(result, :fixed_source_control)
    return true
end

function validate_fixed_source_load_flow_rows!(result::DeckParseResult)
    constraints = result.fixed_source_constraint_rows
    controls = result.fixed_source_control_rows
    if isempty(constraints)
        isempty(controls) || add_issue!(
            result.validation,
            missing_data("FIX SOURCE", "control card has no preceding constraint rows"),
        )
        return result
    end
    length(controls) == 1 || add_issue!(
        result.validation,
        missing_data("FIX SOURCE", "constraints require exactly one trailing control card"),
    )
    source_rows = Dict(row.node => row for row in result.over5a_source_rows)
    constrained_nodes = Set{Symbol}()
    for constraint in constraints
        node_names = filter(name -> !isempty(String(name)), collect(constraint.source_node_names))
        length(node_names) in (1, 2, 3) || add_issue!(
            result.validation,
            invalid_value(
                "line $(constraint.line_no)",
                "FIX SOURCE constraint must name one scalar source or two or three phase sources",
            ),
        )
        allunique(node_names) || add_issue!(
            result.validation,
            invalid_value(
                "line $(constraint.line_no)",
                "FIX SOURCE constraint source names must be unique",
            ),
        )
        constraint.minimum_voltage <= constraint.maximum_voltage || add_issue!(
            result.validation,
            invalid_value(
                "line $(constraint.line_no)",
                "FIX SOURCE minimum voltage exceeds maximum voltage",
            ),
        )
        constraint.minimum_angle_deg <= constraint.maximum_angle_deg || add_issue!(
            result.validation,
            invalid_value(
                "line $(constraint.line_no)",
                "FIX SOURCE minimum angle exceeds maximum angle",
            ),
        )
        constraint.voltage_peak === missing || constraint.voltage_peak > 0.0 || add_issue!(
            result.validation,
            invalid_value(
                "line $(constraint.line_no)",
                "FIX SOURCE P/V target voltage must be positive",
            ),
        )
        for node_name in node_names
            source_row = get(source_rows, node_name, nothing)
            if source_row === nothing ||
               !(source_row.iform in (11, 14)) ||
               source_row.node_value <= 0
                add_issue!(
                    result.validation,
                    invalid_value(
                        "line $(constraint.line_no)",
                        "FIX SOURCE node $node_name must name a positive-node constant or sinusoidal voltage source",
                    ),
                )
            elseif source_row.iform == 11 &&
                   constraint.constraint_kind != :active_power_voltage
                add_issue!(
                    result.validation,
                    invalid_value(
                        "line $(constraint.line_no)",
                        "constant-source FIX SOURCE constraints require active power and voltage",
                    ),
                )
            end
            node_name in constrained_nodes && add_issue!(
                result.validation,
                invalid_value(
                    "line $(constraint.line_no)",
                    "FIX SOURCE node $node_name occurs in more than one constraint group",
                ),
            )
            push!(constrained_nodes, node_name)
        end
    end
    source_domains = unique(
        row.iform
        for row in values(source_rows)
        if row.iform in (11, 14) && row.node_value > 0
    )
    length(source_domains) <= 1 || add_issue!(
        result.validation,
        invalid_value(
            "FIX SOURCE",
            "constant and sinusoidal positive-node source domains cannot be mixed",
        ),
    )
    return result
end

function fixed_float_or_default!(result::DeckParseResult, image::AbstractString,
                                 line_no::Int, first_col::Int, last_col::Int,
                                 field::AbstractString, default::Real)
    raw = fixed_field(image, first_col, last_col)
    isempty(raw) && return Float64(default)
    value = tryparse_deck_float(raw)
    if value === nothing
        add_issue!(result.validation,
                   invalid_value("line $line_no",
                                 "$field=$(raw): expected fixed-field Float64 in columns $first_col-$last_col"))
        return nothing
    end
    return value
end

function fixed_int_or_default!(result::DeckParseResult, image::AbstractString,
                               line_no::Int, first_col::Int, last_col::Int,
                               field::AbstractString, default::Int)
    raw = fixed_field(image, first_col, last_col)
    isempty(raw) && return default
    value = tryparse(Int, raw)
    if value === nothing
        add_issue!(result.validation,
                   invalid_value("line $line_no",
                                 "$field=$(raw): expected fixed-field Int in columns $first_col-$last_col"))
        return nothing
    end
    return value
end

const CONTROL_SYSTEM_NUMBER_RE =
    r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][+-]?\d+)?"

function control_system_number_values(text::AbstractString)
    values = Float64[]
    for found in eachmatch(CONTROL_SYSTEM_NUMBER_RE, String(text))
        value = tryparse_deck_float(found.match)
        value === nothing || push!(values, value)
    end
    return values
end

control_system_optional_value(values::AbstractVector{Float64}, index::Int) =
    index <= length(values) ? values[index] : missing

function control_system_card_type(line::AbstractString)
    stripped = strip(String(line))
    isempty(stripped) && return nothing
    width = min(lastindex(stripped), 2)
    return tryparse(Int, stripped[1:width])
end

function control_system_fixed_name(line::AbstractString)
    name = fixed_field(fixed_image(line), 3, 8)
    return Symbol(name)
end

function control_system_signal_term(token::AbstractString)
    text = strip(String(token))
    isempty(text) && return nothing
    sign = first(text)
    (sign == '+' || sign == '-') || return nothing
    name_start = nextind(text, firstindex(text))
    name_start <= lastindex(text) || return nothing
    name_text = strip(text[name_start:lastindex(text)])
    isempty(name_text) && return nothing
    return DeckControlSystemSignalTerm(Symbol(name_text), sign == '-' ? -1 : 1)
end

function control_system_signal_terms(tokens)
    terms = DeckControlSystemSignalTerm[]
    for token in tokens
        term = control_system_signal_term(token)
        term === nothing && return nothing
        push!(terms, term)
    end
    return terms
end

function control_system_first_signed_signal(line::AbstractString)
    found = match(r"[+-]([A-Z0-9]+)", strip(String(line)))
    found === nothing && return missing
    return Symbol(found.captures[1])
end

function control_system_device_input_terms(line::AbstractString)
    image = fixed_image(line)
    terms = DeckControlSystemSignalTerm[]
    for first_col in 11:8:43
        sign = image[first_col]
        name = fixed_field(image, first_col + 1, first_col + 6)
        isempty(name) && continue
        polarity = sign == '-' ? -1 : sign == '*' ? 9 : sign == '+' ? 1 : 0
        polarity == 0 && continue
        push!(terms, DeckControlSystemSignalTerm(Symbol(name), polarity))
    end
    return terms
end

function control_system_optional_fixed_name(
    image::AbstractString,
    first_col::Int,
    last_col::Int,
)
    name = fixed_field(image, first_col, last_col)
    return isempty(name) ? missing : Symbol(name)
end

function control_system_tail_signal_names(line::AbstractString)
    text = String(line)
    tail = lastindex(text) >= 51 ? text[51:lastindex(text)] : ""
    without_numbers = replace(tail, r"[+-]?(?:\d+\.\d*|\.\d+)(?:[EeDd][+-]?\d+)?" => " ")
    compact = replace(uppercase(without_numbers), r"[^A-Z0-9]" => "")
    isempty(compact) && return Symbol[]
    names = Symbol[]
    for start in 1:6:lastindex(compact)
        stop = min(lastindex(compact), start + 5)
        name = strip(compact[start:stop])
        isempty(name) || push!(names, Symbol(name))
    end
    return names
end

function control_system_parameter_values(line::AbstractString)
    image = fixed_image(line)
    return Float64[
        something(fixed_float_value(image, first_col, first_col + 5), 0.0)
        for first_col in (51, 57, 63)
    ]
end

function control_system_output_signal_names(line::AbstractString, request_type::Int)
    request_type == 33 || return Symbol[]
    text = String(line)
    fields = lastindex(text) >= 3 ? text[3:lastindex(text)] : ""
    isempty(fields) && return Symbol[]
    names = Symbol[]
    for start in 1:6:lastindex(fields)
        stop = min(lastindex(fields), start + 5)
        name = strip(fields[start:stop])
        isempty(name) || push!(names, Symbol(name))
    end
    return names
end

function record_control_system_card!(result::DeckParseResult, kind::Symbol)
    record_card!(result, :fixed_field)
    record_card!(result, :control_system_card_input)
    record_card!(result, kind)
    return true
end

function record_control_system_card_blocker!(result::DeckParseResult, blocker::Symbol)
    record_fixed_blocker!(result, :control_system_card_blocked, blocker)
    record_card!(result, :control_system_card_input)
    return true
end

function parse_control_system_function_card!(
    result::DeckParseResult,
    line::AbstractString,
    tokens,
    line_no::Int,
)::Bool
    image = fixed_image(line)
    order_text = strip(fixed_field(image, 1, 2))
    order = isempty(order_text) ? 0 : tryparse(Int, order_text)
    order === nothing && return false
    if !(0 <= order <= 7)
        record_control_system_card_blocker!(result, :control_system_function_order)
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "control-system function order must be between zero and seven",
            ),
        )
        return true
    end
    terms = DeckControlSystemSignalTerm[]
    for sign_col in (11, 19, 27, 35, 43)
        signal_text = strip(fixed_field(image, sign_col + 1, sign_col + 6))
        isempty(signal_text) && continue
        sign = fixed_field(image, sign_col, sign_col)
        if sign == "+"
            push!(terms, DeckControlSystemSignalTerm(Symbol(signal_text), 1))
        elseif sign == "-"
            push!(terms, DeckControlSystemSignalTerm(Symbol(signal_text), -1))
        else
            return false
        end
    end
    if isempty(terms)
        length(tokens) >= 2 || return false
        terms = control_system_signal_terms(tokens[2:end])
    end
    terms === nothing && return false
    gain = something(fixed_float_value(image, 51, 56), 1.0)
    gain == 0.0 && (gain = 1.0)
    lower_limit = something(fixed_float_value(image, 57, 62), missing)
    upper_limit = something(fixed_float_value(image, 63, 68), missing)
    lower_limit_signal = control_system_optional_fixed_name(image, 69, 74)
    upper_limit_signal = control_system_optional_fixed_name(image, 75, 80)
    if lower_limit === 0.0 && upper_limit === 0.0 &&
       lower_limit_signal === missing && upper_limit_signal === missing
        lower_limit = missing
        upper_limit = missing
    end
    push!(
        result.control_system_function_rows,
        DeckControlSystemFunctionRow(
            line_no,
            control_system_fixed_name(line),
            terms,
            order,
            Float64(gain),
            order == 0 ? [1.0] : Float64[],
            order == 0 ? [1.0] : Float64[],
            lower_limit,
            upper_limit,
            lower_limit_signal,
            upper_limit_signal,
            order == 0,
            String(line),
        ),
    )
    return record_control_system_card!(result, :control_system_function_card)
end

function control_system_function_coefficients_pending(result::DeckParseResult)
    isempty(result.control_system_function_rows) && return false
    row = last(result.control_system_function_rows)
    return row.order > 0 && !row.coefficients_complete
end

function parse_control_system_function_coefficient_card!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
)
    row = last(result.control_system_function_rows)
    values = control_system_number_values(line)
    required = row.order + 1
    if length(values) < required
        record_control_system_card_blocker!(
            result,
            :control_system_function_coefficient_parse_error,
        )
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "order-$(row.order) control-system function requires $required coefficients",
            ),
        )
        return true
    end
    coefficients = Float64.(values[1:required])
    if isempty(row.numerator_coefficients)
        append!(row.numerator_coefficients, coefficients)
        return record_control_system_card!(
            result,
            :control_system_function_numerator_card,
        )
    end
    append!(row.denominator_coefficients, coefficients)
    row.coefficients_complete = true
    return record_control_system_card!(
        result,
        :control_system_function_denominator_card,
    )
end

function parse_control_system_expression_card!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    group_type::Int,
)::Bool
    image = fixed_image(line)
    fixed_field(image, 11, 11) == "=" || return false
    expression = strip(fixed_field(image, 12, 80))
    if isempty(expression)
        record_control_system_card_blocker!(result, :control_system_expression_empty)
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "control-system expression is empty"),
        )
        return true
    end
    push!(
        result.control_system_expression_rows,
        DeckControlSystemExpressionRow(
            line_no,
            group_type,
            control_system_fixed_name(line),
            expression,
            String(line),
        ),
    )
    return record_control_system_card!(result, :control_system_expression_card)
end

function parse_control_system_source_card!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    source_type::Int,
)::Bool
    text = String(line)
    values = control_system_number_values(lastindex(text) >= 9 ? text[9:lastindex(text)] : "")
    image = fixed_image(line)
    activation_start_time_s = something(fixed_float_value(image, 61, 70), 0.0)
    stop_value = fixed_float_value(image, 71, 80)
    activation_stop_time_s = stop_value === nothing || stop_value == 0.0 ?
        Inf : Float64(stop_value)
    push!(
        result.control_system_source_rows,
        DeckControlSystemSourceRow(
            line_no,
            source_type,
            control_system_fixed_name(line),
            control_system_optional_value(values, 1),
            control_system_optional_value(values, 2),
            control_system_optional_value(values, 3),
            activation_start_time_s,
            activation_stop_time_s,
            values,
            String(line),
        ),
    )
    return record_control_system_card!(result, :control_system_source_card)
end

function parse_control_system_device_card!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    group_type::Int,
)::Bool
    image = fixed_image(line)
    initial_issues = length(result.validation.issues)
    device_type = fixed_int_field!(
        result,
        image,
        line_no,
        9,
        10,
        "control_system_device_type",
    )
    if device_type === nothing || length(result.validation.issues) != initial_issues
        record_control_system_card_blocker!(result, :control_system_device_card_parse_error)
        return true
    end
    input_terms = control_system_device_input_terms(line)
    first_input = isempty(input_terms) ? missing : first(input_terms).name
    control_signal = control_system_optional_fixed_name(image, 69, 74)
    reference_signal = control_system_optional_fixed_name(image, 75, 80)
    tail_signal_names = Symbol[]
    control_signal === missing || push!(tail_signal_names, control_signal)
    reference_signal === missing || push!(tail_signal_names, reference_signal)
    push!(
        result.control_system_device_rows,
        DeckControlSystemDeviceRow(
            line_no,
            group_type,
            control_system_fixed_name(line),
            device_type,
            input_terms,
            first_input,
            tail_signal_names,
            control_signal,
            reference_signal,
            control_system_parameter_values(line),
            Float64[],
            Float64[],
            !(device_type in (55, 56, 57)),
            String(line),
        ),
    )
    return record_control_system_card!(result, :control_system_device_card)
end

function control_system_device_table_pending(result::DeckParseResult)
    isempty(result.control_system_device_rows) && return false
    row = last(result.control_system_device_rows)
    return row.device_type in (55, 56, 57) && !row.table_complete
end

function parse_control_system_device_table_card!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
)
    row = last(result.control_system_device_rows)
    image = fixed_image(line)
    input_value = fixed_float_value(image, 1, 16)
    if input_value === nothing
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "expected type-$(row.device_type) control-device table value in columns 1-16",
            ),
        )
        return record_control_system_card_blocker!(
            result,
            :control_system_device_table_parse_error,
        )
    end
    if input_value == 9999.0
        row.table_complete = true
        return record_control_system_card!(result, :control_system_device_table_end)
    end
    output_value = fixed_float_value(image, 17, 32)
    if row.device_type == 56 && output_value === nothing
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "type-56 point table requires an output value in columns 17-32",
            ),
        )
        return record_control_system_card_blocker!(
            result,
            :control_system_device_table_parse_error,
        )
    end
    push!(row.table_input_values, Float64(input_value))
    push!(row.table_output_values, something(output_value, 0.0))
    return record_control_system_card!(result, :control_system_device_table_value)
end

function parse_control_system_output_request_card!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    request_type::Int,
)::Bool
    push!(
        result.control_system_output_request_rows,
        DeckControlSystemOutputRequestRow(
            line_no,
            request_type,
            request_type == 1,
            control_system_output_signal_names(line, request_type),
            String(line),
        ),
    )
    return record_control_system_card!(result, :control_system_output_request_card)
end

function parse_control_system_hybrid_card!(
    result::DeckParseResult,
    line::AbstractString,
    tokens,
    line_no::Int,
)::Bool
    control_system_device_table_pending(result) &&
        return parse_control_system_device_table_card!(result, line, line_no)
    control_system_function_coefficients_pending(result) &&
        return parse_control_system_function_coefficient_card!(
            result,
            line,
            line_no,
        )
    row_type = control_system_card_type(line)
    if row_type !== nothing
        if (11 <= row_type <= 24) || (90 <= row_type <= 93)
            return parse_control_system_source_card!(result, line, line_no, row_type)
        elseif row_type in (88, 98, 99) &&
               parse_control_system_expression_card!(result, line, line_no, row_type)
            return true
        elseif row_type in (88, 98, 99)
            return parse_control_system_device_card!(result, line, line_no, row_type)
        elseif row_type in (1, 33)
            return parse_control_system_output_request_card!(result, line, line_no, row_type)
        end
    elseif parse_control_system_function_card!(result, line, tokens, line_no)
        return true
    end
    record_control_system_card_blocker!(result, :control_system_card_unsupported)
    add_issue!(
        result.validation,
        unknown_field("line $line_no", "Unsupported control-system fixed-field row in TACS HYBRID section"),
    )
    return true
end

function line_constants_scaled_float_text(raw::AbstractString, decimal_places::Int)
    text = strip(String(raw))
    isempty(text) && return nothing
    if occursin('.', text) || occursin('E', uppercase(text)) || occursin('D', uppercase(text))
        return tryparse_deck_float(text)
    end
    value = tryparse_deck_float(text)
    value === nothing && return nothing
    return value * 10.0^(-decimal_places)
end

function line_constants_float_field!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    first_col::Int,
    last_col::Int,
    decimal_places::Int,
    field::AbstractString,
)
    raw = fixed_field(image, first_col, last_col)
    isempty(raw) && return nothing
    value = line_constants_scaled_float_text(raw, decimal_places)
    if value === nothing
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "$field=$(raw): expected line-constants fixed-field numeric value in columns $first_col-$last_col",
            ),
        )
    end
    return value
end

function line_constants_float_or_default!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    first_col::Int,
    last_col::Int,
    decimal_places::Int,
    field::AbstractString,
    default::Real,
)
    value = line_constants_float_field!(
        result,
        image,
        line_no,
        first_col,
        last_col,
        decimal_places,
        field,
    )
    return value === nothing ? Float64(default) : Float64(value)
end

function line_constants_int_or_default!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    first_col::Int,
    last_col::Int,
    field::AbstractString,
    default::Int,
)
    raw = fixed_field(image, first_col, last_col)
    isempty(raw) && return default
    value = tryparse(Int, raw)
    if value === nothing
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "$field=$(raw): expected line-constants fixed-field Int in columns $first_col-$last_col",
            ),
        )
        return default
    end
    return value
end

line_constants_blank(image::AbstractString, first_col::Int, last_col::Int)::Bool =
    isempty(fixed_field(image, first_col, last_col))

function line_constants_inherited_int(
    result::DeckParseResult,
    raw_value::Int,
    image::AbstractString,
    first_col::Int,
    last_col::Int,
    field::Symbol,
)
    if line_constants_blank(image, first_col, last_col) &&
       !isempty(result.line_constants_conductor_cards)
        return Int(getfield(result.line_constants_conductor_cards[end], field))
    end
    return raw_value
end

function line_constants_inherited_float(
    result::DeckParseResult,
    raw_value::Float64,
    image::AbstractString,
    first_col::Int,
    last_col::Int,
    field::Symbol,
)
    if line_constants_blank(image, first_col, last_col) &&
       !isempty(result.line_constants_conductor_cards)
        return Float64(getfield(result.line_constants_conductor_cards[end], field))
    end
    return raw_value
end

function line_constants_option_card_kind(tokens)
    isempty(tokens) && return nothing
    first = normalized_deck_token(tokens[1])
    if first == "metric"
        return :line_constants_metric_units
    elseif first == "english"
        return :line_constants_english_units
    elseif first == "frequency"
        return :line_constants_frequency_loop_printout
    elseif first == "branch"
        return :line_constants_branch_names
    elseif first == "change"
        return :line_constants_change_case
    elseif first == "special" && length(tokens) >= 2 &&
           normalized_deck_token(tokens[2]) == "double"
        return :line_constants_special_double_circuit_transposed
    elseif first == "untransposed"
        return :line_constants_untransposed_modeling
    elseif first == "transposed"
        return :line_constants_transposed_modeling
    end
    return nothing
end

function parse_line_constants_conductor_card!(
    result::DeckParseResult,
    line::AbstractString,
    tokens,
    line_no::Int,
)::Bool
    option_kind = line_constants_option_card_kind(tokens)
    if option_kind !== nothing
        record_control_card!(result, option_kind, tokens, line_no)
        return true
    end

    image = fixed_image(line)
    initial_issues = length(result.validation.issues)
    phase_number = line_constants_int_or_default!(
        result,
        image,
        line_no,
        1,
        3,
        "line_constants_phase_number",
        0,
    )
    skin_effect_type = line_constants_float_or_default!(
        result,
        image,
        line_no,
        4,
        8,
        4,
        "line_constants_skin_effect_type",
        0.0,
    )
    resistance = line_constants_float_or_default!(
        result,
        image,
        line_no,
        9,
        16,
        5,
        "line_constants_resistance",
        0.0,
    )
    reactance_type = line_constants_int_or_default!(
        result,
        image,
        line_no,
        17,
        18,
        "line_constants_reactance_type",
        0,
    )
    reactance_or_gmr = line_constants_float_or_default!(
        result,
        image,
        line_no,
        19,
        26,
        5,
        "line_constants_reactance_or_gmr",
        0.0,
    )
    diameter = line_constants_float_or_default!(
        result,
        image,
        line_no,
        27,
        34,
        5,
        "line_constants_diameter",
        0.0,
    )
    horizontal = line_constants_float_or_default!(
        result,
        image,
        line_no,
        35,
        42,
        3,
        "line_constants_horizontal_position",
        0.0,
    )
    tower_height = line_constants_float_or_default!(
        result,
        image,
        line_no,
        43,
        50,
        3,
        "line_constants_tower_height",
        0.0,
    )
    midspan_height = line_constants_float_or_default!(
        result,
        image,
        line_no,
        51,
        58,
        3,
        "line_constants_midspan_height",
        0.0,
    )
    bundle_spacing = line_constants_float_or_default!(
        result,
        image,
        line_no,
        59,
        66,
        5,
        "line_constants_bundle_spacing",
        0.0,
    )
    bundle_angle = line_constants_float_or_default!(
        result,
        image,
        line_no,
        67,
        72,
        2,
        "line_constants_bundle_angle",
        0.0,
    )
    conductor_name = fixed_field(image, 73, 78)
    bundle_count = line_constants_int_or_default!(
        result,
        image,
        line_no,
        79,
        80,
        "line_constants_bundle_count",
        0,
    )

    if length(result.validation.issues) != initial_issues
        record_fixed_blocker!(
            result,
            :line_constants_conductor_card_blocked,
            :line_constants_conductor_card_parse_error,
        )
        return true
    end

    phase_number = line_constants_inherited_int(
        result,
        phase_number,
        image,
        1,
        3,
        :phase_number,
    )
    skin_effect_type = line_constants_inherited_float(
        result,
        skin_effect_type,
        image,
        4,
        8,
        :skin_effect_type,
    )
    resistance = line_constants_inherited_float(
        result,
        resistance,
        image,
        9,
        16,
        :resistance_ohm_per_mile,
    )
    reactance_type = line_constants_inherited_int(
        result,
        reactance_type,
        image,
        17,
        18,
        :reactance_type,
    )
    reactance_or_gmr = line_constants_inherited_float(
        result,
        reactance_or_gmr,
        image,
        19,
        26,
        :reactance_or_gmr,
    )
    diameter = line_constants_inherited_float(
        result,
        diameter,
        image,
        27,
        34,
        :diameter_inches,
    )
    if line_constants_blank(image, 43, 50) && tower_height == 0.0
        tower_height = midspan_height
    end
    if line_constants_blank(image, 51, 58) && midspan_height == 0.0
        midspan_height = tower_height
    end
    average_height = (tower_height + 2.0 * midspan_height) / 3.0

    push!(
        result.line_constants_conductor_cards,
        DeckLineConstantsConductorCard(
            line_no,
            phase_number,
            skin_effect_type,
            resistance,
            reactance_type,
            reactance_or_gmr,
            diameter,
            horizontal,
            tower_height,
            midspan_height,
            average_height,
            bundle_spacing,
            bundle_angle,
            conductor_name,
            bundle_count,
        ),
    )
    record_card!(result, :fixed_field)
    record_card!(result, :line_constants_conductor_card)
    record_card!(result, Symbol("line_constants_conductor_phase_", phase_number))
    if bundle_count > 1
        record_card!(result, :line_constants_bundled_conductor_card)
    end
    for _ in 1:_line_constants_conductor_bundle_count(result.line_constants_conductor_cards[end])
        record_card!(result, :line_constants_physical_conductor)
    end
    return true
end

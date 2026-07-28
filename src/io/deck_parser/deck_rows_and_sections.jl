
function deck_output_width_request_rows(result::DeckParseResult)
    return DeckOutputWidthRequestRow[
        DeckOutputWidthRequestRow(row.line_no, row.column_width, row.raw_text)
        for row in result.output_width_request_rows
    ]
end

function deck_peak_voltage_monitor_request_rows(result::DeckParseResult)
    return DeckPeakVoltageMonitorRequestRow[
        DeckPeakVoltageMonitorRequestRow(row.line_no, row.raw_text)
        for row in result.peak_voltage_monitor_request_rows
    ]
end

function deck_diagnostic_print_request_rows(result::DeckParseResult)
    return DeckDiagnosticPrintRequestRow[
        DeckDiagnosticPrintRequestRow(row.line_no, row.loop_print_controls, row.raw_text)
        for row in result.diagnostic_print_request_rows
    ]
end

function deck_tacs_warning_limit_request_rows(result::DeckParseResult)
    return DeckTACSWarningLimitRequestRow[
        DeckTACSWarningLimitRequestRow(
            row.line_no,
            row.warning_limit,
            row.begin_time_s,
            row.raw_text,
        )
        for row in result.tacs_warning_limit_request_rows
    ]
end

function deck_plot_file_request_rows(result::DeckParseResult)
    return DeckPlotFileRequestRow[
        DeckPlotFileRequestRow(row.line_no, row.plot_file_mode, row.raw_text)
        for row in result.plot_file_request_rows
    ]
end

function deck_switch_logic_request_rows(result::DeckParseResult)
    return DeckSwitchLogicRequestRow[
        DeckSwitchLogicRequestRow(row.line_no, row.control_value, row.raw_text)
        for row in result.switch_logic_request_rows
    ]
end

function deck_simulation_control_request_rows(result::DeckParseResult)
    return DeckSimulationControlRequestRow[
        DeckSimulationControlRequestRow(
            row.line_no,
            row.request_kind,
            row.numeric_value,
            row.integer_value,
            row.raw_text,
        )
        for row in result.simulation_control_request_rows
    ]
end

function deck_printout_frequency_change_rows(result::DeckParseResult)
    return DeckPrintoutFrequencyChangeRow[
        DeckPrintoutFrequencyChangeRow(
            row.request_line_no,
            row.payload_line_no,
            row.change_steps,
            row.multipliers,
            row.active_pair_count,
            row.provided_value_count,
            row.request_text,
            row.payload_text,
        )
        for row in result.printout_frequency_change_rows
    ]
end

function deck_study_option_request_rows(result::DeckParseResult)
    return DeckStudyOptionRequestRow[
        DeckStudyOptionRequestRow(
            row.line_no,
            row.request_kind,
            copy(row.numeric_values),
            copy(row.integer_values),
            copy(row.text_values),
            row.raw_text,
        )
        for row in result.study_option_request_rows
    ]
end

function deck_fixed_source_constraint_rows(result::DeckParseResult)
    return copy(result.fixed_source_constraint_rows)
end

function deck_fixed_source_control_rows(result::DeckParseResult)
    return copy(result.fixed_source_control_rows)
end

function deck_time_horizon_control_rows(result::DeckParseResult)
    return DeckTimeHorizonControlRow[
        DeckTimeHorizonControlRow(
            row.line_no,
            row.time_step_s,
            row.stop_time_s,
            row.inductance_frequency_hz,
            row.capacitance_frequency_hz,
            row.epsilon,
            row.matrix_tolerance,
            row.start_time_s,
            row.raw_text,
        )
        for row in result.time_horizon_control_rows
    ]
end

function deck_output_schedule_control_rows(result::DeckParseResult)
    return DeckOutputScheduleControlRow[
        DeckOutputScheduleControlRow(
            row.line_no,
            row.print_interval_steps,
            row.plot_interval_steps,
            row.network_print_enabled,
            row.steady_state_print_enabled,
            row.extrema_print_enabled,
            row.terminal_conditions_punch_enabled,
            row.restart_snapshot_enabled,
            row.plot_file_retention_mode,
            row.energization_count,
            row.special_request_control_word,
            row.raw_text,
        )
        for row in result.output_schedule_control_rows
    ]
end

function deck_case_boundary_rows(result::DeckParseResult)
    return DeckCaseBoundaryRow[
        DeckCaseBoundaryRow(row.line_no, row.boundary_kind, row.raw_text)
        for row in result.case_boundary_rows
    ]
end

function _control_system_signal_terms_copy(terms::AbstractVector{DeckControlSystemSignalTerm})
    return DeckControlSystemSignalTerm[
        DeckControlSystemSignalTerm(term.name, term.polarity) for term in terms
    ]
end

function deck_control_system_function_rows(result::DeckParseResult)
    return DeckControlSystemFunctionRow[
        DeckControlSystemFunctionRow(
            row.line_no,
            row.name,
            _control_system_signal_terms_copy(row.input_terms),
            row.order,
            row.gain,
            copy(row.numerator_coefficients),
            copy(row.denominator_coefficients),
            row.lower_limit,
            row.upper_limit,
            row.lower_limit_signal,
            row.upper_limit_signal,
            row.coefficients_complete,
            row.raw_text,
        )
        for row in result.control_system_function_rows
    ]
end

function deck_control_system_expression_rows(result::DeckParseResult)
    return DeckControlSystemExpressionRow[
        DeckControlSystemExpressionRow(
            row.line_no,
            row.group_type,
            row.name,
            row.expression,
            row.raw_text,
        )
        for row in result.control_system_expression_rows
    ]
end

function deck_control_system_source_rows(result::DeckParseResult)
    return DeckControlSystemSourceRow[
        DeckControlSystemSourceRow(
            row.line_no,
            row.source_type,
            row.name,
            row.amplitude,
            row.delay_or_time_constant,
            row.phase_or_width,
            row.activation_start_time_s,
            row.activation_stop_time_s,
            copy(row.numeric_values),
            row.raw_text,
        )
        for row in result.control_system_source_rows
    ]
end

function deck_control_system_device_rows(result::DeckParseResult)
    return DeckControlSystemDeviceRow[
        DeckControlSystemDeviceRow(
            row.line_no,
            row.group_type,
            row.name,
            row.device_type,
            copy(row.input_terms),
            row.first_input,
            copy(row.tail_signal_names),
            row.control_signal,
            row.reference_signal,
            copy(row.parameter_values),
            copy(row.table_input_values),
            copy(row.table_output_values),
            row.table_complete,
            row.raw_text,
        )
        for row in result.control_system_device_rows
    ]
end

function deck_control_system_output_request_rows(result::DeckParseResult)
    return DeckControlSystemOutputRequestRow[
        DeckControlSystemOutputRequestRow(
            row.line_no,
            row.request_type,
            row.all_signals,
            copy(row.signal_names),
            row.raw_text,
        )
        for row in result.control_system_output_request_rows
    ]
end

function deck_control_system_switch_rows(result::DeckParseResult)
    return DeckControlSystemSwitchRow[
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
            row.delayed_arc,
            row.parameter_source_kind,
            row.parameter_reference_index,
            row.parameter_reference_line_no,
            row.layout_kind,
            row.event_output_code,
            row.output_code,
            row.raw_text,
        )
        for row in result.control_system_switch_rows
    ]
end

function _control_system_output_signal_lookup(
    rows::AbstractVector{DeckControlSystemOutputRequestRow},
)
    lookup = Dict{Symbol,Tuple{Int,Int,Int}}()
    linear_index = 0
    for (request_index, row) in enumerate(rows)
        for signal_index in eachindex(row.signal_names)
            linear_index += 1
            signal = row.signal_names[signal_index]
            haskey(lookup, signal) || (lookup[signal] = (request_index, signal_index, linear_index))
        end
    end
    return lookup
end

_deck_node_index_or_missing(result::DeckParseResult, node::Symbol) =
    get(result.node_map, node, missing)

function deck_control_system_switch_coupling_rows(result::DeckParseResult)
    output_lookup = _control_system_output_signal_lookup(result.control_system_output_request_rows)
    couplings = DeckControlSystemSwitchCouplingRow[]
    for row in result.control_system_switch_rows
        output_position = get(output_lookup, row.control_signal, nothing)
        push!(
            couplings,
            DeckControlSystemSwitchCouplingRow(
                row.line_no,
                row.switch_type,
                row.from_node,
                row.to_node,
                _deck_node_index_or_missing(result, row.from_node),
                _deck_node_index_or_missing(result, row.to_node),
                row.control_signal,
                row.gate_signal,
                row.clamp_signal,
                output_position === nothing ? missing : output_position[1],
                output_position === nothing ? missing : output_position[2],
                output_position === nothing ? missing : output_position[3],
                row.output_code,
                row.initial_state,
            ),
        )
    end
    return couplings
end

function deck_node_initial_condition_rows(result::DeckParseResult)
    return DeckNodeInitialConditionRow[
        DeckNodeInitialConditionRow(
            row.line_no,
            row.condition_kind,
            row.node,
            row.node_index,
            row.reference_node,
            row.reference_node_index,
            row.real_value,
            row.imaginary_value,
            row.raw_text,
        )
        for row in result.node_initial_condition_rows
    ]
end

_sequence_numeric_value(value::Missing) = 0.0
_sequence_numeric_value(value::Real) = Float64(value)

function _coupled_phase_pi_row_groups(
    rows::AbstractVector{DeckCoupledPhasePiSectionRow},
)
    groups = Vector{Vector{DeckCoupledPhasePiSectionRow}}()
    index = 1
    while index <= length(rows)
        row = rows[index]
        if row.phase_index != 1
            index += 1
            continue
        end
        candidate = DeckCoupledPhasePiSectionRow[row]
        next_index = index + 1
        while next_index <= length(rows)
            next_row = rows[next_index]
            next_row.phase_index == length(candidate) + 1 || break
            push!(candidate, next_row)
            next_index += 1
        end
        push!(groups, candidate)
        index = next_index
    end
    return groups
end

function _coupled_phase_pi_section_reference(
    accepted_sections::AbstractVector{CoupledLumpedPhasePiSection},
    row::DeckCoupledPhasePiSectionRow,
)
    row.reference_from_node_value isa Missing && return nothing
    row.reference_to_node_value isa Missing && return nothing
    for section in Iterators.reverse(accepted_sections)
        isempty(section.from_node_indices) && continue
        if section.from_node_indices[1] == row.reference_from_node_value &&
           section.to_node_indices[1] == row.reference_to_node_value
            return section
        end
    end
    return nothing
end

function _copy_coupled_phase_pi_section(
    source::CoupledLumpedPhasePiSection,
    rows::AbstractVector{DeckCoupledPhasePiSectionRow},
    group_index::Int,
)
    phase_count = length(rows)
    phase_count == source.phase_count ||
        throw(ArgumentError("coupled phase PI COPY row count does not match reference section"))
    return coupled_lumped_phase_pi_section(
        Symbol("coupled_lumped_phase_pi_section_", group_index),
        [row.phase_index for row in rows],
        [row.from_node for row in rows],
        [row.to_node for row in rows],
        [row.from_node_value for row in rows],
        [row.to_node_value for row in rows],
        source.phase_resistance_matrix,
        source.phase_inductance_matrix,
        source.phase_capacitance_matrix;
        line_numbers = [row.line_no for row in rows],
    )
end

function _coupled_phase_pi_section_from_rows(
    result::DeckParseResult,
    rows::AbstractVector{DeckCoupledPhasePiSectionRow},
    group_index::Int,
    accepted_sections::AbstractVector{CoupledLumpedPhasePiSection},
)
    phase_count = length(rows)
    phase_count > 0 || throw(ArgumentError("coupled phase PI section group is empty"))
    [row.phase_index for row in rows] == collect(1:phase_count) ||
        throw(ArgumentError("coupled phase PI rows must be contiguous from phase 1"))
    first_row = first(rows)
    if first_row.reference_kind != :none
        source = _coupled_phase_pi_section_reference(accepted_sections, first_row)
        source === nothing &&
            throw(ArgumentError("coupled phase PI COPY reference does not match a prior section"))
        return _copy_coupled_phase_pi_section(source, rows, group_index)
    end
    resistance = zeros(Float64, phase_count, phase_count)
    inductance = zeros(Float64, phase_count, phase_count)
    capacitance = zeros(Float64, phase_count, phase_count)
    for row in rows
        expected = row.phase_index
        length(row.raw_resistance_values) >= expected ||
            throw(ArgumentError("coupled phase PI row $(row.line_no) is missing resistance values"))
        length(row.raw_inductance_values) >= expected ||
            throw(ArgumentError("coupled phase PI row $(row.line_no) is missing inductance values"))
        length(row.raw_capacitance_values) >= expected ||
            throw(ArgumentError("coupled phase PI row $(row.line_no) is missing capacitance values"))
        row_index = row.phase_index
        for column in 1:row_index
            resistance_value = row.raw_resistance_values[column]
            inductance_value = fixed_card_branch_timestep_inductance(
                result,
                row.raw_inductance_values[column],
            )
            capacitance_value = fixed_card_branch_timestep_capacitance(
                result,
                row.raw_capacitance_values[column],
            )
            resistance[row_index, column] = resistance_value
            resistance[column, row_index] = resistance_value
            inductance[row_index, column] = inductance_value
            inductance[column, row_index] = inductance_value
            capacitance[row_index, column] = capacitance_value
            capacitance[column, row_index] = capacitance_value
        end
    end
    return coupled_lumped_phase_pi_section(
        Symbol("coupled_lumped_phase_pi_section_", group_index),
        [row.phase_index for row in rows],
        [row.from_node for row in rows],
        [row.to_node for row in rows],
        [row.from_node_value for row in rows],
        [row.to_node_value for row in rows],
        resistance,
        inductance,
        capacitance;
        line_numbers = [row.line_no for row in rows],
    )
end

function deck_coupled_lumped_phase_pi_sections(result::DeckParseResult)
    cascaded_row_indices = Set(
        Iterators.flatten(
            row.source_section_row_indices for row in result.cascaded_pi_request_rows
        ),
    )
    rows = [
        row
        for (index, row) in enumerate(result.coupled_phase_pi_section_rows)
        if !(index in cascaded_row_indices)
    ]
    groups = _coupled_phase_pi_row_groups(rows)
    sections = CoupledLumpedPhasePiSection[]
    for (index, group) in enumerate(groups)
        push!(sections, _coupled_phase_pi_section_from_rows(result, group, index, sections))
    end
    return sections
end

deck_cascaded_pi_request_rows(result::DeckParseResult) =
    copy(result.cascaded_pi_request_rows)

function _cascaded_pi_phase_ordered_section(
    section::CoupledLumpedPhasePiSection,
    phase_map::Vector{Int},
    name::Symbol;
    resistance_matrix::AbstractMatrix{<:Real}=section.phase_resistance_matrix,
    inductance_matrix::AbstractMatrix{<:Real}=section.phase_inductance_matrix,
    capacitance_matrix::AbstractMatrix{<:Real}=section.phase_capacitance_matrix,
)
    phase_count = section.phase_count
    sort(phase_map) == collect(1:phase_count) ||
        throw(ArgumentError("cascade phase order must be a complete permutation"))
    return coupled_lumped_phase_pi_section(
        name,
        collect(1:phase_count),
        section.from_nodes,
        section.to_nodes,
        section.from_node_indices,
        section.to_node_indices,
        resistance_matrix[phase_map, phase_map],
        inductance_matrix[phase_map, phase_map],
        capacitance_matrix[phase_map, phase_map];
        line_numbers = section.line_numbers,
    )
end

function _cascaded_pi_explicit_section_matrices(
    result::DeckParseResult,
    block::DeckCascadedPiBlock,
)
    block.explicit_resistance_values === nothing && return nothing
    resistance = copy(block.explicit_resistance_values)
    inductance = map(
        value -> fixed_card_branch_timestep_inductance(result, value),
        block.explicit_inductance_values,
    )
    capacitance = map(
        value -> fixed_card_branch_timestep_capacitance(result, value),
        block.explicit_capacitance_values,
    )
    return (; resistance, inductance, capacitance)
end

function _cascaded_pi_local_terminal(terminal::Int, phase_map::Vector{Int})
    terminal == 0 && return 0
    terminal < 0 && return terminal
    source_phase = abs(terminal)
    local_phase = findfirst(==(source_phase), phase_map)
    local_phase === nothing &&
        throw(ArgumentError("cascade shunt terminal is absent from its phase order"))
    return local_phase
end

function _cascaded_pi_model_block(
    result::DeckParseResult,
    source_section::CoupledLumpedPhasePiSection,
    block::DeckCascadedPiBlock,
    block_index::Int,
)
    explicit = _cascaded_pi_explicit_section_matrices(result, block)
    section = if explicit === nothing
        _cascaded_pi_phase_ordered_section(
            source_section,
            block.phase_map,
            Symbol("cascaded_phase_pi_block_", block_index, "_section"),
        )
    else
        _cascaded_pi_phase_ordered_section(
            source_section,
            block.phase_map,
            Symbol("cascaded_phase_pi_block_", block_index, "_explicit_section");
            resistance_matrix = explicit.resistance,
            inductance_matrix = explicit.inductance,
            capacitance_matrix = explicit.capacitance,
        )
    end
    inverse_phase_map = invperm(block.phase_map)
    series_impedances = CascadedPiSeriesImpedance[]
    for row in block.series_impedances
        push!(
            series_impedances,
            cascaded_pi_series_impedance(
                inverse_phase_map[row.phase_index],
                row.raw_resistance_ohm,
                fixed_card_branch_timestep_inductance(
                    result,
                    row.raw_inductance_value,
                ),
                fixed_card_branch_timestep_capacitance(
                    result,
                    row.raw_capacitance_value,
                );
                open_circuit = row.open_circuit,
            ),
        )
    end
    shunt_impedances = CascadedPiShuntImpedance[]
    for row in block.shunt_impedances
        push!(
            shunt_impedances,
            cascaded_pi_shunt_impedance(
                _cascaded_pi_local_terminal(row.from_terminal, block.phase_map),
                _cascaded_pi_local_terminal(row.to_terminal, block.phase_map),
                row.raw_resistance_ohm,
                fixed_card_branch_timestep_inductance(
                    result,
                    row.raw_inductance_value,
                ),
                fixed_card_branch_timestep_capacitance(
                    result,
                    row.raw_capacitance_value,
                ),
            ),
        )
    end
    return cascaded_phase_pi_block(
        section,
        block.multiplicity,
        block.section_scale;
        series_impedances = series_impedances,
        shunt_impedances = shunt_impedances,
    )
end

function _deck_cascaded_phase_pi_equivalent(
    result::DeckParseResult,
    request::DeckCascadedPiRequestRow,
    request_index::Int,
)
    rows = result.coupled_phase_pi_section_rows[request.source_section_row_indices]
    all(row -> row.reference_kind == :none, rows) ||
        throw(ArgumentError("cascaded PI source section must own direct matrix values"))
    section = _coupled_phase_pi_section_from_rows(
        result,
        rows,
        request_index,
        CoupledLumpedPhasePiSection[],
    )
    blocks = CascadedPhasePiBlock[
        _cascaded_pi_model_block(result, section, block, block_index)
        for (block_index, block) in enumerate(request.blocks)
    ]
    return cascaded_phase_pi_equivalent(
        request.name,
        blocks,
        request.frequency_hz,
    )
end

function deck_cascaded_phase_pi_equivalents(result::DeckParseResult)
    return CascadedPhasePiEquivalent[
        _deck_cascaded_phase_pi_equivalent(result, request, request_index)
        for (request_index, request) in enumerate(result.cascaded_pi_request_rows)
    ]
end

function _coupled_lumped_sequence_row_groups(rows::AbstractVector{DeckCoupledLineRow})
    groups = Vector{Vector{DeckCoupledLineRow}}()
    index = 1
    while index <= length(rows)
        row = rows[index]
        if row.line_kind != :mutual_source_equivalent
            index += 1
            continue
        end
        if index + 2 <= length(rows)
            candidate = rows[index:(index + 2)]
            if all(candidate_row -> candidate_row.line_kind == :mutual_source_equivalent,
                   candidate) &&
               [candidate_row.phase_index for candidate_row in candidate] == [1, 2, 3]
                push!(groups, collect(candidate))
                index += 3
                continue
            end
        end
        index += 1
    end
    return groups
end

function _coupled_lumped_sequence_impedance_from_rows(
    rows::AbstractVector{DeckCoupledLineRow},
    group_index::Int,
)
    length(rows) == 3 || throw(ArgumentError("coupled lumped sequence group must contain three rows"))
    third_row_resistance = _sequence_numeric_value(rows[3].sequence_resistance)
    third_row_inductance = _sequence_numeric_value(rows[3].sequence_inductance)
    if third_row_resistance != 0.0 || third_row_inductance != 0.0
        throw(ArgumentError("explicit third-row coupled R-L matrix input is not translated yet"))
    end
    return coupled_lumped_sequence_impedance(
        Symbol("coupled_lumped_sequence_impedance_", group_index),
        [row.phase_index for row in rows],
        [row.from_node for row in rows],
        [row.to_node for row in rows],
        [row.from_node_value for row in rows],
        [row.to_node_value for row in rows],
        _sequence_numeric_value(rows[1].sequence_resistance),
        _sequence_numeric_value(rows[2].sequence_resistance),
        _sequence_numeric_value(rows[1].sequence_inductance),
        _sequence_numeric_value(rows[2].sequence_inductance);
        line_numbers = [row.line_no for row in rows],
    )
end

function deck_coupled_lumped_sequence_impedances(result::DeckParseResult)
    groups = _coupled_lumped_sequence_row_groups(result.coupled_line_rows)
    return CoupledLumpedSequenceImpedance[
        _coupled_lumped_sequence_impedance_from_rows(group, index)
        for (index, group) in enumerate(groups)
    ]
end

function _distributed_transposed_line_row_groups(rows::AbstractVector{DeckCoupledLineRow})
    groups = Vector{Vector{DeckCoupledLineRow}}()
    index = 1
    while index <= length(rows)
        row = rows[index]
        if row.line_kind != :distributed_transmission_line ||
           row.sampled_frequency_data_requested
            index += 1
            continue
        end
        if index + 2 <= length(rows)
            candidate = rows[index:(index + 2)]
            if all(candidate_row ->
                       candidate_row.line_kind == :distributed_transmission_line &&
                       !candidate_row.sampled_frequency_data_requested,
                   candidate) &&
               [candidate_row.phase_index for candidate_row in candidate] == [1, 2, 3]
                push!(groups, collect(candidate))
                index += 3
                continue
            end
        end
        index += 1
    end
    return groups
end

function deck_fixed_time_horizon_options(result::DeckParseResult)
    for row in result.time_horizon_control_rows
        return (
            line_no = row.line_no,
            dt_s = row.time_step_s,
            tmax_s = row.stop_time_s,
            x_frequency_hz = row.inductance_frequency_hz,
            c_frequency_hz = row.capacitance_frequency_hz,
            epsilon = row.epsilon,
            tolerance = row.matrix_tolerance,
            start_time_s = row.start_time_s,
        )
    end
    return (
        line_no = 0,
        dt_s = 0.0,
        tmax_s = 0.0,
        x_frequency_hz = 0.0,
        c_frequency_hz = 0.0,
        epsilon = 0.0,
        tolerance = 0.0,
        start_time_s = 0.0,
    )
end

function fixed_card_branch_timestep_inductance(
    result::DeckParseResult,
    raw_inductance::Float64,
)
    raw_inductance == 0.0 && return 0.0
    options = deck_fixed_time_horizon_options(result)
    x_frequency_hz = Float64(options.x_frequency_hz)
    return x_frequency_hz > 0.0 ?
        raw_inductance / (2.0 * pi * x_frequency_hz) :
        raw_inductance / 1000.0
end

function fixed_card_branch_timestep_capacitance(
    result::DeckParseResult,
    raw_capacitance::Float64,
    ;
    legacy_microfarad_units::Bool = true,
)
    raw_capacitance == 0.0 && return 0.0
    legacy_microfarad_units || return raw_capacitance
    options = deck_fixed_time_horizon_options(result)
    c_frequency_hz = Float64(options.c_frequency_hz)
    capacitance = raw_capacitance / 1.0e6
    return c_frequency_hz > 0.0 ?
        capacitance / (2.0 * pi * c_frequency_hz) :
        capacitance
end

function _integer_control_value(value)
    value === nothing && return 0
    return Int(round(Float64(value)))
end

function deck_output_schedule_options(result::DeckParseResult)
    for row in result.output_schedule_control_rows
        return (
            line_no = row.line_no,
            print_interval_steps = row.print_interval_steps,
            plot_interval_steps = row.plot_interval_steps,
            network_print_enabled = row.network_print_enabled,
            steady_state_print_enabled = row.steady_state_print_enabled,
            extrema_print_enabled = row.extrema_print_enabled,
            terminal_conditions_punch_enabled =
                row.terminal_conditions_punch_enabled,
            restart_snapshot_enabled = row.restart_snapshot_enabled,
            plot_file_retention_mode = row.plot_file_retention_mode,
            energization_count = row.energization_count,
            special_request_control_word = row.special_request_control_word,
        )
    end
    return (
        line_no = 0,
        print_interval_steps = 0,
        plot_interval_steps = 0,
        network_print_enabled = false,
        steady_state_print_enabled = false,
        extrema_print_enabled = false,
        terminal_conditions_punch_enabled = false,
        restart_snapshot_enabled = false,
        plot_file_retention_mode = 0,
        energization_count = 0,
        special_request_control_word = 0,
    )
end

function _deck_nominal_steady_state_frequency_hz(result::DeckParseResult)
    options = deck_fixed_time_horizon_options(result)
    options.x_frequency_hz > 0.0 && return options.x_frequency_hz
    for row in result.over5a_source_rows
        row.sfreq > 0.0 && return Float64(row.sfreq) / (2.0 * pi)
    end
    return 60.0
end

function _required_line_numeric_value(value, row::DeckCoupledLineRow, field::Symbol)
    ismissing(value) &&
        throw(ArgumentError("distributed line row $(row.line_no) missing $field"))
    return Float64(value)
end

function _distributed_line_length(rows::AbstractVector{DeckCoupledLineRow})
    length_value = _required_line_numeric_value(rows[1].line_length, rows[1], :line_length)
    for row in rows[2:3]
        ismissing(row.line_length) && continue
        candidate = Float64(row.line_length)
        candidate == length_value ||
            throw(ArgumentError("distributed transposed line rows must share line length"))
    end
    return length_value
end

function _line_row_terminal_key(rows::AbstractVector{DeckCoupledLineRow})
    return Tuple((row.from_node, row.to_node) for row in rows)
end

function _line_row_reference_key(rows::AbstractVector{DeckCoupledLineRow})
    return Tuple((row.reference_from_node, row.reference_to_node) for row in rows)
end

function _kc_lee_untransposed_line_row_groups(rows::AbstractVector{DeckCoupledLineRow})
    groups = Vector{Vector{DeckCoupledLineRow}}()
    index = 1
    while index <= length(rows)
        row = rows[index]
        if row.line_kind != :kc_lee_untransposed_line
            index += 1
            continue
        end
        group = DeckCoupledLineRow[]
        expected_phase = 1
        while index <= length(rows)
            candidate = rows[index]
            candidate.line_kind == :kc_lee_untransposed_line || break
            candidate.phase_index == expected_phase || break
            push!(group, candidate)
            index += 1
            expected_phase += 1
        end
        if isempty(group)
            index += 1
        else
            push!(groups, group)
        end
    end
    return groups
end

function _kc_lee_group_at(rows::AbstractVector{DeckCoupledLineRow}, index::Int)
    index <= length(rows) || return DeckCoupledLineRow[]
    rows[index].line_kind == :kc_lee_untransposed_line || return DeckCoupledLineRow[]
    group = DeckCoupledLineRow[]
    expected_phase = 1
    while index <= length(rows)
        row = rows[index]
        row.line_kind == :kc_lee_untransposed_line || break
        row.phase_index == expected_phase || break
        push!(group, row)
        index += 1
        expected_phase += 1
    end
    return group
end

function _complete_kc_lee_modal_parameter_group(rows::AbstractVector{DeckCoupledLineRow})
    isempty(rows) && return false
    return all(
        row -> !ismissing(row.raw_resistance) &&
               !ismissing(row.raw_inductance) &&
               !ismissing(row.raw_capacitance) &&
               !ismissing(row.line_length),
        rows,
    )
end

function _complete_kc_lee_reference_group(rows::AbstractVector{DeckCoupledLineRow})
    isempty(rows) && return false
    return all(
        row -> !ismissing(row.reference_from_node) &&
               !ismissing(row.reference_to_node),
        rows,
    )
end

function _kc_lee_modal_transform_rows_between(
    result::DeckParseResult,
    first_line_no::Int,
    next_line_no::Int,
)
    return DeckLineModalTransformRow[
        row for row in result.line_modal_transform_rows
        if first_line_no < row.line_no < next_line_no
    ]
end

function _kc_lee_real_modal_to_phase_transform(
    rows::AbstractVector{DeckLineModalTransformRow},
    phase_count::Int,
)
    expected_row_count = 2 * phase_count
    length(rows) == expected_row_count ||
        throw(ArgumentError("modal transform requires $expected_row_count real/imaginary rows"))
    all(row -> length(row.values) == phase_count, rows) ||
        throw(ArgumentError("modal transform rows must each contain $phase_count values"))
    modal_to_phase = zeros(Float64, phase_count, phase_count)
    for phase_index in 1:phase_count
        real_row = rows[2 * phase_index - 1]
        imaginary_row = rows[2 * phase_index]
        maximum(abs.(imaginary_row.values)) <= 1.0e-12 ||
            throw(ArgumentError("complex modal transforms are not yet executable in the real-time line owner"))
        modal_to_phase[phase_index, :] .= real_row.values
    end
    return modal_to_phase
end

function _kc_lee_modal_parameter_data(
    result::DeckParseResult,
    rows::AbstractVector{DeckCoupledLineRow},
    transform_rows::AbstractVector{DeckLineModalTransformRow},
)
    phase_count = length(rows)
    phases = [row.phase_index for row in rows]
    phases == collect(1:phase_count) ||
        throw(ArgumentError("modal line rows must have contiguous phase indices"))
    line_length = abs(_required_line_numeric_value(rows[1].line_length, rows[1], :line_length))
    line_length > 0.0 || throw(ArgumentError("modal line length must be positive"))
    signed_characteristic = Float64[]
    total_resistance = Float64[]
    propagation_times = Float64[]
    for row in rows
        row_length = abs(_required_line_numeric_value(row.line_length, row, :line_length))
        abs(row_length - line_length) <= max(1.0e-12, 64.0 * eps(Float64) * line_length) ||
            throw(ArgumentError("modal line rows must share line length"))
        resistance_per_length =
            _required_line_numeric_value(row.raw_resistance, row, :raw_resistance)
        characteristic =
            _required_line_numeric_value(row.raw_inductance, row, :raw_inductance)
        speed = _required_line_numeric_value(row.raw_capacitance, row, :raw_capacitance)
        speed > 0.0 || throw(ArgumentError("modal propagation speed must be positive"))
        characteristic != 0.0 ||
            throw(ArgumentError("modal characteristic impedance must be nonzero"))
        punched_constant_parameter_row =
            !ismissing(row.line_length) && Float64(row.line_length) < 0.0
        push!(
            signed_characteristic,
            resistance_per_length == 0.0 || punched_constant_parameter_row ?
                -abs(characteristic) :
                characteristic,
        )
        push!(total_resistance, resistance_per_length * line_length)
        push!(propagation_times, line_length / speed)
    end
    return (
        signed_characteristic_impedances = signed_characteristic,
        total_resistances = total_resistance,
        propagation_times_s = propagation_times,
        modal_to_phase_transform =
            _kc_lee_real_modal_to_phase_transform(transform_rows, phase_count),
    )
end

function _kc_lee_modal_state_from_group(
    rows::AbstractVector{DeckCoupledLineRow},
    data,
    group_index::Int,
)
    phase_count = length(rows)
    return distributed_modal_line_branch_state(
        Symbol("distributed_modal_line_branch_state_", group_index),
        collect(1:phase_count),
        [row.from_node for row in rows],
        [row.to_node for row in rows],
        [row.from_node_value for row in rows],
        [row.to_node_value for row in rows],
        data.signed_characteristic_impedances,
        data.total_resistances,
        data.propagation_times_s,
        data.modal_to_phase_transform;
        line_numbers = [row.line_no for row in rows],
        modal_admittance_denominator = 1.0,
    )
end

function _single_phase_distributed_line_rows(rows::AbstractVector{DeckCoupledLineRow})
    selected = DeckCoupledLineRow[]
    index = 1
    while index <= length(rows)
        row = rows[index]
        if row.line_kind != :distributed_transmission_line ||
           row.sampled_frequency_data_requested
            index += 1
            continue
        end
        if index + 2 <= length(rows)
            candidate = rows[index:(index + 2)]
            if all(candidate_row ->
                       candidate_row.line_kind == :distributed_transmission_line &&
                       !candidate_row.sampled_frequency_data_requested,
                   candidate) &&
               [candidate_row.phase_index for candidate_row in candidate] == [1, 2, 3]
                index += 3
                continue
            end
        end
        if row.phase_index == 1
            push!(selected, row)
        end
        index += 1
    end
    return selected
end

function _single_phase_distributed_line_modal_state_from_row(
    row::DeckCoupledLineRow,
    state_index::Int,
)
    line_length = abs(_required_line_numeric_value(row.line_length, row, :line_length))
    line_length > 0.0 || throw(ArgumentError("single-phase distributed line length must be positive"))
    resistance_per_length =
        _required_line_numeric_value(row.raw_resistance, row, :raw_resistance)
    characteristic =
        _required_line_numeric_value(row.raw_inductance, row, :raw_inductance)
    speed = _required_line_numeric_value(row.raw_capacitance, row, :raw_capacitance)
    characteristic != 0.0 ||
        throw(ArgumentError("single-phase distributed line characteristic impedance must be nonzero"))
    speed > 0.0 ||
        throw(ArgumentError("single-phase distributed line propagation speed must be positive"))
    signed_characteristic =
        resistance_per_length == 0.0 || (!ismissing(row.line_length) && Float64(row.line_length) < 0.0) ?
        -abs(characteristic) :
        characteristic
    return distributed_modal_line_branch_state(
        Symbol("single_phase_distributed_line_branch_state_", state_index),
        [1],
        [row.from_node],
        [row.to_node],
        [row.from_node_value],
        [row.to_node_value],
        [signed_characteristic],
        [resistance_per_length * line_length],
        [line_length / speed],
        reshape([1.0], 1, 1);
        line_numbers = [row.line_no],
        modal_admittance_denominator = 1.0,
    )
end

function _deck_single_phase_distributed_line_modal_branch_states(result::DeckParseResult)
    rows = _single_phase_distributed_line_rows(result.coupled_line_rows)
    return DistributedTransposedLineModalBranchState[
        _single_phase_distributed_line_modal_state_from_row(row, index)
        for (index, row) in enumerate(rows)
    ]
end

function deck_sampled_frequency_line_coefficients(
    result::DeckParseResult,
    timestep_s::Real,
)
    coefficients = SampledLineWeightingCoefficients[]
    configured_epsilon = deck_fixed_time_horizon_options(result).epsilon
    tail_relative_tolerance = configured_epsilon > 0.0 ? configured_epsilon : 1.0e-8
    for row in result.sampled_frequency_line_rows
        push!(
            coefficients,
            sampled_line_weighting_coefficients(
                line_weighting_samples(
                    row.propagation_time_s,
                    row.propagation_amplitude,
                ),
                line_weighting_samples(
                    row.admittance_time_s,
                    row.admittance_amplitude,
                ),
                timestep_s,
                row.characteristic_impedance_ohm;
                propagation_peak_index = row.propagation_peak_index,
                admittance_rise_index = row.admittance_rise_index,
                propagation_cutoff_fraction = row.propagation_cutoff_fraction,
                admittance_cutoff_fraction = row.admittance_cutoff_fraction,
                total_resistance_ohm = row.total_resistance_ohm,
                maximum_tail_iterations = row.maximum_tail_iterations,
                tail_relative_tolerance = tail_relative_tolerance,
            ),
        )
    end
    return coefficients
end

function deck_sampled_frequency_line_elements(
    result::DeckParseResult,
    timestep_s::Real,
)
    coefficients = deck_sampled_frequency_line_coefficients(result, timestep_s)
    elements = Any[]
    row_index = 1
    while row_index <= length(result.sampled_frequency_line_rows)
        row = result.sampled_frequency_line_rows[row_index]
        branch = result.coupled_line_rows[row.branch_row_index]
        if branch.phase_index == 1 && row_index + 2 <= length(result.sampled_frequency_line_rows)
            group_rows = result.sampled_frequency_line_rows[row_index:(row_index + 2)]
            group_branches = DeckCoupledLineRow[
                result.coupled_line_rows[group_row.branch_row_index]
                for group_row in group_rows
            ]
            if [group_branch.phase_index for group_branch in group_branches] == [1, 2, 3]
                loss_factors = Float64[
                    _sampled_frequency_line_loss_factor(group_branch)
                    for group_branch in group_branches
                ]
                push!(
                    elements,
                    sampled_frequency_dependent_line_group(
                        [group_row.from_node_index for group_row in group_rows],
                        [group_row.to_node_index for group_row in group_rows],
                        coefficients[row_index:(row_index + 2)];
                        loss_factors = loss_factors,
                    ),
                )
                row_index += 3
                continue
            end
        end
        branch.phase_index == 1 || throw(ArgumentError(
            "sampled frequency-dependent line phase rows must be grouped as [-1, -2, -3]",
        ))
        loss = _sampled_frequency_line_loss_factor(branch)
        push!(
            elements,
            sampled_frequency_dependent_line(
                row.from_node_index,
                row.to_node_index,
                coefficients[row_index];
                loss_factor = loss,
            ),
        )
        row_index += 1
    end
    return elements
end

function _sampled_frequency_line_loss_factor(branch::DeckCoupledLineRow)
    resistance_per_length = ismissing(branch.raw_resistance) ?
        0.0 : abs(Float64(branch.raw_resistance))
    line_length = ismissing(branch.line_length) ?
        0.0 : abs(Float64(branch.line_length))
    loss = max(resistance_per_length * line_length, 1.0e-5)
    loss < 1.0 || throw(ArgumentError(
        "sampled frequency-dependent line loss coefficient must be below one",
    ))
    return loss
end

function deck_sampled_frequency_line_element_names(result::DeckParseResult)
    names = Symbol[]
    row_index = 1
    group_index = 0
    while row_index <= length(result.sampled_frequency_line_rows)
        row = result.sampled_frequency_line_rows[row_index]
        branch = result.coupled_line_rows[row.branch_row_index]
        if branch.phase_index == 1 && row_index + 2 <= length(result.sampled_frequency_line_rows)
            group_rows = result.sampled_frequency_line_rows[row_index:(row_index + 2)]
            phases = Int[
                result.coupled_line_rows[group_row.branch_row_index].phase_index
                for group_row in group_rows
            ]
            if phases == [1, 2, 3]
                group_index += 1
                push!(names, Symbol("sampled_frequency_line_group_", group_index))
                row_index += 3
                continue
            end
        end
        push!(names, row.name)
        row_index += 1
    end
    return names
end

function _terminal_surge_phase_one_branch(
    result::DeckParseResult,
    row::DeckCoupledLineRow,
)
    isempty(result.over2_branch_rows) && return nothing
    for candidate in Iterators.reverse(result.over2_branch_rows)
        candidate.line_no < row.line_no || continue
        candidate.branch_type == 1 || return nothing
        candidate.branch_kind == :conductance || return nothing
        candidate.source_kind == :scalar || return nothing
        candidate.to_node_value == 0 || return nothing
        candidate.raw_resistance > 0.0 || return nothing
        candidate.raw_inductance == 0.0 || return nothing
        candidate.raw_capacitance == 0.0 || return nothing
        return candidate
    end
    return nothing
end

function _terminal_surge_impedance_admittance_from_row!(
    result::DeckParseResult,
    row::DeckCoupledLineRow,
    group_index::Int,
)
    phase_one = _terminal_surge_phase_one_branch(result, row)
    phase_one === nothing &&
        throw(ArgumentError("terminal surge type-2 row requires a preceding type-1 scalar surge impedance"))
    phase_two_self = _required_line_numeric_value(row.raw_inductance, row, :raw_inductance)
    mutual = _required_line_numeric_value(row.raw_resistance, row, :raw_resistance)
    phase_two_self > 0.0 ||
        throw(ArgumentError("terminal surge type-2 self impedance must be positive"))
    impedance = [
        phase_one.raw_resistance mutual
        mutual phase_two_self
    ]
    preexisting = zeros(Float64, 2, 2)
    preexisting[1, 1] = phase_one.conductance
    return terminal_surge_impedance_admittance(
        Symbol("terminal_surge_impedance_admittance_", group_index),
        [phase_one.from_node, row.from_node],
        [phase_one.from_node_value, row.from_node_value],
        impedance;
        line_numbers = [phase_one.line_no, row.line_no],
        preexisting_admittance_matrix = preexisting,
    )
end

function accept_terminal_surge_impedance_row!(
    result::DeckParseResult,
    row::DeckCoupledLineRow,
    group_index::Int,
)
    admittance = _terminal_surge_impedance_admittance_from_row!(result, row, group_index)
    push!(result.elements, admittance)
    push!(result.element_line_numbers, first(admittance.line_numbers))
    push!(result.element_names, admittance.name)
    record_card!(result, :fixed_card_terminal_surge_impedance_admittance)
    record_card!(result, :fixed_card_terminal_surge_impedance_coupled_correction)
    return admittance
end

function accept_single_phase_distributed_line_row!(
    result::DeckParseResult,
    row::DeckCoupledLineRow,
    state_index::Int,
)
    state = _single_phase_distributed_line_modal_state_from_row(row, state_index)
    distributed_transposed_line_history_state(
        state;
        timestep_s = deck_fixed_time_horizon_options(result).dt_s,
        steady_state_frequency_hz = _deck_nominal_steady_state_frequency_hz(result),
        history_storage_start_index = 1,
        initialized_from_steady_state = false,
        name = Symbol("single_phase_distributed_line_history_state_", state_index),
    )
    admittance = distributed_transposed_line_companion_admittance(
        state;
        name = Symbol("single_phase_distributed_line_companion_admittance_", state_index),
    )
    push!(result.elements, admittance)
    push!(result.element_line_numbers, first(admittance.line_numbers))
    push!(result.element_names, admittance.name)
    record_card!(result, :fixed_card_single_phase_distributed_line_modal_branch_state)
    record_card!(result, :fixed_card_single_phase_distributed_line_history_state)
    record_card!(result, :fixed_card_single_phase_distributed_line_modal_timestep_update)
    record_card!(result, :fixed_card_single_phase_distributed_line_phase_current_injection)
    record_card!(result, :fixed_card_single_phase_distributed_line_companion_admittance)
    return admittance
end

function _deck_kc_lee_modal_branch_states(result::DeckParseResult)
    groups = _kc_lee_untransposed_line_row_groups(result.coupled_line_rows)
    isempty(groups) && return DistributedTransposedLineModalBranchState[]
    source_data = Dict{Any,Any}()
    states = DistributedTransposedLineModalBranchState[]
    for (group_index, group) in enumerate(groups)
        next_line_no =
            group_index == length(groups) ? typemax(Int) : first(groups[group_index + 1]).line_no
        if _complete_kc_lee_modal_parameter_group(group)
            transform_rows = _kc_lee_modal_transform_rows_between(
                result,
                last(group).line_no,
                next_line_no,
            )
            data = _kc_lee_modal_parameter_data(result, group, transform_rows)
            source_data[_line_row_terminal_key(group)] = data
        elseif _complete_kc_lee_reference_group(group)
            reference_key = _line_row_reference_key(group)
            haskey(source_data, reference_key) ||
                throw(ArgumentError("modal line copy group references an unknown source span"))
            data = source_data[reference_key]
        else
            throw(ArgumentError("modal line group requires parameters or reference terminals"))
        end
        push!(states, _kc_lee_modal_state_from_group(group, data, group_index))
    end
    return states
end

function accept_kc_lee_modal_line_groups!(result::DeckParseResult)
    states = _deck_kc_lee_modal_branch_states(result)
    isempty(states) && return result
    for (index, state) in enumerate(states)
        admittance = distributed_transposed_line_companion_admittance(
            state;
            name = Symbol("distributed_line_companion_admittance_", index),
        )
        push!(result.elements, admittance)
        push!(result.element_line_numbers, first(admittance.line_numbers))
        push!(result.element_names, admittance.name)
        record_card!(result, :fixed_card_kc_lee_modal_branch_state)
        record_card!(result, :fixed_card_kc_lee_history_state)
        record_card!(result, :fixed_card_kc_lee_modal_timestep_update)
        record_card!(result, :fixed_card_kc_lee_phase_current_injection)
        record_card!(result, :fixed_card_kc_lee_companion_admittance)
    end
    return result
end

function _distributed_transposed_line_constants_from_rows(
    result::DeckParseResult,
    rows::AbstractVector{DeckCoupledLineRow},
    group_index::Int,
)
    length(rows) == 3 || throw(ArgumentError("distributed transposed line group must contain three rows"))
    third_sequence_values = (
        rows[3].raw_resistance,
        rows[3].raw_inductance,
        rows[3].raw_capacitance,
    )
    explicit_third_sequence =
        any(value -> !ismissing(value) && Float64(value) != 0.0, third_sequence_values)
    options = deck_fixed_time_horizon_options(result)
    return distributed_transposed_line_constants(
        Symbol("distributed_transposed_line_constants_", group_index),
        [row.phase_index for row in rows],
        [row.from_node for row in rows],
        [row.to_node for row in rows],
        [row.from_node_value for row in rows],
        [row.to_node_value for row in rows],
        _required_line_numeric_value(rows[1].raw_resistance, rows[1], :raw_resistance),
        _required_line_numeric_value(rows[2].raw_resistance, rows[2], :raw_resistance),
        _required_line_numeric_value(rows[1].raw_inductance, rows[1], :raw_inductance),
        _required_line_numeric_value(rows[2].raw_inductance, rows[2], :raw_inductance),
        _required_line_numeric_value(rows[1].raw_capacitance, rows[1], :raw_capacitance),
        _required_line_numeric_value(rows[2].raw_capacitance, rows[2], :raw_capacitance),
        _distributed_line_length(rows);
        x_frequency_hz = options.x_frequency_hz,
        c_frequency_hz = options.c_frequency_hz,
        line_numbers = [row.line_no for row in rows],
        third_sequence_resistance_per_length = explicit_third_sequence ?
                                               _required_line_numeric_value(
            rows[3].raw_resistance,
            rows[3],
            :raw_resistance,
        ) : nothing,
        third_sequence_inductance_input = explicit_third_sequence ?
                                          _required_line_numeric_value(
            rows[3].raw_inductance,
            rows[3],
            :raw_inductance,
        ) : nothing,
        third_sequence_capacitance_input = explicit_third_sequence ?
                                           _required_line_numeric_value(
            rows[3].raw_capacitance,
            rows[3],
            :raw_capacitance,
        ) : nothing,
    )
end

function deck_distributed_transposed_line_constants(result::DeckParseResult)
    groups = _distributed_transposed_line_row_groups(result.coupled_line_rows)
    return DistributedTransposedLineConstants[
        _distributed_transposed_line_constants_from_rows(result, group, index)
        for (index, group) in enumerate(groups)
    ]
end

function deck_distributed_transposed_line_modal_branch_states(result::DeckParseResult)
    states = DistributedTransposedLineModalBranchState[
        distributed_transposed_line_modal_branch_state(
            constants;
            name = Symbol("distributed_transposed_line_modal_branch_state_", index),
        )
        for (index, constants) in enumerate(deck_distributed_transposed_line_constants(result))
    ]
    append!(states, _deck_single_phase_distributed_line_modal_branch_states(result))
    append!(states, _deck_kc_lee_modal_branch_states(result))
    return states
end

function deck_distributed_transposed_line_steady_state_pi_equivalents(
    result::DeckParseResult;
    first_storage_index::Integer=1,
    steady_state_frequency_hz::Union{Nothing,Real}=nothing,
)
    frequency = steady_state_frequency_hz === nothing ?
        _deck_nominal_steady_state_frequency_hz(result) :
        Float64(steady_state_frequency_hz)
    storage_index = Int(first_storage_index)
    equivalents = DistributedTransposedLineSteadyStatePiEquivalent[]
    for (index, constants) in enumerate(deck_distributed_transposed_line_constants(result))
        equivalent = distributed_transposed_line_steady_state_pi_equivalent(
            constants;
            steady_state_frequency_hz = frequency,
            storage_start_index = storage_index,
            name = Symbol("distributed_transposed_line_steady_state_pi_equivalent_", index),
        )
        push!(equivalents, equivalent)
        storage_index = equivalent.storage_end_index + 1
    end
    return equivalents
end

function deck_distributed_transposed_line_history_states(
    result::DeckParseResult;
    first_history_storage_index::Integer=1,
    steady_state_frequency_hz::Union{Nothing,Real}=nothing,
    initialized_from_steady_state::Bool=false,
)
    modal_states = deck_distributed_transposed_line_modal_branch_states(result)
    isempty(modal_states) && return DistributedTransposedLineHistoryState[]
    options = deck_fixed_time_horizon_options(result)
    options.dt_s > 0.0 ||
        throw(ArgumentError("distributed transposed line history states require a positive fixed timestep"))
    frequency = steady_state_frequency_hz === nothing ?
        _deck_nominal_steady_state_frequency_hz(result) :
        Float64(steady_state_frequency_hz)
    storage_index = Int(first_history_storage_index)
    histories = DistributedTransposedLineHistoryState[]
    for (index, state) in enumerate(modal_states)
        history = distributed_transposed_line_history_state(
            state;
            timestep_s = options.dt_s,
            steady_state_frequency_hz = frequency,
            history_storage_start_index = storage_index,
            initialized_from_steady_state = initialized_from_steady_state,
            name = Symbol("distributed_transposed_line_history_state_", index),
        )
        push!(histories, history)
        storage_index = history.next_history_storage_index
    end
    return histories
end

function deck_distributed_transposed_line_companion_admittances(result::DeckParseResult)
    return DistributedTransposedLineCompanionAdmittance[
        distributed_transposed_line_companion_admittance(
            state;
            name = Symbol("distributed_line_companion_admittance_", index),
        )
        for (index, state) in enumerate(deck_distributed_transposed_line_modal_branch_states(result))
    ]
end

function deck_over2_branch_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.over2_branch_rows)
end

function deck_zinc_oxide_nonlinear_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.zinc_oxide_nonlinear_rows)
end

function deck_zinc_oxide_initialization_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.zinc_oxide_initialization_rows)
end

function deck_zinc_oxide_breakpoint_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.zinc_oxide_breakpoint_rows)
end

function deck_nonlinear_resistance_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.nonlinear_resistance_rows)
end

function deck_nonlinear_resistance_initialization_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.nonlinear_resistance_initialization_rows)
end

function deck_nonlinear_resistance_point_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.nonlinear_resistance_point_rows)
end

function deck_triggered_timed_resistance_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.triggered_timed_resistance_rows)
end

function deck_triggered_timed_resistance_point_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.triggered_timed_resistance_point_rows)
end

function deck_switching_nonlinear_resistor_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.switching_nonlinear_resistor_rows)
end

function deck_switching_nonlinear_resistor_point_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.switching_nonlinear_resistor_point_rows)
end

function deck_piecewise_nonlinear_inductor_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.piecewise_nonlinear_inductor_rows)
end

function deck_piecewise_nonlinear_inductor_point_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.piecewise_nonlinear_inductor_point_rows)
end

function deck_hysteretic_inductor_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.hysteretic_inductor_rows)
end

function deck_hysteretic_inductor_point_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.hysteretic_inductor_point_rows)
end

function deck_arrester_nonlinear_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.arrester_nonlinear_rows)
end

function deck_arrester_constant_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.arrester_constant_rows)
end

deck_over2_branch_names(result::DeckParseResult) =
    Symbol[row.name for row in deck_over2_branch_rows(result)]

deck_over2_branch_line_numbers(result::DeckParseResult) =
    Int[row.line_no for row in deck_over2_branch_rows(result)]

deck_over2_branch_kinds(result::DeckParseResult) =
    Symbol[row.branch_kind for row in deck_over2_branch_rows(result)]

deck_over2_branch_layout_kinds(result::DeckParseResult) =
    Symbol[row.layout_kind for row in deck_over2_branch_rows(result)]

deck_over2_branch_source_kinds(result::DeckParseResult) =
    Symbol[row.source_kind for row in deck_over2_branch_rows(result)]

deck_over2_branch_reference_kinds(result::DeckParseResult) =
    Symbol[row.reference_kind for row in deck_over2_branch_rows(result)]

deck_over2_branch_reference_names(result::DeckParseResult) =
    Symbol[row.reference_name for row in deck_over2_branch_rows(result)]

deck_over2_branch_reference_line_numbers(result::DeckParseResult) =
    Int[row.reference_line_no for row in deck_over2_branch_rows(result)]

deck_over2_branch_from_node_names(result::DeckParseResult) =
    Symbol[row.from_node for row in deck_over2_branch_rows(result)]

deck_over2_branch_to_node_names(result::DeckParseResult) =
    Symbol[row.to_node for row in deck_over2_branch_rows(result)]

deck_over2_branch_from_node_indices(result::DeckParseResult) =
    Int[row.from_node_value for row in deck_over2_branch_rows(result)]

deck_over2_branch_to_node_indices(result::DeckParseResult) =
    Int[row.to_node_value for row in deck_over2_branch_rows(result)]

deck_over2_branch_raw_resistance_values(result::DeckParseResult) =
    Float64[row.raw_resistance for row in deck_over2_branch_rows(result)]

deck_over2_branch_raw_inductance_values(result::DeckParseResult) =
    Float64[row.raw_inductance for row in deck_over2_branch_rows(result)]

deck_over2_branch_raw_capacitance_values(result::DeckParseResult) =
    Float64[row.raw_capacitance for row in deck_over2_branch_rows(result)]

deck_over2_branch_conductance_values(result::DeckParseResult) =
    Float64[row.conductance for row in deck_over2_branch_rows(result)]

deck_over2_branch_resistance_values(result::DeckParseResult) =
    Float64[row.resistance for row in deck_over2_branch_rows(result)]

deck_over2_branch_inductance_values(result::DeckParseResult) =
    Float64[row.inductance for row in deck_over2_branch_rows(result)]

deck_over2_branch_capacitance_values(result::DeckParseResult) =
    Float64[row.capacitance for row in deck_over2_branch_rows(result)]

deck_over2_branch_output_codes(result::DeckParseResult) =
    Int[row.output_code for row in deck_over2_branch_rows(result)]

function deck_model_summary(result::DeckParseResult)
    source_inputs = deck_over5a_source_update_inputs(result)
    return DeckModelSummary(
        result.source,
        deck_node_names(result),
        copy(result.element_names),
        node_count(result),
        length(result.elements),
        deck_branch_names(result),
        deck_branch_kinds(result),
        deck_branch_from_node_names(result),
        deck_branch_to_node_names(result),
        deck_branch_from_node_indices(result),
        deck_branch_to_node_indices(result),
        deck_branch_conductance_values(result),
        deck_branch_resistance_values(result),
        deck_branch_inductance_values(result),
        deck_branch_capacitance_values(result),
        deck_branch_previous_current_values(result),
        deck_branch_previous_voltage_values(result),
        deck_branch_line_numbers(result),
        deck_branch_count(result),
        deck_bergeron_line_names(result),
        deck_bergeron_line_line_numbers(result),
        deck_bergeron_line_from_node_names(result),
        deck_bergeron_line_to_node_names(result),
        deck_bergeron_line_from_node_indices(result),
        deck_bergeron_line_to_node_indices(result),
        deck_bergeron_line_surge_impedance_values(result),
        deck_bergeron_line_surge_admittance_values(result),
        deck_bergeron_line_travel_time_s_values(result),
        deck_bergeron_line_dt_s_values(result),
        deck_bergeron_line_attenuation_values(result),
        deck_bergeron_line_delay_step_counts(result),
        deck_bergeron_line_write_indices(result),
        deck_bergeron_line_history_current_from_values(result),
        deck_bergeron_line_history_current_to_values(result),
        deck_bergeron_line_terminal_voltage_from_values(result),
        deck_bergeron_line_terminal_voltage_to_values(result),
        deck_bergeron_line_terminal_current_from_values(result),
        deck_bergeron_line_terminal_current_to_values(result),
        deck_bergeron_line_traveling_wave_from_values(result),
        deck_bergeron_line_traveling_wave_to_values(result),
        length(result.bergeron_line_rows),
        deck_coupled_lumped_sequence_impedances(result),
        length(deck_coupled_lumped_sequence_impedances(result)),
        deck_coupled_lumped_phase_pi_sections(result),
        length(deck_coupled_lumped_phase_pi_sections(result)),
        deck_distributed_transposed_line_constants(result),
        length(deck_distributed_transposed_line_constants(result)),
        deck_distributed_transposed_line_modal_branch_states(result),
        length(deck_distributed_transposed_line_modal_branch_states(result)),
        deck_distributed_transposed_line_steady_state_pi_equivalents(result),
        length(deck_distributed_transposed_line_steady_state_pi_equivalents(result)),
        deck_distributed_transposed_line_history_states(result),
        length(deck_distributed_transposed_line_history_states(result)),
        deck_distributed_transposed_line_companion_admittances(result),
        length(deck_distributed_transposed_line_companion_admittances(result)),
        deck_over2_branch_names(result),
        deck_over2_branch_line_numbers(result),
        deck_over2_branch_kinds(result),
        deck_over2_branch_layout_kinds(result),
        deck_over2_branch_source_kinds(result),
        deck_over2_branch_reference_kinds(result),
        deck_over2_branch_reference_names(result),
        deck_over2_branch_reference_line_numbers(result),
        deck_over2_branch_from_node_names(result),
        deck_over2_branch_to_node_names(result),
        deck_over2_branch_from_node_indices(result),
        deck_over2_branch_to_node_indices(result),
        deck_over2_branch_raw_resistance_values(result),
        deck_over2_branch_raw_inductance_values(result),
        deck_over2_branch_raw_capacitance_values(result),
        deck_over2_branch_conductance_values(result),
        deck_over2_branch_resistance_values(result),
        deck_over2_branch_inductance_values(result),
        deck_over2_branch_capacitance_values(result),
        deck_over2_branch_output_codes(result),
        length(result.over2_branch_rows),
        deck_over15_output_request_names(result),
        deck_over15_output_request_output_kinds(result),
        deck_over15_output_request_request_kinds(result),
        deck_over15_output_request_layout_kinds(result),
        deck_over15_output_request_line_numbers(result),
        deck_over15_output_request_output_codes(result),
        deck_over15_output_request_node_names(result),
        deck_over15_output_request_node_indices(result),
        deck_over15_output_request_branch_names(result),
        deck_over15_output_request_branch_indices(result),
        length(result.over15_output_request_rows),
        deck_over16_output_channel_names(result),
        deck_over16_output_node_names(result),
        deck_over16_output_node_indices(result),
        deck_over16_output_channel_line_numbers(result),
        length(result.over16_output_channels),
        deck_over16_branch_voltage_output_names(result),
        deck_over16_branch_voltage_branch_names(result),
        deck_over16_branch_voltage_branch_indices(result),
        deck_over16_branch_voltage_output_line_numbers(result),
        length(result.over16_branch_voltage_outputs),
        deck_over16_branch_current_output_names(result),
        deck_over16_branch_current_branch_names(result),
        deck_over16_branch_current_branch_indices(result),
        deck_over16_branch_current_output_line_numbers(result),
        length(result.over16_branch_current_outputs),
        deck_over16_branch_power_output_names(result),
        deck_over16_branch_power_branch_names(result),
        deck_over16_branch_power_branch_indices(result),
        deck_over16_branch_power_output_line_numbers(result),
        length(result.over16_branch_power_outputs),
        copy(source_inputs.names),
        copy(source_inputs.nodes),
        copy(source_inputs.node_values),
        copy(source_inputs.iform_values),
        copy(source_inputs.line_numbers),
        copy(source_inputs.layout_kinds),
        copy(source_inputs.tstart_values),
        copy(source_inputs.tstop_values),
        copy(source_inputs.crest_values),
        copy(source_inputs.time1_values),
        copy(source_inputs.time2_values),
        copy(source_inputs.sfreq_values),
        source_inputs.kconst,
        deck_over16_source_card_kinds(result),
        deck_over16_source_card_values(result),
        deck_over16_source_card_provided_value_counts(result),
        deck_over16_source_card_line_numbers(result),
        length(result.over16_source_card_rows),
        deck_over16_source_interpolation_values(result),
        deck_over16_source_interpolation_provided_value_counts(result),
        deck_over16_source_interpolation_line_numbers(result),
        length(result.over16_source_interpolation_rows),
        deck_over16_source_tacs_override_positions(result),
        deck_over16_source_tacs_override_xtcs_indices(result),
        deck_over16_source_tacs_override_line_numbers(result),
        length(result.over16_source_tacs_override_rows),
        deck_over16_source_analytic_values(result),
        deck_over16_source_analytic_provided_value_counts(result),
        deck_over16_source_analytic_line_numbers(result),
        length(result.over16_source_analytic_rows),
        deck_over5_switch_names(result),
        deck_over5_switch_line_numbers(result),
        deck_over5_switch_from_node_names(result),
        deck_over5_switch_to_node_names(result),
        deck_over5_switch_from_node_indices(result),
        deck_over5_switch_to_node_indices(result),
        deck_over5_switch_layout_kinds(result),
        deck_over5_switch_raw_close_time_s_values(result),
        deck_over5_switch_raw_open_time_s_values(result),
        deck_over5_switch_close_time_s_values(result),
        deck_over5_switch_open_time_s_values(result),
        deck_over5_switch_initially_closed_flags(result),
        deck_over5_switch_measuring_flags(result),
        deck_over5_switch_closed_markers(result),
        deck_over5_switch_marker_texts(result),
        deck_over5_switch_on_conductance_values(result),
        deck_over5_switch_off_conductance_values(result),
        deck_over5_switch_output_codes(result),
        length(result.over5_switch_rows),
        deck_time_switch_names(result),
        deck_time_switch_line_numbers(result),
        deck_time_switch_from_node_names(result),
        deck_time_switch_to_node_names(result),
        deck_time_switch_close_time_s_values(result),
        deck_time_switch_open_time_s_values(result),
        deck_time_switch_initially_closed_flags(result),
        deck_time_switch_on_conductance_values(result),
        deck_time_switch_off_conductance_values(result),
        deck_time_switch_count(result),
        deck_control_card_kinds(result),
        deck_control_card_labels(result),
        deck_control_card_line_numbers(result),
        deck_control_card_tokens(result),
        length(result.control_cards),
        copy(result.card_counts),
    )
end

deck_control_card_kinds(result::DeckParseResult) = [card.kind for card in result.control_cards]

deck_control_card_labels(result::DeckParseResult) = [card.label for card in result.control_cards]

deck_control_card_line_numbers(result::DeckParseResult) = Int[card.line_no for card in result.control_cards]

deck_control_card_tokens(result::DeckParseResult) =
    Vector{String}[copy(card.tokens) for card in result.control_cards]

function deck_over5a_source_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.over5a_source_rows)
end

deck_over5a_source_names(result::DeckParseResult) =
    Symbol[row.name for row in deck_over5a_source_rows(result)]

deck_over5a_source_node_names(result::DeckParseResult) =
    Symbol[row.node for row in deck_over5a_source_rows(result)]

deck_over5a_source_node_values(result::DeckParseResult) =
    Int[row.node_value for row in deck_over5a_source_rows(result)]

deck_over5a_source_iform_values(result::DeckParseResult) =
    Int[row.iform for row in deck_over5a_source_rows(result)]

deck_over5a_source_layout_kinds(result::DeckParseResult) =
    Symbol[row.layout_kind for row in deck_over5a_source_rows(result)]

function deck_over5a_source_line_numbers(result::DeckParseResult)
    assert_deck_valid!(result)
    return Int[row.line_no for row in deck_over5a_source_rows(result)]
end

function deck_over5a_source_crest_values(rows::Vector{DeckOVER5ASourceRow})
    values = [row.crest for row in rows]
    required_length = length(rows)
    for (index, row) in enumerate(rows)
        if row.iform == 16
            required_length = max(required_length, index + 2)
        elseif row.iform == 17
            required_length = max(required_length, index + 1)
        end
    end
    while length(values) < required_length
        push!(values, 0.0)
    end
    return values
end

deck_over5a_source_crest_values(result::DeckParseResult) =
    deck_over5a_source_crest_values(deck_over5a_source_rows(result))

deck_over5a_source_time1_values(result::DeckParseResult) =
    Float64[row.time1 for row in deck_over5a_source_rows(result)]

deck_over5a_source_time2_values(result::DeckParseResult) =
    Float64[row.time2 for row in deck_over5a_source_rows(result)]

deck_over5a_source_sfreq_values(result::DeckParseResult) =
    Float64[row.sfreq for row in deck_over5a_source_rows(result)]

deck_over5a_source_tstart_values(result::DeckParseResult) =
    Float64[row.tstart for row in deck_over5a_source_rows(result)]

deck_over5a_source_tstop_values(result::DeckParseResult) =
    Float64[row.tstop for row in deck_over5a_source_rows(result)]

function deck_over5a_source_update_inputs(result::DeckParseResult)
    rows = deck_over5a_source_rows(result)
    return DeckOVER5ASourceUpdateInputs(
        [row.name for row in rows],
        [row.node for row in rows],
        [row.node_value for row in rows],
        [row.iform for row in rows],
        [row.line_no for row in rows],
        [row.layout_kind for row in rows],
        [row.tstart for row in rows],
        [row.tstop for row in rows],
        deck_over5a_source_crest_values(rows),
        [row.time1 for row in rows],
        [row.time2 for row in rows],
        [row.sfreq for row in rows],
        length(rows),
    )
end

function deck_over5_switch_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.over5_switch_rows)
end

deck_over5_switch_names(result::DeckParseResult) =
    Symbol[row.name for row in deck_over5_switch_rows(result)]

deck_over5_switch_line_numbers(result::DeckParseResult) =
    Int[row.line_no for row in deck_over5_switch_rows(result)]

deck_over5_switch_from_node_names(result::DeckParseResult) =
    Symbol[row.from_node for row in deck_over5_switch_rows(result)]

deck_over5_switch_to_node_names(result::DeckParseResult) =
    Symbol[row.to_node for row in deck_over5_switch_rows(result)]

deck_over5_switch_from_node_indices(result::DeckParseResult) =
    Int[row.from_node_value for row in deck_over5_switch_rows(result)]

deck_over5_switch_to_node_indices(result::DeckParseResult) =
    Int[row.to_node_value for row in deck_over5_switch_rows(result)]

deck_over5_switch_layout_kinds(result::DeckParseResult) =
    Symbol[row.layout_kind for row in deck_over5_switch_rows(result)]

deck_over5_switch_raw_close_time_s_values(result::DeckParseResult) =
    Float64[row.raw_close_time_s for row in deck_over5_switch_rows(result)]

deck_over5_switch_raw_open_time_s_values(result::DeckParseResult) =
    Float64[row.raw_open_time_s for row in deck_over5_switch_rows(result)]

deck_over5_switch_close_time_s_values(result::DeckParseResult) =
    Float64[row.close_time_s for row in deck_over5_switch_rows(result)]

deck_over5_switch_open_time_s_values(result::DeckParseResult) =
    Float64[row.open_time_s for row in deck_over5_switch_rows(result)]

deck_over5_switch_initially_closed_flags(result::DeckParseResult) =
    Bool[row.initially_closed for row in deck_over5_switch_rows(result)]

deck_over5_switch_measuring_flags(result::DeckParseResult) =
    Bool[row.measuring for row in deck_over5_switch_rows(result)]

deck_over5_switch_closed_markers(result::DeckParseResult) =
    String[row.closed_marker for row in deck_over5_switch_rows(result)]

deck_over5_switch_marker_texts(result::DeckParseResult) =
    String[row.marker_text for row in deck_over5_switch_rows(result)]

deck_over5_switch_on_conductance_values(result::DeckParseResult) =
    Float64[row.on_conductance for row in deck_over5_switch_rows(result)]

deck_over5_switch_off_conductance_values(result::DeckParseResult) =
    Float64[row.off_conductance for row in deck_over5_switch_rows(result)]

deck_over5_switch_output_codes(result::DeckParseResult) =
    Int[row.output_code for row in deck_over5_switch_rows(result)]

function deck_over15_output_request_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.over15_output_request_rows)
end

deck_over15_output_request_names(result::DeckParseResult) =
    Symbol[row.name for row in deck_over15_output_request_rows(result)]

deck_over15_output_request_output_kinds(result::DeckParseResult) =
    Symbol[row.output_kind for row in deck_over15_output_request_rows(result)]

deck_over15_output_request_request_kinds(result::DeckParseResult) =
    Symbol[row.request_kind for row in deck_over15_output_request_rows(result)]

deck_over15_output_request_layout_kinds(result::DeckParseResult) =
    Symbol[row.layout_kind for row in deck_over15_output_request_rows(result)]

deck_over15_output_request_line_numbers(result::DeckParseResult) =
    Int[row.line_no for row in deck_over15_output_request_rows(result)]

deck_over15_output_request_output_codes(result::DeckParseResult) =
    Int[row.output_code for row in deck_over15_output_request_rows(result)]

deck_over15_output_request_node_names(result::DeckParseResult) =
    Symbol[row.node for row in deck_over15_output_request_rows(result)]

deck_over15_output_request_node_indices(result::DeckParseResult) =
    Int[row.node_value for row in deck_over15_output_request_rows(result)]

deck_over15_output_request_branch_names(result::DeckParseResult) =
    Symbol[row.branch for row in deck_over15_output_request_rows(result)]

deck_over15_output_request_branch_indices(result::DeckParseResult) =
    Int[row.branch_value for row in deck_over15_output_request_rows(result)]

function deck_time_switch_indices(result::DeckParseResult)
    assert_deck_valid!(result)
    return Int[index for (index, element) in enumerate(result.elements) if element isa TimeSwitch]
end

deck_time_switch_count(result::DeckParseResult) = length(deck_time_switch_indices(result))

deck_time_switch_names(result::DeckParseResult) =
    Symbol[result.element_names[index] for index in deck_time_switch_indices(result)]

deck_time_switch_line_numbers(result::DeckParseResult) =
    Int[result.element_line_numbers[index] for index in deck_time_switch_indices(result)]

deck_time_switch_from_node_indices(result::DeckParseResult) =
    Int[result.elements[index].a for index in deck_time_switch_indices(result)]

deck_time_switch_to_node_indices(result::DeckParseResult) =
    Int[result.elements[index].b for index in deck_time_switch_indices(result)]

function _deck_node_names_for_indices(result::DeckParseResult, indices::Vector{Int})
    node_names = deck_node_names(result)
    return Symbol[index == 0 ? :ground : node_names[index] for index in indices]
end

deck_branch_from_node_names(result::DeckParseResult) =
    _deck_node_names_for_indices(result, deck_branch_from_node_indices(result))

deck_branch_to_node_names(result::DeckParseResult) =
    _deck_node_names_for_indices(result, deck_branch_to_node_indices(result))

deck_time_switch_from_node_names(result::DeckParseResult) =
    _deck_node_names_for_indices(result, deck_time_switch_from_node_indices(result))

deck_time_switch_to_node_names(result::DeckParseResult) =
    _deck_node_names_for_indices(result, deck_time_switch_to_node_indices(result))

deck_time_switch_close_time_s_values(result::DeckParseResult) =
    Float64[Float64(result.elements[index].close_time_s) for index in deck_time_switch_indices(result)]

deck_time_switch_open_time_s_values(result::DeckParseResult) =
    Float64[Float64(result.elements[index].open_time_s) for index in deck_time_switch_indices(result)]

deck_time_switch_initially_closed_flags(result::DeckParseResult) =
    Bool[result.elements[index].initially_closed for index in deck_time_switch_indices(result)]

deck_time_switch_on_conductance_values(result::DeckParseResult) =
    Float64[Float64(result.elements[index].on_conductance) for index in deck_time_switch_indices(result)]

deck_time_switch_off_conductance_values(result::DeckParseResult) =
    Float64[Float64(result.elements[index].off_conductance) for index in deck_time_switch_indices(result)]

deck_over16_output_channel_names(result::DeckParseResult) =
    [channel.name for channel in result.over16_output_channels]

deck_over16_output_node_names(result::DeckParseResult) =
    [channel.node for channel in result.over16_output_channels]

deck_over16_output_channel_line_numbers(result::DeckParseResult) =
    Int[channel.line_no for channel in result.over16_output_channels]

function deck_over16_output_node_indices(result::DeckParseResult)
    assert_deck_valid!(result)
    return [result.node_map[channel.node] for channel in result.over16_output_channels]
end

deck_over16_branch_voltage_output_names(result::DeckParseResult) =
    [request.name for request in result.over16_branch_voltage_outputs]

deck_over16_branch_voltage_branch_names(result::DeckParseResult) =
    [request.branch for request in result.over16_branch_voltage_outputs]

deck_over16_branch_voltage_output_line_numbers(result::DeckParseResult) =
    Int[request.line_no for request in result.over16_branch_voltage_outputs]

function deck_over16_branch_voltage_branch_indices(result::DeckParseResult)
    assert_deck_valid!(result)
    return [findfirst(==(request.branch), result.element_names)::Int
            for request in result.over16_branch_voltage_outputs]
end

deck_over16_branch_current_output_names(result::DeckParseResult) =
    [request.name for request in result.over16_branch_current_outputs]

deck_over16_branch_current_branch_names(result::DeckParseResult) =
    [request.branch for request in result.over16_branch_current_outputs]

deck_over16_branch_current_output_line_numbers(result::DeckParseResult) =
    Int[request.line_no for request in result.over16_branch_current_outputs]

function deck_over16_branch_current_branch_indices(result::DeckParseResult)
    assert_deck_valid!(result)
    return [findfirst(==(request.branch), result.element_names)::Int
            for request in result.over16_branch_current_outputs]
end

deck_over16_branch_power_output_names(result::DeckParseResult) =
    [request.name for request in result.over16_branch_power_outputs]

deck_over16_branch_power_branch_names(result::DeckParseResult) =
    [request.branch for request in result.over16_branch_power_outputs]

deck_over16_branch_power_output_line_numbers(result::DeckParseResult) =
    Int[request.line_no for request in result.over16_branch_power_outputs]

function deck_over16_branch_power_branch_indices(result::DeckParseResult)
    assert_deck_valid!(result)
    return [findfirst(==(request.branch), result.element_names)::Int
            for request in result.over16_branch_power_outputs]
end

function deck_over16_source_card_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return [
        DeckOVER16SourceCardRow(row.kind, copy(row.values), row.provided_value_count, row.line_no)
        for row in result.over16_source_card_rows
    ]
end

deck_over16_source_card_kinds(result::DeckParseResult) =
    [row.kind for row in result.over16_source_card_rows]

function deck_over16_source_card_values(result::DeckParseResult)
    assert_deck_valid!(result)
    return [copy(row.values) for row in result.over16_source_card_rows]
end

deck_over16_source_card_provided_value_counts(result::DeckParseResult) =
    Int[row.provided_value_count for row in result.over16_source_card_rows]

deck_over16_source_card_line_numbers(result::DeckParseResult) =
    Int[row.line_no for row in result.over16_source_card_rows]

function deck_over16_source_interpolation_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return [
        DeckOVER16SourceInterpolationRow(copy(row.values), row.provided_value_count, row.line_no)
        for row in result.over16_source_interpolation_rows
    ]
end

function deck_over16_source_interpolation_values(result::DeckParseResult)
    assert_deck_valid!(result)
    return [copy(row.values) for row in result.over16_source_interpolation_rows]
end

deck_over16_source_interpolation_provided_value_counts(result::DeckParseResult) =
    Int[row.provided_value_count for row in result.over16_source_interpolation_rows]

deck_over16_source_interpolation_line_numbers(result::DeckParseResult) =
    Int[row.line_no for row in result.over16_source_interpolation_rows]

function deck_over16_source_tacs_override_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return [
        DeckOVER16SourceTACSOverrideRow(row.position, row.xtcs_index, row.line_no)
        for row in result.over16_source_tacs_override_rows
    ]
end

deck_over16_source_tacs_override_positions(result::DeckParseResult) =
    [row.position for row in result.over16_source_tacs_override_rows]

deck_over16_source_tacs_override_xtcs_indices(result::DeckParseResult) =
    [row.xtcs_index for row in result.over16_source_tacs_override_rows]

deck_over16_source_tacs_override_line_numbers(result::DeckParseResult) =
    Int[row.line_no for row in result.over16_source_tacs_override_rows]

function deck_over16_source_analytic_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return [
        DeckOVER16SourceAnalyticRow(copy(row.values), row.provided_value_count, row.line_no)
        for row in result.over16_source_analytic_rows
    ]
end

function deck_over16_source_analytic_values(result::DeckParseResult)
    assert_deck_valid!(result)
    return [copy(row.values) for row in result.over16_source_analytic_rows]
end

deck_over16_source_analytic_provided_value_counts(result::DeckParseResult) =
    Int[row.provided_value_count for row in result.over16_source_analytic_rows]

deck_over16_source_analytic_line_numbers(result::DeckParseResult) =
    Int[row.line_no for row in result.over16_source_analytic_rows]

function deck_element_kind(element)
    if element isa ConductanceBranch
        return :conductance
    elseif element isa SeriesRLBranch
        return :series_rl
    elseif element isa SeriesRLCBranch
        return :series_rlc
    elseif element isa CapacitorBranch
        return :capacitor
    elseif element isa CurrentInjection
        return :current
    elseif element isa TheveninSource
        return :source
    elseif element isa IdealSwitch
        return :switch
    elseif element isa TimeSwitch
        return :time_switch
    elseif element isa BergeronLine
        return :bergeron_line
    elseif element isa DistributedTransposedLineCompanionAdmittance
        return :distributed_transposed_line_companion_admittance
    elseif element isa TerminalSurgeImpedanceAdmittance
        return :terminal_surge_impedance_admittance
    end
    return Symbol(nameof(typeof(element)))
end

function deck_asset_tables(result::DeckParseResult)
    assert_deck_valid!(result)
    buses = [
        Dict{Symbol,Any}(
            :id => name,
            :index => index,
            :source => result.source,
        )
        for (index, name) in enumerate(deck_node_names(result))
    ]
    elements = [
        Dict{Symbol,Any}(
            :id => result.element_names[index],
            :kind => deck_element_kind(element),
            :index => index,
            :line_no => result.element_line_numbers[index],
            :source => result.source,
        )
        for (index, element) in enumerate(result.elements)
    ]
    bergeron_lines = [
        Dict{Symbol,Any}(
            :id => row.name,
            :kind => :bergeron_line,
            :index => index,
            :line_no => row.line_no,
            :from_node => row.from_node,
            :to_node => row.to_node,
            :from_node_index => row.from_node_value,
            :to_node_index => row.to_node_value,
            :surge_impedance => row.surge_impedance,
            :surge_admittance => row.surge_admittance,
            :travel_time_s => row.travel_time_s,
            :dt_s => row.dt_s,
            :attenuation => row.attenuation,
            :delay_steps => row.delay_steps,
            :write_index => row.write_index,
            :history_current_from => row.history_current_from,
            :history_current_to => row.history_current_to,
            :terminal_voltage_from => row.terminal_voltage_from,
            :terminal_voltage_to => row.terminal_voltage_to,
            :terminal_current_from => row.terminal_current_from,
            :terminal_current_to => row.terminal_current_to,
            :traveling_wave_from => row.traveling_wave_from,
            :traveling_wave_to => row.traveling_wave_to,
            :source => result.source,
        )
        for (index, row) in enumerate(result.bergeron_line_rows)
    ]
    coupled_lumped_sequences = [
        Dict{Symbol,Any}(
            :id => impedance.name,
            :kind => :coupled_lumped_sequence_impedance,
            :index => index,
            :line_numbers => copy(impedance.line_numbers),
            :phase_indices => copy(impedance.phase_indices),
            :from_nodes => copy(impedance.from_nodes),
            :to_nodes => copy(impedance.to_nodes),
            :from_node_indices => copy(impedance.from_node_indices),
            :to_node_indices => copy(impedance.to_node_indices),
            :zero_sequence_resistance => impedance.zero_sequence_resistance,
            :positive_sequence_resistance => impedance.positive_sequence_resistance,
            :zero_sequence_inductance => impedance.zero_sequence_inductance,
            :positive_sequence_inductance => impedance.positive_sequence_inductance,
            :phase_resistance_matrix => copy(impedance.phase_resistance_matrix),
            :phase_inductance_matrix => copy(impedance.phase_inductance_matrix),
            :phase_capacitance_matrix => copy(impedance.phase_capacitance_matrix),
            :stored_resistance_values => copy(impedance.stored_resistance_values),
            :stored_inductance_values => copy(impedance.stored_inductance_values),
            :stored_capacitance_values => copy(impedance.stored_capacitance_values),
            :source => result.source,
        )
        for (index, impedance) in enumerate(deck_coupled_lumped_sequence_impedances(result))
    ]
    distributed_transposed_lines = [
        Dict{Symbol,Any}(
            :id => constants.name,
            :kind => :distributed_transposed_line_constants,
            :index => index,
            :line_numbers => copy(constants.line_numbers),
            :phase_indices => copy(constants.phase_indices),
            :from_nodes => copy(constants.from_nodes),
            :to_nodes => copy(constants.to_nodes),
            :from_node_indices => copy(constants.from_node_indices),
            :to_node_indices => copy(constants.to_node_indices),
            :line_length => constants.line_length,
            :x_frequency_hz => constants.x_frequency_hz,
            :c_frequency_hz => constants.c_frequency_hz,
            :sequence_resistance_per_length => copy(constants.sequence_resistance_per_length),
            :sequence_inductance_input_values => copy(constants.sequence_inductance_input_values),
            :sequence_capacitance_input_values => copy(constants.sequence_capacitance_input_values),
            :sequence_inductance_h_per_length => copy(constants.sequence_inductance_h_per_length),
            :sequence_capacitance_f_per_length => copy(constants.sequence_capacitance_f_per_length),
            :sequence_characteristic_impedances => copy(constants.sequence_characteristic_impedances),
            :sequence_total_resistances => copy(constants.sequence_total_resistances),
            :sequence_propagation_times_s => copy(constants.sequence_propagation_times_s),
            :ci_values => copy(constants.ci_values),
            :ck_values => copy(constants.ck_values),
            :cik_values => copy(constants.cik_values),
            :phase_resistance_matrix => copy(constants.phase_resistance_matrix),
            :phase_inductance_matrix => copy(constants.phase_inductance_matrix),
            :phase_capacitance_matrix => copy(constants.phase_capacitance_matrix),
            :stored_resistance_values => copy(constants.stored_resistance_values),
            :stored_inductance_values => copy(constants.stored_inductance_values),
            :stored_capacitance_values => copy(constants.stored_capacitance_values),
            :source => result.source,
        )
        for (index, constants) in enumerate(deck_distributed_transposed_line_constants(result))
    ]
    line_constants_conductor_cards = [
        Dict{Symbol,Any}(
            :id => Symbol("line_constants_conductor_card_", index),
            :kind => :line_constants_conductor_card,
            :index => index,
            :line_no => row.line_no,
            :phase_number => row.phase_number,
            :skin_effect_type => row.skin_effect_type,
            :resistance_ohm_per_mile => row.resistance_ohm_per_mile,
            :reactance_type => row.reactance_type,
            :reactance_or_gmr => row.reactance_or_gmr,
            :diameter_inches => row.diameter_inches,
            :horizontal_ft => row.horizontal_ft,
            :tower_height_ft => row.tower_height_ft,
            :midspan_height_ft => row.midspan_height_ft,
            :average_height_ft => row.average_height_ft,
            :bundle_spacing_inches => row.bundle_spacing_inches,
            :bundle_angle_deg => row.bundle_angle_deg,
            :conductor_name => row.conductor_name,
            :bundle_conductor_count => row.bundle_conductor_count,
            :source => result.source,
        )
        for (index, row) in enumerate(result.line_constants_conductor_cards)
    ]
    line_constants_physical_conductors = [
        Dict{Symbol,Any}(
            :id => Symbol("line_constants_physical_conductor_", index),
            :kind => :line_constants_physical_conductor,
            :index => index,
            :source_card_index => row.source_card_index,
            :line_no => row.line_no,
            :bundle_ordinal => row.bundle_ordinal,
            :phase_number => row.phase_number,
            :skin_effect_type => row.skin_effect_type,
            :resistance_ohm_per_mile => row.resistance_ohm_per_mile,
            :reactance_type => row.reactance_type,
            :reactance_or_gmr => row.reactance_or_gmr,
            :diameter_inches => row.diameter_inches,
            :horizontal_ft => row.horizontal_ft,
            :average_height_ft => row.average_height_ft,
            :conductor_name => row.conductor_name,
            :source => result.source,
        )
        for (index, row) in enumerate(deck_line_constants_physical_conductors(result))
    ]
    line_constants_frequency_cards = [
        Dict{Symbol,Any}(
            :id => Symbol("line_constants_frequency_card_", index),
            :kind => :line_constants_frequency_card,
            :index => index,
            :line_no => row.line_no,
            :earth_resistivity_ohm_m => row.earth_resistivity_ohm_m,
            :frequency_hz => row.frequency_hz,
            :carson_correction_factor => row.carson_correction_factor,
            :capacitance_print_flags => row.capacitance_print_flags,
            :impedance_print_flags => row.impedance_print_flags,
            :matrix_output_selector => row.matrix_output_selector,
            :distance_miles => row.distance_miles,
            :punch_request => row.punch_request,
            :alternate_punch_flags => row.alternate_punch_flags,
            :frequency_decade_count => row.frequency_decade_count,
            :points_per_decade => row.points_per_decade,
            :line_model_punch_request => row.line_model_punch_request,
            :modal_output_flag => row.modal_output_flag,
            :transform_output_flag => row.transform_output_flag,
            :conductance_mho_per_mile => row.conductance_mho_per_mile,
            :source => result.source,
        )
        for (index, row) in enumerate(result.line_constants_frequency_cards)
    ]
    cable_constants_cases = [
        Dict{Symbol,Any}(
            :id => Symbol("cable_constants_case_", index),
            :kind => :cable_constants_case,
            :index => index,
            :line_no => row.line_no,
            :cable_kind_code => row.cable_kind_code,
            :surface_position_code => row.surface_position_code,
            :phase_count => row.phase_count,
            :earth_model_code => row.earth_model_code,
            :modal_output_flag => row.modal_output_flag,
            :impedance_output_flag => row.impedance_output_flag,
            :admittance_output_flag => row.admittance_output_flag,
            :pipe_count => row.pipe_count,
            :grounding_selector => row.grounding_selector,
            :layer_counts => copy(row.layer_counts),
            :boundary_radii_m => copy(row.boundary_radii_m),
            :resistivity_ohm_m => copy(row.resistivity_ohm_m),
            :conductor_relative_permeability => copy(row.conductor_relative_permeability),
            :insulation_relative_permeability => copy(row.insulation_relative_permeability),
            :insulation_relative_permittivity => copy(row.insulation_relative_permittivity),
            :depths_m => copy(row.depths_m),
            :horizontal_positions_m => copy(row.horizontal_positions_m),
            :frequency_card_count => length(row.frequency_cards),
            :source => result.source,
        )
        for (index, row) in enumerate(result.cable_constants_cases)
    ]
    power_frequency_requests = [
        Dict{Symbol,Any}(
            :id => Symbol("power_frequency_request_", index),
            :kind => :power_frequency_request,
            :index => index,
            :line_no => row.line_no,
            :frequency_hz => row.frequency_hz,
            :source => result.source,
        )
        for (index, row) in enumerate(result.power_frequency_request_rows)
    ]
    universal_machine_dimension_requests = [
        Dict{Symbol,Any}(
            :id => Symbol("universal_machine_dimension_request_", index),
            :kind => :universal_machine_dimension_request,
            :index => index,
            :line_no => row.line_no,
            :coil_table_size => row.coil_table_size,
            :machine_table_size => row.machine_table_size,
            :output_table_size => row.output_table_size,
            :bus_table_size => row.bus_table_size,
            :source => result.source,
        )
        for (index, row) in enumerate(result.universal_machine_dimension_request_rows)
    ]
    universal_machine_sections = [
        Dict{Symbol,Any}(
            :id => Symbol("universal_machine_section_", index),
            :kind => :universal_machine_section,
            :index => index,
            :line_no => row.line_no,
            :machine_count => row.machine_count,
            :source => result.source,
        )
        for (index, row) in enumerate(result.universal_machine_section_rows)
    ]
    universal_machine_definitions = [
        Dict{Symbol,Any}(
            :id => Symbol("universal_machine_definition_", index),
            :kind => :universal_machine_definition,
            :index => index,
            :line_no => row.line_no,
            :machine_index => row.machine_index,
            :card_index => row.card_index,
            :machine_type => row.machine_type,
            :value1 => row.value1,
            :value2 => row.value2,
            :mechanical_damping_coefficient => row.mechanical_damping_coefficient,
            :speed_convergence_margin => row.speed_convergence_margin,
            :node => row.node,
            :source => result.source,
        )
        for (index, row) in enumerate(result.universal_machine_definition_rows)
    ]
    universal_machine_coils = [
        Dict{Symbol,Any}(
            :id => Symbol("universal_machine_coil_", index),
            :kind => :universal_machine_coil,
            :index => index,
            :line_no => row.line_no,
            :machine_index => row.machine_index,
            :coil_index => row.coil_index,
            :resistance => row.resistance,
            :inductance => row.inductance,
            :terminal_node => row.terminal_node,
            :output_flag => row.output_flag,
            :source => result.source,
        )
        for (index, row) in enumerate(result.universal_machine_coil_rows)
    ]
    universal_machine_terminals = [
        Dict{Symbol,Any}(
            :id => Symbol("universal_machine_terminal_", index),
            :kind => :universal_machine_terminal,
            :index => index,
            :line_no => row.line_no,
            :machine_index => row.machine_index,
            :terminal_index => row.terminal_index,
            :terminal_node => row.terminal_node,
            :reference_node => row.reference_node,
            :terminal_node_value => row.terminal_node_value,
            :reference_node_value => row.reference_node_value,
            :source => result.source,
        )
        for (index, row) in enumerate(result.universal_machine_terminal_rows)
    ]
    universal_machine_generated_branches = [
        Dict{Symbol,Any}(
            :id => Symbol("universal_machine_generated_branch_", index),
            :kind => :universal_machine_generated_branch,
            :index => index,
            :line_no => row.line_no,
            :machine_index => row.machine_index,
            :branch_index => row.branch_index,
            :from_node => row.from_node,
            :to_node => row.to_node,
            :from_node_value => row.from_node_value,
            :to_node_value => row.to_node_value,
            :reactance => row.reactance,
            :source => result.source,
        )
        for (index, row) in enumerate(result.universal_machine_generated_branch_rows)
    ]
    universal_machine_speed_capacitors = [
        Dict{Symbol,Any}(
            :id => Symbol("universal_machine_speed_capacitor_", index),
            :kind => :universal_machine_speed_capacitor,
            :index => index,
            :line_no => row.line_no,
            :machine_index => row.machine_index,
            :capacitor_node => row.capacitor_node,
            :mass_node => row.mass_node,
            :capacitor_node_value => row.capacitor_node_value,
            :mass_node_value => row.mass_node_value,
            :resistance => row.resistance,
            :capacitance => row.capacitance,
            :source => result.source,
        )
        for (index, row) in enumerate(result.universal_machine_speed_capacitor_rows)
    ]
    universal_machine_node_summaries = [
        Dict{Symbol,Any}(
            :id => Symbol("universal_machine_node_summary_", index),
            :kind => :universal_machine_node_summary,
            :index => index,
            :line_no => row.line_no,
            :machine_index => row.machine_index,
            :mass_node => row.mass_node,
            :mass_node_value => row.mass_node_value,
            :mechanical_slack_source => row.mechanical_slack_source,
            :mechanical_slack_source_index => row.mechanical_slack_source_index,
            :field_slack_source => row.field_slack_source,
            :field_slack_source_index => row.field_slack_source_index,
            :source => result.source,
        )
        for (index, row) in enumerate(result.universal_machine_node_summary_rows)
    ]
    universal_machine_output_summaries = [
        Dict{Symbol,Any}(
            :id => Symbol("universal_machine_output_summary_", index),
            :kind => :universal_machine_output_summary,
            :index => index,
            :line_no => row.line_no,
            :machine_count => row.machine_count,
            :output_count => row.output_count,
            :tacs_transfer_count => row.tacs_transfer_count,
            :source => result.source,
        )
        for (index, row) in enumerate(result.universal_machine_output_summary_rows)
    ]
    synchronous_machine_terminal_voltages = [
        Dict{Symbol,Any}(
            :id => Symbol("synchronous_machine_terminal_voltage_", index),
            :kind => :synchronous_machine_terminal_voltage,
            :index => index,
            :line_no => row.line_no,
            :machine_index => row.machine_index,
            :source_type => row.source_type,
            :phase_index => row.phase_index,
            :terminal_node => row.terminal_node,
            :terminal_node_value => row.terminal_node_value,
            :peak_terminal_voltage => row.peak_terminal_voltage,
            :frequency_hz => row.frequency_hz,
            :angle_deg => row.angle_deg,
            :source => result.source,
        )
        for (index, row) in enumerate(result.synchronous_machine_terminal_voltage_rows)
    ]
    synchronous_machine_tolerances = [
        Dict{Symbol,Any}(
            :id => Symbol("synchronous_machine_tolerance_", index),
            :kind => :synchronous_machine_tolerance,
            :index => index,
            :line_no => row.line_no,
            :machine_index => row.machine_index,
            :values => copy(row.values),
            :source => result.source,
        )
        for (index, row) in enumerate(result.synchronous_machine_tolerance_rows)
    ]
    synchronous_machine_parameter_fittings = [
        Dict{Symbol,Any}(
            :id => Symbol("synchronous_machine_parameter_fitting_", index),
            :kind => :synchronous_machine_parameter_fitting,
            :index => index,
            :line_no => row.line_no,
            :machine_index => row.machine_index,
            :value => row.value,
            :source => result.source,
        )
        for (index, row) in enumerate(result.synchronous_machine_parameter_fitting_rows)
    ]
    synchronous_machine_model_parameters = [
        Dict{Symbol,Any}(
            :id => Symbol("synchronous_machine_model_parameter_", index),
            :kind => row.parameter_kind,
            :index => index,
            :line_no => row.line_no,
            :machine_index => row.machine_index,
            :values => copy(row.values),
            :positional_values => copy(row.positional_values),
            :source => result.source,
        )
        for (index, row) in enumerate(result.synchronous_machine_model_parameter_rows)
    ]
    synchronous_machine_masses = [
        Dict{Symbol,Any}(
            :id => Symbol("synchronous_machine_mass_", index),
            :kind => :synchronous_machine_mass,
            :index => index,
            :line_no => row.line_no,
            :machine_index => row.machine_index,
            :mass_index => row.mass_index,
            :values => copy(row.values),
            :torque_fraction => row.torque_fraction,
            :inertia => row.inertia,
            :speed_deviation_damping => row.speed_deviation_damping,
            :mutual_damping => row.mutual_damping,
            :shaft_stiffness => row.shaft_stiffness,
            :absolute_speed_damping => row.absolute_speed_damping,
            :source => result.source,
        )
        for (index, row) in enumerate(result.synchronous_machine_mass_rows)
    ]
    synchronous_machine_output_requests = [
        Dict{Symbol,Any}(
            :id => Symbol("synchronous_machine_output_request_", index),
            :kind => :synchronous_machine_output_request,
            :index => index,
            :line_no => row.line_no,
            :machine_index => row.machine_index,
            :group_index => row.group_index,
            :output_codes => copy(row.output_codes),
            :dynamic_output_count => row.dynamic_output_count,
            :source => result.source,
        )
        for (index, row) in enumerate(result.synchronous_machine_output_request_rows)
    ]
    synchronous_machine_control_interfaces = [
        Dict{Symbol,Any}(
            :id => Symbol("synchronous_machine_control_interface_", index),
            :kind => row.coupling_kind,
            :index => index,
            :line_no => row.line_no,
            :machine_index => row.machine_index,
            :interface_code => row.interface_code,
            :direction => row.direction,
            :signal_name => row.signal_name,
            :variable_index => row.variable_index,
            :source => result.source,
        )
        for (index, row) in enumerate(result.synchronous_machine_control_interface_rows)
    ]
    synchronous_machine_output_summaries = [
        Dict{Symbol,Any}(
            :id => Symbol("synchronous_machine_output_summary_", index),
            :kind => :synchronous_machine_output_summary,
            :index => index,
            :line_no => row.line_no,
            :machine_count => row.machine_count,
            :terminal_count => row.terminal_count,
            :mass_count => row.mass_count,
            :output_count => row.output_count,
            :source => result.source,
        )
        for (index, row) in enumerate(result.synchronous_machine_output_summary_rows)
    ]
    tacs_dimension_requests = [
        Dict{Symbol,Any}(
            :id => Symbol("tacs_dimension_request_", index),
            :kind => row.allocation_kind,
            :index => index,
            :line_no => row.line_no,
            :payload_line_no => row.payload_line_no,
            :request_values => collect(row.request_values),
            :list_sizes => collect(row.list_sizes),
            :provided_value_count => row.provided_value_count,
            :source => result.source,
        )
        for (index, row) in enumerate(result.tacs_dimension_request_rows)
    ]
    output_width_requests = [
        Dict{Symbol,Any}(
            :id => Symbol("output_width_request_", index),
            :kind => :output_width_request,
            :index => index,
            :line_no => row.line_no,
            :column_width => row.column_width,
            :source => result.source,
        )
        for (index, row) in enumerate(result.output_width_request_rows)
    ]
    peak_voltage_monitor_requests = [
        Dict{Symbol,Any}(
            :id => Symbol("peak_voltage_monitor_request_", index),
            :kind => :peak_voltage_monitor_request,
            :index => index,
            :line_no => row.line_no,
            :source => result.source,
        )
        for (index, row) in enumerate(result.peak_voltage_monitor_request_rows)
    ]
    diagnostic_print_requests = [
        Dict{Symbol,Any}(
            :id => Symbol("diagnostic_print_request_", index),
            :kind => :diagnostic_print_request,
            :index => index,
            :line_no => row.line_no,
            :loop_print_controls => collect(row.loop_print_controls),
            :source => result.source,
        )
        for (index, row) in enumerate(result.diagnostic_print_request_rows)
    ]
    tacs_warning_limit_requests = [
        Dict{Symbol,Any}(
            :id => Symbol("tacs_warning_limit_request_", index),
            :kind => :tacs_warning_limit_request,
            :index => index,
            :line_no => row.line_no,
            :warning_limit => row.warning_limit,
            :begin_time_s => row.begin_time_s,
            :source => result.source,
        )
        for (index, row) in enumerate(result.tacs_warning_limit_request_rows)
    ]
    plot_file_requests = [
        Dict{Symbol,Any}(
            :id => Symbol("plot_file_request_", index),
            :kind => :plot_file_request,
            :index => index,
            :line_no => row.line_no,
            :plot_file_mode => row.plot_file_mode,
            :source => result.source,
        )
        for (index, row) in enumerate(result.plot_file_request_rows)
    ]
    switch_logic_requests = [
        Dict{Symbol,Any}(
            :id => Symbol("switch_logic_request_", index),
            :kind => :switch_logic_request,
            :index => index,
            :line_no => row.line_no,
            :control_value => row.control_value,
            :source => result.source,
        )
        for (index, row) in enumerate(result.switch_logic_request_rows)
    ]
    simulation_control_requests = [
        Dict{Symbol,Any}(
            :id => Symbol("simulation_control_request_", index),
            :kind => row.request_kind,
            :index => index,
            :line_no => row.line_no,
            :numeric_value => row.numeric_value,
            :integer_value => row.integer_value,
            :source => result.source,
        )
        for (index, row) in enumerate(result.simulation_control_request_rows)
    ]
    printout_frequency_changes = [
        Dict{Symbol,Any}(
            :id => Symbol("printout_frequency_change_", index),
            :kind => :printout_frequency_change,
            :index => index,
            :request_line_no => row.request_line_no,
            :payload_line_no => row.payload_line_no,
            :change_steps => collect(row.change_steps),
            :multipliers => collect(row.multipliers),
            :active_pair_count => row.active_pair_count,
            :provided_value_count => row.provided_value_count,
            :source => result.source,
        )
        for (index, row) in enumerate(result.printout_frequency_change_rows)
    ]
    study_option_requests = [
        Dict{Symbol,Any}(
            :id => Symbol("study_option_request_", index),
            :kind => row.request_kind,
            :index => index,
            :line_no => row.line_no,
            :numeric_values => copy(row.numeric_values),
            :integer_values => copy(row.integer_values),
            :text_values => copy(row.text_values),
            :source => result.source,
        )
        for (index, row) in enumerate(result.study_option_request_rows)
    ]
    time_horizon_controls = [
        Dict{Symbol,Any}(
            :id => Symbol("time_horizon_control_", index),
            :kind => :time_horizon_control,
            :index => index,
            :line_no => row.line_no,
            :time_step_s => row.time_step_s,
            :stop_time_s => row.stop_time_s,
            :inductance_frequency_hz => row.inductance_frequency_hz,
            :capacitance_frequency_hz => row.capacitance_frequency_hz,
            :epsilon => row.epsilon,
            :matrix_tolerance => row.matrix_tolerance,
            :start_time_s => row.start_time_s,
            :source => result.source,
        )
        for (index, row) in enumerate(result.time_horizon_control_rows)
    ]
    output_schedule_controls = [
        Dict{Symbol,Any}(
            :id => Symbol("output_schedule_control_", index),
            :kind => :output_schedule_control,
            :index => index,
            :line_no => row.line_no,
            :print_interval_steps => row.print_interval_steps,
            :plot_interval_steps => row.plot_interval_steps,
            :network_print_enabled => row.network_print_enabled,
            :steady_state_print_enabled => row.steady_state_print_enabled,
            :extrema_print_enabled => row.extrema_print_enabled,
            :terminal_conditions_punch_enabled =>
                row.terminal_conditions_punch_enabled,
            :restart_snapshot_enabled => row.restart_snapshot_enabled,
            :plot_file_retention_mode => row.plot_file_retention_mode,
            :energization_count => row.energization_count,
            :special_request_control_word => row.special_request_control_word,
            :source => result.source,
        )
        for (index, row) in enumerate(result.output_schedule_control_rows)
    ]
    case_boundaries = [
        Dict{Symbol,Any}(
            :id => Symbol("case_boundary_", index),
            :kind => row.boundary_kind,
            :index => index,
            :line_no => row.line_no,
            :source => result.source,
        )
        for (index, row) in enumerate(result.case_boundary_rows)
    ]
    control_system_functions = [
        Dict{Symbol,Any}(
            :id => row.name,
            :kind => :control_system_function,
            :index => index,
            :line_no => row.line_no,
            :input_signal_names => Symbol[term.name for term in row.input_terms],
            :input_signal_polarities => Int[term.polarity for term in row.input_terms],
            :order => row.order,
            :gain => row.gain,
            :source => result.source,
        )
        for (index, row) in enumerate(result.control_system_function_rows)
    ]
    control_system_expressions = [
        Dict{Symbol,Any}(
            :id => row.name,
            :kind => :control_system_expression,
            :index => index,
            :line_no => row.line_no,
            :group_type => row.group_type,
            :expression => row.expression,
            :source => result.source,
        )
        for (index, row) in enumerate(result.control_system_expression_rows)
    ]
    control_system_sources = [
        Dict{Symbol,Any}(
            :id => row.name,
            :kind => :control_system_source,
            :index => index,
            :line_no => row.line_no,
            :source_type => row.source_type,
            :amplitude => row.amplitude,
            :delay_or_time_constant => row.delay_or_time_constant,
            :phase_or_width => row.phase_or_width,
            :numeric_values => copy(row.numeric_values),
            :source => result.source,
        )
        for (index, row) in enumerate(result.control_system_source_rows)
    ]
    control_system_devices = [
        Dict{Symbol,Any}(
            :id => row.name,
            :kind => :control_system_device,
            :index => index,
            :line_no => row.line_no,
            :group_type => row.group_type,
            :device_type => row.device_type,
            :first_input => row.first_input,
            :tail_signal_names => copy(row.tail_signal_names),
            :parameter_values => copy(row.parameter_values),
            :source => result.source,
        )
        for (index, row) in enumerate(result.control_system_device_rows)
    ]
    control_system_output_requests = [
        Dict{Symbol,Any}(
            :id => Symbol("control_system_output_request_", index),
            :kind => :control_system_output_request,
            :index => index,
            :line_no => row.line_no,
            :request_type => row.request_type,
            :all_signals => row.all_signals,
            :signal_names => copy(row.signal_names),
            :source => result.source,
        )
        for (index, row) in enumerate(result.control_system_output_request_rows)
    ]
    control_system_switches = [
        Dict{Symbol,Any}(
            :id => Symbol("control_system_switch_", index),
            :kind => :control_system_switch,
            :index => index,
            :line_no => row.line_no,
            :switch_type => row.switch_type,
            :from_node => row.from_node,
            :to_node => row.to_node,
            :ignition_voltage => row.ignition_voltage,
            :holding_current => row.holding_current,
            :deionization_time_s => row.deionization_time_s,
            :initial_state => row.initial_state,
            :control_signal => row.control_signal,
            :gate_signal => row.gate_signal,
            :clamp_signal => row.clamp_signal,
            :switching_delay_or_model => row.switching_delay_or_model,
            :delayed_arc => row.delayed_arc,
            :parameter_source_kind => row.parameter_source_kind,
            :parameter_reference_index => row.parameter_reference_index,
            :parameter_reference_line_no => row.parameter_reference_line_no,
            :layout_kind => row.layout_kind,
            :event_output_code => row.event_output_code,
            :output_code => row.output_code,
            :source => result.source,
        )
        for (index, row) in enumerate(result.control_system_switch_rows)
    ]
    control_system_switch_couplings = [
        Dict{Symbol,Any}(
            :id => Symbol("control_system_switch_coupling_", index),
            :kind => :control_system_switch_coupling,
            :index => index,
            :line_no => row.line_no,
            :switch_type => row.switch_type,
            :from_node => row.from_node,
            :to_node => row.to_node,
            :from_node_index => row.from_node_index,
            :to_node_index => row.to_node_index,
            :control_signal => row.control_signal,
            :gate_signal => row.gate_signal,
            :clamp_signal => row.clamp_signal,
            :control_output_request_index => row.control_output_request_index,
            :control_output_signal_index => row.control_output_signal_index,
            :control_output_linear_index => row.control_output_linear_index,
            :output_code => row.output_code,
            :initial_state => row.initial_state,
            :source => result.source,
        )
        for (index, row) in enumerate(deck_control_system_switch_coupling_rows(result))
    ]
    node_initial_conditions = [
        Dict{Symbol,Any}(
            :id => Symbol("node_initial_condition_", index),
            :kind => row.condition_kind,
            :index => index,
            :line_no => row.line_no,
            :node => row.node,
            :node_index => row.node_index,
            :reference_node => row.reference_node,
            :reference_node_index => row.reference_node_index,
            :real_value => row.real_value,
            :imaginary_value => row.imaginary_value,
            :source => result.source,
        )
        for (index, row) in enumerate(deck_node_initial_condition_rows(result))
    ]
    distributed_transposed_modal_states = [
        Dict{Symbol,Any}(
            :id => state.name,
            :kind => :distributed_transposed_line_modal_branch_state,
            :index => index,
            :line_numbers => copy(state.line_numbers),
            :phase_indices => copy(state.phase_indices),
            :modal_sequence_indices => copy(state.modal_sequence_indices),
            :continuation_copy_source_indices => copy(state.continuation_copy_source_indices),
            :from_nodes => copy(state.from_nodes),
            :to_nodes => copy(state.to_nodes),
            :from_node_indices => copy(state.from_node_indices),
            :to_node_indices => copy(state.to_node_indices),
            :group_length => state.group_length,
            :modal_signed_characteristic_impedances =>
                copy(state.modal_signed_characteristic_impedances),
            :modal_total_resistances => copy(state.modal_total_resistances),
            :modal_propagation_times_s => copy(state.modal_propagation_times_s),
            :source => result.source,
        )
        for (index, state) in
            enumerate(deck_distributed_transposed_line_modal_branch_states(result))
    ]
    distributed_transposed_pi_equivalents = [
        Dict{Symbol,Any}(
            :id => equivalent.name,
            :kind => :distributed_transposed_line_steady_state_pi_equivalent,
            :index => index,
            :line_numbers => copy(equivalent.line_numbers),
            :phase_indices => copy(equivalent.phase_indices),
            :modal_sequence_indices => copy(equivalent.modal_sequence_indices),
            :steady_state_frequency_hz => equivalent.steady_state_frequency_hz,
            :angular_frequency_rad_s => equivalent.angular_frequency_rad_s,
            :storage_start_index => equivalent.storage_start_index,
            :storage_end_index => equivalent.storage_end_index,
            :storage_row_indices => copy(equivalent.storage_row_indices),
            :phase_row_indices => copy(equivalent.phase_row_indices),
            :phase_column_indices => copy(equivalent.phase_column_indices),
            :phase_series_resistance_values =>
                copy(equivalent.phase_series_resistance_values),
            :phase_series_reactance_values =>
                copy(equivalent.phase_series_reactance_values),
            :phase_shunt_conductance_values =>
                copy(equivalent.phase_shunt_conductance_values),
            :phase_shunt_susceptance_values =>
                copy(equivalent.phase_shunt_susceptance_values),
            :source => result.source,
        )
        for (index, equivalent) in
            enumerate(deck_distributed_transposed_line_steady_state_pi_equivalents(result))
    ]
    distributed_transposed_history_states = [
        Dict{Symbol,Any}(
            :id => history.name,
            :kind => :distributed_transposed_line_history_state,
            :index => index,
            :line_numbers => copy(history.line_numbers),
            :phase_indices => copy(history.phase_indices),
            :modal_sequence_indices => copy(history.modal_sequence_indices),
            :timestep_s => history.timestep_s,
            :steady_state_frequency_hz => history.steady_state_frequency_hz,
            :angular_step_rad => history.angular_step_rad,
            :history_storage_start_index => history.history_storage_start_index,
            :next_history_storage_index => history.next_history_storage_index,
            :storage_start_indices => copy(history.storage_start_indices),
            :storage_end_indices => copy(history.storage_end_indices),
            :storage_lengths => copy(history.storage_lengths),
            :history_sample_counts => copy(history.history_sample_counts),
            :initialized_from_steady_state => history.initialized_from_steady_state,
            :source => result.source,
        )
        for (index, history) in enumerate(deck_distributed_transposed_line_history_states(result))
    ]
    distributed_transposed_companion_admittances = [
        Dict{Symbol,Any}(
            :id => admittance.name,
            :kind => :distributed_transposed_line_companion_admittance,
            :index => index,
            :line_numbers => copy(admittance.line_numbers),
            :modal_sequence_indices => copy(admittance.modal_sequence_indices),
            :from_node_indices => copy(admittance.from_node_indices),
            :to_node_indices => copy(admittance.to_node_indices),
            :modal_companion_admittance_values =>
                copy(admittance.modal_companion_admittance_values),
            :phase_companion_admittance_matrix =>
                copy(admittance.phase_companion_admittance_matrix),
            :source => result.source,
        )
        for (index, admittance) in
            enumerate(deck_distributed_transposed_line_companion_admittances(result))
    ]
    branch_indices = _deck_scalar_branch_indices(result)
    branch_elements = [
        Dict{Symbol,Any}(
            :id => result.element_names[element_index],
            :kind => deck_element_kind(result.elements[element_index]),
            :index => index,
            :element_index => element_index,
            :from_node => deck_branch_from_node_names(result)[index],
            :to_node => deck_branch_to_node_names(result)[index],
            :from_node_index => deck_branch_from_node_indices(result)[index],
            :to_node_index => deck_branch_to_node_indices(result)[index],
            :conductance => deck_branch_conductance_values(result)[index],
            :resistance => deck_branch_resistance_values(result)[index],
            :inductance => deck_branch_inductance_values(result)[index],
            :capacitance => deck_branch_capacitance_values(result)[index],
            :previous_current => deck_branch_previous_current_values(result)[index],
            :previous_voltage => deck_branch_previous_voltage_values(result)[index],
            :line_no => deck_branch_line_numbers(result)[index],
            :source => result.source,
        )
        for (index, element_index) in enumerate(branch_indices)
    ]
    over2_branches = [
        Dict{Symbol,Any}(
            :id => row.name,
            :kind => :over2_branch,
            :index => index,
            :line_no => row.line_no,
            :branch_type => row.branch_type,
            :branch_kind => row.branch_kind,
            :layout_kind => row.layout_kind,
            :source_kind => row.source_kind,
            :reference_kind => row.reference_kind,
            :reference_name => row.reference_name,
            :reference_line_no => row.reference_line_no,
            :from_node => row.from_node,
            :to_node => row.to_node,
            :from_node_index => row.from_node_value,
            :to_node_index => row.to_node_value,
            :raw_resistance => row.raw_resistance,
            :raw_inductance => row.raw_inductance,
            :raw_capacitance => row.raw_capacitance,
            :conductance => row.conductance,
            :resistance => row.resistance,
            :inductance => row.inductance,
            :capacitance => row.capacitance,
            :output_code => row.output_code,
            :source => result.source,
        )
        for (index, row) in enumerate(result.over2_branch_rows)
    ]
    zinc_oxide_nonlinear_rows = [
        Dict{Symbol,Any}(
            :id => row.name,
            :kind => :zinc_oxide_nonlinear_row,
            :index => index,
            :line_no => row.line_no,
            :nonlinear_type => row.nonlinear_type,
            :from_node => row.from_node,
            :to_node => row.to_node,
            :from_node_index => row.from_node_index,
            :to_node_index => row.to_node_index,
            :raw_resistance => row.raw_resistance,
            :raw_inductance => row.raw_inductance,
            :raw_capacitance_marker => row.raw_capacitance_marker,
            :output_code => row.output_code,
            :first_characteristic_index => row.first_characteristic_index,
            :source_kind => row.source_kind,
            :reference_kind => row.reference_kind,
            :reference_index => row.reference_index,
            :reference_name => row.reference_name,
            :reference_line_no => row.reference_line_no,
            :source => result.source,
        )
        for (index, row) in enumerate(result.zinc_oxide_nonlinear_rows)
    ]
    zinc_oxide_initialization_rows = [
        Dict{Symbol,Any}(
            :id => Symbol("zinc_oxide_initialization_", index),
            :kind => :zinc_oxide_initialization,
            :index => index,
            :nonlinear_row_index => row.nonlinear_row_index,
            :line_no => row.line_no,
            :reference_voltage => row.reference_voltage,
            :gap_voltage => row.gap_voltage,
            :initial_voltage => row.initial_voltage,
            :source => result.source,
        )
        for (index, row) in enumerate(result.zinc_oxide_initialization_rows)
    ]
    zinc_oxide_breakpoints = [
        Dict{Symbol,Any}(
            :id => Symbol("zinc_oxide_breakpoint_", index),
            :kind => :zinc_oxide_breakpoint,
            :index => index,
            :nonlinear_row_index => row.nonlinear_row_index,
            :breakpoint_index => row.breakpoint_index,
            :line_no => row.line_no,
            :current_coefficient => row.current_coefficient,
            :voltage_exponent => row.voltage_exponent,
            :voltage => row.voltage,
            :source => result.source,
        )
        for (index, row) in enumerate(result.zinc_oxide_breakpoint_rows)
    ]
    nonlinear_resistance_rows = [
        Dict{Symbol,Any}(
            :id => row.name,
            :kind => :nonlinear_resistance_row,
            :index => index,
            :line_no => row.line_no,
            :nonlinear_type => row.nonlinear_type,
            :element_kind => row.element_kind,
            :from_node => row.from_node,
            :to_node => row.to_node,
            :from_node_index => row.from_node_index,
            :to_node_index => row.to_node_index,
            :steady_state_reference => row.steady_state_reference,
            :secondary_reference => row.secondary_reference,
            :table_marker => row.table_marker,
            :output_code => row.output_code,
            :first_characteristic_index => row.first_characteristic_index,
            :source_kind => row.source_kind,
            :reference_kind => row.reference_kind,
            :reference_index => row.reference_index,
            :reference_name => row.reference_name,
            :reference_line_no => row.reference_line_no,
            :source => result.source,
        )
        for (index, row) in enumerate(result.nonlinear_resistance_rows)
    ]
    nonlinear_resistance_initializations = [
        Dict{Symbol,Any}(
            :id => Symbol("nonlinear_resistance_initialization_", index),
            :kind => :nonlinear_resistance_initialization,
            :index => index,
            :nonlinear_row_index => row.nonlinear_row_index,
            :line_no => row.line_no,
            :table_voltage_offset => row.table_voltage_offset,
            :gap_voltage => row.gap_voltage,
            :initial_voltage => row.initial_voltage,
            :source => result.source,
        )
        for (index, row) in enumerate(result.nonlinear_resistance_initialization_rows)
    ]
    nonlinear_resistance_points = [
        Dict{Symbol,Any}(
            :id => Symbol("nonlinear_resistance_point_", index),
            :kind => :nonlinear_resistance_point,
            :index => index,
            :nonlinear_row_index => row.nonlinear_row_index,
            :point_index => row.point_index,
            :line_no => row.line_no,
            :ordinate_value => row.ordinate_value,
            :coordinate_value => row.coordinate_value,
            :source => result.source,
        )
        for (index, row) in enumerate(result.nonlinear_resistance_point_rows)
    ]
    triggered_timed_resistance_rows = [
        Dict{Symbol,Any}(
            :id => row.name,
            :kind => :triggered_timed_resistance,
            :index => index,
            :line_no => row.line_no,
            :from_node => row.from_node,
            :to_node => row.to_node,
            :from_node_index => row.from_node_index,
            :to_node_index => row.to_node_index,
            :trigger_voltage_v => row.trigger_voltage_v,
            :arm_time_s => row.arm_time_s,
            :output_code => row.output_code,
            :source => result.source,
        )
        for (index, row) in enumerate(result.triggered_timed_resistance_rows)
    ]
    triggered_timed_resistance_points = [
        Dict{Symbol,Any}(
            :id => Symbol("triggered_timed_resistance_point_", index),
            :kind => :triggered_timed_resistance_point,
            :index => index,
            :resistance_row_index => row.resistance_row_index,
            :point_index => row.point_index,
            :line_no => row.line_no,
            :elapsed_time_s => row.elapsed_time_s,
            :resistance_ohm => row.resistance_ohm,
            :source => result.source,
        )
        for (index, row) in enumerate(result.triggered_timed_resistance_point_rows)
    ]
    switching_nonlinear_resistor_rows = [
        Dict{Symbol,Any}(
            :id => row.name,
            :kind => :switching_nonlinear_resistor,
            :index => index,
            :line_no => row.line_no,
            :from_node => row.from_node,
            :to_node => row.to_node,
            :from_node_index => row.from_node_index,
            :to_node_index => row.to_node_index,
            :turn_on_voltage => row.turn_on_voltage,
            :minimum_on_time_s => row.minimum_on_time_s,
            :activation_segment_count => row.activation_segment_count,
            :turn_off_voltage => row.turn_off_voltage,
            :output_code => row.output_code,
            :single_flash => row.single_flash,
            :source_kind => row.source_kind,
            :reference_kind => row.reference_kind,
            :reference_index => row.reference_index,
            :reference_name => row.reference_name,
            :reference_line_no => row.reference_line_no,
            :source => result.source,
        )
        for (index, row) in enumerate(result.switching_nonlinear_resistor_rows)
    ]
    switching_nonlinear_resistor_points = [
        Dict{Symbol,Any}(
            :id => Symbol("switching_nonlinear_resistor_point_", index),
            :kind => :switching_nonlinear_resistor_point,
            :index => index,
            :resistor_row_index => row.resistor_row_index,
            :point_index => row.point_index,
            :line_no => row.line_no,
            :current_a => row.current_a,
            :voltage_v => row.voltage_v,
            :source => result.source,
        )
        for (index, row) in enumerate(result.switching_nonlinear_resistor_point_rows)
    ]
    arrester_nonlinear_rows = [
        Dict{Symbol,Any}(
            :id => row.name,
            :kind => :arrester_nonlinear_row,
            :index => index,
            :line_no => row.line_no,
            :nonlinear_type => row.nonlinear_type,
            :from_node => row.from_node,
            :to_node => row.to_node,
            :from_node_index => row.from_node_index,
            :to_node_index => row.to_node_index,
            :flashover_voltage => row.flashover_voltage,
            :voltage_division_factor => row.voltage_division_factor,
            :current_division_factor => row.current_division_factor,
            :output_code => row.output_code,
            :first_constant_index => row.first_constant_index,
            :source_kind => row.source_kind,
            :reference_kind => row.reference_kind,
            :reference_index => row.reference_index,
            :reference_name => row.reference_name,
            :reference_line_no => row.reference_line_no,
            :source => result.source,
        )
        for (index, row) in enumerate(result.arrester_nonlinear_rows)
    ]
    arrester_constants = [
        Dict{Symbol,Any}(
            :id => Symbol("arrester_constants_", index),
            :kind => :arrester_constants,
            :index => index,
            :nonlinear_row_index => row.nonlinear_row_index,
            :first_constant_index => row.first_constant_index,
            :line_no => row.line_no,
            :values => copy(row.values),
            :value_count => length(row.values),
            :source => result.source,
        )
        for (index, row) in enumerate(result.arrester_constant_rows)
    ]
    over5a_sources = [
        Dict{Symbol,Any}(
            :id => row.name,
            :kind => :over5a_source,
            :index => index,
            :node => row.node,
            :node_index => row.node_value,
            :iform => row.iform,
            :layout_kind => row.layout_kind,
            :crest => row.crest,
            :time1 => row.time1,
            :time2 => row.time2,
            :sfreq => row.sfreq,
            :tstart => row.tstart,
            :tstop => row.tstop,
            :line_no => row.line_no,
            :source => result.source,
        )
        for (index, row) in enumerate(result.over5a_source_rows)
    ]
    over5_switches = [
        Dict{Symbol,Any}(
            :id => row.name,
            :kind => :over5_switch,
            :index => index,
            :line_no => row.line_no,
            :layout_kind => row.layout_kind,
            :name => row.name,
            :from_node => row.from_node,
            :to_node => row.to_node,
            :from_node_index => row.from_node_value,
            :to_node_index => row.to_node_value,
            :switch_type => row.switch_type,
            :output_code => row.output_code,
            :raw_close_time_s => row.raw_close_time_s,
            :raw_open_time_s => row.raw_open_time_s,
            :close_time_s => row.close_time_s,
            :open_time_s => row.open_time_s,
            :initially_closed => row.initially_closed,
            :measuring => row.measuring,
            :closed_marker => row.closed_marker,
            :marker_text => row.marker_text,
            :on_conductance => row.on_conductance,
            :off_conductance => row.off_conductance,
            :source => result.source,
        )
        for (index, row) in enumerate(result.over5_switch_rows)
    ]
    over15_output_requests = [
        Dict{Symbol,Any}(
            :id => row.name,
            :kind => row.output_kind,
            :index => index,
            :request_kind => row.request_kind,
            :layout_kind => row.layout_kind,
            :line_no => row.line_no,
            :output_code => row.output_code,
            :node => row.node,
            :node_index => row.node_value,
            :branch => row.branch,
            :branch_index => row.branch_value,
            :source => result.source,
        )
        for (index, row) in enumerate(result.over15_output_request_rows)
    ]
    over16_outputs = [
        Dict{Symbol,Any}(
            :id => channel.name,
            :kind => :over16_voltage_output,
            :index => index,
            :node => channel.node,
            :node_index => result.node_map[channel.node],
            :line_no => channel.line_no,
            :source => result.source,
        )
        for (index, channel) in enumerate(result.over16_output_channels)
    ]
    over16_branch_voltage_outputs = [
        Dict{Symbol,Any}(
            :id => request.name,
            :kind => :over16_branch_voltage_output,
            :index => index,
            :branch => request.branch,
            :branch_index => findfirst(==(request.branch), result.element_names),
            :line_no => request.line_no,
            :source => result.source,
        )
        for (index, request) in enumerate(result.over16_branch_voltage_outputs)
    ]
    over16_branch_current_outputs = [
        Dict{Symbol,Any}(
            :id => request.name,
            :kind => :over16_branch_current_output,
            :index => index,
            :branch => request.branch,
            :branch_index => findfirst(==(request.branch), result.element_names),
            :line_no => request.line_no,
            :source => result.source,
        )
        for (index, request) in enumerate(result.over16_branch_current_outputs)
    ]
    over16_branch_power_outputs = [
        Dict{Symbol,Any}(
            :id => request.name,
            :kind => :over16_branch_power_energy_output,
            :index => index,
            :branch => request.branch,
            :branch_index => findfirst(==(request.branch), result.element_names),
            :line_no => request.line_no,
            :source => result.source,
        )
        for (index, request) in enumerate(result.over16_branch_power_outputs)
    ]
    over16_source_cards = [
        Dict{Symbol,Any}(
            :id => Symbol("source_card_", index),
            :kind => row.kind,
            :index => index,
            :values => copy(row.values),
            :value_count => length(row.values),
            :provided_value_count => row.provided_value_count,
            :line_no => row.line_no,
            :source => result.source,
        )
        for (index, row) in enumerate(result.over16_source_card_rows)
    ]
    over16_source_interpolations = [
        Dict{Symbol,Any}(
            :id => Symbol("source_interpolation_", index),
            :kind => :over16_source_interpolation,
            :index => index,
            :values => copy(row.values),
            :value_count => length(row.values),
            :provided_value_count => row.provided_value_count,
            :line_no => row.line_no,
            :source => result.source,
        )
        for (index, row) in enumerate(result.over16_source_interpolation_rows)
    ]
    over16_source_tacs_overrides = [
        Dict{Symbol,Any}(
            :id => Symbol("source_tacs_override_", index),
            :kind => :over16_source_tacs_override,
            :index => index,
            :position => row.position,
            :xtcs_index => row.xtcs_index,
            :line_no => row.line_no,
            :source => result.source,
        )
        for (index, row) in enumerate(result.over16_source_tacs_override_rows)
    ]
    over16_source_analytics = [
        Dict{Symbol,Any}(
            :id => Symbol("source_analytic_", index),
            :kind => :over16_source_analytic,
            :index => index,
            :values => copy(row.values),
            :value_count => length(row.values),
            :provided_value_count => row.provided_value_count,
            :line_no => row.line_no,
            :source => result.source,
        )
        for (index, row) in enumerate(result.over16_source_analytic_rows)
    ]
    time_switch_indices = deck_time_switch_indices(result)
    node_names = deck_node_names(result)
    time_switches = Vector{Dict{Symbol,Any}}()
    for (index, element_index) in enumerate(time_switch_indices)
        switch = result.elements[element_index]
        push!(
            time_switches,
            Dict{Symbol,Any}(
                :id => result.element_names[element_index],
                :kind => :time_switch,
                :index => index,
                :element_index => element_index,
                :line_no => result.element_line_numbers[element_index],
                :from_node => node_names[switch.a],
                :to_node => node_names[switch.b],
                :from_node_index => switch.a,
                :to_node_index => switch.b,
                :close_time_s => switch.close_time_s,
                :open_time_s => switch.open_time_s,
                :initially_closed => switch.initially_closed,
                :on_conductance => switch.on_conductance,
                :off_conductance => switch.off_conductance,
                :source => result.source,
            ),
        )
    end
    deck_control_cards = [
        Dict{Symbol,Any}(
            :id => Symbol("control_", index),
            :kind => card.kind,
            :label => card.label,
            :line_no => card.line_no,
            :tokens => copy(card.tokens),
            :source => result.source,
        )
        for (index, card) in enumerate(result.control_cards)
    ]
    return Dict{Symbol,Vector{Dict{Symbol,Any}}}(
        :buses => buses,
        :emt_elements => elements,
        :bergeron_lines => bergeron_lines,
        :coupled_lumped_sequence_impedances => coupled_lumped_sequences,
        :distributed_transposed_line_constants => distributed_transposed_lines,
        :line_constants_conductor_cards => line_constants_conductor_cards,
        :line_constants_physical_conductors => line_constants_physical_conductors,
        :line_constants_frequency_cards => line_constants_frequency_cards,
        :cable_constants_cases => cable_constants_cases,
        :power_frequency_requests => power_frequency_requests,
        :universal_machine_dimension_requests => universal_machine_dimension_requests,
        :universal_machine_sections => universal_machine_sections,
        :universal_machine_definitions => universal_machine_definitions,
        :universal_machine_coils => universal_machine_coils,
        :universal_machine_terminals => universal_machine_terminals,
        :universal_machine_generated_branches => universal_machine_generated_branches,
        :universal_machine_speed_capacitors => universal_machine_speed_capacitors,
        :universal_machine_node_summaries => universal_machine_node_summaries,
        :universal_machine_output_summaries => universal_machine_output_summaries,
        :synchronous_machine_terminal_voltages =>
            synchronous_machine_terminal_voltages,
        :synchronous_machine_tolerances => synchronous_machine_tolerances,
        :synchronous_machine_parameter_fittings =>
            synchronous_machine_parameter_fittings,
        :synchronous_machine_model_parameters =>
            synchronous_machine_model_parameters,
        :synchronous_machine_masses => synchronous_machine_masses,
        :synchronous_machine_output_requests => synchronous_machine_output_requests,
        :synchronous_machine_control_interfaces => synchronous_machine_control_interfaces,
        :synchronous_machine_output_summaries => synchronous_machine_output_summaries,
        :tacs_dimension_requests => tacs_dimension_requests,
        :output_width_requests => output_width_requests,
        :peak_voltage_monitor_requests => peak_voltage_monitor_requests,
        :diagnostic_print_requests => diagnostic_print_requests,
        :tacs_warning_limit_requests => tacs_warning_limit_requests,
        :plot_file_requests => plot_file_requests,
        :switch_logic_requests => switch_logic_requests,
        :simulation_control_requests => simulation_control_requests,
        :printout_frequency_changes => printout_frequency_changes,
        :study_option_requests => study_option_requests,
        :time_horizon_controls => time_horizon_controls,
        :output_schedule_controls => output_schedule_controls,
        :case_boundaries => case_boundaries,
        :control_system_functions => control_system_functions,
        :control_system_expressions => control_system_expressions,
        :control_system_sources => control_system_sources,
        :control_system_devices => control_system_devices,
        :control_system_output_requests => control_system_output_requests,
        :control_system_switches => control_system_switches,
        :control_system_switch_couplings => control_system_switch_couplings,
        :node_initial_conditions => node_initial_conditions,
        :distributed_transposed_line_modal_branch_states =>
            distributed_transposed_modal_states,
        :distributed_transposed_line_steady_state_pi_equivalents =>
            distributed_transposed_pi_equivalents,
        :distributed_transposed_line_history_states =>
            distributed_transposed_history_states,
        :distributed_transposed_line_companion_admittances =>
            distributed_transposed_companion_admittances,
        :branch_elements => branch_elements,
        :over2_branches => over2_branches,
        :zinc_oxide_nonlinear_rows => zinc_oxide_nonlinear_rows,
        :zinc_oxide_initializations => zinc_oxide_initialization_rows,
        :zinc_oxide_breakpoints => zinc_oxide_breakpoints,
        :nonlinear_resistance_rows => nonlinear_resistance_rows,
        :nonlinear_resistance_initializations => nonlinear_resistance_initializations,
        :nonlinear_resistance_points => nonlinear_resistance_points,
        :triggered_timed_resistance_rows => triggered_timed_resistance_rows,
        :triggered_timed_resistance_points => triggered_timed_resistance_points,
        :switching_nonlinear_resistor_rows => switching_nonlinear_resistor_rows,
        :switching_nonlinear_resistor_points => switching_nonlinear_resistor_points,
        :arrester_nonlinear_rows => arrester_nonlinear_rows,
        :arrester_constants => arrester_constants,
        :over5a_sources => over5a_sources,
        :over5_switches => over5_switches,
        :over15_output_requests => over15_output_requests,
        :over16_outputs => over16_outputs,
        :over16_branch_voltage_outputs => over16_branch_voltage_outputs,
        :over16_branch_current_outputs => over16_branch_current_outputs,
        :over16_branch_power_outputs => over16_branch_power_outputs,
        :over16_source_cards => over16_source_cards,
        :over16_source_interpolations => over16_source_interpolations,
        :over16_source_tacs_overrides => over16_source_tacs_overrides,
        :over16_source_analytics => over16_source_analytics,
        :time_switches => time_switches,
        :deck_control_cards => deck_control_cards,
    )
end

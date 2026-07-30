
function parse_deck_file(path::AbstractString)
    return parse_deck_lines(readlines(path); source=String(path))
end

function parse_deck_lines(lines; source::AbstractString="deck")
    result = DeckParseResult(String(source), Any[], Int[], Dict{Symbol,Int}(), Symbol[],
                             DeckBergeronLineRow[],
                             DeckCoupledLineRow[],
                             DeckLineModalTransformRow[],
                             DeckCoupledPhasePiSectionRow[],
                             DeckCascadedPiRequestRow[],
                             DeckSampledFrequencyLineRow[],
                             DeckSemlyenLineGroupRow[],
                             DeckRationalLineGroupRow[],
                             DeckGeneratorEquivalentRow[],
                             DeckLineConstantsConductorCard[],
                             DeckLineConstantsFrequencyCard[],
                             DeckCableConstantsCase[],
                             DeckPowerFrequencyRequestRow[],
                             DeckUniversalMachineDimensionRequestRow[],
                             DeckUniversalMachineSectionRow[],
                             DeckUniversalMachineDefinitionRow[],
                             DeckUniversalMachineCoilRow[],
                             DeckUniversalMachineTerminalRow[],
                             DeckUniversalMachineGeneratedBranchRow[],
                             DeckUniversalMachineSpeedCapacitorRow[],
                             DeckUniversalMachineNodeSummaryRow[],
                             DeckUniversalMachineOutputSummaryRow[],
                             DeckSynchronousMachineTerminalVoltageRow[],
                             DeckSynchronousMachineToleranceRow[],
                             DeckSynchronousMachineParameterFittingRow[],
                             DeckSynchronousMachineModelParameterRow[],
                             DeckSynchronousMachineMassRow[],
                             DeckSynchronousMachineOutputRequestRow[],
                             DeckSynchronousMachineControlInterfaceRow[],
                             DeckSynchronousMachineOutputSummaryRow[],
                             DeckTACSDimensionRequestRow[],
                             DeckOutputWidthRequestRow[],
                             DeckPeakVoltageMonitorRequestRow[],
                             DeckDiagnosticPrintRequestRow[],
                             DeckTACSWarningLimitRequestRow[],
                             DeckPlotFileRequestRow[],
                             DeckSwitchLogicRequestRow[],
                             DeckSimulationControlRequestRow[],
                             DeckPrintoutFrequencyChangeRow[],
                             DeckStudyOptionRequestRow[],
                             DeckFixedSourceConstraintRow[],
                             DeckFixedSourceControlRow[],
                             DeckTimeHorizonControlRow[],
                             DeckOutputScheduleControlRow[],
                             DeckCaseBoundaryRow[],
                             DeckControlSystemFunctionRow[],
                             DeckControlSystemExpressionRow[],
                             DeckControlSystemSourceRow[],
                             DeckControlSystemDeviceRow[],
                             DeckControlSystemOutputRequestRow[],
                             DeckControlSystemSwitchRow[],
                             DeckNodeInitialConditionRow[],
                             DeckOVER2BranchRow[],
                             DeckZincOxideNonlinearRow[],
                             DeckZincOxideInitializationRow[],
                             DeckZincOxideBreakpointRow[],
                             DeckNonlinearResistanceRow[],
                             DeckNonlinearResistanceInitializationRow[],
                             DeckNonlinearResistancePointRow[],
                             DeckTriggeredTimedResistanceRow[],
                             DeckTriggeredTimedResistancePointRow[],
                             DeckSwitchingNonlinearResistorRow[],
                             DeckSwitchingNonlinearResistorPointRow[],
                             DeckPiecewiseNonlinearInductorRow[],
                             DeckPiecewiseNonlinearInductorPointRow[],
                             DeckHystereticInductorRow[],
                             DeckHystereticInductorPointRow[],
                             DeckArresterNonlinearRow[],
                             DeckArresterConstantRow[],
                             DeckOVER5ASourceRow[],
                             DeckOVER5SwitchRow[],
                             DeckOVER15OutputRequestRow[],
                             DeckOVER16OutputChannel[],
                             DeckOVER16BranchVoltageOutput[],
                             DeckOVER16BranchCurrentOutput[],
                             DeckOVER16BranchPowerOutput[],
                             DeckOVER16SourceCardRow[],
                             DeckOVER16SourceInterpolationRow[],
                             DeckOVER16SourceTACSOverrideRow[],
                             DeckOVER16SourceAnalyticRow[],
                             DeckControlCard[],
                             Dict{Symbol,Int}(),
                             Dict{Symbol,Vector{String}}(),
                             validation_result(source=String(source)))
    active_section = nothing
    fixed_miscellaneous_control_count = 0
    fixed_saturated_transformer_intake_active = false
    pending_tacs_dimension_request = nothing
    pending_printout_frequency_change_request = nothing
    pending_zinc_oxide_table = nothing
    pending_nonlinear_resistance_table = nothing
    pending_triggered_timed_resistance_table = nothing
    pending_switching_nonlinear_resistor_table = nothing
    pending_piecewise_nonlinear_inductor_table = nothing
    pending_hysteretic_inductor_table = nothing
    pending_arrester_constants = nothing
    pending_universal_machine_data = nothing
    pending_synchronous_machine_data = nothing
    pending_cable_constants = nothing
    pending_generator_equivalent = nothing
    pending_dc_simulator_source = nothing
    pending_cascaded_pi = nothing
    pending_sampled_frequency_line = nothing
    pending_semlyen_line = nothing
    pending_rational_frequency_line = nothing
    semlyen_transform_imaginary_parts = true
    fixed_branch_vintage_mode = 0
    fixed_card_kc_lee_transform_active = false
    fixed_source_load_flow_requested = false
    abort_data_case_line = nothing
    for (line_no, line) in deck_logical_records!(result, lines)
        tokens = deck_tokens(line)
        if !isempty(result.control_system_switch_rows)
            delayed_arc_row = last(result.control_system_switch_rows)
            if delayed_arc_row.switching_delay_or_model !== missing &&
               delayed_arc_row.switching_delay_or_model == 7777.0 &&
               delayed_arc_row.delayed_arc === nothing &&
               line_no > delayed_arc_row.line_no
                parse_delayed_arc_switch_continuation!(
                    result,
                    length(result.control_system_switch_rows),
                    line,
                    line_no,
                )
                continue
            end
        end
        if pending_rational_frequency_line !== nothing
            pending_rational_frequency_line = parse_rational_frequency_line_card!(
                result,
                pending_rational_frequency_line,
                line,
                line_no,
            )
            continue
        elseif pending_semlyen_line !== nothing
            pending_semlyen_line = parse_semlyen_line_card!(
                result,
                pending_semlyen_line,
                line,
                line_no,
            )
            continue
        elseif pending_sampled_frequency_line !== nothing
            if pending_sampled_frequency_line.configuration === nothing
                pending_sampled_frequency_line =
                    parse_sampled_frequency_line_configuration!(
                        result,
                        pending_sampled_frequency_line,
                        line,
                        line_no,
                    )
            else
                pending_sampled_frequency_line =
                    parse_sampled_frequency_line_points!(
                        result,
                        pending_sampled_frequency_line,
                        line,
                        line_no,
                    )
            end
            continue
        elseif pending_cascaded_pi !== nothing &&
               pending_cascaded_pi.configuration === nothing &&
               coupled_phase_pi_numeric_continuation_row(fixed_image(line))
            append_coupled_phase_pi_continuation_row!(
                result,
                fixed_image(line),
                line_no,
                length(result.validation.issues),
            )
            continue
        elseif pending_cascaded_pi !== nothing &&
           pending_cascaded_pi.configuration === nothing &&
           pending_cascaded_pi.phase_count > 0 &&
           cascaded_pi_source_section_complete(result, pending_cascaded_pi)
            pending_cascaded_pi = parse_cascaded_pi_configuration!(
                result,
                pending_cascaded_pi,
                line,
                line_no,
            )
            continue
        elseif pending_cascaded_pi !== nothing &&
               pending_cascaded_pi.configuration !== nothing
            pending_cascaded_pi = parse_cascaded_pi_detail!(
                result,
                pending_cascaded_pi,
                line,
                line_no,
            )
            continue
        end
        semlyen_marker = semlyen_transform_imaginary_parts_marker(tokens)
        if semlyen_marker !== nothing
            semlyen_transform_imaginary_parts = semlyen_marker
            record_card!(
                result,
                semlyen_marker ?
                    :semlyen_complex_modal_transform_enabled :
                    :semlyen_real_modal_transform_enabled,
            )
            continue
        end
        if abort_data_case_line !== nothing && !begin_new_data_case_card(tokens)
            record_case_boundary!(result, :aborted_case_discarded_card, tokens, line_no)
            continue
        elseif abort_data_case_line !== nothing
            abort_data_case_line = nothing
        end
        if pending_generator_equivalent !== nothing
            finished = parse_generator_equivalent_card!(
                result,
                pending_generator_equivalent,
                line,
                line_no,
            )
            finished && (pending_generator_equivalent = nothing)
            continue
        end
        if pending_dc_simulator_source !== nothing
            parse_dc_simulator_secondary_card!(
                result,
                pending_dc_simulator_source,
                line,
                line_no,
            )
            pending_dc_simulator_source = nothing
            continue
        end
        if pending_tacs_dimension_request !== nothing
            parse_tacs_dimension_request_payload!(
                result,
                pending_tacs_dimension_request.allocation_kind,
                pending_tacs_dimension_request.line_no,
                pending_tacs_dimension_request.request_text,
                tokens,
                line_no,
            )
            pending_tacs_dimension_request = nothing
            continue
        end
        if pending_printout_frequency_change_request !== nothing
            parse_printout_frequency_change_payload!(
                result,
                pending_printout_frequency_change_request.line_no,
                pending_printout_frequency_change_request.request_text,
                tokens,
                line_no,
            )
            pending_printout_frequency_change_request = nothing
            continue
        end
        if pending_zinc_oxide_table !== nothing
            if zinc_oxide_table_sentinel(line)
                if !pending_zinc_oxide_table.initialization_read
                    add_issue!(
                        result.validation,
                        missing_data(
                            "line $line_no",
                            "expected zinc-oxide initialization row before table sentinel",
                        ),
                    )
                elseif pending_zinc_oxide_table.breakpoint_count == 0
                    add_issue!(
                        result.validation,
                        missing_data(
                            "line $line_no",
                            "expected at least one zinc-oxide breakpoint row before table sentinel",
                        ),
                    )
                else
                    record_card!(result, :fixed_field)
                    record_card!(result, :zinc_oxide_table_end)
                end
                pending_zinc_oxide_table = nothing
                continue
            end
            marker_before_table = fixed_card_section_marker(tokens)
            named_marker_before_table = fixed_card_named_section_marker(tokens)
            begins_new_card =
                !isempty(tokens) && normalized_deck_token(tokens[1]) == "begin"
            if marker_before_table !== nothing || named_marker_before_table !== nothing ||
               begins_new_card
                add_issue!(
                    result.validation,
                    missing_data(
                        "line $line_no",
                        "expected zinc-oxide table sentinel 9999. before the next deck section",
                    ),
                )
                pending_zinc_oxide_table = nothing
            elseif !pending_zinc_oxide_table.initialization_read
                parse_zinc_oxide_initialization_row!(
                    result,
                    line,
                    line_no,
                    pending_zinc_oxide_table.row_index,
                )
                pending_zinc_oxide_table = (
                    row_index = pending_zinc_oxide_table.row_index,
                    initialization_read = true,
                    breakpoint_count = pending_zinc_oxide_table.breakpoint_count,
                )
                continue
            else
                issues_before = length(result.validation.issues)
                parse_zinc_oxide_breakpoint_row!(
                    result,
                    line,
                    line_no,
                    pending_zinc_oxide_table.row_index,
                    pending_zinc_oxide_table.breakpoint_count + 1,
                )
                breakpoint_count =
                    length(result.validation.issues) == issues_before ?
                    pending_zinc_oxide_table.breakpoint_count + 1 :
                    pending_zinc_oxide_table.breakpoint_count
                pending_zinc_oxide_table = (
                    row_index = pending_zinc_oxide_table.row_index,
                    initialization_read = pending_zinc_oxide_table.initialization_read,
                    breakpoint_count = breakpoint_count,
                )
                continue
            end
        end
        if pending_piecewise_nonlinear_inductor_table !== nothing
            if nonlinear_resistance_table_sentinel(line)
                if pending_piecewise_nonlinear_inductor_table.point_count < 2
                    add_issue!(
                        result.validation,
                        missing_data(
                            "line $line_no",
                            "expected at least two nonlinear-inductor current/flux points before table sentinel",
                        ),
                    )
                else
                    record_card!(result, :fixed_field)
                    record_card!(result, :piecewise_nonlinear_inductor_table_end)
                end
                pending_piecewise_nonlinear_inductor_table = nothing
                continue
            end
            marker_before_table = fixed_card_section_marker(tokens)
            named_marker_before_table = fixed_card_named_section_marker(tokens)
            begins_new_card =
                !isempty(tokens) && normalized_deck_token(tokens[1]) == "begin"
            if marker_before_table !== nothing || named_marker_before_table !== nothing ||
               begins_new_card
                add_issue!(
                    result.validation,
                    missing_data(
                        "line $line_no",
                        "expected nonlinear-inductor table sentinel 9999. before the next deck section",
                    ),
                )
                pending_piecewise_nonlinear_inductor_table = nothing
            else
                issues_before = length(result.validation.issues)
                parse_piecewise_nonlinear_inductor_point_row!(
                    result,
                    line,
                    line_no,
                    pending_piecewise_nonlinear_inductor_table.row_index,
                    pending_piecewise_nonlinear_inductor_table.point_count + 1,
                )
                point_count =
                    length(result.validation.issues) == issues_before ?
                    pending_piecewise_nonlinear_inductor_table.point_count + 1 :
                    pending_piecewise_nonlinear_inductor_table.point_count
                pending_piecewise_nonlinear_inductor_table = (
                    row_index = pending_piecewise_nonlinear_inductor_table.row_index,
                    point_count = point_count,
                )
                continue
            end
        end
        if pending_hysteretic_inductor_table !== nothing
            if nonlinear_resistance_table_sentinel(line)
                if pending_hysteretic_inductor_table.point_count < 3
                    add_issue!(
                        result.validation,
                        missing_data(
                            "line $line_no",
                            "expected at least three type-96 current/flux points before table sentinel",
                        ),
                    )
                else
                    record_card!(result, :fixed_field)
                    record_card!(result, :hysteretic_inductor_table_end)
                end
                pending_hysteretic_inductor_table = nothing
                continue
            end
            marker_before_table = fixed_card_section_marker(tokens)
            named_marker_before_table = fixed_card_named_section_marker(tokens)
            begins_new_card =
                !isempty(tokens) && normalized_deck_token(tokens[1]) == "begin"
            if marker_before_table !== nothing || named_marker_before_table !== nothing ||
               begins_new_card
                add_issue!(
                    result.validation,
                    missing_data(
                        "line $line_no",
                        "expected type-96 table sentinel 9999. before the next deck section",
                    ),
                )
                pending_hysteretic_inductor_table = nothing
            else
                issues_before = length(result.validation.issues)
                parse_hysteretic_inductor_point_row!(
                    result,
                    line,
                    line_no,
                    pending_hysteretic_inductor_table.row_index,
                    pending_hysteretic_inductor_table.point_count + 1,
                )
                point_count =
                    length(result.validation.issues) == issues_before ?
                    pending_hysteretic_inductor_table.point_count + 1 :
                    pending_hysteretic_inductor_table.point_count
                pending_hysteretic_inductor_table = (
                    row_index = pending_hysteretic_inductor_table.row_index,
                    point_count = point_count,
                )
                continue
            end
        end
        if pending_nonlinear_resistance_table !== nothing
            if nonlinear_resistance_table_sentinel(line)
                if !pending_nonlinear_resistance_table.initialization_read
                    add_issue!(
                        result.validation,
                        missing_data(
                            "line $line_no",
                            "expected nonlinear resistance initialization row before table sentinel",
                        ),
                    )
                elseif pending_nonlinear_resistance_table.point_count < 2
                    add_issue!(
                        result.validation,
                        missing_data(
                            "line $line_no",
                            "expected at least two nonlinear resistance points before table sentinel",
                        ),
                    )
                else
                    record_card!(result, :fixed_field)
                    record_card!(result, :nonlinear_resistance_table_end)
                end
                pending_nonlinear_resistance_table = nothing
                continue
            end
            marker_before_table = fixed_card_section_marker(tokens)
            named_marker_before_table = fixed_card_named_section_marker(tokens)
            begins_new_card =
                !isempty(tokens) && normalized_deck_token(tokens[1]) == "begin"
            if marker_before_table !== nothing || named_marker_before_table !== nothing ||
               begins_new_card
                add_issue!(
                    result.validation,
                    missing_data(
                        "line $line_no",
                        "expected nonlinear resistance table sentinel 9999. before the next deck section",
                    ),
                )
                pending_nonlinear_resistance_table = nothing
            elseif !pending_nonlinear_resistance_table.initialization_read
                parse_nonlinear_resistance_initialization_row!(
                    result,
                    line,
                    line_no,
                    pending_nonlinear_resistance_table.row_index,
                )
                pending_nonlinear_resistance_table = (
                    row_index = pending_nonlinear_resistance_table.row_index,
                    initialization_read = true,
                    point_count = pending_nonlinear_resistance_table.point_count,
                )
                continue
            else
                issues_before = length(result.validation.issues)
                parse_nonlinear_resistance_point_row!(
                    result,
                    line,
                    line_no,
                    pending_nonlinear_resistance_table.row_index,
                    pending_nonlinear_resistance_table.point_count + 1,
                )
                point_count =
                    length(result.validation.issues) == issues_before ?
                    pending_nonlinear_resistance_table.point_count + 1 :
                    pending_nonlinear_resistance_table.point_count
                pending_nonlinear_resistance_table = (
                    row_index = pending_nonlinear_resistance_table.row_index,
                    initialization_read = pending_nonlinear_resistance_table.initialization_read,
                    point_count = point_count,
                )
                continue
            end
        end
        if pending_triggered_timed_resistance_table !== nothing
            if triggered_timed_resistance_table_sentinel(line)
                if pending_triggered_timed_resistance_table.point_count < 1
                    add_issue!(
                        result.validation,
                        missing_data(
                            "line $line_no",
                            "expected at least one timed-resistance schedule point before table sentinel",
                        ),
                    )
                else
                    record_card!(result, :fixed_field)
                    record_card!(result, :triggered_timed_resistance_table_end)
                end
                pending_triggered_timed_resistance_table = nothing
                continue
            end
            marker_before_table = fixed_card_section_marker(tokens)
            named_marker_before_table = fixed_card_named_section_marker(tokens)
            begins_new_card =
                !isempty(tokens) && normalized_deck_token(tokens[1]) == "begin"
            if marker_before_table !== nothing || named_marker_before_table !== nothing ||
               begins_new_card
                add_issue!(
                    result.validation,
                    missing_data(
                        "line $line_no",
                        "expected timed-resistance table sentinel 9999. before the next deck section",
                    ),
                )
                pending_triggered_timed_resistance_table = nothing
            else
                issues_before = length(result.validation.issues)
                parse_triggered_timed_resistance_point_row!(
                    result,
                    line,
                    line_no,
                    pending_triggered_timed_resistance_table.row_index,
                    pending_triggered_timed_resistance_table.point_count + 1,
                )
                point_count =
                    length(result.validation.issues) == issues_before ?
                    pending_triggered_timed_resistance_table.point_count + 1 :
                    pending_triggered_timed_resistance_table.point_count
                pending_triggered_timed_resistance_table = (
                    row_index = pending_triggered_timed_resistance_table.row_index,
                    point_count = point_count,
                )
                continue
            end
        end
        if pending_switching_nonlinear_resistor_table !== nothing
            sentinel = switching_nonlinear_resistor_table_sentinel(line)
            if sentinel !== nothing
                if pending_switching_nonlinear_resistor_table.point_count < 1
                    add_issue!(
                        result.validation,
                        missing_data(
                            "line $line_no",
                            "expected at least one switching nonlinear resistor V-I point before table sentinel",
                        ),
                    )
                else
                    finish_switching_nonlinear_resistor_table!(
                        result,
                        pending_switching_nonlinear_resistor_table.row_index,
                        sentinel.single_flash,
                    )
                    record_card!(result, :fixed_field)
                    record_card!(result, :switching_nonlinear_resistor_table_end)
                end
                pending_switching_nonlinear_resistor_table = nothing
                continue
            end
            marker_before_table = fixed_card_section_marker(tokens)
            named_marker_before_table = fixed_card_named_section_marker(tokens)
            begins_new_card =
                !isempty(tokens) && normalized_deck_token(tokens[1]) == "begin"
            if marker_before_table !== nothing || named_marker_before_table !== nothing ||
               begins_new_card
                add_issue!(
                    result.validation,
                    missing_data(
                        "line $line_no",
                        "expected switching nonlinear resistor table sentinel 9999. before the next deck section",
                    ),
                )
                pending_switching_nonlinear_resistor_table = nothing
            else
                issues_before = length(result.validation.issues)
                parse_switching_nonlinear_resistor_point_row!(
                    result,
                    line,
                    line_no,
                    pending_switching_nonlinear_resistor_table.row_index,
                    pending_switching_nonlinear_resistor_table.point_count + 1,
                )
                point_count =
                    length(result.validation.issues) == issues_before ?
                    pending_switching_nonlinear_resistor_table.point_count + 1 :
                    pending_switching_nonlinear_resistor_table.point_count
                pending_switching_nonlinear_resistor_table = (
                    row_index = pending_switching_nonlinear_resistor_table.row_index,
                    point_count = point_count,
                )
                continue
            end
        end
        if pending_synchronous_machine_data !== nothing
            finished = parse_synchronous_machine_data_card!(
                result,
                pending_synchronous_machine_data,
                line,
                tokens,
                line_no,
            )
            if finished
                pending_synchronous_machine_data = nothing
                active_section = :blank_source
            end
            continue
        elseif active_section == :blank_source &&
               synchronous_machine_terminal_voltage_card(line)
            pending_synchronous_machine_data = SynchronousMachineDataParseState(
                length(result.synchronous_machine_output_summary_rows) + 1,
            )
            parse_synchronous_machine_data_card!(
                result,
                pending_synchronous_machine_data,
                line,
                tokens,
                line_no,
            )
            active_section = :blank_source
            continue
        end
        if pending_cable_constants !== nothing
            finished = parse_cable_constants_card!(
                result,
                pending_cable_constants,
                line,
                tokens,
                line_no,
            )
            if finished
                pending_cable_constants = nothing
                active_section = nothing
            end
            continue
        end
        if active_section == :blank_source &&
           dc_simulator_primary_card_candidate(line)
            pending_dc_simulator_source =
                parse_dc_simulator_primary_card(result, line, line_no)
            continue
        end
        output_end_sentinels_before = get(result.card_counts, :bpa_fixed_output_end_sentinel, 0)
        short_section_marker = fixed_card_short_section_marker(active_section, tokens)
        if short_section_marker !== nothing
            record_control_card!(result, short_section_marker.kind, tokens, line_no)
            active_section =
                short_section_marker.kind == :fixed_card_source_section_end &&
                fixed_source_load_flow_requested ?
                :fixed_source_load_flow : short_section_marker.next_section
            fixed_saturated_transformer_intake_active = false
            fixed_branch_vintage_mode = 0
            fixed_card_kc_lee_transform_active = false
            continue
        end
        section_marker = fixed_card_section_marker(tokens)
        if section_marker !== nothing
            if section_marker.kind == :fixed_card_case_terminator
                record_case_boundary!(result, :run_termination, tokens, line_no)
            end
            record_control_card!(result, section_marker.kind, tokens, line_no)
            if active_section == :control_system_hybrid &&
               section_marker.kind == :fixed_card_tacs_section_end
                if control_system_function_coefficients_pending(result)
                    row = last(result.control_system_function_rows)
                    record_control_system_card_blocker!(
                        result,
                        :control_system_function_coefficients_incomplete,
                    )
                    add_issue!(
                        result.validation,
                        missing_data(
                            "line $(row.line_no)",
                            "order-$(row.order) control-system function $(row.name) " *
                            "requires numerator and denominator coefficient cards",
                        ),
                    )
                end
                record_card!(result, :control_system_card_input)
                record_card!(result, :control_system_section_end)
            end
            active_section =
                section_marker.kind == :fixed_card_source_section_end &&
                fixed_source_load_flow_requested ?
                :fixed_source_load_flow : section_marker.next_section
            fixed_saturated_transformer_intake_active = false
            fixed_branch_vintage_mode = 0
            fixed_card_kc_lee_transform_active = false
            continue
        end
        if active_section == :fixed_source_load_flow
            if parse_fixed_source_load_flow_card!(result, line, line_no)
                active_section = :blank_initial_condition_or_output
                fixed_source_load_flow_requested = false
            end
            continue
        end
        named_section_marker = fixed_card_named_section_marker(tokens)
        if named_section_marker !== nothing
            record_control_card!(result, named_section_marker.kind, tokens, line_no)
            if named_section_marker.kind == :control_system_hybrid_section
                record_card!(result, :control_system_card_input)
                record_card!(result, :control_system_section_start)
            elseif named_section_marker.kind == :cable_constants_section
                pending_cable_constants = CableConstantsParseState(line_no)
                record_card!(result, :cable_constants_section_start)
            end
            active_section = named_section_marker.next_section
            fixed_saturated_transformer_intake_active = false
            fixed_branch_vintage_mode = 0
            fixed_card_kc_lee_transform_active = false
            continue
        end
        if active_section === nothing &&
           fixed_card_miscellaneous_control_candidate(tokens)
            parse_fixed_card_miscellaneous_control!(result, tokens, line_no)
            fixed_miscellaneous_control_count += 1
            if fixed_miscellaneous_control_count >= 2
                active_section = :blank_branch
            end
            continue
        end
        if active_section == :blank_branch &&
           coupled_lumped_numeric_continuation_row(result, fixed_image(line))
            append_coupled_lumped_continuation_row!(
                result,
                fixed_image(line),
                line_no,
                length(result.validation.issues),
            )
            continue
        end
        if active_section == :blank_branch &&
           fixed_miscellaneous_control_count == 2 &&
           isempty(result.elements) &&
           fixed_card_miscellaneous_control_candidate(tokens)
            parse_fixed_card_miscellaneous_control!(result, tokens, line_no)
            fixed_miscellaneous_control_count += 1
            continue
        end
        tacs_dimension_marker = tacs_dimension_request_marker(tokens)
        if tacs_dimension_marker !== nothing
            record_tacs_dimension_request_marker!(
                result,
                tacs_dimension_marker.allocation_kind,
                tokens,
                line_no,
            )
            pending_tacs_dimension_request = (
                allocation_kind = tacs_dimension_marker.allocation_kind,
                line_no = line_no,
                request_text = join(token_strings(tokens), " "),
            )
            active_section = next_deck_section(active_section, tokens)
            continue
        end
        if printout_frequency_change_request_marker(tokens)
            record_printout_frequency_change_request_marker!(result, tokens, line_no)
            pending_printout_frequency_change_request = (
                line_no = line_no,
                request_text = join(token_strings(tokens), " "),
            )
            active_section = next_deck_section(active_section, tokens)
            continue
        end
        if abort_data_case_card(tokens)
            parse_control_card!(result, tokens, line_no)
            abort_data_case_line = line_no
            active_section = nothing
            fixed_miscellaneous_control_count = 0
            fixed_saturated_transformer_intake_active = false
            pending_zinc_oxide_table = nothing
            pending_nonlinear_resistance_table = nothing
            pending_triggered_timed_resistance_table = nothing
            pending_switching_nonlinear_resistor_table = nothing
            pending_piecewise_nonlinear_inductor_table = nothing
            pending_arrester_constants = nothing
            pending_universal_machine_data = nothing
            pending_synchronous_machine_data = nothing
            pending_cable_constants = nothing
            pending_dc_simulator_source = nothing
            pending_cascaded_pi = nothing
            pending_sampled_frequency_line = nothing
            fixed_branch_vintage_mode = 0
            fixed_card_kc_lee_transform_active = false
            continue
        end
        if universal_machine_data_section_marker(tokens)
            pending_universal_machine_data = UniversalMachineDataParseState()
            record_control_card!(result, :universal_machine_data_section, tokens, line_no)
            record_card!(result, :universal_machine_data_section)
            active_section = :blank_universal_machine
            fixed_saturated_transformer_intake_active = false
            continue
        end
        if pending_universal_machine_data !== nothing
            finished = parse_universal_machine_data_card!(
                result,
                pending_universal_machine_data,
                line,
                tokens,
                line_no,
            )
            if finished
                pending_universal_machine_data = nothing
                active_section = :blank_source
            else
                active_section = :blank_universal_machine
            end
            continue
        end
        if pending_arrester_constants !== nothing
            marker_before_constants = fixed_card_section_marker(tokens)
            named_marker_before_constants = fixed_card_named_section_marker(tokens)
            begins_new_card =
                !isempty(tokens) && normalized_deck_token(tokens[1]) == "begin"
            if marker_before_constants !== nothing ||
               named_marker_before_constants !== nothing ||
               begins_new_card
                add_issue!(
                    result.validation,
                    missing_data(
                        "line $line_no",
                        "expected arrester constants before the next deck section",
                    ),
                )
                pending_arrester_constants = nothing
            else
                issues_before = length(result.validation.issues)
                values_read = parse_arrester_constant_row!(
                    result,
                    line,
                    line_no,
                    pending_arrester_constants.row_index,
                    pending_arrester_constants.next_constant_index,
                )
                if length(result.validation.issues) == issues_before
                    next_constant_index =
                        pending_arrester_constants.next_constant_index + values_read
                    pending_arrester_constants = next_constant_index > 18 ? nothing : (
                        row_index = pending_arrester_constants.row_index,
                        next_constant_index = next_constant_index,
                    )
                end
                continue
            end
        end
        if active_section == :blank_branch && !fixed_saturated_transformer_intake_active &&
           fixed_card_saturated_transformer_header_card(tokens)
            fixed_saturated_transformer_intake_active = true
            record_card!(result, :fixed_field)
            record_card!(result, :fixed_card_saturated_transformer_intake)
            continue
        end
        if active_section == :blank_branch && generator_equivalent_header_card(line)
            pending_generator_equivalent =
                start_generator_equivalent_parse!(result, line, line_no)
            continue
        end
        if active_section == :blank_branch && cascaded_pi_header_card(line)
            pending_cascaded_pi = parse_cascaded_pi_header!(result, line, line_no)
            continue
        end
        if active_section == :blank_branch && semlyen_line_header_card(line)
            pending_semlyen_line = start_semlyen_line_parse!(
                result,
                line,
                line_no,
                semlyen_transform_imaginary_parts,
            )
            continue
        end
        if active_section == :blank_branch && rational_frequency_line_header_card(line)
            pending_rational_frequency_line = start_rational_frequency_line_parse!(
                result,
                line,
                line_no,
            )
            continue
        end
        if active_section == :blank_branch && fixed_saturated_transformer_intake_active
            record_card!(result, :fixed_field)
            record_card!(result, :fixed_card_saturated_transformer_intake)
            continue
        end
        vintage_mode = fixed_card_branch_vintage_mode(tokens)
        if active_section == :blank_branch && vintage_mode !== nothing
            fixed_branch_vintage_mode = vintage_mode
            fixed_card_kc_lee_transform_active =
                vintage_mode == 0 &&
                get(result.card_counts, :fixed_card_kc_lee_untransposed_line_row, 0) > 0
            record_control_card!(
                result,
                vintage_mode == 0 ?
                :fixed_card_branch_vintage_mode_disabled :
                :fixed_card_branch_vintage_mode_enabled,
                tokens,
                line_no,
            )
            record_card!(result, :fixed_card_branch_vintage_mode)
            continue
        end
        if active_section == :blank_branch && fixed_card_kc_lee_transform_active &&
           fixed_card_kc_lee_transform_row_candidate(tokens)
            parse_fixed_card_kc_lee_transform_row!(result, tokens, line_no)
            continue
        end
        zinc_oxide_rows_before = length(result.zinc_oxide_nonlinear_rows)
        nonlinear_resistance_rows_before = length(result.nonlinear_resistance_rows)
        triggered_timed_resistance_rows_before =
            length(result.triggered_timed_resistance_rows)
        switching_resistor_rows_before = length(result.switching_nonlinear_resistor_rows)
        piecewise_nonlinear_inductor_rows_before =
            length(result.piecewise_nonlinear_inductor_rows)
        hysteretic_inductor_rows_before = length(result.hysteretic_inductor_rows)
        arrester_rows_before = length(result.arrester_nonlinear_rows)
        kc_lee_rows_before =
            get(result.card_counts, :fixed_card_kc_lee_untransposed_line_row, 0)
        coupled_line_rows_before = length(result.coupled_line_rows)
        if parse_fixed_field_section_card!(
            result,
            active_section,
            line,
            tokens,
            line_no;
            branch_vintage_mode = fixed_branch_vintage_mode,
        )
            if active_section == :blank_branch &&
               length(result.zinc_oxide_nonlinear_rows) > zinc_oxide_rows_before &&
               result.zinc_oxide_nonlinear_rows[end].source_kind != :copy_reference
                pending_zinc_oxide_table = (
                    row_index = length(result.zinc_oxide_nonlinear_rows),
                    initialization_read = false,
                    breakpoint_count = 0,
                )
            end
            if active_section == :blank_branch &&
               length(result.nonlinear_resistance_rows) > nonlinear_resistance_rows_before &&
               result.nonlinear_resistance_rows[end].source_kind != :copy_reference
                pending_nonlinear_resistance_table = (
                    row_index = length(result.nonlinear_resistance_rows),
                    initialization_read = false,
                    point_count = 0,
                )
            end
            if active_section == :blank_branch &&
               length(result.triggered_timed_resistance_rows) >
               triggered_timed_resistance_rows_before &&
               result.triggered_timed_resistance_rows[end].source_kind != :copy_reference
                pending_triggered_timed_resistance_table = (
                    row_index = length(result.triggered_timed_resistance_rows),
                    point_count = 0,
                )
            end
            if active_section == :blank_branch &&
               length(result.switching_nonlinear_resistor_rows) >
               switching_resistor_rows_before &&
               result.switching_nonlinear_resistor_rows[end].source_kind != :copy_reference
                pending_switching_nonlinear_resistor_table = (
                    row_index = length(result.switching_nonlinear_resistor_rows),
                    point_count = 0,
                )
            end
            if active_section == :blank_branch &&
               length(result.piecewise_nonlinear_inductor_rows) >
               piecewise_nonlinear_inductor_rows_before &&
               result.piecewise_nonlinear_inductor_rows[end].source_kind != :copy_reference
                pending_piecewise_nonlinear_inductor_table = (
                    row_index = length(result.piecewise_nonlinear_inductor_rows),
                    point_count = 0,
                )
            end
            if active_section == :blank_branch &&
               length(result.hysteretic_inductor_rows) > hysteretic_inductor_rows_before &&
               result.hysteretic_inductor_rows[end].source_kind != :copy_reference
                pending_hysteretic_inductor_table = (
                    row_index = length(result.hysteretic_inductor_rows),
                    point_count = 0,
                )
            end
            if active_section == :blank_branch &&
               length(result.arrester_nonlinear_rows) > arrester_rows_before &&
               result.arrester_nonlinear_rows[end].source_kind != :copy_reference
                pending_arrester_constants = (
                    row_index = length(result.arrester_nonlinear_rows),
                    next_constant_index = 1,
                )
            end
            if active_section == :blank_branch &&
               get(result.card_counts, :fixed_card_kc_lee_untransposed_line_row, 0) >
               kc_lee_rows_before
                fixed_card_kc_lee_transform_active = false
            end
            if active_section == :blank_branch &&
               length(result.coupled_line_rows) > coupled_line_rows_before &&
               result.coupled_line_rows[end].sampled_frequency_data_requested
                pending_sampled_frequency_line = start_sampled_frequency_line_parse(
                    result,
                    length(result.coupled_line_rows),
                )
            end
            if active_section == :blank_output &&
               get(result.card_counts, :bpa_fixed_output_end_sentinel, 0) > output_end_sentinels_before
                active_section = nothing
            else
                active_section = next_deck_section(active_section, tokens)
            end
            continue
        end
        parse_deck_card!(result, tokens, line_no)
        study_option_request_kind(tokens) == :fixed_source_declaration &&
            (fixed_source_load_flow_requested = true)
        active_section = next_deck_section(active_section, tokens)
        if !isempty(tokens) && normalized_deck_token(tokens[1]) == "begin"
            fixed_miscellaneous_control_count = 0
            fixed_saturated_transformer_intake_active = false
            pending_zinc_oxide_table = nothing
            pending_nonlinear_resistance_table = nothing
            pending_triggered_timed_resistance_table = nothing
            pending_switching_nonlinear_resistor_table = nothing
            pending_piecewise_nonlinear_inductor_table = nothing
            pending_hysteretic_inductor_table = nothing
            pending_arrester_constants = nothing
            pending_universal_machine_data = nothing
            pending_synchronous_machine_data = nothing
            pending_dc_simulator_source = nothing
            pending_semlyen_line = nothing
            pending_rational_frequency_line = nothing
            fixed_branch_vintage_mode = 0
            fixed_card_kc_lee_transform_active = false
        end
    end
    if pending_generator_equivalent !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "end of deck",
                "frequency-dependent generator-equivalent card sequence ended before both modal sentinels",
            ),
        )
    end
    if pending_dc_simulator_source !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "line $(pending_dc_simulator_source.line_no)",
                "expected the second DC-simulator source card",
            ),
        )
    end
    if pending_tacs_dimension_request !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "line $(pending_tacs_dimension_request.line_no)",
                "expected numeric card after $(pending_tacs_dimension_request.allocation_kind) TACS dimension request",
            ),
        )
    end
    if pending_printout_frequency_change_request !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "line $(pending_printout_frequency_change_request.line_no)",
                "expected numeric card after printout frequency change request",
            ),
        )
    end
    if pending_zinc_oxide_table !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "end of deck",
                "expected zinc-oxide table sentinel 9999. after type-92 row",
            ),
        )
    end
    if pending_nonlinear_resistance_table !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "end of deck",
                "expected nonlinear resistance table sentinel 9999. after type-91/type-92 row",
            ),
        )
    end
    if pending_triggered_timed_resistance_table !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "end of deck",
                "expected timed-resistance table sentinel 9999. after the resistor row",
            ),
        )
    end
    if pending_switching_nonlinear_resistor_table !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "end of deck",
                "expected switching nonlinear resistor table sentinel 9999. after the resistor row",
            ),
        )
    end
    if pending_piecewise_nonlinear_inductor_table !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "end of deck",
                "expected nonlinear-inductor table sentinel 9999. after the inductor row",
            ),
        )
    end
    if pending_hysteretic_inductor_table !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "end of deck",
                "expected type-96 table sentinel 9999. after hysteretic-inductor row",
            ),
        )
    end
    if pending_arrester_constants !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "end of deck",
                "expected 18 arrester constants after type-94 row",
            ),
        )
    end
    if pending_universal_machine_data !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "end of deck",
                "expected BLANK card terminating all universal-machine data",
            ),
        )
    end
    if pending_synchronous_machine_data !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "end of deck",
                "expected FINISH terminating type-59 synchronous-machine data",
            ),
        )
    end
    if pending_cable_constants !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "end of deck",
                "expected BLANK CARD TERMINATING CABLE CONSTANTS CASES",
            ),
        )
    end
    if pending_cascaded_pi !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "end of deck",
                "expected a cascaded PI section configuration and STOP CASCADE terminator",
            ),
        )
    end
    if pending_sampled_frequency_line !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "end of deck",
                "expected sampled frequency-dependent line configuration and weighting points",
            ),
        )
    end
    if pending_semlyen_line !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "end of deck",
                "expected the remaining Semlyen modes and both modal transformation matrices",
            ),
        )
    end
    if pending_rational_frequency_line !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "end of deck",
                "expected the remaining rational line modes and complex modal transform",
            ),
        )
    end
    if abort_data_case_line !== nothing
        add_issue!(
            result.validation,
            missing_data(
                "line $(abort_data_case_line)",
                "expected BEGIN NEW DATA CASE after ABORT DATA CASE request",
            ),
        )
    end
    for row in result.control_system_switch_rows
        if row.switching_delay_or_model !== missing &&
           row.switching_delay_or_model == 7777.0 &&
           row.delayed_arc === nothing
            add_issue!(
                result.validation,
                missing_data(
                    "line $(row.line_no)",
                    "expected four-field delayed-arc continuation after the switch row",
                ),
            )
        end
    end
    validate_over5a_controlled_source_rows!(result)
    validate_fixed_source_load_flow_rows!(result)
    validate_coupled_line_rows!(result)
    validate_coupled_phase_pi_section_rows!(result)
    validate_cascaded_pi_requests!(result)
    validate_steady_state_frequency_networks!(result)
    validate_over16_output_channels!(result)
    validate_over16_branch_voltage_outputs!(result)
    validate_over16_branch_current_outputs!(result)
    validate_over16_branch_power_outputs!(result)
    return result
end

function _fixed_column_c_capacitance_card(line::AbstractString)::Bool
    return isempty(strip(fixed_field(line, 1, 2))) &&
           uppercase(strip(fixed_field(line, 3, 8))) == "C" &&
           fixed_float_value(line, 39, 44) !== nothing
end

function strip_deck_line(raw_line)::String
    line = String(raw_line)
    stripped = lstrip(line)
    if isempty(stripped)
        return ""
    end
    first_char = first(stripped)
    if first_char == '$' && !startswith(uppercase(stripped), "\$VINTAGE")
        return ""
    end
    if uppercase(first_char) == 'C' &&
       (length(stripped) == 1 || isspace(stripped[nextind(stripped, firstindex(stripped))])) &&
       !_fixed_column_c_capacitance_card(line)
        return ""
    end
    line = split(line, '#'; limit=2)[1]
    line = split(line, '!'; limit=2)[1]
    return rstrip(line)
end

deck_tokens(line::AbstractString) = split(replace(String(line), ',' => ' '))

const EXPLICIT_DECK_CARD_TOKENS = Set((
    "begin", "end", "section", "blank", "eldc", "last", "adc", "abort",
    "branch", "bpa_branch", "bus", "node",
    "source", "bpa_source", "current", "current_source", "source_current",
    "bpa_current_source", "conductance", "resistor", "rl", "inductor",
    "capacitor", "breqiv", "breqiv3", "bergeron_line", "line",
    "switch", "bpa_switch", "time_switch", "output", "bpa_output",
    "power", "absolute", "relative", "aumd", "atd", "rtd", "ow1", "ow8",
    "pvm", "peak", "adp", "alternate", "twl", "tacs", "custom", "msl",
    "modify",
    "obc", "cs", "mdc", "rte", "bpvs", "todr", "usst", "omit",
    "miscellaneous", "redefine", "time", "user", "printout",
    "fr", "file", "request", "cz", "convert", "hs", "hauer", "setup",
    "lmfs", "model", "freq", "scan",
    "tl", "type99", "ff", "free", "format", "pph", "plotter", "plpi", "printer", "mvo", "mode",
    "asu", "analytic", "lopo", "limit", "lbu", "linear", "an", "auto",
    "rb", "renumber", "fs", "frequency", "d", "diagnostic", "ui", "hr",
    "high", "ao", "average", "sos", "zo", "zinc", "fxs", "fix",
    "over16_output", "over16_voltage_output",
    "over16_branch_voltage", "over16_branch_voltage_output",
    "over16_branch_current", "over16_branch_current_output",
    "over16_branch_power", "over16_branch_power_output",
    "over16_branch_energy", "over16_branch_energy_output",
    "over16_source_card",
    "over16_source_interpolation", "over16_source_interpolation_values",
    "over16_source_interp", "over16_source_interp_values",
    "over16_source_tacs_override", "over16_vstacs_override",
    "over16_source_analytic", "over16_source_analytic_values",
))

const BPA_FIXED_SOURCE_CONDUCTANCE = 1.0e12
const BPA_FIXED_ACCEPTED_SOURCE_TYPES = Tuple(1:18)
const BPA_FIXED_VOLTBC_SOURCE_TYPES = Tuple(1:10)
const BPA_FIXED_ANALYTIC_SOURCE_TYPES = (11, 12, 13, 14, 15)
const BPA_FIXED_CONTROLLED_SOURCE_TYPES = (16, 17, 18)
const BPA_FIXED_INCREMENTAL_IFORM_SOURCE_TYPES = BPA_FIXED_VOLTBC_SOURCE_TYPES
const BPA_FIXED_INFINITE_STOP_SOURCE_TYPES = Tuple(1:15)

const BPA_FIXED_UNSUPPORTED_SECTION_BLOCKERS = Dict{Symbol,Tuple{Symbol,String}}(
    :blank_tacs => (
        :bpa_fixed_tacs_blocked,
        "Unsupported BPA fixed-field BLANK TACS row: TACS/NTACS/ELEC execution requires RC5 state ownership and oracle coverage",
    ),
    :blank_machine => (
        :bpa_fixed_machine_blocked,
        "Unsupported BPA fixed-field BLANK MACHINE row: machine deck rows require RC6 machine solve/coupling ownership and oracle coverage",
    ),
    :blank_synchronous_machine => (
        :bpa_fixed_machine_blocked,
        "Unsupported BPA fixed-field BLANK SYNCHRONOUS_MACHINE row: machine deck rows require RC6 machine solve/coupling ownership and oracle coverage",
    ),
    :blank_universal_machine => (
        :bpa_fixed_machine_blocked,
        "Unsupported BPA fixed-field BLANK UNIVERSAL_MACHINE row: machine deck rows require RC6 machine solve/coupling ownership and oracle coverage",
    ),
    :blank_line => (
        :bpa_fixed_line_transformer_network_blocked,
        "Unsupported BPA fixed-field BLANK LINE row: distributed, multiphase, and frequency-dependent line rows require RC7 line ownership and oracle coverage",
    ),
    :blank_transformer => (
        :bpa_fixed_line_transformer_network_blocked,
        "Unsupported BPA fixed-field BLANK TRANSFORMER row: transformer/network rows require translated branch or line ownership and oracle coverage",
    ),
    :blank_network => (
        :bpa_fixed_line_transformer_network_blocked,
        "Unsupported BPA fixed-field BLANK NETWORK row: network rows require translated topology/admittance ownership and oracle coverage",
    ),
    :blank_frequency => (
        :bpa_fixed_line_transformer_network_blocked,
        "Unsupported BPA fixed-field BLANK FREQUENCY row: frequency-scan and frequency-dependent line rows require RC7 ownership and oracle coverage",
    ),
    :blank_nonlinear => (
        :bpa_fixed_nonlinear_blocked,
        "Unsupported BPA fixed-field BLANK NONLINEAR row: nonlinear table/card rows require RC3 nonlinear ownership and oracle coverage",
    ),
    :blank_plot => (
        :bpa_fixed_report_plot_blocked,
        "Unsupported BPA fixed-field BLANK PLOT row: report/plot rows require RC8 report and plot ownership and oracle coverage",
    ),
    :blank_request => (
        :bpa_fixed_report_plot_blocked,
        "Unsupported BPA fixed-field BLANK REQUEST row: report/request rows require RC8 report ownership and oracle coverage",
    ),
    :blank_monitor => (
        :bpa_fixed_report_plot_blocked,
        "Unsupported BPA fixed-field BLANK MONITOR row: monitor/report rows require RC8 report ownership and oracle coverage",
    ),
    :blank_statistics => (
        :bpa_fixed_report_plot_blocked,
        "Unsupported BPA fixed-field BLANK STATISTICS row: statistical output rows require translated source/switch/report ownership and oracle coverage",
    ),
    :blank_load => (
        :bpa_fixed_deck_build_blocked,
        "Unsupported BPA fixed-field BLANK LOAD row: load rows require translated shared-model or branch/source ownership and BUILD* oracle coverage",
    ),
    :blank_data => (
        :bpa_fixed_deck_build_blocked,
        "Unsupported BPA fixed-field BLANK DATA row: general data rows require DATAIN/TREAD/BUILD* ownership and oracle coverage",
    ),
    :blank_misc => (
        :bpa_fixed_deck_build_blocked,
        "Unsupported BPA fixed-field BLANK MISC row: miscellaneous data rows require DATAIN/TREAD/BUILD* ownership and oracle coverage",
    ),
    :blank_miscellaneous => (
        :bpa_fixed_deck_build_blocked,
        "Unsupported BPA fixed-field BLANK MISCELLANEOUS row: miscellaneous data rows require DATAIN/TREAD/BUILD* ownership and oracle coverage",
    ),
    :blank_measurement => (
        :bpa_fixed_report_plot_blocked,
        "Unsupported BPA fixed-field BLANK MEASUREMENT row: measurement/output rows require translated report ownership and oracle coverage",
    ),
)

function bpa_fixed_source_type_accepted(source_type::Int)
    return source_type in BPA_FIXED_ACCEPTED_SOURCE_TYPES || source_type >= 60
end

function bpa_fixed_source_blocker(source_type::Int)
    if source_type == 19
        return (:bpa_fixed_source_blocked_tacs_machine_or_statistics,
                "Malformed OVER5A type-19 row: a universal-machine section marker is the standalone fixed-field value 19 or the explicit free-field form 19 UM")
    end
    return (:bpa_fixed_source_blocked_other,
            "Unsupported OVER5A fixed-field source type $source_type")
end

function next_deck_section(current_section, tokens)
    isempty(tokens) && return current_section
    card = normalized_deck_token(tokens[1])
    if card == "blank" && length(tokens) >= 2
        category = normalized_deck_token(tokens[2])
        kind = get(BPA_BLANK_SECTION_KINDS, category, nothing)
        return kind === nothing ? current_section : kind
    elseif card == "begin"
        return nothing
    elseif card == "end"
        return nothing
    end
    return current_section
end

function fixed_card_section_marker(tokens)
    isempty(tokens) && return nothing
    normalized_deck_token(tokens[1]) == "blank" || return nothing
    length(tokens) == 1 && return (kind = :fixed_card_case_terminator, next_section = nothing)
    all_tail = Set(normalized_deck_token(token) for token in tokens[2:end])
    if ("terminating" in all_tail || "ending" in all_tail) && "case" in all_tail
        return (kind = :fixed_card_case_terminator, next_section = nothing)
    end
    length(tokens) >= 3 || return nothing
    normalized_deck_token(tokens[2]) == "card" || return nothing
    tail = Set(normalized_deck_token(token) for token in tokens[3:end])
    if "tacs" in tail && "data" in tail
        return (kind = :fixed_card_tacs_section_end, next_section = :blank_branch)
    elseif "initial" in tail && "conditions" in tail
        return (kind = :fixed_card_initial_conditions_section_end, next_section = :blank_output)
    elseif "conductor" in tail || "conductors" in tail
        return (kind = :fixed_card_conductor_section_end, next_section = :line_constants_frequency)
    elseif "frequency" in tail
        return (kind = :fixed_card_frequency_section_end, next_section = :line_constants_conductor)
    elseif "line" in tail && "constants" in tail && "cases" in tail
        return (kind = :fixed_card_line_constants_section_end, next_section = nothing)
    elseif "branch" in tail || "branches" in tail
        return (kind = :fixed_card_branch_section_end, next_section = :blank_switch)
    elseif "switch" in tail || "switches" in tail
        return (kind = :fixed_card_switch_section_end, next_section = :blank_source)
    elseif "source" in tail || "sources" in tail
        return (kind = :fixed_card_source_section_end, next_section = :blank_initial_condition_or_output)
    elseif "output" in tail || "outputs" in tail || "requests" in tail
        return (kind = :fixed_card_output_section_end, next_section = :blank_plot)
    elseif "plot" in tail || "plots" in tail
        return (kind = :fixed_card_plot_section_end, next_section = nothing)
    end
    return nothing
end

function fixed_card_short_section_marker(current_section, tokens)
    isempty(tokens) && return nothing
    normalized_deck_token(tokens[1]) == "blank" || return nothing
    length(tokens) == 2 || return nothing
    category = normalized_deck_token(tokens[2])
    if current_section == :blank_branch && category == "branch"
        return (kind = :fixed_card_branch_section_end, next_section = :blank_switch)
    elseif current_section == :blank_switch && category == "switch"
        return (kind = :fixed_card_switch_section_end, next_section = :blank_source)
    elseif current_section == :blank_source && category == "source"
        return (
            kind = :fixed_card_source_section_end,
            next_section = :blank_initial_condition_or_output,
        )
    elseif current_section in (:blank_initial_condition_or_output, :blank_output) &&
           category == "output"
        return (kind = :fixed_card_output_section_end, next_section = :blank_plot)
    elseif current_section == :blank_plot && category == "plot"
        return (kind = :fixed_card_plot_section_end, next_section = nothing)
    end
    return nothing
end

function fixed_card_named_section_marker(tokens)
    isempty(tokens) && return nothing
    card = normalized_deck_token(tokens[1])
    if card == "tacs" && length(tokens) == 2 &&
       normalized_deck_token(tokens[2]) == "hybrid"
        return (kind = :control_system_hybrid_section, next_section = :control_system_hybrid)
    elseif card == "tacs" && length(tokens) == 2 &&
           normalized_deck_token(tokens[2]) == "data"
        return (kind = :blank_tacs, next_section = :blank_tacs)
    elseif card == "line" && length(tokens) >= 2 &&
           normalized_deck_token(tokens[2]) == "constants"
        return (kind = :blank_line, next_section = :line_constants_conductor)
    elseif card == "cable" && length(tokens) >= 2 &&
           normalized_deck_token(tokens[2]) == "constants"
        return (kind = :cable_constants_section, next_section = :cable_constants)
    end
    return nothing
end

function fixed_card_miscellaneous_control_candidate(tokens)::Bool
    length(tokens) >= 2 || return false
    all(token -> tryparse_deck_float(String(token)) !== nothing, tokens[1:min(length(tokens), 4)])
end

function tacs_dimension_request_marker(tokens)
    isempty(tokens) && return nothing
    first_token = compact_deck_keyword(tokens[1])
    if first_token == "atd"
        return (allocation_kind = :absolute_tacs_dimension_request,)
    elseif first_token == "rtd"
        return (allocation_kind = :relative_tacs_dimension_request,)
    end
    length(tokens) >= 3 || return nothing
    third_token = compact_deck_keyword(tokens[3])
    third_token == "dimensions" || return nothing
    second_token = compact_deck_keyword(tokens[2])
    second_token == "tacs" || return nothing
    if first_token == "absolute"
        return (allocation_kind = :absolute_tacs_dimension_request,)
    elseif first_token == "relative"
        return (allocation_kind = :relative_tacs_dimension_request,)
    end
    return nothing
end

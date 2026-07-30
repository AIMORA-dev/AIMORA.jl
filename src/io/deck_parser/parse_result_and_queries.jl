
using ..Branches:
    CapacitorBranch,
    ConductanceBranch,
    CurrentInjection,
    GeneratorEquivalentModalBranch,
    IdealTransformerVoltageConstraint,
    SeriesRLBranch,
    SeriesRLCBranch,
    TheveninSource,
    single_phase_breqiv_history_injection,
    three_phase_breqiv_history_injection,
    generator_equivalent_history_injection
using ..Nodal: NodalSystem
using ..Sources:
    analytic_current_injection_source,
    analytic_thevenin_source,
    sinusoidal_thevenin_source,
    sinusoidal_value
using ..Switches: IdealSwitch, TimeSwitch
using ..Lines:
    BergeronLine,
    CascadedPiSeriesImpedance,
    CascadedPiShuntImpedance,
    CascadedPhasePiBlock,
    CoupledLumpedSequenceImpedance,
    CoupledLumpedPhasePiSection,
    CascadedPhasePiEquivalent,
    DistributedTransposedLineConstants,
    DistributedTransposedLineModalBranchState,
    DistributedTransposedLineSteadyStatePiEquivalent,
    DistributedTransposedLineHistoryState,
    DistributedTransposedLineCompanionAdmittance,
    TerminalSurgeImpedanceAdmittance,
    LineWeightingSamples,
    SampledLineWeightingCoefficients,
    SampledFrequencyDependentLine,
    SampledFrequencyDependentLineGroup,
    SemlyenModeParameters,
    SemlyenFrequencyDependentLine,
    ComplexModalBergeronLine,
    LineModalTransform,
    PoleResidueTransfer,
    PoleResidueReductionResult,
    RationalLineModeConversion,
    coupled_lumped_matrix_impedance,
    coupled_lumped_sequence_impedance,
    coupled_lumped_phase_pi_section,
    cascaded_phase_pi_equivalent,
    cascaded_pi_series_impedance,
    cascaded_pi_shunt_impedance,
    cascaded_phase_pi_block,
    distributed_transposed_line_constants,
    distributed_transposed_line_modal_branch_state,
    distributed_modal_line_branch_state,
    distributed_transposed_line_steady_state_pi_equivalent,
    distributed_transposed_line_history_state,
    distributed_transposed_line_companion_admittance,
    terminal_surge_impedance_admittance,
    line_weighting_samples,
    sampled_line_weighting_coefficients,
    sampled_frequency_dependent_line,
    sampled_frequency_dependent_line_group,
    semlyen_frequency_dependent_line,
    rational_frequency_dependent_mode_parameters,
    pole_residue_transfer_for_timestep,
    semlyen_rational_terms,
    semlyen_line_physical_checks,
    line_history_currents,
    line_surge_admittance,
    line_terminal_currents,
    line_terminal_voltages,
    line_traveling_waves
using ..ValidationCore:
    ValidationResult,
    add_issue!,
    assert_valid!,
    invalid_value,
    missing_data,
    unknown_field,
    validation_result

export DeckParseResult,
       DeckExecutableCase,
       DeckCaseSequence,
       DeckControlCard,
       DeckBergeronLineRow,
       DeckCoupledLineRow,
       DeckLineModalTransformRow,
       DeckCoupledPhasePiSectionRow,
       DeckCascadedPiSeriesImpedanceRow,
       DeckCascadedPiShuntImpedanceRow,
       DeckCascadedPiBlock,
       DeckCascadedPiRequestRow,
       DeckSampledFrequencyLineRow,
       DeckSemlyenModeRow,
       DeckSemlyenLineGroupRow,
       DeckRationalLineModeRow,
       DeckRationalLineGroupRow,
       DeckGeneratorEquivalentModalBranchRow,
       DeckGeneratorEquivalentRow,
       DeckLineConstantsConductorCard,
       DeckLineConstantsFrequencyCard,
       DeckLineConstantsPhysicalConductor,
       DeckCableConstantsFrequencyCard,
       DeckCableConstantsCase,
       DeckPowerFrequencyRequestRow,
       DeckUniversalMachineDimensionRequestRow,
       DeckUniversalMachineSectionRow,
       DeckUniversalMachineDefinitionRow,
       DeckUniversalMachineCoilRow,
       DeckUniversalMachineTerminalRow,
       DeckUniversalMachineGeneratedBranchRow,
       DeckUniversalMachineSpeedCapacitorRow,
       DeckUniversalMachineNodeSummaryRow,
       DeckUniversalMachineOutputSummaryRow,
       DeckSynchronousMachineTerminalVoltageRow,
       DeckSynchronousMachineToleranceRow,
       DeckSynchronousMachineParameterFittingRow,
       DeckSynchronousMachineModelParameterRow,
       DeckSynchronousMachineMassRow,
       DeckSynchronousMachineOutputRequestRow,
       DeckSynchronousMachineControlInterfaceRow,
       DeckSynchronousMachineOutputSummaryRow,
       DeckTACSDimensionRequestRow,
       DeckOutputWidthRequestRow,
       DeckPeakVoltageMonitorRequestRow,
       DeckDiagnosticPrintRequestRow,
       DeckTACSWarningLimitRequestRow,
       DeckPlotFileRequestRow,
       DeckSwitchLogicRequestRow,
       DeckSimulationControlRequestRow,
       DeckPrintoutFrequencyChangeRow,
       DeckStudyOptionRequestRow,
       DeckFixedSourceConstraintRow,
       DeckFixedSourceControlRow,
       DeckTimeHorizonControlRow,
       DeckOutputScheduleControlRow,
       DeckCaseBoundaryRow,
       DeckRestartRequest,
       DeckRestartSwitchMutation,
       DeckRestartControlSourceMutation,
       DeckControlSystemSignalTerm,
       DeckControlSystemFunctionRow,
       DeckControlSystemExpressionRow,
       DeckControlSystemSourceRow,
       DeckControlSystemDeviceRow,
       DeckControlSystemOutputRequestRow,
       DeckDelayedArcSwitchParameters,
       DeckControlSystemSwitchRow,
       DeckControlSystemSwitchCouplingRow,
       DeckNodeInitialConditionRow,
       DeckZincOxideBreakpointRow,
       DeckZincOxideInitializationRow,
       DeckZincOxideNonlinearRow,
       DeckNonlinearResistanceInitializationRow,
       DeckNonlinearResistancePointRow,
       DeckNonlinearResistanceRow,
       DeckTriggeredTimedResistancePointRow,
       DeckTriggeredTimedResistanceRow,
       DeckSwitchingNonlinearResistorPointRow,
       DeckSwitchingNonlinearResistorRow,
       DeckPiecewiseNonlinearInductorPointRow,
       DeckPiecewiseNonlinearInductorRow,
       DeckHystereticInductorPointRow,
       DeckHystereticInductorRow,
       DeckOVER2BranchRow,
       DeckOVER5ASourceRow,
       DeckOVER5SwitchRow,
       DeckOVER5ASourceUpdateInputs,
       DeckOVER15OutputRequestRow,
       DeckOVER16OutputChannel,
       DeckOVER16BranchVoltageOutput,
       DeckOVER16BranchCurrentOutput,
       DeckOVER16BranchPowerOutput,
       DeckOVER16SourceCardRow,
       DeckOVER16SourceInterpolationRow,
       DeckOVER16SourceTACSOverrideRow,
       DeckOVER16SourceAnalyticRow,
       DeckSaturatedTransformerBreakpointRow,
       DeckSaturatedTransformerCopyRow,
       DeckSaturatedTransformerHeaderRow,
       DeckSaturatedTransformerIntake,
       DeckSaturatedTransformerWindingRow,
       DeckTransformerBranchShuntCapacitanceRow,
       DeckModelSummary,
       assert_deck_valid!,
       deck_asset_tables,
       deck_branch_capacitance_values,
       deck_branch_conductance_values,
       deck_branch_count,
       deck_branch_from_node_indices,
       deck_branch_from_node_names,
       deck_branch_inductance_values,
       deck_branch_kinds,
       deck_branch_line_numbers,
       deck_branch_names,
       deck_branch_previous_current_values,
       deck_branch_previous_voltage_values,
       deck_branch_resistance_values,
       deck_branch_to_node_indices,
       deck_branch_to_node_names,
       deck_control_card_kinds,
       deck_control_card_labels,
       deck_control_card_line_numbers,
       deck_control_card_tokens,
       deck_control_system_device_rows,
       deck_control_system_function_rows,
       deck_control_system_expression_rows,
       deck_control_system_output_request_rows,
       deck_control_system_source_rows,
       deck_control_system_switch_coupling_rows,
       deck_control_system_switch_rows,
       deck_node_initial_condition_rows,
       deck_fixed_time_horizon_options,
       deck_output_schedule_options,
       deck_bergeron_line_attenuation_values,
       deck_bergeron_line_delay_step_counts,
       deck_bergeron_line_dt_s_values,
       deck_bergeron_line_from_node_indices,
       deck_bergeron_line_from_node_names,
       deck_bergeron_line_history_current_from_values,
       deck_bergeron_line_history_current_to_values,
       deck_bergeron_line_line_numbers,
       deck_bergeron_line_names,
       deck_bergeron_line_rows,
       deck_bergeron_line_surge_admittance_values,
       deck_bergeron_line_surge_impedance_values,
       deck_bergeron_line_terminal_current_from_values,
       deck_bergeron_line_terminal_current_to_values,
       deck_bergeron_line_terminal_voltage_from_values,
       deck_bergeron_line_terminal_voltage_to_values,
       deck_bergeron_line_to_node_indices,
       deck_bergeron_line_to_node_names,
       deck_bergeron_line_travel_time_s_values,
       deck_bergeron_line_traveling_wave_from_values,
       deck_bergeron_line_traveling_wave_to_values,
       deck_bergeron_line_write_indices,
       deck_coupled_line_from_node_indices,
       deck_coupled_line_from_node_names,
       deck_coupled_line_kinds,
       deck_coupled_line_line_numbers,
       deck_coupled_line_names,
       deck_coupled_line_phase_indices,
       deck_coupled_line_rows,
       deck_line_modal_transform_rows,
       deck_coupled_line_sequence_inductance_values,
       deck_coupled_line_sequence_resistance_values,
       deck_coupled_line_to_node_indices,
       deck_coupled_line_to_node_names,
       deck_coupled_line_types,
       deck_coupled_lumped_sequence_impedances,
       deck_coupled_phase_pi_section_rows,
       deck_coupled_lumped_phase_pi_sections,
       deck_cascaded_pi_request_rows,
       deck_cascaded_phase_pi_equivalents,
       deck_sampled_frequency_line_rows,
       deck_sampled_frequency_line_coefficients,
       deck_sampled_frequency_line_elements,
       deck_sampled_frequency_line_element_names,
       deck_semlyen_line_groups,
       deck_rational_frequency_line_groups,
       deck_rational_frequency_line_elements,
       deck_rational_frequency_line_element_names,
       deck_semlyen_line_elements,
       deck_semlyen_line_element_names,
       deck_generator_equivalent_rows,
       deck_distributed_transposed_line_constants,
       deck_distributed_transposed_line_modal_branch_states,
       deck_distributed_transposed_line_steady_state_pi_equivalents,
       deck_distributed_transposed_line_history_states,
       deck_distributed_transposed_line_companion_admittances,
       deck_line_constants_conductor_cards,
       deck_line_constants_frequency_cards,
       deck_line_constants_physical_conductors,
       deck_cable_constants_cases,
       deck_power_frequency_request_rows,
       deck_universal_machine_dimension_request_rows,
       deck_universal_machine_section_rows,
       deck_universal_machine_definition_rows,
       deck_universal_machine_coil_rows,
       deck_universal_machine_terminal_rows,
       deck_universal_machine_generated_branch_rows,
       deck_universal_machine_speed_capacitor_rows,
       deck_universal_machine_node_summary_rows,
       deck_universal_machine_output_summary_rows,
       deck_synchronous_machine_terminal_voltage_rows,
       deck_synchronous_machine_tolerance_rows,
       deck_synchronous_machine_parameter_fitting_rows,
       deck_synchronous_machine_model_parameter_rows,
       deck_synchronous_machine_mass_rows,
       deck_synchronous_machine_output_request_rows,
       deck_synchronous_machine_control_interface_rows,
       deck_synchronous_machine_output_summary_rows,
       deck_tacs_dimension_request_rows,
       deck_output_width_request_rows,
       deck_peak_voltage_monitor_request_rows,
       deck_diagnostic_print_request_rows,
       deck_tacs_warning_limit_request_rows,
       deck_plot_file_request_rows,
       deck_switch_logic_request_rows,
       deck_simulation_control_request_rows,
       deck_printout_frequency_change_rows,
       deck_study_option_request_rows,
       deck_fixed_source_constraint_rows,
       deck_fixed_source_control_rows,
       deck_time_horizon_control_rows,
       deck_output_schedule_control_rows,
       deck_case_boundary_rows,
       parse_emt_restart_request,
       deck_model_summary,
       deck_nodal_system,
       deck_node_names,
       deck_time_switch_close_time_s_values,
       deck_time_switch_count,
       deck_time_switch_from_node_indices,
       deck_time_switch_from_node_names,
       deck_time_switch_initially_closed_flags,
       deck_time_switch_line_numbers,
       deck_time_switch_off_conductance_values,
       deck_time_switch_on_conductance_values,
       deck_time_switch_names,
       deck_time_switch_open_time_s_values,
       deck_time_switch_to_node_indices,
       deck_time_switch_to_node_names,
       deck_over16_branch_voltage_branch_indices,
       deck_over16_branch_voltage_branch_names,
       deck_over16_branch_voltage_output_line_numbers,
       deck_over16_branch_voltage_output_names,
       deck_over16_branch_current_branch_indices,
       deck_over16_branch_current_branch_names,
       deck_over16_branch_current_output_line_numbers,
       deck_over16_branch_current_output_names,
       deck_over2_branch_capacitance_values,
       deck_over2_branch_conductance_values,
       deck_over2_branch_from_node_indices,
       deck_over2_branch_from_node_names,
       deck_over2_branch_inductance_values,
       deck_over2_branch_kinds,
       deck_over2_branch_layout_kinds,
       deck_over2_branch_line_numbers,
       deck_over2_branch_names,
       deck_over2_branch_output_codes,
       deck_over2_branch_raw_capacitance_values,
       deck_over2_branch_raw_inductance_values,
       deck_over2_branch_raw_resistance_values,
       deck_over2_branch_reference_kinds,
       deck_over2_branch_reference_line_numbers,
       deck_over2_branch_reference_names,
       deck_over2_branch_resistance_values,
       deck_over2_branch_rows,
       deck_over2_branch_source_kinds,
       deck_over2_branch_to_node_indices,
       deck_over2_branch_to_node_names,
       deck_arrester_constant_rows,
       deck_arrester_nonlinear_rows,
       deck_nonlinear_resistance_initialization_rows,
       deck_nonlinear_resistance_point_rows,
       deck_nonlinear_resistance_rows,
       deck_triggered_timed_resistance_point_rows,
       deck_triggered_timed_resistance_rows,
       deck_switching_nonlinear_resistor_point_rows,
       deck_switching_nonlinear_resistor_rows,
       deck_piecewise_nonlinear_inductor_point_rows,
       deck_piecewise_nonlinear_inductor_rows,
       deck_hysteretic_inductor_point_rows,
       deck_hysteretic_inductor_rows,
       deck_zinc_oxide_breakpoint_rows,
       deck_zinc_oxide_initialization_rows,
       deck_zinc_oxide_nonlinear_rows,
       deck_over5a_source_crest_values,
       deck_over5a_source_iform_values,
       deck_over5a_source_layout_kinds,
       deck_over5a_source_line_numbers,
       deck_over5a_source_names,
       deck_over5a_source_node_names,
       deck_over5a_source_node_values,
       deck_over5a_source_rows,
       deck_over5a_source_sfreq_values,
       deck_over5a_source_time1_values,
       deck_over5a_source_time2_values,
       deck_over5a_source_tstart_values,
       deck_over5a_source_tstop_values,
       deck_over5a_source_update_inputs,
       deck_over5_switch_closed_markers,
       deck_over5_switch_critical_current_values,
       deck_over5_switch_random_opening_standard_deviation_s_values,
       deck_over5_switch_type_values,
       deck_over5_switch_close_time_s_values,
       deck_over5_switch_from_node_indices,
       deck_over5_switch_from_node_names,
       deck_over5_switch_initially_closed_flags,
       deck_over5_switch_layout_kinds,
       deck_over5_switch_line_numbers,
       deck_over5_switch_marker_texts,
       deck_over5_switch_measuring_flags,
       deck_over5_switch_names,
       deck_over5_switch_off_conductance_values,
       deck_over5_switch_on_conductance_values,
       deck_over5_switch_open_time_s_values,
       deck_over5_switch_output_codes,
       deck_over5_switch_raw_close_time_s_values,
       deck_over5_switch_raw_open_time_s_values,
       deck_over5_switch_rows,
       deck_over5_switch_to_node_indices,
       deck_over5_switch_to_node_names,
       deck_over15_output_request_branch_indices,
       deck_over15_output_request_branch_names,
       deck_over15_output_request_layout_kinds,
       deck_over15_output_request_line_numbers,
       deck_over15_output_request_names,
       deck_over15_output_request_node_indices,
       deck_over15_output_request_node_names,
       deck_over15_output_request_output_codes,
       deck_over15_output_request_output_kinds,
       deck_over15_output_request_request_kinds,
       deck_over15_output_request_rows,
       deck_over16_branch_power_branch_indices,
       deck_over16_branch_power_branch_names,
       deck_over16_branch_power_output_line_numbers,
       deck_over16_branch_power_output_names,
       deck_over16_source_card_kinds,
       deck_over16_source_card_rows,
       deck_over16_source_card_values,
       deck_over16_source_interpolation_rows,
       deck_over16_source_interpolation_values,
       deck_over16_source_tacs_override_rows,
       deck_over16_source_tacs_override_positions,
       deck_over16_source_tacs_override_xtcs_indices,
       deck_over16_source_analytic_rows,
       deck_over16_source_analytic_values,
       deck_over16_source_analytic_line_numbers,
       deck_over16_source_analytic_provided_value_counts,
       deck_over16_output_channel_line_numbers,
       deck_over16_output_channel_names,
       deck_over16_output_node_indices,
       deck_over16_output_node_names,
       deck_over16_source_card_line_numbers,
       deck_over16_source_card_provided_value_counts,
       deck_over16_source_interpolation_line_numbers,
       deck_over16_source_interpolation_provided_value_counts,
       deck_over16_source_tacs_override_line_numbers,
       node_count,
       parse_deck_file,
       parse_deck_lines,
       parse_deck_case_sequence,
       parse_deck_file_sequence,
       parse_saturated_transformer_branch_section_intake_file,
       parse_saturated_transformer_branch_section_intake_lines,
       parse_saturated_transformer_branch_section_shunt_capacitance_rows,
       parse_saturated_transformer_intake_lines,
       parse_saturated_transformer_intake_file

struct DeckOVER16OutputChannel
    name::Symbol
    node::Symbol
    line_no::Int
end

struct DeckOVER16BranchVoltageOutput
    name::Symbol
    branch::Symbol
    line_no::Int
end

struct DeckOVER16BranchCurrentOutput
    name::Symbol
    branch::Symbol
    line_no::Int
end

struct DeckOVER16BranchPowerOutput
    name::Symbol
    branch::Symbol
    line_no::Int
end

struct DeckOVER15OutputRequestRow
    name::Symbol
    output_kind::Symbol
    request_kind::Symbol
    layout_kind::Symbol
    line_no::Int
    output_code::Int
    node::Symbol
    node_value::Int
    branch::Symbol
    branch_value::Int
end

struct DeckOVER16SourceCardRow
    kind::Symbol
    values::Vector{Float64}
    provided_value_count::Int
    line_no::Int
end

struct DeckOVER16SourceInterpolationRow
    values::Vector{Float64}
    provided_value_count::Int
    line_no::Int
end

struct DeckOVER16SourceTACSOverrideRow
    position::Int
    xtcs_index::Int
    line_no::Int
end

struct DeckOVER16SourceAnalyticRow
    values::Vector{Float64}
    provided_value_count::Int
    line_no::Int
end

struct DeckSaturatedTransformerHeaderRow
    name::Symbol
    reference_name::Symbol
    line_no::Int
    initial_current::Union{Missing,Float64}
    initial_flux::Union{Missing,Float64}
    magnetizing_resistance::Union{Missing,Float64}
end

struct DeckSaturatedTransformerBreakpointRow
    transformer_name::Symbol
    line_no::Int
    current::Float64
    flux::Float64
end

struct DeckSaturatedTransformerWindingRow
    transformer_name::Symbol
    line_no::Int
    winding_number::Int
    from_node::Symbol
    to_node::Symbol
    resistance::Union{Missing,Float64}
    inductance::Union{Missing,Float64}
    turns::Union{Missing,Float64}
    inherited_parameters::Bool
end

struct DeckSaturatedTransformerCopyRow
    reference_name::Symbol
    target_name::Symbol
    line_no::Int
end

struct DeckTransformerBranchShuntCapacitanceRow
    from_node::Symbol
    line_no::Int
    capacitance::Float64
end

struct DeckSaturatedTransformerIntake
    source::String
    transformers::Vector{DeckSaturatedTransformerHeaderRow}
    breakpoints::Vector{DeckSaturatedTransformerBreakpointRow}
    windings::Vector{DeckSaturatedTransformerWindingRow}
    copies::Vector{DeckSaturatedTransformerCopyRow}
    card_counts::Dict{Symbol,Int}
    validation::ValidationResult
end

struct DeckBergeronLineRow
    name::Symbol
    from_node::Symbol
    to_node::Symbol
    from_node_value::Int
    to_node_value::Int
    line_no::Int
    surge_impedance::Float64
    surge_admittance::Float64
    travel_time_s::Float64
    dt_s::Float64
    attenuation::Float64
    delay_steps::Int
    write_index::Int
    history_current_from::Float64
    history_current_to::Float64
    terminal_voltage_from::Float64
    terminal_voltage_to::Float64
    terminal_current_from::Float64
    terminal_current_to::Float64
    traveling_wave_from::Float64
    traveling_wave_to::Float64
end

struct DeckCoupledLineRow
    name::Symbol
    from_node::Symbol
    to_node::Symbol
    from_node_value::Int
    to_node_value::Int
    reference_from_node::Union{Missing,Symbol}
    reference_to_node::Union{Missing,Symbol}
    reference_from_node_value::Union{Missing,Int}
    reference_to_node_value::Union{Missing,Int}
    line_no::Int
    line_type::Int
    line_kind::Symbol
    phase_index::Int
    sequence_resistance::Union{Missing,Float64}
    sequence_inductance::Union{Missing,Float64}
    triangular_resistance_values::Vector{Float64}
    triangular_inductance_values::Vector{Float64}
    raw_resistance::Union{Missing,Float64}
    raw_inductance::Union{Missing,Float64}
    raw_capacitance::Union{Missing,Float64}
    line_length::Union{Missing,Float64}
    sampled_frequency_data_requested::Bool
end

struct DeckLineModalTransformRow
    line_no::Int
    values::Vector{Float64}
    raw_text::String
end

struct DeckCoupledPhasePiSectionRow
    name::Symbol
    from_node::Symbol
    to_node::Symbol
    from_node_value::Int
    to_node_value::Int
    line_no::Int
    branch_type::Int
    phase_index::Int
    reference_kind::Symbol
    reference_from_node::Union{Missing,Symbol}
    reference_to_node::Union{Missing,Symbol}
    reference_from_node_value::Union{Missing,Int}
    reference_to_node_value::Union{Missing,Int}
    raw_resistance_values::Vector{Float64}
    raw_inductance_values::Vector{Float64}
    raw_capacitance_values::Vector{Float64}
    continuation_line_numbers::Vector{Int}
    output_code::Int
end

struct DeckCascadedPiSeriesImpedanceRow
    line_no::Int
    phase_index::Int
    raw_resistance_ohm::Float64
    raw_inductance_value::Float64
    raw_capacitance_value::Float64
    open_circuit::Bool
end

struct DeckCascadedPiShuntImpedanceRow
    line_no::Int
    from_terminal::Int
    to_terminal::Int
    raw_resistance_ohm::Float64
    raw_inductance_value::Float64
    raw_capacitance_value::Float64
end

struct DeckCascadedPiBlock
    configuration_line_no::Int
    section_scale::Float64
    multiplicity::Int
    phase_map::Vector{Int}
    series_impedances::Vector{DeckCascadedPiSeriesImpedanceRow}
    shunt_impedances::Vector{DeckCascadedPiShuntImpedanceRow}
    explicit_resistance_values::Union{Nothing,Matrix{Float64}}
    explicit_inductance_values::Union{Nothing,Matrix{Float64}}
    explicit_capacitance_values::Union{Nothing,Matrix{Float64}}
    detail_line_numbers::Vector{Int}
end

struct DeckCascadedPiRequestRow
    name::Symbol
    header_line_no::Int
    configuration_line_no::Int
    terminator_line_no::Int
    phase_count::Int
    frequency_hz::Float64
    section_scale::Float64
    section_count::Int
    phase_map::Vector{Int}
    source_section_row_indices::Vector{Int}
    blocks::Vector{DeckCascadedPiBlock}
end

struct DeckSampledFrequencyLineRow
    name::Symbol
    branch_row_index::Int
    branch_line_no::Int
    configuration_line_no::Int
    sample_line_numbers::Vector{Int}
    from_node::Symbol
    to_node::Symbol
    from_node_index::Int
    to_node_index::Int
    propagation_peak_index::Int
    admittance_rise_index::Int
    characteristic_impedance_ohm::Float64
    propagation_cutoff_fraction::Float64
    admittance_cutoff_fraction::Float64
    total_resistance_ohm::Float64
    maximum_tail_iterations::Int
    propagation_time_s::Vector{Float64}
    propagation_amplitude::Vector{Float64}
    admittance_time_s::Vector{Float64}
    admittance_amplitude::Vector{Float64}
end

struct DeckSemlyenModeRow
    header_line_no::Int
    detail_line_numbers::Vector{Int}
    from_node::Symbol
    to_node::Symbol
    from_node_index::Int
    to_node_index::Int
    mode_index::Int
    parameters::SemlyenModeParameters
end

struct DeckSemlyenLineGroupRow
    name::Symbol
    phase_count::Int
    modes::Vector{DeckSemlyenModeRow}
    voltage_modal_to_phase::Matrix{ComplexF64}
    current_modal_to_phase::Matrix{ComplexF64}
    transform_line_numbers::Vector{Int}
    imaginary_transform_parts_supplied::Bool
end

struct DeckRationalLineModeRow
    header_line_no::Int
    detail_line_numbers::Vector{Int}
    from_node::Symbol
    to_node::Symbol
    from_node_index::Int
    to_node_index::Int
    mode_index::Int
    conversion::RationalLineModeConversion
    characteristic_impedance_reduction::PoleResidueReductionResult
    propagation_reduction::PoleResidueReductionResult
end

struct DeckRationalLineGroupRow
    name::Symbol
    phase_count::Int
    modes::Vector{DeckRationalLineModeRow}
    voltage_modal_to_phase::Matrix{ComplexF64}
    current_modal_to_phase::Matrix{ComplexF64}
    transform_line_numbers::Vector{Int}
    source_format::Symbol
end

struct DeckGeneratorEquivalentModalBranchRow
    line_no::Int
    sequence_kind::Symbol
    raw_resistance::Float64
    raw_inductance::Float64
    raw_capacitance::Float64
    raw_damping_resistance::Float64
    branch::GeneratorEquivalentModalBranch
end

struct DeckGeneratorEquivalentRow
    name::Symbol
    header_line_numbers::Vector{Int}
    marker_line_no::Int
    from_node_names::Vector{Symbol}
    to_node_names::Vector{Symbol}
    from_node_indices::Vector{Int}
    to_node_indices::Vector{Int}
    zero_mode_branches::Vector{DeckGeneratorEquivalentModalBranchRow}
    positive_mode_branches::Vector{DeckGeneratorEquivalentModalBranchRow}
    output_codes::Vector{Int}
end

struct DeckLineConstantsConductorCard
    line_no::Int
    phase_number::Int
    skin_effect_type::Float64
    resistance_ohm_per_mile::Float64
    reactance_type::Int
    reactance_or_gmr::Float64
    diameter_inches::Float64
    horizontal_ft::Float64
    tower_height_ft::Float64
    midspan_height_ft::Float64
    average_height_ft::Float64
    bundle_spacing_inches::Float64
    bundle_angle_deg::Float64
    conductor_name::String
    bundle_conductor_count::Int
end

struct DeckLineConstantsFrequencyCard
    line_no::Int
    earth_resistivity_ohm_m::Float64
    frequency_hz::Float64
    carson_correction_factor::Float64
    capacitance_print_flags::NTuple{6,Int}
    impedance_print_flags::NTuple{6,Int}
    matrix_output_selector::Int
    distance_miles::Float64
    punch_request::Int
    alternate_punch_flags::NTuple{5,Int}
    frequency_decade_count::Int
    points_per_decade::Int
    line_model_punch_request::Int
    modal_output_flag::Int
    transform_output_flag::Int
    conductance_mho_per_mile::Float64
end

struct DeckLineConstantsPhysicalConductor
    source_card_index::Int
    line_no::Int
    bundle_ordinal::Int
    phase_number::Int
    skin_effect_type::Float64
    resistance_ohm_per_mile::Float64
    reactance_type::Int
    reactance_or_gmr::Float64
    diameter_inches::Float64
    horizontal_ft::Float64
    average_height_ft::Float64
    conductor_name::String
end

struct DeckCableConstantsFrequencyCard
    line_no::Int
    earth_resistivity_ohm_m::Float64
    start_frequency_hz::Float64
    decade_count::Int
    points_per_decade::Int
    distance_m::Float64
    card_output_flag::Int
    transform_flag::Int
end

struct DeckCableConstantsCase
    line_no::Int
    cable_kind_code::Int
    surface_position_code::Int
    phase_count::Int
    earth_model_code::Int
    modal_output_flag::Int
    impedance_output_flag::Int
    admittance_output_flag::Int
    pipe_count::Int
    grounding_selector::Int
    pipe_radii_m::Vector{Float64}
    pipe_resistivity_ohm_m::Float64
    pipe_relative_permeability::Float64
    pipe_inner_insulator_relative_permittivity::Float64
    pipe_outer_insulator_relative_permittivity::Float64
    cable_to_pipe_center_distances_m::Vector{Float64}
    cable_to_pipe_angles_rad::Vector{Float64}
    pipe_depths_m::Vector{Float64}
    pipe_horizontal_positions_m::Vector{Float64}
    layer_counts::Vector{Int}
    boundary_radii_m::Matrix{Float64}
    resistivity_ohm_m::Matrix{Float64}
    conductor_relative_permeability::Matrix{Float64}
    insulation_relative_permeability::Matrix{Float64}
    insulation_relative_permittivity::Matrix{Float64}
    depths_m::Vector{Float64}
    horizontal_positions_m::Vector{Float64}
    frequency_cards::Vector{DeckCableConstantsFrequencyCard}
end

struct DeckPowerFrequencyRequestRow
    line_no::Int
    frequency_hz::Float64
    raw_text::String
end

struct DeckUniversalMachineDimensionRequestRow
    line_no::Int
    coil_table_size::Int
    machine_table_size::Int
    output_table_size::Int
    bus_table_size::Int
    raw_text::String
end

struct DeckUniversalMachineSectionRow
    line_no::Int
    machine_count::Int
    input_layout::Symbol
    parameter_basis::Symbol
    remanent_flux_enabled::Bool
    initialization_mode::Symbol
    detailed_machine_input::Bool
    maximum_shaft_mass_count::Int
    terminal_coupling::Symbol
    raw_text::String
end

struct DeckUniversalMachineDefinitionRow
    line_no::Int
    machine_index::Int
    card_index::Int
    machine_type::Int
    d_axis_coil_count::Int
    q_axis_coil_count::Int
    torque_output_flag::Int
    speed_output_flag::Int
    angle_output_flag::Int
    pole_pair_count::Int
    rotor_mass::Union{Missing,Float64}
    saturation_mode::Int
    saturated_inductance::Union{Missing,Float64}
    saturation_flux::Union{Missing,Float64}
    remanent_flux::Union{Missing,Float64}
    value1::Union{Missing,Float64}
    value2::Union{Missing,Float64}
    mechanical_damping_coefficient::Union{Missing,Float64}
    speed_convergence_margin::Union{Missing,Float64}
    node::Union{Missing,Symbol}
    raw_text::String
end

struct DeckUniversalMachineCoilRow
    line_no::Int
    machine_index::Int
    coil_index::Int
    resistance::Float64
    inductance::Float64
    terminal_node::Union{Missing,Symbol}
    control_signal::Union{Missing,Symbol}
    output_flag::Int
    initial_history_current::Float64
    raw_text::String
end

struct DeckUniversalMachineTerminalRow
    line_no::Int
    machine_index::Int
    terminal_index::Int
    terminal_node::Symbol
    reference_node::Symbol
    terminal_node_value::Int
    reference_node_value::Int
    raw_text::String
end

struct DeckUniversalMachineGeneratedBranchRow
    line_no::Int
    machine_index::Int
    branch_index::Int
    from_node::Symbol
    to_node::Symbol
    from_node_value::Int
    to_node_value::Int
    reactance::Union{Missing,Float64}
    raw_text::String
end

struct DeckUniversalMachineSpeedCapacitorRow
    line_no::Int
    machine_index::Int
    capacitor_node::Symbol
    mass_node::Symbol
    capacitor_node_value::Int
    mass_node_value::Int
    resistance::Float64
    capacitance::Float64
    raw_text::String
end

struct DeckUniversalMachineNodeSummaryRow
    line_no::Int
    machine_index::Int
    mass_node::Symbol
    mass_node_value::Int
    mechanical_slack_source::Union{Missing,Symbol}
    mechanical_slack_source_index::Int
    field_slack_source::Union{Missing,Symbol}
    field_slack_source_index::Int
    raw_text::String
end

struct DeckUniversalMachineOutputSummaryRow
    line_no::Int
    machine_count::Int
    output_count::Int
    tacs_transfer_count::Int
    raw_text::String
end

struct DeckSynchronousMachineTerminalVoltageRow
    line_no::Int
    machine_index::Int
    source_group_index::Int
    source_type::Int
    phase_index::Int
    terminal_node::Symbol
    terminal_node_value::Int
    source_node_value::Int
    peak_terminal_voltage::Union{Missing,Float64}
    frequency_hz::Union{Missing,Float64}
    angle_deg::Union{Missing,Float64}
    raw_text::String
end

struct DeckSynchronousMachineToleranceRow
    line_no::Int
    machine_index::Int
    values::Vector{Float64}
    raw_text::String
end

struct DeckSynchronousMachineParameterFittingRow
    line_no::Int
    machine_index::Int
    value::Union{Missing,Float64}
    raw_text::String
end

struct DeckSynchronousMachineModelParameterRow
    line_no::Int
    machine_index::Int
    parameter_kind::Symbol
    values::Vector{Float64}
    positional_values::Vector{Union{Missing,Float64}}
    raw_text::String
end

struct DeckSynchronousMachineMassRow
    line_no::Int
    machine_index::Int
    mass_index::Int
    values::Vector{Float64}
    torque_fraction::Float64
    inertia::Float64
    speed_deviation_damping::Float64
    mutual_damping::Float64
    shaft_stiffness::Float64
    absolute_speed_damping::Float64
    raw_text::String
end

struct DeckSynchronousMachineOutputRequestRow
    line_no::Int
    machine_index::Int
    group_index::Int
    output_codes::Vector{Int}
    dynamic_output_count::Int
    raw_text::String
end

struct DeckSynchronousMachineControlInterfaceRow
    line_no::Int
    machine_index::Int
    interface_code::Int
    direction::Symbol
    coupling_kind::Symbol
    signal_name::Symbol
    variable_index::Int
    raw_text::String
end

struct DeckSynchronousMachineOutputSummaryRow
    line_no::Int
    machine_count::Int
    terminal_count::Int
    mass_count::Int
    output_count::Int
    raw_text::String
end

struct DeckTACSDimensionRequestRow
    line_no::Int
    payload_line_no::Int
    allocation_kind::Symbol
    request_values::NTuple{8,Float64}
    list_sizes::NTuple{8,Int}
    provided_value_count::Int
    request_text::String
    payload_text::String
end

struct DeckOutputWidthRequestRow
    line_no::Int
    column_width::Int
    raw_text::String
end

struct DeckPeakVoltageMonitorRequestRow
    line_no::Int
    raw_text::String
end

struct DeckDiagnosticPrintRequestRow
    line_no::Int
    loop_print_controls::NTuple{4,Int}
    raw_text::String
end

struct DeckTACSWarningLimitRequestRow
    line_no::Int
    warning_limit::Int
    begin_time_s::Float64
    raw_text::String
end

struct DeckPlotFileRequestRow
    line_no::Int
    plot_file_mode::Int
    raw_text::String
end

struct DeckSwitchLogicRequestRow
    line_no::Int
    control_value::Int
    raw_text::String
end

struct DeckSimulationControlRequestRow
    line_no::Int
    request_kind::Symbol
    numeric_value::Union{Missing,Float64}
    integer_value::Union{Missing,Int}
    raw_text::String
end

struct DeckPrintoutFrequencyChangeRow
    request_line_no::Int
    payload_line_no::Int
    change_steps::NTuple{5,Int}
    multipliers::NTuple{5,Int}
    active_pair_count::Int
    provided_value_count::Int
    request_text::String
    payload_text::String
end

struct DeckStudyOptionRequestRow
    line_no::Int
    request_kind::Symbol
    numeric_values::Vector{Float64}
    integer_values::Vector{Int}
    text_values::Vector{String}
    raw_text::String
end

struct DeckFixedSourceConstraintRow
    line_no::Int
    constraint_kind::Symbol
    source_node_names::NTuple{3,Symbol}
    active_power::Union{Missing,Float64}
    reactive_power::Union{Missing,Float64}
    voltage_peak::Union{Missing,Float64}
    angle_deg::Union{Missing,Float64}
    minimum_voltage::Float64
    maximum_voltage::Float64
    minimum_angle_deg::Float64
    maximum_angle_deg::Float64
    raw_text::String
end

struct DeckFixedSourceControlRow
    line_no::Int
    print_voltage_changes::Bool
    maximum_iterations::Int
    voltage_change_report_interval::Int
    print_final_sources::Bool
    relative_power_tolerance::Float64
    voltage_correction_factor::Float64
    angle_correction_factor::Float64
    raw_text::String
end

struct DeckTimeHorizonControlRow
    line_no::Int
    time_step_s::Float64
    stop_time_s::Float64
    inductance_frequency_hz::Float64
    capacitance_frequency_hz::Float64
    epsilon::Float64
    matrix_tolerance::Float64
    start_time_s::Float64
    raw_text::String
end

struct DeckOutputScheduleControlRow
    line_no::Int
    print_interval_steps::Int
    plot_interval_steps::Int
    network_print_enabled::Bool
    steady_state_print_enabled::Bool
    extrema_print_enabled::Bool
    terminal_conditions_punch_enabled::Bool
    restart_snapshot_enabled::Bool
    plot_file_retention_mode::Int
    energization_count::Int
    special_request_control_word::Int
    raw_text::String
end

struct DeckCaseBoundaryRow
    line_no::Int
    boundary_kind::Symbol
    raw_text::String
end

struct DeckControlSystemSignalTerm
    name::Symbol
    polarity::Int
end

mutable struct DeckControlSystemFunctionRow
    line_no::Int
    name::Symbol
    input_terms::Vector{DeckControlSystemSignalTerm}
    order::Int
    gain::Float64
    numerator_coefficients::Vector{Float64}
    denominator_coefficients::Vector{Float64}
    lower_limit::Union{Missing,Float64}
    upper_limit::Union{Missing,Float64}
    lower_limit_signal::Union{Missing,Symbol}
    upper_limit_signal::Union{Missing,Symbol}
    coefficients_complete::Bool
    raw_text::String
end

struct DeckControlSystemExpressionRow
    line_no::Int
    group_type::Int
    name::Symbol
    expression::String
    raw_text::String
end

struct DeckControlSystemSourceRow
    line_no::Int
    source_type::Int
    name::Symbol
    amplitude::Union{Missing,Float64}
    delay_or_time_constant::Union{Missing,Float64}
    phase_or_width::Union{Missing,Float64}
    activation_start_time_s::Float64
    activation_stop_time_s::Float64
    numeric_values::Vector{Float64}
    raw_text::String
end

mutable struct DeckControlSystemDeviceRow
    line_no::Int
    group_type::Int
    name::Symbol
    device_type::Int
    input_terms::Vector{DeckControlSystemSignalTerm}
    first_input::Union{Missing,Symbol}
    tail_signal_names::Vector{Symbol}
    control_signal::Union{Missing,Symbol}
    reference_signal::Union{Missing,Symbol}
    parameter_values::Vector{Float64}
    table_input_values::Vector{Float64}
    table_output_values::Vector{Float64}
    table_complete::Bool
    raw_text::String
end

struct DeckControlSystemOutputRequestRow
    line_no::Int
    request_type::Int
    all_signals::Bool
    signal_names::Vector{Symbol}
    raw_text::String
end

struct DeckDelayedArcSwitchParameters
    line_no::Int
    current_coefficient::Float64
    current_exponent::Float64
    time_scale_s::Float64
    cutoff_current_a::Float64
    raw_text::String
end

struct DeckControlSystemSwitchRow
    line_no::Int
    switch_type::Int
    from_node::Symbol
    to_node::Symbol
    ignition_voltage::Union{Missing,Float64}
    holding_current::Union{Missing,Float64}
    deionization_time_s::Union{Missing,Float64}
    initial_state::String
    control_signal::Symbol
    gate_signal::Union{Missing,Symbol}
    clamp_signal::Union{Missing,Symbol}
    switching_delay_or_model::Union{Missing,Float64}
    delayed_arc::Union{Nothing,DeckDelayedArcSwitchParameters}
    parameter_source_kind::Symbol
    parameter_reference_index::Int
    parameter_reference_line_no::Int
    layout_kind::Symbol
    event_output_code::Int
    output_code::Int
    raw_text::String
end

struct DeckControlSystemSwitchCouplingRow
    line_no::Int
    switch_type::Int
    from_node::Symbol
    to_node::Symbol
    from_node_index::Union{Missing,Int}
    to_node_index::Union{Missing,Int}
    control_signal::Symbol
    gate_signal::Union{Missing,Symbol}
    clamp_signal::Union{Missing,Symbol}
    control_output_request_index::Union{Missing,Int}
    control_output_signal_index::Union{Missing,Int}
    control_output_linear_index::Union{Missing,Int}
    output_code::Int
    initial_state::String
end

struct DeckNodeInitialConditionRow
    line_no::Int
    condition_kind::Symbol
    node::Symbol
    node_index::Int
    reference_node::Union{Missing,Symbol}
    reference_node_index::Union{Missing,Int}
    real_value::Float64
    imaginary_value::Float64
    raw_text::String
end

struct DeckControlCard
    kind::Symbol
    label::String
    line_no::Int
    tokens::Vector{String}
end

struct DeckOVER2BranchRow
    name::Symbol
    from_node::Symbol
    to_node::Symbol
    from_node_value::Int
    to_node_value::Int
    line_no::Int
    branch_type::Int
    branch_kind::Symbol
    layout_kind::Symbol
    source_kind::Symbol
    reference_kind::Symbol
    reference_name::Symbol
    reference_line_no::Int
    raw_resistance::Float64
    raw_inductance::Float64
    raw_capacitance::Float64
    conductance::Float64
    resistance::Float64
    inductance::Float64
    capacitance::Float64
    output_code::Int
end

struct DeckZincOxideNonlinearRow
    name::Symbol
    from_node::Symbol
    to_node::Symbol
    from_node_index::Int
    to_node_index::Int
    line_no::Int
    nonlinear_type::Int
    raw_resistance::Float64
    raw_inductance::Float64
    raw_capacitance_marker::Float64
    output_code::Int
    first_characteristic_index::Int
    source_kind::Symbol
    reference_kind::Symbol
    reference_index::Int
    reference_name::Symbol
    reference_line_no::Int
    raw_text::String
end

struct DeckZincOxideInitializationRow
    nonlinear_row_index::Int
    line_no::Int
    reference_voltage::Float64
    gap_voltage::Float64
    initial_voltage::Float64
    raw_text::String
end

struct DeckZincOxideBreakpointRow
    nonlinear_row_index::Int
    breakpoint_index::Int
    line_no::Int
    current_coefficient::Float64
    voltage_exponent::Float64
    voltage::Float64
    raw_text::String
end

struct DeckNonlinearResistanceRow
    name::Symbol
    from_node::Symbol
    to_node::Symbol
    from_node_index::Int
    to_node_index::Int
    line_no::Int
    nonlinear_type::Int
    element_kind::Symbol
    steady_state_reference::Float64
    secondary_reference::Float64
    table_marker::Float64
    output_code::Int
    first_characteristic_index::Int
    source_kind::Symbol
    reference_kind::Symbol
    reference_index::Int
    reference_name::Symbol
    reference_line_no::Int
    raw_text::String
end

struct DeckNonlinearResistanceInitializationRow
    nonlinear_row_index::Int
    line_no::Int
    table_voltage_offset::Float64
    gap_voltage::Float64
    initial_voltage::Float64
    raw_text::String
end

struct DeckNonlinearResistancePointRow
    nonlinear_row_index::Int
    point_index::Int
    line_no::Int
    ordinate_value::Float64
    coordinate_value::Float64
    raw_text::String
end

"""A passive resistance schedule armed by time and triggered by branch voltage."""
struct DeckTriggeredTimedResistanceRow
    name::Symbol
    from_node::Symbol
    to_node::Symbol
    from_node_index::Int
    to_node_index::Int
    line_no::Int
    trigger_voltage_v::Float64
    arm_time_s::Float64
    output_code::Int
    source_kind::Symbol
    reference_kind::Symbol
    reference_index::Int
    reference_name::Symbol
    reference_line_no::Int
    raw_text::String
end

struct DeckTriggeredTimedResistancePointRow
    resistance_row_index::Int
    point_index::Int
    line_no::Int
    elapsed_time_s::Float64
    resistance_ohm::Float64
    raw_text::String
end

"""A symmetric piecewise-linear resistor with voltage-triggered switching."""
struct DeckSwitchingNonlinearResistorRow
    name::Symbol
    from_node::Symbol
    to_node::Symbol
    from_node_index::Int
    to_node_index::Int
    line_no::Int
    turn_on_voltage::Float64
    minimum_on_time_s::Float64
    activation_segment_count::Int
    turn_off_voltage::Float64
    output_code::Int
    single_flash::Bool
    source_kind::Symbol
    reference_kind::Symbol
    reference_index::Int
    reference_name::Symbol
    reference_line_no::Int
    raw_text::String
end

struct DeckSwitchingNonlinearResistorPointRow
    resistor_row_index::Int
    point_index::Int
    line_no::Int
    current_a::Float64
    voltage_v::Float64
    raw_text::String
end

"""A passive piecewise-linear flux-current inductor with explicit initial state."""
struct DeckPiecewiseNonlinearInductorRow
    name::Symbol
    from_node::Symbol
    to_node::Symbol
    from_node_index::Int
    to_node_index::Int
    line_no::Int
    steady_state_current_a::Float64
    steady_state_flux_wb::Float64
    output_code::Int
    source_kind::Symbol
    reference_kind::Symbol
    reference_index::Int
    reference_name::Symbol
    reference_line_no::Int
    raw_text::String
end

struct DeckPiecewiseNonlinearInductorPointRow
    inductor_row_index::Int
    point_index::Int
    line_no::Int
    current_a::Float64
    flux_wb::Float64
    raw_text::String
end

struct DeckHystereticInductorRow
    name::Symbol
    from_node::Symbol
    to_node::Symbol
    from_node_index::Int
    to_node_index::Int
    line_no::Int
    steady_state_flux_Wb::Float64
    steady_state_current_A::Float64
    residual_flux_Wb::Float64
    output_code::Int
    first_point_index::Int
    source_kind::Symbol
    reference_kind::Symbol
    reference_index::Int
    reference_name::Symbol
    reference_line_no::Int
    raw_text::String
end

struct DeckHystereticInductorPointRow
    hysteretic_inductor_row_index::Int
    point_index::Int
    line_no::Int
    current_A::Float64
    flux_Wb::Float64
    raw_text::String
end

struct DeckArresterNonlinearRow
    name::Symbol
    from_node::Symbol
    to_node::Symbol
    from_node_index::Int
    to_node_index::Int
    line_no::Int
    nonlinear_type::Int
    flashover_voltage::Float64
    voltage_division_factor::Float64
    current_division_factor::Float64
    output_code::Int
    first_constant_index::Int
    source_kind::Symbol
    reference_kind::Symbol
    reference_index::Int
    reference_name::Symbol
    reference_line_no::Int
    raw_text::String
end

struct DeckArresterConstantRow
    nonlinear_row_index::Int
    first_constant_index::Int
    line_no::Int
    values::Vector{Float64}
    raw_text::String
end

struct DeckOVER5ASourceRow
    name::Symbol
    node::Symbol
    node_value::Int
    line_no::Int
    layout_kind::Symbol
    iform::Int
    crest::Float64
    time1::Float64
    time2::Float64
    sfreq::Float64
    tstart::Float64
    tstop::Float64
end

struct DeckOVER5SwitchRow
    name::Symbol
    from_node::Symbol
    to_node::Symbol
    from_node_value::Int
    to_node_value::Int
    line_no::Int
    switch_type::Int
    output_code::Int
    layout_kind::Symbol
    raw_close_time_s::Float64
    raw_open_time_s::Float64
    close_time_s::Float64
    open_time_s::Float64
    initially_closed::Bool
    measuring::Bool
    closed_marker::String
    marker_text::String
    critical_current_a::Float64
    random_opening_standard_deviation_s::Float64
    on_conductance::Float64
    off_conductance::Float64
end

struct DeckOVER5ASourceUpdateInputs
    names::Vector{Symbol}
    nodes::Vector{Symbol}
    node_values::Vector{Int}
    iform_values::Vector{Int}
    line_numbers::Vector{Int}
    layout_kinds::Vector{Symbol}
    tstart_values::Vector{Float64}
    tstop_values::Vector{Float64}
    crest_values::Vector{Float64}
    time1_values::Vector{Float64}
    time2_values::Vector{Float64}
    sfreq_values::Vector{Float64}
    kconst::Int
end

struct DeckParseResult
    source::String
    elements::Vector{Any}
    element_line_numbers::Vector{Int}
    node_map::Dict{Symbol,Int}
    element_names::Vector{Symbol}
    bergeron_line_rows::Vector{DeckBergeronLineRow}
    coupled_line_rows::Vector{DeckCoupledLineRow}
    line_modal_transform_rows::Vector{DeckLineModalTransformRow}
    coupled_phase_pi_section_rows::Vector{DeckCoupledPhasePiSectionRow}
    cascaded_pi_request_rows::Vector{DeckCascadedPiRequestRow}
    sampled_frequency_line_rows::Vector{DeckSampledFrequencyLineRow}
    semlyen_line_groups::Vector{DeckSemlyenLineGroupRow}
    rational_frequency_line_groups::Vector{DeckRationalLineGroupRow}
    generator_equivalent_rows::Vector{DeckGeneratorEquivalentRow}
    line_constants_conductor_cards::Vector{DeckLineConstantsConductorCard}
    line_constants_frequency_cards::Vector{DeckLineConstantsFrequencyCard}
    cable_constants_cases::Vector{DeckCableConstantsCase}
    power_frequency_request_rows::Vector{DeckPowerFrequencyRequestRow}
    universal_machine_dimension_request_rows::Vector{DeckUniversalMachineDimensionRequestRow}
    universal_machine_section_rows::Vector{DeckUniversalMachineSectionRow}
    universal_machine_definition_rows::Vector{DeckUniversalMachineDefinitionRow}
    universal_machine_coil_rows::Vector{DeckUniversalMachineCoilRow}
    universal_machine_terminal_rows::Vector{DeckUniversalMachineTerminalRow}
    universal_machine_generated_branch_rows::Vector{DeckUniversalMachineGeneratedBranchRow}
    universal_machine_speed_capacitor_rows::Vector{DeckUniversalMachineSpeedCapacitorRow}
    universal_machine_node_summary_rows::Vector{DeckUniversalMachineNodeSummaryRow}
    universal_machine_output_summary_rows::Vector{DeckUniversalMachineOutputSummaryRow}
    synchronous_machine_terminal_voltage_rows::Vector{DeckSynchronousMachineTerminalVoltageRow}
    synchronous_machine_tolerance_rows::Vector{DeckSynchronousMachineToleranceRow}
    synchronous_machine_parameter_fitting_rows::Vector{DeckSynchronousMachineParameterFittingRow}
    synchronous_machine_model_parameter_rows::Vector{DeckSynchronousMachineModelParameterRow}
    synchronous_machine_mass_rows::Vector{DeckSynchronousMachineMassRow}
    synchronous_machine_output_request_rows::Vector{DeckSynchronousMachineOutputRequestRow}
    synchronous_machine_control_interface_rows::Vector{DeckSynchronousMachineControlInterfaceRow}
    synchronous_machine_output_summary_rows::Vector{DeckSynchronousMachineOutputSummaryRow}
    tacs_dimension_request_rows::Vector{DeckTACSDimensionRequestRow}
    output_width_request_rows::Vector{DeckOutputWidthRequestRow}
    peak_voltage_monitor_request_rows::Vector{DeckPeakVoltageMonitorRequestRow}
    diagnostic_print_request_rows::Vector{DeckDiagnosticPrintRequestRow}
    tacs_warning_limit_request_rows::Vector{DeckTACSWarningLimitRequestRow}
    plot_file_request_rows::Vector{DeckPlotFileRequestRow}
    switch_logic_request_rows::Vector{DeckSwitchLogicRequestRow}
    simulation_control_request_rows::Vector{DeckSimulationControlRequestRow}
    printout_frequency_change_rows::Vector{DeckPrintoutFrequencyChangeRow}
    study_option_request_rows::Vector{DeckStudyOptionRequestRow}
    fixed_source_constraint_rows::Vector{DeckFixedSourceConstraintRow}
    fixed_source_control_rows::Vector{DeckFixedSourceControlRow}
    time_horizon_control_rows::Vector{DeckTimeHorizonControlRow}
    output_schedule_control_rows::Vector{DeckOutputScheduleControlRow}
    case_boundary_rows::Vector{DeckCaseBoundaryRow}
    control_system_function_rows::Vector{DeckControlSystemFunctionRow}
    control_system_expression_rows::Vector{DeckControlSystemExpressionRow}
    control_system_source_rows::Vector{DeckControlSystemSourceRow}
    control_system_device_rows::Vector{DeckControlSystemDeviceRow}
    control_system_output_request_rows::Vector{DeckControlSystemOutputRequestRow}
    control_system_switch_rows::Vector{DeckControlSystemSwitchRow}
    node_initial_condition_rows::Vector{DeckNodeInitialConditionRow}
    over2_branch_rows::Vector{DeckOVER2BranchRow}
    zinc_oxide_nonlinear_rows::Vector{DeckZincOxideNonlinearRow}
    zinc_oxide_initialization_rows::Vector{DeckZincOxideInitializationRow}
    zinc_oxide_breakpoint_rows::Vector{DeckZincOxideBreakpointRow}
    nonlinear_resistance_rows::Vector{DeckNonlinearResistanceRow}
    nonlinear_resistance_initialization_rows::Vector{DeckNonlinearResistanceInitializationRow}
    nonlinear_resistance_point_rows::Vector{DeckNonlinearResistancePointRow}
    triggered_timed_resistance_rows::Vector{DeckTriggeredTimedResistanceRow}
    triggered_timed_resistance_point_rows::Vector{DeckTriggeredTimedResistancePointRow}
    switching_nonlinear_resistor_rows::Vector{DeckSwitchingNonlinearResistorRow}
    switching_nonlinear_resistor_point_rows::Vector{DeckSwitchingNonlinearResistorPointRow}
    piecewise_nonlinear_inductor_rows::Vector{DeckPiecewiseNonlinearInductorRow}
    piecewise_nonlinear_inductor_point_rows::Vector{DeckPiecewiseNonlinearInductorPointRow}
    hysteretic_inductor_rows::Vector{DeckHystereticInductorRow}
    hysteretic_inductor_point_rows::Vector{DeckHystereticInductorPointRow}
    arrester_nonlinear_rows::Vector{DeckArresterNonlinearRow}
    arrester_constant_rows::Vector{DeckArresterConstantRow}
    over5a_source_rows::Vector{DeckOVER5ASourceRow}
    over5_switch_rows::Vector{DeckOVER5SwitchRow}
    over15_output_request_rows::Vector{DeckOVER15OutputRequestRow}
    over16_output_channels::Vector{DeckOVER16OutputChannel}
    over16_branch_voltage_outputs::Vector{DeckOVER16BranchVoltageOutput}
    over16_branch_current_outputs::Vector{DeckOVER16BranchCurrentOutput}
    over16_branch_power_outputs::Vector{DeckOVER16BranchPowerOutput}
    over16_source_card_rows::Vector{DeckOVER16SourceCardRow}
    over16_source_interpolation_rows::Vector{DeckOVER16SourceInterpolationRow}
    over16_source_tacs_override_rows::Vector{DeckOVER16SourceTACSOverrideRow}
    over16_source_analytic_rows::Vector{DeckOVER16SourceAnalyticRow}
    control_cards::Vector{DeckControlCard}
    card_counts::Dict{Symbol,Int}
    pending_fixed_owner_names::Dict{Symbol,Vector{String}}
    validation::ValidationResult
end

struct DeckModelSummary
    source::String
    node_names::Vector{Symbol}
    element_names::Vector{Symbol}
    node_count::Int
    element_count::Int
    branch_names::Vector{Symbol}
    branch_kinds::Vector{Symbol}
    branch_from_node_names::Vector{Symbol}
    branch_to_node_names::Vector{Symbol}
    branch_from_node_indices::Vector{Int}
    branch_to_node_indices::Vector{Int}
    branch_conductance_values::Vector{Float64}
    branch_resistance_values::Vector{Float64}
    branch_inductance_values::Vector{Float64}
    branch_capacitance_values::Vector{Float64}
    branch_previous_current_values::Vector{Float64}
    branch_previous_voltage_values::Vector{Float64}
    branch_line_numbers::Vector{Int}
    branch_count::Int
    bergeron_line_names::Vector{Symbol}
    bergeron_line_line_numbers::Vector{Int}
    bergeron_line_from_node_names::Vector{Symbol}
    bergeron_line_to_node_names::Vector{Symbol}
    bergeron_line_from_node_indices::Vector{Int}
    bergeron_line_to_node_indices::Vector{Int}
    bergeron_line_surge_impedance_values::Vector{Float64}
    bergeron_line_surge_admittance_values::Vector{Float64}
    bergeron_line_travel_time_s_values::Vector{Float64}
    bergeron_line_dt_s_values::Vector{Float64}
    bergeron_line_attenuation_values::Vector{Float64}
    bergeron_line_delay_step_counts::Vector{Int}
    bergeron_line_write_indices::Vector{Int}
    bergeron_line_history_current_from_values::Vector{Float64}
    bergeron_line_history_current_to_values::Vector{Float64}
    bergeron_line_terminal_voltage_from_values::Vector{Float64}
    bergeron_line_terminal_voltage_to_values::Vector{Float64}
    bergeron_line_terminal_current_from_values::Vector{Float64}
    bergeron_line_terminal_current_to_values::Vector{Float64}
    bergeron_line_traveling_wave_from_values::Vector{Float64}
    bergeron_line_traveling_wave_to_values::Vector{Float64}
    bergeron_line_count::Int
    coupled_lumped_sequence_impedances::Vector{CoupledLumpedSequenceImpedance}
    coupled_lumped_sequence_count::Int
    coupled_lumped_phase_pi_sections::Vector{CoupledLumpedPhasePiSection}
    coupled_lumped_phase_pi_section_count::Int
    distributed_transposed_line_constants::Vector{DistributedTransposedLineConstants}
    distributed_transposed_line_count::Int
    distributed_transposed_line_modal_branch_states::Vector{DistributedTransposedLineModalBranchState}
    distributed_transposed_line_modal_branch_state_count::Int
    distributed_transposed_line_steady_state_pi_equivalents::Vector{DistributedTransposedLineSteadyStatePiEquivalent}
    distributed_transposed_line_steady_state_pi_equivalent_count::Int
    distributed_transposed_line_history_states::Vector{DistributedTransposedLineHistoryState}
    distributed_transposed_line_history_state_count::Int
    distributed_transposed_line_companion_admittances::Vector{DistributedTransposedLineCompanionAdmittance}
    distributed_transposed_line_companion_admittance_count::Int
    over2_branch_names::Vector{Symbol}
    over2_branch_line_numbers::Vector{Int}
    over2_branch_kinds::Vector{Symbol}
    over2_branch_layout_kinds::Vector{Symbol}
    over2_branch_source_kinds::Vector{Symbol}
    over2_branch_reference_kinds::Vector{Symbol}
    over2_branch_reference_names::Vector{Symbol}
    over2_branch_reference_line_numbers::Vector{Int}
    over2_branch_from_node_names::Vector{Symbol}
    over2_branch_to_node_names::Vector{Symbol}
    over2_branch_from_node_indices::Vector{Int}
    over2_branch_to_node_indices::Vector{Int}
    over2_branch_raw_resistance_values::Vector{Float64}
    over2_branch_raw_inductance_values::Vector{Float64}
    over2_branch_raw_capacitance_values::Vector{Float64}
    over2_branch_conductance_values::Vector{Float64}
    over2_branch_resistance_values::Vector{Float64}
    over2_branch_inductance_values::Vector{Float64}
    over2_branch_capacitance_values::Vector{Float64}
    over2_branch_output_codes::Vector{Int}
    over2_branch_count::Int
    over15_output_request_names::Vector{Symbol}
    over15_output_request_output_kinds::Vector{Symbol}
    over15_output_request_request_kinds::Vector{Symbol}
    over15_output_request_layout_kinds::Vector{Symbol}
    over15_output_request_line_numbers::Vector{Int}
    over15_output_request_output_codes::Vector{Int}
    over15_output_request_node_names::Vector{Symbol}
    over15_output_request_node_indices::Vector{Int}
    over15_output_request_branch_names::Vector{Symbol}
    over15_output_request_branch_indices::Vector{Int}
    over15_output_request_count::Int
    over16_output_channel_names::Vector{Symbol}
    over16_output_node_names::Vector{Symbol}
    over16_output_node_indices::Vector{Int}
    over16_output_channel_line_numbers::Vector{Int}
    over16_output_channel_count::Int
    over16_branch_voltage_output_names::Vector{Symbol}
    over16_branch_voltage_branch_names::Vector{Symbol}
    over16_branch_voltage_branch_indices::Vector{Int}
    over16_branch_voltage_output_line_numbers::Vector{Int}
    over16_branch_voltage_output_count::Int
    over16_branch_current_output_names::Vector{Symbol}
    over16_branch_current_branch_names::Vector{Symbol}
    over16_branch_current_branch_indices::Vector{Int}
    over16_branch_current_output_line_numbers::Vector{Int}
    over16_branch_current_output_count::Int
    over16_branch_power_output_names::Vector{Symbol}
    over16_branch_power_branch_names::Vector{Symbol}
    over16_branch_power_branch_indices::Vector{Int}
    over16_branch_power_output_line_numbers::Vector{Int}
    over16_branch_power_output_count::Int
    over5a_source_names::Vector{Symbol}
    over5a_source_node_names::Vector{Symbol}
    over5a_source_node_values::Vector{Int}
    over5a_source_iform_values::Vector{Int}
    over5a_source_line_numbers::Vector{Int}
    over5a_source_layout_kinds::Vector{Symbol}
    over5a_source_tstart_values::Vector{Float64}
    over5a_source_tstop_values::Vector{Float64}
    over5a_source_crest_values::Vector{Float64}
    over5a_source_time1_values::Vector{Float64}
    over5a_source_time2_values::Vector{Float64}
    over5a_source_sfreq_values::Vector{Float64}
    over5a_source_count::Int
    over16_source_card_kinds::Vector{Symbol}
    over16_source_card_values::Vector{Vector{Float64}}
    over16_source_card_provided_value_counts::Vector{Int}
    over16_source_card_line_numbers::Vector{Int}
    over16_source_card_count::Int
    over16_source_interpolation_values::Vector{Vector{Float64}}
    over16_source_interpolation_provided_value_counts::Vector{Int}
    over16_source_interpolation_line_numbers::Vector{Int}
    over16_source_interpolation_count::Int
    over16_source_tacs_override_positions::Vector{Int}
    over16_source_tacs_override_xtcs_indices::Vector{Int}
    over16_source_tacs_override_line_numbers::Vector{Int}
    over16_source_tacs_override_count::Int
    over16_source_analytic_values::Vector{Vector{Float64}}
    over16_source_analytic_provided_value_counts::Vector{Int}
    over16_source_analytic_line_numbers::Vector{Int}
    over16_source_analytic_count::Int
    over5_switch_names::Vector{Symbol}
    over5_switch_line_numbers::Vector{Int}
    over5_switch_from_node_names::Vector{Symbol}
    over5_switch_to_node_names::Vector{Symbol}
    over5_switch_from_node_indices::Vector{Int}
    over5_switch_to_node_indices::Vector{Int}
    over5_switch_layout_kinds::Vector{Symbol}
    over5_switch_raw_close_time_s_values::Vector{Float64}
    over5_switch_raw_open_time_s_values::Vector{Float64}
    over5_switch_close_time_s_values::Vector{Float64}
    over5_switch_open_time_s_values::Vector{Float64}
    over5_switch_initially_closed_flags::Vector{Bool}
    over5_switch_measuring_flags::Vector{Bool}
    over5_switch_closed_markers::Vector{String}
    over5_switch_marker_texts::Vector{String}
    over5_switch_on_conductance_values::Vector{Float64}
    over5_switch_off_conductance_values::Vector{Float64}
    over5_switch_output_codes::Vector{Int}
    over5_switch_count::Int
    time_switch_names::Vector{Symbol}
    time_switch_line_numbers::Vector{Int}
    time_switch_from_node_names::Vector{Symbol}
    time_switch_to_node_names::Vector{Symbol}
    time_switch_close_time_s_values::Vector{Float64}
    time_switch_open_time_s_values::Vector{Float64}
    time_switch_initially_closed_flags::Vector{Bool}
    time_switch_on_conductance_values::Vector{Float64}
    time_switch_off_conductance_values::Vector{Float64}
    time_switch_count::Int
    control_card_kinds::Vector{Symbol}
    control_card_labels::Vector{String}
    control_card_line_numbers::Vector{Int}
    control_card_tokens::Vector{Vector{String}}
    control_card_count::Int
    card_counts::Dict{Symbol,Int}
end

node_count(result::DeckParseResult)::Int = length(result.node_map)

assert_deck_valid!(result::DeckParseResult) = assert_valid!(result.validation)

function deck_nodal_system(result::DeckParseResult)
    assert_deck_valid!(result)
    return NodalSystem(node_count(result), Tuple(result.elements))
end

function deck_node_names(result::DeckParseResult)
    names = Vector{Symbol}(undef, length(result.node_map))
    for (name, index) in result.node_map
        index > 0 || throw(ArgumentError("node indices must be positive"))
        index <= length(names) || throw(ArgumentError("node index $index exceeds node count"))
        names[index] = name
    end
    return names
end

function _deck_scalar_branch_indices(result::DeckParseResult)
    assert_deck_valid!(result)
    return Int[
        index for (index, element) in enumerate(result.elements)
        if element isa ConductanceBranch ||
           element isa SeriesRLBranch ||
           element isa SeriesRLCBranch ||
           element isa CapacitorBranch
    ]
end

deck_branch_count(result::DeckParseResult) = length(_deck_scalar_branch_indices(result))

deck_branch_names(result::DeckParseResult) =
    Symbol[result.element_names[index] for index in _deck_scalar_branch_indices(result)]

deck_branch_kinds(result::DeckParseResult) =
    Symbol[deck_element_kind(result.elements[index]) for index in _deck_scalar_branch_indices(result)]

deck_branch_line_numbers(result::DeckParseResult) =
    Int[result.element_line_numbers[index] for index in _deck_scalar_branch_indices(result)]

_deck_branch_from_node(element::Union{ConductanceBranch,SeriesRLBranch,SeriesRLCBranch,CapacitorBranch}) = element.a
_deck_branch_to_node(element::Union{ConductanceBranch,SeriesRLBranch,SeriesRLCBranch,CapacitorBranch}) = element.b

deck_branch_from_node_indices(result::DeckParseResult) =
    Int[_deck_branch_from_node(result.elements[index]) for index in _deck_scalar_branch_indices(result)]

deck_branch_to_node_indices(result::DeckParseResult) =
    Int[_deck_branch_to_node(result.elements[index]) for index in _deck_scalar_branch_indices(result)]

_deck_branch_conductance_value(element::ConductanceBranch) = element.g
_deck_branch_conductance_value(::Union{SeriesRLBranch,SeriesRLCBranch,CapacitorBranch}) = 0.0

_deck_branch_resistance_value(element::ConductanceBranch) =
    iszero(element.g) ? Inf : inv(element.g)
_deck_branch_resistance_value(element::SeriesRLBranch) = element.r
_deck_branch_resistance_value(element::SeriesRLCBranch) = element.r
_deck_branch_resistance_value(::CapacitorBranch) = 0.0

_deck_branch_inductance_value(::ConductanceBranch) = 0.0
_deck_branch_inductance_value(element::SeriesRLBranch) = element.l
_deck_branch_inductance_value(element::SeriesRLCBranch) = element.l
_deck_branch_inductance_value(::CapacitorBranch) = 0.0

_deck_branch_capacitance_value(::Union{ConductanceBranch,SeriesRLBranch}) = 0.0
_deck_branch_capacitance_value(element::SeriesRLCBranch) = element.c
_deck_branch_capacitance_value(element::CapacitorBranch) = element.c

_deck_branch_previous_current_value(::ConductanceBranch) = 0.0
_deck_branch_previous_current_value(element::Union{SeriesRLBranch,SeriesRLCBranch,CapacitorBranch}) = element.i_prev

_deck_branch_previous_voltage_value(::ConductanceBranch) = 0.0
_deck_branch_previous_voltage_value(element::Union{SeriesRLBranch,SeriesRLCBranch,CapacitorBranch}) = element.v_prev

deck_branch_conductance_values(result::DeckParseResult) =
    Float64[_deck_branch_conductance_value(result.elements[index]) for index in _deck_scalar_branch_indices(result)]

deck_branch_resistance_values(result::DeckParseResult) =
    Float64[_deck_branch_resistance_value(result.elements[index]) for index in _deck_scalar_branch_indices(result)]

deck_branch_inductance_values(result::DeckParseResult) =
    Float64[_deck_branch_inductance_value(result.elements[index]) for index in _deck_scalar_branch_indices(result)]

deck_branch_capacitance_values(result::DeckParseResult) =
    Float64[_deck_branch_capacitance_value(result.elements[index]) for index in _deck_scalar_branch_indices(result)]

deck_branch_previous_current_values(result::DeckParseResult) =
    Float64[_deck_branch_previous_current_value(result.elements[index]) for index in _deck_scalar_branch_indices(result)]

deck_branch_previous_voltage_values(result::DeckParseResult) =
    Float64[_deck_branch_previous_voltage_value(result.elements[index]) for index in _deck_scalar_branch_indices(result)]

function deck_bergeron_line_rows(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.bergeron_line_rows)
end

deck_bergeron_line_names(result::DeckParseResult) =
    Symbol[row.name for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_line_numbers(result::DeckParseResult) =
    Int[row.line_no for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_from_node_names(result::DeckParseResult) =
    Symbol[row.from_node for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_to_node_names(result::DeckParseResult) =
    Symbol[row.to_node for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_from_node_indices(result::DeckParseResult) =
    Int[row.from_node_value for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_to_node_indices(result::DeckParseResult) =
    Int[row.to_node_value for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_surge_impedance_values(result::DeckParseResult) =
    Float64[row.surge_impedance for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_surge_admittance_values(result::DeckParseResult) =
    Float64[row.surge_admittance for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_travel_time_s_values(result::DeckParseResult) =
    Float64[row.travel_time_s for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_dt_s_values(result::DeckParseResult) =
    Float64[row.dt_s for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_attenuation_values(result::DeckParseResult) =
    Float64[row.attenuation for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_delay_step_counts(result::DeckParseResult) =
    Int[row.delay_steps for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_write_indices(result::DeckParseResult) =
    Int[row.write_index for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_history_current_from_values(result::DeckParseResult) =
    Float64[row.history_current_from for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_history_current_to_values(result::DeckParseResult) =
    Float64[row.history_current_to for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_terminal_voltage_from_values(result::DeckParseResult) =
    Float64[row.terminal_voltage_from for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_terminal_voltage_to_values(result::DeckParseResult) =
    Float64[row.terminal_voltage_to for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_terminal_current_from_values(result::DeckParseResult) =
    Float64[row.terminal_current_from for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_terminal_current_to_values(result::DeckParseResult) =
    Float64[row.terminal_current_to for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_traveling_wave_from_values(result::DeckParseResult) =
    Float64[row.traveling_wave_from for row in deck_bergeron_line_rows(result)]

deck_bergeron_line_traveling_wave_to_values(result::DeckParseResult) =
    Float64[row.traveling_wave_to for row in deck_bergeron_line_rows(result)]

deck_coupled_line_rows(result::DeckParseResult) = copy(result.coupled_line_rows)

deck_sampled_frequency_line_rows(result::DeckParseResult) =
    copy(result.sampled_frequency_line_rows)

deck_semlyen_line_groups(result::DeckParseResult) = copy(result.semlyen_line_groups)

deck_rational_frequency_line_groups(result::DeckParseResult) =
    copy(result.rational_frequency_line_groups)

deck_line_modal_transform_rows(result::DeckParseResult) =
    DeckLineModalTransformRow[
        DeckLineModalTransformRow(row.line_no, copy(row.values), row.raw_text)
        for row in result.line_modal_transform_rows
    ]

deck_coupled_line_names(result::DeckParseResult) =
    Symbol[row.name for row in deck_coupled_line_rows(result)]

deck_coupled_line_line_numbers(result::DeckParseResult) =
    Int[row.line_no for row in deck_coupled_line_rows(result)]

deck_coupled_line_types(result::DeckParseResult) =
    Int[row.line_type for row in deck_coupled_line_rows(result)]

deck_coupled_line_kinds(result::DeckParseResult) =
    Symbol[row.line_kind for row in deck_coupled_line_rows(result)]

deck_coupled_line_phase_indices(result::DeckParseResult) =
    Int[row.phase_index for row in deck_coupled_line_rows(result)]

deck_coupled_line_sequence_resistance_values(result::DeckParseResult) =
    Union{Missing,Float64}[row.sequence_resistance for row in deck_coupled_line_rows(result)]

deck_coupled_line_sequence_inductance_values(result::DeckParseResult) =
    Union{Missing,Float64}[row.sequence_inductance for row in deck_coupled_line_rows(result)]

deck_coupled_line_from_node_names(result::DeckParseResult) =
    Symbol[row.from_node for row in deck_coupled_line_rows(result)]

deck_coupled_line_to_node_names(result::DeckParseResult) =
    Symbol[row.to_node for row in deck_coupled_line_rows(result)]

deck_coupled_line_from_node_indices(result::DeckParseResult) =
    Int[row.from_node_value for row in deck_coupled_line_rows(result)]

deck_coupled_line_to_node_indices(result::DeckParseResult) =
    Int[row.to_node_value for row in deck_coupled_line_rows(result)]

deck_coupled_phase_pi_section_rows(result::DeckParseResult) =
    copy(result.coupled_phase_pi_section_rows)

deck_generator_equivalent_rows(result::DeckParseResult) =
    copy(result.generator_equivalent_rows)

function deck_line_constants_conductor_cards(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.line_constants_conductor_cards)
end

function deck_line_constants_frequency_cards(result::DeckParseResult)
    assert_deck_valid!(result)
    return copy(result.line_constants_frequency_cards)
end

function deck_cable_constants_cases(result::DeckParseResult)
    assert_deck_valid!(result)
    return DeckCableConstantsCase[
        DeckCableConstantsCase(
            row.line_no,
            row.cable_kind_code,
            row.surface_position_code,
            row.phase_count,
            row.earth_model_code,
            row.modal_output_flag,
            row.impedance_output_flag,
            row.admittance_output_flag,
            row.pipe_count,
            row.grounding_selector,
            copy(row.pipe_radii_m),
            row.pipe_resistivity_ohm_m,
            row.pipe_relative_permeability,
            row.pipe_inner_insulator_relative_permittivity,
            row.pipe_outer_insulator_relative_permittivity,
            copy(row.cable_to_pipe_center_distances_m),
            copy(row.cable_to_pipe_angles_rad),
            copy(row.pipe_depths_m),
            copy(row.pipe_horizontal_positions_m),
            copy(row.layer_counts),
            copy(row.boundary_radii_m),
            copy(row.resistivity_ohm_m),
            copy(row.conductor_relative_permeability),
            copy(row.insulation_relative_permeability),
            copy(row.insulation_relative_permittivity),
            copy(row.depths_m),
            copy(row.horizontal_positions_m),
            DeckCableConstantsFrequencyCard[
                DeckCableConstantsFrequencyCard(
                    card.line_no,
                    card.earth_resistivity_ohm_m,
                    card.start_frequency_hz,
                    card.decade_count,
                    card.points_per_decade,
                    card.distance_m,
                    card.card_output_flag,
                    card.transform_flag,
                ) for card in row.frequency_cards
            ],
        ) for row in result.cable_constants_cases
    ]
end

function _line_constants_conductor_bundle_count(card::DeckLineConstantsConductorCard)
    return max(card.bundle_conductor_count, 1)
end

function _line_constants_physical_conductor_rows(
    cards::AbstractVector{DeckLineConstantsConductorCard},
)
    expanded = DeckLineConstantsPhysicalConductor[]
    for (source_card_index, card) in enumerate(cards)
        bundle_count = _line_constants_conductor_bundle_count(card)
        bundle_radius_ft = 0.0
        bundle_step_rad = 0.0
        if bundle_count > 1
            angle = (pi - 2.0 * pi / bundle_count) / 2.0
            bundle_radius_ft = card.bundle_spacing_inches / (24.0 * cos(angle))
            bundle_step_rad = 2.0 * pi / bundle_count
        end
        reference_angle_rad = deg2rad(card.bundle_angle_deg)
        for bundle_ordinal in 1:bundle_count
            horizontal = card.horizontal_ft
            height = card.average_height_ft
            if bundle_count > 1
                offset_angle = reference_angle_rad - bundle_step_rad * bundle_ordinal
                horizontal += bundle_radius_ft * cos(offset_angle)
                height += bundle_radius_ft * sin(offset_angle)
            end
            push!(
                expanded,
                DeckLineConstantsPhysicalConductor(
                    source_card_index,
                    card.line_no,
                    bundle_ordinal,
                    card.phase_number,
                    card.skin_effect_type,
                    card.resistance_ohm_per_mile,
                    card.reactance_type,
                    card.reactance_or_gmr,
                    card.diameter_inches,
                    horizontal,
                    height,
                    card.conductor_name,
                ),
            )
        end
    end

    phase_rows = [row for row in expanded if row.phase_number > 0]
    grounded_rows = [row for row in expanded if row.phase_number <= 0]
    ordered = DeckLineConstantsPhysicalConductor[]
    max_bundle = isempty(phase_rows) ? 0 : maximum(row.bundle_ordinal for row in phase_rows)
    for bundle_ordinal in 1:max_bundle
        for row in phase_rows
            row.bundle_ordinal == bundle_ordinal && push!(ordered, row)
        end
    end
    append!(ordered, grounded_rows)
    return ordered
end

const TACS_DIMENSION_VALUE_COUNT = 8
const DEFAULT_TACS_TOTAL_STORAGE_CELLS = 90000
const DEFAULT_TACS_NUMERIC_BYTE_COUNTS = (2, 2, 2, 1)

function deck_line_constants_physical_conductors(result::DeckParseResult)
    assert_deck_valid!(result)
    return _line_constants_physical_conductor_rows(result.line_constants_conductor_cards)
end

function deck_power_frequency_request_rows(result::DeckParseResult)
    return DeckPowerFrequencyRequestRow[
        DeckPowerFrequencyRequestRow(row.line_no, row.frequency_hz, row.raw_text)
        for row in result.power_frequency_request_rows
    ]
end

function deck_universal_machine_dimension_request_rows(result::DeckParseResult)
    return DeckUniversalMachineDimensionRequestRow[
        DeckUniversalMachineDimensionRequestRow(
            row.line_no,
            row.coil_table_size,
            row.machine_table_size,
            row.output_table_size,
            row.bus_table_size,
            row.raw_text,
        )
        for row in result.universal_machine_dimension_request_rows
    ]
end

function deck_universal_machine_section_rows(result::DeckParseResult)
    return DeckUniversalMachineSectionRow[
        DeckUniversalMachineSectionRow(
            row.line_no,
            row.machine_count,
            row.input_layout,
            row.parameter_basis,
            row.remanent_flux_enabled,
            row.initialization_mode,
            row.detailed_machine_input,
            row.maximum_shaft_mass_count,
            row.terminal_coupling,
            row.raw_text,
        )
        for row in result.universal_machine_section_rows
    ]
end

function deck_universal_machine_definition_rows(result::DeckParseResult)
    return DeckUniversalMachineDefinitionRow[
        DeckUniversalMachineDefinitionRow(
            row.line_no,
            row.machine_index,
            row.card_index,
            row.machine_type,
            row.d_axis_coil_count,
            row.q_axis_coil_count,
            row.torque_output_flag,
            row.speed_output_flag,
            row.angle_output_flag,
            row.pole_pair_count,
            row.rotor_mass,
            row.saturation_mode,
            row.saturated_inductance,
            row.saturation_flux,
            row.remanent_flux,
            row.value1,
            row.value2,
            row.mechanical_damping_coefficient,
            row.speed_convergence_margin,
            row.node,
            row.raw_text,
        )
        for row in result.universal_machine_definition_rows
    ]
end

function deck_universal_machine_coil_rows(result::DeckParseResult)
    return DeckUniversalMachineCoilRow[
        DeckUniversalMachineCoilRow(
            row.line_no,
            row.machine_index,
            row.coil_index,
            row.resistance,
            row.inductance,
            row.terminal_node,
            row.control_signal,
            row.output_flag,
            row.initial_history_current,
            row.raw_text,
        )
        for row in result.universal_machine_coil_rows
    ]
end

function deck_universal_machine_terminal_rows(result::DeckParseResult)
    return DeckUniversalMachineTerminalRow[
        DeckUniversalMachineTerminalRow(
            row.line_no,
            row.machine_index,
            row.terminal_index,
            row.terminal_node,
            row.reference_node,
            row.terminal_node_value,
            row.reference_node_value,
            row.raw_text,
        )
        for row in result.universal_machine_terminal_rows
    ]
end

function deck_universal_machine_generated_branch_rows(result::DeckParseResult)
    return DeckUniversalMachineGeneratedBranchRow[
        DeckUniversalMachineGeneratedBranchRow(
            row.line_no,
            row.machine_index,
            row.branch_index,
            row.from_node,
            row.to_node,
            row.from_node_value,
            row.to_node_value,
            row.reactance,
            row.raw_text,
        )
        for row in result.universal_machine_generated_branch_rows
    ]
end

function deck_universal_machine_speed_capacitor_rows(result::DeckParseResult)
    return DeckUniversalMachineSpeedCapacitorRow[
        DeckUniversalMachineSpeedCapacitorRow(
            row.line_no,
            row.machine_index,
            row.capacitor_node,
            row.mass_node,
            row.capacitor_node_value,
            row.mass_node_value,
            row.resistance,
            row.capacitance,
            row.raw_text,
        )
        for row in result.universal_machine_speed_capacitor_rows
    ]
end

function deck_universal_machine_node_summary_rows(result::DeckParseResult)
    return DeckUniversalMachineNodeSummaryRow[
        DeckUniversalMachineNodeSummaryRow(
            row.line_no,
            row.machine_index,
            row.mass_node,
            row.mass_node_value,
            row.mechanical_slack_source,
            row.mechanical_slack_source_index,
            row.field_slack_source,
            row.field_slack_source_index,
            row.raw_text,
        )
        for row in result.universal_machine_node_summary_rows
    ]
end

function deck_universal_machine_output_summary_rows(result::DeckParseResult)
    return DeckUniversalMachineOutputSummaryRow[
        DeckUniversalMachineOutputSummaryRow(
            row.line_no,
            row.machine_count,
            row.output_count,
            row.tacs_transfer_count,
            row.raw_text,
        )
        for row in result.universal_machine_output_summary_rows
    ]
end

function deck_synchronous_machine_terminal_voltage_rows(result::DeckParseResult)
    return DeckSynchronousMachineTerminalVoltageRow[
        DeckSynchronousMachineTerminalVoltageRow(
            row.line_no,
            row.machine_index,
            row.source_group_index,
            row.source_type,
            row.phase_index,
            row.terminal_node,
            row.terminal_node_value,
            row.source_node_value,
            row.peak_terminal_voltage,
            row.frequency_hz,
            row.angle_deg,
            row.raw_text,
        )
        for row in result.synchronous_machine_terminal_voltage_rows
    ]
end

function deck_synchronous_machine_tolerance_rows(result::DeckParseResult)
    return DeckSynchronousMachineToleranceRow[
        DeckSynchronousMachineToleranceRow(
            row.line_no,
            row.machine_index,
            copy(row.values),
            row.raw_text,
        )
        for row in result.synchronous_machine_tolerance_rows
    ]
end

function deck_synchronous_machine_parameter_fitting_rows(result::DeckParseResult)
    return DeckSynchronousMachineParameterFittingRow[
        DeckSynchronousMachineParameterFittingRow(
            row.line_no,
            row.machine_index,
            row.value,
            row.raw_text,
        )
        for row in result.synchronous_machine_parameter_fitting_rows
    ]
end

function deck_synchronous_machine_model_parameter_rows(result::DeckParseResult)
    return DeckSynchronousMachineModelParameterRow[
        DeckSynchronousMachineModelParameterRow(
            row.line_no,
            row.machine_index,
            row.parameter_kind,
            copy(row.values),
            copy(row.positional_values),
            row.raw_text,
        )
        for row in result.synchronous_machine_model_parameter_rows
    ]
end

function deck_synchronous_machine_mass_rows(result::DeckParseResult)
    return DeckSynchronousMachineMassRow[
        DeckSynchronousMachineMassRow(
            row.line_no,
            row.machine_index,
            row.mass_index,
            copy(row.values),
            row.torque_fraction,
            row.inertia,
            row.speed_deviation_damping,
            row.mutual_damping,
            row.shaft_stiffness,
            row.absolute_speed_damping,
            row.raw_text,
        )
        for row in result.synchronous_machine_mass_rows
    ]
end

function deck_synchronous_machine_output_request_rows(result::DeckParseResult)
    return DeckSynchronousMachineOutputRequestRow[
        DeckSynchronousMachineOutputRequestRow(
            row.line_no,
            row.machine_index,
            row.group_index,
            copy(row.output_codes),
            row.dynamic_output_count,
            row.raw_text,
        )
        for row in result.synchronous_machine_output_request_rows
    ]
end

function deck_synchronous_machine_control_interface_rows(result::DeckParseResult)
    return DeckSynchronousMachineControlInterfaceRow[
        DeckSynchronousMachineControlInterfaceRow(
            row.line_no,
            row.machine_index,
            row.interface_code,
            row.direction,
            row.coupling_kind,
            row.signal_name,
            row.variable_index,
            row.raw_text,
        )
        for row in result.synchronous_machine_control_interface_rows
    ]
end

function deck_synchronous_machine_output_summary_rows(result::DeckParseResult)
    return DeckSynchronousMachineOutputSummaryRow[
        DeckSynchronousMachineOutputSummaryRow(
            row.line_no,
            row.machine_count,
            row.terminal_count,
            row.mass_count,
            row.output_count,
            row.raw_text,
        )
        for row in result.synchronous_machine_output_summary_rows
    ]
end

function deck_tacs_dimension_request_rows(result::DeckParseResult)
    return DeckTACSDimensionRequestRow[
        DeckTACSDimensionRequestRow(
            row.line_no,
            row.payload_line_no,
            row.allocation_kind,
            row.request_values,
            row.list_sizes,
            row.provided_value_count,
            row.request_text,
            row.payload_text,
        )
        for row in result.tacs_dimension_request_rows
    ]
end

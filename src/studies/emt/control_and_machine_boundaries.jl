
function _deck_with_prescribed_machine_source(
    parsed::DeckParser.DeckParseResult,
    source_node::Int,
    source_row,
    crest::Float64,
    angular_frequency_rad_s::Float64,
)
    runtime_parsed = deepcopy(parsed)
    source_indices = findall(runtime_parsed.elements) do element
        element isa TheveninSource && element.node == source_node
    end
    length(source_indices) == 1 ||
        throw(ArgumentError("machine initialization requires one prescribed source element"))
    source_index = only(source_indices)
    source_element = runtime_parsed.elements[source_index]
    runtime_parsed.elements[source_index] = TheveninSource(
        source_element.node,
        source_element.g,
        time_s -> _machine_drive_source_value(
            source_row,
            crest,
            time_s;
            angular_frequency_rad_s = angular_frequency_rad_s,
        ),
    )
    return runtime_parsed
end

struct ManualDirectMachineInitialVoltageBoundary
    machine_index::Int
    machine_type::Int
    power_terminal_nodes::NTuple{3,Int}
    power_terminal_voltages::NTuple{3,Float64}
end

function _manual_direct_machine_initial_voltage_boundary(
    parsed::DeckParser.DeckParseResult,
    machine_index::Int,
)
    card1 = _deck_universal_machine_definition(parsed, machine_index, 1)
    card1.machine_type in (8, 9, 10, 11, 12) || throw(ArgumentError(
        "manual direct-machine voltage boundary requires machine type 8 through 12",
    ))
    _deck_universal_machine_initialization_mode(parsed) == :manual ||
        throw(ArgumentError("manual direct-machine voltage boundary requires manual initialization"))
    network_nodes = _deck_universal_machine_network_nodes(parsed, machine_index)
    return ManualDirectMachineInitialVoltageBoundary(
        machine_index,
        card1.machine_type,
        Tuple(Int.(network_nodes.power)),
        (0.0, 0.0, 0.0),
    )
end

function _connected_nodal_component_nodes(
    system::NodalSystem,
    seed_nodes::AbstractVector{<:Integer},
)
    included = falses(system.node_count)
    for node in seed_nodes
        1 <= node <= system.node_count || throw(ArgumentError(
            "connected-component seed node $node is outside the nodal system",
        ))
        included[node] = true
    end

    changed = true
    while changed
        changed = false
        for element in system.elements
            fields = fieldnames(typeof(element))
            (:a in fields && :b in fields) || continue
            terminal_nodes = Int[]
            for field in (:a, :b)
                value = getfield(element, field)
                if value isa Integer
                    value != 0 && push!(terminal_nodes, Int(value))
                elseif value isa AbstractVector{<:Integer}
                    append!(terminal_nodes, Int.(filter(!=(0), value)))
                end
            end
            any(node -> included[node], terminal_nodes) || continue
            for node in terminal_nodes
                included[node] && continue
                included[node] = true
                changed = true
            end
        end
    end
    return findall(included)
end

function run_deck_universal_machine_horizon(
    parsed::DeckParser.DeckParseResult;
    machine_index::Int=1,
    time_step_s::Real=deck_fixed_step_horizon(parsed).dt_s,
)
    DeckParser.assert_deck_valid!(parsed)
    horizon = deck_fixed_step_horizon(parsed)
    horizon.step_count > 0 ||
        throw(ArgumentError("machine network horizon requires at least one dynamic step"))
    dt_s = Float64(time_step_s)
    dt_s == horizon.dt_s ||
        throw(ArgumentError("machine network horizon time step must match the deck time step"))
    network_nodes = _deck_universal_machine_network_nodes(parsed, machine_index)
    machine_section = _deck_universal_machine_section(parsed)
    predicted_current_coupling =
        machine_section.terminal_coupling == :predicted_current
    card1 = _deck_universal_machine_definition(parsed, machine_index, 1)
    automatic_direct_machine_initialization =
        card1.machine_type in (8, 9, 10, 11, 12) &&
        _deck_universal_machine_initialization_mode(parsed) == :automatic ?
        _deck_direct_machine_automatic_initialization(
            parsed,
            machine_index,
        ) : nothing
    single_phase_initialization =
        card1.machine_type in (6, 7) &&
        _deck_universal_machine_initialization_mode(parsed) == :automatic ?
        _deck_single_phase_induction_initialization(parsed, machine_index) : nothing
    fixed_source_synchronous_steady_state =
        card1.machine_type in (1, 2) &&
        !isempty(DeckParser.deck_fixed_source_constraint_rows(parsed)) ?
        _deck_fixed_source_synchronous_steady_state(parsed) : nothing
    steady_state = automatic_direct_machine_initialization !== nothing ?
        automatic_direct_machine_initialization.steady_state :
        single_phase_initialization !== nothing ?
        single_phase_initialization.steady_state :
        fixed_source_synchronous_steady_state === nothing ?
        deck_steady_state_voltage_phasors(parsed) :
        fixed_source_synchronous_steady_state
    state = automatic_direct_machine_initialization !== nothing ?
        automatic_direct_machine_initialization.state :
        single_phase_initialization === nothing ?
        deck_coupled_dq_machine_initial_state(
            parsed;
            machine_index = machine_index,
            steady_state = steady_state,
        ) : single_phase_initialization.state
    parameters = deck_coupled_dq_machine_parameters(
        parsed;
        machine_index = machine_index,
        time_step_s = dt_s,
    )
    normalized_machine_coordinates =
        machine_section.parameter_basis == :power_frequency_normalized
    power_axis_count = parameters.machine_type in (1, 3, 4) ? 3.0 :
        parameters.machine_type in (2, 5) ? 2.0 : 1.0
    speed_rad_s_per_network_unit = normalized_machine_coordinates ?
        parameters.synchronous_electrical_speed_rad_s : 1.0
    torque_network_units_per_Nm = normalized_machine_coordinates ?
        parameters.synchronous_electrical_speed_rad_s / power_axis_count : 1.0
    impedance_rad_s_per_Nm_per_network_unit =
        speed_rad_s_per_network_unit * torque_network_units_per_Nm
    synchronous_field = parameters.machine_type in (1, 2)
    single_phase_machine = parameters.machine_type in (6, 7)
    manually_initialized_direct_machine =
        parameters.machine_type in (8, 9, 10, 11, 12) &&
        _deck_universal_machine_initialization_mode(parsed) == :manual
    manually_initialized_machine =
        _deck_universal_machine_initialization_mode(parsed) == :manual
    external_field = parameters.machine_type in (1, 2, 6, 7, 8)
    excitation_nodes = external_field ? Int[network_nodes.field] : Int[]
    excitation_sources = external_field ? [network_nodes.field_source] : Any[]
    if parameters.machine_type == 7
        push!(excitation_nodes, Int(network_nodes.q_axis_field))
        push!(excitation_sources, network_nodes.q_axis_field_source)
    end
    card4 = manually_initialized_machine ? nothing :
        _deck_universal_machine_definition(parsed, machine_index, 4)
    slip = card4 === nothing || card4.value1 === missing ? 0.0 :
        Float64(card4.value1) / 100.0
    excitation_source_frequencies_hz = if single_phase_machine
        fill(
            abs(slip) * _deck_steady_state_frequency_hz(parsed),
            length(excitation_sources),
        )
    elseif external_field
        Float64[source.sfreq / (2.0 * pi) for source in excitation_sources]
    else
        Float64[]
    end
    excitation_source_angular_frequencies_rad_s =
        2.0 * pi .* excitation_source_frequencies_hz
    desired_excitation_node_voltages = Float64[
        _machine_drive_source_value(
            source,
            source.crest,
            0.0;
            angular_frequency_rad_s = angular_frequency,
        ) for (source, angular_frequency) in zip(
            excitation_sources,
            excitation_source_angular_frequencies_rad_s,
        )
    ]
    synchronous_field && !manually_initialized_machine &&
        (desired_excitation_node_voltages[1] =
            state.current_values[4] / parameters.coil_conductances[4])
    initialized_excitation_crests = copy(desired_excitation_node_voltages)
    if synchronous_field && !manually_initialized_machine
        unit_field_at_initial_step =
            _machine_drive_source_value(
                network_nodes.field_source,
                1.0,
                0.0;
                angular_frequency_rad_s =
                    excitation_source_angular_frequencies_rad_s[1],
            )
        abs(unit_field_at_initial_step) > eps(Float64) ||
            throw(ArgumentError("wound-field excitation source is zero at initialization"))
        initialized_excitation_crests[1] =
            desired_excitation_node_voltages[1] / unit_field_at_initial_step
    elseif external_field
        initialized_excitation_crests .= Float64[source.crest for source in excitation_sources]
    end
    # The horizon is replayable from one parsed deck. Isolate its mutable
    # branch companions without changing the mutation contract of generic
    # step contexts used by report and oracle owners.
    runtime_parsed = deepcopy(parsed)
    power_terminal_norton_conductances = zeros(Float64, 3)
    if predicted_current_coupling
        if parameters.machine_type in (1, 3, 4)
            power_terminal_norton_conductances .=
                parameters.coil_conductances[2]
        elseif parameters.machine_type in (2, 5)
            power_terminal_norton_conductances[2:3] .=
                parameters.coil_conductances[2:3]
        else
            power_terminal_norton_conductances[3] =
                parameters.coil_conductances[3]
        end
        terminal_rows = sort!(
            [
                row for row in DeckParser.deck_universal_machine_terminal_rows(runtime_parsed)
                if row.machine_index == machine_index
            ];
            by = row -> row.terminal_index,
        )
        for position in network_nodes.active_power_positions
            phase_norton_conductance =
                power_terminal_norton_conductances[position]
            phase_norton_conductance > 0.0 ||
                throw(ArgumentError("predicted-current coupling requires positive power-terminal conductance"))
            row = terminal_rows[position]
            push!(
                runtime_parsed.elements,
                ConductanceBranch(
                    row.terminal_node_value,
                    row.reference_node_value,
                    phase_norton_conductance,
                ),
            )
            push!(runtime_parsed.element_line_numbers, row.line_no)
            push!(
                runtime_parsed.element_names,
                Symbol("machine_terminal_norton_conductance_", machine_index, "_", position),
            )
        end
    end
    for index in eachindex(excitation_sources)
        runtime_parsed = _deck_with_prescribed_machine_source(
            runtime_parsed,
            excitation_nodes[index],
            excitation_sources[index],
            initialized_excitation_crests[index],
            excitation_source_angular_frequencies_rad_s[index],
        )
    end
    context = initialize_step_context(
        runtime_parsed;
        dt_s = dt_s,
        t_end_s = horizon.t_end_s,
        recorded_step_indices = [0, horizon.step_count],
    )
    if !manually_initialized_machine
        _seed_steady_state_network_state!(context, steady_state)
        if automatic_direct_machine_initialization !== nothing &&
           parameters.machine_type in (9, 10, 11, 12)
            _seed_direct_machine_power_leakage_currents!(
                context,
                runtime_parsed,
                [machine_index],
                [
                    automatic_direct_machine_initialization.armature_current_injection,
                ],
            )
        end
        _seed_machine_mechanical_inertia!(
            context.system,
            network_nodes.mechanical,
            state.mechanical_speed_rad_s / speed_rad_s_per_network_unit,
        )
    end

    call_count = horizon.step_count + 1
    coil_count = length(parameters.coil_conductances)
    outputs = zeros(Float64, coil_count + 3, call_count)
    currents = zeros(Float64, coil_count, call_count)
    histories = zeros(Float64, coil_count, call_count)
    substitution_count = if parameters.machine_type in (1, 2)
        5
    elseif parameters.machine_type == 6
        3
    elseif parameters.machine_type == 7
        4
    elseif parameters.machine_type == 8
        3
    elseif parameters.machine_type == 9
        2
    elseif parameters.machine_type == 10
        2
    elseif parameters.machine_type == 11
        2
    elseif parameters.machine_type == 12
        2
    elseif parameters.machine_type in (3, 5)
        4
    else
        7
    end
    substitutions = zeros(Float64, substitution_count, call_count)
    predicted_terminal_currents = zeros(Float64, 3, call_count)
    power_terminal_network_currents = zeros(Float64, 3, call_count)
    terminal_prediction = MachineTerminalPredictionState()
    d_flux = zeros(Float64, call_count)
    q_flux = zeros(Float64, call_count)
    torque = zeros(Float64, call_count)
    speed = zeros(Float64, call_count)
    angle = zeros(Float64, call_count)
    iterations = zeros(Int, call_count)
    times = Float64[step * dt_s for step in 0:horizon.step_count]
    terminal_voltages = zeros(Float64, 3, call_count)
    rotor_thevenin = zeros(Float64, 3, 3, call_count)
    speed_thevenin = zeros(Float64, call_count)
    torque_impedance = zeros(Float64, call_count)
    drive_source_values = zeros(Float64, call_count)
    excitation_source_matrix = zeros(
        Float64,
        length(excitation_sources),
        call_count,
    )
    compensated_voltages = zeros(Float64, context.system.node_count, call_count)
    runtime_output_values = zeros(
        Float64,
        length(context.output_channel_names),
        call_count,
    )

    !isempty(excitation_sources) && !manually_initialized_direct_machine &&
        (excitation_source_matrix[:, 1] .= desired_excitation_node_voltages)

    function record_machine_result!(index, result)
        outputs[:, index] .= result.output_values
        currents[:, index] .= result.current_values
        histories[:, index] .= result.history_currents
        substitutions[:, index] .= result.current_substitution_values
        if predicted_current_coupling
            substitutions[1:(end - 1), index] .= 0.0
        end
        d_flux[index] = result.d_axis_flux
        q_flux[index] = result.q_axis_flux
        torque[index] = result.generated_torque
        speed[index] = result.mechanical_speed_rad_s
        angle[index] = result.mechanical_angle_rad
        iterations[index] = result.iteration_count
        return result
    end

    active_power_positions = network_nodes.active_power_positions
    active_power_count = length(active_power_positions)
    stator_count = coil_count - 3
    single_phase_excitation_impedances = single_phase_machine ?
        Float64[
            _deck_single_phase_excitation_resistance_ohm(parsed, node)
            for node in excitation_nodes
        ] : Float64[]
    stored_power_voltages(values) = Float64[
        node > 0 ? values[node] : 0.0 for node in network_nodes.power
    ]
    initial_power_voltages = if manually_initialized_machine
        parameters.machine_type in (8, 9, 10, 11, 12) ?
            collect(
                _manual_direct_machine_initial_voltage_boundary(
                    parsed,
                    machine_index,
                ).power_terminal_voltages,
            ) :
            stored_power_voltages(steady_state.node_voltage_values)
    else
        stored_power_voltages(steady_state.node_voltage_values)
    end
    initial_machine_power_voltages = copy(initial_power_voltages)
    initial_machine_rotor_thevenin = zeros(Float64, 3, 3)
    if automatic_direct_machine_initialization !== nothing &&
       parameters.machine_type in (9, 10, 11, 12)
        active_position = only(network_nodes.active_power_positions)
        initial_machine_power_voltages[active_position] =
            automatic_direct_machine_initialization.runtime_power_terminal_open_circuit_voltage
        initial_machine_rotor_thevenin[active_position, active_position] =
            automatic_direct_machine_initialization.runtime_power_terminal_thevenin_impedance
    end
    initial_stator_voltages = manually_initialized_direct_machine ?
        zeros(Float64, stator_count) : external_field ?
        vcat(
            -desired_excitation_node_voltages,
            zeros(Float64, stator_count - length(excitation_sources)),
        ) :
        zeros(Float64, stator_count)

    initial_prediction_history = predicted_current_coupling ?
        copy(state.history_currents[1:3]) : Float64[]
    retained_automatic_histories =
        automatic_direct_machine_initialization !== nothing &&
        parameters.machine_type in (9, 10, 11, 12) ?
        copy(state.history_currents) : Float64[]
    initial_result = coupled_dq_machine_step!(
        state,
        parameters;
        power_terminal_voltages = initial_machine_power_voltages,
        rotor_thevenin_matrix = initial_machine_rotor_thevenin,
        mechanical_speed_thevenin_rad_s = state.mechanical_speed_rad_s,
        generated_torque_impedance = 0.0,
        stator_terminal_voltages = initial_stator_voltages,
        stator_thevenin_matrix = zeros(stator_count, stator_count),
        initial_step = true,
    )
    isempty(retained_automatic_histories) ||
        (state.history_currents .= retained_automatic_histories)
    record_machine_result!(1, initial_result)
    if predicted_current_coupling
        histories[1:3, 1] .= initial_prediction_history
        predicted_terminal_currents[:, 1] .= predict_machine_terminal_currents!(
            terminal_prediction,
            state,
            parameters;
            time_s = 0.0,
        )
        predicted_terminal_currents[
            setdiff(1:3, network_nodes.active_power_positions),
            1,
        ] .= 0.0
        state.history_currents[1] = 0.0
        state.history_currents[2] =
            terminal_prediction.previous_d_axis_internal_flux_Wb
        state.history_currents[3] =
            terminal_prediction.previous_q_axis_internal_flux_Wb
    end
    terminal_voltages[:, 1] .= initial_power_voltages
    speed_thevenin[1] = manually_initialized_machine ?
        0.0 : state.mechanical_speed_rad_s
    if manually_initialized_machine
        fill!(@view(compensated_voltages[:, 1]), 0.0)
        for (position, node) in pairs(network_nodes.power)
            node > 0 &&
                (compensated_voltages[node, 1] = initial_power_voltages[position])
        end
        for index in eachindex(excitation_nodes)
            compensated_voltages[excitation_nodes[index], 1] =
                desired_excitation_node_voltages[index]
        end
    else
        compensated_voltages[:, 1] .= steady_state.node_voltage_values
        compensated_voltages[network_nodes.mechanical, 1] =
            state.mechanical_speed_rad_s / speed_rad_s_per_network_unit
        for index in eachindex(excitation_nodes)
            compensated_voltages[excitation_nodes[index], 1] =
                desired_excitation_node_voltages[index]
        end
    end
    active_power_nodes = network_nodes.power[active_power_positions]
    power_compensation_nodes =
        predicted_current_coupling ? Int[] : active_power_nodes
    compensation_nodes = vcat(
        power_compensation_nodes,
        network_nodes.mechanical,
        excitation_nodes,
    )
    calibration_nodes = vcat(compensation_nodes, network_nodes.drive)
    mechanical_position = length(power_compensation_nodes) + 1
    drive_position = length(calibration_nodes)
    initialized_drive_crest = if manually_initialized_machine
        network_nodes.drive_source.crest
    else
        calibration = solve_step_with_compensated_current_injections!(
            deepcopy(context.system),
            dt_s,
            dt_s,
            zeros(Float64, context.system.node_count),
            calibration_nodes,
            (_voltage, _impedance) -> zeros(Float64, length(calibration_nodes)),
        )
        machine_torque_current =
            -initial_result.generated_torque * torque_network_units_per_Nm
        drive_transfer =
            calibration.compensation_impedance[mechanical_position, drive_position]
        abs(drive_transfer) > eps(Float64) ||
            throw(ArgumentError("mechanical drive source is disconnected from the machine node"))
        original_drive_at_first_step = _machine_drive_source_value(
            network_nodes.drive_source,
            network_nodes.drive_source.crest,
            dt_s,
        )
        drive_delta = (
            state.mechanical_speed_rad_s / speed_rad_s_per_network_unit -
            calibration.open_circuit_voltage[network_nodes.mechanical] -
            calibration.compensation_impedance[
                mechanical_position,
                mechanical_position,
            ] * machine_torque_current
        ) / drive_transfer
        unit_drive_at_first_step =
            _machine_drive_source_value(network_nodes.drive_source, 1.0, dt_s)
        abs(unit_drive_at_first_step) > eps(Float64) ||
            throw(ArgumentError("mechanical drive source is zero at the first dynamic step"))
        (original_drive_at_first_step + drive_delta) / unit_drive_at_first_step
    end
    initial_drive_source_network_value = _machine_drive_source_value(
        network_nodes.drive_source,
        initialized_drive_crest,
        0.0,
    )
    drive_source_values[1] =
        initial_drive_source_network_value / torque_network_units_per_Nm
    initial_base_injections = zeros(Float64, context.system.node_count)
    initial_base_injections[network_nodes.drive] =
        initial_drive_source_network_value -
        _machine_drive_source_value(
            network_nodes.drive_source,
            network_nodes.drive_source.crest,
            0.0,
        )
    initial_power_currents = predicted_current_coupling ?
        Float64[] :
        initial_result.output_values[3 .+ active_power_positions]
    initial_compensation_currents = external_field ?
        vcat(
            initial_power_currents,
            -initial_result.generated_torque * torque_network_units_per_Nm,
            initial_result.output_values[7:(6 + length(excitation_sources))],
        ) :
        vcat(
            initial_power_currents,
            -initial_result.generated_torque * torque_network_units_per_Nm,
        )
    if !manually_initialized_machine
        initial_network_result = solve_step_with_compensated_current_injections!(
            deepcopy(context.system),
            0.0,
            dt_s,
            initial_base_injections,
            compensation_nodes,
            (_voltage, _impedance) -> initial_compensation_currents,
        )
        initial_dynamic_nodes = _connected_nodal_component_nodes(
            context.system,
            vcat(network_nodes.mechanical, network_nodes.drive, excitation_nodes),
        )
        compensated_voltages[initial_dynamic_nodes, 1] .=
            initial_network_result.voltage[initial_dynamic_nodes]
    end
    context.step_index = 0
    context.t_s = 0.0
    initial_voltage_values = @view compensated_voltages[:, 1]
    _update_deck_power_energy_state!(context, initial_voltage_values)
    _record_context_outputs!(runtime_output_values, 1, context, initial_voltage_values)

    network_solve_count = 0
    network_correction_count = 0
    for index in 2:call_count
        time_s = times[index]
        base_injections = zeros(Float64, context.system.node_count)
        base_injections[network_nodes.drive] =
            _machine_drive_source_value(
                network_nodes.drive_source,
                initialized_drive_crest,
                time_s,
            ) -
            _machine_drive_source_value(
                network_nodes.drive_source,
                network_nodes.drive_source.crest,
                time_s,
            )
        if predicted_current_coupling
            for position in active_power_positions
                base_injections[network_nodes.power[position]] +=
                    predicted_terminal_currents[position, index - 1]
            end
        end
        drive_source_values[index] = _machine_drive_source_value(
            network_nodes.drive_source,
            initialized_drive_crest,
            time_s,
        ) / torque_network_units_per_Nm
        for source_index in eachindex(excitation_sources)
            excitation_source_matrix[source_index, index] =
                _machine_drive_source_value(
                    excitation_sources[source_index],
                    initialized_excitation_crests[source_index],
                    time_s,
                    angular_frequency_rad_s =
                        excitation_source_angular_frequencies_rad_s[source_index],
                )
        end
        result_ref = Ref{Any}()
        power_prediction_history = predicted_current_coupling ?
            copy(state.history_currents[1:3]) : Float64[]
        compensation_current = function (open_circuit_voltage, impedance)
            power_terminal_voltages =
                stored_power_voltages(open_circuit_voltage)
            prescribed_power_currents = if predicted_current_coupling
                values = zeros(Float64, 3)
                for position in active_power_positions
                    values[position] =
                        predicted_terminal_currents[position, index - 1] -
                        power_terminal_norton_conductances[position] *
                        power_terminal_voltages[position]
                end
                values
            else
                Float64[]
            end
            excitation_voltages = single_phase_machine ?
                zeros(Float64, stator_count) : external_field ?
                vcat(
                    -open_circuit_voltage[network_nodes.field],
                    zeros(Float64, stator_count - 1),
                ) : zeros(Float64, stator_count)
            excitation_thevenin = if single_phase_machine
                Matrix(Diagonal(single_phase_excitation_impedances))
            elseif external_field
                zeros(Float64, stator_count, stator_count)
            else
                zeros(Float64, stator_count, stator_count)
            end
            rotor_thevenin_matrix = zeros(Float64, 3, 3)
            if !single_phase_machine && !predicted_current_coupling
                rotor_thevenin_matrix[active_power_positions, active_power_positions] .=
                    impedance[
                        1:length(power_compensation_nodes),
                        1:length(power_compensation_nodes),
                    ]
            end
            result = coupled_dq_machine_step!(
                state,
                parameters;
                power_terminal_voltages = power_terminal_voltages,
                rotor_thevenin_matrix = rotor_thevenin_matrix,
                mechanical_speed_thevenin_rad_s =
                    open_circuit_voltage[network_nodes.mechanical] *
                    speed_rad_s_per_network_unit,
                generated_torque_impedance =
                    impedance[mechanical_position, mechanical_position] *
                    impedance_rad_s_per_Nm_per_network_unit,
                stator_terminal_voltages = excitation_voltages,
                stator_thevenin_matrix = excitation_thevenin,
                prescribed_power_terminal_currents =
                    prescribed_power_currents,
                initial_step = false,
            )
            result_ref[] = result
            power_currents = predicted_current_coupling ?
                Float64[] :
                result.output_values[3 .+ active_power_positions]
            return external_field ?
                vcat(
                    power_currents,
                    -result.generated_torque * torque_network_units_per_Nm,
                    result.output_values[7:(6 + length(excitation_sources))],
                ) :
                vcat(
                    power_currents,
                    -result.generated_torque * torque_network_units_per_Nm,
                )
        end
        network_result = solve_step_with_compensated_current_injections!(
            context.system,
            time_s,
            dt_s,
            base_injections,
            compensation_nodes,
            compensation_current,
        )
        result = result_ref[]
        record_machine_result!(index, result)
        if predicted_current_coupling
            histories[1:3, index] .= power_prediction_history
            predicted_terminal_currents[:, index] .=
                predict_machine_terminal_currents!(
                    terminal_prediction,
                    state,
                    parameters;
                    time_s = time_s,
                )
            predicted_terminal_currents[
                setdiff(1:3, active_power_positions),
                index,
            ] .= 0.0
            state.history_currents[1] = 0.0
            state.history_currents[2] =
                terminal_prediction.previous_d_axis_internal_flux_Wb
            state.history_currents[3] =
                terminal_prediction.previous_q_axis_internal_flux_Wb
        end
        terminal_voltages[:, index] .= stored_power_voltages(
            network_result.open_circuit_voltage,
        )
        if predicted_current_coupling
            for position in active_power_positions
                power_terminal_network_currents[position, index] =
                    predicted_terminal_currents[position, index - 1] -
                    power_terminal_norton_conductances[position] *
                    terminal_voltages[position, index]
            end
        end
        if !single_phase_machine && !predicted_current_coupling
            rotor_thevenin[active_power_positions, active_power_positions, index] .=
                network_result.compensation_impedance[
                    1:length(power_compensation_nodes),
                    1:length(power_compensation_nodes),
                ]
        end
        speed_thevenin[index] =
            network_result.open_circuit_voltage[network_nodes.mechanical] *
            speed_rad_s_per_network_unit
        torque_impedance[index] =
            network_result.compensation_impedance[
                mechanical_position,
                mechanical_position,
            ] * impedance_rad_s_per_Nm_per_network_unit
        compensated_voltages[:, index] .= network_result.voltage
        if single_phase_machine
            for source_index in eachindex(excitation_nodes)
                compensated_voltages[excitation_nodes[source_index], index] =
                    result.output_values[6 + source_index] *
                    single_phase_excitation_impedances[source_index]
            end
        end
        context.step_index = index - 1
        context.t_s = time_s
        voltage_values = @view compensated_voltages[:, index]
        _update_deck_power_energy_state!(context, voltage_values)
        _record_context_outputs!(runtime_output_values, index, context, voltage_values)
        network_solve_count += 1
        network_correction_count += 1
    end

    report_output_names = unique(vcat(
        _deck_requested_electrical_output_names(parsed),
        _deck_control_system_trace_output_names(parsed),
    ))
    report_output_indices = Int[]
    for name in report_output_names
        index = findfirst(==(name), context.output_channel_names)
        index === nothing && throw(ArgumentError(
            "deck report output channel $name is missing from the native machine runtime",
        ))
        push!(report_output_indices, index)
    end
    report_output_values = runtime_output_values[report_output_indices, :]
    control_system_step_count =
        context.control_system_runtime === nothing ? 0 :
        context.control_system_runtime.executed_step_count

    return DeckUniversalMachineHorizon(
        parsed.source,
        parameters.machine_type,
        machine_section.input_layout,
        machine_section.parameter_basis,
        machine_section.remanent_flux_enabled,
        machine_section.initialization_mode,
        machine_section.maximum_shaft_mass_count,
        machine_section.terminal_coupling,
        times,
        outputs,
        currents,
        histories,
        substitutions,
        predicted_terminal_currents,
        power_terminal_network_currents,
        terminal_prediction.update_count,
        d_flux,
        q_flux,
        torque,
        speed,
        angle,
        iterations,
        state.call_count,
        call_count,
        call_count,
        ordered_node_names(parsed.node_map),
        terminal_voltages,
        rotor_thevenin,
        speed_thevenin,
        torque_impedance,
        parameters.series_path_leakage_inductance_h,
        parameters.effective_armature_leakage_inductance_h,
        parameters.effective_compound_field_leakage_inductance_h,
        compensated_voltages,
        network_solve_count,
        network_correction_count,
        drive_source_values,
        isempty(excitation_sources) ? zeros(Float64, call_count) :
            vec(excitation_source_matrix[1, :]),
        isempty(excitation_source_frequencies_hz) ? 0.0 :
            excitation_source_frequencies_hz[1],
        length(excitation_sources) < 2 ? zeros(Float64, call_count) :
            vec(excitation_source_matrix[2, :]),
        length(excitation_source_frequencies_hz) < 2 ? 0.0 :
            excitation_source_frequencies_hz[2],
        (1 + length(excitation_sources)) * call_count,
        call_count,
        network_solve_count == call_count - 1 &&
            (
                !predicted_current_coupling ||
                terminal_prediction.update_count == call_count
            ),
        Symbol[],
        report_output_names,
        report_output_values,
        control_system_step_count,
    )
end

function _universal_machine_output_headers(horizon::DeckUniversalMachineHorizon)
    excitation_count = size(horizon.output_values, 1) - 6
    excitation_headers = horizon.machine_type in (1, 2) ?
        ["field_current_A", "q_axis_damper_current_A"] :
        horizon.machine_type == 6 ? ["field_current_A"] :
        horizon.machine_type == 7 ?
        ["d_axis_rotor_current_A", "q_axis_rotor_current_A"] :
        horizon.machine_type == 8 ? ["separate_field_current_A"] :
        horizon.machine_type == 9 ?
        ["compound_field_current_A", "series_field_current_A"] :
        horizon.machine_type == 10 ?
        ["inactive_compound_field_current_A", "series_field_current_A"] :
        horizon.machine_type == 11 ?
        ["parallel_compound_field_current_A", "series_field_current_A"] :
        horizon.machine_type == 12 ?
        ["shunt_field_current_A", "inactive_series_field_current_A"] :
        excitation_count == 3 ?
        ["stator_current_$(index)_A" for index in 1:excitation_count] :
        ["shorted_rotor_axis_current_$(index)_A" for index in 1:excitation_count]
    angle_header = horizon.machine_type in (1, 2) ?
        "synchronous_torque_angle_rad" : "mechanical_angle_rad"
    power_headers = horizon.machine_type in (6, 7, 8, 9, 10, 11, 12) ?
        [
            "power_zero_axis_current_A",
            "orthogonal_dummy_current_A",
            horizon.machine_type in (8, 9, 10, 11, 12) ? "armature_current_A" :
            "single_phase_stator_current_A",
        ] :
        horizon.machine_type in (2, 5) ?
        [
            "power_zero_axis_current_A",
            "armature_axis_1_current_A",
            "armature_axis_2_current_A",
        ] :
        ["phase_a_current_A", "phase_b_current_A", "phase_c_current_A"]
    return vcat([
        "torque_Nm", "mechanical_speed_rad_s", angle_header,
    ], power_headers, excitation_headers)
end

function _deck_universal_machine_trace(
    parsed::DeckParser.DeckParseResult,
    horizon::DeckUniversalMachineHorizon,
)
    length(horizon.time_s) == size(horizon.compensated_voltage_values, 2) ||
        throw(ArgumentError("universal-machine voltage horizon is incomplete"))
    length(horizon.time_s) == size(horizon.output_values, 2) ||
        throw(ArgumentError("universal-machine output horizon is incomplete"))
    voltage_extrema = _sampled_trace_extrema(
        horizon.compensated_voltage_values,
        horizon.time_s,
    )
    combined_output_values = vcat(
        horizon.report_output_values,
        horizon.output_values,
    )
    output_extrema = _sampled_trace_extrema(combined_output_values, horizon.time_s)
    timing = deck_fixed_step_horizon(parsed)
    return DeckEMTTrace(
        horizon.source,
        timing.dt_s,
        last(horizon.time_s),
        copy(parsed.node_map),
        copy(horizon.node_names),
        copy(parsed.element_names),
        copy(horizon.time_s),
        copy(horizon.compensated_voltage_values),
        vcat(
            copy(horizon.report_output_names),
            Symbol.(_universal_machine_output_headers(horizon)),
        ),
        Int[],
        combined_output_values,
        voltage_extrema.maximum_values,
        voltage_extrema.maximum_times_s,
        voltage_extrema.minimum_values,
        voltage_extrema.minimum_times_s,
        output_extrema.maximum_values,
        output_extrema.maximum_times_s,
        output_extrema.minimum_values,
        output_extrema.minimum_times_s,
    )
end

function write_universal_machine_horizon_csv(
    path::AbstractString,
    horizon::DeckUniversalMachineHorizon,
)
    output_headers = _universal_machine_output_headers(horizon)
    prediction_headers =
        horizon.requested_terminal_coupling == :predicted_current ?
        [
            "predicted_phase_a_current_injection_A",
            "predicted_phase_b_current_injection_A",
            "predicted_phase_c_current_injection_A",
            "accepted_phase_a_terminal_network_current_A",
            "accepted_phase_b_terminal_network_current_A",
            "accepted_phase_c_terminal_network_current_A",
        ] : String[]
    series_path_headers = horizon.machine_type in (9, 10, 11, 12) ? [
        "series_path_leakage_inductance_H",
        "effective_armature_leakage_inductance_H",
        "effective_compound_field_leakage_inductance_H",
    ] : String[]
    headers = vcat(["time_s"], output_headers, prediction_headers, [
        "d_axis_flux_Wb", "q_axis_flux_Wb", "iteration_count",
        "phase_a_thevenin_voltage_V", "phase_b_thevenin_voltage_V",
        "phase_c_thevenin_voltage_V", "mechanical_speed_thevenin_rad_s",
        "generated_torque_impedance_rad_s_per_Nm",
    ], series_path_headers)
    open(path, "w") do io
        println(io, join(headers, ','))
        for index in eachindex(horizon.time_s)
            values = Any[
                horizon.time_s[index],
                horizon.output_values[:, index]...,
            ]
            horizon.requested_terminal_coupling == :predicted_current &&
                append!(
                    values,
                    (
                        horizon.predicted_terminal_current_injections_A[:, index]...,
                        horizon.power_terminal_network_currents_A[:, index]...,
                    ),
                )
            append!(
                values,
                (
                    horizon.d_axis_flux[index],
                    horizon.q_axis_flux[index],
                    horizon.iteration_counts[index],
                    horizon.power_terminal_voltages[:, index]...,
                    horizon.mechanical_speed_thevenin_rad_s[index],
                    horizon.generated_torque_impedance[index],
                ),
            )
            if horizon.machine_type in (9, 10, 11, 12)
                append!(
                    values,
                    (
                        horizon.series_path_leakage_inductance_h,
                        horizon.effective_armature_leakage_inductance_h,
                        horizon.effective_compound_field_leakage_inductance_h,
                    ),
                )
            end
            println(io, join(values, ','))
        end
    end
    return path
end

function deck_output_step_indices(
    parsed::DeckParser.DeckParseResult,
    dt_s::Float64,
    t_end_s::Float64;
    schedule::Symbol = :print_and_plot,
)
    step_count = fixed_step_count(dt_s, t_end_s)
    schedule == :all_steps && return collect(0:step_count)
    schedule in (:print, :plot, :print_and_plot) ||
        throw(ArgumentError("unsupported deck output schedule $schedule"))
    options = DeckParser.deck_output_schedule_options(parsed)
    intervals = Int[]
    schedule in (:print, :print_and_plot) &&
        options.print_interval_steps > 0 &&
        push!(intervals, Int(options.print_interval_steps))
    schedule in (:plot, :print_and_plot) &&
        options.plot_interval_steps > 0 &&
        push!(intervals, Int(options.plot_interval_steps))
    isempty(intervals) &&
        throw(ArgumentError("deck does not define a positive $schedule output interval"))
    steps = Set{Int}((0, step_count))
    for interval in intervals
        for step in interval:interval:step_count
            push!(steps, step)
        end
    end
    return sort!(collect(steps))
end

function _element_trace_name(element_names::AbstractVector{Symbol}, index::Int)
    return index <= length(element_names) ?
        element_names[index] :
        Symbol("element_", index)
end

function _trace_output_channel_names(elements::Tuple, element_names::AbstractVector{Symbol})
    names = Symbol[]
    _trace_output_channel_names!(names, elements, element_names, 1)
    return names
end

_trace_output_channel_names!(names::Vector{Symbol}, elements::Tuple{}, element_names, index::Int) = names

function _trace_output_channel_names!(names::Vector{Symbol}, elements::Tuple, element_names, index::Int)
    element = first(elements)
    trace_output_channel_names!(names, _element_trace_name(element_names, index), element)
    return _trace_output_channel_names!(names, Base.tail(elements), element_names, index + 1)
end

_record_trace_outputs!(output::AbstractMatrix{Float64}, sample::Int, elements::Tuple{}, voltage, first_channel::Int) = first_channel

function _record_trace_outputs!(output::AbstractMatrix{Float64}, sample::Int, elements::Tuple, voltage, first_channel::Int)
    next_channel = trace_output_values!(output, first_channel, sample, first(elements), voltage)
    return _record_trace_outputs!(output, sample, Base.tail(elements), voltage, next_channel)
end

function _record_context_outputs!(
    output::AbstractMatrix{Float64},
    sample::Int,
    context::EMTStepContext,
    voltage,
)
    next_channel =
        _record_trace_outputs!(output, sample, context.system.elements, voltage, 1)
    for (offset, node_index) in enumerate(context.output_node_indices)
        output[next_channel + offset - 1, sample] = voltage[node_index]
    end
    next_channel += length(context.output_node_indices)
    for branch_index in context.branch_voltage_output_branch_indices
        snapshot = _deck_branch_output_snapshot(context, branch_index, voltage)
        output[next_channel, sample] =
            branch_index in context.branch_power_output_branch_indices ?
            snapshot.branch_voltage *
            _deck_branch_current_value(context, branch_index, voltage) :
            snapshot.branch_voltage
        next_channel += 1
    end
    switch_currents =
        isempty(context.switch_voltage_output_switch_indices) &&
        isempty(context.switch_current_output_switch_indices) ?
        context.switch_current_step_values :
        _deck_time_switch_report_current_values!(context, voltage, context.t_s)
    if !isempty(context.switch_voltage_output_switch_indices)
        switch_voltages = context.switch_voltage_step_values
        switch_powers =
            _deck_time_switch_report_power_values!(context, voltage, switch_currents)
        for switch_index in context.switch_voltage_output_switch_indices
            output[next_channel, sample] =
                context.deck_over5_switch_output_codes[switch_index] > 3 ?
                switch_powers[switch_index] :
                switch_voltages[switch_index]
            next_channel += 1
        end
    end
    if !isempty(context.switch_current_output_switch_indices)
        for switch_index in context.switch_current_output_switch_indices
            output[next_channel, sample] =
                context.deck_over5_switch_output_codes[switch_index] > 3 ?
                context.switch_energy_values[switch_index] :
                switch_currents[switch_index]
            next_channel += 1
        end
    end
    for branch_index in context.branch_current_output_branch_indices
        output[next_channel, sample] =
            branch_index in context.branch_power_output_branch_indices ?
            context.branch_energy_values[branch_index] :
            _deck_branch_current_value(context, branch_index, voltage)
        next_channel += 1
    end
    return _record_control_system_network_outputs!(
        output,
        sample,
        context,
        voltage,
        next_channel,
    )
end

function _update_deck_power_energy_state!(
    context::EMTStepContext,
    voltage::AbstractVector{Float64},
)
    for branch_index in context.branch_power_output_branch_indices
        snapshot = _deck_branch_output_snapshot(context, branch_index, voltage)
        power = snapshot.branch_voltage *
                _deck_branch_current_value(context, branch_index, voltage)
        if context.step_index > 0
            context.branch_energy_values[branch_index] +=
                0.5 * context.dt_s *
                (context.branch_previous_power_values[branch_index] + power)
        end
        context.branch_previous_power_values[branch_index] = power
    end
    if any(>(3), context.deck_over5_switch_output_codes)
        switch_currents =
            _deck_time_switch_report_current_values!(context, voltage, context.t_s)
        switch_powers =
            _deck_time_switch_report_power_values!(context, voltage, switch_currents)
        for switch_index in eachindex(context.deck_over5_switch_output_codes)
            context.deck_over5_switch_output_codes[switch_index] > 3 || continue
            power = switch_powers[switch_index]
            if context.step_index > 0
                context.switch_energy_values[switch_index] +=
                    0.5 * context.dt_s *
                    (context.switch_previous_power_values[switch_index] + power)
            end
            context.switch_previous_power_values[switch_index] = power
        end
    end
    return context
end

function initialize_step_context(
    lines;
    dt_s::Float64 = 20e-6,
    t_end_s::Float64 = 0.0,
    source::AbstractString = "deck",
    recorded_step_indices = nothing,
    source_signal_provider::AbstractSourceSignalProvider = IdentitySourceSignalProvider(),
)
    parsed = DeckParser.parse_deck_lines(lines; source = source)
    return initialize_step_context(
        parsed;
        dt_s = dt_s,
        t_end_s = t_end_s,
        recorded_step_indices = recorded_step_indices,
        source_signal_provider = source_signal_provider,
    )
end

function _shift_time_switch_event(
    time_s::Real,
    delay_s::Float64,
    horizon_s::Float64=Inf,
)
    event_time = Float64(time_s)
    (!isfinite(event_time) || delay_s == 0.0) && return event_time
    event_time < 0.0 && return event_time
    if event_time > horizon_s
        return event_time + delay_s
    end
    detection_time = 0.0
    while detection_time < event_time
        detection_time += delay_s
    end
    return detection_time + delay_s
end

function _delay_deck_time_switch_events!(
    elements::Vector,
    delay_s::Float64,
    horizon_s::Float64=Inf,
)
    delay_s >= 0.0 || throw(ArgumentError("time-switch event delay must be nonnegative"))
    delay_s == 0.0 && return elements
    for index in eachindex(elements)
        element = elements[index]
        element isa TimeSwitch || continue
        elements[index] = TimeSwitch(
            element.a,
            element.b;
            close_time_s =
                _shift_time_switch_event(element.close_time_s, delay_s, horizon_s),
            open_time_s =
                _shift_time_switch_event(element.open_time_s, delay_s, horizon_s),
            initially_closed = element.initially_closed,
            on_conductance = element.on_conductance,
            off_conductance = element.off_conductance,
        )
    end
    return elements
end

function _deck_current_zero_switch_names(parsed::DeckParser.DeckParseResult)
    t_end_s = Float64(DeckParser.deck_fixed_time_horizon_options(parsed).tmax_s)
    return Set(
        row.name for row in DeckParser.deck_over5_switch_rows(parsed)
        if 0.0 <= row.open_time_s <= t_end_s
    )
end

function _convert_deck_current_zero_switches!(
    elements::Vector,
    element_names::Vector{Symbol},
    parsed::DeckParser.DeckParseResult,
)
    current_zero_names = _deck_current_zero_switch_names(parsed)
    isempty(current_zero_names) && return elements
    open_request_times = Dict(
        row.name => Float64(row.open_time_s)
        for row in DeckParser.deck_over5_switch_rows(parsed)
        if row.name in current_zero_names
    )
    critical_currents = Dict(
        row.name => Float64(row.critical_current_a)
        for row in DeckParser.deck_over5_switch_rows(parsed)
        if row.name in current_zero_names
    )
    for index in eachindex(elements, element_names)
        element_names[index] in current_zero_names || continue
        element = elements[index]
        element isa TimeSwitch || throw(ArgumentError(
            "current-zero switch row must own a time-switch model",
        ))
        current_zero = CurrentZeroSwitch(
            element;
            critical_current_a = critical_currents[element_names[index]],
        )
        current_zero.open_request_time_s = open_request_times[element_names[index]]
        elements[index] = current_zero
    end
    return elements
end

function _deck_runtime_output_context_kwargs(
    parsed::DeckParser.DeckParseResult;
    time_switch_event_delay_s::Float64 = 0.0,
    event_horizon_s::Float64 = Inf,
)
    time_switch_event_delay_s >= 0.0 ||
        throw(ArgumentError("time-switch event delay must be nonnegative"))
    return (
        deck_output_channel_names = DeckParser.deck_over16_output_channel_names(parsed),
        deck_output_node_indices = DeckParser.deck_over16_output_node_indices(parsed),
        deck_branch_voltage_output_names =
            DeckParser.deck_over16_branch_voltage_output_names(parsed),
        deck_branch_voltage_output_branch_indices =
            DeckParser.deck_over16_branch_voltage_branch_indices(parsed),
        deck_branch_current_output_names =
            DeckParser.deck_over16_branch_current_output_names(parsed),
        deck_branch_current_output_branch_indices =
            DeckParser.deck_over16_branch_current_branch_indices(parsed),
        deck_branch_power_output_branch_indices =
            DeckParser.deck_over16_branch_power_branch_indices(parsed),
        deck_time_switch_names = DeckParser.deck_time_switch_names(parsed),
        deck_time_switch_from_node_indices =
            DeckParser.deck_time_switch_from_node_indices(parsed),
        deck_time_switch_to_node_indices =
            DeckParser.deck_time_switch_to_node_indices(parsed),
        deck_time_switch_close_time_s_values =
            _shift_time_switch_event.(
                DeckParser.deck_time_switch_close_time_s_values(parsed),
                time_switch_event_delay_s,
                event_horizon_s,
            ),
        deck_time_switch_open_time_s_values =
            _shift_time_switch_event.(
                DeckParser.deck_time_switch_open_time_s_values(parsed),
                time_switch_event_delay_s,
                event_horizon_s,
            ),
        deck_time_switch_initially_closed_flags =
            DeckParser.deck_time_switch_initially_closed_flags(parsed),
        deck_time_switch_on_conductance_values =
            DeckParser.deck_time_switch_on_conductance_values(parsed),
        deck_time_switch_off_conductance_values =
            DeckParser.deck_time_switch_off_conductance_values(parsed),
        deck_over5_switch_output_codes =
            DeckParser.deck_over5_switch_output_codes(parsed),
    )
end

function initialize_step_context(
    parsed::DeckParser.DeckParseResult;
    dt_s::Float64 = 20e-6,
    t_end_s::Float64 = 0.0,
    recorded_step_indices = nothing,
    time_switch_event_delay_s::Float64 = 0.0,
    current_zero_switching::Bool = false,
    source_signal_provider::AbstractSourceSignalProvider = IdentitySourceSignalProvider(),
)
    DeckParser.assert_deck_valid!(parsed)
    elements = Any[parsed.elements...]
    element_names = copy(parsed.element_names)
    sampled_frequency_lines =
        DeckParser.deck_sampled_frequency_line_elements(parsed, dt_s)
    append!(elements, sampled_frequency_lines)
    append!(
        element_names,
        DeckParser.deck_sampled_frequency_line_element_names(parsed),
    )
    semlyen_lines = DeckParser.deck_semlyen_line_elements(parsed, dt_s)
    append!(elements, semlyen_lines)
    append!(element_names, DeckParser.deck_semlyen_line_element_names(parsed))
    rational_frequency_lines =
        DeckParser.deck_rational_frequency_line_elements(parsed, dt_s)
    append!(elements, rational_frequency_lines)
    append!(
        element_names,
        DeckParser.deck_rational_frequency_line_element_names(parsed),
    )
    nonlinear_inductor_rows =
        DeckParser.deck_piecewise_nonlinear_inductor_rows(parsed)
    if any(
        row -> row.nonlinear_type == PSEUDO_NONLINEAR_INDUCTOR_TYPE,
        nonlinear_inductor_rows,
    )
        nonlinear_inductor_config =
            deck_piecewise_nonlinear_inductor_current_config(
                parsed;
                delta2 = dt_s / 2.0,
            )
        nodal_nonlinear_inductor_config = merge(
            nonlinear_inductor_config,
            (
                nonlinear_from_nodes =
                    [row.from_node_index for row in nonlinear_inductor_rows],
                nonlinear_to_nodes =
                    [row.to_node_index for row in nonlinear_inductor_rows],
            ),
        )
        nonlinear_slope_branches =
            saturated_transformer_nonlinear_slope_branches(
                nodal_nonlinear_inductor_config,
            )
        append!(elements, nonlinear_slope_branches.elements)
        append!(element_names, nonlinear_slope_branches.element_names)
    end
    _append_switching_nonlinear_resistor_safety_shunts!(
        elements,
        element_names,
        parsed,
    )
    _delay_deck_time_switch_events!(elements, time_switch_event_delay_s, t_end_s)
    current_zero_switching &&
        _convert_deck_current_zero_switches!(elements, element_names, parsed)
    source_function_runtime, control_system_runtime =
        _append_dynamic_source_and_control_elements!(
        elements,
        element_names,
        parsed,
        dt_s,
        source_signal_provider,
    )
    system = NodalSystem(maximum(values(parsed.node_map); init = 0), elements)
    return initialize_step_context(
        system;
        node_map = parsed.node_map,
        element_names = element_names,
        source_function_runtime = source_function_runtime,
        control_system_runtime = control_system_runtime,
        _deck_runtime_output_context_kwargs(
            parsed;
            time_switch_event_delay_s = time_switch_event_delay_s,
            event_horizon_s = t_end_s,
        )...,
        dt_s = dt_s,
        t_end_s = t_end_s,
        source = parsed.source,
        recorded_step_indices = recorded_step_indices,
    )
end

function _control_system_device_row_executable(row::DeckParser.DeckControlSystemDeviceRow)
    50 <= row.device_type <= 67 || return false
    length(row.parameter_values) == 3 || return false
    input_required = row.device_type in (
        50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 62, 64, 65, 66, 67,
    )
    input_required && isempty(row.input_terms) && return false
    row.device_type == 60 && length(row.input_terms) < 3 && return false
    if row.device_type in (55, 56, 57)
        row.table_complete && !isempty(row.table_input_values) || return false
    end
    if row.device_type == 56
        length(row.table_output_values) == length(row.table_input_values) || return false
    elseif row.device_type == 53
        maximum_delay = row.control_signal === missing ?
            row.parameter_values[2] : row.parameter_values[3]
        maximum_delay >= 0.0 || return false
    elseif row.device_type == 66
        row.parameter_values[1] > 0.0 || return false
    elseif row.device_type == 67
        row.parameter_values[1] != 0.0 || return false
    end
    return true
end

function deck_over16_boundary_plan(parsed::DeckParser.DeckParseResult)
    DeckParser.assert_deck_valid!(parsed)
    source_inputs = DeckParser.deck_over5a_source_update_inputs(parsed)
    control_device_rows = DeckParser.deck_control_system_device_rows(parsed)
    full_control_device_execution =
        !isempty(control_device_rows) &&
        all(_control_system_device_row_executable, control_device_rows) &&
        all(
            row -> row.coefficients_complete,
            DeckParser.deck_control_system_function_rows(parsed),
        )
    deferred_effects = Symbol[
        :full_bpa_deck_grammar,
        :full_init_step_calc_elec_orchestration,
        :full_distributed_line_card_grammar,
        :multiphase_modal_transform,
        :frequency_dependent_line_fitting,
        :report_file_writers,
        :external_bpa_executable_waveform_comparison,
    ]
    full_control_device_execution || insert!(deferred_effects, 7, :tacs_execution)
    return DeckOVER16BoundaryPlan(
        :deck_owned_over16_boundary_plan,
        DeckParser.deck_branch_names(parsed),
        DeckParser.deck_branch_kinds(parsed),
        DeckParser.deck_branch_from_node_names(parsed),
        DeckParser.deck_branch_to_node_names(parsed),
        DeckParser.deck_branch_from_node_indices(parsed),
        DeckParser.deck_branch_to_node_indices(parsed),
        DeckParser.deck_branch_conductance_values(parsed),
        DeckParser.deck_branch_resistance_values(parsed),
        DeckParser.deck_branch_inductance_values(parsed),
        DeckParser.deck_branch_capacitance_values(parsed),
        DeckParser.deck_branch_previous_current_values(parsed),
        DeckParser.deck_branch_previous_voltage_values(parsed),
        DeckParser.deck_branch_line_numbers(parsed),
        DeckParser.deck_branch_count(parsed),
        DeckParser.deck_bergeron_line_names(parsed),
        DeckParser.deck_bergeron_line_line_numbers(parsed),
        DeckParser.deck_bergeron_line_from_node_names(parsed),
        DeckParser.deck_bergeron_line_to_node_names(parsed),
        DeckParser.deck_bergeron_line_from_node_indices(parsed),
        DeckParser.deck_bergeron_line_to_node_indices(parsed),
        DeckParser.deck_bergeron_line_surge_impedance_values(parsed),
        DeckParser.deck_bergeron_line_surge_admittance_values(parsed),
        DeckParser.deck_bergeron_line_travel_time_s_values(parsed),
        DeckParser.deck_bergeron_line_dt_s_values(parsed),
        DeckParser.deck_bergeron_line_attenuation_values(parsed),
        DeckParser.deck_bergeron_line_delay_step_counts(parsed),
        DeckParser.deck_bergeron_line_write_indices(parsed),
        DeckParser.deck_bergeron_line_history_current_from_values(parsed),
        DeckParser.deck_bergeron_line_history_current_to_values(parsed),
        DeckParser.deck_bergeron_line_terminal_voltage_from_values(parsed),
        DeckParser.deck_bergeron_line_terminal_voltage_to_values(parsed),
        DeckParser.deck_bergeron_line_terminal_current_from_values(parsed),
        DeckParser.deck_bergeron_line_terminal_current_to_values(parsed),
        DeckParser.deck_bergeron_line_traveling_wave_from_values(parsed),
        DeckParser.deck_bergeron_line_traveling_wave_to_values(parsed),
        length(DeckParser.deck_bergeron_line_rows(parsed)),
        DeckParser.deck_over2_branch_names(parsed),
        DeckParser.deck_over2_branch_line_numbers(parsed),
        DeckParser.deck_over2_branch_kinds(parsed),
        DeckParser.deck_over2_branch_layout_kinds(parsed),
        DeckParser.deck_over2_branch_source_kinds(parsed),
        DeckParser.deck_over2_branch_reference_kinds(parsed),
        DeckParser.deck_over2_branch_reference_names(parsed),
        DeckParser.deck_over2_branch_reference_line_numbers(parsed),
        DeckParser.deck_over2_branch_from_node_names(parsed),
        DeckParser.deck_over2_branch_to_node_names(parsed),
        DeckParser.deck_over2_branch_from_node_indices(parsed),
        DeckParser.deck_over2_branch_to_node_indices(parsed),
        DeckParser.deck_over2_branch_raw_resistance_values(parsed),
        DeckParser.deck_over2_branch_raw_inductance_values(parsed),
        DeckParser.deck_over2_branch_raw_capacitance_values(parsed),
        DeckParser.deck_over2_branch_conductance_values(parsed),
        DeckParser.deck_over2_branch_resistance_values(parsed),
        DeckParser.deck_over2_branch_inductance_values(parsed),
        DeckParser.deck_over2_branch_capacitance_values(parsed),
        DeckParser.deck_over2_branch_output_codes(parsed),
        length(DeckParser.deck_over2_branch_rows(parsed)),
        DeckParser.deck_over16_output_channel_names(parsed),
        DeckParser.deck_over16_output_node_names(parsed),
        DeckParser.deck_over16_output_node_indices(parsed),
        DeckParser.deck_over16_output_channel_line_numbers(parsed),
        DeckParser.deck_over16_branch_voltage_output_names(parsed),
        DeckParser.deck_over16_branch_voltage_branch_names(parsed),
        DeckParser.deck_over16_branch_voltage_branch_indices(parsed),
        DeckParser.deck_over16_branch_voltage_output_line_numbers(parsed),
        DeckParser.deck_over16_branch_current_output_names(parsed),
        DeckParser.deck_over16_branch_current_branch_names(parsed),
        DeckParser.deck_over16_branch_current_branch_indices(parsed),
        DeckParser.deck_over16_branch_current_output_line_numbers(parsed),
        DeckParser.deck_over16_branch_power_output_names(parsed),
        DeckParser.deck_over16_branch_power_branch_names(parsed),
        DeckParser.deck_over16_branch_power_branch_indices(parsed),
        DeckParser.deck_over16_branch_power_output_line_numbers(parsed),
        DeckParser.deck_over15_output_request_names(parsed),
        DeckParser.deck_over15_output_request_output_kinds(parsed),
        DeckParser.deck_over15_output_request_request_kinds(parsed),
        DeckParser.deck_over15_output_request_layout_kinds(parsed),
        DeckParser.deck_over15_output_request_line_numbers(parsed),
        DeckParser.deck_over15_output_request_output_codes(parsed),
        DeckParser.deck_over15_output_request_node_names(parsed),
        DeckParser.deck_over15_output_request_node_indices(parsed),
        DeckParser.deck_over15_output_request_branch_names(parsed),
        DeckParser.deck_over15_output_request_branch_indices(parsed),
        DeckParser.deck_time_switch_names(parsed),
        DeckParser.deck_time_switch_line_numbers(parsed),
        DeckParser.deck_time_switch_from_node_names(parsed),
        DeckParser.deck_time_switch_to_node_names(parsed),
        DeckParser.deck_time_switch_from_node_indices(parsed),
        DeckParser.deck_time_switch_to_node_indices(parsed),
        DeckParser.deck_time_switch_close_time_s_values(parsed),
        DeckParser.deck_time_switch_open_time_s_values(parsed),
        DeckParser.deck_time_switch_initially_closed_flags(parsed),
        DeckParser.deck_time_switch_on_conductance_values(parsed),
        DeckParser.deck_time_switch_off_conductance_values(parsed),
        DeckParser.deck_time_switch_count(parsed),
        DeckParser.deck_over5_switch_names(parsed),
        DeckParser.deck_over5_switch_line_numbers(parsed),
        DeckParser.deck_over5_switch_from_node_names(parsed),
        DeckParser.deck_over5_switch_to_node_names(parsed),
        DeckParser.deck_over5_switch_from_node_indices(parsed),
        DeckParser.deck_over5_switch_to_node_indices(parsed),
        DeckParser.deck_over5_switch_layout_kinds(parsed),
        DeckParser.deck_over5_switch_raw_close_time_s_values(parsed),
        DeckParser.deck_over5_switch_raw_open_time_s_values(parsed),
        DeckParser.deck_over5_switch_close_time_s_values(parsed),
        DeckParser.deck_over5_switch_open_time_s_values(parsed),
        DeckParser.deck_over5_switch_initially_closed_flags(parsed),
        DeckParser.deck_over5_switch_measuring_flags(parsed),
        DeckParser.deck_over5_switch_closed_markers(parsed),
        DeckParser.deck_over5_switch_marker_texts(parsed),
        DeckParser.deck_over5_switch_type_values(parsed),
        DeckParser.deck_over5_switch_critical_current_values(parsed),
        DeckParser.deck_over5_switch_random_opening_standard_deviation_s_values(parsed),
        DeckParser.deck_over5_switch_on_conductance_values(parsed),
        DeckParser.deck_over5_switch_off_conductance_values(parsed),
        DeckParser.deck_over5_switch_output_codes(parsed),
        length(DeckParser.deck_over5_switch_rows(parsed)),
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
        DeckParser.deck_over16_source_card_kinds(parsed),
        DeckParser.deck_over16_source_card_values(parsed),
        DeckParser.deck_over16_source_card_provided_value_counts(parsed),
        DeckParser.deck_over16_source_card_line_numbers(parsed),
        length(DeckParser.deck_over16_source_card_kinds(parsed)),
        DeckParser.deck_over16_source_interpolation_values(parsed),
        DeckParser.deck_over16_source_interpolation_provided_value_counts(parsed),
        DeckParser.deck_over16_source_interpolation_line_numbers(parsed),
        length(DeckParser.deck_over16_source_interpolation_values(parsed)),
        DeckParser.deck_over16_source_tacs_override_positions(parsed),
        DeckParser.deck_over16_source_tacs_override_xtcs_indices(parsed),
        DeckParser.deck_over16_source_tacs_override_line_numbers(parsed),
        length(DeckParser.deck_over16_source_tacs_override_positions(parsed)),
        DeckParser.deck_over16_source_analytic_values(parsed),
        DeckParser.deck_over16_source_analytic_provided_value_counts(parsed),
        DeckParser.deck_over16_source_analytic_line_numbers(parsed),
        length(DeckParser.deck_over16_source_analytic_values(parsed)),
        (
            :over16_output_report_labels_1003_1658,
            :over16_branch_power_energy_labels_1600_1642,
            :over16_output_writer_labels_1643_1654,
            :over16_post_extrema_timestep_advance_labels_1647_1661,
            :over16_source_post_advance_time_labels_1661_1249,
            :over16_source_card_labels_1945_1247_11247,
            :over16_source_supplied_signal_labels_11247_11248,
            :over5a_fixed_source_labels_332_333_349_362_6304_6358,
            :over16_source_row_labels_1249_1300,
            :over16_source_row_labels_7362_7368_1258,
            :over5_switch_time_labels_3483_3506_218_220_216,
            :over5_switch_measuring_labels_218_70218,
            :over2_branch_labels_1_8331_64114,
            :over2_branch_output_code_labels_144_54208,
            :over2_branch_copy_reference_labels_175_177_21696_54111,
            :over2_over10_over13_bergeron_line_delay_owner,
        ),
        Tuple(deferred_effects),
    )
end

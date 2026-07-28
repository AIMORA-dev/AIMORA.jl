struct DeckDirectMachineFleetHorizon
    source::String
    machine_types::Vector{Int}
    time_s::Vector{Float64}
    output_values::Array{Float64,3}
    current_values::Array{Float64,3}
    history_currents::Array{Float64,3}
    current_substitution_values::Array{Float64,3}
    d_axis_flux::Matrix{Float64}
    q_axis_flux::Matrix{Float64}
    generated_torque::Matrix{Float64}
    mechanical_speed_rad_s::Matrix{Float64}
    mechanical_angle_rad::Matrix{Float64}
    iteration_count::Matrix{Int}
    call_counts::Vector{Int}
    node_names::Vector{Symbol}
    terminal_voltage_values::Array{Float64,3}
    mechanical_speed_thevenin_rad_s::Matrix{Float64}
    compensation_impedances::Array{Float64,3}
    cross_machine_impedance::Vector{Float64}
    electrical_cross_impedance::Vector{Float64}
    mechanical_cross_impedance::Vector{Float64}
    coupling_iteration_counts::Vector{Int}
    coupling_residuals::Vector{Float64}
    compensated_voltage_values::Matrix{Float64}
    coil_control_signal_names::Matrix{Union{Missing,Symbol}}
    coil_control_voltages::Array{Float64,3}
    control_output_names::Vector{Symbol}
    control_output_values::Matrix{Float64}
    control_system_execution_count::Int
    initialized_drive_crests::Vector{Float64}
    network_solve_count::Int
    network_correction_count::Int
    terminal_rhs_mutation_count::Int
    network_control_mutation_count::Int
    report_mutation_count::Int
    complete_machine_path::Bool
    deferred_effects::Vector{Symbol}
end

function _direct_machine_fleet_cross_impedance(
    impedance::AbstractMatrix{<:Real},
    machine_count::Int,
)
    maximum_cross = 0.0
    for left in 1:machine_count, right in 1:machine_count
        left == right && continue
        for left_position in (left, machine_count + left),
            right_position in (right, machine_count + right)
            maximum_cross = max(
                maximum_cross,
                abs(Float64(impedance[left_position, right_position])),
            )
        end
    end
    return maximum_cross
end

function _direct_machine_fleet_cross_impedance_components(
    impedance::AbstractMatrix{<:Real},
    machine_count::Int,
)
    electrical = 0.0
    mechanical = 0.0
    for left in 1:machine_count, right in 1:machine_count
        left == right && continue
        electrical = max(electrical, abs(Float64(impedance[left, right])))
        mechanical = max(
            mechanical,
            abs(Float64(impedance[machine_count + left, machine_count + right])),
        )
    end
    return (; electrical, mechanical)
end

function run_deck_direct_machine_fleet_horizon(
    parsed::DeckParser.DeckParseResult;
    time_step_s::Real=deck_fixed_step_horizon(parsed).dt_s,
    coupling_absolute_tolerance::Real=1.0e-10,
    coupling_relative_tolerance::Real=1.0e-10,
    coupling_max_iterations::Int=50,
    coupling_relaxation::Real=0.7,
)
    DeckParser.assert_deck_valid!(parsed)
    sections = DeckParser.deck_universal_machine_section_rows(parsed)
    length(sections) == 1 ||
        throw(ArgumentError("direct-machine fleet requires one universal-machine section"))
    section = only(sections)
    section.initialization_mode == :manual ||
        throw(ArgumentError("direct-machine fleet requires manual initialization"))
    machine_count = section.machine_count
    machine_count > 1 ||
        throw(ArgumentError("direct-machine fleet requires at least two machines"))

    horizon = deck_fixed_step_horizon(parsed)
    dt_s = Float64(time_step_s)
    dt_s == horizon.dt_s ||
        throw(ArgumentError("direct-machine fleet time step must match the deck time step"))
    absolute_tolerance = Float64(coupling_absolute_tolerance)
    relative_tolerance = Float64(coupling_relative_tolerance)
    relaxation = Float64(coupling_relaxation)
    absolute_tolerance >= 0.0 && isfinite(absolute_tolerance) ||
        throw(ArgumentError("coupling_absolute_tolerance must be finite and nonnegative"))
    relative_tolerance >= 0.0 && isfinite(relative_tolerance) ||
        throw(ArgumentError("coupling_relative_tolerance must be finite and nonnegative"))
    coupling_max_iterations > 0 ||
        throw(ArgumentError("coupling_max_iterations must be positive"))
    0.0 < relaxation <= 1.0 ||
        throw(ArgumentError("coupling_relaxation must be in (0, 1]"))
    parameters = [
        deck_direct_current_machine_parameters(
            parsed;
            machine_index,
            time_step_s = dt_s,
        ) for machine_index in 1:machine_count
    ]
    machine_types = [parameter.machine_type for parameter in parameters]
    all(machine_type -> machine_type in (9, 10), machine_types) ||
        throw(ArgumentError("direct-machine fleet currently owns machine types 9 and 10"))
    all(parameter -> length(parameter.coil_conductances) == 5, parameters) ||
        throw(ArgumentError("direct-machine fleet requires five-coil owners"))

    steady_state = deck_steady_state_voltage_phasors(parsed)
    states = [
        deck_direct_current_machine_initial_state(
            parsed;
            machine_index,
            steady_state,
        ) for machine_index in 1:machine_count
    ]
    network_nodes = [
        _deck_universal_machine_network_nodes(parsed, machine_index)
        for machine_index in 1:machine_count
    ]
    all(nodes -> length(nodes.active_power_positions) == 1, network_nodes) ||
        throw(ArgumentError("direct-machine fleet requires one active armature terminal per machine"))

    context = initialize_step_context(
        parsed;
        dt_s,
        t_end_s = horizon.t_end_s,
        recorded_step_indices = [0, horizon.step_count],
    )
    _seed_steady_state_network_state!(context, steady_state)
    for machine_index in 1:machine_count
        _seed_machine_mechanical_inertia!(
            context.system,
            network_nodes[machine_index].mechanical,
            states[machine_index].mechanical_speed_rad_s,
        )
    end

    call_count = horizon.step_count + 1
    times = Float64[step * dt_s for step in 0:horizon.step_count]
    outputs = zeros(Float64, 8, machine_count, call_count)
    currents = zeros(Float64, 5, machine_count, call_count)
    histories = zeros(Float64, 5, machine_count, call_count)
    substitutions = zeros(Float64, 2, machine_count, call_count)
    d_flux = zeros(Float64, machine_count, call_count)
    q_flux = zeros(Float64, machine_count, call_count)
    torque = zeros(Float64, machine_count, call_count)
    speed = zeros(Float64, machine_count, call_count)
    angle = zeros(Float64, machine_count, call_count)
    iterations = zeros(Int, machine_count, call_count)
    terminal_voltages = zeros(Float64, 3, machine_count, call_count)
    mechanical_speed_thevenin = zeros(Float64, machine_count, call_count)
    compensation_impedances = zeros(
        Float64,
        2 * machine_count,
        2 * machine_count,
        call_count,
    )
    cross_machine_impedance = zeros(Float64, call_count)
    electrical_cross_impedance = zeros(Float64, call_count)
    mechanical_cross_impedance = zeros(Float64, call_count)
    coupling_iteration_counts = ones(Int, call_count)
    coupling_residuals = zeros(Float64, call_count)
    compensated_voltages = zeros(Float64, context.system.node_count, call_count)
    coil_control_signal_names = Matrix{Union{Missing,Symbol}}(
        missing,
        5,
        machine_count,
    )
    for row in DeckParser.deck_universal_machine_coil_rows(parsed)
        row.machine_index <= machine_count && row.coil_index <= 5 || continue
        coil_control_signal_names[row.coil_index, row.machine_index] = row.control_signal
    end
    controlled_coil_count = count(value -> !ismissing(value), coil_control_signal_names)
    control_runtime = controlled_coil_count == 0 ? nothing :
                      _deck_control_system_network_runtime(parsed, dt_s)
    controlled_coil_count == 0 || control_runtime !== nothing ||
        throw(ArgumentError("controlled universal-machine coils require a TACS section"))
    control_runtime === nothing ||
        _initialize_control_system_network_steady_state!(
            control_runtime,
            context.system.v,
        )
    control_output_names = sort!(
        unique(Symbol[
            signal_name for signal_name in coil_control_signal_names
            if !ismissing(signal_name)
        ]);
        by = String,
    )
    coil_control_voltages = zeros(Float64, 5, machine_count, call_count)
    control_output_values = zeros(Float64, length(control_output_names), call_count)
    control_system_execution_count = 0

    function execute_controls!(sample_index::Int)
        control_runtime === nothing && return
        if sample_index > 1
            _advance_control_system_network_runtime!(
                control_runtime,
                sample_index - 1,
                times[sample_index],
            )
            control_system_execution_count += 1
        end
        for machine_index in 1:machine_count, coil_index in 1:5
            signal_name = coil_control_signal_names[coil_index, machine_index]
            ismissing(signal_name) && continue
            haskey(control_runtime.state.values, signal_name) ||
                throw(ArgumentError("missing TACS coil-control signal $signal_name"))
            coil_control_voltages[coil_index, machine_index, sample_index] =
                control_runtime.state.values[signal_name]
        end
        for (output_index, output_name) in enumerate(control_output_names)
            control_output_values[output_index, sample_index] =
                control_runtime.state.values[output_name]
        end
        return
    end

    function record_result!(machine_index::Int, sample_index::Int, result)
        outputs[:, machine_index, sample_index] .= result.output_values
        currents[:, machine_index, sample_index] .= result.current_values
        histories[:, machine_index, sample_index] .= result.history_currents
        substitutions[:, machine_index, sample_index] .=
            result.current_substitution_values
        d_flux[machine_index, sample_index] = result.d_axis_flux
        q_flux[machine_index, sample_index] = result.q_axis_flux
        torque[machine_index, sample_index] = result.generated_torque
        speed[machine_index, sample_index] = result.mechanical_speed_rad_s
        angle[machine_index, sample_index] = result.mechanical_angle_rad
        iterations[machine_index, sample_index] = result.iteration_count
        return result
    end

    execute_controls!(1)
    initial_results = Vector{Any}(undef, machine_count)
    for machine_index in 1:machine_count
        result = coupled_dq_machine_step!(
            states[machine_index],
            parameters[machine_index];
            power_terminal_voltages = zeros(3),
            rotor_thevenin_matrix = zeros(3, 3),
            mechanical_speed_thevenin_rad_s =
                states[machine_index].mechanical_speed_rad_s,
            generated_torque_impedance = 0.0,
            stator_terminal_voltages = zeros(2),
            stator_thevenin_matrix = zeros(2, 2),
            coil_control_voltages = coil_control_voltages[:, machine_index, 1],
            initial_step = true,
        )
        initial_results[machine_index] = result
        record_result!(machine_index, 1, result)
        mechanical_node = network_nodes[machine_index].mechanical
        compensated_voltages[mechanical_node, 1] =
            steady_state.node_voltage_values[mechanical_node]
        mechanical_speed_thevenin[machine_index, 1] =
            steady_state.node_voltage_values[mechanical_node]
    end

    power_nodes = Int[
        nodes.power[only(nodes.active_power_positions)] for nodes in network_nodes
    ]
    mechanical_nodes = Int[nodes.mechanical for nodes in network_nodes]
    drive_nodes = Int[nodes.drive for nodes in network_nodes]
    compensation_nodes = vcat(power_nodes, mechanical_nodes)
    calibration_nodes = vcat(compensation_nodes, drive_nodes)
    calibration = solve_step_with_compensated_current_injections!(
        deepcopy(context.system),
        dt_s,
        dt_s,
        zeros(Float64, context.system.node_count),
        calibration_nodes,
        (_voltage, _impedance) -> zeros(Float64, length(calibration_nodes)),
    )
    mechanical_positions = (machine_count + 1):(2 * machine_count)
    drive_positions = (2 * machine_count + 1):(3 * machine_count)
    torque_currents = -Float64[result.generated_torque for result in initial_results]
    desired_speeds = Float64[
        steady_state.node_voltage_values[node] for node in mechanical_nodes
    ]
    drive_transfer = calibration.compensation_impedance[
        mechanical_positions,
        drive_positions,
    ]
    drive_residual = desired_speeds .-
        calibration.open_circuit_voltage[mechanical_nodes] .-
        calibration.compensation_impedance[
            mechanical_positions,
            mechanical_positions,
        ] * torque_currents
    drive_deltas = drive_transfer \ drive_residual
    initialized_drive_crests = zeros(Float64, machine_count)
    for machine_index in 1:machine_count
        source = network_nodes[machine_index].drive_source
        original_value = _machine_drive_source_value(source, source.crest, dt_s)
        unit_value = _machine_drive_source_value(source, 1.0, dt_s)
        abs(unit_value) > eps(Float64) ||
            throw(ArgumentError("machine $machine_index drive source is zero at the first dynamic step"))
        initialized_drive_crests[machine_index] =
            (original_value + drive_deltas[machine_index]) / unit_value
    end

    network_solve_count = 0
    previous_compensation_currents = vcat(
        Float64[result.output_values[6] for result in initial_results],
        -Float64[result.generated_torque for result in initial_results],
    )
    for sample_index in 2:call_count
        time_s = times[sample_index]
        execute_controls!(sample_index)
        base_injections = zeros(Float64, context.system.node_count)
        for machine_index in 1:machine_count
            source = network_nodes[machine_index].drive_source
            base_injections[drive_nodes[machine_index]] =
                _machine_drive_source_value(
                    source,
                    initialized_drive_crests[machine_index],
                    time_s,
                ) - _machine_drive_source_value(source, source.crest, time_s)
        end
        results = Vector{Any}(undef, machine_count)
        compensation_current = function (open_circuit_voltage, impedance)
            cross_impedance =
                _direct_machine_fleet_cross_impedance(impedance, machine_count)
            cross_machine_impedance[sample_index] = cross_impedance
            components = _direct_machine_fleet_cross_impedance_components(
                impedance,
                machine_count,
            )
            electrical_cross_impedance[sample_index] = components.electrical
            mechanical_cross_impedance[sample_index] = components.mechanical
            for machine_index in 1:machine_count
                mechanical_position = machine_count + machine_index
                max(
                    abs(impedance[machine_index, mechanical_position]),
                    abs(impedance[mechanical_position, machine_index]),
                ) <= 1.0e-12 || throw(ArgumentError(
                    "machine $machine_index has unsupported within-owner electrical/mechanical network transfer",
                ))
            end

            base_states = deepcopy.(states)
            function evaluate_coupled_currents(guess)
                trial_states = deepcopy.(base_states)
                trial_results = Vector{Any}(undef, machine_count)
                trial_currents = zeros(Float64, 2 * machine_count)
                for machine_index in 1:machine_count
                    nodes = network_nodes[machine_index]
                    active_position = only(nodes.active_power_positions)
                    mechanical_position = machine_count + machine_index
                    local_positions = (machine_index, mechanical_position)
                    power_thevenin = open_circuit_voltage[nodes.power[active_position]]
                    speed_thevenin = open_circuit_voltage[nodes.mechanical]
                    for position in eachindex(guess)
                        position in local_positions && continue
                        power_thevenin += impedance[machine_index, position] * guess[position]
                        speed_thevenin += impedance[mechanical_position, position] * guess[position]
                    end
                    stored_power_voltages = zeros(Float64, 3)
                    stored_power_voltages[active_position] = power_thevenin
                    rotor_thevenin = zeros(Float64, 3, 3)
                    rotor_thevenin[active_position, active_position] =
                        impedance[machine_index, machine_index]
                    result = coupled_dq_machine_step!(
                        trial_states[machine_index],
                        parameters[machine_index];
                        power_terminal_voltages = stored_power_voltages,
                        rotor_thevenin_matrix = rotor_thevenin,
                        mechanical_speed_thevenin_rad_s = speed_thevenin,
                        generated_torque_impedance =
                            impedance[mechanical_position, mechanical_position],
                        stator_terminal_voltages = zeros(2),
                        stator_thevenin_matrix = zeros(2, 2),
                        coil_control_voltages =
                            coil_control_voltages[:, machine_index, sample_index],
                        initial_step = false,
                    )
                    trial_results[machine_index] = result
                    trial_currents[machine_index] =
                        result.output_values[3 + active_position]
                    trial_currents[mechanical_position] = -result.generated_torque
                end
                return trial_currents, trial_states, trial_results
            end

            if cross_impedance <= 1.0e-12
                converged_currents, converged_states, converged_results =
                    evaluate_coupled_currents(zeros(Float64, 2 * machine_count))
                coupling_iteration_counts[sample_index] = 1
                coupling_residuals[sample_index] = 0.0
            else
                guess = copy(previous_compensation_currents)
                converged_currents = similar(guess)
                converged_states = deepcopy.(base_states)
                converged_results = Vector{Any}(undef, machine_count)
                converged = false
                for iteration in 1:coupling_max_iterations
                    next_currents, trial_states, trial_results =
                        evaluate_coupled_currents(guess)
                    residual = maximum(abs.(next_currents .- guess); init = 0.0)
                    allowance = absolute_tolerance + relative_tolerance *
                        maximum(abs, next_currents; init = 0.0)
                    coupling_iteration_counts[sample_index] = iteration
                    coupling_residuals[sample_index] = residual
                    if residual <= allowance
                        converged_currents .= next_currents
                        converged_states = trial_states
                        converged_results = trial_results
                        converged = true
                        break
                    end
                    guess .= (1.0 - relaxation) .* guess .+ relaxation .* next_currents
                end
                converged || throw(ArgumentError(
                    "direct-machine fleet current coupling failed to converge at t=$time_s; residual=$(coupling_residuals[sample_index])",
                ))
            end
            for machine_index in 1:machine_count
                states[machine_index] = converged_states[machine_index]
                results[machine_index] = converged_results[machine_index]
            end
            previous_compensation_currents .= converged_currents
            return converged_currents
        end
        network_result = solve_step_with_compensated_current_injections!(
            context.system,
            time_s,
            dt_s,
            base_injections,
            compensation_nodes,
            compensation_current,
        )
        for machine_index in 1:machine_count
            record_result!(machine_index, sample_index, results[machine_index])
            nodes = network_nodes[machine_index]
            terminal_voltages[:, machine_index, sample_index] .= Float64[
                node > 0 ? network_result.open_circuit_voltage[node] : 0.0
                for node in nodes.power
            ]
            mechanical_speed_thevenin[machine_index, sample_index] =
                network_result.open_circuit_voltage[nodes.mechanical]
        end
        compensation_impedances[:, :, sample_index] .=
            network_result.compensation_impedance
        compensated_voltages[:, sample_index] .= network_result.voltage
        network_solve_count += 1
    end

    call_counts = Int[state.call_count for state in states]
    complete_controls = control_runtime === nothing ||
        control_system_execution_count == call_count - 1
    complete_path = all(==(call_count), call_counts) &&
        network_solve_count == call_count - 1 && complete_controls
    return DeckDirectMachineFleetHorizon(
        parsed.source,
        machine_types,
        times,
        outputs,
        currents,
        histories,
        substitutions,
        d_flux,
        q_flux,
        torque,
        speed,
        angle,
        iterations,
        call_counts,
        ordered_node_names(parsed.node_map),
        terminal_voltages,
        mechanical_speed_thevenin,
        compensation_impedances,
        cross_machine_impedance,
        electrical_cross_impedance,
        mechanical_cross_impedance,
        coupling_iteration_counts,
        coupling_residuals,
        compensated_voltages,
        coil_control_signal_names,
        coil_control_voltages,
        control_output_names,
        control_output_values,
        control_system_execution_count,
        initialized_drive_crests,
        network_solve_count,
        network_solve_count,
        machine_count * call_count,
        machine_count * call_count,
        machine_count * call_count,
        complete_path,
        Symbol[],
    )
end

function write_direct_machine_fleet_csv(
    path::AbstractString,
    horizon::DeckDirectMachineFleetHorizon,
)
    machine_count = length(horizon.machine_types)
    headers = String[
        "time_s",
        "coupling_iterations",
        "coupling_residual_A",
        "cross_machine_impedance_ohm",
        "electrical_cross_impedance_ohm",
        "mechanical_cross_impedance_rad_s_per_Nm",
    ]
    for machine_index in 1:machine_count
        prefix = "machine_$(machine_index)_"
        append!(headers, [
            prefix * "type",
            prefix * "torque_Nm",
            prefix * "mechanical_speed_rad_s",
            prefix * "mechanical_angle_rad",
            prefix * "d_axis_flux_Wb",
            prefix * "q_axis_flux_Wb",
            [prefix * "output_$(output_index)" for output_index in 1:8]...,
        ])
    end
    append!(headers, ["tacs_$(name)_V" for name in horizon.control_output_names])
    open(path, "w") do io
        println(io, join(headers, ','))
        for sample_index in eachindex(horizon.time_s)
            values = Any[
                horizon.time_s[sample_index],
                horizon.coupling_iteration_counts[sample_index],
                horizon.coupling_residuals[sample_index],
                horizon.cross_machine_impedance[sample_index],
                horizon.electrical_cross_impedance[sample_index],
                horizon.mechanical_cross_impedance[sample_index],
            ]
            for machine_index in 1:machine_count
                append!(values, (
                    horizon.machine_types[machine_index],
                    horizon.generated_torque[machine_index, sample_index],
                    horizon.mechanical_speed_rad_s[machine_index, sample_index],
                    horizon.mechanical_angle_rad[machine_index, sample_index],
                    horizon.d_axis_flux[machine_index, sample_index],
                    horizon.q_axis_flux[machine_index, sample_index],
                ))
                append!(values, horizon.output_values[:, machine_index, sample_index])
            end
            append!(values, horizon.control_output_values[:, sample_index])
            println(io, join(values, ','))
        end
    end
    return path
end

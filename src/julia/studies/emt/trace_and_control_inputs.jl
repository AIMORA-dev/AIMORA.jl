
function deck_trace(context::EMTStepContext)
    context.step_index > context.step_count ||
        throw(ArgumentError("fixed-step EMT context has not completed all samples"))
    context.trace_write_index == length(context.recorded_step_indices) + 1 ||
        throw(ArgumentError("fixed-step EMT context did not record all requested samples"))
    return DeckEMTTrace(
        context.source,
        context.dt_s,
        context.t_end_s,
        copy(context.node_map),
        copy(context.node_names),
        copy(context.element_names),
        copy(context.time_s),
        copy(context.voltage_pu),
        copy(context.output_channel_names),
        copy(context.output_node_indices),
        copy(context.output_pu),
        copy(context.node_maximum_values),
        copy(context.node_maximum_times_s),
        copy(context.node_minimum_values),
        copy(context.node_minimum_times_s),
        copy(context.output_maximum_values),
        copy(context.output_maximum_times_s),
        copy(context.output_minimum_values),
        copy(context.output_minimum_times_s),
    )
end

function _deck_synchronous_machine_model_values(
    parsed::DeckParser.DeckParseResult,
    kind::Symbol,
    count::Int,
    machine_index::Int=1,
)
    rows = [
        row for row in DeckParser.deck_synchronous_machine_model_parameter_rows(parsed)
        if row.machine_index == machine_index && row.parameter_kind == kind
    ]
    length(rows) == 1 || throw(ArgumentError(
        "synchronous-machine deck requires exactly one $(kind) parameter row",
    ))
    values = zeros(Float64, count)
    for index in 1:min(count, length(only(rows).positional_values))
        value = only(rows).positional_values[index]
        value isa Missing || (values[index] = Float64(value))
    end
    return values
end

function deck_synchronous_machine_parameters(
    parsed::DeckParser.DeckParseResult;
    machine_index::Int=1,
)
    DeckParser.assert_deck_valid!(parsed)
    machine_index > 0 || throw(ArgumentError("machine_index must be positive"))
    definition = _deck_synchronous_machine_model_values(
        parsed,
        :definition,
        11,
        machine_index,
    )
    model_option = _deck_synchronous_machine_model_values(
        parsed,
        :model_option,
        7,
        machine_index,
    )
    reactance = _deck_synchronous_machine_model_values(
        parsed,
        :electrical_reactance,
        8,
        machine_index,
    )
    time_constant =
        _deck_synchronous_machine_model_values(
            parsed,
            :electrical_time_constant,
            7,
            machine_index,
        )
    numask = round(Int, definition[1])
    mass_rows = sort(
        [
            row for row in DeckParser.deck_synchronous_machine_mass_rows(parsed)
            if row.machine_index == machine_index
        ];
        by = row -> row.mass_index,
    )
    length(mass_rows) == numask || throw(ArgumentError(
        "synchronous-machine mass row count does not match the definition card",
    ))
    [row.mass_index for row in mass_rows] == collect(1:numask) ||
        throw(ArgumentError("synchronous-machine mass indices must be consecutive"))
    masses = SynchronousMachineRotorMass[
        SynchronousMachineRotorMass(
            row.torque_fraction,
            row.inertia,
            row.speed_deviation_damping,
            row.mutual_damping,
            row.shaft_stiffness,
            row.absolute_speed_damping,
        )
        for row in mass_rows
    ]
    fitting_rows = [
        row for row in DeckParser.deck_synchronous_machine_parameter_fitting_rows(parsed)
        if row.machine_index == machine_index
    ]
    fitting_value =
        isempty(fitting_rows) || only(fitting_rows).value isa Missing ? -2.0 :
        Float64(only(fitting_rows).value)
    delta_connected = any(
        row -> row.machine_index == machine_index && row.parameter_kind == :delta_connection,
        DeckParser.deck_synchronous_machine_model_parameter_rows(parsed),
    )
    return SynchronousMachineParameters(
        round(Int, definition[2]),
        round(Int, definition[3]),
        round(Int, definition[4]),
        definition[7],
        definition[8],
        definition[9],
        definition[10],
        definition[11],
        Tuple(model_option),
        reactance[1],
        reactance[2],
        reactance[3],
        reactance[4],
        reactance[5],
        reactance[6],
        reactance[7],
        reactance[8],
        time_constant[1],
        time_constant[2],
        time_constant[3],
        time_constant[4],
        time_constant[5],
        time_constant[6],
        time_constant[7],
        fitting_value,
        delta_connected,
        masses,
    )
end

function _deck_synchronous_machine_runtime_context(
    parsed::DeckParser.DeckParseResult,
    dt_s::Float64,
    t_end_s::Float64;
    saturated_transformer_branch_runtime_enabled::Bool,
    coupled_lumped_sequence_history_enabled::Bool,
    recorded_step_indices=nothing,
)
    transformer_intake =
        saturated_transformer_branch_runtime_enabled ?
        _deck_runtime_saturated_transformer_intake(parsed) :
        nothing
    return transformer_intake === nothing ?
           initialize_step_context(
               parsed;
               dt_s = dt_s,
               t_end_s = t_end_s,
               recorded_step_indices,
           ) :
           saturated_transformer_branch_augmented_step_context(
               parsed,
               transformer_intake;
               dt_s = dt_s,
               t_end_s = t_end_s,
               transformer_branch_shunt_capacitance_rows =
                   _deck_transformer_branch_shunt_capacitance_rows(parsed, nothing),
               include_coupled_lumped_sequence_history =
                   coupled_lumped_sequence_history_enabled,
               current_zero_switching = true,
               recorded_step_indices,
           )
end

function _deck_synchronous_machine_terminal_network_admittance(
    context,
    terminal_nodes::Vector{Int},
)
    Base.@nospecialize context
    admittance = zeros(Float64, context.system.node_count, context.system.node_count)
    rhs = zeros(Float64, context.system.node_count)
    for element in context.system.elements
        Nodal.stamp!(admittance, rhs, element, 0.0, context.dt_s)
    end
    return admittance[terminal_nodes, terminal_nodes]
end

function _deck_synchronous_machine_declared_rotor_reference(
    parsed::DeckParser.DeckParseResult,
    parameters::SynchronousMachineParameters,
    terminal_rows,
    sample,
    machine_index::Int,
)
    sample.time_zero_ground_fault && sample.external_excitation_port_initialization ||
        return nothing
    any(
        row -> row.machine_index == machine_index &&
               row.coupling_kind == :external_field_voltage_input,
        DeckParser.deck_synchronous_machine_control_interface_rows(parsed),
    ) || return nothing
    reference_index = findfirst(
        row -> !(row.peak_terminal_voltage isa Missing),
        terminal_rows,
    )
    reference_index === nothing && return nothing
    reference_row = terminal_rows[reference_index]
    peak_voltage = Float64(reference_row.peak_terminal_voltage)
    base_angle_deg = reference_row.angle_deg isa Missing ? 0.0 :
        Float64(reference_row.angle_deg)
    terminal_phasors = ComplexF64[]
    for row in terminal_rows
        local_peak = row.peak_terminal_voltage isa Missing ? peak_voltage :
            Float64(row.peak_terminal_voltage)
        local_angle = _synchronous_machine_terminal_phase_angle_deg(
            row,
            base_angle_deg,
        )
        push!(terminal_phasors, local_peak * cis(deg2rad(local_angle)))
    end
    phase_rotation = cis(2.0 * pi / 3.0)
    positive_sequence = (
        terminal_phasors[1] +
        phase_rotation * terminal_phasors[2] +
        phase_rotation^2 * terminal_phasors[3]
    ) / 3.0
    abs(positive_sequence) > eps(Float64) || throw(ArgumentError(
        "declared synchronous-machine winding reference must be nonzero",
    ))
    pole_pairs = parameters.pole_count / 2.0
    pole_pairs > 0.0 || throw(ArgumentError(
        "synchronous-machine pole count must define a positive pole-pair count",
    ))
    # Type-59 declares phase-A terminal voltage 30 degrees behind the Park
    # rotor reference used by SOLVUM. Retain that reference when the physical
    # terminal phasor collapses to zero under the time-zero fault topology.
    electrical_reference = angle(positive_sequence) + pi / 6.0
    abs(electrical_reference) <= 32.0 * eps(Float64) &&
        (electrical_reference = 0.0)
    return electrical_reference / pole_pairs
end

function _deck_synchronous_machine_initial_state(
    parsed::DeckParser.DeckParseResult,
    context,
    sample;
    machine_index::Int,
    phase_current_phasors=nothing,
)
    Base.@nospecialize context sample
    terminal_rows = sort(
        [
            row for row in DeckParser.deck_synchronous_machine_terminal_voltage_rows(parsed)
            if row.machine_index == machine_index
        ];
        by = row -> row.phase_index,
    )
    length(terminal_rows) == 3 || throw(ArgumentError(
        "synchronous-machine initialization requires three terminal phases",
    ))
    terminal_nodes = Int[row.terminal_node_value for row in terminal_rows]
    tolerances = [
        row for row in DeckParser.deck_synchronous_machine_tolerance_rows(parsed)
        if row.machine_index == machine_index
    ]
    tolerance_values = isempty(tolerances) ? Float64[] : only(tolerances).values
    initialization_tolerance = get(tolerance_values, 1, 1.0e-9)
    omega_tolerance = get(tolerance_values, 2, 1.0e-4)
    speed_tolerance = get(tolerance_values, 3, 1.0e-5)
    max_iterations = round(Int, get(tolerance_values, 4, 10.0))
    machine_current_phasors = phase_current_phasors === nothing ?
        sample.node_current_phasors[terminal_nodes] :
        ComplexF64.(phase_current_phasors)
    length(machine_current_phasors) == 3 || throw(ArgumentError(
        "synchronous-machine initialization requires three current phasors",
    ))
    parameters = deck_synchronous_machine_parameters(parsed; machine_index)
    declared_rotor_reference = _deck_synchronous_machine_declared_rotor_reference(
        parsed,
        parameters,
        terminal_rows,
        sample,
        machine_index,
    )
    initializer = Base.inferencebarrier(synchronous_machine_initial_state)
    return initializer(
        parameters,
        sample.node_voltage_phasors[terminal_nodes],
        machine_current_phasors;
        time_step_s = context.dt_s,
        frequency_hz = sample.steady_state_frequency_hz,
        damping_control = initialization_tolerance,
        initialization_tolerance,
        omega_tolerance,
        speed_tolerance,
        max_iterations,
        terminal_network_admittance =
            _deck_synchronous_machine_terminal_network_admittance(context, terminal_nodes),
        initial_mechanical_angle_rad = declared_rotor_reference,
    )
end

Base.@noinline function deck_synchronous_machine_initial_state(
    parsed::DeckParser.DeckParseResult;
    machine_index::Int=1,
    time_step_s::Real=deck_fixed_step_horizon(parsed).dt_s,
    dynamic_step_count::Int=deck_fixed_step_horizon(parsed).step_count,
    saturated_transformer_branch_runtime_enabled::Bool=true,
    coupled_lumped_sequence_history_enabled::Bool=true,
)
    horizon = deck_fixed_step_horizon(parsed)
    dt_s = Float64(time_step_s)
    dt_s == horizon.dt_s ||
        throw(ArgumentError("synchronous-machine time step must match the deck time step"))
    1 <= dynamic_step_count <= horizon.step_count ||
        throw(ArgumentError("dynamic_step_count must be within the parsed deck horizon"))
    context_builder = Base.inferencebarrier(_deck_synchronous_machine_runtime_context)
    context = context_builder(
        parsed,
        dt_s,
        dynamic_step_count * dt_s;
        saturated_transformer_branch_runtime_enabled,
        coupled_lumped_sequence_history_enabled,
        recorded_step_indices = [0],
    )
    sample_builder =
        Base.inferencebarrier(_deck_synchronous_machine_network_initial_sample)
    sample = sample_builder(parsed, context)
    sample === nothing &&
        throw(ArgumentError("synchronous-machine terminal initial voltages are missing"))
    return _deck_synchronous_machine_initial_state(
        parsed,
        context,
        sample;
        machine_index,
    )
end

const SYNCHRONOUS_MACHINE_ELECTRICAL_OUTPUT_KINDS = Dict(
    1 => :d_axis_current,
    2 => :q_axis_current,
    3 => :zero_sequence_current,
    4 => :field_current,
    5 => :d_axis_damper_current,
    6 => :generator_field_current,
    7 => :q_axis_damper_current,
    8 => :phase_a_current,
    9 => :phase_b_current,
    10 => :phase_c_current,
    11 => :mechanical_power,
    12 => :field_current_magnitude,
    13 => :field_current_angle,
    14 => :electromagnetic_torque,
    15 => :exciter_torque,
    16 => :field_current_angle,
)

function _deck_synchronous_machine_channel_name(
    machine_index::Int,
    kind::Symbol,
    index::Int,
)
    suffix = index == 0 ? "" : string("_", index)
    return Symbol("synchronous_machine_", machine_index, "_", kind, suffix)
end

function _deck_synchronous_machine_dynamic_output_spec(
    parsed::DeckParser.DeckParseResult;
    machine_index::Union{Nothing,Int}=nothing,
)
    rows = [
        row for row in DeckParser.deck_synchronous_machine_output_request_rows(parsed)
        if machine_index === nothing || row.machine_index == machine_index
    ]
    isempty(rows) && return nothing
    summaries = [
        row for row in DeckParser.deck_synchronous_machine_output_summary_rows(parsed)
        if machine_index === nothing || row.machine_count == machine_index
    ]
    mass_count = maximum((row.mass_count for row in summaries); init = 0)
    output_names = Symbol[]
    electrical_codes = Int[]
    angle_indices = Int[]
    speed_indices = Int[]
    shaft_indices = Int[]
    for row in rows
        if row.group_index == 1
            for code in row.output_codes
                kind = get(SYNCHRONOUS_MACHINE_ELECTRICAL_OUTPUT_KINDS, code, nothing)
                kind === nothing && continue
                push!(
                    output_names,
                    _deck_synchronous_machine_channel_name(row.machine_index, kind, 0),
                )
                push!(electrical_codes, code)
            end
        elseif row.group_index == 2
            for index in row.output_codes
                push!(
                    output_names,
                    _deck_synchronous_machine_channel_name(
                        row.machine_index,
                        :mechanical_angle,
                        index,
                    ),
                )
                push!(angle_indices, index)
            end
        elseif row.group_index == 3
            for index in row.output_codes
                push!(
                    output_names,
                    _deck_synchronous_machine_channel_name(
                        row.machine_index,
                        :mechanical_speed,
                        index,
                    ),
                )
                push!(speed_indices, index)
            end
        elseif row.group_index == 4
            for index in row.output_codes
                push!(
                    output_names,
                    _deck_synchronous_machine_channel_name(
                        row.machine_index,
                        :shaft_torque,
                        index,
                    ),
                )
                push!(shaft_indices, index)
            end
        end
    end
    isempty(output_names) && return nothing
    mechanical_required_count = maximum(
        vcat(angle_indices, speed_indices, (index + 1 for index in shaft_indices)...);
        init = 0,
    )
    numask = max(mass_count, mechanical_required_count, 1)
    cu_count = max(
        maximum((code for code in electrical_codes if 1 <= code <= 7); init = 0),
        7,
    )
    histq_count = 2 * numask
    shp_count = isempty(shaft_indices) ? 0 : 5 * numask + maximum(shaft_indices)
    return (
        output_names = output_names,
        electrical_codes = electrical_codes,
        angle_indices = angle_indices,
        speed_indices = speed_indices,
        shaft_indices = shaft_indices,
        numask = numask,
        cu_count = cu_count,
        histq_count = histq_count,
        shp_count = shp_count,
    )
end

function _deck_synchronous_machine_output_values(
    spec,
    time_s::AbstractVector{<:Real},
    dt_s::Real,
)
    output_count = length(spec.output_names)
    values = Matrix{Float64}(undef, output_count, length(time_s))
    output_state = OVER16OutputReportState(
        vsmout_values = zeros(output_count),
        msmout = 0,
        writer = OVER16OutputWriterControlState(isplot = 0, isprin = 0, istep = 0),
        peaknd_values = [0.0, 0.0, 0.0],
    )
    output_kwargs = (
        cu_values = zeros(Float64, spec.cu_count),
        d6 = 0.0,
        d7 = 0.0,
        d8 = 0.0,
        q3 = 0.0,
        sf4 = 0.0,
        sf5 = 0.0,
        cd = 0.0,
        cexc = 0.0,
        angle_history_indices = spec.angle_indices,
        speed_history_indices = spec.speed_indices,
        shaft_segment_indices = spec.shaft_indices,
        histq_values = zeros(Float64, spec.histq_count),
        shp_values = zeros(Float64, spec.shp_count),
        iu = 0,
        numask = spec.numask,
        n22 = 2 * spec.numask,
        acb = 0.0,
        cz = 1.0,
        omega = 0.0,
        radeg = 1.0,
    )
    effective_dt = Float64(dt_s)
    if effective_dt <= 0.0 && length(time_s) >= 2
        effective_dt = abs(Float64(time_s[2]) - Float64(time_s[1]))
    end
    effective_dt > 0.0 || (effective_dt = 1.0)
    for (sample_index, time_value) in enumerate(time_s)
        update = over16_output_report_update!(
            output_state,
            [0.0],
            Float64(time_value),
            effective_dt,
            sample_index - 1;
            iaverg = 0,
            branch_output_count = 0,
            vsmout_output_codes = spec.electrical_codes,
            vsmout_output_kwargs = output_kwargs,
            iplot = 0,
            m4plot = 0,
            iout = 0,
            kbase = 0,
        )
        length(update.output_volti_values) == output_count + 1 ||
            throw(ArgumentError("synchronous-machine output owner returned an unexpected output count"))
        values[:, sample_index] .= update.output_volti_values[2:end]
    end
    return values
end

function _deck_synchronous_machine_output_values_from_matrix(
    spec,
    time_s::AbstractVector{<:Real},
    output_values,
)
    values = Matrix{Float64}(output_values)
    expected_rows = length(spec.output_names)
    size(values, 1) == expected_rows ||
        throw(ArgumentError("synchronous-machine output values must have one row per parsed output channel"))
    size(values, 2) >= length(time_s) ||
        throw(ArgumentError("synchronous-machine output values must cover every trace sample"))
    return copy(values[:, 1:length(time_s)])
end

function _append_deck_synchronous_machine_outputs(
    trace::DeckEMTTrace,
    parsed::DeckParser.DeckParseResult,
    output_values = nothing,
)
    spec = _deck_synchronous_machine_dynamic_output_spec(parsed)
    spec === nothing && return trace
    values =
        output_values === nothing ?
        _deck_synchronous_machine_output_values(spec, trace.time_s, trace.dt_s) :
        _deck_synchronous_machine_output_values_from_matrix(
            spec,
            trace.time_s,
            output_values,
        )
    appended_extrema = _sampled_trace_extrema(values, trace.time_s)
    return DeckEMTTrace(
        trace.source,
        trace.dt_s,
        trace.t_end_s,
        copy(trace.node_map),
        copy(trace.node_names),
        copy(trace.element_names),
        copy(trace.time_s),
        copy(trace.voltage_pu),
        vcat(copy(trace.output_channel_names), spec.output_names),
        copy(trace.output_node_indices),
        vcat(copy(trace.output_pu), values),
        copy(trace.node_maximum_values),
        copy(trace.node_maximum_times_s),
        copy(trace.node_minimum_values),
        copy(trace.node_minimum_times_s),
        vcat(copy(trace.output_maximum_values), appended_extrema.maximum_values),
        vcat(copy(trace.output_maximum_times_s), appended_extrema.maximum_times_s),
        vcat(copy(trace.output_minimum_values), appended_extrema.minimum_values),
        vcat(copy(trace.output_minimum_times_s), appended_extrema.minimum_times_s),
    )
end

function _synchronous_machine_dynamic_output_values(
    spec,
    result,
    time_s::Float64,
    electrical_coefficients::AbstractVector{<:Real},
    electrical_speed_rad_s::Float64,
)
    preview = over16_synchronous_machine_output_preview(
        spec.electrical_codes;
        cu_values = result.cu_values,
        d6 = result.phase_a_current_residual,
        d7 = result.phase_b_current_residual,
        d8 = result.phase_c_current_residual,
        q3 = result.mechanical_power,
        cd = result.electromagnetic_torque,
        cexc = result.exciter_torque,
        angle_history_indices = spec.angle_indices,
        speed_history_indices = spec.speed_indices,
        shaft_segment_indices = spec.shaft_indices,
        histq_values = result.histq_values,
        shp_values = result.shp_values,
        iu = 0,
        numask = spec.numask,
        n22 = 2 * spec.numask,
        acb = electrical_speed_rad_s * time_s,
        cz = Float64(electrical_coefficients[26]),
        omega = electrical_speed_rad_s,
        radeg = 180.0 / pi,
    )
    preview.output_count == length(spec.output_names) ||
        throw(ArgumentError("synchronous-machine output count does not match the deck request"))
    return preview.output_values
end

function _deck_synchronous_machine_control_runtime(
    parsed::DeckParser.DeckParseResult,
    dt_s::Float64,
    numask::Int,
)
    interface_rows = DeckParser.deck_synchronous_machine_control_interface_rows(parsed)
    isempty(interface_rows) && return nothing
    all(row -> row.machine_index == 1, interface_rows) || throw(ArgumentError(
        "the single-machine horizon cannot consume control rows for another machine",
    ))
    field_rows = [
        row for row in interface_rows
        if row.coupling_kind == :field_voltage_multiplier
    ]
    external_field_rows = [
        row for row in interface_rows
        if row.coupling_kind == :external_field_voltage_input
    ]
    length(field_rows) + length(external_field_rows) <= 1 || throw(ArgumentError(
        "a synchronous machine may have only one field-voltage control signal",
    ))
    total_torque_rows = [
        row for row in interface_rows
        if row.coupling_kind == :total_applied_torque_input
    ]
    length(total_torque_rows) <= 1 || throw(ArgumentError(
        "a synchronous machine may have only one total applied-torque control signal",
    ))
    torque_rows = [
        row for row in interface_rows
        if row.coupling_kind == :mechanical_torque_input
    ]
    torque_signals = Union{Nothing,Symbol}[nothing for _ in 1:numask]
    for row in torque_rows
        1 <= row.variable_index <= numask || throw(ArgumentError(
            "mechanical-torque control mass index $(row.variable_index) is outside 1:$numask",
        ))
        torque_signals[row.variable_index] === nothing || throw(ArgumentError(
            "rotor mass $(row.variable_index) has more than one mechanical-torque control signal",
        ))
        torque_signals[row.variable_index] = row.signal_name
    end
    transfer_rows = [row for row in interface_rows if row.direction == :machine_to_control]
    signal_values = _deck_control_system_source_values(
        parsed,
        Dict{Symbol,Float64}(),
        NamedTuple[],
        0.0,
    )
    for row in DeckParser.deck_control_system_function_rows(parsed)
        get!(signal_values, row.name, 0.0)
    end
    for row in DeckParser.deck_control_system_expression_rows(parsed)
        get!(signal_values, row.name, 0.0)
    end
    for row in DeckParser.deck_control_system_device_rows(parsed)
        get!(signal_values, row.name, 0.0)
    end
    signals = ConstantControlSignal[
        ConstantControlSignal(name, signal_values[name])
        for name in sort!(collect(keys(signal_values)); by = String)
    ]
    output_names = _deck_control_system_output_names(parsed, NamedTuple[])
    state = ControlSystemExecutionState(
        signals,
        _deck_control_system_assignments(parsed, NamedTuple[]);
        functions = _deck_control_system_functions(parsed, NamedTuple[]),
        devices = _deck_control_system_devices(parsed),
        output_names,
        deltat_s = dt_s,
        frequency_hz = _deck_control_system_frequency_hz(parsed),
    )
    for row in interface_rows
        haskey(state.values, row.signal_name) || throw(ArgumentError(
            "unresolved synchronous-machine control signal $(row.signal_name)",
        ))
    end
    expression_rows = DeckParser.deck_control_system_expression_rows(parsed)
    control_expressions = if isempty(expression_rows)
        ControlExpressionRuntime[]
    else
        advance_control_system_state!(state, 0, 0.0; execute_devices = false)
        _deck_control_system_expression_runtimes(parsed, state)
    end
    supplemental = _control_system_supplemental_device_runtime(parsed, state, dt_s)
    _advance_control_system_owner!(
        state,
        supplemental,
        control_expressions,
        0,
        0.0,
    )
    return (
        state = state,
        supplemental = supplemental,
        control_expressions = control_expressions,
        field_signal = isempty(field_rows) ? nothing : only(field_rows).signal_name,
        external_field_voltage_signal =
            isempty(external_field_rows) ? nothing :
            only(external_field_rows).signal_name,
        total_applied_torque_signal =
            isempty(total_torque_rows) ? nothing : only(total_torque_rows).signal_name,
        mechanical_torque_signals = torque_signals,
        transfer_rows = transfer_rows,
        interface_rows = interface_rows,
        interface_state = SynchronousMachineTACSInterfaceState(
            zeros(Float64, length(transfer_rows)),
        ),
    )
end

function _deck_synchronous_machine_fleet_control_runtime(
    parsed::DeckParser.DeckParseResult,
    dt_s::Float64,
    numasks::Vector{Int},
)
    interface_rows = DeckParser.deck_synchronous_machine_control_interface_rows(parsed)
    isempty(interface_rows) && return nothing
    machine_count = length(numasks)
    all(row -> 1 <= row.machine_index <= machine_count, interface_rows) ||
        throw(ArgumentError("synchronous-machine control row has an invalid machine index"))

    signal_values = _deck_control_system_source_values(
        parsed,
        Dict{Symbol,Float64}(),
        NamedTuple[],
        0.0,
    )
    for row in DeckParser.deck_control_system_function_rows(parsed)
        get!(signal_values, row.name, 0.0)
    end
    for row in DeckParser.deck_control_system_expression_rows(parsed)
        get!(signal_values, row.name, 0.0)
    end
    for row in DeckParser.deck_control_system_device_rows(parsed)
        get!(signal_values, row.name, 0.0)
    end
    signals = ConstantControlSignal[
        ConstantControlSignal(name, signal_values[name])
        for name in sort!(collect(keys(signal_values)); by = String)
    ]
    output_names = _deck_control_system_output_names(parsed, NamedTuple[])
    state = ControlSystemExecutionState(
        signals,
        _deck_control_system_assignments(parsed, NamedTuple[]);
        functions = _deck_control_system_functions(parsed, NamedTuple[]),
        devices = _deck_control_system_devices(parsed),
        output_names,
        deltat_s = dt_s,
        frequency_hz = _deck_control_system_frequency_hz(parsed),
    )
    for row in interface_rows
        haskey(state.values, row.signal_name) || throw(ArgumentError(
            "unresolved synchronous-machine control signal $(row.signal_name)",
        ))
    end
    expression_rows = DeckParser.deck_control_system_expression_rows(parsed)
    control_expressions = if isempty(expression_rows)
        ControlExpressionRuntime[]
    else
        advance_control_system_state!(state, 0, 0.0; execute_devices = false)
        _deck_control_system_expression_runtimes(parsed, state)
    end
    supplemental = _control_system_supplemental_device_runtime(parsed, state, dt_s)
    _advance_control_system_owner!(
        state,
        supplemental,
        control_expressions,
        0,
        0.0,
    )

    machine_interfaces = NamedTuple[]
    for machine_index in 1:machine_count
        rows = [row for row in interface_rows if row.machine_index == machine_index]
        field_rows = [
            row for row in rows if row.coupling_kind == :field_voltage_multiplier
        ]
        external_field_rows = [
            row for row in rows
            if row.coupling_kind == :external_field_voltage_input
        ]
        length(field_rows) + length(external_field_rows) <= 1 || throw(ArgumentError(
            "a synchronous machine may have only one field-voltage control signal",
        ))
        total_torque_rows = [
            row for row in rows
            if row.coupling_kind == :total_applied_torque_input
        ]
        length(total_torque_rows) <= 1 || throw(ArgumentError(
            "a synchronous machine may have only one total applied-torque control signal",
        ))
        torque_signals = Union{Nothing,Symbol}[
            nothing for _ in 1:numasks[machine_index]
        ]
        for row in rows
            row.coupling_kind == :mechanical_torque_input || continue
            1 <= row.variable_index <= numasks[machine_index] ||
                throw(ArgumentError("mechanical-torque control mass index is out of range"))
            torque_signals[row.variable_index] === nothing || throw(ArgumentError(
                "a rotor mass has more than one mechanical-torque control signal",
            ))
            torque_signals[row.variable_index] = row.signal_name
        end
        transfer_rows = [row for row in rows if row.direction == :machine_to_control]
        push!(machine_interfaces, (
            state = state,
            field_signal = isempty(field_rows) ? nothing : only(field_rows).signal_name,
            external_field_voltage_signal =
                isempty(external_field_rows) ? nothing :
                only(external_field_rows).signal_name,
            total_applied_torque_signal =
                isempty(total_torque_rows) ? nothing : only(total_torque_rows).signal_name,
            mechanical_torque_signals = torque_signals,
            transfer_rows = transfer_rows,
            interface_rows = rows,
            interface_state = SynchronousMachineTACSInterfaceState(
                zeros(Float64, length(transfer_rows)),
            ),
        ))
    end
    return (
        state = state,
        supplemental = supplemental,
        control_expressions = control_expressions,
        interface_rows = interface_rows,
        machine_interfaces = machine_interfaces,
    )
end

function _synchronous_machine_control_transfer_request(
    row::DeckParser.DeckSynchronousMachineControlInterfaceRow,
    numask::Int,
    nloce::Int,
)
    request = row.variable_index
    if row.interface_code == 73 && nloce == 0
        request <= 17 || throw(ArgumentError(
            "machine-state control request must not exceed electrical variable 17",
        ))
        request != 15 || throw(ArgumentError(
            "machine-state control request 15 requires an exciter mass",
        ))
        return -request
    end
    request <= 3 * numask - 1 || throw(ArgumentError(
        "rotor-mass control request exceeds the machine history/shaft range",
    ))
    return request
end

function _synchronous_machine_control_transfer!(
    runtime,
    result,
    electrical_coefficients::AbstractVector{<:Real},
    numask::Int,
    nloce::Int,
    ;
    publish::Bool=true,
)
    Base.@nospecialize runtime result
    isempty(runtime.transfer_rows) && return Float64[]
    transfer_inputs = (
        cu_values = result.cu_values,
        cv1 = result.d_axis_voltage,
        cv2 = result.q_axis_voltage,
        cv3 = result.zero_sequence_voltage,
        a3 = result.d_axis_current,
        a4 = result.q_axis_current,
        c1 = result.d_axis_rotor_current_1,
        c2 = result.d_axis_rotor_current_2,
        c3 = result.q_axis_rotor_current_1,
        c4 = result.q_axis_rotor_current_2,
        q3 = result.mechanical_power,
        cd = result.electromagnetic_torque,
        cexc = result.exciter_torque,
        ac1 = result.q_axis_internal_voltage,
        ac2 = result.d_axis_internal_voltage,
        acde = electrical_coefficients[29],
        acdf = electrical_coefficients[30],
        elp_i26_19 = electrical_coefficients[19],
        elp_i26_21 = electrical_coefficients[21],
        elp_i75_4 = electrical_coefficients[31],
        histq_values = result.histq_values,
        shp_values = result.shp_values,
        iu = 0,
        numask = numask,
        n22 = 2 * numask,
    )
    generic_indices = findall(runtime.transfer_rows) do row
        row.coupling_kind in (:machine_state_output, :rotor_mass_output)
    end
    request_values = zeros(Float64, length(runtime.transfer_rows))
    if !isempty(generic_indices)
        requests = Int[
            _synchronous_machine_control_transfer_request(
                runtime.transfer_rows[index],
                numask,
                nloce,
            )
            for index in generic_indices
        ]
        generic_state = SynchronousMachineTACSInterfaceState(
            zeros(Float64, length(generic_indices)),
        )
        transfer = synchronous_machine_tacs_transfer_update!(
            generic_state,
            requests;
            transfer_inputs...,
        )
        request_values[generic_indices] .= transfer.request_values
    end
    for (index, row) in enumerate(runtime.transfer_rows)
        if row.coupling_kind == :exciter_voltage_output
            request_values[index] = result.external_field_voltage_control_executed ?
                Float64(result.external_field_voltage_input_pu) :
                Float64(result.field_voltage_multiplier)
        elseif row.coupling_kind == :exciter_current_output
            result.exciter_port_sensor_closed || throw(ArgumentError(
                "synchronous-machine exciter-current transfer requires a closed sensor",
            ))
            request_values[index] = Float64(result.exciter_port_current_a)
        elseif !(row.coupling_kind in (:machine_state_output, :rotor_mass_output))
            throw(ArgumentError(
                "unsupported synchronous-machine control transfer $(row.coupling_kind)",
            ))
        end
    end
    runtime.interface_state.etac_values .= request_values
    runtime.interface_state.lmset = length(request_values)
    runtime.interface_state.transfer_count = length(request_values)
    runtime.interface_state.transfer_mutated = !isempty(request_values)
    if publish
        for (index, row) in enumerate(runtime.transfer_rows)
            runtime.state.values[row.signal_name] = request_values[index]
        end
    end
    return request_values
end

function _advance_synchronous_machine_control!(runtime, step_index::Int, time_s::Float64)
    Base.@nospecialize runtime
    return _advance_control_system_owner!(
        runtime.state,
        runtime.supplemental,
        runtime.control_expressions,
        step_index,
        time_s,
    )
end

function _append_synchronous_machine_control_outputs(
    trace::DeckEMTTrace,
    names::Vector{Symbol},
    values::Matrix{Float64},
)
    isempty(names) && return trace
    size(values) == (length(names), length(trace.time_s)) || throw(ArgumentError(
        "control output matrix must contain one row per name and trace sample",
    ))
    output_names = copy(trace.output_channel_names)
    output_values = copy(trace.output_pu)
    maxima = copy(trace.output_maximum_values)
    maxima_times = copy(trace.output_maximum_times_s)
    minima = copy(trace.output_minimum_values)
    minima_times = copy(trace.output_minimum_times_s)
    extrema = _sampled_trace_extrema(values, trace.time_s)
    for (row_index, name) in enumerate(names)
        output_index = findfirst(==(name), output_names)
        if output_index === nothing
            push!(output_names, name)
            output_values = vcat(output_values, reshape(values[row_index, :], 1, :))
            push!(maxima, extrema.maximum_values[row_index])
            push!(maxima_times, extrema.maximum_times_s[row_index])
            push!(minima, extrema.minimum_values[row_index])
            push!(minima_times, extrema.minimum_times_s[row_index])
        else
            output_values[output_index, :] .= values[row_index, :]
            maxima[output_index] = extrema.maximum_values[row_index]
            maxima_times[output_index] = extrema.maximum_times_s[row_index]
            minima[output_index] = extrema.minimum_values[row_index]
            minima_times[output_index] = extrema.minimum_times_s[row_index]
        end
    end
    return DeckEMTTrace(
        trace.source,
        trace.dt_s,
        trace.t_end_s,
        copy(trace.node_map),
        copy(trace.node_names),
        copy(trace.element_names),
        copy(trace.time_s),
        copy(trace.voltage_pu),
        output_names,
        copy(trace.output_node_indices),
        output_values,
        copy(trace.node_maximum_values),
        copy(trace.node_maximum_times_s),
        copy(trace.node_minimum_values),
        copy(trace.node_minimum_times_s),
        maxima,
        maxima_times,
        minima,
        minima_times,
    )
end

function _deck_output_channel_metadata_entry(row, branch_names, branch_from, branch_to)
    if row.output_kind == :node_voltage
        return DeckOutputChannelMetadata(
            row.name,
            :node_voltage,
            "V",
            String(row.node),
            "",
        )
    elseif row.output_kind in (:branch_voltage, :branch_current)
        branch_index = findfirst(==(row.branch), branch_names)
        branch_index === nothing && throw(ArgumentError(
            "output channel $(row.name) refers to unknown branch $(row.branch)",
        ))
        return DeckOutputChannelMetadata(
            row.name,
            row.output_kind,
            row.output_kind == :branch_voltage ? "V" : "A",
            String(branch_from[branch_index]),
            branch_to[branch_index] == :ground ? "" : String(branch_to[branch_index]),
        )
    end
    throw(ArgumentError("unsupported physical output quantity $(row.output_kind)"))
end

function _record_deck_output_metadata!(
    metadata::Dict{Symbol,DeckOutputChannelMetadata},
    row::DeckOutputChannelMetadata,
)
    if haskey(metadata, row.name)
        metadata[row.name] == row || throw(ArgumentError(
            "conflicting physical output metadata for channel $(row.name)",
        ))
        return metadata
    end
    metadata[row.name] = row
    return metadata
end

function _record_deck_branch_output_metadata!(
    metadata::Dict{Symbol,DeckOutputChannelMetadata},
    names,
    branch_indices,
    quantity::Symbol,
    branch_from,
    branch_to,
)
    for (name, branch_index) in zip(names, branch_indices)
        1 <= branch_index <= length(branch_from) || throw(ArgumentError(
            "output channel $name refers to branch index $branch_index outside the deck model",
        ))
        lower = branch_to[branch_index] == :ground ? "" : String(branch_to[branch_index])
        _record_deck_output_metadata!(
            metadata,
            DeckOutputChannelMetadata(
                name,
                quantity,
                quantity == :branch_voltage ? "V" : "A",
                String(branch_from[branch_index]),
                lower,
            ),
        )
    end
    return metadata
end

function _deck_control_output_unit(name::Symbol)
    name == :TIMEX && return "s"
    name == :ISTEP && return "count"
    name == :DELTAT && return "s"
    name == :FREQHZ && return "Hz"
    name == :OMEGAR && return "rad/s"
    return "1"
end

function _deck_synchronous_machine_output_unit(name::Symbol)
    text = String(name)
    occursin("_current", text) && return "A"
    occursin("_power", text) && return "W"
    occursin("_torque", text) && return "N*m"
    occursin("_angle", text) && return "rad"
    occursin("_speed", text) && return "rad/s"
    return "1"
end

function _record_deck_derived_output_metadata!(
    metadata::Dict{Symbol,DeckOutputChannelMetadata},
    parsed::DeckParser.DeckParseResult,
    trace::DeckEMTTrace,
)
    for row in DeckParser.deck_control_system_switch_coupling_rows(parsed)
        signal_code = min(row.output_code, 3)
        signal_code in (1, 3) || continue
        row.from_node_index === missing && continue
        row.to_node_index === missing && continue
        name = Symbol(String(row.from_node), ".", String(row.to_node))
        _record_deck_output_metadata!(
            metadata,
            DeckOutputChannelMetadata(
                name,
                :branch_current,
                "A",
                String(row.from_node),
                row.to_node == :ground ? "" : String(row.to_node),
            ),
        )
    end
    control_names = Set(_deck_control_system_trace_output_names(parsed))
    for name in trace.output_channel_names
        haskey(metadata, name) && continue
        if name in control_names
            _record_deck_output_metadata!(
                metadata,
                DeckOutputChannelMetadata(
                    name,
                    :tacs,
                    _deck_control_output_unit(name),
                    "TACS",
                    String(name),
                ),
            )
        elseif startswith(String(name), "synchronous_machine_")
            _record_deck_output_metadata!(
                metadata,
                DeckOutputChannelMetadata(
                    name,
                    :synchronous_machine,
                    _deck_synchronous_machine_output_unit(name),
                    "SYNCHRONOUS_MACHINE",
                    String(name),
                ),
            )
        end
    end
    return metadata
end

function _deck_output_channel_metadata_by_name(
    parsed::DeckParser.DeckParseResult,
    trace::DeckEMTTrace,
)
    DeckParser.assert_deck_valid!(parsed)
    branch_names = DeckParser.deck_branch_names(parsed)
    branch_from = DeckParser.deck_branch_from_node_names(parsed)
    branch_to = DeckParser.deck_branch_to_node_names(parsed)
    metadata = Dict{Symbol,DeckOutputChannelMetadata}()
    for row in DeckParser.deck_over15_output_request_rows(parsed)
        _record_deck_output_metadata!(
            metadata,
            _deck_output_channel_metadata_entry(row, branch_names, branch_from, branch_to),
        )
    end
    output_channel_names = DeckParser.deck_over16_output_channel_names(parsed)
    output_node_names = DeckParser.deck_over16_output_node_names(parsed)
    for (name, node) in zip(output_channel_names, output_node_names)
        _record_deck_output_metadata!(
            metadata,
            DeckOutputChannelMetadata(
                name,
                :node_voltage,
                "V",
                String(node),
                "",
            ),
        )
    end
    _record_deck_branch_output_metadata!(
        metadata,
        DeckParser.deck_over16_branch_voltage_output_names(parsed),
        DeckParser.deck_over16_branch_voltage_branch_indices(parsed),
        :branch_voltage,
        branch_from,
        branch_to,
    )
    _record_deck_branch_output_metadata!(
        metadata,
        DeckParser.deck_over16_branch_current_output_names(parsed),
        DeckParser.deck_over16_branch_current_branch_indices(parsed),
        :branch_current,
        branch_from,
        branch_to,
    )

    switch_names = DeckParser.deck_over5_switch_names(parsed)
    switch_from = DeckParser.deck_over5_switch_from_node_names(parsed)
    switch_to = DeckParser.deck_over5_switch_to_node_names(parsed)
    switch_codes = DeckParser.deck_over5_switch_output_codes(parsed)
    for index in eachindex(switch_names, switch_from, switch_to, switch_codes)
        signal_code = _switch_signal_output_code(switch_codes[index])
        lower = switch_to[index] == :ground ? "" : String(switch_to[index])
        if signal_code >= 2
            name = Symbol("switch_voltage_", String(switch_names[index]))
            _record_deck_output_metadata!(
                metadata,
                DeckOutputChannelMetadata(
                    name,
                    :branch_voltage,
                    "V",
                    String(switch_from[index]),
                    lower,
                ),
            )
        end
        if signal_code in (1, 3)
            name = Symbol("switch_current_", String(switch_names[index]))
            _record_deck_output_metadata!(
                metadata,
                DeckOutputChannelMetadata(
                    name,
                    :branch_current,
                    "A",
                    String(switch_from[index]),
                    lower,
                ),
            )
        end
    end
    _record_deck_derived_output_metadata!(metadata, parsed, trace)
    return metadata
end

function deck_report_output_trace(
    parsed::DeckParser.DeckParseResult,
    trace::DeckEMTTrace,
)
    metadata = _deck_output_channel_metadata_by_name(parsed, trace)
    indices = Int[
        index for (index, name) in enumerate(trace.output_channel_names)
        if haskey(metadata, name)
    ]
    isempty(indices) && throw(ArgumentError(
        "trace contains no deck-requested report or binary-plot channels",
    ))
    return DeckEMTTrace(
        trace.source,
        trace.dt_s,
        trace.t_end_s,
        copy(trace.node_map),
        copy(trace.node_names),
        copy(trace.element_names),
        copy(trace.time_s),
        copy(trace.voltage_pu),
        copy(trace.output_channel_names[indices]),
        copy(trace.output_node_indices),
        copy(trace.output_pu[indices, :]),
        copy(trace.node_maximum_values),
        copy(trace.node_maximum_times_s),
        copy(trace.node_minimum_values),
        copy(trace.node_minimum_times_s),
        copy(trace.output_maximum_values[indices]),
        copy(trace.output_maximum_times_s[indices]),
        copy(trace.output_minimum_values[indices]),
        copy(trace.output_minimum_times_s[indices]),
    )
end

function deck_output_channel_metadata(
    parsed::DeckParser.DeckParseResult,
    trace::DeckEMTTrace,
)
    metadata = _deck_output_channel_metadata_by_name(parsed, trace)
    ordered = DeckOutputChannelMetadata[]
    sizehint!(ordered, length(trace.output_channel_names))
    for name in trace.output_channel_names
        haskey(metadata, name) || throw(ArgumentError(
            "trace output channel $name has no typed physical deck metadata",
        ))
        push!(ordered, metadata[name])
    end
    length(ordered) == size(trace.output_pu, 1) || throw(ArgumentError(
        "trace output metadata count must match the output value rows",
    ))
    return ordered
end

function run_deck_synchronous_machine_horizon(
    parsed::DeckParser.DeckParseResult,
    state::SynchronousMachineDynamicState;
    numask::Int,
    nlocg::Int,
    nloce::Int=0,
    time_step_s::Real=deck_fixed_step_horizon(parsed).dt_s,
    dynamic_step_count::Int=deck_fixed_step_horizon(parsed).step_count,
    angle_half_step_inverse::Real,
    speed_tolerance::Real,
    omega_tolerance::Real=Inf,
    speed_floor::Real=1.0e-12,
    max_iterations::Int=1,
    damping_ratio::Real,
    rotor_angle_extrapolation_interval::Real,
    speed_voltage_factor::Real,
    electrical_speed_rad_s::Real,
    electrical_angle_increment::Real,
    saturated_transformer_branch_runtime_enabled::Bool=true,
    coupled_lumped_sequence_history_enabled::Bool=true,
    recorded_step_indices=nothing,
)
    DeckParser.assert_deck_valid!(parsed)
    horizon = deck_fixed_step_horizon(parsed)
    horizon.step_count > 0 ||
        throw(ArgumentError("synchronous-machine horizon requires at least one dynamic step"))
    dt_s = Float64(time_step_s)
    dt_s == horizon.dt_s ||
        throw(ArgumentError("synchronous-machine time step must match the deck time step"))
    1 <= dynamic_step_count <= horizon.step_count || throw(ArgumentError(
        "dynamic_step_count must be within the parsed deck horizon",
    ))
    runtime_t_end_s = dynamic_step_count * dt_s
    spec = _deck_synchronous_machine_dynamic_output_spec(parsed)
    spec === nothing &&
        throw(ArgumentError("synchronous-machine horizon requires parsed output requests"))
    terminal_rows = sort(
        DeckParser.deck_synchronous_machine_terminal_voltage_rows(parsed);
        by = row -> row.phase_index,
    )
    length(terminal_rows) == 3 ||
        throw(ArgumentError("synchronous-machine horizon requires three terminal phases"))
    source_types = unique(row.source_type for row in terminal_rows)
    length(source_types) == 1 || throw(ArgumentError(
        "all terminal phases of a synchronous machine must use one source type",
    ))
    source_type = only(source_types)
    terminal_nodes = Int[row.terminal_node_value for row in terminal_rows]
    all(node -> 1 <= node <= length(parsed.node_map), terminal_nodes) ||
        throw(ArgumentError("synchronous-machine terminal node is outside the deck network"))
    delta_connected = _deck_synchronous_machine_delta_connected(parsed, 1)
    delta_reference_phase_admittance = delta_connected ?
        state.electrical_coefficients[27] - state.electrical_coefficients[28] : nothing

    control_runtime = _deck_synchronous_machine_control_runtime(parsed, dt_s, numask)
    control_output_names = control_runtime === nothing ? Symbol[] :
        copy(control_runtime.state.output_names)
    control_transfer_count = control_runtime === nothing ? 0 :
        length(control_runtime.transfer_rows)

    context = _deck_synchronous_machine_runtime_context(
        parsed,
        dt_s,
        runtime_t_end_s;
        saturated_transformer_branch_runtime_enabled,
        coupled_lumped_sequence_history_enabled,
        recorded_step_indices,
    )
    initial_sample = _initial_voltage_sample_for_context(
        _deck_synchronous_machine_network_initial_sample(parsed, context),
        context.system.node_count,
    )
    initial_sample === nothing &&
        throw(ArgumentError("synchronous-machine terminal initial voltages are missing"))
    _seed_steady_state_network_state!(context, initial_sample)
    initial_output_values = _steady_state_initial_output_values(context)
    _apply_steady_state_initial_sample!(context, initial_sample, initial_output_values)
    context.step_index = 1
    context.t_s = min(context.step_index, context.step_count) * context.dt_s

    call_count = dynamic_step_count + 1
    recorded_steps = _recorded_trace_step_indices(
        dynamic_step_count,
        recorded_step_indices,
    )
    recorded_times = Float64[step * dt_s for step in recorded_steps]
    sample_count = length(recorded_steps)
    terminal_voltages = zeros(Float64, 3, sample_count)
    terminal_currents = zeros(Float64, 3, sample_count)
    terminal_open_circuit_voltages = zeros(Float64, 3, sample_count)
    terminal_impedances = zeros(Float64, 3, 3, sample_count)
    terminal_admittances = zeros(Float64, 3, 3, sample_count)
    machine_outputs = zeros(Float64, length(spec.output_names), sample_count)
    control_outputs = zeros(Float64, length(control_output_names), sample_count)
    machine_control_transfers = zeros(Float64, control_transfer_count, sample_count)
    field_voltage_multipliers = ones(Float64, sample_count)
    external_field_voltage_inputs_pu = zeros(Float64, sample_count)
    mechanical_torque_multipliers = ones(Float64, numask, sample_count)
    total_applied_torques = fill(NaN, sample_count)
    mechanical_histories = zeros(Float64, length(state.equation_state.histq_values), sample_count)
    iterations = zeros(Int, sample_count)
    next_record_position = Ref(1)
    control_system_step_count = control_runtime === nothing ? 0 : 1
    machine_control_input_count = 0
    machine_control_transfer_count = 0

    function advance_machine!(step_index::Int, phase_voltages)
        winding_voltages = _synchronous_machine_winding_voltages(
            phase_voltages,
            delta_connected,
        )
        field_voltage_multiplier =
            control_runtime === nothing || control_runtime.field_signal === nothing ?
            nothing : control_runtime.state.values[control_runtime.field_signal]
        external_field_voltage_input_pu =
            control_runtime === nothing ||
            control_runtime.external_field_voltage_signal === nothing ?
            nothing : control_runtime.state.values[
                control_runtime.external_field_voltage_signal
            ]
        total_applied_torque =
            control_runtime === nothing ||
            control_runtime.total_applied_torque_signal === nothing ?
            nothing : control_runtime.state.values[
                control_runtime.total_applied_torque_signal
            ]
        torque_multipliers = ones(Float64, numask)
        if control_runtime !== nothing
            for mass_index in 1:numask
                signal = control_runtime.mechanical_torque_signals[mass_index]
                signal === nothing ||
                    (torque_multipliers[mass_index] = control_runtime.state.values[signal])
            end
        end
        mass_torque_control_active = control_runtime !== nothing &&
            any(!isnothing, control_runtime.mechanical_torque_signals)
        result = synchronous_machine_dynamic_step!(
            state;
            phase_voltages = winding_voltages,
            numask = numask,
            nlocg = nlocg,
            nloce = nloce,
            delta2 = dt_s / 2.0,
            angle_half_step_inverse = angle_half_step_inverse,
            speed_tolerance = speed_tolerance,
            omega_tolerance = omega_tolerance,
            speed_floor = speed_floor,
            max_iterations = max_iterations,
            machine_sequence_index = 1,
            timestep_index = step_index,
            active_execution_chain = 0,
            damping_ratio = damping_ratio,
            rotor_angle_extrapolation_interval = rotor_angle_extrapolation_interval,
            speed_voltage_factor = speed_voltage_factor,
            electrical_angle_increment = electrical_angle_increment,
            field_voltage_multiplier = field_voltage_multiplier,
            external_field_voltage_input_pu = external_field_voltage_input_pu,
            mechanical_torque_multipliers =
                mass_torque_control_active ? torque_multipliers : nothing,
            total_applied_torque = total_applied_torque,
        )
        if control_runtime !== nothing
            machine_control_input_count +=
                (control_runtime.field_signal === nothing ? 0 : 1) +
                (control_runtime.external_field_voltage_signal === nothing ? 0 : 1) +
                (control_runtime.total_applied_torque_signal === nothing ? 0 : 1) +
                count(!isnothing, control_runtime.mechanical_torque_signals)
        end
        transfer_values = control_runtime === nothing ? Float64[] :
            _synchronous_machine_control_transfer!(
                control_runtime,
                result,
                state.electrical_coefficients,
                numask,
                nloce;
                publish = step_index > 0,
            )
        if control_runtime !== nothing && step_index > 0
            _advance_synchronous_machine_control!(
                control_runtime,
                step_index,
                step_index * dt_s,
            )
            control_system_step_count += 1
            machine_control_transfer_count += length(transfer_values)
        end
        position = next_record_position[]
        if position <= sample_count && recorded_steps[position] == step_index
            terminal_voltages[:, position] .= winding_voltages
            terminal_currents[:, position] .= (
                result.phase_a_current_residual,
                result.phase_b_current_residual,
                result.phase_c_current_residual,
            )
            terminal_admittances[:, :, position] .=
                _deck_synchronous_machine_terminal_admittance(
                    parsed,
                    state,
                    1,
                    delta_reference_phase_admittance,
                )
            machine_outputs[:, position] .=
                _synchronous_machine_dynamic_output_values(
                    spec,
                    result,
                    step_index * dt_s,
                    state.electrical_coefficients,
                    Float64(electrical_speed_rad_s),
                )
            if control_runtime !== nothing
                control_outputs[:, position] .=
                    control_system_output_values(control_runtime.state)
                machine_control_transfers[:, position] .= transfer_values
                field_voltage_multipliers[position] =
                    field_voltage_multiplier === nothing ? 1.0 :
                    Float64(field_voltage_multiplier)
                external_field_voltage_inputs_pu[position] =
                    external_field_voltage_input_pu === nothing ? 0.0 :
                    Float64(external_field_voltage_input_pu)
                mechanical_torque_multipliers[:, position] .= torque_multipliers
                total_applied_torques[position] =
                    total_applied_torque === nothing ? NaN :
                    Float64(total_applied_torque)
            end
            mechanical_histories[:, position] .= state.equation_state.histq_values
            iterations[position] = result.iteration_count
            next_record_position[] = position + 1
        end
        return result
    end

    advance_machine!(0, initial_sample.node_voltage_values[terminal_nodes])
    if first(recorded_steps) == 0
        terminal_open_circuit_voltages[:, 1] .= terminal_voltages[:, 1]
    end
    network_solve_count = 0
    for step_index in 1:dynamic_step_count
        _advance_source_function_network!(context)
        machine_history = _synchronous_machine_terminal_currents(
            state.current_history,
            delta_connected,
        )
        raw_terminal_admittance =
            _deck_synchronous_machine_terminal_admittance(parsed, state, 1)
        terminal_admittance =
            _deck_synchronous_machine_terminal_admittance(
                parsed,
                state,
                1,
                delta_reference_phase_admittance,
            )
        compensation_current = function (open_circuit_voltage, terminal_impedance)
            open_terminal_voltage = open_circuit_voltage[terminal_nodes]
            return (I + terminal_admittance * terminal_impedance) \
                   (-terminal_admittance * open_terminal_voltage - machine_history)
        end
        network_result = solve_step_with_compensated_current_injections!(
            context.system,
            step_index * dt_s,
            dt_s,
            zeros(Float64, context.system.node_count),
            terminal_nodes,
            compensation_current,
            switch_time_s = (step_index - 1) * dt_s,
        )
        result = advance_machine!(
            step_index,
            network_result.voltage[terminal_nodes],
        )
        terminal_current_residuals = _synchronous_machine_terminal_currents(
            Float64[
                result.phase_a_current_residual,
                result.phase_b_current_residual,
                result.phase_c_current_residual,
            ],
            delta_connected,
        )
        # The reported winding current uses the full machine companion, while
        # the network current uses the PAST baseline plus the delta UPDATE increment.
        terminal_current_residuals .+=
            (raw_terminal_admittance - terminal_admittance) *
            network_result.voltage[terminal_nodes]
        maximum(abs, terminal_current_residuals .- network_result.compensation_current) <=
            1.0e-8 + 1.0e-12 * maximum(abs, network_result.compensation_current) ||
            throw(ArgumentError("synchronous-machine current compensation did not converge"))
        position = next_record_position[] - 1
        if 1 <= position <= sample_count && recorded_steps[position] == step_index
            terminal_open_circuit_voltages[:, position] .=
                _synchronous_machine_winding_voltages(
                    network_result.open_circuit_voltage[terminal_nodes],
                    delta_connected,
                )
            terminal_impedances[:, :, position] .=
                network_result.compensation_impedance
        end
        record_step!(context, network_result.voltage)
        network_solve_count += 1
    end
    next_record_position[] == sample_count + 1 || throw(ArgumentError(
        "synchronous-machine horizon did not record every requested step",
    ))
    trace = _append_deck_synchronous_machine_outputs(
        deck_trace(context),
        parsed,
        machine_outputs,
    )
    trace = _append_synchronous_machine_control_outputs(
        trace,
        control_output_names,
        control_outputs,
    )
    complete_machine_control_coupling =
        control_runtime !== nothing &&
        control_system_step_count == dynamic_step_count + 1 &&
        machine_control_input_count ==
            (dynamic_step_count + 1) * count(
                row -> row.direction == :control_to_machine,
                control_runtime.interface_rows,
            ) &&
        machine_control_transfer_count == dynamic_step_count * control_transfer_count
    deferred_effects = Symbol[]
    if control_runtime !== nothing && !complete_machine_control_coupling
        push!(deferred_effects, :incomplete_machine_control_coupling)
    end
    return DeckSynchronousMachineHorizon(
        parsed.source,
        trace,
        recorded_times,
        source_type,
        terminal_voltages,
        terminal_currents,
        terminal_open_circuit_voltages,
        terminal_impedances,
        terminal_admittances,
        machine_outputs,
        iterations,
        state,
        network_solve_count,
        state.call_count,
        3 * call_count,
        length(spec.output_names) * call_count,
        control_output_names,
        control_outputs,
        machine_control_transfers,
        field_voltage_multipliers,
        external_field_voltage_inputs_pu,
        mechanical_torque_multipliers,
        total_applied_torques,
        mechanical_histories,
        control_system_step_count,
        machine_control_input_count,
        machine_control_transfer_count,
        complete_machine_control_coupling,
        deferred_effects,
    )
end

function run_deck_synchronous_machine_horizon(
    parsed::DeckParser.DeckParseResult;
    time_step_s::Real=deck_fixed_step_horizon(parsed).dt_s,
    dynamic_step_count::Int=deck_fixed_step_horizon(parsed).step_count,
    saturated_transformer_branch_runtime_enabled::Bool=true,
    coupled_lumped_sequence_history_enabled::Bool=true,
    recorded_step_indices=nothing,
)
    initialization = deck_synchronous_machine_initial_state(
        parsed;
        time_step_s,
        dynamic_step_count,
        saturated_transformer_branch_runtime_enabled,
        coupled_lumped_sequence_history_enabled,
    )
    return run_deck_synchronous_machine_horizon(
        parsed,
        initialization.state;
        numask = initialization.numask,
        nlocg = initialization.nlocg,
        nloce = initialization.nloce,
        time_step_s,
        dynamic_step_count,
        angle_half_step_inverse = initialization.angle_half_step_inverse,
        speed_tolerance = initialization.speed_tolerance,
        omega_tolerance = initialization.omega_tolerance,
        speed_floor = initialization.speed_floor,
        max_iterations = initialization.max_iterations,
        damping_ratio = initialization.damping_ratio,
        rotor_angle_extrapolation_interval =
            initialization.rotor_angle_extrapolation_interval,
        speed_voltage_factor = initialization.speed_voltage_factor,
        electrical_speed_rad_s = initialization.electrical_speed_rad_s,
        electrical_angle_increment = initialization.electrical_angle_increment,
        saturated_transformer_branch_runtime_enabled,
        coupled_lumped_sequence_history_enabled,
        recorded_step_indices,
    )
end

function run_deck_synchronous_machine_fleet_horizon(
    parsed::DeckParser.DeckParseResult;
    time_step_s::Real=deck_fixed_step_horizon(parsed).dt_s,
    dynamic_step_count::Int=deck_fixed_step_horizon(parsed).step_count,
    saturated_transformer_branch_runtime_enabled::Bool=true,
    coupled_lumped_sequence_history_enabled::Bool=true,
    recorded_step_indices=nothing,
)
    DeckParser.assert_deck_valid!(parsed)
    summaries = DeckParser.deck_synchronous_machine_output_summary_rows(parsed)
    machine_count = maximum((row.machine_count for row in summaries); init = 0)
    machine_count > 1 || throw(ArgumentError(
        "synchronous-machine fleet horizon requires at least two machines",
    ))
    horizon = deck_fixed_step_horizon(parsed)
    dt_s = Float64(time_step_s)
    dt_s == horizon.dt_s || throw(ArgumentError(
        "synchronous-machine fleet time step must match the deck time step",
    ))
    1 <= dynamic_step_count <= horizon.step_count || throw(ArgumentError(
        "dynamic_step_count must be within the parsed deck horizon",
    ))

    terminal_node_indices = zeros(Int, 3, machine_count)
    source_types = zeros(Int, machine_count)
    source_group_indices = zeros(Int, machine_count)
    for machine_index in 1:machine_count
        rows = sort(
            [
                row for row in DeckParser.deck_synchronous_machine_terminal_voltage_rows(parsed)
                if row.machine_index == machine_index
            ];
            by = row -> row.phase_index,
        )
        length(rows) == 3 || throw(ArgumentError(
            "each synchronous machine requires three terminal phases",
        ))
        machine_source_types = unique(row.source_type for row in rows)
        length(machine_source_types) == 1 || throw(ArgumentError(
            "all terminal phases of a synchronous machine must use one source type",
        ))
        source_types[machine_index] = only(machine_source_types)
        machine_source_groups = unique(row.source_group_index for row in rows)
        length(machine_source_groups) == 1 || throw(ArgumentError(
            "all terminal phases of a synchronous machine must use one source group",
        ))
        source_group_indices[machine_index] = only(machine_source_groups)
        terminal_node_indices[:, machine_index] .=
            Int[row.terminal_node_value for row in rows]
    end
    flat_terminal_nodes = vec(terminal_node_indices)
    unique_terminal_nodes = unique(flat_terminal_nodes)
    terminal_node_positions = Dict(
        node => position for (position, node) in enumerate(unique_terminal_nodes)
    )
    terminal_node_incidence = zeros(
        Float64,
        length(flat_terminal_nodes),
        length(unique_terminal_nodes),
    )
    for (row, node) in enumerate(flat_terminal_nodes)
        terminal_node_incidence[row, terminal_node_positions[node]] = 1.0
    end
    delta_connections = Bool[
        _deck_synchronous_machine_delta_connected(parsed, machine_index)
        for machine_index in 1:machine_count
    ]

    context = _deck_synchronous_machine_runtime_context(
        parsed,
        dt_s,
        dynamic_step_count * dt_s;
        saturated_transformer_branch_runtime_enabled,
        coupled_lumped_sequence_history_enabled,
        recorded_step_indices,
    )
    network_initial_sample = _deck_synchronous_machine_network_initial_sample(parsed, context)
    network_initial_sample === nothing && throw(ArgumentError(
        "synchronous-machine fleet terminal initial voltages are missing",
    ))
    initial_sample = _initial_voltage_sample_for_context(
        network_initial_sample,
        context.system.node_count,
    )

    source_group_counts = Dict(
        group => count(==(group), source_group_indices)
        for group in unique(source_group_indices)
    )
    for group in unique(source_group_indices)
        members = findall(==(group), source_group_indices)
        reference_nodes = terminal_node_indices[:, first(members)]
        all(member -> terminal_node_indices[:, member] == reference_nodes, members) ||
            throw(ArgumentError(
                "parallel synchronous-machine parts must inherit one terminal source",
            ))
    end
    initializations = [
        _deck_synchronous_machine_initial_state(
            parsed,
            context,
            network_initial_sample;
            machine_index,
            phase_current_phasors =
                network_initial_sample.node_current_phasors[
                    terminal_node_indices[:, machine_index]
                ] ./ source_group_counts[source_group_indices[machine_index]],
        )
        for machine_index in 1:machine_count
    ]
    states = SynchronousMachineDynamicState[
        initialization.state for initialization in initializations
    ]
    numasks = Int[initialization.numask for initialization in initializations]
    specs = [
        _deck_synchronous_machine_dynamic_output_spec(parsed; machine_index)
        for machine_index in 1:machine_count
    ]
    any(isnothing, specs) && throw(ArgumentError(
        "every synchronous machine requires parsed dynamic output requests",
    ))
    machine_output_names = reduce(
        vcat,
        (spec.output_names for spec in specs);
        init = Symbol[],
    )
    output_offsets = cumsum(vcat(0, [length(spec.output_names) for spec in specs]))

    control_runtime = _deck_synchronous_machine_fleet_control_runtime(
        parsed,
        dt_s,
        numasks,
    )
    control_output_names = control_runtime === nothing ? Symbol[] :
        copy(control_runtime.state.output_names)
    transfer_counts = control_runtime === nothing ? zeros(Int, machine_count) :
        Int[
            length(interface.transfer_rows)
            for interface in control_runtime.machine_interfaces
        ]
    transfer_offsets = cumsum(vcat(0, transfer_counts))
    transfer_count = sum(transfer_counts)

    _seed_steady_state_network_state!(context, initial_sample)
    _apply_steady_state_initial_sample!(
        context,
        initial_sample,
        _steady_state_initial_output_values(context),
    )
    context.step_index = 1
    context.t_s = min(context.step_index, context.step_count) * context.dt_s

    recorded_steps = _recorded_trace_step_indices(
        dynamic_step_count,
        recorded_step_indices,
    )
    recorded_times = Float64[step * dt_s for step in recorded_steps]
    sample_count = length(recorded_steps)
    terminal_voltages = zeros(Float64, 3, machine_count, sample_count)
    terminal_currents = zeros(Float64, 3, machine_count, sample_count)
    terminal_open_circuit_voltages = zeros(Float64, 3, machine_count, sample_count)
    terminal_impedances = zeros(Float64, 3 * machine_count, 3 * machine_count, sample_count)
    terminal_admittances = zeros(Float64, 3 * machine_count, 3 * machine_count, sample_count)
    machine_outputs = zeros(Float64, length(machine_output_names), sample_count)
    mechanical_histories = [
        zeros(Float64, length(state.equation_state.histq_values), sample_count)
        for state in states
    ]
    electrical_coefficients = [
        zeros(Float64, length(state.electrical_coefficients), sample_count)
        for state in states
    ]
    saturation_enabled = Bool[state.saturation_enabled for state in states]
    delta_reference_phase_admittances = Float64[
        state.electrical_coefficients[27] - state.electrical_coefficients[28]
        for state in states
    ]
    d_axis_saturation_regions = zeros(Int, machine_count, sample_count)
    q_axis_saturation_regions = zeros(Int, machine_count, sample_count)
    saturation_refactor_counts = zeros(Int, machine_count, sample_count)
    mechanical_torque_multipliers = [
        ones(Float64, numasks[machine_index], sample_count)
        for machine_index in 1:machine_count
    ]
    total_applied_torques = fill(NaN, machine_count, sample_count)
    field_voltage_multipliers = ones(Float64, machine_count, sample_count)
    external_field_voltage_inputs_pu = zeros(Float64, machine_count, sample_count)
    control_outputs = zeros(Float64, length(control_output_names), sample_count)
    machine_control_transfers = zeros(Float64, transfer_count, sample_count)
    next_record_position = Ref(1)
    control_system_step_count = control_runtime === nothing ? 0 : 1
    machine_control_input_count = 0
    machine_control_transfer_count = 0

    function fleet_terminal_admittance(delta_network_mutation::Bool=true)
        admittance = zeros(Float64, 3 * machine_count, 3 * machine_count)
        for machine_index in 1:machine_count
            block = (3 * (machine_index - 1) + 1):(3 * machine_index)
            admittance[block, block] .= _deck_synchronous_machine_terminal_admittance(
                parsed,
                states[machine_index],
                machine_index,
                delta_network_mutation ?
                    delta_reference_phase_admittances[machine_index] : nothing,
            )
        end
        return admittance
    end

    function advance_fleet!(step_index::Int, phase_voltages::Matrix{Float64})
        results = NamedTuple[]
        transfer_values = zeros(Float64, transfer_count)
        step_field_multipliers = ones(Float64, machine_count)
        step_external_field_voltages_pu = zeros(Float64, machine_count)
        step_torque_multipliers = [ones(Float64, count) for count in numasks]
        step_total_applied_torques = fill(NaN, machine_count)
        for machine_index in 1:machine_count
            initialization = initializations[machine_index]
            interface = control_runtime === nothing ? nothing :
                control_runtime.machine_interfaces[machine_index]
            field_multiplier =
                interface === nothing || interface.field_signal === nothing ?
                nothing : interface.state.values[interface.field_signal]
            field_multiplier === nothing ||
                (step_field_multipliers[machine_index] = field_multiplier)
            external_field_voltage_input_pu =
                interface === nothing ||
                interface.external_field_voltage_signal === nothing ?
                nothing : interface.state.values[
                    interface.external_field_voltage_signal
                ]
            external_field_voltage_input_pu === nothing ||
                (step_external_field_voltages_pu[machine_index] =
                    external_field_voltage_input_pu)
            mass_torque_control_active = interface !== nothing &&
                any(!isnothing, interface.mechanical_torque_signals)
            if interface !== nothing
                total_signal = interface.total_applied_torque_signal
                total_signal === nothing ||
                    (step_total_applied_torques[machine_index] =
                        interface.state.values[total_signal])
                for mass_index in 1:numasks[machine_index]
                    signal = interface.mechanical_torque_signals[mass_index]
                    signal === nothing ||
                        (step_torque_multipliers[machine_index][mass_index] =
                            interface.state.values[signal])
                end
                machine_control_input_count +=
                    (interface.field_signal === nothing ? 0 : 1) +
                    (interface.external_field_voltage_signal === nothing ? 0 : 1) +
                    (interface.total_applied_torque_signal === nothing ? 0 : 1) +
                    count(!isnothing, interface.mechanical_torque_signals)
            end
            winding_voltages = _synchronous_machine_winding_voltages(
                phase_voltages[:, machine_index],
                delta_connections[machine_index],
            )
            result = synchronous_machine_dynamic_step!(
                states[machine_index];
                phase_voltages = winding_voltages,
                numask = initialization.numask,
                nlocg = initialization.nlocg,
                nloce = initialization.nloce,
                delta2 = dt_s / 2.0,
                angle_half_step_inverse = initialization.angle_half_step_inverse,
                speed_tolerance = initialization.speed_tolerance,
                omega_tolerance = initialization.omega_tolerance,
                speed_floor = initialization.speed_floor,
                max_iterations = initialization.max_iterations,
                machine_sequence_index = machine_index,
                timestep_index = step_index,
                active_execution_chain = 0,
                damping_ratio = initialization.damping_ratio,
                rotor_angle_extrapolation_interval =
                    initialization.rotor_angle_extrapolation_interval,
                speed_voltage_factor = initialization.speed_voltage_factor,
                electrical_angle_increment = initialization.electrical_angle_increment,
                field_voltage_multiplier = field_multiplier,
                external_field_voltage_input_pu =
                    external_field_voltage_input_pu,
                mechanical_torque_multipliers =
                    mass_torque_control_active ?
                    step_torque_multipliers[machine_index] : nothing,
                total_applied_torque =
                    isnan(step_total_applied_torques[machine_index]) ? nothing :
                    step_total_applied_torques[machine_index],
            )
            push!(results, result)
            if interface !== nothing
                local_transfer = _synchronous_machine_control_transfer!(
                    interface,
                    result,
                    states[machine_index].electrical_coefficients,
                    initialization.numask,
                    initialization.nloce;
                    publish = step_index > 0,
                )
                local_range =
                    (transfer_offsets[machine_index] + 1):transfer_offsets[machine_index + 1]
                transfer_values[local_range] .= local_transfer
            end
        end
        if control_runtime !== nothing && step_index > 0
            _advance_synchronous_machine_control!(
                control_runtime,
                step_index,
                step_index * dt_s,
            )
            control_system_step_count += 1
            machine_control_transfer_count += transfer_count
        end

        position = next_record_position[]
        if position <= sample_count && recorded_steps[position] == step_index
            terminal_admittances[:, :, position] .= fleet_terminal_admittance()
            for machine_index in 1:machine_count
                result = results[machine_index]
                terminal_voltages[:, machine_index, position] .=
                    _synchronous_machine_winding_voltages(
                        phase_voltages[:, machine_index],
                        delta_connections[machine_index],
                    )
                terminal_currents[:, machine_index, position] .= (
                    result.phase_a_current_residual,
                    result.phase_b_current_residual,
                    result.phase_c_current_residual,
                )
                output_range =
                    (output_offsets[machine_index] + 1):output_offsets[machine_index + 1]
                machine_outputs[output_range, position] .=
                    _synchronous_machine_dynamic_output_values(
                        specs[machine_index],
                        result,
                        step_index * dt_s,
                        states[machine_index].electrical_coefficients,
                        initializations[machine_index].electrical_speed_rad_s,
                    )
                mechanical_histories[machine_index][:, position] .=
                    states[machine_index].equation_state.histq_values
                electrical_coefficients[machine_index][:, position] .=
                    states[machine_index].electrical_coefficients
                d_axis_saturation_regions[machine_index, position] =
                    states[machine_index].d_axis_saturation_region
                q_axis_saturation_regions[machine_index, position] =
                    states[machine_index].q_axis_saturation_region
                saturation_refactor_counts[machine_index, position] =
                    states[machine_index].saturation_refactor_count
                mechanical_torque_multipliers[machine_index][:, position] .=
                    step_torque_multipliers[machine_index]
                total_applied_torques[machine_index, position] =
                    step_total_applied_torques[machine_index]
                field_voltage_multipliers[machine_index, position] =
                    step_field_multipliers[machine_index]
                external_field_voltage_inputs_pu[machine_index, position] =
                    step_external_field_voltages_pu[machine_index]
            end
            if control_runtime !== nothing
                control_outputs[:, position] .=
                    control_system_output_values(control_runtime.state)
                machine_control_transfers[:, position] .= transfer_values
            end
            next_record_position[] = position + 1
        end
        return results
    end

    initial_terminal_voltages = reshape(
        initial_sample.node_voltage_values[flat_terminal_nodes],
        3,
        machine_count,
    )
    advance_fleet!(0, initial_terminal_voltages)
    if first(recorded_steps) == 0
        terminal_open_circuit_voltages[:, :, 1] .= terminal_voltages[:, :, 1]
    end

    network_solve_count = 0
    for step_index in 1:dynamic_step_count
        _advance_source_function_network!(context)
        terminal_admittance = fleet_terminal_admittance()
        raw_terminal_admittance = fleet_terminal_admittance(false)
        machine_history = reduce(
            vcat,
            (
                _synchronous_machine_terminal_currents(
                    states[machine_index].current_history,
                    delta_connections[machine_index],
                )
                for machine_index in 1:machine_count
            ),
        )
        network_terminal_admittance =
            transpose(terminal_node_incidence) * terminal_admittance *
            terminal_node_incidence
        network_machine_history =
            transpose(terminal_node_incidence) * machine_history
        compensation_current = function (open_circuit_voltage, terminal_impedance)
            open_terminal_voltage = open_circuit_voltage[unique_terminal_nodes]
            return (I + network_terminal_admittance * terminal_impedance) \
                   (-network_terminal_admittance * open_terminal_voltage -
                    network_machine_history)
        end
        network_result = solve_step_with_compensated_current_injections!(
            context.system,
            step_index * dt_s,
            dt_s,
            zeros(Float64, context.system.node_count),
            unique_terminal_nodes,
            compensation_current,
            switch_time_s = (step_index - 1) * dt_s,
        )
        results = advance_fleet!(
            step_index,
            reshape(
                terminal_node_incidence * network_result.voltage[unique_terminal_nodes],
                3,
                machine_count,
            ),
        )
        residual_currents = reduce(
            vcat,
            (
                _synchronous_machine_terminal_currents(
                    Float64[
                        results[machine_index].phase_a_current_residual,
                        results[machine_index].phase_b_current_residual,
                        results[machine_index].phase_c_current_residual,
                    ],
                    delta_connections[machine_index],
                )
                for machine_index in 1:machine_count
            ),
        )
        # Reconcile full winding residuals with the delta-adjusted network companion.
        terminal_phase_voltages =
            terminal_node_incidence * network_result.voltage[unique_terminal_nodes]
        residual_currents .+=
            (raw_terminal_admittance - terminal_admittance) * terminal_phase_voltages
        network_residual_currents =
            transpose(terminal_node_incidence) * residual_currents
        maximum(abs, network_residual_currents .- network_result.compensation_current) <=
            1.0e-8 + 1.0e-12 * maximum(abs, network_result.compensation_current) ||
            throw(ArgumentError("synchronous-machine fleet compensation did not converge"))
        position = next_record_position[] - 1
        if 1 <= position <= sample_count && recorded_steps[position] == step_index
            open_terminal_voltages = reshape(
                terminal_node_incidence *
                network_result.open_circuit_voltage[unique_terminal_nodes],
                3,
                machine_count,
            )
            for machine_index in 1:machine_count
                terminal_open_circuit_voltages[:, machine_index, position] .=
                    _synchronous_machine_winding_voltages(
                        open_terminal_voltages[:, machine_index],
                        delta_connections[machine_index],
                    )
            end
            terminal_impedances[:, :, position] .=
                terminal_node_incidence * network_result.compensation_impedance *
                transpose(terminal_node_incidence)
        end
        record_step!(context, network_result.voltage)
        network_solve_count += 1
    end
    next_record_position[] == sample_count + 1 || throw(ArgumentError(
        "synchronous-machine fleet horizon did not record every requested step",
    ))

    trace = _append_synchronous_machine_control_outputs(
        deck_trace(context),
        machine_output_names,
        machine_outputs,
    )
    trace = _append_synchronous_machine_control_outputs(
        trace,
        control_output_names,
        control_outputs,
    )
    complete_machine_network_coupling =
        network_solve_count == dynamic_step_count &&
        all(state -> state.call_count == dynamic_step_count + 1, states)
    complete_machine_control_coupling = control_runtime === nothing || (
        control_system_step_count == dynamic_step_count + 1 &&
        machine_control_input_count ==
            (dynamic_step_count + 1) * count(
                row -> row.direction == :control_to_machine,
                control_runtime.interface_rows,
            ) &&
        machine_control_transfer_count == dynamic_step_count * transfer_count
    )
    deferred_effects = Symbol[]
    complete_machine_network_coupling ||
        push!(deferred_effects, :incomplete_machine_fleet_network_coupling)
    complete_machine_control_coupling ||
        push!(deferred_effects, :incomplete_machine_control_coupling)
    call_count = sum(state.call_count for state in states)
    return DeckSynchronousMachineFleetHorizon(
        parsed.source,
        trace,
        recorded_times,
        source_types,
        terminal_node_indices,
        terminal_voltages,
        terminal_currents,
        terminal_open_circuit_voltages,
        terminal_impedances,
        terminal_admittances,
        machine_output_names,
        machine_outputs,
        mechanical_histories,
        electrical_coefficients,
        saturation_enabled,
        d_axis_saturation_regions,
        q_axis_saturation_regions,
        saturation_refactor_counts,
        states,
        network_solve_count,
        call_count,
        3 * call_count,
        length(machine_output_names) * (dynamic_step_count + 1),
        control_output_names,
        control_outputs,
        machine_control_transfers,
        field_voltage_multipliers,
        external_field_voltage_inputs_pu,
        mechanical_torque_multipliers,
        total_applied_torques,
        control_system_step_count,
        machine_control_input_count,
        machine_control_transfer_count,
        complete_machine_network_coupling,
        complete_machine_control_coupling,
        deferred_effects,
    )
end

function _deck_control_system_frequency_hz(parsed::DeckParser.DeckParseResult)
    for row in parsed.over5a_source_rows
        if isfinite(row.sfreq) && row.sfreq != 0.0
            return abs(row.sfreq) / (2.0 * pi)
        end
    end
    return 0.0
end

function _control_system_signal_index(
    name::Symbol,
    ordinary_slots::Dict{Symbol,Int},
    device_slots::Dict{Symbol,Int},
)
    haskey(device_slots, name) && return device_slots[name]
    haskey(ordinary_slots, name) && return ordinary_slots[name]
    throw(ArgumentError("unresolved control-system signal $name"))
end

function _control_system_optional_signal_index(
    name::Union{Missing,Symbol},
    ordinary_slots::Dict{Symbol,Int},
    device_slots::Dict{Symbol,Int};
    counter_name::Bool=false,
)
    name === missing && return 0
    counter_name && name == :COUNTR && return -9999
    return _control_system_signal_index(name, ordinary_slots, device_slots)
end

function _append_control_system_device_parameters!(
    values::Vector{Float64},
    row::DeckParser.DeckControlSystemDeviceRow,
    dt_s::Float64,
)
    length(row.parameter_values) == 3 ||
        throw(ArgumentError("control-system device parameters require three fixed fields"))
    p1, p2, p3 = row.parameter_values
    base = length(values) + 1
    table_start = 0
    table_end = 0
    transport_pointer = 0
    initialize_transport = false
    initialize_rms = false
    initial_integer = 0
    rms_sample_count = 0
    reference_is_value =
        row.reference_signal !== missing && row.reference_signal == :VALUE

    if row.device_type == 50
        append!(values, (p1, p2, p3, p2 == 0.0 ? 0.5 : p2))
    elseif row.device_type in (51, 52)
        append!(values, (p1, p2, p3 == 0.0 ? 2.0 : p3))
    elseif row.device_type == 53
        maximum_delay = row.control_signal === missing ? p2 : p3
        maximum_delay >= 0.0 ||
            throw(ArgumentError("type-53 maximum delay must be nonnegative"))
        history_count = max(1, ceil(Int, maximum_delay / dt_s - 10.0e-12))
        history_start = base + 4
        delay_output_initial = reference_is_value ? 0.0 : p1
        history_values = zeros(Float64, history_count)
        reference_is_value && (history_values[1] = p1)
        append!(
            values,
            (
                Float64(history_start),
                p2,
                Float64(history_count),
                delay_output_initial,
            ),
        )
        append!(values, history_values)
        transport_pointer = history_start
        initialize_transport = !reference_is_value
    elseif row.device_type == 54
        append!(values, (p1 == 0.0 ? -9999.0 : p1, p2, p3 == 0.0 ? -9999.0 : p3))
    elseif row.device_type == 55
        row.table_complete && !isempty(row.table_input_values) ||
            throw(ArgumentError("type-55 digitizer requires a complete input table"))
        append!(values, (p1, p2, p3))
        table_start = length(values) + 1
        append!(values, row.table_input_values)
        table_end = length(values)
    elseif row.device_type == 56
        row.table_complete && !isempty(row.table_input_values) ||
            throw(ArgumentError("type-56 point nonlinearity requires a complete table"))
        length(row.table_output_values) == length(row.table_input_values) ||
            throw(ArgumentError("type-56 input/output table lengths must match"))
        append!(values, (p1, p2, p3))
        table_start = length(values) + 1
        for index in eachindex(row.table_input_values)
            push!(values, row.table_input_values[index], row.table_output_values[index])
        end
        table_end = length(values) - 1
    elseif row.device_type == 57
        row.table_complete && !isempty(row.table_input_values) ||
            throw(ArgumentError("type-57 time switch requires a complete time table"))
        append!(values, (p1, 0.0, 0.0))
        table_start = length(values) + 1
        append!(values, row.table_input_values)
        table_end = length(values)
    elseif row.device_type == 58
        denominator = p1 == 0.0 ? 1.0 : p1
        dynamic_coefficient = 2.0 / dt_s * p3 / denominator
        static_coefficient = p2 / denominator
        append!(
            values,
            (
                static_coefficient,
                static_coefficient + dynamic_coefficient,
                static_coefficient - dynamic_coefficient,
            ),
        )
    elseif row.device_type == 59
        append!(values, ((p1 == 0.0 ? 1.0 : p1) / dt_s, p2, p3))
    elseif row.device_type == 60
        append!(values, (p1, p2))
    elseif row.device_type == 62
        append!(values, (0.0, p2, p3))
    elseif row.device_type == 66
        scale = p1 * dt_s
        scale > 0.0 || throw(ArgumentError("type-66 RMS period scale must be positive"))
        sample_count = max(1, round(Int, inv(scale)))
        initial_value = reference_is_value ? p2 : 0.0
        append!(values, (scale,))
        append!(values, fill(initial_value, sample_count))
        initial_integer = 0
        rms_sample_count = sample_count
        initialize_rms = !reference_is_value
    elseif row.device_type == 67
        p1 != 0.0 || throw(ArgumentError("type-67 limiter gain must be nonzero"))
        append!(values, (p1, p2 * p1, p3 * p1, 0.0))
    else
        append!(values, (p1, p2, p3))
    end

    return (
        base = base,
        table_start = table_start,
        table_end = table_end,
        transport_pointer = transport_pointer,
        initialize_transport = initialize_transport,
        initialize_rms = initialize_rms,
        initial_integer = initial_integer,
        rms_sample_count = rms_sample_count,
    )
end

function _control_system_supplemental_device_runtime(
    parsed::DeckParser.DeckParseResult,
    state::ControlSystemExecutionState,
    dt_s::Float64,
)
    parsed_rows = DeckParser.deck_control_system_device_rows(parsed)
    isempty(parsed_rows) && return nothing
    all(row -> 50 <= row.device_type <= 67, parsed_rows) || return nothing

    device_names = Symbol[row.name for row in parsed_rows]
    length(unique(device_names)) == length(device_names) ||
        throw(ArgumentError("control-system supplemental device names must be unique"))
    device_name_set = Set(device_names)
    ordinary_name_set = Set{Symbol}(
        name for name in keys(state.values) if !(name in device_name_set)
    )
    for function_row in state.functions
        push!(ordinary_name_set, function_row.output_name)
        for term in function_row.input_terms
            term.name in device_name_set || push!(ordinary_name_set, term.name)
        end
    end
    for row in parsed_rows
        for term in row.input_terms
            term.name in device_name_set || push!(ordinary_name_set, term.name)
        end
        for name in (row.control_signal, row.reference_signal)
            name === missing && continue
            name in (:COUNTR, :VALUE) && continue
            name in device_name_set || push!(ordinary_name_set, name)
        end
    end
    ordinary_names = sort!(collect(ordinary_name_set); by = String)
    for name in ordinary_names
        get!(state.values, name, 0.0)
    end
    ordinary_slots = Dict(name => index for (index, name) in enumerate(ordinary_names))
    ordinary_count = length(ordinary_names)
    device_slots = Dict(
        name => ordinary_count + index for (index, name) in enumerate(device_names)
    )

    parsup_values = Float64[]
    rows = OVER16CSUPDeviceRow[]
    transport_pointers = zeros(Int, length(parsed_rows))
    initialize_transport = falses(length(parsed_rows))
    initialize_rms = falses(length(parsed_rows))
    device_integer_values = zeros(Int, length(parsed_rows))
    for (index, parsed_row) in enumerate(parsed_rows)
        parameter_layout = _append_control_system_device_parameters!(
            parsup_values,
            parsed_row,
            dt_s,
        )
        # TACS1A stores supplemental signal terms in reverse card order.
        input_terms = OVER16CSUPDeviceInputTerm[
            over16_csup_device_input_term(
                term.polarity == 0 ? 0 : _control_system_signal_index(
                    term.name,
                    ordinary_slots,
                    device_slots,
                );
                scale = Float64(term.polarity),
            )
            for term in reverse(parsed_row.input_terms)
        ]
        control_index = parsed_row.device_type == 66 ? parameter_layout.rms_sample_count :
            parsed_row.device_type in (50, 55, 56, 57) ? 0 :
            _control_system_optional_signal_index(
                parsed_row.control_signal,
                ordinary_slots,
                device_slots;
                counter_name = parsed_row.device_type == 58,
            )
        reference_index = if parsed_row.device_type == 53
            parameter_layout.transport_pointer
        elseif parsed_row.device_type in (50, 55, 56, 57, 66)
            0
        else
            _control_system_optional_signal_index(
                parsed_row.reference_signal,
                ordinary_slots,
                device_slots,
            )
        end
        push!(
            rows,
            over16_csup_device_row(
                index,
                parsed_row.device_type,
                parameter_layout.base;
                input_terms = input_terms,
                control_index = control_index,
                reference_index = reference_index,
                table_start_index = parameter_layout.table_start,
                table_end_index = parameter_layout.table_end,
            ),
        )
        transport_pointers[index] = parameter_layout.transport_pointer
        initialize_transport[index] = parameter_layout.initialize_transport
        initialize_rms[index] = parameter_layout.initialize_rms
        device_integer_values[index] = parameter_layout.initial_integer
    end

    xtcs_values = zeros(Float64, ordinary_count + length(rows))
    for (index, name) in enumerate(ordinary_names)
        xtcs_values[index] = state.values[name]
    end
    output_slots = Int[ordinary_count + index for index in eachindex(rows)]
    for index in eachindex(device_names)
        xtcs_values[output_slots[index]] = get(state.values, device_names[index], 0.0)
    end
    runtime = ControlSystemSupplementalDeviceRuntime(
        OVER16CSUPState(
            xtcs_values;
            parsup_values = parsup_values,
            device_integer_values = device_integer_values,
        ),
        rows,
        [index == length(rows) ? 0 : index + 1 for index in eachindex(rows)],
        ordinary_names,
        device_names,
        output_slots,
        transport_pointers,
        initialize_transport,
        initialize_rms,
        false,
        0,
    )
    _initialize_control_system_device_histories!(runtime)
    for index in eachindex(device_names)
        state.values[device_names[index]] = runtime.state.xtcs_values[output_slots[index]]
    end
    return runtime
end

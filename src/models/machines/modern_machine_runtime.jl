export ModernMachineAcceptedState,
       ModernMachineRuntime,
       ModernMachineRuntimeSnapshot,
       ModernMachineEventKind,
       MachineFieldVoltageEvent,
       MachineRotorVoltageEvent,
       MachineMechanicalTorqueEvent,
       MachineVoltageReferenceEvent,
       MachineSpeedReferenceEvent,
       MachineControlEnableEvent,
       ModernMachineEvent,
       modern_machine_runtime,
       modern_machine_runtime_snapshot,
       restore_modern_machine_runtime_snapshot!,
       set_machine_port_inputs!,
       schedule_machine_event!,
       advance_modern_machine!,
       modern_machine_runtime_diagnostics

mutable struct ModernMachineAcceptedState
    flux_wb::Vector{Float64}
    winding_current_a::Vector{Float64}
    terminal_voltage_v::Vector{Float64}
    terminal_current_a::Vector{Float64}
    field_voltage_v::Float64
    rotor_voltage_dq_v::NTuple{2,Float64}
    electromagnetic_torque_nm::Float64
    terminal_power_w::Float64
    field_power_w::Float64
    rotor_port_power_w::Float64
    copper_loss_w::Float64
    magnetic_coenergy_j::Float64
    supplied_electrical_energy_j::Float64
    supplied_mechanical_energy_j::Float64
    dissipated_energy_j::Float64
    maximum_energy_residual_j::Float64
    maximum_energy_quadrature_defect_j::Float64
    maximum_kcl_residual_a::Float64
    maximum_flux_residual_wb::Float64
    maximum_nonlinear_iterations::Int
    accepted_time_s::Float64
    accepted_step_count::Int
    event_count::Int
    rollback_count::Int
end

mutable struct _ModernMachineCandidate
    flux_wb::Vector{Float64}
    winding_current_a::Vector{Float64}
    current_flux_jacobian_per_h::Matrix{Float64}
    terminal_voltage_v::Vector{Float64}
    terminal_current_a::Vector{Float64}
    terminal_jacobian_s::Matrix{Float64}
    terminal_flux_sensitivity_wb_per_v::Matrix{Float64}
    control_state::MachineControlState
    shaft_state::MachineShaftState
    shaft_workspace::_MachineShaftWorkspace
    electromagnetic_workspace::_MachineElectromagneticWorkspace
    field_voltage_v::Float64
    rotor_voltage_dq_v::NTuple{2,Float64}
    mechanical_torque_nm::Float64
    electromagnetic_torque_nm::Float64
    terminal_power_w::Float64
    field_power_w::Float64
    rotor_port_power_w::Float64
    copper_loss_w::Float64
    magnetic_coenergy_j::Float64
    flux_residual_wb::Float64
    nonlinear_iterations::Int
    differential_inductance_margin_h::Float64
    evaluated::Bool
end

@enum ModernMachineEventKind begin
    MachineFieldVoltageEvent
    MachineRotorVoltageEvent
    MachineMechanicalTorqueEvent
    MachineVoltageReferenceEvent
    MachineSpeedReferenceEvent
    MachineControlEnableEvent
end

struct ModernMachineEvent
    id::Symbol
    time_s::Float64
    kind::ModernMachineEventKind
    values::NTuple{2,Float64}
    enabled_value::Bool
    priority::Int

    function ModernMachineEvent(
        id::Symbol,
        time_s::Real,
        kind::ModernMachineEventKind;
        value::Real=0.0,
        values=(Float64(value), 0.0),
        enabled_value::Bool=true,
        priority::Integer=0,
    )
        time = Float64(time_s)
        event_values = Tuple(Float64.(values))
        isempty(String(id)) && throw(ArgumentError("machine event id must not be empty"))
        isfinite(time) && time >= 0.0 || throw(ArgumentError(
            "machine event time must be finite and nonnegative",
        ))
        length(event_values) == 2 && all(isfinite, event_values) || throw(ArgumentError(
            "machine event values must contain two finite values",
        ))
        return new(id, time, kind, event_values, enabled_value, Int(priority))
    end
end

mutable struct ModernMachineRuntime <: AbstractNonlinearCurrentDevice
    preparation::ModernMachinePreparation
    terminal_nodes::NTuple{4,Int}
    accepted_state::ModernMachineAcceptedState
    control_state::MachineControlState
    shaft_state::MachineShaftState
    inputs::MachinePortInputs
    candidate::_ModernMachineCandidate
    events::Vector{ModernMachineEvent}
    next_event_index::Int
    candidate_time_s::Float64
    candidate_step_s::Float64
    companion_method::Symbol
    prepared::Bool
    preparation_count::Int
    trial_evaluation_count::Int
    accepted_evaluation_count::Int
    rejected_trial_count::Int
    initial_stored_energy_j::Float64
end

struct ModernMachineRuntimeSnapshot
    schema::Symbol
    deterministic_signature_sha256::String
    accepted_state::ModernMachineAcceptedState
    control_state::MachineControlState
    shaft_state::MachineShaftState
    inputs::MachinePortInputs
    next_event_index::Int
    preparation_count::Int
    trial_evaluation_count::Int
    accepted_evaluation_count::Int
    rejected_trial_count::Int
end

function _machine_electromagnetic_mass_index(specification::ModernMachineSpecification)
    return _machine_shaft_mass_index(specification, specification.electromagnetic_mass)
end

function _machine_initial_accepted_state(
    preparation::ModernMachinePreparation,
    shaft_state::MachineShaftState,
)
    specification = preparation.specification
    mass_index = _machine_electromagnetic_mass_index(specification)
    electrical_angle = specification.pole_pairs * shaft_state.angle_rad[mass_index]
    terminal_voltage = vcat(collect(specification.initial_phase_voltage_v), 0.0)
    evaluation = machine_electromagnetic_evaluation(
        preparation,
        preparation.initial_flux_wb,
    )
    output_map = _machine_terminal_output_matrix(preparation, electrical_angle)
    terminal_current = output_map * evaluation.current_a
    field_current = preparation.layout.field_index === nothing ? 0.0 :
        evaluation.current_a[preparation.layout.field_index]
    rotor_power = 0.0
    for (index, rotor_d_index) in enumerate(preparation.layout.rotor_d_indices)
        preparation.specification.rotor_branches[index].terminal_exposed || continue
        rotor_power += specification.initial_rotor_voltage_dq_v[1] *
            evaluation.current_a[rotor_d_index]
        rotor_power += specification.initial_rotor_voltage_dq_v[2] *
            evaluation.current_a[preparation.layout.rotor_q_indices[index]]
    end
    return ModernMachineAcceptedState(
        copy(preparation.initial_flux_wb),
        evaluation.current_a,
        terminal_voltage,
        terminal_current,
        specification.initial_field_voltage_v,
        specification.initial_rotor_voltage_dq_v,
        evaluation.electromagnetic_torque_nm,
        dot(terminal_voltage, terminal_current),
        specification.initial_field_voltage_v * field_current,
        rotor_power,
        evaluation.copper_loss_w,
        evaluation.magnetic_coenergy_j,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        abs(sum(terminal_current)),
        0.0,
        0,
        0.0,
        0,
        0,
        0,
    )
end

function _machine_candidate(
    preparation::ModernMachinePreparation,
    control_state::MachineControlState,
    shaft_state::MachineShaftState,
)
    state_count = length(preparation.initial_flux_wb)
    return _ModernMachineCandidate(
        copy(preparation.initial_flux_wb),
        zeros(state_count),
        zeros(state_count, state_count),
        zeros(4),
        zeros(4),
        zeros(4, 4),
        zeros(state_count, 4),
        deepcopy(control_state),
        deepcopy(shaft_state),
        _machine_shaft_workspace(preparation.specification),
        _machine_electromagnetic_workspace(preparation),
        0.0,
        (0.0, 0.0),
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0,
        0.0,
        false,
    )
end

function modern_machine_runtime(
    preparation::ModernMachinePreparation,
    terminal_nodes=(1, 2, 3, 0);
    events=ModernMachineEvent[],
)
    nodes = Tuple(Int.(terminal_nodes))
    length(nodes) == 4 || throw(DimensionMismatch(
        "modern machine runtime requires four terminal nodes",
    ))
    all(>=(0), nodes) && all(>(0), nodes[1:3]) || throw(ArgumentError(
        "modern machine phase nodes must be positive and neutral may be ground",
    ))
    length(unique(nodes)) == 4 || throw(ArgumentError(
        "modern machine terminal nodes must be distinct",
    ))
    ordered_events = sort!(ModernMachineEvent[events...]; by=event -> (
        event.time_s,
        event.priority,
        String(event.id),
    ))
    length(unique(getfield.(ordered_events, :id))) == length(ordered_events) ||
        throw(ArgumentError("modern machine event ids must be unique"))
    control_state = machine_control_state(preparation.specification)
    shaft_state = machine_shaft_state(preparation.specification)
    accepted_state = _machine_initial_accepted_state(preparation, shaft_state)
    initial_stored_energy = accepted_state.magnetic_coenergy_j +
        shaft_state.kinetic_energy_j + shaft_state.elastic_energy_j
    return ModernMachineRuntime(
        preparation,
        nodes,
        accepted_state,
        control_state,
        shaft_state,
        machine_port_inputs(preparation.specification),
        _machine_candidate(preparation, control_state, shaft_state),
        ordered_events,
        1,
        accepted_state.accepted_time_s,
        preparation.specification.settings.timestep_s,
        :trapezoidal,
        false,
        0,
        0,
        0,
        0,
        initial_stored_energy,
    )
end

nonlinear_terminal_nodes(runtime::ModernMachineRuntime) = runtime.terminal_nodes
nonlinear_device_formulation(::ModernMachineRuntime) = PhysicalConstitutiveCurrent
nonlinear_device_provenance(runtime::ModernMachineRuntime) =
    runtime.preparation.specification.provenance

function set_machine_port_inputs!(
    runtime::ModernMachineRuntime;
    field_voltage_v::Real=runtime.inputs.field_voltage_v,
    rotor_voltage_dq_v=runtime.inputs.rotor_voltage_dq_v,
    mechanical_torque_nm::Real=runtime.inputs.mechanical_torque_nm,
    voltage_reference_v::Real=runtime.inputs.voltage_reference_v,
    speed_reference_rad_s::Real=runtime.inputs.speed_reference_rad_s,
    control_enabled::Bool=runtime.inputs.control_enabled,
)
    values = Float64.((
        field_voltage_v,
        rotor_voltage_dq_v...,
        mechanical_torque_nm,
        voltage_reference_v,
        speed_reference_rad_s,
    ))
    all(isfinite, values) || throw(ArgumentError("machine port inputs must be finite"))
    runtime.inputs.field_voltage_v = values[1]
    runtime.inputs.rotor_voltage_dq_v = (values[2], values[3])
    runtime.inputs.mechanical_torque_nm = values[4]
    runtime.inputs.voltage_reference_v = values[5]
    runtime.inputs.speed_reference_rad_s = values[6]
    runtime.inputs.control_enabled = control_enabled
    return runtime
end

function schedule_machine_event!(runtime::ModernMachineRuntime, event::ModernMachineEvent)
    event.time_s + runtime.preparation.specification.settings.event_time_tolerance_s >=
        runtime.accepted_state.accepted_time_s || throw(ArgumentError(
            "machine event cannot be scheduled before accepted time",
        ))
    event.id in getfield.(runtime.events, :id) && throw(ArgumentError(
        "machine event id is already scheduled",
    ))
    push!(runtime.events, event)
    sort!(runtime.events; by=item -> (item.time_s, item.priority, String(item.id)))
    runtime.next_event_index = findfirst(
        item -> item.time_s + runtime.preparation.specification.settings.event_time_tolerance_s >=
            runtime.accepted_state.accepted_time_s,
        runtime.events,
    )
    runtime.next_event_index === nothing && (runtime.next_event_index = length(runtime.events) + 1)
    return runtime
end

function _apply_machine_event!(runtime::ModernMachineRuntime, event::ModernMachineEvent)
    if event.kind === MachineFieldVoltageEvent
        runtime.inputs.field_voltage_v = event.values[1]
    elseif event.kind === MachineRotorVoltageEvent
        runtime.inputs.rotor_voltage_dq_v = event.values
    elseif event.kind === MachineMechanicalTorqueEvent
        runtime.inputs.mechanical_torque_nm = event.values[1]
    elseif event.kind === MachineVoltageReferenceEvent
        runtime.inputs.voltage_reference_v = event.values[1]
    elseif event.kind === MachineSpeedReferenceEvent
        runtime.inputs.speed_reference_rad_s = event.values[1]
    elseif event.kind === MachineControlEnableEvent
        runtime.inputs.control_enabled = event.enabled_value
    else
        error("unsupported modern machine event kind")
    end
    runtime.accepted_state.event_count += 1
    runtime.next_event_index += 1
    return runtime
end

function _apply_due_machine_events!(runtime::ModernMachineRuntime, time_s::Float64)
    tolerance = runtime.preparation.specification.settings.event_time_tolerance_s
    applied = 0
    while true
        event = _machine_next_event(runtime)
        event === nothing && break
        abs(event.time_s - time_s) <= tolerance || break
        _apply_machine_event!(runtime, event)
        applied += 1
    end
    return applied
end

function _machine_next_event(runtime::ModernMachineRuntime)
    runtime.next_event_index <= length(runtime.events) || return nothing
    return runtime.events[runtime.next_event_index]
end

function nonlinear_device_event_surfaces(runtime::ModernMachineRuntime)
    event = _machine_next_event(runtime)
    event === nothing && return ()
    surface = NonlinearDeviceEventSurface(
        event.id,
        (_runtime, time_s) -> time_s - event.time_s,
        (machine_runtime, event_time_s) ->
            _apply_due_machine_events!(machine_runtime, event_time_s);
        direction=:rising,
        priority=event.priority,
        topology_invalidating=false,
        candidate_time=(_runtime -> event.time_s),
    )
    return (surface,)
end

function prepare_nonlinear_device_step!(
    runtime::ModernMachineRuntime,
    time_s::Float64,
    step_s::Float64,
    companion_method::Symbol,
)
    isfinite(time_s) || throw(ArgumentError("machine candidate time must be finite"))
    isfinite(step_s) && step_s > 0.0 || throw(ArgumentError(
        "machine candidate step must be finite and positive",
    ))
    normalized_companion_method = companion_method in (
        :TrapezoidalCompanion,
        :trapezoidal,
    ) ? :trapezoidal : companion_method in (
        :BackwardEulerCompanion,
        :backward_euler,
    ) ? :backward_euler : companion_method
    normalized_companion_method in (:trapezoidal, :backward_euler) ||
        _modern_machine_refusal(
            :unsupported_companion,
            :prepare_step,
            runtime.preparation.specification,
            "modern machine execution requires trapezoidal stepping or the event-localization backward-Euler companion";
            diagnostics=(companion_method=companion_method,),
        )
    expected_step = runtime.preparation.specification.settings.timestep_s
    step_tolerance = max(64.0 * eps(Float64) * expected_step, 1.0e-15)
    step_valid = normalized_companion_method === :trapezoidal ?
        abs(step_s - expected_step) <= step_tolerance :
        step_s <= expected_step + step_tolerance
    step_valid ||
        _modern_machine_refusal(
            :timestep_identity_mismatch,
            :prepare_step,
            runtime.preparation.specification,
            "candidate timestep is incompatible with the prepared fixed-step identity";
            diagnostics=(
                candidate_step_s=step_s,
                prepared_step_s=expected_step,
                companion_method=normalized_companion_method,
            ),
        )
    runtime.candidate_time_s = time_s
    runtime.candidate_step_s = step_s
    runtime.companion_method = normalized_companion_method
    runtime.prepared = true
    runtime.preparation_count += 1
    candidate = runtime.candidate
    _copy_machine_control_state!(candidate.control_state, runtime.control_state)
    _copy_machine_shaft_state!(candidate.shaft_state, runtime.shaft_state)
    mass_index = _machine_electromagnetic_mass_index(runtime.preparation.specification)
    update_machine_control_state!(
        candidate.control_state,
        runtime.preparation.specification,
        runtime.inputs,
        runtime.accepted_state.terminal_voltage_v,
        runtime.shaft_state.speed_rad_s[mass_index],
        time_s,
    )
    candidate.field_voltage_v = runtime.inputs.control_enabled ?
        candidate.control_state.field_voltage_v : runtime.inputs.field_voltage_v
    candidate.mechanical_torque_nm = runtime.inputs.control_enabled ?
        candidate.control_state.mechanical_torque_nm : runtime.inputs.mechanical_torque_nm
    candidate.rotor_voltage_dq_v = runtime.inputs.rotor_voltage_dq_v
    candidate.evaluated = false
    return runtime
end

function _solve_machine_electrical_candidate!(
    runtime::ModernMachineRuntime,
    terminal_voltage_v::AbstractVector{<:Real},
)
    runtime.prepared || throw(ArgumentError(
        "machine nonlinear step must be prepared before trial evaluation",
    ))
    preparation = runtime.preparation
    accepted = runtime.accepted_state
    candidate = runtime.candidate
    specification = preparation.specification
    settings = specification.settings
    mass_index = _machine_electromagnetic_mass_index(specification)
    mechanical_angle = runtime.shaft_state.angle_rad[mass_index]
    mechanical_speed = runtime.shaft_state.speed_rad_s[mass_index]
    electrical_angle = specification.pole_pairs * mechanical_angle
    workspace = candidate.electromagnetic_workspace
    speed_matrix = _machine_electrical_speed_matrix!(
        workspace.speed_matrix_per_s,
        preparation,
        mechanical_speed,
    )
    previous_voltage = _machine_voltage_vector!(
        workspace.previous_voltage_v,
        preparation,
        accepted.terminal_voltage_v,
        electrical_angle,
        accepted.field_voltage_v,
        accepted.rotor_voltage_dq_v,
        workspace,
    )
    candidate_voltage = _machine_voltage_vector!(
        workspace.candidate_voltage_v,
        preparation,
        terminal_voltage_v,
        electrical_angle,
        candidate.field_voltage_v,
        candidate.rotor_voltage_dq_v,
        workspace,
    )
    previous_derivative = workspace.previous_flux_derivative_wb_per_s
    _machine_flux_derivative!(
        previous_derivative,
        candidate.winding_current_a,
        candidate.current_flux_jacobian_per_h,
        preparation,
        accepted.flux_wb,
        previous_voltage,
        speed_matrix,
        workspace,
    )
    step = runtime.candidate_step_s
    candidate_weight = runtime.companion_method === :trapezoidal ? 0.5 : 1.0
    previous_weight = runtime.companion_method === :trapezoidal ? 0.5 : 0.0
    for index in eachindex(candidate.flux_wb)
        candidate.flux_wb[index] =
            accepted.flux_wb[index] + step * previous_derivative[index]
    end
    residual_norm = Inf
    evaluation = _machine_coenergy_current_hessian!(
        candidate.winding_current_a,
        candidate.current_flux_jacobian_per_h,
        workspace.displaced_flux_wb,
        preparation,
        candidate.flux_wb,
    )
    tangent = workspace.tangent
    residual = workspace.residual_wb
    derivative = workspace.candidate_flux_derivative_wb_per_s
    iterations = 0
    for iteration in 1:settings.maximum_nonlinear_iterations
        iterations = iteration
        evaluation = _machine_flux_derivative!(
            derivative,
            candidate.winding_current_a,
            candidate.current_flux_jacobian_per_h,
            preparation,
            candidate.flux_wb,
            candidate_voltage,
            speed_matrix,
            workspace,
        )
        residual_norm = 0.0
        scale = 1.0
        for index in eachindex(residual)
            residual[index] = candidate.flux_wb[index] - accepted.flux_wb[index] -
                step * (
                    previous_weight * previous_derivative[index] +
                    candidate_weight * derivative[index]
                )
            residual_norm = max(residual_norm, abs(residual[index]))
            scale = max(scale, abs(candidate.flux_wb[index]), abs(accepted.flux_wb[index]))
        end
        for column in axes(tangent, 2)
            for row in axes(tangent, 1)
                tangent[row, column] = (row == column ? 1.0 : 0.0) -
                    candidate_weight * step * (
                        -preparation.resistance_ohm[row] *
                            candidate.current_flux_jacobian_per_h[row, column] +
                        speed_matrix[row, column]
                    )
            end
        end
        residual_norm <= settings.nonlinear_tolerance * scale && break
        reciprocal_condition = _factor_machine_tangent!(
            tangent,
            workspace.tangent_pivots,
            workspace.inverse_column,
        )
        reciprocal_condition > 64.0 * eps(Float64) || _modern_machine_refusal(
            :singular_electrical_tangent,
            :trial_evaluation,
            specification,
            "machine trapezoidal electrical tangent is singular";
            diagnostics=(reciprocal_condition=reciprocal_condition,),
        )
        _solve_machine_lu_vector!(tangent, workspace.tangent_pivots, residual)
        for index in eachindex(candidate.flux_wb)
            candidate.flux_wb[index] -= residual[index]
        end
    end
    scale = 1.0
    for value in candidate.flux_wb
        scale = max(scale, abs(value))
    end
    residual_norm <= settings.nonlinear_tolerance * scale || _modern_machine_refusal(
        :electrical_nonconvergence,
        :trial_evaluation,
        specification,
        "machine trapezoidal electrical state did not converge";
        diagnostics=(iterations=iterations, residual_wb=residual_norm),
    )
    # Re-evaluate at the converged flux so the returned current and exact tangent agree.
    evaluation = _machine_flux_derivative!(
        derivative,
        candidate.winding_current_a,
        candidate.current_flux_jacobian_per_h,
        preparation,
        candidate.flux_wb,
        candidate_voltage,
        speed_matrix,
        workspace,
    )
    input_map = _machine_terminal_input_matrix!(
        workspace.terminal_input_map,
        preparation,
        electrical_angle,
        workspace,
    )
    output_map = _machine_terminal_output_matrix!(
        workspace.terminal_output_map,
        preparation,
        electrical_angle,
        workspace,
    )
    copyto!(candidate.terminal_voltage_v, terminal_voltage_v)
    _machine_matrix_vector_product!(
        candidate.terminal_current_a,
        output_map,
        candidate.winding_current_a,
    )
    for column in axes(tangent, 2)
        for row in axes(tangent, 1)
            tangent[row, column] = (row == column ? 1.0 : 0.0) -
                candidate_weight * step * (
                    -preparation.resistance_ohm[row] *
                        candidate.current_flux_jacobian_per_h[row, column] +
                    speed_matrix[row, column]
                )
        end
    end
    reciprocal_condition = _factor_machine_tangent!(
        tangent,
        workspace.tangent_pivots,
        workspace.inverse_column,
    )
    reciprocal_condition > 64.0 * eps(Float64) || _modern_machine_refusal(
        :singular_electrical_tangent,
        :trial_evaluation,
        specification,
        "machine converged terminal sensitivity tangent is singular";
        diagnostics=(reciprocal_condition=reciprocal_condition,),
    )
    for column in axes(candidate.terminal_flux_sensitivity_wb_per_v, 2)
        for row in axes(candidate.terminal_flux_sensitivity_wb_per_v, 1)
            candidate.terminal_flux_sensitivity_wb_per_v[row, column] =
                candidate_weight * step * input_map[row, column]
        end
    end
    _solve_machine_lu_matrix!(
        tangent,
        workspace.tangent_pivots,
        candidate.terminal_flux_sensitivity_wb_per_v,
    )
    _machine_matrix_matrix_product!(
        workspace.terminal_jacobian_product,
        candidate.current_flux_jacobian_per_h,
        candidate.terminal_flux_sensitivity_wb_per_v,
    )
    _machine_matrix_matrix_product!(
        candidate.terminal_jacobian_s,
        output_map,
        workspace.terminal_jacobian_product,
    )
    candidate.electromagnetic_torque_nm = evaluation.electromagnetic_torque_nm
    candidate.magnetic_coenergy_j = evaluation.magnetic_coenergy_j
    candidate.copper_loss_w = evaluation.copper_loss_w
    candidate.terminal_power_w = dot(
        candidate.terminal_voltage_v,
        candidate.terminal_current_a,
    )
    candidate.field_power_w = preparation.layout.field_index === nothing ? 0.0 :
        candidate.field_voltage_v *
            candidate.winding_current_a[preparation.layout.field_index]
    candidate.rotor_port_power_w = 0.0
    for (index, rotor_d_index) in enumerate(preparation.layout.rotor_d_indices)
        specification.rotor_branches[index].terminal_exposed || continue
        candidate.rotor_port_power_w += candidate.rotor_voltage_dq_v[1] *
            candidate.winding_current_a[rotor_d_index]
        candidate.rotor_port_power_w += candidate.rotor_voltage_dq_v[2] *
            candidate.winding_current_a[preparation.layout.rotor_q_indices[index]]
    end
    candidate.flux_residual_wb = residual_norm
    candidate.nonlinear_iterations = iterations
    candidate.differential_inductance_margin_h =
        evaluation.differential_inductance_margin_h
    candidate.evaluated = true
    runtime.trial_evaluation_count += 1
    return candidate
end

function nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    runtime::ModernMachineRuntime,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::Float64,
)
    length(terminal_current_a) >= 4 || throw(DimensionMismatch(
        "machine current workspace must contain four terminals",
    ))
    size(terminal_jacobian_s, 1) >= 4 && size(terminal_jacobian_s, 2) >= 4 ||
        throw(DimensionMismatch("machine Jacobian workspace must be at least 4x4"))
    isfinite(time_s) || throw(ArgumentError("machine trial time must be finite"))
    candidate = _solve_machine_electrical_candidate!(runtime, terminal_voltage_v)
    terminal_current_a[1:4] .= candidate.terminal_current_a
    terminal_jacobian_s[1:4, 1:4] .= candidate.terminal_jacobian_s
    return nothing
end

function _maximum_absolute_difference(left, right)
    axes(left) == axes(right) || return Inf
    maximum_difference = 0.0
    for index in eachindex(left, right)
        maximum_difference = max(
            maximum_difference,
            abs(Float64(left[index]) - Float64(right[index])),
        )
    end
    return maximum_difference
end

function _machine_accept_energy!(runtime::ModernMachineRuntime)
    accepted = runtime.accepted_state
    candidate = runtime.candidate
    old_stored = accepted.magnetic_coenergy_j + runtime.shaft_state.kinetic_energy_j +
        runtime.shaft_state.elastic_energy_j
    step = runtime.candidate_step_s
    controls_active = runtime.inputs.control_enabled &&
        runtime.preparation.specification.controls.enabled
    previous_mechanical_torque_nm = controls_active ?
        runtime.control_state.mechanical_torque_nm : runtime.inputs.mechanical_torque_nm
    _advance_machine_shaft!(
        candidate.shaft_state,
        runtime.preparation.specification,
        Float64(accepted.electromagnetic_torque_nm),
        Float64(candidate.electromagnetic_torque_nm),
        Float64(previous_mechanical_torque_nm),
        Float64(candidate.mechanical_torque_nm),
        step,
        candidate.shaft_workspace,
    )
    new_stored = candidate.magnetic_coenergy_j + candidate.shaft_state.kinetic_energy_j +
        candidate.shaft_state.elastic_energy_j
    electrical_power_old = accepted.terminal_power_w + accepted.field_power_w +
        accepted.rotor_port_power_w
    electrical_power_new = candidate.terminal_power_w + candidate.field_power_w +
        candidate.rotor_port_power_w
    mechanical_power_old = previous_mechanical_torque_nm *
        runtime.shaft_state.speed_rad_s[end]
    mechanical_power_new = candidate.mechanical_torque_nm *
        candidate.shaft_state.speed_rad_s[end]
    supplied_electrical = 0.5 * step * (electrical_power_old + electrical_power_new)
    supplied_mechanical = 0.5 * step * (mechanical_power_old + mechanical_power_new)
    dissipated = 0.5 * step * (accepted.copper_loss_w + candidate.copper_loss_w) +
        step * candidate.shaft_state.damping_loss_w
    energy_quadrature_defect = (new_stored - old_stored) -
        (supplied_electrical + supplied_mechanical - dissipated)
    companion_flux_work = 0.0
    for index in eachindex(accepted.flux_wb, candidate.flux_wb)
        average_current = 0.5 * (
            accepted.winding_current_a[index] + candidate.winding_current_a[index]
        )
        companion_flux_work += average_current *
            (candidate.flux_wb[index] - accepted.flux_wb[index])
    end
    layout = runtime.preparation.layout
    d_index = layout.stator_d_index
    q_index = layout.stator_q_index
    average_flux_d = 0.5 * (accepted.flux_wb[d_index] + candidate.flux_wb[d_index])
    average_flux_q = 0.5 * (accepted.flux_wb[q_index] + candidate.flux_wb[q_index])
    average_current_d = 0.5 * (
        accepted.winding_current_a[d_index] + candidate.winding_current_a[d_index]
    )
    average_current_q = 0.5 * (
        accepted.winding_current_a[q_index] + candidate.winding_current_a[q_index]
    )
    mass_index = _machine_electromagnetic_mass_index(
        runtime.preparation.specification,
    )
    conversion_work = step * runtime.preparation.specification.pole_pairs *
        runtime.shaft_state.speed_rad_s[mass_index] *
        (average_flux_d * average_current_q - average_flux_q * average_current_d)
    companion_net_electrical_work = companion_flux_work + conversion_work
    companion_shaft_dissipation = step * candidate.shaft_state.damping_loss_w
    energy_residual = (new_stored - old_stored) - (
        companion_net_electrical_work + supplied_mechanical -
        companion_shaft_dissipation
    )
    all(isfinite, (
        supplied_electrical,
        supplied_mechanical,
        dissipated,
        energy_residual,
        energy_quadrature_defect,
        companion_net_electrical_work,
        conversion_work,
    )) || _modern_machine_refusal(
        :nonfinite_energy,
        :accept,
        runtime.preparation.specification,
        "machine energy accounting became nonfinite",
    )
    dissipated >= -runtime.preparation.specification.settings.energy_tolerance_j ||
        _modern_machine_refusal(
            :active_dissipation,
            :accept,
            runtime.preparation.specification,
            "machine passive loss accounting became active";
            diagnostics=(dissipated_energy_j=dissipated,),
        )
    accepted.supplied_electrical_energy_j += supplied_electrical
    accepted.supplied_mechanical_energy_j += supplied_mechanical
    accepted.dissipated_energy_j += max(dissipated, 0.0)
    accepted.maximum_energy_residual_j = max(
        accepted.maximum_energy_residual_j,
        abs(energy_residual),
    )
    accepted.maximum_energy_quadrature_defect_j = max(
        accepted.maximum_energy_quadrature_defect_j,
        abs(energy_quadrature_defect),
    )
    return nothing
end

function accept_nonlinear_device_state!(
    runtime::ModernMachineRuntime,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    time_s::Float64,
)
    candidate = runtime.candidate
    if !candidate.evaluated || candidate.terminal_voltage_v != terminal_voltage_v
        candidate = _solve_machine_electrical_candidate!(runtime, terminal_voltage_v)
    end
    _maximum_absolute_difference(candidate.terminal_current_a, terminal_current_a) <=
        max(1.0e-9, 1.0e-9 * maximum(abs, terminal_current_a; init=0.0)) ||
        throw(ArgumentError("accepted machine current differs from its converged trial"))
    _maximum_absolute_difference(candidate.terminal_jacobian_s, terminal_jacobian_s) <=
        max(1.0e-9, 1.0e-9 * maximum(abs, terminal_jacobian_s; init=0.0)) ||
        throw(ArgumentError("accepted machine Jacobian differs from its converged trial"))
    _machine_accept_energy!(runtime)
    accepted = runtime.accepted_state
    accepted.flux_wb .= candidate.flux_wb
    accepted.winding_current_a .= candidate.winding_current_a
    accepted.terminal_voltage_v .= candidate.terminal_voltage_v
    accepted.terminal_current_a .= candidate.terminal_current_a
    accepted.field_voltage_v = candidate.field_voltage_v
    accepted.rotor_voltage_dq_v = candidate.rotor_voltage_dq_v
    accepted.electromagnetic_torque_nm = candidate.electromagnetic_torque_nm
    accepted.terminal_power_w = candidate.terminal_power_w
    accepted.field_power_w = candidate.field_power_w
    accepted.rotor_port_power_w = candidate.rotor_port_power_w
    accepted.copper_loss_w = candidate.copper_loss_w
    accepted.magnetic_coenergy_j = candidate.magnetic_coenergy_j
    accepted.maximum_kcl_residual_a = max(
        accepted.maximum_kcl_residual_a,
        abs(sum(candidate.terminal_current_a)),
    )
    accepted.maximum_flux_residual_wb = max(
        accepted.maximum_flux_residual_wb,
        candidate.flux_residual_wb,
    )
    accepted.maximum_nonlinear_iterations = max(
        accepted.maximum_nonlinear_iterations,
        candidate.nonlinear_iterations,
    )
    accepted.accepted_time_s = time_s
    accepted.accepted_step_count += 1
    _copy_machine_control_state!(runtime.control_state, candidate.control_state)
    _copy_machine_shaft_state!(runtime.shaft_state, candidate.shaft_state)
    runtime.accepted_evaluation_count += 1
    return runtime
end

function accept_nonlinear_device_state!(
    runtime::ModernMachineRuntime,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    time_s::Float64,
)
    candidate = runtime.candidate
    return accept_nonlinear_device_state!(
        runtime,
        terminal_voltage_v,
        terminal_current_a,
        candidate.terminal_jacobian_s,
        time_s,
    )
end

function finish_nonlinear_device_step!(runtime::ModernMachineRuntime)
    runtime.prepared = false
    runtime.candidate.evaluated = false
    return runtime
end

function advance_modern_machine!(
    runtime::ModernMachineRuntime,
    terminal_voltage_v;
    time_s::Real=runtime.accepted_state.accepted_time_s + runtime.candidate_step_s,
)
    time = Float64(time_s)
    step = time - runtime.accepted_state.accepted_time_s
    abs(step - runtime.preparation.specification.settings.timestep_s) <=
        max(64.0 * eps(Float64) * abs(step), 1.0e-15) || throw(ArgumentError(
            "machine standalone advance requires exactly one prepared timestep",
        ))
    event = _machine_next_event(runtime)
    if event !== nothing && event.time_s < time -
        runtime.preparation.specification.settings.event_time_tolerance_s
        throw(ArgumentError(
            "machine standalone advance crossed an event boundary; advance exactly to the event first",
        ))
    end
    if event !== nothing && abs(event.time_s - time) <=
        runtime.preparation.specification.settings.event_time_tolerance_s
        _apply_due_machine_events!(runtime, time)
    end
    terminal_voltage = terminal_voltage_v isa AbstractVector{Float64} ?
        terminal_voltage_v : Float64.(terminal_voltage_v)
    prepare_nonlinear_device_step!(runtime, time, step, :trapezoidal)
    candidate = _solve_machine_electrical_candidate!(runtime, terminal_voltage)
    accept_nonlinear_device_state!(
        runtime,
        candidate.terminal_voltage_v,
        candidate.terminal_current_a,
        candidate.terminal_jacobian_s,
        time,
    )
    finish_nonlinear_device_step!(runtime)
    return runtime
end

function modern_machine_runtime_snapshot(runtime::ModernMachineRuntime)
    runtime.prepared && throw(ArgumentError(
        "machine snapshot is unavailable during an unaccepted trial",
    ))
    return ModernMachineRuntimeSnapshot(
        :aimora_modern_machine_snapshot_v1,
        runtime.preparation.deterministic_signature_sha256,
        deepcopy(runtime.accepted_state),
        deepcopy(runtime.control_state),
        deepcopy(runtime.shaft_state),
        deepcopy(runtime.inputs),
        runtime.next_event_index,
        runtime.preparation_count,
        runtime.trial_evaluation_count,
        runtime.accepted_evaluation_count,
        runtime.rejected_trial_count,
    )
end

function restore_modern_machine_runtime_snapshot!(
    runtime::ModernMachineRuntime,
    snapshot::ModernMachineRuntimeSnapshot,
)
    snapshot.schema === :aimora_modern_machine_snapshot_v1 || throw(ArgumentError(
        "machine snapshot schema is unsupported",
    ))
    snapshot.deterministic_signature_sha256 ==
        runtime.preparation.deterministic_signature_sha256 ||
        _modern_machine_refusal(
            :snapshot_identity_mismatch,
            :restore,
            runtime.preparation.specification,
            "machine snapshot identity differs from the prepared family/topology/settings",
        )
    runtime.accepted_state = deepcopy(snapshot.accepted_state)
    runtime.control_state = deepcopy(snapshot.control_state)
    runtime.shaft_state = deepcopy(snapshot.shaft_state)
    runtime.inputs = deepcopy(snapshot.inputs)
    runtime.next_event_index = snapshot.next_event_index
    runtime.preparation_count = snapshot.preparation_count
    runtime.trial_evaluation_count = snapshot.trial_evaluation_count
    runtime.accepted_evaluation_count = snapshot.accepted_evaluation_count
    runtime.rejected_trial_count = snapshot.rejected_trial_count
    runtime.prepared = false
    runtime.candidate.evaluated = false
    return runtime
end

function modern_machine_runtime_diagnostics(runtime::ModernMachineRuntime)
    accepted = runtime.accepted_state
    specification = runtime.preparation.specification
    return (
        machine=specification.id,
        family=_MODERN_MACHINE_FAMILY_IDS[specification.family],
        operating_mode=_MODERN_MACHINE_MODE_IDS[specification.operating_mode],
        signature_sha256=runtime.preparation.deterministic_signature_sha256,
        accepted_time_s=accepted.accepted_time_s,
        accepted_step_count=accepted.accepted_step_count,
        event_count=accepted.event_count,
        preparation_count=runtime.preparation_count,
        trial_evaluation_count=runtime.trial_evaluation_count,
        accepted_evaluation_count=runtime.accepted_evaluation_count,
        rejected_trial_count=runtime.rejected_trial_count,
        maximum_nonlinear_iterations=accepted.maximum_nonlinear_iterations,
        maximum_flux_residual_wb=accepted.maximum_flux_residual_wb,
        maximum_kcl_residual_a=accepted.maximum_kcl_residual_a,
        maximum_energy_residual_j=accepted.maximum_energy_residual_j,
        maximum_energy_quadrature_defect_j=
            accepted.maximum_energy_quadrature_defect_j,
        maximum_angular_momentum_residual_nms=
            runtime.shaft_state.maximum_angular_momentum_residual_nms,
        magnetic_coenergy_j=accepted.magnetic_coenergy_j,
        kinetic_energy_j=runtime.shaft_state.kinetic_energy_j,
        elastic_energy_j=runtime.shaft_state.elastic_energy_j,
        dissipated_energy_j=accepted.dissipated_energy_j,
        control_sample_count=runtime.control_state.sample_count,
        field_limited=runtime.control_state.field_limited,
        torque_limited=runtime.control_state.torque_limited,
        next_event_index=runtime.next_event_index,
    )
end

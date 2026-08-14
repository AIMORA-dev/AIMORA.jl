function _rebuild_transformer_terminal_fault_admittance!(
    state::TransformerApparatusEventState,
)
    fill!(state.terminal_fault_admittance_s, 0.0)
    for id in sort!(collect(keys(state.active_terminal_faults)); by=String)
        indices, conductance_s = state.active_terminal_faults[id]
        for index in indices
            state.terminal_fault_admittance_s[index, index] += conductance_s
        end
    end
    for id in sort!(collect(keys(state.active_grounding_paths)); by=String)
        indices, conductance_s = state.active_grounding_paths[id]
        for index in indices
            state.terminal_fault_admittance_s[index, index] += conductance_s
        end
    end
    for id in sort!(collect(keys(state.active_winding_faults)); by=String)
        from_index, to_index, conductance_s = state.active_winding_faults[id]
        state.terminal_fault_admittance_s[from_index, from_index] += conductance_s
        state.terminal_fault_admittance_s[to_index, to_index] += conductance_s
        state.terminal_fault_admittance_s[from_index, to_index] -= conductance_s
        state.terminal_fault_admittance_s[to_index, from_index] -= conductance_s
    end
    return state
end

function _validate_transformer_event_command(
    runtime::TransformerApparatusRuntime,
    command::TransformerApparatusEventCommand,
)
    state = runtime.event_state
    command.id in state.applied_event_ids && throw(ArgumentError(
        "transformer event identity $(command.id) has already been applied",
    ))
    any(existing -> existing.id == command.id, state.scheduled_commands) &&
        throw(ArgumentError("transformer event identity $(command.id) is already scheduled"))
    accepted_time = runtime.accepted_state.accepted_time_s
    tolerance = 64.0 * eps(Float64) * max(abs(accepted_time), abs(command.time_s), 1.0)
    command.time_s > accepted_time + tolerance || throw(ArgumentError(
        "transformer event time must make forward progress from accepted state",
    ))
    terminal_count = length(runtime.terminal_nodes)
    if command.kind in (
        TransformerBreakerOpenEvent,
        TransformerBreakerCloseEvent,
        TransformerTerminalFaultApplyEvent,
        TransformerGroundingApplyEvent,
        TransformerWindingFaultApplyEvent,
    )
        all(index -> index <= terminal_count, command.target_indices) ||
            throw(BoundsError(runtime.terminal_nodes, command.target_indices))
    end
    if command.kind in (TransformerTapChangeEvent, TransformerPhaseShiftChangeEvent)
        transform = something(command.terminal_transform)
        size(transform) == (terminal_count, terminal_count) || throw(DimensionMismatch(
            "transformer event terminal transform must match the complete terminal order",
        ))
        condition_number = cond(transform)
        isfinite(condition_number) && condition_number <= 1.0e12 ||
            throw(ArgumentError(
                "transformer event terminal transform must be nonsingular and well-conditioned",
            ))
        command.kind === TransformerPhaseShiftChangeEvent &&
            runtime.preparation.specification.phase_count != 3 && throw(ArgumentError(
                "transformer phase-shift events require a three-phase apparatus",
            ))
    end
    if command.kind === TransformerInternalFaultApplyEvent
        runtime.candidate isa _TransformerNetworkCandidate ||
            _transformer_refusal(
                :unsupported_internal_fault,
                :event_schedule,
                runtime.preparation.specification.id,
                runtime.preparation.specification.tier,
                "internal faults require an explicitly represented grey- or white-box node",
            )
        internal = runtime.candidate.internal_node_indices
        all(index -> index in internal, command.target_indices) ||
            _transformer_refusal(
                :unrepresented_internal_fault,
                :event_schedule,
                runtime.preparation.specification.id,
                runtime.preparation.specification.tier,
                "internal fault targets must be explicitly represented internal nodes";
                diagnostics=(targets=Tuple(command.target_indices),),
            )
    end
    if command.reference_id !== nothing
        reference = something(command.reference_id)
        reference_index = findfirst(
            existing -> existing.id == reference,
            state.scheduled_commands,
        )
        reference_index !== nothing || reference in state.applied_event_ids || throw(ArgumentError(
                "transformer clear event references an unknown apply event $reference",
            ))
        if reference_index !== nothing
            referenced_command = state.scheduled_commands[something(reference_index)]
            expected_apply_kind = Dict(
                TransformerTerminalFaultClearEvent => TransformerTerminalFaultApplyEvent,
                TransformerWindingFaultClearEvent => TransformerWindingFaultApplyEvent,
                TransformerGroundingClearEvent => TransformerGroundingApplyEvent,
                TransformerInternalFaultClearEvent => TransformerInternalFaultApplyEvent,
            )[command.kind]
            referenced_command.kind === expected_apply_kind || throw(ArgumentError(
                "transformer clear event references an incompatible apply-event kind",
            ))
            command.time_s >= referenced_command.time_s || throw(ArgumentError(
                "transformer clear event cannot precede its apply event",
            ))
        end
    end
    return command
end

"""Schedule one future apparatus command without mutating physical state."""
function queue_transformer_apparatus_event!(
    runtime::TransformerApparatusRuntime,
    command::TransformerApparatusEventCommand,
)
    runtime.prepared && throw(ArgumentError(
        "transformer events cannot be scheduled during an active trial step",
    ))
    _validate_transformer_event_command(runtime, command)
    push!(runtime.event_state.scheduled_commands, command)
    sort!(runtime.event_state.scheduled_commands; by=row -> (
        row.time_s,
        row.priority,
        String(row.id),
    ))
    return runtime
end

function _transformer_command_is_applied(
    runtime::TransformerApparatusRuntime,
    command::TransformerApparatusEventCommand,
)
    return command.id in runtime.event_state.applied_event_ids
end

function _transformer_event_value(
    runtime::TransformerApparatusRuntime,
    command::TransformerApparatusEventCommand,
    _time_s::Real,
)
    _transformer_command_is_applied(runtime, command) && return nothing
    return nothing
end

function _transformer_event_candidate_time(
    runtime::TransformerApparatusRuntime,
    command::TransformerApparatusEventCommand,
)
    _transformer_command_is_applied(runtime, command) && return nothing
    return command.time_s
end

function _remove_transformer_active_event!(dictionary, reference_id, label)
    haskey(dictionary, reference_id) || throw(ArgumentError(
        "transformer $label clear event references an inactive apply event $reference_id",
    ))
    delete!(dictionary, reference_id)
    return dictionary
end

function _apply_transformer_apparatus_command!(
    runtime::TransformerApparatusRuntime,
    command::TransformerApparatusEventCommand,
    time_s::Real,
)
    runtime.prepared && throw(ArgumentError(
        "transformer event transition requires an accepted trial boundary",
    ))
    time = Float64(time_s)
    tolerance = 64.0 * eps(Float64) * max(abs(time), abs(command.time_s), 1.0)
    abs(time - command.time_s) <= tolerance || throw(ArgumentError(
        "transformer event transition time does not match its scheduled boundary",
    ))
    _transformer_command_is_applied(runtime, command) && throw(ArgumentError(
        "transformer event $(command.id) cannot be applied twice",
    ))
    state = runtime.event_state
    kind = command.kind
    if kind === TransformerEnergizeEvent
        fill!(state.terminal_closed, true)
    elseif kind === TransformerDeenergizeEvent
        fill!(state.terminal_closed, false)
    elseif kind === TransformerBreakerOpenEvent
        state.terminal_closed[command.target_indices] .= false
    elseif kind === TransformerBreakerCloseEvent
        state.terminal_closed[command.target_indices] .= true
    elseif kind === TransformerTerminalFaultApplyEvent
        state.active_terminal_faults[command.id] =
            (copy(command.target_indices), command.conductance_s)
    elseif kind === TransformerTerminalFaultClearEvent
        _remove_transformer_active_event!(
            state.active_terminal_faults,
            something(command.reference_id),
            "terminal fault",
        )
    elseif kind === TransformerWindingFaultApplyEvent
        state.active_winding_faults[command.id] = (
            command.target_indices[1],
            command.target_indices[2],
            command.conductance_s,
        )
    elseif kind === TransformerWindingFaultClearEvent
        _remove_transformer_active_event!(
            state.active_winding_faults,
            something(command.reference_id),
            "winding fault",
        )
    elseif kind === TransformerGroundingApplyEvent
        state.active_grounding_paths[command.id] =
            (copy(command.target_indices), command.conductance_s)
    elseif kind === TransformerGroundingClearEvent
        _remove_transformer_active_event!(
            state.active_grounding_paths,
            something(command.reference_id),
            "grounding",
        )
    elseif kind in (TransformerTapChangeEvent, TransformerPhaseShiftChangeEvent)
        state.terminal_transform .= something(command.terminal_transform)
        kind === TransformerTapChangeEvent ?
            (state.tap_change_count += 1) :
            (state.phase_shift_change_count += 1)
    elseif kind === TransformerInternalFaultApplyEvent
        state.active_internal_faults[command.id] = (
            command.target_indices[1],
            command.target_indices[2],
            command.conductance_s,
        )
    elseif kind === TransformerInternalFaultClearEvent
        _remove_transformer_active_event!(
            state.active_internal_faults,
            something(command.reference_id),
            "internal fault",
        )
    else
        error("unreachable transformer event kind")
    end
    _rebuild_transformer_terminal_fault_admittance!(state)
    state.previous_fault_current_a .=
        state.terminal_fault_admittance_s * state.previous_network_voltage_v
    push!(state.applied_event_ids, command.id)
    push!(state.occurrences, TransformerApparatusEventOccurrence(
        command.id,
        command.kind,
        time,
        command.priority,
        command.topology_invalidating,
    ))
    command.topology_invalidating && (state.topology_transition_count += 1)
    return runtime
end

function nonlinear_device_event_surfaces(runtime::TransformerApparatusRuntime)
    surfaces = NonlinearDeviceEventSurface[]
    accepted_time = runtime.accepted_state.accepted_time_s
    for command in runtime.event_state.scheduled_commands
        _transformer_command_is_applied(runtime, command) && continue
        tolerance = 64.0 * eps(Float64) * max(
            abs(accepted_time),
            abs(command.time_s),
            1.0,
        )
        command.time_s >= accepted_time - tolerance || throw(ArgumentError(
            "transformer event $(command.id) was missed before the accepted time",
        ))
        push!(surfaces, NonlinearDeviceEventSurface(
            command.id,
            (device, time_s) -> _transformer_event_value(device, command, time_s),
            (device, time_s) ->
                _apply_transformer_apparatus_command!(device, command, time_s);
            direction=:rising,
            priority=command.priority,
            topology_invalidating=command.topology_invalidating,
            candidate_time=device -> _transformer_event_candidate_time(device, command),
        ))
    end
    return Tuple(surfaces)
end

function _transformer_open_terminal_correction(jacobian, residual)
    isempty(residual) && return Float64[]
    decomposition = svd(jacobian)
    largest = maximum(decomposition.S; init=0.0)
    threshold = max(size(jacobian)...) * eps(Float64) * max(largest, 1.0)
    inverse_values = map(value -> value > threshold ? inv(value) : 0.0, decomposition.S)
    correction = -(decomposition.V * Diagonal(inverse_values) *
        transpose(decomposition.U) * residual)
    remaining = jacobian * correction + residual
    maximum(abs, remaining; init=0.0) <=
        1.0e-9 * max(maximum(abs, residual; init=0.0), 1.0) ||
        throw(ArgumentError(
            "transformer open-circuit terminal constraint is inconsistent",
        ))
    return correction
end

function _transformer_event_trial!(
    runtime::TransformerApparatusRuntime,
    network_terminal_voltage_v::AbstractVector{Float64},
)
    state = runtime.event_state
    apparatus_voltage = runtime.candidate_apparatus_terminal_voltage_v
    if all(state.terminal_closed)
        mul!(apparatus_voltage, state.terminal_transform, network_terminal_voltage_v)
        current, jacobian = _transformer_raw_trial!(runtime, apparatus_voltage)
        network_current = runtime.candidate_network_terminal_current_a
        network_jacobian = runtime.candidate_network_terminal_jacobian_s
        mul!(
            network_current,
            state.terminal_fault_admittance_s,
            network_terminal_voltage_v,
        )
        mul!(network_current, transpose(state.terminal_transform), current, 1.0, 1.0)
        mul!(network_jacobian, transpose(state.terminal_transform), jacobian)
        mul!(jacobian, network_jacobian, state.terminal_transform)
        network_jacobian .= jacobian .+ state.terminal_fault_admittance_s
        return network_current, network_jacobian
    end
    closed = findall(state.terminal_closed)
    open = findall(.!state.terminal_closed)
    apparatus_voltage .= runtime.accepted_state.terminal_voltage_v
    if !isempty(closed)
        apparatus_voltage[closed] .=
            state.terminal_transform[closed, :] * network_terminal_voltage_v
    end
    current = nothing
    jacobian = nothing
    if !isempty(open)
        converged = false
        for _iteration in 1:runtime.preparation.specification.settings.maximum_local_nonlinear_iterations
            current, jacobian = _transformer_raw_trial!(runtime, apparatus_voltage)
            residual = current[open]
            scale = max(maximum(abs, current; init=0.0), 1.0)
            if maximum(abs, residual; init=0.0) <=
               runtime.preparation.specification.settings.nonlinear_residual_relative_tolerance *
               scale
                converged = true
                break
            end
            apparatus_voltage[open] .+=
                _transformer_open_terminal_correction(jacobian[open, open], residual)
        end
        converged || _transformer_refusal(
            :open_circuit_nonconvergence,
            :trial,
            runtime.preparation.specification.id,
            runtime.preparation.specification.tier,
            "transformer open terminal voltage did not satisfy zero terminal current",
        )
    end
    current, jacobian = _transformer_raw_trial!(runtime, apparatus_voltage)
    network_current = runtime.candidate_network_terminal_current_a
    network_jacobian = runtime.candidate_network_terminal_jacobian_s
    network_current .= state.terminal_fault_admittance_s * network_terminal_voltage_v
    network_jacobian .= state.terminal_fault_admittance_s
    if !isempty(closed)
        closed_transform = state.terminal_transform[closed, :]
        reduced_jacobian = if isempty(open)
            jacobian[closed, closed]
        else
            jacobian[closed, closed] .+
                jacobian[closed, open] *
                (_transformer_open_terminal_correction(
                    jacobian[open, open],
                    jacobian[open, closed],
                ))
        end
        network_current .+= transpose(closed_transform) * current[closed]
        network_jacobian .+=
            transpose(closed_transform) * reduced_jacobian * closed_transform
    end
    return network_current, network_jacobian
end

function _accept_transformer_event_state!(
    runtime::TransformerApparatusRuntime,
    network_voltage,
    network_current,
)
    state = runtime.event_state
    fault_current = runtime.candidate_apparatus_terminal_voltage_v
    mul!(fault_current, state.terminal_fault_admittance_s, network_voltage)
    external_fault_energy_increment = runtime.candidate_step_s *
        _transformer_endpoint_inner_product(
            state.previous_network_voltage_v,
            network_voltage,
            state.previous_fault_current_a,
            fault_current,
            runtime.companion_method === :backward_euler,
        )
    state.event_energy_j += external_fault_energy_increment
    state.external_fault_energy_j += external_fault_energy_increment
    state.event_energy_j >= -512.0 * eps(Float64) || throw(ArgumentError(
        "transformer passive fault/grounding energy became negative",
    ))
    state.previous_network_voltage_v .= network_voltage
    state.previous_fault_current_a .= fault_current
    runtime.candidate_network_terminal_current_a .= network_current
    accepted_state = runtime.accepted_state
    balance_residual = runtime.initial_stored_energy_j +
        accepted_state.supplied_energy_j + state.external_fault_energy_j -
        _transformer_total_stored_energy(accepted_state) -
        _transformer_total_dissipated_energy(accepted_state) - state.event_energy_j
    tolerance = runtime.preparation.specification.settings.energy_absolute_tolerance_j +
        1024.0 * eps(Float64) * max(
            abs(accepted_state.supplied_energy_j),
            abs(runtime.initial_stored_energy_j),
            abs(state.event_energy_j),
            1.0,
        )
    unexplained_residual = balance_residual - state.numerical_dissipation_energy_j
    energy_admissible = runtime.companion_method === :trapezoidal ?
        abs(unexplained_residual) <= tolerance : unexplained_residual >= -tolerance
    energy_admissible || throw(ArgumentError(
        "transformer cumulative energy balance became active: " *
        "unexplained_residual_j=$(unexplained_residual), " *
        "balance_residual_j=$(balance_residual), " *
        "numerical_dissipation_j=$(state.numerical_dissipation_energy_j), " *
        "initial_stored_energy_j=$(runtime.initial_stored_energy_j), " *
        "supplied_energy_j=$(accepted_state.supplied_energy_j), " *
        "stored_energy_j=$(_transformer_total_stored_energy(accepted_state)), " *
        "device_dissipation_j=$(_transformer_total_dissipated_energy(accepted_state)), " *
        "external_fault_energy_j=$(state.external_fault_energy_j), " *
        "internal_fault_energy_j=$(state.internal_fault_energy_j), " *
        "tolerance_j=$(tolerance)",
    ))
    runtime.companion_method === :backward_euler &&
        (state.numerical_dissipation_energy_j = max(
            state.numerical_dissipation_energy_j,
            balance_residual,
        ))
    state.maximum_energy_balance_residual_j = max(
        state.maximum_energy_balance_residual_j,
        abs(balance_residual - state.numerical_dissipation_energy_j),
    )
    return state
end

function transformer_apparatus_event_state(runtime::TransformerApparatusRuntime)
    state = runtime.event_state
    return (
        terminal_closed=copy(state.terminal_closed),
        terminal_transform=copy(state.terminal_transform),
        active_terminal_faults=deepcopy(state.active_terminal_faults),
        active_winding_faults=deepcopy(state.active_winding_faults),
        active_grounding_paths=deepcopy(state.active_grounding_paths),
        active_internal_faults=deepcopy(state.active_internal_faults),
        applied_event_ids=copy(state.applied_event_ids),
        topology_transition_count=state.topology_transition_count,
        tap_change_count=state.tap_change_count,
        phase_shift_change_count=state.phase_shift_change_count,
        event_energy_j=state.event_energy_j,
        external_fault_energy_j=state.external_fault_energy_j,
        internal_fault_energy_j=state.internal_fault_energy_j,
        numerical_dissipation_energy_j=state.numerical_dissipation_energy_j,
        maximum_energy_balance_residual_j=state.maximum_energy_balance_residual_j,
    )
end

transformer_apparatus_event_occurrences(runtime::TransformerApparatusRuntime) =
    copy(runtime.event_state.occurrences)

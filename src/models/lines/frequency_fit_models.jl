
function _validate_distributed_modal_runtime_inputs(
    state::DistributedTransposedLineModalBranchState,
    history::DistributedTransposedLineHistoryState,
)
    state.group_length == history.phase_count ||
        throw(ArgumentError("modal branch state and history state phase counts must match"))
    state.modal_sequence_indices == history.modal_sequence_indices ||
        throw(ArgumentError("modal branch state and history state sequence order must match"))
    length(history.modal_history_from_values) == state.group_length ||
        throw(ArgumentError("from-side modal history count must match modal branch state"))
    length(history.modal_history_to_values) == state.group_length ||
        throw(ArgumentError("to-side modal history count must match modal branch state"))
    length(history.history_read_indices) == state.group_length ||
        throw(ArgumentError("history read cursor count must match modal branch state"))
    return nothing
end

function _resolve_modal_history_interpolation_factors(
    state::DistributedTransposedLineModalBranchState,
    history::DistributedTransposedLineHistoryState,
    override::Union{Nothing,AbstractVector{<:Real}},
)
    if override === nothing
        factors = isempty(history.modal_history_interpolation_factors) ?
            _distributed_history_interpolation_factors!(
                Float64[],
                state,
                history.timestep_s,
            ) :
            copy(history.modal_history_interpolation_factors)
    else
        factors = _copy_real_modal_values(
            "modal_history_interpolation_factors",
            override,
            state.group_length,
        )
    end
    length(factors) == state.group_length ||
        throw(ArgumentError("modal_history_interpolation_factors count must match modal state"))
    all(value -> 0.0 <= value <= 1.0, factors) ||
        throw(ArgumentError("modal_history_interpolation_factors must be between zero and one"))
    return factors
end

function _phase_current_injection_from_modal(
    modal_currents::AbstractVector{Float64},
    modal_transform::AbstractMatrix{Float64},
)
    return Vector{Float64}(modal_transform * modal_currents)
end

function _packed_history_offset(
    history::DistributedTransposedLineHistoryState,
    slot::Int,
)
    offset = slot - history.history_storage_start_index + 1
    offset > 0 ||
        throw(ArgumentError("packed history slot precedes history storage start index"))
    return offset
end

function _ensure_packed_history_slot!(
    values::Vector{Float64},
    history::DistributedTransposedLineHistoryState,
    slot::Int,
)
    offset = _packed_history_offset(history, slot)
    if offset > length(values)
        old_length = length(values)
        resize!(values, offset)
        for index in (old_length + 1):offset
            values[index] = 0.0
        end
    end
    return offset
end

function _packed_history_value(
    values::Vector{Float64},
    history::DistributedTransposedLineHistoryState,
    slot::Int,
)
    return values[_ensure_packed_history_slot!(values, history, slot)]
end

function _set_packed_history_value!(
    values::Vector{Float64},
    history::DistributedTransposedLineHistoryState,
    slot::Int,
    value::Float64,
)
    values[_ensure_packed_history_slot!(values, history, slot)] = value
    return value
end

function distributed_transposed_line_modal_timestep_update!(
    state::DistributedTransposedLineModalBranchState,
    history::DistributedTransposedLineHistoryState,
    terminal_voltage_from::AbstractVector{<:Real},
    terminal_voltage_to::AbstractVector{<:Real};
    modal_history_interpolation_factors::Union{Nothing,AbstractVector{<:Real}}=nothing,
    name::Symbol=Symbol(string(history.name), "_modal_timestep_update"),
    collect_diagnostics::Bool=true,
)
    _validate_distributed_modal_runtime_inputs(state, history)
    phase_count = state.group_length
    phase_from = _copy_real_modal_values(
        "terminal_voltage_from",
        terminal_voltage_from,
        phase_count,
    )
    phase_to = _copy_real_modal_values(
        "terminal_voltage_to",
        terminal_voltage_to,
        phase_count,
    )
    interpolation_factors = _resolve_modal_history_interpolation_factors(
        state,
        history,
        modal_history_interpolation_factors,
    )
    modal_to_phase_transform = state.modal_to_phase_transform_matrix
    phase_to_modal_transform = state.phase_to_modal_transform_matrix
    modal_voltage_from = Vector{Float64}(phase_to_modal_transform * phase_from)
    modal_voltage_to = Vector{Float64}(phase_to_modal_transform * phase_to)
    norton_coefficients = distributed_transposed_line_norton_coefficients(state)
    if !isempty(history.packed_history_from_values)
        isempty(history.packed_history_to_values) &&
            throw(ArgumentError("packed from/to history arrays must both be present"))
        length(unique(history.history_read_indices)) == 1 ||
            throw(ArgumentError("packed distributed-line histories require one shared read cursor"))
        modal_current_from = zeros(Float64, phase_count)
        modal_current_to = zeros(Float64, phase_count)
        modal_injection_from = zeros(Float64, phase_count)
        modal_injection_to = zeros(Float64, phase_count)
        read_before = copy(history.history_read_indices)
        read_after = similar(read_before)
        write_indices = similar(read_before)
        next_write_indices = similar(read_before)
        cursor = read_before[1] + 1
        for mode_index in 1:phase_count
            base_slot = history.storage_start_indices[mode_index]
            delay_count =
                history.storage_end_indices[mode_index] -
                history.storage_start_indices[mode_index] - 2
            delay_count >= 0 ||
                throw(ArgumentError("packed distributed-line delay count must be nonnegative"))
            read_slot = base_slot + cursor
            second_read_slot = read_slot + 1
            write_slot = read_slot + delay_count
            next_write_slot = write_slot + 1
            write_indices[mode_index] = write_slot
            next_write_indices[mode_index] = next_write_slot

            modal_admittance =
                norton_coefficients.modal_history_source_admittance_values[mode_index]
            modal_damping =
                norton_coefficients.modal_history_damping_values[mode_index]
            old_from_write = _packed_history_value(
                history.packed_history_from_values,
                history,
                write_slot,
            )
            old_to_write = _packed_history_value(
                history.packed_history_to_values,
                history,
                write_slot,
            )
            modal_injection_from[mode_index] = old_from_write
            modal_injection_to[mode_index] = old_to_write
            _set_packed_history_value!(
                history.packed_history_to_values,
                history,
                write_slot,
                modal_admittance * modal_voltage_from[mode_index] -
                modal_damping * old_from_write,
            )
            _set_packed_history_value!(
                history.packed_history_from_values,
                history,
                write_slot,
                modal_admittance * modal_voltage_to[mode_index] -
                modal_damping * old_to_write,
            )

            factor = interpolation_factors[mode_index]
            incident_from =
                _packed_history_value(
                    history.packed_history_from_values,
                    history,
                    read_slot,
                ) * factor +
                _packed_history_value(
                    history.packed_history_from_values,
                    history,
                    second_read_slot,
                ) * (1.0 - factor)
            incident_to =
                _packed_history_value(
                    history.packed_history_to_values,
                    history,
                    read_slot,
                ) * factor +
                _packed_history_value(
                    history.packed_history_to_values,
                    history,
                    second_read_slot,
                ) * (1.0 - factor)
            if state.modal_signed_characteristic_impedances[mode_index] >= 0.0
                total_incident = incident_from + incident_to
                reflected_difference = (incident_from - incident_to) * modal_damping
                incident_from = 0.5 * (total_incident + reflected_difference)
                incident_to = 0.5 * (total_incident - reflected_difference)
            end
            _set_packed_history_value!(
                history.packed_history_from_values,
                history,
                next_write_slot,
                incident_from,
            )
            _set_packed_history_value!(
                history.packed_history_to_values,
                history,
                next_write_slot,
                incident_to,
            )
            modal_current_from[mode_index] = incident_from
            modal_current_to[mode_index] = incident_to
            read_after[mode_index] = cursor
        end
        history.history_read_indices .= read_after
        phase_current_from =
            _phase_current_injection_from_modal(modal_injection_from, modal_to_phase_transform)
        phase_current_to =
            _phase_current_injection_from_modal(modal_injection_to, modal_to_phase_transform)
        return DistributedTransposedLineModalTimestepUpdate(
            name,
            copy(state.phase_indices),
            copy(state.line_numbers),
            copy(state.modal_sequence_indices),
            modal_voltage_from,
            modal_voltage_to,
            modal_current_from,
            modal_current_to,
            phase_current_from,
            phase_current_to,
            read_before,
            read_after,
            write_indices,
            next_write_indices,
            interpolation_factors,
            phase_count,
        )
    end
    modal_current_from = zeros(Float64, phase_count)
    modal_current_to = zeros(Float64, phase_count)
    read_before = copy(history.history_read_indices)
    read_after = similar(read_before)
    write_indices = similar(read_before)
    next_write_indices = similar(read_before)

    for mode_index in 1:phase_count
        from_history = history.modal_history_from_values[mode_index]
        to_history = history.modal_history_to_values[mode_index]
        sample_count = history.history_sample_counts[mode_index]
        length(from_history) == sample_count && length(to_history) == sample_count ||
            throw(ArgumentError("modal history lengths must match history_sample_counts"))
        sample_count >= 2 ||
            throw(ArgumentError("modal history update requires at least two samples"))

        read_index = history.history_read_indices[mode_index]
        1 <= read_index <= sample_count ||
            throw(ArgumentError("history_read_indices entries must reference modal histories"))
        second_read_index = read_index == sample_count ? 1 : read_index + 1
        write_index = read_index == 1 ? sample_count : read_index - 1
        next_write_index = write_index == sample_count ? 1 : write_index + 1
        write_indices[mode_index] = write_index
        next_write_indices[mode_index] = next_write_index

        modal_admittance =
            norton_coefficients.modal_history_source_admittance_values[mode_index]
        modal_damping =
            norton_coefficients.modal_history_damping_values[mode_index]
        old_from_write = from_history[write_index]
        old_to_write = to_history[write_index]
        to_history[write_index] =
            modal_admittance * modal_voltage_from[mode_index] -
            modal_damping * old_from_write
        from_history[write_index] =
            modal_admittance * modal_voltage_to[mode_index] -
            modal_damping * old_to_write

        factor = interpolation_factors[mode_index]
        incident_from =
            from_history[read_index] * factor +
            from_history[second_read_index] * (1.0 - factor)
        incident_to =
            to_history[read_index] * factor +
            to_history[second_read_index] * (1.0 - factor)
        if state.modal_signed_characteristic_impedances[mode_index] >= 0.0
            total_incident = incident_from + incident_to
            reflected_difference = (incident_from - incident_to) * modal_damping
            incident_from = 0.5 * (total_incident + reflected_difference)
            incident_to = 0.5 * (total_incident - reflected_difference)
        end
        from_history[next_write_index] = incident_from
        to_history[next_write_index] = incident_to
        modal_current_from[mode_index] = incident_from
        modal_current_to[mode_index] = incident_to
        history.history_read_indices[mode_index] = second_read_index
        read_after[mode_index] = second_read_index
    end

    collect_diagnostics || return nothing
    phase_current_from =
        _phase_current_injection_from_modal(modal_current_from, modal_to_phase_transform)
    phase_current_to =
        _phase_current_injection_from_modal(modal_current_to, modal_to_phase_transform)
    return DistributedTransposedLineModalTimestepUpdate(
        name,
        copy(state.phase_indices),
        copy(state.line_numbers),
        copy(state.modal_sequence_indices),
        modal_voltage_from,
        modal_voltage_to,
        modal_current_from,
        modal_current_to,
        phase_current_from,
        phase_current_to,
        read_before,
        read_after,
        write_indices,
        next_write_indices,
        interpolation_factors,
        phase_count,
    )
end

function distributed_transposed_line_history_current_injection!(
    rhs::AbstractVector{Float64},
    state::DistributedTransposedLineModalBranchState,
    history::DistributedTransposedLineHistoryState;
    name::Symbol=Symbol(string(history.name), "_history_current_injection"),
    collect_diagnostics::Bool=true,
)
    _validate_distributed_modal_runtime_inputs(state, history)
    phase_count = state.group_length
    modal_current_from = zeros(Float64, phase_count)
    modal_current_to = zeros(Float64, phase_count)

    if !isempty(history.packed_history_from_values)
        isempty(history.packed_history_to_values) &&
            throw(ArgumentError("packed from/to history arrays must both be present"))
        length(unique(history.history_read_indices)) == 1 ||
            throw(ArgumentError("packed distributed-line histories require one shared read cursor"))
        cursor = history.history_read_indices[1] + 1
        for mode_index in 1:phase_count
            base_slot = history.storage_start_indices[mode_index]
            delay_count =
                history.storage_end_indices[mode_index] -
                history.storage_start_indices[mode_index] - 2
            delay_count >= 0 ||
                throw(ArgumentError("packed distributed-line delay count must be nonnegative"))
            write_slot = base_slot + cursor + delay_count
            modal_current_from[mode_index] = _packed_history_value(
                history.packed_history_from_values,
                history,
                write_slot,
            )
            modal_current_to[mode_index] = _packed_history_value(
                history.packed_history_to_values,
                history,
                write_slot,
            )
        end
    else
        for mode_index in 1:phase_count
            from_history = history.modal_history_from_values[mode_index]
            to_history = history.modal_history_to_values[mode_index]
            sample_count = history.history_sample_counts[mode_index]
            length(from_history) == sample_count && length(to_history) == sample_count ||
                throw(ArgumentError("modal history lengths must match history_sample_counts"))
            read_index = history.history_read_indices[mode_index]
            1 <= read_index <= sample_count ||
                throw(ArgumentError("history_read_indices entries must reference modal histories"))
            write_index = read_index == 1 ? sample_count : read_index - 1
            modal_current_from[mode_index] = from_history[write_index]
            modal_current_to[mode_index] = to_history[write_index]
        end
    end

    phase_current_from =
        _phase_current_injection_from_modal(
            modal_current_from,
            state.modal_to_phase_transform_matrix,
        )
    phase_current_to =
        _phase_current_injection_from_modal(
            modal_current_to,
            state.modal_to_phase_transform_matrix,
        )
    from_nodes = _distributed_line_phase_node_indices(
        "from_node_indices",
        state.from_node_indices,
        phase_count,
        length(rhs),
    )
    to_nodes = _distributed_line_phase_node_indices(
        "to_node_indices",
        state.to_node_indices,
        phase_count,
        length(rhs),
    )
    rhs_before = collect_diagnostics ? copy(rhs) : Float64[]
    rhs_update_count = 0
    for phase_index in 1:phase_count
        from_node = from_nodes[phase_index]
        if from_node != 0
            rhs[from_node] += phase_current_from[phase_index]
            rhs_update_count += 1
        end
        to_node = to_nodes[phase_index]
        if to_node != 0
            rhs[to_node] += phase_current_to[phase_index]
            rhs_update_count += 1
        end
    end
    collect_diagnostics || return nothing
    return DistributedTransposedLinePhaseCurrentInjection(
        name,
        from_nodes,
        to_nodes,
        phase_current_from,
        phase_current_to,
        rhs_before,
        copy(rhs),
        rhs_update_count,
        phase_count,
    )
end

function distributed_transposed_line_modal_timestep_update(
    state::DistributedTransposedLineModalBranchState,
    history::DistributedTransposedLineHistoryState,
    terminal_voltage_from::AbstractVector{<:Real},
    terminal_voltage_to::AbstractVector{<:Real};
    modal_history_interpolation_factors::Union{Nothing,AbstractVector{<:Real}}=nothing,
    name::Symbol=Symbol(string(history.name), "_modal_timestep_update"),
)
    copied_history = DistributedTransposedLineHistoryState(
        history.name,
        copy(history.phase_indices),
        copy(history.line_numbers),
        copy(history.modal_sequence_indices),
        history.timestep_s,
        history.steady_state_frequency_hz,
        history.angular_step_rad,
        history.history_storage_start_index,
        history.next_history_storage_index,
        copy(history.storage_start_indices),
        copy(history.storage_end_indices),
        copy(history.storage_lengths),
        copy(history.history_sample_counts),
        copy(history.modal_history_interpolation_factors),
        copy(history.history_read_indices),
        copy(history.outgoing_wave_from_real_values),
        copy(history.outgoing_wave_from_imag_values),
        copy(history.outgoing_wave_to_real_values),
        copy(history.outgoing_wave_to_imag_values),
        [copy(values) for values in history.modal_history_from_values],
        [copy(values) for values in history.modal_history_to_values],
        copy(history.packed_history_from_values),
        copy(history.packed_history_to_values),
        history.initialized_from_steady_state,
        history.phase_count,
    )
    return distributed_transposed_line_modal_timestep_update!(
        state,
        copied_history,
        terminal_voltage_from,
        terminal_voltage_to;
        modal_history_interpolation_factors = modal_history_interpolation_factors,
        name = name,
    )
end

function _distributed_line_phase_node_indices(
    name::AbstractString,
    indices::AbstractVector{<:Integer},
    phase_count::Int,
    rhs_length::Int,
)
    length(indices) == phase_count ||
        throw(ArgumentError("$name must have $phase_count entries"))
    converted = Int.(indices)
    all(index -> 0 <= index <= rhs_length, converted) ||
        throw(ArgumentError("$name entries must be zero or valid RHS node indices"))
    return converted
end

function distributed_transposed_line_phase_current_injection!(
    rhs::AbstractVector{Float64},
    update::DistributedTransposedLineModalTimestepUpdate,
    from_node_indices::AbstractVector{<:Integer},
    to_node_indices::AbstractVector{<:Integer};
    name::Symbol=Symbol(string(update.name), "_phase_current_injection"),
)
    phase_count = update.phase_count
    length(update.phase_current_injection_from_values) == phase_count ||
        throw(ArgumentError("from-side phase current injection count must match phase count"))
    length(update.phase_current_injection_to_values) == phase_count ||
        throw(ArgumentError("to-side phase current injection count must match phase count"))
    from_nodes = _distributed_line_phase_node_indices(
        "from_node_indices",
        from_node_indices,
        phase_count,
        length(rhs),
    )
    to_nodes = _distributed_line_phase_node_indices(
        "to_node_indices",
        to_node_indices,
        phase_count,
        length(rhs),
    )
    rhs_before = copy(rhs)
    rhs_update_count = 0
    for phase_index in 1:phase_count
        from_node = from_nodes[phase_index]
        if from_node != 0
            rhs[from_node] += update.phase_current_injection_from_values[phase_index]
            rhs_update_count += 1
        end
        to_node = to_nodes[phase_index]
        if to_node != 0
            rhs[to_node] += update.phase_current_injection_to_values[phase_index]
            rhs_update_count += 1
        end
    end
    return DistributedTransposedLinePhaseCurrentInjection(
        name,
        from_nodes,
        to_nodes,
        copy(update.phase_current_injection_from_values),
        copy(update.phase_current_injection_to_values),
        rhs_before,
        copy(rhs),
        rhs_update_count,
        phase_count,
    )
end

function distributed_transposed_line_phase_current_injection(
    rhs_values::AbstractVector{<:Real},
    update::DistributedTransposedLineModalTimestepUpdate,
    from_node_indices::AbstractVector{<:Integer},
    to_node_indices::AbstractVector{<:Integer};
    name::Symbol=Symbol(string(update.name), "_phase_current_injection"),
)
    rhs = _copy_real_modal_values("rhs_values", rhs_values, length(rhs_values))
    return distributed_transposed_line_phase_current_injection!(
        rhs,
        update,
        from_node_indices,
        to_node_indices;
        name = name,
    )
end

function _distributed_transposed_line_norton_coefficient(
    signed_characteristic_impedance::Float64,
    total_resistance::Float64,
    modal_admittance_denominator::Float64,
)
    modal_admittance_denominator > 0.0 ||
        throw(ArgumentError("modal admittance denominator must be positive"))
    isfinite(signed_characteristic_impedance) &&
        signed_characteristic_impedance != 0.0 ||
        throw(ArgumentError("modal characteristic impedance must be finite and nonzero"))
    isfinite(total_resistance) ||
        throw(ArgumentError("modal total resistance must be finite"))
    characteristic_impedance = abs(signed_characteristic_impedance)
    resistance = abs(total_resistance)
    if signed_characteristic_impedance < 0.0
        return (
            companion_admittance =
                2.0 / (characteristic_impedance * modal_admittance_denominator),
            history_source_admittance =
                2.0 / (characteristic_impedance * modal_admittance_denominator),
            history_damping = exp(-0.5 * resistance / characteristic_impedance),
        )
    end
    shifted_impedance = characteristic_impedance + 0.25 * resistance
    shifted_impedance > 0.0 ||
        throw(ArgumentError("distributed-line shifted impedance must be positive"))
    source_impedance =
        characteristic_impedance +
        0.5 * resistance +
        resistance * resistance / (16.0 * characteristic_impedance)
    source_impedance > 0.0 ||
        throw(ArgumentError("distributed-line source impedance must be positive"))
    damping = (shifted_impedance - 0.5 * resistance) / shifted_impedance
    damping > 0.0 ||
        throw(ArgumentError("distributed-line history damping must be positive"))
    return (
        companion_admittance =
            2.0 / (shifted_impedance * modal_admittance_denominator),
        history_source_admittance =
            2.0 / (source_impedance * modal_admittance_denominator),
        history_damping = damping,
    )
end

function distributed_transposed_line_norton_coefficients(
    state::DistributedTransposedLineModalBranchState,
)
    admittances = Float64[]
    history_source_admittances = Float64[]
    damping = Float64[]
    for index in 1:state.group_length
        coefficient = _distributed_transposed_line_norton_coefficient(
            state.modal_signed_characteristic_impedances[index],
            state.modal_total_resistances[index],
            state.modal_admittance_denominator,
        )
        push!(admittances, coefficient.companion_admittance)
        push!(history_source_admittances, coefficient.history_source_admittance)
        push!(damping, coefficient.history_damping)
    end
    return (
        source = :distributed_transposed_line_norton_coefficients,
        outcome = :norton_history_coefficients,
        modal_companion_admittance_values = admittances,
        modal_history_source_admittance_values = history_source_admittances,
        modal_history_damping_values = damping,
        phase_count = state.phase_count,
    )
end

function distributed_transposed_line_companion_admittance(
    state::DistributedTransposedLineModalBranchState;
    name::Symbol=Symbol(string(state.name), "_companion_admittance"),
)
    state.group_length == state.phase_count ||
        throw(ArgumentError("distributed transposed line companion admittance requires one modal row per phase"))
    modal_values =
        distributed_transposed_line_norton_coefficients(state).
        modal_companion_admittance_values
    all(value -> isfinite(value) && value > 0.0, modal_values) ||
        throw(ArgumentError("modal companion admittance values must be finite and positive"))
    phase_matrix =
        state.modal_to_phase_transform_matrix *
        Diagonal(modal_values) *
        state.phase_to_modal_transform_matrix
    return DistributedTransposedLineCompanionAdmittance(
        name,
        copy(state.from_node_indices),
        copy(state.to_node_indices),
        copy(state.line_numbers),
        copy(state.modal_sequence_indices),
        Float64.(modal_values),
        Matrix{Float64}(phase_matrix),
        state.phase_count,
    )
end

function distributed_transposed_line_companion_admittance_stamp!(
    y::AbstractMatrix{Float64},
    admittance::DistributedTransposedLineCompanionAdmittance,
)
    phase_count = admittance.phase_count
    size(admittance.phase_companion_admittance_matrix, 1) == phase_count &&
        size(admittance.phase_companion_admittance_matrix, 2) == phase_count ||
        throw(ArgumentError("distributed transposed line phase admittance matrix size mismatch"))
    size(y, 1) == size(y, 2) ||
        throw(ArgumentError("nodal admittance matrix must be square"))
    node_limit = size(y, 1)
    all(node -> 0 <= node <= node_limit, admittance.from_node_indices) ||
        throw(ArgumentError("from_node_indices entries must be zero or valid matrix node indices"))
    all(node -> 0 <= node <= node_limit, admittance.to_node_indices) ||
        throw(ArgumentError("to_node_indices entries must be zero or valid matrix node indices"))
    for col in 1:phase_count
        for row in 1:phase_count
            value =
                DISTRIBUTED_LINE_ENDPOINT_NORTON_SHARE *
                admittance.phase_companion_admittance_matrix[row, col]
            from_row = admittance.from_node_indices[row]
            from_col = admittance.from_node_indices[col]
            to_row = admittance.to_node_indices[row]
            to_col = admittance.to_node_indices[col]
            from_row != 0 && from_col != 0 && (y[from_row, from_col] += value)
            to_row != 0 && to_col != 0 && (y[to_row, to_col] += value)
        end
    end
    return y
end

function stamp!(
    y,
    rhs,
    admittance::DistributedTransposedLineCompanionAdmittance,
    t::Float64,
    dt::Float64,
)
    distributed_transposed_line_companion_admittance_stamp!(y, admittance)
    return nothing
end

function update!(admittance::DistributedTransposedLineCompanionAdmittance, v, dt::Float64)
    return nothing
end

function terminal_surge_impedance_admittance_stamp!(
    y::AbstractMatrix{Float64},
    admittance::TerminalSurgeImpedanceAdmittance,
)
    phase_count = admittance.phase_count
    size(admittance.stamped_admittance_matrix, 1) == phase_count &&
        size(admittance.stamped_admittance_matrix, 2) == phase_count ||
        throw(ArgumentError("terminal surge stamped admittance matrix size mismatch"))
    size(y, 1) == size(y, 2) ||
        throw(ArgumentError("nodal admittance matrix must be square"))
    node_limit = size(y, 1)
    all(node -> 0 <= node <= node_limit, admittance.node_indices) ||
        throw(ArgumentError("terminal surge node indices must be zero or valid matrix node indices"))
    for col in 1:phase_count
        for row in 1:phase_count
            stamp_admittance_entry!(
                y,
                admittance.node_indices[row],
                0,
                admittance.node_indices[col],
                0,
                admittance.stamped_admittance_matrix[row, col],
            )
        end
    end
    return y
end

function stamp!(
    y,
    rhs,
    admittance::TerminalSurgeImpedanceAdmittance,
    t::Float64,
    dt::Float64,
)
    terminal_surge_impedance_admittance_stamp!(y, admittance)
    return nothing
end

function update!(admittance::TerminalSurgeImpedanceAdmittance, v, dt::Float64)
    return nothing
end

function _transposed_line_modal_transform()
    return [
        inv(sqrt(3.0)) inv(sqrt(2.0)) inv(sqrt(6.0))
        inv(sqrt(3.0)) -inv(sqrt(2.0)) inv(sqrt(6.0))
        inv(sqrt(3.0)) 0.0 -2.0 * inv(sqrt(6.0))
    ]
end

function _steady_state_storage_conversion_factors(
    x_frequency_hz::Float64,
    c_frequency_hz::Float64,
    steady_state_frequency_hz::Float64,
)
    steady_state_frequency_hz > 0.0 ||
        throw(ArgumentError("steady_state_frequency_hz must be positive"))
    inductive_scale = 1000.0 / (2.0 * pi)
    capacitive_scale = 1000.0 * inductive_scale
    if x_frequency_hz > 0.0
        inductive_scale = x_frequency_hz * 1.0e-3 * 2.0 * pi * inductive_scale
    end
    if c_frequency_hz > 0.0
        capacitive_scale = c_frequency_hz * 2.0 * pi * capacitive_scale
    end
    return (
        series_reactance_scale = inductive_scale / steady_state_frequency_hz,
        shunt_susceptance_scale = capacitive_scale / steady_state_frequency_hz,
    )
end

function _distributed_line_modal_pi_term(
    signed_characteristic_impedance::Float64,
    total_resistance::Float64,
    propagation_time_s::Float64,
    angular_frequency_rad_s::Float64,
)
    angle = angular_frequency_rad_s * propagation_time_s
    surge_impedance = signed_characteristic_impedance
    resistance = total_resistance
    if resistance < 0.0
        resistance = -resistance
        denr = angle * surge_impedance
        denx = angle / surge_impedance
        quarter_pi = 0.5 * pi
        magnitude = sqrt(resistance^2 + denr^2)
        phase = atan(denr, resistance)
        half_angle = 0.5 * (phase - quarter_pi)
        denominator_magnitude = sqrt(magnitude / denx)
        rzero = denominator_magnitude * cos(half_angle)
        xzero = denominator_magnitude * sin(half_angle)
        rotated_angle = half_angle + quarter_pi
        numerator_magnitude = sqrt(magnitude * denx)
        exponent_real = numerator_magnitude * cos(rotated_angle)
        exponent_imag = numerator_magnitude * sin(rotated_angle)
        positive = exp(exponent_real)
        negative = inv(positive)
        cosh_real = 0.5 * (positive + negative)
        sinh_real = cosh_real - negative
        sin_imag = sin(exponent_imag)
        cos_imag = cos(exponent_imag)
        sinhgr = sinh_real * cos_imag
        sinhgi = cosh_real * sin_imag
        denr = rzero * sinhgr - xzero * sinhgi
        denx = rzero * sinhgi + xzero * sinhgr
        numerator_real = cosh_real * cos_imag - 1.0
        numerator_imag = sinh_real * sin_imag
    elseif surge_impedance >= 0.0
        ratio = resistance / surge_impedance
        half_angle = 0.5 * angle
        sine = sin(half_angle)
        cosine = cos(half_angle)
        ratio_squared = ratio^2
        sine_squared = sine^2
        denr = (1.0 - sine_squared * ((48.0 + ratio_squared) / 32.0)) * resistance
        ratio_squared_over_eight = ratio_squared / 8.0
        denx = (3.0 * ratio_squared_over_eight + 2.0) * sine * cosine *
            surge_impedance
        numerator_real = (-2.0 - ratio_squared_over_eight) * sine_squared
        numerator_imag = ratio * sine * cosine
    else
        half_decay = 0.5 * (resistance / surge_impedance)
        sine = sin(angle)
        cosine = cos(angle)
        positive = exp(-half_decay)
        negative = inv(positive)
        even_part = positive + negative
        odd_part = positive - negative
        numerator_real = cosine * even_part - 2.0
        numerator_imag = sine * odd_part
        denr = cosine * odd_part
        denx = sine * even_part
        denominator_norm = denr^2 + denx^2
        scale = -0.5 * surge_impedance
        series_resistance = scale * denr
        series_reactance = scale * denx
        denominator = -surge_impedance * denominator_norm
        shunt_conductance =
            (numerator_real * denr + numerator_imag * denx) / denominator
        shunt_susceptance =
            (numerator_imag * denr - numerator_real * denx) / denominator
        return (
            series_impedance = complex(series_resistance, series_reactance),
            shunt_admittance = complex(shunt_conductance, shunt_susceptance),
        )
    end
    denominator = denr^2 + denx^2
    shunt_conductance = (numerator_real * denr + numerator_imag * denx) / denominator
    shunt_susceptance = (numerator_imag * denr - numerator_real * denx) / denominator
    return (
        series_impedance = complex(denr, denx),
        shunt_admittance = complex(shunt_conductance, shunt_susceptance),
    )
end

function _lower_triangular_complex_values(matrix::AbstractMatrix{ComplexF64})
    row_count, column_count = size(matrix)
    row_count == column_count || throw(ArgumentError("matrix must be square"))
    row_indices = Int[]
    column_indices = Int[]
    values = ComplexF64[]
    for row in 1:row_count
        for column in 1:row
            push!(row_indices, row)
            push!(column_indices, column)
            push!(values, matrix[row, column])
        end
    end
    return (rows = row_indices, columns = column_indices, values = values)
end

function distributed_transposed_line_steady_state_pi_equivalent(
    constants::DistributedTransposedLineConstants;
    steady_state_frequency_hz::Real,
    storage_start_index::Integer=1,
    name::Symbol=Symbol(string(constants.name), "_steady_state_pi_equivalent"),
)
    constants.phase_indices == [1, 2, 3] ||
        throw(ArgumentError("distributed transposed line pi equivalent requires phase indices [1, 2, 3]"))
    start_index = Int(storage_start_index)
    start_index > 0 || throw(ArgumentError("storage_start_index must be positive"))
    steady_frequency =
        _positive_finite_impedance_value("steady_state_frequency_hz", steady_state_frequency_hz)
    angular_frequency = 2.0 * pi * steady_frequency
    state = distributed_transposed_line_modal_branch_state(constants)
    modal_terms = [
        _distributed_line_modal_pi_term(
            state.modal_signed_characteristic_impedances[index],
            state.modal_total_resistances[index],
            state.modal_propagation_times_s[index],
            angular_frequency,
        )
        for index in 1:state.group_length
    ]
    modal_transform = _transposed_line_modal_transform()
    series_admittances = inv.(ComplexF64[term.series_impedance for term in modal_terms])
    shunt_admittances = ComplexF64[term.shunt_admittance for term in modal_terms]
    phase_series_admittance =
        modal_transform * Diagonal(series_admittances) * transpose(modal_transform)
    phase_series_impedance = inv(phase_series_admittance)
    phase_shunt_admittance =
        modal_transform * Diagonal(shunt_admittances) * transpose(modal_transform)
    conversion = _steady_state_storage_conversion_factors(
        constants.x_frequency_hz,
        constants.c_frequency_hz,
        steady_frequency,
    )
    series_storage = _lower_triangular_complex_values(phase_series_impedance)
    shunt_storage = _lower_triangular_complex_values(phase_shunt_admittance)
    storage_row_count = length(series_storage.values)
    storage_indices = collect(start_index:(start_index + storage_row_count - 1))
    return DistributedTransposedLineSteadyStatePiEquivalent(
        name,
        copy(constants.phase_indices),
        copy(constants.line_numbers),
        copy(state.modal_sequence_indices),
        steady_frequency,
        angular_frequency,
        constants.x_frequency_hz,
        constants.c_frequency_hz,
        start_index,
        last(storage_indices),
        storage_indices,
        series_storage.rows,
        series_storage.columns,
        real.(series_storage.values),
        imag.(series_storage.values) .* conversion.series_reactance_scale,
        real.(shunt_storage.values),
        imag.(shunt_storage.values) .* conversion.shunt_susceptance_scale,
        Matrix{ComplexF64}(phase_series_impedance),
        Matrix{ComplexF64}(phase_shunt_admittance),
        constants.phase_count,
    )
end

struct LineFrequencyPoint
    frequency_hz::Float64
    characteristic_impedance::ComplexF64
    propagation_constant::ComplexF64
    propagation_factor::ComplexF64
end

struct SemlyenMartiFrequencyScan
    start_frequency_hz::Float64
    decade_count::Int
    intervals_per_decade::Int
    frequency_count::Int
    frequencies_hz::Vector{Float64}
    decade_base_frequencies_hz::Vector{Float64}
    decade_indices::Vector{Int}
    interval_indices::Vector{Int}
    logarithmic_ratio::Float64
end

struct CableFrequencyScanLoopSchedule
    start_frequency_hz::Float64
    steady_state_frequency_hz::Float64
    decade_count::Int
    points_per_decade::Int
    frequency_count::Int
    frequencies_hz::Vector{Float64}
    decade_base_frequencies_hz::Vector{Float64}
    decade_indices::Vector{Int}
    interval_indices::Vector{Int}
    output_indices::Vector{Int}
    logarithmic_ratio::Float64
    initial_earth_resistivity_ohm_m::Float64
    final_earth_resistivity_ohm_m::Float64
    earth_resistivity_values_ohm_m::Vector{Float64}
    initial_distance_m::Float64
    final_distance_m::Float64
    distance_values_m::Vector{Float64}
    retained_distance_m::Float64
    retained_earth_resistivity_ohm_m::Float64
    voltage_control_values::Vector{Float64}
    resistivity_history_values::Vector{Float64}
    resistivity_history_frequencies_hz::Vector{Float64}
    unit_distance_record_count::Int
    unit_distance_km::Float64
    card_output_flag::Int
    transform_flag::Int
    alteration_mode::Int
    modal_output_enabled::Bool
    retry_requested::Bool
    final_retry_flag::Int
    output_row_count::Int
    loop_executed::Bool
end

struct LineModalTransform
    phase_to_modal::Matrix{ComplexF64}
    modal_to_phase::Matrix{ComplexF64}
end

mutable struct LineModeUnwindState
    sequence::Vector{Int}
    previous_sequence::Vector{Int}
    current_metrics::Matrix{Float64}
    previous_metrics::Matrix{Float64}
    second_previous_metrics::Matrix{Float64}
    update_count::Int
end

struct LineModalSolution
    transform::LineModalTransform
    eigenvalues::Vector{ComplexF64}
    propagation_roots::Vector{ComplexF64}
    modal_eigen_order::Vector{Int}
    mode_sequence::Vector{Int}
    mode_metrics::Matrix{Float64}
    ordered_mode_metrics::Matrix{Float64}
    eigenvector_residual_max_abs_error::Float64
    inverse_product_max_abs_error::Float64
    update_count::Int
end

struct LineModalSolutionScan
    frequencies_hz::Vector{Float64}
    angular_frequencies_rad_s::Vector{Float64}
    yz_matrices::Vector{Matrix{ComplexF64}}
    solutions::Vector{LineModalSolution}
    modal_eigenvalues::Vector{Vector{ComplexF64}}
    modal_eigen_orders::Vector{Vector{Int}}
    mode_sequences::Vector{Vector{Int}}
    frequency_order::Vector{Int}
    mode_order_update_count::Int
    frequency_count::Int
    mode_count::Int
    eigenvector_residual_max_abs_error::Float64
    inverse_product_max_abs_error::Float64
end

struct CableGeometryConductor
    radius_m::Float64
    horizontal_position_m::Float64
    depth_m::Float64
    resistivity_ohm_m::Float64
    relative_permittivity::Float64
    relative_permeability::Float64
end

struct CableGeometryConstants
    conductors::Vector{CableGeometryConductor}
    phase_conductor_counts::Vector{Int}
    grounded_conductor_count::Int
    radius_m::Vector{Float64}
    horizontal_position_m::Vector{Float64}
    depth_m::Vector{Float64}
    direct_distance_m::Matrix{Float64}
    image_distance_m::Matrix{Float64}
    angle_rad::Matrix{Float64}
    potential_log_matrix::Matrix{Float64}
    resistivity_ohm_m::Vector{Float64}
    conductivity_s_per_m::Vector{Float64}
    relative_permittivity::Vector{Float64}
    permittivity_f_per_m::Vector{Float64}
    relative_permeability::Vector{Float64}
    permeability_h_per_m::Vector{Float64}
    conductor_count::Int
    phase_count::Int
end

struct CableElectrostaticAdmittance
    frequency_hz::Float64
    angular_frequency_rad_s::Float64
    relative_permittivity::Float64
    potential_coefficient_matrix::Matrix{Float64}
    capacitance_matrix_f_per_m::Matrix{Float64}
    shunt_admittance_matrix_s_per_m::Matrix{ComplexF64}
    conductor_count::Int
    phase_count::Int
end

struct CablePhaseElectrostaticAdmittance
    frequency_hz::Float64
    angular_frequency_rad_s::Float64
    relative_permittivity::Float64
    conductor_to_phase_matrix::Matrix{Float64}
    phase_capacitance_matrix_f_per_m::Matrix{Float64}
    phase_shunt_admittance_matrix_s_per_m::Matrix{ComplexF64}
    conductor_admittance::CableElectrostaticAdmittance
    conductor_count::Int
    phase_count::Int
end

struct CablePiSectionReportState
    punch_requested::Bool
    total_length_m::Float64
    section_length_m::Float64
    section_count::Int
    crossbonded::Bool
    modeling_kind::Symbol
    sheath_grounding_resistance_ohm::Float64
    report_executed::Bool
end

struct CablePipeSheathDerivedState
    cable_kind::Symbol
    surface_position_kind::Symbol
    phase_count::Int
    conductor_count::Int
    pipe_count::Int
    selected_grounded_conductor_count::Int
    matrix_width_conductor_count::Int
    admittance_conductor_count::Int
    pi_reduction_conductor_count::Int
    outer_boundary_radii_m::Vector{Float64}
    conductor_depths_m::Vector{Float64}
    conductor_distances_m::Vector{Float64}
    conductor_pipe_center_distances_m::Vector{Float64}
    conductor_angles_rad::Vector{Float64}
    pipe_return_included::Bool
    pipe_radii_m::Vector{Float64}
    pipe_resistivity_ohm_m::Float64
    pipe_relative_permeability::Float64
    pipe_inner_insulator_relative_permittivity::Float64
    pipe_outer_insulator_relative_permittivity::Float64
    pipe_outer_insulator_defaulted::Bool
    layer_wave_speeds_m_per_s::Matrix{Float64}
    core_inner_diffusion_factors::Vector{Float64}
    core_outer_diffusion_factors::Vector{Float64}
    sheath_inner_diffusion_factors::Vector{Float64}
    sheath_outer_diffusion_factors::Vector{Float64}
    armor_inner_diffusion_factors::Vector{Float64}
    armor_outer_diffusion_factors::Vector{Float64}
    core_insulation_log_ratios::Vector{Float64}
    sheath_insulation_log_ratios::Vector{Float64}
    armor_insulation_log_ratios::Vector{Float64}
    pipe_inner_conductor_limit_m::Float64
    pipe_radius_allows_inner_conductors::Bool
    layer_count_order_valid::Bool
    no_ungrounded_conductor_stop_required::Bool
    crossbond_layout_validation_passed::Bool
    crossbond_layout_validation_error::Symbol
    pi_model_validation_passed::Bool
    pi_model_validation_error::Symbol
    pi_section_report::CablePiSectionReportState
    derived_state_executed::Bool
end

struct CableConductorSeriesImpedance
    frequency_hz::Float64
    angular_frequency_rad_s::Float64
    conductor_external_impedance_matrix_ohm_per_m::Matrix{ComplexF64}
    conductor_internal_impedance_ohm_per_m::Vector{ComplexF64}
    conductor_skin_effect_delta_impedance_ohm_per_m::Vector{ComplexF64}
    conductor_earth_return_impedance_matrix_ohm_per_m::Matrix{ComplexF64}
    conductor_series_impedance_matrix_ohm_per_m::Matrix{ComplexF64}
    active_conductor_series_impedance_matrix_ohm_per_m::Matrix{ComplexF64}
    conductor_to_phase_average_matrix::Matrix{Float64}
    phase_series_impedance_unreduced_matrix_ohm_per_m::Matrix{ComplexF64}
    phase_series_impedance_matrix_ohm_per_m::Matrix{ComplexF64}
    bounded_skin_effect_executed::Bool
    skin_effect_internal_impedance_executed::Bool
    bounded_earth_return_executed::Bool
    earth_return_impedance_executed::Bool
    grounded_conductor_reduction_executed::Bool
    earth_resistivity_ohm_m::Float64
    grounded_conductor_count::Int
    conductor_count::Int
    phase_count::Int
end

struct CablePhaseLineConstants
    frequency_hz::Float64
    angular_frequency_rad_s::Float64
    line_length_m::Float64
    phase_series_impedance_matrix_ohm_per_m::Matrix{ComplexF64}
    phase_shunt_admittance_matrix_s_per_m::Matrix{ComplexF64}
    phase_zy_matrix_per_m2::Matrix{ComplexF64}
    modal_solution::LineModalSolution
    modal_series_impedance_ohm_per_m::Vector{ComplexF64}
    modal_shunt_admittance_s_per_m::Vector{ComplexF64}
    frequency_points::Vector{LineFrequencyPoint}
    modal_diagonalization_max_abs_error::Float64
    propagation_root_max_abs_error::Float64
    phase_admittance::CablePhaseElectrostaticAdmittance
    phase_count::Int
end

struct CableFrequencyScanSeriesImpedance
    frequencies_hz::Vector{Float64}
    angular_frequencies_rad_s::Vector{Float64}
    frequency_rows::Vector{CableConductorSeriesImpedance}
    phase_series_impedance_matrices_ohm_per_m::Vector{Matrix{ComplexF64}}
    frequency_order::Vector{Int}
    frequency_count::Int
    phase_count::Int
end

struct CableFrequencyScanLineConstants
    frequencies_hz::Vector{Float64}
    angular_frequencies_rad_s::Vector{Float64}
    line_length_m::Float64
    frequency_constants::Vector{CablePhaseLineConstants}
    sample_rows::Vector{Vector{LineFrequencyPoint}}
    frequency_order::Vector{Int}
    frequency_count::Int
    phase_count::Int
    skin_effect_internal_impedance_executed::Bool
    earth_return_impedance_executed::Bool
    grounded_conductor_reduction_executed::Bool
end

struct CableFrequencyScanGeneratedLineConstants
    series_impedance::CableFrequencyScanSeriesImpedance
    line_constants::CableFrequencyScanLineConstants
end

struct FrequencyDependentLineModalResponse
    phase_voltage::Vector{ComplexF64}
    modal_voltage::Vector{ComplexF64}
    modal_admittance::Vector{ComplexF64}
    modal_current::Vector{ComplexF64}
    propagated_modal_current::Vector{ComplexF64}
    phase_current::Vector{ComplexF64}
    propagated_phase_current::Vector{ComplexF64}
end

mutable struct FrequencyDependentLineModalState
    transform::LineModalTransform
    frequency_points::Vector{LineFrequencyPoint}
    phase_voltage::Vector{ComplexF64}
    modal_voltage::Vector{ComplexF64}
    modal_admittance::Vector{ComplexF64}
    modal_current::Vector{ComplexF64}
    propagated_modal_current::Vector{ComplexF64}
    phase_current::Vector{ComplexF64}
    propagated_phase_current::Vector{ComplexF64}
    update_count::Int
end

mutable struct FrequencyDependentLineRuntimeState
    modal_state::FrequencyDependentLineModalState
    from_phase_voltage::Vector{ComplexF64}
    to_phase_voltage::Vector{ComplexF64}
    phase_voltage_difference::Vector{ComplexF64}
    previous_sending_phase_current::Vector{ComplexF64}
    previous_receiving_phase_current::Vector{ComplexF64}
    sending_phase_current::Vector{ComplexF64}
    receiving_phase_current::Vector{ComplexF64}
    update_count::Int
end

struct LineFrequencySampleFitResult
    frequency_hz::Float64
    frequency_points::Vector{LineFrequencyPoint}
    lower_frequency_hz::Float64
    upper_frequency_hz::Float64
    interpolation_weight::Float64
    exact_frequency_match::Bool
    line_length::Float64
    mode_count::Int
end

struct LineRecursiveConvolutionFitResult
    sample_frequencies_hz::Vector{Float64}
    pole_decay::Matrix{ComplexF64}
    residue::Matrix{ComplexF64}
    modal_response_samples::Matrix{ComplexF64}
    fitted_modal_response::Matrix{ComplexF64}
    residual_modal_response::Matrix{ComplexF64}
    max_abs_error::Float64
    dt_s::Float64
    mode_count::Int
    term_count::Int
end

struct LineStepResponseExponentialFitResult
    angular_frequencies_rad_s::Vector{Float64}
    frequency_response_values::Vector{Float64}
    step_response_values::Vector{Float64}
    final_value::Float64
    fit_span_s::Float64
    time_start_s::Float64
    time_step_s::Float64
    first_amplitude::Float64
    first_time_constant_s::Float64
    second_amplitude::Float64
    second_time_constant_s::Float64
    delay_s::Float64
    steady_state_shift::ComplexF64
    normalized_square_error::Float64
    iteration_count::Int
    variable_count::Int
    fit_executed::Bool
    steady_state_adjustment_executed::Bool
end

struct SemlyenLineExponentialConvolutionCoefficients
    pole::ComplexF64
    residue::ComplexF64
    dt_s::Float64
    decay::ComplexF64
    current_gain::ComplexF64
    delayed_gain::ComplexF64
    taylor_series_used::Bool
    taylor_iteration_count::Int
end

struct SemlyenLineHarmonicHistoryUpdate
    pole::ComplexF64
    residue::ComplexF64
    angular_frequency_rad_s::Float64
    transfer_response::ComplexF64
    quadrature_transfer_response::ComplexF64
    from_voltage_phasor::ComplexF64
    to_voltage_phasor::ComplexF64
    previous_from_history::ComplexF64
    previous_to_history::ComplexF64
    from_increment::ComplexF64
    to_increment::ComplexF64
    from_history::ComplexF64
    to_history::ComplexF64
    history_mutated::Bool
end

struct CableFrequencyScanModalResponseSamples
    sample_frequencies_hz::Vector{Float64}
    modal_response_samples::Matrix{ComplexF64}
    response_kind::Symbol
    mode_count::Int
    frequency_count::Int
end

mutable struct FrequencyDependentLineSampleRuntimeState
    sample_frequencies_hz::Vector{Float64}
    sample_rows::Vector{Vector{LineFrequencyPoint}}
    line_length::Float64
    current_fit::LineFrequencySampleFitResult
    runtime_state::FrequencyDependentLineRuntimeState
    previous_frequency_hz::Float64
    current_frequency_hz::Float64
    previous_interpolation_weight::Float64
    current_interpolation_weight::Float64
    frequency_update_count::Int
end

mutable struct FrequencyDependentLineModalSampleRuntimeState
    sample_frequencies_hz::Vector{Float64}
    sample_rows::Vector{Vector{LineFrequencyPoint}}
    modal_yz_matrices::Vector{Matrix{ComplexF64}}
    line_length::Float64
    mode_order_state::LineModeUnwindState
    current_modal_solution::LineModalSolution
    sample_runtime_state::FrequencyDependentLineSampleRuntimeState
    previous_frequency_hz::Float64
    current_frequency_hz::Float64
    modal_solution_update_count::Int
end

mutable struct CableFrequencyScanRuntimeState
    scan::CableFrequencyScanLineConstants
    modal_sample_state::FrequencyDependentLineModalSampleRuntimeState
    update_count::Int
end

mutable struct FrequencyDependentLineRecursiveConvolutionState
    modal_sample_state::FrequencyDependentLineModalSampleRuntimeState
    pole_decay::Matrix{ComplexF64}
    residue::Matrix{ComplexF64}
    modal_history_state::Matrix{ComplexF64}
    previous_modal_history_state::Matrix{ComplexF64}
    convolution_modal_current::Vector{ComplexF64}
    previous_convolution_modal_current::Vector{ComplexF64}
    convolution_phase_current::Vector{ComplexF64}
    previous_convolution_phase_current::Vector{ComplexF64}
    sending_phase_current::Vector{ComplexF64}
    receiving_phase_current::Vector{ComplexF64}
    skin_effect_internal_impedance_executed::Bool
    earth_return_impedance_executed::Bool
    frequency_dependent_fitting_executed::Bool
    frequency_loop_executed::Bool
    pipe_sheath_side_effects_executed::Bool
    update_count::Int
end

mutable struct CableFrequencyScanRecursiveConvolutionState
    scan::CableFrequencyScanLineConstants
    recursive_state::FrequencyDependentLineRecursiveConvolutionState
    update_count::Int
end

function _complex_square_matrix(real_part::AbstractMatrix, imag_part::AbstractMatrix, label::AbstractString)
    size(real_part) == size(imag_part) ||
        throw(ArgumentError("$label real and imaginary matrices must have identical dimensions"))
    rows, cols = size(real_part)
    rows == cols && rows > 0 ||
        throw(ArgumentError("$label matrix must be nonempty and square"))
    matrix = complex.(Float64.(real_part), Float64.(imag_part))
    all(value -> isfinite(real(value)) && isfinite(imag(value)), matrix) ||
        throw(ArgumentError("$label matrix entries must be finite"))
    return matrix
end

function LineModalTransform(
    phase_to_modal_re::AbstractMatrix,
    phase_to_modal_im::AbstractMatrix,
    modal_to_phase_re::AbstractMatrix,
    modal_to_phase_im::AbstractMatrix,
)
    phase_to_modal = _complex_square_matrix(phase_to_modal_re, phase_to_modal_im, "phase_to_modal")
    modal_to_phase = _complex_square_matrix(modal_to_phase_re, modal_to_phase_im, "modal_to_phase")
    size(phase_to_modal) == size(modal_to_phase) ||
        throw(ArgumentError("phase_to_modal and modal_to_phase matrices must have identical dimensions"))
    return LineModalTransform(phase_to_modal, modal_to_phase)
end

LineModalTransform(phase_to_modal::AbstractMatrix{<:Complex}, modal_to_phase::AbstractMatrix{<:Complex}) =
    LineModalTransform(real.(phase_to_modal), imag.(phase_to_modal), real.(modal_to_phase), imag.(modal_to_phase))

function LineModeUnwindState(mode_count::Integer)
    count = Int(mode_count)
    0 < count <= 15 || throw(ArgumentError("line mode count must be between 1 and 15"))
    return LineModeUnwindState(
        collect(1:count),
        collect(1:count),
        zeros(Float64, 4, count),
        zeros(Float64, 4, count),
        zeros(Float64, 4, count),
        0,
    )
end

function _line_modal_dimension(transform::LineModalTransform)::Int
    n = size(transform.phase_to_modal, 1)
    size(transform.phase_to_modal, 2) == n &&
        size(transform.modal_to_phase, 1) == n &&
        size(transform.modal_to_phase, 2) == n ||
        throw(ArgumentError("line modal transform matrices must have identical dimensions"))
    return n
end

function _checked_line_mode_metrics(metrics::AbstractMatrix)
    rows, cols = size(metrics)
    rows == 4 && 0 < cols <= 15 ||
        throw(ArgumentError("line mode metrics must be a 4 x mode_count matrix with at most 15 modes"))
    checked = Matrix{Float64}(undef, rows, cols)
    for col in 1:cols, row in 1:rows
        value = Float64(metrics[row, col])
        isfinite(value) || throw(ArgumentError("line mode metrics must be finite"))
        checked[row, col] = value
    end
    return checked
end

function _check_line_mode_state!(state::LineModeUnwindState, mode_count::Int)
    length(state.sequence) == mode_count &&
        length(state.previous_sequence) == mode_count &&
        size(state.current_metrics) == (4, mode_count) &&
        size(state.previous_metrics) == (4, mode_count) &&
        size(state.second_previous_metrics) == (4, mode_count) ||
        throw(ArgumentError("line mode unwind state dimensions must match the mode metrics"))
    if all(==(0), state.sequence)
        state.sequence .= 1:mode_count
    end
    sort(state.sequence) == collect(1:mode_count) ||
        throw(ArgumentError("line mode sequence must be a permutation of 1:mode_count"))
    return state
end

function _reset_line_mode_unwind_state!(state::LineModeUnwindState)
    mode_count = length(state.sequence)
    state.sequence .= 1:mode_count
    state.previous_sequence .= 1:mode_count
    fill!(state.current_metrics, 0.0)
    fill!(state.previous_metrics, 0.0)
    fill!(state.second_previous_metrics, 0.0)
    return state
end

function _line_ordered_mode_metrics(metrics::AbstractMatrix{Float64}, sequence::AbstractVector{Int})
    ordered = Matrix{Float64}(undef, size(metrics))
    for (target_col, source_col) in pairs(sequence)
        ordered[:, target_col] .= metrics[:, source_col]
    end
    return ordered
end

_line_crossed_metric(current::AbstractMatrix{Float64}, previous::AbstractMatrix{Float64}, row::Int, left::Int, right::Int) =
    (previous[row, left] - previous[row, right]) * (current[row, left] - current[row, right]) <= 0.0

function _line_mode_continuity_kept(
    current::AbstractMatrix{Float64},
    previous::AbstractMatrix{Float64},
    second_previous::AbstractMatrix{Float64},
    row::Int,
    left::Int,
    right::Int,
)
    left_den = previous[row, left] - second_previous[row, left]
    right_den = previous[row, right] - second_previous[row, right]
    left_den == 0.0 && (left_den = 1.0e-12)
    right_den == 0.0 && (right_den = 1.0e-12)
    left_same = (current[row, left] - previous[row, left]) / left_den
    right_same = (current[row, right] - previous[row, right]) / right_den
    left_swapped = (current[row, right] - previous[row, left]) / left_den
    right_swapped = (current[row, left] - previous[row, right]) / right_den
    return abs(left_same - 1.0) < abs(left_swapped - 1.0) &&
        abs(right_same - 1.0) < abs(right_swapped - 1.0)
end

function _smooth_line_mode_angle_magnitude!(
    sequence::Vector{Int},
    raw_metrics::Matrix{Float64},
    previous_metrics::Matrix{Float64},
    second_previous_metrics::Matrix{Float64},
)
    mode_count = length(sequence)
    iteration = 0
    while iteration <= mode_count * (mode_count - 1)
        current = _line_ordered_mode_metrics(raw_metrics, sequence)
        changed = false
        for left in 1:(mode_count - 1), right in (left + 1):mode_count
            _line_crossed_metric(current, previous_metrics, 2, left, right) ||
                _line_crossed_metric(current, previous_metrics, 1, left, right) ||
                continue
            _line_mode_continuity_kept(current, previous_metrics, second_previous_metrics, 2, left, right) &&
                continue
            _line_mode_continuity_kept(current, previous_metrics, second_previous_metrics, 1, left, right) &&
                continue
            sequence[left], sequence[right] = sequence[right], sequence[left]
            changed = true
            current = _line_ordered_mode_metrics(raw_metrics, sequence)
        end
        iteration += 1
        changed || break
    end
    return sequence
end

function _line_mode_db_rejects_swap(
    current::AbstractMatrix{Float64},
    previous::AbstractMatrix{Float64},
    second_previous::AbstractMatrix{Float64},
    left::Int,
    right::Int,
)
    left_den = previous[4, left] - second_previous[4, left]
    right_den = previous[4, right] - second_previous[4, right]
    left_den == 0.0 && (left_den = 1.0e-12)
    right_den == 0.0 && (right_den = 1.0e-12)
    left_swapped = (current[4, right] - previous[4, left]) / left_den
    right_swapped = (current[4, left] - previous[4, right]) / right_den
    return abs(left_swapped - 1.0) > 1.0 && abs(right_swapped - 1.0) > 1.0
end

function _smooth_line_mode_velocity_db!(
    sequence::Vector{Int},
    raw_metrics::Matrix{Float64},
    previous_metrics::Matrix{Float64},
    second_previous_metrics::Matrix{Float64},
)
    mode_count = length(sequence)
    iteration = 0
    while iteration <= mode_count * (mode_count - 1)
        current = _line_ordered_mode_metrics(raw_metrics, sequence)
        start_mode = 0
        for mode in 1:mode_count
            denominator = previous_metrics[3, mode] - second_previous_metrics[3, mode]
            denominator == 0.0 && (denominator = 1.0e-12)
            slope = (current[3, mode] - previous_metrics[3, mode]) / denominator
            if abs(slope - 1.0) > 0.5
                start_mode = mode
                break
            end
        end
        start_mode == 0 && break
        changed = false
        for left in start_mode:(mode_count - 1), right in (left + 1):mode_count
            _line_crossed_metric(current, previous_metrics, 3, left, right) ||
                _line_crossed_metric(current, previous_metrics, 4, left, right) ||
                continue
            _line_mode_continuity_kept(current, previous_metrics, second_previous_metrics, 3, left, right) &&
                continue
            _line_mode_db_rejects_swap(current, previous_metrics, second_previous_metrics, left, right) &&
                continue
            sequence[left], sequence[right] = sequence[right], sequence[left]
            changed = true
            current = _line_ordered_mode_metrics(raw_metrics, sequence)
        end
        iteration += 1
        changed || break
    end
    return sequence
end

function line_mode_unwind!(
    state::LineModeUnwindState,
    metrics::AbstractMatrix;
    ntol::Integer = 2,
    nrp::Integer = 0,
)
    raw_metrics = _checked_line_mode_metrics(metrics)
    mode_count = size(raw_metrics, 2)
    _check_line_mode_state!(state, mode_count)
    ntol_int = Int(ntol)
    nrp_int = Int(nrp)
    if ntol_int <= 1 && nrp_int != 0
        _reset_line_mode_unwind_state!(state)
    end

    before_sequence = copy(state.sequence)
    state.previous_sequence .= state.sequence
    state.second_previous_metrics .= state.previous_metrics
    state.previous_metrics .= state.current_metrics
    state.current_metrics .= _line_ordered_mode_metrics(raw_metrics, state.sequence)

    if ntol_int >= 2
        _smooth_line_mode_angle_magnitude!(
            state.sequence,
            raw_metrics,
            state.previous_metrics,
            state.second_previous_metrics,
        )
        _smooth_line_mode_velocity_db!(
            state.sequence,
            raw_metrics,
            state.previous_metrics,
            state.second_previous_metrics,
        )
        state.current_metrics .= _line_ordered_mode_metrics(raw_metrics, state.sequence)
    end
    state.update_count += 1
    sequence_mutated = state.sequence != before_sequence
    return (
        source = :line_mode_unwind,
        outcome = :state_mutation,
        fortran_files = (:OVER44_FOR, :OVER47_FOR),
        fortran_routines = (:MODAL, :UNWIND),
        common_regions = (:LINEMODEL,),
        sequence = copy(state.sequence),
        previous_sequence = before_sequence,
        ordered_metrics = copy(state.current_metrics),
        previous_metrics = copy(state.previous_metrics),
        second_previous_metrics = copy(state.second_previous_metrics),
        update_count = state.update_count,
        sequence_mutated = sequence_mutated,
        mode_order_state_mutated = sequence_mutated || state.update_count > 0,
        legacy_fortran_in_loop = false,
        full_modal_eigenvector_generation_executed = false,
        deferred_calls = (:DCEIGN_bulk_oracle, :frequency_dependent_fitting, :line_timestep_bulk_oracle),
    )
end

function line_mode_unwind(metrics::AbstractMatrix; ntol::Integer = 2, nrp::Integer = 0)
    raw_metrics = _checked_line_mode_metrics(metrics)
    state = LineModeUnwindState(size(raw_metrics, 2))
    return line_mode_unwind!(state, raw_metrics; ntol = ntol, nrp = nrp)
end

function _checked_line_frequency_points(
    frequency_points::AbstractVector{LineFrequencyPoint},
    mode_count::Int,
)
    length(frequency_points) == mode_count ||
        throw(ArgumentError("frequency point count must match the modal transform dimension"))
    points = LineFrequencyPoint[]
    sizehint!(points, mode_count)
    for point in frequency_points
        isfinite(point.frequency_hz) && point.frequency_hz >= 0.0 ||
            throw(ArgumentError("line frequency values must be finite and nonnegative"))
        zc = ComplexF64(point.characteristic_impedance)
        zc != 0.0 + 0.0im ||
            throw(ArgumentError("line characteristic impedance must be nonzero"))
        all(isfinite, (real(zc), imag(zc))) ||
            throw(ArgumentError("line characteristic impedance values must be finite"))
        gamma = ComplexF64(point.propagation_constant)
        all(isfinite, (real(gamma), imag(gamma))) ||
            throw(ArgumentError("line propagation constants must be finite"))
        factor = ComplexF64(point.propagation_factor)
        all(isfinite, (real(factor), imag(factor))) ||
            throw(ArgumentError("line propagation factors must be finite"))
        push!(points, LineFrequencyPoint(point.frequency_hz, zc, gamma, factor))
    end
    return points
end

function _line_frequency_row_tolerance(frequency_hz::Float64)
    return max(1.0e-12, 64.0 * eps(Float64) * max(1.0, abs(frequency_hz)))
end

function _checked_line_frequency_sample_rows_with_order(sample_rows::AbstractVector)
    isempty(sample_rows) && throw(ArgumentError("line frequency sample rows must be nonempty"))
    first_row = sample_rows[1]
    first_row isa AbstractVector{LineFrequencyPoint} ||
        throw(ArgumentError("line frequency sample rows must contain LineFrequencyPoint vectors"))
    mode_count = length(first_row)
    mode_count > 0 || throw(ArgumentError("line frequency sample rows must contain at least one mode"))
    frequencies = Vector{Float64}(undef, length(sample_rows))
    rows = Vector{Vector{LineFrequencyPoint}}(undef, length(sample_rows))
    for row_index in eachindex(sample_rows)
        row = sample_rows[row_index]
        row isa AbstractVector{LineFrequencyPoint} ||
            throw(ArgumentError("line frequency sample rows must contain LineFrequencyPoint vectors"))
        points = _checked_line_frequency_points(row, mode_count)
        frequency = points[1].frequency_hz
        tolerance = _line_frequency_row_tolerance(frequency)
        for point in points
            abs(point.frequency_hz - frequency) <= tolerance ||
                throw(ArgumentError("all modes in a line frequency sample row must share one frequency"))
        end
        frequencies[row_index] = frequency
        rows[row_index] = points
    end
    order = sortperm(frequencies)
    sorted_frequencies = frequencies[order]
    sorted_rows = rows[order]
    for idx in 2:length(sorted_frequencies)
        abs(sorted_frequencies[idx] - sorted_frequencies[idx - 1]) >
            _line_frequency_row_tolerance(sorted_frequencies[idx]) ||
            throw(ArgumentError("line frequency sample rows must not contain duplicate frequencies"))
    end
    return sorted_frequencies, sorted_rows, mode_count, order
end

function _checked_line_frequency_sample_rows(sample_rows::AbstractVector)
    frequencies, rows, mode_count, _ = _checked_line_frequency_sample_rows_with_order(sample_rows)
    return frequencies, rows, mode_count
end

function _line_frequency_interpolation_weight(lower_frequency::Float64, upper_frequency::Float64, target_frequency::Float64)
    if lower_frequency > 0.0 && target_frequency > 0.0 && upper_frequency > 0.0
        return (log(target_frequency) - log(lower_frequency)) / (log(upper_frequency) - log(lower_frequency))
    end
    return (target_frequency - lower_frequency) / (upper_frequency - lower_frequency)
end

_line_complex_lerp(lower::ComplexF64, upper::ComplexF64, weight::Float64) =
    lower + weight * (upper - lower)

function _checked_line_target_frequency(target_frequency_hz::Real)
    target_frequency = Float64(target_frequency_hz)
    isfinite(target_frequency) && target_frequency >= 0.0 ||
        throw(ArgumentError("target_frequency_hz must be finite and nonnegative"))
    return target_frequency
end

function _checked_line_length(line_length::Real)
    line_length_value = Float64(line_length)
    isfinite(line_length_value) && line_length_value >= 0.0 ||
        throw(ArgumentError("line_length must be finite and nonnegative"))
    return line_length_value
end

function _line_frequency_sample_fit_from_sorted(
    frequencies::AbstractVector{Float64},
    rows,
    mode_count::Int,
    target_frequency::Float64,
    line_length_value::Float64,
)
    exact_index = findfirst(
        frequency -> abs(frequency - target_frequency) <= _line_frequency_row_tolerance(frequency),
        frequencies,
    )
    if exact_index !== nothing
        points = LineFrequencyPoint[
            LineFrequencyPoint(
                target_frequency,
                point.characteristic_impedance,
                point.propagation_constant,
                exp(-point.propagation_constant * line_length_value),
            )
            for point in rows[exact_index]
        ]
        return LineFrequencySampleFitResult(
            target_frequency,
            points,
            frequencies[exact_index],
            frequencies[exact_index],
            0.0,
            true,
            line_length_value,
            mode_count,
        )
    end
    first(frequencies) < target_frequency < last(frequencies) ||
        throw(ArgumentError("target_frequency_hz must be within the sampled frequency range"))
    upper_index = findfirst(frequency -> frequency > target_frequency, frequencies)
    lower_index = upper_index - 1
    lower_frequency = frequencies[lower_index]
    upper_frequency = frequencies[upper_index]
    weight = _line_frequency_interpolation_weight(lower_frequency, upper_frequency, target_frequency)
    points = LineFrequencyPoint[]
    sizehint!(points, mode_count)
    for mode in 1:mode_count
        lower_point = rows[lower_index][mode]
        upper_point = rows[upper_index][mode]
        characteristic_impedance = _line_complex_lerp(
            lower_point.characteristic_impedance,
            upper_point.characteristic_impedance,
            weight,
        )
        propagation_constant = _line_complex_lerp(
            lower_point.propagation_constant,
            upper_point.propagation_constant,
            weight,
        )
        push!(
            points,
            LineFrequencyPoint(
                target_frequency,
                characteristic_impedance,
                propagation_constant,
                exp(-propagation_constant * line_length_value),
            ),
        )
    end
    return LineFrequencySampleFitResult(
        target_frequency,
        points,
        lower_frequency,
        upper_frequency,
        weight,
        false,
        line_length_value,
        mode_count,
    )
end


function universal_machine_stator_excitation_current_solution_preview(;
    coil_conductances::AbstractVector{<:Real},
    coil_predictor_factors::AbstractVector{<:Real},
    coil_input_voltages::AbstractVector{<:Real},
    coil_history_currents::AbstractVector{<:Real},
    axis_kinds,
    stator_thevenin_matrix,
    d_axis_flux::Real,
    q_axis_flux::Real,
    half_step_interval_s::Real,
    three_phase_induction_terminal_coil::Bool = false,
)
    coil_count = length(coil_conductances)
    1 <= coil_count <= 3 ||
        throw(ArgumentError("stator/excitation current solve requires one to three active coils"))
    conductances = _machine_real_vector("coil_conductances", coil_conductances, coil_count)
    predictor_factors = _machine_real_vector("coil_predictor_factors", coil_predictor_factors, coil_count)
    input_voltages = _machine_real_vector("coil_input_voltages", coil_input_voltages, coil_count)
    history_currents = _machine_real_vector("coil_history_currents", coil_history_currents, coil_count)
    axes = _machine_symbol_vector("axis_kinds", axis_kinds, coil_count)
    all(kind -> kind in (:d_axis, :q_axis, :uncoupled), axes) ||
        throw(ArgumentError("axis_kinds entries must be :d_axis, :q_axis, or :uncoupled"))
    thevenin_matrix = _machine_real_matrix(
        "stator_thevenin_matrix",
        stator_thevenin_matrix,
        coil_count,
        coil_count,
    )

    half_step = Float64(half_step_interval_s)
    half_step > 0.0 || throw(ArgumentError("half_step_interval_s must be positive"))
    d_flux = Float64(d_axis_flux)
    q_flux = Float64(q_axis_flux)
    isfinite(d_flux) && isfinite(q_flux) ||
        throw(ArgumentError("stator/excitation flux inputs must be finite"))

    coil_admittances = conductances .* predictor_factors
    rhs = -coil_admittances .* input_voltages .+ history_currents
    for index in 1:coil_count
        if axes[index] == :d_axis
            rhs[index] -= d_flux * coil_admittances[index] / half_step
        elseif axes[index] == :q_axis
            rhs[index] -= q_flux * coil_admittances[index] / half_step
        end
    end
    if three_phase_induction_terminal_coil
        rhs[coil_count] = -coil_admittances[coil_count] * input_voltages[coil_count] +
                          history_currents[coil_count]
    end

    matrix = zeros(Float64, coil_count, coil_count)
    for row in 1:coil_count, col in 1:coil_count
        matrix[row, col] = coil_admittances[row] * thevenin_matrix[row, col]
    end
    for index in 1:coil_count
        matrix[index, index] += 1.0
    end

    currents = _universal_machine_active_coil_solution(matrix, rhs)
    return (
        source = :universal_machine_stator_excitation_current_solution,
        axis_kinds = axes,
        coil_admittances = coil_admittances,
        solution_matrix = matrix,
        right_hand_side = rhs,
        current_values = currents,
        three_phase_induction_terminal_coil = three_phase_induction_terminal_coil,
        stator_excitation_solution_executed = true,
        stator_excitation_solution_mutated = false,
        complete_machine_equation_solution = false,
        full_machine_equation_solution = false,
        deferred_effects = (
            :special_direct_machine_coupling,
            :history_vector_update,
            :flux_saturation_iteration,
            :mechanical_speed_iteration,
            :multi_machine_orchestration,
            :tacs_transfer,
            :output_report_mutation,
        ),
    )
end

function universal_machine_stator_excitation_current_solution!(
    state::UniversalMachineStatorExcitationCurrentState;
    kwargs...,
)
    preview = universal_machine_stator_excitation_current_solution_preview(; kwargs...)
    empty!(state.current_values)
    empty!(state.right_hand_side)
    empty!(state.axis_kinds)
    append!(state.current_values, preview.current_values)
    append!(state.right_hand_side, preview.right_hand_side)
    append!(state.axis_kinds, preview.axis_kinds)
    state.solution_matrix = copy(preview.solution_matrix)
    state.stator_excitation_solution_mutated = true
    return merge(
        preview,
        (
            stator_excitation_solution_mutated = state.stator_excitation_solution_mutated,
            current_values = copy(state.current_values),
            solution_matrix = copy(state.solution_matrix),
            right_hand_side = copy(state.right_hand_side),
            axis_kinds = copy(state.axis_kinds),
        ),
    )
end

function universal_machine_rotor_current_solution_preview(;
    coil_conductances::AbstractVector{<:Real},
    coil_predictor_factors::AbstractVector{<:Real},
    coil_input_voltages::AbstractVector{<:Real},
    coil_history_currents::AbstractVector{<:Real},
    coil_resistances::AbstractVector{<:Real},
    rotor_thevenin_matrix,
    d_axis_flux::Real,
    q_axis_flux::Real,
    electrical_speed_rad_s::Real,
    half_step_interval_s::Real,
    direct_axis_coupling_scale::Real = 1.0,
    quadrature_axis_extra_impedance::Real = 0.0,
    cage_coupling_fraction::Real = 0.0,
    rotor_damping_resistance::Real = 0.0,
)
    conductances = _machine_real_vector("coil_conductances", coil_conductances, 3)
    predictor_factors = _machine_real_vector("coil_predictor_factors", coil_predictor_factors, 3)
    input_voltages = _machine_real_vector("coil_input_voltages", coil_input_voltages, 3)
    history_currents = _machine_real_vector("coil_history_currents", coil_history_currents, 3)
    resistances = _machine_real_vector("coil_resistances", coil_resistances, 3)
    thevenin_matrix = _machine_real_matrix("rotor_thevenin_matrix", rotor_thevenin_matrix, 3, 3)

    half_step = Float64(half_step_interval_s)
    half_step > 0.0 || throw(ArgumentError("half_step_interval_s must be positive"))
    speed = Float64(electrical_speed_rad_s)
    d_flux = Float64(d_axis_flux)
    q_flux = Float64(q_axis_flux)
    direct_scale = Float64(direct_axis_coupling_scale)
    q_extra = Float64(quadrature_axis_extra_impedance)
    cage_fraction = Float64(cage_coupling_fraction)
    damping_resistance = Float64(rotor_damping_resistance)
    scalars = (speed, d_flux, q_flux, direct_scale, q_extra, cage_fraction, damping_resistance)
    all(isfinite, scalars) ||
        throw(ArgumentError("rotor-current equation scalar inputs must be finite"))

    coil_admittances = conductances .* predictor_factors
    quadrature_impedance = q_extra + cage_fraction * damping_resistance
    matrix = zeros(Float64, 3, 3)
    for row in 1:3, col in 1:3
        matrix[row, col] = coil_admittances[row] * thevenin_matrix[row, col]
    end
    matrix[1, 1] += 1.0
    matrix[2, 2] += 1.0
    matrix[3, 3] += coil_admittances[3] * quadrature_impedance + 1.0
    matrix[2, 3] += direct_scale * speed * coil_admittances[2] * resistances[3]
    matrix[3, 2] -= direct_scale * speed * coil_admittances[3] * resistances[2]

    rhs = -coil_admittances .* input_voltages .+ history_currents
    rhs[2] -= coil_admittances[2] * (d_flux / half_step + q_flux * speed)
    rhs[3] += coil_admittances[3] * (d_flux * speed - q_flux / half_step)
    rhs[3] -= coil_admittances[3] * d_flux * cage_fraction / half_step

    currents = _universal_machine_three_phase_solution(matrix, rhs)
    return (
        source = :universal_machine_rotor_current_solution,
        coil_admittances = coil_admittances,
        solution_matrix = matrix,
        right_hand_side = rhs,
        current_values = currents,
        d_axis_current = currents[2],
        q_axis_current = currents[3],
        zero_sequence_current = currents[1],
        quadrature_axis_impedance = quadrature_impedance,
        rotor_current_solution_executed = true,
        rotor_current_solution_mutated = false,
        complete_machine_equation_solution = false,
        full_machine_equation_solution = false,
        deferred_effects = (
            :stator_excitation_current_solution,
            :flux_saturation_iteration,
            :mechanical_speed_iteration,
            :multi_machine_orchestration,
            :tacs_transfer,
            :output_report_mutation,
        ),
    )
end

function universal_machine_rotor_current_solution!(
    state::UniversalMachineRotorCurrentState;
    kwargs...,
)
    preview = universal_machine_rotor_current_solution_preview(; kwargs...)
    empty!(state.current_values)
    empty!(state.right_hand_side)
    append!(state.current_values, preview.current_values)
    append!(state.right_hand_side, preview.right_hand_side)
    state.solution_matrix = copy(preview.solution_matrix)
    state.rotor_current_solution_mutated = true
    return merge(
        preview,
        (
            rotor_current_solution_mutated = state.rotor_current_solution_mutated,
            current_values = copy(state.current_values),
            solution_matrix = copy(state.solution_matrix),
            right_hand_side = copy(state.right_hand_side),
        ),
    )
end

function universal_machine_direct_coupling_preview(;
    current_values::AbstractVector{<:Real},
    voltage_coupling_fraction::Real = 0.0,
    resistance_coupling_fraction::Real = 0.0,
    quadrature_axis_feedback::Real = 0.0,
    direct_axis_feedback::Real = 0.0,
    transfer_current_enabled::Bool = true,
)
    currents = _machine_float_vector(current_values)
    length(currents) >= 5 ||
        throw(ArgumentError("direct-machine coupling requires at least five universal-machine current values"))
    voltage_fraction = Float64(voltage_coupling_fraction)
    resistance_fraction = Float64(resistance_coupling_fraction)
    q_feedback = Float64(quadrature_axis_feedback)
    d_feedback = Float64(direct_axis_feedback)
    all(isfinite, (voltage_fraction, resistance_fraction, q_feedback, d_feedback)) ||
        throw(ArgumentError("direct-machine coupling scalar inputs must be finite"))

    coupling_fraction = voltage_fraction + resistance_fraction
    coupled = copy(currents)
    if coupling_fraction != 0.0
        denominator = 1.0 - q_feedback * d_feedback
        denominator != 0.0 ||
            throw(ArgumentError("direct-machine coupling feedback denominator is singular"))
        q_axis_current = currents[3]
        direct_transfer_current = currents[4]
        coupled[3] = (q_axis_current + q_feedback * direct_transfer_current) / denominator
        coupled[4] = (d_feedback * q_axis_current + direct_transfer_current) / denominator
        if transfer_current_enabled
            coupled[5] = coupling_fraction * coupled[3] - resistance_fraction * coupled[4]
        end
    end

    return (
        source = :universal_machine_direct_coupling,
        current_values = coupled,
        coupling_fraction = coupling_fraction,
        voltage_coupling_fraction = voltage_fraction,
        resistance_coupling_fraction = resistance_fraction,
        transfer_current_enabled = transfer_current_enabled,
        direct_coupling_executed = true,
        direct_coupling_mutated = false,
        complete_machine_equation_solution = false,
        full_machine_equation_solution = false,
        deferred_effects = (
            :history_vector_update,
            :flux_saturation_iteration,
            :mechanical_speed_iteration,
            :multi_machine_orchestration,
            :tacs_transfer,
            :output_report_mutation,
        ),
    )
end

function universal_machine_direct_coupling!(
    state::UniversalMachineDirectCouplingState;
    kwargs...,
)
    preview = universal_machine_direct_coupling_preview(; kwargs...)
    empty!(state.current_values)
    append!(state.current_values, preview.current_values)
    state.coupling_fraction = preview.coupling_fraction
    state.voltage_coupling_fraction = preview.voltage_coupling_fraction
    state.resistance_coupling_fraction = preview.resistance_coupling_fraction
    state.direct_coupling_mutated = true
    return merge(
        preview,
        (
            direct_coupling_mutated = state.direct_coupling_mutated,
            current_values = copy(state.current_values),
        ),
    )
end

function universal_machine_history_current_update_preview(;
    current_values::AbstractVector{<:Real},
    input_voltages::AbstractVector{<:Real},
    coil_conductances::AbstractVector{<:Real},
    coil_predictor_factors::AbstractVector{<:Real},
    coil_resistances::AbstractVector{<:Real},
    axis_kinds,
    rotor_thevenin_matrix,
    stator_thevenin_matrix,
    d_axis_flux::Real,
    q_axis_flux::Real,
    electrical_speed_rad_s::Real,
    half_step_interval_s::Real,
    three_phase_induction_terminal_coil::Bool = false,
    quadrature_axis_extra_impedance::Real = 0.0,
    direct_machine_coupling_fraction::Real = 0.0,
    direct_machine_voltage_coupling_fraction::Real = 0.0,
    direct_machine_resistance_coupling_fraction::Real = 0.0,
    direct_machine_rotor_resistance::Real = 0.0,
    speed_coupling_scale::Real = 1.0,
    direct_machine_stator_history_coupling_enabled::Bool = true,
)
    coil_count = length(current_values)
    4 <= coil_count <= 6 ||
        throw(ArgumentError("history-current update requires three rotor coils plus one to three stator/excitation coils"))
    stator_count = coil_count - 3
    currents = _machine_real_vector("current_values", current_values, coil_count)
    voltages = _machine_real_vector("input_voltages", input_voltages, coil_count)
    conductances = _machine_real_vector("coil_conductances", coil_conductances, coil_count)
    predictor_factors = _machine_real_vector("coil_predictor_factors", coil_predictor_factors, coil_count)
    resistances = _machine_real_vector("coil_resistances", coil_resistances, coil_count)
    axes = _machine_symbol_vector("axis_kinds", axis_kinds, stator_count)
    all(kind -> kind in (:d_axis, :q_axis, :uncoupled), axes) ||
        throw(ArgumentError("axis_kinds entries must be :d_axis, :q_axis, or :uncoupled"))
    rotor_thevenin = _machine_real_matrix("rotor_thevenin_matrix", rotor_thevenin_matrix, 3, 3)
    stator_thevenin =
        _machine_real_matrix("stator_thevenin_matrix", stator_thevenin_matrix, stator_count, stator_count)

    half_step = Float64(half_step_interval_s)
    half_step > 0.0 || throw(ArgumentError("half_step_interval_s must be positive"))
    d_flux = Float64(d_axis_flux)
    q_flux = Float64(q_axis_flux)
    speed = Float64(electrical_speed_rad_s)
    extra_impedance = Float64(quadrature_axis_extra_impedance)
    existing_coupling_fraction = Float64(direct_machine_coupling_fraction)
    voltage_coupling_fraction = Float64(direct_machine_voltage_coupling_fraction)
    resistance_coupling_fraction = Float64(direct_machine_resistance_coupling_fraction)
    damping_resistance = Float64(direct_machine_rotor_resistance)
    coupling_scale = Float64(speed_coupling_scale)
    coupling_fraction =
        existing_coupling_fraction + voltage_coupling_fraction + resistance_coupling_fraction
    all(
        isfinite,
        (
            d_flux,
            q_flux,
            speed,
            extra_impedance,
            existing_coupling_fraction,
            voltage_coupling_fraction,
            resistance_coupling_fraction,
            damping_resistance,
            coupling_scale,
        ),
    ) || throw(ArgumentError("history-current scalar inputs must be finite"))
    if voltage_coupling_fraction != 0.0 || resistance_coupling_fraction != 0.0
        coil_count >= 5 ||
            throw(ArgumentError("direct-machine history coupling requires at least five coils"))
    end

    coil_admittances = conductances .* predictor_factors
    history = zeros(Float64, coil_count)
    rotor_matrix = zeros(Float64, 3, 3)
    for row in 1:3, col in 1:3
        rotor_matrix[row, col] = -coil_admittances[row] * rotor_thevenin[row, col]
    end

    d_axis_admittance = coil_admittances[2]
    q_axis_admittance = coil_admittances[3]
    rotor_matrix[3, 3] -= (extra_impedance + coupling_fraction * damping_resistance) * q_axis_admittance
    for index in 1:3
        rotor_matrix[index, index] +=
            (conductances[index] * resistances[index] / half_step - 1.0) * predictor_factors[index]
    end
    rotor_matrix[2, 3] -= coupling_scale * speed * d_axis_admittance * resistances[3]
    rotor_matrix[3, 2] += coupling_scale * speed * q_axis_admittance * resistances[2]

    for row in 1:3
        value = 0.0
        for col in 1:3
            value += rotor_matrix[row, col] * currents[col]
        end
        history[row] = value - coil_admittances[row] * voltages[row]
    end
    history[2] += d_axis_admittance * (d_flux / half_step - q_flux * speed)
    history[3] += q_axis_admittance * (d_flux * speed + q_flux / half_step)
    history[3] += d_flux * coupling_fraction * q_axis_admittance / half_step
    if voltage_coupling_fraction != 0.0
        history[3] += voltage_coupling_fraction * q_axis_admittance * extra_impedance * currents[5]
    end
    if resistance_coupling_fraction != 0.0
        resistance_history_impedance = extra_impedance + damping_resistance - resistances[5] / half_step
        history[3] +=
            resistance_coupling_fraction * q_axis_admittance * resistance_history_impedance *
            currents[5]
    end

    for stator_index in 1:stator_count
        coil_index = stator_index + 3
        admittance = coil_admittances[coil_index]
        history[coil_index] =
            -admittance * voltages[coil_index] +
            (-predictor_factors[coil_index] + admittance * resistances[coil_index] / half_step) *
            currents[coil_index]
        if axes[stator_index] == :d_axis
            history[coil_index] += d_flux * admittance / half_step
        elseif axes[stator_index] == :q_axis
            history[coil_index] += q_flux * admittance / half_step
        end
    end
    if direct_machine_stator_history_coupling_enabled &&
       (voltage_coupling_fraction != 0.0 || resistance_coupling_fraction != 0.0)
        direct_axis_stator_admittance = coil_admittances[4]
        if resistance_coupling_fraction != 0.0
            history[4] -= resistance_coupling_fraction * d_flux * direct_axis_stator_admittance / half_step
            resistance_history_impedance =
                extra_impedance + damping_resistance - resistances[5] / half_step
            history[4] +=
                resistance_coupling_fraction * direct_axis_stator_admittance *
                resistance_history_impedance * currents[3]
            history[4] -=
                resistance_coupling_fraction * (extra_impedance + damping_resistance) *
                direct_axis_stator_admittance * currents[4] / half_step
        end
        if voltage_coupling_fraction != 0.0
            history[4] +=
                voltage_coupling_fraction * direct_axis_stator_admittance * extra_impedance *
                currents[3]
            history[4] -=
                voltage_coupling_fraction * extra_impedance * direct_axis_stator_admittance *
                currents[4] / half_step
        end
    end
    if three_phase_induction_terminal_coil
        admittance = coil_admittances[coil_count]
        history[coil_count] =
            (admittance * resistances[coil_count] / half_step - predictor_factors[coil_count]) *
            currents[coil_count] - admittance * voltages[coil_count]
    end

    stator_drop = zeros(Float64, stator_count)
    for row in 1:stator_count
        value = 0.0
        for col in 1:stator_count
            value += stator_thevenin[row, col] * currents[col + 3]
        end
        stator_drop[row] = value
        history[row + 3] -= coil_admittances[row + 3] * value
    end

    return (
        source = :universal_machine_history_current_update,
        history_currents = history,
        rotor_history_matrix = rotor_matrix,
        stator_thevenin_drop = stator_drop,
        coil_admittances = coil_admittances,
        axis_kinds = axes,
        three_phase_induction_terminal_coil = three_phase_induction_terminal_coil,
        history_current_update_executed = true,
        history_current_update_mutated = false,
        complete_machine_equation_solution = false,
        full_machine_equation_solution = false,
        deferred_effects = (
            :special_direct_machine_coupling,
            :three_phase_induction_stator_branch,
            :flux_saturation_iteration,
            :mechanical_speed_iteration,
            :multi_machine_orchestration,
            :tacs_transfer,
            :output_report_mutation,
        ),
    )
end

function universal_machine_history_current_update!(
    state::UniversalMachineHistoryCurrentState;
    kwargs...,
)
    preview = universal_machine_history_current_update_preview(; kwargs...)
    empty!(state.history_currents)
    empty!(state.stator_thevenin_drop)
    empty!(state.axis_kinds)
    append!(state.history_currents, preview.history_currents)
    append!(state.stator_thevenin_drop, preview.stator_thevenin_drop)
    append!(state.axis_kinds, preview.axis_kinds)
    state.rotor_history_matrix = copy(preview.rotor_history_matrix)
    state.history_current_update_mutated = true
    return merge(
        preview,
        (
            history_current_update_mutated = state.history_current_update_mutated,
            history_currents = copy(state.history_currents),
            rotor_history_matrix = copy(state.rotor_history_matrix),
            stator_thevenin_drop = copy(state.stator_thevenin_drop),
            axis_kinds = copy(state.axis_kinds),
        ),
    )
end

function universal_machine_mechanical_iteration_preview(;
    current_values::AbstractVector{<:Real},
    d_axis_flux::Real,
    q_axis_flux::Real,
    d_axis_reactance::Real,
    q_axis_reactance::Real,
    pole_pair_count::Real,
    mechanical_speed_rad_s::Real,
    previous_mechanical_speed_rad_s::Real,
    mechanical_angle_rad::Real,
    half_step_interval_s::Real,
    speed_tolerance::Real,
    torque_normalization::Real = 1.0,
    normalize_generated_torque::Bool = false,
    tacs_torque_controlled::Bool = false,
    mechanical_torque_input::Real = 0.0,
    damping_coefficient::Real = 0.0,
    inertia::Real = 1.0,
    speed_history::Real = 0.0,
    mechanical_speed_thevenin::Real = 0.0,
    generated_torque_impedance::Real = 0.0,
    shared_generated_torques::AbstractVector{<:Real} = Float64[],
    shared_torque_impedances::AbstractVector{<:Real} = Float64[],
    electrical_angle_scale::Real = 1.0,
    max_iterations::Int = 51,
    start_mode::Bool = true,
)
    length(current_values) >= 3 ||
        throw(ArgumentError("mechanical iteration requires at least d/q universal-machine currents"))
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    length(shared_generated_torques) == length(shared_torque_impedances) ||
        throw(ArgumentError("shared torque vectors must have matching lengths"))
    currents = _machine_real_vector("current_values", current_values, length(current_values))
    shared_torques =
        _machine_real_vector("shared_generated_torques", shared_generated_torques, length(shared_generated_torques))
    shared_impedances =
        _machine_real_vector("shared_torque_impedances", shared_torque_impedances, length(shared_torque_impedances))

    d_flux = Float64(d_axis_flux)
    q_flux = Float64(q_axis_flux)
    d_reactance = Float64(d_axis_reactance)
    q_reactance = Float64(q_axis_reactance)
    pole_pairs = Float64(pole_pair_count)
    speed = Float64(mechanical_speed_rad_s)
    previous_speed = Float64(previous_mechanical_speed_rad_s)
    angle = Float64(mechanical_angle_rad)
    half_step = Float64(half_step_interval_s)
    tolerance = Float64(speed_tolerance)
    normalization = Float64(torque_normalization)
    torque_input = Float64(mechanical_torque_input)
    damping = Float64(damping_coefficient)
    rotor_inertia = Float64(inertia)
    history_speed = Float64(speed_history)
    thevenin_speed = Float64(mechanical_speed_thevenin)
    torque_impedance = Float64(generated_torque_impedance)
    angle_scale = Float64(electrical_angle_scale)
    scalars = (
        d_flux,
        q_flux,
        d_reactance,
        q_reactance,
        pole_pairs,
        speed,
        previous_speed,
        angle,
        half_step,
        tolerance,
        normalization,
        torque_input,
        damping,
        rotor_inertia,
        history_speed,
        thevenin_speed,
        torque_impedance,
        angle_scale,
    )
    all(isfinite, scalars) || throw(ArgumentError("mechanical iteration scalar inputs must be finite"))
    half_step > 0.0 || throw(ArgumentError("half_step_interval_s must be positive"))
    tolerance >= 0.0 || throw(ArgumentError("speed_tolerance must be nonnegative"))
    pole_pairs != 0.0 || throw(ArgumentError("pole_pair_count must be nonzero"))
    if normalize_generated_torque
        normalization != 0.0 || throw(ArgumentError("torque_normalization must be nonzero"))
    end
    if tacs_torque_controlled
        rotor_inertia != 0.0 || throw(ArgumentError("inertia must be nonzero for torque-input iteration"))
    end

    d_axis_current = currents[2]
    q_axis_current = currents[3]
    reluctance_torque = (d_reactance - q_reactance) * d_axis_current * q_axis_current
    generated_torque =
        (reluctance_torque - d_axis_current * q_flux + q_axis_current * d_flux) * pole_pairs
    if normalize_generated_torque
        generated_torque /= normalization
    end

    predicted_speed = 2.0 * speed - previous_speed
    old_speed = speed
    solved_speed = speed
    final_speed = speed
    torque_increment = 0.0
    damping_torque = 0.0
    angle_candidate = angle
    electrical_angle = pole_pairs * angle_scale * angle
    electrical_speed = pole_pairs * predicted_speed
    speed_error = Inf
    iteration_count = 0
    converged = !start_mode
    failed = false

    if start_mode
        for iteration in 1:max_iterations
            iteration_count = iteration
            angle_candidate = angle + half_step * (old_speed + predicted_speed)
            electrical_angle = pole_pairs * angle_scale * angle_candidate
            electrical_speed = pole_pairs * predicted_speed
            if tacs_torque_controlled
                damping_torque = damping * solved_speed
                torque_increment =
                    (torque_input - generated_torque - damping_torque) * half_step / rotor_inertia
                solved_speed = history_speed + torque_increment
            else
                solved_speed = thevenin_speed - torque_impedance * generated_torque
                for index in eachindex(shared_torques)
                    solved_speed -= shared_impedances[index] * shared_torques[index]
                end
                torque_increment = 0.0
                damping_torque = 0.0
            end

            speed_error = solved_speed - predicted_speed
            if abs(speed_error) <= tolerance
                final_speed = predicted_speed
                converged = true
                break
            end
            predicted_speed = solved_speed
            final_speed = solved_speed
        end
        failed = !converged
    end

    speed_history_next = final_speed + torque_increment
    return (
        source = :universal_machine_mechanical_iteration,
        generated_torque = generated_torque,
        reluctance_torque = reluctance_torque * pole_pairs,
        d_axis_current = d_axis_current,
        q_axis_current = q_axis_current,
        predicted_speed_rad_s = predicted_speed,
        solved_speed_rad_s = solved_speed,
        mechanical_speed_rad_s = final_speed,
        previous_mechanical_speed_rad_s = speed,
        input_previous_mechanical_speed_rad_s = previous_speed,
        mechanical_angle_rad = angle_candidate,
        electrical_angle_rad = electrical_angle,
        electrical_speed_rad_s = electrical_speed,
        torque_increment = torque_increment,
        damping_torque = damping_torque,
        speed_history_next = speed_history_next,
        speed_error_rad_s = speed_error,
        iteration_count = iteration_count,
        converged = converged,
        mechanical_iteration_failed = failed,
        mechanical_iteration_executed = true,
        mechanical_iteration_mutated = false,
        complete_machine_equation_solution = false,
        full_machine_equation_solution = false,
        deferred_effects = (
            :flux_saturation_iteration,
            :multi_machine_orchestration,
            :tacs_transfer,
            :output_report_mutation,
            :full_machine_equation_solution,
        ),
    )
end

function universal_machine_mechanical_iteration!(
    state::UniversalMachineMechanicalIterationState;
    kwargs...,
)
    preview = universal_machine_mechanical_iteration_preview(; kwargs...)
    state.previous_mechanical_speed_rad_s = preview.previous_mechanical_speed_rad_s
    state.mechanical_speed_rad_s = preview.mechanical_speed_rad_s
    state.mechanical_angle_rad = preview.mechanical_angle_rad
    state.predicted_speed_rad_s = preview.predicted_speed_rad_s
    state.solved_speed_rad_s = preview.solved_speed_rad_s
    state.generated_torque = preview.generated_torque
    state.torque_increment = preview.torque_increment
    state.iteration_count = preview.iteration_count
    state.converged = preview.converged
    state.mechanical_iteration_mutated = true
    return merge(
        preview,
        (
            mechanical_iteration_mutated = state.mechanical_iteration_mutated,
            mechanical_speed_rad_s = state.mechanical_speed_rad_s,
            mechanical_angle_rad = state.mechanical_angle_rad,
            predicted_speed_rad_s = state.predicted_speed_rad_s,
            solved_speed_rad_s = state.solved_speed_rad_s,
        ),
    )
end

function universal_machine_flux_saturation_preview(;
    current_values::AbstractVector{<:Real},
    d_axis_stator_current_indices::AbstractVector{<:Integer} = Int[],
    q_axis_stator_current_indices::AbstractVector{<:Integer} = Int[],
    d_axis_subtract_current_index::Int = 0,
    d_axis_unsaturated_inductance::Real,
    d_axis_saturated_inductance::Real,
    q_axis_unsaturated_inductance::Real,
    q_axis_saturated_inductance::Real,
    d_axis_saturation_flux::Real,
    q_axis_saturation_flux::Real,
    d_axis_remanent_flux::Real = 0.0,
    q_axis_remanent_flux::Real = 0.0,
    d_axis_saturation_mode::Int = 0,
    q_axis_saturation_mode::Int = 0,
)
    length(current_values) >= 3 ||
        throw(ArgumentError("flux saturation requires at least d/q universal-machine currents"))
    currents = _machine_real_vector("current_values", current_values, length(current_values))
    d_unsaturated = Float64(d_axis_unsaturated_inductance)
    d_saturated = Float64(d_axis_saturated_inductance)
    q_unsaturated = Float64(q_axis_unsaturated_inductance)
    q_saturated = Float64(q_axis_saturated_inductance)
    d_saturation_flux = Float64(d_axis_saturation_flux)
    q_saturation_flux = Float64(q_axis_saturation_flux)
    d_remanent = Float64(d_axis_remanent_flux)
    q_remanent = Float64(q_axis_remanent_flux)
    all(
        isfinite,
        (
            d_unsaturated,
            d_saturated,
            q_unsaturated,
            q_saturated,
            d_saturation_flux,
            q_saturation_flux,
            d_remanent,
            q_remanent,
        ),
    ) || throw(ArgumentError("flux saturation scalar inputs must be finite"))
    d_axis_saturation_mode >= 0 || throw(ArgumentError("d_axis_saturation_mode must be nonnegative"))
    q_axis_saturation_mode >= 0 || throw(ArgumentError("q_axis_saturation_mode must be nonnegative"))

    axis_totals = _universal_machine_axis_current_totals(
        currents,
        d_axis_stator_current_indices,
        q_axis_stator_current_indices,
        d_axis_subtract_current_index,
    )
    d_current = axis_totals.d_axis_current
    q_current = axis_totals.q_axis_current

    d_offset = 0.0
    q_offset = 0.0
    if d_axis_saturation_mode == 5
        current_magnitude = hypot(d_current, q_current)
        if current_magnitude != 0.0 && d_unsaturated != 0.0
            d_axis_projection = abs(d_current * d_saturation_flux / current_magnitude)
            q_axis_projection = abs(q_current * q_saturation_flux / current_magnitude)
            saturation_ratio = d_saturated / d_unsaturated
            d_offset = d_axis_projection + saturation_ratio * (d_remanent - d_axis_projection)
            q_offset = q_axis_projection + saturation_ratio * (q_remanent - q_axis_projection)
        end
    end

    function axis_flux(
        current::Float64,
        unsaturated::Float64,
        saturated::Float64,
        saturation_flux::Float64,
        remanent::Float64,
        saturation_mode::Int,
        offset::Float64,
    )
        if saturation_mode != 0 && unsaturated != 0.0
            threshold = saturation_flux / unsaturated
            if current > threshold || current < -threshold
                return (current >= 0.0 ? offset : -offset) + saturated * current
            end
        end
        flux = unsaturated * current + remanent
        current < 0.0 && (flux -= 2.0 * remanent)
        return flux
    end

    d_flux = axis_flux(
        d_current,
        d_unsaturated,
        d_saturated,
        d_saturation_flux,
        d_remanent,
        d_axis_saturation_mode,
        d_offset,
    )
    q_flux = axis_flux(
        q_current,
        q_unsaturated,
        q_saturated,
        q_saturation_flux,
        q_remanent,
        q_axis_saturation_mode,
        q_offset,
    )

    return (
        source = :universal_machine_flux_saturation,
        d_axis_current = d_current,
        q_axis_current = q_current,
        d_axis_flux = d_flux,
        q_axis_flux = q_flux,
        d_axis_saturation_offset = d_offset,
        q_axis_saturation_offset = q_offset,
        d_axis_saturated = d_axis_saturation_mode != 0 &&
                           d_unsaturated != 0.0 &&
                           (d_current > d_saturation_flux / d_unsaturated ||
                            d_current < -d_saturation_flux / d_unsaturated),
        q_axis_saturated = q_axis_saturation_mode != 0 &&
                           q_unsaturated != 0.0 &&
                           (q_current > q_saturation_flux / q_unsaturated ||
                            q_current < -q_saturation_flux / q_unsaturated),
        total_saturation_mode = d_axis_saturation_mode == 5,
        flux_saturation_executed = true,
        flux_saturation_mutated = false,
        complete_machine_equation_solution = false,
        full_machine_equation_solution = false,
        deferred_effects = (
            :flux_saturation_jacobian_iteration,
            :multi_machine_orchestration,
            :tacs_transfer,
            :output_report_mutation,
            :full_machine_equation_solution,
        ),
    )
end

function universal_machine_flux_saturation!(
    state::UniversalMachineFluxSaturationState;
    kwargs...,
)
    preview = universal_machine_flux_saturation_preview(; kwargs...)
    state.d_axis_flux = preview.d_axis_flux
    state.q_axis_flux = preview.q_axis_flux
    state.d_axis_current = preview.d_axis_current
    state.q_axis_current = preview.q_axis_current
    state.d_axis_saturation_offset = preview.d_axis_saturation_offset
    state.q_axis_saturation_offset = preview.q_axis_saturation_offset
    state.flux_saturation_mutated = true
    return merge(
        preview,
        (
            flux_saturation_mutated = state.flux_saturation_mutated,
            d_axis_flux = state.d_axis_flux,
            q_axis_flux = state.q_axis_flux,
            d_axis_current = state.d_axis_current,
            q_axis_current = state.q_axis_current,
            d_axis_saturation_offset = state.d_axis_saturation_offset,
            q_axis_saturation_offset = state.q_axis_saturation_offset,
        ),
    )
end

function _universal_machine_current_solution_from_flux(;
    coil_conductances::Vector{Float64},
    coil_predictor_factors::Vector{Float64},
    coil_input_voltages::Vector{Float64},
    coil_history_currents::Vector{Float64},
    coil_resistances::Vector{Float64},
    axis_kinds::Vector{Symbol},
    rotor_thevenin_matrix::Matrix{Float64},
    stator_thevenin_matrix::Matrix{Float64},
    d_axis_flux::Float64,
    q_axis_flux::Float64,
    electrical_speed_rad_s::Float64,
    half_step_interval_s::Float64,
    direct_axis_coupling_scale::Float64,
    quadrature_axis_extra_impedance::Float64,
    cage_coupling_fraction::Float64,
    rotor_damping_resistance::Float64,
    direct_machine_voltage_coupling_fraction::Float64,
    direct_machine_resistance_coupling_fraction::Float64,
    quadrature_axis_feedback::Float64,
    direct_axis_feedback::Float64,
    direct_machine_transfer_current_enabled::Bool,
)
    stator_count = length(axis_kinds)
    rotor = universal_machine_rotor_current_solution_preview(
        coil_conductances = coil_conductances[1:3],
        coil_predictor_factors = coil_predictor_factors[1:3],
        coil_input_voltages = coil_input_voltages[1:3],
        coil_history_currents = coil_history_currents[1:3],
        coil_resistances = coil_resistances[1:3],
        rotor_thevenin_matrix = rotor_thevenin_matrix,
        d_axis_flux = d_axis_flux,
        q_axis_flux = q_axis_flux,
        electrical_speed_rad_s = electrical_speed_rad_s,
        half_step_interval_s = half_step_interval_s,
        direct_axis_coupling_scale = direct_axis_coupling_scale,
        quadrature_axis_extra_impedance = quadrature_axis_extra_impedance,
        cage_coupling_fraction = cage_coupling_fraction,
        rotor_damping_resistance = rotor_damping_resistance,
    )
    stator = universal_machine_stator_excitation_current_solution_preview(
        coil_conductances = coil_conductances[4:end],
        coil_predictor_factors = coil_predictor_factors[4:end],
        coil_input_voltages = coil_input_voltages[4:end],
        coil_history_currents = coil_history_currents[4:end],
        axis_kinds = axis_kinds,
        stator_thevenin_matrix = stator_thevenin_matrix,
        d_axis_flux = d_axis_flux,
        q_axis_flux = q_axis_flux,
        half_step_interval_s = half_step_interval_s,
        three_phase_induction_terminal_coil = stator_count == 3,
    )
    currents = vcat(rotor.current_values, stator.current_values)
    if direct_machine_voltage_coupling_fraction != 0.0 ||
       direct_machine_resistance_coupling_fraction != 0.0 ||
       quadrature_axis_feedback != 0.0 ||
       direct_axis_feedback != 0.0
        coupled = universal_machine_direct_coupling_preview(
            current_values = currents,
            voltage_coupling_fraction = direct_machine_voltage_coupling_fraction,
            resistance_coupling_fraction = direct_machine_resistance_coupling_fraction,
            quadrature_axis_feedback = quadrature_axis_feedback,
            direct_axis_feedback = direct_axis_feedback,
            transfer_current_enabled = direct_machine_transfer_current_enabled,
        )
        currents = coupled.current_values
    end
    return currents
end

function _universal_machine_flux_region_solution(;
    base_d_axis_current::Float64,
    base_q_axis_current::Float64,
    d_axis_d_flux_sensitivity::Float64,
    q_axis_d_flux_sensitivity::Float64,
    d_axis_q_flux_sensitivity::Float64,
    q_axis_q_flux_sensitivity::Float64,
    d_axis_unsaturated_inductance::Float64,
    d_axis_saturated_inductance::Float64,
    q_axis_unsaturated_inductance::Float64,
    q_axis_saturated_inductance::Float64,
    d_axis_saturation_flux::Float64,
    q_axis_saturation_flux::Float64,
    d_axis_remanent_flux::Float64,
    q_axis_remanent_flux::Float64,
    d_axis_saturation_mode::Int,
    q_axis_saturation_mode::Int,
)
    saturation_request = 1 + d_axis_saturation_mode + 2 * q_axis_saturation_mode
    d_region = 0.0
    q_region = 0.0
    signed_remanent_scale = 0.0
    d_current_sign = 0.0
    q_current_sign = 0.0
    d_threshold_current = Inf
    q_threshold_current = Inf
    d_offset = 0.0
    q_offset = 0.0
    decision_code = 1
    remanent_clip_enabled = true

    final_for_region = function (d_region::Float64, q_region::Float64, remanent_clip::Bool)
        d_axis_offset = d_current_sign *
                        ((1.0 - d_region) * signed_remanent_scale * d_axis_remanent_flux +
                         d_region * d_offset)
        q_axis_offset = q_current_sign *
                        ((1.0 - q_region) * signed_remanent_scale * q_axis_remanent_flux +
                         q_region * q_offset)
        d_rhs = d_axis_d_flux_sensitivity * d_axis_offset +
                d_axis_q_flux_sensitivity * q_axis_offset +
                base_d_axis_current
        q_rhs = q_axis_q_flux_sensitivity * q_axis_offset +
                q_axis_d_flux_sensitivity * d_axis_offset +
                base_q_axis_current
        d_inductance = d_axis_unsaturated_inductance * (1.0 - d_region) +
                       d_axis_saturated_inductance * d_region
        q_inductance = q_axis_unsaturated_inductance * (1.0 - q_region) +
                       q_axis_saturated_inductance * q_region
        d_jacobian = 1.0 - d_axis_d_flux_sensitivity * d_inductance
        q_jacobian = 1.0 - q_axis_q_flux_sensitivity * q_inductance
        d_q_coupling = -d_axis_q_flux_sensitivity * q_inductance
        q_d_coupling = -q_axis_d_flux_sensitivity * d_inductance
        determinant = d_jacobian * q_jacobian - d_q_coupling * q_d_coupling
        determinant != 0.0 ||
            throw(ArgumentError("flux/current Jacobian determinant is zero"))
        d_current = (q_jacobian * d_rhs - d_q_coupling * q_rhs) / determinant
        q_current = (d_jacobian * q_rhs - q_d_coupling * d_rhs) / determinant
        d_flux = d_axis_unsaturated_inductance * (1.0 - d_region) * d_current
        q_flux = q_axis_unsaturated_inductance * (1.0 - q_region) * q_current
        if saturation_request != 1
            d_flux += d_axis_saturated_inductance * d_region * d_current
            d_flux += d_current_sign * (1.0 - d_region) * d_axis_remanent_flux
            d_flux += d_region * d_current_sign * d_offset
            q_flux += q_axis_saturated_inductance * q_region * q_current
            q_flux += q_current_sign * (1.0 - q_region) * q_axis_remanent_flux
            q_flux += q_region * q_current_sign * q_offset
        end
        if !remanent_clip
            if d_region == 0.0 && d_flux * d_flux <= d_axis_remanent_flux * d_axis_remanent_flux
                d_flux = d_current_sign * d_axis_remanent_flux
            end
            if q_region == 0.0 && q_flux * q_flux <= q_axis_remanent_flux * q_axis_remanent_flux
                q_flux = q_current_sign * q_axis_remanent_flux
            end
        end
        return (;
            d_axis_current = d_current,
            q_axis_current = q_current,
            d_axis_flux = d_flux,
            q_axis_flux = q_flux,
            d_region = d_region,
            q_region = q_region,
            d_threshold_measure = d_current_sign * d_current,
            q_threshold_measure = q_current_sign * q_current,
        )
    end

    result = final_for_region(d_region, q_region, remanent_clip_enabled)
    for iteration in 1:16
        if decision_code >= 5
            if saturation_request == 2
                if result.d_threshold_measure > d_threshold_current
                    return merge(result, (; iteration_count = iteration))
                end
                decision_code = 2
                remanent_clip_enabled = false
                d_region = 0.0
                q_region = 0.0
            elseif saturation_request == 3
                if result.q_threshold_measure > q_threshold_current
                    return merge(result, (; iteration_count = iteration))
                end
                decision_code = 2
                remanent_clip_enabled = false
                d_region = 0.0
                q_region = 0.0
            elseif saturation_request >= 4
                if saturation_request == 5
                    remanent_clip_enabled = false
                    if result.d_threshold_measure > d_threshold_current ||
                       result.q_threshold_measure > q_threshold_current
                        return merge(result, (; iteration_count = iteration))
                    end
                    decision_code = 2
                    d_region = 0.0
                    q_region = 0.0
                elseif result.d_threshold_measure <= d_threshold_current &&
                       result.q_threshold_measure <= q_threshold_current
                    decision_code = 2
                    remanent_clip_enabled = false
                    d_region = 0.0
                    q_region = 0.0
                elseif result.d_threshold_measure > d_threshold_current &&
                       result.q_threshold_measure > q_threshold_current
                    return merge(result, (; iteration_count = iteration))
                elseif result.d_threshold_measure > d_threshold_current
                    decision_code = 4
                    remanent_clip_enabled = false
                    d_region = 1.0
                    q_region = 0.0
                else
                    decision_code = 3
                    remanent_clip_enabled = false
                    d_region = 0.0
                    q_region = 1.0
                end
            else
                return merge(result, (; iteration_count = iteration))
            end
        elseif decision_code == 1
            if saturation_request == 1
                return merge(result, (; iteration_count = iteration))
            end
            d_current_sign = result.d_axis_current < 0.0 ? -1.0 : 1.0
            q_current_sign = result.q_axis_current < 0.0 ? -1.0 : 1.0
            d_axis_unsaturated_inductance != 0.0 ||
                throw(ArgumentError("d-axis unsaturated inductance must be nonzero for flux iteration"))
            q_axis_unsaturated_inductance != 0.0 ||
                throw(ArgumentError("q-axis unsaturated inductance must be nonzero for flux iteration"))
            d_threshold_current =
                (d_axis_saturation_flux - d_axis_remanent_flux) / d_axis_unsaturated_inductance
            q_threshold_current =
                (q_axis_saturation_flux - q_axis_remanent_flux) / q_axis_unsaturated_inductance
            signed_remanent_scale = 1.0
            decision_code = 5
            if d_axis_saturation_mode == 5
                current_magnitude = hypot(result.d_axis_current, result.q_axis_current)
                if current_magnitude == 0.0
                    d_projection = d_axis_saturation_flux
                    q_projection = d_projection
                else
                    d_threshold_current *= d_current_sign * result.d_axis_current / current_magnitude
                    q_threshold_current *= q_current_sign * result.q_axis_current / current_magnitude
                    d_projection = d_current_sign * result.d_axis_current * d_axis_saturation_flux /
                                   current_magnitude
                    q_projection = q_current_sign * result.q_axis_current * d_axis_saturation_flux /
                                   current_magnitude
                end
                saturation_ratio = d_axis_saturated_inductance / d_axis_unsaturated_inductance
                d_offset = d_projection + saturation_ratio * (d_axis_remanent_flux - d_projection)
                q_offset = q_projection + saturation_ratio * (q_axis_remanent_flux - q_projection)
                saturation_request = 5
            end
            if saturation_request == 2
                d_region = 1.0
                q_region = 0.0
            elseif saturation_request == 3
                d_region = 0.0
                q_region = 1.0
            elseif saturation_request >= 4
                d_region = 1.0
                q_region = 1.0
            else
                d_region = 0.0
                q_region = 0.0
            end
        elseif decision_code == 2
            return merge(result, (; iteration_count = iteration))
        elseif decision_code == 3
            if result.q_threshold_measure > q_threshold_current
                return merge(result, (; iteration_count = iteration))
            end
            decision_code = 2
            remanent_clip_enabled = false
            d_region = 0.0
            q_region = 0.0
        elseif decision_code == 4
            if result.d_threshold_measure > d_threshold_current
                return merge(result, (; iteration_count = iteration))
            end
            decision_code = 2
            remanent_clip_enabled = false
            d_region = 0.0
            q_region = 0.0
        end
        result = final_for_region(d_region, q_region, remanent_clip_enabled)
    end
    throw(ArgumentError("flux/current Jacobian saturation decision did not converge"))
end

function universal_machine_flux_current_jacobian_preview(;
    coil_conductances::AbstractVector{<:Real},
    coil_predictor_factors::AbstractVector{<:Real},
    coil_input_voltages::AbstractVector{<:Real},
    coil_history_currents::AbstractVector{<:Real},
    coil_resistances::AbstractVector{<:Real},
    axis_kinds,
    rotor_thevenin_matrix,
    stator_thevenin_matrix,
    d_axis_stator_current_indices::AbstractVector{<:Integer} = Int[],
    q_axis_stator_current_indices::AbstractVector{<:Integer} = Int[],
    d_axis_subtract_current_index::Int = 0,
    d_axis_unsaturated_inductance::Real,
    d_axis_saturated_inductance::Real,
    q_axis_unsaturated_inductance::Real,
    q_axis_saturated_inductance::Real,
    d_axis_saturation_flux::Real,
    q_axis_saturation_flux::Real,
    d_axis_remanent_flux::Real = 0.0,
    q_axis_remanent_flux::Real = 0.0,
    d_axis_saturation_mode::Int = 0,
    q_axis_saturation_mode::Int = 0,
    electrical_speed_rad_s::Real,
    half_step_interval_s::Real,
    direct_axis_coupling_scale::Real = 1.0,
    quadrature_axis_extra_impedance::Real = 0.0,
    cage_coupling_fraction::Real = 0.0,
    rotor_damping_resistance::Real = 0.0,
    direct_machine_voltage_coupling_fraction::Real = 0.0,
    direct_machine_resistance_coupling_fraction::Real = 0.0,
    quadrature_axis_feedback::Real = 0.0,
    direct_axis_feedback::Real = 0.0,
    direct_machine_transfer_current_enabled::Bool = true,
)
    coil_count = length(coil_conductances)
    4 <= coil_count <= 6 ||
        throw(ArgumentError("flux/current Jacobian requires three rotor coils plus one to three stator coils"))
    stator_count = coil_count - 3
    conductances = _machine_real_vector("coil_conductances", coil_conductances, coil_count)
    predictor_factors = _machine_real_vector("coil_predictor_factors", coil_predictor_factors, coil_count)
    input_voltages = _machine_real_vector("coil_input_voltages", coil_input_voltages, coil_count)
    history_currents = _machine_real_vector("coil_history_currents", coil_history_currents, coil_count)
    resistances = _machine_real_vector("coil_resistances", coil_resistances, coil_count)
    axes = _machine_symbol_vector("axis_kinds", axis_kinds, stator_count)
    all(kind -> kind in (:d_axis, :q_axis, :uncoupled), axes) ||
        throw(ArgumentError("axis_kinds entries must be :d_axis, :q_axis, or :uncoupled"))
    rotor_thevenin = _machine_real_matrix("rotor_thevenin_matrix", rotor_thevenin_matrix, 3, 3)
    stator_thevenin =
        _machine_real_matrix("stator_thevenin_matrix", stator_thevenin_matrix, stator_count, stator_count)

    half_step = Float64(half_step_interval_s)
    half_step > 0.0 || throw(ArgumentError("half_step_interval_s must be positive"))
    speed = Float64(electrical_speed_rad_s)
    direct_scale = Float64(direct_axis_coupling_scale)
    q_extra = Float64(quadrature_axis_extra_impedance)
    cage_fraction = Float64(cage_coupling_fraction)
    damping_resistance = Float64(rotor_damping_resistance)
    voltage_fraction = Float64(direct_machine_voltage_coupling_fraction)
    resistance_fraction = Float64(direct_machine_resistance_coupling_fraction)
    q_feedback = Float64(quadrature_axis_feedback)
    d_feedback = Float64(direct_axis_feedback)
    d_unsaturated = Float64(d_axis_unsaturated_inductance)
    d_saturated = Float64(d_axis_saturated_inductance)
    q_unsaturated = Float64(q_axis_unsaturated_inductance)
    q_saturated = Float64(q_axis_saturated_inductance)
    d_saturation_flux = Float64(d_axis_saturation_flux)
    q_saturation_flux = Float64(q_axis_saturation_flux)
    d_remanent = Float64(d_axis_remanent_flux)
    q_remanent = Float64(q_axis_remanent_flux)
    all(
        isfinite,
        (
            speed,
            direct_scale,
            q_extra,
            cage_fraction,
            damping_resistance,
            voltage_fraction,
            resistance_fraction,
            q_feedback,
            d_feedback,
            d_unsaturated,
            d_saturated,
            q_unsaturated,
            q_saturated,
            d_saturation_flux,
            q_saturation_flux,
            d_remanent,
            q_remanent,
        ),
    ) || throw(ArgumentError("flux/current Jacobian scalar inputs must be finite"))
    d_axis_saturation_mode >= 0 || throw(ArgumentError("d_axis_saturation_mode must be nonnegative"))
    q_axis_saturation_mode >= 0 || throw(ArgumentError("q_axis_saturation_mode must be nonnegative"))

    default_d_indices = [index + 3 for (index, kind) in pairs(axes) if kind == :d_axis]
    default_q_indices = [index + 3 for (index, kind) in pairs(axes) if kind == :q_axis]
    d_indices = isempty(d_axis_stator_current_indices) ?
                default_d_indices :
                Int[Int(index) for index in d_axis_stator_current_indices]
    q_indices = isempty(q_axis_stator_current_indices) ?
                default_q_indices :
                Int[Int(index) for index in q_axis_stator_current_indices]

    current_solution = function (d_flux::Float64, q_flux::Float64)
        return _universal_machine_current_solution_from_flux(
            coil_conductances = conductances,
            coil_predictor_factors = predictor_factors,
            coil_input_voltages = input_voltages,
            coil_history_currents = history_currents,
            coil_resistances = resistances,
            axis_kinds = axes,
            rotor_thevenin_matrix = rotor_thevenin,
            stator_thevenin_matrix = stator_thevenin,
            d_axis_flux = d_flux,
            q_axis_flux = q_flux,
            electrical_speed_rad_s = speed,
            half_step_interval_s = half_step,
            direct_axis_coupling_scale = direct_scale,
            quadrature_axis_extra_impedance = q_extra,
            cage_coupling_fraction = cage_fraction,
            rotor_damping_resistance = damping_resistance,
            direct_machine_voltage_coupling_fraction = voltage_fraction,
            direct_machine_resistance_coupling_fraction = resistance_fraction,
            quadrature_axis_feedback = q_feedback,
            direct_axis_feedback = d_feedback,
            direct_machine_transfer_current_enabled = direct_machine_transfer_current_enabled,
        )
    end

    base_currents = current_solution(0.0, 0.0)
    base_axes = _universal_machine_axis_current_totals(
        base_currents,
        d_indices,
        q_indices,
        d_axis_subtract_current_index,
    )
    d_probe_currents = current_solution(1.0, 0.0)
    d_probe_axes = _universal_machine_axis_current_totals(
        d_probe_currents,
        d_indices,
        q_indices,
        d_axis_subtract_current_index,
    )
    q_probe_currents = current_solution(0.0, 1.0)
    q_probe_axes = _universal_machine_axis_current_totals(
        q_probe_currents,
        d_indices,
        q_indices,
        d_axis_subtract_current_index,
    )

    region = _universal_machine_flux_region_solution(
        base_d_axis_current = base_axes.d_axis_current,
        base_q_axis_current = base_axes.q_axis_current,
        d_axis_d_flux_sensitivity = d_probe_axes.d_axis_current - base_axes.d_axis_current,
        q_axis_d_flux_sensitivity = d_probe_axes.q_axis_current - base_axes.q_axis_current,
        d_axis_q_flux_sensitivity = q_probe_axes.d_axis_current - base_axes.d_axis_current,
        q_axis_q_flux_sensitivity = q_probe_axes.q_axis_current - base_axes.q_axis_current,
        d_axis_unsaturated_inductance = d_unsaturated,
        d_axis_saturated_inductance = d_saturated,
        q_axis_unsaturated_inductance = q_unsaturated,
        q_axis_saturated_inductance = q_saturated,
        d_axis_saturation_flux = d_saturation_flux,
        q_axis_saturation_flux = q_saturation_flux,
        d_axis_remanent_flux = d_remanent,
        q_axis_remanent_flux = q_remanent,
        d_axis_saturation_mode = d_axis_saturation_mode,
        q_axis_saturation_mode = q_axis_saturation_mode,
    )
    final_currents = current_solution(region.d_axis_flux, region.q_axis_flux)
    final_axes = _universal_machine_axis_current_totals(
        final_currents,
        d_indices,
        q_indices,
        d_axis_subtract_current_index,
    )

    return (
        source = :universal_machine_flux_current_jacobian,
        current_values = final_currents,
        d_axis_flux = region.d_axis_flux,
        q_axis_flux = region.q_axis_flux,
        d_axis_current = final_axes.d_axis_current,
        q_axis_current = final_axes.q_axis_current,
        region_d_axis_current = region.d_axis_current,
        region_q_axis_current = region.q_axis_current,
        base_d_axis_current = base_axes.d_axis_current,
        base_q_axis_current = base_axes.q_axis_current,
        d_axis_d_flux_sensitivity = d_probe_axes.d_axis_current - base_axes.d_axis_current,
        q_axis_d_flux_sensitivity = d_probe_axes.q_axis_current - base_axes.q_axis_current,
        d_axis_q_flux_sensitivity = q_probe_axes.d_axis_current - base_axes.d_axis_current,
        q_axis_q_flux_sensitivity = q_probe_axes.q_axis_current - base_axes.q_axis_current,
        d_axis_saturated = region.d_region != 0.0,
        q_axis_saturated = region.q_region != 0.0,
        iteration_count = region.iteration_count,
        flux_current_jacobian_executed = true,
        flux_current_jacobian_mutated = false,
        complete_machine_equation_solution = false,
        full_machine_equation_solution = false,
        deferred_effects = (
            :multi_machine_orchestration,
            :tacs_transfer,
            :output_report_mutation,
            :full_machine_equation_solution,
        ),
    )
end

function universal_machine_flux_current_jacobian!(
    state::UniversalMachineFluxCurrentJacobianState;
    kwargs...,
)
    preview = universal_machine_flux_current_jacobian_preview(; kwargs...)
    empty!(state.current_values)
    append!(state.current_values, preview.current_values)
    state.d_axis_flux = preview.d_axis_flux
    state.q_axis_flux = preview.q_axis_flux
    state.d_axis_current = preview.d_axis_current
    state.q_axis_current = preview.q_axis_current
    state.base_d_axis_current = preview.base_d_axis_current
    state.base_q_axis_current = preview.base_q_axis_current
    state.d_axis_d_flux_sensitivity = preview.d_axis_d_flux_sensitivity
    state.q_axis_d_flux_sensitivity = preview.q_axis_d_flux_sensitivity
    state.d_axis_q_flux_sensitivity = preview.d_axis_q_flux_sensitivity
    state.q_axis_q_flux_sensitivity = preview.q_axis_q_flux_sensitivity
    state.d_axis_saturated = preview.d_axis_saturated
    state.q_axis_saturated = preview.q_axis_saturated
    state.iteration_count = preview.iteration_count
    state.flux_current_jacobian_mutated = true
    return merge(
        preview,
        (
            current_values = copy(state.current_values),
            d_axis_flux = state.d_axis_flux,
            q_axis_flux = state.q_axis_flux,
            d_axis_current = state.d_axis_current,
            q_axis_current = state.q_axis_current,
            base_d_axis_current = state.base_d_axis_current,
            base_q_axis_current = state.base_q_axis_current,
            d_axis_d_flux_sensitivity = state.d_axis_d_flux_sensitivity,
            q_axis_d_flux_sensitivity = state.q_axis_d_flux_sensitivity,
            d_axis_q_flux_sensitivity = state.d_axis_q_flux_sensitivity,
            q_axis_q_flux_sensitivity = state.q_axis_q_flux_sensitivity,
            d_axis_saturated = state.d_axis_saturated,
            q_axis_saturated = state.q_axis_saturated,
            iteration_count = state.iteration_count,
            flux_current_jacobian_mutated = state.flux_current_jacobian_mutated,
        ),
    )
end

function _induction_machine_clarke_transform()
    inverse_sqrt_three = inv(sqrt(3.0))
    inverse_sqrt_two = inv(sqrt(2.0))
    return [
        inverse_sqrt_three inverse_sqrt_three inverse_sqrt_three
        sqrt(2.0) * inverse_sqrt_three -inverse_sqrt_two * inverse_sqrt_three -inverse_sqrt_two * inverse_sqrt_three
        0.0 -inverse_sqrt_two inverse_sqrt_two
    ]
end

function _coupled_dq_power_base_transform(machine_type::Int)
    machine_type in (2, 5) && return [
        0.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0
    ]
    machine_type in (6, 8, 9, 10, 11, 12) && return [
        0.0 0.0 0.0
        0.0 0.0 0.0
        0.0 0.0 1.0
    ]
    machine_type == 7 && return [
        0.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 -1.0
    ]
    return _induction_machine_clarke_transform()
end

function _induction_machine_power_transform(
    electrical_angle_rad::Float64,
    machine_type::Int=3,
)
    base_transform = _coupled_dq_power_base_transform(machine_type)
    machine_type in (6, 8, 9, 10, 11, 12) && return base_transform
    rotation = [
        1.0 0.0 0.0
        0.0 cos(electrical_angle_rad) -sin(electrical_angle_rad)
        0.0 sin(electrical_angle_rad) cos(electrical_angle_rad)
    ]
    return rotation * base_transform
end

coupled_dq_power_transform(electrical_angle_rad::Real; machine_type::Integer=3) =
    _induction_machine_power_transform(
        Float64(electrical_angle_rad),
        Int(machine_type),
    )

function predict_machine_terminal_currents!(
    prediction::MachineTerminalPredictionState,
    state::InductionMachineState,
    parameters::InductionMachineParameters;
    time_s::Real,
)
    length(state.current_values) >= 3 ||
        throw(ArgumentError("terminal-current prediction requires zero-, d-, and q-axis currents"))
    time = Float64(time_s)
    isfinite(time) || throw(ArgumentError("terminal-current prediction time must be finite"))
    interval_s = 2.0 * parameters.half_step_interval_s
    interval_s > 0.0 ||
        throw(ArgumentError("terminal-current prediction interval must be positive"))
    reference_speed = parameters.synchronous_electrical_speed_rad_s
    axis_scale = parameters.machine_type in (6, 8, 9, 10, 11, 12) ? 0.0 : 1.0
    reference_angle = reference_speed * time
    rotor_electrical_angle =
        parameters.pole_pair_count * axis_scale * state.mechanical_angle_rad
    stationary_transform_angle =
        (reference_angle - rotor_electrical_angle) * axis_scale
    integration_interval = interval_s
    next_reference_angle =
        (reference_angle + interval_s * reference_speed) * axis_scale
    cosine = cos(stationary_transform_angle)
    sine = sin(stationary_transform_angle)
    d_axis_flux =
        cosine * state.d_axis_flux - sine * state.q_axis_flux
    q_axis_flux =
        sine * state.d_axis_flux + cosine * state.q_axis_flux
    d_axis_current =
        cosine * state.current_values[2] - sine * state.current_values[3]
    q_axis_current =
        sine * state.current_values[2] + cosine * state.current_values[3]
    d_axis_internal_flux =
        d_axis_flux + parameters.coil_reactances[2] * d_axis_current
    q_axis_internal_flux =
        q_axis_flux + parameters.coil_reactances[3] * q_axis_current

    previous_d_axis_flux = prediction.initialized ?
        prediction.previous_d_axis_flux_Wb : d_axis_flux
    previous_q_axis_flux = prediction.initialized ?
        prediction.previous_q_axis_flux_Wb : q_axis_flux
    previous_d_axis_internal_flux = prediction.initialized ?
        prediction.previous_d_axis_internal_flux_Wb : d_axis_internal_flux
    previous_q_axis_internal_flux = prediction.initialized ?
        prediction.previous_q_axis_internal_flux_Wb : q_axis_internal_flux
    predicted_d_axis_internal_flux =
        2.0 * d_axis_internal_flux - previous_d_axis_internal_flux
    predicted_q_axis_internal_flux =
        2.0 * q_axis_internal_flux - previous_q_axis_internal_flux
    stationary_drive = Float64[
        0.0,
        (
            (previous_d_axis_flux - d_axis_flux) / integration_interval -
            reference_speed * predicted_q_axis_internal_flux
        ) * parameters.coil_conductances[2],
        (
            (previous_q_axis_flux - q_axis_flux) / integration_interval +
            reference_speed * predicted_d_axis_internal_flux
        ) * parameters.coil_conductances[3],
    ]
    phase_current_injections =
        transpose(
            _induction_machine_power_transform(
                next_reference_angle,
                parameters.machine_type,
            ),
        ) * stationary_drive
    all(isfinite, phase_current_injections) ||
        throw(ArgumentError("terminal-current prediction produced a nonfinite injection"))

    prediction.previous_d_axis_flux_Wb = d_axis_flux
    prediction.previous_q_axis_flux_Wb = q_axis_flux
    prediction.previous_d_axis_internal_flux_Wb = d_axis_internal_flux
    prediction.previous_q_axis_internal_flux_Wb = q_axis_internal_flux
    prediction.phase_current_injections_A .= phase_current_injections
    prediction.update_count += 1
    prediction.initialized = true
    return copy(prediction.phase_current_injections_A)
end

function _induction_machine_stator_transform()
    shifted_rows = _induction_machine_clarke_transform()[[2, 3, 1], :]
    return shifted_rows[:, [2, 3, 1]]
end

function induction_machine_axis_fluxes(
    parameters::InductionMachineParameters,
    currents::AbstractVector{<:Real},
)
    current_values = _machine_real_vector(
        "induction-machine current_values",
        currents,
        length(parameters.coil_conductances),
    )
    d_start = 4
    q_start = d_start + parameters.d_axis_coil_count
    flux = universal_machine_flux_saturation_preview(
        current_values = current_values,
        d_axis_stator_current_indices = collect(
            d_start:(d_start + parameters.d_axis_coil_count - 1),
        ),
        q_axis_stator_current_indices = collect(
            q_start:(q_start + parameters.q_axis_coil_count - 1),
        ),
        d_axis_unsaturated_inductance = parameters.d_axis_unsaturated_inductance,
        d_axis_saturated_inductance = parameters.d_axis_saturated_inductance,
        q_axis_unsaturated_inductance = parameters.q_axis_unsaturated_inductance,
        q_axis_saturated_inductance = parameters.q_axis_saturated_inductance,
        d_axis_saturation_flux = parameters.d_axis_saturation_flux,
        q_axis_saturation_flux = parameters.q_axis_saturation_flux,
        d_axis_remanent_flux = parameters.d_axis_remanent_flux,
        q_axis_remanent_flux = parameters.q_axis_remanent_flux,
        d_axis_saturation_mode = parameters.d_axis_saturation_mode,
        q_axis_saturation_mode = parameters.q_axis_saturation_mode,
        d_axis_subtract_current_index = parameters.machine_type == 12 ? 5 : 0,
    )
    return flux.d_axis_flux, flux.q_axis_flux
end

function _induction_machine_generated_torque(
    parameters::InductionMachineParameters,
    currents::Vector{Float64},
    d_axis_flux::Float64,
    q_axis_flux::Float64,
)
    reluctance =
        (parameters.coil_reactances[2] - parameters.coil_reactances[3]) *
        currents[2] * currents[3]
    return (
        reluctance - currents[2] * q_axis_flux + currents[3] * d_axis_flux
    ) * parameters.pole_pair_count
end

function _induction_machine_voltage_context(
    parameters::InductionMachineParameters,
    mechanical_angle_rad::Float64,
    power_terminal_voltages::AbstractVector{<:Real},
    stator_terminal_voltages::AbstractVector{<:Real},
    rotor_thevenin_matrix,
    stator_thevenin_matrix,
)
    power_voltages = _machine_real_vector(
        "power_terminal_voltages",
        power_terminal_voltages,
        3,
    )
    rotor_thevenin = _machine_real_matrix(
        "rotor_thevenin_matrix",
        rotor_thevenin_matrix,
        3,
        3,
    )
    power_transform = _induction_machine_power_transform(
        parameters.pole_pair_count * mechanical_angle_rad,
        parameters.machine_type,
    )
    stator_count = length(parameters.coil_conductances) - 3
    stator_transform = if parameters.machine_type == 4
        _induction_machine_stator_transform()
    else
        Matrix{Float64}(I, stator_count, stator_count)
    end
    transformed_power_voltages = power_transform * power_voltages
    direct_machine_scalar_thevenin = parameters.machine_type in (9, 10, 11, 12) ?
        (power_transform * rotor_thevenin * transpose(power_transform))[3, 3] : 0.0
    direct_machine_series_resistance =
        parameters.machine_type == 11 && parameters.coil_conductances[5] != 0.0 ?
        inv(parameters.coil_conductances[5]) : 0.0
    stator_voltages = if parameters.machine_type == 4
        _machine_real_vector("stator_terminal_voltages", stator_terminal_voltages, 3)
    elseif parameters.machine_type in (1, 2, 6, 7, 8)
        _machine_real_vector(
            "excitation_terminal_voltages",
            stator_terminal_voltages,
            stator_count,
        )
    elseif parameters.machine_type in (9, 10, 11, 12)
        [-transformed_power_voltages[3], 0.0]
    else
        zeros(Float64, stator_count)
    end
    stator_thevenin = if parameters.machine_type == 4
        _machine_real_matrix("stator_thevenin_matrix", stator_thevenin_matrix, 3, 3)
    elseif parameters.machine_type in (1, 2, 6, 7, 8)
        _machine_real_matrix(
            "excitation_thevenin_matrix",
            stator_thevenin_matrix,
            stator_count,
            stator_count,
        )
    elseif parameters.machine_type in (9, 10, 11, 12)
        [direct_machine_scalar_thevenin + direct_machine_series_resistance 0.0; 0.0 0.0]
    else
        zeros(Float64, stator_count, stator_count)
    end
    transformed_rotor_thevenin =
        power_transform * rotor_thevenin * transpose(power_transform)
    parameters.machine_type in (9, 10, 11, 12) && (transformed_rotor_thevenin[3, 3] = 0.0)
    return (
        power_transform = power_transform,
        stator_transform = stator_transform,
        input_voltages = vcat(
            transformed_power_voltages,
            stator_transform * stator_voltages,
        ),
        rotor_thevenin_matrix = transformed_rotor_thevenin,
        stator_thevenin_matrix =
            stator_transform * stator_thevenin * transpose(stator_transform),
        direct_machine_scalar_thevenin = direct_machine_scalar_thevenin,
    )
end

function _wound_field_synchronous_torque_angle(
    parameters::InductionMachineParameters,
    state::InductionMachineState,
)
    elapsed_time_s = state.call_count * (2.0 * parameters.half_step_interval_s)
    synchronous_reference_angle =
        pi / 2.0 + parameters.synchronous_electrical_speed_rad_s * elapsed_time_s
    return rem2pi(
        parameters.pole_pair_count * state.mechanical_angle_rad -
        synchronous_reference_angle,
        RoundNearest,
    )
end

function induction_machine_step!(
    state::InductionMachineState,
    parameters::InductionMachineParameters;
    power_terminal_voltages::AbstractVector{<:Real},
    rotor_thevenin_matrix,
    mechanical_speed_thevenin_rad_s::Real,
    generated_torque_impedance::Real,
    stator_terminal_voltages::AbstractVector{<:Real}=zeros(3),
    stator_thevenin_matrix=zeros(3, 3),
    coil_control_voltages::AbstractVector{<:Real}=Float64[],
    prescribed_power_terminal_currents::AbstractVector{<:Real}=Float64[],
    initial_step::Bool=false,
    max_speed_iterations::Int=51,
)
    max_speed_iterations > 0 ||
        throw(ArgumentError("max_speed_iterations must be positive"))
    speed_thevenin = Float64(mechanical_speed_thevenin_rad_s)
    torque_impedance = Float64(generated_torque_impedance)
    all(isfinite, (speed_thevenin, torque_impedance)) ||
        throw(ArgumentError("mechanical Thevenin inputs must be finite"))

    old_speed = state.mechanical_speed_rad_s
    old_angle = state.mechanical_angle_rad
    predicted_speed = 2.0 * old_speed - state.previous_mechanical_speed_rad_s
    currents = copy(state.current_values)
    d_axis_flux = state.d_axis_flux
    q_axis_flux = state.q_axis_flux
    generated_torque = state.generated_torque
    iteration_count = 0
    converged = initial_step
    speed_residual_rad_s = 0.0
    coil_count = length(parameters.coil_conductances)
    control_voltages = isempty(coil_control_voltages) ?
        zeros(Float64, coil_count) :
        _machine_real_vector("coil_control_voltages", coil_control_voltages, coil_count)
    prescribed_phase_currents = isempty(prescribed_power_terminal_currents) ?
        nothing :
        _machine_real_vector(
            "prescribed_power_terminal_currents",
            prescribed_power_terminal_currents,
            3,
        )
    length(state.current_values) == coil_count && length(state.history_currents) == coil_count ||
        throw(ArgumentError("coupled d/q machine parameter and state coil counts differ"))
    stator_axes = if parameters.machine_type in (11, 12)
        # OVER16 cancels the ordinary coil-4 d-axis flux voltage with the
        # parallel-field resistance-coupling term.  Coil 4 still contributes
        # to total d-axis flux, but its individual voltage equation is uncoupled.
        [:uncoupled, :d_axis]
    else
        vcat(
            fill(:d_axis, parameters.d_axis_coil_count),
            fill(:q_axis, parameters.q_axis_coil_count),
            parameters.machine_type == 4 ? [:uncoupled] : Symbol[],
        )
    end
    voltage_context = _induction_machine_voltage_context(
        parameters,
        old_angle,
        power_terminal_voltages,
        stator_terminal_voltages,
        rotor_thevenin_matrix,
        stator_thevenin_matrix,
    )
    controlled_inputs = copy(voltage_context.input_voltages)
    controlled_inputs[4:end] .-= control_voltages[4:end]
    voltage_context = merge(voltage_context, (; input_voltages = controlled_inputs))

    if initial_step
        d_axis_flux, q_axis_flux = induction_machine_axis_fluxes(parameters, currents)
        generated_torque = _induction_machine_generated_torque(
            parameters,
            currents,
            d_axis_flux,
            q_axis_flux,
        )
    else
        for iteration in 1:max_speed_iterations
            iteration_count = iteration
            angle = old_angle + parameters.half_step_interval_s *
                                (old_speed + predicted_speed)
            voltage_context = _induction_machine_voltage_context(
                parameters,
                angle,
                power_terminal_voltages,
                stator_terminal_voltages,
                rotor_thevenin_matrix,
                stator_thevenin_matrix,
            )
            controlled_inputs = copy(voltage_context.input_voltages)
            controlled_inputs[4:end] .-= control_voltages[4:end]
            voltage_context = merge(voltage_context, (; input_voltages = controlled_inputs))
            direct_machine = parameters.machine_type in (9, 10, 11, 12)
            resistance_coupled_direct_machine = parameters.machine_type in (11, 12)
            scalar_thevenin = voltage_context.direct_machine_scalar_thevenin
            series_resistance = parameters.machine_type in (9, 10, 11) &&
                parameters.coil_conductances[5] != 0.0 ?
                inv(parameters.coil_conductances[5]) : 0.0
            q_axis_admittance =
                parameters.coil_conductances[3] * parameters.coil_predictor_factors[3]
            compound_field_admittance =
                parameters.coil_conductances[4] * parameters.coil_predictor_factors[4]
            direct_feedback_impedance = scalar_thevenin + (
                resistance_coupled_direct_machine ?
                (parameters.series_path_leakage_inductance_h + series_resistance) /
                parameters.half_step_interval_s : 0.0
            )
            solve_conductances = parameters.coil_conductances
            solve_histories = state.history_currents
            if prescribed_phase_currents !== nothing
                prescribed_axis_currents =
                    voltage_context.power_transform * prescribed_phase_currents
                solve_conductances = copy(parameters.coil_conductances)
                solve_conductances[1:3] .= 0.0
                solve_histories = copy(state.history_currents)
                solve_histories[1:3] .= prescribed_axis_currents
            end
            flux_current = universal_machine_flux_current_jacobian_preview(
                coil_conductances = solve_conductances,
                coil_predictor_factors = parameters.coil_predictor_factors,
                coil_input_voltages = voltage_context.input_voltages,
                coil_history_currents = solve_histories,
                # The extracted equations use shifted leakage reactances in these terms.
                coil_resistances = parameters.coil_reactances,
                axis_kinds = stator_axes,
                rotor_thevenin_matrix = voltage_context.rotor_thevenin_matrix,
                stator_thevenin_matrix = voltage_context.stator_thevenin_matrix,
                d_axis_unsaturated_inductance =
                    parameters.d_axis_unsaturated_inductance,
                d_axis_saturated_inductance =
                    parameters.d_axis_saturated_inductance,
                q_axis_unsaturated_inductance =
                    parameters.q_axis_unsaturated_inductance,
                q_axis_saturated_inductance =
                    parameters.q_axis_saturated_inductance,
                d_axis_saturation_flux = parameters.d_axis_saturation_flux,
                q_axis_saturation_flux = parameters.q_axis_saturation_flux,
                d_axis_remanent_flux = parameters.d_axis_remanent_flux,
                q_axis_remanent_flux = parameters.q_axis_remanent_flux,
                d_axis_saturation_mode = parameters.d_axis_saturation_mode,
                q_axis_saturation_mode = parameters.q_axis_saturation_mode,
                d_axis_stator_current_indices = parameters.machine_type in (11, 12) ?
                    [4, 5] : Int[],
                d_axis_subtract_current_index = parameters.machine_type == 12 ? 5 : 0,
                electrical_speed_rad_s = parameters.pole_pair_count * predicted_speed,
                half_step_interval_s = parameters.half_step_interval_s,
                quadrature_axis_extra_impedance =
                    direct_machine ? scalar_thevenin : 0.0,
                cage_coupling_fraction = direct_machine ? 1.0 : 0.0,
                rotor_damping_resistance = series_resistance,
                direct_machine_voltage_coupling_fraction =
                    direct_machine && !resistance_coupled_direct_machine ? 1.0 : 0.0,
                direct_machine_resistance_coupling_fraction =
                    resistance_coupled_direct_machine ? 1.0 : 0.0,
                quadrature_axis_feedback =
                    direct_machine ? q_axis_admittance * direct_feedback_impedance : 0.0,
                direct_axis_feedback = direct_machine ?
                    compound_field_admittance * direct_feedback_impedance : 0.0,
                direct_machine_transfer_current_enabled = parameters.machine_type != 12,
            )
            currents = flux_current.current_values
            d_axis_flux = flux_current.d_axis_flux
            q_axis_flux = flux_current.q_axis_flux
            generated_torque = _induction_machine_generated_torque(
                parameters,
                currents,
                d_axis_flux,
                q_axis_flux,
            )
            solved_speed = speed_thevenin - torque_impedance * generated_torque
            speed_residual_rad_s = abs(solved_speed - predicted_speed)
            if speed_residual_rad_s <= parameters.speed_tolerance_rad_s
                state.mechanical_angle_rad = angle
                converged = true
                break
            end
            predicted_speed = solved_speed
        end
        converged || throw(
            MachineConvergenceError(
                parameters.machine_type,
                state.call_count + 1,
                iteration_count,
                speed_residual_rad_s,
                parameters.speed_tolerance_rad_s,
            ),
        )
        state.previous_mechanical_speed_rad_s = old_speed
        state.mechanical_speed_rad_s = predicted_speed
    end

    history_rotor_thevenin =
        initial_step ? zeros(3, 3) : voltage_context.rotor_thevenin_matrix
    history_stator_thevenin =
        initial_step ? zeros(length(stator_axes), length(stator_axes)) :
        voltage_context.stator_thevenin_matrix
    direct_machine = parameters.machine_type in (9, 10, 11, 12)
    resistance_coupled_direct_machine = parameters.machine_type in (11, 12)
    scalar_thevenin = voltage_context.direct_machine_scalar_thevenin
    series_resistance = parameters.machine_type in (9, 10, 11) &&
        parameters.coil_conductances[5] != 0.0 ?
        inv(parameters.coil_conductances[5]) : 0.0
    history = universal_machine_history_current_update_preview(
        current_values = currents,
        input_voltages = voltage_context.input_voltages,
        coil_conductances = parameters.coil_conductances,
        coil_predictor_factors = parameters.coil_predictor_factors,
        coil_resistances = parameters.coil_reactances,
        axis_kinds = stator_axes,
        rotor_thevenin_matrix = history_rotor_thevenin,
        stator_thevenin_matrix = history_stator_thevenin,
        d_axis_flux = d_axis_flux,
        q_axis_flux = q_axis_flux,
        electrical_speed_rad_s =
            parameters.pole_pair_count * state.mechanical_speed_rad_s,
        half_step_interval_s = parameters.half_step_interval_s,
        three_phase_induction_terminal_coil = parameters.machine_type == 4,
        quadrature_axis_extra_impedance = direct_machine ? scalar_thevenin : 0.0,
        direct_machine_voltage_coupling_fraction =
            direct_machine && !resistance_coupled_direct_machine ? 1.0 : 0.0,
        direct_machine_resistance_coupling_fraction =
            resistance_coupled_direct_machine ? 1.0 : 0.0,
        direct_machine_rotor_resistance = series_resistance,
        direct_machine_stator_history_coupling_enabled = !direct_machine,
    )

    power_currents = transpose(voltage_context.power_transform) * currents[1:3]
    parameters.machine_type in (9, 10, 11, 12) && (power_currents[3] -= currents[4])
    stator_currents = transpose(voltage_context.stator_transform) * currents[4:coil_count]
    reported_angle = parameters.machine_type in (1, 2) ?
        _wound_field_synchronous_torque_angle(parameters, state) :
        state.mechanical_angle_rad
    outputs = vcat(
        generated_torque,
        state.mechanical_speed_rad_s,
        reported_angle,
        power_currents,
        -stator_currents,
    )
    state.current_values .= currents
    state.history_currents .= history.history_currents
    state.d_axis_flux = d_axis_flux
    state.q_axis_flux = q_axis_flux
    state.generated_torque = generated_torque
    state.output_values .= outputs
    state.call_count += 1

    return (
        source = :coupled_dq_machine_step,
        machine_type = parameters.machine_type,
        current_values = copy(state.current_values),
        history_currents = copy(state.history_currents),
        d_axis_flux = state.d_axis_flux,
        q_axis_flux = state.q_axis_flux,
        generated_torque = state.generated_torque,
        mechanical_speed_rad_s = state.mechanical_speed_rad_s,
        mechanical_angle_rad = reported_angle,
        rotor_mechanical_angle_rad = state.mechanical_angle_rad,
        power_terminal_currents = power_currents,
        stator_terminal_currents = stator_currents,
        current_substitution_values = parameters.machine_type == 4 ?
            vcat(-power_currents, stator_currents, generated_torque) :
            parameters.machine_type in (1, 2) ?
            vcat(-power_currents, stator_currents[1], generated_torque) :
            parameters.machine_type in (6, 8) ?
            [
                -power_currents[3],
                stator_currents[1],
                generated_torque,
            ] :
            parameters.machine_type in (9, 10, 11, 12) ?
            [
                -power_currents[3],
                generated_torque,
            ] :
            parameters.machine_type == 7 ?
            [
                -power_currents[3],
                stator_currents[1],
                stator_currents[2],
                generated_torque,
            ] :
            vcat(-power_currents, generated_torque),
        output_values = copy(state.output_values),
        coil_control_voltages = copy(control_voltages),
        iteration_count = iteration_count,
        converged = converged,
        call_count = state.call_count,
        state_mutated = true,
        terminal_rhs_mutated = true,
        output_mutated = true,
        full_coupled_dq_machine_equation_step = true,
        full_induction_machine_equation_step = parameters.machine_type in (3, 4, 5, 6, 7),
        full_wound_field_synchronous_equation_step = parameters.machine_type in (1, 2),
        full_two_phase_synchronous_equation_step = parameters.machine_type == 2,
        full_single_phase_machine_equation_step = parameters.machine_type in (6, 7),
        full_two_phase_rotor_machine_equation_step = parameters.machine_type == 7,
        full_separately_excited_dc_machine_equation_step = parameters.machine_type == 8,
        full_series_compound_dc_machine_equation_step = parameters.machine_type == 9,
        full_series_field_dc_machine_equation_step = parameters.machine_type == 10,
        full_parallel_compound_dc_machine_equation_step = parameters.machine_type == 11,
        full_self_excited_shunt_dc_machine_equation_step = parameters.machine_type == 12,
        full_type4_machine_equation_step = parameters.machine_type == 4,
    )
end

coupled_dq_machine_step!(
    state::InductionMachineState,
    parameters::InductionMachineParameters;
    kwargs...,
) = induction_machine_step!(state, parameters; kwargs...)

function universal_machine_type4_step!(
    state::InductionMachineState,
    parameters::InductionMachineParameters;
    kwargs...,
)
    parameters.machine_type == 4 ||
        throw(ArgumentError("universal_machine_type4_step! requires machine type 4"))
    return induction_machine_step!(state, parameters; kwargs...)
end

function universal_machine_shared_torque_prediction_preview(;
    current_torques::AbstractVector{<:Real},
    previous_torque_history::AbstractVector{<:Real},
    shared_network_counts::AbstractVector{<:Integer},
    substitution_slot_indices::AbstractVector{<:Integer},
    current_substitution_values::AbstractVector{<:Real} = Float64[],
    initial_step::Bool = false,
)
    machine_count = length(current_torques)
    machine_count == length(previous_torque_history) == length(shared_network_counts) ==
        length(substitution_slot_indices) ||
        throw(ArgumentError("shared torque vectors must have matching machine counts"))
    torques = _machine_real_vector("current_torques", current_torques, machine_count)
    history = _machine_real_vector("previous_torque_history", previous_torque_history, machine_count)
    share_counts = _machine_int_vector("shared_network_counts", shared_network_counts, machine_count)
    substitution_indices = _machine_int_vector(
        "substitution_slot_indices",
        substitution_slot_indices,
        machine_count,
    )
    all(count -> count >= 0, share_counts) ||
        throw(ArgumentError("shared_network_counts entries must be nonnegative"))
    maximum_slot = maximum(max(index, 0) for index in substitution_indices; init = 0)
    substitutions =
        isempty(current_substitution_values) ?
        zeros(Float64, maximum_slot) :
        _machine_real_vector("current_substitution_values", current_substitution_values,
            length(current_substitution_values))
    length(substitutions) >= maximum_slot ||
        throw(ArgumentError("current_substitution_values is shorter than the requested substitution slot"))

    predicted = copy(history)
    next_history = copy(history)
    active_indices = Int[]
    written_indices = Int[]
    for machine_index in 1:machine_count
        share_counts[machine_index] == 0 && continue
        push!(active_indices, machine_index)
        predicted[machine_index] =
            initial_step ? torques[machine_index] : 2.0 * torques[machine_index] - history[machine_index]
        next_history[machine_index] = torques[machine_index]
        slot = substitution_indices[machine_index]
        if slot > 0
            substitutions[slot] = torques[machine_index]
            push!(written_indices, slot)
        end
    end

    return (
        source = :universal_machine_shared_torque_prediction,
        predicted_torques = predicted,
        torque_history = next_history,
        current_substitution_values = substitutions,
        active_machine_indices = active_indices,
        written_substitution_indices = written_indices,
        shared_machine_count = length(active_indices),
        substitution_write_count = length(written_indices),
        shared_torque_prediction_executed = true,
        shared_torque_prediction_mutated = false,
        complete_machine_equation_solution = false,
        full_machine_equation_solution = false,
        deferred_effects = (
            :tacs_transfer,
            :output_report_mutation,
            :full_machine_equation_solution,
        ),
    )
end

function universal_machine_shared_torque_prediction!(
    state::UniversalMachineSharedTorquePredictionState;
    kwargs...,
)
    preview = universal_machine_shared_torque_prediction_preview(; kwargs...)
    empty!(state.predicted_torques)
    empty!(state.torque_history)
    empty!(state.current_substitution_values)
    empty!(state.active_machine_indices)
    empty!(state.written_substitution_indices)
    append!(state.predicted_torques, preview.predicted_torques)
    append!(state.torque_history, preview.torque_history)
    append!(state.current_substitution_values, preview.current_substitution_values)
    append!(state.active_machine_indices, preview.active_machine_indices)
    append!(state.written_substitution_indices, preview.written_substitution_indices)
    state.shared_torque_prediction_mutated = true
    return merge(
        preview,
        (
            predicted_torques = copy(state.predicted_torques),
            torque_history = copy(state.torque_history),
            current_substitution_values = copy(state.current_substitution_values),
            active_machine_indices = copy(state.active_machine_indices),
            written_substitution_indices = copy(state.written_substitution_indices),
            shared_torque_prediction_mutated = state.shared_torque_prediction_mutated,
        ),
    )
end

function universal_machine_postsolve_update_preview(;
    coil_parameters::AbstractVector{<:Real},
    source_crests::AbstractVector{<:Real}=Float64[],
    current_values::AbstractVector{<:Real},
    prediction_values::AbstractVector{<:Real},
    history_values::AbstractVector{<:Real},
    generated_torque::Real=0.0,
    shared_network_active::Bool=false,
    type59_exciter_active::Bool=false,
    stored_input_mode::Int=0,
    machine_type::Int=0,
    input_mode::Int=0,
    connection_index::Int=0,
    exciter_voltage_node::Int=0,
    exciter_current_node::Int=0,
    node_values::AbstractVector{<:Real}=Float64[],
    output_scale::Real=1.0,
    power_terminal_prediction_enabled::Bool=false,
    machine_index::Int=1,
    machine_count::Int=1,
    d_axis_flux::Real=0.0,
    q_axis_flux::Real=0.0,
    d_axis_reactance::Real=0.0,
    q_axis_reactance::Real=0.0,
    coil_conductances::AbstractVector{<:Real}=Float64[],
    tau_rad::Real=0.0,
    theta_electric_rad::Real=0.0,
    reference_speed_rad_s::Real=0.0,
    time_step_s::Real=0.0,
    integration_step_s::Real=1.0,
    con1::Real=0.0,
    con4::Real=0.0,
    initial_step::Bool=false,
    previous_d_axis_flux_memory::Real=0.0,
    previous_q_axis_flux_memory::Real=0.0,
    history_start::Int=1,
    prediction_loop_marker::Int=0,
    emit_prediction_report_text::Bool=false,
)
    machine_index > 0 || throw(ArgumentError("machine_index must be positive"))
    machine_count >= machine_index || throw(ArgumentError("machine_count must include machine_index"))
    history_start > 0 || throw(ArgumentError("history_start must be positive"))
    params = _machine_float_vector(coil_parameters)
    length(params) >= 6 || throw(ArgumentError("coil_parameters must contain at least six entries"))
    crests = _machine_float_vector(source_crests)
    currents = _machine_float_vector(current_values)
    length(currents) >= 4 || throw(ArgumentError("current_values must contain zero/d/q and exciter current entries"))
    predictions = _machine_float_vector(prediction_values)
    histories = _machine_float_vector(history_values)
    length(histories) >= history_start + 2 ||
        throw(ArgumentError("history_values must contain the three power-coil history entries"))
    gpar =
        isempty(coil_conductances) ?
        ones(Float64, 3) :
        _machine_real_vector("coil_conductances", coil_conductances, length(coil_conductances))
    length(gpar) >= 3 || throw(ArgumentError("coil_conductances must contain at least three entries"))
    nodes = _machine_float_vector(node_values)
    torque = Float64(generated_torque)
    scale = Float64(output_scale)
    d_flux = Float64(d_axis_flux)
    q_flux = Float64(q_axis_flux)
    d_reactance = Float64(d_axis_reactance)
    q_reactance = Float64(q_axis_reactance)
    tau = Float64(tau_rad)
    theta = Float64(theta_electric_rad)
    reference_speed = Float64(reference_speed_rad_s)
    time_step = Float64(time_step_s)
    integration_step = Float64(integration_step_s)
    c1 = Float64(con1)
    c4 = Float64(con4)
    previous_d_memory = Float64(previous_d_axis_flux_memory)
    previous_q_memory = Float64(previous_q_axis_flux_memory)
    scalars = (
        torque,
        scale,
        d_flux,
        q_flux,
        d_reactance,
        q_reactance,
        tau,
        theta,
        reference_speed,
        time_step,
        integration_step,
        c1,
        c4,
        previous_d_memory,
        previous_q_memory,
    )
    all(isfinite, scalars) || throw(ArgumentError("postsolve scalar inputs must be finite"))
    integration_step != 0.0 || throw(ArgumentError("integration_step_s must be nonzero"))

    mutated_machine_type = machine_type
    mutated_input_mode = input_mode
    final_prediction_loop_marker = prediction_loop_marker
    shared_torque_stored = false
    exciter_torque_injected = false
    prediction_executed = false
    drive_start_index = 0
    previous_current_index = 0
    drive_values = Float64[]

    if shared_network_active
        params[3] = torque
        shared_torque_stored = true
    end

    if type59_exciter_active
        mutated_machine_type = 13
        mutated_input_mode = stored_input_mode
        if exciter_voltage_node != 0
            connection_index > 0 || throw(ArgumentError("connection_index must be positive for exciter injection"))
            _machine_check_index(connection_index, crests, "source_crests")
            _machine_check_index(exciter_voltage_node, nodes, "node_values")
            _machine_check_index(exciter_current_node, nodes, "node_values")
            nodes[exciter_voltage_node] != 0.0 ||
                throw(ArgumentError("exciter voltage node value must be nonzero"))
            exciter_current = currents[4] * scale
            crests[connection_index] = -nodes[exciter_current_node] * exciter_current / nodes[exciter_voltage_node]
            params[4] = -connection_index
            params[5] = exciter_voltage_node
            params[6] = exciter_current_node
            exciter_torque_injected = true
        end
    end

    if power_terminal_prediction_enabled
        previous_current_index = 3 * (machine_index - 1) + 1
        drive_start_index = 3 * (machine_count + machine_index - 1) + 1
        length(predictions) >= drive_start_index + 2 ||
            throw(ArgumentError("prediction_values must contain previous-current and next-step drive slots"))
        angle_scale = 1.0 - c1 - c4
        transform_angle = (tau - theta) * angle_scale
        if mutated_input_mode == 1
            transform_angle *= reference_speed
        end
        cos_angle = cos(transform_angle)
        sin_angle = sin(transform_angle)
        original_d_flux = d_flux
        original_q_flux = q_flux
        original_d_current = currents[2]
        original_q_current = currents[3]
        d_flux = cos_angle * original_d_flux - sin_angle * original_q_flux
        q_flux = sin_angle * original_d_flux + cos_angle * original_q_flux
        currents[2] = cos_angle * original_d_current - sin_angle * original_q_current
        currents[3] = sin_angle * original_d_current + cos_angle * original_q_current
        d_axis_internal_flux = d_flux + d_reactance * currents[2]
        q_axis_internal_flux = q_flux + q_reactance * currents[3]
        d_memory = initial_step ? d_flux : previous_d_memory
        q_memory = initial_step ? q_flux : previous_q_memory
        if initial_step
            predictions[previous_current_index] = 0.0
            histories[history_start] = 0.0
            histories[history_start + 1] = d_axis_internal_flux
            histories[history_start + 2] = q_axis_internal_flux
        end
        d_axis_predicted_flux = 2.0 * d_axis_internal_flux - histories[history_start + 1]
        q_axis_predicted_flux = 2.0 * q_axis_internal_flux - histories[history_start + 2]
        d_derivative = (d_memory - d_flux) / integration_step
        q_derivative = (q_memory - q_flux) / integration_step
        drive_values = [
            0.0,
            (d_derivative - reference_speed * q_axis_predicted_flux) * gpar[2],
            (q_derivative + reference_speed * d_axis_predicted_flux) * gpar[3],
        ]
        predictions[previous_current_index] = currents[1]
        histories[history_start + 1] = d_axis_internal_flux
        histories[history_start + 2] = q_axis_internal_flux
        next_theta = (tau + time_step * reference_speed) * angle_scale
        if mutated_input_mode == 1
            next_theta *= reference_speed
        end
        theta = next_theta
        for offset in 0:2
            predictions[drive_start_index + offset] = drive_values[offset + 1]
        end
        if machine_index == machine_count
            final_prediction_loop_marker = 3
        end
        prediction_executed = true
    end

    prediction_report_text_formatting = emit_prediction_report_text && prediction_executed
    prediction_report_text_lines =
        prediction_report_text_formatting ?
        _universal_machine_prediction_report_text_lines(
            machine_index,
            drive_start_index,
            drive_start_index + 2,
            predictions[drive_start_index:(drive_start_index + 2)],
        ) :
        String[]

    return (
        source = :universal_machine_postsolve_update,
        fortran_labels =
            prediction_report_text_formatting ?
            (16100, 18000, 18020, 18030, 18100, 18120, 18130, 18132) :
            (16100, 18000, 18020, 18030, 18100, 18120),
        coil_parameters = params,
        source_crests = crests,
        current_values = currents,
        prediction_values = predictions,
        history_values = histories,
        d_axis_flux = d_flux,
        q_axis_flux = q_flux,
        theta_electric_rad = theta,
        machine_type = mutated_machine_type,
        input_mode = mutated_input_mode,
        final_prediction_loop_marker = final_prediction_loop_marker,
        previous_current_index = previous_current_index,
        drive_start_index = drive_start_index,
        drive_values = drive_values,
        prediction_report_text_lines = prediction_report_text_lines,
        prediction_report_text_line_count = length(prediction_report_text_lines),
        shared_torque_stored = shared_torque_stored,
        exciter_torque_injected = exciter_torque_injected,
        power_terminal_prediction_executed = prediction_executed,
        postsolve_update_mutated = false,
        prediction_report_text_mutated = false,
        prediction_report_text_formatting = prediction_report_text_formatting,
        complete_output_report_mutation = false,
        complete_machine_equation_solution = false,
        full_machine_equation_solution = false,
        deferred_effects = (
            :complete_output_report_mutation,
            :complete_machine_equation_solution,
            :full_machine_equation_solution,
        ),
    )
end

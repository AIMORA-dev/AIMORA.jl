export MachineControlState,
       MachineShaftState,
       MachinePortInputs,
       machine_control_state,
       machine_shaft_state,
       machine_port_inputs,
       update_machine_control_state!,
       machine_shaft_derivative,
       advance_machine_shaft!

mutable struct MachineControlState
    sensed_voltage_v::Float64
    sampled_speed_rad_s::Float64
    excitation_state_v::Float64
    governor_state_nm::Float64
    stabilizer_washout_state::Float64
    stabilizer_lead_lag_state::Float64
    field_voltage_v::Float64
    mechanical_torque_nm::Float64
    field_limited::Bool
    torque_limited::Bool
    next_task_time_s::Float64
    sample_count::Int
end

mutable struct MachineShaftState
    angle_rad::Vector{Float64}
    speed_rad_s::Vector{Float64}
    coupling_torque_nm::Vector{Float64}
    kinetic_energy_j::Float64
    elastic_energy_j::Float64
    damping_loss_w::Float64
    maximum_angular_momentum_residual_nms::Float64
end

mutable struct _MachineShaftWorkspace
    old_angle_rad::Vector{Float64}
    old_speed_rad_s::Vector{Float64}
    predicted_angle_rad::Vector{Float64}
    predicted_speed_rad_s::Vector{Float64}
    first_acceleration_rad_s2::Vector{Float64}
    candidate_acceleration_rad_s2::Vector{Float64}
    first_coupling_torque_nm::Vector{Float64}
    candidate_coupling_torque_nm::Vector{Float64}
    mass_torque_nm::Vector{Float64}
end

mutable struct MachinePortInputs
    field_voltage_v::Float64
    rotor_voltage_dq_v::NTuple{2,Float64}
    mechanical_torque_nm::Float64
    voltage_reference_v::Float64
    speed_reference_rad_s::Float64
    control_enabled::Bool
end

function _copy_machine_control_state!(
    destination::MachineControlState,
    source::MachineControlState,
)
    destination.sensed_voltage_v = source.sensed_voltage_v
    destination.sampled_speed_rad_s = source.sampled_speed_rad_s
    destination.excitation_state_v = source.excitation_state_v
    destination.governor_state_nm = source.governor_state_nm
    destination.stabilizer_washout_state = source.stabilizer_washout_state
    destination.stabilizer_lead_lag_state = source.stabilizer_lead_lag_state
    destination.field_voltage_v = source.field_voltage_v
    destination.mechanical_torque_nm = source.mechanical_torque_nm
    destination.field_limited = source.field_limited
    destination.torque_limited = source.torque_limited
    destination.next_task_time_s = source.next_task_time_s
    destination.sample_count = source.sample_count
    return destination
end

function _copy_machine_shaft_state!(
    destination::MachineShaftState,
    source::MachineShaftState,
)
    length(destination.angle_rad) == length(source.angle_rad) ||
        throw(DimensionMismatch("machine shaft angle state lengths differ"))
    length(destination.speed_rad_s) == length(source.speed_rad_s) ||
        throw(DimensionMismatch("machine shaft speed state lengths differ"))
    length(destination.coupling_torque_nm) == length(source.coupling_torque_nm) ||
        throw(DimensionMismatch("machine shaft coupling state lengths differ"))
    copyto!(destination.angle_rad, source.angle_rad)
    copyto!(destination.speed_rad_s, source.speed_rad_s)
    copyto!(destination.coupling_torque_nm, source.coupling_torque_nm)
    destination.kinetic_energy_j = source.kinetic_energy_j
    destination.elastic_energy_j = source.elastic_energy_j
    destination.damping_loss_w = source.damping_loss_w
    destination.maximum_angular_momentum_residual_nms =
        source.maximum_angular_momentum_residual_nms
    return destination
end

function machine_port_inputs(specification::ModernMachineSpecification)
    controls = specification.controls
    return MachinePortInputs(
        specification.initial_field_voltage_v,
        specification.initial_rotor_voltage_dq_v,
        specification.initial_mechanical_torque_nm,
        controls.voltage_reference_v,
        controls.speed_reference_rad_s,
        controls.enabled,
    )
end

function machine_control_state(specification::ModernMachineSpecification)
    controls = specification.controls
    field_voltage = clamp(
        specification.initial_field_voltage_v,
        controls.field_voltage_min_v,
        controls.field_voltage_max_v,
    )
    torque = clamp(
        specification.initial_mechanical_torque_nm,
        controls.torque_min_nm,
        controls.torque_max_nm,
    )
    return MachineControlState(
        0.0,
        first(specification.shaft_masses).initial_speed_rad_s,
        field_voltage,
        torque,
        0.0,
        0.0,
        field_voltage,
        torque,
        field_voltage != specification.initial_field_voltage_v,
        torque != specification.initial_mechanical_torque_nm,
        controls.task_phase_s,
        0,
    )
end

function _machine_shaft_mass_index(
    specification::ModernMachineSpecification,
    mass_id::Symbol,
)
    for index in eachindex(specification.shaft_masses)
        specification.shaft_masses[index].id === mass_id && return index
    end
    error("machine shaft mass disappeared from its validated topology")
end

function _machine_shaft_energy(
    specification::ModernMachineSpecification,
    angle_rad::Vector{Float64},
    speed_rad_s::Vector{Float64},
)
    kinetic = sum(
        0.5 * mass.inertia_kg_m2 * speed_rad_s[index]^2 for
        (index, mass) in enumerate(specification.shaft_masses)
    )
    elastic = 0.0
    for coupling in specification.shaft_couplings
        left = _machine_shaft_mass_index(specification, coupling.left_mass)
        right = _machine_shaft_mass_index(specification, coupling.right_mass)
        twist = angle_rad[left] - angle_rad[right] - coupling.initial_twist_rad
        elastic += 0.5 * coupling.stiffness_nm_per_rad * twist^2
    end
    return kinetic, elastic
end

function machine_shaft_state(specification::ModernMachineSpecification)
    angle = getfield.(specification.shaft_masses, :initial_angle_rad)
    speed = getfield.(specification.shaft_masses, :initial_speed_rad_s)
    kinetic, elastic = _machine_shaft_energy(specification, angle, speed)
    return MachineShaftState(
        angle,
        speed,
        zeros(length(specification.shaft_couplings)),
        kinetic,
        elastic,
        0.0,
        0.0,
    )
end

function _machine_shaft_workspace(specification::ModernMachineSpecification)
    mass_count = length(specification.shaft_masses)
    coupling_count = length(specification.shaft_couplings)
    return _MachineShaftWorkspace(
        zeros(mass_count),
        zeros(mass_count),
        zeros(mass_count),
        zeros(mass_count),
        zeros(mass_count),
        zeros(mass_count),
        zeros(coupling_count),
        zeros(coupling_count),
        zeros(mass_count),
    )
end

function _machine_terminal_rms_v(terminal_voltage_v::AbstractVector{<:Real})
    length(terminal_voltage_v) == 4 || throw(DimensionMismatch(
        "machine terminal voltage must contain three phases and neutral",
    ))
    neutral = Float64(terminal_voltage_v[4])
    squared_voltage_sum = 0.0
    for index in 1:3
        squared_voltage_sum += abs2(Float64(terminal_voltage_v[index]) - neutral)
    end
    return sqrt(squared_voltage_sum / 3.0)
end

function update_machine_control_state!(
    state::MachineControlState,
    specification::ModernMachineSpecification,
    inputs::MachinePortInputs,
    terminal_voltage_v::AbstractVector{<:Real},
    mechanical_speed_rad_s::Real,
    accepted_time_s::Real,
)
    controls = specification.controls
    time = Float64(accepted_time_s)
    speed = Float64(mechanical_speed_rad_s)
    isfinite(time) && isfinite(speed) || throw(ArgumentError(
        "machine control time and speed must be finite",
    ))
    if !inputs.control_enabled || !controls.enabled
        state.field_voltage_v = inputs.field_voltage_v
        state.mechanical_torque_nm = inputs.mechanical_torque_nm
        return false
    end
    tolerance = max(64.0 * eps(Float64) * max(abs(time), 1.0), 1.0e-14)
    time + tolerance >= state.next_task_time_s || return false
    sample_interval = controls.task_period_s
    voltage = _machine_terminal_rms_v(terminal_voltage_v)
    speed_error = inputs.speed_reference_rad_s - speed
    previous_sampled_speed = state.sampled_speed_rad_s
    speed_derivative = (speed - previous_sampled_speed) / sample_interval
    washout_decay = exp(-sample_interval / controls.stabilizer_washout_s)
    state.stabilizer_washout_state = washout_decay * state.stabilizer_washout_state +
        (1.0 - washout_decay) * speed_derivative
    lead_input = controls.stabilizer_gain * state.stabilizer_washout_state
    lag_decay = exp(-sample_interval / controls.stabilizer_lag_s)
    lead_fraction = controls.stabilizer_lead_s / controls.stabilizer_lag_s
    state.stabilizer_lead_lag_state = lag_decay * state.stabilizer_lead_lag_state +
        (1.0 - lag_decay) * lead_fraction * lead_input
    excitation_target = inputs.field_voltage_v + controls.excitation_gain * (
        inputs.voltage_reference_v - voltage + state.stabilizer_lead_lag_state
    )
    excitation_decay = exp(-sample_interval / controls.excitation_time_constant_s)
    unconstrained_excitation = excitation_decay * state.excitation_state_v +
        (1.0 - excitation_decay) * excitation_target
    field_voltage = clamp(
        unconstrained_excitation,
        controls.field_voltage_min_v,
        controls.field_voltage_max_v,
    )
    # Tracking anti-windup: the stored state is the physical limited output.
    state.excitation_state_v = field_voltage
    droop_torque = controls.governor_droop_rad_s_per_nm == 0.0 ? 0.0 :
        speed_error / controls.governor_droop_rad_s_per_nm
    torque_target = inputs.mechanical_torque_nm + droop_torque
    governor_decay = exp(-sample_interval / controls.governor_time_constant_s)
    unconstrained_torque = governor_decay * state.governor_state_nm +
        (1.0 - governor_decay) * torque_target
    torque = clamp(
        unconstrained_torque,
        controls.torque_min_nm,
        controls.torque_max_nm,
    )
    state.governor_state_nm = torque
    state.sensed_voltage_v = voltage
    state.sampled_speed_rad_s = speed
    state.field_voltage_v = field_voltage
    state.mechanical_torque_nm = torque
    state.field_limited = field_voltage != unconstrained_excitation
    state.torque_limited = torque != unconstrained_torque
    state.sample_count += 1
    while state.next_task_time_s <= time + tolerance
        state.next_task_time_s += controls.task_period_s
    end
    return true
end

function machine_shaft_derivative(
    specification::ModernMachineSpecification,
    angle_rad::Vector{Float64},
    speed_rad_s::Vector{Float64},
    electromagnetic_torque_nm::Float64,
    mechanical_torque_nm::Float64,
)
    mass_count = length(specification.shaft_masses)
    length(angle_rad) == mass_count && length(speed_rad_s) == mass_count ||
        throw(DimensionMismatch("machine shaft state does not match its specification"))
    torque = zeros(mass_count)
    coupling_torque = zeros(length(specification.shaft_couplings))
    electromagnetic_index = _machine_shaft_mass_index(
        specification,
        specification.electromagnetic_mass,
    )
    torque[electromagnetic_index] += electromagnetic_torque_nm
    torque[end] += mechanical_torque_nm
    damping_loss = 0.0
    for (index, mass) in enumerate(specification.shaft_masses)
        damping_torque = mass.damping_nm_s_per_rad * speed_rad_s[index]
        torque[index] -= damping_torque
        damping_loss += mass.damping_nm_s_per_rad * speed_rad_s[index]^2
    end
    for (index, coupling) in enumerate(specification.shaft_couplings)
        left = _machine_shaft_mass_index(specification, coupling.left_mass)
        right = _machine_shaft_mass_index(specification, coupling.right_mass)
        relative_angle = angle_rad[left] - angle_rad[right] - coupling.initial_twist_rad
        relative_speed = speed_rad_s[left] - speed_rad_s[right]
        branch_torque = coupling.stiffness_nm_per_rad * relative_angle +
            coupling.damping_nm_s_per_rad * relative_speed
        coupling_torque[index] = branch_torque
        torque[left] -= branch_torque
        torque[right] += branch_torque
        damping_loss += coupling.damping_nm_s_per_rad * relative_speed^2
    end
    acceleration = Float64[
        torque[index] / specification.shaft_masses[index].inertia_kg_m2 for
        index in 1:mass_count
    ]
    return speed_rad_s, acceleration, coupling_torque, damping_loss
end

function _machine_shaft_derivative!(
    acceleration_rad_s2::Vector{Float64},
    coupling_torque_nm::Vector{Float64},
    mass_torque_nm::Vector{Float64},
    specification::ModernMachineSpecification,
    angle_rad::Vector{Float64},
    speed_rad_s::Vector{Float64},
    electromagnetic_torque_nm::Float64,
    mechanical_torque_nm::Float64,
)
    mass_count = length(specification.shaft_masses)
    length(angle_rad) == mass_count && length(speed_rad_s) == mass_count &&
        length(acceleration_rad_s2) == mass_count &&
        length(mass_torque_nm) == mass_count || throw(DimensionMismatch(
            "machine shaft workspace does not match its specification",
        ))
    length(coupling_torque_nm) == length(specification.shaft_couplings) ||
        throw(DimensionMismatch("machine shaft coupling workspace has invalid length"))
    fill!(mass_torque_nm, 0.0)
    electromagnetic_index = _machine_shaft_mass_index(
        specification,
        specification.electromagnetic_mass,
    )
    mass_torque_nm[electromagnetic_index] += electromagnetic_torque_nm
    mass_torque_nm[end] += mechanical_torque_nm
    damping_loss_w = 0.0
    for index in eachindex(specification.shaft_masses)
        mass = specification.shaft_masses[index]
        damping_torque = mass.damping_nm_s_per_rad * speed_rad_s[index]
        mass_torque_nm[index] -= damping_torque
        damping_loss_w += mass.damping_nm_s_per_rad * speed_rad_s[index]^2
    end
    for index in eachindex(specification.shaft_couplings)
        coupling = specification.shaft_couplings[index]
        left = _machine_shaft_mass_index(specification, coupling.left_mass)
        right = _machine_shaft_mass_index(specification, coupling.right_mass)
        relative_angle = angle_rad[left] - angle_rad[right] - coupling.initial_twist_rad
        relative_speed = speed_rad_s[left] - speed_rad_s[right]
        branch_torque = coupling.stiffness_nm_per_rad * relative_angle +
            coupling.damping_nm_s_per_rad * relative_speed
        coupling_torque_nm[index] = branch_torque
        mass_torque_nm[left] -= branch_torque
        mass_torque_nm[right] += branch_torque
        damping_loss_w += coupling.damping_nm_s_per_rad * relative_speed^2
    end
    for index in eachindex(acceleration_rad_s2)
        acceleration_rad_s2[index] = mass_torque_nm[index] /
            specification.shaft_masses[index].inertia_kg_m2
    end
    return damping_loss_w
end

function _advance_machine_shaft!(
    state::MachineShaftState,
    specification::ModernMachineSpecification,
    previous_electromagnetic_torque_nm::Float64,
    candidate_electromagnetic_torque_nm::Float64,
    previous_mechanical_torque_nm::Float64,
    candidate_mechanical_torque_nm::Float64,
    step_s::Float64,
    workspace::_MachineShaftWorkspace,
)
    copyto!(workspace.old_angle_rad, state.angle_rad)
    copyto!(workspace.old_speed_rad_s, state.speed_rad_s)
    first_damping_loss_w = _machine_shaft_derivative!(
        workspace.first_acceleration_rad_s2,
        workspace.first_coupling_torque_nm,
        workspace.mass_torque_nm,
        specification,
        workspace.old_angle_rad,
        workspace.old_speed_rad_s,
        previous_electromagnetic_torque_nm,
        previous_mechanical_torque_nm,
    )
    for index in eachindex(workspace.predicted_angle_rad)
        workspace.predicted_angle_rad[index] = workspace.old_angle_rad[index] +
            step_s * workspace.old_speed_rad_s[index]
        workspace.predicted_speed_rad_s[index] = workspace.old_speed_rad_s[index] +
            step_s * workspace.first_acceleration_rad_s2[index]
    end
    candidate_damping_loss_w = _machine_shaft_derivative!(
        workspace.candidate_acceleration_rad_s2,
        workspace.candidate_coupling_torque_nm,
        workspace.mass_torque_nm,
        specification,
        workspace.predicted_angle_rad,
        workspace.predicted_speed_rad_s,
        candidate_electromagnetic_torque_nm,
        candidate_mechanical_torque_nm,
    )
    for index in eachindex(state.angle_rad)
        state.angle_rad[index] = workspace.old_angle_rad[index] + 0.5 * step_s * (
            workspace.old_speed_rad_s[index] + workspace.predicted_speed_rad_s[index]
        )
        state.speed_rad_s[index] = workspace.old_speed_rad_s[index] + 0.5 * step_s * (
            workspace.first_acceleration_rad_s2[index] +
            workspace.candidate_acceleration_rad_s2[index]
        )
    end
    for index in eachindex(state.coupling_torque_nm)
        state.coupling_torque_nm[index] = 0.5 * (
            workspace.first_coupling_torque_nm[index] +
            workspace.candidate_coupling_torque_nm[index]
        )
    end
    state.damping_loss_w = 0.5 * (first_damping_loss_w + candidate_damping_loss_w)
    kinetic, elastic = _machine_shaft_energy(
        specification,
        state.angle_rad,
        state.speed_rad_s,
    )
    old_momentum = 0.0
    new_momentum = 0.0
    damping_momentum_change = 0.0
    for index in eachindex(specification.shaft_masses)
        mass = specification.shaft_masses[index]
        old_momentum += mass.inertia_kg_m2 * workspace.old_speed_rad_s[index]
        new_momentum += mass.inertia_kg_m2 * state.speed_rad_s[index]
        damping_momentum_change += mass.damping_nm_s_per_rad *
            (workspace.old_speed_rad_s[index] + state.speed_rad_s[index])
    end
    expected_change = 0.5 * step_s * (
        previous_electromagnetic_torque_nm + candidate_electromagnetic_torque_nm +
        previous_mechanical_torque_nm + candidate_mechanical_torque_nm -
        damping_momentum_change
    )
    residual = abs((new_momentum - old_momentum) - expected_change)
    state.maximum_angular_momentum_residual_nms = max(
        state.maximum_angular_momentum_residual_nms,
        residual,
    )
    state.kinetic_energy_j = kinetic
    state.elastic_energy_j = elastic
    return state
end

function advance_machine_shaft!(
    state::MachineShaftState,
    specification::ModernMachineSpecification,
    previous_electromagnetic_torque_nm::Real,
    candidate_electromagnetic_torque_nm::Real,
    previous_mechanical_torque_nm::Real,
    candidate_mechanical_torque_nm::Real,
    step_s::Real,
)
    step = Float64(step_s)
    previous_torque_e = Float64(previous_electromagnetic_torque_nm)
    candidate_torque_e = Float64(candidate_electromagnetic_torque_nm)
    previous_torque_m = Float64(previous_mechanical_torque_nm)
    candidate_torque_m = Float64(candidate_mechanical_torque_nm)
    isfinite(step) && step > 0.0 || throw(ArgumentError(
        "machine shaft step must be finite and positive",
    ))
    old_angle = copy(state.angle_rad)
    old_speed = copy(state.speed_rad_s)
    derivative_0 = machine_shaft_derivative(
        specification,
        old_angle,
        old_speed,
        previous_torque_e,
        previous_torque_m,
    )
    predicted_angle = old_angle .+ step .* derivative_0[1]
    predicted_speed = old_speed .+ step .* derivative_0[2]
    derivative_1 = machine_shaft_derivative(
        specification,
        predicted_angle,
        predicted_speed,
        candidate_torque_e,
        candidate_torque_m,
    )
    state.angle_rad .= old_angle .+ 0.5 * step .* (derivative_0[1] .+ derivative_1[1])
    state.speed_rad_s .= old_speed .+ 0.5 * step .* (derivative_0[2] .+ derivative_1[2])
    state.coupling_torque_nm .= 0.5 .* (derivative_0[3] .+ derivative_1[3])
    state.damping_loss_w = 0.5 * (derivative_0[4] + derivative_1[4])
    kinetic, elastic = _machine_shaft_energy(
        specification,
        state.angle_rad,
        state.speed_rad_s,
    )
    old_momentum = sum(
        specification.shaft_masses[index].inertia_kg_m2 * old_speed[index] for
        index in eachindex(old_speed)
    )
    new_momentum = sum(
        specification.shaft_masses[index].inertia_kg_m2 * state.speed_rad_s[index] for
        index in eachindex(state.speed_rad_s)
    )
    expected_change = 0.5 * step * (
        previous_torque_e + candidate_torque_e +
        previous_torque_m + candidate_torque_m - sum(
            mass.damping_nm_s_per_rad *
                (old_speed[index] + state.speed_rad_s[index]) for
            (index, mass) in enumerate(specification.shaft_masses)
        )
    )
    residual = abs((new_momentum - old_momentum) - expected_change)
    state.maximum_angular_momentum_residual_nms = max(
        state.maximum_angular_momentum_residual_nms,
        residual,
    )
    state.kinetic_energy_j = kinetic
    state.elastic_energy_j = elastic
    return state
end

function advance_machine_shaft!(
    state::MachineShaftState,
    specification::ModernMachineSpecification,
    electromagnetic_torque_nm::Real,
    mechanical_torque_nm::Real,
    step_s::Real,
)
    return advance_machine_shaft!(
        state,
        specification,
        electromagnetic_torque_nm,
        electromagnetic_torque_nm,
        mechanical_torque_nm,
        mechanical_torque_nm,
        step_s,
    )
end

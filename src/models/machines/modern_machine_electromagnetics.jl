export MachineElectromagneticEvaluation,
       machine_phase_transform,
       machine_phase_to_rotor,
       machine_rotor_to_phase,
       machine_coenergy_current_hessian,
       machine_electromagnetic_evaluation,
       machine_electromagnetic_torque_nm,
       machine_magnetic_coenergy_j

struct MachineElectromagneticEvaluation
    current_a::Vector{Float64}
    current_flux_jacobian_per_h::Matrix{Float64}
    electromagnetic_torque_nm::Float64
    magnetic_coenergy_j::Float64
    copper_loss_w::Float64
    differential_inductance_margin_h::Float64
end

mutable struct _MachineElectromagneticWorkspace
    displaced_flux_wb::Vector{Float64}
    speed_matrix_per_s::Matrix{Float64}
    previous_voltage_v::Vector{Float64}
    candidate_voltage_v::Vector{Float64}
    previous_flux_derivative_wb_per_s::Vector{Float64}
    candidate_flux_derivative_wb_per_s::Vector{Float64}
    residual_wb::Vector{Float64}
    tangent::Matrix{Float64}
    terminal_input_map::Matrix{Float64}
    terminal_output_map::Matrix{Float64}
    terminal_jacobian_product::Matrix{Float64}
    phase_transform::Matrix{Float64}
    phase_voltage_v::Vector{Float64}
    rotor_voltage_v::Vector{Float64}
    phase_terminal_map::Matrix{Float64}
    tangent_pivots::Vector{Int}
    inverse_column::Vector{Float64}
end

function _machine_electromagnetic_workspace(preparation::ModernMachinePreparation)
    state_count = length(preparation.initial_flux_wb)
    return _MachineElectromagneticWorkspace(
        zeros(state_count),
        zeros(state_count, state_count),
        zeros(state_count),
        zeros(state_count),
        zeros(state_count),
        zeros(state_count),
        zeros(state_count),
        zeros(state_count, state_count),
        zeros(state_count, 4),
        zeros(4, state_count),
        zeros(state_count, 4),
        zeros(3, 3),
        zeros(3),
        zeros(3),
        zeros(3, 4),
        zeros(Int, state_count),
        zeros(state_count),
    )
end

function _solve_machine_lu_vector!(
    factors::Matrix{Float64},
    pivots::Vector{Int},
    right_hand_side::Vector{Float64},
)
    state_count = length(right_hand_side)
    for column in 1:state_count
        pivot = pivots[column]
        if pivot != column
            right_hand_side[column], right_hand_side[pivot] =
                right_hand_side[pivot], right_hand_side[column]
        end
    end
    for row in 2:state_count
        value = right_hand_side[row]
        for column in 1:(row - 1)
            value -= factors[row, column] * right_hand_side[column]
        end
        right_hand_side[row] = value
    end
    for row in state_count:-1:1
        value = right_hand_side[row]
        for column in (row + 1):state_count
            value -= factors[row, column] * right_hand_side[column]
        end
        right_hand_side[row] = value / factors[row, row]
    end
    return right_hand_side
end

function _solve_machine_lu_matrix!(
    factors::Matrix{Float64},
    pivots::Vector{Int},
    right_hand_side::Matrix{Float64},
)
    state_count = size(right_hand_side, 1)
    for factor_column in 1:state_count
        pivot = pivots[factor_column]
        if pivot != factor_column
            for right_column in axes(right_hand_side, 2)
                right_hand_side[factor_column, right_column],
                    right_hand_side[pivot, right_column] =
                    right_hand_side[pivot, right_column],
                    right_hand_side[factor_column, right_column]
            end
        end
    end
    for right_column in axes(right_hand_side, 2)
        for row in 2:state_count
            value = right_hand_side[row, right_column]
            for factor_column in 1:(row - 1)
                value -= factors[row, factor_column] *
                    right_hand_side[factor_column, right_column]
            end
            right_hand_side[row, right_column] = value
        end
        for row in state_count:-1:1
            value = right_hand_side[row, right_column]
            for factor_column in (row + 1):state_count
                value -= factors[row, factor_column] *
                    right_hand_side[factor_column, right_column]
            end
            right_hand_side[row, right_column] = value / factors[row, row]
        end
    end
    return right_hand_side
end

function _factor_machine_tangent!(
    tangent::Matrix{Float64},
    pivots::Vector{Int},
    inverse_column::Vector{Float64},
)
    state_count = size(tangent, 1)
    size(tangent, 2) == state_count && length(pivots) == state_count &&
        length(inverse_column) == state_count || throw(DimensionMismatch(
            "machine tangent factorization workspace has invalid dimensions",
        ))
    one_norm = 0.0
    infinity_norm = 0.0
    minimum_diagonal_dominance = Inf
    maximum_entry = 0.0
    for column in 1:state_count
        column_sum = 0.0
        for row in 1:state_count
            entry = abs(tangent[row, column])
            column_sum += entry
            maximum_entry = max(maximum_entry, entry)
        end
        one_norm = max(one_norm, column_sum)
    end
    for row in 1:state_count
        row_sum = 0.0
        off_diagonal_sum = 0.0
        for column in 1:state_count
            entry = abs(tangent[row, column])
            row_sum += entry
            row == column || (off_diagonal_sum += entry)
        end
        infinity_norm = max(infinity_norm, row_sum)
        minimum_diagonal_dominance = min(
            minimum_diagonal_dominance,
            abs(tangent[row, row]) - off_diagonal_sum,
        )
    end
    maximum_entry > 0.0 || return 0.0
    pivot_floor = 64.0 * eps(Float64) * maximum_entry
    for column in 1:state_count
        pivot_row = column
        pivot_magnitude = abs(tangent[column, column])
        for row in (column + 1):state_count
            candidate_magnitude = abs(tangent[row, column])
            if candidate_magnitude > pivot_magnitude
                pivot_magnitude = candidate_magnitude
                pivot_row = row
            end
        end
        pivot_magnitude > pivot_floor || return 0.0
        pivots[column] = pivot_row
        if pivot_row != column
            for entry_column in 1:state_count
                tangent[column, entry_column], tangent[pivot_row, entry_column] =
                    tangent[pivot_row, entry_column], tangent[column, entry_column]
            end
        end
        pivot_value = tangent[column, column]
        for row in (column + 1):state_count
            tangent[row, column] /= pivot_value
            multiplier = tangent[row, column]
            for entry_column in (column + 1):state_count
                tangent[row, entry_column] -= multiplier * tangent[column, entry_column]
            end
        end
    end
    if minimum_diagonal_dominance > 0.0
        return minimum_diagonal_dominance / infinity_norm
    end
    inverse_norm = 0.0
    for column in 1:state_count
        fill!(inverse_column, 0.0)
        inverse_column[column] = 1.0
        _solve_machine_lu_vector!(tangent, pivots, inverse_column)
        column_sum = 0.0
        for value in inverse_column
            column_sum += abs(value)
        end
        inverse_norm = max(inverse_norm, column_sum)
    end
    return inv(one_norm * inverse_norm)
end

function _machine_phase_transform!(transform::Matrix{Float64}, electrical_angle_rad::Float64)
    size(transform) == (3, 3) || throw(DimensionMismatch(
        "machine phase-transform workspace must be 3x3",
    ))
    phase_shift = 2.0 * pi / 3.0
    scale = sqrt(2.0 / 3.0)
    zero_scale = inv(sqrt(3.0))
    for column in 1:3
        transform[1, column] = zero_scale
    end
    transform[2, 1] = scale * cos(electrical_angle_rad)
    transform[2, 2] = scale * cos(electrical_angle_rad - phase_shift)
    transform[2, 3] = scale * cos(electrical_angle_rad + phase_shift)
    transform[3, 1] = -scale * sin(electrical_angle_rad)
    transform[3, 2] = -scale * sin(electrical_angle_rad - phase_shift)
    transform[3, 3] = -scale * sin(electrical_angle_rad + phase_shift)
    return transform
end

"""Power-invariant zero/d/q transform with q positive against increasing electrical angle."""
function machine_phase_transform(electrical_angle_rad::Real)
    angle = Float64(electrical_angle_rad)
    isfinite(angle) || throw(ArgumentError("machine electrical angle must be finite"))
    phase_shift = 2.0 * pi / 3.0
    scale = sqrt(2.0 / 3.0)
    zero_scale = inv(sqrt(3.0))
    return [
        zero_scale zero_scale zero_scale
        scale * cos(angle) scale * cos(angle - phase_shift) scale * cos(angle + phase_shift)
        -scale * sin(angle) -scale * sin(angle - phase_shift) -scale * sin(angle + phase_shift)
    ]
end

function machine_phase_to_rotor(phase_values, electrical_angle_rad::Real)
    values = phase_values isa AbstractVector{Float64} ?
        phase_values : Float64.(phase_values)
    length(values) == 3 || throw(DimensionMismatch(
        "machine phase vector must contain exactly three values",
    ))
    all(isfinite, values) || throw(ArgumentError("machine phase values must be finite"))
    return machine_phase_transform(electrical_angle_rad) * values
end

function machine_rotor_to_phase(rotor_values, electrical_angle_rad::Real)
    values = rotor_values isa AbstractVector{Float64} ?
        rotor_values : Float64.(rotor_values)
    length(values) == 3 || throw(DimensionMismatch(
        "machine zero/d/q vector must contain exactly three values",
    ))
    all(isfinite, values) || throw(ArgumentError("machine rotor values must be finite"))
    return transpose(machine_phase_transform(electrical_angle_rad)) * values
end

function _machine_saturation_terms(
    law::MachineMagneticCoenergyLaw,
    d_axis_flux_wb::Float64,
    q_axis_flux_wb::Float64,
)
    radial = law.radial_coefficient_per_wb2_h
    cross = law.cross_coefficient_per_wb2_h
    d2 = d_axis_flux_wb^2
    q2 = q_axis_flux_wb^2
    radius2 = d2 + q2
    radius = sqrt(radius2)
    radius <= law.maximum_flux_wb || throw(ModernMachineRefusal(
        :flux_domain_exceeded,
        :electromagnetic_evaluation,
        :unknown,
        WoundFieldSynchronousMachine,
        "machine flux magnitude exceeds the declared coenergy domain",
        (flux_magnitude_wb=radius, maximum_flux_wb=law.maximum_flux_wb),
    ))
    current_d = radial * radius2 * d_axis_flux_wb +
        cross * d_axis_flux_wb * q2
    current_q = radial * radius2 * q_axis_flux_wb +
        cross * q_axis_flux_wb * d2
    hessian_dd = radial * (3.0 * d2 + q2) + cross * q2
    hessian_qq = radial * (d2 + 3.0 * q2) + cross * d2
    hessian_dq = 2.0 * (radial + cross) * d_axis_flux_wb * q_axis_flux_wb
    coenergy = 0.25 * radial * radius2^2 + 0.5 * cross * d2 * q2
    return current_d, current_q, hessian_dd, hessian_qq, hessian_dq, coenergy
end

"""Return current, symmetric differential current/flux map, and convex coenergy."""
function machine_coenergy_current_hessian(
    preparation::ModernMachinePreparation,
    flux_wb,
)
    flux = flux_wb isa AbstractVector{Float64} ? flux_wb : Float64.(flux_wb)
    length(flux) == length(preparation.initial_flux_wb) || throw(DimensionMismatch(
        "machine flux vector does not match its preparation",
    ))
    all(isfinite, flux) || throw(ArgumentError("machine flux vector must be finite"))
    displaced_flux = flux - preparation.permanent_flux_offset_wb
    current = preparation.inverse_inductance_per_h * displaced_flux
    linear_coenergy = 0.5 * dot(displaced_flux, current)
    hessian = copy(preparation.inverse_inductance_per_h)
    layout = preparation.layout
    d_index = layout.stator_d_index
    q_index = layout.stator_q_index
    d_axis_flux = displaced_flux[d_index]
    q_axis_flux = displaced_flux[q_index]
    saturation = _machine_saturation_terms(
        preparation.specification.saturation,
        d_axis_flux,
        q_axis_flux,
    )
    current[d_index] += saturation[1]
    current[q_index] += saturation[2]
    hessian[d_index, d_index] += saturation[3]
    hessian[q_index, q_index] += saturation[4]
    hessian[d_index, q_index] += saturation[5]
    hessian[q_index, d_index] += saturation[5]
    all(isfinite, hessian) || throw(ArgumentError(
        "machine differential current/flux map must remain finite",
    ))
    # The prepared linear map is positive definite. The admitted quartic
    # saturation law has a positive-semidefinite Hessian because its
    # nonnegative cross coefficient is bounded by twice its radial coefficient.
    # The inverse maximum absolute row sum is therefore a conservative lower
    # bound on differential inductance, without a trial-time eigendecomposition.
    maximum_absolute_row_sum = 0.0
    for row in axes(hessian, 1)
        absolute_row_sum = 0.0
        for column in axes(hessian, 2)
            absolute_row_sum += abs(hessian[row, column])
        end
        maximum_absolute_row_sum = max(maximum_absolute_row_sum, absolute_row_sum)
    end
    maximum_absolute_row_sum > 0.0 || throw(ArgumentError(
        "machine differential current/flux map must remain positive definite",
    ))
    differential_inductance_margin = inv(maximum_absolute_row_sum)
    return current, hessian, linear_coenergy + saturation[6],
        differential_inductance_margin
end

function machine_electromagnetic_torque_nm(
    preparation::ModernMachinePreparation,
    flux_wb,
    current_a,
)
    flux = flux_wb isa AbstractVector{Float64} ? flux_wb : Float64.(flux_wb)
    current = current_a isa AbstractVector{Float64} ? current_a : Float64.(current_a)
    layout = preparation.layout
    return preparation.specification.pole_pairs * (
        flux[layout.stator_d_index] * current[layout.stator_q_index] -
        flux[layout.stator_q_index] * current[layout.stator_d_index]
    )
end

function machine_magnetic_coenergy_j(
    preparation::ModernMachinePreparation,
    flux_wb,
)
    return machine_coenergy_current_hessian(preparation, flux_wb)[3]
end

function machine_electromagnetic_evaluation(
    preparation::ModernMachinePreparation,
    flux_wb,
)
    current, hessian, coenergy, inductance_margin =
        machine_coenergy_current_hessian(preparation, flux_wb)
    torque = machine_electromagnetic_torque_nm(preparation, flux_wb, current)
    copper_loss = dot(current, preparation.resistance_ohm .* current)
    copper_loss >= -64.0 * eps(Float64) || throw(ArgumentError(
        "machine winding loss must be passive",
    ))
    return MachineElectromagneticEvaluation(
        current,
        hessian,
        torque,
        coenergy,
        max(copper_loss, 0.0),
        inductance_margin,
    )
end

function _machine_coenergy_current_hessian!(
    current_a::Vector{Float64},
    hessian_per_h::Matrix{Float64},
    displaced_flux_wb::Vector{Float64},
    preparation::ModernMachinePreparation,
    flux_wb::Vector{Float64},
)
    state_count = length(preparation.initial_flux_wb)
    length(flux_wb) == state_count && length(current_a) == state_count &&
        length(displaced_flux_wb) == state_count &&
        size(hessian_per_h) == (state_count, state_count) ||
        throw(DimensionMismatch("machine electromagnetic workspace has invalid dimensions"))
    for index in eachindex(flux_wb)
        isfinite(flux_wb[index]) || throw(ArgumentError(
            "machine flux vector must be finite",
        ))
        displaced_flux_wb[index] =
            flux_wb[index] - preparation.permanent_flux_offset_wb[index]
    end
    _machine_matrix_vector_product!(
        current_a,
        preparation.inverse_inductance_per_h,
        displaced_flux_wb,
    )
    linear_coenergy_j = 0.5 * dot(displaced_flux_wb, current_a)
    copyto!(hessian_per_h, preparation.inverse_inductance_per_h)
    layout = preparation.layout
    d_index = layout.stator_d_index
    q_index = layout.stator_q_index
    saturation = _machine_saturation_terms(
        preparation.specification.saturation,
        displaced_flux_wb[d_index],
        displaced_flux_wb[q_index],
    )
    current_a[d_index] += saturation[1]
    current_a[q_index] += saturation[2]
    hessian_per_h[d_index, d_index] += saturation[3]
    hessian_per_h[q_index, q_index] += saturation[4]
    hessian_per_h[d_index, q_index] += saturation[5]
    hessian_per_h[q_index, d_index] += saturation[5]
    maximum_absolute_row_sum = 0.0
    for row in axes(hessian_per_h, 1)
        absolute_row_sum = 0.0
        for column in axes(hessian_per_h, 2)
            entry = hessian_per_h[row, column]
            isfinite(entry) || throw(ArgumentError(
                "machine differential current/flux map must remain finite",
            ))
            absolute_row_sum += abs(entry)
        end
        maximum_absolute_row_sum = max(maximum_absolute_row_sum, absolute_row_sum)
    end
    maximum_absolute_row_sum > 0.0 || throw(ArgumentError(
        "machine differential current/flux map must remain positive definite",
    ))
    electromagnetic_torque_nm = preparation.specification.pole_pairs * (
        flux_wb[d_index] * current_a[q_index] -
        flux_wb[q_index] * current_a[d_index]
    )
    copper_loss_w = 0.0
    for index in eachindex(current_a)
        copper_loss_w += preparation.resistance_ohm[index] * current_a[index]^2
    end
    copper_loss_w >= -64.0 * eps(Float64) || throw(ArgumentError(
        "machine winding loss must be passive",
    ))
    return (
        electromagnetic_torque_nm=electromagnetic_torque_nm,
        magnetic_coenergy_j=linear_coenergy_j + saturation[6],
        copper_loss_w=max(copper_loss_w, 0.0),
        differential_inductance_margin_h=inv(maximum_absolute_row_sum),
    )
end

function _machine_matrix_vector_product!(
    output::Vector{Float64},
    matrix::Matrix{Float64},
    input::Vector{Float64},
)
    size(matrix, 1) == length(output) && size(matrix, 2) == length(input) ||
        throw(DimensionMismatch("machine matrix-vector product dimensions differ"))
    @inbounds for row in axes(matrix, 1)
        value = 0.0
        for column in axes(matrix, 2)
            value += matrix[row, column] * input[column]
        end
        output[row] = value
    end
    return output
end

function _machine_matrix_matrix_product!(
    output::Matrix{Float64},
    left::Matrix{Float64},
    right::Matrix{Float64},
)
    size(left, 1) == size(output, 1) &&
        size(right, 2) == size(output, 2) &&
        size(left, 2) == size(right, 1) || throw(DimensionMismatch(
        "machine matrix product dimensions differ",
    ))
    @inbounds for column in axes(output, 2)
        for row in axes(output, 1)
            value = 0.0
            for inner in axes(left, 2)
                value += left[row, inner] * right[inner, column]
            end
            output[row, column] = value
        end
    end
    return output
end

function _machine_electrical_speed_matrix!(
    speed_matrix_per_s::Matrix{Float64},
    preparation::ModernMachinePreparation,
    mechanical_speed_rad_s::Float64,
)
    fill!(speed_matrix_per_s, 0.0)
    electrical_speed_rad_s =
        preparation.specification.pole_pairs * mechanical_speed_rad_s
    layout = preparation.layout
    speed_matrix_per_s[layout.stator_d_index, layout.stator_q_index] =
        electrical_speed_rad_s
    speed_matrix_per_s[layout.stator_q_index, layout.stator_d_index] =
        -electrical_speed_rad_s
    return speed_matrix_per_s
end

function _machine_voltage_vector!(
    voltage_v::Vector{Float64},
    preparation::ModernMachinePreparation,
    terminal_voltage_v::AbstractVector{<:Real},
    electrical_angle_rad::Float64,
    field_voltage_v::Float64,
    rotor_voltage_dq_v::NTuple{2,Float64},
    workspace::_MachineElectromagneticWorkspace,
)
    length(terminal_voltage_v) == 4 || throw(DimensionMismatch(
        "modern machine runtime requires phase-a, phase-b, phase-c, and neutral voltages",
    ))
    fill!(workspace.phase_voltage_v, 0.0)
    for terminal in eachindex(terminal_voltage_v)
        terminal_value = Float64(terminal_voltage_v[terminal])
        isfinite(terminal_value) || throw(ArgumentError(
            "modern machine terminal voltage must be finite",
        ))
        for phase in axes(preparation.terminal_voltage_map, 1)
            workspace.phase_voltage_v[phase] +=
                preparation.terminal_voltage_map[phase, terminal] * terminal_value
        end
    end
    _machine_phase_transform!(workspace.phase_transform, electrical_angle_rad)
    _machine_matrix_vector_product!(
        workspace.rotor_voltage_v,
        workspace.phase_transform,
        workspace.phase_voltage_v,
    )
    fill!(voltage_v, 0.0)
    layout = preparation.layout
    voltage_v[layout.zero_index] = workspace.rotor_voltage_v[1]
    voltage_v[layout.stator_d_index] = workspace.rotor_voltage_v[2]
    voltage_v[layout.stator_q_index] = workspace.rotor_voltage_v[3]
    layout.field_index === nothing || (voltage_v[layout.field_index] = field_voltage_v)
    for branch_index in eachindex(layout.rotor_d_indices)
        preparation.specification.rotor_branches[branch_index].terminal_exposed || continue
        voltage_v[layout.rotor_d_indices[branch_index]] = rotor_voltage_dq_v[1]
        voltage_v[layout.rotor_q_indices[branch_index]] = rotor_voltage_dq_v[2]
    end
    return voltage_v
end

function _machine_terminal_input_matrix!(
    input_map::Matrix{Float64},
    preparation::ModernMachinePreparation,
    electrical_angle_rad::Float64,
    workspace::_MachineElectromagneticWorkspace,
)
    _machine_phase_transform!(workspace.phase_transform, electrical_angle_rad)
    _machine_matrix_matrix_product!(
        workspace.phase_terminal_map,
        workspace.phase_transform,
        preparation.terminal_voltage_map,
    )
    fill!(input_map, 0.0)
    layout = preparation.layout
    for terminal in axes(input_map, 2)
        input_map[layout.zero_index, terminal] = workspace.phase_terminal_map[1, terminal]
        input_map[layout.stator_d_index, terminal] =
            workspace.phase_terminal_map[2, terminal]
        input_map[layout.stator_q_index, terminal] =
            workspace.phase_terminal_map[3, terminal]
    end
    return input_map
end

function _machine_terminal_output_matrix!(
    output_map::Matrix{Float64},
    preparation::ModernMachinePreparation,
    electrical_angle_rad::Float64,
    workspace::_MachineElectromagneticWorkspace,
)
    _machine_phase_transform!(workspace.phase_transform, electrical_angle_rad)
    fill!(output_map, 0.0)
    layout = preparation.layout
    state_indices = (layout.zero_index, layout.stator_d_index, layout.stator_q_index)
    for axis in 1:3
        state_index = state_indices[axis]
        for terminal in axes(output_map, 1)
            coefficient = 0.0
            for phase in 1:3
                coefficient += preparation.terminal_current_map[terminal, phase] *
                    workspace.phase_transform[axis, phase]
            end
            output_map[terminal, state_index] = coefficient
        end
    end
    return output_map
end

function _machine_flux_derivative!(
    derivative_wb_per_s::Vector{Float64},
    current_a::Vector{Float64},
    current_flux_jacobian_per_h::Matrix{Float64},
    preparation::ModernMachinePreparation,
    flux_wb::Vector{Float64},
    voltage_v::Vector{Float64},
    speed_matrix_per_s::Matrix{Float64},
    workspace::_MachineElectromagneticWorkspace,
)
    evaluation = _machine_coenergy_current_hessian!(
        current_a,
        current_flux_jacobian_per_h,
        workspace.displaced_flux_wb,
        preparation,
        flux_wb,
    )
    _machine_matrix_vector_product!(derivative_wb_per_s, speed_matrix_per_s, flux_wb)
    for index in eachindex(derivative_wb_per_s)
        derivative_wb_per_s[index] += voltage_v[index] -
            preparation.resistance_ohm[index] * current_a[index]
    end
    return evaluation
end

function _machine_electrical_speed_matrix(
    preparation::ModernMachinePreparation,
    mechanical_speed_rad_s::Float64,
)
    state_count = length(preparation.initial_flux_wb)
    matrix = zeros(state_count, state_count)
    electrical_speed = preparation.specification.pole_pairs * mechanical_speed_rad_s
    layout = preparation.layout
    matrix[layout.stator_d_index, layout.stator_q_index] = electrical_speed
    matrix[layout.stator_q_index, layout.stator_d_index] = -electrical_speed
    return matrix
end

function _machine_voltage_vector(
    preparation::ModernMachinePreparation,
    terminal_voltage_v::AbstractVector{<:Real},
    electrical_angle_rad::Float64,
    field_voltage_v::Float64,
    rotor_voltage_dq_v::NTuple{2,Float64},
)
    length(terminal_voltage_v) == 4 || throw(DimensionMismatch(
        "modern machine runtime requires phase-a, phase-b, phase-c, and neutral voltages",
    ))
    terminal_voltage = terminal_voltage_v isa AbstractVector{Float64} ?
        terminal_voltage_v : Float64.(terminal_voltage_v)
    all(isfinite, terminal_voltage) || throw(ArgumentError(
        "modern machine terminal voltage must be finite",
    ))
    phase_voltage = preparation.terminal_voltage_map * terminal_voltage
    rotor_voltage = machine_phase_to_rotor(phase_voltage, electrical_angle_rad)
    voltage = zeros(length(preparation.initial_flux_wb))
    layout = preparation.layout
    voltage[layout.zero_index] = rotor_voltage[1]
    voltage[layout.stator_d_index] = rotor_voltage[2]
    voltage[layout.stator_q_index] = rotor_voltage[3]
    layout.field_index === nothing || (voltage[layout.field_index] = field_voltage_v)
    if !isempty(layout.rotor_d_indices)
        exposed = findfirst(getfield.(preparation.specification.rotor_branches, :terminal_exposed))
        if exposed !== nothing
            voltage[layout.rotor_d_indices[exposed]] = rotor_voltage_dq_v[1]
            voltage[layout.rotor_q_indices[exposed]] = rotor_voltage_dq_v[2]
        end
    end
    return voltage
end

function _machine_terminal_input_matrix(
    preparation::ModernMachinePreparation,
    electrical_angle_rad::Float64,
)
    state_count = length(preparation.initial_flux_wb)
    matrix = zeros(state_count, 4)
    transform = machine_phase_transform(electrical_angle_rad)
    phase_terminal_map = transform * preparation.terminal_voltage_map
    layout = preparation.layout
    matrix[layout.zero_index, :] .= phase_terminal_map[1, :]
    matrix[layout.stator_d_index, :] .= phase_terminal_map[2, :]
    matrix[layout.stator_q_index, :] .= phase_terminal_map[3, :]
    return matrix
end

function _machine_terminal_output_matrix(
    preparation::ModernMachinePreparation,
    electrical_angle_rad::Float64,
)
    axis_selector = zeros(3, length(preparation.initial_flux_wb))
    layout = preparation.layout
    axis_selector[1, layout.zero_index] = 1.0
    axis_selector[2, layout.stator_d_index] = 1.0
    axis_selector[3, layout.stator_q_index] = 1.0
    return preparation.terminal_current_map *
        transpose(machine_phase_transform(electrical_angle_rad)) * axis_selector
end

function _machine_flux_derivative(
    preparation::ModernMachinePreparation,
    flux_wb::Vector{Float64},
    voltage_v::Vector{Float64},
    speed_matrix::Matrix{Float64},
)
    evaluation = machine_electromagnetic_evaluation(preparation, flux_wb)
    derivative = voltage_v - preparation.resistance_ohm .* evaluation.current_a +
        speed_matrix * flux_wb
    return derivative, evaluation
end

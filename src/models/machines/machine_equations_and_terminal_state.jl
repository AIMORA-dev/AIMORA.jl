
function universal_machine_postsolve_update!(
    state::UniversalMachinePostsolveState;
    kwargs...,
)
    preview = universal_machine_postsolve_update_preview(
        ;
        coil_parameters = state.coil_parameters,
        source_crests = state.source_crests,
        current_values = state.current_values,
        prediction_values = state.prediction_values,
        history_values = state.history_values,
        d_axis_flux = state.d_axis_flux,
        q_axis_flux = state.q_axis_flux,
        theta_electric_rad = state.theta_electric_rad,
        machine_type = state.machine_type,
        input_mode = state.input_mode,
        prediction_loop_marker = state.prediction_loop_marker,
        kwargs...,
    )
    empty!(state.coil_parameters)
    empty!(state.source_crests)
    empty!(state.current_values)
    empty!(state.prediction_values)
    empty!(state.history_values)
    empty!(state.prediction_report_text_lines)
    append!(state.coil_parameters, preview.coil_parameters)
    append!(state.source_crests, preview.source_crests)
    append!(state.current_values, preview.current_values)
    append!(state.prediction_values, preview.prediction_values)
    append!(state.history_values, preview.history_values)
    append!(state.prediction_report_text_lines, preview.prediction_report_text_lines)
    state.d_axis_flux = preview.d_axis_flux
    state.q_axis_flux = preview.q_axis_flux
    state.theta_electric_rad = preview.theta_electric_rad
    state.machine_type = preview.machine_type
    state.input_mode = preview.input_mode
    state.prediction_loop_marker = preview.final_prediction_loop_marker
    state.postsolve_update_mutated = true
    state.prediction_report_text_mutated = preview.prediction_report_text_formatting
    return merge(
        preview,
        (
            coil_parameters = copy(state.coil_parameters),
            source_crests = copy(state.source_crests),
            current_values = copy(state.current_values),
            prediction_values = copy(state.prediction_values),
            history_values = copy(state.history_values),
            prediction_report_text_lines = copy(state.prediction_report_text_lines),
            prediction_report_text_line_count = length(state.prediction_report_text_lines),
            d_axis_flux = state.d_axis_flux,
            q_axis_flux = state.q_axis_flux,
            theta_electric_rad = state.theta_electric_rad,
            machine_type = state.machine_type,
            input_mode = state.input_mode,
            final_prediction_loop_marker = state.prediction_loop_marker,
            postsolve_update_mutated = state.postsolve_update_mutated,
            prediction_report_text_mutated = state.prediction_report_text_mutated,
        ),
    )
end

function _machine_matrix_values(name::String, values)
    matrix = Matrix{Float64}(values)
    all(isfinite, matrix) || throw(ArgumentError("$name entries must be finite"))
    return matrix
end

const _UNIVERSAL_MACHINE_REPORT_SEPARATOR = "0-------------------------------------------------"

function _machine_fortran_i(value::Int, width::Int)
    text = string(value)
    length(text) <= width || return repeat("*", width)
    return lpad(text, width)
end

function _machine_fortran_e14_5(value::Float64)
    isfinite(value) || throw(ArgumentError("Fortran E14.5 field requires a finite value"))
    negative = signbit(value)
    magnitude = abs(value)
    exponent = 0
    mantissa = 0.0
    if magnitude != 0.0
        exponent = floor(Int, log10(magnitude)) + 1
        mantissa = magnitude / 10.0^exponent
    end
    mantissa_text = @sprintf("%.5f", mantissa)
    if startswith(mantissa_text, "1.")
        exponent += 1
        mantissa_text = "0.10000"
    end
    exponent_text = lpad(string(abs(exponent)), 2, '0')
    body = string(
        negative ? "-" : "",
        mantissa_text,
        "E",
        exponent < 0 ? "-" : "+",
        exponent_text,
    )
    length(body) <= 14 || return repeat("*", 14)
    return lpad(body, 14)
end

function _machine_fortran_i_list(values::AbstractVector{Int}, width::Int)
    return join((_machine_fortran_i(value, width) for value in values), "")
end

function _machine_fortran_e14_5_list(values)
    return join((_machine_fortran_e14_5(Float64(value)) for value in values), "")
end

function _universal_machine_prediction_report_text_lines(
    machine_index::Int,
    prediction_start_index::Int,
    prediction_end_index::Int,
    prediction_values,
)
    values = _machine_real_vector("prediction_values", prediction_values, 3)
    return String[
        "",
        " ************************************** PREDICTED A,B,C POWER VOLTAGES/RESISTANCE FOR NEXT TIME-STEP OF UM NUMBER" *
        _machine_fortran_i(machine_index, 4) *
        " :",
        "      UMCURP(" *
        _machine_fortran_i(prediction_start_index, 3) *
        ":" *
        _machine_fortran_i(prediction_end_index, 3) *
        ") :" *
        "   " *
        _machine_fortran_e14_5(values[1]) *
        "   " *
        _machine_fortran_e14_5(values[2]) *
        "   " *
        _machine_fortran_e14_5(values[3]),
    ]
end

function _universal_machine_matrix_report_lines(
    header::String,
    continuation_indent::Int,
    values::AbstractVector{Float64},
)
    length(values) == 9 || throw(ArgumentError("machine report matrix text requires nine values"))
    matrix = reshape(values, 3, 3)
    lines = String[
        header * "   " *
        _machine_fortran_e14_5_list([matrix[1, 1], matrix[1, 2], matrix[1, 3]]),
    ]
    for row in 2:3
        push!(
            lines,
            repeat(" ", continuation_indent) *
            _machine_fortran_e14_5_list([matrix[row, 1], matrix[row, 2], matrix[row, 3]]),
        )
    end
    return lines
end

function _universal_machine_report_text_lines(report)
    lines = String[]
    report.report_block_executed || return lines

    if report.initial_header_candidate == 0 && report.diagnostic_level >= 1
        push!(lines, _UNIVERSAL_MACHINE_REPORT_SEPARATOR)
    end
    if report.time_header_written
        push!(lines, _UNIVERSAL_MACHINE_REPORT_SEPARATOR)
        push!(
            lines,
            "0TIME =" *
            _machine_fortran_e14_5(report.time_s) *
            "***************************  NCOMP =" *
            _machine_fortran_i(report.component_count, 3),
        )
    end
    report.test_output_header_written && push!(lines, "0*** TEST OUTPUT :")
    if report.current_1_3_written
        push!(
            lines,
            "0UMCUR(1)=" *
            _machine_fortran_e14_5_list(report.current_report_values[1:1]) *
            "  UMCUR(2)=" *
            _machine_fortran_e14_5_list(report.current_report_values[2:2]) *
            "  UMCUR(3)=" *
            _machine_fortran_e14_5_list(report.current_report_values[3:3]),
        )
    end
    if report.current_4_6_written
        push!(
            lines,
            " UMCUR(4)=" *
            _machine_fortran_e14_5_list(report.current_report_values[4:4]) *
            "  UMCUR(5)=" *
            _machine_fortran_e14_5_list(report.current_report_values[5:5]) *
            "  UMCUR(6)=" *
            _machine_fortran_e14_5_list(report.current_report_values[6:6]),
        )
    end
    if report.flux_iteration_written
        push!(
            lines,
            "     FLXD=" *
            _machine_fortran_e14_5(report.flux_report_values[1]) *
            "      FLXQ=" *
            _machine_fortran_e14_5(report.flux_report_values[2]) *
            "    NITROM=" *
            _machine_fortran_i(round(Int, report.flux_report_values[3]), 3),
        )
    end
    if report.node_table_written
        length(report.first_terminal_report_nodes) <= 15 ||
            throw(ArgumentError("bounded universal-machine report text supports up to fifteen terminal nodes"))
        push!(lines, "0NODVO1 :   " * _machine_fortran_i_list(report.first_terminal_report_nodes, 4))
        push!(lines, " NODVO2 :   " * _machine_fortran_i_list(report.second_terminal_report_nodes, 4))
        push!(lines, " NODOM  :   " * _machine_fortran_i(report.mechanical_report_node, 4))
    end
    if report.thevenin_scalar_written
        push!(lines, "0ZTHEVM = " * _machine_fortran_e14_5(report.mechanical_thevenin_resistance))
    end
    if report.rotor_matrix_written
        append!(lines, _universal_machine_matrix_report_lines(
            "0ZTHEVR(3,3):",
            16,
            report.rotor_thevenin_values,
        ))
    end
    if report.stator_matrix_written
        append!(lines, _universal_machine_matrix_report_lines(
            "0ZTHS3(3,3):",
            15,
            report.stator_thevenin_values,
        ))
    end
    return lines
end

function universal_machine_report_schedule_preview(;
    start_step_index::Int,
    machine_index::Int,
    machine_count::Int,
    last_report_step::Int=0,
    output_interval_steps::Int=0,
    diagnostic_level::Int=0,
    time_s::Real=0.0,
    component_count::Int=0,
    current_values::AbstractVector{<:Real}=Float64[],
    include_extra_power_coil_currents::Bool=false,
    d_axis_flux::Real=0.0,
    q_axis_flux::Real=0.0,
    iteration_count::Int=0,
    first_terminal_nodes::AbstractVector{<:Integer}=Int[],
    second_terminal_nodes::AbstractVector{<:Integer}=Int[],
    mechanical_node::Int=0,
    mechanical_thevenin_resistance::Real=0.0,
    rotor_thevenin_matrix::AbstractMatrix{<:Real}=zeros(0, 0),
    stator_thevenin_matrix::AbstractMatrix{<:Real}=zeros(0, 0),
    emit_text::Bool=false,
)
    machine_index > 0 || throw(ArgumentError("machine_index must be positive"))
    machine_count >= machine_index || throw(ArgumentError("machine_count must include machine_index"))
    output_interval_steps >= 0 || throw(ArgumentError("output_interval_steps must be nonnegative"))
    diagnostic_level >= 0 || throw(ArgumentError("diagnostic_level must be nonnegative"))
    currents = _machine_float_vector(current_values)
    if diagnostic_level >= 3
        length(currents) >= 3 || throw(ArgumentError("current_values must contain at least three entries"))
        if include_extra_power_coil_currents
            length(currents) >= 6 || throw(ArgumentError("current_values must contain six entries"))
        end
    end
    node1 = Int[Int(value) for value in first_terminal_nodes]
    node2 = Int[Int(value) for value in second_terminal_nodes]
    length(node1) == length(node2) ||
        throw(ArgumentError("terminal node vectors must have matching lengths"))
    rotor_matrix = _machine_matrix_values("rotor_thevenin_matrix", rotor_thevenin_matrix)
    stator_matrix = _machine_matrix_values("stator_thevenin_matrix", stator_thevenin_matrix)
    if diagnostic_level >= 2
        isempty(rotor_matrix) || size(rotor_matrix) == (3, 3) ||
            throw(ArgumentError("rotor_thevenin_matrix must be 3x3 when supplied"))
        isempty(stator_matrix) || size(stator_matrix) == (3, 3) ||
            throw(ArgumentError("stator_thevenin_matrix must be 3x3 when supplied"))
    end
    time = Float64(time_s)
    d_flux = Float64(d_axis_flux)
    q_flux = Float64(q_axis_flux)
    zthevm = Float64(mechanical_thevenin_resistance)
    all(isfinite, (time, d_flux, q_flux, zthevm)) ||
        throw(ArgumentError("report schedule scalar inputs must be finite"))

    initial_header_candidate = start_step_index + machine_index - 1
    report_skipped_to_postsolve = false
    report_block_executed = true
    next_due_step = last_report_step + output_interval_steps
    separator_write_count = 0
    time_header_written = false

    if initial_header_candidate == 0
        diagnostic_level >= 1 && (separator_write_count += 1)
        if diagnostic_level >= 2
            separator_write_count += 1
            time_header_written = true
        end
    elseif start_step_index != 0
        if output_interval_steps == 0 || start_step_index != next_due_step
            report_skipped_to_postsolve = true
            report_block_executed = false
        elseif machine_index == 1 && diagnostic_level >= 2
            separator_write_count += 1
            time_header_written = true
        end
    end

    final_last_report_step =
        report_block_executed && machine_index == machine_count ?
        start_step_index : last_report_step
    last_report_step_mutated = final_last_report_step != last_report_step

    test_output_header_written = report_block_executed && diagnostic_level >= 2
    current_1_3_written = report_block_executed && diagnostic_level >= 3
    current_4_6_written =
        current_1_3_written && include_extra_power_coil_currents
    flux_iteration_written = report_block_executed && diagnostic_level >= 2
    node_table_written =
        report_block_executed && start_step_index <= 1 && diagnostic_level > 3
    thevenin_scalar_written = report_block_executed && diagnostic_level >= 2
    rotor_matrix_written = report_block_executed && diagnostic_level >= 2
    stator_matrix_written = report_block_executed && diagnostic_level >= 2

    current_report_values = zeros(Float64, 6)
    current_1_3_written && (current_report_values[1:3] .= currents[1:3])
    current_4_6_written && (current_report_values[4:6] .= currents[4:6])
    flux_report_values =
        flux_iteration_written ? [d_flux, q_flux, Float64(iteration_count)] : zeros(Float64, 3)
    first_terminal_report_nodes = node_table_written ? copy(node1) : zeros(Int, length(node1))
    second_terminal_report_nodes = node_table_written ? copy(node2) : zeros(Int, length(node2))
    mechanical_report_node = node_table_written ? mechanical_node : 0
    rotor_matrix_values =
        rotor_matrix_written && !isempty(rotor_matrix) ? vec(copy(rotor_matrix)) : zeros(Float64, 9)
    stator_matrix_values =
        stator_matrix_written && !isempty(stator_matrix) ? vec(copy(stator_matrix)) : zeros(Float64, 9)

    report_context = (
        initial_header_candidate = initial_header_candidate,
        start_step_index = start_step_index,
        machine_index = machine_index,
        diagnostic_level = diagnostic_level,
        time_s = time,
        component_count = component_count,
        report_block_executed = report_block_executed,
        time_header_written = time_header_written,
        test_output_header_written = test_output_header_written,
        current_1_3_written = current_1_3_written,
        current_4_6_written = current_4_6_written,
        current_report_values = current_report_values,
        flux_iteration_written = flux_iteration_written,
        flux_report_values = flux_report_values,
        node_table_written = node_table_written,
        first_terminal_report_nodes = first_terminal_report_nodes,
        second_terminal_report_nodes = second_terminal_report_nodes,
        mechanical_report_node = mechanical_report_node,
        thevenin_scalar_written = thevenin_scalar_written,
        mechanical_thevenin_resistance = thevenin_scalar_written ? zthevm : 0.0,
        rotor_matrix_written = rotor_matrix_written,
        rotor_thevenin_values = rotor_matrix_values,
        stator_matrix_written = stator_matrix_written,
        stator_thevenin_values = stator_matrix_values,
    )
    report_text_lines = emit_text ? _universal_machine_report_text_lines(report_context) : String[]

    return (
        source = :universal_machine_report_schedule,
        fortran_labels = (
            16000, 16005, 16008, 16012, 16014, 16016, 16020, 16052, 16054,
            16055, 16056, 16057, 16060, 16061, 16062, 16077, 16078, 16080,
            16082,
        ),
        start_step_index = start_step_index,
        machine_index = machine_index,
        machine_count = machine_count,
        last_report_step_before = last_report_step,
        output_interval_steps = output_interval_steps,
        next_due_step = next_due_step,
        final_last_report_step = final_last_report_step,
        last_report_step_mutated = last_report_step_mutated,
        diagnostic_level = diagnostic_level,
        report_block_executed = report_block_executed,
        report_skipped_to_postsolve = report_skipped_to_postsolve,
        separator_write_count = separator_write_count,
        time_header_written = time_header_written,
        time_header_values = time_header_written ? [time, Float64(component_count)] : zeros(Float64, 2),
        test_output_header_written = test_output_header_written,
        current_1_3_written = current_1_3_written,
        current_4_6_written = current_4_6_written,
        current_report_values = current_report_values,
        flux_iteration_written = flux_iteration_written,
        flux_report_values = flux_report_values,
        node_table_written = node_table_written,
        first_terminal_report_nodes = first_terminal_report_nodes,
        second_terminal_report_nodes = second_terminal_report_nodes,
        mechanical_report_node = mechanical_report_node,
        thevenin_scalar_written = thevenin_scalar_written,
        mechanical_thevenin_resistance = thevenin_scalar_written ? zthevm : 0.0,
        rotor_matrix_written = rotor_matrix_written,
        rotor_thevenin_values = rotor_matrix_values,
        stator_matrix_written = stator_matrix_written,
        stator_thevenin_values = stator_matrix_values,
        report_text_lines = report_text_lines,
        report_text_line_count = length(report_text_lines),
        report_schedule_mutated = false,
        report_text_mutated = false,
        exact_text_formatting = emit_text,
        complete_machine_equation_solution = false,
        full_machine_equation_solution = false,
        deferred_effects = (
            :complete_machine_equation_solution,
            :full_machine_equation_solution,
        ),
    )
end

function universal_machine_report_schedule_update!(
    state::UniversalMachineReportScheduleState;
    kwargs...,
)
    preview = universal_machine_report_schedule_preview(
        ;
        last_report_step = state.last_report_step,
        kwargs...,
    )
    state.last_report_step = preview.final_last_report_step
    empty!(state.report_text_lines)
    append!(state.report_text_lines, preview.report_text_lines)
    state.report_schedule_mutated = preview.last_report_step_mutated
    state.report_text_mutated = preview.exact_text_formatting
    return merge(
        preview,
        (
            final_last_report_step = state.last_report_step,
            report_text_lines = copy(state.report_text_lines),
            report_text_line_count = length(state.report_text_lines),
            report_schedule_mutated = state.report_schedule_mutated,
            report_text_mutated = state.report_text_mutated,
        ),
    )
end

function _machine_uniform_rows_per_call(call_indices::Vector{Int})
    isempty(call_indices) && return 0
    counts = Dict{Int,Int}()
    for call_index in call_indices
        call_index > 0 || throw(ArgumentError("call indices must be positive"))
        counts[call_index] = get(counts, call_index, 0) + 1
    end
    rows_per_call = first(values(counts))
    all(==(rows_per_call), values(counts)) ||
        throw(ArgumentError("terminal-current rows must have a uniform row count per call"))
    return rows_per_call
end

function _machine_terminal_output_slot_count(output_slots::Vector{Int})
    slots = sort!(unique(slot for slot in output_slots if slot > 0))
    return length(slots)
end

function _machine_time_bounds(times_s::Vector{Float64})
    isempty(times_s) && return (0.0, 0.0)
    for index in 2:length(times_s)
        times_s[index] + eps(Float64) >= times_s[index - 1] ||
            throw(ArgumentError("terminal-current times must be nondecreasing"))
    end
    return (first(times_s), last(times_s))
end

function _machine_symbol_vector(name::String, values, expected_count::Int)
    length(values) == expected_count ||
        throw(ArgumentError("$name must contain $expected_count entries"))
    return Symbol[Symbol(value) for value in values]
end

function _machine_bool_vector(name::String, values, expected_count::Int)
    length(values) == expected_count ||
        throw(ArgumentError("$name must contain $expected_count entries"))
    return Bool[Bool(value) for value in values]
end

function synchronous_machine_equation_step_preview(
    cu_values::AbstractVector{<:Real},
    histq_values::AbstractVector{<:Real},
    shp_values::AbstractVector{<:Real};
    phase_voltages::AbstractVector{<:Real},
    current_history::AbstractVector{<:Real},
    electrical_coefficients::AbstractVector{<:Real},
    numask::Int,
    nlocg::Int,
    nloce::Int = 0,
    delta2::Real,
    angle_half_step_inverse::Real,
    speed_tolerance::Real,
    omega_tolerance::Real = Inf,
    speed_floor::Real = 1.0e-12,
    max_iterations::Int = 1,
    machine_sequence_index::Int = 0,
    timestep_index::Int = 0,
    active_execution_chain::Int = 0,
    mechanical_torque_multipliers::Union{Nothing,AbstractVector{<:Real}}=nothing,
    total_applied_torque::Union{Nothing,Real}=nothing,
    applied_torque_distribution::Union{Nothing,AbstractVector{<:Real}}=nothing,
    tenm6::Real = 1.0e-6,
    asqrt3::Real = inv(sqrt(3.0)),
    sqrt3::Real = sqrt(3.0),
    sqrt32::Real = sqrt(3.0) / 2.0,
)
    length(phase_voltages) == 3 || throw(ArgumentError("phase_voltages must contain three phase voltages"))
    length(current_history) == 3 || throw(ArgumentError("current_history must contain three phase-current history values"))
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    Float64(delta2) > 0.0 || throw(ArgumentError("delta2 must be positive"))
    Float64(angle_half_step_inverse) != 0.0 ||
        throw(ArgumentError("angle_half_step_inverse must be nonzero"))
    Float64(speed_tolerance) >= 0.0 || throw(ArgumentError("speed_tolerance must be nonnegative"))
    Float64(omega_tolerance) >= 0.0 || throw(ArgumentError("omega_tolerance must be nonnegative"))
    Float64(speed_floor) > 0.0 || throw(ArgumentError("speed_floor must be positive"))
    nlocg > 0 || throw(ArgumentError("nlocg must identify the rotor angle history index"))
    machine_sequence_index >= 0 || throw(ArgumentError("machine_sequence_index must be nonnegative"))
    timestep_index >= 0 || throw(ArgumentError("timestep_index must be nonnegative"))
    active_execution_chain >= 0 || throw(ArgumentError("active_execution_chain must be nonnegative"))

    cu = _machine_float_vector(cu_values)
    histq = _machine_float_vector(histq_values)
    shp = _machine_float_vector(shp_values)
    elp = _machine_float_vector(electrical_coefficients)
    indices = _synchronous_machine_equation_check_storage(cu, histq, shp, elp, numask, nlocg, nloce)
    torque_multipliers = mechanical_torque_multipliers === nothing ?
        ones(Float64, numask) : Float64.(mechanical_torque_multipliers)
    length(torque_multipliers) == numask || throw(ArgumentError(
        "mechanical_torque_multipliers must contain one value per rotor mass",
    ))
    all(isfinite, torque_multipliers) || throw(ArgumentError(
        "mechanical_torque_multipliers must be finite",
    ))
    total_torque = total_applied_torque === nothing ? nothing :
        Float64(total_applied_torque)
    total_torque === nothing || isfinite(total_torque) || throw(ArgumentError(
        "total_applied_torque must be finite",
    ))
    torque_distribution_scale = 1.0
    if total_torque !== nothing
        applied_torque_range = (indices.n27 + 1):(indices.n27 + numask)
        base_torque = sum(@view histq[applied_torque_range])
        if abs(base_torque) > eps(Float64)
            torque_distribution_scale = total_torque / base_torque
        else
            applied_torque_distribution === nothing && throw(ArgumentError(
                "zero-base total applied-torque control requires an explicit rotor-mass distribution",
            ))
            distribution_weights = Float64.(applied_torque_distribution)
            length(distribution_weights) == numask || throw(ArgumentError(
                "applied_torque_distribution must contain one value per rotor mass",
            ))
            all(isfinite, distribution_weights) || throw(ArgumentError(
                "applied_torque_distribution entries must be finite",
            ))
            weight_sum = sum(distribution_weights)
            abs(weight_sum) > eps(Float64) || throw(ArgumentError(
                "total applied-torque control requires declared rotor-mass torque fractions",
            ))
            histq[applied_torque_range] .=
                total_torque .* distribution_weights ./ weight_sum
        end
    end

    v1 = Float64(phase_voltages[1])
    v2 = Float64(phase_voltages[2])
    v3 = Float64(phase_voltages[3])
    a1 = elp[27]
    a2 = elp[28]
    d6 = -a1 * v1 - a2 * (v2 + v3) - Float64(current_history[1])
    d7 = -a1 * v2 - a2 * (v1 + v3) - Float64(current_history[2])
    d8 = -a1 * v3 - a2 * (v1 + v2) - Float64(current_history[3])
    cv3 = (v1 + v2 + v3) * Float64(asqrt3)
    a5 = -(cv3 - cu[3]) * elp[17]
    cz = elp[26]
    czt = cz * Float64(tenm6)
    acde = elp[29]
    acdf = elp[30]
    q3 = -cu[11] * cu[12] * czt
    etot = cu[20]
    sum_value = cu[21]
    dsped = histq[indices.ksg]
    spdd = Inf
    iteration_count = 0
    converged = false

    a3 = a4 = cv1 = cv2 = c1 = c2 = c3 = c4 = ac1 = ac2 = cd = cexc = spdn = dang = d2 = 0.0
    for iteration in 1:max_iterations
        iteration_count = iteration
        etot *= Float64(angle_half_step_inverse)
        sum_value *= Float64(angle_half_step_inverse)
        au = etot * 0.5
        av = sum_value * Float64(sqrt32)
        tsd = av - au
        tsc = -au - av
        a3 = d6 * etot + d7 * tsd + d8 * tsc
        cv1 = v1 * etot + v2 * tsd + v3 * tsc
        au *= Float64(sqrt3)
        av = sum_value * 0.5
        tsd = au - av
        tsc = -av - au
        a4 = d6 * sum_value + d7 * tsc + d8 * tsd
        cv2 = v1 * sum_value + v2 * tsc + v3 * tsd
        c1 = elp[32] * a3 + elp[36] * cu[4] + elp[37] * cu[5]
        c2 = elp[33] * a3 + elp[37] * cu[4] + elp[39] * cu[5]
        c3 = elp[34] * a4 + elp[40] * cu[6] + elp[41] * cu[7]
        c4 = elp[35] * a4 + elp[41] * cu[6] + elp[43] * cu[7]
        ac2 = elp[1] * a3 + elp[2] * c1 + elp[4] * c2
        ac1 = elp[9] * a4 + elp[10] * c3 + elp[12] * c4
        cd = (ac2 * acde * a4 - ac1 * acdf * a3) * czt
        cexc = nloce == 0 ? 0.0 : histq[indices.ksex]

        jt = indices.n27
        for ik in indices.ikw:indices.ikp
            jt += 1
            mass_index = ik - indices.ikw + 1
            histq[ik] =
                histq[jt] * torque_distribution_scale *
                torque_multipliers[mass_index] / histq[ik]
        end
        histq[indices.ksg] -= cd
        if nloce != 0
            cexc = (q3 * c1) / cexc
            histq[indices.ksex] -= cexc
        end

        jt = indices.n26
        for ik in indices.ikw:indices.ikp
            jt += 1
            histq[jt] = -histq[ik]
            histq[ik] = histq[ik] - histq[ik + indices.num2]
        end
        _machine_bansol_segment!(shp, 0, histq, indices.ikw - 1, numask)
        spdn = histq[indices.ksg]
        dsped_denominator = abs(dsped) <= Float64(speed_floor) ? Float64(speed_floor) : dsped
        spdd = abs((spdn - dsped) / dsped_denominator)
        if spdd <= Float64(speed_tolerance)
            converged = true
            break
        end
        dsped = spdn
        dang = (histq[nlocg + indices.num2] + Float64(delta2) * spdn) * cz
        etot = cos(dang)
        sum_value = sin(dang)
    end

    if !converged && spdd > Float64(omega_tolerance)
        nonconvergence_report_iteration_index = iteration_count + 1
        return (
            source = :synchronous_machine_equation_step,
            cu_values = cu,
            histq_values = histq,
            shp_values = shp,
            phase_voltages = Float64.(phase_voltages),
            current_history = Float64.(current_history),
            phase_a_current_residual = d6,
            phase_b_current_residual = d7,
            phase_c_current_residual = d8,
            iteration_count = iteration_count,
            converged = false,
            failed_nonconvergence = true,
            nonconvergence_report_iteration_index = nonconvergence_report_iteration_index,
            machine_iteration_warning = false,
            machine_nonconvergence_stop = true,
            max_relative_speed_error = spdd,
            stop_status_code = 206,
            stop_previous_execution_chain = active_execution_chain,
            stop_machine_sequence_index = machine_sequence_index,
            stop_timestep_index = timestep_index,
            stop_relative_speed_error = spdd,
            kill_code = 104,
            machine_equation_step_executed = true,
            mechanical_torque_multipliers = copy(torque_multipliers),
            total_applied_torque = total_torque,
            mechanical_torque_control_executed =
                mechanical_torque_multipliers !== nothing ||
                total_torque !== nothing,
            complete_machine_equation_solution = false,
            full_solvum_execution = false,
            deferred_effects = (
                :nonconverged_machine_stop_report,
                :multi_machine_solvum_orchestration,
                :tacs_transfer,
                :output_report_mutation,
            ),
        )
    end

    kc = indices.n22 + 5 * numask
    for ka in 1:numask
        jt = ka
        ik = indices.ikv + ka
        kb = indices.n26 + ka
        kd = kc + ka
        histq[ik + indices.num2] = histq[kb] - shp[kd]
        histq[jt] = histq[jt + indices.num2] + Float64(delta2) * histq[ik]
    end
    cu[1] = a3
    cu[2] = a4
    cu[3] = a5
    cu[4] = c1
    cu[5] = c2
    cu[6] = c3
    cu[7] = c4
    d2 = spdn * cz
    q3 = cu[11] * cu[12]
    dang = histq[nlocg] * cz

    return (
        source = :synchronous_machine_equation_step,
        cu_values = cu,
        histq_values = histq,
        shp_values = shp,
        phase_voltages = Float64.(phase_voltages),
        current_history = Float64.(current_history),
        phase_a_current_residual = d6,
        phase_b_current_residual = d7,
        phase_c_current_residual = d8,
        d_axis_current = a3,
        q_axis_current = a4,
        zero_sequence_current = a5,
        d_axis_voltage = cv1,
        q_axis_voltage = cv2,
        zero_sequence_voltage = cv3,
        d_axis_rotor_current_1 = c1,
        d_axis_rotor_current_2 = c2,
        q_axis_rotor_current_1 = c3,
        q_axis_rotor_current_2 = c4,
        d_axis_internal_voltage = ac2,
        q_axis_internal_voltage = ac1,
        electromagnetic_torque = cd,
        exciter_torque = cexc,
        mechanical_speed = spdn,
        mechanical_speed_error = spdd,
        mechanical_angle = dang,
        scaled_mechanical_speed = d2,
        mechanical_power = q3,
        iteration_count = iteration_count,
        converged = converged,
        failed_nonconvergence = false,
        nonconvergence_report_iteration_index = converged ? iteration_count : iteration_count + 1,
        machine_iteration_warning = !converged,
        machine_nonconvergence_stop = false,
        stop_status_code = 0,
        stop_previous_execution_chain = 0,
        stop_machine_sequence_index = 0,
        stop_timestep_index = 0,
        stop_relative_speed_error = 0.0,
        kill_code = 0,
        machine_equation_step_executed = true,
        mechanical_torque_multipliers = copy(torque_multipliers),
        total_applied_torque = total_torque,
        mechanical_torque_control_executed =
            mechanical_torque_multipliers !== nothing || total_torque !== nothing,
        machine_equation_state_mutated = false,
        complete_machine_equation_solution = false,
        full_solvum_execution = false,
        deferred_effects = (
            :multi_machine_solvum_orchestration,
            :tacs_transfer,
            :output_report_mutation,
            :nonconverged_machine_stop_report,
            :full_machine_equation_solution,
        ),
    )
end

function synchronous_machine_equation_step_preview(
    state::SynchronousMachineEquationState;
    kwargs...,
)
    return synchronous_machine_equation_step_preview(
        state.cu_values,
        state.histq_values,
        state.shp_values;
        kwargs...,
    )
end

function synchronous_machine_equation_step!(
    state::SynchronousMachineEquationState;
    kwargs...,
)
    preview = synchronous_machine_equation_step_preview(state; kwargs...)
    empty!(state.cu_values)
    empty!(state.histq_values)
    empty!(state.shp_values)
    append!(state.cu_values, preview.cu_values)
    append!(state.histq_values, preview.histq_values)
    append!(state.shp_values, preview.shp_values)
    state.iteration_count = preview.iteration_count
    state.converged = preview.converged
    state.equation_mutated = preview.machine_equation_step_executed
    return merge(
        preview,
        (
            machine_equation_state_mutated = state.equation_mutated,
            machine_equation_step_executed = state.equation_mutated,
        ),
    )
end

function synchronous_machine_companion_update!(
    state::SynchronousMachineEquationState,
    current_history::Vector{Float64},
    electrical_coefficients::Vector{Float64},
    equation_result;
    numask::Int,
    nlocg::Int,
    delta2::Real,
    damping_ratio::Real,
    rotor_angle_extrapolation_interval::Real,
    speed_voltage_factor::Real,
    electrical_angle_increment::Real,
    field_voltage_multiplier::Union{Nothing,Real}=nothing,
    external_field_voltage_input_pu::Union{Nothing,Real}=nothing,
    previous_external_field_voltage_input_pu::Union{Nothing,Real}=nothing,
    saturation_enabled::Bool=false,
    d_axis_saturation_region::Int=0,
    q_axis_saturation_region::Int=0,
)
    length(current_history) == 3 ||
        throw(ArgumentError("current_history must contain three phase values"))
    equation_result.failed_nonconvergence &&
        throw(ArgumentError("cannot update companion history after a failed machine solve"))
    cu = state.cu_values
    histq = state.histq_values
    shp = state.shp_values
    elp = electrical_coefficients
    indices = _synchronous_machine_equation_check_storage(
        cu,
        histq,
        shp,
        elp,
        numask,
        nlocg,
        0,
    )
    delta2_value = Float64(delta2)
    damping = Float64(damping_ratio)
    extrapolation_interval = Float64(rotor_angle_extrapolation_interval)
    speed_factor = Float64(speed_voltage_factor)
    angle_increment = Float64(electrical_angle_increment)
    (external_field_voltage_input_pu === nothing) ==
    (previous_external_field_voltage_input_pu === nothing) ||
        throw(ArgumentError(
            "current and previous external field-voltage inputs must be supplied together",
        ))
    field_voltage_multiplier !== nothing &&
    external_field_voltage_input_pu !== nothing &&
        throw(ArgumentError(
            "field-voltage multiplier and external field-voltage input are mutually exclusive",
        ))
    external_field_voltage_input_pu === nothing ||
        isfinite(Float64(external_field_voltage_input_pu)) ||
        throw(ArgumentError("external_field_voltage_input_pu must be finite"))

    accumulated_coupling = 0.0
    angle_index = 1
    speed_index = indices.ikv
    predicted_index = speed_index + indices.num2
    matrix_index = indices.n22 + 6 * numask
    previous_angle = histq[angle_index]
    previous_speed = histq[speed_index + 1]
    if numask != 1
        for mass_index in angle_index:(indices.ikv - 1)
            speed_index += 1
            predicted_index += 1
            matrix_index += 4
            next_angle = histq[mass_index + 1]
            next_speed = histq[speed_index + 1]
            diagonal_angle = shp[matrix_index - 2]
            diagonal_speed = shp[matrix_index]
            histq[predicted_index] +=
                shp[matrix_index - 3] * previous_angle +
                shp[matrix_index - 1] * previous_speed +
                diagonal_angle * next_angle +
                diagonal_speed * next_speed +
                accumulated_coupling
            accumulated_coupling =
                diagonal_angle * previous_angle + diagonal_speed * previous_speed
            previous_angle = next_angle
            previous_speed = next_speed
        end
        accumulated_coupling += shp[matrix_index + 1] * previous_angle
    end
    histq[predicted_index + 1] +=
        shp[matrix_index + 3] * previous_speed + accumulated_coupling

    matrix_index = indices.n22
    speed_index = indices.ikv
    predicted_index = 1 + indices.num2
    for angle_index in 1:indices.ikv
        speed_index += 1
        histq[predicted_index] =
            histq[angle_index] + delta2_value * histq[speed_index]
        matrix_index += 1
        previous_matrix_speed = shp[matrix_index]
        shp[matrix_index] = histq[speed_index]
        histq[speed_index] = 2.0 * histq[speed_index] - previous_matrix_speed
        predicted_index += 1
    end

    mechanical_angle = equation_result.mechanical_angle
    scaled_speed = equation_result.scaled_mechanical_speed
    alpha =
        9.0 * (cu[22] - mechanical_angle) + cu[23] +
        extrapolation_interval * (scaled_speed + cu[24])
    cu[23] = cu[22]
    cu[22] = mechanical_angle
    cu[24] = scaled_speed

    d_axis_current = equation_result.d_axis_current
    q_axis_current = equation_result.q_axis_current
    d_axis_rotor_current_1 = equation_result.d_axis_rotor_current_1
    d_axis_rotor_current_2 = equation_result.d_axis_rotor_current_2
    q_axis_rotor_current_1 = equation_result.q_axis_rotor_current_1
    q_axis_rotor_current_2 = equation_result.q_axis_rotor_current_2
    saturation = if saturation_enabled
        _synchronous_machine_saturation_update!(
            elp;
            d_axis_current,
            q_axis_current,
            d_axis_rotor_current_1,
            d_axis_rotor_current_2,
            q_axis_rotor_current_1,
            q_axis_rotor_current_2,
            d_axis_region = d_axis_saturation_region,
            q_axis_region = q_axis_saturation_region,
            damping_ratio = damping,
        )
    else
        (
            flux_magnitude = 0.0,
            d_axis_region = d_axis_saturation_region,
            q_axis_region = q_axis_saturation_region,
            d_axis_factor = 1.0,
            q_axis_factor = 1.0,
            refactorized = false,
        )
    end
    d_axis_flux_voltage =
        equation_result.d_axis_internal_voltage * elp[29] +
        elp[19] * d_axis_current
    q_axis_flux_voltage =
        -equation_result.q_axis_internal_voltage * elp[30] -
        elp[19] * q_axis_current
    predicted_q_axis_current =
        (2.5 * q_axis_current - 1.5 * cu[14] + cu[9]) * 0.5
    cu[14] = cu[9]
    cu[9] = q_axis_current

    field_voltage = cu[11]
    if external_field_voltage_input_pu !== nothing
        current_external_voltage = Float64(external_field_voltage_input_pu)
        previous_external_voltage =
            Float64(previous_external_field_voltage_input_pu)
        all(isfinite, (current_external_voltage, previous_external_voltage)) ||
            throw(ArgumentError("external field-voltage inputs must be finite"))
        # The alternative UMDATA port is oriented into the field winding. The
        # trapezoidal companion consumes the interval-average applied voltage
        # through its present and damped-history terms. Its initial previous
        # value is zero, so a step input contributes half its magnitude on the
        # first history update and its full magnitude thereafter.
        field_voltage = -0.5 * (1.0 + damping) * (
            previous_external_voltage + current_external_voltage
        )
    elseif field_voltage_multiplier !== nothing
        cu[12] = Float64(field_voltage_multiplier)
        field_voltage *= cu[12]
    end
    d_axis_history =
        elp[44] * d_axis_current + elp[45] * d_axis_rotor_current_1 +
        elp[46] * d_axis_rotor_current_2 +
        (q_axis_flux_voltage * scaled_speed - equation_result.d_axis_voltage) * damping
    q_axis_history =
        elp[47] * q_axis_current + elp[48] * q_axis_rotor_current_1 +
        elp[49] * q_axis_rotor_current_2 +
        (d_axis_flux_voltage * scaled_speed - equation_result.q_axis_voltage) * damping
    zero_sequence_history =
        elp[18] * cu[3] - equation_result.zero_sequence_voltage * damping
    rotor_history = Float64[
        elp[50] * d_axis_current + elp[51] * d_axis_rotor_current_1 +
        elp[52] * d_axis_rotor_current_2,
        elp[53] * d_axis_current + elp[54] * d_axis_rotor_current_1 +
        elp[55] * d_axis_rotor_current_2,
        elp[56] * q_axis_current + elp[57] * q_axis_rotor_current_1 +
        elp[58] * q_axis_rotor_current_2,
        elp[59] * q_axis_current + elp[60] * q_axis_rotor_current_1 +
        elp[61] * q_axis_rotor_current_2,
    ]
    rotor_history[1] -= equation_result.mechanical_power * damping + field_voltage
    d_axis_rotor_correction = -(elp[70] * rotor_history[1] + elp[71] * rotor_history[2])
    q_axis_rotor_correction = -(elp[72] * rotor_history[3] + elp[73] * rotor_history[4])
    cu[1] = d_axis_history + d_axis_rotor_correction
    cu[2] = q_axis_history + q_axis_rotor_correction
    cu[3] = zero_sequence_history
    cu[4:7] .= rotor_history

    q_axis_rotor_correction -= elp[79] * predicted_q_axis_current
    rotor_offset = cu[19] - mechanical_angle
    rotor_cosine = cos(rotor_offset)
    rotor_sine = sin(rotor_offset)
    synchronous_d_axis_voltage =
        rotor_cosine * d_axis_flux_voltage + rotor_sine * q_axis_flux_voltage
    synchronous_q_axis_voltage =
        rotor_sine * d_axis_flux_voltage - rotor_cosine * q_axis_flux_voltage
    predicted_d_axis_speed_voltage =
        (2.5 * synchronous_d_axis_voltage - 1.5 * cu[15] + cu[16]) * speed_factor
    predicted_q_axis_speed_voltage =
        (2.5 * synchronous_q_axis_voltage - 1.5 * cu[17] + cu[18]) * speed_factor
    cu[15] = cu[16]
    cu[16] = synchronous_d_axis_voltage
    cu[17] = cu[18]
    cu[18] = synchronous_q_axis_voltage

    alpha_cosine = cos(alpha)
    alpha_sine = sin(alpha)
    cu[19] += angle_increment
    cu[20] = alpha_cosine
    cu[21] = alpha_sine
    reference_cosine = cos(cu[19])
    reference_sine = sin(cu[19])
    transform_d = reference_cosine * alpha_cosine + reference_sine * alpha_sine
    transform_q = reference_sine * alpha_cosine - reference_cosine * alpha_sine
    corrected_d_axis_history =
        transform_d * d_axis_rotor_correction -
        transform_q * q_axis_rotor_correction
    corrected_q_axis_history =
        transform_q * d_axis_rotor_correction +
        transform_d * q_axis_rotor_correction
    phase_d_axis_history =
        (rotor_cosine * d_axis_history - rotor_sine * q_axis_history -
         predicted_q_axis_speed_voltage + corrected_d_axis_history) * elp[80]
    phase_q_axis_history =
        (rotor_sine * d_axis_history + rotor_cosine * q_axis_history +
         predicted_d_axis_speed_voltage + corrected_q_axis_history) * elp[80]
    phase_zero_sequence_history = zero_sequence_history * elp[81]
    current_history[1] =
        -phase_d_axis_history * reference_cosine -
        phase_q_axis_history * reference_sine -
        phase_zero_sequence_history
    phase_b_d_axis = -0.5 * reference_cosine + (sqrt(3.0) / 2.0) * reference_sine
    phase_c_d_axis = -0.5 * reference_cosine - (sqrt(3.0) / 2.0) * reference_sine
    phase_b_q_axis = -0.5 * reference_sine - (sqrt(3.0) / 2.0) * reference_cosine
    phase_c_q_axis = -0.5 * reference_sine + (sqrt(3.0) / 2.0) * reference_cosine
    current_history[2] =
        -phase_b_d_axis * phase_d_axis_history -
        phase_b_q_axis * phase_q_axis_history -
        phase_zero_sequence_history
    current_history[3] =
        -phase_c_d_axis * phase_d_axis_history -
        phase_c_q_axis * phase_q_axis_history -
        phase_zero_sequence_history

    return (
        source = :synchronous_machine_companion_update,
        current_history = copy(current_history),
        companion_history_mutated = true,
        saturation_flux_magnitude = saturation.flux_magnitude,
        d_axis_saturation_region = saturation.d_axis_region,
        q_axis_saturation_region = saturation.q_axis_region,
        d_axis_saturation_factor = saturation.d_axis_factor,
        q_axis_saturation_factor = saturation.q_axis_factor,
        saturation_matrix_refactorized = saturation.refactorized,
    )
end

function synchronous_machine_dynamic_step!(
    state::SynchronousMachineDynamicState;
    phase_voltages::AbstractVector{<:Real},
    numask::Int,
    nlocg::Int,
    delta2::Real,
    angle_half_step_inverse::Real,
    speed_tolerance::Real,
    damping_ratio::Real,
    rotor_angle_extrapolation_interval::Real,
    speed_voltage_factor::Real,
    electrical_angle_increment::Real,
    field_voltage_multiplier::Union{Nothing,Real}=nothing,
    external_field_voltage_input_pu::Union{Nothing,Real}=nothing,
    mechanical_torque_multipliers::Union{Nothing,AbstractVector{<:Real}}=nothing,
    total_applied_torque::Union{Nothing,Real}=nothing,
    kwargs...,
)
    field_voltage_multiplier !== nothing &&
    external_field_voltage_input_pu !== nothing &&
        throw(ArgumentError(
            "field-voltage multiplier and external field-voltage input are mutually exclusive",
        ))
    external_field_voltage_input_pu === nothing ||
        isfinite(Float64(external_field_voltage_input_pu)) ||
        throw(ArgumentError("external_field_voltage_input_pu must be finite"))
    equation_result = synchronous_machine_equation_step!(
        state.equation_state;
        phase_voltages = phase_voltages,
        current_history = state.current_history,
        electrical_coefficients = state.electrical_coefficients,
        numask = numask,
        nlocg = nlocg,
        delta2 = delta2,
        angle_half_step_inverse = angle_half_step_inverse,
        speed_tolerance = speed_tolerance,
        mechanical_torque_multipliers = mechanical_torque_multipliers,
        total_applied_torque = total_applied_torque,
        applied_torque_distribution = state.applied_torque_distribution,
        kwargs...,
    )
    previous_external_field_voltage_input_pu =
        external_field_voltage_input_pu === nothing ? nothing :
        state.exciter_port_state.voltage_input_pu
    companion_result = synchronous_machine_companion_update!(
        state.equation_state,
        state.current_history,
        state.electrical_coefficients,
        equation_result;
        numask = numask,
        nlocg = nlocg,
        delta2 = delta2,
        damping_ratio = damping_ratio,
        rotor_angle_extrapolation_interval = rotor_angle_extrapolation_interval,
        speed_voltage_factor = speed_voltage_factor,
        electrical_angle_increment = electrical_angle_increment,
        saturation_enabled = state.saturation_enabled,
        d_axis_saturation_region = state.d_axis_saturation_region,
        q_axis_saturation_region = state.q_axis_saturation_region,
        field_voltage_multiplier = field_voltage_multiplier,
        external_field_voltage_input_pu = external_field_voltage_input_pu,
        previous_external_field_voltage_input_pu =
            previous_external_field_voltage_input_pu,
    )
    state.call_count += 1
    state.companion_history_mutated = companion_result.companion_history_mutated
    state.d_axis_saturation_region = companion_result.d_axis_saturation_region
    state.q_axis_saturation_region = companion_result.q_axis_saturation_region
    companion_result.saturation_matrix_refactorized &&
        (state.saturation_refactor_count += 1)
    referred_field_current_a = equation_result.d_axis_rotor_current_1
    field_winding_current_a =
        referred_field_current_a /
        state.exciter_port_state.current_reduction_factor
    exciter_port_result = synchronous_machine_exciter_port_update!(
        state.exciter_port_state,
        field_winding_current_a,
        external_field_voltage_input_pu = external_field_voltage_input_pu,
    )
    return merge(
        equation_result,
        (
            companion_current_history = companion_result.current_history,
            companion_history_mutated = state.companion_history_mutated,
            saturation_enabled = state.saturation_enabled,
            saturation_flux_magnitude = companion_result.saturation_flux_magnitude,
            d_axis_saturation_region = state.d_axis_saturation_region,
            q_axis_saturation_region = state.q_axis_saturation_region,
            d_axis_saturation_factor = companion_result.d_axis_saturation_factor,
            q_axis_saturation_factor = companion_result.q_axis_saturation_factor,
            saturation_matrix_refactorized =
                companion_result.saturation_matrix_refactorized,
            saturation_refactor_count = state.saturation_refactor_count,
            dynamic_call_count = state.call_count,
            exciter_port_field_winding_current_a =
                exciter_port_result.field_winding_current_a,
            exciter_port_current_reduction_factor =
                exciter_port_result.current_reduction_factor,
            exciter_port_current_a = exciter_port_result.terminal_current_a,
            exciter_port_sensor_closed = exciter_port_result.sensor_closed,
            exciter_port_update_count = exciter_port_result.update_count,
            exciter_port_current_mutated = exciter_port_result.current_mutated,
            field_voltage_multiplier =
                field_voltage_multiplier === nothing ? 1.0 :
                Float64(field_voltage_multiplier),
            external_field_voltage_input_pu =
                external_field_voltage_input_pu === nothing ? 0.0 :
                Float64(external_field_voltage_input_pu),
            previous_external_field_voltage_input_pu =
                previous_external_field_voltage_input_pu === nothing ? 0.0 :
                Float64(previous_external_field_voltage_input_pu),
            field_voltage_control_executed =
                field_voltage_multiplier !== nothing ||
                external_field_voltage_input_pu !== nothing,
            external_field_voltage_control_executed =
                external_field_voltage_input_pu !== nothing,
            exciter_port_voltage_input_pu =
                exciter_port_result.voltage_input_pu,
            exciter_port_previous_voltage_input_pu =
                exciter_port_result.previous_voltage_input_pu,
            exciter_port_voltage_update_count =
                exciter_port_result.voltage_update_count,
            mechanical_torque_multipliers =
                mechanical_torque_multipliers === nothing ?
                ones(Float64, numask) : Float64.(mechanical_torque_multipliers),
            total_applied_torque = total_applied_torque === nothing ? nothing :
                Float64(total_applied_torque),
            mechanical_torque_control_executed =
                mechanical_torque_multipliers !== nothing ||
                total_applied_torque !== nothing,
            complete_machine_equation_solution = true,
            deferred_effects = (:tacs_transfer, :output_report_mutation),
        ),
    )
end

function _machine_source_stage_counts(source_stage_names::Vector{Symbol})
    return (
        created = count(==(:created), source_stage_names),
        changed = count(==(:changed), source_stage_names),
        slack = count(==(:slack), source_stage_names),
    )
end

function machine_network_coupling_state_preview(
    header_call_indices::AbstractVector{<:Integer};
    machine_counts::AbstractVector{<:Integer},
    coupling_row_counts::AbstractVector{<:Integer},
    output_counts::AbstractVector{<:Integer},
    source_constant_starts::AbstractVector{<:Integer},
    branch_counts::AbstractVector{<:Integer},
    vector_call_indices::AbstractVector{<:Integer},
    vector_names,
    vector_indices::AbstractVector{<:Integer},
    vector_values::AbstractVector{<:Real},
    source_call_indices::AbstractVector{<:Integer},
    source_stage_names,
    source_kinds,
    source_indices::AbstractVector{<:Integer},
    source_has_nodes::AbstractVector{Bool},
    source_nodes::AbstractVector{<:Integer},
    source_constant_indices::AbstractVector{<:Integer},
    source_frequencies_hz::AbstractVector{<:Real},
    source_has_crests::AbstractVector{Bool},
    source_crests::AbstractVector{<:Real},
    source_has_time1::AbstractVector{Bool},
    source_time1_s::AbstractVector{<:Real},
    source_has_start_times::AbstractVector{Bool},
    source_start_times_s::AbstractVector{<:Real},
    source_has_stop_times::AbstractVector{Bool},
    source_stop_times_s::AbstractVector{<:Real},
    exit_call_indices::AbstractVector{<:Integer},
    exit_start_indices::AbstractVector{<:Integer},
    exit_source_constant_starts::AbstractVector{<:Integer},
    exit_branch_counts::AbstractVector{<:Integer},
    exit_reference_frequencies_rad_s::AbstractVector{<:Real},
    exit_twopi_values::AbstractVector{<:Real},
)
    header_count = length(header_call_indices)
    vector_count = length(vector_call_indices)
    source_count = length(source_call_indices)
    exit_count = length(exit_call_indices)

    headers = _machine_int_vector("header_call_indices", header_call_indices, header_count)
    machines = _machine_int_vector("machine_counts", machine_counts, header_count)
    coupling_rows = _machine_int_vector("coupling_row_counts", coupling_row_counts, header_count)
    outputs = _machine_int_vector("output_counts", output_counts, header_count)
    source_constants =
        _machine_int_vector("source_constant_starts", source_constant_starts, header_count)
    branches = _machine_int_vector("branch_counts", branch_counts, header_count)

    vector_calls = _machine_int_vector("vector_call_indices", vector_call_indices, vector_count)
    vector_name_values = _machine_symbol_vector("vector_names", vector_names, vector_count)
    vector_index_values = _machine_int_vector("vector_indices", vector_indices, vector_count)
    vector_value_values = _machine_real_vector("vector_values", vector_values, vector_count)

    source_calls = _machine_int_vector("source_call_indices", source_call_indices, source_count)
    source_stage_values = _machine_symbol_vector("source_stage_names", source_stage_names, source_count)
    source_kind_values = _machine_symbol_vector("source_kinds", source_kinds, source_count)
    source_index_values = _machine_int_vector("source_indices", source_indices, source_count)
    source_node_flags = _machine_bool_vector("source_has_nodes", source_has_nodes, source_count)
    source_node_values = _machine_int_vector("source_nodes", source_nodes, source_count)
    source_constant_values =
        _machine_int_vector("source_constant_indices", source_constant_indices, source_count)
    source_frequency_values =
        _machine_real_vector("source_frequencies_hz", source_frequencies_hz, source_count)
    source_crest_flags = _machine_bool_vector("source_has_crests", source_has_crests, source_count)
    source_crest_values = _machine_real_vector("source_crests", source_crests, source_count)
    source_time1_flags = _machine_bool_vector("source_has_time1", source_has_time1, source_count)
    source_time1_values = _machine_real_vector("source_time1_s", source_time1_s, source_count)
    source_start_flags =
        _machine_bool_vector("source_has_start_times", source_has_start_times, source_count)
    source_start_values =
        _machine_real_vector("source_start_times_s", source_start_times_s, source_count)
    source_stop_flags =
        _machine_bool_vector("source_has_stop_times", source_has_stop_times, source_count)
    source_stop_values = _machine_real_vector("source_stop_times_s", source_stop_times_s, source_count)

    exit_calls = _machine_int_vector("exit_call_indices", exit_call_indices, exit_count)
    exit_starts = _machine_int_vector("exit_start_indices", exit_start_indices, exit_count)
    exit_source_constants =
        _machine_int_vector("exit_source_constant_starts", exit_source_constant_starts, exit_count)
    exit_branches = _machine_int_vector("exit_branch_counts", exit_branch_counts, exit_count)
    exit_reference_frequencies =
        _machine_real_vector("exit_reference_frequencies_rad_s", exit_reference_frequencies_rad_s, exit_count)
    exit_twopi = _machine_real_vector("exit_twopi_values", exit_twopi_values, exit_count)

    header_count > 0 || throw(ArgumentError("header_call_indices must be nonempty"))
    vector_count > 0 || throw(ArgumentError("vector_call_indices must be nonempty"))
    source_count > 0 || throw(ArgumentError("source_call_indices must be nonempty"))
    exit_count > 0 || throw(ArgumentError("exit_call_indices must be nonempty"))
    all(>(0), headers) || throw(ArgumentError("header call indices must be positive"))
    all(>(0), vector_calls) || throw(ArgumentError("vector call indices must be positive"))
    all(>(0), source_calls) || throw(ArgumentError("source call indices must be positive"))
    all(>(0), exit_calls) || throw(ArgumentError("exit call indices must be positive"))
    all(>(0), machines) || throw(ArgumentError("machine counts must be positive"))
    all(>(0), coupling_rows) || throw(ArgumentError("coupling row counts must be positive"))
    all(>(0), outputs) || throw(ArgumentError("output counts must be positive"))
    all(>=(0), branches) || throw(ArgumentError("branch counts must be nonnegative"))
    all(pair -> !pair[1] || pair[2] != 0, zip(source_node_flags, source_node_values)) ||
        throw(ArgumentError("flagged source nodes must be nonzero"))
    source_stage_count = _machine_source_stage_counts(source_stage_values)

    return (
        source = :machine_network_coupling_state,
        header_call_indices = headers,
        machine_counts = machines,
        coupling_row_counts = coupling_rows,
        output_counts = outputs,
        source_constant_starts = source_constants,
        branch_counts = branches,
        vector_call_indices = vector_calls,
        vector_names = vector_name_values,
        vector_indices = vector_index_values,
        vector_values = vector_value_values,
        source_call_indices = source_calls,
        source_stage_names = source_stage_values,
        source_kinds = source_kind_values,
        source_indices = source_index_values,
        source_has_nodes = source_node_flags,
        source_nodes = source_node_values,
        source_constant_indices = source_constant_values,
        source_frequencies_hz = source_frequency_values,
        source_has_crests = source_crest_flags,
        source_crests = source_crest_values,
        source_has_time1 = source_time1_flags,
        source_time1_s = source_time1_values,
        source_has_start_times = source_start_flags,
        source_start_times_s = source_start_values,
        source_has_stop_times = source_stop_flags,
        source_stop_times_s = source_stop_values,
        exit_call_indices = exit_calls,
        exit_start_indices = exit_starts,
        exit_source_constant_starts = exit_source_constants,
        exit_branch_counts = exit_branches,
        exit_reference_frequencies_rad_s = exit_reference_frequencies,
        exit_twopi_values = exit_twopi,
        header_count = header_count,
        vector_row_count = vector_count,
        created_source_count = source_stage_count.created,
        changed_source_count = source_stage_count.changed,
        slack_source_count = source_stage_count.slack,
        exit_count = exit_count,
        max_machine_count = maximum(machines),
        max_output_count = maximum(outputs),
        max_coupling_row_count = maximum(coupling_rows),
        initial_source_constant_start = first(source_constants),
        final_source_constant_start = last(exit_source_constants),
        max_branch_count = maximum(branches),
        machine_network_rows_mutated = false,
        machine_network_coupling_replayed = false,
        machine_equation_solution_executed = false,
        saturation_equivalence_executed = false,
        deferred_effects = (
            :machine_equation_solution,
            :machine_saturation_equivalence,
            :tacs_machine_interface,
        ),
    )
end

function machine_network_coupling_state_update!(
    state::MachineNetworkCouplingState,
    header_call_indices::AbstractVector{<:Integer};
    kwargs...,
)
    preview = machine_network_coupling_state_preview(header_call_indices; kwargs...)
    for field in (
        :header_call_indices,
        :machine_counts,
        :coupling_row_counts,
        :output_counts,
        :source_constant_starts,
        :branch_counts,
        :vector_call_indices,
        :vector_names,
        :vector_indices,
        :vector_values,
        :source_call_indices,
        :source_stage_names,
        :source_kinds,
        :source_indices,
        :source_has_nodes,
        :source_nodes,
        :source_constant_indices,
        :source_frequencies_hz,
        :source_has_crests,
        :source_crests,
        :source_has_time1,
        :source_time1_s,
        :source_has_start_times,
        :source_start_times_s,
        :source_has_stop_times,
        :source_stop_times_s,
        :exit_call_indices,
        :exit_start_indices,
        :exit_source_constant_starts,
        :exit_branch_counts,
        :exit_reference_frequencies_rad_s,
        :exit_twopi_values,
    )
        empty!(getfield(state, field))
        append!(getfield(state, field), getproperty(preview, field))
    end
    state.header_count = preview.header_count
    state.vector_row_count = preview.vector_row_count
    state.created_source_count = preview.created_source_count
    state.changed_source_count = preview.changed_source_count
    state.slack_source_count = preview.slack_source_count
    state.exit_count = preview.exit_count
    state.machine_network_rows_mutated = preview.header_count > 0
    return merge(
        preview,
        (
            machine_network_rows_mutated = state.machine_network_rows_mutated,
            machine_network_coupling_replayed = state.machine_network_rows_mutated,
        ),
    )
end

function machine_terminal_current_state_preview(
    call_indices::AbstractVector{<:Integer};
    row_indices::AbstractVector{<:Integer},
    times_s::AbstractVector{<:Real},
    phase_nodes::AbstractVector{<:Integer},
    terminal_nodes::AbstractVector{<:Integer},
    output_slots::AbstractVector{<:Integer},
    coupling_flags::AbstractVector{<:Integer},
    terminal_currents::AbstractVector{<:Real},
    output_values::AbstractVector{<:Real},
)
    row_count = length(call_indices)
    calls = _machine_int_vector("call_indices", call_indices, row_count)
    rows = _machine_int_vector("row_indices", row_indices, row_count)
    times = _machine_real_vector("times_s", times_s, row_count)
    phases = _machine_int_vector("phase_nodes", phase_nodes, row_count)
    terminals = _machine_int_vector("terminal_nodes", terminal_nodes, row_count)
    slots = _machine_int_vector("output_slots", output_slots, row_count)
    flags = _machine_int_vector("coupling_flags", coupling_flags, row_count)
    currents = _machine_real_vector("terminal_currents", terminal_currents, row_count)
    outputs = _machine_real_vector("output_values", output_values, row_count)

    all(>(0), rows) || throw(ArgumentError("row indices must be positive"))
    all(>=(0), phases) || throw(ArgumentError("phase nodes must be nonnegative"))
    all(>=(0), terminals) || throw(ArgumentError("terminal nodes must be nonnegative"))
    all(>=(0), slots) || throw(ArgumentError("output slots must be nonnegative"))
    all(>=(0), flags) || throw(ArgumentError("coupling flags must be nonnegative"))
    rows_per_call = _machine_uniform_rows_per_call(calls)
    first_time, last_time = _machine_time_bounds(times)

    return (
        source = :machine_terminal_current_state,
        call_indices = calls,
        row_indices = rows,
        times_s = times,
        phase_nodes = phases,
        terminal_nodes = terminals,
        output_slots = slots,
        coupling_flags = flags,
        terminal_currents = currents,
        output_values = outputs,
        row_count = row_count,
        call_count = length(unique(calls)),
        rows_per_call = rows_per_call,
        terminal_output_slot_count = _machine_terminal_output_slot_count(slots),
        first_time_s = first_time,
        last_time_s = last_time,
        terminal_rows_mutated = false,
        terminal_current_state_replayed = false,
        machine_equation_solution_executed = false,
        network_coupling_executed = false,
        deferred_effects = (
            :machine_equation_solution,
            :network_coupling,
            :saturation_equivalence,
        ),
    )
end

function machine_terminal_current_state_update!(
    state::MachineTerminalCurrentState,
    call_indices::AbstractVector{<:Integer};
    kwargs...,
)
    preview = machine_terminal_current_state_preview(call_indices; kwargs...)
    empty!(state.call_indices)
    empty!(state.row_indices)
    empty!(state.times_s)
    empty!(state.phase_nodes)
    empty!(state.terminal_nodes)
    empty!(state.output_slots)
    empty!(state.coupling_flags)
    empty!(state.terminal_currents)
    empty!(state.output_values)
    append!(state.call_indices, preview.call_indices)
    append!(state.row_indices, preview.row_indices)
    append!(state.times_s, preview.times_s)
    append!(state.phase_nodes, preview.phase_nodes)
    append!(state.terminal_nodes, preview.terminal_nodes)
    append!(state.output_slots, preview.output_slots)
    append!(state.coupling_flags, preview.coupling_flags)
    append!(state.terminal_currents, preview.terminal_currents)
    append!(state.output_values, preview.output_values)
    state.row_count = preview.row_count
    state.call_count = preview.call_count
    state.rows_per_call = preview.rows_per_call
    state.terminal_output_slot_count = preview.terminal_output_slot_count
    state.first_time_s = preview.first_time_s
    state.last_time_s = preview.last_time_s
    state.terminal_rows_mutated = preview.row_count > 0
    return merge(
        preview,
        (
            terminal_rows_mutated = state.terminal_rows_mutated,
            terminal_current_state_replayed = state.terminal_rows_mutated,
        ),
    )
end

function _synchronous_machine_tacs_field_components(
    a3::Real,
    a4::Real,
    c1::Real,
    c2::Real,
    c3::Real,
    c4::Real,
    elp_i26_21::Real,
    elp_i75_4::Real,
)
    sf4 = Float64(a3) * Float64(elp_i26_21) + Float64(c1) + Float64(c2)
    sf5 = (Float64(a4) * Float64(elp_i26_21) + Float64(c3) + Float64(c4)) * Float64(elp_i75_4)
    return sf4, sf5
end

function synchronous_machine_tacs_transfer_preview(
    ismtac_requests::AbstractVector{Int};
    cu_values::AbstractVector{<:Real}=Float64[],
    cv1::Real=0.0,
    cv2::Real=0.0,
    cv3::Real=0.0,
    a3::Real=0.0,
    a4::Real=0.0,
    c1::Real=0.0,
    c2::Real=0.0,
    c3::Real=0.0,
    c4::Real=0.0,
    q3::Real=0.0,
    cd::Real=0.0,
    cexc::Real=0.0,
    ac1::Real=0.0,
    ac2::Real=0.0,
    acde::Real=0.0,
    acdf::Real=0.0,
    elp_i26_19::Real=0.0,
    elp_i26_21::Real=0.0,
    elp_i75_4::Real=1.0,
    histq_values::AbstractVector{<:Real}=Float64[],
    shp_values::AbstractVector{<:Real}=Float64[],
    iu::Int=0,
    numask::Int=0,
    n22::Int=0,
    lmset::Int=0,
    cu_offset::Int=0,
)
    iu >= 0 || throw(ArgumentError("iu must be nonnegative"))
    numask >= 0 || throw(ArgumentError("numask must be nonnegative"))
    n22 >= 0 || throw(ArgumentError("n22 must be nonnegative"))
    lmset >= 0 || throw(ArgumentError("lmset must be nonnegative"))
    cu_offset >= 0 || throw(ArgumentError("cu_offset must be nonnegative"))

    cu = _machine_float_vector(cu_values)
    histq = _machine_float_vector(histq_values)
    shp = _machine_float_vector(shp_values)
    request_values = Float64[]
    request_kinds = Symbol[]
    source_labels = Int[]
    source_indices = Int[]

    for request in ismtac_requests
        request != 0 || throw(ArgumentError("ISMTAC request entries must be nonzero"))
        if request < 0
            code = -request
            value =
                if 1 <= code <= 7
                    cu_index = cu_offset + code
                    cu_index <= length(cu) ||
                        throw(ArgumentError("cu_values must contain every requested negative ISMTAC CU code"))
                    push!(request_kinds, :cu)
                    push!(source_labels, 1680)
                    cu[cu_index]
                elseif code == 8
                    push!(request_kinds, :cv1)
                    push!(source_labels, 1680)
                    Float64(cv1)
                elseif code == 9
                    push!(request_kinds, :cv2)
                    push!(source_labels, 1681)
                    Float64(cv2)
                elseif code == 10
                    push!(request_kinds, :cv3)
                    push!(source_labels, 1682)
                    Float64(cv3)
                elseif code == 11
                    push!(request_kinds, :mechanical_power)
                    push!(source_labels, 1683)
                    Float64(q3)
                elseif code == 12 || code == 13
                    sf4, sf5 = _synchronous_machine_tacs_field_components(
                        a3,
                        a4,
                        c1,
                        c2,
                        c3,
                        c4,
                        elp_i26_21,
                        elp_i75_4,
                    )
                    if code == 12
                        push!(request_kinds, :field_current_magnitude)
                        push!(source_labels, 1684)
                        sqrt(sf4 * sf4 + sf5 * sf5)
                    else
                        push!(request_kinds, :field_current_angle)
                        push!(source_labels, 1685)
                        atan(sf5, sf4)
                    end
                elseif code == 14
                    push!(request_kinds, :electromagnetic_torque)
                    push!(source_labels, 1686)
                    Float64(cd)
                elseif code == 15
                    push!(request_kinds, :exciter_torque)
                    push!(source_labels, 1687)
                    Float64(cexc)
                elseif code == 16
                    push!(request_kinds, :internal_voltage_d)
                    push!(source_labels, 1688)
                    Float64(ac2) * Float64(acde) + Float64(elp_i26_19) * Float64(a3)
                elseif code == 17
                    push!(request_kinds, :internal_voltage_q)
                    push!(source_labels, 1689)
                    Float64(ac1) * Float64(acdf) + Float64(elp_i26_19) * Float64(a4)
                else
                    throw(ArgumentError("negative ISMTAC request codes must be in the range -17:-1"))
                end
            push!(request_values, value)
            push!(source_indices, code)
        else
            numask > 0 || throw(ArgumentError("numask must be positive for positive ISMTAC machine-history requests"))
            num2 = 2 * numask
            if request <= num2
                histq_index = iu + request
                _machine_check_index(histq_index, histq, "histq_values")
                push!(request_values, histq[histq_index])
                push!(request_kinds, :mechanical_history)
                push!(source_labels, 1696)
                push!(source_indices, request)
            else
                relative = request - num2
                relative > 0 || throw(ArgumentError("positive ISMTAC shaft request must exceed 2*numask"))
                angle_index = iu + relative
                speed_index = angle_index + numask
                shp_index = n22 + request
                _machine_check_index(angle_index, histq, "histq_values")
                _machine_check_index(angle_index + 1, histq, "histq_values")
                _machine_check_index(speed_index, histq, "histq_values")
                _machine_check_index(speed_index + 1, histq, "histq_values")
                _machine_check_index(shp_index, shp, "shp_values")
                _machine_check_index(shp_index + numask, shp, "shp_values")
                push!(
                    request_values,
                    shp[shp_index + numask] * (histq[angle_index] - histq[angle_index + 1]) +
                    shp[shp_index] * (histq[speed_index] - histq[speed_index + 1]),
                )
                push!(request_kinds, :shaft_torque)
                push!(source_labels, 1697)
                push!(source_indices, request)
            end
        end
    end

    transfer_count = length(request_values)
    storage_indices =
        transfer_count == 0 ? Int[] : collect((lmset + 1):(lmset + transfer_count))
    return (
        source = :over16_synchronous_machine_tacs_transfer,
        fortran_labels = (1680, 1681, 1682, 1683, 1684, 1685, 1686, 1687, 1688, 1689, 1696, 1697),
        ismtac_requests = collect(ismtac_requests),
        transfer_count = transfer_count,
        lmset = lmset,
        final_lmset = lmset + transfer_count,
        etac_storage_indices = storage_indices,
        request_values = request_values,
        request_kinds = request_kinds,
        source_labels = source_labels,
        source_indices = source_indices,
        cu_offset = cu_offset,
        etac_mutated = false,
        machine_tacs_transfer_mutated = false,
        machine_solve_executed = false,
        tacs_executed = false,
        solvum_executed = false,
        deferred_effects = (
            :full_machine_solve,
            :complete_machine_matrix_integration,
            :remaining_saturable_machine_paths,
            :tacs_execution,
            :network_timestep_coupling,
        ),
    )
end

function synchronous_machine_tacs_transfer_update!(
    state::SynchronousMachineTACSInterfaceState,
    ismtac_requests::AbstractVector{Int};
    kwargs...,
)
    preview = synchronous_machine_tacs_transfer_preview(
        ismtac_requests;
        kwargs...,
        lmset = state.lmset,
    )
    if preview.transfer_count > 0
        last_index = state.lmset + preview.transfer_count
        length(state.etac_values) >= last_index ||
            throw(ArgumentError("etac_values must contain every ETAC output storage index"))
        for i in 1:preview.transfer_count
            state.etac_values[state.lmset + i] = preview.request_values[i]
        end
    end
    state.lmset = preview.final_lmset
    state.transfer_count = preview.transfer_count
    state.transfer_mutated = preview.transfer_count > 0
    return merge(
        preview,
        (
            etac_values = copy(state.etac_values),
            final_lmset = state.lmset,
            etac_mutated = state.transfer_mutated,
            machine_tacs_transfer_mutated = state.transfer_mutated,
        ),
    )
end

function _past_machine_tacs_kwargs(
    past_state::PASTMachineHistoryState,
    kwargs::NamedTuple,
)
    haskey(kwargs, :histq_values) &&
        throw(ArgumentError("histq_values must come from the PAST machine-history state"))
    haskey(kwargs, :shp_values) &&
        throw(ArgumentError("shp_values must come from the PAST machine-history state"))
    return merge(kwargs, (histq_values = past_state.histq, shp_values = past_state.shp))
end

function synchronous_machine_tacs_transfer_from_past_preview(
    past_state::PASTMachineHistoryState,
    ismtac_requests::AbstractVector{Int};
    kwargs...,
)
    preview = synchronous_machine_tacs_transfer_preview(
        ismtac_requests;
        _past_machine_tacs_kwargs(past_state, (; kwargs...))...,
    )
    return merge(
        preview,
        (
            source = :over16_synchronous_machine_tacs_transfer_from_past,
            past_machine_history_consumed = true,
            past_state_updated = past_state.past_state_updated,
            past_state_mutated = false,
            machine_solve_executed = false,
        ),
    )
end

function synchronous_machine_tacs_transfer_from_past_update!(
    state::SynchronousMachineTACSInterfaceState,
    past_state::PASTMachineHistoryState,
    ismtac_requests::AbstractVector{Int};
    kwargs...,
)
    preview = synchronous_machine_tacs_transfer_from_past_preview(
        past_state,
        ismtac_requests;
        kwargs...,
        lmset = state.lmset,
    )
    if preview.transfer_count > 0
        last_index = state.lmset + preview.transfer_count
        length(state.etac_values) >= last_index ||
            throw(ArgumentError("etac_values must contain every ETAC output storage index"))
        for i in 1:preview.transfer_count
            state.etac_values[state.lmset + i] = preview.request_values[i]
        end
    end
    state.lmset = preview.final_lmset
    state.transfer_count = preview.transfer_count
    state.transfer_mutated = preview.transfer_count > 0
    return merge(
        preview,
        (
            etac_values = copy(state.etac_values),
            final_lmset = state.lmset,
            etac_mutated = state.transfer_mutated,
            machine_tacs_transfer_mutated = state.transfer_mutated,
        ),
    )
end

const _SYNCHRONOUS_MACHINE_TACS_REQUESTED_PASS = 6644
const _SYNCHRONOUS_MACHINE_TACS_COMPLETED_PASS = 7766

function _machine_table_int(value::Real, name::String)
    number = Float64(value)
    isfinite(number) || throw(ArgumentError("$name must be finite"))
    isinteger(number) || throw(ArgumentError("$name must be an integer-valued table entry"))
    return Int(number)
end

function _machine_output_table_bounds(table::AbstractVector{Float64}, output_count::Int)
    output_count >= 0 || throw(ArgumentError("output_count must be nonnegative"))
    length(table) >= output_count + 3 ||
        throw(ArgumentError("machine_output_table must contain the TACS interface header"))
    interface_requested = table[output_count + 1] == -9999.0
    body_start = output_count + 4
    last_entry = _machine_table_int(table[output_count + 3], "last machine-output TACS table entry")
    last_entry <= length(table) ||
        throw(ArgumentError("machine_output_table is shorter than the advertised TACS table"))
    return (; interface_requested, body_start, last_transfer_row = last_entry - 2)
end

function _synchronous_machine_initialize_rotor_mass_history!(
    table::Vector{Float64},
    bounds;
    node_values::Vector{Float64},
    switch_positions::Vector{Float64},
    mechanical_angles_rad::Vector{Float64},
    pole_pair_counts::Vector{Float64},
    generator_mass_nodes::Vector{Int},
)
    initialized_rows = Int[]
    source_labels = Int[]
    machine_index = 0
    row = bounds.body_start
    last_scan_row = bounds.last_transfer_row - 2
    while row <= last_scan_row
        code = _machine_table_int(table[row], "machine-output TACS code")
        if code <= -300
            row + 4 <= length(table) ||
                throw(ArgumentError("rotor-mass TACS block header is incomplete"))
            machine_index += 1
            _machine_check_index(machine_index, mechanical_angles_rad, "mechanical_angles_rad")
            _machine_check_index(machine_index, pole_pair_counts, "pole_pair_counts")
            _machine_check_index(machine_index, generator_mass_nodes, "generator_mass_nodes")
            entry_count = _machine_table_int(table[row + 3], "rotor-mass TACS entry count")
            mass_count = _machine_table_int(table[row + 4], "rotor-mass count")
            mass_count > 0 || throw(ArgumentError("rotor-mass count must be positive"))
            entry_count >= 5 * mass_count ||
                throw(ArgumentError("rotor-mass TACS entry count is shorter than the mass rows"))
            row + entry_count - 1 <= length(table) ||
                throw(ArgumentError("rotor-mass TACS block exceeds machine_output_table storage"))
            pole_pairs = pole_pair_counts[machine_index]
            pole_pairs != 0.0 || throw(ArgumentError("pole_pair_counts entries must be nonzero"))
            generator_angle = mechanical_angles_rad[machine_index] + 2.0 * pi / (4.0 * pole_pairs)
            first_node = _machine_table_int(table[row + 1], "first rotor-mass node index")
            _machine_check_index(first_node, node_values, "node_values")
            table[row + 3] = mass_count == 1 ? generator_angle : 0.0
            table[row + 4] = node_values[first_node]
            push!(initialized_rows, row)
            push!(source_labels, 17010)

            generator_mass_index = mass_count == 1 ? 1 : 0
            generator_shift = 0.0
            switch_start = -code - 300
            for mass_index in 2:mass_count
                mass_row = row + (mass_index - 1) * 5
                previous_row = row + (mass_index - 2) * 5
                node = _machine_table_int(table[mass_row + 1], "rotor-mass node index")
                previous_node = _machine_table_int(table[previous_row + 1], "previous rotor-mass node index")
                switch_index = switch_start + mass_index - 1
                _machine_check_index(node, node_values, "node_values")
                _machine_check_index(previous_node, node_values, "node_values")
                _machine_check_index(switch_index, switch_positions, "switch_positions")
                table[mass_row + 3] =
                    table[previous_row + 3] +
                    table[mass_row + 3] * switch_positions[switch_index] +
                    table[mass_row + 4] * (node_values[previous_node] - node_values[node])
                if node == generator_mass_nodes[machine_index]
                    generator_mass_index = mass_index
                    generator_shift = table[mass_row + 3] - generator_angle
                end
                table[mass_row + 4] = node_values[node]
                push!(initialized_rows, mass_row)
                push!(source_labels, 17015)
            end
            generator_mass_index != 0 ||
                throw(ArgumentError("generator_mass_nodes must identify one rotor mass in each multi-mass block"))
            if mass_count > 1
                for mass_index in 1:mass_count
                    mass_row = row + (mass_index - 1) * 5
                    table[mass_row + 3] -= generator_shift
                    if mass_index == generator_mass_index
                        table[mass_row + 3] = generator_angle
                    end
                    push!(source_labels, 17030)
                end
            end
            row += entry_count
        else
            row += 1
        end
    end
    return (; initialized_rows, initialization_source_labels = source_labels)
end

function _synchronous_machine_transfer_rotor_mass_tacs!(
    table::Vector{Float64},
    tacs::Vector{Float64},
    bounds;
    node_values::Vector{Float64},
    switch_positions::Vector{Float64},
    time_step_s::Float64,
)
    transfer_codes = Int[]
    transfer_kinds = Symbol[]
    tacs_indices = Int[]
    tacs_written_values = Float64[]
    table_history_indices = Int[]
    table_history_values = Float64[]
    source_labels = Int[]
    row = bounds.body_start
    half_step = time_step_s / 2.0
    while row <= bounds.last_transfer_row
        row + 2 <= length(table) ||
            throw(ArgumentError("machine-output TACS transfer row is incomplete"))
        code = _machine_table_int(table[row], "machine-output TACS code")
        node_or_switch = _machine_table_int(table[row + 1], "machine-output TACS source index")
        tacs_index = _machine_table_int(table[row + 2], "machine-output TACS destination index")
        if code == -1
            _machine_check_index(node_or_switch, switch_positions, "switch_positions")
            _machine_check_index(tacs_index, tacs, "tacs_values")
            tacs[tacs_index] = -switch_positions[node_or_switch]
            push!(transfer_kinds, :shaft_torque)
            push!(source_labels, 17070)
            row += 3
        elseif code == -2
            _machine_check_index(node_or_switch, node_values, "node_values")
            _machine_check_index(tacs_index, tacs, "tacs_values")
            tacs[tacs_index] = node_values[node_or_switch]
            push!(transfer_kinds, :mass_speed)
            push!(source_labels, 17100)
            row += 3
        elseif code <= -300 || code == -3
            row + 4 <= length(table) ||
                throw(ArgumentError("rotor-mass angle transfer row is incomplete"))
            if tacs_index != 0
                _machine_check_index(node_or_switch, node_values, "node_values")
                _machine_check_index(tacs_index, tacs, "tacs_values")
                tacs[tacs_index] = table[row + 3] + table[row + 4] * node_values[node_or_switch] * half_step
                table[row + 3] = tacs[tacs_index]
                table[row + 4] = half_step * node_values[node_or_switch]
                push!(tacs_indices, tacs_index)
                push!(tacs_written_values, tacs[tacs_index])
                push!(table_history_indices, row + 3)
                push!(table_history_values, table[row + 3])
                push!(table_history_indices, row + 4)
                push!(table_history_values, table[row + 4])
            end
            push!(transfer_kinds, :mass_angle)
            push!(source_labels, 17114)
            push!(transfer_codes, code)
            row += 5
            continue
        elseif code == -299
            row + 4 <= length(table) ||
                throw(ArgumentError("rotor-mass skip row is incomplete"))
            push!(transfer_kinds, :mass_angle_skip)
            push!(source_labels, 17116)
            push!(transfer_codes, code)
            row += 5
            continue
        elseif code == -4
            _machine_check_index(node_or_switch, switch_positions, "switch_positions")
            _machine_check_index(tacs_index, tacs, "tacs_values")
            tacs[tacs_index] = -switch_positions[node_or_switch]
            push!(transfer_kinds, :exciter_current)
            push!(source_labels, 17120)
            row += 3
        elseif code == -5
            _machine_check_index(node_or_switch, node_values, "node_values")
            _machine_check_index(tacs_index, tacs, "tacs_values")
            tacs[tacs_index] = node_values[node_or_switch]
            push!(transfer_kinds, :exciter_voltage)
            push!(source_labels, 17130)
            row += 3
        else
            throw(ArgumentError("unsupported machine-output TACS transfer code $code"))
        end
        push!(transfer_codes, code)
        push!(tacs_indices, tacs_index)
        push!(tacs_written_values, tacs[tacs_index])
    end
    return (;
        transfer_codes,
        transfer_kinds,
        tacs_indices,
        tacs_written_values,
        table_history_indices,
        table_history_values,
        transfer_source_labels = source_labels,
    )
end

function synchronous_machine_rotor_mass_tacs_preview(
    machine_output_table::AbstractVector{<:Real};
    output_count::Int=0,
    tacs_values::AbstractVector{<:Real}=Float64[],
    node_values::AbstractVector{<:Real}=Float64[],
    switch_positions::AbstractVector{<:Real}=Float64[],
    mechanical_angles_rad::AbstractVector{<:Real}=Float64[],
    pole_pair_counts::AbstractVector{<:Real}=Float64[],
    generator_mass_nodes::AbstractVector{<:Integer}=Int[],
    time_step_s::Real=0.0,
    initial_step::Bool=false,
    transfer_pass_marker::Int=0,
)
    table = _machine_float_vector(machine_output_table)
    tacs = _machine_float_vector(tacs_values)
    nodes = _machine_float_vector(node_values)
    switches = _machine_float_vector(switch_positions)
    angles = _machine_float_vector(mechanical_angles_rad)
    pole_pairs = _machine_float_vector(pole_pair_counts)
    generator_nodes = Int[Int(node) for node in generator_mass_nodes]
    time_step = Float64(time_step_s)
    time_step >= 0.0 || throw(ArgumentError("time_step_s must be nonnegative"))
    bounds = _machine_output_table_bounds(table, output_count)
    final_marker = transfer_pass_marker
    mass_result = (initialized_rows = Int[], initialization_source_labels = Int[])
    transfer_result = (
        transfer_codes = Int[],
        transfer_kinds = Symbol[],
        tacs_indices = Int[],
        tacs_written_values = Float64[],
        table_history_indices = Int[],
        table_history_values = Float64[],
        transfer_source_labels = Int[],
    )
    if bounds.interface_requested
        if initial_step && table[output_count + 2] != 0.0
            mass_result = _synchronous_machine_initialize_rotor_mass_history!(
                table,
                bounds;
                node_values = nodes,
                switch_positions = switches,
                mechanical_angles_rad = angles,
                pole_pair_counts = pole_pairs,
                generator_mass_nodes = generator_nodes,
            )
        end
        transfer_allowed = initial_step || transfer_pass_marker == _SYNCHRONOUS_MACHINE_TACS_REQUESTED_PASS
        if transfer_allowed
            transfer_result = _synchronous_machine_transfer_rotor_mass_tacs!(
                table,
                tacs,
                bounds;
                node_values = nodes,
                switch_positions = switches,
                time_step_s = time_step,
            )
            if !initial_step && transfer_pass_marker == _SYNCHRONOUS_MACHINE_TACS_REQUESTED_PASS
                final_marker = _SYNCHRONOUS_MACHINE_TACS_COMPLETED_PASS
            end
        end
    end
    return (
        source = :over16_synchronous_machine_rotor_mass_tacs_transfer,
        fortran_labels = (
            17000, 17010, 17015, 17020, 17030, 17040, 17042, 17050, 17060,
            17070, 17100, 17110, 17112, 17114, 17116, 17120, 17130, 17140,
        ),
        machine_output_table = table,
        tacs_values = tacs,
        initialized_mass_rows = mass_result.initialized_rows,
        initialization_source_labels = mass_result.initialization_source_labels,
        transfer_codes = transfer_result.transfer_codes,
        transfer_kinds = transfer_result.transfer_kinds,
        tacs_indices = transfer_result.tacs_indices,
        tacs_written_values = transfer_result.tacs_written_values,
        table_history_indices = transfer_result.table_history_indices,
        table_history_values = transfer_result.table_history_values,
        transfer_source_labels = transfer_result.transfer_source_labels,
        transfer_pass_marker = transfer_pass_marker,
        final_transfer_pass_marker = final_marker,
        mass_history_initialized = false,
        tacs_transfer_mutated = false,
        machine_output_table_mutated = false,
        complete_tacs_execution = false,
        complete_output_report_mutation = false,
        full_machine_equation_solution = false,
        deferred_effects = (
            :complete_tacs_execution,
            :complete_output_report_mutation,
            :full_machine_equation_solution,
        ),
    )
end

function synchronous_machine_rotor_mass_tacs_update!(
    state::SynchronousMachineRotorMassTACSState;
    kwargs...,
)
    preview = synchronous_machine_rotor_mass_tacs_preview(
        state.machine_output_table;
        tacs_values = state.tacs_values,
        transfer_pass_marker = state.transfer_pass_marker,
        kwargs...,
    )
    empty!(state.machine_output_table)
    append!(state.machine_output_table, preview.machine_output_table)
    empty!(state.tacs_values)
    append!(state.tacs_values, preview.tacs_values)
    state.transfer_pass_marker = preview.final_transfer_pass_marker
    state.mass_history_initialized = !isempty(preview.initialized_mass_rows)
    state.transfer_mutated = !isempty(preview.transfer_codes)
    return merge(
        preview,
        (
            machine_output_table = copy(state.machine_output_table),
            tacs_values = copy(state.tacs_values),
            final_transfer_pass_marker = state.transfer_pass_marker,
            mass_history_initialized = state.mass_history_initialized,
            tacs_transfer_mutated = state.transfer_mutated,
            machine_output_table_mutated = state.mass_history_initialized || state.transfer_mutated,
        ),
    )
end

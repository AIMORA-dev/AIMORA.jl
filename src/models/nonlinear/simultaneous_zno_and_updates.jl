function _float_vector_values_differ(
    first::AbstractVector{<:Real},
    second::AbstractVector{<:Real},
)
    length(first) == length(second) || return true
    for index in eachindex(first, second)
        Float64(first[index]) == Float64(second[index]) || return true
    end
    return false
end

struct NonlinearCurrentStepResult{
    R<:AbstractVector{Float64},
    N<:AbstractVector{Float64},
    S<:AbstractVector{Float64},
    H<:AbstractVector{Float64},
}
    rhs::R
    nonlinear_sparse_ykm::N
    saturated_transformer_sparse_ykm::S
    hysteretic_inductor_admittance_deltas::H
    saturated_transformer_sparse_update_count::Int
    saturated_transformer_sparse_retriangularization_request_count::Int
    current_update_count::Int
    nonlinear_current_compensation_applied::Bool
    nonlinear_current_state_mutated::Bool
    rhs_mutated::Bool
    anonl_mutated::Bool
    curr_mutated::Bool
    cursub_mutated::Bool
    vzero_mutated::Bool
    vnonl_mutated::Bool
    ilast_mutated::Bool
    cchar_mutated::Bool
    vchar_mutated::Bool
    gslope_mutated::Bool
    current_update_count_mutated::Bool
end

function over16_nonlinear_current_compensation_update!(
    state::OVER16NonlinearInverseColumnState,
    voltages::AbstractVector{<:Real},
    rhs::AbstractVector{<:Real},
    nonlinear_types::AbstractVector{Int},
    nonlinear_from_nodes::AbstractVector{Int},
    nonlinear_to_nodes::AbstractVector{Int},
    nonlinear_admittance_nodes::AbstractVector{Int},
    nonlinear_table_end_indices::AbstractVector{Int},
    nonlinear_subsystem_indices::AbstractVector{Int},
    subsystem_begin_indices::AbstractVector{Int},
    subsystem_owner_rows::AbstractVector{Int},
    subsystem_simultaneous_flags::AbstractVector{Int},
    cchar::AbstractVector{<:Real},
    vchar::AbstractVector{<:Real};
    result_mode::Val{D}=Val(true),
    vnonl::AbstractVector{<:Real}=state.vnonl,
    gslope::AbstractVector{<:Real}=Float64[],
    kwargs...,
) where {D}
    cchar_input = isempty(state.cchar) ? cchar : state.cchar
    vchar_input = isempty(state.vchar) ? vchar : state.vchar
    gslope_input = isempty(state.gslope) ? gslope : state.gslope
    network_response_kwargs = NamedTuple()
    haskey(kwargs, :network_current_response_columns) ||
        (network_response_kwargs = merge(
            network_response_kwargs,
            (network_current_response_columns = state.znonl,),
        ))
    haskey(kwargs, :network_response_node_count) ||
        (network_response_kwargs = merge(
            network_response_kwargs,
            (network_response_node_count = state.ntot,),
        ))
    haskey(kwargs, :network_response_component_count) ||
        (network_response_kwargs = merge(
            network_response_kwargs,
            (network_response_component_count = state.ncomp,),
        ))
    preview = over16_nonlinear_current_compensation_update(
        voltages,
        rhs,
        nonlinear_types,
        nonlinear_from_nodes,
        nonlinear_to_nodes,
        nonlinear_admittance_nodes,
        nonlinear_table_end_indices,
        nonlinear_subsystem_indices,
        subsystem_begin_indices,
        subsystem_owner_rows,
        subsystem_simultaneous_flags,
        cchar_input,
        vchar_input;
        initial_anonl = state.anonl,
        initial_vzero = state.vzero,
        initial_ilast = state.ilast,
        initial_curr = state.curr,
        initial_cursub = state.cursub,
        vnonl = vnonl,
        gslope = gslope_input,
        network_response_kwargs...,
        kwargs...,
    )

    rhs_mutated = _float_vector_values_differ(preview.rhs, rhs)
    anonl_mutated = preview.anonl != state.anonl
    curr_mutated = preview.curr != state.curr
    cursub_mutated = preview.cursub != state.cursub
    vzero_mutated = preview.vzero != state.vzero
    vnonl_mutated = preview.vnonl != state.vnonl
    ilast_mutated = preview.ilast != state.ilast
    cchar_mutated = preview.cchar != state.cchar
    vchar_mutated = preview.vchar != state.vchar
    gslope_mutated = preview.gslope != state.gslope

    resize!(state.anonl, length(preview.anonl))
    state.anonl .= preview.anonl
    resize!(state.curr, length(preview.curr))
    state.curr .= preview.curr
    resize!(state.cursub, length(preview.cursub))
    state.cursub .= preview.cursub
    resize!(state.vzero, length(preview.vzero))
    state.vzero .= preview.vzero
    resize!(state.vnonl, length(preview.vnonl))
    state.vnonl .= preview.vnonl
    resize!(state.ilast, length(preview.ilast))
    state.ilast .= preview.ilast
    resize!(state.cchar, length(preview.cchar))
    state.cchar .= preview.cchar
    resize!(state.vchar, length(preview.vchar))
    state.vchar .= preview.vchar
    resize!(state.gslope, length(preview.gslope))
    state.gslope .= preview.gslope
    state.current_update_count += 1

    # The update counter is advanced on every accepted call, so the aggregate
    # state mutation contract remains true even when all numerical vectors are
    # unchanged.
    state_mutated = true
    step_result = NonlinearCurrentStepResult(
        preview.rhs,
        preview.nonlinear_sparse_ykm,
        preview.saturated_transformer_sparse_ykm,
        preview.hysteretic_inductor_admittance_deltas,
        preview.saturated_transformer_sparse_update_count,
        preview.saturated_transformer_sparse_retriangularization_request_count,
        state.current_update_count,
        preview.nonlinear_current_compensation_applied,
        state_mutated,
        rhs_mutated,
        anonl_mutated,
        curr_mutated,
        cursub_mutated,
        vzero_mutated,
        vnonl_mutated,
        ilast_mutated,
        cchar_mutated,
        vchar_mutated,
        gslope_mutated,
        true,
    )
    D === false && return step_result
    D === true || throw(ArgumentError("result_mode must be Val(true) or Val(false)"))
    return merge(
        preview,
        (
            current_update_count = step_result.current_update_count,
            rhs_mutated = step_result.rhs_mutated,
            anonl_mutated = step_result.anonl_mutated,
            curr_mutated = step_result.curr_mutated,
            cursub_mutated = step_result.cursub_mutated,
            vzero_mutated = step_result.vzero_mutated,
            vnonl_mutated = step_result.vnonl_mutated,
            ilast_mutated = step_result.ilast_mutated,
            cchar_mutated = step_result.cchar_mutated,
            vchar_mutated = step_result.vchar_mutated,
            gslope_mutated = step_result.gslope_mutated,
            current_update_count_mutated = step_result.current_update_count_mutated,
            nonlinear_current_state_mutated = step_result.nonlinear_current_state_mutated,
        ),
    )
end

function _over16_append_block_diagonal(
    existing::AbstractMatrix{<:Real},
    block::AbstractMatrix{<:Real},
)
    existing_rows, existing_cols = size(existing)
    block_rows, block_cols = size(block)
    existing_rows == existing_cols ||
        throw(ArgumentError("existing simultaneous ZNO diagnostic matrix must be square"))
    block_rows == block_cols ||
        throw(ArgumentError("simultaneous ZNO diagnostic block must be square"))
    existing_rows == 0 && return Matrix{Float64}(block)
    combined = zeros(Float64, existing_rows + block_rows, existing_cols + block_cols)
    combined[1:existing_rows, 1:existing_cols] .= existing
    combined[
        existing_rows + 1:existing_rows + block_rows,
        existing_cols + 1:existing_cols + block_cols,
    ] .= block
    return combined
end

function _over16_type94_arrester_update!(
    cchar::Vector{Float64},
    vchar::Vector{Float64},
    a_start::Int,
    b_start::Int,
    source_resistance_term::Float64,
    source_voltage_term::Float64,
    current::Float64,
    delta2::Float64,
    deltat::Float64,
    epsiln::Float64,
)
    a_start >= 1 && a_start + 17 <= length(cchar) ||
        throw(ArgumentError("type-94 cchar state must include A(1:18)"))
    b_start >= 1 && b_start + 10 <= length(vchar) ||
        throw(ArgumentError("type-94 vchar state must include B(1:11)"))
    a(offset) = cchar[a_start + offset - 1]
    b(offset) = vchar[b_start + offset - 1]
    function setb!(offset::Int, value::Float64)
        vchar[b_start + offset - 1] = value
        return value
    end
    a(3) != 0.0 || throw(ArgumentError("type-94 valve block inductance A(3) must not be zero"))
    b(10) != 0.0 || throw(ArgumentError("type-94 voltage division factor B(10) must not be zero"))
    b(11) != 0.0 || throw(ArgumentError("type-94 current division factor B(11) must not be zero"))

    if b(6) <= 0.0
        for offset in 2:9
            setb!(offset, 0.0)
        end
        setb!(6, 1.0)
    end

    sign = 1.0
    setb!(1, b(1) / b(10))
    setb!(3, b(3) / b(11))
    if b(1) <= 0.0
        setb!(1, -b(1))
        b(3) < 0.0 && setb!(3, -b(3))
        sign = -1.0
    end

    ylb = delta2 / a(3)
    ylb != 0.0 || throw(ArgumentError("type-94 valve block admittance must not be zero"))
    cb = b(3) + ylb * b(8)
    rb = b(3) > epsiln ? a(2) * a(1) * b(3) ^ (a(2) - 1.0) : 0.0
    be = (1.0 - a(2)) * a(1) * b(3) ^ a(2)

    region = b(6)
    if region < 3.5
        if b(3) - b(4) < 0.0 && b(3) <= a(18)
            region = 4.0
        elseif region <= 1.5
            if b(2) >= a(7)
                region = b(5) < a(12) ? 2.0 : 3.0
            end
        elseif region <= 2.5 && b(5) >= a(12)
            region = 3.0
        end
    end
    setb!(6, region)

    if region < 1.5
        f0 = b(7) * b(3)
        dfdv = 0.0
        dfdi = b(7)
    elseif region < 2.5
        f1 = a(8) * (a(9) - b(2))
        f2 = a(10) + a(11) * b(3)
        f0 = f1 / f2
        dfdv = -a(8) / f2
        dfdi = -a(11) * f1 / f2^2
    elseif region < 3.5
        wx = a(14) + a(15) * b(5)
        f1 = a(8) * (a(9) - b(2))
        f2 = a(10) + a(11) * b(3)
        f0 = f1 / f2 + wx * (a(13) - b(2))
        dfdv = -a(8) / f2 - wx
        dfdi = a(11) * f1 / f2^2
    else
        f1 = a(8) * (a(16) - b(2))
        f2 = a(10) + a(17) * b(3)^2
        f0 = f1 / f2
        dfdv = -a(8) / f2
        dfdi = -2.0 * a(17) * f1 * b(3) / f2^2
    end
    dfdi == 0.0 && (dfdi = epsiln)
    gi = b(3) + (dfdv * b(2) - f0) / dfdi
    yg = (-dfdv + 2.0 / deltat) / dfdi
    yg != 0.0 || throw(ArgumentError("type-94 gap admittance must not be zero"))
    cg = b(2) * (-dfdv - 2.0 / deltat) / dfdi - b(3) + gi
    art = 1.0 / ylb + rb + 1.0 / yg
    vblock = -cb / ylb
    vgap = -(cg + gi) / yg
    avt = vblock + be + vgap
    adjusted_svt = sign < 0.0 ? -source_voltage_term : source_voltage_term
    denominator = art * b(10) / b(11) - source_resistance_term
    denominator != 0.0 || throw(ArgumentError("type-94 arrester correction denominator must not be zero"))
    carst = (adjusted_svt - avt * b(10)) / denominator
    corrected_current = carst / b(11)
    if carst <= 0.0
        setb!(6, 0.0)
        setb!(3, 0.0)
        setb!(1, 0.0)
        return 0.0
    end

    vgap += corrected_current / yg
    if b(6) <= 1.5
        cip = (corrected_current + b(3)) * 0.5
        corrected_current > a(6) && (cip = a(6))
        setb!(7, b(7) + a(5) * cip * deltat)
    elseif b(6) <= 3.5
        setb!(5, b(5) + (vgap * corrected_current + b(2) * b(3)) * delta2)
    end
    setb!(4, b(3))
    setb!(3, carst)
    setb!(2, vgap)
    setb!(1, b(10) * (avt + art * corrected_current))
    setb!(8, vblock + corrected_current / ylb)
    if sign < 0.0
        setb!(3, -b(3))
        setb!(1, -b(1))
        carst = -carst
    end
    setb!(9, b(9) + b(1) * b(3) * deltat)
    return carst
end

function _over16_simultaneous_zno_record_chain(
    subnetwork_next_indices::AbstractVector{Int},
    subsystem_begin_index::Int,
)
    1 <= subsystem_begin_index <= length(subnetwork_next_indices) ||
        throw(ArgumentError("subsystem_begin_index must address subnetwork_next_indices"))
    records = Int[]
    seen = Set{Int}()
    record = subsystem_begin_index
    while true
        record in seen && throw(ArgumentError("subnetwork records must close at the subsystem head"))
        push!(seen, record)
        push!(records, record)
        next_record = subnetwork_next_indices[record]
        1 <= next_record <= length(subnetwork_next_indices) ||
            throw(ArgumentError("subnetwork_next_indices entries must address subnetwork records"))
        next_record == subsystem_begin_index && break
        record = next_record
    end
    return records
end

function _over16_simultaneous_zno_thevenin_matrix(
    znonl::AbstractVector{<:Real},
    subnetwork_from_nodes::AbstractVector{Int},
    subnetwork_to_nodes::AbstractVector{Int},
    records::Vector{Int},
    ntot::Int,
)
    element_count = length(records)
    zthevenin = zeros(Float64, element_count, element_count)
    for (component, _) in enumerate(records)
        offset = (component - 1) * ntot
        for (row, record) in enumerate(records)
            from_node = subnetwork_from_nodes[record]
            to_node = subnetwork_to_nodes[record]
            zthevenin[row, component] = znonl[from_node + offset] - znonl[to_node + offset]
        end
    end
    return zthevenin
end

function _over16_simultaneous_zno_column_status(zthevenin::AbstractMatrix{Float64})
    element_count = size(zthevenin, 2)
    singular = [all(iszero, view(zthevenin, :, column)) for column in 1:element_count]
    dependent = zeros(Int, element_count)
    for column in 1:(element_count - 1)
        singular[column] && continue
        for candidate in (column + 1):element_count
            (singular[candidate] || dependent[candidate] != 0) && continue
            if zthevenin[:, column] == zthevenin[:, candidate]
                dependent[candidate] = column
            end
        end
    end
    return singular, dependent
end

function _over16_zno_current_and_derivative(
    voltage::Float64,
    reference_voltage::Float64,
    last_segment::Int,
    table_start::Int,
    table_end::Int,
    cchar::Vector{Float64},
    gslope::Vector{Float64},
    vchar::Vector{Float64},
)
    reference_voltage != 0.0 ||
        throw(ArgumentError("simultaneous ZNO ANONL reference voltage must not be zero"))
    1 <= table_start <= table_end <= length(cchar) ||
        throw(ArgumentError("simultaneous ZNO table range must address cchar"))
    table_end <= length(vchar) ||
        throw(ArgumentError("simultaneous ZNO table range must address vchar"))
    table_end <= length(gslope) ||
        throw(ArgumentError("simultaneous ZNO table range must address gslope"))
    abs(last_segment) >= table_start ||
        throw(ArgumentError("simultaneous ZNO ilast must be inside or beyond table range"))

    scaled_voltage = abs(voltage / reference_voltage)
    if last_segment < 0
        scan_start = -last_segment + 1
        scan_end = table_end
    else
        scan_start = table_start
        scan_end = last_segment
    end
    scan_start += 1
    scan_start <= scan_end + 1 ||
        throw(ArgumentError("simultaneous ZNO segment scan range is invalid"))

    search_steps = 0
    segment = scan_end
    for candidate in scan_start:scan_end
        search_steps += 1
        if vchar[candidate] >= scaled_voltage
            segment = candidate - 1
            break
        end
    end
    segment = clamp(segment, table_start, table_end)
    coefficient = cchar[segment]
    if segment < scan_start
        current = coefficient * voltage
        derivative = coefficient
    else
        exponent = gslope[segment]
        magnitude = coefficient * scaled_voltage^exponent
        current = voltage < 0.0 ? -magnitude : magnitude
        derivative = voltage == 0.0 ? coefficient : exponent * current / voltage
    end
    return current, derivative, search_steps
end

function _over16_simultaneous_nonlinear_current_and_derivative(
    element_type::Int,
    voltage::Float64,
    reference_voltage::Float64,
    last_segment::Int,
    table_start::Int,
    table_end::Int,
    cchar::Vector{Float64},
    gslope::Vector{Float64},
    vchar::Vector{Float64},
    time_state::Float64,
    t::Float64,
    epsiln::Float64,
)
    if element_type == 1
        current, derivative, search_steps = _over16_zno_current_and_derivative(
            voltage,
            reference_voltage,
            last_segment,
            table_start,
            table_end,
            cchar,
            gslope,
            vchar,
        )
        return current, derivative, last_segment, search_steps
    end
    element_type in (2, 3) ||
        throw(ArgumentError("simultaneous nonlinear element type must be 1, 2, or 3"))
    if last_segment > 0
        return 0.0, 0.0, last_segment, 0
    end
    1 <= table_start <= table_end <= length(cchar) ||
        throw(ArgumentError("simultaneous nonlinear table range must address cchar"))
    table_end <= length(vchar) ||
        throw(ArgumentError("simultaneous nonlinear table range must address vchar"))
    table_end <= length(gslope) ||
        throw(ArgumentError("simultaneous nonlinear table range must address gslope"))

    table_coordinate = element_type == 3 ? t + time_state : voltage
    scan_coordinate =
        element_type == 2 && vchar[table_start] == 0.0 ?
        abs(table_coordinate) :
        table_coordinate
    segment, stored_segment, search_steps = _over16_piecewise_segment_index(
        scan_coordinate,
        -last_segment,
        table_start,
        table_end,
        vchar,
    )
    slope = gslope[segment]
    value = slope * scan_coordinate + cchar[segment]
    if scan_coordinate * table_coordinate < 0.0
        value = -value
    end
    updated_last = -stored_segment
    if element_type == 2
        return value, slope, updated_last, search_steps
    end
    resistance = value == 0.0 ? epsiln : value
    conductance = 1.0 / resistance
    return conductance * voltage, conductance, updated_last, search_steps
end

function _over16_piecewise_segment_index(
    coordinate::Float64,
    last_segment::Int,
    table_start::Int,
    table_end::Int,
    vchar::Vector{Float64},
)
    last_segment >= table_start - 1 ||
        throw(ArgumentError("simultaneous nonlinear ilast must be at or immediately below the table"))
    search_segment = max(last_segment, table_start)
    if coordinate > vchar[search_segment]
        search_steps = 0
        for candidate in search_segment:table_end
            search_steps += 1
            if coordinate <= vchar[candidate]
                stored_segment = candidate - 1
                return max(stored_segment, table_start), stored_segment, search_steps
            end
        end
        return table_end, table_end, search_steps
    end
    search_steps = 0
    segment = search_segment
    for _ in table_start:search_segment
        search_steps += 1
        if coordinate >= vchar[segment]
            return segment, segment, search_steps
        end
        segment -= 1
    end
    stored_segment = segment
    return max(stored_segment, table_start), stored_segment, search_steps
end

function _over16_check_simultaneous_zno_vector_lengths(
    name::String,
    values::AbstractVector,
    required_length::Int,
)
    length(values) >= required_length ||
        throw(ArgumentError("$name length must cover subnetwork records"))
    return nothing
end

function _check_over16_simultaneous_zno_inputs(
    voltages::AbstractVector{<:Real},
    rhs::AbstractVector{<:Real},
    cchar::AbstractVector{<:Real},
    vchar::AbstractVector{<:Real},
    subnetwork_next_indices::AbstractVector{Int},
    subnetwork_from_nodes::AbstractVector{Int},
    subnetwork_to_nodes::AbstractVector{Int},
    subnetwork_nonlinear_indices::AbstractVector{Int},
    subnetwork_element_types::AbstractVector{Int},
    subsystem_begin_index::Int,
    nonlinear_admittance_start_indices::AbstractVector{Int},
    nonlinear_table_end_indices::AbstractVector{Int},
    initial_anonl::AbstractVector{<:Real},
    initial_vzero::AbstractVector{<:Real},
    initial_ilast::AbstractVector{Int},
    initial_curr::AbstractVector{<:Real},
    initial_cursub::AbstractVector{<:Real},
    vnonl::AbstractVector{<:Real},
    gslope::AbstractVector{<:Real},
    gap_status_values::AbstractVector{<:Real},
    ntot::Int,
)
    length(voltages) == ntot || throw(ArgumentError("voltages length must match ntot"))
    length(rhs) == ntot || throw(ArgumentError("rhs length must match ntot"))
    _check_over16_finite_vector("voltages", voltages)
    _check_over16_finite_vector("rhs", rhs)
    _check_over16_finite_vector("cchar", cchar)
    _check_over16_finite_vector("vchar", vchar)
    isempty(gslope) || _check_over16_finite_vector("gslope", gslope)
    isempty(subnetwork_next_indices) &&
        throw(ArgumentError("subnetwork_next_indices must not be empty"))
    1 <= subsystem_begin_index <= length(subnetwork_next_indices) ||
        throw(ArgumentError("subsystem_begin_index must address subnetwork_next_indices"))
    required_length = length(subnetwork_next_indices)
    _over16_check_simultaneous_zno_vector_lengths(
        "subnetwork_from_nodes",
        subnetwork_from_nodes,
        required_length,
    )
    _over16_check_simultaneous_zno_vector_lengths(
        "subnetwork_to_nodes",
        subnetwork_to_nodes,
        required_length,
    )
    _over16_check_simultaneous_zno_vector_lengths(
        "subnetwork_nonlinear_indices",
        subnetwork_nonlinear_indices,
        required_length,
    )
    _over16_check_simultaneous_zno_vector_lengths(
        "subnetwork_element_types",
        subnetwork_element_types,
        required_length,
    )
    _check_over16_finite_vector("initial_anonl", initial_anonl)
    _check_over16_vzero_vector(initial_vzero)
    isempty(initial_curr) || _check_over16_finite_vector("initial_curr", initial_curr)
    isempty(initial_cursub) || _check_over16_finite_vector("initial_cursub", initial_cursub)
    isempty(vnonl) || all(value -> isfinite(Float64(value)) || Float64(value) == Inf, vnonl) ||
        throw(ArgumentError("vnonl entries must be finite or Inf"))
    count = length(initial_anonl)
    length(initial_vzero) == count ||
        throw(ArgumentError("initial_vzero length must match initial_anonl"))
    length(initial_ilast) == count ||
        throw(ArgumentError("initial_ilast length must match initial_anonl"))
    isempty(initial_curr) || length(initial_curr) == count ||
        throw(ArgumentError("initial_curr length must match initial_anonl"))
    isempty(vnonl) || length(vnonl) == count ||
        throw(ArgumentError("vnonl length must match initial_anonl"))
    isempty(gap_status_values) || length(gap_status_values) == count ||
        throw(ArgumentError("gap_status_values length must match initial_anonl"))
    isempty(gap_status_values) ||
        _check_over16_finite_vector("gap_status_values", gap_status_values)
    length(nonlinear_admittance_start_indices) == count ||
        throw(ArgumentError("nonlinear_admittance_start_indices length must match initial_anonl"))
    length(nonlinear_table_end_indices) == count ||
        throw(ArgumentError("nonlinear_table_end_indices length must match initial_anonl"))

    records = _over16_simultaneous_zno_record_chain(
        subnetwork_next_indices,
        subsystem_begin_index,
    )
    for record in records
        from_node = subnetwork_from_nodes[record]
        to_node = subnetwork_to_nodes[record]
        1 <= from_node <= ntot ||
            throw(ArgumentError("subnetwork_from_nodes entries must address voltages"))
        1 <= to_node <= ntot ||
            throw(ArgumentError("subnetwork_to_nodes entries must address voltages"))
        nonlinear_index = subnetwork_nonlinear_indices[record]
        1 <= nonlinear_index <= count ||
            throw(ArgumentError("subnetwork_nonlinear_indices entries must address nonlinear state"))
        element_type = subnetwork_element_types[record]
        element_type in (1, 2, 3) ||
            throw(ArgumentError("simultaneous nonlinear element type must be 1, 2, or 3"))
        table_start = nonlinear_admittance_start_indices[nonlinear_index]
        table_end = nonlinear_table_end_indices[nonlinear_index]
        1 <= table_start <= table_end <= length(cchar) ||
            throw(ArgumentError("nonlinear ZNO table range must address cchar"))
        table_end <= length(vchar) ||
            throw(ArgumentError("nonlinear ZNO table range must address vchar"))
        if element_type in (2, 3)
            !isempty(gslope) ||
                throw(ArgumentError("gslope must be provided for simultaneous type-2/type-3 elements"))
            table_end <= length(gslope) ||
                throw(ArgumentError("nonlinear type-2/type-3 table range must address gslope"))
        end
        initial_anonl[nonlinear_index] != 0.0 ||
            throw(ArgumentError("initial_anonl entries used by ZNO must not be zero"))
        abs(initial_ilast[nonlinear_index]) >= table_start ||
            throw(ArgumentError("initial_ilast entries used by ZNO must be inside the table"))
    end
    return nothing
end

function _check_over16_nonlinear_dimensions(ntot::Int, ncomp::Int)
    ntot > 0 || throw(ArgumentError("ntot must be positive"))
    ncomp > 0 || throw(ArgumentError("ncomp must be positive"))
    return nothing
end

function _check_over16_znonl_layout(znonl::AbstractVector{<:Real}, ntot::Int, ncomp::Int)
    length(znonl) == ntot * ncomp ||
        throw(ArgumentError("znonl length must equal ntot * ncomp"))
    _check_over16_finite_vector("znonl", znonl)
    return nothing
end

function _check_over16_finite_vector(name::String, values::AbstractVector{<:Real})
    for value in values
        isfinite(Float64(value)) || throw(ArgumentError("$name entries must be finite"))
    end
    return nothing
end

function _check_over16_vzero_vector(values::AbstractVector{<:Real})
    for value in values
        !isnan(Float64(value)) || throw(ArgumentError("vzero entries must not be NaN"))
    end
    return nothing
end

function _check_over16_vnonl_vector(values::AbstractVector{<:Real})
    for value in values
        candidate = Float64(value)
        (isfinite(candidate) || candidate == Inf) ||
            throw(ArgumentError("vnonl entries must be finite or Inf"))
    end
    return nothing
end

function _over16_check_source_column_vector(
    name::String,
    values::AbstractVector{Int},
    expected::Int,
)
    length(values) == expected || throw(ArgumentError("$name length must be $expected"))
    return nothing
end

function _over16_check_nonlinear_source_metadata(
    nonlinear_types::AbstractVector{Int},
    nonlinear_admittance_nodes::AbstractVector{Int},
    initial_vzero::AbstractVector{<:Real},
    initial_ilast::AbstractVector{Int},
)
    count = length(nonlinear_types)
    length(nonlinear_admittance_nodes) == count ||
        throw(ArgumentError("nonlinear_admittance_nodes length must match nonlinear_types"))
    isempty(initial_vzero) || length(initial_vzero) == count ||
        throw(ArgumentError("initial_vzero length must match nonlinear_types"))
    isempty(initial_ilast) || length(initial_ilast) == count ||
        throw(ArgumentError("initial_ilast length must match nonlinear_types"))
    isempty(initial_vzero) || _check_over16_vzero_vector(initial_vzero)
    return nothing
end

function _over16_source_column_representative(node::Int, kode::AbstractVector{Int})
    target = node
    redirects = 0
    seen = Set{Int}()
    while kode[target] >= target
        target in seen && throw(ArgumentError("kode source-column chain must not repeat"))
        push!(seen, target)
        target = kode[target]
        1 <= target <= length(kode) ||
            throw(ArgumentError("kode source-column target must be within ntot"))
        redirects += 1
    end
    return target, redirects
end

function _check_over16_nonlinear_sparse_workspace(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int},
    ntot::Int,
    partition_boundary::Int,
    iupper::Int,
)
    length(km) == length(ykm) ||
        throw(ArgumentError("km and ykm lengths must match"))
    length(kk) == ntot || throw(ArgumentError("kk length must match ntot"))
    1 <= partition_boundary <= ntot ||
        throw(ArgumentError("partition_boundary must be within ntot"))
    1 <= iupper <= length(km) ||
        throw(ArgumentError("iupper must be within km/ykm length"))
    _check_over16_finite_vector("ykm", ykm)

    index = 1
    while index <= iupper
        marker = km[index]
        marker < 0 || throw(ArgumentError("each sparse row must start with a negative diagonal marker"))
        row = abs(marker)
        1 <= row <= partition_boundary ||
            throw(ArgumentError("sparse row diagonal must be within partition_boundary"))
        row_end = abs(kk[row])
        index <= row_end <= iupper ||
            throw(ArgumentError("kk row end must cover the sparse row within iupper"))
        for entry in index:row_end
            column_marker = km[entry]
            column_marker != 0 || throw(ArgumentError("km entries must not be zero"))
            column = abs(column_marker)
            1 <= column <= ntot ||
                throw(ArgumentError("km entries must be within ntot"))
            if entry > index
                column_marker > 0 ||
                    throw(ArgumentError("off-diagonal km entries must be positive"))
            end
        end
        index = row_end + 1
    end
    index == iupper + 1 ||
        throw(ArgumentError("sparse rows must consume entries through iupper"))
    return nothing
end

function _check_over16_kode(kode::AbstractVector{Int}, ntot::Int)
    length(kode) == ntot || throw(ArgumentError("kode length must match ntot"))
    for value in kode
        value == 0 || 1 <= value <= ntot ||
            throw(ArgumentError("kode entries must be zero or within ntot"))
    end
    return nothing
end

function _over16_nonlinear_difference_inputs(
    nonlinear_types::AbstractVector{Int},
    nonlinear_from_nodes::AbstractVector{Int},
    nonlinear_to_nodes::AbstractVector{Int},
    nonlinear_source_flags::AbstractVector{Int},
    initial_anonl::AbstractVector{<:Real},
    ntot::Int,
)
    count = length(nonlinear_types)
    length(nonlinear_from_nodes) == count ||
        throw(ArgumentError("nonlinear_from_nodes length must match nonlinear_types"))
    length(nonlinear_to_nodes) == count ||
        throw(ArgumentError("nonlinear_to_nodes length must match nonlinear_types"))
    flags = isempty(nonlinear_source_flags) ? zeros(Int, count) : collect(nonlinear_source_flags)
    length(flags) == count ||
        throw(ArgumentError("nonlinear_source_flags length must match nonlinear_types"))
    anonl = isempty(initial_anonl) ? zeros(Float64, count) : Float64.(initial_anonl)
    length(anonl) == count ||
        throw(ArgumentError("initial_anonl length must match nonlinear_types"))
    _check_over16_finite_vector("initial_anonl", anonl)

    for index in eachindex(nonlinear_types)
        nonlinear_type = nonlinear_types[index]
        if nonlinear_type < 0 || nonlinear_type > 920
            continue
        end
        from_node = nonlinear_from_nodes[index]
        to_node = abs(nonlinear_to_nodes[index])
        1 <= from_node <= ntot ||
            throw(ArgumentError("nonlinear_from_nodes entries must be within ntot"))
        1 <= to_node <= ntot ||
            throw(ArgumentError("nonlinear_to_nodes entries must be within ntot"))
    end
    return anonl, flags
end

function _over16_check_nonlinear_current_inputs(
    voltages::AbstractVector{<:Real},
    rhs::AbstractVector{<:Real},
    nonlinear_types::AbstractVector{Int},
    nonlinear_from_nodes::AbstractVector{Int},
    nonlinear_to_nodes::AbstractVector{Int},
    nonlinear_admittance_nodes::AbstractVector{Int},
    nonlinear_table_end_indices::AbstractVector{Int},
    nonlinear_subsystem_indices::AbstractVector{Int},
    subsystem_begin_indices::AbstractVector{Int},
    subsystem_owner_rows::AbstractVector{Int},
    subsystem_simultaneous_flags::AbstractVector{Int},
    cchar::AbstractVector{<:Real},
    vchar::AbstractVector{<:Real},
    initial_anonl::AbstractVector{<:Real},
    initial_vzero::AbstractVector{<:Real},
    initial_ilast::AbstractVector{Int},
    initial_curr::AbstractVector{<:Real},
    initial_cursub::AbstractVector{<:Real},
    vnonl::AbstractVector{<:Real},
)
    count = length(nonlinear_types)
    length(nonlinear_from_nodes) == count ||
        throw(ArgumentError("nonlinear_from_nodes length must match nonlinear_types"))
    length(nonlinear_to_nodes) == count ||
        throw(ArgumentError("nonlinear_to_nodes length must match nonlinear_types"))
    length(nonlinear_admittance_nodes) == count ||
        throw(ArgumentError("nonlinear_admittance_nodes length must match nonlinear_types"))
    length(nonlinear_table_end_indices) == count ||
        throw(ArgumentError("nonlinear_table_end_indices length must match nonlinear_types"))
    length(nonlinear_subsystem_indices) == count ||
        throw(ArgumentError("nonlinear_subsystem_indices length must match nonlinear_types"))
    length(initial_anonl) == count ||
        throw(ArgumentError("initial_anonl length must match nonlinear_types"))
    length(initial_vzero) == count ||
        throw(ArgumentError("initial_vzero length must match nonlinear_types"))
    length(initial_ilast) == count ||
        throw(ArgumentError("initial_ilast length must match nonlinear_types"))
    isempty(initial_curr) || length(initial_curr) == count ||
        throw(ArgumentError("initial_curr length must match nonlinear_types"))
    isempty(vnonl) || length(vnonl) == count ||
        throw(ArgumentError("vnonl length must match nonlinear_types"))
    !isempty(subsystem_begin_indices) ||
        throw(ArgumentError("subsystem_begin_indices must not be empty"))
    length(subsystem_owner_rows) == length(subsystem_simultaneous_flags) ||
        throw(ArgumentError("subsystem owner rows and flags lengths must match"))
    length(voltages) == length(rhs) ||
        throw(ArgumentError("voltages and rhs lengths must match"))
    _check_over16_finite_vector("voltages", voltages)
    _check_over16_finite_vector("rhs", rhs)
    _check_over16_finite_vector("cchar", cchar)
    _check_over16_finite_vector("vchar", vchar)
    _check_over16_finite_vector("initial_anonl", initial_anonl)
    _check_over16_vzero_vector(initial_vzero)
    isempty(initial_curr) || _check_over16_finite_vector("initial_curr", initial_curr)
    isempty(initial_cursub) || _check_over16_finite_vector("initial_cursub", initial_cursub)
    isempty(vnonl) || _check_over16_vnonl_vector(vnonl)
    for subsystem_index in nonlinear_subsystem_indices
        1 <= subsystem_index <= length(subsystem_begin_indices) ||
            throw(ArgumentError("nonlinear_subsystem_indices entries must address subsystem_begin_indices"))
    end
    return nothing
end

function _over16_required_cursub_length(
    nonlinear_subsystem_indices::AbstractVector{Int},
    subsystem_begin_indices::AbstractVector{Int},
)
    required = 0
    for subsystem_index in nonlinear_subsystem_indices
        head = subsystem_begin_indices[subsystem_index]
        head <= 0 && continue
        required = max(required, div(head, 5) + 1)
    end
    return required
end

function _over16_interpolated_table_current(
    h1::Float64,
    h2::Float64,
    lower_index::Int,
    cchar::Vector{Float64},
)
    denominator = abs(h1) + abs(h2)
    denominator > 0.0 ||
        return cchar[lower_index]
    delta = abs(h1) * (cchar[lower_index + 1] - cchar[lower_index]) / denominator
    return cchar[lower_index] + delta
end

function _over16_find_increasing_nonlinear_current(
    voltage_delta::Float64,
    slope::Float64,
    last::Int,
    table_end::Int,
    cchar::Vector{Float64},
    vchar::Vector{Float64},
)
    h1 = voltage_delta + slope * cchar[last] - vchar[last]
    steps = 0
    while true
        last < table_end ||
            throw(ArgumentError("nonlinear current interpolation exceeded upper table limit"))
        h2 = voltage_delta + slope * cchar[last + 1] - vchar[last + 1]
        steps += 1
        if h1 * h2 <= 0.0
            return _over16_interpolated_table_current(h1, h2, last, cchar), last, steps
        end
        last += 1
        h1 = h2
    end
end

function _over16_find_decreasing_nonlinear_current(
    voltage_delta::Float64,
    slope::Float64,
    last::Int,
    table_start::Int,
    cchar::Vector{Float64},
    vchar::Vector{Float64},
)
    h2 = voltage_delta + slope * cchar[last + 1] - vchar[last + 1]
    steps = 0
    while true
        last >= table_start ||
            throw(ArgumentError("nonlinear current interpolation exceeded lower table limit"))
        h1 = voltage_delta + slope * cchar[last] - vchar[last]
        steps += 1
        if h1 * h2 <= 0.0
            return _over16_interpolated_table_current(h1, h2, last, cchar), last, steps
        end
        h2 = h1
        last -= 1
    end
end

function _over16_nonlinear_column_index(row::Int, component::Int, ntot::Int)
    return row + (component - 1) * ntot
end

function _over16_copy_kode_group!(
    znonl::Vector{Float64},
    kode::AbstractVector{Int},
    source_row::Int,
    ntot::Int,
    ncomp::Int,
)
    first_target = kode[source_row]
    first_target == 0 && return 0
    target = first_target
    copied = 0
    seen = Set{Int}()
    while true
        1 <= target <= ntot ||
            throw(ArgumentError("kode group target must be within ntot"))
        for component in 1:ncomp
            znonl[_over16_nonlinear_column_index(target, component, ntot)] =
                znonl[_over16_nonlinear_column_index(source_row, component, ntot)]
        end
        copied += 1
        kode[target] == source_row && return copied
        target in seen && throw(ArgumentError("kode group must close back to source row"))
        push!(seen, target)
        target = kode[target]
        target != 0 || throw(ArgumentError("kode group must close back to source row"))
    end
end

mutable struct DiodeValveSwitch <: EMTElement
    a::Int
    b::Int
    threshold_v::Float64
    holding_current::Float64
    on_conductance::Float64
    off_conductance::Float64
    closed::Bool
    last_voltage::Float64
    last_current::Float64
    last_conductance::Float64
end

function DiodeValveSwitch(
    a::Int,
    b::Int;
    threshold_v::Real=0.0,
    holding_current::Real=0.0,
    on_conductance::Real=1.0e9,
    off_conductance::Real=0.0,
    initially_closed::Bool=false,
)
    a >= 0 && b >= 0 || throw(ArgumentError("diode nodes must be nonnegative"))
    threshold = Float64(threshold_v)
    holding = Float64(holding_current)
    on = Float64(on_conductance)
    off = Float64(off_conductance)
    isfinite(threshold) || throw(ArgumentError("threshold_v must be finite"))
    isfinite(holding) && holding >= 0.0 ||
        throw(ArgumentError("holding_current must be finite and nonnegative"))
    isfinite(on) && isfinite(off) && on >= 0.0 && off >= 0.0 ||
        throw(ArgumentError("diode conductances must be finite and nonnegative"))
    on >= off || throw(ArgumentError("on_conductance must be at least off_conductance"))
    conductance = initially_closed ? on : off
    return DiodeValveSwitch(a, b, threshold, holding, on, off, initially_closed, 0.0, 0.0, conductance)
end

function diode_next_closed(s::DiodeValveSwitch, voltage::Real, current::Real)::Bool
    return s.closed ? Float64(current) >= s.holding_current : Float64(voltage) >= s.threshold_v
end

diode_conductance(s::DiodeValveSwitch)::Float64 = s.closed ? s.on_conductance : s.off_conductance

function stamp!(
    y::AbstractMatrix{Float64},
    _rhs::AbstractVector{Float64},
    s::DiodeValveSwitch,
    _t::Float64,
    _dt::Float64,
)
    conductance = diode_conductance(s)
    s.last_conductance = conductance
    stamp_conductance!(y, s.a, s.b, conductance)
    return nothing
end

function update!(s::DiodeValveSwitch, voltages::AbstractVector{Float64}, _dt::Float64)
    voltage = Branches.branch_voltage(voltages, s.a, s.b)
    current = s.last_conductance * voltage
    s.last_voltage = voltage
    s.last_current = current
    s.closed = diode_next_closed(s, voltage, current)
    s.last_conductance = diode_conductance(s)
    return nothing
end

mutable struct SaturableInductorBranch <: EMTElement
    a::Int
    b::Int
    r::Float64
    l_unsaturated::Float64
    l_saturated::Float64
    current_threshold::Float64
    i_prev::Float64
    v_prev::Float64
    i_last::Float64
    last_inductance::Float64
    saturated::Bool
end

function SaturableInductorBranch(
    a::Int,
    b::Int,
    r::Real,
    l_unsaturated::Real,
    l_saturated::Real,
    current_threshold::Real;
    i_prev::Real=0.0,
    v_prev::Real=0.0,
)
    a >= 0 && b >= 0 || throw(ArgumentError("inductor nodes must be nonnegative"))
    resistance = Float64(r)
    l_u = Float64(l_unsaturated)
    l_s = Float64(l_saturated)
    threshold = Float64(current_threshold)
    i0 = Float64(i_prev)
    v0 = Float64(v_prev)
    isfinite(resistance) && resistance >= 0.0 ||
        throw(ArgumentError("r must be finite and nonnegative"))
    isfinite(l_u) && isfinite(l_s) && l_u > 0.0 && l_s > 0.0 ||
        throw(ArgumentError("inductances must be finite and positive"))
    isfinite(threshold) && threshold >= 0.0 ||
        throw(ArgumentError("current_threshold must be finite and nonnegative"))
    isfinite(i0) && isfinite(v0) ||
        throw(ArgumentError("initial current and voltage must be finite"))
    l0 = abs(i0) >= threshold ? l_s : l_u
    return SaturableInductorBranch(a, b, resistance, l_u, l_s, threshold, i0, v0, i0, l0, l0 == l_s)
end

function effective_inductance(b::SaturableInductorBranch, current::Real)::Float64
    return abs(Float64(current)) >= b.current_threshold ? b.l_saturated : b.l_unsaturated
end

effective_inductance(b::SaturableInductorBranch)::Float64 = effective_inductance(b, b.i_prev)

function saturable_inductor_companion(b::SaturableInductorBranch, dt::Float64)
    inductance = effective_inductance(b)
    g, ih = Companion.series_rl_companion(b.r, inductance, b.i_prev, b.v_prev, dt)
    return g, ih, inductance
end

function stamp!(
    y::AbstractMatrix{Float64},
    rhs::AbstractVector{Float64},
    b::SaturableInductorBranch,
    _t::Float64,
    dt::Float64,
)
    g, ih, inductance = saturable_inductor_companion(b, dt)
    b.last_inductance = inductance
    stamp_conductance!(y, b.a, b.b, g)
    stamp_history_current!(rhs, b.a, b.b, ih)
    return nothing
end

function update!(b::SaturableInductorBranch, voltages::AbstractVector{Float64}, dt::Float64)
    g, ih, inductance = saturable_inductor_companion(b, dt)
    voltage = Branches.branch_voltage(voltages, b.a, b.b)
    current = g * voltage + ih
    b.v_prev = voltage
    b.i_prev = current
    b.i_last = current
    b.last_inductance = inductance
    b.saturated = abs(current) >= b.current_threshold
    return nothing
end

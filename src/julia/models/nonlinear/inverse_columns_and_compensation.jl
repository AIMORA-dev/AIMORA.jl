
function saturated_transformer_sparse_admittance_update(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kks::AbstractVector{Int},
    from_node::Int,
    to_node::Int,
    partition_boundary::Int,
    admittance_delta::Real,
)
    length(km) == length(ykm) ||
        throw(ArgumentError("km and ykm lengths must match"))
    !isempty(km) || throw(ArgumentError("km must not be empty"))
    !isempty(kks) || throw(ArgumentError("kks must not be empty"))
    1 <= from_node <= length(kks) ||
        throw(ArgumentError("from_node must address kks"))
    1 <= to_node <= length(kks) ||
        throw(ArgumentError("to_node must address kks"))
    from_node != to_node || throw(ArgumentError("saturated transformer branch must not be a self-loop"))
    1 <= partition_boundary <= length(kks) ||
        throw(ArgumentError("partition_boundary must address kks"))
    _check_over16_finite_vector("ykm", ykm)
    delta = Float64(admittance_delta)
    isfinite(delta) || throw(ArgumentError("admittance_delta must be finite"))
    for pointer in kks
        1 <= pointer <= length(km) + 1 ||
            throw(ArgumentError("kks entries must point one past a sparse row"))
    end

    work = Float64.(ykm)
    lower_node = min(from_node, to_node)
    upper_node = max(from_node, to_node)
    positive_indices = Int[]
    negative_indices = Int[]
    partition_skip_count = 0
    row_scan_count = 0
    retriangularization_required = false

    if lower_node > partition_boundary
        partition_skip_count += 1
    else
        if lower_node != 1
            row_scan_count += 1
            index = kks[lower_node]
            found_coupling = false
            while index > 1
                index -= 1
                node = abs(km[index])
                if node == lower_node
                    work[index] += delta
                    push!(positive_indices, index)
                elseif node == upper_node
                    work[index] -= delta
                    push!(negative_indices, index)
                    found_coupling = true
                    break
                end
            end
            found_coupling ||
                throw(ArgumentError("lower sparse row does not contain upper branch node"))
        end

        if upper_node > partition_boundary
            partition_skip_count += 1
        else
            row_scan_count += 1
            index = kks[upper_node]
            found_diagonal = false
            while index > 1
                index -= 1
                node = abs(km[index])
                if node == upper_node
                    work[index] += delta
                    push!(positive_indices, index)
                    found_diagonal = true
                    break
                elseif node == lower_node
                    work[index] -= delta
                    push!(negative_indices, index)
                end
            end
            found_diagonal ||
                throw(ArgumentError("upper sparse row does not contain its diagonal"))
        end
        retriangularization_required = true
    end

    return (
        source = :saturated_transformer_sparse_admittance_update,
        outcome = :state_mutation,
        fortran_files = (:OVER16_FOR,),
        fortran_routines = (:SUBTS3,),
        fortran_labels = OVER16_SATURATED_TRANSFORMER_ADMITTANCE_LABELS,
        km = collect(km),
        ykm = work,
        kks = collect(kks),
        from_node = from_node,
        to_node = to_node,
        partition_boundary = partition_boundary,
        admittance_delta = delta,
        positive_update_indices = positive_indices,
        negative_update_indices = negative_indices,
        positive_update_count = length(positive_indices),
        negative_update_count = length(negative_indices),
        partition_skip_count = partition_skip_count,
        row_scan_count = row_scan_count,
        retriangularization_required = retriangularization_required,
        sparse_admittance_mutated = work != Float64.(ykm),
        mutation_order = (
            :lower_sparse_row_scan,
            :upper_sparse_row_scan,
            :retriangularization_request,
        ),
        deferred_calls = retriangularization_required ?
            [:sparse_factor_update, :retriangularization_execution, :bulk_last14_oracle] :
            [:bulk_last14_oracle],
        replacement_ready = false,
    )
end

function _over16_integer_table_pointer(value::Float64, label::AbstractString)
    rounded = round(Int, value)
    abs(value - rounded) <= 1.0e-9 ||
        throw(ArgumentError("$label must contain an integer table pointer"))
    return rounded
end

function _over16_first_table_index(
    start_index::Int,
    stop_index::Int,
    predicate::Function,
)
    for candidate in start_index:stop_index
        predicate(candidate) && return candidate
    end
    return stop_index + 1
end

function hysteretic_inductor_current_update(
    cchar::AbstractVector{<:Real},
    vchar::AbstractVector{<:Real},
    gslope::AbstractVector{<:Real};
    hysteresis_state_start_index::Int,
    major_loop_start_index::Int,
    stored_flux::Real,
    branch_voltage::Real,
    delta2::Real,
    flzero::Real=0.0,
)
    n6 = hysteresis_state_start_index
    n9 = major_loop_start_index
    n6 >= 1 || throw(ArgumentError("hysteresis_state_start_index must be positive"))
    n9 >= 1 || throw(ArgumentError("major_loop_start_index must be positive"))
    n6 + 5 <= length(cchar) ||
        throw(ArgumentError("hysteretic inductor state must cover CCHAR(N6:N6+5)"))
    n6 + 5 <= length(vchar) ||
        throw(ArgumentError("hysteretic inductor state must cover VCHAR(N6:N6+5)"))
    n6 + 5 <= length(gslope) ||
        throw(ArgumentError("hysteretic inductor state must cover GSLOPE(N6:N6+5)"))
    cchar_values = Float64.(cchar)
    vchar_values = Float64.(vchar)
    gslope_values = Float64.(gslope)
    _check_over16_finite_vector("cchar", cchar_values)
    _check_over16_finite_vector("vchar", vchar_values)
    _check_over16_finite_vector("gslope", gslope_values)
    n7 = _over16_integer_table_pointer(cchar_values[n6], "hysteretic inductor point count")
    n7 > 0 || throw(ArgumentError("hysteretic inductor point count must be positive"))
    n10 = n9 + n7 - 1
    n12 = n6 + 2
    maximum_required = n10 + n7 + 2
    n10 <= length(cchar_values) ||
        throw(ArgumentError("hysteretic inductor major-loop current table is incomplete"))
    maximum_required <= length(cchar_values) ||
        throw(ArgumentError("hysteretic inductor CCHAR extension table is incomplete"))
    maximum_required <= length(vchar_values) ||
        throw(ArgumentError("hysteretic inductor VCHAR extension table is incomplete"))
    maximum_required <= length(gslope_values) ||
        throw(ArgumentError("hysteretic inductor GSLOPE extension table is incomplete"))
    delta = Float64(delta2)
    voltage = Float64(branch_voltage)
    flux = Float64(stored_flux)
    zero = Float64(flzero)
    isfinite(delta) && delta > 0.0 ||
        throw(ArgumentError("delta2 must be finite and positive"))
    isfinite(voltage) || throw(ArgumentError("branch_voltage must be finite"))
    isfinite(flux) || throw(ArgumentError("stored_flux must be finite"))
    isfinite(zero) && zero >= 0.0 || throw(ArgumentError("flzero must be finite and nonnegative"))

    flux += voltage * delta
    cchar_values[n6 + 3] = voltage * gslope_values[n6 + 1] + gslope_values[n6]
    branch_current = cchar_values[n6 + 3]
    if cchar_values[n6 + 4] < 0.0 && voltage > zero
        cchar_values[n6 + 4] = 0.0
    end

    terminal_branch = :upper
    reversal_count = 0
    trajectory_limit_count = 0
    d9 = 0.0

    function set_trace_index!(index::Int)
        cchar_values[n12] = Float64(index)
        return index
    end

    function downer_final!()
        n14 = _over16_integer_table_pointer(cchar_values[n12], "hysteretic inductor trace index")
        n13 = n14 + n7 + 1
        1 <= n14 <= length(gslope_values) ||
            throw(ArgumentError("hysteretic inductor trace index must address GSLOPE"))
        n13 <= length(gslope_values) ||
            throw(ArgumentError("hysteretic inductor trace extension must address GSLOPE"))
        if cchar_values[n6 + 5] == 1.0
            d13 = vchar_values[n13]
            d14 = cchar_values[n13]
        else
            denominator = gslope_values[n14] * (1.0 + vchar_values[n6])
            denominator != 0.0 ||
                throw(ArgumentError("hysteretic inductor downer denominator must not be zero"))
            d13 = 1.0 / denominator
            d14 = (gslope_values[n13] - gslope_values[n14] * vchar_values[n6 + 1]) * d13
        end
        return d13, d14
    end

    function upper_final!()
        n14 = _over16_integer_table_pointer(cchar_values[n12], "hysteretic inductor trace index")
        n13 = n14 + n7 + 1
        1 <= n14 <= length(gslope_values) ||
            throw(ArgumentError("hysteretic inductor trace index must address GSLOPE"))
        n13 <= length(gslope_values) ||
            throw(ArgumentError("hysteretic inductor trace extension must address GSLOPE"))
        if cchar_values[n6 + 5] == 1.0
            d13 = vchar_values[n13]
            d14 = -cchar_values[n13]
        else
            denominator = gslope_values[n14] * (1.0 - vchar_values[n6])
            denominator != 0.0 ||
                throw(ArgumentError("hysteretic inductor upper denominator must not be zero"))
            d13 = 1.0 / denominator
            d14 = (gslope_values[n14] * vchar_values[n6 + 1] - gslope_values[n13]) * d13
        end
        return d13, d14
    end

    function rebuild_minor_loop!(d6::Float64, d9::Float64, d11::Float64, d13::Float64)
        d12 = d11
        if !(vchar_values[n6 + 4] < vchar_values[n6 + 5] - zero) &&
           !(vchar_values[n6 + 4] > vchar_values[n6 + 5] + zero)
            if cchar_values[n6 + 1] == 1.0
                vchar_values[n6 + 5] = vchar_values[n10]
                gslope_values[n6 + 5] = cchar_values[n10]
            else
                vchar_values[n6 + 5] = vchar_values[n9]
                gslope_values[n6 + 5] = cchar_values[n9]
            end
            d13 = 0.0
        end
        denominator = vchar_values[n10] - vchar_values[n6 + 4]
        d15 = d12 * (vchar_values[n10] - vchar_values[n6 + 5]) / denominator
        if cchar_values[n6 + 1] == -1.0
            denominator = -vchar_values[n10] - vchar_values[n6 + 4]
            d15 = d12 * (-vchar_values[n10] - vchar_values[n6 + 5]) / denominator
        end
        if d13 > d15
            d13 = d15
            trajectory_limit_count += 1
        end
        minor_denominator = vchar_values[n6 + 4] - vchar_values[n6 + 5]
        minor_denominator != 0.0 ||
            throw(ArgumentError("hysteretic inductor minor-loop denominator must not be zero"))
        vchar_values[n6] = (d12 - d13) / minor_denominator
        vchar_values[n6 + 1] = d12 - vchar_values[n6] * vchar_values[n6 + 4]
        return cchar_values[n6 + 1] == 1.0 ? upper_final!() : downer_final!()
    end

    if (flux + zero < vchar_values[n6 + 2] && cchar_values[n6 + 1] == 1.0) ||
       (flux - zero > vchar_values[n6 + 2] && cchar_values[n6 + 1] == -1.0)
        reversal_count += 1
        cchar_values[n6 + 1] = -cchar_values[n6 + 1]
        cchar_values[n6 + 4] += 1.0
        cchar_values[n6 + 5] = 0.0
        if cchar_values[n6 + 1] == 1.0
            if vchar_values[n6 + 2] <= -vchar_values[n10]
                vchar_values[n6 + 5] = vchar_values[n10]
                gslope_values[n6 + 5] = cchar_values[n10]
                vchar_values[n6] = 0.0
                vchar_values[n6 + 1] = 0.0
                vchar_values[n6 + 4] = vchar_values[n9]
                gslope_values[n6 + 4] = cchar_values[n9]
                set_trace_index!(n9)
                terminal_branch = :upper
            else
                if cchar_values[n6 + 4] <= 1.0
                    vchar_values[n6 + 4] = vchar_values[n10]
                    gslope_values[n6 + 4] = cchar_values[n10]
                    d6 = 0.0
                else
                    index = _over16_first_table_index(
                        n9,
                        n10,
                        candidate -> !(gslope_values[n6 + 4] > cchar_values[candidate]),
                    )
                    set_trace_index!(index)
                    n13 = index + n7 + 1
                    d9 = vchar_values[n13] * gslope_values[n6 + 4] + cchar_values[n13]
                    d6 = vchar_values[n6 + 4] - d9
                end
                index = _over16_first_table_index(
                    n9,
                    n10,
                    candidate -> !(vchar_values[n6 + 3] > cchar_values[candidate]),
                )
                set_trace_index!(index)
                n13 = index + n7 + 1
                d10 = vchar_values[n13] * vchar_values[n6 + 3] + cchar_values[n13]
                d11 = vchar_values[n6 + 2] - d10
                n14 = _over16_first_table_index(
                    n9,
                    n10,
                    candidate -> !(gslope_values[n6 + 5] > cchar_values[candidate]),
                )
                n14_extension = n14 + n7 + 1
                d10 = vchar_values[n14_extension] * gslope_values[n6 + 5] +
                    cchar_values[n14_extension]
                d13 = vchar_values[n6 + 5] - d10
                if vchar_values[n6 + 2] < vchar_values[n6 + 5]
                    d14 = 0.0
                    abs(d13) > zero && (d14 = d11 * d6 / d13)
                    d9 += d14
                    vchar_values[n6 + 5] = d9
                    gslope_values[n6 + 5] = gslope_values[n6 + 4]
                    d13 = d14
                else
                    vchar_values[n6 + 5] = vchar_values[n6 + 4]
                    gslope_values[n6 + 5] = gslope_values[n6 + 4]
                    d13 = d6
                end
                vchar_values[n6 + 4] = vchar_values[n6 + 2]
                gslope_values[n6 + 4] = vchar_values[n6 + 3]
                d13, d14 = rebuild_minor_loop!(d6, d9, d11, d13)
                terminal_branch = cchar_values[n6 + 1] == 1.0 ? :upper : :downer
            end
        else
            if vchar_values[n6 + 2] >= vchar_values[n10]
                vchar_values[n6 + 5] = vchar_values[n9]
                gslope_values[n6 + 5] = cchar_values[n9]
                vchar_values[n6] = 0.0
                vchar_values[n6 + 1] = 0.0
                vchar_values[n6 + 4] = vchar_values[n10]
                gslope_values[n6 + 4] = cchar_values[n10]
                set_trace_index!(n9)
                terminal_branch = :downer
            else
                if cchar_values[n6 + 4] <= 1.0
                    vchar_values[n6 + 4] = vchar_values[n9]
                    gslope_values[n6 + 4] = cchar_values[n9]
                    d6 = 0.0
                else
                    index = _over16_first_table_index(
                        n9,
                        n10,
                        candidate -> !(-gslope_values[n6 + 4] > cchar_values[candidate]),
                    )
                    set_trace_index!(index)
                    n13 = index + n7 + 1
                    d9 = vchar_values[n13] * gslope_values[n6 + 4] - cchar_values[n13]
                    d6 = d9 - vchar_values[n6 + 4]
                end
                index = _over16_first_table_index(
                    n9,
                    n10,
                    candidate -> !(-vchar_values[n6 + 3] > cchar_values[candidate]),
                )
                set_trace_index!(index)
                n13 = index + n7 + 1
                d10 = vchar_values[n13] * vchar_values[n6 + 3] - cchar_values[n13]
                d11 = d10 - vchar_values[n6 + 2]
                n14 = _over16_first_table_index(
                    n9,
                    n10,
                    candidate -> !(-gslope_values[n6 + 5] > cchar_values[candidate]),
                )
                n14_extension = n14 + n7 + 1
                d10 = vchar_values[n14_extension] * gslope_values[n6 + 5] -
                    cchar_values[n14_extension]
                d13 = d10 - vchar_values[n6 + 5]
                if vchar_values[n6 + 2] > vchar_values[n6 + 5]
                    d14 = 0.0
                    abs(d13) > zero && (d14 = d11 * d6 / d13)
                    d9 -= d14
                    vchar_values[n6 + 5] = d9
                    gslope_values[n6 + 5] = gslope_values[n6 + 4]
                    d13 = d14
                else
                    vchar_values[n6 + 5] = vchar_values[n6 + 4]
                    gslope_values[n6 + 5] = gslope_values[n6 + 4]
                    d13 = d6
                end
                vchar_values[n6 + 4] = vchar_values[n6 + 2]
                gslope_values[n6 + 4] = vchar_values[n6 + 3]
                d13, d14 = rebuild_minor_loop!(d6, d9, d11, d13)
                terminal_branch = cchar_values[n6 + 1] == 1.0 ? :upper : :downer
            end
        end
    else
        if abs(flux) >= vchar_values[n10]
            vchar_values[n6] = 0.0
            vchar_values[n6 + 1] = 0.0
        end
        d7 = vchar_values[n6] * flux + vchar_values[n6 + 1]
        if d7 < 0.0
            d7 = 0.0
            vchar_values[n6] = 0.0
            vchar_values[n6 + 1] = 0.0
        end
        if cchar_values[n6 + 1] == 1.0
            if cchar_values[n6 + 5] == 1.0
                index = _over16_first_table_index(
                    n9,
                    n10,
                    candidate -> !(-flux > vchar_values[candidate]),
                )
                set_trace_index!(index)
            else
                d8 = flux - d7
                index = _over16_first_table_index(
                    n9,
                    n10,
                    candidate -> !(-branch_current > cchar_values[candidate]),
                )
                set_trace_index!(index)
                n13 = index + n7 + 1
                d9 = vchar_values[n13] * branch_current - cchar_values[n13]
                d10 = d9 - d8
                if d7 <= d10 + zero
                    index = _over16_first_table_index(
                        n9,
                        n10,
                        candidate -> !(d8 > vchar_values[candidate]),
                    )
                    set_trace_index!(index)
                else
                    cchar_values[n6 + 5] = 1.0
                    index = _over16_first_table_index(
                        n9,
                        n10,
                        candidate -> !(-flux > vchar_values[candidate]),
                    )
                    set_trace_index!(index)
                end
            end
            terminal_branch = :upper
        else
            if cchar_values[n6 + 5] == 1.0
                index = _over16_first_table_index(
                    n9,
                    n10,
                    candidate -> !(flux > vchar_values[candidate]),
                )
                set_trace_index!(index)
            else
                d8 = flux + d7
                index = _over16_first_table_index(
                    n9,
                    n10,
                    candidate -> !(branch_current > cchar_values[candidate]),
                )
                set_trace_index!(index)
                n13 = index + n7 + 1
                d9 = vchar_values[n13] * branch_current + cchar_values[n13]
                d10 = d8 - d9
                if d7 <= d10 + zero
                    index = _over16_first_table_index(
                        n9,
                        n10,
                        candidate -> !(-d8 > vchar_values[candidate]),
                    )
                    set_trace_index!(index)
                else
                    cchar_values[n6 + 5] = 1.0
                    index = _over16_first_table_index(
                        n9,
                        n10,
                        candidate -> !(flux > vchar_values[candidate]),
                    )
                    set_trace_index!(index)
                end
            end
            terminal_branch = :downer
        end
        d13, d14 = terminal_branch == :upper ? upper_final!() : downer_final!()
    end

    d13 != 0.0 || throw(ArgumentError("hysteretic inductor trajectory denominator must not be zero"))
    if gslope_values[n12] == cchar_values[n6 + 1] &&
       gslope_values[n6 + 3] == cchar_values[n12]
        admittance_delta = 0.0
    else
        d15 = delta / d13
        admittance_delta = d15 - gslope_values[n6 + 1]
        gslope_values[n6 + 1] = d15
    end
    d16 = (flux - d14 + delta * voltage) / d13
    source_current_delta = d16 - gslope_values[n6]
    gslope_values[n6] = d16
    gslope_values[n12] = cchar_values[n6 + 1]
    gslope_values[n6 + 3] = cchar_values[n12]
    vchar_values[n6 + 3] = cchar_values[n6 + 3]
    vchar_values[n12] = flux
    stored_flux_after_step = flux + voltage * delta

    return (
        source = :hysteretic_inductor_current_update,
        outcome = :state_mutation,
        fortran_files = (:OVER16_FOR,),
        fortran_routines = (:SUBTS3,),
        fortran_labels = HYSTERETIC_INDUCTOR_REFERENCE_LABELS,
        cchar = cchar_values,
        vchar = vchar_values,
        gslope = gslope_values,
        current = branch_current,
        source_current_delta = source_current_delta,
        companion_current = gslope_values[n6],
        admittance_delta = admittance_delta,
        stored_flux = stored_flux_after_step,
        trace_index = _over16_integer_table_pointer(cchar_values[n12], "hysteretic inductor trace index"),
        reversal_count = reversal_count,
        trajectory_limit_count = trajectory_limit_count,
        terminal_branch = terminal_branch,
        mutation_order = (
            :flux_predictor,
            :current_prediction,
            :hysteresis_trace_selection,
            :companion_admittance_update,
            :source_current_delta,
            :stored_flux_update,
        ),
        deferred_calls = admittance_delta == 0.0 ?
            [:bulk_last14_oracle] :
            [:sparse_admittance_restamp, :bulk_last14_oracle],
        replacement_ready = false,
    )
end

function over16_nonlinear_source_column_assembly(
    ntot::Int,
    ncomp::Int,
    kode::AbstractVector{Int},
    source_begin_indices::AbstractVector{Int},
    source_next_indices::AbstractVector{Int},
    source_from_nodes::AbstractVector{Int},
    source_to_nodes::AbstractVector{Int},
    source_activity_flags::AbstractVector{Int};
    nonlinear_types::AbstractVector{Int}=Int[],
    nonlinear_admittance_nodes::AbstractVector{Int}=Int[],
    initial_vzero::AbstractVector{<:Real}=Float64[],
    initial_ilast::AbstractVector{Int}=Int[],
    num99::Int=0,
    fltinf::Real=Inf,
)
    _check_over16_nonlinear_dimensions(ntot, ncomp)
    _check_over16_kode(kode, ntot)
    link_count = length(source_next_indices)
    _over16_check_source_column_vector("source_from_nodes", source_from_nodes, link_count)
    _over16_check_source_column_vector("source_to_nodes", source_to_nodes, link_count)
    _over16_check_source_column_vector("source_activity_flags", source_activity_flags, ntot)
    num99 >= 0 || throw(ArgumentError("num99 must be nonnegative"))
    infinity = Float64(fltinf)
    infinity > 0.0 || throw(ArgumentError("fltinf must be positive"))

    znonl = zeros(Float64, ntot * ncomp)
    nonlinear_count = length(nonlinear_types)
    _over16_check_nonlinear_source_metadata(
        nonlinear_types,
        nonlinear_admittance_nodes,
        initial_vzero,
        initial_ilast,
    )
    vzero = isempty(initial_vzero) ? zeros(Float64, nonlinear_count) : Float64.(initial_vzero)
    ilast = isempty(initial_ilast) ? zeros(Int, nonlinear_count) : collect(initial_ilast)

    nonlinear_initialized_count = 0
    if nonlinear_count != num99
        for index in eachindex(nonlinear_types)
            nonlinear_type = nonlinear_types[index]
            if nonlinear_type < 0 || nonlinear_type > 920
                continue
            end
            node = nonlinear_admittance_nodes[index]
            vzero[index] = node < 0 ? 0.0 : -infinity
            ilast[index] = abs(node)
            nonlinear_initialized_count += 1
        end
    end

    source_link_count = 0
    source_column_count = 0
    active_terminal_count = 0
    kode_redirect_count = 0
    reference_node_skip_count = 0
    inactive_terminal_skip_count = 0
    for head in source_begin_indices
        head <= 0 && continue
        1 <= head <= link_count ||
            throw(ArgumentError("source_begin_indices entries must be within source link count"))
        offset = 0
        link = head
        seen_links = Set{Int}()
        while true
            link in seen_links &&
                throw(ArgumentError("source links must close only at their head"))
            push!(seen_links, link)
            source_link_count += 1
            component = div(offset, ntot) + 1
            component <= ncomp ||
                throw(ArgumentError("source link count exceeds ncomp"))
            for (terminal_node, sign) in (
                (source_from_nodes[link], -1.0),
                (source_to_nodes[link], 1.0),
            )
                1 <= terminal_node <= ntot ||
                    throw(ArgumentError("source terminal nodes must be within ntot"))
                if source_activity_flags[terminal_node] == 0
                    inactive_terminal_skip_count += 1
                    continue
                end
                representative, redirects =
                    _over16_source_column_representative(terminal_node, kode)
                kode_redirect_count += redirects
                if kode[representative] == 1
                    reference_node_skip_count += 1
                    continue
                end
                znonl[representative + offset] = sign
                active_terminal_count += 1
                source_column_count += 1
            end
            next_link = source_next_indices[link]
            1 <= next_link <= link_count ||
                throw(ArgumentError("source_next_indices entries must be within source link count"))
            next_link == head && break
            link = next_link
            offset += ntot
        end
    end

    return (
        source = :over16_nonlinear_source_column_assembly,
        outcome = :state_mutation,
        fortran_files = (:OVER16_FOR,),
        fortran_routines = (:SUBTS1,),
        fortran_labels = OVER16_NONLINEAR_SOURCE_COLUMN_LABELS,
        ntot = ntot,
        ncomp = ncomp,
        num99 = num99,
        znonl = znonl,
        vzero = vzero,
        ilast = ilast,
        nonlinear_initialized_count = nonlinear_initialized_count,
        source_link_count = source_link_count,
        active_terminal_count = active_terminal_count,
        source_column_count = source_column_count,
        kode_redirect_count = kode_redirect_count,
        reference_node_skip_count = reference_node_skip_count,
        inactive_terminal_skip_count = inactive_terminal_skip_count,
        nonlinear_source_columns_built = true,
        nonlinear_source_column_state_mutated = false,
        znonl_mutated = false,
        vzero_mutated = false,
        ilast_mutated = false,
        mutation_order = (:znonl_zero_fill, :nonlinear_vzero_ilast_init, :source_column_write),
        deferred_calls = [:nonlinear_table_card_mutation, :bulk_last14_oracle],
        tacs_executed = false,
        solvum_executed = false,
        replacement_ready = false,
    )
end

function over16_nonlinear_source_column_assembly!(
    state::OVER16NonlinearInverseColumnState,
    kode::AbstractVector{Int},
    source_begin_indices::AbstractVector{Int},
    source_next_indices::AbstractVector{Int},
    source_from_nodes::AbstractVector{Int},
    source_to_nodes::AbstractVector{Int},
    source_activity_flags::AbstractVector{Int};
    kwargs...,
)
    znonl_before = copy(state.znonl)
    vzero_before = copy(state.vzero)
    ilast_before = copy(state.ilast)
    count_before = state.source_column_update_count
    preview = over16_nonlinear_source_column_assembly(
        state.ntot,
        state.ncomp,
        kode,
        source_begin_indices,
        source_next_indices,
        source_from_nodes,
        source_to_nodes,
        source_activity_flags;
        initial_vzero = state.vzero,
        initial_ilast = state.ilast,
        kwargs...,
    )

    resize!(state.znonl, length(preview.znonl))
    state.znonl .= preview.znonl
    resize!(state.vzero, length(preview.vzero))
    state.vzero .= preview.vzero
    resize!(state.ilast, length(preview.ilast))
    state.ilast .= preview.ilast
    state.source_column_update_count += 1

    znonl_mutated = state.znonl != znonl_before
    vzero_mutated = state.vzero != vzero_before
    ilast_mutated = state.ilast != ilast_before
    count_mutated = state.source_column_update_count != count_before
    state_mutated = znonl_mutated || vzero_mutated || ilast_mutated || count_mutated
    return merge(
        preview,
        (
            znonl = copy(state.znonl),
            vzero = copy(state.vzero),
            ilast = copy(state.ilast),
            source_column_update_count = state.source_column_update_count,
            znonl_mutated = znonl_mutated,
            vzero_mutated = vzero_mutated,
            ilast_mutated = ilast_mutated,
            source_column_update_count_mutated = count_mutated,
            nonlinear_source_column_state_mutated = state_mutated,
        ),
    )
end

function over16_nonlinear_inverse_column_solution(
    znonl::AbstractVector{<:Real},
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int},
    kode::AbstractVector{Int};
    ntot::Int,
    ncomp::Int,
    partition_boundary::Int,
    iupper::Int=length(km),
    nonlinear_types::AbstractVector{Int}=Int[],
    nonlinear_from_nodes::AbstractVector{Int}=Int[],
    nonlinear_to_nodes::AbstractVector{Int}=Int[],
    nonlinear_source_flags::AbstractVector{Int}=Int[],
    initial_anonl::AbstractVector{<:Real}=Float64[],
    delta2::Real=1.0,
)
    _check_over16_nonlinear_dimensions(ntot, ncomp)
    _check_over16_znonl_layout(znonl, ntot, ncomp)
    _check_over16_nonlinear_sparse_workspace(km, ykm, kk, ntot, partition_boundary, iupper)
    _check_over16_kode(kode, ntot)
    delta = Float64(delta2)
    isfinite(delta) || throw(ArgumentError("delta2 must be finite"))

    work = Float64.(znonl)
    scratch = zeros(Float64, ncomp)
    yvalues = Float64.(ykm)
    anonl, flags = _over16_nonlinear_difference_inputs(
        nonlinear_types,
        nonlinear_from_nodes,
        nonlinear_to_nodes,
        nonlinear_source_flags,
        initial_anonl,
        ntot,
    )

    ii = 1
    forward_row_count = 0
    forward_elimination_count = 0
    forward_partition_skip_count = 0
    while ii <= iupper
        l = abs(km[ii])
        for component in 1:ncomp
            index = _over16_nonlinear_column_index(l, component, ntot)
            scratch[component] = work[index]
            work[index] *= yvalues[ii]
        end
        forward_row_count += 1
        row_end = abs(kk[l])
        while true
            ii += 1
            ii > row_end && break
            k = km[ii]
            if k <= partition_boundary
                target = k
                for component in 1:ncomp
                    work[target] -= scratch[component] * yvalues[ii]
                    target += ntot
                end
                forward_elimination_count += 1
            else
                ii = row_end + 1
                forward_partition_skip_count += 1
                break
            end
        end
    end

    backward_accumulation_count = 0
    kode_group_copy_count = 0
    while ii != 1
        fill!(scratch, 0.0)
        while true
            ii -= 1
            ii >= 1 || throw(ArgumentError("sparse workspace back substitution underflowed"))
            k = km[ii]
            if k >= 0
                target = k
                for component in 1:ncomp
                    scratch[component] -= work[target] * yvalues[ii]
                    target += ntot
                end
                backward_accumulation_count += 1
                continue
            end

            l = abs(k)
            for component in 1:ncomp
                work[_over16_nonlinear_column_index(l, component, ntot)] += scratch[component]
            end
            kode_group_copy_count += _over16_copy_kode_group!(work, kode, l, ntot, ncomp)
            break
        end
    end

    difference_count = 0
    for index in eachindex(nonlinear_types)
        nonlinear_type = nonlinear_types[index]
        if nonlinear_type < 0 || nonlinear_type > 920
            continue
        end
        from_node = nonlinear_from_nodes[index]
        to_node = abs(nonlinear_to_nodes[index])
        value = work[from_node] - work[to_node]
        if nonlinear_type != 94 && flags[index] <= 0
            value *= delta
        end
        anonl[index] = value
        difference_count += 1
    end

    return (
        source = :over16_nonlinear_inverse_column_solution,
        outcome = :state_mutation,
        fortran_files = (:OVER16_FOR,),
        fortran_routines = (:SUBTS1,),
        fortran_labels = OVER16_NONLINEAR_INVERSE_COLUMN_LABELS,
        ntot = ntot,
        ncomp = ncomp,
        partition_boundary = partition_boundary,
        iupper = iupper,
        znonl = work,
        anonl = anonl,
        voltbc = copy(scratch),
        forward_row_count = forward_row_count,
        forward_elimination_count = forward_elimination_count,
        forward_partition_skip_count = forward_partition_skip_count,
        backward_accumulation_count = backward_accumulation_count,
        kode_group_copy_count = kode_group_copy_count,
        nonlinear_difference_count = difference_count,
        nonlinear_inverse_columns_built = true,
        nonlinear_inverse_column_state_mutated = false,
        znonl_mutated = false,
        anonl_mutated = false,
        voltbc_mutated = false,
        mutation_order = (:forward_sparse_solve, :backward_sparse_solve, :kode_group_copy, :anonl_difference),
        deferred_calls = [:znonl_source_column_assembly, :nonlinear_table_card_mutation, :bulk_last14_oracle],
        tacs_executed = false,
        solvum_executed = false,
        replacement_ready = false,
    )
end

function over16_nonlinear_inverse_column_solution!(
    state::OVER16NonlinearInverseColumnState,
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int},
    kode::AbstractVector{Int};
    kwargs...,
)
    znonl_before = copy(state.znonl)
    anonl_before = copy(state.anonl)
    voltbc_before = copy(state.voltbc)
    update_count_before = state.update_count
    preview = over16_nonlinear_inverse_column_solution(
        state.znonl,
        km,
        ykm,
        kk,
        kode;
        ntot = state.ntot,
        ncomp = state.ncomp,
        initial_anonl = state.anonl,
        kwargs...,
    )

    resize!(state.znonl, length(preview.znonl))
    state.znonl .= preview.znonl
    resize!(state.anonl, length(preview.anonl))
    state.anonl .= preview.anonl
    resize!(state.voltbc, length(preview.voltbc))
    state.voltbc .= preview.voltbc
    state.update_count += 1

    znonl_mutated = state.znonl != znonl_before
    anonl_mutated = state.anonl != anonl_before
    voltbc_mutated = state.voltbc != voltbc_before
    update_count_mutated = state.update_count != update_count_before
    state_mutated = znonl_mutated || anonl_mutated || voltbc_mutated || update_count_mutated
    return merge(
        preview,
        (
            znonl = copy(state.znonl),
            anonl = copy(state.anonl),
            voltbc = copy(state.voltbc),
            update_count = state.update_count,
            znonl_mutated = znonl_mutated,
            anonl_mutated = anonl_mutated,
            voltbc_mutated = voltbc_mutated,
            update_count_mutated = update_count_mutated,
            nonlinear_inverse_column_state_mutated = state_mutated,
        ),
    )
end

function over16_simultaneous_zno_solution(
    znonl::AbstractVector{<:Real},
    voltages::AbstractVector{<:Real},
    rhs::AbstractVector{<:Real},
    cchar::AbstractVector{<:Real},
    vchar::AbstractVector{<:Real},
    subnetwork_next_indices::AbstractVector{Int},
    subnetwork_from_nodes::AbstractVector{Int},
    subnetwork_to_nodes::AbstractVector{Int},
    subnetwork_nonlinear_indices::AbstractVector{Int},
    subnetwork_element_types::AbstractVector{Int};
    ntot::Int,
    ncomp::Int,
    subsystem_begin_index::Int,
    nonlinear_admittance_start_indices::AbstractVector{Int},
    nonlinear_table_end_indices::AbstractVector{Int},
    initial_anonl::AbstractVector{<:Real},
    initial_vzero::AbstractVector{<:Real},
    initial_ilast::AbstractVector{Int},
    initial_curr::AbstractVector{<:Real}=Float64[],
    initial_cursub::AbstractVector{<:Real}=Float64[],
    vnonl::AbstractVector{<:Real}=Float64[],
    gslope::AbstractVector{<:Real}=Float64[],
    gap_status_values::AbstractVector{<:Real}=Float64[],
    t::Real=0.0,
    deltat::Real=0.0,
    epszno::Real=1.0e-10,
    znolim::Tuple{<:Real,<:Real}=(0.2, 1.5),
    max_iterations::Int=20,
    fltinf::Real=Inf,
    epsiln::Real=1.0e-12,
)
    _check_over16_nonlinear_dimensions(ntot, ncomp)
    _check_over16_znonl_layout(znonl, ntot, ncomp)
    _check_over16_simultaneous_zno_inputs(
        voltages,
        rhs,
        cchar,
        vchar,
        subnetwork_next_indices,
        subnetwork_from_nodes,
        subnetwork_to_nodes,
        subnetwork_nonlinear_indices,
        subnetwork_element_types,
        subsystem_begin_index,
        nonlinear_admittance_start_indices,
        nonlinear_table_end_indices,
        initial_anonl,
        initial_vzero,
        initial_ilast,
        initial_curr,
        initial_cursub,
        vnonl,
        gslope,
        gap_status_values,
        ntot,
    )
    time = Float64(t)
    isfinite(time) || throw(ArgumentError("t must be finite"))
    dt = Float64(deltat)
    isfinite(dt) || throw(ArgumentError("deltat must be finite"))
    tolerance = Float64(epszno)
    isfinite(tolerance) && tolerance > 0.0 ||
        throw(ArgumentError("epszno must be finite and positive"))
    voltage_limit = (Float64(znolim[1]), Float64(znolim[2]))
    all(value -> isfinite(value) && value > 0.0, voltage_limit) ||
        throw(ArgumentError("znolim entries must be finite and positive"))
    max_iterations >= 0 || throw(ArgumentError("max_iterations must be nonnegative"))
    infinity = Float64(fltinf)
    infinity > 0.0 || throw(ArgumentError("fltinf must be positive"))
    epsilon = Float64(epsiln)
    isfinite(epsilon) && epsilon > 0.0 ||
        throw(ArgumentError("epsiln must be finite and positive"))

    records = _over16_simultaneous_zno_record_chain(
        subnetwork_next_indices,
        subsystem_begin_index,
    )
    element_count = length(records)
    zthevenin = _over16_simultaneous_zno_thevenin_matrix(
        znonl,
        subnetwork_from_nodes,
        subnetwork_to_nodes,
        records,
        ntot,
    )
    singular, dependent = _over16_simultaneous_zno_column_status(zthevenin)
    independent_columns = [index for index in 1:element_count if !singular[index] && dependent[index] == 0]
    independent_count = length(independent_columns)
    dependent_count = count(!=(0), dependent)
    singular_count = count(identity, singular)
    piecewise_resistance_count = count(
        record -> subnetwork_element_types[record] == 2,
        records,
    )
    time_varying_resistance_count = count(
        record -> subnetwork_element_types[record] == 3,
        records,
    )

    curr = isempty(initial_curr) ? zeros(Float64, length(initial_anonl)) : Float64.(initial_curr)
    cursub = isempty(initial_cursub) ?
        zeros(Float64, maximum(records; init = 0) ÷ 5 + 1) :
        Float64.(initial_cursub)
    rhs_values = Float64.(rhs)
    vzero = Float64.(initial_vzero)
    ilast = collect(initial_ilast)
    anonl = Float64.(initial_anonl)
    cchar_values = Float64.(cchar)
    vchar_values = Float64.(vchar)
    gslope_values = isempty(gslope) ? zeros(Float64, length(cchar_values)) : Float64.(gslope)
    nonlinear_gap_limits =
        isempty(vnonl) ? fill(infinity, length(anonl)) : Float64.(vnonl)
    nonlinear_gap_status =
        isempty(gap_status_values) ? zeros(Float64, length(anonl)) : Float64.(gap_status_values)
    voltage_values = Float64.(voltages)

    element_voltage = zeros(Float64, element_count)
    thevenin_voltage = zeros(Float64, element_count)
    previous_current = copy(curr)
    for (position, record) in enumerate(records)
        nonlinear_index = subnetwork_nonlinear_indices[record]
        from_node = subnetwork_from_nodes[record]
        to_node = subnetwork_to_nodes[record]
        thevenin_voltage[position] = voltage_values[from_node] - voltage_values[to_node]
        element_voltage[position] = vzero[nonlinear_index]
    end

    independent_position = Dict(column => position for (position, column) in enumerate(independent_columns))
    base_inverse = independent_count == 0 ?
        zeros(Float64, 0, 0) :
        inv(zthevenin[independent_columns, independent_columns])
    max_voltage_correction = Inf
    previous_voltage_correction = Inf
    iteration_count = 0
    zno_segment_search_count = 0
    correction_scale_count = 0
    voltage_limit_count = 0
    converged = independent_count == 0

    while independent_count > 0
        iteration_count += 1
        jacobian = copy(base_inverse)
        residual = zeros(Float64, independent_count)
        for (element_position, record) in enumerate(records)
            singular[element_position] && continue
            representative = dependent[element_position] == 0 ?
                element_position : dependent[element_position]
            row = independent_position[representative]
            nonlinear_index = subnetwork_nonlinear_indices[record]
            current, derivative, updated_last, search_steps =
                _over16_simultaneous_nonlinear_current_and_derivative(
                subnetwork_element_types[record],
                element_voltage[element_position],
                anonl[nonlinear_index],
                ilast[nonlinear_index],
                nonlinear_admittance_start_indices[nonlinear_index],
                nonlinear_table_end_indices[nonlinear_index],
                cchar_values,
                gslope_values,
                vchar_values,
                nonlinear_gap_limits[nonlinear_index],
                time,
                epsilon,
            )
            curr[nonlinear_index] = current
            ilast[nonlinear_index] = updated_last
            zno_segment_search_count += search_steps
            jacobian[row, row] -= derivative
            residual[row] += current
        end
        independent_voltage_delta =
            element_voltage[independent_columns] .- thevenin_voltage[independent_columns]
        # ZINCOX stores the inverse column-major but walks it sequentially in this residual loop.
        residual .-= transpose(base_inverse) * independent_voltage_delta
        if max_voltage_correction <= tolerance
            converged = true
            break
        end

        correction = jacobian \ residual
        max_voltage_correction = 0.0
        for (local_index, element_position) in enumerate(independent_columns)
            record = records[element_position]
            nonlinear_index = subnetwork_nonlinear_indices[record]
            scaled = abs(correction[local_index] / anonl[nonlinear_index])
            max_voltage_correction = max(max_voltage_correction, scaled)
        end
        if max_voltage_correction > voltage_limit[1]
            scale = voltage_limit[1] / max_voltage_correction
            correction .*= scale
            correction_scale_count += 1
        end
        if iteration_count >= 4 && isfinite(previous_voltage_correction) &&
                previous_voltage_correction != 0.0
            relative_change =
                abs(max_voltage_correction - previous_voltage_correction) /
                abs(previous_voltage_correction)
            if relative_change <= 10.0 * epsilon
                correction .*= 0.1
            end
        end
        for (local_index, element_position) in enumerate(independent_columns)
            record = records[element_position]
            nonlinear_index = subnetwork_nonlinear_indices[record]
            element_voltage[element_position] += correction[local_index]
            if subnetwork_element_types[record] == 1
                limit = voltage_limit[2] * anonl[nonlinear_index]
                if abs(element_voltage[element_position]) >= limit
                    element_voltage[element_position] = copysign(limit, element_voltage[element_position])
                    voltage_limit_count += 1
                end
            end
        end
        for element_position in 1:element_count
            source = dependent[element_position]
            source != 0 && (element_voltage[element_position] = element_voltage[source])
        end
        previous_voltage_correction = max_voltage_correction
        iteration_count > max_iterations && break
    end

    rhs_update_count = 0
    cursub_update_count = 0
    for (element_position, record) in enumerate(records)
        nonlinear_index = subnetwork_nonlinear_indices[record]
        if singular[element_position]
            current, _, updated_last, search_steps =
                _over16_simultaneous_nonlinear_current_and_derivative(
                subnetwork_element_types[record],
                thevenin_voltage[element_position],
                anonl[nonlinear_index],
                ilast[nonlinear_index],
                nonlinear_admittance_start_indices[nonlinear_index],
                nonlinear_table_end_indices[nonlinear_index],
                cchar_values,
                gslope_values,
                vchar_values,
                nonlinear_gap_limits[nonlinear_index],
                time,
                epsilon,
            )
            curr[nonlinear_index] = current
            ilast[nonlinear_index] = updated_last
            element_voltage[element_position] = thevenin_voltage[element_position]
            zno_segment_search_count += search_steps
        end

        gap_limit = nonlinear_gap_limits[nonlinear_index]
        if gap_limit != infinity
            ils = ilast[nonlinear_index]
            element_type = subnetwork_element_types[record]
            if ils >= 0
                if abs(element_voltage[element_position]) > gap_limit
                    ils = -ils
                end
                if element_type == 3 && ils < 0
                    nonlinear_gap_limits[nonlinear_index] = -(time + dt)
                end
            elseif element_type != 3
                if curr[nonlinear_index] * previous_current[nonlinear_index] < 0.0
                    ils = -ils
                end
                if element_type != 1 && ils >= 0
                    gap = nonlinear_gap_status[nonlinear_index]
                    if gap < 0.0
                        ils = -ils
                    elseif gap > 0.0
                        nonlinear_gap_limits[nonlinear_index] = infinity
                    end
                end
            end
            ilast[nonlinear_index] = ils
        end
        vzero[nonlinear_index] = element_voltage[element_position]
        cursub_index = div(record, 5) + 1
        cursub_index <= length(cursub) ||
            throw(ArgumentError("initial_cursub length must cover subnetwork record heads"))
        cursub[cursub_index] = curr[nonlinear_index]
        from_node = subnetwork_from_nodes[record]
        to_node = subnetwork_to_nodes[record]
        rhs_values[from_node] -= curr[nonlinear_index]
        rhs_values[to_node] += curr[nonlinear_index]
        cursub_update_count += 1
        rhs_update_count += 2
    end

    return (
        source = :over16_simultaneous_zno_solution,
        outcome = :state_mutation,
        fortran_files = (:OVER16_FOR,),
        fortran_routines = (:ZINCOX,),
        fortran_labels = OVER16_SIMULTANEOUS_ZNO_LABELS,
        records = Tuple(records),
        independent_columns = Tuple(independent_columns),
        dependent_columns = Tuple(dependent),
        singular_columns = Tuple(findall(identity, singular)),
        zthevenin = zthevenin,
        thevenin_voltage = thevenin_voltage,
        element_voltage = element_voltage,
        rhs = rhs_values,
        curr = curr,
        cursub = cursub,
        vzero = vzero,
        ilast = ilast,
        vnonl = nonlinear_gap_limits,
        element_count = element_count,
        independent_count = independent_count,
        dependent_count = dependent_count,
        singular_count = singular_count,
        piecewise_resistance_count = piecewise_resistance_count,
        time_varying_resistance_count = time_varying_resistance_count,
        iteration_count = iteration_count,
        converged = converged,
        max_voltage_correction = isfinite(max_voltage_correction) ? max_voltage_correction : 0.0,
        correction_scale_count = correction_scale_count,
        voltage_limit_count = voltage_limit_count,
        zno_segment_search_count = zno_segment_search_count,
        rhs_update_count = rhs_update_count,
        cursub_update_count = cursub_update_count,
        current_sign_change_count =
            count(index -> curr[index] * previous_current[index] < 0.0, eachindex(curr)),
        simultaneous_zno_solution_applied = true,
        nonlinear_current_state_mutated = false,
        rhs_mutated = false,
        curr_mutated = false,
        cursub_mutated = false,
        vzero_mutated = false,
        ilast_mutated = false,
        mutation_order = (
            :extract_zthevenin,
            :classify_singular_dependent_columns,
            :newton_current_balance,
            :voltage_limit,
            :final_current_gap_state,
            :rhs_current_injection,
        ),
        deferred_calls = [:full_last14_card_execution],
        tacs_executed = false,
        solvum_executed = false,
        replacement_ready = false,
    )
end

function over16_simultaneous_zno_solution!(
    state::OVER16NonlinearInverseColumnState,
    voltages::AbstractVector{<:Real},
    rhs::AbstractVector{<:Real},
    cchar::AbstractVector{<:Real},
    vchar::AbstractVector{<:Real},
    subnetwork_next_indices::AbstractVector{Int},
    subnetwork_from_nodes::AbstractVector{Int},
    subnetwork_to_nodes::AbstractVector{Int},
    subnetwork_nonlinear_indices::AbstractVector{Int},
    subnetwork_element_types::AbstractVector{Int};
    kwargs...,
)
    rhs_before = Float64.(rhs)
    curr_before = copy(state.curr)
    cursub_before = copy(state.cursub)
    vzero_before = copy(state.vzero)
    vnonl_before = copy(state.vnonl)
    ilast_before = copy(state.ilast)
    count_before = state.current_update_count
    owner_kwargs = haskey(kwargs, :vnonl) ? (; kwargs...) : (; vnonl = state.vnonl, kwargs...)
    preview = over16_simultaneous_zno_solution(
        state.znonl,
        voltages,
        rhs,
        cchar,
        isempty(state.vchar) ? vchar : state.vchar,
        subnetwork_next_indices,
        subnetwork_from_nodes,
        subnetwork_to_nodes,
        subnetwork_nonlinear_indices,
        subnetwork_element_types;
        ntot = state.ntot,
        ncomp = state.ncomp,
        initial_anonl = state.anonl,
        initial_vzero = state.vzero,
        initial_ilast = state.ilast,
        initial_curr = state.curr,
        initial_cursub = state.cursub,
        owner_kwargs...,
    )

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
    state.current_update_count += 1

    rhs_mutated = preview.rhs != rhs_before
    curr_mutated = state.curr != curr_before
    cursub_mutated = state.cursub != cursub_before
    vzero_mutated = state.vzero != vzero_before
    vnonl_mutated = state.vnonl != vnonl_before
    ilast_mutated = state.ilast != ilast_before
    count_mutated = state.current_update_count != count_before
    state_mutated =
        curr_mutated ||
        cursub_mutated ||
        vzero_mutated ||
        vnonl_mutated ||
        ilast_mutated ||
        count_mutated
    return merge(
        preview,
        (
            curr = copy(state.curr),
            cursub = copy(state.cursub),
            vzero = copy(state.vzero),
            vnonl = copy(state.vnonl),
            ilast = copy(state.ilast),
            current_update_count = state.current_update_count,
            rhs_mutated = rhs_mutated,
            curr_mutated = curr_mutated,
            cursub_mutated = cursub_mutated,
            vzero_mutated = vzero_mutated,
            vnonl_mutated = vnonl_mutated,
            ilast_mutated = ilast_mutated,
            current_update_count_mutated = count_mutated,
            nonlinear_current_state_mutated = state_mutated,
        ),
    )
end

function _over16_node_voltage(voltages::AbstractVector{Float64}, node::Int)
    node == 0 && return 0.0
    return voltages[node]
end

function _over16_add_rhs_delta!(rhs::Vector{Float64}, node::Int, delta::Float64)
    node == 0 && return 0
    rhs[node] += delta
    return 1
end

function over16_nonlinear_current_compensation_update(
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
    initial_anonl::AbstractVector{<:Real},
    initial_vzero::AbstractVector{<:Real},
    initial_ilast::AbstractVector{Int},
    initial_curr::AbstractVector{<:Real}=Float64[],
    initial_cursub::AbstractVector{<:Real}=Float64[],
    vnonl::AbstractVector{<:Real}=Float64[],
    gslope::AbstractVector{<:Real}=Float64[],
    delta2::Real=1.0,
    deltat::Real=2.0 * Float64(delta2),
    t::Real=0.0,
    epsiln::Real=1.0e-12,
    flzero::Real=0.0,
    fltinf::Real=1.0e99,
    num99::Int=0,
    minimum_on_time_values::AbstractVector{<:Real}=Float64[],
    timed_resistance_arm_time_values::AbstractVector{<:Real}=Float64[],
    rearm_time_state_indices::AbstractVector{Int}=Int[],
    single_flash_flags::AbstractVector{Bool}=Bool[],
    simultaneous_zno_config::Union{Nothing,NamedTuple}=nothing,
    nonlinear_sparse_config::Union{Nothing,NamedTuple}=nothing,
    saturated_transformer_sparse_config::Union{Nothing,NamedTuple}=nothing,
    network_current_response_columns::AbstractVector{<:Real}=Float64[],
    network_response_node_count::Int=0,
    network_response_component_count::Int=0,
    complete_nonlinear_source_loop::Bool=false,
)
    count = length(nonlinear_types)
    _over16_check_nonlinear_current_inputs(
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
        cchar,
        vchar,
        initial_anonl,
        initial_vzero,
        initial_ilast,
        initial_curr,
        initial_cursub,
        vnonl,
    )
    num99 >= 0 || throw(ArgumentError("num99 must be nonnegative"))
    delta = Float64(delta2)
    isfinite(delta) || throw(ArgumentError("delta2 must be finite"))
    dt = Float64(deltat)
    isfinite(dt) && dt > 0.0 || throw(ArgumentError("deltat must be finite and positive"))
    epsilon = Float64(epsiln)
    isfinite(epsilon) && epsilon > 0.0 || throw(ArgumentError("epsiln must be finite and positive"))
    flux_zero = Float64(flzero)
    isfinite(flux_zero) && flux_zero >= 0.0 ||
        throw(ArgumentError("flzero must be finite and nonnegative"))
    time_value = Float64(t)
    isfinite(time_value) || throw(ArgumentError("t must be finite"))
    infinity_value = Float64(fltinf)
    infinity_value > 0.0 || throw(ArgumentError("fltinf must be positive"))
    minimum_on_times = isempty(minimum_on_time_values) ?
        zeros(Float64, count) : Float64.(minimum_on_time_values)
    timed_resistance_arm_times = isempty(timed_resistance_arm_time_values) ?
        zeros(Float64, count) : Float64.(timed_resistance_arm_time_values)
    rearm_indices = isempty(rearm_time_state_indices) ?
        zeros(Int, count) : collect(rearm_time_state_indices)
    one_shot_flags = isempty(single_flash_flags) ?
        falses(count) : collect(single_flash_flags)
    length(minimum_on_times) == count &&
        length(timed_resistance_arm_times) == count &&
        length(rearm_indices) == count &&
        length(one_shot_flags) == count ||
        throw(ArgumentError("timed and switching resistor state arrays must match nonlinear count"))

    rhs_values = Float64.(rhs)
    anonl = Float64.(initial_anonl)
    vzero = Float64.(initial_vzero)
    ilast = collect(initial_ilast)
    curr = isempty(initial_curr) ? zeros(Float64, count) : Float64.(initial_curr)
    cursub = isempty(initial_cursub) ?
        zeros(Float64, _over16_required_cursub_length(nonlinear_subsystem_indices, subsystem_begin_indices)) :
        Float64.(initial_cursub)
    nonlinear_voltages = isempty(vnonl) ? zeros(Float64, count) : Float64.(vnonl)
    voltage_values = Float64.(voltages)
    cchar_values = Float64.(cchar)
    vchar_values = Float64.(vchar)
    gslope_values = Float64.(gslope)
    active_sparse_config = nonlinear_sparse_config === nothing ?
        saturated_transformer_sparse_config :
        nonlinear_sparse_config

    processed_count = 0
    interpolation_count = 0
    increasing_count = 0
    decreasing_count = 0
    unchanged_count = 0
    negative_type_skip_count = 0
    missing_subsystem_skip_count = 0
    simultaneous_zno_deferred_count = 0
    type94_processed_count = 0
    type94_flashover_count = 0
    type94_cleared_count = 0
    type94_skip_count = 0
    saturated_transformer_update_count = 0
    saturated_transformer_segment_change_count = 0
    saturated_transformer_polarity_warning_count = 0
    saturated_transformer_polarity_reversal_count = 0
    saturated_transformer_finit_deltas = zeros(Float64, count)
    saturated_transformer_admittance_deltas = zeros(Float64, count)
    saturated_transformer_sparse_ykm =
        active_sparse_config === nothing ?
        Float64[] :
        Float64.(get(active_sparse_config, :ykm, Float64[]))
    saturated_transformer_sparse_from_nodes =
        active_sparse_config === nothing ?
        Int[] :
        Int.(get(
            active_sparse_config,
            :from_nodes,
            nonlinear_from_nodes,
        ))
    saturated_transformer_sparse_to_nodes =
        active_sparse_config === nothing ?
        Int[] :
        Int.(get(
            active_sparse_config,
            :to_nodes,
            nonlinear_to_nodes,
        ))
    if active_sparse_config !== nothing
        length(saturated_transformer_sparse_from_nodes) == count ||
            throw(ArgumentError("nonlinear sparse from_nodes length must match nonlinear count"))
        length(saturated_transformer_sparse_to_nodes) == count ||
            throw(ArgumentError("nonlinear sparse to_nodes length must match nonlinear count"))
    end
    saturated_transformer_sparse_update_count = 0
    saturated_transformer_sparse_positive_update_count = 0
    saturated_transformer_sparse_negative_update_count = 0
    saturated_transformer_sparse_partition_skip_count = 0
    saturated_transformer_sparse_retriangularization_request_count = 0
    hysteretic_inductor_update_count = 0
    hysteretic_inductor_rhs_update_count = 0
    hysteretic_inductor_sparse_update_count = 0
    hysteretic_inductor_reversal_count = 0
    hysteretic_inductor_trajectory_limit_count = 0
    hysteretic_inductor_admittance_deltas = zeros(Float64, count)
    hysteretic_inductor_source_current_deltas = zeros(Float64, count)
    switching_resistor_update_count = 0
    switching_resistor_activation_count = 0
    switching_resistor_deactivation_count = 0
    switching_resistor_segment_change_count = 0
    switching_resistor_polarity_reversal_count = 0
    switching_resistor_reported_currents = zeros(Float64, count)
    switching_resistor_accepted_currents = zeros(Float64, count)
    switching_resistor_admittance_deltas = zeros(Float64, count)
    switching_resistor_companion_current_deltas = zeros(Float64, count)
    timed_resistance_update_count = 0
    timed_resistance_activation_count = 0
    timed_resistance_segment_change_count = 0
    timed_resistance_accepted_currents = zeros(Float64, count)
    timed_resistance_admittance_deltas = zeros(Float64, count)
    piecewise_nonlinear_inductor_update_count = 0
    piecewise_nonlinear_inductor_segment_change_count = 0
    piecewise_nonlinear_inductor_iteration_count = 0
    piecewise_nonlinear_inductor_accepted_currents = zeros(Float64, count)
    piecewise_nonlinear_inductor_accepted_fluxes = zeros(Float64, count)
    piecewise_nonlinear_inductor_accepted_voltages = zeros(Float64, count)
    piecewise_nonlinear_inductor_network_response = zeros(Float64, 0, 0)
    rhs_update_count = 0
    cursub_update_count = 0
    table_limit_count = 0
    arrester_state_update_count = 0
    simultaneous_zno_solution_count = 0
    simultaneous_zno_element_count = 0
    simultaneous_piecewise_resistance_count = 0
    simultaneous_time_varying_resistance_count = 0
    simultaneous_zno_independent_count = 0
    simultaneous_zno_dependent_count = 0
    simultaneous_zno_singular_count = 0
    simultaneous_zno_iteration_count = 0
    simultaneous_zno_converged = true
    simultaneous_zno_correction_scale_count = 0
    simultaneous_zno_voltage_limit_count = 0
    simultaneous_zno_segment_search_count = 0
    simultaneous_zno_max_voltage_correction = 0.0
    simultaneous_zno_element_voltage = Float64[]
    simultaneous_zno_zthevenin = zeros(Float64, 0, 0)
    simultaneous_zno_dispatch_slot_indices = Int[]
    simultaneous_zno_dispatch_next_slots = Int[]
    simultaneous_zno_dispatch_from_nodes = Int[]
    simultaneous_zno_dispatch_to_nodes = Int[]
    simultaneous_zno_dispatch_nonlinear_indices = Int[]
    simultaneous_zno_dispatch_type_codes = Int[]
    simultaneous_zno_dispatch_nonlad = Int[]
    simultaneous_zno_dispatch_nonle = Int[]
    simultaneous_zno_dispatch_anonl = Float64[]
    simultaneous_zno_dispatch_initial_vzero = Float64[]
    simultaneous_zno_dispatch_initial_ilast = Int[]
    simultaneous_zno_dispatch_initial_curr = Float64[]
    simultaneous_zno_dispatch_initial_vnonl = Float64[]
    simultaneous_zno_dispatch_gap_status = Float64[]

    piecewise_nonlinear_inductor_indices =
        findall(==(PIECEWISE_NONLINEAR_INDUCTOR_TYPE), nonlinear_types)
    if !isempty(piecewise_nonlinear_inductor_indices)
        network_response_node_count > 0 || throw(ArgumentError(
            "true nonlinear inductors require solved network response columns",
        ))
        network_response_component_count >= maximum(piecewise_nonlinear_inductor_indices) ||
            throw(ArgumentError("network response does not cover every nonlinear inductor"))
        length(network_current_response_columns) ==
            network_response_node_count * network_response_component_count ||
            throw(ArgumentError("network response column dimensions are inconsistent"))
        owner_count = length(piecewise_nonlinear_inductor_indices)
        response = zeros(Float64, owner_count, owner_count)
        measured_voltage = zeros(Float64, owner_count)
        for (row_position, nonlinear_index) in
            enumerate(piecewise_nonlinear_inductor_indices)
            from_node = nonlinear_from_nodes[nonlinear_index]
            to_node = abs(nonlinear_to_nodes[nonlinear_index])
            measured_voltage[row_position] = voltage_values[from_node] -
                _over16_node_voltage(voltage_values, to_node)
            for (column_position, response_index) in
                enumerate(piecewise_nonlinear_inductor_indices)
                offset = (response_index - 1) * network_response_node_count
                from_response = from_node == 0 ? 0.0 :
                    Float64(network_current_response_columns[offset + from_node])
                to_response = to_node == 0 ? 0.0 :
                    Float64(network_current_response_columns[offset + to_node])
                response[row_position, column_position] = from_response - to_response
            end
        end
        piecewise_nonlinear_inductor_network_response = response
        update = piecewise_nonlinear_inductor_network_step(
            nonlinear_voltages[piecewise_nonlinear_inductor_indices],
            measured_voltage,
            ilast[piecewise_nonlinear_inductor_indices],
            nonlinear_admittance_nodes[piecewise_nonlinear_inductor_indices],
            abs.(nonlinear_table_end_indices[piecewise_nonlinear_inductor_indices]),
            cchar_values,
            vchar_values,
            response;
            half_timestep_s = delta,
        )
        for (position, nonlinear_index) in
            enumerate(piecewise_nonlinear_inductor_indices)
            curr[nonlinear_index] = update.accepted_current_a[position]
            vzero[nonlinear_index] = update.accepted_flux_wb[position]
            nonlinear_voltages[nonlinear_index] = update.predictor_flux_wb[position]
            ilast[nonlinear_index] = update.active_segments[position]
            anonl[nonlinear_index] = delta * response[position, position]
            piecewise_nonlinear_inductor_accepted_currents[nonlinear_index] =
                update.accepted_current_a[position]
            piecewise_nonlinear_inductor_accepted_fluxes[nonlinear_index] =
                update.accepted_flux_wb[position]
            piecewise_nonlinear_inductor_accepted_voltages[nonlinear_index] =
                update.accepted_branch_voltage_v[position]
        end
        piecewise_nonlinear_inductor_update_count = owner_count
        piecewise_nonlinear_inductor_segment_change_count = update.segment_change_count
        piecewise_nonlinear_inductor_iteration_count = update.iteration_count
    end

    if count != num99
        for index in eachindex(nonlinear_types)
            nonlinear_type = nonlinear_types[index]
            if nonlinear_type < 0 &&
               nonlinear_type != SATURATED_TRANSFORMER_NONLINEAR_TYPE &&
               nonlinear_type != HYSTERETIC_INDUCTOR_NONLINEAR_TYPE &&
               nonlinear_type != TRIGGERED_TIMED_RESISTANCE_TYPE &&
               nonlinear_type != SWITCHING_NONLINEAR_RESISTOR_TYPE
                negative_type_skip_count += 1
                continue
            end
            subsystem_index = nonlinear_subsystem_indices[index]
            head = subsystem_begin_indices[subsystem_index]
            if head <= 0
                missing_subsystem_skip_count += 1
                continue
            end
            1 <= head <= length(subsystem_simultaneous_flags) ||
                throw(ArgumentError("subsystem_begin_indices entries must address subsystem flags"))
            if subsystem_simultaneous_flags[head] != 0
                if subsystem_owner_rows[head] == index
                    if simultaneous_zno_config === nothing
                        simultaneous_zno_deferred_count += 1
                    else
                        subnetwork_next_indices = get(
                            simultaneous_zno_config,
                            :subnetwork_next_indices,
                            Int[],
                        )
                        subnetwork_from_nodes = get(
                            simultaneous_zno_config,
                            :subnetwork_from_nodes,
                            Int[],
                        )
                        subnetwork_to_nodes = get(
                            simultaneous_zno_config,
                            :subnetwork_to_nodes,
                            Int[],
                        )
                        subnetwork_nonlinear_indices = get(
                            simultaneous_zno_config,
                            :subnetwork_nonlinear_indices,
                            Int[],
                        )
                        subnetwork_element_types = get(
                            simultaneous_zno_config,
                            :subnetwork_element_types,
                            Int[],
                        )
                        gap_status_values = get(
                            simultaneous_zno_config,
                            :gap_status_values,
                            Float64[],
                        )
                        for slot_index in eachindex(subnetwork_nonlinear_indices)
                            nonlinear_index = subnetwork_nonlinear_indices[slot_index]
                            if !(
                                nonlinear_index != 0 ||
                                subnetwork_element_types[slot_index] != 0 ||
                                subnetwork_from_nodes[slot_index] != 0 ||
                                subnetwork_to_nodes[slot_index] != 0 ||
                                subnetwork_next_indices[slot_index] != 0
                            )
                                continue
                            end
                            1 <= nonlinear_index <= count || continue
                            push!(simultaneous_zno_dispatch_slot_indices, slot_index)
                            push!(
                                simultaneous_zno_dispatch_next_slots,
                                subnetwork_next_indices[slot_index],
                            )
                            push!(
                                simultaneous_zno_dispatch_from_nodes,
                                subnetwork_from_nodes[slot_index],
                            )
                            push!(
                                simultaneous_zno_dispatch_to_nodes,
                                subnetwork_to_nodes[slot_index],
                            )
                            push!(
                                simultaneous_zno_dispatch_nonlinear_indices,
                                nonlinear_index,
                            )
                            push!(
                                simultaneous_zno_dispatch_type_codes,
                                subnetwork_element_types[slot_index],
                            )
                            push!(
                                simultaneous_zno_dispatch_nonlad,
                                nonlinear_admittance_nodes[nonlinear_index],
                            )
                            push!(
                                simultaneous_zno_dispatch_nonle,
                                nonlinear_table_end_indices[nonlinear_index],
                            )
                            push!(simultaneous_zno_dispatch_anonl, anonl[nonlinear_index])
                            push!(
                                simultaneous_zno_dispatch_initial_vzero,
                                vzero[nonlinear_index],
                            )
                            push!(
                                simultaneous_zno_dispatch_initial_ilast,
                                ilast[nonlinear_index],
                            )
                            push!(
                                simultaneous_zno_dispatch_initial_curr,
                                curr[nonlinear_index],
                            )
                            push!(
                                simultaneous_zno_dispatch_initial_vnonl,
                                nonlinear_voltages[nonlinear_index],
                            )
                            push!(
                                simultaneous_zno_dispatch_gap_status,
                                nonlinear_index <= length(gap_status_values) ?
                                    gap_status_values[nonlinear_index] :
                                    0.0,
                            )
                        end
                        zno_result = over16_simultaneous_zno_solution(
                            get(simultaneous_zno_config, :znonl, Float64[]),
                            voltage_values,
                            rhs_values,
                            cchar_values,
                            vchar_values,
                            subnetwork_next_indices,
                            subnetwork_from_nodes,
                            subnetwork_to_nodes,
                            subnetwork_nonlinear_indices,
                            subnetwork_element_types;
                            ntot = length(voltage_values),
                            ncomp = get(simultaneous_zno_config, :ncomp, 1),
                            subsystem_begin_index = get(
                                simultaneous_zno_config,
                                :subsystem_begin_index,
                                head,
                            ),
                            nonlinear_admittance_start_indices = nonlinear_admittance_nodes,
                            nonlinear_table_end_indices = nonlinear_table_end_indices,
                            initial_anonl = anonl,
                            initial_vzero = vzero,
                            initial_ilast = ilast,
                            initial_curr = curr,
                            initial_cursub = cursub,
                            vnonl = nonlinear_voltages,
                            gslope = get(simultaneous_zno_config, :gslope, Float64[]),
                            gap_status_values = get(
                                simultaneous_zno_config,
                                :gap_status_values,
                                Float64[],
                            ),
                            t = get(simultaneous_zno_config, :t, 0.0),
                            deltat = get(simultaneous_zno_config, :deltat, 0.0),
                            epszno = get(simultaneous_zno_config, :epszno, 1.0e-10),
                            znolim = get(simultaneous_zno_config, :znolim, (0.2, 1.5)),
                            max_iterations = get(
                                simultaneous_zno_config,
                                :max_iterations,
                                20,
                            ),
                            fltinf = get(simultaneous_zno_config, :fltinf, Inf),
                            epsiln = epsilon,
                        )
                        rhs_values = zno_result.rhs
                        curr = zno_result.curr
                        cursub = zno_result.cursub
                        vzero = zno_result.vzero
                        ilast = zno_result.ilast
                        nonlinear_voltages = zno_result.vnonl
                        processed_count += zno_result.element_count
                        rhs_update_count += zno_result.rhs_update_count
                        cursub_update_count += zno_result.cursub_update_count
                        simultaneous_zno_solution_count += 1
                        simultaneous_zno_element_count += zno_result.element_count
                        simultaneous_piecewise_resistance_count +=
                            zno_result.piecewise_resistance_count
                        simultaneous_time_varying_resistance_count +=
                            zno_result.time_varying_resistance_count
                        simultaneous_zno_independent_count +=
                            zno_result.independent_count
                        simultaneous_zno_dependent_count +=
                            zno_result.dependent_count
                        simultaneous_zno_singular_count +=
                            zno_result.singular_count
                        simultaneous_zno_iteration_count +=
                            zno_result.iteration_count
                        simultaneous_zno_converged &=
                            zno_result.converged
                        simultaneous_zno_correction_scale_count +=
                            zno_result.correction_scale_count
                        simultaneous_zno_voltage_limit_count +=
                            zno_result.voltage_limit_count
                        simultaneous_zno_segment_search_count +=
                            zno_result.zno_segment_search_count
                        simultaneous_zno_max_voltage_correction = max(
                            simultaneous_zno_max_voltage_correction,
                            zno_result.max_voltage_correction,
                        )
                        append!(
                            simultaneous_zno_element_voltage,
                            zno_result.element_voltage,
                        )
                        simultaneous_zno_zthevenin =
                            _over16_append_block_diagonal(
                                simultaneous_zno_zthevenin,
                                zno_result.zthevenin,
                            )
                    end
                end
                continue
            end
            from_node = nonlinear_from_nodes[index]
            to_node = abs(nonlinear_to_nodes[index])
            1 <= from_node <= length(voltage_values) ||
                throw(ArgumentError("nonlinear_from_nodes entries must address voltages"))
            0 <= to_node <= length(voltage_values) ||
                throw(ArgumentError("nonlinear_to_nodes entries must address voltages"))
            voltage_delta =
                voltage_values[from_node] - _over16_node_voltage(voltage_values, to_node)

            if nonlinear_type == PIECEWISE_NONLINEAR_INDUCTOR_TYPE
                piecewise_nonlinear_inductor_update_count > 0 ||
                    throw(ArgumentError("true nonlinear-inductor update was not prepared"))
            elseif nonlinear_type == TRIGGERED_TIMED_RESISTANCE_TYPE
                table_start = nonlinear_admittance_nodes[index]
                table_end = abs(nonlinear_table_end_indices[index])
                1 <= table_start <= table_end <= length(cchar_values) ||
                    throw(ArgumentError("timed-resistance table range must address cchar"))
                table_end <= length(gslope_values) ||
                    throw(ArgumentError("timed-resistance table range must address gslope"))
                update = triggered_timed_resistance_step(
                    round(Int, curr[index]),
                    anonl[index],
                    voltage_delta,
                    table_start,
                    table_end,
                    cchar_values,
                    gslope_values;
                    trigger_voltage_v = nonlinear_voltages[index],
                    arm_time_s = timed_resistance_arm_times[index],
                    time_s = time_value,
                )
                curr[index] = update.active_segment
                anonl[index] = update.activation_time_s
                timed_resistance_accepted_currents[index] = update.accepted_current_a
                timed_resistance_admittance_deltas[index] = update.admittance_delta_s
                timed_resistance_update_count += 1
                update.activated && (timed_resistance_activation_count += 1)
                timed_resistance_segment_change_count += update.segment_change_count
                if update.admittance_delta_s != 0.0 && active_sparse_config !== nothing
                    restamp = saturated_transformer_sparse_admittance_update(
                        get(active_sparse_config, :km, Int[]),
                        saturated_transformer_sparse_ykm,
                        get(active_sparse_config, :kks, Int[]),
                        saturated_transformer_sparse_from_nodes[index],
                        saturated_transformer_sparse_to_nodes[index],
                        get(active_sparse_config, :partition_boundary,
                            length(get(active_sparse_config, :kks, Int[]))),
                        update.admittance_delta_s,
                    )
                    saturated_transformer_sparse_ykm = restamp.ykm
                    saturated_transformer_sparse_update_count += 1
                    saturated_transformer_sparse_positive_update_count +=
                        restamp.positive_update_count
                    saturated_transformer_sparse_negative_update_count +=
                        restamp.negative_update_count
                    saturated_transformer_sparse_partition_skip_count +=
                        restamp.partition_skip_count
                    restamp.retriangularization_required &&
                        (saturated_transformer_sparse_retriangularization_request_count += 1)
                end
                processed_count += 1
                continue
            elseif nonlinear_type == SWITCHING_NONLINEAR_RESISTOR_TYPE
                table_start = nonlinear_admittance_nodes[index]
                table_end = abs(nonlinear_table_end_indices[index])
                1 <= table_start <= table_end <= length(cchar_values) ||
                    throw(ArgumentError("switching resistor table range must address cchar"))
                table_end <= length(vchar_values) && table_end <= length(gslope_values) ||
                    throw(ArgumentError("switching resistor table range must address vchar and gslope"))
                rearm_index = rearm_indices[index]
                1 <= rearm_index <= length(vchar_values) ||
                    throw(ArgumentError("switching resistor rearm state must address vchar"))
                cursub_index = div(head, 5) + 1
                cursub_index <= length(cursub) ||
                    throw(ArgumentError("cursub length must cover switching resistor subsystem heads"))
                activation_segment_count = round(Int, anonl[index])
                update = switching_nonlinear_resistor_step(
                    curr[index],
                    cursub[cursub_index],
                    voltage_delta,
                    table_start,
                    table_end,
                    cchar_values,
                    gslope_values,
                    vchar_values;
                    turn_on_voltage_v = nonlinear_voltages[index],
                    turn_off_voltage_v = vzero[index],
                    activation_segment_count = activation_segment_count,
                    minimum_on_time_s = minimum_on_times[index],
                    rearm_time_s = vchar_values[rearm_index],
                    time_s = time_value,
                    single_flash = one_shot_flags[index],
                    infinity = infinity_value,
                    voltage_tolerance = flux_zero,
                )
                curr[index] = update.current_segment
                cursub[cursub_index] = update.companion_current_a
                nonlinear_voltages[index] = update.turn_on_voltage_v
                vchar_values[rearm_index] = update.rearm_time_s
                switching_resistor_reported_currents[index] = update.reported_current_a
                switching_resistor_accepted_currents[index] = update.accepted_current_a
                switching_resistor_admittance_deltas[index] = update.admittance_delta_s
                switching_resistor_companion_current_deltas[index] =
                    update.companion_current_delta_a
                if update.companion_current_a != 0.0
                    rhs_update_count += _over16_add_rhs_delta!(
                        rhs_values,
                        from_node,
                        -update.companion_current_a,
                    )
                    rhs_update_count += _over16_add_rhs_delta!(
                        rhs_values,
                        to_node,
                        update.companion_current_a,
                    )
                end
                cursub_update_count += 1
                switching_resistor_update_count += 1
                update.activated && (switching_resistor_activation_count += 1)
                update.deactivated && (switching_resistor_deactivation_count += 1)
                update.segment_changed && (switching_resistor_segment_change_count += 1)
                update.polarity_reversed && (switching_resistor_polarity_reversal_count += 1)
                if update.admittance_delta_s != 0.0 && active_sparse_config !== nothing
                    restamp = saturated_transformer_sparse_admittance_update(
                        get(active_sparse_config, :km, Int[]),
                        saturated_transformer_sparse_ykm,
                        get(active_sparse_config, :kks, Int[]),
                        saturated_transformer_sparse_from_nodes[index],
                        saturated_transformer_sparse_to_nodes[index],
                        get(active_sparse_config, :partition_boundary,
                            length(get(active_sparse_config, :kks, Int[]))),
                        update.admittance_delta_s,
                    )
                    saturated_transformer_sparse_ykm = restamp.ykm
                    saturated_transformer_sparse_update_count += 1
                    saturated_transformer_sparse_positive_update_count +=
                        restamp.positive_update_count
                    saturated_transformer_sparse_negative_update_count +=
                        restamp.negative_update_count
                    saturated_transformer_sparse_partition_skip_count +=
                        restamp.partition_skip_count
                    restamp.retriangularization_required &&
                        (saturated_transformer_sparse_retriangularization_request_count += 1)
                end
                processed_count += 1
                continue
            elseif nonlinear_type == HYSTERETIC_INDUCTOR_NONLINEAR_TYPE
                update = hysteretic_inductor_current_update(
                    cchar_values,
                    vchar_values,
                    gslope_values;
                    hysteresis_state_start_index = nonlinear_admittance_nodes[index],
                    major_loop_start_index = ilast[index],
                    stored_flux = nonlinear_voltages[index],
                    branch_voltage = voltage_delta,
                    delta2 = delta,
                    flzero = flux_zero,
                )
                cchar_values = update.cchar
                vchar_values = update.vchar
                gslope_values = update.gslope
                curr[index] = update.current
                anonl[index] = update.companion_current
                nonlinear_voltages[index] = update.stored_flux
                hysteretic_inductor_admittance_deltas[index] = update.admittance_delta
                hysteretic_inductor_source_current_deltas[index] = update.source_current_delta
                rhs_update_count +=
                    _over16_add_rhs_delta!(rhs_values, from_node, -update.source_current_delta)
                rhs_update_count +=
                    _over16_add_rhs_delta!(rhs_values, to_node, update.source_current_delta)
                hysteretic_inductor_rhs_update_count += 1
                hysteretic_inductor_update_count += 1
                hysteretic_inductor_reversal_count += update.reversal_count
                hysteretic_inductor_trajectory_limit_count += update.trajectory_limit_count
                if update.admittance_delta != 0.0 && active_sparse_config !== nothing
                    restamp = saturated_transformer_sparse_admittance_update(
                        get(active_sparse_config, :km, Int[]),
                        saturated_transformer_sparse_ykm,
                        get(active_sparse_config, :kks, Int[]),
                        saturated_transformer_sparse_from_nodes[index],
                        saturated_transformer_sparse_to_nodes[index],
                        get(
                            active_sparse_config,
                            :partition_boundary,
                            length(get(active_sparse_config, :kks, Int[])),
                        ),
                        update.admittance_delta,
                    )
                    saturated_transformer_sparse_ykm = restamp.ykm
                    saturated_transformer_sparse_update_count += 1
                    saturated_transformer_sparse_positive_update_count +=
                        restamp.positive_update_count
                    saturated_transformer_sparse_negative_update_count +=
                        restamp.negative_update_count
                    saturated_transformer_sparse_partition_skip_count +=
                        restamp.partition_skip_count
                    restamp.retriangularization_required &&
                        (saturated_transformer_sparse_retriangularization_request_count += 1)
                    hysteretic_inductor_sparse_update_count += 1
                end
                processed_count += 1
                continue
            elseif nonlinear_type == SATURATED_TRANSFORMER_NONLINEAR_TYPE
                table_start = nonlinear_admittance_nodes[index]
                table_end = abs(nonlinear_table_end_indices[index])
                1 <= table_start <= table_end <= length(cchar_values) ||
                    throw(ArgumentError("saturated transformer table range must address cchar"))
                table_end <= length(vchar_values) ||
                    throw(ArgumentError("saturated transformer table range must address vchar"))
                table_end <= length(gslope_values) ||
                    throw(ArgumentError("gslope must cover saturated transformer table ranges"))
                segment_update = saturated_transformer_segment_update(
                    curr[index],
                    anonl[index],
                    nonlinear_voltages[index],
                    voltage_delta,
                    table_start,
                    table_end,
                    cchar_values,
                    gslope_values,
                    vchar_values;
                    delta2 = delta,
                )
                curr[index] = segment_update.current_segment
                anonl[index] = segment_update.companion_current
                nonlinear_voltages[index] = segment_update.stored_voltage
                saturated_transformer_finit_deltas[index] =
                    segment_update.companion_current_delta
                saturated_transformer_admittance_deltas[index] =
                    segment_update.admittance_delta
                rhs_update_count += _over16_add_rhs_delta!(
                    rhs_values,
                    from_node,
                    segment_update.finitial_from_delta,
                )
                rhs_update_count += _over16_add_rhs_delta!(
                    rhs_values,
                    to_node,
                    segment_update.finitial_to_delta,
                )
                saturated_transformer_update_count += 1
                segment_update.segment_changed &&
                    (saturated_transformer_segment_change_count += 1)
                segment_update.polarity_mismatch_warning &&
                    (saturated_transformer_polarity_warning_count += 1)
                segment_update.polarity_reversed &&
                    (saturated_transformer_polarity_reversal_count += 1)
                if segment_update.segment_changed &&
                        active_sparse_config !== nothing
                    restamp = saturated_transformer_sparse_admittance_update(
                        get(active_sparse_config, :km, Int[]),
                        saturated_transformer_sparse_ykm,
                        get(active_sparse_config, :kks, Int[]),
                        saturated_transformer_sparse_from_nodes[index],
                        saturated_transformer_sparse_to_nodes[index],
                        get(
                            active_sparse_config,
                            :partition_boundary,
                            length(get(active_sparse_config, :kks, Int[])),
                        ),
                        segment_update.admittance_delta,
                    )
                    saturated_transformer_sparse_ykm = restamp.ykm
                    saturated_transformer_sparse_update_count += 1
                    saturated_transformer_sparse_positive_update_count +=
                        restamp.positive_update_count
                    saturated_transformer_sparse_negative_update_count +=
                        restamp.negative_update_count
                    saturated_transformer_sparse_partition_skip_count +=
                        restamp.partition_skip_count
                    restamp.retriangularization_required &&
                        (saturated_transformer_sparse_retriangularization_request_count += 1)
                end
            elseif nonlinear_type == 94
                a_start = abs(nonlinear_admittance_nodes[index])
                b_start = nonlinear_table_end_indices[index]
                b_start > 0 || throw(ArgumentError("type-94 nonlinear_table_end_indices entries must be positive"))
                vzero[index] = voltage_delta
                b6_index = b_start + 5
                b6_index <= length(vchar_values) ||
                    throw(ArgumentError("type-94 vchar state must include B(6)"))
                should_call_arrest =
                    vchar_values[b6_index] != 0.0 ||
                    abs(nonlinear_voltages[index]) <= abs(voltage_delta)
                if should_call_arrest
                    vchar_values[b_start] = voltage_delta
                    prior_region = vchar_values[b6_index]
                    curr[index] = _over16_type94_arrester_update!(
                        cchar_values,
                        vchar_values,
                        a_start,
                        b_start,
                        anonl[index],
                        voltage_delta,
                        curr[index],
                        delta,
                        dt,
                        epsilon,
                    )
                    type94_processed_count += 1
                    arrester_state_update_count += 1
                    if prior_region == 0.0
                        type94_flashover_count += 1
                    end
                    if vchar_values[b6_index] == 0.0
                        type94_cleared_count += 1
                    end
                else
                    type94_skip_count += 1
                end
            else
                table_start = nonlinear_admittance_nodes[index]
                table_end = abs(nonlinear_table_end_indices[index])
                1 <= table_start < table_end <= length(cchar_values) ||
                    throw(ArgumentError("nonlinear table range must address cchar/vchar"))
                table_end <= length(vchar_values) ||
                    throw(ArgumentError("nonlinear table range must address vchar"))
                table_start <= ilast[index] < table_end ||
                    throw(ArgumentError("ilast must lie inside the nonlinear table interval"))

                slope = anonl[index]
                last = ilast[index]
                voltage_delta = voltage_delta * delta + nonlinear_voltages[index]
                difference = voltage_delta - vzero[index]
                vzero[index] = voltage_delta
                if difference > 0.0
                    current, last, steps = _over16_find_increasing_nonlinear_current(
                        voltage_delta,
                        slope,
                        last,
                        table_end,
                        cchar_values,
                        vchar_values,
                    )
                    curr[index] = current
                    increasing_count += 1
                    interpolation_count += 1
                    table_limit_count += steps
                elseif difference < 0.0
                    current, last, steps = _over16_find_decreasing_nonlinear_current(
                        voltage_delta,
                        slope,
                        last,
                        table_start,
                        cchar_values,
                        vchar_values,
                    )
                    curr[index] = current
                    decreasing_count += 1
                    interpolation_count += 1
                    table_limit_count += steps
                else
                    unchanged_count += 1
                end
                ilast[index] = last
            end

            if nonlinear_type != SATURATED_TRANSFORMER_NONLINEAR_TYPE
                cursub_index = div(head, 5) + 1
                cursub_index <= length(cursub) ||
                    throw(ArgumentError("cursub length must cover nonlinear subsystem heads"))
                cursub[cursub_index] = curr[index]
                rhs_update_count += _over16_add_rhs_delta!(rhs_values, from_node, -curr[index])
                rhs_update_count += _over16_add_rhs_delta!(rhs_values, to_node, curr[index])
                cursub_update_count += 1
            end
            processed_count += 1
        end
    end

    return (
        source = :over16_nonlinear_current_compensation_update,
        outcome = :state_mutation,
        fortran_files = (:OVER16_FOR,),
        fortran_routines = (:SUBTS3,),
        fortran_labels = OVER16_NONLINEAR_CURRENT_COMPENSATION_LABELS,
        rhs = rhs_values,
        anonl = anonl,
        curr = curr,
        cursub = cursub,
        vzero = vzero,
        ilast = ilast,
        vnonl = nonlinear_voltages,
        cchar = cchar_values,
        vchar = vchar_values,
        gslope = gslope_values,
        processed_nonlinear_current_count = processed_count,
        nonlinear_current_interpolation_count = interpolation_count,
        increasing_current_count = increasing_count,
        decreasing_current_count = decreasing_count,
        unchanged_current_count = unchanged_count,
        negative_type_skip_count = negative_type_skip_count,
        missing_subsystem_skip_count = missing_subsystem_skip_count,
        simultaneous_zno_deferred_count = simultaneous_zno_deferred_count,
        simultaneous_zno_solution_count = simultaneous_zno_solution_count,
        simultaneous_zno_element_count = simultaneous_zno_element_count,
        simultaneous_piecewise_resistance_count = simultaneous_piecewise_resistance_count,
        simultaneous_time_varying_resistance_count = simultaneous_time_varying_resistance_count,
        simultaneous_zno_independent_count = simultaneous_zno_independent_count,
        simultaneous_zno_dependent_count = simultaneous_zno_dependent_count,
        simultaneous_zno_singular_count = simultaneous_zno_singular_count,
        simultaneous_zno_iteration_count = simultaneous_zno_iteration_count,
        simultaneous_zno_converged = simultaneous_zno_converged,
        simultaneous_zno_correction_scale_count = simultaneous_zno_correction_scale_count,
        simultaneous_zno_voltage_limit_count = simultaneous_zno_voltage_limit_count,
        simultaneous_zno_segment_search_count = simultaneous_zno_segment_search_count,
        simultaneous_zno_max_voltage_correction = simultaneous_zno_max_voltage_correction,
        simultaneous_zno_element_voltage = simultaneous_zno_element_voltage,
        simultaneous_zno_zthevenin = simultaneous_zno_zthevenin,
        simultaneous_zno_dispatch_slot_indices = simultaneous_zno_dispatch_slot_indices,
        simultaneous_zno_dispatch_next_slots = simultaneous_zno_dispatch_next_slots,
        simultaneous_zno_dispatch_from_nodes = simultaneous_zno_dispatch_from_nodes,
        simultaneous_zno_dispatch_to_nodes = simultaneous_zno_dispatch_to_nodes,
        simultaneous_zno_dispatch_nonlinear_indices = simultaneous_zno_dispatch_nonlinear_indices,
        simultaneous_zno_dispatch_type_codes = simultaneous_zno_dispatch_type_codes,
        simultaneous_zno_dispatch_nonlad = simultaneous_zno_dispatch_nonlad,
        simultaneous_zno_dispatch_nonle = simultaneous_zno_dispatch_nonle,
        simultaneous_zno_dispatch_anonl = simultaneous_zno_dispatch_anonl,
        simultaneous_zno_dispatch_initial_vzero = simultaneous_zno_dispatch_initial_vzero,
        simultaneous_zno_dispatch_initial_ilast = simultaneous_zno_dispatch_initial_ilast,
        simultaneous_zno_dispatch_initial_curr = simultaneous_zno_dispatch_initial_curr,
        simultaneous_zno_dispatch_initial_vnonl = simultaneous_zno_dispatch_initial_vnonl,
        simultaneous_zno_dispatch_gap_status = simultaneous_zno_dispatch_gap_status,
        simultaneous_zno_mutation_order = simultaneous_zno_solution_count > 0 ?
            (
                :subsystem_owner_scan,
                :simultaneous_zinc_oxide_solution,
                :state_vector_commit,
                :rhs_current_injection,
            ) :
            Symbol[],
        type94_processed_count = type94_processed_count,
        type94_flashover_count = type94_flashover_count,
        type94_cleared_count = type94_cleared_count,
        type94_skip_count = type94_skip_count,
        saturated_transformer_update_count = saturated_transformer_update_count,
        saturated_transformer_segment_change_count = saturated_transformer_segment_change_count,
        saturated_transformer_polarity_warning_count = saturated_transformer_polarity_warning_count,
        saturated_transformer_polarity_reversal_count = saturated_transformer_polarity_reversal_count,
        saturated_transformer_finit_deltas = saturated_transformer_finit_deltas,
        saturated_transformer_admittance_deltas = saturated_transformer_admittance_deltas,
        saturated_transformer_sparse_ykm = saturated_transformer_sparse_ykm,
        nonlinear_sparse_ykm = saturated_transformer_sparse_ykm,
        saturated_transformer_sparse_from_nodes = saturated_transformer_sparse_from_nodes,
        saturated_transformer_sparse_to_nodes = saturated_transformer_sparse_to_nodes,
        saturated_transformer_sparse_update_count = saturated_transformer_sparse_update_count,
        nonlinear_sparse_update_count = saturated_transformer_sparse_update_count,
        saturated_transformer_sparse_positive_update_count = saturated_transformer_sparse_positive_update_count,
        saturated_transformer_sparse_negative_update_count = saturated_transformer_sparse_negative_update_count,
        saturated_transformer_sparse_partition_skip_count = saturated_transformer_sparse_partition_skip_count,
        saturated_transformer_sparse_retriangularization_request_count =
            saturated_transformer_sparse_retriangularization_request_count,
        nonlinear_sparse_retriangularization_request_count =
            saturated_transformer_sparse_retriangularization_request_count,
        hysteretic_inductor_update_count = hysteretic_inductor_update_count,
        hysteretic_inductor_rhs_update_count = hysteretic_inductor_rhs_update_count,
        hysteretic_inductor_sparse_update_count = hysteretic_inductor_sparse_update_count,
        hysteretic_inductor_reversal_count = hysteretic_inductor_reversal_count,
        hysteretic_inductor_trajectory_limit_count = hysteretic_inductor_trajectory_limit_count,
        hysteretic_inductor_admittance_deltas = hysteretic_inductor_admittance_deltas,
        hysteretic_inductor_source_current_deltas = hysteretic_inductor_source_current_deltas,
        switching_resistor_update_count = switching_resistor_update_count,
        switching_resistor_activation_count = switching_resistor_activation_count,
        switching_resistor_deactivation_count = switching_resistor_deactivation_count,
        switching_resistor_segment_change_count = switching_resistor_segment_change_count,
        switching_resistor_polarity_reversal_count = switching_resistor_polarity_reversal_count,
        switching_resistor_reported_currents = switching_resistor_reported_currents,
        switching_resistor_accepted_currents = switching_resistor_accepted_currents,
        switching_resistor_admittance_deltas = switching_resistor_admittance_deltas,
        switching_resistor_companion_current_deltas =
            switching_resistor_companion_current_deltas,
        switching_resistor_deferred_effects = (),
        timed_resistance_update_count = timed_resistance_update_count,
        timed_resistance_activation_count = timed_resistance_activation_count,
        timed_resistance_segment_change_count = timed_resistance_segment_change_count,
        timed_resistance_accepted_currents = timed_resistance_accepted_currents,
        timed_resistance_admittance_deltas = timed_resistance_admittance_deltas,
        timed_resistance_deferred_effects = (),
        piecewise_nonlinear_inductor_update_count =
            piecewise_nonlinear_inductor_update_count,
        piecewise_nonlinear_inductor_segment_change_count =
            piecewise_nonlinear_inductor_segment_change_count,
        piecewise_nonlinear_inductor_iteration_count =
            piecewise_nonlinear_inductor_iteration_count,
        piecewise_nonlinear_inductor_accepted_currents =
            piecewise_nonlinear_inductor_accepted_currents,
        piecewise_nonlinear_inductor_accepted_fluxes =
            piecewise_nonlinear_inductor_accepted_fluxes,
        piecewise_nonlinear_inductor_accepted_voltages =
            piecewise_nonlinear_inductor_accepted_voltages,
        piecewise_nonlinear_inductor_network_response =
            piecewise_nonlinear_inductor_network_response,
        piecewise_nonlinear_inductor_deferred_effects = (),
        rhs_update_count = rhs_update_count,
        cursub_update_count = cursub_update_count,
        table_search_step_count = table_limit_count,
        arrester_state_update_count = arrester_state_update_count,
        nonlinear_current_compensation_applied = true,
        nonlinear_current_state_mutated = false,
        rhs_mutated = false,
        curr_mutated = false,
        cursub_mutated = false,
        vzero_mutated = false,
        ilast_mutated = false,
        vchar_mutated = false,
        mutation_order = piecewise_nonlinear_inductor_update_count > 0 ?
            (
                :voltage_difference,
                :table_current_or_arrester_update,
                :saturated_transformer_segment_update,
                :piecewise_nonlinear_inductor_coupled_solution,
                :hysteretic_inductor_update,
                :cursub_write,
                :rhs_current_injection,
            ) :
            (
                :voltage_difference,
                :table_current_or_arrester_update,
                :saturated_transformer_segment_update,
                :hysteretic_inductor_update,
                :cursub_write,
                :rhs_current_injection,
            ),
        deferred_calls = complete_nonlinear_source_loop ?
            Symbol[] :
            simultaneous_zno_deferred_count > 0 ?
            [:zincox_simultaneous_zno_solution, :bulk_last14_oracle] :
            simultaneous_zno_solution_count > 0 ?
            [:full_last14_card_execution, :bulk_last14_oracle] :
            [:bulk_last14_oracle],
        tacs_executed = false,
        solvum_executed = false,
        replacement_ready = complete_nonlinear_source_loop,
    )
end

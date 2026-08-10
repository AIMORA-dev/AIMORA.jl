
function coupled_lumped_sequence_augmented_step_context(
    parsed::DeckParser.DeckParseResult;
    dt_s::Float64 = 20e-6,
    t_end_s::Float64 = 0.0,
    time_switch_event_delay_s::Float64 = 0.0,
    current_zero_switching::Bool = false,
    recorded_step_indices = nothing,
    source_signal_provider::AbstractSourceSignalProvider = IdentitySourceSignalProvider(),
)
    DeckParser.assert_deck_valid!(parsed)
    source_equivalent = coupled_lumped_sequence_history_injection_elements(parsed)
    deck_elements = Any[parsed.elements...]
    deck_element_names = copy(parsed.element_names)
    _append_switching_nonlinear_resistor_safety_shunts!(
        deck_elements,
        deck_element_names,
        parsed,
    )
    _delay_deck_time_switch_events!(
        deck_elements,
        time_switch_event_delay_s,
        t_end_s,
    )
    current_zero_switching &&
        _convert_deck_current_zero_switches!(deck_elements, deck_element_names, parsed)
    elements = Any[deck_elements...; source_equivalent.elements...]
    element_names = Symbol[deck_element_names...; source_equivalent.element_names...]
    source_function_runtime, control_system_runtime =
        _append_dynamic_source_and_control_elements!(
        elements,
        element_names,
        parsed,
        dt_s,
        source_signal_provider,
    )
    node_count = maximum(values(parsed.node_map); init = 0)
    return initialize_step_context(
        NodalSystem(node_count, elements);
        node_map = parsed.node_map,
        element_names = element_names,
        source_function_runtime = source_function_runtime,
        control_system_runtime = control_system_runtime,
        _deck_runtime_output_context_kwargs(
            parsed;
            time_switch_event_delay_s = time_switch_event_delay_s,
            event_horizon_s = t_end_s,
        )...,
        dt_s = dt_s,
        t_end_s = t_end_s,
        source = parsed.source,
        recorded_step_indices = recorded_step_indices,
    )
end

function deck_saturated_transformer_nonlinear_current_config(
    parsed::DeckParser.DeckParseResult,
    saturated_transformer_intake;
    delta2::Real = 1.0,
    winding_number::Int = 1,
    saturated_transformer_sparse_config::Union{Nothing,NamedTuple} = nothing,
)
    DeckParser.assert_deck_valid!(parsed)
    arrays = saturated_transformer_nonlinear_arrays(saturated_transformer_intake)
    physical_node_map = _saturated_transformer_physical_node_map(parsed, arrays)
    voltage_source_node_names, terminal_node_names =
        _deck_saturated_transformer_winding_node_names(arrays, winding_number)
    voltage_source_node_indices = [
        _deck_saturated_transformer_node_index(physical_node_map, name)
        for name in voltage_source_node_names
    ]
    terminal_node_indices = [
        _deck_saturated_transformer_node_index(physical_node_map, name)
        for name in terminal_node_names
    ]
    branch_assembly = saturated_transformer_winding_branch_assembly(
        arrays,
        physical_node_map;
        nonlinear_winding_number = winding_number,
    )
    resolved_sparse_config = _deck_saturated_transformer_sparse_config(
        saturated_transformer_sparse_config,
        branch_assembly.nonlinear_from_node_indices,
        branch_assembly.nonlinear_to_node_indices,
    )
    current_config = saturated_transformer_winding_current_config(
        arrays,
        branch_assembly.nonlinear_from_node_indices,
        physical_node_map;
        delta2 = delta2,
        winding_number = winding_number,
        saturated_transformer_sparse_config = resolved_sparse_config,
    )
    required_node_count = maximum(
        vcat(
            current_config.nonlinear_from_nodes,
            current_config.nonlinear_to_nodes,
            branch_assembly.internal_top_node_indices,
        );
        init = 0,
    )
    safety_shunts = _saturated_transformer_safety_shunts(
        parsed,
        merge(
            current_config,
            (saturated_transformer_branch_assembly = branch_assembly,),
        ),
    )
    return merge(
        current_config,
        (
            saturated_transformer_branch_assembly = branch_assembly,
            saturated_transformer_internal_top_node_indices =
                branch_assembly.internal_top_node_indices,
            saturated_transformer_internal_top_voltage_source_nodes =
                voltage_source_node_indices,
            saturated_transformer_required_node_count = required_node_count,
            saturated_transformer_terminal_node_indices = terminal_node_indices,
            safety_shunt_resistor_names = [shunt.owner_name for shunt in safety_shunts],
            safety_shunt_resistance_ohm_values =
                [shunt.resistance_ohm for shunt in safety_shunts],
            safety_shunt_conductance_s_values =
                [shunt.conductance_s for shunt in safety_shunts],
        ),
    )
end

function _deck_zinc_oxide_ordered_rows(
    nonlinear_rows::Vector{DeckParser.DeckZincOxideNonlinearRow},
    initialization_rows::Vector{DeckParser.DeckZincOxideInitializationRow},
    breakpoint_rows::Vector{DeckParser.DeckZincOxideBreakpointRow},
)
    count = length(nonlinear_rows)
    initializations = Vector{DeckParser.DeckZincOxideInitializationRow}(undef, count)
    seen_initialization = falses(count)
    breakpoints = [
        DeckParser.DeckZincOxideBreakpointRow[] for _ in 1:count
    ]
    for row in initialization_rows
        1 <= row.nonlinear_row_index <= count ||
            throw(ArgumentError("zinc-oxide initialization row index is outside nonlinear rows"))
        !seen_initialization[row.nonlinear_row_index] ||
            throw(ArgumentError("zinc-oxide nonlinear row has duplicate initialization rows"))
        initializations[row.nonlinear_row_index] = row
        seen_initialization[row.nonlinear_row_index] = true
    end
    for row in breakpoint_rows
        1 <= row.nonlinear_row_index <= count ||
            throw(ArgumentError("zinc-oxide breakpoint row index is outside nonlinear rows"))
        push!(breakpoints[row.nonlinear_row_index], row)
    end
    for index in 1:count
        row = nonlinear_rows[index]
        row.source_kind == :copy_reference || continue
        reference_index = row.reference_index
        1 <= reference_index < index ||
            throw(ArgumentError("zinc-oxide COPY rows must reference a prior nonlinear row"))
        seen_initialization[reference_index] ||
            throw(ArgumentError("zinc-oxide COPY reference row has no initialization"))
        initializations[index] = initializations[reference_index]
        breakpoints[index] = breakpoints[reference_index]
        seen_initialization[index] = true
    end
    all(seen_initialization) ||
        throw(ArgumentError("every zinc-oxide nonlinear row requires one initialization row"))
    for rows in breakpoints
        isempty(rows) && throw(ArgumentError("every zinc-oxide nonlinear row requires at least one breakpoint"))
        sort!(rows; by = row -> row.breakpoint_index)
        for (index, row) in enumerate(rows)
            row.breakpoint_index == index ||
                throw(ArgumentError("zinc-oxide breakpoint indices must be consecutive"))
        end
    end
    return initializations, breakpoints
end

function _deck_zinc_oxide_initial_slope(
    initialization::DeckParser.DeckZincOxideInitializationRow,
    first_breakpoint::DeckParser.DeckZincOxideBreakpointRow,
)
    reference_voltage = Float64(initialization.reference_voltage)
    reference_voltage > 0.0 ||
        throw(ArgumentError("zinc-oxide reference voltage must be positive"))
    coefficient = Float64(first_breakpoint.current_coefficient)
    exponent = Float64(first_breakpoint.voltage_exponent)
    voltage = Float64(first_breakpoint.voltage)
    coefficient > 0.0 ||
        throw(ArgumentError("zinc-oxide breakpoint current coefficient must be positive"))
    voltage > 0.0 ||
        throw(ArgumentError("zinc-oxide breakpoint voltage must be positive"))
    return coefficient * voltage^exponent / (voltage * reference_voltage)
end

function _deck_zinc_oxide_gap_limit(
    initialization::DeckParser.DeckZincOxideInitializationRow,
    fltinf::Float64,
)
    reference_voltage = Float64(initialization.reference_voltage)
    gap_voltage = Float64(initialization.gap_voltage)
    gap_voltage <= 0.0 && return fltinf
    return gap_voltage * reference_voltage
end

_deck_reference_shifted_runtime_node(node_index::Int, deck_node_count::Int) =
    node_index == 0 ? 1 : deck_node_count + 2 - node_index

_deck_zinc_oxide_runtime_node(node_index::Int, deck_node_count::Int) =
    _deck_reference_shifted_runtime_node(node_index, deck_node_count)

function deck_zinc_oxide_nonlinear_current_config(
    parsed::DeckParser.DeckParseResult;
    delta2::Real = 1.0,
    fltinf::Real = 1.0e99,
    epszno::Real = 1.0e-10,
    znolim::Tuple{<:Real,<:Real} = (0.2, 1.5),
    max_iterations::Int = 20,
    znonl::AbstractVector{<:Real} = Float64[],
    ncomp::Union{Nothing,Int} = nothing,
    t::Real = 0.0,
    deltat::Real = 2.0 * Float64(delta2),
)
    DeckParser.assert_deck_valid!(parsed)
    nonlinear_rows = DeckParser.deck_zinc_oxide_nonlinear_rows(parsed)
    count = length(nonlinear_rows)
    count > 0 || throw(ArgumentError("deck has no zinc-oxide nonlinear rows"))
    initialization_rows, breakpoint_groups = _deck_zinc_oxide_ordered_rows(
        nonlinear_rows,
        DeckParser.deck_zinc_oxide_initialization_rows(parsed),
        DeckParser.deck_zinc_oxide_breakpoint_rows(parsed),
    )
    delta = Float64(delta2)
    isfinite(delta) && delta > 0.0 ||
        throw(ArgumentError("delta2 must be finite and positive"))
    infinity = Float64(fltinf)
    infinity > 0.0 || throw(ArgumentError("fltinf must be positive"))
    tolerance = Float64(epszno)
    isfinite(tolerance) && tolerance > 0.0 ||
        throw(ArgumentError("epszno must be finite and positive"))
    voltage_limit = (Float64(znolim[1]), Float64(znolim[2]))
    all(value -> isfinite(value) && value > 0.0, voltage_limit) ||
        throw(ArgumentError("znolim entries must be finite and positive"))
    max_iterations >= 0 || throw(ArgumentError("max_iterations must be nonnegative"))
    component_count = ncomp === nothing ? count : ncomp
    component_count > 0 || throw(ArgumentError("ncomp must be positive"))

    cchar = Float64[]
    gslope = Float64[]
    vchar = Float64[]
    table_start_indices = zeros(Int, count)
    table_end_indices = zeros(Int, count)
    initial_reference_voltages = zeros(Float64, count)
    initial_residual_voltages = zeros(Float64, count)
    initial_gap_limits = zeros(Float64, count)
    initial_table_indices = zeros(Int, count)
    for index in eachindex(nonlinear_rows)
        row = nonlinear_rows[index]
        if row.source_kind == :copy_reference
            reference_index = row.reference_index
            1 <= reference_index < index ||
                throw(ArgumentError("zinc-oxide COPY rows must reference a prior nonlinear row"))
            table_start_indices[index] = table_start_indices[reference_index]
            table_end_indices[index] = table_end_indices[reference_index]
            initial_reference_voltages[index] = initial_reference_voltages[reference_index]
            initial_residual_voltages[index] = initial_residual_voltages[reference_index]
            initial_gap_limits[index] = initial_gap_limits[reference_index]
            initial_table_indices[index] = initial_table_indices[reference_index]
            continue
        end
        initialization = initialization_rows[index]
        row_breakpoints = breakpoint_groups[index]
        table_start = length(cchar) + 1
        push!(cchar, _deck_zinc_oxide_initial_slope(initialization, first(row_breakpoints)))
        push!(gslope, 0.0)
        push!(vchar, 0.0)
        for breakpoint in row_breakpoints
            push!(cchar, Float64(breakpoint.current_coefficient))
            push!(gslope, Float64(breakpoint.voltage_exponent))
            push!(vchar, Float64(breakpoint.voltage))
        end
        table_end = length(cchar)
        table_start_indices[index] = table_start
        table_end_indices[index] = table_end
        initial_reference_voltages[index] = Float64(initialization.reference_voltage)
        initial_residual_voltages[index] = Float64(initialization.initial_voltage)
        initial_gap_limits[index] = _deck_zinc_oxide_gap_limit(initialization, infinity)
        initial_table_indices[index] = table_end
    end

    deck_node_count = length(parsed.node_map)
    shifted_from_nodes = [
        _deck_zinc_oxide_runtime_node(row.from_node_index, deck_node_count)
        for row in nonlinear_rows
    ]
    shifted_to_nodes = [
        _deck_zinc_oxide_runtime_node(row.to_node_index, deck_node_count)
        for row in nonlinear_rows
    ]
    deck_to_runtime_node_indices = [
        _deck_zinc_oxide_runtime_node(index, deck_node_count)
        for index in 1:deck_node_count
    ]
    required_node_count = max(
        deck_node_count + 1,
        maximum(vcat(shifted_from_nodes, shifted_to_nodes); init = 1),
    )
    source_activity_flags = zeros(Int, required_node_count)
    for node in vcat(shifted_from_nodes, shifted_to_nodes)
        if node > 1
            source_activity_flags[node] = 1
        end
    end
    last_slot = 1 + 5 * (count - 1)
    subnetwork_next_indices = zeros(Int, last_slot)
    subnetwork_from_nodes = zeros(Int, last_slot)
    subnetwork_to_nodes = zeros(Int, last_slot)
    subnetwork_nonlinear_indices = zeros(Int, last_slot)
    subnetwork_element_types = zeros(Int, last_slot)
    for index in 1:count
        slot = 1 + 5 * (index - 1)
        subnetwork_next_indices[slot] = index == count ? 1 : slot + 5
        subnetwork_from_nodes[slot] = shifted_from_nodes[index]
        subnetwork_to_nodes[slot] = shifted_to_nodes[index]
        subnetwork_nonlinear_indices[slot] = index
        subnetwork_element_types[slot] = 1
    end
    subsystem_owner_rows = zeros(Int, last_slot)
    subsystem_simultaneous_flags = zeros(Int, last_slot)
    subsystem_owner_rows[1] = 1
    subsystem_simultaneous_flags[1] = 1

    base_config = (
        source = :deck_zinc_oxide_nonlinear_current_config,
        outcome = :nonlinear_current_config,
        nonlinear_types = fill(921, count),
        nonlinear_from_nodes = shifted_from_nodes,
        nonlinear_to_nodes = shifted_to_nodes,
        nonlinear_admittance_nodes = table_start_indices,
        nonlinear_table_end_indices = table_end_indices,
        nonlinear_subsystem_indices = ones(Int, count),
        subsystem_begin_indices = [1],
        subsystem_owner_rows = subsystem_owner_rows,
        subsystem_simultaneous_flags = subsystem_simultaneous_flags,
        cchar = cchar,
        vchar = vchar,
        gslope = gslope,
        delta2 = delta,
        epszno = tolerance,
        znolim = voltage_limit,
        max_iterations = max_iterations,
        ncomp = component_count,
        fltinf = infinity,
        t = Float64(t),
        deltat = Float64(deltat),
        deck_owned_voltage_context = true,
        reference_shifted_voltage_context = true,
        deck_to_runtime_node_indices = deck_to_runtime_node_indices,
        reference_node_index = 1,
        nonlinear_required_node_count = required_node_count,
        initialize_nonlinear_state = true,
        initial_companion_current_values = initial_reference_voltages,
        initial_characteristic_current_values = initial_residual_voltages,
        initial_stored_voltage_values = initial_gap_limits,
        initial_current_segment_values = zeros(Float64, count),
        initial_table_index_values = initial_table_indices,
        initial_cursub_values = zeros(Float64, count + 1),
        nonlinear_current_segments = zeros(Float64, count),
        use_state_nonlinear_inverse_columns = true,
        nonlinear_inverse_config = (
            require_fortran_sparse_factor_result = true,
            ntot = required_node_count,
            ncomp = component_count,
            source_begin_indices = [1],
            source_next_indices = subnetwork_next_indices,
            source_from_nodes = subnetwork_from_nodes,
            source_to_nodes = subnetwork_to_nodes,
            source_activity_flags = source_activity_flags,
            nonlinear_types = fill(921, count),
            nonlinear_admittance_nodes = table_start_indices,
            nonlinear_from_nodes = shifted_from_nodes,
            nonlinear_to_nodes = shifted_to_nodes,
            nonlinear_source_flags = zeros(Int, count),
            partition_boundary = max(required_node_count - 1, 1),
            delta2 = delta,
            fltinf = infinity,
        ),
        nonlinear_source_begin_indices = [1],
        nonlinear_source_next_indices = subnetwork_next_indices,
        nonlinear_source_from_nodes = subnetwork_from_nodes,
        nonlinear_source_to_nodes = subnetwork_to_nodes,
        nonlinear_source_activity_flags = source_activity_flags,
        nonlinear_source_flags = zeros(Int, count),
        nonlinear_inverse_partition_boundary = max(required_node_count - 1, 1),
        nonlinear_deck_from_nodes = [row.from_node_index for row in nonlinear_rows],
        nonlinear_deck_to_nodes = [row.to_node_index for row in nonlinear_rows],
        nonlinear_deck_from_node_names = [row.from_node for row in nonlinear_rows],
        nonlinear_deck_to_node_names = [row.to_node for row in nonlinear_rows],
        subnetwork_next_indices = subnetwork_next_indices,
        subnetwork_from_nodes = subnetwork_from_nodes,
        subnetwork_to_nodes = subnetwork_to_nodes,
        subnetwork_nonlinear_indices = subnetwork_nonlinear_indices,
        subnetwork_element_types = subnetwork_element_types,
        table_entry_count = length(cchar),
        nonlinear_row_count = count,
        nonlinear_reference_indices = [row.reference_index for row in nonlinear_rows],
        nonlinear_source_kinds = [row.source_kind for row in nonlinear_rows],
        nonlinear_output_codes = [row.output_code for row in nonlinear_rows],
        nonlinear_owner_names = [row.name for row in nonlinear_rows],
        nonlinear_owner_line_numbers = [row.line_no for row in nonlinear_rows],
        fortran_files = (:OVER2_FOR, :OVER15_FOR, :OVER16_FOR),
        fortran_labels = (185, 21, 29, 31, 37, 38, 39, 182, 184, 269, 7234, 3522, 196),
        mutation_order = (
            :type_92_nonlinear_row,
            :initial_reference_gap_state,
            :characteristic_breakpoint_table,
            :copy_reference_table_pointer,
            :initial_linear_slope_backfill,
            :subnetwork_record_chain,
            :source_column_config,
            :nonlinear_state_seed,
        ),
        deferred_calls = (),
        complete_nonlinear_source_loop = true,
        replacement_ready = true,
    )
    isempty(znonl) && return base_config
    znonl_values = Float64.(znonl)
    length(znonl_values) == required_node_count * component_count ||
        throw(ArgumentError("znonl length must equal nonlinear_required_node_count * ncomp"))
    return merge(
        base_config,
        (
            simultaneous_zno_config = (
                znonl = znonl_values,
                ncomp = component_count,
                gslope = gslope,
                gap_status_values = zeros(Float64, count),
                t = Float64(t),
                deltat = Float64(deltat),
                subnetwork_next_indices = subnetwork_next_indices,
                subnetwork_from_nodes = subnetwork_from_nodes,
                subnetwork_to_nodes = subnetwork_to_nodes,
                subnetwork_nonlinear_indices = subnetwork_nonlinear_indices,
                subnetwork_element_types = subnetwork_element_types,
                subsystem_begin_index = 1,
                epszno = tolerance,
                znolim = voltage_limit,
                max_iterations = max_iterations,
                fltinf = infinity,
            ),
        ),
    )
end

function _deck_nonlinear_resistance_ordered_rows(
    nonlinear_rows::Vector{DeckParser.DeckNonlinearResistanceRow},
    initialization_rows::Vector{DeckParser.DeckNonlinearResistanceInitializationRow},
    point_rows::Vector{DeckParser.DeckNonlinearResistancePointRow},
)
    count = length(nonlinear_rows)
    initializations =
        Vector{DeckParser.DeckNonlinearResistanceInitializationRow}(undef, count)
    seen_initialization = falses(count)
    point_groups = [
        DeckParser.DeckNonlinearResistancePointRow[] for _ in 1:count
    ]
    for row in initialization_rows
        1 <= row.nonlinear_row_index <= count ||
            throw(ArgumentError("nonlinear resistance initialization row index is outside nonlinear rows"))
        !seen_initialization[row.nonlinear_row_index] ||
            throw(ArgumentError("nonlinear resistance row has duplicate initialization rows"))
        initializations[row.nonlinear_row_index] = row
        seen_initialization[row.nonlinear_row_index] = true
    end
    for row in point_rows
        1 <= row.nonlinear_row_index <= count ||
            throw(ArgumentError("nonlinear resistance point row index is outside nonlinear rows"))
        push!(point_groups[row.nonlinear_row_index], row)
    end
    for index in 1:count
        row = nonlinear_rows[index]
        row.source_kind == :copy_reference || continue
        reference_index = row.reference_index
        1 <= reference_index < index ||
            throw(ArgumentError("nonlinear resistance COPY rows must reference a prior nonlinear row"))
        nonlinear_rows[reference_index].nonlinear_type == row.nonlinear_type ||
            throw(ArgumentError("nonlinear resistance COPY reference type must match the current row"))
        seen_initialization[reference_index] ||
            throw(ArgumentError("nonlinear resistance COPY reference row has no initialization"))
        initializations[index] = initializations[reference_index]
        point_groups[index] = point_groups[reference_index]
        seen_initialization[index] = true
    end
    all(seen_initialization) ||
        throw(ArgumentError("every nonlinear resistance row requires one initialization row"))
    for rows in point_groups
        length(rows) >= 2 ||
            throw(ArgumentError("every nonlinear resistance row requires at least two points"))
        sort!(rows; by = row -> row.point_index)
        for (index, row) in enumerate(rows)
            row.point_index == index ||
                throw(ArgumentError("nonlinear resistance point indices must be consecutive"))
        end
    end
    return initializations, point_groups
end

function _deck_nonlinear_resistance_segments(
    row::DeckParser.DeckNonlinearResistanceRow,
    initialization::DeckParser.DeckNonlinearResistanceInitializationRow,
    points::Vector{DeckParser.DeckNonlinearResistancePointRow},
    infinity::Float64,
)
    coordinate_offset =
        row.element_kind == :piecewise_resistance ?
        Float64(initialization.table_voltage_offset) :
        0.0
    coordinates = [
        Float64(point.coordinate_value) + coordinate_offset * Float64(point.ordinate_value)
        for point in points
    ]
    ordinates = [Float64(point.ordinate_value) for point in points]
    reference_scale = maximum(abs.([first(coordinates), last(coordinates)]); init = 0.0)
    if row.element_kind == :piecewise_resistance
        raw_coordinates = [Float64(point.coordinate_value) for point in points]
        for index in 2:length(raw_coordinates)
            raw_coordinates[index] > raw_coordinates[index - 1] ||
                throw(ArgumentError("piecewise resistance point coordinates must be increasing"))
        end
        if coordinates[1] > 0.0 && ordinates[1] > 0.0
            pushfirst!(coordinates, 0.0)
            pushfirst!(ordinates, 0.0)
        end
    end

    segment_count = length(coordinates) - 1
    segment_count >= 1 ||
        throw(ArgumentError("nonlinear resistance table must produce at least one segment"))
    cchar = Vector{Float64}(undef, segment_count)
    gslope = Vector{Float64}(undef, segment_count)
    vchar = Vector{Float64}(undef, segment_count)
    for segment in 1:segment_count
        left_coordinate = coordinates[segment]
        right_coordinate = coordinates[segment + 1]
        left_ordinate = ordinates[segment]
        right_ordinate = ordinates[segment + 1]
        denominator = right_coordinate - left_coordinate
        slope =
            denominator == 0.0 && segment == segment_count ?
            0.0 :
            (right_ordinate - left_ordinate) / denominator
        intercept = right_ordinate - right_coordinate * slope
        cchar[segment] = intercept
        gslope[segment] = slope
        vchar[segment] = left_coordinate
    end
    if row.element_kind == :piecewise_resistance && vchar[1] < 0.0
        vchar[1] = -infinity
    end
    return cchar, gslope, vchar, reference_scale
end

function deck_nonlinear_resistance_current_config(
    parsed::DeckParser.DeckParseResult;
    delta2::Real = 1.0,
    fltinf::Real = 1.0e99,
    epszno::Real = 1.0e-10,
    znolim::Tuple{<:Real,<:Real} = (0.2, 1.5),
    max_iterations::Int = 20,
    znonl::AbstractVector{<:Real} = Float64[],
    ncomp::Union{Nothing,Int} = nothing,
    t::Real = 0.0,
    deltat::Real = 2.0 * Float64(delta2),
)
    DeckParser.assert_deck_valid!(parsed)
    nonlinear_rows = DeckParser.deck_nonlinear_resistance_rows(parsed)
    count = length(nonlinear_rows)
    count > 0 || throw(ArgumentError("deck has no nonlinear resistance rows"))
    initialization_rows, point_groups = _deck_nonlinear_resistance_ordered_rows(
        nonlinear_rows,
        DeckParser.deck_nonlinear_resistance_initialization_rows(parsed),
        DeckParser.deck_nonlinear_resistance_point_rows(parsed),
    )
    delta = Float64(delta2)
    isfinite(delta) && delta > 0.0 ||
        throw(ArgumentError("delta2 must be finite and positive"))
    infinity = Float64(fltinf)
    infinity > 0.0 || throw(ArgumentError("fltinf must be positive"))
    tolerance = Float64(epszno)
    isfinite(tolerance) && tolerance > 0.0 ||
        throw(ArgumentError("epszno must be finite and positive"))
    voltage_limit = (Float64(znolim[1]), Float64(znolim[2]))
    all(value -> isfinite(value) && value > 0.0, voltage_limit) ||
        throw(ArgumentError("znolim entries must be finite and positive"))
    max_iterations >= 0 || throw(ArgumentError("max_iterations must be nonnegative"))
    component_count = ncomp === nothing ? count : ncomp
    component_count > 0 || throw(ArgumentError("ncomp must be positive"))

    cchar = Float64[]
    gslope = Float64[]
    vchar = Float64[]
    table_start_indices = zeros(Int, count)
    table_end_indices = zeros(Int, count)
    fortran_initial_table_indices = zeros(Int, count)
    runtime_initial_table_indices = zeros(Int, count)
    initial_reference_values = zeros(Float64, count)
    initial_voltage_values = zeros(Float64, count)
    initial_gap_values = zeros(Float64, count)
    gap_status_values = zeros(Float64, count)
    nonlinear_types = zeros(Int, count)
    subnetwork_element_type_values = zeros(Int, count)
    segment_counts = zeros(Int, count)
    for index in eachindex(nonlinear_rows)
        row = nonlinear_rows[index]
        if row.source_kind == :copy_reference
            reference_index = row.reference_index
            1 <= reference_index < index ||
                throw(ArgumentError("nonlinear resistance COPY rows must reference a prior nonlinear row"))
            nonlinear_rows[reference_index].nonlinear_type == row.nonlinear_type ||
                throw(ArgumentError("nonlinear resistance COPY reference type must match the current row"))
            table_start_indices[index] = table_start_indices[reference_index]
            table_end_indices[index] = table_end_indices[reference_index]
            segment_counts[index] = segment_counts[reference_index]
            initial_voltage_values[index] = initial_voltage_values[reference_index]
            initial_gap_values[index] = initial_gap_values[reference_index]
            initial_reference_values[index] = initial_reference_values[reference_index]
            fortran_initial_table_indices[index] = fortran_initial_table_indices[reference_index]
            runtime_initial_table_indices[index] = runtime_initial_table_indices[reference_index]
            gap_status_values[index] = gap_status_values[reference_index]
            nonlinear_types[index] = nonlinear_types[reference_index]
            subnetwork_element_type_values[index] =
                subnetwork_element_type_values[reference_index]
            continue
        end
        initialization = initialization_rows[index]
        segments, slopes, boundaries, reference_scale = _deck_nonlinear_resistance_segments(
            row,
            initialization,
            point_groups[index],
            infinity,
        )
        table_start = length(cchar) + 1
        append!(cchar, segments)
        append!(gslope, slopes)
        append!(vchar, boundaries)
        table_end = length(cchar)
        table_start_indices[index] = table_start
        table_end_indices[index] = table_end
        segment_counts[index] = length(segments)
        initial_voltage_values[index] = Float64(initialization.initial_voltage)
        if row.element_kind == :time_varying_resistance
            start_voltage = Float64(initialization.table_voltage_offset)
            initial_reference_values[index] = 2.0 * start_voltage
            initial_gap_values[index] = start_voltage
            fortran_initial_table_indices[index] = 1
            runtime_initial_table_indices[index] = table_start
            gap_status_values[index] = -1.0
            nonlinear_types[index] = 923
            subnetwork_element_type_values[index] = 3
        else
            gap_voltage = Float64(initialization.gap_voltage)
            if gap_voltage < 0.0
                gap_voltage = infinity
            end
            initial_gap_values[index] = gap_voltage
            initial_reference_values[index] = reference_scale
            local_initial = gap_voltage == infinity || gap_voltage == 0.0 ? -1 : 1
            fortran_initial_table_indices[index] = local_initial
            runtime_initial_table_indices[index] = local_initial < 0 ? -table_start : table_start
            gap_status_values[index] = row.steady_state_reference
            nonlinear_types[index] = 922
            subnetwork_element_type_values[index] = 2
        end
    end

    deck_node_count = length(parsed.node_map)
    shifted_from_nodes = [
        _deck_reference_shifted_runtime_node(row.from_node_index, deck_node_count)
        for row in nonlinear_rows
    ]
    shifted_to_nodes = [
        _deck_reference_shifted_runtime_node(row.to_node_index, deck_node_count)
        for row in nonlinear_rows
    ]
    deck_to_runtime_node_indices = [
        _deck_reference_shifted_runtime_node(index, deck_node_count)
        for index in 1:deck_node_count
    ]
    required_node_count = max(
        deck_node_count + 1,
        maximum(vcat(shifted_from_nodes, shifted_to_nodes); init = 1),
    )
    source_activity_flags = zeros(Int, required_node_count)
    for node in vcat(shifted_from_nodes, shifted_to_nodes)
        node > 1 && (source_activity_flags[node] = 1)
    end
    last_slot = 1 + 5 * (count - 1)
    subnetwork_next_indices = zeros(Int, last_slot)
    subnetwork_from_nodes = zeros(Int, last_slot)
    subnetwork_to_nodes = zeros(Int, last_slot)
    subnetwork_nonlinear_indices = zeros(Int, last_slot)
    subnetwork_element_types = zeros(Int, last_slot)
    for index in 1:count
        slot = 1 + 5 * (index - 1)
        subnetwork_next_indices[slot] = index == count ? 1 : slot + 5
        subnetwork_from_nodes[slot] = shifted_from_nodes[index]
        subnetwork_to_nodes[slot] = shifted_to_nodes[index]
        subnetwork_nonlinear_indices[slot] = index
        subnetwork_element_types[slot] = subnetwork_element_type_values[index]
    end
    subsystem_owner_rows = zeros(Int, last_slot)
    subsystem_simultaneous_flags = zeros(Int, last_slot)
    subsystem_owner_rows[1] = 1
    subsystem_simultaneous_flags[1] = 1

    base_config = (
        source = :deck_nonlinear_resistance_current_config,
        outcome = :nonlinear_current_config,
        nonlinear_types = nonlinear_types,
        nonlinear_from_nodes = shifted_from_nodes,
        nonlinear_to_nodes = shifted_to_nodes,
        nonlinear_admittance_nodes = table_start_indices,
        nonlinear_table_end_indices = table_end_indices,
        nonlinear_subsystem_indices = ones(Int, count),
        subsystem_begin_indices = [1],
        subsystem_owner_rows = subsystem_owner_rows,
        subsystem_simultaneous_flags = subsystem_simultaneous_flags,
        cchar = cchar,
        vchar = vchar,
        gslope = gslope,
        delta2 = delta,
        epszno = tolerance,
        znolim = voltage_limit,
        max_iterations = max_iterations,
        ncomp = component_count,
        fltinf = infinity,
        t = Float64(t),
        deltat = Float64(deltat),
        deck_owned_voltage_context = true,
        reference_shifted_voltage_context = true,
        deck_to_runtime_node_indices = deck_to_runtime_node_indices,
        reference_node_index = 1,
        nonlinear_required_node_count = required_node_count,
        initialize_nonlinear_state = true,
        initial_companion_current_values = initial_reference_values,
        initial_characteristic_current_values = initial_voltage_values,
        initial_stored_voltage_values = initial_gap_values,
        initial_current_segment_values = zeros(Float64, count),
        initial_table_index_values = runtime_initial_table_indices,
        initial_cursub_values = zeros(Float64, count + 1),
        nonlinear_current_segments = zeros(Float64, count),
        use_state_nonlinear_inverse_columns = true,
        nonlinear_inverse_config = (
            require_fortran_sparse_factor_result = true,
            ntot = required_node_count,
            ncomp = component_count,
            source_begin_indices = [1],
            source_next_indices = subnetwork_next_indices,
            source_from_nodes = subnetwork_from_nodes,
            source_to_nodes = subnetwork_to_nodes,
            source_activity_flags = source_activity_flags,
            nonlinear_types = nonlinear_types,
            nonlinear_admittance_nodes = table_start_indices,
            nonlinear_from_nodes = shifted_from_nodes,
            nonlinear_to_nodes = shifted_to_nodes,
            nonlinear_source_flags = zeros(Int, count),
            partition_boundary = max(required_node_count - 1, 1),
            delta2 = delta,
            fltinf = infinity,
        ),
        nonlinear_source_begin_indices = [1],
        nonlinear_source_next_indices = subnetwork_next_indices,
        nonlinear_source_from_nodes = subnetwork_from_nodes,
        nonlinear_source_to_nodes = subnetwork_to_nodes,
        nonlinear_source_activity_flags = source_activity_flags,
        nonlinear_source_flags = zeros(Int, count),
        nonlinear_inverse_partition_boundary = max(required_node_count - 1, 1),
        nonlinear_deck_from_nodes = [row.from_node_index for row in nonlinear_rows],
        nonlinear_deck_to_nodes = [row.to_node_index for row in nonlinear_rows],
        nonlinear_deck_from_node_names = [row.from_node for row in nonlinear_rows],
        nonlinear_deck_to_node_names = [row.to_node for row in nonlinear_rows],
        fortran_initial_table_index_values = fortran_initial_table_indices,
        fortran_gap_status_values = gap_status_values,
        subnetwork_next_indices = subnetwork_next_indices,
        subnetwork_from_nodes = subnetwork_from_nodes,
        subnetwork_to_nodes = subnetwork_to_nodes,
        subnetwork_nonlinear_indices = subnetwork_nonlinear_indices,
        subnetwork_element_types = subnetwork_element_types,
        table_entry_count = length(cchar),
        nonlinear_row_count = count,
        nonlinear_resistance_segment_counts = segment_counts,
        nonlinear_resistance_element_kinds = [row.element_kind for row in nonlinear_rows],
        nonlinear_reference_indices = [row.reference_index for row in nonlinear_rows],
        nonlinear_source_kinds = [row.source_kind for row in nonlinear_rows],
        nonlinear_output_codes = [row.output_code for row in nonlinear_rows],
        nonlinear_owner_names = [row.name for row in nonlinear_rows],
        nonlinear_owner_line_numbers = [row.line_no for row in nonlinear_rows],
        fortran_files = (:OVER2_FOR, :OVER16_FOR),
        fortran_labels = (21, 40, 45, 46, 47, 49, 51, 52, 53, 54, 55, 70, 74, 75, 76, 78, 182, 184, 269),
        mutation_order = (
            :nonlinear_resistance_row,
            :initial_gap_state,
            :characteristic_point_table,
            :copy_reference_table_pointer,
            :optional_zero_point_insertion,
            :slope_intercept_table,
            :subnetwork_record_chain,
            :nonlinear_state_seed,
        ),
        deferred_calls = (),
        complete_nonlinear_source_loop = true,
        replacement_ready = true,
    )
    isempty(znonl) && return base_config
    znonl_values = Float64.(znonl)
    length(znonl_values) == required_node_count * component_count ||
        throw(ArgumentError("znonl length must equal nonlinear_required_node_count * ncomp"))
    return merge(
        base_config,
        (
            simultaneous_zno_config = (
                znonl = znonl_values,
                ncomp = component_count,
                gslope = gslope,
                gap_status_values = gap_status_values,
                t = Float64(t),
                deltat = Float64(deltat),
                subnetwork_next_indices = subnetwork_next_indices,
                subnetwork_from_nodes = subnetwork_from_nodes,
                subnetwork_to_nodes = subnetwork_to_nodes,
                subnetwork_nonlinear_indices = subnetwork_nonlinear_indices,
                subnetwork_element_types = subnetwork_element_types,
                subsystem_begin_index = 1,
                epszno = tolerance,
                znolim = voltage_limit,
                max_iterations = max_iterations,
                fltinf = infinity,
            ),
        ),
    )
end

function _deck_arrester_ordered_constants(
    nonlinear_rows::Vector{DeckParser.DeckArresterNonlinearRow},
    constant_rows::Vector{DeckParser.DeckArresterConstantRow},
)
    count = length(nonlinear_rows)
    values = [Float64[] for _ in 1:count]
    next_index = ones(Int, count)
    for row in constant_rows
        1 <= row.nonlinear_row_index <= count ||
            throw(ArgumentError("arrester constant row index is outside nonlinear rows"))
        row.first_constant_index == next_index[row.nonlinear_row_index] ||
            throw(ArgumentError("arrester constant rows must be consecutive"))
        append!(values[row.nonlinear_row_index], row.values)
        next_index[row.nonlinear_row_index] += length(row.values)
    end
    for (index, row) in enumerate(nonlinear_rows)
        if row.source_kind == :copy_reference
            1 <= row.reference_index < index ||
                throw(ArgumentError("arrester copy rows must reference a prior row"))
            isempty(values[index]) ||
                throw(ArgumentError("arrester copy rows must not carry duplicate constants"))
            values[index] = copy(values[row.reference_index])
        else
            length(values[index]) == 18 ||
                throw(ArgumentError("every direct arrester nonlinear row requires 18 constants"))
        end
    end
    return values
end

function _deck_extend_float_vector!(values::Vector{Float64}, target_length::Int)
    target_length >= 0 || throw(ArgumentError("target length must be nonnegative"))
    old_length = length(values)
    old_length >= target_length && return values
    resize!(values, target_length)
    for index in (old_length + 1):target_length
        values[index] = 0.0
    end
    return values
end

function deck_arrester_nonlinear_current_config(
    parsed::DeckParser.DeckParseResult;
    delta2::Real = 1.0,
    ncomp::Union{Nothing,Int} = nothing,
    t::Real = 0.0,
    deltat::Real = 2.0 * Float64(delta2),
)
    DeckParser.assert_deck_valid!(parsed)
    nonlinear_rows = DeckParser.deck_arrester_nonlinear_rows(parsed)
    count = length(nonlinear_rows)
    count > 0 || throw(ArgumentError("deck has no arrester nonlinear rows"))
    constant_groups =
        _deck_arrester_ordered_constants(
            nonlinear_rows,
            DeckParser.deck_arrester_constant_rows(parsed),
        )
    delta = Float64(delta2)
    isfinite(delta) && delta > 0.0 ||
        throw(ArgumentError("delta2 must be finite and positive"))
    component_count = ncomp === nothing ? count : ncomp
    component_count > 0 || throw(ArgumentError("ncomp must be positive"))

    cchar = Float64[]
    vchar = Float64[]
    constant_start_indices = zeros(Int, count)
    state_start_indices = zeros(Int, count)
    initial_state_indices = zeros(Int, count)
    initial_flashover_voltage_values = zeros(Float64, count)
    character_index = 0
    for index in eachindex(nonlinear_rows)
        row = nonlinear_rows[index]
        state_start = character_index + 1
        if row.source_kind == :copy_reference
            1 <= row.reference_index < index ||
                throw(ArgumentError("arrester copy rows must reference a prior row"))
            constant_start_indices[index] = constant_start_indices[row.reference_index]
            character_index += 11
        else
            constant_start_indices[index] = state_start
            _deck_extend_float_vector!(cchar, constant_start_indices[index] - 1)
            append!(cchar, constant_groups[index])
            character_index += 18
        end
        _deck_extend_float_vector!(vchar, state_start - 1)
        append!(vchar, zeros(Float64, 11))
        vchar[state_start + 9] = row.voltage_division_factor
        vchar[state_start + 10] = row.current_division_factor
        state_start_indices[index] = state_start
        initial_state_indices[index] = state_start
        initial_flashover_voltage_values[index] = row.flashover_voltage
    end

    deck_node_count = length(parsed.node_map)
    shifted_from_nodes = [
        _deck_reference_shifted_runtime_node(row.from_node_index, deck_node_count)
        for row in nonlinear_rows
    ]
    shifted_to_nodes = [
        _deck_reference_shifted_runtime_node(row.to_node_index, deck_node_count)
        for row in nonlinear_rows
    ]
    deck_to_runtime_node_indices = [
        _deck_reference_shifted_runtime_node(index, deck_node_count)
        for index in 1:deck_node_count
    ]
    required_node_count = max(
        deck_node_count + 1,
        maximum(vcat(shifted_from_nodes, shifted_to_nodes); init = 1),
    )
    source_activity_flags = zeros(Int, required_node_count)
    for node in vcat(shifted_from_nodes, shifted_to_nodes)
        node > 1 && (source_activity_flags[node] = 1)
    end
    last_slot = 1 + 5 * (count - 1)
    source_next_indices = zeros(Int, last_slot)
    source_from_nodes = zeros(Int, last_slot)
    source_to_nodes = zeros(Int, last_slot)
    source_begin_indices = zeros(Int, count)
    subsystem_begin_indices = zeros(Int, count)
    subsystem_owner_rows = zeros(Int, last_slot)
    subsystem_simultaneous_flags = zeros(Int, last_slot)
    for index in 1:count
        slot = 1 + 5 * (index - 1)
        source_next_indices[slot] = slot
        source_from_nodes[slot] = shifted_from_nodes[index]
        source_to_nodes[slot] = shifted_to_nodes[index]
        source_begin_indices[index] = slot
        subsystem_begin_indices[index] = slot
        subsystem_owner_rows[slot] = index
    end

    return (
        source = :deck_arrester_nonlinear_current_config,
        outcome = :nonlinear_current_config,
        nonlinear_types = fill(94, count),
        nonlinear_from_nodes = shifted_from_nodes,
        nonlinear_to_nodes = shifted_to_nodes,
        nonlinear_admittance_nodes = constant_start_indices,
        nonlinear_table_end_indices = state_start_indices,
        nonlinear_subsystem_indices = collect(1:count),
        subsystem_begin_indices = subsystem_begin_indices,
        subsystem_owner_rows = subsystem_owner_rows,
        subsystem_simultaneous_flags = subsystem_simultaneous_flags,
        cchar = cchar,
        vchar = vchar,
        gslope = zeros(Float64, length(cchar)),
        delta2 = delta,
        ncomp = component_count,
        t = Float64(t),
        deltat = Float64(deltat),
        deck_owned_voltage_context = true,
        reference_shifted_voltage_context = true,
        deck_to_runtime_node_indices = deck_to_runtime_node_indices,
        reference_node_index = 1,
        nonlinear_required_node_count = required_node_count,
        initialize_nonlinear_state = true,
        initial_companion_current_values = zeros(Float64, count),
        initial_characteristic_current_values = zeros(Float64, count),
        initial_stored_voltage_values = initial_flashover_voltage_values,
        initial_runtime_voltage_values = zeros(Float64, count),
        initial_current_segment_values = zeros(Float64, count),
        initial_table_index_values = initial_state_indices,
        initial_cursub_values = zeros(Float64, count + 1),
        use_state_nonlinear_inverse_columns = true,
        nonlinear_inverse_config = (
            require_fortran_sparse_factor_result = true,
            ntot = required_node_count,
            ncomp = component_count,
            source_begin_indices = source_begin_indices,
            source_next_indices = source_next_indices,
            source_from_nodes = source_from_nodes,
            source_to_nodes = source_to_nodes,
            source_activity_flags = source_activity_flags,
            nonlinear_types = fill(94, count),
            nonlinear_admittance_nodes = constant_start_indices,
            nonlinear_from_nodes = shifted_from_nodes,
            nonlinear_to_nodes = shifted_to_nodes,
            nonlinear_source_flags = zeros(Int, count),
            partition_boundary = max(required_node_count - 1, 1),
            delta2 = delta,
            fltinf = Inf,
        ),
        nonlinear_source_begin_indices = source_begin_indices,
        nonlinear_source_next_indices = source_next_indices,
        nonlinear_source_from_nodes = source_from_nodes,
        nonlinear_source_to_nodes = source_to_nodes,
        nonlinear_source_activity_flags = source_activity_flags,
        nonlinear_source_flags = zeros(Int, count),
        nonlinear_inverse_partition_boundary = max(required_node_count - 1, 1),
        nonlinear_deck_from_nodes = [row.from_node_index for row in nonlinear_rows],
        nonlinear_deck_to_nodes = [row.to_node_index for row in nonlinear_rows],
        nonlinear_deck_from_node_names = [row.from_node for row in nonlinear_rows],
        nonlinear_deck_to_node_names = [row.to_node for row in nonlinear_rows],
        nonlinear_reference_indices = [row.reference_index for row in nonlinear_rows],
        nonlinear_source_kinds = [row.source_kind for row in nonlinear_rows],
        nonlinear_output_codes = [row.output_code for row in nonlinear_rows],
        nonlinear_owner_names = [row.name for row in nonlinear_rows],
        nonlinear_owner_line_numbers = [row.line_no for row in nonlinear_rows],
        fortran_nonlinear_admittance_nodes = -constant_start_indices,
        fortran_nonlinear_state_start_indices = state_start_indices,
        table_entry_count = length(cchar),
        state_entry_count = length(vchar),
        nonlinear_row_count = count,
        fortran_files = (:OVER2_FOR, :OVER16_FOR),
        fortran_labels = (185, 270, 271, 4271, 182, 184, 269, 274, 3271, 1521, 1561, 8259),
        mutation_order = (
            :arrester_nonlinear_row,
            :state_zero_fill,
            :division_factor_defaults,
            :constant_table,
            :copy_reference_constant_pointer,
            :nonlinear_state_seed,
        ),
        deferred_calls = (),
        complete_nonlinear_source_loop = true,
        replacement_ready = true,
    )
end

function _deck_saturated_transformer_voltage_context(
    nonlinear_current_config::NamedTuple,
    voltage::AbstractVector{Float64},
)
    required_node_count = get(
        nonlinear_current_config,
        :saturated_transformer_required_node_count,
        length(voltage),
    )
    values = copy(voltage)
    if length(values) < required_node_count
        resize!(values, required_node_count)
        for index in (length(voltage) + 1):required_node_count
            values[index] = 0.0
        end
    end
    internal_nodes = get(
        nonlinear_current_config,
        :saturated_transformer_internal_top_node_indices,
        Int[],
    )
    source_nodes = get(
        nonlinear_current_config,
        :saturated_transformer_internal_top_voltage_source_nodes,
        Int[],
    )
    length(internal_nodes) == length(source_nodes) ||
        throw(ArgumentError("saturated transformer internal voltage source length mismatch"))
    for (internal_node, source_node) in zip(internal_nodes, source_nodes)
        internal_node <= length(voltage) && continue
        values[internal_node] = source_node == 0 ? 0.0 : voltage[source_node]
    end
    return values
end

function _deck_saturated_transformer_runtime_node(
    node_index::Int,
    deck_to_runtime_node_indices::AbstractVector{<:Integer},
)
    node_index == 0 && return 1
    if 1 <= node_index <= length(deck_to_runtime_node_indices)
        return Int(deck_to_runtime_node_indices[node_index])
    end
    return node_index + 1
end

function _deck_saturated_transformer_fill_internal_voltage!(
    values::Vector{Float64},
    nonlinear_current_config::NamedTuple,
)
    internal_nodes = Int.(get(
        nonlinear_current_config,
        :saturated_transformer_internal_top_node_indices,
        Int[],
    ))
    source_nodes = Int.(get(
        nonlinear_current_config,
        :saturated_transformer_internal_top_voltage_source_nodes,
        Int[],
    ))
    length(internal_nodes) == length(source_nodes) ||
        throw(ArgumentError("saturated transformer internal voltage source length mismatch"))
    for (internal_node, source_node) in zip(internal_nodes, source_nodes)
        1 <= internal_node <= length(values) ||
            throw(ArgumentError("saturated transformer internal voltage node is outside node count"))
        source_node == 0 && (values[internal_node] = 0.0; continue)
        1 <= source_node <= length(values) ||
            throw(ArgumentError("saturated transformer source voltage node is outside node count"))
        values[internal_node] = values[source_node]
    end
    return values
end

function _deck_nonlinear_voltage_context(
    nonlinear_current_config::NamedTuple,
    voltage::AbstractVector{Float64},
)
    if get(nonlinear_current_config, :reference_shifted_voltage_context, false)
        required_node_count = Int(
            get(
                nonlinear_current_config,
                :nonlinear_required_node_count,
                length(voltage) + 1,
            ),
        )
        required_node_count >= 1 ||
            throw(ArgumentError("reference-shifted nonlinear voltage context requires at least one node"))
        values = zeros(Float64, required_node_count)
        runtime_indices = Int.(get(
            nonlinear_current_config,
            :deck_to_runtime_node_indices,
            Int[],
        ))
        if isempty(runtime_indices)
            copied_count = min(length(voltage), required_node_count - 1)
            for index in 1:copied_count
                values[index + 1] = voltage[index]
            end
        else
            for index in 1:min(length(voltage), length(runtime_indices))
                runtime_index = runtime_indices[index]
                1 <= runtime_index <= required_node_count ||
                    throw(ArgumentError("deck nonlinear runtime voltage index is outside node count"))
                values[runtime_index] = voltage[index]
            end
        end
        return _deck_saturated_transformer_fill_internal_voltage!(
            values,
            nonlinear_current_config,
        )
    end
    return _deck_saturated_transformer_voltage_context(nonlinear_current_config, voltage)
end

function _deck_saturated_transformer_source_context!(
    over16_state::OVER16AcceptedTimestepState,
    required_node_count::Int,
)
    if required_node_count <= length(over16_state.source.f_values)
        _ensure_float_vector_length!(
            over16_state.source.nonlinear_current_compensation_values,
            required_node_count,
        )
        return nothing
    end
    previous_count = length(over16_state.source.f_values)
    resize!(over16_state.source.f_values, required_node_count)
    resize!(over16_state.source.e_values, required_node_count)
    resize!(over16_state.source.nonlinear_current_compensation_values, required_node_count)
    for index in (previous_count + 1):required_node_count
        over16_state.source.f_values[index] = 0.0
        over16_state.source.e_values[index] = 0.0
        over16_state.source.nonlinear_current_compensation_values[index] = 0.0
    end
    return nothing
end

function _over16_step_sparse_switch_state_flow_result(update)
    return hasproperty(update, :over16_sparse_switch_state_flow_result) ?
        update.over16_sparse_switch_state_flow_result : nothing
end

function _over16_step_pass_count(update)
    flow = _over16_step_sparse_switch_state_flow_result(update)
    return flow === nothing ? 1 : flow.pass_count
end

function _over16_step_mutation_count(
    update,
    accepted_field::Symbol,
    flow_count_field::Union{Nothing,Symbol}=nothing,
)
    flow = _over16_step_sparse_switch_state_flow_result(update)
    if flow !== nothing && flow_count_field !== nothing
        return getproperty(flow, flow_count_field)
    end
    return getproperty(update.over16_update, accepted_field) ? 1 : 0
end

function _over16_step_pass_updates(update)
    flow = _over16_step_sparse_switch_state_flow_result(update)
    return flow === nothing ? (update.over16_update,) : flow.pass_updates
end

function _over16_boundary_pass_updates(boundary_updates)
    updates = Any[]
    for step_update in boundary_updates
        append!(updates, _over16_step_pass_updates(step_update))
    end
    return updates
end

function _over16_boundary_pass_mutation_count(boundary_updates, field::Symbol)
    return count(update -> getproperty(update, field), _over16_boundary_pass_updates(boundary_updates))
end

function _over16_boundary_pass_bool_count(pass_updates, field::Symbol)
    return count(pass_updates) do update
        hasproperty(update, field) && Bool(getproperty(update, field))
    end
end

function _over16_boundary_pass_float_values(pass_updates, field::Symbol)
    values = Float64[]
    for update in pass_updates
        hasproperty(update, field) || continue
        append!(values, Float64.(getproperty(update, field)))
    end
    return values
end

const _OVER16_OWNER_PHASE_RESULT_FIELDS = (
    (:switch_scan_result, :switch_scan),
    (:switch_operation_result, :switch_operation),
    (:switch_status_result, :switch_status),
    (:switch_order_result, :switch_order),
    (:switch_admittance_result, :switch_topology_admittance),
    (:switch_retriangularization_result, :switch_retriangularization),
    (:switch_sparse_factor_result, :switch_sparse_factor_workspace),
    (:switch_fortran_sparse_factor_result, :switch_fortran_sparse_factor_workspace),
    (:nonlinear_source_column_result, :nonlinear_source_column),
    (:nonlinear_inverse_result, :nonlinear_inverse_column),
    (:switch_network_solution_result, :switch_network_solution),
    (:nonlinear_current_result, :nonlinear_current_compensation),
    (:switch_current_result, :switch_current),
    (:switch_post_current_result, :switch_post_current),
    (:switch_post_current_operation_queue_result, :switch_post_current_operation_queue),
    (:switch_bvalue_result, :switch_bvalue),
    (:switch_alteration_result, :switch_alteration),
    (:tacs_result, :tacs_utility),
    (:ntacs3_result, :ntacs3_utility),
    (:csup_pre_solve_result, :csup_pre_solve),
    (:tacs_linear_result, :tacs_linear_solve),
    (:tacs_post_solve_result, :tacs_post_solve),
    (:csup_termination_result, :csup_termination),
    (:output_report_result, :output_report),
    (:post_extrema_result, :post_extrema),
    (:source_result, :source_update),
)

function _over16_update_owner_phase_names(update)
    phase_names = Symbol[]
    for (field, phase_name) in _OVER16_OWNER_PHASE_RESULT_FIELDS
        if hasproperty(update, field) && getproperty(update, field) !== nothing
            push!(phase_names, phase_name)
        end
    end
    return phase_names
end

function _push_over16_owner_phase_trace!(
    step_indices::Vector{Int},
    pass_indices::Vector{Int},
    phase_names::Vector{Symbol},
    step_index::Int,
    pass_index::Int,
    update,
)
    for phase_name in _over16_update_owner_phase_names(update)
        push!(step_indices, step_index)
        push!(pass_indices, pass_index)
        push!(phase_names, phase_name)
    end
    return nothing
end

function _over16_owner_phase_trace_values(boundary_run)
    step_indices = Int[]
    pass_indices = Int[]
    phase_names = Symbol[]
    pass_lengths = Int[]
    for step_update in boundary_run.over16_updates
        step_index = hasproperty(step_update, :step_index) ? step_update.step_index : 0
        flow = _over16_step_sparse_switch_state_flow_result(step_update)
        if flow === nothing
            length_before = length(phase_names)
            _push_over16_owner_phase_trace!(
                step_indices,
                pass_indices,
                phase_names,
                step_index,
                1,
                step_update.over16_update,
            )
            push!(pass_lengths, length(phase_names) - length_before)
        else
            for (pass_index, pass_update) in enumerate(flow.pass_updates)
                length_before = length(phase_names)
                _push_over16_owner_phase_trace!(
                    step_indices,
                    pass_indices,
                    phase_names,
                    step_index,
                    pass_index,
                    pass_update,
                )
                push!(pass_lengths, length(phase_names) - length_before)
            end
        end
    end
    return (
        step_indices = step_indices,
        pass_indices = pass_indices,
        phase_names = phase_names,
        pass_lengths = pass_lengths,
    )
end

function _over16_owner_phase_family(phase_name::Symbol)
    phase_text = String(phase_name)
    if startswith(phase_text, "switch_")
        return :switch
    elseif startswith(phase_text, "nonlinear_")
        return :nonlinear
    elseif phase_name == :output_report
        return :output
    elseif phase_name == :post_extrema
        return :post_extrema
    elseif phase_name == :source_update
        return :source
    elseif phase_name in (
        :tacs_utility,
        :ntacs3_utility,
        :csup_pre_solve,
        :tacs_linear_solve,
        :tacs_post_solve,
        :csup_termination,
    )
        return :tacs_elec
    end
    return :other
end

function _over16_owner_family_trace_values(owner_phase_trace)
    phase_families = Symbol[
        _over16_owner_phase_family(phase_name)
        for phase_name in owner_phase_trace.phase_names
    ]
    family_names = Symbol[]
    family_counts = Int[]
    for family in phase_families
        index = findfirst(==(family), family_names)
        if index === nothing
            push!(family_names, family)
            push!(family_counts, 1)
        else
            family_counts[index] += 1
        end
    end
    pass_family_counts = Int[]
    offset = 0
    for pass_length in owner_phase_trace.pass_lengths
        pass_families = Symbol[]
        for family in phase_families[(offset + 1):(offset + pass_length)]
            family in pass_families || push!(pass_families, family)
        end
        push!(pass_family_counts, length(pass_families))
        offset += pass_length
    end
    return (
        phase_families = phase_families,
        family_names = family_names,
        family_counts = family_counts,
        pass_family_counts = pass_family_counts,
        cross_family_pass_count = count(>(1), pass_family_counts),
    )
end

function _over16_sparse_switch_state_flow_repeat_suppressed_config_keys(boundary_run)
    suppressed_keys = Symbol[]
    for flow in _over16_sparse_switch_state_flow_results(boundary_run)
        hasproperty(flow, :repeat_pass_suppressed_config_keys) || continue
        for key in flow.repeat_pass_suppressed_config_keys
            key in suppressed_keys || push!(suppressed_keys, key)
        end
    end
    return suppressed_keys
end

function _deck_over16_output_step_config(
    plan::DeckOVER16BoundaryPlan,
    context::EMTStepContext,
    ::DeckStepConfigFeatures{O,BV,BC,BP,S,SC,ST,SA,SX},
) where {O,BV,BC,BP,S,SC,ST,SA,SX}
    O || return NamedTuple()
    kwargs = (
        ivolt = 0,
        iaverg = 0,
        selected_node_indices = plan.output_node_indices,
    )
    config = (
        output_report_config = (
            t = context.t_s,
            deltat = context.dt_s,
            istep = context.step_index,
            kwargs = kwargs,
        ),
    )
    branch_voltage_config =
        BV ?
        (
            deck_branch_voltage_output_config = (
                output_names = plan.branch_voltage_output_names,
                branch_names = plan.branch_voltage_branch_names,
                branch_indices = plan.branch_voltage_branch_indices,
            ),
        ) : NamedTuple()
    branch_current_config =
        BC ?
        (
            deck_branch_current_output_config = (
                output_names = plan.branch_current_output_names,
                branch_names = plan.branch_current_branch_names,
                branch_indices = plan.branch_current_branch_indices,
            ),
        ) : NamedTuple()
    branch_power_config =
        BP ?
        (
            deck_branch_power_output_config = (
                output_names = plan.branch_power_output_names,
                branch_names = plan.branch_power_branch_names,
                branch_indices = plan.branch_power_branch_indices,
            ),
        ) : NamedTuple()
    return (
        output_report_config = merge(
            config.output_report_config,
            branch_voltage_config,
            branch_current_config,
            branch_power_config,
        ),
    )
end

function _deck_over16_post_advance_source_time(context::EMTStepContext)
    # SUBTS3 advances T at label 1661 before source-row labels 1249-1300.
    return _deck_over16_source_update_time(context, context.dt_s)
end

function _deck_over16_source_update_time(
    context::EMTStepContext,
    source_time_offset_s::Float64,
)
    return context.t_s + source_time_offset_s
end

function _deck_over16_source_card_step_kwargs(plan::DeckOVER16BoundaryPlan,
                                              context::EMTStepContext)
    runtime = context.source_function_runtime
    row_index = runtime === nothing ?
        context.step_index + 1 :
        runtime.next_input_row_index
    if row_index > plan.source_card_row_count
        row_index <= plan.source_interpolation_row_count ||
            return plan.source_card_row_count == 0 ? NamedTuple() : (nchain = 16,)
        return (
            interpolated_values = copy(plan.source_interpolation_values[row_index]),
        )
    end
    kind = plan.source_card_kinds[row_index]
    card_values = plan.source_card_values[row_index]
    card_kwargs = (
        nchain = 17,
        kolbeg = kind == :free_field ? 1 : 0,
        card_values = card_values,
    )
    identity_values =
        !isempty(card_values) && card_values[1] == 9999.0 ?
        zeros(Float64, 10) :
        copy(card_values)
    row_index > plan.source_interpolation_row_count && return merge(
        card_kwargs,
        (
            interpolated_values = identity_values,
            source_signal_provider_applied = false,
        ),
    )
    return merge(
        card_kwargs,
        (
            interpolated_values = copy(plan.source_interpolation_values[row_index]),
            source_signal_provider_applied = true,
        ),
    )
end

function _deck_over16_source_tacs_override_kwargs(plan::DeckOVER16BoundaryPlan)
    plan.source_tacs_override_count == 0 && return NamedTuple()
    nstacs = maximum(plan.source_tacs_override_positions)
    vstacs_indices = zeros(Int, nstacs)
    for (position, xtcs_index) in
        zip(plan.source_tacs_override_positions, plan.source_tacs_override_xtcs_indices)
        vstacs_indices[position] = xtcs_index
    end
    return (
        vstacs_indices = vstacs_indices,
        nstacs = nstacs,
    )
end

function _deck_over16_source_analytic_kwargs(plan::DeckOVER16BoundaryPlan,
                                             context::EMTStepContext)
    runtime = context.source_function_runtime
    row_index = runtime === nothing ?
        context.step_index + 1 :
        runtime.next_input_row_index
    row_index <= plan.source_analytic_row_count || return NamedTuple()
    return (
        kanal = 1,
        analytic_values = copy(plan.source_analytic_values[row_index]),
    )
end

function _deck_source_card_step_kwargs(
    plan::DeckOVER16BoundaryPlan,
    context::EMTStepContext,
    ::Val{SC},
) where {SC}
    return SC ? _deck_over16_source_card_step_kwargs(plan, context) : NamedTuple()
end

function _deck_source_tacs_override_kwargs(
    plan::DeckOVER16BoundaryPlan,
    ::Val{ST},
) where {ST}
    return ST ? _deck_over16_source_tacs_override_kwargs(plan) : NamedTuple()
end

function _deck_source_analytic_kwargs(
    plan::DeckOVER16BoundaryPlan,
    context::EMTStepContext,
    ::Val{SA},
) where {SA}
    return SA ? _deck_over16_source_analytic_kwargs(plan, context) : NamedTuple()
end

function _deck_over5a_source_step_config(
    plan::DeckOVER16BoundaryPlan,
    context::EMTStepContext,
    over16_state::OVER16AcceptedTimestepState,
    ::DeckStepConfigFeatures{O,BV,BC,BP,S,SC,ST,SA,SX},
    source_time_offset_s::Float64,
) where {O,BV,BC,BP,S,SC,ST,SA,SX}
    S || return NamedTuple()
    source_card_kwargs =
        _deck_source_card_step_kwargs(plan, context, Val(SC))
    source_time_s = _deck_over16_source_update_time(context, source_time_offset_s)
    source_runtime = context.source_function_runtime
    if SX && source_runtime !== nothing &&
       source_signal_interpolation_active(source_runtime.signal_provider) &&
       !get(source_card_kwargs, :source_signal_provider_applied, false) &&
       !(over16_state.source.source_card.iread != 0 &&
         get(source_card_kwargs, :nchain, 16) == 16)
        input_values = over16_state.source.source_card.iread == 0 ?
            zeros(Float64, 10) :
            Float64.(get(
                source_card_kwargs,
                :card_values,
                over16_state.source.source_card.voltbc_values,
            ))
        if !isempty(input_values) && input_values[1] == 9999.0
            fill!(input_values, 0.0)
        end
        source_card_kwargs = merge(
            source_card_kwargs,
            (
                interpolated_values = source_signal_values(
                    source_runtime.signal_provider,
                    input_values,
                    source_time_s,
                ),
                source_signal_provider_applied = true,
            ),
        )
    end
    source_tacs_override_kwargs =
        _deck_source_tacs_override_kwargs(plan, Val(ST))
    source_analytic_kwargs =
        _deck_source_analytic_kwargs(plan, context, Val(SA))
    control_signal_values = if ST && context.control_system_runtime !== nothing
        referenced_indices = Int[]
        append!(referenced_indices, plan.source_tacs_override_xtcs_indices)
        for index in eachindex(plan.source_iform_values)
            source_type = abs(plan.source_iform_values[index])
            if source_type == 17 || source_type >= 60
                push!(referenced_indices, round(Int, plan.source_sfreq_values[index]))
            end
        end
        _control_system_signal_slot_values(
            context;
            maximum_index = isempty(referenced_indices) ? 0 : maximum(referenced_indices),
        )
    else
        over16_state.tacs.xtcs_values
    end
    if SA && isempty(source_analytic_kwargs) && source_runtime !== nothing &&
       source_signal_analytic_active(source_runtime.signal_provider)
        analytic_program = _source_signal_program_analytic_kwargs(
            source_runtime.signal_provider,
            over16_state.source.source_card,
            source_card_kwargs,
            merge(
                source_tacs_override_kwargs,
                (xtcs_values = control_signal_values,),
            ),
            source_time_s,
        )
        source_analytic_kwargs = analytic_program.kwargs
    elseif SA && isempty(source_analytic_kwargs) && source_runtime !== nothing &&
           source_runtime.internal_analytic_requested
        throw(ArgumentError(
            "analytic source usage requires a typed SourceSignalProgram",
        ))
    end
    source_kwargs = merge(
        (
            kconst = plan.source_row_count,
            crest_values = plan.source_crest_values,
            time1_values = plan.source_time1_values,
            time2_values = plan.source_time2_values,
            sfreq_values = plan.source_sfreq_values,
            xtcs_values = control_signal_values,
            delta2 = context.dt_s / 2.0,
        ),
        source_card_kwargs,
        source_tacs_override_kwargs,
        source_analytic_kwargs,
    )
    return (
        source_config = (
            node_values = plan.source_node_values,
            iform_values = plan.source_iform_values,
            tstart_values = plan.source_tstart_values,
            tstop_values = plan.source_tstop_values,
            t = source_time_s,
            constraint_t = context.t_s,
            kwargs = source_kwargs,
        ),
    )
end

function _deck_over16_post_extrema_begmax_values(over16_state::OVER16AcceptedTimestepState)
    value_count = max(over16_state.post_extrema.maxout + 1, 3)
    values = fill(Inf, value_count)
    values[1] = 0.0
    values[2] = 0.0
    return values
end

function _deck_over16_post_extrema_step_config(over16_state::OVER16AcceptedTimestepState)
    return (
        post_extrema_config = (
            kwargs = (
                begmax_values = _deck_over16_post_extrema_begmax_values(over16_state),
            ),
        ),
    )
end

function _deck_nonlinear_inverse_step_config(
    nonlinear_current_config::NamedTuple,
    over16_state::OVER16AcceptedTimestepState,
    ;
    allow_pending_fortran_sparse_factor::Bool=false,
)
    haskey(nonlinear_current_config, :nonlinear_inverse_config) || return NamedTuple()
    inverse_config = nonlinear_current_config.nonlinear_inverse_config
    sparse_workspace_ready =
        over16_state.switch_fortran_sparse_factor.workspace_update_count > 0 ||
        allow_pending_fortran_sparse_factor ||
        !get(inverse_config, :require_fortran_sparse_factor_result, true)
    sparse_workspace_ready || return NamedTuple()
    required_ntot = Int(get(inverse_config, :ntot, over16_state.nonlinear_inverse.ntot))
    if length(over16_state.switch_topology.kode) != required_ntot
        if all(iszero, over16_state.switch_topology.kode) && required_ntot > 0
            resize!(over16_state.switch_topology.kode, required_ntot)
            fill!(over16_state.switch_topology.kode, 0)
        else
            return NamedTuple()
        end
    end

    prepared = merge(
        inverse_config,
        (
            require_fortran_sparse_factor_result = false,
            partition_boundary = length(over16_state.switch_fortran_sparse_factor.kk),
        ),
    )
    seed_initial_state =
        over16_state.nonlinear_inverse.source_column_update_count == 0 &&
        over16_state.nonlinear_inverse.update_count == 0 &&
        over16_state.nonlinear_inverse.current_update_count == 0
    if seed_initial_state
        prepared = merge(
            prepared,
            (
                anonl = get(
                    nonlinear_current_config,
                    :initial_companion_current_values,
                    Float64[],
                ),
                vzero = get(
                    nonlinear_current_config,
                    :initial_characteristic_current_values,
                    Float64[],
                ),
                vnonl = get(
                    nonlinear_current_config,
                    :initial_runtime_voltage_values,
                    get(
                        nonlinear_current_config,
                        :initial_stored_voltage_values,
                        Float64[],
                    ),
                ),
                ilast = get(
                    nonlinear_current_config,
                    :initial_table_index_values,
                    Int[],
                ),
                curr = get(
                    nonlinear_current_config,
                    :initial_current_segment_values,
                    Float64[],
                ),
                cursub = get(
                    nonlinear_current_config,
                    :initial_cursub_values,
                    Float64[],
                ),
                vchar = get(nonlinear_current_config, :vchar, Float64[]),
            ),
        )
    end
    return (nonlinear_inverse_config = prepared,)
end

function _deck_nonlinear_current_step_config(
    nonlinear_current_config,
    over16_state::OVER16AcceptedTimestepState,
    context::EMTStepContext,
    ;
    allow_pending_fortran_sparse_factor::Bool=false,
)
    nonlinear_current_config === nothing && return NamedTuple()
    return merge(
        _deck_nonlinear_inverse_step_config(
            nonlinear_current_config,
            over16_state;
            allow_pending_fortran_sparse_factor =
                allow_pending_fortran_sparse_factor,
        ),
        (
            nonlinear_current_config = merge(
                nonlinear_current_config,
                (
                    deck_owned_voltage_context = true,
                    initialize_nonlinear_state = true,
                    t = context.t_s,
                    deltat = context.dt_s,
                ),
            ),
        ),
    )
end

function _reference_shifted_node(index::Int)
    index >= 0 || throw(ArgumentError("node group endpoint must be nonnegative"))
    return index + 1
end

function _deck_time_switch_closed_mask!(
    closed_mask::AbstractVector{Bool},
    plan::DeckOVER16BoundaryPlan,
    time_s::Float64,
)
    length(closed_mask) == plan.switch_count || throw(ArgumentError(
        "time-switch closed-mask length must match the boundary plan",
    ))
    @inbounds for index in eachindex(closed_mask)
        closed_mask[index] = _deck_time_switch_closed_at(
            plan.switch_initially_closed_flags[index],
            plan.switch_close_time_s_values[index],
            plan.switch_open_time_s_values[index],
            time_s,
        )
    end
    return closed_mask
end

function DeckTimeSwitchStepWorkspace(
    plan::DeckOVER16BoundaryPlan,
    node_count::Int,
)
    switch_count = plan.switch_count
    shifted_from_nodes = Vector{Int}(undef, switch_count)
    shifted_to_nodes = Vector{Int}(undef, switch_count)
    @inbounds for index in 1:switch_count
        shifted_from_nodes[index] =
            _reference_shifted_node(plan.switch_from_node_indices[index])
        shifted_to_nodes[index] =
            _reference_shifted_node(plan.switch_to_node_indices[index])
    end
    operation_queue = Int[]
    grouped_from_nodes = Int[]
    grouped_to_nodes = Int[]
    sizehint!(operation_queue, switch_count)
    sizehint!(grouped_from_nodes, switch_count)
    sizehint!(grouped_to_nodes, switch_count)
    return DeckTimeSwitchStepWorkspace(
        falses(switch_count),
        falses(switch_count),
        shifted_from_nodes,
        shifted_to_nodes,
        operation_queue,
        zeros(Float64, switch_count),
        zeros(Int, switch_count),
        fill(Inf, switch_count),
        zeros(Float64, switch_count),
        zeros(Float64, switch_count),
        zeros(Int, node_count),
        grouped_from_nodes,
        grouped_to_nodes,
    )
end

function _deck_time_switch_reference_shifted_nodes(plan::DeckOVER16BoundaryPlan)
    shifted_from_nodes = Vector{Int}(undef, plan.switch_count)
    shifted_to_nodes = Vector{Int}(undef, plan.switch_count)
    @inbounds for index in 1:plan.switch_count
        shifted_from_nodes[index] =
            _reference_shifted_node(plan.switch_from_node_indices[index])
        shifted_to_nodes[index] =
            _reference_shifted_node(plan.switch_to_node_indices[index])
    end
    return shifted_from_nodes, shifted_to_nodes
end

_deck_time_switch_position(closed::Bool) = closed ? 2 : 5

function _deck_time_switch_operation_queue!(
    modswt::Vector{Int},
    previous_closed_mask::AbstractVector{Bool},
    requested_closed_mask::AbstractVector{Bool},
)
    length(previous_closed_mask) == length(requested_closed_mask) ||
        throw(ArgumentError("time-switch closed-mask lengths must match"))
    empty!(modswt)
    for index in eachindex(previous_closed_mask, requested_closed_mask)
        if !previous_closed_mask[index] && requested_closed_mask[index]
            push!(modswt, index)
        elseif previous_closed_mask[index] && !requested_closed_mask[index]
            push!(modswt, -index)
        end
    end
    return modswt
end

function _sync_deck_time_switch_over16_state!(
    over16_state::OVER16AcceptedTimestepState,
    initial_closed_mask::AbstractVector{Bool},
    node_count::Int,
)
    switch_count = length(initial_closed_mask)
    initialize_deck_switch_state =
        length(over16_state.switch_topology.closed_mask) != switch_count ||
        (
            over16_state.switch_operation.accumulated_operation_count == 0 &&
            over16_state.switch_admittance.retriangularization_count == 0 &&
            over16_state.switch_fortran_sparse_factor.workspace_update_count == 0
        )
    if initialize_deck_switch_state
        resize!(over16_state.switch_topology.closed_mask, switch_count)
        over16_state.switch_topology.closed_mask .= initial_closed_mask
        over16_state.switch_topology.closed_switch_count = count(identity, initial_closed_mask)
        over16_state.switch_topology.first_group_head = 0
    end
    length(over16_state.switch_topology.nextsw) == switch_count ||
        resize!(over16_state.switch_topology.nextsw, switch_count)
    length(over16_state.switch_topology.kode) == node_count ||
        resize!(over16_state.switch_topology.kode, node_count)
    fill!(over16_state.switch_topology.nextsw, 0)
    fill!(over16_state.switch_topology.kode, 0)

    if initialize_deck_switch_state ||
       length(over16_state.switch_scan.positions) != switch_count
        resize!(over16_state.switch_scan.positions, switch_count)
        @inbounds for index in eachindex(initial_closed_mask)
            over16_state.switch_scan.positions[index] =
                _deck_time_switch_position(initial_closed_mask[index])
        end
    end
    length(over16_state.switch_scan.elapsed_open_times) == switch_count || begin
        resize!(over16_state.switch_scan.elapsed_open_times, switch_count)
        fill!(over16_state.switch_scan.elapsed_open_times, 0.0)
    end

    over16_state.switch_operation.closed_switch_count =
        count(identity, over16_state.switch_topology.closed_mask)

    if switch_count == 0
        if !isempty(over16_state.switch_admittance.base_admittance)
            over16_state.switch_admittance = OVER16SwitchAdmittanceState(
                zeros(Float64, 0, 0),
            )
        end
        if !isempty(over16_state.switch_retriangularization.factor)
            over16_state.switch_retriangularization =
                OVER16SwitchRetriangularizationState(zeros(Float64, 0, 0))
        end
        isempty(over16_state.switch_sparse_factor.kk) ||
            (over16_state.switch_sparse_factor =
                OVER16SwitchSparseFactorWorkspaceState(Int[], Float64[], Int[]))
        isempty(over16_state.switch_fortran_sparse_factor.kk) ||
            (over16_state.switch_fortran_sparse_factor =
                OVER16FortranSparseFactorWorkspaceState(Int[], Float64[], Int[]))
        return over16_state
    end

    base_admittance_size = size(over16_state.switch_admittance.base_admittance)
    active_admittance_size = size(over16_state.switch_admittance.admittance)
    if base_admittance_size != (node_count, node_count) ||
       active_admittance_size != (node_count, node_count)
        over16_state.switch_admittance =
            OVER16SwitchAdmittanceState(node_count; switch_count = switch_count)
    elseif length(over16_state.switch_admittance.switch_conductances) != switch_count
        resize!(over16_state.switch_admittance.switch_conductances, switch_count)
        fill!(over16_state.switch_admittance.switch_conductances, 0.0)
    end

    size(over16_state.switch_retriangularization.factor, 1) == node_count &&
    size(over16_state.switch_retriangularization.factor, 2) == node_count ||
        (over16_state.switch_retriangularization =
            OVER16SwitchRetriangularizationState(node_count))
    length(over16_state.switch_sparse_factor.kk) == node_count ||
        (over16_state.switch_sparse_factor = OVER16SwitchSparseFactorWorkspaceState(node_count))
    length(over16_state.switch_fortran_sparse_factor.kk) == node_count ||
        (over16_state.switch_fortran_sparse_factor =
            OVER16FortranSparseFactorWorkspaceState(node_count))

    length(over16_state.switch_current.rhs) == node_count || begin
        resize!(over16_state.switch_current.rhs, node_count)
        fill!(over16_state.switch_current.rhs, 0.0)
    end
    length(over16_state.switch_current.network_solution) == node_count || begin
        resize!(over16_state.switch_current.network_solution, node_count)
        fill!(over16_state.switch_current.network_solution, 0.0)
    end
    length(over16_state.switch_current.switch_currents) == switch_count || begin
        resize!(over16_state.switch_current.switch_currents, switch_count)
        fill!(over16_state.switch_current.switch_currents, 0.0)
    end
    length(over16_state.switch_current.current_products) == switch_count || begin
        resize!(over16_state.switch_current.current_products, switch_count)
        fill!(over16_state.switch_current.current_products, 0.0)
    end

    if initialize_deck_switch_state ||
       length(over16_state.switch_post_current.positions) != switch_count
        resize!(over16_state.switch_post_current.positions, switch_count)
        @inbounds for index in eachindex(over16_state.switch_topology.closed_mask)
            over16_state.switch_post_current.positions[index] =
                _deck_time_switch_position(
                    over16_state.switch_topology.closed_mask[index],
                )
        end
    end
    length(over16_state.switch_post_current.switch_currents) == switch_count || begin
        resize!(over16_state.switch_post_current.switch_currents, switch_count)
        fill!(over16_state.switch_post_current.switch_currents, 0.0)
    end
    length(over16_state.switch_post_current.energies) == switch_count || begin
        resize!(over16_state.switch_post_current.energies, switch_count)
        fill!(over16_state.switch_post_current.energies, 0.0)
    end

    return nothing
end

function _deck_configure_current_extinction_step!(
    ::Tuple{},
    _element_names::Vector{Symbol},
    _target_name::Symbol,
    _time_s::Float64,
    _switch_index::Int,
    _requested_closed_mask::Vector{Bool},
    _workspace::DeckTimeSwitchStepWorkspace,
    _element_index::Int = 1,
)
    return nothing
end

function _deck_configure_current_extinction_step!(
    elements::Tuple,
    element_names::Vector{Symbol},
    target_name::Symbol,
    time_s::Float64,
    switch_index::Int,
    requested_closed_mask::Vector{Bool},
    workspace::DeckTimeSwitchStepWorkspace,
    element_index::Int = 1,
)
    element = first(elements)
    if @inbounds element_names[element_index] == target_name
        if element isa CurrentZeroSwitch
            prepare_current_zero_switch!(element, time_s)
            requested_closed_mask[switch_index] = element.closed
            workspace.open_times[switch_index] = element.open_request_time_s
            workspace.critical_currents[switch_index] = element.critical_current_a
            workspace.delay_times[switch_index] = element.open_delay_time_s
        elseif element isa TimeSwitch && element.current_extinction !== nothing
            prepare_current_zero_switch!(element, time_s)
            requested_closed_mask[switch_index] =
                element.current_extinction.closed
            workspace.open_times[switch_index] = element.open_time_s
            workspace.critical_currents[switch_index] =
                element.current_extinction.critical_current_a
            workspace.delay_times[switch_index] =
                element.current_extinction.not_before_time_s
        end
        return nothing
    end
    return _deck_configure_current_extinction_step!(
        Base.tail(elements),
        element_names,
        target_name,
        time_s,
        switch_index,
        requested_closed_mask,
        workspace,
        element_index + 1,
    )
end


function _deck_configure_current_extinction_step!(
    elements::NodalElementSequence,
    element_names::Vector{Symbol},
    target_name::Symbol,
    time_s::Float64,
    switch_index::Int,
    requested_closed_mask::Vector{Bool},
    workspace::DeckTimeSwitchStepWorkspace,
    element_index::Int=1,
)
    for index in element_index:length(elements)
        @inbounds element_names[index] == target_name || continue
        element = elements[index]
        if element isa CurrentZeroSwitch
            prepare_current_zero_switch!(element, time_s)
            requested_closed_mask[switch_index] = element.closed
            workspace.open_times[switch_index] = element.open_request_time_s
            workspace.critical_currents[switch_index] = element.critical_current_a
            workspace.delay_times[switch_index] = element.open_delay_time_s
        elseif element isa TimeSwitch && element.current_extinction !== nothing
            prepare_current_zero_switch!(element, time_s)
            requested_closed_mask[switch_index] =
                element.current_extinction.closed
            workspace.open_times[switch_index] = element.open_time_s
            workspace.critical_currents[switch_index] =
                element.current_extinction.critical_current_a
            workspace.delay_times[switch_index] =
                element.current_extinction.not_before_time_s
        end
        return nothing
    end
    return nothing
end

function _deck_time_switch_sparse_state_config!(
    plan::DeckOVER16BoundaryPlan,
    context::EMTStepContext,
    over16_state::OVER16AcceptedTimestepState,
    workspace::DeckTimeSwitchStepWorkspace,
)
    node_count = context.system.node_count + 1
    _sync_deck_time_switch_over16_state!(
        over16_state,
        plan.switch_initially_closed_flags,
        node_count,
    )
    requested_closed_mask = _deck_time_switch_closed_mask!(
        workspace.requested_closed_mask,
        plan,
        context.t_s,
    )
    fill!(workspace.open_times, Inf)
    fill!(workspace.critical_currents, 0.0)
    fill!(workspace.delay_times, 0.0)
    @inbounds for index in eachindex(requested_closed_mask)
        _deck_configure_current_extinction_step!(
            context.system.elements,
            context.element_names,
            plan.switch_names[index],
            context.t_s,
            index,
            requested_closed_mask,
            workspace,
        )
    end
    previous_closed_mask = workspace.previous_closed_mask
    copyto!(previous_closed_mask, over16_state.switch_topology.closed_mask)
    modswt = _deck_time_switch_operation_queue!(
        workspace.operation_queue,
        previous_closed_mask,
        requested_closed_mask,
    )
    @inbounds for index in eachindex(requested_closed_mask)
        position = _deck_time_switch_position(requested_closed_mask[index])
        over16_state.switch_scan.positions[index] = position
        over16_state.switch_post_current.positions[index] = position
    end
    empty!(over16_state.switch_operation.modswt)
    append!(over16_state.switch_operation.modswt, modswt)
    over16_state.switch_operation.closed_switch_count =
        count(identity, previous_closed_mask)

    shifted_from_nodes = workspace.shifted_from_nodes
    shifted_to_nodes = workspace.shifted_to_nodes
    rebuild_sparse_factor =
        !isempty(modswt) ||
        over16_state.switch_admittance.retriangularization_count == 0 ||
        over16_state.switch_fortran_sparse_factor.workspace_update_count == 0
    fill!(workspace.source_voltage_differences, 0.0)
    fill!(workspace.source_indices, 0)

    return (
        switch_base_admittance_from_step_context = true,
        switch_operation_enabled = !isempty(modswt),
        switch_admittance_config = (
            from_nodes = shifted_from_nodes,
            to_nodes = shifted_to_nodes,
            modswt = modswt,
            partition_boundary = node_count,
            node_count = node_count,
            reference_node = 1,
            closed_conductances = plan.switch_on_conductance_values,
            open_conductances = plan.switch_off_conductance_values,
            request_retriangularization = rebuild_sparse_factor,
        ),
        switch_retriangularization_config = (
            enabled = false,
            require_retriangularization_request = false,
        ),
        switch_fortran_sparse_factor_config = (
            enabled = rebuild_sparse_factor,
            partition_boundary = node_count,
            first_factor_row = 2,
            require_retriangularization_result = false,
        ),
        switch_network_solution_config = (
            require_fortran_sparse_factor_result = false,
            partition_boundary = node_count,
            first_factor_row = 2,
        ),
        switch_current_config = (
            from_nodes = shifted_from_nodes,
            to_nodes = shifted_to_nodes,
            use_network_solution = true,
        ),
        switch_post_current_config = (
            source_voltage_differences = workspace.source_voltage_differences,
            source_indices = workspace.source_indices,
            open_times = workspace.open_times,
            critical_currents = workspace.critical_currents,
            delay_times = workspace.delay_times,
            t = context.t_s,
            dt = context.dt_s,
        ),
        switch_bvalue_enabled = true,
        switch_alteration_config = NamedTuple(),
        sparse_switch_state_flow_enabled = true,
    )
end

function _deck_time_switch_sparse_node_group_config(
    plan::DeckOVER16BoundaryPlan,
    context::EMTStepContext,
    workspace::DeckTimeSwitchStepWorkspace,
)
    closed_mask = workspace.requested_closed_mask
    enabled = any(closed_mask)
    shifted_from_nodes = workspace.shifted_from_nodes
    shifted_to_nodes = workspace.shifted_to_nodes
    node_count = context.system.node_count + 1
    node_group_successors = workspace.node_group_successors
    grouped_from_nodes = workspace.grouped_from_nodes
    grouped_to_nodes = workspace.grouped_to_nodes
    empty!(grouped_from_nodes)
    empty!(grouped_to_nodes)
    if enabled
        ordering = over16_switch_simple_ordering(
            shifted_from_nodes,
            shifted_to_nodes,
            closed_mask,
            node_count;
            node_count = node_count,
            reference_node = 1,
        )
        copyto!(node_group_successors, ordering.kode)
        @inbounds for index in eachindex(closed_mask)
            closed_mask[index] || continue
            push!(grouped_from_nodes, plan.switch_from_node_indices[index])
            push!(grouped_to_nodes, plan.switch_to_node_indices[index])
        end
    else
        fill!(node_group_successors, 0)
    end
    return (
        sparse_node_group_config = (
            enabled = enabled,
            node_group_successors = node_group_successors,
            grouped_switch_from_nodes = grouped_from_nodes,
            grouped_switch_to_nodes = grouped_to_nodes,
            grouped_switch_count = length(grouped_from_nodes),
        ),
    )
end

function _deck_current_zero_sparse_node_group_config(
    plan::DeckOVER16BoundaryPlan,
    context::EMTStepContext,
)
    plan.switch_count == 0 && return NamedTuple()
    closed_mask = falses(plan.switch_count)
    for row in 1:plan.switch_count
        element_index = findfirst(==(plan.switch_names[row]), context.element_names)
        element_index === nothing && throw(ArgumentError(
            "deck switch $(plan.switch_names[row]) is missing from the timestep context",
        ))
        element = context.system.elements[element_index]
        current_extinction_enabled(element) &&
            prepare_current_zero_switch!(element, context.t_s)
        element isa Union{IdealSwitch,TimeSwitch,CurrentZeroSwitch} ||
            throw(ArgumentError("deck switch runtime element has an unsupported type"))
        closed_mask[row] = switch_closed(element, context.t_s)
    end
    any(closed_mask) || return NamedTuple()
    shifted_from_nodes, shifted_to_nodes =
        _deck_time_switch_reference_shifted_nodes(plan)
    node_count = context.system.node_count + 1
    ordering = over16_switch_simple_ordering(
        shifted_from_nodes,
        shifted_to_nodes,
        closed_mask,
        node_count;
        node_count = node_count,
        reference_node = 1,
    )
    return (
        sparse_node_group_config = (
            node_group_successors = ordering.kode,
            grouped_switch_from_nodes = plan.switch_from_node_indices[closed_mask],
            grouped_switch_to_nodes = plan.switch_to_node_indices[closed_mask],
            grouped_switch_count = count(identity, closed_mask),
        ),
    )
end

function _deck_over16_step_config(plan::DeckOVER16BoundaryPlan,
                                  context::EMTStepContext,
                                  over16_state::OVER16AcceptedTimestepState,
                                  nonlinear_current_config,
                                  features::DeckStepConfigFeatures,
                                  time_switch_workspace::DeckTimeSwitchStepWorkspace,
                                  ::Val{P},
                                  ::Val{T},
                                  source_time_offset_s::Float64) where {P,T}
    time_switch_sparse_state_config =
        T ?
        _deck_time_switch_sparse_state_config!(
            plan,
            context,
            over16_state,
            time_switch_workspace,
        ) :
        NamedTuple()
    return merge(
        _deck_over16_output_step_config(plan, context, features),
        _deck_over5a_source_step_config(
            plan,
            context,
            over16_state,
            features,
            source_time_offset_s,
        ),
        P ?
            _deck_over16_post_extrema_step_config(over16_state) :
            NamedTuple(),
        _deck_nonlinear_current_step_config(
            nonlinear_current_config,
            over16_state,
            context,
            allow_pending_fortran_sparse_factor = T,
        ),
        time_switch_sparse_state_config,
        T ?
            _deck_time_switch_sparse_node_group_config(
                plan,
                context,
                time_switch_workspace,
            ) :
            _deck_current_zero_sparse_node_group_config(plan, context),
    )
end

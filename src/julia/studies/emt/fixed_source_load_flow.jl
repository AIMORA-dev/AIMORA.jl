struct FixedSourceLoadFlowResult
    source::Symbol
    outcome::Symbol
    converged::Bool
    iteration_count::Int
    maximum_iterations::Int
    relative_power_tolerance::Float64
    network_topology_kinds::Vector{Symbol}
    initial_switch_names::Vector{Symbol}
    initial_switch_from_node_names::Vector{Symbol}
    initial_switch_to_node_names::Vector{Symbol}
    initial_switch_closed_flags::Vector{Bool}
    initial_switch_conductances::Vector{Float64}
    source_row_indices::Vector{Int}
    source_node_names::Vector{Symbol}
    source_node_indices::Vector{Int}
    source_voltage_phasors::Vector{ComplexF64}
    source_active_powers::Vector{Float64}
    source_reactive_powers::Vector{Float64}
    constraint_kinds::Vector{Symbol}
    constraint_source_row_indices::Vector{Vector{Int}}
    constraint_source_node_names::Vector{Vector{Symbol}}
    constraint_target_active_powers::Vector{Union{Missing,Float64}}
    constraint_target_reactive_powers::Vector{Union{Missing,Float64}}
    constraint_target_voltage_peaks::Vector{Union{Missing,Float64}}
    constraint_target_angles_deg::Vector{Union{Missing,Float64}}
    constraint_active_powers::Vector{Float64}
    constraint_reactive_powers::Vector{Float64}
    constraint_active_power_mismatches::Vector{Float64}
    constraint_reactive_power_mismatches::Vector{Float64}
    node_voltage_phasors::Vector{ComplexF64}
end

struct FixedSourceInitialSwitchBoundary
    name::Symbol
    from_node_name::Symbol
    to_node_name::Symbol
    from_node::Int
    to_node::Int
    initially_closed::Bool
    on_conductance::Float64
    off_conductance::Float64
    control_signal::Union{Nothing,Symbol}
end

function _fixed_source_initial_switch_boundaries(
    parsed::DeckParser.DeckParseResult,
)
    boundaries = FixedSourceInitialSwitchBoundary[]
    for row in DeckParser.deck_over5_switch_rows(parsed)
        push!(boundaries, FixedSourceInitialSwitchBoundary(
            row.name,
            row.from_node,
            row.to_node,
            Int(row.from_node_value),
            Int(row.to_node_value),
            row.initially_closed,
            Float64(row.on_conductance),
            Float64(row.off_conductance),
            nothing,
        ))
    end
    for (index, row) in
        enumerate(DeckParser.deck_control_system_switch_coupling_rows(parsed))
        row.from_node_index !== missing && row.to_node_index !== missing ||
            throw(ArgumentError(
                "FIX SOURCE controlled switch on line $(row.line_no) has unresolved nodes",
            ))
        push!(boundaries, FixedSourceInitialSwitchBoundary(
            _control_system_switch_element_name(row, index),
            row.from_node,
            row.to_node,
            Int(row.from_node_index),
            Int(row.to_node_index),
            _control_system_switch_initially_closed(row.initial_state),
            1.0e9,
            0.0,
            row.control_signal,
        ))
    end
    return boundaries
end

function _fixed_source_control(parsed::DeckParser.DeckParseResult)
    rows = DeckParser.deck_fixed_source_control_rows(parsed)
    length(rows) == 1 || throw(ArgumentError(
        "FIX SOURCE load flow requires exactly one trailing control card",
    ))
    return only(rows)
end

function _fixed_source_supported_topology!(parsed::DeckParser.DeckParseResult)
    definitions = DeckParser.deck_universal_machine_definition_rows(parsed)
    card1_rows = [row for row in definitions if row.card_index == 1]
    all(row -> row.machine_type in 1:12, card1_rows) || throw(ArgumentError(
        "FIX SOURCE load flow requires an accepted universal-machine type",
    ))
    return parsed
end

function _fixed_source_network_topology_kinds(parsed::DeckParser.DeckParseResult)
    kinds = Symbol[]
    isempty(DeckParser.deck_over2_branch_rows(parsed)) || push!(kinds, :scalar_branch)
    isempty(DeckParser.deck_coupled_lumped_sequence_impedances(parsed)) ||
        push!(kinds, :coupled_lumped_sequence_impedance)
    isempty(DeckParser.deck_coupled_lumped_phase_pi_sections(parsed)) ||
        push!(kinds, :coupled_lumped_phase_pi_section)
    isempty(DeckParser.deck_cascaded_phase_pi_equivalents(parsed)) ||
        push!(kinds, :cascaded_phase_pi_equivalent)
    isempty(DeckParser.deck_over5_switch_rows(parsed)) ||
        push!(kinds, :initial_switch_topology)
    isempty(DeckParser.deck_control_system_switch_rows(parsed)) ||
        push!(kinds, :controlled_switch_initial_topology)
    isempty(DeckParser.deck_generator_equivalent_rows(parsed)) ||
        push!(kinds, :frequency_dependent_generator_equivalent)
    isempty(DeckParser.deck_distributed_transposed_line_modal_branch_states(parsed)) ||
        push!(kinds, :distributed_transposed_line)
    isempty(DeckParser.deck_semlyen_line_groups(parsed)) ||
        push!(kinds, :semlyen_frequency_dependent_line)
    isempty(DeckParser.deck_rational_frequency_line_groups(parsed)) ||
        push!(kinds, :rational_frequency_dependent_line)
    any(
        row -> row.card_index == 1 && row.machine_type in (3, 4, 5, 6, 7),
        DeckParser.deck_universal_machine_definition_rows(parsed),
    ) && push!(kinds, :induction_machine_steady_state_equivalent)
    any(
        row -> row.card_index == 1 && row.machine_type == 1,
        DeckParser.deck_universal_machine_definition_rows(parsed),
    ) && push!(kinds, :wound_field_synchronous_machine_voltage_boundary)
    any(
        row -> row.card_index == 1 && row.machine_type == 2,
        DeckParser.deck_universal_machine_definition_rows(parsed),
    ) && push!(kinds, :two_phase_synchronous_machine_voltage_boundary)
    any(
        row -> row.card_index == 1 && row.machine_type == 8,
        DeckParser.deck_universal_machine_definition_rows(parsed),
    ) && push!(kinds, :separately_excited_dc_machine_voltage_boundary)
    any(
        row -> row.card_index == 1 && row.machine_type in (9, 10, 11, 12),
        DeckParser.deck_universal_machine_definition_rows(parsed),
    ) && push!(kinds, :manual_direct_current_machine_voltage_boundary)
    return kinds
end

struct FixedSourceMachineVoltageBoundary{N}
    machine_index::Int
    machine_type::Int
    terminal_nodes::NTuple{N,Int}
    terminal_voltage_phasors::NTuple{N,ComplexF64}
end

function _fixed_source_prescribed_machine_voltage_boundary(
    parsed::DeckParser.DeckParseResult,
    machine_index::Int,
    card1,
)
    card4 = _deck_universal_machine_definition(parsed, machine_index, 4)
    card4.value1 === missing && throw(ArgumentError(
        "type-$(card1.machine_type) FIX SOURCE load flow requires a requested terminal voltage",
    ))
    card1.machine_type in (1, 2) && card4.value2 === missing && throw(ArgumentError(
        "type-$(card1.machine_type) FIX SOURCE load flow requires a requested terminal angle",
    ))
    machine_branches = sort!(
        [
            row for row in DeckParser.deck_universal_machine_generated_branch_rows(parsed)
            if row.machine_index == machine_index
        ];
        by = row -> row.branch_index,
    )
    expected_branch_count = card1.machine_type == 8 ? 2 : 3
    length(machine_branches) == expected_branch_count || throw(ArgumentError(
        "type-$(card1.machine_type) FIX SOURCE load flow requires $expected_branch_count generated coil branches",
    ))
    voltage_peak = Float64(card4.value1)
    voltage_peak > 0.0 || throw(ArgumentError(
        "type-$(card1.machine_type) FIX SOURCE requested terminal voltage must be positive",
    ))
    base_angle_deg = card4.value2 === missing ? 0.0 : Float64(card4.value2)
    branch_positions, phase_offsets_deg = if card1.machine_type == 1
        ((1, 2, 3), (0.0, -120.0, 120.0))
    elseif card1.machine_type == 2
        ((2, 3), (0.0, 90.0))
    else
        ((1,), (0.0,))
    end
    terminal_nodes = ntuple(length(branch_positions)) do position
        branch = machine_branches[branch_positions[position]]
        branch.to_node_value > 0 || throw(ArgumentError(
            "type-$(card1.machine_type) FIX SOURCE load flow requires physical power terminals",
        ))
        Int(branch.to_node_value)
    end
    terminal_phasors = ntuple(length(branch_positions)) do position
        voltage_peak * cis(deg2rad(base_angle_deg + phase_offsets_deg[position]))
    end
    return FixedSourceMachineVoltageBoundary(
        machine_index,
        card1.machine_type,
        terminal_nodes,
        terminal_phasors,
    )
end

function _fixed_source_manual_direct_machine_voltage_boundary(
    parsed::DeckParser.DeckParseResult,
    machine_index::Int,
)
    initial = _manual_direct_machine_initial_voltage_boundary(parsed, machine_index)
    active_positions = findall(>(0), initial.power_terminal_nodes)
    active_positions == [3] || throw(ArgumentError(
        "type-$(initial.machine_type) FIX SOURCE load flow requires one physical armature terminal in slot 3",
    ))
    position = only(active_positions)
    return FixedSourceMachineVoltageBoundary(
        machine_index,
        initial.machine_type,
        (initial.power_terminal_nodes[position],),
        (ComplexF64(initial.power_terminal_voltages[position]),),
    )
end

function _fixed_source_machine_voltage_boundary(
    parsed::DeckParser.DeckParseResult,
    machine_index::Int,
)
    card1 = _deck_universal_machine_definition(parsed, machine_index, 1)
    if card1.machine_type in (1, 2) ||
       (card1.machine_type == 8 &&
        _deck_universal_machine_initialization_mode(parsed) == :automatic)
        return _fixed_source_prescribed_machine_voltage_boundary(
            parsed,
            machine_index,
            card1,
        )
    elseif card1.machine_type in (8, 9, 10, 11, 12)
        return _fixed_source_manual_direct_machine_voltage_boundary(parsed, machine_index)
    end
    return nothing
end

function _fixed_source_machine_voltage_boundaries(
    parsed::DeckParser.DeckParseResult,
)
    definitions = DeckParser.deck_universal_machine_definition_rows(parsed)
    nodes = Int[]
    phasors = ComplexF64[]
    for machine_index in sort(unique(row.machine_index for row in definitions))
        boundary = _fixed_source_machine_voltage_boundary(parsed, machine_index)
        boundary === nothing && continue
        append!(nodes, boundary.terminal_nodes)
        append!(phasors, boundary.terminal_voltage_phasors)
    end
    allunique(nodes) || throw(ArgumentError(
        "FIX SOURCE machine voltage boundaries must own distinct terminals",
    ))
    return nodes, phasors
end

function _fixed_source_rows(parsed::DeckParser.DeckParseResult)
    all_rows = DeckParser.deck_over5a_source_rows(parsed)
    isempty(all_rows) && throw(ArgumentError(
        "FIX SOURCE load flow requires constant or sinusoidal voltage sources",
    ))
    source_row_indices = findall(
        row -> row.iform in (11, 14) && row.node_value > 0,
        all_rows,
    )
    isempty(source_row_indices) && throw(ArgumentError(
        "FIX SOURCE load flow requires positive-node constant or sinusoidal voltage sources",
    ))
    rows = all_rows[source_row_indices]
    ignored_rows = setdiff(collect(eachindex(all_rows)), source_row_indices)
    machine_definitions = DeckParser.deck_universal_machine_definition_rows(parsed)
    all(
        index -> !isempty(machine_definitions) &&
                 all_rows[index].iform in (11, 14) &&
                 all_rows[index].node_value < 0,
        ignored_rows,
    ) || throw(ArgumentError(
        "FIX SOURCE load flow found a non-machine source outside its positive-node electrical source set",
    ))
    nodes = Int[row.node_value for row in rows]
    allunique(nodes) || throw(ArgumentError(
        "FIX SOURCE load flow requires one voltage source per source node",
    ))
    source_types = unique(row.iform for row in rows)
    length(source_types) == 1 || throw(ArgumentError(
        "FIX SOURCE load flow cannot mix constant and sinusoidal source domains",
    ))
    if only(source_types) == 14
        frequencies = Float64[row.sfreq for row in rows]
        all(frequency -> isapprox(frequency, first(frequencies); atol = 1.0e-12), frequencies) ||
            throw(ArgumentError("FIX SOURCE load flow requires one steady-state source frequency"))
    end
    return rows, source_row_indices
end

function _fixed_source_constant_source_domain(parsed::DeckParser.DeckParseResult)
    rows, _ = _fixed_source_rows(parsed)
    return only(unique(row.iform for row in rows)) == 11
end

function _stamp_fixed_source_dc_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
)
    for row in DeckParser.deck_over2_branch_rows(parsed)
        admittance = if row.branch_kind == :conductance
            complex(Float64(row.conductance), 0.0)
        elseif row.branch_kind == :series_rl && row.resistance > 0.0
            complex(inv(Float64(row.resistance)), 0.0)
        else
            continue
        end
        _stamp_complex_branch_admittance!(
            matrix,
            row.from_node_value,
            row.to_node_value,
            admittance,
        )
    end
    for row in DeckParser.deck_universal_machine_speed_capacitor_rows(parsed)
        row.resistance > 0.0 && row.capacitance == 0.0 || continue
        _stamp_complex_branch_admittance!(
            matrix,
            row.capacitor_node_value,
            row.mass_node_value,
            complex(inv(Float64(row.resistance)), 0.0),
        )
    end
    return matrix
end

function _fixed_source_network_admittance(parsed::DeckParser.DeckParseResult)
    node_count = maximum(values(parsed.node_map); init = 0)
    matrix = zeros(ComplexF64, node_count, node_count)
    frequency_partition = DeckParser.deck_steady_state_frequency_partition(parsed)
    isempty(frequency_partition.unsupported_topology_kinds) || throw(ArgumentError(
        "FIX SOURCE load flow found unsupported steady-state topology",
    ))
    constant_source_domain = _fixed_source_constant_source_domain(parsed)
    if constant_source_domain
        isempty(DeckParser.deck_coupled_lumped_sequence_impedances(parsed)) &&
        isempty(DeckParser.deck_coupled_lumped_phase_pi_sections(parsed)) &&
        isempty(DeckParser.deck_cascaded_phase_pi_equivalents(parsed)) &&
        isempty(DeckParser.deck_generator_equivalent_rows(parsed)) &&
        isempty(DeckParser.deck_distributed_transposed_line_modal_branch_states(parsed)) &&
        isempty(DeckParser.deck_semlyen_line_groups(parsed)) &&
        isempty(DeckParser.deck_rational_frequency_line_groups(parsed)) &&
        !any(
            row -> row.card_index == 1 && row.machine_type in (3, 4, 5, 6, 7),
            DeckParser.deck_universal_machine_definition_rows(parsed),
        ) || throw(ArgumentError(
            "constant-source FIX SOURCE load flow requires scalar DC network topology",
        ))
        _stamp_fixed_source_dc_admittance!(matrix, parsed)
    else
        _stamp_deck_branch_steady_state_admittance!(matrix, parsed, frequency_partition)
    end
    _stamp_deck_open_switch_steady_state_admittance!(matrix, parsed)
    for boundary in _fixed_source_initial_switch_boundaries(parsed)
        boundary.control_signal === nothing && continue
        if boundary.initially_closed
            boundary.from_node != 0 && boundary.to_node != 0 && continue
            conductance = boundary.on_conductance
        else
            conductance = boundary.off_conductance
        end
        _stamp_complex_branch_admittance!(
            matrix,
            boundary.from_node,
            boundary.to_node,
            complex(conductance, 0.0),
        )
    end
    if !constant_source_domain
        _stamp_coupled_lumped_sequence_steady_state_admittance!(matrix, parsed)
        _stamp_generator_equivalent_steady_state_admittance!(matrix, parsed)
        _stamp_coupled_lumped_phase_pi_steady_state_admittance!(matrix, parsed)
        _stamp_cascaded_phase_pi_steady_state_admittance!(matrix, parsed)
        _stamp_distributed_line_steady_state_admittance!(matrix, parsed)
        _stamp_semlyen_line_steady_state_admittance!(matrix, parsed)
        _stamp_deck_induction_machine_steady_state_admittance!(matrix, parsed)
    end
    return matrix
end

function _fixed_source_node_representatives(
    parsed::DeckParser.DeckParseResult,
    node_count::Int,
)
    representatives = _steady_state_closed_switch_representatives(parsed, node_count)

    function representative(node::Int)
        root = node
        while representatives[root] != root
            root = representatives[root]
        end
        while representatives[node] != node
            parent = representatives[node]
            representatives[node] = root
            node = parent
        end
        return root
    end

    function union_nodes(left::Int, right::Int)
        (left == 0 || right == 0) && throw(ArgumentError(
            "constant-source FIX SOURCE does not yet support a grounded ideal inductor",
        ))
        left_root = representative(left)
        right_root = representative(right)
        left_root == right_root && return nothing
        retained = max(left_root, right_root)
        replaced = min(left_root, right_root)
        representatives[replaced] = retained
        return nothing
    end

    for boundary in _fixed_source_initial_switch_boundaries(parsed)
        boundary.control_signal === nothing && continue
        boundary.initially_closed || continue
        boundary.from_node == 0 || boundary.to_node == 0 ||
            union_nodes(boundary.from_node, boundary.to_node)
    end
    for node in eachindex(representatives)
        representatives[node] = representative(node)
    end
    _fixed_source_constant_source_domain(parsed) || return representatives

    for row in DeckParser.deck_over2_branch_rows(parsed)
        row.branch_kind == :series_rl && iszero(row.resistance) || continue
        union_nodes(Int(row.from_node_value), Int(row.to_node_value))
    end
    for row in DeckParser.deck_universal_machine_generated_branch_rows(parsed)
        row.reactance === missing && continue
        union_nodes(Int(row.from_node_value), Int(row.to_node_value))
    end
    for node in eachindex(representatives)
        representatives[node] = representative(node)
    end
    return representatives
end

function _deck_fixed_source_synchronous_steady_state(
    parsed::DeckParser.DeckParseResult,
)
    result = deck_fixed_source_load_flow(parsed)
    result.converged || throw(ArgumentError(
        "synchronous-machine initialization requires a converged FIX SOURCE solution",
    ))
    admittance = _fixed_source_network_admittance(parsed)
    frequency_partition = DeckParser.deck_steady_state_frequency_partition(parsed)
    phasors = copy(result.node_voltage_phasors)
    load_flow_currents = admittance * phasors
    for branch in DeckParser.deck_universal_machine_generated_branch_rows(parsed)
        card1 = _deck_universal_machine_definition(parsed, branch.machine_index, 1)
        card1.machine_type in (1, 2) || continue
        branch.reactance === missing && throw(ArgumentError(
            "synchronous-machine power-coil leakage is missing",
        ))
        physical_node = Int(branch.to_node_value)
        internal_node = Int(branch.from_node_value)
        frequency_hz = DeckParser.deck_node_steady_state_frequency_hz(
            frequency_partition,
            physical_node,
            internal_node,
            _deck_steady_state_frequency_hz(parsed),
        )
        leakage_inductance_h = DeckParser.fixed_card_branch_timestep_inductance(
            parsed,
            Float64(branch.reactance),
        )
        machine_terminal_current = -load_flow_currents[physical_node]
        phasors[internal_node] = phasors[physical_node] -
            im * (2.0 * pi * frequency_hz * leakage_inductance_h) *
            machine_terminal_current
    end
    node_current_phasors = admittance * phasors
    return (
        source = :fixed_source_synchronous_machine_steady_state,
        outcome = :steady_state_initial_voltage_sample,
        steady_state_frequency_hz = _deck_steady_state_frequency_hz(parsed),
        node_steady_state_frequencies_hz =
            copy(frequency_partition.node_frequencies_hz),
        timestep_s = DeckParser.deck_fixed_time_horizon_options(parsed).dt_s,
        node_map = Dict{Symbol,Int}(parsed.node_map),
        node_names = ordered_node_names(parsed.node_map),
        node_voltage_phasors = phasors,
        node_voltage_values = real.(phasors),
        node_current_phasors = node_current_phasors,
        fixed_source_load_flow = result,
    )
end

mutable struct FixedSourceVoltageWorkspace{F}
    factor::F
    unknown_admittance::Matrix{ComplexF64}
    boundary_admittance::Matrix{ComplexF64}
    fixed_indices::Vector{Int}
    unknown_indices::Vector{Int}
    node_reduced_indices::Vector{Int}
    fixed_phasors::Vector{ComplexF64}
    reduced_voltage::Vector{ComplexF64}
    unknown_rhs::Vector{ComplexF64}
    unknown_solution::Vector{ComplexF64}
    residual::Vector{ComplexF64}
    node_voltage::Vector{ComplexF64}
    check_residual::Bool
end

function _fixed_source_factor(admittance::Matrix{ComplexF64})
    isempty(admittance) && return nothing, false
    try
        return lu(admittance), false
    catch error
        error isa SingularException || rethrow()
    end
    return qr(admittance, ColumnNorm()), true
end

function _fixed_source_voltage_workspace(
    admittance::Matrix{ComplexF64},
    fixed_nodes::Vector{Int},
    node_representatives::Vector{Int},
)
    node_count = size(admittance, 1)
    length(node_representatives) == node_count || throw(ArgumentError(
        "FIX SOURCE switch representatives must cover every network node",
    ))
    fixed_representatives = node_representatives[fixed_nodes]
    allunique(fixed_representatives) || throw(ArgumentError(
        "FIX SOURCE does not support multiple voltage sources in one closed-switch group",
    ))

    active_representatives = sort(unique(node_representatives))
    reduced_index = zeros(Int, maximum(active_representatives; init = 0))
    for (index, representative) in enumerate(active_representatives)
        reduced_index[representative] = index
    end
    node_reduced_indices = Int[
        reduced_index[node_representatives[node]] for node in 1:node_count
    ]
    reduced_admittance = zeros(
        ComplexF64,
        length(active_representatives),
        length(active_representatives),
    )
    for row in 1:node_count, column in 1:node_count
        reduced_admittance[
            node_reduced_indices[row],
            node_reduced_indices[column],
        ] += admittance[row, column]
    end

    fixed_indices = Int[reduced_index[representative] for representative in fixed_representatives]
    unknown_indices = setdiff(collect(eachindex(active_representatives)), fixed_indices)
    unknown_admittance = reduced_admittance[unknown_indices, unknown_indices]
    boundary_admittance = reduced_admittance[unknown_indices, fixed_indices]
    factor, check_residual = _fixed_source_factor(unknown_admittance)
    unknown_count = length(unknown_indices)
    return FixedSourceVoltageWorkspace(
        factor,
        unknown_admittance,
        boundary_admittance,
        fixed_indices,
        unknown_indices,
        node_reduced_indices,
        zeros(ComplexF64, length(fixed_indices)),
        zeros(ComplexF64, length(active_representatives)),
        zeros(ComplexF64, unknown_count),
        zeros(ComplexF64, unknown_count),
        zeros(ComplexF64, unknown_count),
        zeros(ComplexF64, node_count),
        check_residual,
    )
end

function _fixed_source_voltage_solution!(
    workspace::FixedSourceVoltageWorkspace,
    source_phasors::Vector{ComplexF64},
    machine_phasors::Vector{ComplexF64},
)
    fixed_count = length(workspace.fixed_phasors)
    length(source_phasors) + length(machine_phasors) == fixed_count ||
        throw(ArgumentError("FIX SOURCE phasor count must match fixed network nodes"))
    copyto!(workspace.fixed_phasors, 1, source_phasors, 1, length(source_phasors))
    copyto!(
        workspace.fixed_phasors,
        length(source_phasors) + 1,
        machine_phasors,
        1,
        length(machine_phasors),
    )
    @inbounds for index in eachindex(workspace.fixed_indices)
        workspace.reduced_voltage[workspace.fixed_indices[index]] =
            workspace.fixed_phasors[index]
    end
    if !isempty(workspace.unknown_indices)
        mul!(
            workspace.unknown_rhs,
            workspace.boundary_admittance,
            workspace.fixed_phasors,
            -1.0,
            0.0,
        )
        copyto!(workspace.unknown_solution, workspace.unknown_rhs)
        ldiv!(workspace.factor, workspace.unknown_solution)
        if workspace.check_residual
            mul!(
                workspace.residual,
                workspace.unknown_admittance,
                workspace.unknown_solution,
            )
            @inbounds for index in eachindex(workspace.residual)
                workspace.residual[index] -= workspace.unknown_rhs[index]
            end
            residual_scale = max(norm(workspace.unknown_rhs, Inf), 1.0)
            norm(workspace.residual, Inf) <= 1.0e-10 * residual_scale ||
                throw(ArgumentError(
                    "rank-deficient steady-state network equations are inconsistent",
                ))
        end
        @inbounds for index in eachindex(workspace.unknown_indices)
            workspace.reduced_voltage[workspace.unknown_indices[index]] =
                workspace.unknown_solution[index]
        end
    end
    @inbounds for node in eachindex(workspace.node_voltage)
        workspace.node_voltage[node] =
            workspace.reduced_voltage[workspace.node_reduced_indices[node]]
    end
    return workspace.node_voltage
end

function _fixed_source_indexed_sum(values::Vector{Float64}, indices::Vector{Int})
    total = 0.0
    @inbounds for index in indices
        total += values[index]
    end
    return total
end

function _fixed_source_constraint_source_rows(constraint, source_rows)
    indices = Int[]
    for node_name in constraint.source_node_names
        isempty(String(node_name)) && continue
        index = findfirst(row -> row.node == node_name, source_rows)
        index === nothing && throw(ArgumentError(
            "FIX SOURCE node $(node_name) does not name an owned voltage source",
        ))
        push!(indices, index)
    end
    isempty(indices) && throw(ArgumentError("FIX SOURCE constraint has no source nodes"))
    return indices
end

function _fixed_source_bounded_angle(
    angle_deg::Float64,
    correction_deg::Float64,
    minimum_deg::Float64,
    maximum_deg::Float64,
)
    correction = correction_deg
    while !(minimum_deg <= angle_deg - correction <= maximum_deg) &&
          abs(correction) > eps(Float64)
        correction *= 0.1
    end
    return minimum_deg <= angle_deg - correction <= maximum_deg ?
        angle_deg - correction : angle_deg
end

function deck_fixed_source_load_flow(parsed::DeckParser.DeckParseResult)
    DeckParser.assert_deck_valid!(parsed)
    constraints = DeckParser.deck_fixed_source_constraint_rows(parsed)
    isempty(constraints) && throw(ArgumentError("deck does not request FIX SOURCE load flow"))
    _fixed_source_supported_topology!(parsed)
    control = _fixed_source_control(parsed)
    source_rows, source_row_indices = _fixed_source_rows(parsed)
    source_nodes = Int[row.node_value for row in source_rows]
    source_phasors = ComplexF64[_deck_source_voltage_phasor(row) for row in source_rows]
    machine_voltage_nodes, machine_voltage_phasors =
        _fixed_source_machine_voltage_boundaries(parsed)
    fixed_nodes = vcat(source_nodes, machine_voltage_nodes)
    allunique(fixed_nodes) || throw(ArgumentError(
        "FIX SOURCE deck sources and machine voltage boundaries must own distinct nodes",
    ))
    switch_boundaries = _fixed_source_initial_switch_boundaries(parsed)
    constraint_source_rows = [
        _fixed_source_constraint_source_rows(constraint, source_rows)
        for constraint in constraints
    ]
    for (constraint, rows) in zip(constraints, constraint_source_rows)
        any(source_rows[index].iform == 11 for index in rows) &&
        constraint.constraint_kind != :active_power_voltage && throw(ArgumentError(
            "constant-source FIX SOURCE constraints require active power and voltage",
        ))
    end
    admittance = _fixed_source_network_admittance(parsed)
    node_representatives = _fixed_source_node_representatives(
        parsed,
        size(admittance, 1),
    )
    voltage_workspace = _fixed_source_voltage_workspace(
        admittance,
        fixed_nodes,
        node_representatives,
    )
    power_targets = Float64[]
    for constraint in constraints
        constraint.active_power === missing || push!(power_targets, abs(constraint.active_power))
        constraint.reactive_power === missing || push!(power_targets, abs(constraint.reactive_power))
    end
    power_scale = maximum(power_targets; init = 0.0)
    power_scale = max(power_scale, DeckParser.deck_fixed_time_horizon_options(parsed).tolerance)

    node_voltage_phasors = voltage_workspace.node_voltage
    node_currents = similar(node_voltage_phasors)
    source_active_powers = zeros(Float64, length(source_rows))
    source_reactive_powers = zeros(Float64, length(source_rows))
    constraint_active_powers = zeros(Float64, length(constraints))
    constraint_reactive_powers = zeros(Float64, length(constraints))
    active_mismatches = zeros(Float64, length(constraints))
    reactive_mismatches = zeros(Float64, length(constraints))
    converged = false
    completed_iterations = 0

    for iteration in 0:control.maximum_iterations
        _fixed_source_voltage_solution!(
            voltage_workspace,
            source_phasors,
            machine_voltage_phasors,
        )
        mul!(node_currents, admittance, node_voltage_phasors)
        for (source_index, node) in enumerate(source_nodes)
            representative = node_representatives[node]
            group_current = zero(ComplexF64)
            @inbounds for group_node in eachindex(node_representatives)
                node_representatives[group_node] == representative || continue
                group_current += node_currents[group_node]
            end
            power_factor = source_rows[source_index].iform == 11 ? 1.0 : 0.5
            source_power = power_factor * node_voltage_phasors[node] * conj(group_current)
            source_active_powers[source_index] = real(source_power)
            source_reactive_powers[source_index] = imag(source_power)
        end
        converged = true
        for (constraint_index, constraint) in enumerate(constraints)
            rows = constraint_source_rows[constraint_index]
            active_power = _fixed_source_indexed_sum(source_active_powers, rows)
            reactive_power = _fixed_source_indexed_sum(source_reactive_powers, rows)
            constraint_active_powers[constraint_index] = active_power
            constraint_reactive_powers[constraint_index] = reactive_power
            active_mismatch = constraint.constraint_kind == :angle_reactive_power ? 0.0 :
                (active_power - Float64(constraint.active_power)) / power_scale
            reactive_mismatch = constraint.constraint_kind == :active_power_voltage ? 0.0 :
                (reactive_power - Float64(constraint.reactive_power)) / power_scale
            active_mismatches[constraint_index] = active_mismatch
            reactive_mismatches[constraint_index] = reactive_mismatch
            converged &= abs(active_mismatch) <= control.relative_power_tolerance &&
                         abs(reactive_mismatch) <= control.relative_power_tolerance
        end
        completed_iterations = iteration
        (converged || iteration == control.maximum_iterations) && break

        iteration_scale = ((control.maximum_iterations - iteration) /
                           control.maximum_iterations)^2
        for (constraint_index, constraint) in enumerate(constraints)
            angle_correction_deg = clamp(
                active_mismatches[constraint_index] *
                iteration_scale * control.angle_correction_factor,
                -1.0,
                1.0,
            )
            voltage_correction = clamp(
                reactive_mismatches[constraint_index] *
                iteration_scale * control.voltage_correction_factor,
                -0.01,
                0.01,
            )
            for (phase_index, source_index) in
                enumerate(constraint_source_rows[constraint_index])
                source_phasor = source_phasors[source_index]
                constant_source = source_rows[source_index].iform == 11
                angle_deg = rad2deg(angle(source_phasor))
                magnitude = abs(source_phasor)
                if constant_source
                    angle_deg = 0.0
                elseif constraint.constraint_kind == :angle_reactive_power
                    angle_deg = Float64(constraint.angle_deg) + (phase_index == 2 ? -120.0 :
                        phase_index == 3 ? 120.0 : 0.0)
                else
                    angle_deg = _fixed_source_bounded_angle(
                        angle_deg,
                        angle_correction_deg,
                        constraint.minimum_angle_deg,
                        constraint.maximum_angle_deg,
                    )
                end
                if constraint.constraint_kind == :active_power_voltage
                    magnitude = Float64(constraint.voltage_peak)
                else
                    magnitude = clamp(
                        magnitude * (1.0 - voltage_correction),
                        constraint.minimum_voltage,
                        constraint.maximum_voltage,
                    )
                end
                source_phasors[source_index] = magnitude * cis(deg2rad(angle_deg))
            end
        end
    end

    return FixedSourceLoadFlowResult(
        :deck_fixed_source_load_flow,
        :steady_state_source_constraint_solution,
        converged,
        completed_iterations,
        control.maximum_iterations,
        control.relative_power_tolerance,
        _fixed_source_network_topology_kinds(parsed),
        Symbol[boundary.name for boundary in switch_boundaries],
        Symbol[boundary.from_node_name for boundary in switch_boundaries],
        Symbol[boundary.to_node_name for boundary in switch_boundaries],
        Bool[boundary.initially_closed for boundary in switch_boundaries],
        Float64[
            boundary.initially_closed ? boundary.on_conductance : boundary.off_conductance
            for boundary in switch_boundaries
        ],
        source_row_indices,
        Symbol[row.node for row in source_rows],
        source_nodes,
        source_phasors,
        source_active_powers,
        source_reactive_powers,
        Symbol[constraint.constraint_kind for constraint in constraints],
        [copy(indices) for indices in constraint_source_rows],
        [Symbol[source_rows[index].node for index in indices] for indices in constraint_source_rows],
        Union{Missing,Float64}[constraint.active_power for constraint in constraints],
        Union{Missing,Float64}[constraint.reactive_power for constraint in constraints],
        Union{Missing,Float64}[constraint.voltage_peak for constraint in constraints],
        Union{Missing,Float64}[constraint.angle_deg for constraint in constraints],
        constraint_active_powers,
        constraint_reactive_powers,
        active_mismatches,
        reactive_mismatches,
        node_voltage_phasors,
    )
end

function _fixed_source_runtime_row(row, phasor::ComplexF64)
    crest = row.iform == 11 ? real(phasor) : abs(phasor)
    phase = row.iform == 11 ? row.time1 : angle(phasor)
    return DeckParser.DeckOVER5ASourceRow(
        row.name,
        row.node,
        row.node_value,
        row.line_no,
        row.layout_kind,
        row.iform,
        crest,
        phase,
        row.time2,
        row.sfreq,
        row.tstart,
        row.tstop,
    )
end

function apply_deck_fixed_source_load_flow(parsed::DeckParser.DeckParseResult)
    result = deck_fixed_source_load_flow(parsed)
    runtime_deck = deepcopy(parsed)
    source_rows = runtime_deck.over5a_source_rows
    constrained_source_indices = Set(Iterators.flatten(result.constraint_source_row_indices))
    for (source_index, source_row_index) in enumerate(result.source_row_indices)
        source_index in constrained_source_indices || continue
        row = source_rows[source_row_index]
        phasor = result.source_voltage_phasors[source_index]
        runtime_row = _fixed_source_runtime_row(row, phasor)
        source_rows[source_row_index] = runtime_row
        element_indices = findall(runtime_deck.elements) do element
            element isa TheveninSource && element.node == row.node_value
        end
        length(element_indices) == 1 || throw(ArgumentError(
            "FIX SOURCE runtime requires one live voltage-source element per source node",
        ))
        element_index = only(element_indices)
        element = runtime_deck.elements[element_index]
        runtime_deck.elements[element_index] = TheveninSource(
            element.node,
            element.g,
            AnalyticSourceSignal(
                runtime_row.iform,
                runtime_row.crest,
                runtime_row.time1,
                runtime_row.sfreq,
                runtime_row.tstart,
                runtime_row.tstop,
            ),
        )
    end
    return (deck = runtime_deck, load_flow = result)
end

_fixed_source_report_value(::Missing) = ""
_fixed_source_report_value(value::Float64) = @sprintf("%.12g", value)

function write_fixed_source_load_flow_report(io::IO, result::FixedSourceLoadFlowResult)
    println(io, "FIX SOURCE LOAD FLOW")
    println(io, "converged = $(result.converged)")
    println(io, "iterations = $(result.iteration_count)")
    println(io, "network_topologies = $(join(String.(result.network_topology_kinds), '+'))")
    println(io, "switch,name,from_node,to_node,initial_state,conductance_siemens")
    for index in eachindex(result.initial_switch_names)
        println(io, join((
            string(index),
            String(result.initial_switch_names[index]),
            String(result.initial_switch_from_node_names[index]),
            String(result.initial_switch_to_node_names[index]),
            result.initial_switch_closed_flags[index] ? "closed" : "open",
            _fixed_source_report_value(result.initial_switch_conductances[index]),
        ), ','))
    end
    println(
        io,
        "constraint,kind,nodes,target_active_power,target_reactive_power," *
        "target_voltage_peak,target_angle_deg,active_power,reactive_power",
    )
    for index in eachindex(result.constraint_kinds)
        println(io, join((
            string(index),
            String(result.constraint_kinds[index]),
            join(String.(result.constraint_source_node_names[index]), '+'),
            _fixed_source_report_value(result.constraint_target_active_powers[index]),
            _fixed_source_report_value(result.constraint_target_reactive_powers[index]),
            _fixed_source_report_value(result.constraint_target_voltage_peaks[index]),
            _fixed_source_report_value(result.constraint_target_angles_deg[index]),
            _fixed_source_report_value(result.constraint_active_powers[index]),
            _fixed_source_report_value(result.constraint_reactive_powers[index]),
        ), ','))
    end
    println(io, "node,voltage_peak,angle_deg,active_power,reactive_power")
    for index in eachindex(result.source_node_names)
        @printf(
            io,
            "%s,%.12g,%.12g,%.12g,%.12g\n",
            String(result.source_node_names[index]),
            abs(result.source_voltage_phasors[index]),
            rad2deg(angle(result.source_voltage_phasors[index])),
            result.source_active_powers[index],
            result.source_reactive_powers[index],
        )
    end
    return io
end

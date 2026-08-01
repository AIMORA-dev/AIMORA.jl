export DeckSteadyStateFrequencyPartition,
       DeckMixedSteadyStateFrequencyConflict,
       deck_steady_state_frequency_partition,
       deck_node_steady_state_frequency_hz

struct DeckSteadyStateFrequencyPartition
    node_source_row_indices::Vector{Int}
    node_frequencies_hz::Vector{Float64}
    active_source_row_indices::Vector{Int}
    source_frequencies_hz::Vector{Union{Nothing,Float64}}
    source_successor_indices::Vector{Int}
    source_groups::Vector{Vector{Int}}
    subnetwork_frequencies_hz::Vector{Float64}
    subnetwork_node_indices::Vector{Vector{Int}}
    inactive_node_indices::Vector{Int}
    scalar_branch_count::Int
    initially_closed_switch_count::Int
    topology_kinds::Vector{Symbol}
end

struct DeckMixedSteadyStateFrequencyConflict <: Exception
    connection_kind::Symbol
    connection_index::Int
    from_node_index::Int
    to_node_index::Int
    from_source_row_index::Int
    to_source_row_index::Int
    from_frequency_hz::Float64
    to_frequency_hz::Float64
end

function Base.showerror(io::IO, error::DeckMixedSteadyStateFrequencyConflict)
    print(
        io,
        "mixed steady-state source frequencies ",
        error.from_frequency_hz,
        " Hz (source row ",
        error.from_source_row_index,
        ") and ",
        error.to_frequency_hz,
        " Hz (source row ",
        error.to_source_row_index,
        ") meet across ",
        error.connection_kind,
        " ",
        error.connection_index,
        " joining nodes ",
        error.from_node_index,
        " and ",
        error.to_node_index,
    )
end

struct _DeckSteadyStateFrequencyEdge
    kind::Symbol
    index::Int
    from_node_index::Int
    to_node_index::Int
end

function _active_steady_state_frequency_source(row::DeckOVER5ASourceRow)
    return abs(Int(row.iform)) == 14 &&
           (Float64(row.tstart) == 5432.0 || Float64(row.tstart) < 0.0)
end

_source_frequency_hz(row::DeckOVER5ASourceRow) = Float64(row.sfreq) / (2.0 * pi)

function _steady_state_frequency_edges(result::DeckParseResult)
    edges = _DeckSteadyStateFrequencyEdge[]
    for (index, row) in enumerate(result.over2_branch_rows)
        from_node = abs(Int(row.from_node_value))
        to_node = abs(Int(row.to_node_value))
        (from_node == 0 || to_node == 0) && continue
        push!(edges, _DeckSteadyStateFrequencyEdge(:branch, index, from_node, to_node))
    end
    for (index, row) in enumerate(result.over5_switch_rows)
        row.initially_closed || continue
        from_node = abs(Int(row.from_node_value))
        to_node = abs(Int(row.to_node_value))
        (from_node == 0 || to_node == 0) && continue
        push!(edges, _DeckSteadyStateFrequencyEdge(:switch, index, from_node, to_node))
    end
    coupled_group_index = 0
    coupled_group_nodes = Int[]
    for row in result.coupled_line_rows
        if row.phase_index == 1 && !isempty(coupled_group_nodes)
            coupled_group_index += 1
            _append_steady_state_frequency_clique_edges!(
                edges,
                :coupled_line,
                coupled_group_index,
                coupled_group_nodes,
            )
            empty!(coupled_group_nodes)
        end
        append!(
            coupled_group_nodes,
            (
                abs(Int(row.from_node_value)),
                abs(Int(row.to_node_value)),
                ismissing(row.reference_from_node_value) ?
                    0 : abs(Int(row.reference_from_node_value)),
                ismissing(row.reference_to_node_value) ?
                    0 : abs(Int(row.reference_to_node_value)),
            ),
        )
    end
    if !isempty(coupled_group_nodes)
        coupled_group_index += 1
        _append_steady_state_frequency_clique_edges!(
            edges,
            :coupled_line,
            coupled_group_index,
            coupled_group_nodes,
        )
    end
    phase_pi_group_index = 0
    phase_pi_group_nodes = Int[]
    for row in result.coupled_phase_pi_section_rows
        if row.phase_index == 1 && !isempty(phase_pi_group_nodes)
            phase_pi_group_index += 1
            _append_steady_state_frequency_clique_edges!(
                edges,
                :coupled_phase_pi,
                phase_pi_group_index,
                phase_pi_group_nodes,
            )
            empty!(phase_pi_group_nodes)
        end
        append!(
            phase_pi_group_nodes,
            (
                abs(Int(row.from_node_value)),
                abs(Int(row.to_node_value)),
                ismissing(row.reference_from_node_value) ?
                    0 : abs(Int(row.reference_from_node_value)),
                ismissing(row.reference_to_node_value) ?
                    0 : abs(Int(row.reference_to_node_value)),
            ),
        )
    end
    if !isempty(phase_pi_group_nodes)
        phase_pi_group_index += 1
        _append_steady_state_frequency_clique_edges!(
            edges,
            :coupled_phase_pi,
            phase_pi_group_index,
            phase_pi_group_nodes,
        )
    end
    for (index, row) in enumerate(result.bergeron_line_rows)
        from_node = abs(Int(row.from_node_value))
        to_node = abs(Int(row.to_node_value))
        (from_node == 0 || to_node == 0) && continue
        push!(edges, _DeckSteadyStateFrequencyEdge(
            :bergeron_line,
            index,
            from_node,
            to_node,
        ))
    end
    for (index, row) in enumerate(result.generator_equivalent_rows)
        _append_steady_state_frequency_clique_edges!(
            edges,
            :generator_equivalent,
            index,
            vcat(row.from_node_indices, row.to_node_indices),
        )
    end
    for (index, element) in enumerate(result.elements)
        hasproperty(element, :from_nodes) && hasproperty(element, :to_nodes) ||
            continue
        _append_steady_state_frequency_clique_edges!(
            edges,
            :multiphase_element,
            index,
            vcat(getproperty(element, :from_nodes), getproperty(element, :to_nodes)),
        )
    end
    return edges
end

function _append_steady_state_frequency_clique_edges!(
    edges::Vector{_DeckSteadyStateFrequencyEdge},
    kind::Symbol,
    index::Int,
    node_values,
)
    nodes = sort!(unique(Int[abs(Int(node)) for node in node_values if Int(node) != 0]))
    isempty(nodes) && return edges
    anchor = first(nodes)
    for node in Iterators.drop(nodes, 1)
        push!(edges, _DeckSteadyStateFrequencyEdge(kind, index, anchor, node))
    end
    return edges
end

function _frequency_conflict(
    edge::_DeckSteadyStateFrequencyEdge,
    left_source::Int,
    right_source::Int,
    source_frequencies_hz::Vector{Union{Nothing,Float64}},
)
    left_frequency = source_frequencies_hz[left_source]
    right_frequency = source_frequencies_hz[right_source]
    left_frequency === nothing && error("active source frequency is missing")
    right_frequency === nothing && error("active source frequency is missing")
    return DeckMixedSteadyStateFrequencyConflict(
        edge.kind,
        edge.index,
        edge.from_node_index,
        edge.to_node_index,
        left_source,
        right_source,
        left_frequency,
        right_frequency,
    )
end

function _source_group_find!(parents::Vector{Int}, source_index::Int)
    root = source_index
    while parents[root] != root
        root = parents[root]
    end
    while parents[source_index] != source_index
        parent = parents[source_index]
        parents[source_index] = root
        source_index = parent
    end
    return root
end

function _source_group_union!(parents::Vector{Int}, left::Int, right::Int)
    left_root = _source_group_find!(parents, left)
    right_root = _source_group_find!(parents, right)
    left_root == right_root && return false
    if left_root < right_root
        parents[right_root] = left_root
    else
        parents[left_root] = right_root
    end
    return true
end

function _source_groups_and_successors!(
    parents::Vector{Int},
    active_source_row_indices::Vector{Int},
)
    members_by_root = Dict{Int,Vector{Int}}()
    for source_index in active_source_row_indices
        root = _source_group_find!(parents, source_index)
        push!(get!(members_by_root, root, Int[]), source_index)
    end
    source_groups = sort!(collect(values(members_by_root)); by = first)
    source_successors = collect(eachindex(parents))
    for members in source_groups
        sort!(members)
        for index in eachindex(members)
            source_successors[members[index]] = members[mod1(index + 1, length(members))]
        end
    end
    return source_groups, source_successors
end

function deck_steady_state_frequency_partition(result::DeckParseResult)
    node_count = maximum(values(result.node_map); init = 0)
    source_count = length(result.over5a_source_rows)
    node_sources = zeros(Int, node_count)
    source_frequencies_hz = Union{Nothing,Float64}[nothing for _ in 1:source_count]
    active_sources = Int[]
    parents = collect(1:source_count)

    for (source_index, row) in enumerate(result.over5a_source_rows)
        _active_steady_state_frequency_source(row) || continue
        frequency_hz = _source_frequency_hz(row)
        isfinite(frequency_hz) && frequency_hz > 0.0 || throw(ArgumentError(
            "steady-state sinusoidal source row $source_index requires a finite positive frequency",
        ))
        source_frequencies_hz[source_index] = frequency_hz
        push!(active_sources, source_index)
        node_index = abs(Int(row.node_value))
        1 <= node_index <= node_count || throw(ArgumentError(
            "steady-state source row $source_index references a node outside the deck network",
        ))
        prior_source = node_sources[node_index]
        if prior_source == 0
            node_sources[node_index] = source_index
        else
            prior_frequency = source_frequencies_hz[prior_source]
            prior_frequency === nothing && error("active source frequency is missing")
            if frequency_hz != prior_frequency
                throw(
                    DeckMixedSteadyStateFrequencyConflict(
                        :source_node,
                        source_index,
                        node_index,
                        node_index,
                        prior_source,
                        source_index,
                        prior_frequency,
                        frequency_hz,
                    ),
                )
            end
            _source_group_union!(parents, prior_source, source_index)
        end
    end

    active_frequencies = unique(
        source_frequencies_hz[index]::Float64 for index in active_sources
    )
    edges = _steady_state_frequency_edges(result)
    if length(active_frequencies) > 1
        changed = true
        while changed
            changed = false
            for edge in edges
                left_source = node_sources[edge.from_node_index]
                right_source = node_sources[edge.to_node_index]
                if left_source == 0 && right_source == 0
                    continue
                elseif left_source == 0
                    node_sources[edge.from_node_index] = right_source
                    changed = true
                    continue
                elseif right_source == 0
                    node_sources[edge.to_node_index] = left_source
                    changed = true
                    continue
                end
                left_frequency = source_frequencies_hz[left_source]
                right_frequency = source_frequencies_hz[right_source]
                left_frequency === nothing && error("active source frequency is missing")
                right_frequency === nothing && error("active source frequency is missing")
                left_frequency == right_frequency || throw(
                    _frequency_conflict(
                        edge,
                        left_source,
                        right_source,
                        source_frequencies_hz,
                    ),
                )
                _source_group_union!(parents, left_source, right_source)
            end
        end
    end

    first_active_source = isempty(active_sources) ? 0 : first(active_sources)
    inactive_nodes = Int[]
    if first_active_source != 0
        for node_index in eachindex(node_sources)
            if node_sources[node_index] == 0
                node_sources[node_index] = -first_active_source
                push!(inactive_nodes, node_index)
            end
        end
    else
        append!(inactive_nodes, eachindex(node_sources))
    end

    node_frequencies_hz = zeros(Float64, node_count)
    for node_index in eachindex(node_sources)
        source_index = abs(node_sources[node_index])
        source_index == 0 && continue
        frequency_hz = source_frequencies_hz[source_index]
        frequency_hz === nothing && error("node frequency source is missing")
        node_frequencies_hz[node_index] = frequency_hz
    end

    source_groups, source_successors =
        _source_groups_and_successors!(parents, active_sources)
    topology_kinds = unique(edge.kind for edge in edges)
    positive_source_roots = sort!(unique(Int[
        _source_group_find!(parents, source_index)
        for source_index in node_sources
        if source_index > 0
    ]))
    subnetwork_frequencies_hz = Float64[]
    subnetwork_node_indices = Vector{Int}[]
    for source_root in positive_source_roots
        frequency_hz = source_frequencies_hz[source_root]
        frequency_hz === nothing && error("subnetwork frequency source is missing")
        push!(subnetwork_frequencies_hz, frequency_hz)
        push!(
            subnetwork_node_indices,
            findall(
                source_index ->
                    source_index > 0 &&
                    _source_group_find!(parents, source_index) == source_root,
                node_sources,
            ),
        )
    end

    return DeckSteadyStateFrequencyPartition(
        node_sources,
        node_frequencies_hz,
        active_sources,
        source_frequencies_hz,
        source_successors,
        source_groups,
        subnetwork_frequencies_hz,
        subnetwork_node_indices,
        inactive_nodes,
        length(result.over2_branch_rows),
        count(row -> row.initially_closed, result.over5_switch_rows),
        topology_kinds,
    )
end

function deck_node_steady_state_frequency_hz(
    partition::DeckSteadyStateFrequencyPartition,
    from_node_index::Integer,
    to_node_index::Integer,
    default_frequency_hz::Real,
)
    from_node = abs(Int(from_node_index))
    to_node = abs(Int(to_node_index))
    frequency_hz = nothing
    for node_index in (from_node, to_node)
        node_index == 0 && continue
        1 <= node_index <= length(partition.node_frequencies_hz) || throw(ArgumentError(
            "steady-state branch endpoint is outside the frequency partition",
        ))
        node_frequency = partition.node_frequencies_hz[node_index]
        node_frequency == 0.0 && continue
        if frequency_hz === nothing
            frequency_hz = node_frequency
        elseif node_frequency != frequency_hz
            throw(ArgumentError("steady-state branch endpoints have different frequencies"))
        end
    end
    return frequency_hz === nothing ? Float64(default_frequency_hz) : frequency_hz
end

function validate_steady_state_frequency_networks!(result::DeckParseResult)
    partition = try
        deck_steady_state_frequency_partition(result)
    catch error
        if error isa DeckMixedSteadyStateFrequencyConflict
            record_card!(result, :mixed_frequency_subnetwork_conflict)
            add_issue!(
                result.validation,
                invalid_value("steady-state frequency topology", sprint(showerror, error)),
            )
            return result
        end
        rethrow()
    end
    distinct_frequencies = unique(
        filter(>(0.0), partition.node_frequencies_hz),
    )
    if length(distinct_frequencies) > 1
        record_card!(result, :mixed_frequency_subnetwork_partition)
        record_card!(result, :mixed_frequency_subnetwork_scalar_branch_topology)
        partition.initially_closed_switch_count > 0 &&
            record_card!(result, :mixed_frequency_subnetwork_switch_topology)
    end
    return result
end

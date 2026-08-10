struct _EMTInitializationRefusal <: Exception
    code::Symbol
    owner::Symbol
    quantity::Symbol
    message::String
    context::NamedTuple
end

function Base.showerror(io::IO, refusal::_EMTInitializationRefusal)
    print(io, refusal.message)
end

function _empty_emt_initialization_topology_report()
    return EMTInitializationTopologyReport(
        0,
        0,
        Vector{Vector{Int}}(),
        Vector{Vector{Int}}(),
        BitVector(),
        Vector{Vector{Int}}(),
        0,
        Inf,
        Inf,
        Inf,
        :not_evaluated,
    )
end

function _emt_initialization_topology_report(diagnostics, node_count::Int)
    return EMTInitializationTopologyReport(
        node_count,
        diagnostics.reduced_node_count,
        copy.(diagnostics.switch_node_groups),
        copy.(diagnostics.connected_components),
        copy(diagnostics.referenced_components),
        copy.(diagnostics.unreferenced_components),
        diagnostics.numerical_rank,
        diagnostics.condition_estimate,
        diagnostics.maximum_residual_a,
        diagnostics.relative_residual,
        diagnostics.classification,
    )
end

function _emt_initialization_classification_failure(topology::EMTInitializationTopologyReport)
    code = topology.classification in (
        :islanded,
        :nonunique,
        :infeasible,
        :ill_conditioned,
    ) ? topology.classification : :network_equilibrium
    return _EMTInitializationRefusal(
        code,
        :network_topology,
        :nodal_voltage,
        "EMT initialization network classification $(topology.classification): " *
        "rank $(topology.numerical_rank)/$(topology.reduced_node_count), " *
        "condition $(topology.condition_estimate), residual " *
        "$(topology.maximum_residual_a) A",
        (
            node_count=topology.node_count,
            reduced_node_count=topology.reduced_node_count,
            unreferenced_components=copy.(topology.unreferenced_components),
        ),
    )
end

function _emt_admittance_symmetry_error(
    admittance::AbstractMatrix{ComplexF64},
)
    error = 0.0
    for column in axes(admittance, 2), row in axes(admittance, 1)
        error = max(
            error,
            abs(admittance[row, column] - admittance[column, row]),
        )
    end
    return error
end

function _emt_minimum_dissipative_eigenvalue(
    admittance::AbstractMatrix{ComplexF64},
)
    row_count, column_count = size(admittance)
    row_count == column_count || throw(DimensionMismatch(
        "dissipative admittance diagnostic requires a square matrix",
    ))
    if all(value -> iszero(imag(value)), admittance)
        dissipative = Matrix{Float64}(undef, row_count, column_count)
        for column in 1:column_count, row in 1:row_count
            dissipative[row, column] = 0.5 * (
                real(admittance[row, column]) +
                real(admittance[column, row])
            )
        end
        return minimum(eigvals!(Symmetric(dissipative)); init=0.0)
    end
    dissipative = Matrix{ComplexF64}(undef, row_count, column_count)
    for column in 1:column_count, row in 1:row_count
        dissipative[row, column] = 0.5 * (
            admittance[row, column] + conj(admittance[column, row])
        )
    end
    return minimum(eigvals!(Hermitian(dissipative)); init=0.0)
end

function _emt_initialization_frequency_point_evidence(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
    frequency_hz::Float64,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition;
    frequency_assignment::Symbol,
    apply_operating_constraints::Bool,
    passive_conductance_network::Bool=false,
)
    node_count = maximum(values(parsed.node_map); init=0)
    node_count > 0 || throw(_EMTInitializationRefusal(
        :missing_network,
        :network_topology,
        :node,
        "EMT initialization requires at least one network node",
        (source=parsed.source,),
    ))
    admittance, rhs, representatives = _deck_steady_state_nodal_equations(
        parsed,
        node_count;
        frequency_partition,
        formulation=request.formulation,
        default_frequency_hz=frequency_hz,
    )
    fixed_node_phasors = apply_operating_constraints ?
        _emt_operating_point_voltage_constraints(
            parsed,
            request,
            frequency_hz,
            representatives,
            frequency_partition,
        ) : Dict{Int,ComplexF64}()
    tolerances = request.tolerances
    natural_diagnostics = if passive_conductance_network &&
                             !isempty(fixed_node_phasors)
        _solve_grouped_harmonic_linear_system(
            admittance,
            rhs,
            representatives;
            current_absolute_a=tolerances.current_absolute_a,
            current_relative=tolerances.current_relative,
            rank_relative_threshold_multiplier=
                tolerances.rank_relative_threshold_multiplier,
            maximum_condition_estimate=tolerances.maximum_condition_estimate,
            passive_conductance_network,
        )
    else
        nothing
    end
    redundant_operating_constraints =
        natural_diagnostics !== nothing &&
        natural_diagnostics.solution !== nothing &&
        all(fixed_node_phasors) do (representative, target)
            natural = natural_diagnostics.solution[representative]
            allowance = tolerances.voltage_absolute_v +
                tolerances.voltage_relative * max(abs(natural), abs(target))
            abs(natural - target) <= allowance
        end
    diagnostics = if redundant_operating_constraints
        natural_diagnostics
    elseif isempty(fixed_node_phasors)
        _solve_grouped_harmonic_linear_system(
            admittance,
            rhs,
            representatives;
            current_absolute_a=tolerances.current_absolute_a,
            current_relative=tolerances.current_relative,
            rank_relative_threshold_multiplier=
                tolerances.rank_relative_threshold_multiplier,
            maximum_condition_estimate=tolerances.maximum_condition_estimate,
            passive_conductance_network,
        )
    else
        _solve_grouped_constrained_harmonic_linear_system(
            admittance,
            rhs,
            representatives,
            fixed_node_phasors;
            current_absolute_a=tolerances.current_absolute_a,
            current_relative=tolerances.current_relative,
            rank_relative_threshold_multiplier=
                tolerances.rank_relative_threshold_multiplier,
            maximum_condition_estimate=tolerances.maximum_condition_estimate,
            passive_conductance_network,
        )
    end
    topology = _emt_initialization_topology_report(diagnostics, node_count)
    symmetry_error = _emt_admittance_symmetry_error(admittance)
    minimum_dissipative_eigenvalue =
        passive_conductance_network && topology.classification === :unique ?
        0.0 : _emt_minimum_dissipative_eigenvalue(admittance)
    physical_angular_frequency = 2.0 * pi * frequency_hz
    reactive_angular_frequency = _emt_reactive_angular_frequency(
        request.formulation,
        physical_angular_frequency,
    )
    node_voltage_phasors = diagnostics.solution === nothing ? ComplexF64[] :
        ComplexF64.(diagnostics.solution)
    source_injection_phasors = ComplexF64.(rhs)
    operating_constraint_current_phasors =
        hasproperty(diagnostics, :constraint_reaction_current_phasors) ?
        ComplexF64.(diagnostics.constraint_reaction_current_phasors) :
        zeros(ComplexF64, node_count)
    if request.time_origin_s != 0.0
        for node in 1:node_count
            rotation = cis(
                2.0 * pi *
                frequency_partition.node_frequencies_hz[node] *
                request.time_origin_s,
            )
            node_voltage_phasors[node] *= rotation
            source_injection_phasors[node] *= rotation
            operating_constraint_current_phasors[node] *= rotation
        end
    end
    matrix_scale = max(norm(admittance, Inf), 1.0)
    symmetry_passed = symmetry_error <= 1.0e-11 * matrix_scale
    passivity_floor = -1.0e-11 * matrix_scale
    passed = topology.classification === :unique &&
        symmetry_passed &&
        minimum_dissipative_eigenvalue >= passivity_floor &&
        !isempty(node_voltage_phasors)
    point = EMTInitializationFrequencyPoint(
        _emt_harmonic_formulation_symbol(request.formulation),
        frequency_assignment,
        frequency_hz,
        reactive_angular_frequency,
        copy(frequency_partition.node_frequencies_hz),
        copy(frequency_partition.node_source_row_indices),
        copy(frequency_partition.source_successor_indices),
        length(frequency_partition.subnetwork_node_indices),
        node_voltage_phasors,
        source_injection_phasors,
        operating_constraint_current_phasors,
        topology,
        symmetry_error,
        minimum_dissipative_eigenvalue,
        passed,
    )
    return (;
        point,
        natural_equilibrium_confirmed=redundant_operating_constraints ||
            isempty(fixed_node_phasors),
    )
end

function _emt_initialization_frequency_point(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
    frequency_hz::Float64,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition;
    frequency_assignment::Symbol,
    apply_operating_constraints::Bool,
    passive_conductance_network::Bool=false,
)
    return _emt_initialization_frequency_point_evidence(
        parsed,
        request,
        frequency_hz,
        frequency_partition;
        frequency_assignment,
        apply_operating_constraints,
        passive_conductance_network,
    ).point
end

function _emt_initialization_frequency_point(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
    frequency_hz::Float64,
)
    return _emt_initialization_frequency_point(
        parsed,
        request,
        frequency_hz,
        _network_frequency_scan_partition(parsed, frequency_hz);
        frequency_assignment=:uniform_scan,
        apply_operating_constraints=false,
        passive_conductance_network=
            _emt_frequency_invariant_passive_conductance_network(parsed),
    )
end

function _emt_operating_point_quantity_target(
    quantity::OperatingPointQuantity,
    operating_point::EMTOperatingPoint,
)
    quantity.quantity === :node_voltage_peak_phasor || throw(
        _EMTInitializationRefusal(
            :unsupported_quantity,
            :operating_point_mapping,
            quantity.quantity,
            "only explicit peak node-voltage phasors are admitted for EMT operating-point constraints",
            (asset=quantity.asset, phase=quantity.phase),
        ),
    )
    quantity.unit in ("V", "kV", "MV", "pu") || throw(
        _EMTInitializationRefusal(
            :unknown_unit,
            :operating_point_mapping,
            quantity.quantity,
            "unsupported operating-point voltage unit $(quantity.unit)",
            (asset=quantity.asset, unit=quantity.unit),
        ),
    )
    quantity.basis in (
        "peak_node_to_ground",
        "peak_phase_to_ground",
        "absolute_si_peak",
        "per_unit_peak_node_to_ground",
    ) || throw(_EMTInitializationRefusal(
        :unknown_basis,
        :operating_point_mapping,
        quantity.quantity,
        "unsupported operating-point voltage basis $(quantity.basis)",
        (asset=quantity.asset, basis=quantity.basis),
    ))
    quantity.orientation == "node_to_ground" || throw(
        _EMTInitializationRefusal(
            :unknown_orientation,
            :operating_point_mapping,
            quantity.quantity,
            "operating-point node voltage must use node_to_ground orientation",
            (asset=quantity.asset, orientation=quantity.orientation),
        ),
    )
    quantity.phase === :not_applicable ||
        quantity.phase in operating_point.phase_order || throw(
        _EMTInitializationRefusal(
            :invalid_phase,
            :operating_point_mapping,
            quantity.quantity,
            "operating-point phase is absent from the declared phase order",
            (asset=quantity.asset, phase=quantity.phase),
        ),
    )
    quantity.provenance.units == quantity.unit || throw(
        _EMTInitializationRefusal(
            :provenance_mismatch,
            :operating_point_mapping,
            quantity.quantity,
            "operating-point quantity and provenance units do not match",
            (
                asset=quantity.asset,
                quantity_unit=quantity.unit,
                provenance_units=quantity.provenance.units,
            ),
        ),
    )
    quantity.provenance.nature in (
        PhysicalModelParameter,
        ScalingBasisParameter,
    ) || throw(_EMTInitializationRefusal(
        :wrong_parameter_nature,
        :operating_point_mapping,
        quantity.quantity,
        "an imported operating point cannot be numerical-policy data",
        (asset=quantity.asset, nature=quantity.provenance.nature),
    ))
    return quantity.orientation_sign * quantity.scale_to_si * quantity.value
end

function _emt_operating_point_voltage_constraints(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
    frequency_hz::Float64,
    representatives::AbstractVector{<:Integer},
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition,
)
    operating_point = request.operating_point
    operating_point === nothing && return Dict{Int,ComplexF64}()
    frequency_hz == request.frequency_hz || return Dict{Int,ComplexF64}()
    operating_point isa EMTOperatingPoint || throw(_EMTInitializationRefusal(
        :unsupported_schema,
        :operating_point_mapping,
        :schema,
        "operating_point must be an EMTOperatingPoint",
        (runtime_type=string(typeof(operating_point)),),
    ))
    _validate_emt_operating_point_signatures(operating_point, request)
    fixed_by_representative = Dict{Int,Tuple{ComplexF64,Float64,Symbol}}()
    for quantity in operating_point.quantities
        node = get(parsed.node_map, quantity.asset, 0)
        node > 0 || throw(_EMTInitializationRefusal(
            :missing_asset,
            :operating_point_mapping,
            quantity.quantity,
            "operating-point asset $(quantity.asset) is not a network node",
            (asset=quantity.asset,),
        ))
        node_frequency_hz = frequency_partition.node_frequencies_hz[node]
        isapprox(
            node_frequency_hz,
            operating_point.frequency_hz;
            atol=1.0e-12,
            rtol=1.0e-12,
        ) || throw(_EMTInitializationRefusal(
            :frequency_mapping_mismatch,
            :operating_point_mapping,
            quantity.quantity,
            "operating-point quantity $(quantity.asset) belongs to a different frequency subnetwork",
            (
                asset=quantity.asset,
                operating_point_frequency_hz=operating_point.frequency_hz,
                node_frequency_hz,
            ),
        ))
        target_at_origin = _emt_operating_point_quantity_target(
            quantity,
            operating_point,
        )
        target_reference = target_at_origin /
            cis(2.0 * pi * node_frequency_hz * request.time_origin_s)
        representative = Int(representatives[node])
        uncertainty = quantity.absolute_uncertainty * quantity.scale_to_si
        if haskey(fixed_by_representative, representative)
            previous, previous_uncertainty, previous_asset =
                fixed_by_representative[representative]
            allowance = request.tolerances.voltage_absolute_v +
                request.tolerances.voltage_relative *
                    max(abs(previous), abs(target_reference)) +
                previous_uncertainty + uncertainty
            abs(previous - target_reference) <= allowance || throw(
                _EMTInitializationRefusal(
                    :mode_inconsistent,
                    :operating_point_mapping,
                    quantity.quantity,
                    "closed-switch node group has inconsistent operating-point voltages",
                    (
                        first_asset=previous_asset,
                        second_asset=quantity.asset,
                        residual_v=abs(previous - target_reference),
                        allowance_v=allowance,
                    ),
                ),
            )
        else
            fixed_by_representative[representative] = (
                target_reference,
                uncertainty,
                quantity.asset,
            )
        end
    end
    return Dict{Int,ComplexF64}(
        representative => value[1]
        for (representative, value) in fixed_by_representative
    )
end

function _emt_initialization_scan(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
)
    frequencies = sort!(unique(vcat(
        request.frequency_hz,
        request.frequency_grid_hz,
    )))
    declared_partition = DeckParser.deck_steady_state_frequency_partition(parsed)
    declared_source_frequencies = sort!(unique(Float64[
        declared_partition.source_frequencies_hz[source_index]
        for source_index in declared_partition.active_source_row_indices
    ]))
    mixed_declared_subnetworks = length(declared_source_frequencies) > 1
    primary_partition = mixed_declared_subnetworks ? declared_partition :
        _network_frequency_scan_partition(parsed, request.frequency_hz)
    invariant_named_network = !mixed_declared_subnetworks &&
        _emt_frequency_invariant_named_network(parsed)
    passive_conductance_network = invariant_named_network &&
        _emt_frequency_invariant_passive_conductance_network(parsed)
    primary_evidence = _emt_initialization_frequency_point_evidence(
        parsed,
        request,
        request.frequency_hz,
        primary_partition;
        frequency_assignment=:initial_operating_point,
        apply_operating_constraints=true,
        passive_conductance_network,
    )
    primary_point = primary_evidence.point
    invariant_source_frequencies_hz = invariant_named_network ? Float64[
        element.value.frequency
        for element in parsed.elements
        if element isa Union{TheveninSource,CurrentInjection}
    ] : Float64[]
    points = EMTInitializationFrequencyPoint[primary_point]
    invariant_unforced_point = nothing
    for frequency_hz in frequencies
        !mixed_declared_subnetworks && frequency_hz == request.frequency_hz &&
            continue
        source_active = any(invariant_source_frequencies_hz) do source_frequency_hz
            isapprox(
                frequency_hz,
                source_frequency_hz;
                atol=1.0e-12,
                rtol=1.0e-12,
            )
        end
        if invariant_named_network &&
           primary_evidence.natural_equilibrium_confirmed &&
           !source_active
            if invariant_unforced_point === nothing
                invariant_unforced_point =
                    _emt_frequency_invariant_unforced_scan_point(
                        primary_point,
                        request,
                        frequency_hz,
                    )
                push!(points, invariant_unforced_point)
            else
                push!(
                    points,
                    _emt_frequency_invariant_scan_point(
                        invariant_unforced_point,
                        request,
                        frequency_hz,
                    ),
                )
            end
            continue
        end
        frequency_partition = _network_frequency_scan_partition(
            parsed,
            frequency_hz,
        )
        point = _emt_initialization_frequency_point(
            parsed,
            request,
            frequency_hz,
            frequency_partition;
            frequency_assignment=:uniform_scan,
            apply_operating_constraints=false,
            passive_conductance_network,
        )
        push!(
            points,
            point,
        )
    end
    return points
end

function _emt_frequency_invariant_named_network(
    parsed::DeckParser.DeckParseResult,
)
    isempty(parsed.elements) && return false
    all(
        element -> element isa Union{
            ConductanceBranch,
            TheveninSource,
            CurrentInjection,
        },
        parsed.elements,
    ) || return false
    owners = _basic_harmonic_element_owners(parsed)
    return length(owners) == length(parsed.elements) &&
        all(owner -> owner === :named_basic_element, owners)
end

function _emt_frequency_invariant_passive_conductance_network(
    parsed::DeckParser.DeckParseResult,
)
    _emt_frequency_invariant_named_network(parsed) || return false
    return all(parsed.elements) do element
        if element isa ConductanceBranch
            isfinite(element.g) && element.g >= 0.0
        elseif element isa TheveninSource
            isfinite(element.g) && element.g >= 0.0
        else
            element isa CurrentInjection
        end
    end
end

function _emt_frequency_invariant_scan_point(
    template::EMTInitializationFrequencyPoint,
    request::EMTInitializationRequest,
    frequency_hz::Float64,
)
    node_count = length(template.node_voltage_phasors)
    node_physical_frequencies_hz = fill(frequency_hz, node_count)
    node_voltage_phasors = template.node_voltage_phasors
    source_injection_phasors = template.source_injection_phasors
    operating_constraint_current_phasors =
        template.operating_constraint_current_phasors
    if request.time_origin_s != 0.0
        node_voltage_phasors = copy(node_voltage_phasors)
        source_injection_phasors = copy(source_injection_phasors)
        operating_constraint_current_phasors =
            copy(operating_constraint_current_phasors)
        for node in 1:node_count
            rotation = cis(
                2.0 * pi *
                (node_physical_frequencies_hz[node] -
                 template.node_physical_frequencies_hz[node]) *
                request.time_origin_s,
            )
            node_voltage_phasors[node] *= rotation
            source_injection_phasors[node] *= rotation
            operating_constraint_current_phasors[node] *= rotation
        end
    end
    return EMTInitializationFrequencyPoint(
        template.formulation,
        :uniform_scan,
        frequency_hz,
        _emt_reactive_angular_frequency(
            request.formulation,
            2.0 * pi * frequency_hz,
        ),
        node_physical_frequencies_hz,
        template.node_frequency_source_row_indices,
        template.source_frequency_successor_indices,
        template.frequency_subnetwork_count,
        node_voltage_phasors,
        source_injection_phasors,
        operating_constraint_current_phasors,
        template.topology,
        template.admittance_symmetry_max_abs_error,
        template.minimum_dissipative_eigenvalue_s,
        template.passed,
    )
end

function _emt_frequency_invariant_unforced_scan_point(
    template::EMTInitializationFrequencyPoint,
    request::EMTInitializationRequest,
    frequency_hz::Float64,
)
    node_count = length(template.node_voltage_phasors)
    return EMTInitializationFrequencyPoint(
        template.formulation,
        :uniform_scan,
        frequency_hz,
        _emt_reactive_angular_frequency(
            request.formulation,
            2.0 * pi * frequency_hz,
        ),
        fill(frequency_hz, node_count),
        template.node_frequency_source_row_indices,
        template.source_frequency_successor_indices,
        template.frequency_subnetwork_count,
        zeros(ComplexF64, node_count),
        zeros(ComplexF64, node_count),
        zeros(ComplexF64, node_count),
        template.topology,
        template.admittance_symmetry_max_abs_error,
        template.minimum_dissipative_eigenvalue_s,
        template.passed,
    )
end

function _validate_emt_initialization_source_frequency(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
)
    deck_source_frequencies = Float64[
        abs(Float64(row.sfreq)) / (2.0 * pi)
        for row in DeckParser.deck_over5a_source_rows(parsed)
        if abs(Int(row.iform)) == 14 &&
           (Float64(row.tstart) == 5432.0 || Float64(row.tstart) < 0.0)
    ]
    for element in parsed.elements
        element isa Union{TheveninSource,CurrentInjection} || continue
        element.value isa SinusoidalSourceSignal || continue
        signal = element.value
        if signal.frequency > 0.0 && abs(signal.offset) >
           64.0 * eps(Float64) * max(abs(signal.amplitude), 1.0)
            throw(_EMTInitializationRefusal(
                :unsupported_source_spectrum,
                :source_state,
                :frequency_hz,
                "a named source with simultaneous DC offset and AC amplitude requires a multifrequency initializer",
                (
                    source_frequency_hz=signal.frequency,
                    offset=signal.offset,
                    amplitude=signal.amplitude,
                ),
            ))
        end
        signal.frequency > 0.0 && abs(signal.amplitude) > 0.0 &&
            isapprox(
                signal.frequency,
                request.frequency_hz;
                atol=1.0e-12,
                rtol=1.0e-12,
            ) || throw(_EMTInitializationRefusal(
                :source_frequency_mismatch,
                :source_state,
                :frequency_hz,
                "a named sinusoidal source does not match the primary initialization frequency",
                (
                    requested_frequency_hz=request.frequency_hz,
                    source_frequency_hz=signal.frequency,
                ),
            ))
    end
    isempty(deck_source_frequencies) && return request
    unique_frequencies = sort!(unique(deck_source_frequencies))
    any(
        frequency -> isapprox(
            frequency,
            request.frequency_hz;
            atol=1.0e-12,
            rtol=1.0e-12,
        ),
        unique_frequencies,
    ) || throw(_EMTInitializationRefusal(
        :source_frequency_mismatch,
        :source_state,
        :frequency_hz,
        "the primary initialization frequency must match one declared source subnetwork",
        (
            requested_frequency_hz=request.frequency_hz,
            source_frequencies_hz=unique_frequencies,
        ),
    ))
    return request
end

function _emt_initialization_primary_point(
    points::Vector{EMTInitializationFrequencyPoint},
    frequency_hz::Float64,
)
    index = findfirst(
        point -> point.frequency_assignment === :initial_operating_point &&
            point.physical_frequency_hz == frequency_hz,
        points,
    )
    index === nothing && error("initialization scan omitted its primary frequency")
    return points[index]
end

function _validate_emt_operating_point_signatures(
    operating_point::EMTOperatingPoint,
    request::EMTInitializationRequest,
)
    fields = (
        (:project_signature, operating_point.project_signature,
            request.project_signature),
        (:settings_signature, operating_point.settings_signature,
            request.settings_signature),
        (:model_signature, operating_point.model_signature,
            request.model_signature),
    )
    for (field, source, target) in fields
        source == target || throw(_EMTInitializationRefusal(
            :stale_signature,
            :operating_point_mapping,
            field,
            "operating-point $field does not match the initialization request",
            (source_signature=source, target_signature=target),
        ))
    end
    isapprox(
        operating_point.frequency_hz,
        request.frequency_hz;
        atol=0.0,
        rtol=64.0 * eps(Float64),
    ) || throw(_EMTInitializationRefusal(
        :frequency_mismatch,
        :operating_point_mapping,
        :frequency_hz,
        "operating-point frequency does not match the initialization frequency",
        (
            source_frequency_hz=operating_point.frequency_hz,
            target_frequency_hz=request.frequency_hz,
        ),
    ))
    operating_point.time_origin_s == request.time_origin_s || throw(
        _EMTInitializationRefusal(
            :time_origin_mismatch,
            :operating_point_mapping,
            :time_origin_s,
            "operating-point time origin does not match the initialization request",
            (
                source_time_origin_s=operating_point.time_origin_s,
                target_time_origin_s=request.time_origin_s,
            ),
        ),
    )
    return operating_point
end

function _emt_operating_point_mappings(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
    point::EMTInitializationFrequencyPoint,
)
    operating_point = request.operating_point
    operating_point === nothing && return OperatingPointMappingRecord[]
    operating_point isa EMTOperatingPoint || throw(_EMTInitializationRefusal(
        :unsupported_schema,
        :operating_point_mapping,
        :schema,
        "operating_point must be an EMTOperatingPoint",
        (runtime_type=string(typeof(operating_point)),),
    ))
    _validate_emt_operating_point_signatures(operating_point, request)
    mappings = OperatingPointMappingRecord[]
    for quantity in operating_point.quantities
        target_value = _emt_operating_point_quantity_target(
            quantity,
            operating_point,
        )
        node_index = get(parsed.node_map, quantity.asset, 0)
        node_index > 0 || throw(_EMTInitializationRefusal(
            :missing_asset,
            :operating_point_mapping,
            quantity.quantity,
            "operating-point asset $(quantity.asset) is not a network node",
            (asset=quantity.asset,),
        ))
        required_value = point.node_voltage_phasors[node_index]
        uncertainty = quantity.absolute_uncertainty * quantity.scale_to_si
        residual = abs(target_value - required_value)
        allowance = request.tolerances.voltage_absolute_v +
            request.tolerances.voltage_relative * abs(required_value) +
            uncertainty
        passed = isfinite(residual) && residual <= allowance
        switch_group_index = findfirst(
            group -> node_index in group,
            point.topology.switch_node_groups,
        )
        reaction_nodes = switch_group_index === nothing ? Int[node_index] :
            point.topology.switch_node_groups[switch_group_index]
        constraint_current = sum(
            point.operating_constraint_current_phasors[reaction_nodes];
            init=0.0 + 0.0im,
        )
        push!(
            mappings,
            OperatingPointMappingRecord(
                quantity.asset,
                quantity.quantity,
                quantity.phase,
                quantity.value,
                target_value,
                quantity.unit,
                "V",
                quantity.basis,
                quantity.orientation,
                quantity.scale_to_si,
                quantity.orientation_sign,
                uncertainty,
                residual,
                constraint_current,
                passed,
            ),
        )
        passed || throw(_EMTInitializationRefusal(
            :infeasible_operating_target,
            :operating_point_mapping,
            quantity.quantity,
            "operating-point quantity $(quantity.asset)/$(quantity.quantity) " *
            "does not close the EMT equilibrium",
            (asset=quantity.asset, residual=residual, allowance=allowance),
        ))
    end
    return mappings
end

function _emt_initial_voltage_sample(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
    point::EMTInitializationFrequencyPoint,
)
    output_node_indices = DeckParser.deck_over16_output_node_indices(parsed)
    values = real.(point.node_voltage_phasors)
    return (
        source=:consistent_emt_initialization,
        outcome=:steady_state_initial_voltage_sample,
        harmonic_formulation=point.formulation,
        exact_discrete_histories=
            request.formulation isa TimestepMatchedFormulation,
        steady_state_frequency_hz=point.physical_frequency_hz,
        reactive_angular_frequency_rad_s=
            point.reactive_angular_frequency_rad_s,
        time_origin_s=request.time_origin_s,
        timestep_s=request.formulation isa TimestepMatchedFormulation ?
            request.formulation.timestep_s :
            DeckParser.deck_fixed_time_horizon_options(parsed).dt_s,
        node_steady_state_frequencies_hz=
            copy(point.node_physical_frequencies_hz),
        node_frequency_source_row_indices=
            copy(point.node_frequency_source_row_indices),
        source_frequency_successor_indices=
            copy(point.source_frequency_successor_indices),
        steady_state_frequency_subnetwork_count=
            point.frequency_subnetwork_count,
        node_map=Dict{Symbol,Int}(parsed.node_map),
        node_names=ordered_node_names(parsed.node_map),
        node_voltage_phasors=copy(point.node_voltage_phasors),
        node_voltage_values=values,
        output_node_indices=output_node_indices,
        output_voltage_values=Float64[
            node == 0 ? 0.0 : values[node]
            for node in output_node_indices
        ],
        source_row_count=length(DeckParser.deck_over5a_source_rows(parsed)),
        distributed_transposed_line_count=length(
            DeckParser.deck_distributed_transposed_line_modal_branch_states(parsed),
        ),
        saturated_transformer_branch_count=0,
        saturated_transformer_ideal_branch_count=0,
        saturated_transformer_linearized_nonlinear_branch_count=0,
    )
end

function _shift_emt_initialization_time_origin!(
    prepared::PreparedEMTStudy,
    time_origin_s::Float64,
)
    time_origin_s == 0.0 && return prepared
    runtime = prepared.runtime_template
    plan = runtime.plan
    for source_index in eachindex(plan.source_iform_values)
        abs(plan.source_iform_values[source_index]) == 14 || continue
        plan.source_time1_values[source_index] +=
            plan.source_sfreq_values[source_index] * time_origin_s
    end
    for signal in runtime.context.analytic_source_signals
        signal.source_type == 14 || continue
        signal.time1_s += signal.angular_frequency_or_rate * time_origin_s
    end
    for element in runtime.context.system.elements
        element isa Union{TheveninSource,CurrentInjection} || continue
        element.value isa SinusoidalSourceSignal || continue
        signal = element.value
        signal.phase += 2.0 * pi * signal.frequency * time_origin_s
    end
    return prepared
end

function _append_emt_initialization_state!(
    records::Vector{EMTInitializationStateRecord},
    owner::Symbol,
    state_family::Symbol,
    instance_count::Integer,
    initialization_basis::Symbol,
)
    instance_count == 0 && return records
    push!(
        records,
        EMTInitializationStateRecord(
            owner,
            state_family,
            instance_count,
            initialization_basis,
        ),
    )
    return records
end

function _emt_initialization_state_inventory(prepared::PreparedEMTStudy)
    context = prepared.runtime_template.context
    parsed = prepared.parsed
    elements = context.system.elements
    records = EMTInitializationStateRecord[]
    _append_emt_initialization_state!(
        records,
        :network_voltage,
        :algebraic,
        context.system.node_count,
        :harmonic_network_equilibrium,
    )
    _append_emt_initialization_state!(
        records,
        :network_topology,
        :discrete,
        1,
        :ranked_connected_component_classification,
    )
    source_element_count = count(
        element -> element isa Union{TheveninSource,CurrentInjection},
        elements,
    )
    source_program_count = context.source_function_runtime === nothing ? 0 :
        max(context.source_function_runtime.plan.source_row_count, 1)
    source_count = max(
        source_element_count,
        length(context.analytic_source_signals),
        source_program_count,
    )
    _append_emt_initialization_state!(
        records,
        :source_state,
        :algebraic,
        source_count,
        :phasor_time_mapping,
    )
    history_kinds = context.electromagnetic_history_plan.kinds
    step_configs = prepared.runtime_template.step_configs
    algebraic_named_network =
        _emt_frequency_invariant_named_network(parsed) &&
        isempty(history_kinds) &&
        context.control_system_runtime === nothing &&
        !(step_configs isa DynamicDeckStepConfigProvider &&
          step_configs.nonlinear_current_config !== nothing) &&
        context.deck_time_switch_count == 0 &&
        isempty(context.series_rlc_alterations)
    if algebraic_named_network
        _append_emt_initialization_state!(
            records,
            :energy_accumulator,
            :algebraic,
            length(context.branch_energy_values) +
                length(context.switch_energy_values),
            :time_zero_energy_reference,
        )
        _append_emt_initialization_state!(
            records,
            :output_cursor,
            :discrete,
            1,
            :time_zero_output_epoch,
        )
        _append_emt_initialization_state!(
            records,
            :checkpoint_continuation_state,
            :checkpoint,
            1,
            :complete_prepared_runtime_before_first_advance,
        )
        sort!(records; by=record -> (String(record.state_family), String(record.owner)))
        return records
    end
    for (kind, owner, basis) in (
        (SERIES_RL_HISTORY, :series_rl_history, :exact_companion_recurrence),
        (SERIES_RLC_HISTORY, :series_rlc_history, :exact_companion_recurrence),
        (CAPACITOR_HISTORY, :capacitor_history, :exact_companion_recurrence),
        (COUPLED_INDUCTIVE_HISTORY, :coupled_inductive_history, :exact_coupled_recurrence),
        (COUPLED_SERIES_RL_HISTORY, :coupled_series_rl_history, :exact_coupled_recurrence),
        (BREQIV_HISTORY, :breqiv_history, :frequency_domain_history_mapping),
    )
        _append_emt_initialization_state!(
            records,
            owner,
            :delayed_history,
            count(==(kind), history_kinds),
            basis,
        )
    end
    _append_emt_initialization_state!(
        records,
        :complex_modal_line_history,
        :delayed_history,
        count(element -> element isa ComplexModalBergeronLine, elements),
        :traveling_wave_harmonic_prehistory,
    )
    _append_emt_initialization_state!(
        records,
        :frequency_dependent_line_history,
        :delayed_history,
        count(element -> element isa Union{
            SemlyenFrequencyDependentLine,
            SampledFrequencyDependentLine,
            SampledFrequencyDependentLineGroup,
        }, elements),
        :recursive_convolution_harmonic_prehistory,
    )
    _append_emt_initialization_state!(
        records,
        :distributed_line_history,
        :delayed_history,
        length(DeckParser.deck_distributed_transposed_line_modal_branch_states(parsed)),
        :modal_traveling_wave_harmonic_prehistory,
    )
    control_runtime = context.control_system_runtime
    if control_runtime !== nothing
        control_state_count = length(control_runtime.state.values) +
            length(control_runtime.state.function_states)
        _append_emt_initialization_state!(
            records,
            :control_state,
            :discrete,
            max(control_state_count, 1),
            :control_network_steady_state,
        )
        _append_emt_initialization_state!(
            records,
            :control_frequency_history,
            :delayed_history,
            length(control_runtime.frequency_initializations),
            :control_frequency_response,
        )
    end
    if step_configs isa DynamicDeckStepConfigProvider &&
       step_configs.nonlinear_current_config !== nothing
        nonlinear_types = Int.(
            step_configs.nonlinear_current_config.nonlinear_types,
        )
        saturated_transformer_count = get(
            step_configs.nonlinear_current_config,
            :saturated_transformer_residual_flux_initialized,
            false,
        ) ? length(get(
            step_configs.nonlinear_current_config,
            :saturated_transformer_internal_top_node_indices,
            Int[],
        )) : 0
        _append_emt_initialization_state!(
            records,
            :pseudo_nonlinear_inductor_state,
            :continuous,
            count(==(PSEUDO_NONLINEAR_INDUCTOR_TYPE), nonlinear_types),
            :harmonic_flux_companion_prehistory,
        )
        _append_emt_initialization_state!(
            records,
            :saturated_transformer_magnetic_state,
            :continuous,
            saturated_transformer_count,
            :declared_flux_characteristic_equilibrium,
        )
        _append_emt_initialization_state!(
            records,
            :piecewise_nonlinear_inductor_state,
            :continuous,
            count(==(PIECEWISE_NONLINEAR_INDUCTOR_TYPE), nonlinear_types),
            :harmonic_characteristic_equilibrium,
        )
        _append_emt_initialization_state!(
            records,
            :hysteretic_magnetic_state,
            :continuous,
            count(==(HYSTERETIC_INDUCTOR_NONLINEAR_TYPE), nonlinear_types),
            :model_owned_hysteresis_equilibrium,
        )
    end
    switch_count = count(element -> element isa Union{
        IdealSwitch,
        TimeSwitch,
        CurrentZeroSwitch,
        TACSControlledSwitch,
    }, elements)
    _append_emt_initialization_state!(
        records,
        :switch_mode,
        :discrete,
        switch_count,
        :declared_initial_topology,
    )
    _append_emt_initialization_state!(
        records,
        :switch_event_state,
        :scheduler,
        context.deck_time_switch_count,
        :initial_event_surface_classification,
    )
    _append_emt_initialization_state!(
        records,
        :energy_accumulator,
        :algebraic,
        length(context.branch_energy_values) + length(context.switch_energy_values),
        :time_zero_energy_reference,
    )
    _append_emt_initialization_state!(
        records,
        :series_rlc_alteration_schedule,
        :scheduler,
        length(context.series_rlc_alterations),
        :ordered_future_event_cursor,
    )
    _append_emt_initialization_state!(
        records,
        :output_cursor,
        :discrete,
        1,
        :time_zero_output_epoch,
    )
    _append_emt_initialization_state!(
        records,
        :checkpoint_continuation_state,
        :checkpoint,
        1,
        :complete_prepared_runtime_before_first_advance,
    )
    sort!(records; by=record -> (String(record.state_family), String(record.owner)))
    return records
end

function _emt_initialized_state_owners(
    records::Vector{EMTInitializationStateRecord},
)
    return sort!(unique(record.owner for record in records); by=String)
end

function _emt_unsupported_initialization_owners(
    parsed::DeckParser.DeckParseResult,
)
    elements = parsed.elements
    unsupported = Symbol[]
    any(element -> element isa PowerSemiconductorSwitch ||
        element isa PowerSemiconductorBridgeLeg, elements) &&
        push!(unsupported, :switch_detailed_periodic_prehistory)
    isempty(DeckParser.deck_universal_machine_definition_rows(parsed)) ||
        push!(unsupported, :universal_machine_state)
    isempty(DeckParser.deck_synchronous_machine_terminal_voltage_rows(parsed)) ||
        push!(unsupported, :synchronous_machine_state)
    nonlinear_periodic_count = sum(length, (
        DeckParser.deck_zinc_oxide_nonlinear_rows(parsed),
        DeckParser.deck_nonlinear_resistance_rows(parsed),
        DeckParser.deck_triggered_timed_resistance_rows(parsed),
        DeckParser.deck_switching_nonlinear_resistor_rows(parsed),
        DeckParser.deck_arrester_nonlinear_rows(parsed),
    ))
    nonlinear_periodic_count == 0 ||
        push!(unsupported, :nonlinear_periodic_state)
    any(
        row -> row.condition_kind == :node_voltage_initial_condition,
        DeckParser.deck_node_initial_condition_rows(parsed),
    ) && push!(unsupported, :declared_node_initial_condition)
    any(
        row -> row.random_opening_standard_deviation_s > 0.0,
        DeckParser.deck_over5_switch_rows(parsed),
    ) && push!(unsupported, :stochastic_switch_state)
    return sort!(unique(unsupported); by=String)
end

function _emt_initialization_probe_runtime(prepared::PreparedEMTStudy)
    runtime = deepcopy(prepared.runtime_template)
    context = runtime.context
    context.step_count >= 1 || throw(_EMTInitializationRefusal(
        :missing_probe_horizon,
        :no_artificial_transient,
        :time_window,
        "initialization requires at least one timestep for its coupled transient probe",
        (step_count=context.step_count,),
    ))
    context.t_end_s = context.dt_s
    context.step_count = 1
    context.time_s = Vector{Float64}(undef, 2)
    context.voltage_pu = Matrix{Float64}(undef, context.system.node_count, 2)
    context.output_pu = Matrix{Float64}(
        undef,
        length(context.output_channel_names),
        2,
    )
    context.recorded_step_indices = Int[0, 1]
    context.trace_write_index = 1
    return runtime
end

function _emt_initialization_probe_metric(
    prepared::PreparedEMTStudy,
    point::EMTInitializationFrequencyPoint,
    request::EMTInitializationRequest,
)
    runtime = _emt_initialization_probe_runtime(prepared)
    run = _run_prepared_dynamic_deck!(runtime; collect_run_diagnostics=true)
    trace = run.trace
    size(trace.voltage_pu, 2) >= 2 || throw(_EMTInitializationRefusal(
        :incomplete_probe,
        :no_artificial_transient,
        :node_voltage,
        "initialization transient probe did not produce two voltage samples",
        (sample_count=size(trace.voltage_pu, 2),),
    ))
    timestep_s = runtime.context.dt_s
    rotations = ComplexF64[
        cis(2.0 * pi * frequency_hz * timestep_s)
        for frequency_hz in point.node_physical_frequencies_hz
    ]
    expected = real.(point.node_voltage_phasors .* rotations)
    actual = Float64.(trace.voltage_pu[:, 2])
    error = actual - expected
    reference_scale = max(norm(expected), request.tolerances.voltage_absolute_v)
    raw_normalized_rms = norm(error) / reference_scale
    raw_scaled_discontinuity = maximum(abs, error; init=0.0) /
        max(maximum(abs, expected; init=0.0), request.tolerances.voltage_absolute_v)
    physical_to_trapezoidal_warping = if request.formulation isa
                                         PhysicalFrequencyFormulation
        maximum(point.node_physical_frequencies_hz; init=0.0) do frequency_hz
            frequency_step = pi * frequency_hz * timestep_s
            frequency_step == 0.0 ? 0.0 :
                abs(tan(frequency_step) / frequency_step - 1.0)
        end
    else
        0.0
    end
    physical_frequency_probe_bound =
        2.0 * physical_to_trapezoidal_warping /
        max(1.0 - physical_to_trapezoidal_warping, eps(Float64))
    maximum_physical_step_angle = 2.0 * pi *
        maximum(point.node_physical_frequencies_hz; init=0.0) * timestep_s
    physical_frequency_probe_bound = max(
        physical_frequency_probe_bound,
        0.5 * maximum_physical_step_angle^2,
    )
    normalized_rms = max(
        0.0,
        raw_normalized_rms - physical_frequency_probe_bound,
    )
    scaled_discontinuity = max(
        0.0,
        raw_scaled_discontinuity - physical_frequency_probe_bound,
    )
    threshold = request.tolerances.no_artificial_transient_normalized_rms
    passed = isfinite(normalized_rms) && isfinite(scaled_discontinuity) &&
        normalized_rms <= threshold &&
        scaled_discontinuity <=
            request.tolerances.first_step_scaled_discontinuity
    return NoArtificialTransientMetric(
        :node_voltage,
        "V",
        request.time_origin_s,
        request.time_origin_s + timestep_s,
        normalized_rms,
        scaled_discontinuity,
        max(raw_normalized_rms, raw_scaled_discontinuity),
        0.0,
        threshold,
        passed,
    )
end

function _emt_initialization_residuals(
    point::EMTInitializationFrequencyPoint,
    mappings::Vector{OperatingPointMappingRecord},
    request::EMTInitializationRequest,
    fixed_source_load_flow::Union{Nothing,FixedSourceLoadFlowResult}=nothing,
)
    tolerances = request.tolerances
    current_scale = max(
        norm(point.source_injection_phasors, Inf),
        tolerances.current_absolute_a,
    )
    current_allowance = tolerances.current_absolute_a +
        tolerances.current_relative * current_scale
    residuals = EMTInitializationResidual[
        EMTInitializationResidual(
            :network_equilibrium,
            :nodal_kcl,
            "A",
            point.topology.maximum_residual_a,
            tolerances.current_absolute_a,
            tolerances.current_relative,
            current_scale,
            0.0,
            point.topology.maximum_residual_a / current_allowance,
            point.topology.maximum_residual_a <= current_allowance,
        ),
    ]
    for mapping in mappings
        allowance = tolerances.voltage_absolute_v +
            tolerances.voltage_relative * abs(mapping.target_value_si) +
            mapping.absolute_uncertainty_si
        push!(
            residuals,
            EMTInitializationResidual(
                :operating_point_mapping,
                mapping.quantity,
                mapping.target_unit,
                mapping.residual,
                tolerances.voltage_absolute_v,
                tolerances.voltage_relative,
                abs(mapping.target_value_si),
                mapping.absolute_uncertainty_si,
                mapping.residual / allowance,
                mapping.passed,
            ),
        )
    end
    if fixed_source_load_flow !== nothing
        for constraint_index in eachindex(
            fixed_source_load_flow.constraint_kinds,
        )
            for (quantity, unit, target, actual) in (
                (
                    :active_power,
                    "W",
                    fixed_source_load_flow.constraint_target_active_powers[
                        constraint_index
                    ],
                    fixed_source_load_flow.constraint_active_powers[
                        constraint_index
                    ],
                ),
                (
                    :reactive_power,
                    "var",
                    fixed_source_load_flow.constraint_target_reactive_powers[
                        constraint_index
                    ],
                    fixed_source_load_flow.constraint_reactive_powers[
                        constraint_index
                    ],
                ),
            )
                target === missing && continue
                reference_scale = max(abs(Float64(target)), 1.0)
                allowance = tolerances.power_absolute_w +
                    tolerances.power_relative * reference_scale
                residual = abs(Float64(actual) - Float64(target))
                push!(
                    residuals,
                    EMTInitializationResidual(
                        :fixed_source_operating_point,
                        quantity,
                        unit,
                        residual,
                        tolerances.power_absolute_w,
                        tolerances.power_relative,
                        reference_scale,
                        0.0,
                        residual / allowance,
                        residual <= allowance,
                    ),
                )
            end
        end
    end
    return residuals
end

function _emt_pseudo_nonlinear_initialization_residuals(
    prepared::PreparedEMTStudy,
    point::EMTInitializationFrequencyPoint,
    request::EMTInitializationRequest,
)
    step_configs = prepared.runtime_template.step_configs
    step_configs isa DynamicDeckStepConfigProvider ||
        return EMTInitializationResidual[]
    config = step_configs.nonlinear_current_config
    config === nothing && return EMTInitializationResidual[]
    nonlinear_types = Int.(config.nonlinear_types)
    nonlinear_indices = findall(
        ==(PSEUDO_NONLINEAR_INDUCTOR_TYPE),
        nonlinear_types,
    )
    isempty(nonlinear_indices) && return EMTInitializationResidual[]
    from_nodes = Int.(config.nonlinear_deck_from_nodes)
    to_nodes = Int.(config.nonlinear_deck_to_nodes)
    stored_fluxes = Float64.(config.initial_stored_voltage_values)
    companion_currents = Float64.(config.initial_companion_current_values)
    table_indices = Int.(config.initial_table_index_values)
    slopes = Float64.(config.gslope)
    delta2 = Float64(config.delta2)
    flux_residuals = Float64[]
    flux_scales = Float64[]
    current_residuals = Float64[]
    current_scales = Float64[]
    for index in nonlinear_indices
        from_node = from_nodes[index]
        to_node = to_nodes[index]
        branch_phasor =
            _node_voltage_phasor(point.node_voltage_phasors, from_node) -
            _node_voltage_phasor(point.node_voltage_phasors, to_node)
        endpoint_nodes = filter(!=(0), (from_node, to_node))
        frequency_hz = point.node_physical_frequencies_hz[first(endpoint_nodes)]
        all(
            node -> point.node_physical_frequencies_hz[node] == frequency_hz,
            endpoint_nodes,
        ) || throw(_EMTInitializationRefusal(
            :frequency_mapping_mismatch,
            :pseudo_nonlinear_inductor_state,
            :frequency_hz,
            "pseudo-nonlinear inductor endpoints belong to different frequency subnetworks",
            (owner_index=index, from_node, to_node),
        ))
        angular_frequency = 2.0 * pi * frequency_hz
        angular_frequency > 0.0 || throw(_EMTInitializationRefusal(
            :unsupported_dc_state,
            :pseudo_nonlinear_inductor_state,
            :flux,
            "pseudo-nonlinear inductor prehistory requires a positive harmonic frequency",
            (owner_index=index, frequency_hz),
        ))
        expected_flux = imag(branch_phasor) / angular_frequency -
            delta2 * real(branch_phasor)
        actual_flux = stored_fluxes[index]
        push!(flux_residuals, abs(actual_flux - expected_flux))
        push!(flux_scales, abs(expected_flux))
        table_index = table_indices[index]
        1 <= table_index <= length(slopes) || throw(
            _EMTInitializationRefusal(
                :invalid_characteristic_state,
                :pseudo_nonlinear_inductor_state,
                :table_index,
                "pseudo-nonlinear inductor initialization selected an invalid characteristic segment",
                (owner_index=index, table_index),
            ),
        )
        expected_current = expected_flux * slopes[table_index] / delta2
        actual_current = companion_currents[index]
        push!(current_residuals, abs(actual_current - expected_current))
        push!(current_scales, abs(expected_current))
    end
    tolerances = request.tolerances
    flux_residual = maximum(flux_residuals)
    flux_scale = maximum(flux_scales; init=0.0)
    flux_allowance = tolerances.flux_absolute_wb
    current_residual = maximum(current_residuals)
    current_scale = maximum(current_scales; init=0.0)
    current_allowance = tolerances.current_absolute_a +
        tolerances.current_relative * current_scale
    return EMTInitializationResidual[
        EMTInitializationResidual(
            :pseudo_nonlinear_inductor_state,
            :flux_history_recurrence,
            "Wb",
            flux_residual,
            tolerances.flux_absolute_wb,
            0.0,
            flux_scale,
            0.0,
            flux_residual / flux_allowance,
            flux_residual <= flux_allowance,
        ),
        EMTInitializationResidual(
            :pseudo_nonlinear_inductor_state,
            :companion_current_recurrence,
            "A",
            current_residual,
            tolerances.current_absolute_a,
            tolerances.current_relative,
            current_scale,
            0.0,
            current_residual / current_allowance,
            current_residual <= current_allowance,
        ),
    ]
end

function _emt_saturated_transformer_initialization_residuals(
    prepared::PreparedEMTStudy,
    point::EMTInitializationFrequencyPoint,
    request::EMTInitializationRequest,
)
    step_configs = prepared.runtime_template.step_configs
    step_configs isa DynamicDeckStepConfigProvider ||
        return EMTInitializationResidual[]
    config = step_configs.nonlinear_current_config
    config === nothing && return EMTInitializationResidual[]
    get(
        config,
        :saturated_transformer_residual_flux_initialized,
        false,
    ) || return EMTInitializationResidual[]
    count = length(get(
        config,
        :saturated_transformer_internal_top_node_indices,
        Int[],
    ))
    count > 0 || return EMTInitializationResidual[]
    from_nodes = Int.(config.nonlinear_from_nodes[1:count])
    to_nodes = Int.(config.nonlinear_to_nodes[1:count])
    declared_currents = Float64.(
        config.nonlinear_steady_state_current_values[1:count],
    )
    declared_fluxes = Float64.(
        config.nonlinear_steady_state_flux_values[1:count],
    )
    stored_fluxes = Float64.(config.initial_stored_voltage_values[1:count])
    companion_currents = Float64.(
        config.initial_companion_current_values[1:count],
    )
    table_indices = Int.(config.initial_table_index_values[1:count])
    slopes = Float64.(config.gslope)
    delta2 = Float64(config.delta2)
    flux_errors = Float64[]
    current_errors = Float64[]
    current_scales = Float64[]
    residual_fluxes = Float64[]
    for index in 1:count
        branch_phasor =
            _node_voltage_phasor(point.node_voltage_phasors, from_nodes[index]) -
            _node_voltage_phasor(point.node_voltage_phasors, to_nodes[index])
        branch_voltage = real(branch_phasor)
        reconstructed_flux = stored_fluxes[index] + delta2 * branch_voltage
        push!(flux_errors, abs(reconstructed_flux - declared_fluxes[index]))
        table_index = table_indices[index]
        1 <= table_index <= length(slopes) || throw(
            _EMTInitializationRefusal(
                :invalid_characteristic_state,
                :saturated_transformer_magnetic_state,
                :table_index,
                "saturated-transformer initialization selected an invalid characteristic segment",
                (transformer_index=index, table_index),
            ),
        )
        reconstructed_current = companion_currents[index] +
            slopes[table_index] * branch_voltage
        push!(current_errors, abs(reconstructed_current - declared_currents[index]))
        push!(current_scales, abs(declared_currents[index]))
        harmonic_flux = point.physical_frequency_hz == 0.0 ? 0.0 :
            imag(branch_phasor) / (2.0 * pi * point.physical_frequency_hz)
        push!(residual_fluxes, declared_fluxes[index] - harmonic_flux)
    end
    tolerances = request.tolerances
    flux_error = maximum(flux_errors; init=0.0)
    current_error = maximum(current_errors; init=0.0)
    current_scale = maximum(current_scales; init=0.0)
    current_allowance = tolerances.current_absolute_a +
        tolerances.current_relative * current_scale
    maximum_residual_flux = maximum(abs, residual_fluxes; init=0.0)
    return EMTInitializationResidual[
        EMTInitializationResidual(
            :saturated_transformer_magnetic_state,
            :declared_time_zero_flux,
            "Wb-turn",
            flux_error,
            tolerances.flux_absolute_wb,
            0.0,
            maximum(abs, declared_fluxes; init=0.0),
            0.0,
            flux_error / tolerances.flux_absolute_wb,
            flux_error <= tolerances.flux_absolute_wb,
        ),
        EMTInitializationResidual(
            :saturated_transformer_magnetic_state,
            :characteristic_current_equilibrium,
            "A",
            current_error,
            tolerances.current_absolute_a,
            tolerances.current_relative,
            current_scale,
            0.0,
            current_error / current_allowance,
            current_error <= current_allowance,
        ),
        EMTInitializationResidual(
            :saturated_transformer_magnetic_state,
            :residual_flux_magnitude,
            "Wb-turn",
            maximum_residual_flux,
            max(maximum_residual_flux, tolerances.flux_absolute_wb),
            0.0,
            maximum_residual_flux,
            0.0,
            maximum_residual_flux == 0.0 ? 0.0 : 1.0,
            true,
        ),
    ]
end

function _emt_piecewise_nonlinear_initialization_residuals(
    prepared::PreparedEMTStudy,
    point::EMTInitializationFrequencyPoint,
    request::EMTInitializationRequest,
)
    step_configs = prepared.runtime_template.step_configs
    step_configs isa DynamicDeckStepConfigProvider ||
        return EMTInitializationResidual[]
    config = step_configs.nonlinear_current_config
    config === nothing && return EMTInitializationResidual[]
    nonlinear_types = Int.(config.nonlinear_types)
    indices = findall(==(PIECEWISE_NONLINEAR_INDUCTOR_TYPE), nonlinear_types)
    isempty(indices) && return EMTInitializationResidual[]
    from_nodes = Int.(config.nonlinear_deck_from_nodes)
    to_nodes = Int.(config.nonlinear_deck_to_nodes)
    table_starts = Int.(config.nonlinear_admittance_nodes)
    table_ends = Int.(config.nonlinear_table_end_indices)
    currents_a = Float64.(config.cchar)
    fluxes_wb = Float64.(config.vchar)
    initialized_fluxes = Float64.(
        config.piecewise_nonlinear_inductor_initial_flux_values,
    )
    predictor_fluxes = Float64.(
        config.piecewise_nonlinear_inductor_initial_predictor_flux_values,
    )
    initialized_currents = Float64.(
        config.piecewise_nonlinear_inductor_initial_current_values,
    )
    initialized_segments = Int.(
        config.piecewise_nonlinear_inductor_initial_segment_values,
    )
    delta2 = Float64(config.delta2)
    characteristic_residuals = Float64[]
    characteristic_scales = Float64[]
    predictor_residuals = Float64[]
    predictor_scales = Float64[]
    current_residuals = Float64[]
    current_scales = Float64[]
    segment_residuals = Float64[]
    for index in indices
        from_node = from_nodes[index]
        to_node = to_nodes[index]
        branch_phasor =
            _node_voltage_phasor(point.node_voltage_phasors, from_node) -
            _node_voltage_phasor(point.node_voltage_phasors, to_node)
        frequency_hz = _nonlinear_initial_frequency_hz(
            prepared.runtime_template.steady_state_initial_sample,
            from_node,
            to_node,
        )
        reactive_angular_frequency = _steady_state_reactive_angular_frequency(
            prepared.runtime_template.steady_state_initial_sample,
            frequency_hz,
        )
        expected_flux_wb = imag(branch_phasor) / reactive_angular_frequency
        expected_predictor_flux_wb = expected_flux_wb + delta2 * real(branch_phasor)
        state = _piecewise_nonlinear_inductor_characteristic_state(
            initialized_currents[index],
            initialized_fluxes[index],
            table_starts[index],
            table_ends[index],
            currents_a,
            fluxes_wb;
            flux_tolerance_wb=Float64(config.flzero),
        )
        push!(
            characteristic_residuals,
            abs(initialized_fluxes[index] - state.flux_wb),
        )
        push!(characteristic_scales, abs(state.flux_wb))
        push!(
            predictor_residuals,
            abs(predictor_fluxes[index] - expected_predictor_flux_wb),
        )
        push!(predictor_scales, abs(expected_predictor_flux_wb))
        expected_current_a = begin
            declared_current_a = Float64(
                config.nonlinear_steady_state_current_values[index],
            )
            declared_flux_wb = Float64(
                config.nonlinear_steady_state_flux_values[index],
            )
            if declared_current_a == 0.0 && declared_flux_wb == 0.0
                0.0
            else
                secant_inductance_h = declared_flux_wb / declared_current_a
                real(
                    branch_phasor /
                    complex(0.0, reactive_angular_frequency * secant_inductance_h),
                )
            end
        end
        push!(current_residuals, abs(initialized_currents[index] - expected_current_a))
        push!(current_scales, abs(expected_current_a))
        push!(segment_residuals, abs(initialized_segments[index] - state.segment))
    end
    tolerances = request.tolerances
    characteristic_residual = maximum(characteristic_residuals; init=0.0)
    characteristic_scale = maximum(characteristic_scales; init=0.0)
    predictor_residual = maximum(predictor_residuals; init=0.0)
    predictor_scale = maximum(predictor_scales; init=0.0)
    current_residual = maximum(current_residuals; init=0.0)
    current_scale = maximum(current_scales; init=0.0)
    segment_residual = maximum(segment_residuals; init=0.0)
    flux_allowance = tolerances.flux_absolute_wb
    current_allowance = tolerances.current_absolute_a +
        tolerances.current_relative * current_scale
    return EMTInitializationResidual[
        EMTInitializationResidual(
            :piecewise_nonlinear_inductor_state,
            :characteristic_flux_equilibrium,
            "Wb",
            characteristic_residual,
            tolerances.flux_absolute_wb,
            0.0,
            characteristic_scale,
            0.0,
            characteristic_residual / flux_allowance,
            characteristic_residual <= flux_allowance,
        ),
        EMTInitializationResidual(
            :piecewise_nonlinear_inductor_state,
            :flux_half_step_recurrence,
            "Wb",
            predictor_residual,
            tolerances.flux_absolute_wb,
            0.0,
            predictor_scale,
            0.0,
            predictor_residual / flux_allowance,
            predictor_residual <= flux_allowance,
        ),
        EMTInitializationResidual(
            :piecewise_nonlinear_inductor_state,
            :harmonic_current_equilibrium,
            "A",
            current_residual,
            tolerances.current_absolute_a,
            tolerances.current_relative,
            current_scale,
            0.0,
            current_residual / current_allowance,
            current_residual <= current_allowance,
        ),
        EMTInitializationResidual(
            :piecewise_nonlinear_inductor_state,
            :active_characteristic_segment,
            "index",
            segment_residual,
            0.0,
            0.0,
            1.0,
            0.0,
            segment_residual == 0.0 ? 0.0 : Inf,
            segment_residual == 0.0,
        ),
    ]
end

function _emt_hysteretic_initialization_residuals(
    prepared::PreparedEMTStudy,
    point::EMTInitializationFrequencyPoint,
    request::EMTInitializationRequest,
)
    step_configs = prepared.runtime_template.step_configs
    step_configs isa DynamicDeckStepConfigProvider ||
        return EMTInitializationResidual[]
    config = step_configs.nonlinear_current_config
    config === nothing && return EMTInitializationResidual[]
    nonlinear_types = Int.(config.nonlinear_types)
    indices = findall(==(HYSTERETIC_INDUCTOR_NONLINEAR_TYPE), nonlinear_types)
    isempty(indices) && return EMTInitializationResidual[]
    from_nodes = Int.(config.nonlinear_deck_from_nodes)
    to_nodes = Int.(config.nonlinear_deck_to_nodes)
    state_starts = Int.(config.nonlinear_admittance_nodes)
    time_zero_fluxes = Float64.(config.hysteretic_initial_flux_values)
    runtime_fluxes = Float64.(config.initial_runtime_voltage_values)
    currents = Float64.(config.hysteretic_initial_current_values)
    companion_currents = Float64.(config.initial_companion_current_values)
    slopes = Float64.(config.gslope)
    delta2 = Float64(config.delta2)
    flux_residuals = Float64[]
    flux_scales = Float64[]
    current_residuals = Float64[]
    current_scales = Float64[]
    for index in indices
        branch_phasor =
            _node_voltage_phasor(point.node_voltage_phasors, from_nodes[index]) -
            _node_voltage_phasor(point.node_voltage_phasors, to_nodes[index])
        branch_voltage_v = real(branch_phasor)
        reconstructed_flux_wb =
            runtime_fluxes[index] + delta2 * branch_voltage_v
        push!(
            flux_residuals,
            abs(reconstructed_flux_wb - time_zero_fluxes[index]),
        )
        push!(flux_scales, abs(time_zero_fluxes[index]))
        companion_admittance_s = slopes[state_starts[index] + 1]
        reconstructed_current_a =
            companion_currents[index] + companion_admittance_s * branch_voltage_v
        push!(
            current_residuals,
            abs(reconstructed_current_a - currents[index]),
        )
        push!(current_scales, abs(currents[index]))
    end
    tolerances = request.tolerances
    flux_residual = maximum(flux_residuals; init=0.0)
    flux_scale = maximum(flux_scales; init=0.0)
    flux_allowance = tolerances.flux_absolute_wb
    current_residual = maximum(current_residuals; init=0.0)
    current_scale = maximum(current_scales; init=0.0)
    current_allowance = tolerances.current_absolute_a +
        tolerances.current_relative * current_scale
    return EMTInitializationResidual[
        EMTInitializationResidual(
            :hysteretic_magnetic_state,
            :flux_half_step_recurrence,
            "Wb",
            flux_residual,
            tolerances.flux_absolute_wb,
            0.0,
            flux_scale,
            0.0,
            flux_residual / flux_allowance,
            flux_residual <= flux_allowance,
        ),
        EMTInitializationResidual(
            :hysteretic_magnetic_state,
            :companion_current_equilibrium,
            "A",
            current_residual,
            tolerances.current_absolute_a,
            tolerances.current_relative,
            current_scale,
            0.0,
            current_residual / current_allowance,
            current_residual <= current_allowance,
        ),
    ]
end

struct _EMTInitializationDigestWriter
    context::SHA.SHA2_256_CTX
    pending_bytes::Vector{UInt8}
    symbol_identifiers::Dict{Symbol,UInt64}
    type_identifiers::Dict{Tuple{DataType,Bool},UInt64}
end

const _EMT_INITIALIZATION_DIGEST_BUFFER_BYTES = 16 * 1024

_EMTInitializationDigestWriter() =
    _EMTInitializationDigestWriter(
        SHA.SHA2_256_CTX(),
        UInt8[],
        Dict{Symbol,UInt64}(),
        Dict{Tuple{DataType,Bool},UInt64}(),
    )

function _flush_emt_initialization_digest!(
    writer::_EMTInitializationDigestWriter,
)
    isempty(writer.pending_bytes) && return writer
    SHA.update!(writer.context, writer.pending_bytes)
    empty!(writer.pending_bytes)
    return writer
end

function _update_emt_initialization_digest!(
    writer::_EMTInitializationDigestWriter,
    bytes::NTuple{N,UInt8},
) where {N}
    pending_bytes = writer.pending_bytes
    for byte in bytes
        push!(pending_bytes, byte)
    end
    length(pending_bytes) >=
        _EMT_INITIALIZATION_DIGEST_BUFFER_BYTES &&
        _flush_emt_initialization_digest!(writer)
    return writer
end

function _update_emt_initialization_digest!(
    writer::_EMTInitializationDigestWriter,
    bytes::AbstractVector{UInt8},
)
    if length(bytes) >= _EMT_INITIALIZATION_DIGEST_BUFFER_BYTES
        _flush_emt_initialization_digest!(writer)
        SHA.update!(writer.context, bytes)
    else
        append!(writer.pending_bytes, bytes)
        length(writer.pending_bytes) >=
            _EMT_INITIALIZATION_DIGEST_BUFFER_BYTES &&
            _flush_emt_initialization_digest!(writer)
    end
    return writer
end

function _write_emt_initialization_digest_tag!(
    writer::_EMTInitializationDigestWriter,
    tag::UInt8,
)
    return _update_emt_initialization_digest!(writer, (tag,))
end

function _write_emt_initialization_digest_uint64!(
    writer::_EMTInitializationDigestWriter,
    value::UInt64,
)
    bytes = ntuple(
        index -> UInt8((value >> (8 * (index - 1))) & 0xff),
        Val(8),
    )
    return _update_emt_initialization_digest!(writer, bytes)
end

function _write_emt_initialization_digest_text!(
    writer::_EMTInitializationDigestWriter,
    value::AbstractString,
)
    bytes = codeunits(value)
    _write_emt_initialization_digest_uint64!(writer, UInt64(length(bytes)))
    return _update_emt_initialization_digest!(writer, bytes)
end

function _write_emt_initialization_digest_symbol!(
    writer::_EMTInitializationDigestWriter,
    value::Symbol,
)
    identifier = get(writer.symbol_identifiers, value, UInt64(0))
    if identifier != 0
        _write_emt_initialization_digest_tag!(writer, 0x00)
        return _write_emt_initialization_digest_uint64!(writer, identifier)
    end
    identifier = UInt64(length(writer.symbol_identifiers) + 1)
    writer.symbol_identifiers[value] = identifier
    _write_emt_initialization_digest_tag!(writer, 0x01)
    _write_emt_initialization_digest_uint64!(writer, identifier)
    return _write_emt_initialization_digest_text!(writer, String(value))
end

function _write_emt_initialization_digest_type!(
    writer::_EMTInitializationDigestWriter,
    value::DataType;
    include_parameters::Bool=true,
)
    key = (value, include_parameters)
    identifier = get(writer.type_identifiers, key, UInt64(0))
    if identifier != 0
        _write_emt_initialization_digest_tag!(writer, 0x00)
        return _write_emt_initialization_digest_uint64!(writer, identifier)
    end
    identifier = UInt64(length(writer.type_identifiers) + 1)
    writer.type_identifiers[key] = identifier
    _write_emt_initialization_digest_tag!(writer, 0x01)
    _write_emt_initialization_digest_uint64!(writer, identifier)
    _write_emt_initialization_digest_text!(writer, string(parentmodule(value)))
    _write_emt_initialization_digest_symbol!(writer, nameof(value))
    include_parameters || return writer
    _write_emt_initialization_digest_uint64!(writer, UInt64(length(value.parameters)))
    for parameter in value.parameters
        if parameter isa DataType
            _write_emt_initialization_digest_type!(
                writer,
                parameter;
                include_parameters=false,
            )
        elseif parameter isa UnionAll
            _write_emt_initialization_digest_text!(writer, string(parameter))
        else
            _write_canonical_emt_initialization_state!(
                writer,
                parameter,
                IdDict{Any,Nothing}(),
            )
        end
    end
    return writer
end

function _emt_initialization_bulk_digest_eltype(::Type{T}) where {T}
    return T <: Union{
        Bool,
        Int8,
        Int16,
        Int32,
        Int64,
        Int128,
        UInt8,
        UInt16,
        UInt32,
        UInt64,
        UInt128,
        Float16,
        Float32,
        Float64,
        Complex{Float16},
        Complex{Float32},
        Complex{Float64},
    }
end

function _write_emt_initialization_digest_array_values!(
    writer::_EMTInitializationDigestWriter,
    value::AbstractArray,
    active_objects::IdDict{Any,Nothing},
)
    if value isa Array &&
       _emt_initialization_bulk_digest_eltype(eltype(value)) &&
       ENDIAN_BOM == 0x04030201
        _write_emt_initialization_digest_tag!(writer, 0x01)
        _update_emt_initialization_digest!(
            writer,
            reinterpret(UInt8, vec(value)),
        )
        return writer
    end
    _write_emt_initialization_digest_tag!(writer, 0x00)
    for index in eachindex(value)
        isassigned(value, index) || throw(ArgumentError(
            "deterministic EMT initialization state cannot contain an unassigned array entry",
        ))
        _write_canonical_emt_initialization_state!(
            writer,
            value[index],
            active_objects,
        )
    end
    return writer
end

function _write_emt_initialization_tuple_values_loop!(
    writer::_EMTInitializationDigestWriter,
    value::Tuple,
    active_objects::IdDict{Any,Nothing},
)
    for entry in value
        _write_canonical_emt_initialization_state!(
            writer,
            entry,
            active_objects,
        )
    end
    return writer
end

function _write_emt_initialization_homogeneous_tuple_values!(
    writer::_EMTInitializationDigestWriter,
    value::Tuple{Vararg{T}},
    active_objects::IdDict{Any,Nothing},
) where {T}
    for index in eachindex(value)
        _write_canonical_emt_initialization_state!(
            writer,
            value[index],
            active_objects,
        )
    end
    return writer
end

@generated function _write_emt_initialization_tuple_values!(
    writer::_EMTInitializationDigestWriter,
    value::T,
    active_objects::IdDict{Any,Nothing},
) where {T<:Tuple}
    field_count = fieldcount(T)
    if field_count > 32
        field_types = fieldtypes(T)
        if all(==(first(field_types)), field_types)
            return :(_write_emt_initialization_homogeneous_tuple_values!(
                writer,
                value,
                active_objects,
            ))
        elseif all(==(field_types[2]), field_types[2:end])
            return quote
                _write_canonical_emt_initialization_state!(
                    writer,
                    getfield(value, 1),
                    active_objects,
                )
                _write_emt_initialization_homogeneous_tuple_values!(
                    writer,
                    Base.tail(value),
                    active_objects,
                )
            end
        end
        return :(_write_emt_initialization_tuple_values_loop!(
            writer,
            value,
            active_objects,
        ))
    end
    expressions = [
        :(_write_canonical_emt_initialization_state!(
            writer,
            getfield(value, $index),
            active_objects,
        ))
        for index in 1:field_count
    ]
    return quote
        $(expressions...)
        writer
    end
end

@generated function _write_emt_initialization_named_tuple_fields!(
    writer::_EMTInitializationDigestWriter,
    value::T,
    active_objects::IdDict{Any,Nothing},
) where {T<:NamedTuple}
    names = fieldnames(T)
    expressions = Expr[]
    for (index, name) in enumerate(names)
        push!(
            expressions,
            :(_write_emt_initialization_digest_symbol!(
                writer,
                $(QuoteNode(name)),
            )),
            :(_write_canonical_emt_initialization_state!(
                writer,
                getfield(value, $index),
                active_objects,
            )),
        )
    end
    return quote
        $(expressions...)
        writer
    end
end

@generated function _write_emt_initialization_struct_fields!(
    writer::_EMTInitializationDigestWriter,
    value::T,
    active_objects::IdDict{Any,Nothing},
) where {T}
    names = fieldnames(T)
    expressions = Expr[]
    for (index, name) in enumerate(names)
        push!(
            expressions,
            :(_write_emt_initialization_digest_symbol!(
                writer,
                $(QuoteNode(name)),
            )),
            :(_write_canonical_emt_initialization_state!(
                writer,
                getfield(value, $index),
                active_objects,
            )),
        )
    end
    return quote
        $(expressions...)
        writer
    end
end

@generated function _write_emt_initialization_record_values!(
    writer::_EMTInitializationDigestWriter,
    value::T,
    active_objects::IdDict{Any,Nothing},
) where {T}
    expressions = [
        :(_write_canonical_emt_initialization_state!(
            writer,
            getfield(value, $index),
            active_objects,
        ))
        for index in 1:fieldcount(T)
    ]
    return quote
        $(expressions...)
        writer
    end
end

function _write_emt_initialization_record_schema!(
    writer::_EMTInitializationDigestWriter,
    ::Type{T},
) where {T}
    names = fieldnames(T)
    types = fieldtypes(T)
    _write_emt_initialization_digest_uint64!(writer, UInt64(length(names)))
    for (name, field_type) in zip(names, types)
        _write_emt_initialization_digest_symbol!(writer, name)
        if field_type isa DataType
            _write_emt_initialization_digest_type!(writer, field_type)
        else
            _write_emt_initialization_digest_text!(writer, string(field_type))
        end
    end
    return writer
end

function _write_emt_initialization_nodal_element_batch!(
    writer::_EMTInitializationDigestWriter,
    batch::AbstractVector{T},
    active_objects::IdDict{Any,Nothing},
) where {T}
    _write_emt_initialization_digest_type!(writer, T)
    _write_emt_initialization_digest_uint64!(writer, UInt64(length(batch)))
    compact_records = isbitstype(T) && fieldcount(T) > 0
    _write_emt_initialization_digest_tag!(
        writer,
        compact_records ? 0x01 : 0x00,
    )
    if compact_records
        _write_emt_initialization_record_schema!(writer, T)
        for element in batch
            _write_emt_initialization_record_values!(
                writer,
                element,
                active_objects,
            )
        end
    else
        for element in batch
            _write_canonical_emt_initialization_state!(
                writer,
                element,
                active_objects,
            )
        end
    end
    return writer
end

_write_emt_initialization_nodal_element_batches!(
    writer::_EMTInitializationDigestWriter,
    ::Tuple{},
    _active_objects::IdDict{Any,Nothing},
) = writer

function _write_emt_initialization_nodal_element_batches!(
    writer::_EMTInitializationDigestWriter,
    batches::Tuple,
    active_objects::IdDict{Any,Nothing},
)
    _write_emt_initialization_nodal_element_batch!(
        writer,
        first(batches),
        active_objects,
    )
    return _write_emt_initialization_nodal_element_batches!(
        writer,
        Base.tail(batches),
        active_objects,
    )
end

function _write_emt_initialization_nodal_elements!(
    writer::_EMTInitializationDigestWriter,
    elements::NodalElementSequence,
    active_objects::IdDict{Any,Nothing},
)
    batches = elements.contiguous_type_batches
    sum(length, batches; init=0) == length(elements) || throw(ArgumentError(
        "deterministic EMT initialization state found an inconsistent nodal element sequence",
    ))
    _write_emt_initialization_digest_tag!(writer, 0x19)
    _write_emt_initialization_digest_uint64!(writer, UInt64(length(elements)))
    _write_emt_initialization_digest_uint64!(writer, UInt64(length(batches)))
    return _write_emt_initialization_nodal_element_batches!(
        writer,
        batches,
        active_objects,
    )
end

function _emt_initialization_arrays_bitwise_equal(
    first::Array{T},
    second::Array{T},
) where {T}
    size(first) == size(second) || return false
    isbitstype(T) || return false
    return reinterpret(UInt8, vec(first)) == reinterpret(UInt8, vec(second))
end

@generated function _write_emt_initialization_nodal_system_fields!(
    writer::_EMTInitializationDigestWriter,
    system::T,
    active_objects::IdDict{Any,Nothing},
) where {T<:NodalSystem}
    names = fieldnames(T)
    admittance_index = findfirst(==(:y), names)
    factor_index = findfirst(==(:y_factor), names)
    if isnothing(admittance_index) || isnothing(factor_index)
        return :(throw(ArgumentError(
            "deterministic EMT initialization requires nodal admittance and factor workspaces",
        )))
    end
    expressions = Expr[]
    for (index, name) in enumerate(names)
        push!(
            expressions,
            :(_write_emt_initialization_digest_symbol!(
                writer,
                $(QuoteNode(name)),
            )),
        )
        if index == factor_index
            push!(expressions, quote
                factor_matches_admittance =
                    _emt_initialization_arrays_bitwise_equal(
                        getfield(system, $admittance_index),
                        getfield(system, $factor_index),
                    )
                _write_emt_initialization_digest_tag!(
                    writer,
                    factor_matches_admittance ? 0x01 : 0x00,
                )
                factor_matches_admittance ||
                    _write_canonical_emt_initialization_state!(
                        writer,
                        getfield(system, $factor_index),
                        active_objects,
                    )
            end)
        else
            push!(expressions, :(_write_canonical_emt_initialization_state!(
                writer,
                getfield(system, $index),
                active_objects,
            )))
        end
    end
    return quote
        $(expressions...)
        writer
    end
end

function _write_emt_initialization_nodal_system!(
    writer::_EMTInitializationDigestWriter,
    system::NodalSystem,
    active_objects::IdDict{Any,Nothing},
)
    _write_emt_initialization_digest_tag!(writer, 0x18)
    _write_emt_initialization_digest_type!(writer, typeof(system))
    _write_emt_initialization_digest_uint64!(
        writer,
        UInt64(fieldcount(typeof(system))),
    )
    return _write_emt_initialization_nodal_system_fields!(
        writer,
        system,
        active_objects,
    )
end

function _emt_initialization_state_sort_key(value)
    return Tuple(_emt_initialization_state_digest(value))
end

function _write_canonical_emt_initialization_state!(
    writer::_EMTInitializationDigestWriter,
    value,
    active_objects::IdDict{Any,Nothing},
)
    if value === nothing
        return _write_emt_initialization_digest_tag!(writer, 0x00)
    elseif value === missing
        return _write_emt_initialization_digest_tag!(writer, 0x01)
    elseif value isa Bool
        _write_emt_initialization_digest_tag!(writer, 0x02)
        return _write_emt_initialization_digest_tag!(writer, value ? 0x01 : 0x00)
    elseif value isa Signed && sizeof(value) <= 8
        _write_emt_initialization_digest_tag!(writer, 0x03)
        return _write_emt_initialization_digest_uint64!(
            writer,
            reinterpret(UInt64, Int64(value)),
        )
    elseif value isa Unsigned && sizeof(value) <= 8
        _write_emt_initialization_digest_tag!(writer, 0x04)
        return _write_emt_initialization_digest_uint64!(writer, UInt64(value))
    elseif value isa Float64
        _write_emt_initialization_digest_tag!(writer, 0x05)
        return _write_emt_initialization_digest_uint64!(
            writer,
            reinterpret(UInt64, value),
        )
    elseif value isa Float32
        _write_emt_initialization_digest_tag!(writer, 0x06)
        return _write_emt_initialization_digest_uint64!(
            writer,
            UInt64(reinterpret(UInt32, value)),
        )
    elseif value isa Float16
        _write_emt_initialization_digest_tag!(writer, 0x07)
        return _write_emt_initialization_digest_uint64!(
            writer,
            UInt64(reinterpret(UInt16, value)),
        )
    elseif value isa AbstractFloat
        _write_emt_initialization_digest_tag!(writer, 0x08)
        return _write_emt_initialization_digest_text!(writer, repr(value))
    elseif value isa Char
        _write_emt_initialization_digest_tag!(writer, 0x09)
        return _write_emt_initialization_digest_uint64!(writer, UInt64(value))
    elseif value isa Symbol
        _write_emt_initialization_digest_tag!(writer, 0x0a)
        return _write_emt_initialization_digest_symbol!(writer, value)
    elseif value isa AbstractString
        _write_emt_initialization_digest_tag!(writer, 0x0b)
        return _write_emt_initialization_digest_text!(writer, value)
    elseif value isa Enum
        _write_emt_initialization_digest_tag!(writer, 0x0c)
        _write_emt_initialization_digest_type!(writer, typeof(value))
        return _write_emt_initialization_digest_uint64!(
            writer,
            UInt64(Integer(value)),
        )
    elseif value isa DataType
        _write_emt_initialization_digest_tag!(writer, 0x0d)
        return _write_emt_initialization_digest_type!(writer, value)
    elseif value isa Type
        _write_emt_initialization_digest_tag!(writer, 0x17)
        return _write_emt_initialization_digest_text!(writer, string(value))
    elseif value isa Module
        _write_emt_initialization_digest_tag!(writer, 0x0e)
        return _write_emt_initialization_digest_text!(writer, string(value))
    elseif value isa Complex
        _write_emt_initialization_digest_tag!(writer, 0x0f)
        _write_canonical_emt_initialization_state!(writer, real(value), active_objects)
        return _write_canonical_emt_initialization_state!(
            writer,
            imag(value),
            active_objects,
        )
    elseif value isa Ptr
        value == C_NULL || throw(ArgumentError(
            "deterministic EMT initialization state cannot contain a live pointer",
        ))
        return _write_emt_initialization_digest_tag!(writer, 0x10)
    end

    tracked = ismutabletype(typeof(value))
    if tracked
        haskey(active_objects, value) && throw(ArgumentError(
            "deterministic EMT initialization state cannot contain a reference cycle",
        ))
        active_objects[value] = nothing
    end
    try
        if value isa NodalSystem
            _write_emt_initialization_nodal_system!(
                writer,
                value,
                active_objects,
            )
        elseif value isa NodalElementSequence
            _write_emt_initialization_nodal_elements!(
                writer,
                value,
                active_objects,
            )
        elseif value isa NamedTuple
            _write_emt_initialization_digest_tag!(writer, 0x11)
            names = keys(value)
            _write_emt_initialization_digest_uint64!(writer, UInt64(length(names)))
            _write_emt_initialization_named_tuple_fields!(
                writer,
                value,
                active_objects,
            )
        elseif value isa Tuple
            _write_emt_initialization_digest_tag!(writer, 0x12)
            _write_emt_initialization_digest_uint64!(writer, UInt64(length(value)))
            _write_emt_initialization_tuple_values!(
                writer,
                value,
                active_objects,
            )
        elseif value isa AbstractArray
            _write_emt_initialization_digest_tag!(writer, 0x13)
            element_type = eltype(value)
            if element_type isa DataType
                _write_emt_initialization_digest_type!(
                    writer,
                    element_type;
                    include_parameters=isempty(value),
                )
            else
                _write_emt_initialization_digest_text!(writer, string(element_type))
            end
            _write_emt_initialization_digest_uint64!(writer, UInt64(ndims(value)))
            for dimension in size(value)
                _write_emt_initialization_digest_uint64!(writer, UInt64(dimension))
            end
            _write_emt_initialization_digest_array_values!(
                writer,
                value,
                active_objects,
            )
        elseif value isa AbstractDict
            _write_emt_initialization_digest_tag!(writer, 0x14)
            _write_emt_initialization_digest_uint64!(writer, UInt64(length(value)))
            if keytype(value) === Symbol
                ordered_keys = sort!(collect(keys(value)))
                for key in ordered_keys
                    _write_canonical_emt_initialization_state!(
                        writer,
                        key,
                        active_objects,
                    )
                    _write_canonical_emt_initialization_state!(
                        writer,
                        value[key],
                        active_objects,
                    )
                end
            else
                entries = [
                    (_emt_initialization_state_sort_key(key), key, entry)
                    for (key, entry) in value
                ]
                sort!(entries; by=first)
                for (_, key, entry) in entries
                    _write_canonical_emt_initialization_state!(
                        writer,
                        key,
                        active_objects,
                    )
                    _write_canonical_emt_initialization_state!(
                        writer,
                        entry,
                        active_objects,
                    )
                end
            end
        elseif value isa AbstractSet
            _write_emt_initialization_digest_tag!(writer, 0x15)
            entries = [
                (_emt_initialization_state_sort_key(entry), entry)
                for entry in value
            ]
            sort!(entries; by=first)
            _write_emt_initialization_digest_uint64!(writer, UInt64(length(entries)))
            for (_, entry) in entries
                _write_canonical_emt_initialization_state!(writer, entry, active_objects)
            end
        else
            _write_emt_initialization_digest_tag!(writer, 0x16)
            _write_emt_initialization_digest_type!(writer, typeof(value))
            names = fieldnames(typeof(value))
            _write_emt_initialization_digest_uint64!(writer, UInt64(length(names)))
            if isempty(names)
                _write_emt_initialization_digest_text!(writer, repr(value))
            else
                _write_emt_initialization_struct_fields!(
                    writer,
                    value,
                    active_objects,
                )
            end
        end
    finally
        tracked && delete!(active_objects, value)
    end
    return writer
end

function _emt_initialization_state_digest(value)
    writer = _EMTInitializationDigestWriter()
    _write_canonical_emt_initialization_state!(
        writer,
        value,
        IdDict{Any,Nothing}(),
    )
    _flush_emt_initialization_digest!(writer)
    return SHA.digest!(writer.context)
end

function _emt_initialization_state_signature(
    prepared::Union{
        PreparedEMTStudy,
        PreparedMachineEMTStudy,
        PreparedAverageValueGridFollowingEMTStudy,
    },
    request::EMTInitializationRequest,
    mappings::Vector{OperatingPointMappingRecord},
)
    io = IOBuffer()
    write(io, "aimora.emt.initialization.state.v4\n")
    for value in (
        request.project_signature,
        request.settings_signature,
        request.model_signature,
        String(_emt_harmonic_formulation_symbol(request.formulation)),
        repr(request.frequency_hz),
        repr(request.time_origin_s),
    )
        write(io, value, '\n')
    end
    if request.operating_point isa EMTOperatingPoint
        write(
            io,
            String(request.operating_point.source_representation),
            '\n',
            request.operating_point.source_state_signature,
            '\n',
        )
    end
    accepted_state = if prepared isa PreparedEMTStudy
        prepared.runtime_template
    elseif prepared isa PreparedMachineEMTStudy
        (
            machine_family=prepared.machine_family,
            initialization_state=prepared.initialization_state,
        )
    else
        (
            network=prepared.network.runtime_template,
            converter_equilibrium=prepared.equilibrium,
            converter_ownership=prepared.converter,
        )
    end
    write(
        io,
        "complete_accepted_state_sha256=",
        bytes2hex(_emt_initialization_state_digest(accepted_state)),
        '\n',
    )
    for mapping in mappings
        write(
            io,
            String(mapping.asset),
            '/',
            String(mapping.quantity),
            '/',
            String(mapping.phase),
            '=',
            bitstring(real(mapping.target_value_si)),
            ',',
            bitstring(imag(mapping.target_value_si)),
            '\n',
        )
    end
    return bytes2hex(sha256(take!(io)))
end

function _emt_machine_frequency_point(
    sample,
    request::EMTInitializationRequest,
    ;
    frequency_assignment::Symbol=:model_owned_machine_equilibrium,
)
    node_voltage_phasors = ComplexF64.(sample.node_voltage_phasors)
    node_count = length(node_voltage_phasors)
    frequency_hz = Float64(sample.steady_state_frequency_hz)
    isapprox(
        frequency_hz,
        request.frequency_hz;
        atol=1.0e-12,
        rtol=1.0e-12,
    ) || throw(_EMTInitializationRefusal(
        :source_frequency_mismatch,
        :machine_equilibrium,
        :frequency_hz,
        "machine-owned steady-state frequency does not match the initialization request",
        (requested_frequency_hz=request.frequency_hz, machine_frequency_hz=frequency_hz),
    ))
    diagnostics = sample.topology_diagnostics
    topology = _emt_initialization_topology_report(diagnostics, node_count)
    topology.classification === :unique || throw(
        _emt_initialization_classification_failure(topology),
    )
    topology.condition_estimate <= request.tolerances.maximum_condition_estimate ||
        throw(_EMTInitializationRefusal(
            :ill_conditioned,
            :network_topology,
            :nodal_voltage,
            "machine-owned harmonic network exceeds the requested condition limit",
            (
                condition_estimate=topology.condition_estimate,
                maximum_condition_estimate=
                    request.tolerances.maximum_condition_estimate,
            ),
        ))
    admittance = ComplexF64.(sample.steady_state_admittance)
    source_injections = ComplexF64.(
        sample.steady_state_source_injection_phasors,
    )
    size(admittance) == (node_count, node_count) || throw(
        _EMTInitializationRefusal(
            :incomplete_model_state,
            :machine_network_equilibrium,
            :admittance,
            "machine-owned steady-state sample omitted its complete nodal admittance",
            (; node_count, admittance_size=size(admittance)),
        ),
    )
    length(source_injections) == node_count || throw(
        _EMTInitializationRefusal(
            :incomplete_model_state,
            :machine_network_equilibrium,
            :source_injection,
            "machine-owned steady-state sample omitted its complete source injection",
            (; node_count, source_injection_count=length(source_injections)),
        ),
    )
    symmetry_error = maximum(abs, admittance - transpose(admittance); init=0.0)
    dissipative = Hermitian(0.5 .* (admittance + adjoint(admittance)))
    minimum_dissipative_eigenvalue = minimum(eigvals(dissipative); init=0.0)
    node_frequencies = hasproperty(sample, :node_steady_state_frequencies_hz) ?
        Float64.(sample.node_steady_state_frequencies_hz) :
        fill(frequency_hz, node_count)
    source_rows = hasproperty(sample, :node_frequency_source_row_indices) ?
        Int.(sample.node_frequency_source_row_indices) : zeros(Int, node_count)
    if length(node_frequencies) < node_count
        append!(
            node_frequencies,
            fill(frequency_hz, node_count - length(node_frequencies)),
        )
    end
    if length(source_rows) < node_count
        append!(source_rows, zeros(Int, node_count - length(source_rows)))
    end
    length(node_frequencies) == node_count && length(source_rows) == node_count ||
        throw(_EMTInitializationRefusal(
            :incomplete_model_state,
            :machine_network_equilibrium,
            :frequency_partition,
            "steady-state frequency ownership does not cover the augmented network",
            (;
                node_count,
                frequency_count=length(node_frequencies),
                source_owner_count=length(source_rows),
            ),
        ))
    successors = hasproperty(sample, :source_frequency_successor_indices) ?
        Int.(sample.source_frequency_successor_indices) : Int[]
    subnetwork_count = hasproperty(sample, :steady_state_frequency_subnetwork_count) ?
        Int(sample.steady_state_frequency_subnetwork_count) : 1
    return EMTInitializationFrequencyPoint(
        :physical_frequency,
        frequency_assignment,
        frequency_hz,
        2.0 * pi * frequency_hz,
        node_frequencies,
        source_rows,
        successors,
        subnetwork_count,
        node_voltage_phasors,
        source_injections,
        zeros(ComplexF64, node_count),
        topology,
        symmetry_error,
        minimum_dissipative_eigenvalue,
        true,
    )
end

function _emt_model_initialization_residual(
    owner::Symbol,
    quantity::Symbol,
    unit::AbstractString,
    value::Real,
    absolute_tolerance::Real,
    relative_tolerance::Real,
    reference_scale::Real,
)
    residual = Float64(value)
    absolute = Float64(absolute_tolerance)
    relative = Float64(relative_tolerance)
    scale = Float64(reference_scale)
    allowance = absolute + relative * scale
    scaled = allowance > 0.0 ? residual / allowance :
        (residual == 0.0 ? 0.0 : Inf)
    return EMTInitializationResidual(
        owner,
        quantity,
        String(unit),
        residual,
        absolute,
        relative,
        scale,
        0.0,
        scaled,
        isfinite(residual) && residual <= allowance,
    )
end

function _emt_machine_network_residual(
    point::EMTInitializationFrequencyPoint,
    request::EMTInitializationRequest,
)
    tolerance = request.tolerances.current_relative
    return _emt_model_initialization_residual(
        :machine_network_equilibrium,
        :scaled_nodal_backward_error,
        "1",
        point.topology.relative_residual,
        tolerance,
        0.0,
        1.0,
    )
end

function _emt_validate_machine_initialization_request(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
)
    request.formulation isa PhysicalFrequencyFormulation || throw(
        _EMTInitializationRefusal(
            :unsupported_formulation,
            :machine_equilibrium,
            :harmonic_formulation,
            "the admitted machine initializer owns a physical-frequency equilibrium; timestep-matched machine recurrence has not been proven",
            (requested_formulation=_emt_harmonic_formulation_symbol(request.formulation),),
        ),
    )
    all(frequency -> isapprox(
        frequency,
        request.frequency_hz;
        atol=1.0e-12,
        rtol=1.0e-12,
    ), request.frequency_grid_hz) || throw(_EMTInitializationRefusal(
        :unsupported_frequency_scan,
        :machine_equilibrium,
        :frequency_grid_hz,
        "model-owned machine initialization accepts one operating frequency; whole-network scan points must be requested from the network scan owner",
        (frequency_grid_hz=copy(request.frequency_grid_hz),),
    ))
    request.operating_point === nothing || throw(_EMTInitializationRefusal(
        :unsupported_operating_point_mapping,
        :machine_equilibrium,
        :operating_point,
        "machine operating-point import requires machine-owned current, torque, flux, and control quantities rather than voltage-only constraints",
        (operating_point_type=string(typeof(request.operating_point)),),
    ))
    request.time_origin_s == 0.0 || throw(_EMTInitializationRefusal(
        :unsupported_time_origin,
        :machine_equilibrium,
        :time_origin_s,
        "machine-owned initialization currently requires the deck time origin",
        (time_origin_s=request.time_origin_s,),
    ))
    _validate_emt_initialization_source_frequency(parsed, request)
    return request
end

function _emt_machine_state_inventory(
    parsed::DeckParser.DeckParseResult,
    machine_family::Symbol,
    electrical_state_count::Int,
    mechanical_state_count::Int,
)
    records = EMTInitializationStateRecord[]
    _append_emt_initialization_state!(
        records,
        :network_voltage,
        :algebraic,
        maximum(values(parsed.node_map); init=0),
        :model_owned_harmonic_network_equilibrium,
    )
    _append_emt_initialization_state!(
        records,
        :network_topology,
        :discrete,
        1,
        :ranked_connected_component_classification,
    )
    _append_emt_initialization_state!(
        records,
        Symbol(machine_family, :_electrical_state),
        :continuous,
        electrical_state_count,
        :model_owned_machine_equilibrium,
    )
    _append_emt_initialization_state!(
        records,
        Symbol(machine_family, :_mechanical_state),
        :continuous,
        mechanical_state_count,
        :model_owned_torque_speed_equilibrium,
    )
    switch_count = length(DeckParser.deck_over5_switch_rows(parsed)) +
        length(DeckParser.deck_control_system_switch_coupling_rows(parsed))
    _append_emt_initialization_state!(
        records,
        :switch_mode,
        :discrete,
        switch_count,
        :declared_initial_topology,
    )
    _append_emt_initialization_state!(
        records,
        :switch_event_state,
        :scheduler,
        length(DeckParser.deck_over5_switch_rows(parsed)),
        :initial_event_surface_classification,
    )
    control_count = length(
        DeckParser.deck_synchronous_machine_control_interface_rows(parsed),
    )
    _append_emt_initialization_state!(
        records,
        :machine_control_state,
        :discrete,
        control_count,
        :declared_control_equilibrium,
    )
    transformer_count = get(
        parsed.card_counts,
        :fixed_card_saturated_transformer_intake,
        0,
    )
    _append_emt_initialization_state!(
        records,
        :machine_terminal_transformer_state,
        :continuous,
        transformer_count,
        :coupled_transformer_branch_equilibrium,
    )
    _append_emt_initialization_state!(
        records,
        :output_cursor,
        :discrete,
        1,
        :time_zero_output_epoch,
    )
    sort!(records; by=record -> (String(record.state_family), String(record.owner)))
    return records
end

function _emt_synchronous_machine_preparation(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
)
    machine_indices = sort!(unique(
        row.machine_index for row in
        DeckParser.deck_synchronous_machine_terminal_voltage_rows(parsed)
    ))
    length(machine_indices) == 1 || throw(_EMTInitializationRefusal(
        :unsupported_machine_fleet_initialization,
        :synchronous_machine_state,
        :machine_count,
        "the cohesive synchronous-machine initializer currently requires one admitted machine",
        (machine_indices=copy(machine_indices),),
    ))
    machine_index = only(machine_indices)
    timestep_s = DeckParser.deck_fixed_time_horizon_options(parsed).dt_s
    context = _deck_synchronous_machine_runtime_context(
        parsed,
        timestep_s,
        timestep_s;
        saturated_transformer_branch_runtime_enabled=true,
        coupled_lumped_sequence_history_enabled=true,
        recorded_step_indices=[0],
    )
    sample = _deck_synchronous_machine_network_initial_sample(
        parsed,
        context;
        strict_topology_classification=true,
    )
    sample === nothing && throw(_EMTInitializationRefusal(
        :missing_model_state,
        :synchronous_machine_state,
        :terminal_voltage,
        "synchronous-machine terminal equilibrium is missing",
        (machine_index=machine_index,),
    ))
    point = _emt_machine_frequency_point(sample, request)
    initialization = _deck_synchronous_machine_initial_state(
        parsed,
        context,
        sample;
        machine_index,
    )
    horizon = run_deck_synchronous_machine_horizon(
        parsed,
        deepcopy(initialization.state);
        numask=initialization.numask,
        nlocg=initialization.nlocg,
        nloce=initialization.nloce,
        time_step_s=timestep_s,
        dynamic_step_count=1,
        angle_half_step_inverse=initialization.angle_half_step_inverse,
        speed_tolerance=initialization.speed_tolerance,
        omega_tolerance=initialization.omega_tolerance,
        speed_floor=initialization.speed_floor,
        max_iterations=initialization.max_iterations,
        damping_ratio=initialization.damping_ratio,
        rotor_angle_extrapolation_interval=
            initialization.rotor_angle_extrapolation_interval,
        speed_voltage_factor=initialization.speed_voltage_factor,
        electrical_speed_rad_s=initialization.electrical_speed_rad_s,
        electrical_angle_increment=initialization.electrical_angle_increment,
        saturated_transformer_branch_runtime_enabled=true,
        coupled_lumped_sequence_history_enabled=true,
        recorded_step_indices=[0, 1],
    )
    terminal_rows = sort!(
        [
            row for row in DeckParser.deck_synchronous_machine_terminal_voltage_rows(parsed)
            if row.machine_index == machine_index
        ];
        by=row -> row.phase_index,
    )
    terminal_nodes = Int[row.terminal_node_value for row in terminal_rows]
    tolerances = request.tolerances
    current_error = maximum(
        abs,
        initialization.phase_current_phasors .-
            sample.node_current_phasors[terminal_nodes];
        init=0.0,
    )
    current_scale = maximum(abs, initialization.phase_current_phasors; init=0.0)
    voltage_error = maximum(
        abs,
        horizon.terminal_voltage_values[:, 1] .-
            real.(sample.node_voltage_phasors[terminal_nodes]);
        init=0.0,
    )
    voltage_scale = maximum(
        abs,
        sample.node_voltage_phasors[terminal_nodes];
        init=0.0,
    )
    terminal_kcl = abs(sum(horizon.terminal_current_values[:, 1]))
    residuals = EMTInitializationResidual[
        _emt_machine_network_residual(point, request),
        _emt_model_initialization_residual(
            :synchronous_machine_state,
            :terminal_current_equilibrium,
            "A",
            current_error,
            tolerances.current_absolute_a,
            tolerances.current_relative,
            current_scale,
        ),
        _emt_model_initialization_residual(
            :synchronous_machine_state,
            :terminal_voltage_equilibrium,
            "V",
            voltage_error,
            tolerances.voltage_absolute_v,
            tolerances.voltage_relative,
            voltage_scale,
        ),
        _emt_model_initialization_residual(
            :synchronous_machine_state,
            :terminal_zero_sequence_kcl,
            "A",
            terminal_kcl,
            tolerances.current_absolute_a,
            tolerances.current_relative,
            maximum(abs, horizon.terminal_current_values[:, 1]; init=0.0),
        ),
    ]
    maximum_scaled = maximum(residual.scaled_value for residual in residuals)
    metric = NoArtificialTransientMetric(
        :synchronous_machine_time_zero_equilibrium,
        "1",
        request.time_origin_s,
        request.time_origin_s,
        maximum_scaled <= 1.0 ? 0.0 : maximum_scaled - 1.0,
        maximum_scaled <= 1.0 ? 0.0 : maximum_scaled - 1.0,
        0.0,
        0.0,
        request.tolerances.no_artificial_transient_normalized_rms,
        all(residual -> residual.passed, residuals),
    )
    accepted_state = (
        machine_initialization=deepcopy(initialization),
        network_voltage_phasors=copy(sample.node_voltage_phasors),
        terminal_voltage_values=copy(horizon.terminal_voltage_values[:, 1]),
        terminal_current_values=copy(horizon.terminal_current_values[:, 1]),
        machine_output_values=copy(horizon.machine_output_values[:, 1]),
        mechanical_history_values=copy(horizon.mechanical_history_values[:, 1]),
        control_output_names=copy(horizon.control_output_names),
        control_output_values=copy(horizon.control_output_values[:, 1]),
        switch_node_groups=copy.(point.topology.switch_node_groups),
        output_cursor=0,
    )
    prepared = PreparedMachineEMTStudy(
        :synchronous_machine,
        accepted_state,
        horizon,
        parsed,
    )
    inventory = _emt_machine_state_inventory(
        parsed,
        :synchronous_machine,
        length(initialization.state.current_history),
        length(initialization.state.equation_state.histq_values),
    )
    warnings = sample.time_zero_ground_fault ? String[
        "The first advance contains the deck-declared switch transition; no-artificial-transient acceptance is therefore evaluated at the complete time-zero machine/network equilibrium.",
    ] : String[]
    return prepared, point, residuals, inventory, metric, warnings
end

function _emt_universal_machine_preparation(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
)
    machine_indices = sort!(unique(
        row.machine_index for row in
        DeckParser.deck_universal_machine_definition_rows(parsed)
        if row.card_index == 1
    ))
    length(machine_indices) == 1 || throw(_EMTInitializationRefusal(
        :unsupported_machine_fleet_initialization,
        :universal_machine_state,
        :machine_count,
        "the cohesive universal-machine initializer currently requires one admitted machine",
        (machine_indices=copy(machine_indices),),
    ))
    machine_index = only(machine_indices)
    card = _deck_universal_machine_definition(parsed, machine_index, 1)
    card.machine_type in 3:12 || throw(_EMTInitializationRefusal(
        :unsupported_machine_family,
        :universal_machine_state,
        :machine_type,
        "universal wound-field synchronous types use a separate admitted initialization owner",
        (machine_type=card.machine_type,),
    ))
    automatic_direct = card.machine_type in 8:12 &&
        _deck_universal_machine_initialization_mode(parsed) == :automatic ?
        deck_direct_current_machine_automatic_initialization(
            parsed;
            machine_index,
        ) : nothing
    single_phase = card.machine_type in (6, 7) &&
        _deck_universal_machine_initialization_mode(parsed) == :automatic ?
        _deck_single_phase_induction_initialization(parsed, machine_index) : nothing
    steady_state = automatic_direct !== nothing ? automatic_direct.steady_state :
        single_phase !== nothing ? single_phase.steady_state :
        deck_steady_state_voltage_phasors(parsed)
    seed_state = automatic_direct !== nothing ? automatic_direct.state :
        single_phase !== nothing ? single_phase.state :
        card.machine_type in 3:7 ?
        deck_induction_machine_initial_state(
            parsed;
            machine_index,
            steady_state,
        ) :
        deck_direct_current_machine_initial_state(
            parsed;
            machine_index,
            steady_state,
        )
    point = _emt_machine_frequency_point(steady_state, request)
    horizon = run_deck_universal_machine_horizon(
        parsed;
        machine_index,
        dynamic_step_count=1,
    )
    time_zero_state = InductionMachineState(
        copy(horizon.current_values[:, 1]),
        copy(horizon.history_currents[:, 1]);
        mechanical_speed_rad_s=horizon.mechanical_speed_rad_s[1],
        previous_mechanical_speed_rad_s=horizon.mechanical_speed_rad_s[1],
        mechanical_angle_rad=horizon.mechanical_angle_rad[1],
    )
    time_zero_state.d_axis_flux = horizon.d_axis_flux[1]
    time_zero_state.q_axis_flux = horizon.q_axis_flux[1]
    time_zero_state.generated_torque = horizon.generated_torque[1]
    time_zero_state.output_values .= horizon.output_values[:, 1]
    time_zero_state.call_count = 1
    tolerances = request.tolerances
    current_seed_error = maximum(
        abs,
        seed_state.current_values .- horizon.current_values[:, 1];
        init=0.0,
    )
    current_scale = maximum(abs, horizon.current_values[:, 1]; init=0.0)
    residuals = EMTInitializationResidual[
        _emt_machine_network_residual(point, request),
        _emt_model_initialization_residual(
            card.machine_type in 3:7 ? :induction_machine_state :
                :direct_current_machine_state,
            :time_zero_current_equilibrium,
            "A",
            current_seed_error,
            tolerances.current_absolute_a,
            tolerances.current_relative,
            current_scale,
        ),
        _emt_model_initialization_residual(
            card.machine_type in 3:7 ? :induction_machine_state :
                :direct_current_machine_state,
            :coupled_runtime_completion,
            "count",
            horizon.complete_induction_machine_path ? 0.0 : 1.0,
            eps(Float64),
            0.0,
            1.0,
        ),
    ]
    if automatic_direct !== nothing &&
       hasproperty(automatic_direct, :armature_kvl_residual)
        push!(
            residuals,
            _emt_model_initialization_residual(
                :direct_current_machine_state,
                :armature_kvl,
                "V",
                abs(Float64(automatic_direct.armature_kvl_residual)),
                tolerances.voltage_absolute_v,
                tolerances.voltage_relative,
                abs(Float64(automatic_direct.requested_armature_voltage)),
            ),
        )
        push!(
            residuals,
            _emt_model_initialization_residual(
                :direct_current_machine_state,
                :field_kvl,
                "V",
                abs(Float64(automatic_direct.field_kvl_residual)),
                tolerances.voltage_absolute_v,
                tolerances.voltage_relative,
                abs(Float64(automatic_direct.field_voltage)),
            ),
        )
    end
    current_envelope_drift = abs(
        norm(horizon.current_values[:, 2]) -
            norm(horizon.current_values[:, 1]),
    ) / max(norm(horizon.current_values[:, 1]), tolerances.current_absolute_a)
    initial_flux_envelope = hypot(horizon.d_axis_flux[1], horizon.q_axis_flux[1])
    flux_envelope_drift = abs(
        hypot(horizon.d_axis_flux[2], horizon.q_axis_flux[2]) -
            initial_flux_envelope,
    ) / max(initial_flux_envelope, tolerances.flux_absolute_wb)
    torque_envelope_drift = abs(
        horizon.generated_torque[2] - horizon.generated_torque[1],
    ) / max(abs(horizon.generated_torque[1]), tolerances.power_absolute_w)
    envelope_drift = maximum((
        current_envelope_drift,
        flux_envelope_drift,
        torque_envelope_drift,
    ))
    frequency_step = pi * request.frequency_hz *
        DeckParser.deck_fixed_time_horizon_options(parsed).dt_s
    physical_to_trapezoidal_warping = frequency_step == 0.0 ? 0.0 :
        abs(tan(frequency_step) / frequency_step - 1.0)
    excess_envelope_drift = max(
        0.0,
        envelope_drift - physical_to_trapezoidal_warping,
    )
    transient_threshold = request.tolerances.no_artificial_transient_normalized_rms
    metric = NoArtificialTransientMetric(
        :machine_periodic_envelope,
        "1",
        request.time_origin_s,
        request.time_origin_s +
            DeckParser.deck_fixed_time_horizon_options(parsed).dt_s,
        excess_envelope_drift,
        excess_envelope_drift,
        envelope_drift,
        0.0,
        transient_threshold,
        isfinite(excess_envelope_drift) &&
            excess_envelope_drift <= transient_threshold,
    )
    accepted_state = (
        machine_state=time_zero_state,
        network_voltage_phasors=copy(steady_state.node_voltage_phasors),
        compensated_voltage_values=copy(horizon.compensated_voltage_values[:, 1]),
        power_terminal_voltage_values=copy(horizon.power_terminal_voltages[:, 1]),
        current_substitution_values=copy(horizon.current_substitution_values[:, 1]),
        drive_source_value=horizon.drive_source_values[1],
        excitation_source_value=horizon.excitation_source_values[1],
        report_output_names=copy(horizon.report_output_names),
        report_output_values=copy(horizon.report_output_values[:, 1]),
        switch_node_groups=copy.(point.topology.switch_node_groups),
        output_cursor=0,
    )
    machine_family = card.machine_type in 3:7 ? :induction_machine :
        :direct_current_machine
    prepared = PreparedMachineEMTStudy(
        machine_family,
        accepted_state,
        horizon,
        parsed,
    )
    inventory = _emt_machine_state_inventory(
        parsed,
        machine_family,
        length(time_zero_state.current_values) +
            length(time_zero_state.history_currents) + 2,
        3,
    )
    return prepared, point, residuals, inventory, metric, String[]
end

function _initialize_model_owned_machine_emt_study(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
)
    points = EMTInitializationFrequencyPoint[]
    residuals = EMTInitializationResidual[]
    state_inventory = EMTInitializationStateRecord[]
    transient_metrics = NoArtificialTransientMetric[]
    try
        working_parsed = deepcopy(parsed)
        _emt_validate_machine_initialization_request(working_parsed, request)
        detailed_synchronous = !isempty(
            DeckParser.deck_synchronous_machine_terminal_voltage_rows(working_parsed),
        )
        universal = !isempty(
            DeckParser.deck_universal_machine_definition_rows(working_parsed),
        )
        detailed_synchronous != universal || throw(_EMTInitializationRefusal(
            :ambiguous_machine_owner,
            :machine_equilibrium,
            :machine_family,
            "one cohesive initialization request must resolve to exactly one machine owner",
            (detailed_synchronous, universal),
        ))
        prepared, point, residuals, state_inventory, metric, warnings =
            detailed_synchronous ?
            _emt_synchronous_machine_preparation(working_parsed, request) :
            _emt_universal_machine_preparation(working_parsed, request)
        push!(points, point)
        push!(transient_metrics, metric)
        all(residual -> residual.passed, residuals) || throw(
            _EMTInitializationRefusal(
                :residual_failure,
                :machine_equilibrium,
                :scaled_residual,
                "one or more machine initialization residuals exceeded their quantity-specific limits",
                (failed_count=count(residual -> !residual.passed, residuals),),
            ),
        )
        metric.passed || throw(_EMTInitializationRefusal(
            :excessive_artificial_transient,
            :machine_equilibrium,
            metric.quantity,
            "machine initialization exceeded its physical-frequency first-step envelope allowance",
            (
                normalized_rms=metric.normalized_rms,
                envelope_drift=metric.low_frequency_envelope_drift,
                threshold=metric.threshold,
            ),
        ))
        initialized_state_owners = _emt_initialized_state_owners(state_inventory)
        mappings = OperatingPointMappingRecord[]
        signature = _emt_initialization_state_signature(
            prepared,
            request,
            mappings,
        )
        report = EMTInitializationReport(
            :accepted,
            :physical_frequency,
            request.frequency_hz,
            request.time_origin_s,
            point.topology,
            points,
            residuals,
            mappings,
            state_inventory,
            initialized_state_owners,
            Symbol[],
            transient_metrics,
            warnings,
            request.project_signature,
            request.settings_signature,
            request.model_signature,
            signature,
        )
        return EMTInitializationResult(prepared, report, nothing)
    catch error
        refusal = error isa _EMTInitializationRefusal ? error :
            _EMTInitializationRefusal(
                :initialization_error,
                :machine_equilibrium,
                :state,
                sprint(showerror, error),
                (exception_type=string(typeof(error)),),
            )
        return _emt_initialization_failure_result(
            request,
            refusal;
            points,
            residuals,
            state_inventory,
            transient_metrics,
        )
    end
end

function _average_value_converter_current_source_phasors(
    parsed::DeckParser.DeckParseResult,
    converter::AverageValueGridFollowingConverterInitialization,
    frequency_hz::Float64,
)
    phasors = ComplexF64[]
    for owner in converter.current_source_names
        element_index = findfirst(==(owner), parsed.element_names)
        element_index === nothing && throw(_EMTInitializationRefusal(
            :missing_model_input,
            :average_value_grid_following_converter,
            :phase_current_source,
            "average-value converter phase-current owner $(owner) is missing from the EMT network",
            (current_source=owner,),
        ))
        element = parsed.elements[element_index]
        element isa CurrentInjection &&
            element.value isa SinusoidalSourceSignal || throw(
                _EMTInitializationRefusal(
                    :unsupported_model_input,
                    :average_value_grid_following_converter,
                    :phase_current_source,
                    "average-value converter phase-current owners must be typed sinusoidal current injections",
                    (current_source=owner, element_type=string(typeof(element))),
                ),
            )
        push!(
            phasors,
            sinusoidal_source_peak_phasor(element.value, frequency_hz),
        )
    end
    return Tuple(phasors)
end

function _average_value_converter_terminal_phasors(
    parsed::DeckParser.DeckParseResult,
    point::EMTInitializationFrequencyPoint,
    converter::AverageValueGridFollowingConverterInitialization,
)
    phasors = ComplexF64[]
    for terminal in converter.terminal_nodes
        node = get(parsed.node_map, terminal, 0)
        node > 0 || throw(_EMTInitializationRefusal(
            :missing_model_input,
            :average_value_grid_following_converter,
            :phase_terminal,
            "average-value converter phase terminal $(terminal) is missing from the EMT network",
            (terminal,),
        ))
        push!(phasors, point.node_voltage_phasors[node])
    end
    return Tuple(phasors)
end

function _average_value_converter_state_values(state::InverterState)
    return (state.id_a, state.iq_a, state.xid, state.xiq)
end

function initialize_emt_study(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
    converter::AverageValueGridFollowingConverterInitialization;
    kwargs...,
)
    network_result = initialize_emt_study(parsed, request; kwargs...)
    initialization_accepted(network_result) || return network_result
    points = copy(network_result.report.frequency_scan)
    mappings = copy(network_result.report.mappings)
    network_residuals = copy(network_result.report.residuals)
    network_metrics = copy(network_result.report.transient_metrics)
    try
        parameters = converter.parameters
        isapprox(
            parameters.f_hz,
            request.frequency_hz;
            atol=1.0e-12,
            rtol=1.0e-12,
        ) || throw(_EMTInitializationRefusal(
            :frequency_mapping_mismatch,
            :average_value_grid_following_converter,
            :frequency_hz,
            "average-value converter frequency must equal the requested EMT operating frequency",
            (converter_frequency_hz=parameters.f_hz,
             request_frequency_hz=request.frequency_hz),
        ))
        primary = _emt_initialization_primary_point(points, request.frequency_hz)
        timestep_s = network_result.prepared.runtime_template.context.dt_s
        equilibrium = try
            grid_following_inverter_equilibrium(
                parameters;
                timestep_s,
            )
        catch error
            throw(_EMTInitializationRefusal(
                :infeasible_model_state,
                :average_value_grid_following_converter,
                :dq_equilibrium,
                sprint(showerror, error),
                (exception_type=string(typeof(error)),),
            ))
        end
        terminal_phasors = _average_value_converter_terminal_phasors(
            parsed,
            primary,
            converter,
        )
        expected_voltage_magnitude_v =
            sqrt(2.0) * parameters.v_ll_rms_v / sqrt(3.0)
        phase_rotations = (
            1.0 + 0.0im,
            cis(-2.0 * pi / 3.0),
            cis(2.0 * pi / 3.0),
        )
        expected_terminal_phasors = ntuple(
            phase -> expected_voltage_magnitude_v * phase_rotations[phase],
            3,
        )
        terminal_voltage_error_v = maximum(
            abs(terminal_phasors[phase] - expected_terminal_phasors[phase])
            for phase in 1:3
        )
        voltage_allowance_v = request.tolerances.voltage_absolute_v +
            request.tolerances.voltage_relative * expected_voltage_magnitude_v
        terminal_voltage_error_v <= voltage_allowance_v || throw(
            _EMTInitializationRefusal(
                :infeasible_model_state,
                :average_value_grid_following_converter,
                :terminal_voltage_phasors,
                "converter terminal phasors do not match the declared balanced grid-voltage basis",
                (error_v=terminal_voltage_error_v,
                 allowance_v=voltage_allowance_v),
            ),
        )
        source_phasors = _average_value_converter_current_source_phasors(
            parsed,
            converter,
            request.frequency_hz,
        )
        current_injection_error_a = maximum(
            abs(
                source_phasors[phase] -
                equilibrium.phase_current_phasors_a[phase],
            )
            for phase in 1:3
        )
        current_scale_a = max(
            maximum(abs, equilibrium.phase_current_phasors_a),
            request.tolerances.current_absolute_a,
        )
        current_allowance_a = request.tolerances.current_absolute_a +
            request.tolerances.current_relative * current_scale_a
        current_injection_error_a <= current_allowance_a || throw(
            _EMTInitializationRefusal(
                :infeasible_model_state,
                :average_value_grid_following_converter,
                :terminal_current_phasors,
                "network current-source phasors do not match the model-owned dq current equilibrium",
                (error_a=current_injection_error_a,
                 allowance_a=current_allowance_a),
            ),
        )
        derivative = equilibrium.derivative
        current_derivative_a_per_s = max(abs(derivative.id_a), abs(derivative.iq_a))
        controller_error_a = max(abs(derivative.xid), abs(derivative.xiq))
        initial_values = _average_value_converter_state_values(equilibrium.state)
        advanced_values = _average_value_converter_state_values(
            equilibrium.one_step_state,
        )
        current_recurrence_error_a = max(
            abs(advanced_values[1] - initial_values[1]),
            abs(advanced_values[2] - initial_values[2]),
        )
        integral_recurrence_error_as = max(
            abs(advanced_values[3] - initial_values[3]),
            abs(advanced_values[4] - initial_values[4]),
        )
        power_row = inverter_row(equilibrium.state, 0.0, parameters)
        active_power_error_w = abs(
            power_row[6] * parameters.s_base_va - equilibrium.active_power_w,
        )
        reactive_power_error_var = abs(
            power_row[7] * parameters.s_base_va - equilibrium.reactive_power_var,
        )
        tolerances = request.tolerances
        converter_residuals = EMTInitializationResidual[
            _emt_model_initialization_residual(
                :average_value_grid_following_converter,
                :balanced_terminal_voltage_mapping,
                "V",
                terminal_voltage_error_v,
                tolerances.voltage_absolute_v,
                tolerances.voltage_relative,
                expected_voltage_magnitude_v,
            ),
            _emt_model_initialization_residual(
                :average_value_grid_following_converter,
                :terminal_current_injection,
                "A",
                current_injection_error_a,
                tolerances.current_absolute_a,
                tolerances.current_relative,
                current_scale_a,
            ),
            _emt_model_initialization_residual(
                :average_value_grid_following_converter,
                :dq_current_derivative,
                "A/s",
                current_derivative_a_per_s,
                tolerances.current_absolute_a / timestep_s,
                tolerances.current_relative,
                current_scale_a / timestep_s,
            ),
            _emt_model_initialization_residual(
                :average_value_grid_following_converter,
                :controller_current_error,
                "A",
                controller_error_a,
                tolerances.current_absolute_a,
                tolerances.current_relative,
                current_scale_a,
            ),
            _emt_model_initialization_residual(
                :average_value_grid_following_converter,
                :one_step_current_recurrence,
                "A",
                current_recurrence_error_a,
                tolerances.current_absolute_a,
                tolerances.current_relative,
                current_scale_a,
            ),
            _emt_model_initialization_residual(
                :average_value_grid_following_converter,
                :one_step_controller_recurrence,
                "A*s",
                integral_recurrence_error_as,
                tolerances.current_absolute_a * timestep_s,
                tolerances.current_relative,
                current_scale_a * timestep_s,
            ),
            _emt_model_initialization_residual(
                :average_value_grid_following_converter,
                :active_power_equilibrium,
                "W",
                active_power_error_w,
                tolerances.power_absolute_w,
                tolerances.power_relative,
                abs(equilibrium.active_power_w),
            ),
            _emt_model_initialization_residual(
                :average_value_grid_following_converter,
                :reactive_power_equilibrium,
                "var",
                reactive_power_error_var,
                tolerances.power_absolute_w,
                tolerances.power_relative,
                abs(equilibrium.reactive_power_var),
            ),
        ]
        all(residual -> residual.passed, converter_residuals) || throw(
            _EMTInitializationRefusal(
                :residual_failure,
                :average_value_grid_following_converter,
                :scaled_residual,
                "average-value converter equilibrium exceeded one or more quantity-specific limits",
                (failed_count=count(residual -> !residual.passed, converter_residuals),),
            ),
        )
        state_inventory = copy(network_result.report.state_inventory)
        _append_emt_initialization_state!(
            state_inventory,
            :converter_filter_current_state,
            :continuous,
            2,
            :balanced_dq_zero_derivative_equilibrium,
        )
        _append_emt_initialization_state!(
            state_inventory,
            :converter_current_controller_state,
            :continuous,
            2,
            :zero_error_integral_bias,
        )
        _append_emt_initialization_state!(
            state_inventory,
            :converter_terminal_current_injection,
            :algebraic,
            3,
            :balanced_positive_sequence_current_phasors,
        )
        _append_emt_initialization_state!(
            state_inventory,
            :converter_dc_voltage_input,
            :algebraic,
            1,
            :declared_constant_dc_boundary,
        )
        sort!(state_inventory; by=record -> (
            String(record.state_family),
            String(record.owner),
        ))
        initial_state_scale = max(
            maximum(abs, initial_values),
            tolerances.current_absolute_a,
        )
        state_recurrence_error = maximum(
            abs(advanced_values[index] - initial_values[index])
            for index in eachindex(initial_values)
        )
        normalized_state_recurrence = state_recurrence_error / initial_state_scale
        converter_metric = NoArtificialTransientMetric(
            :average_value_converter_state,
            "scaled state",
            request.time_origin_s,
            request.time_origin_s + timestep_s,
            normalized_state_recurrence,
            normalized_state_recurrence,
            normalized_state_recurrence,
            0.0,
            tolerances.no_artificial_transient_normalized_rms,
            normalized_state_recurrence <=
                tolerances.no_artificial_transient_normalized_rms,
        )
        converter_metric.passed || throw(_EMTInitializationRefusal(
            :excessive_artificial_transient,
            :average_value_grid_following_converter,
            :state_recurrence,
            "average-value converter state changed during its undisturbed first-step probe",
            (normalized_state_recurrence,
             threshold=tolerances.no_artificial_transient_normalized_rms),
        ))
        prepared = PreparedAverageValueGridFollowingEMTStudy(
            network_result.prepared,
            equilibrium,
            converter,
        )
        signature = _emt_initialization_state_signature(
            prepared,
            request,
            mappings,
        )
        report = EMTInitializationReport(
            :accepted,
            network_result.report.formulation,
            request.frequency_hz,
            request.time_origin_s,
            network_result.report.topology,
            points,
            [network_residuals; converter_residuals],
            mappings,
            state_inventory,
            _emt_initialized_state_owners(state_inventory),
            Symbol[],
            [network_metrics; converter_metric],
            vcat(
                network_result.report.warnings,
                [
                    "The admitted average-value converter has a constant DC-voltage boundary and no switching, PLL, grid-forming, or sampled-control state.",
                ],
            ),
            request.project_signature,
            request.settings_signature,
            request.model_signature,
            signature,
        )
        return EMTInitializationResult(prepared, report, nothing)
    catch error
        refusal = error isa _EMTInitializationRefusal ? error :
            _EMTInitializationRefusal(
                :initialization_error,
                :average_value_grid_following_converter,
                :state,
                sprint(showerror, error),
                (exception_type=string(typeof(error)),),
            )
        return _emt_initialization_failure_result(
            request,
            refusal;
            points,
            mappings,
            residuals=network_residuals,
            transient_metrics=network_metrics,
        )
    end
end

function _emt_initialization_failure_result(
    request::EMTInitializationRequest,
    refusal::_EMTInitializationRefusal;
    points::Vector{EMTInitializationFrequencyPoint}=
        EMTInitializationFrequencyPoint[],
    mappings::Vector{OperatingPointMappingRecord}=
        OperatingPointMappingRecord[],
    residuals::Vector{EMTInitializationResidual}=
        EMTInitializationResidual[],
    initialized_state_owners::Vector{Symbol}=Symbol[],
    state_inventory::Vector{EMTInitializationStateRecord}=
        EMTInitializationStateRecord[],
    unsupported_state_owners::Vector{Symbol}=Symbol[],
    transient_metrics::Vector{NoArtificialTransientMetric}=
        NoArtificialTransientMetric[],
)
    topology = isempty(points) ?
        _empty_emt_initialization_topology_report() :
        last(points).topology
    report = EMTInitializationReport(
        :failed,
        _emt_harmonic_formulation_symbol(request.formulation),
        request.frequency_hz,
        request.time_origin_s,
        topology,
        points,
        residuals,
        mappings,
        state_inventory,
        initialized_state_owners,
        unsupported_state_owners,
        transient_metrics,
        String[],
        request.project_signature,
        request.settings_signature,
        request.model_signature,
        "",
    )
    failure = EMTInitializationFailure(
        refusal.code,
        refusal.owner,
        refusal.quantity,
        refusal.message,
        refusal.context,
    )
    return EMTInitializationResult(nothing, report, failure)
end

function initialize_emt_study(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest;
    timestep_s::Union{Nothing,Real}=nothing,
    t_end_s::Real=timestep_s !== nothing ? Float64(timestep_s) :
        request.formulation isa TimestepMatchedFormulation ?
        request.formulation.timestep_s :
        DeckParser.deck_fixed_time_horizon_options(parsed).dt_s,
    recorded_step_indices=nothing,
    output_schedule::Symbol=:all_steps,
    source_signal_provider::AbstractSourceSignalProvider=
        IdentitySourceSignalProvider(),
)
    points = EMTInitializationFrequencyPoint[]
    mappings = OperatingPointMappingRecord[]
    residuals = EMTInitializationResidual[]
    initialized_state_owners = Symbol[]
    state_inventory = EMTInitializationStateRecord[]
    unsupported_state_owners = Symbol[]
    transient_metrics = NoArtificialTransientMetric[]
    fixed_source_load_flow = nothing
    try
        DeckParser.assert_deck_valid!(parsed)
        if !isempty(DeckParser.deck_synchronous_machine_terminal_voltage_rows(parsed)) ||
           !isempty(DeckParser.deck_universal_machine_definition_rows(parsed))
            return _initialize_model_owned_machine_emt_study(parsed, request)
        end
        working_parsed = deepcopy(parsed)
        saturated_transformer_declared = get(
            working_parsed.card_counts,
            :fixed_card_saturated_transformer_intake,
            0,
        ) > 0
        saturated_transformer_intake =
            _deck_runtime_saturated_transformer_intake(working_parsed)
        saturated_transformer_declared &&
            saturated_transformer_intake === nothing && throw(
                _EMTInitializationRefusal(
                    :missing_source_backed_intake,
                    :saturated_transformer_magnetic_state,
                    :source_path,
                    "saturated-transformer initialization requires its exact source-backed characteristic and winding intake",
                    (source=working_parsed.source,),
                ),
            )
        unsupported_state_owners =
            _emt_unsupported_initialization_owners(working_parsed)
        isempty(unsupported_state_owners) || throw(_EMTInitializationRefusal(
            :unsupported_state_owner,
            :model_state,
            first(unsupported_state_owners),
            "the requested model requires unsupported initialization state: " *
            join(String.(unsupported_state_owners), ", "),
            (unsupported=copy(unsupported_state_owners),),
        ))
        if !isempty(DeckParser.deck_fixed_source_constraint_rows(working_parsed))
            applied = try
                apply_deck_fixed_source_load_flow(
                    working_parsed;
                    relative_power_tolerance=
                        _fixed_source_normalized_power_tolerance(
                            working_parsed,
                            request.tolerances.power_absolute_w,
                            request.tolerances.power_relative,
                        ),
                    maximum_iterations=
                        request.tolerances.operating_point_maximum_iterations,
                )
            catch error
                throw(_EMTInitializationRefusal(
                    :nonconvergent_operating_point,
                    :fixed_source_operating_point,
                    :active_reactive_power,
                    sprint(showerror, error),
                    (exception_type=string(typeof(error)),),
                ))
            end
            working_parsed = applied.deck
            fixed_source_load_flow = applied.load_flow
        end
        _validate_emt_initialization_source_frequency(working_parsed, request)
        transformer_initial_sample = nothing
        if saturated_transformer_intake !== nothing
            request.formulation isa PhysicalFrequencyFormulation || throw(
                _EMTInitializationRefusal(
                    :unsupported_formulation,
                    :saturated_transformer_magnetic_state,
                    :harmonic_formulation,
                    "saturated-transformer residual-flux initialization currently owns a physical-frequency equilibrium",
                    (requested_formulation=
                        _emt_harmonic_formulation_symbol(request.formulation),),
                ),
            )
            request.operating_point === nothing || throw(
                _EMTInitializationRefusal(
                    :unsupported_operating_point_mapping,
                    :saturated_transformer_magnetic_state,
                    :operating_point,
                    "saturated-transformer operating-point import requires explicit winding current and residual-flux quantities",
                    (operating_point_type=string(typeof(request.operating_point)),),
                ),
            )
            request.time_origin_s == 0.0 || throw(
                _EMTInitializationRefusal(
                    :unsupported_time_origin,
                    :saturated_transformer_magnetic_state,
                    :time_origin_s,
                    "declared saturated-transformer residual flux is referenced to the deck time origin",
                    (time_origin_s=request.time_origin_s,),
                ),
            )
            all(frequency -> isapprox(
                frequency,
                request.frequency_hz;
                atol=1.0e-12,
                rtol=1.0e-12,
            ), request.frequency_grid_hz) || throw(
                _EMTInitializationRefusal(
                    :unsupported_frequency_scan,
                    :saturated_transformer_magnetic_state,
                    :frequency_grid_hz,
                    "one residual-flux state cannot be reused across a frequency grid",
                    (frequency_grid_hz=copy(request.frequency_grid_hz),),
                ),
            )
            transformer_initial_sample = deck_steady_state_voltage_phasors(
                working_parsed;
                saturated_transformer_intake,
            )
            points = EMTInitializationFrequencyPoint[
                _emt_machine_frequency_point(
                    transformer_initial_sample,
                    request;
                    frequency_assignment=:initial_operating_point,
                ),
            ]
        else
            points = _emt_initialization_scan(working_parsed, request)
        end
        primary = _emt_initialization_primary_point(points, request.frequency_hz)
        primary.passed || throw(_emt_initialization_classification_failure(
            primary.topology,
        ))
        mappings = _emt_operating_point_mappings(
            working_parsed,
            request,
            primary,
        )
        sample = transformer_initial_sample === nothing ?
            _emt_initial_voltage_sample(working_parsed, request, primary) :
            transformer_initial_sample
        resolved_timestep_s = timestep_s === nothing ?
            request.formulation isa TimestepMatchedFormulation ?
                request.formulation.timestep_s :
                DeckParser.deck_fixed_time_horizon_options(working_parsed).dt_s :
            Float64(timestep_s)
        isfinite(resolved_timestep_s) && resolved_timestep_s > 0.0 || throw(
            _EMTInitializationRefusal(
                :missing_probe_timestep,
                :study_horizon,
                :timestep_s,
                "initialized EMT state requires a finite positive probe timestep",
                (timestep_s=resolved_timestep_s,),
            ),
        )
        horizon = Float64(t_end_s)
        isfinite(horizon) && horizon >= resolved_timestep_s || throw(
            _EMTInitializationRefusal(
                :missing_probe_horizon,
                :study_horizon,
                :t_end_s,
                "initialized EMT horizon must include at least one timestep",
                (timestep_s=resolved_timestep_s, t_end_s=horizon),
            ),
        )
        prepared = try
            prepare_emt_study(
                working_parsed;
                dt_s=resolved_timestep_s,
                t_end_s=horizon,
                initial_voltage_sample=sample,
                saturated_transformer_branch_runtime_enabled=
                    saturated_transformer_intake !== nothing,
                coupled_lumped_sequence_history_enabled=true,
                recorded_step_indices,
                output_schedule,
                source_signal_provider,
            )
        catch error
            message = sprint(showerror, error)
            if saturated_transformer_intake !== nothing &&
               error isa ArgumentError && occursin("characteristic", message)
                throw(_EMTInitializationRefusal(
                    :invalid_characteristic_state,
                    :saturated_transformer_magnetic_state,
                    :declared_current_flux,
                    message,
                    (exception_type=string(typeof(error)),),
                ))
            end
            rethrow()
        end
        _shift_emt_initialization_time_origin!(prepared, request.time_origin_s)
        transient_metric = _emt_initialization_probe_metric(
            prepared,
            primary,
            request,
        )
        push!(transient_metrics, transient_metric)
        transient_metric.passed || throw(_EMTInitializationRefusal(
            :excessive_artificial_transient,
            :no_artificial_transient,
            transient_metric.quantity,
            "initialized state exceeded the no-artificial-transient threshold",
            (
                normalized_rms=transient_metric.normalized_rms,
                maximum_scaled_discontinuity=
                    transient_metric.maximum_scaled_discontinuity,
                threshold=transient_metric.threshold,
            ),
        ))
        residuals = _emt_initialization_residuals(
            primary,
            mappings,
            request,
            fixed_source_load_flow,
        )
        append!(
            residuals,
            _emt_saturated_transformer_initialization_residuals(
                prepared,
                primary,
                request,
            ),
        )
        append!(
            residuals,
            _emt_pseudo_nonlinear_initialization_residuals(
                prepared,
                primary,
                request,
            ),
        )
        append!(
            residuals,
            _emt_piecewise_nonlinear_initialization_residuals(
                prepared,
                primary,
                request,
            ),
        )
        append!(
            residuals,
            _emt_hysteretic_initialization_residuals(
                prepared,
                primary,
                request,
            ),
        )
        all(residual -> residual.passed, residuals) || throw(
            _EMTInitializationRefusal(
                :residual_failure,
                :initialization_residual,
                :scaled_residual,
                "one or more initialization residuals exceeded their quantity-specific limits",
                (failed_count=count(residual -> !residual.passed, residuals),),
            ),
        )
        state_inventory = _emt_initialization_state_inventory(prepared)
        if fixed_source_load_flow !== nothing
            _append_emt_initialization_state!(
                state_inventory,
                :fixed_source_operating_point,
                :algebraic,
                length(fixed_source_load_flow.constraint_kinds),
                :coupled_active_reactive_power_solution,
            )
            sort!(
                state_inventory;
                by=record -> (
                    String(record.state_family),
                    String(record.owner),
                ),
            )
        end
        initialized_state_owners =
            _emt_initialized_state_owners(state_inventory)
        signature = _emt_initialization_state_signature(
            prepared,
            request,
            mappings,
        )
        report = EMTInitializationReport(
            :accepted,
            _emt_harmonic_formulation_symbol(request.formulation),
            request.frequency_hz,
            request.time_origin_s,
            primary.topology,
            points,
            residuals,
            mappings,
            state_inventory,
            initialized_state_owners,
            Symbol[],
            transient_metrics,
            String[],
            request.project_signature,
            request.settings_signature,
            request.model_signature,
            signature,
        )
        return EMTInitializationResult(prepared, report, nothing)
    catch error
        refusal = error isa _EMTInitializationRefusal ? error :
            _EMTInitializationRefusal(
                :initialization_error,
                :initialization_orchestration,
                :state,
                sprint(showerror, error),
                (exception_type=string(typeof(error)),),
            )
        return _emt_initialization_failure_result(
            request,
            refusal;
            points,
            mappings,
            residuals,
            initialized_state_owners,
            state_inventory,
            unsupported_state_owners,
            transient_metrics,
        )
    end
end

function initialize_emt_study(
    lines,
    request::EMTInitializationRequest;
    source::AbstractString="deck",
    kwargs...,
)
    parsed = DeckParser.parse_deck_lines(lines; source)
    return initialize_emt_study(parsed, request; kwargs...)
end

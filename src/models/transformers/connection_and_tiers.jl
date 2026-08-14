export TransformerApparatusTier,
       LowFrequencyTerminalTier,
       BCTRANTerminalTier,
       HybridTransformerTier,
       MagneticEquivalentCircuitTier,
       WidebandBlackBoxTier,
       GreyBoxLadderTier,
       WhiteBoxWindingTier,
       TRANSFORMER_APPARATUS_TIERS,
       ReactorApplication,
       ShuntReactorApplication,
       SeriesReactorApplication,
       NeutralReactorApplication,
       SmoothingReactorApplication,
       ReactorMagneticConstruction,
       AirCoreReactorConstruction,
       IronCoreReactorConstruction,
       ReactorWindingConfiguration,
       SingleCoilReactorWinding,
       SplitWindingReactorWinding,
       MutuallyCoupledReactorWinding,
       ReactorControlMode,
       FixedReactorControl,
       BreakerSwitchedReactorControl,
       ReactorGapModel,
       NoReactorAirGap,
       UniformReactorAirGap,
       EffectiveAreaReactorAirGap,
       ReactorApparatusDefinition,
       TransformerConnectionTopology,
       TransformerSourceRecord,
       transformer_apparatus_contract,
       transformer_coil_voltages!,
       transformer_terminal_currents!,
       transformer_connection_rank,
       transformer_connection_signature,
       transformer_connection_phase_shift

"""Explicit fidelity selected for one transformer or reactor apparatus."""
@enum TransformerApparatusTier begin
    LowFrequencyTerminalTier
    BCTRANTerminalTier
    HybridTransformerTier
    MagneticEquivalentCircuitTier
    WidebandBlackBoxTier
    GreyBoxLadderTier
    WhiteBoxWindingTier
end

const TRANSFORMER_APPARATUS_TIERS = (
    LowFrequencyTerminalTier,
    BCTRANTerminalTier,
    HybridTransformerTier,
    MagneticEquivalentCircuitTier,
    WidebandBlackBoxTier,
    GreyBoxLadderTier,
    WhiteBoxWindingTier,
)

@enum ReactorApplication begin
    ShuntReactorApplication
    SeriesReactorApplication
    NeutralReactorApplication
    SmoothingReactorApplication
end

@enum ReactorMagneticConstruction begin
    AirCoreReactorConstruction
    IronCoreReactorConstruction
end

@enum ReactorWindingConfiguration begin
    SingleCoilReactorWinding
    SplitWindingReactorWinding
    MutuallyCoupledReactorWinding
end

@enum ReactorControlMode begin
    FixedReactorControl
    BreakerSwitchedReactorControl
end

@enum ReactorGapModel begin
    NoReactorAirGap
    UniformReactorAirGap
    EffectiveAreaReactorAirGap
end

"""Explicit physical family and topology declaration for a reactor apparatus."""
struct ReactorApparatusDefinition
    application::ReactorApplication
    magnetic_construction::ReactorMagneticConstruction
    winding_configuration::ReactorWindingConfiguration
    control_mode::ReactorControlMode
    gap_model::ReactorGapModel
    total_air_gap_length_m::Float64
    air_gap_effective_area_factor::Float64

    function ReactorApparatusDefinition(
        application::ReactorApplication,
        magnetic_construction::ReactorMagneticConstruction;
        winding_configuration::ReactorWindingConfiguration=
            SingleCoilReactorWinding,
        control_mode::ReactorControlMode=FixedReactorControl,
        gap_model::ReactorGapModel=NoReactorAirGap,
        total_air_gap_length_m::Real=0.0,
        air_gap_effective_area_factor::Real=1.0,
    )
        gap_length = Float64(total_air_gap_length_m)
        effective_area_factor = Float64(air_gap_effective_area_factor)
        isfinite(gap_length) && gap_length >= 0.0 || throw(ArgumentError(
            "reactor total air-gap length must be finite and nonnegative",
        ))
        isfinite(effective_area_factor) && effective_area_factor >= 1.0 ||
            throw(ArgumentError(
                "reactor air-gap effective-area factor must be finite and at least one",
            ))
        if magnetic_construction === AirCoreReactorConstruction
            gap_model === NoReactorAirGap && gap_length == 0.0 &&
                effective_area_factor == 1.0 || throw(ArgumentError(
                    "air-core reactors cannot declare a ferromagnetic-core air gap",
                ))
        elseif gap_model === NoReactorAirGap
            gap_length == 0.0 && effective_area_factor == 1.0 ||
                throw(ArgumentError(
                    "an ungapped iron-core reactor cannot declare gap geometry",
                ))
        else
            gap_length > 0.0 || throw(ArgumentError(
                "a gapped iron-core reactor requires positive total air-gap length",
            ))
            gap_model === UniformReactorAirGap && effective_area_factor != 1.0 &&
                throw(ArgumentError(
                    "uniform-gap reactor modeling requires unit effective-area factor",
                ))
            gap_model === EffectiveAreaReactorAirGap &&
                effective_area_factor <= 1.0 && throw(ArgumentError(
                    "effective-area gap fringing requires a factor greater than one",
                ))
        end
        return new(
            application,
            magnetic_construction,
            winding_configuration,
            control_mode,
            gap_model,
            gap_length,
            effective_area_factor,
        )
    end
end

const _TRANSFORMER_TIER_IDS = Dict(
    LowFrequencyTerminalTier => :low_frequency_terminal,
    BCTRANTerminalTier => :bctran_terminal,
    HybridTransformerTier => :hybrid_transformer,
    MagneticEquivalentCircuitTier => :magnetic_equivalent_circuit,
    WidebandBlackBoxTier => :wideband_black_box,
    GreyBoxLadderTier => :grey_box_ladder,
    WhiteBoxWindingTier => :white_box_winding,
)

const _TRANSFORMER_TIER_UNSUPPORTED = Dict(
    LowFrequencyTerminalTier => (
        :core_local_flux,
        :hysteresis,
        :remanence,
        :frequency_dependent_winding_loss,
        :internal_winding_voltage,
    ),
    BCTRANTerminalTier => (
        :core_local_flux,
        :hysteresis,
        :remanence,
        :frequency_dependent_winding_loss,
        :internal_winding_voltage,
    ),
    HybridTransformerTier => (
        :black_box_internal_state_as_physical_state,
        :geometry_derived_turn_voltage,
        :three_dimensional_field_solution,
    ),
    MagneticEquivalentCircuitTier => (
        :unrepresented_winding_internal_voltage,
        :three_dimensional_field_solution,
    ),
    WidebandBlackBoxTier => (
        :core_local_flux,
        :winding_internal_voltage,
        :physical_label_for_fitted_state,
    ),
    GreyBoxLadderTier => (
        :unrepresented_turn_voltage,
        :geometry_derived_field_solution,
    ),
    WhiteBoxWindingTier => (
        :three_dimensional_field_solution,
        :insulation_lifetime_prediction,
        :partial_discharge,
    ),
)

function _transformer_contract_fidelity(tier::TransformerApparatusTier)
    tier in (LowFrequencyTerminalTier, BCTRANTerminalTier) &&
        return LegacyDetailed
    tier === WhiteBoxWindingTier && return FieldCoupledDetailed
    return SwitchingStateEquivalent
end

"""Return the public scientific contract for exactly one selected apparatus tier."""
function transformer_apparatus_contract(tier::TransformerApparatusTier)
    tier_id = _TRANSFORMER_TIER_IDS[tier]
    internal_outputs = if tier in (HybridTransformerTier, MagneticEquivalentCircuitTier)
        (
            ContractQuantity(:core_branch_flux_wb; unit="Wb", orientation="declared_magnetic_branch"),
            ContractQuantity(:core_branch_mmf_at; unit="A*turn", orientation="declared_magnetic_branch"),
        )
    elseif tier === GreyBoxLadderTier
        (
            ContractQuantity(:ladder_node_voltage_v; unit="V", orientation="declared_ladder_node_to_reference"),
        )
    elseif tier === WhiteBoxWindingTier
        (
            ContractQuantity(:winding_section_voltage_v; unit="V", orientation="declared_section_start_to_end"),
            ContractQuantity(:winding_section_current_a; unit="A", orientation="positive_along_declared_winding_order"),
        )
    else
        ()
    end
    return ScientificModelContract(
        Symbol(tier_id, :_transformer_apparatus),
        :instantaneous_emt_transformer_and_reactor;
        owner="AIMORA.TransformerApparatus and AIMORA.EMTStudy",
        maturity=:implemented,
        fidelity=_transformer_contract_fidelity(tier),
        validity_domain=ModelValidityDomain(
            Symbol(tier_id, :_fixed_step_domain);
            description="Explicit transformer/reactor tier with ordered physical terminals, declared connection incidence, SI state, fixed-step execution, typed output, and no automatic fallback.",
            bounds=(
                NumericDomainBound(:phase_count; unit="count", lower=1, upper=3),
                NumericDomainBound(:external_terminal_count; unit="count", lower=2, upper=30),
                NumericDomainBound(:rated_power_va; unit="VA", lower=100.0, upper=1.5e9),
                NumericDomainBound(:rated_voltage_v; unit="V", lower=10.0, upper=765.0e3),
                NumericDomainBound(:rated_frequency_hz; unit="Hz", lower=16.7, upper=400.0),
                NumericDomainBound(:timestep_s; unit="s", lower=2.0e-9, upper=100.0e-6),
            ),
            unsupported_phenomena=(
                _TRANSFORMER_TIER_UNSUPPORTED[tier]...,
                :manufacturer_parameter_prediction,
                :protected_standard_conformance,
                :atp_or_pscad_equivalence,
                :certification,
            ),
        ),
        state_inventory=DynamicStateInventory(
            differential=(
                :winding_current,
                :capacitor_charge,
                :passive_response_state,
                :magnetic_flux,
                :loss_energy,
            ),
            algebraic=(
                :terminal_voltage,
                :terminal_current,
                :connection_constraint,
                :magnetic_continuity,
                :instantaneous_power,
            ),
            discrete=(
                :tier_identity,
                :topology_identity,
                :hysteresis_direction,
                :event_mode,
                :factor_signature,
            ),
            delayed_history=(
                :companion_history,
                :rational_history,
                :section_history,
                :energy_history,
            ),
            scheduler=(
                :apparatus_event_boundary,
                :tap_or_switch_task,
                :output_cursor,
            ),
        ),
        inputs=(
            ContractQuantity(:terminal_voltage_v; unit="V", orientation="node_to_declared_reference"),
            ContractQuantity(:timestep_s; unit="s"),
            ContractQuantity(:initialization_frequency_hz; unit="Hz"),
        ),
        outputs=(
            ContractQuantity(:terminal_current_a; unit="A", orientation="positive_into_apparatus"),
            ContractQuantity(:terminal_power_w; unit="W", orientation="positive_into_apparatus"),
            ContractQuantity(:supplied_energy_j; unit="J", orientation="positive_into_apparatus"),
            internal_outputs...,
        ),
        assumptions=(
            "The tier is selected explicitly and unavailable higher-tier data never triggers inferred typical values or lower-tier fallback.",
            "Terminal current is positive into the apparatus and the same signed incidence is used for coil voltage and current injection.",
            "Only outputs represented by the selected tier are available; absent internal quantities remain typed unavailable.",
            "Fixed-step execution is the default and every accepted history advances exactly once.",
        ),
        mutation_order=(
            :validate_identity_topology_and_domain,
            :capture_complete_transaction,
            :apply_due_topology_and_event_commands,
            :assemble_pure_trial_apparatus_state,
            :solve_coupled_network,
            :localize_and_apply_simultaneous_events,
            :verify_kcl_flux_charge_energy_and_passivity,
            :accept_state_history_energy_and_output_once_or_restore,
        ),
    )
end

function _connection_symbols(values, label::AbstractString)
    result = Symbol.(values)
    isempty(result) && throw(ArgumentError("transformer $label must not be empty"))
    length(unique(result)) == length(result) ||
        throw(ArgumentError("transformer $label must be unique and ordered"))
    any(==(Symbol("")), result) &&
        throw(ArgumentError("transformer $label must not contain an empty identity"))
    return result
end

"""Explicit sparse winding connection with `u_coil = A' * e_node` and `i_node = A * i_coil`.

`node_order` omits the reference ground. A coil connected to ground therefore has one
nonzero incidence entry; a floating coil has one `+1` and one `-1`. `coil_winding`
and `coil_phase` bind every column to physical winding and phase identities. A
three-phase clock is verified against explicit node-phase and incidence polarity;
the label never synthesizes incidence or hidden nodes.
"""
struct TransformerConnectionTopology
    node_order::Vector{Symbol}
    node_phase::Vector{Union{Nothing,Symbol}}
    coil_order::Vector{Symbol}
    winding_order::Vector{Symbol}
    phase_order::Vector{Symbol}
    coil_winding::Vector{Symbol}
    coil_phase::Vector{Symbol}
    incidence::Matrix{Float64}
    vector_group::String
    clock_number::Union{Nothing,Int}
    phase_shift_rad::Float64
    grounded_coils::BitVector
    rank::Int
    deterministic_signature_sha256::String
end

function _connection_signature(
    node_order,
    node_phase,
    coil_order,
    winding_order,
    phase_order,
    coil_winding,
    coil_phase,
    incidence,
    vector_group,
    clock_number,
    phase_shift_rad,
    grounded_coils,
)
    io = IOBuffer()
    for values in (
        node_order,
        coil_order,
        winding_order,
        phase_order,
        coil_winding,
        coil_phase,
    )
        println(io, join(String.(values), ','))
    end
    println(io, join(
        (value === nothing ? "neutral" : String(value) for value in node_phase),
        ',',
    ))
    for value in incidence
        println(io, bitstring(value))
    end
    println(io, vector_group)
    println(io, clock_number === nothing ? "no_clock" : clock_number)
    println(io, bitstring(phase_shift_rad))
    println(io, join(Int.(grounded_coils), ','))
    return bytes2hex(sha256(take!(io)))
end

function TransformerConnectionTopology(;
    node_order,
    node_phase=nothing,
    coil_order,
    winding_order,
    phase_order,
    coil_winding,
    coil_phase,
    incidence,
    vector_group::AbstractString,
    clock_number::Union{Nothing,Integer}=nothing,
    phase_shift_rad::Real=0.0,
)
    nodes = _connection_symbols(node_order, "node order")
    coils = _connection_symbols(coil_order, "coil order")
    windings = _connection_symbols(winding_order, "winding order")
    phases = _connection_symbols(phase_order, "phase order")
    node_phase_owner = if node_phase === nothing
        Union{Nothing,Symbol}[nothing for _ in nodes]
    else
        Union{Nothing,Symbol}[
            value === nothing ? nothing : Symbol(value) for value in node_phase
        ]
    end
    length(node_phase_owner) == length(nodes) || throw(DimensionMismatch(
        "transformer node phase ownership must cover every declared node",
    ))
    all(value -> value === nothing || value in phases, node_phase_owner) ||
        throw(ArgumentError("transformer node references an unknown phase"))
    winding_owner = Symbol.(coil_winding)
    phase_owner = Symbol.(coil_phase)
    length(winding_owner) == length(coils) || throw(DimensionMismatch(
        "transformer coil winding ownership must cover every coil",
    ))
    length(phase_owner) == length(coils) || throw(DimensionMismatch(
        "transformer coil phase ownership must cover every coil",
    ))
    all(in(windings), winding_owner) ||
        throw(ArgumentError("transformer coil references an unknown winding"))
    all(in(phases), phase_owner) ||
        throw(ArgumentError("transformer coil references an unknown phase"))
    length(unique(zip(winding_owner, phase_owner))) == length(coils) ||
        throw(ArgumentError("transformer winding-phase ownership must be unique"))
    matrix = Matrix{Float64}(incidence)
    size(matrix) == (length(nodes), length(coils)) || throw(DimensionMismatch(
        "transformer incidence size must be node_count by coil_count",
    ))
    all(isfinite, matrix) ||
        throw(ArgumentError("transformer incidence entries must be finite"))
    incidence_tolerance = 64.0 * eps(Float64)
    for value in matrix
        any(isapprox(value, allowed; atol=incidence_tolerance, rtol=0.0) for allowed in (-1.0, 0.0, 1.0)) ||
            throw(ArgumentError("transformer incidence entries must be -1, 0, or 1"))
    end
    normalized = map(matrix) do value
        abs(value) <= incidence_tolerance ? 0.0 : sign(value)
    end
    grounded = falses(length(coils))
    for coil in axes(normalized, 2)
        nonzero = findall(!iszero, @view normalized[:, coil])
        length(nonzero) in (1, 2) || throw(ArgumentError(
            "transformer coil incidence must connect one node to ground or two distinct nodes",
        ))
        if length(nonzero) == 1
            grounded[coil] = true
        else
            sort(normalized[nonzero, coil]) == [-1.0, 1.0] || throw(ArgumentError(
                "floating transformer coil incidence must contain one +1 and one -1",
            ))
        end
    end
    all(row -> any(!iszero, @view normalized[row, :]), axes(normalized, 1)) ||
        throw(ArgumentError("transformer connection contains an isolated declared node"))
    group = strip(String(vector_group))
    isempty(group) && throw(ArgumentError("transformer vector group must not be empty"))
    clock = clock_number === nothing ? nothing : Int(clock_number)
    clock === nothing || 0 <= clock <= 11 ||
        throw(ArgumentError("transformer clock number must lie from 0 through 11"))
    shift = Float64(phase_shift_rad)
    isfinite(shift) || throw(ArgumentError("transformer phase shift must be finite"))
    if clock !== nothing
        clock_shift = clock * pi / 6.0
        wrapped_difference = mod(shift - clock_shift + pi, 2.0 * pi) - pi
        abs(wrapped_difference) <= 1.0e-10 || throw(ArgumentError(
            "transformer clock number disagrees with the declared phase shift",
        ))
    end
    if length(phases) == 1
        clock in (nothing, 0) || throw(ArgumentError(
            "single-phase transformer connections admit only clock zero",
        ))
    elseif length(phases) == 3
        any(isnothing, node_phase_owner) && all(isnothing, node_phase_owner) &&
            clock !== nothing && throw(ArgumentError(
                "three-phase clock verification requires explicit node phase ownership",
            ))
        if clock !== nothing
            length(windings) >= 2 || throw(ArgumentError(
                "three-phase clock verification requires at least two windings",
            ))
            positive_sequence_factor = cis(2.0 * pi / 3.0)
            phase_excitation = ComplexF64[
                value === nothing ? 0.0 + 0.0im :
                positive_sequence_factor^(-findfirst(==(value), phases) + 1)
                for value in node_phase_owner
            ]
            coil_excitation = transpose(normalized) * phase_excitation
            winding_sequence_voltage = ComplexF64[]
            for winding in windings[1:2]
                phase_voltage = ComplexF64[]
                for phase in phases
                    coil_index = findfirst(
                        index -> winding_owner[index] == winding &&
                            phase_owner[index] == phase,
                        eachindex(coils),
                    )
                    coil_index === nothing && throw(ArgumentError(
                        "three-phase clock verification requires one coil for every winding phase",
                    ))
                    push!(phase_voltage, coil_excitation[coil_index])
                end
                sequence_voltage = sum(
                    positive_sequence_factor^(phase_index - 1) *
                    phase_voltage[phase_index]
                    for phase_index in eachindex(phases)
                ) / 3.0
                abs(sequence_voltage) > 256.0 * eps(Float64) ||
                    throw(ArgumentError(
                        "transformer winding incidence has no positive-sequence voltage",
                    ))
                push!(winding_sequence_voltage, sequence_voltage)
            end
            incidence_shift = angle(
                winding_sequence_voltage[2] / winding_sequence_voltage[1],
            )
            incidence_difference = mod(incidence_shift - shift + pi, 2.0 * pi) - pi
            abs(incidence_difference) <= 1.0e-10 || throw(ArgumentError(
                "transformer clock and phase shift disagree with the explicit incidence",
            ))
        end
    else
        throw(ArgumentError("transformer connection requires one or three ordered phases"))
    end
    matrix_rank = rank(normalized)
    matrix_rank > 0 || throw(ArgumentError("transformer connection incidence has zero rank"))
    signature = _connection_signature(
        nodes,
        node_phase_owner,
        coils,
        windings,
        phases,
        winding_owner,
        phase_owner,
        normalized,
        group,
        clock,
        shift,
        grounded,
    )
    return TransformerConnectionTopology(
        nodes,
        node_phase_owner,
        coils,
        windings,
        phases,
        winding_owner,
        phase_owner,
        normalized,
        group,
        clock,
        shift,
        grounded,
        matrix_rank,
        signature,
    )
end

transformer_connection_rank(topology::TransformerConnectionTopology) = topology.rank
transformer_connection_signature(topology::TransformerConnectionTopology) =
    topology.deterministic_signature_sha256
transformer_connection_phase_shift(topology::TransformerConnectionTopology) =
    topology.phase_shift_rad

function transformer_coil_voltages!(
    coil_voltage_v::AbstractVector{Float64},
    topology::TransformerConnectionTopology,
    node_voltage_v::AbstractVector{Float64},
)
    length(node_voltage_v) == length(topology.node_order) || throw(DimensionMismatch(
        "transformer node voltage count must match the connection topology",
    ))
    length(coil_voltage_v) == length(topology.coil_order) || throw(DimensionMismatch(
        "transformer coil voltage workspace must match the connection topology",
    ))
    mul!(coil_voltage_v, transpose(topology.incidence), node_voltage_v)
    return coil_voltage_v
end

function transformer_terminal_currents!(
    node_current_a::AbstractVector{Float64},
    topology::TransformerConnectionTopology,
    coil_current_a::AbstractVector{Float64},
)
    length(coil_current_a) == length(topology.coil_order) || throw(DimensionMismatch(
        "transformer coil current count must match the connection topology",
    ))
    length(node_current_a) == length(topology.node_order) || throw(DimensionMismatch(
        "transformer node current workspace must match the connection topology",
    ))
    mul!(node_current_a, topology.incidence, coil_current_a)
    return node_current_a
end

"""Content identity and full scientific provenance of one apparatus input source."""
struct TransformerSourceRecord
    id::Symbol
    content_sha256::String
    provenance::ParameterProvenance

    function TransformerSourceRecord(
        id::Symbol,
        content_sha256::AbstractString,
        provenance::ParameterProvenance,
    )
        id == Symbol("") && throw(ArgumentError(
            "transformer source identity must not be empty",
        ))
        signature = lowercase(String(content_sha256))
        occursin(r"^[0-9a-f]{64}$", signature) || throw(ArgumentError(
            "transformer source content identity must be lowercase SHA-256",
        ))
        provenance.nature === PhysicalModelParameter || throw(ArgumentError(
            "transformer source provenance must describe physical model data",
        ))
        return new(id, signature, provenance)
    end
end

export AbstractTransformerResultQuantity,
       TransformerResultAvailable,
       TransformerResultUnavailable,
       transformer_result_quantity_available,
       transformer_result_quantity_value,
       TransformerApparatusEnergyResult,
       TransformerApparatusResidualResult,
       TransformerApparatusResult,
       transformer_apparatus_result

abstract type AbstractTransformerResultQuantity end

"""A quantity physically represented by the explicitly selected apparatus tier."""
struct TransformerResultAvailable{T} <: AbstractTransformerResultQuantity
    value::T
end

"""A quantity deliberately absent from the explicitly selected apparatus tier."""
struct TransformerResultUnavailable <: AbstractTransformerResultQuantity
    quantity::Symbol
    tier::TransformerApparatusTier
    reason::Symbol
end

transformer_result_quantity_available(::TransformerResultAvailable) = true
transformer_result_quantity_available(::TransformerResultUnavailable) = false

transformer_result_quantity_value(quantity::TransformerResultAvailable) = quantity.value

function transformer_result_quantity_value(quantity::TransformerResultUnavailable)
    throw(ArgumentError(
        "transformer result quantity $(quantity.quantity) is unavailable for " *
        "tier $(_TRANSFORMER_TIER_IDS[quantity.tier]): $(quantity.reason)",
    ))
end

struct TransformerApparatusEnergyResult
    initial_stored_energy_j::Float64
    apparatus_supplied_energy_j::Float64
    total_network_supplied_energy_j::Float64
    stored_energy_j::Float64
    apparatus_dissipated_energy_j::Float64
    external_fault_energy_j::Float64
    internal_fault_energy_j::Float64
    total_event_dissipated_energy_j::Float64
    numerical_dissipation_energy_j::Float64
    physical_balance_residual_j::Float64
    unexplained_balance_residual_j::Float64
end

struct TransformerApparatusResidualResult
    maximum_terminal_kcl_residual_a::AbstractTransformerResultQuantity
    maximum_internal_kcl_residual_a::AbstractTransformerResultQuantity
    maximum_magnetic_continuity_residual_wb::AbstractTransformerResultQuantity
    maximum_magnetic_constitutive_residual_at::AbstractTransformerResultQuantity
    maximum_energy_balance_residual_j::Float64
end

"""Immutable, versioned result for one accepted transformer or reactor runtime state.

Every vector-valued output is frozen as a tuple. Quantities not represented by the
selected tier are `TransformerResultUnavailable`; they are never filled with zeros or
inferred from another tier.
"""
struct TransformerApparatusResult
    schema_version::Int
    apparatus::Symbol
    tier::TransformerApparatusTier
    reactor_definition::Union{Nothing,ReactorApparatusDefinition}
    terminal_order::Tuple
    coil_order::Tuple
    source_ids::Tuple
    source_content_sha256::Tuple
    initial_time_s::Float64
    accepted_time_s::Float64
    accepted_step_count::Int
    terminal_voltage_v::Tuple
    terminal_current_a::Tuple
    coil_voltage_v::AbstractTransformerResultQuantity
    coil_current_a::AbstractTransformerResultQuantity
    magnetic_branch_flux_wb::AbstractTransformerResultQuantity
    magnetic_branch_mmf_drop_at::AbstractTransformerResultQuantity
    passive_rational_state::AbstractTransformerResultQuantity
    represented_node_voltage_v::AbstractTransformerResultQuantity
    represented_branch_voltage_v::AbstractTransformerResultQuantity
    represented_branch_current_a::AbstractTransformerResultQuantity
    ladder_node_voltage_v::AbstractTransformerResultQuantity
    winding_section_voltage_v::AbstractTransformerResultQuantity
    winding_section_current_a::AbstractTransformerResultQuantity
    energy::TransformerApparatusEnergyResult
    residuals::TransformerApparatusResidualResult
    event_occurrences::Tuple
    uncertainty::String
    validity_domain::String
    unsupported_outputs::Tuple
    preparation_signature_sha256::String
    snapshot_signature_sha256::String
    deterministic_signature_sha256::String
end

_transformer_result_tuple(values) = Tuple(copy(values))

function _transformer_result_available(values)
    return TransformerResultAvailable(_transformer_result_tuple(values))
end

function _transformer_result_unavailable(
    quantity::Symbol,
    tier::TransformerApparatusTier,
    reason::Symbol,
)
    return TransformerResultUnavailable(quantity, tier, reason)
end

function _transformer_result_energy(runtime::TransformerApparatusRuntime)
    state = runtime.accepted_state
    events = runtime.event_state
    apparatus_supplied = state.supplied_energy_j
    total_network_supplied = apparatus_supplied + events.external_fault_energy_j
    stored = _transformer_total_stored_energy(state)
    apparatus_dissipated = _transformer_total_dissipated_energy(state)
    physical_residual = runtime.initial_stored_energy_j + total_network_supplied -
        stored - apparatus_dissipated - events.event_energy_j
    unexplained_residual = physical_residual - events.numerical_dissipation_energy_j
    return TransformerApparatusEnergyResult(
        runtime.initial_stored_energy_j,
        apparatus_supplied,
        total_network_supplied,
        stored,
        apparatus_dissipated,
        events.external_fault_energy_j,
        events.internal_fault_energy_j,
        events.event_energy_j,
        events.numerical_dissipation_energy_j,
        physical_residual,
        unexplained_residual,
    )
end

function _transformer_result_residuals(runtime::TransformerApparatusRuntime)
    tier = runtime.preparation.specification.tier
    state = runtime.accepted_state
    terminal_kcl = state isa TransformerTerminalMatrixRuntimeState ?
        TransformerResultAvailable(state.maximum_kcl_residual_a) :
        _transformer_result_unavailable(
            :maximum_terminal_kcl_residual_a,
            tier,
            :not_owned_by_selected_runtime_state,
        )
    internal_kcl = state isa TransformerNetworkRuntimeState ?
        TransformerResultAvailable(state.maximum_internal_kcl_residual_a) :
        _transformer_result_unavailable(
            :maximum_internal_kcl_residual_a,
            tier,
            :no_represented_internal_network,
        )
    magnetic_continuity = state isa TransformerMagneticRuntimeState ?
        TransformerResultAvailable(state.maximum_magnetic_continuity_residual_wb) :
        _transformer_result_unavailable(
            :maximum_magnetic_continuity_residual_wb,
            tier,
            :no_explicit_magnetic_graph,
        )
    magnetic_constitutive = state isa TransformerMagneticRuntimeState ?
        TransformerResultAvailable(state.maximum_magnetic_constitutive_residual_at) :
        _transformer_result_unavailable(
            :maximum_magnetic_constitutive_residual_at,
            tier,
            :no_explicit_magnetic_graph,
        )
    return TransformerApparatusResidualResult(
        terminal_kcl,
        internal_kcl,
        magnetic_continuity,
        magnetic_constitutive,
        runtime.event_state.maximum_energy_balance_residual_j,
    )
end

function _transformer_result_internal_quantities(runtime::TransformerApparatusRuntime)
    state = runtime.accepted_state
    tier = runtime.preparation.specification.tier
    unavailable(quantity, reason) =
        _transformer_result_unavailable(quantity, tier, reason)
    coil_voltage = state isa TransformerTerminalMatrixRuntimeState ||
            state isa TransformerMagneticRuntimeState ?
        _transformer_result_available(state.coil_voltage_v) :
        unavailable(:coil_voltage_v, :not_owned_by_selected_tier)
    coil_current = state isa TransformerTerminalMatrixRuntimeState ||
            state isa TransformerMagneticRuntimeState ?
        _transformer_result_available(state.coil_current_a) :
        unavailable(:coil_current_a, :not_owned_by_selected_tier)
    magnetic_flux = state isa TransformerMagneticRuntimeState ?
        _transformer_result_available(state.branch_flux_wb) :
        unavailable(:magnetic_branch_flux_wb, :no_explicit_magnetic_graph)
    magnetic_mmf = state isa TransformerMagneticRuntimeState ?
        _transformer_result_available(state.branch_mmf_drop_at) :
        unavailable(:magnetic_branch_mmf_drop_at, :no_explicit_magnetic_graph)
    rational_state = state isa TransformerWidebandRuntimeState ?
        _transformer_result_available(state.rational_state) :
        unavailable(:passive_rational_state, :no_passive_black_box_state)
    represented_nodes = state isa TransformerNetworkRuntimeState ?
        _transformer_result_available(state.represented_node_voltage_v) :
        unavailable(:represented_node_voltage_v, :no_represented_internal_network)
    represented_branch_voltage = state isa TransformerNetworkRuntimeState ?
        _transformer_result_available(state.branch_voltage_v) :
        unavailable(:represented_branch_voltage_v, :no_represented_internal_network)
    represented_branch_current = state isa TransformerNetworkRuntimeState ?
        _transformer_result_available(state.branch_current_a) :
        unavailable(:represented_branch_current_a, :no_represented_internal_network)
    ladder_nodes = tier === GreyBoxLadderTier ? represented_nodes :
        unavailable(:ladder_node_voltage_v, :selected_tier_is_not_grey_box_ladder)
    section_voltage = tier === WhiteBoxWindingTier ? represented_branch_voltage :
        unavailable(:winding_section_voltage_v, :selected_tier_is_not_white_box_winding)
    section_current = tier === WhiteBoxWindingTier ? represented_branch_current :
        unavailable(:winding_section_current_a, :selected_tier_is_not_white_box_winding)
    return (
        coil_voltage=coil_voltage,
        coil_current=coil_current,
        magnetic_flux=magnetic_flux,
        magnetic_mmf=magnetic_mmf,
        rational_state=rational_state,
        represented_nodes=represented_nodes,
        represented_branch_voltage=represented_branch_voltage,
        represented_branch_current=represented_branch_current,
        ladder_nodes=ladder_nodes,
        section_voltage=section_voltage,
        section_current=section_current,
    )
end

"""Freeze the complete accepted apparatus state into its public typed result contract."""
function transformer_apparatus_result(runtime::TransformerApparatusRuntime)
    runtime.prepared && throw(ArgumentError(
        "transformer result cannot be created during an active trial step",
    ))
    preparation = runtime.preparation
    specification = preparation.specification
    state = runtime.accepted_state
    snapshot = transformer_apparatus_runtime_snapshot(runtime)
    internal = _transformer_result_internal_quantities(runtime)
    source_ids = Tuple(source.id for source in specification.sources)
    source_hashes = Tuple(source.content_sha256 for source in specification.sources)
    unsupported = Tuple(_TRANSFORMER_TIER_UNSUPPORTED[specification.tier])
    signature_io = IOBuffer()
    println(signature_io, "aimora.transformer_apparatus_result.v1")
    println(signature_io, specification.deterministic_signature_sha256)
    println(signature_io, snapshot.deterministic_signature_sha256)
    println(signature_io, join(String.(unsupported), ','))
    signature = bytes2hex(sha256(take!(signature_io)))
    return TransformerApparatusResult(
        1,
        specification.id,
        specification.tier,
        specification.reactor_definition,
        Tuple(preparation.terminal_order),
        Tuple(preparation.coil_order),
        source_ids,
        source_hashes,
        preparation.initial_time_s,
        state.accepted_time_s,
        state.accepted_step_count,
        _transformer_result_tuple(state.terminal_voltage_v),
        _transformer_result_tuple(state.terminal_current_a),
        internal.coil_voltage,
        internal.coil_current,
        internal.magnetic_flux,
        internal.magnetic_mmf,
        internal.rational_state,
        internal.represented_nodes,
        internal.represented_branch_voltage,
        internal.represented_branch_current,
        internal.ladder_nodes,
        internal.section_voltage,
        internal.section_current,
        _transformer_result_energy(runtime),
        _transformer_result_residuals(runtime),
        Tuple(transformer_apparatus_event_occurrences(runtime)),
        specification.uncertainty,
        specification.validity_domain,
        unsupported,
        preparation.preparation_signature_sha256,
        snapshot.deterministic_signature_sha256,
        signature,
    )
end

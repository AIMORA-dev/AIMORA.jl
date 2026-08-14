export TransformerApparatusEventKind,
       TransformerApparatusEventCommand,
       TransformerApparatusEventOccurrence,
       TransformerApparatusEventState,
       TransformerEnergizeEvent,
       TransformerDeenergizeEvent,
       TransformerBreakerOpenEvent,
       TransformerBreakerCloseEvent,
       TransformerTerminalFaultApplyEvent,
       TransformerTerminalFaultClearEvent,
       TransformerWindingFaultApplyEvent,
       TransformerWindingFaultClearEvent,
       TransformerGroundingApplyEvent,
       TransformerGroundingClearEvent,
       TransformerTapChangeEvent,
       TransformerPhaseShiftChangeEvent,
       TransformerInternalFaultApplyEvent,
       TransformerInternalFaultClearEvent,
       queue_transformer_apparatus_event!,
       transformer_apparatus_event_state,
       transformer_apparatus_event_occurrences

@enum TransformerApparatusEventKind begin
    TransformerEnergizeEvent
    TransformerDeenergizeEvent
    TransformerBreakerOpenEvent
    TransformerBreakerCloseEvent
    TransformerTerminalFaultApplyEvent
    TransformerTerminalFaultClearEvent
    TransformerWindingFaultApplyEvent
    TransformerWindingFaultClearEvent
    TransformerGroundingApplyEvent
    TransformerGroundingClearEvent
    TransformerTapChangeEvent
    TransformerPhaseShiftChangeEvent
    TransformerInternalFaultApplyEvent
    TransformerInternalFaultClearEvent
end

const _TRANSFORMER_EVENT_KIND_IDS = Dict(
    TransformerEnergizeEvent => :energize,
    TransformerDeenergizeEvent => :deenergize,
    TransformerBreakerOpenEvent => :breaker_open,
    TransformerBreakerCloseEvent => :breaker_close,
    TransformerTerminalFaultApplyEvent => :terminal_fault_apply,
    TransformerTerminalFaultClearEvent => :terminal_fault_clear,
    TransformerWindingFaultApplyEvent => :winding_fault_apply,
    TransformerWindingFaultClearEvent => :winding_fault_clear,
    TransformerGroundingApplyEvent => :grounding_apply,
    TransformerGroundingClearEvent => :grounding_clear,
    TransformerTapChangeEvent => :tap_change,
    TransformerPhaseShiftChangeEvent => :phase_shift_change,
    TransformerInternalFaultApplyEvent => :internal_fault_apply,
    TransformerInternalFaultClearEvent => :internal_fault_clear,
)

function _transformer_event_priority(kind::TransformerApparatusEventKind)
    kind in (
        TransformerDeenergizeEvent,
        TransformerBreakerOpenEvent,
    ) && return -30
    kind in (
        TransformerTerminalFaultApplyEvent,
        TransformerWindingFaultApplyEvent,
        TransformerInternalFaultApplyEvent,
        TransformerGroundingApplyEvent,
    ) && return -20
    kind in (
        TransformerTerminalFaultClearEvent,
        TransformerWindingFaultClearEvent,
        TransformerInternalFaultClearEvent,
        TransformerGroundingClearEvent,
    ) && return -10
    kind in (TransformerTapChangeEvent, TransformerPhaseShiftChangeEvent) && return 0
    return 10
end

function _transformer_event_topology_invalidating(kind::TransformerApparatusEventKind)
    return kind in (
        TransformerEnergizeEvent,
        TransformerDeenergizeEvent,
        TransformerBreakerOpenEvent,
        TransformerBreakerCloseEvent,
        TransformerTerminalFaultApplyEvent,
        TransformerTerminalFaultClearEvent,
        TransformerWindingFaultApplyEvent,
        TransformerWindingFaultClearEvent,
        TransformerGroundingApplyEvent,
        TransformerGroundingClearEvent,
        TransformerTapChangeEvent,
        TransformerPhaseShiftChangeEvent,
        TransformerInternalFaultApplyEvent,
        TransformerInternalFaultClearEvent,
    )
end

function _transformer_event_requires_apparatus_current_reconstruction(
    kind::TransformerApparatusEventKind,
)
    return kind in (
        TransformerEnergizeEvent,
        TransformerDeenergizeEvent,
        TransformerBreakerOpenEvent,
        TransformerBreakerCloseEvent,
        TransformerTapChangeEvent,
        TransformerPhaseShiftChangeEvent,
        TransformerInternalFaultApplyEvent,
        TransformerInternalFaultClearEvent,
    )
end

"""One exact, solver-neutral apparatus command applied only at a localized accepted boundary.

`target_indices` use the public terminal order except for represented internal faults,
where they use the declared grey-/white-box node order. Fault and grounding
conductance is in siemens. Tap and phase-shift commands replace the complete real
terminal voltage transform; currents use its transpose so ideal transformation
cannot create power. Clear commands identify the corresponding apply command with
`reference_id`.
"""
struct TransformerApparatusEventCommand
    id::Symbol
    kind::TransformerApparatusEventKind
    time_s::Float64
    target_indices::Vector{Int}
    conductance_s::Float64
    terminal_transform::Union{Nothing,Matrix{Float64}}
    reference_id::Union{Nothing,Symbol}
    priority::Int
    topology_invalidating::Bool

    function TransformerApparatusEventCommand(
        id::Symbol,
        kind::TransformerApparatusEventKind,
        time_s::Real;
        target_indices=Int[],
        conductance_s::Real=0.0,
        terminal_transform=nothing,
        reference_id::Union{Nothing,Symbol}=nothing,
    )
        id == Symbol("") && throw(ArgumentError(
            "transformer event identity must not be empty",
        ))
        time = Float64(time_s)
        isfinite(time) && time >= 0.0 || throw(ArgumentError(
            "transformer event time must be finite and nonnegative",
        ))
        targets = Int.(target_indices)
        all(>(0), targets) && length(unique(targets)) == length(targets) ||
            throw(ArgumentError(
                "transformer event target indices must be unique positive integers",
            ))
        conductance = Float64(conductance_s)
        isfinite(conductance) && conductance >= 0.0 || throw(ArgumentError(
            "transformer event conductance must be finite and nonnegative",
        ))
        transform = terminal_transform === nothing ? nothing :
            Matrix{Float64}(terminal_transform)
        transform === nothing || all(isfinite, transform) || throw(ArgumentError(
            "transformer event terminal transform must be finite",
        ))
        clear_kind = kind in (
            TransformerTerminalFaultClearEvent,
            TransformerWindingFaultClearEvent,
            TransformerGroundingClearEvent,
            TransformerInternalFaultClearEvent,
        )
        clear_kind == (reference_id !== nothing) || throw(ArgumentError(
            clear_kind ?
                "transformer clear event requires its apply-event reference identity" :
                "transformer non-clear event cannot carry a clear reference identity",
        ))
        transform_kind = kind in (
            TransformerTapChangeEvent,
            TransformerPhaseShiftChangeEvent,
        )
        transform_kind == (transform !== nothing) || throw(ArgumentError(
            transform_kind ?
                "transformer tap/phase event requires a complete terminal transform" :
                "only transformer tap/phase events may carry a terminal transform",
        ))
        apply_fault_kind = kind in (
            TransformerTerminalFaultApplyEvent,
            TransformerWindingFaultApplyEvent,
            TransformerGroundingApplyEvent,
            TransformerInternalFaultApplyEvent,
        )
        apply_fault_kind == (conductance > 0.0) || throw(ArgumentError(
            apply_fault_kind ?
                "transformer fault/grounding application requires positive conductance" :
                "transformer non-fault event conductance must be zero",
        ))
        kind in (TransformerTerminalFaultApplyEvent, TransformerGroundingApplyEvent) &&
            isempty(targets) && throw(ArgumentError(
                "transformer terminal fault/grounding event requires targets",
            ))
        kind in (
            TransformerWindingFaultApplyEvent,
            TransformerInternalFaultApplyEvent,
        ) && length(targets) != 2 && throw(ArgumentError(
            "transformer winding/internal fault event requires exactly two targets",
        ))
        kind in (TransformerBreakerOpenEvent, TransformerBreakerCloseEvent) &&
            isempty(targets) && throw(ArgumentError(
                "transformer breaker event requires terminal targets",
            ))
        return new(
            id,
            kind,
            time,
            targets,
            conductance,
            transform,
            reference_id,
            _transformer_event_priority(kind),
            _transformer_event_topology_invalidating(kind),
        )
    end
end

struct TransformerApparatusEventOccurrence
    id::Symbol
    kind::TransformerApparatusEventKind
    time_s::Float64
    priority::Int
    topology_invalidating::Bool
end

mutable struct TransformerApparatusEventState
    terminal_closed::BitVector
    terminal_transform::Matrix{Float64}
    active_terminal_faults::Dict{Symbol,Tuple{Vector{Int},Float64}}
    active_winding_faults::Dict{Symbol,Tuple{Int,Int,Float64}}
    active_grounding_paths::Dict{Symbol,Tuple{Vector{Int},Float64}}
    active_internal_faults::Dict{Symbol,Tuple{Int,Int,Float64}}
    terminal_fault_admittance_s::Matrix{Float64}
    scheduled_commands::Vector{TransformerApparatusEventCommand}
    applied_event_ids::Vector{Symbol}
    occurrences::Vector{TransformerApparatusEventOccurrence}
    topology_transition_count::Int
    tap_change_count::Int
    phase_shift_change_count::Int
    event_energy_j::Float64
    external_fault_energy_j::Float64
    internal_fault_energy_j::Float64
    numerical_dissipation_energy_j::Float64
    maximum_energy_balance_residual_j::Float64
    previous_network_voltage_v::Vector{Float64}
    previous_fault_current_a::Vector{Float64}
end

function TransformerApparatusEventState(terminal_count::Integer)
    count = Int(terminal_count)
    count > 0 || throw(ArgumentError(
        "transformer event state requires at least one terminal",
    ))
    return TransformerApparatusEventState(
        trues(count),
        Matrix{Float64}(I, count, count),
        Dict{Symbol,Tuple{Vector{Int},Float64}}(),
        Dict{Symbol,Tuple{Int,Int,Float64}}(),
        Dict{Symbol,Tuple{Vector{Int},Float64}}(),
        Dict{Symbol,Tuple{Int,Int,Float64}}(),
        zeros(count, count),
        TransformerApparatusEventCommand[],
        Symbol[],
        TransformerApparatusEventOccurrence[],
        0,
        0,
        0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        zeros(count),
        zeros(count),
    )
end

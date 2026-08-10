module Switches

import ..Branches: EMTElement,
                   backward_euler_companion_supported,
                   stamp!,
                   stamp_conductance!,
                   update!

export IdealSwitch,
       TimeSwitch,
       CurrentZeroSwitch,
       OVER16SwitchAdmittanceState,
       OVER16SwitchAlterationState,
       OVER16SwitchBValueExportState,
       OVER16SwitchCurrentState,
       OVER16SwitchOperationState,
       OVER16SwitchPostCurrentState,
       OVER16SwitchRetriangularizationState,
       OVER16SwitchScanState,
       OVER16SwitchSparseFactorWorkspaceState,
       OVER16FortranSparseFactorWorkspaceState,
       GroupedSparseNetworkWorkspace,
       OVER16SwitchTopologyState,
       switch_closed,
       switch_conductance,
       configure_current_extinction!,
       current_extinction_enabled,
       prepare_current_zero_switch!,
       apply_current_zero_transition!,
       update_current_zero_switch!,
       over16_switch_margin_history,
       over16_switch_close_critical_current,
       over16_gap_energy_check,
       over16_voltage_open_delay,
       over16_switch_tail_current_injection,
       over16_open_switch_close_decision,
       over16_closed_switch_open_decision,
       over16_switch_clamp_action,
       over16_switch_position_update,
       over16_a8sw_delayed_open_step,
       over16_a8sw_tail_current_table,
       over16_controlled_switch_scan_step,
       over16_controlled_switch_table_scan,
       over16_controlled_switch_table_scan!,
       over16_switch_operation_schedule,
       over16_switch_operation_schedule!,
       over16_switch_status_update,
       over16_switch_status_update!,
       over16_switch_simple_ordering,
       over16_switch_simple_ordering!,
       over16_switch_admittance_update,
       over16_switch_admittance_update!,
       over16_switch_topology_admittance_update!,
       over16_switch_retriangularization_update,
       over16_switch_retriangularization_update!,
       over16_switch_retriangularization_solve,
       over16_switch_sparse_factor_matrix,
       over16_switch_sparse_factor_update,
       over16_switch_sparse_factor_update!,
       over16_switch_sparse_factor_solve,
       sparse_network_admittance_matrix,
       sparse_network_factorization_update,
       sparse_network_factorization_update!,
       grouped_sparse_network_factorization_update!,
       sparse_network_solution,
       grouped_sparse_network_solution!,
       sparse_network_solution_update!,
       over16_fortran_sparse_factor_update,
       over16_fortran_sparse_factor_update!,
       over16_fortran_sparse_factor_solve,
       over16_network_solution_update!,
       over16_switch_current_reconstruction,
       over16_switch_current_reconstruction_table,
       over16_switch_current_reconstruction_table!,
       over16_switch_post_current_transition,
       over16_switch_post_current_transition_table,
       over16_switch_post_current_transition_table!,
       over16_switch_alteration_rebuild_intent,
       over16_switch_alteration_rebuild_update!,
       over16_switch_bvalue_export,
       over16_switch_bvalue_export!

struct IdealSwitch <: EMTElement
    a::Int
    b::Int
    closed::Bool
    on_conductance::Float64
    off_conductance::Float64
end

function IdealSwitch(a::Int, b::Int, closed::Bool; on_conductance::Real=1.0e9,
                     off_conductance::Real=0.0)
    return IdealSwitch(a, b, closed, Float64(on_conductance), Float64(off_conductance))
end

mutable struct CurrentExtinctionState
    not_before_time_s::Float64
    critical_current_a::Float64
    closed::Bool
    opened::Bool
    current_initialized::Bool
    previous_current::Float64
    operation_count::Int
    open_reason::Symbol
    opened_time_s::Float64
end

mutable struct TimeSwitch <: EMTElement
    a::Int
    b::Int
    close_time_s::Float64
    open_time_s::Float64
    initially_closed::Bool
    on_conductance::Float64
    off_conductance::Float64
    current_extinction::Union{Nothing,CurrentExtinctionState}
end

mutable struct CurrentZeroSwitch <: EMTElement
    a::Int
    b::Int
    close_time_s::Float64
    open_request_time_s::Float64
    open_delay_time_s::Float64
    critical_current_a::Float64
    initially_closed::Bool
    on_conductance::Float64
    off_conductance::Float64
    closed::Bool
    opened::Bool
    current_initialized::Bool
    previous_current::Float64
    operation_count::Int
    open_reason::Symbol
end

function CurrentZeroSwitch(
    switch::TimeSwitch;
    open_delay_time_s::Real = 0.0,
    critical_current_a::Real = 0.0,
)
    delay_time = Float64(open_delay_time_s)
    critical_current = Float64(critical_current_a)
    delay_time >= 0.0 && (isfinite(delay_time) || delay_time == Inf) ||
        throw(ArgumentError("current-extinction delay time must be nonnegative"))
    critical_current >= 0.0 && isfinite(critical_current) ||
        throw(ArgumentError("critical current must be finite and nonnegative"))
    initially_active = switch_closed(switch, 0.0)
    return CurrentZeroSwitch(
        switch.a,
        switch.b,
        switch.close_time_s,
        switch.open_time_s,
        delay_time,
        critical_current,
        switch.initially_closed,
        switch.on_conductance,
        switch.off_conductance,
        initially_active,
        false,
        false,
        0.0,
        0,
        :none,
    )
end

mutable struct OVER16SwitchScanState
    positions::Vector{Int}
    elapsed_open_times::Vector{Float64}
    gap_currents::Vector{Float64}
    modswt::Vector{Int}
end

function OVER16SwitchScanState(
    positions::AbstractVector{Int},
    elapsed_open_times::AbstractVector{<:Real};
    gap_currents::AbstractVector{<:Real}=Float64[],
    modswt::AbstractVector{Int}=Int[],
)
    length(positions) == length(elapsed_open_times) ||
        throw(ArgumentError("positions and elapsed_open_times lengths must match"))
    return OVER16SwitchScanState(
        collect(positions),
        Float64.(elapsed_open_times),
        Float64.(gap_currents),
        collect(modswt),
    )
end

mutable struct OVER16SwitchOperationState
    modswt::Vector{Int}
    closed_switch_count::Int
    accumulated_operation_count::Int
end

function OVER16SwitchOperationState(
    modswt::AbstractVector{Int},
    closed_switch_count::Int;
    accumulated_operation_count::Int=0,
)
    closed_switch_count >= 0 ||
        throw(ArgumentError("closed_switch_count must be nonnegative"))
    accumulated_operation_count >= 0 ||
        throw(ArgumentError("accumulated_operation_count must be nonnegative"))
    return OVER16SwitchOperationState(
        collect(modswt),
        closed_switch_count,
        accumulated_operation_count,
    )
end

mutable struct OVER16SwitchTopologyState
    closed_mask::Vector{Bool}
    closed_switch_count::Int
    first_group_head::Int
    nextsw::Vector{Int}
    kode::Vector{Int}
end

function OVER16SwitchTopologyState(
    closed_mask::AbstractVector{Bool};
    closed_switch_count::Int=count(identity, closed_mask),
    first_group_head::Int=0,
    nextsw::AbstractVector{Int}=zeros(Int, length(closed_mask)),
    kode::AbstractVector{Int}=Int[],
)
    switch_count = length(closed_mask)
    closed_switch_count == count(identity, closed_mask) ||
        throw(ArgumentError("closed_switch_count must match closed_mask"))
    0 <= first_group_head <= switch_count ||
        throw(ArgumentError("first_group_head must be between 0 and switch count"))
    length(nextsw) == switch_count ||
        throw(ArgumentError("nextsw length must match closed_mask"))
    return OVER16SwitchTopologyState(
        collect(closed_mask),
        closed_switch_count,
        first_group_head,
        collect(nextsw),
        collect(kode),
    )
end

mutable struct OVER16SwitchAdmittanceState
    base_admittance::Matrix{Float64}
    admittance::Matrix{Float64}
    admittance_workspace::Matrix{Float64}
    switch_conductances::Vector{Float64}
    retriangularization_count::Int
end

struct SwitchOperationStepResult
    processed_modswt::Vector{Int}
    ktrlsw_count::Int
    switch_operation_state_mutated::Bool
end

struct SwitchStatusStepResult
    requires_order_rebuild::Bool
    switch_status_state_mutated::Bool
end

struct SwitchOrderStepResult
    switch_order_state_mutated::Bool
end

struct SwitchAdmittanceStepResult
    should_retriangularize::Bool
    admittance_mutated::Bool
    switch_admittance_state_mutated::Bool
end

struct SwitchTopologyAdmittanceStepResult
    status_result::SwitchStatusStepResult
    order_result::Union{Nothing,SwitchOrderStepResult}
    admittance_result::SwitchAdmittanceStepResult
    switch_topology_state_mutated::Bool
    switch_admittance_state_mutated::Bool
    switch_topology_admittance_state_mutated::Bool
    topology_mutated::Bool
    admittance_mutated::Bool
    should_retriangularize::Bool
end

function OVER16SwitchAdmittanceState(
    base_admittance::AbstractMatrix{<:Real};
    switch_conductances::AbstractVector{<:Real}=Float64[],
    retriangularization_count::Int=0,
)
    size(base_admittance, 1) == size(base_admittance, 2) ||
        throw(ArgumentError("base_admittance must be square"))
    retriangularization_count >= 0 ||
        throw(ArgumentError("retriangularization_count must be nonnegative"))
    base = Float64.(base_admittance)
    for value in base
        isfinite(value) ||
            throw(ArgumentError("base_admittance entries must be finite"))
    end
    conductances = Float64.(switch_conductances)
    for value in conductances
        isfinite(value) && value >= 0.0 ||
            throw(ArgumentError("switch_conductances entries must be finite and nonnegative"))
    end
    return OVER16SwitchAdmittanceState(
        base,
        copy(base),
        copy(base),
        conductances,
        retriangularization_count,
    )
end

function OVER16SwitchAdmittanceState(
    node_count::Int;
    switch_count::Int=0,
    retriangularization_count::Int=0,
)
    node_count > 0 || throw(ArgumentError("node_count must be positive"))
    switch_count >= 0 || throw(ArgumentError("switch_count must be nonnegative"))
    return OVER16SwitchAdmittanceState(
        zeros(Float64, node_count, node_count);
        switch_conductances = zeros(Float64, switch_count),
        retriangularization_count = retriangularization_count,
    )
end

mutable struct OVER16SwitchRetriangularizationState
    factor::Matrix{Float64}
    pivot_values::Vector{Float64}
    factorization_count::Int
end

function OVER16SwitchRetriangularizationState(
    factor::AbstractMatrix{<:Real};
    pivot_values::AbstractVector{<:Real}=Float64[],
    factorization_count::Int=0,
)
    size(factor, 1) == size(factor, 2) ||
        throw(ArgumentError("factor must be square"))
    factorization_count >= 0 ||
        throw(ArgumentError("factorization_count must be nonnegative"))
    dense_factor = Float64.(factor)
    for value in dense_factor
        isfinite(value) ||
            throw(ArgumentError("factor entries must be finite"))
    end
    pivots = isempty(pivot_values) ? zeros(Float64, size(factor, 1)) : Float64.(pivot_values)
    length(pivots) == size(factor, 1) ||
        throw(ArgumentError("pivot_values length must match factor size"))
    for value in pivots
        isfinite(value) ||
            throw(ArgumentError("pivot_values entries must be finite"))
    end
    return OVER16SwitchRetriangularizationState(dense_factor, pivots, factorization_count)
end

function OVER16SwitchRetriangularizationState(
    node_count::Int;
    factorization_count::Int=0,
)
    node_count > 0 || throw(ArgumentError("node_count must be positive"))
    return OVER16SwitchRetriangularizationState(
        zeros(Float64, node_count, node_count);
        pivot_values = zeros(Float64, node_count),
        factorization_count = factorization_count,
    )
end

mutable struct OVER16SwitchSparseFactorWorkspaceState
    km::Vector{Int}
    ykm::Vector{Float64}
    kk::Vector{Int}
    workspace_update_count::Int
end

function OVER16SwitchSparseFactorWorkspaceState(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int};
    workspace_update_count::Int=0,
)
    length(km) == length(ykm) ||
        throw(ArgumentError("km and ykm lengths must match"))
    workspace_update_count >= 0 ||
        throw(ArgumentError("workspace_update_count must be nonnegative"))
    dense_ykm = Float64.(ykm)
    for value in dense_ykm
        isfinite(value) ||
            throw(ArgumentError("ykm entries must be finite"))
    end
    return OVER16SwitchSparseFactorWorkspaceState(
        collect(km),
        dense_ykm,
        collect(kk),
        workspace_update_count,
    )
end

function OVER16SwitchSparseFactorWorkspaceState(
    node_count::Int;
    workspace_update_count::Int=0,
)
    node_count > 0 || throw(ArgumentError("node_count must be positive"))
    return OVER16SwitchSparseFactorWorkspaceState(
        Int[],
        Float64[],
        zeros(Int, node_count);
        workspace_update_count = workspace_update_count,
    )
end

mutable struct OVER16FortranSparseFactorWorkspaceState
    km::Vector{Int}
    ykm::Vector{Float64}
    kk::Vector{Int}
    workspace_update_count::Int
end

"""
Reusable storage for the grouped sparse factorization and solve kernels.

The workspace owns every temporary whose size is determined by the augmented
network. A caller may refactor it whenever admittance or node grouping changes,
then solve any number of right-hand sides without allocating.
"""
mutable struct GroupedSparseNetworkWorkspace
    km::Vector{Int}
    ykm::Vector{Float64}
    kk::Vector{Int}
    pivot_values::Vector{Float64}
    inverse_diagonal_values::Vector{Float64}
    contracted_admittance::Matrix{Float64}
    node_group_representatives::Vector{Int}
    active_rows::Vector{Int}
    row_starts::Vector{Int}
    factor_row::Vector{Float64}
    visit_marks::Vector{Int}
    visit_positions::Vector{Int}
    visit_path::Vector{Int}
    alias_successors::Vector{Int}
    solution::Vector{Float64}
    factorization_count::Int
end

function GroupedSparseNetworkWorkspace(node_count::Int)
    node_count > 1 || throw(ArgumentError(
        "grouped sparse workspace must include a reference and at least one node",
    ))
    km = Int[]
    ykm = Float64[]
    active_rows = Int[]
    pivot_values = Float64[]
    inverse_diagonal_values = Float64[]
    sizehint!(km, node_count * node_count)
    sizehint!(ykm, node_count * node_count)
    sizehint!(active_rows, node_count)
    sizehint!(pivot_values, node_count)
    sizehint!(inverse_diagonal_values, node_count)
    return GroupedSparseNetworkWorkspace(
        km,
        ykm,
        zeros(Int, node_count),
        pivot_values,
        inverse_diagonal_values,
        zeros(Float64, node_count, node_count),
        collect(1:node_count),
        active_rows,
        zeros(Int, node_count),
        zeros(Float64, node_count),
        zeros(Int, node_count),
        zeros(Int, node_count),
        zeros(Int, node_count),
        zeros(Int, node_count),
        zeros(Float64, node_count),
        0,
    )
end

function OVER16FortranSparseFactorWorkspaceState(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int};
    workspace_update_count::Int=0,
)
    length(km) == length(ykm) ||
        throw(ArgumentError("km and ykm lengths must match"))
    workspace_update_count >= 0 ||
        throw(ArgumentError("workspace_update_count must be nonnegative"))
    values = Float64.(ykm)
    for value in values
        isfinite(value) ||
            throw(ArgumentError("ykm entries must be finite"))
    end
    return OVER16FortranSparseFactorWorkspaceState(
        collect(km),
        values,
        collect(kk),
        workspace_update_count,
    )
end

function OVER16FortranSparseFactorWorkspaceState(
    node_count::Int;
    workspace_update_count::Int=0,
)
    node_count > 1 || throw(ArgumentError("node_count must be greater than one"))
    return OVER16FortranSparseFactorWorkspaceState(
        Int[],
        Float64[],
        zeros(Int, node_count);
        workspace_update_count = workspace_update_count,
    )
end

mutable struct OVER16SwitchCurrentState
    rhs::Vector{Float64}
    switch_currents::Vector{Float64}
    current_products::Vector{Float64}
    network_solution::Vector{Float64}
end

function OVER16SwitchCurrentState(
    rhs::AbstractVector{<:Real},
    switch_currents::AbstractVector{<:Real};
    current_products::AbstractVector{<:Real}=zeros(Float64, length(switch_currents)),
    network_solution::AbstractVector{<:Real}=Float64[],
)
    length(current_products) == length(switch_currents) ||
        throw(ArgumentError("current_products length must match switch_currents"))
    isempty(network_solution) || length(network_solution) == length(rhs) ||
        throw(ArgumentError("network_solution length must match rhs when provided"))
    return OVER16SwitchCurrentState(
        Float64.(rhs),
        Float64.(switch_currents),
        Float64.(current_products),
        Float64.(network_solution),
    )
end

mutable struct OVER16SwitchPostCurrentState
    positions::Vector{Int}
    switch_currents::Vector{Float64}
    energies::Vector{Float64}
    modswt::Vector{Int}
end

function OVER16SwitchPostCurrentState(
    positions::AbstractVector{Int},
    switch_currents::AbstractVector{<:Real},
    energies::AbstractVector{<:Real};
    modswt::AbstractVector{Int}=Int[],
)
    switch_count = length(positions)
    length(switch_currents) == switch_count ||
        throw(ArgumentError("switch_currents length must match positions"))
    length(energies) == switch_count ||
        throw(ArgumentError("energies length must match positions"))
    return OVER16SwitchPostCurrentState(
        collect(positions),
        Float64.(switch_currents),
        Float64.(energies),
        collect(modswt),
    )
end

mutable struct OVER16SwitchBValueExportState
    bvalue::Vector{Float64}
    output_count::Int
end

function OVER16SwitchBValueExportState(
    bvalue::AbstractVector{<:Real};
    output_count::Int=length(bvalue),
)
    output_count >= 0 || throw(ArgumentError("output_count must be nonnegative"))
    output_count <= length(bvalue) ||
        throw(ArgumentError("output_count cannot exceed bvalue length"))
    return OVER16SwitchBValueExportState(Float64.(bvalue), output_count)
end

mutable struct OVER16SwitchAlterationState
    ialter::Int
    operation_count::Int
    closed_switch_count::Int
    triangularization_count::Int
    first_group_head::Int
    total_operation_count::Int
    ktrlsw6::Int
end

function OVER16SwitchAlterationState(
    ialter::Int,
    operation_count::Int,
    closed_switch_count::Int,
    triangularization_count::Int;
    first_group_head::Int=0,
    total_operation_count::Int=0,
    ktrlsw6::Int=1,
)
    over16_switch_alteration_rebuild_intent(
        ialter,
        operation_count,
        closed_switch_count,
        triangularization_count,
        first_group_head,
        total_operation_count,
        ktrlsw6,
    )
    return OVER16SwitchAlterationState(
        ialter,
        operation_count,
        closed_switch_count,
        triangularization_count,
        first_group_head,
        total_operation_count,
        ktrlsw6,
    )
end

function TimeSwitch(a::Int, b::Int; close_time_s::Real=Inf, open_time_s::Real=Inf,
                    initially_closed::Bool=false, on_conductance::Real=1.0e9,
                    off_conductance::Real=0.0)
    return TimeSwitch(
        a,
        b,
        Float64(close_time_s),
        Float64(open_time_s),
        initially_closed,
        Float64(on_conductance),
        Float64(off_conductance),
        nothing,
    )
end

switch_closed(s::IdealSwitch, _t::Real)::Bool = s.closed

function switch_closed(s::TimeSwitch, t::Real)::Bool
    s.current_extinction === nothing ||
        return s.current_extinction.closed
    time = Float64(t)
    closed = s.initially_closed
    if isfinite(s.close_time_s) && time >= s.close_time_s
        closed = true
    end
    if isfinite(s.open_time_s) && time >= s.open_time_s
        closed = false
    end
    return closed
end

function switch_conductance(s::Union{IdealSwitch,TimeSwitch}, t::Real)::Float64
    return switch_closed(s, t) ? s.on_conductance : s.off_conductance
end

current_extinction_enabled(::Any) = false
current_extinction_enabled(::CurrentZeroSwitch) = true
current_extinction_enabled(s::TimeSwitch) = s.current_extinction !== nothing
current_extinction_enabled(::Tuple{}) = false
current_extinction_enabled(elements::Tuple) =
    any(current_extinction_enabled, elements)

function configure_current_extinction!(
    switch::TimeSwitch,
    not_before_time_s::Real,
    critical_current_a::Real,
    time_s::Real,
    ;
    currently_closed::Bool = switch_closed(switch, time_s),
)
    not_before = Float64(not_before_time_s)
    critical_current = Float64(critical_current_a)
    time = Float64(time_s)
    not_before >= 0.0 && (isfinite(not_before) || not_before == Inf) ||
        throw(ArgumentError("current-extinction not-before time must be nonnegative"))
    critical_current >= 0.0 && isfinite(critical_current) ||
        throw(ArgumentError("critical current must be finite and nonnegative"))
    if not_before == 0.0 && critical_current == 0.0
        switch.current_extinction = nothing
        return switch
    end
    scheduled_closed =
        currently_closed ||
        (
            isfinite(switch.close_time_s) &&
            time >= switch.close_time_s &&
            time < switch.open_time_s
        )
    previous = switch.current_extinction
    reclosed = !currently_closed && scheduled_closed
    switch.current_extinction = CurrentExtinctionState(
        not_before,
        critical_current,
        scheduled_closed,
        reclosed || previous === nothing ? false : previous.opened,
        reclosed || previous === nothing ? false : previous.current_initialized,
        reclosed || previous === nothing ? 0.0 : previous.previous_current,
        (previous === nothing ? 0 : previous.operation_count) + (reclosed ? 1 : 0),
        reclosed ? :restart_topology_change :
            previous === nothing ? :none : previous.open_reason,
        reclosed || previous === nothing ? Inf : previous.opened_time_s,
    )
    return switch
end

switch_closed(s::CurrentZeroSwitch, _t::Real)::Bool = s.closed

function switch_conductance(s::CurrentZeroSwitch, _t::Real)::Float64
    return s.closed ? s.on_conductance : s.off_conductance
end

function prepare_current_zero_switch!(switch::CurrentZeroSwitch, time_s::Real)
    time = Float64(time_s)
    if !switch.opened && !switch.closed &&
       isfinite(switch.close_time_s) && time >= switch.close_time_s
        switch.closed = true
        switch.operation_count += 1
    end
    return switch
end

function prepare_current_zero_switch!(switch::TimeSwitch, time_s::Real)
    state = switch.current_extinction
    state === nothing && return switch
    time = Float64(time_s)
    if !state.opened && !state.closed &&
       isfinite(switch.close_time_s) && time >= switch.close_time_s
        state.closed = true
        state.operation_count += 1
    end
    return switch
end

function _check_current_zero_transition_reason(reason::Symbol)
    reason in (:current_reversal, :critical_current) || throw(ArgumentError(
        "current-zero transition reason must be :current_reversal or :critical_current",
    ))
    return reason
end

function apply_current_zero_transition!(
    switch::CurrentZeroSwitch,
    reason::Symbol,
    time_s::Real,
)
    _check_current_zero_transition_reason(reason)
    time = Float64(time_s)
    isfinite(time) || throw(ArgumentError("current-zero transition time must be finite"))
    if switch.closed
        switch.closed = false
        switch.opened = true
        switch.operation_count += 1
        switch.open_reason = reason
    end
    return switch
end

function apply_current_zero_transition!(
    switch::TimeSwitch,
    reason::Symbol,
    time_s::Real,
)
    state = switch.current_extinction
    state === nothing && throw(ArgumentError(
        "time switch has no current-extinction transition owner",
    ))
    _check_current_zero_transition_reason(reason)
    time = Float64(time_s)
    isfinite(time) || throw(ArgumentError("current-zero transition time must be finite"))
    if state.closed
        state.closed = false
        state.opened = true
        state.operation_count += 1
        state.open_reason = reason
        state.opened_time_s = time
    end
    return switch
end

function update_current_zero_switch!(
    switch::CurrentZeroSwitch,
    current::Real,
    time_s::Real,
)
    switch.closed || return switch
    current_value = Float64(current)
    time = Float64(time_s)
    current_reversed =
        switch.current_initialized &&
        current_value * switch.previous_current < 0.0
    critical_reached =
        switch.critical_current_a > 0.0 &&
        abs(current_value) < switch.critical_current_a
    if switch.current_initialized &&
       time >= switch.open_request_time_s &&
       time >= switch.open_delay_time_s &&
       (current_reversed || critical_reached)
        apply_current_zero_transition!(
            switch,
            critical_reached ? :critical_current : :current_reversal,
            time,
        )
    end
    switch.previous_current = current_value
    switch.current_initialized = true
    return switch
end

function update_current_zero_switch!(
    switch::TimeSwitch,
    current::Real,
    time_s::Real,
)
    state = switch.current_extinction
    state === nothing && return switch
    state.closed || return switch
    current_value = Float64(current)
    time = Float64(time_s)
    current_reversed =
        state.current_initialized &&
        current_value * state.previous_current < 0.0
    critical_reached =
        state.critical_current_a > 0.0 &&
        abs(current_value) < state.critical_current_a
    if state.current_initialized &&
       time >= switch.open_time_s &&
       time >= state.not_before_time_s &&
       (current_reversed || critical_reached)
        apply_current_zero_transition!(
            switch,
            critical_reached ? :critical_current : :current_reversal,
            time,
        )
    end
    state.previous_current = current_value
    state.current_initialized = true
    return switch
end

function over16_switch_margin_history(
    previous_ck::Float64,
    voltage_difference::Float64,
    delta2::Float64,
    critical_current::Float64,
    branch_multiplier::Float64,
    cik_current::Float64,
    cik_next::Float64,
)
    area = voltage_difference * delta2
    midpoint_ck = previous_ck + area
    updated_ck = midpoint_ck + area
    should_update = !(abs(midpoint_ck) >= critical_current || abs(midpoint_ck) > abs(previous_ck))
    if !should_update
        return updated_ck, cik_current, cik_next, false
    end

    history_current = branch_multiplier * voltage_difference
    return updated_ck, -history_current, history_current + cik_current + cik_next, true
end

function over16_switch_close_critical_current(
    previous_ck::Float64,
    voltage_difference::Float64,
    delta2::Float64,
    open_threshold::Float64,
    branch_multiplier::Float64,
    branch_delay_multiplier::Float64,
    cik_next::Float64,
    critical_current::Float64,
)
    delta2 > 0.0 || throw(ArgumentError("delta2 must be positive"))
    branch_delay_multiplier != 0.0 || throw(ArgumentError("branch_delay_multiplier must be nonzero"))

    history_current = branch_multiplier * voltage_difference
    area = voltage_difference * delta2
    midpoint_ck = previous_ck + area
    updated_ck = midpoint_ck + area
    should_update = !(abs(midpoint_ck) < open_threshold || abs(midpoint_ck) < abs(previous_ck))
    if !should_update
        return updated_ck, -history_current, critical_current, false
    end

    d2 = branch_delay_multiplier / delta2
    multiplier_ratio = branch_multiplier / (delta2 * d2)
    multiplier_ratio != 0.0 || throw(ArgumentError("branch_multiplier must be nonzero"))
    ci1 = cik_next + d2 * area
    d3 = midpoint_ck * (multiplier_ratio + 1.0) - ci1 / d2
    return updated_ck, -history_current, abs(d3 / multiplier_ratio), true
end

function over16_gap_energy_check(previous_energy::Float64, voltage_difference::Float64)
    average_voltage = (voltage_difference + previous_energy) / 2.0
    return average_voltage, voltage_difference, average_voltage >= 0.0
end

function over16_voltage_open_delay(
    voltage_difference::Float64,
    open_threshold::Float64,
    t::Float64,
    delay_offset::Float64,
    current_delay::Float64,
)
    if abs(voltage_difference) < open_threshold
        return current_delay, false
    end
    return t + delay_offset, true
end

function over16_switch_tail_current_injection(
    decay_end_time::Float64,
    t::Float64,
    dt::Float64,
    amplitude::Float64,
    time_constant::Float64,
    cutoff_current::Float64,
)
    dt > 0.0 || throw(ArgumentError("dt must be positive"))
    time_constant != 0.0 || throw(ArgumentError("time_constant must be nonzero"))

    decay_factor = exp((decay_end_time - t - dt) / time_constant)
    current = amplitude * decay_factor
    clear_marker = current <= cutoff_current ||
                   (decay_factor <= 1.0e-4 && cutoff_current == 0.0)
    return current, -current, decay_factor, clear_marker
end

function over16_open_switch_close_decision(
    elapsed_open_time::Float64,
    dt::Float64,
    maximum_open_time::Float64,
    voltage_difference::Float64,
    close_threshold::Float64;
    is_gap::Bool=false,
    control_enabled::Bool=false,
    control_signal::Float64=0.0,
    flzero::Float64=0.0,
    fltinf::Float64=Inf,
)
    dt > 0.0 || throw(ArgumentError("dt must be positive"))

    updated_elapsed = elapsed_open_time + dt
    if updated_elapsed > maximum_open_time
        updated_elapsed = -fltinf
    end

    threshold_voltage = is_gap ? abs(voltage_difference) : voltage_difference
    if threshold_voltage < close_threshold
        return updated_elapsed, false
    end

    if !is_gap && updated_elapsed >= 0.0
        return -fltinf, true
    end

    if !control_enabled
        return updated_elapsed, true
    end
    return updated_elapsed, control_signal > 10.0 * flzero
end

function over16_closed_switch_open_decision(
    switch_current::Float64,
    previous_gap_current::Float64,
    open_threshold::Float64;
    is_gap::Bool=false,
    gap_control_signal::Float64=0.0,
    flzero::Float64=0.0,
)
    if is_gap
        if gap_control_signal > 10.0 * flzero
            return previous_gap_current, false
        end
        updated_gap_current = switch_current
        if switch_current * previous_gap_current < 0.0
            return updated_gap_current, true
        end
        return updated_gap_current, abs(switch_current) < open_threshold
    end
    return previous_gap_current, switch_current < open_threshold
end

function over16_switch_clamp_action(
    switch_type::Int,
    position::Int,
    controlled_node::Int,
    control_index::Int,
    control_signal::Float64,
    flzero::Float64,
)
    m1 = abs(position)
    tolerance = 10.0 * flzero
    if control_index != 0
        if control_signal < -tolerance
            return m1 == 2 ? :open : :skip
        elseif control_signal > tolerance
            return m1 == 2 ? :skip : :close
        end
    elseif switch_type == 8891 || (switch_type == 8890 && controlled_node == 0)
        return :skip
    end

    if switch_type == 8891
        return m1 == 2 ? :open : :skip
    elseif switch_type == 8890 && controlled_node == 0
        return :skip
    elseif m1 == 2
        return :evaluate_closed
    end
    return :evaluate_open
end

function over16_switch_position_update(
    position::Int,
    action::Symbol,
    elapsed_open_time::Float64;
    reset_elapsed::Bool=false,
)
    if action == :open
        magnitude = 5
        modifier = -1
    elseif action == :close
        magnitude = 2
        modifier = 1
    else
        throw(ArgumentError("action must be :open or :close"))
    end
    signed_position = position < 0 ? -magnitude : magnitude
    updated_elapsed = action == :open && reset_elapsed ? 0.0 : elapsed_open_time
    return signed_position, modifier, updated_elapsed
end

function over16_a8sw_delayed_open_step(
    switch_current::Float64,
    previous_current::Float64,
    shape_current::Float64,
    shape_delay::Float64,
    scheduled_open_time::Float64,
    t::Float64,
    dt::Float64,
    coefficient::Float64,
    exponent::Float64,
    time_scale::Float64;
    threshold_marker::Float64=9999.0,
)
    dt > 0.0 || throw(ArgumentError("dt must be positive"))
    time_scale != 0.0 || throw(ArgumentError("time_scale must be nonzero"))

    updated_shape_current = shape_current
    updated_shape_delay = shape_delay
    updated_scheduled_open_time = scheduled_open_time

    if switch_current > 0.0
        delta_current = previous_current - switch_current
        if delta_current < switch_current
            return (
                previous_current = switch_current,
                shape_current = updated_shape_current,
                shape_delay = updated_shape_delay,
                scheduled_open_time = updated_scheduled_open_time,
                open_threshold = nothing,
                opens = false,
            )
        end
        didt = delta_current / dt
        updated_shape_current = didt^exponent * coefficient
        updated_shape_delay = updated_shape_current / (didt * time_scale)
        updated_scheduled_open_time =
            t + switch_current * dt / delta_current + updated_shape_delay * time_scale
    elseif previous_current != 0.0 && previous_current * switch_current <= 0.0
        delta_current = previous_current - switch_current
        didt = delta_current / dt
        updated_shape_current = didt^exponent * coefficient
        updated_shape_delay = updated_shape_current / (didt * time_scale)
        updated_scheduled_open_time =
            t + switch_current * dt / delta_current + updated_shape_delay * time_scale
    end

    if t + dt < updated_scheduled_open_time
        return (
            previous_current = switch_current,
            shape_current = updated_shape_current,
            shape_delay = updated_shape_delay,
            scheduled_open_time = updated_scheduled_open_time,
            open_threshold = nothing,
            opens = false,
        )
    end

    return (
        previous_current = previous_current,
        shape_current = updated_shape_current,
        shape_delay = updated_shape_delay,
        scheduled_open_time = updated_scheduled_open_time,
        open_threshold = threshold_marker,
        opens = true,
    )
end

function over16_a8sw_tail_current_table(
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    switch_types::AbstractVector{Int},
    positions::AbstractVector{Int},
    tail_indices::AbstractVector{Int},
    tail_markers::AbstractVector{<:Real},
    scheduled_open_times::AbstractVector{<:Real},
    amplitudes::AbstractVector{<:Real},
    time_constants::AbstractVector{<:Real},
    cutoff_currents::AbstractVector{<:Real},
    t::Real,
    dt::Real;
    node_count::Int=max(maximum(vcat([1], collect(from_nodes), collect(to_nodes))), 1),
    tail_marker_value::Real=9999.0,
)
    switch_count = length(switch_types)
    _over16_check_length("from_nodes", from_nodes, switch_count)
    _over16_check_length("to_nodes", to_nodes, switch_count)
    _over16_check_length("positions", positions, switch_count)
    _over16_check_length("tail_indices", tail_indices, switch_count)
    _over16_check_length("tail_markers", tail_markers, switch_count)
    _over16_check_length("scheduled_open_times", scheduled_open_times, switch_count)
    _over16_check_length("amplitudes", amplitudes, switch_count)
    _over16_check_length("time_constants", time_constants, switch_count)
    _over16_check_length("cutoff_currents", cutoff_currents, switch_count)
    node_count >= 1 || throw(ArgumentError("node_count must be positive"))

    rhs_assignments = zeros(Float64, node_count)
    tail_currents = zeros(Float64, switch_count)
    decay_factors = zeros(Float64, switch_count)
    updated_tail_markers = Float64.(tail_markers)
    active_rows = Int[]
    cleared_rows = Int[]

    for row in 1:switch_count
        from_node = from_nodes[row]
        to_node = to_nodes[row]
        1 <= from_node <= node_count ||
            throw(ArgumentError("from_nodes entries must be within node_count"))
        1 <= to_node <= node_count ||
            throw(ArgumentError("to_nodes entries must be within node_count"))
        from_node != to_node ||
            throw(ArgumentError("tail-current switch endpoints must be distinct"))

        switch_type = switch_types[row]
        if !(switch_type == 8888 || switch_type == 8890 || switch_type == 8891)
            continue
        end
        if tail_indices[row] <= 0 || abs(positions[row]) != 5 ||
           Float64(tail_markers[row]) != Float64(tail_marker_value)
            continue
        end

        current_from, current_to, decay_factor, clear_marker =
            over16_switch_tail_current_injection(
                Float64(scheduled_open_times[row]),
                Float64(t),
                Float64(dt),
                Float64(amplitudes[row]),
                Float64(time_constants[row]),
                Float64(cutoff_currents[row]),
            )
        rhs_assignments[from_node] = current_from
        rhs_assignments[to_node] = current_to
        tail_currents[row] = current_from
        decay_factors[row] = decay_factor
        push!(active_rows, row)
        if clear_marker
            updated_tail_markers[row] = 0.0
            push!(cleared_rows, row)
        end
    end

    return (
        rhs_assignments = rhs_assignments,
        tail_currents = tail_currents,
        decay_factors = decay_factors,
        updated_tail_markers = updated_tail_markers,
        active_rows = active_rows,
        cleared_rows = cleared_rows,
        rhs_mutated = false,
        topology_mutated = false,
        admittance_mutated = false,
    )
end

function over16_controlled_switch_scan_step(
    switch_type::Int,
    position::Int,
    elapsed_open_time::Float64,
    maximum_open_time::Float64,
    voltage_difference::Float64,
    switch_current::Float64,
    close_threshold::Float64,
    open_threshold::Float64,
    dt::Float64;
    controlled_node::Int=0,
    control_index::Int=0,
    clamp_signal::Float64=0.0,
    close_control_signal::Float64=0.0,
    gap_control_signal::Float64=0.0,
    previous_gap_current::Float64=0.0,
    flzero::Float64=0.0,
    fltinf::Float64=Inf,
)
    scan_action = over16_switch_clamp_action(
        switch_type,
        position,
        controlled_node,
        control_index,
        clamp_signal,
        flzero,
    )

    transition = :none
    updated_elapsed = elapsed_open_time
    updated_gap_current = previous_gap_current
    if scan_action == :skip
        transition = :none
    elseif scan_action == :open || scan_action == :close
        transition = scan_action
    elseif scan_action == :evaluate_closed
        updated_gap_current, opens = over16_closed_switch_open_decision(
            switch_current,
            previous_gap_current,
            open_threshold;
            is_gap = switch_type == 8890,
            gap_control_signal = gap_control_signal,
            flzero = flzero,
        )
        transition = opens ? :open : :none
    elseif scan_action == :evaluate_open
        updated_elapsed, closes = over16_open_switch_close_decision(
            elapsed_open_time,
            dt,
            maximum_open_time,
            voltage_difference,
            close_threshold;
            is_gap = switch_type == 8890,
            control_enabled = controlled_node != 0,
            control_signal = close_control_signal,
            flzero = flzero,
            fltinf = fltinf,
        )
        transition = closes ? :close : :none
    else
        throw(ArgumentError("unsupported switch scan action"))
    end

    if transition == :none
        return (
            scan_action = scan_action,
            transition = transition,
            position = position,
            elapsed_open_time = updated_elapsed,
            gap_current = updated_gap_current,
            modswt_sign = 0,
            altered = false,
        )
    end

    updated_position, modswt_sign, final_elapsed = over16_switch_position_update(
        position,
        transition,
        updated_elapsed;
        reset_elapsed = transition == :open && maximum_open_time != 0.0,
    )
    return (
        scan_action = scan_action,
        transition = transition,
        position = updated_position,
        elapsed_open_time = final_elapsed,
        gap_current = updated_gap_current,
        modswt_sign = modswt_sign,
        altered = true,
    )
end

function _over16_check_length(name::String, values::AbstractVector, expected::Int)
    length(values) == expected || throw(ArgumentError("$name length must be $expected"))
    return nothing
end

function _over16_check_optional_length(name::String, values::AbstractVector, expected::Int)
    (isempty(values) || length(values) == expected) ||
        throw(ArgumentError("$name length must be 0 or $expected"))
    return nothing
end

function _over16_optional_float(values::AbstractVector{Float64}, index::Int, default::Float64)
    return isempty(values) ? default : values[index]
end

function _over16_switch_conductance_vector(
    name::String,
    values::AbstractVector{<:Real},
    switch_count::Int,
    default::Real,
)
    default_value = Float64(default)
    isfinite(default_value) && default_value >= 0.0 ||
        throw(ArgumentError("$(name) default must be finite and nonnegative"))
    if isempty(values)
        return fill(default_value, switch_count)
    end
    _over16_check_length(name, values, switch_count)
    result = Float64.(values)
    for value in result
        isfinite(value) && value >= 0.0 ||
            throw(ArgumentError("$(name) entries must be finite and nonnegative"))
    end
    return result
end

function _over16_optional_int(values::AbstractVector{Int}, index::Int, default::Int)
    return isempty(values) ? default : values[index]
end

function over16_controlled_switch_table_scan(
    switch_types::AbstractVector{Int},
    positions::AbstractVector{Int},
    elapsed_open_times::AbstractVector{Float64},
    maximum_open_times::AbstractVector{Float64},
    voltage_differences::AbstractVector{Float64},
    switch_currents::AbstractVector{Float64},
    close_thresholds::AbstractVector{Float64},
    open_thresholds::AbstractVector{Float64},
    dt::Float64;
    controlled_nodes::AbstractVector{Int}=Int[],
    control_indices::AbstractVector{Int}=Int[],
    clamp_signals::AbstractVector{Float64}=Float64[],
    close_control_signals::AbstractVector{Float64}=Float64[],
    gap_control_signals::AbstractVector{Float64}=Float64[],
    previous_gap_currents::AbstractVector{Float64}=Float64[],
    flzero::Float64=0.0,
    fltinf::Float64=Inf,
    gap_storage_limit::Int=100,
)
    row_count = length(switch_types)
    _over16_check_length("positions", positions, row_count)
    _over16_check_length("elapsed_open_times", elapsed_open_times, row_count)
    _over16_check_length("maximum_open_times", maximum_open_times, row_count)
    _over16_check_length("voltage_differences", voltage_differences, row_count)
    _over16_check_length("switch_currents", switch_currents, row_count)
    _over16_check_length("close_thresholds", close_thresholds, row_count)
    _over16_check_length("open_thresholds", open_thresholds, row_count)
    _over16_check_optional_length("controlled_nodes", controlled_nodes, row_count)
    _over16_check_optional_length("control_indices", control_indices, row_count)
    _over16_check_optional_length("clamp_signals", clamp_signals, row_count)
    _over16_check_optional_length("close_control_signals", close_control_signals, row_count)
    _over16_check_optional_length("gap_control_signals", gap_control_signals, row_count)
    gap_storage_limit >= 0 || throw(ArgumentError("gap_storage_limit must be nonnegative"))

    updated_positions = collect(positions)
    updated_elapsed = collect(elapsed_open_times)
    updated_gap_currents = collect(previous_gap_currents)
    modswt = Int[]
    scan_actions = Vector{Symbol}(undef, row_count)
    transitions = Vector{Symbol}(undef, row_count)
    gap_count = 0

    for row in 1:row_count
        switch_type = switch_types[row]
        if switch_type == 8890
            gap_count += 1
            gap_count <= gap_storage_limit ||
                throw(ArgumentError("OVER16 controlled-switch gap storage overflow"))
            if length(updated_gap_currents) < gap_count
                push!(updated_gap_currents, 0.0)
            end
        end

        if !(switch_type == 8888 || switch_type == 8890 || switch_type == 8891)
            scan_actions[row] = :skip
            transitions[row] = :none
            continue
        end

        previous_gap_current = switch_type == 8890 ? updated_gap_currents[gap_count] : 0.0
        step = over16_controlled_switch_scan_step(
            switch_type,
            updated_positions[row],
            updated_elapsed[row],
            maximum_open_times[row],
            voltage_differences[row],
            switch_currents[row],
            close_thresholds[row],
            open_thresholds[row],
            dt;
            controlled_node = _over16_optional_int(controlled_nodes, row, 0),
            control_index = _over16_optional_int(control_indices, row, 0),
            clamp_signal = _over16_optional_float(clamp_signals, row, 0.0),
            close_control_signal = _over16_optional_float(close_control_signals, row, 0.0),
            gap_control_signal = _over16_optional_float(gap_control_signals, row, 0.0),
            previous_gap_current = previous_gap_current,
            flzero = flzero,
            fltinf = fltinf,
        )

        scan_actions[row] = step.scan_action
        transitions[row] = step.transition
        updated_positions[row] = step.position
        updated_elapsed[row] = step.elapsed_open_time
        if switch_type == 8890
            updated_gap_currents[gap_count] = step.gap_current
        end
        if step.altered
            push!(modswt, step.modswt_sign * row)
        end
    end

    return (
        positions = updated_positions,
        elapsed_open_times = updated_elapsed,
        gap_currents = updated_gap_currents,
        modswt = modswt,
        ktrlsw_count = length(modswt),
        altered = !isempty(modswt),
        gap_count = gap_count,
        scan_actions = scan_actions,
        transitions = transitions,
    )
end

function over16_controlled_switch_table_scan!(
    state::OVER16SwitchScanState,
    switch_types::AbstractVector{Int},
    maximum_open_times::AbstractVector{Float64},
    voltage_differences::AbstractVector{Float64},
    switch_currents::AbstractVector{Float64},
    close_thresholds::AbstractVector{Float64},
    open_thresholds::AbstractVector{Float64},
    dt::Float64;
    controlled_nodes::AbstractVector{Int}=Int[],
    control_indices::AbstractVector{Int}=Int[],
    clamp_signals::AbstractVector{Float64}=Float64[],
    close_control_signals::AbstractVector{Float64}=Float64[],
    gap_control_signals::AbstractVector{Float64}=Float64[],
    flzero::Float64=0.0,
    fltinf::Float64=Inf,
    gap_storage_limit::Int=100,
)
    positions_before = copy(state.positions)
    elapsed_before = copy(state.elapsed_open_times)
    gap_currents_before = copy(state.gap_currents)
    modswt_before = copy(state.modswt)
    preview = over16_controlled_switch_table_scan(
        switch_types,
        state.positions,
        state.elapsed_open_times,
        maximum_open_times,
        voltage_differences,
        switch_currents,
        close_thresholds,
        open_thresholds,
        dt;
        controlled_nodes = controlled_nodes,
        control_indices = control_indices,
        clamp_signals = clamp_signals,
        close_control_signals = close_control_signals,
        gap_control_signals = gap_control_signals,
        previous_gap_currents = state.gap_currents,
        flzero = flzero,
        fltinf = fltinf,
        gap_storage_limit = gap_storage_limit,
    )

    state.positions .= preview.positions
    state.elapsed_open_times .= preview.elapsed_open_times
    resize!(state.gap_currents, length(preview.gap_currents))
    state.gap_currents .= preview.gap_currents
    empty!(state.modswt)
    append!(state.modswt, preview.modswt)

    positions_mutated = state.positions != positions_before
    elapsed_open_times_mutated = state.elapsed_open_times != elapsed_before
    gap_currents_mutated = state.gap_currents != gap_currents_before
    modswt_mutated = state.modswt != modswt_before
    switch_scan_state_mutated =
        positions_mutated || elapsed_open_times_mutated ||
        gap_currents_mutated || modswt_mutated
    return merge(
        preview,
        (
            positions = copy(state.positions),
            elapsed_open_times = copy(state.elapsed_open_times),
            gap_currents = copy(state.gap_currents),
            modswt = copy(state.modswt),
            positions_mutated = positions_mutated,
            elapsed_open_times_mutated = elapsed_open_times_mutated,
            gap_currents_mutated = gap_currents_mutated,
            modswt_mutated = modswt_mutated,
            switch_scan_state_mutated = switch_scan_state_mutated,
            switch_graph_state_mutated = switch_scan_state_mutated,
            topology_mutated = false,
            admittance_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function over16_switch_operation_schedule(
    modswt::AbstractVector{Int},
    closed_switch_count::Int,
    accumulated_operation_count::Int=0,
)
    closed_switch_count >= 0 ||
        throw(ArgumentError("closed_switch_count must be nonnegative"))
    accumulated_operation_count >= 0 ||
        throw(ArgumentError("accumulated_operation_count must be nonnegative"))

    opening_rows = Int[]
    closing_rows = Int[]
    for entry in modswt
        entry != 0 || throw(ArgumentError("MODSWT entries must be nonzero"))
        if entry < 0
            push!(opening_rows, -entry)
        else
            push!(closing_rows, entry)
        end
    end

    operation_count = length(modswt)
    updated_closed_count = closed_switch_count - length(opening_rows) + length(closing_rows)
    updated_closed_count >= 0 ||
        throw(ArgumentError("switch operations produce a negative closed-switch count"))
    return (
        opening_rows = opening_rows,
        closing_rows = closing_rows,
        processed_modswt = vcat(-opening_rows, closing_rows),
        ktrlsw_count = operation_count,
        closed_switch_count = updated_closed_count,
        accumulated_operation_count = accumulated_operation_count + operation_count,
        should_call_switch = operation_count > 0,
    )
end

function over16_switch_operation_schedule!(state::OVER16SwitchOperationState)
    modswt_before = copy(state.modswt)
    closed_count_before = state.closed_switch_count
    accumulated_count_before = state.accumulated_operation_count
    preview = over16_switch_operation_schedule(
        state.modswt,
        state.closed_switch_count,
        state.accumulated_operation_count,
    )

    empty!(state.modswt)
    append!(state.modswt, preview.processed_modswt)
    state.closed_switch_count = preview.closed_switch_count
    state.accumulated_operation_count = preview.accumulated_operation_count

    modswt_mutated = state.modswt != modswt_before
    closed_switch_count_mutated = state.closed_switch_count != closed_count_before
    accumulated_operation_count_mutated =
        state.accumulated_operation_count != accumulated_count_before
    switch_operation_state_mutated =
        modswt_mutated || closed_switch_count_mutated ||
        accumulated_operation_count_mutated
    return merge(
        preview,
        (
            processed_modswt = copy(state.modswt),
            closed_switch_count = state.closed_switch_count,
            accumulated_operation_count = state.accumulated_operation_count,
            modswt_mutated = modswt_mutated,
            closed_switch_count_mutated = closed_switch_count_mutated,
            accumulated_operation_count_mutated = accumulated_operation_count_mutated,
            switch_operation_state_mutated = switch_operation_state_mutated,
            switch_graph_state_mutated = switch_operation_state_mutated,
            topology_mutated = false,
            admittance_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function switch_operation_schedule_lean!(state::OVER16SwitchOperationState)
    state.closed_switch_count >= 0 ||
        throw(ArgumentError("closed_switch_count must be nonnegative"))
    state.accumulated_operation_count >= 0 ||
        throw(ArgumentError("accumulated_operation_count must be nonnegative"))

    opening_count = 0
    closing_count = 0
    @inbounds for entry in state.modswt
        entry != 0 || throw(ArgumentError("MODSWT entries must be nonzero"))
        if entry < 0
            opening_count += 1
        else
            closing_count += 1
        end
    end
    operation_count = opening_count + closing_count
    updated_closed_count =
        state.closed_switch_count - opening_count + closing_count
    updated_closed_count >= 0 ||
        throw(ArgumentError("switch operations produce a negative closed-switch count"))

    # Stable in-place partition: the reference schedule processes openings before
    # closings while retaining the original order within each group.
    next_opening = firstindex(state.modswt)
    @inbounds for source in eachindex(state.modswt)
        if state.modswt[source] < 0
            entry = state.modswt[source]
            for destination in source:-1:(next_opening + 1)
                state.modswt[destination] = state.modswt[destination - 1]
            end
            state.modswt[next_opening] = entry
            next_opening += 1
        end
    end

    state.closed_switch_count = updated_closed_count
    state.accumulated_operation_count += operation_count
    return SwitchOperationStepResult(
        state.modswt,
        operation_count,
        operation_count > 0,
    )
end

function over16_switch_status_update(
    modswt::AbstractVector{Int},
    closed_mask::AbstractVector{Bool},
    closed_switch_count::Int,
    first_group_head::Int=0;
    strict_consistency::Bool=true,
)
    switch_count = length(closed_mask)
    closed_switch_count >= 0 ||
        throw(ArgumentError("closed_switch_count must be nonnegative"))
    closed_switch_count <= switch_count ||
        throw(ArgumentError("closed_switch_count cannot exceed switch count"))
    0 <= first_group_head <= switch_count ||
        throw(ArgumentError("first_group_head must be between 0 and switch count"))
    if strict_consistency && count(identity, closed_mask) != closed_switch_count
        throw(ArgumentError("closed_switch_count must match closed_mask in strict mode"))
    end

    updated_closed_mask = collect(closed_mask)
    updated_closed_count = closed_switch_count
    for entry in modswt
        entry != 0 || throw(ArgumentError("MODSWT entries must be nonzero"))
        row = abs(entry)
        1 <= row <= switch_count ||
            throw(ArgumentError("MODSWT row index out of switch range"))
        if entry < 0
            if strict_consistency && !updated_closed_mask[row]
                throw(ArgumentError("cannot open an already-open switch in strict mode"))
            end
            updated_closed_count -= 1
            updated_closed_mask[row] = false
        else
            if strict_consistency && updated_closed_mask[row]
                throw(ArgumentError("cannot close an already-closed switch in strict mode"))
            end
            updated_closed_count += 1
            updated_closed_mask[row] = true
        end
    end
    updated_closed_count >= 0 ||
        throw(ArgumentError("switch operations produce a negative closed-switch count"))
    updated_closed_count <= switch_count ||
        throw(ArgumentError("switch operations produce a closed-switch count above switch count"))

    return (
        closed_mask = updated_closed_mask,
        closed_switch_count = updated_closed_count,
        first_group_head = updated_closed_count > 0 ? first_group_head : 0,
        clear_nextsw = updated_closed_count == 0,
        clear_kode = updated_closed_count == 0,
        requires_order_rebuild = !isempty(modswt) && updated_closed_count > 0,
        topology_mutated = false,
    )
end

function over16_switch_status_update!(
    state::OVER16SwitchTopologyState,
    modswt::AbstractVector{Int};
    strict_consistency::Bool=true,
)
    closed_mask_before = copy(state.closed_mask)
    closed_count_before = state.closed_switch_count
    first_group_head_before = state.first_group_head
    nextsw_before = copy(state.nextsw)
    kode_before = copy(state.kode)
    preview = over16_switch_status_update(
        modswt,
        state.closed_mask,
        state.closed_switch_count,
        state.first_group_head;
        strict_consistency = strict_consistency,
    )

    state.closed_mask .= preview.closed_mask
    state.closed_switch_count = preview.closed_switch_count
    state.first_group_head = preview.first_group_head
    if preview.clear_nextsw
        fill!(state.nextsw, 0)
    end
    if preview.clear_kode
        fill!(state.kode, 0)
    end

    closed_mask_mutated = state.closed_mask != closed_mask_before
    closed_switch_count_mutated = state.closed_switch_count != closed_count_before
    first_group_head_mutated = state.first_group_head != first_group_head_before
    nextsw_cleared = state.nextsw != nextsw_before
    kode_cleared = state.kode != kode_before
    switch_status_state_mutated =
        closed_mask_mutated || closed_switch_count_mutated ||
        first_group_head_mutated || nextsw_cleared || kode_cleared
    return merge(
        preview,
        (
            closed_mask = copy(state.closed_mask),
            closed_switch_count = state.closed_switch_count,
            first_group_head = state.first_group_head,
            nextsw = copy(state.nextsw),
            kode = copy(state.kode),
            closed_mask_mutated = closed_mask_mutated,
            closed_switch_count_mutated = closed_switch_count_mutated,
            first_group_head_mutated = first_group_head_mutated,
            nextsw_cleared = nextsw_cleared,
            kode_cleared = kode_cleared,
            switch_status_state_mutated = switch_status_state_mutated,
            switch_graph_state_mutated = switch_status_state_mutated,
            topology_mutated = false,
            admittance_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function switch_status_update_lean!(
    state::OVER16SwitchTopologyState,
    modswt::AbstractVector{Int};
    strict_consistency::Bool=true,
)
    switch_count = length(state.closed_mask)
    state.closed_switch_count >= 0 ||
        throw(ArgumentError("closed_switch_count must be nonnegative"))
    state.closed_switch_count <= switch_count ||
        throw(ArgumentError("closed_switch_count cannot exceed switch count"))
    0 <= state.first_group_head <= switch_count ||
        throw(ArgumentError("first_group_head must be between 0 and switch count"))
    if strict_consistency && count(identity, state.closed_mask) != state.closed_switch_count
        throw(ArgumentError("closed_switch_count must match closed_mask in strict mode"))
    end

    updated_closed_count = state.closed_switch_count
    @inbounds for operation_index in eachindex(modswt)
        entry = modswt[operation_index]
        entry != 0 || throw(ArgumentError("MODSWT entries must be nonzero"))
        row = abs(entry)
        1 <= row <= switch_count ||
            throw(ArgumentError("MODSWT row index out of switch range"))
        row_closed = state.closed_mask[row]
        for previous_index in firstindex(modswt):(operation_index - 1)
            previous_entry = modswt[previous_index]
            abs(previous_entry) == row && (row_closed = previous_entry > 0)
        end
        if entry < 0
            strict_consistency && !row_closed &&
                throw(ArgumentError("cannot open an already-open switch in strict mode"))
            updated_closed_count -= 1
        else
            strict_consistency && row_closed &&
                throw(ArgumentError("cannot close an already-closed switch in strict mode"))
            updated_closed_count += 1
        end
    end
    updated_closed_count >= 0 ||
        throw(ArgumentError("switch operations produce a negative closed-switch count"))
    updated_closed_count <= switch_count ||
        throw(ArgumentError("switch operations produce a closed-switch count above switch count"))

    closed_mask_mutated = false
    @inbounds for row in eachindex(state.closed_mask)
        updated = state.closed_mask[row]
        for entry in modswt
            abs(entry) == row && (updated = entry > 0)
        end
        closed_mask_mutated |= updated != state.closed_mask[row]
        state.closed_mask[row] = updated
    end
    closed_switch_count_mutated = updated_closed_count != state.closed_switch_count
    state.closed_switch_count = updated_closed_count

    first_group_head_mutated = false
    nextsw_cleared = false
    kode_cleared = false
    if updated_closed_count == 0
        first_group_head_mutated = state.first_group_head != 0
        state.first_group_head = 0
        @inbounds for index in eachindex(state.nextsw)
            nextsw_cleared |= state.nextsw[index] != 0
            state.nextsw[index] = 0
        end
        @inbounds for index in eachindex(state.kode)
            kode_cleared |= state.kode[index] != 0
            state.kode[index] = 0
        end
    end
    mutated =
        closed_mask_mutated || closed_switch_count_mutated ||
        first_group_head_mutated || nextsw_cleared || kode_cleared
    return SwitchStatusStepResult(
        !isempty(modswt) && updated_closed_count > 0,
        mutated,
    )
end

function _over16_endpoint_is_known(node::Int, partition_boundary::Int, reference_node::Int)
    return node == reference_node || node > partition_boundary
end

function _over16_shares_endpoint(
    a_from::Int,
    a_to::Int,
    b_from::Int,
    b_to::Int,
)
    return a_from == b_from || a_from == b_to || a_to == b_from || a_to == b_to
end

function _over16_insert_kode_chain!(
    kode::Vector{Int},
    from_node::Int,
    to_node::Int,
    node_count::Int,
)
    if kode[from_node] == 0 && kode[to_node] == 0
        kode[from_node] = to_node
        kode[to_node] = from_node
        return nothing
    end

    anchor = from_node
    if kode[from_node] == 0
        anchor = to_node
    end
    insert_node = anchor == from_node ? to_node : from_node
    first_insert = insert_node

    steps = 0
    while kode[anchor] >= anchor
        anchor = kode[anchor]
        steps += 1
        steps <= node_count ||
            throw(ArgumentError("KODE chain search did not terminate"))
    end

    while true
        next_insert = kode[insert_node]
        if !(insert_node <= anchor && insert_node >= kode[anchor])
            old_next = kode[anchor]
            kode[anchor] = insert_node
            kode[insert_node] = old_next
            if insert_node > anchor
                anchor = insert_node
            end
        else
            predecessor = kode[anchor]
            steps = 0
            while !(insert_node < kode[predecessor])
                predecessor = kode[predecessor]
                steps += 1
                steps <= node_count ||
                    throw(ArgumentError("KODE insertion search did not terminate"))
            end
            old_next = kode[predecessor]
            kode[predecessor] = insert_node
            kode[insert_node] = old_next
        end

        if next_insert == 0 || next_insert == first_insert
            return nothing
        end
        insert_node = next_insert
    end
end

function over16_switch_simple_ordering(
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    closed_mask::AbstractVector{Bool},
    partition_boundary::Int;
    expected_closed_switch_count::Union{Nothing,Int}=nothing,
    node_count::Int=max(partition_boundary, maximum(vcat([1], collect(from_nodes), collect(to_nodes)))),
    reference_node::Int=1,
)
    switch_count = length(closed_mask)
    _over16_check_length("from_nodes", from_nodes, switch_count)
    _over16_check_length("to_nodes", to_nodes, switch_count)
    partition_boundary >= 1 ||
        throw(ArgumentError("partition_boundary must be positive"))
    reference_node >= 1 ||
        throw(ArgumentError("reference_node must be positive"))
    node_count >= partition_boundary ||
        throw(ArgumentError("node_count must be at least partition_boundary"))

    for row in 1:switch_count
        from_node = from_nodes[row]
        to_node = to_nodes[row]
        1 <= from_node <= node_count ||
            throw(ArgumentError("from_nodes entries must be within node_count"))
        1 <= to_node <= node_count ||
            throw(ArgumentError("to_nodes entries must be within node_count"))
        from_node != to_node ||
            throw(ArgumentError("closed switch ordering requires distinct endpoints"))
    end

    closed_switch_count = count(identity, closed_mask)
    if expected_closed_switch_count !== nothing &&
       closed_switch_count != expected_closed_switch_count
        throw(ArgumentError("expected_closed_switch_count must match closed_mask"))
    end
    nextsw = zeros(Int, switch_count)
    kode = zeros(Int, node_count)
    kcl_endpoint_indices = zeros(Int, switch_count)
    if closed_switch_count == 0
        return (
            nextsw = nextsw,
            kode = kode,
            first_group_head = 0,
            closed_switch_count = 0,
            ordered_switch_count = 0,
            ordered_rows = Int[],
            kcl_endpoint_indices = kcl_endpoint_indices,
            pass_count = 0,
            topology_mutated = false,
        )
    end

    ordered_switch_count = 0
    previous_count = 0
    previous_switch = 0
    first_switch = 0
    pass_count = 0
    ordered_rows = Int[]

    while ordered_switch_count < closed_switch_count
        for row in 1:switch_count
            if !closed_mask[row] || nextsw[row] != 0
                continue
            end

            from_node = from_nodes[row]
            to_node = to_nodes[row]
            from_blocked =
                _over16_endpoint_is_known(from_node, partition_boundary, reference_node)
            to_blocked =
                _over16_endpoint_is_known(to_node, partition_boundary, reference_node)
            blocked_count = (from_blocked ? 1 : 0) + (to_blocked ? 1 : 0)

            if blocked_count < 2
                for neighbor in 1:switch_count
                    if neighbor == row || !closed_mask[neighbor] || nextsw[neighbor] != 0
                        continue
                    end
                    if !_over16_shares_endpoint(
                        from_node,
                        to_node,
                        from_nodes[neighbor],
                        to_nodes[neighbor],
                    )
                        continue
                    end

                    if !from_blocked &&
                       (from_node == from_nodes[neighbor] || from_node == to_nodes[neighbor])
                        from_blocked = true
                        blocked_count += 1
                    elseif !to_blocked &&
                           (to_node == from_nodes[neighbor] || to_node == to_nodes[neighbor])
                        to_blocked = true
                        blocked_count += 1
                    end

                    if blocked_count == 2
                        break
                    end
                end
            end
            if blocked_count == 2
                continue
            end

            if ordered_switch_count > 0
                nextsw[previous_switch] *= row
            end
            previous_switch = row
            nextsw[row] = from_blocked ? -1 : 1
            kcl_endpoint_indices[row] = from_blocked ? 2 : 1
            ordered_switch_count += 1
            push!(ordered_rows, row)
            if ordered_switch_count == 1
                first_switch = row
            end
            _over16_insert_kode_chain!(kode, from_node, to_node, node_count)
        end

        pass_count += 1
        if ordered_switch_count <= previous_count
            throw(ArgumentError("closed switch ordering did not progress"))
        end
        previous_count = ordered_switch_count
    end

    nextsw[previous_switch] *= first_switch
    return (
        nextsw = nextsw,
        kode = kode,
        first_group_head = first_switch,
        closed_switch_count = closed_switch_count,
        ordered_switch_count = ordered_switch_count,
        ordered_rows = ordered_rows,
        kcl_endpoint_indices = kcl_endpoint_indices,
        pass_count = pass_count,
        topology_mutated = false,
    )
end

function over16_switch_simple_ordering!(
    state::OVER16SwitchTopologyState,
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    partition_boundary::Int;
    node_count::Int=length(state.kode),
    reference_node::Int=1,
)
    node_count > 0 || throw(ArgumentError("node_count must be positive"))
    nextsw_before = copy(state.nextsw)
    kode_before = copy(state.kode)
    first_group_head_before = state.first_group_head
    closed_count_before = state.closed_switch_count
    preview = over16_switch_simple_ordering(
        from_nodes,
        to_nodes,
        state.closed_mask,
        partition_boundary;
        expected_closed_switch_count = state.closed_switch_count,
        node_count = node_count,
        reference_node = reference_node,
    )

    resize!(state.nextsw, length(preview.nextsw))
    state.nextsw .= preview.nextsw
    resize!(state.kode, length(preview.kode))
    state.kode .= preview.kode
    state.first_group_head = preview.first_group_head
    state.closed_switch_count = preview.closed_switch_count

    nextsw_mutated = state.nextsw != nextsw_before
    kode_mutated = state.kode != kode_before
    first_group_head_mutated = state.first_group_head != first_group_head_before
    closed_switch_count_mutated = state.closed_switch_count != closed_count_before
    switch_order_state_mutated =
        nextsw_mutated || kode_mutated || first_group_head_mutated ||
        closed_switch_count_mutated
    return merge(
        preview,
        (
            nextsw = copy(state.nextsw),
            kode = copy(state.kode),
            first_group_head = state.first_group_head,
            closed_switch_count = state.closed_switch_count,
            nextsw_mutated = nextsw_mutated,
            kode_mutated = kode_mutated,
            first_group_head_mutated = first_group_head_mutated,
            closed_switch_count_mutated = closed_switch_count_mutated,
            switch_order_state_mutated = switch_order_state_mutated,
            switch_graph_state_mutated = switch_order_state_mutated,
            topology_mutated = false,
            admittance_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function switch_simple_ordering_lean!(
    state::OVER16SwitchTopologyState,
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    partition_boundary::Int;
    node_count::Int=length(state.kode),
    reference_node::Int=1,
)
    result = over16_switch_simple_ordering!(
        state,
        from_nodes,
        to_nodes,
        partition_boundary;
        node_count = node_count,
        reference_node = reference_node,
    )
    return SwitchOrderStepResult(result.switch_order_state_mutated)
end

function over16_switch_admittance_update(
    base_admittance::AbstractMatrix{<:Real},
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    closed_mask::AbstractVector{Bool};
    closed_conductances::AbstractVector{<:Real}=Float64[],
    open_conductances::AbstractVector{<:Real}=Float64[],
    closed_conductance::Real=1.0e9,
    open_conductance::Real=0.0,
    request_retriangularization::Bool=true,
)
    size(base_admittance, 1) == size(base_admittance, 2) ||
        throw(ArgumentError("base_admittance must be square"))
    node_count = size(base_admittance, 1)
    node_count > 0 || throw(ArgumentError("base_admittance must not be empty"))
    switch_count = length(closed_mask)
    _over16_check_length("from_nodes", from_nodes, switch_count)
    _over16_check_length("to_nodes", to_nodes, switch_count)
    closed_values = _over16_switch_conductance_vector(
        "closed_conductances",
        closed_conductances,
        switch_count,
        closed_conductance,
    )
    open_values = _over16_switch_conductance_vector(
        "open_conductances",
        open_conductances,
        switch_count,
        open_conductance,
    )

    admittance = Float64.(base_admittance)
    for value in admittance
        isfinite(value) ||
            throw(ArgumentError("base_admittance entries must be finite"))
    end
    switch_conductances = zeros(Float64, switch_count)
    closed_rows = Int[]
    open_rows = Int[]
    stamped_rows = Int[]
    for row in 1:switch_count
        from_node = from_nodes[row]
        to_node = to_nodes[row]
        0 <= from_node <= node_count ||
            throw(ArgumentError("from_nodes entries must be between 0 and node count"))
        0 <= to_node <= node_count ||
            throw(ArgumentError("to_nodes entries must be between 0 and node count"))
        from_node != to_node ||
            throw(ArgumentError("switch admittance endpoints must be distinct"))
        from_node != 0 || to_node != 0 ||
            throw(ArgumentError("switch admittance cannot connect ground to ground"))

        is_closed = closed_mask[row]
        conductance = is_closed ? closed_values[row] : open_values[row]
        switch_conductances[row] = conductance
        push!(is_closed ? closed_rows : open_rows, row)
        if conductance != 0.0
            stamp_conductance!(admittance, from_node, to_node, conductance)
            push!(stamped_rows, row)
        end
    end

    deferred_calls = Symbol[]
    if request_retriangularization
        push!(deferred_calls, :last14)
        push!(deferred_calls, :sparse_factor_update)
        push!(deferred_calls, :retriangularization_execution)
    end

    return (
        source = :over16_switch_admittance_update,
        outcome = :state_mutation,
        fortran_files = (:OVER16_FOR,),
        fortran_routines = (:SWITCH, :LAST14),
        fortran_labels = (312, 317, 601, 1010, 3111),
        node_count = node_count,
        switch_count = switch_count,
        closed_switch_count = length(closed_rows),
        closed_rows = closed_rows,
        open_rows = open_rows,
        stamped_rows = stamped_rows,
        base_admittance = Float64.(base_admittance),
        admittance = admittance,
        switch_conductances = switch_conductances,
        should_retriangularize = request_retriangularization,
        should_clear_factor_workspace = request_retriangularization,
        deferred_calls = deferred_calls,
        topology_mutated = false,
        admittance_mutated = false,
        sparse_factor_mutated = false,
        retriangularized = false,
        tacs_executed = false,
        solvum_executed = false,
    )
end

function over16_switch_admittance_update(
    node_count::Int,
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    closed_mask::AbstractVector{Bool};
    kwargs...,
)
    node_count > 0 || throw(ArgumentError("node_count must be positive"))
    return over16_switch_admittance_update(
        zeros(Float64, node_count, node_count),
        from_nodes,
        to_nodes,
        closed_mask;
        kwargs...,
    )
end

function over16_switch_admittance_update!(
    state::OVER16SwitchAdmittanceState,
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    closed_mask::AbstractVector{Bool};
    closed_conductances::AbstractVector{<:Real}=Float64[],
    open_conductances::AbstractVector{<:Real}=Float64[],
    closed_conductance::Real=1.0e9,
    open_conductance::Real=0.0,
    request_retriangularization::Bool=true,
)
    admittance_before = copy(state.admittance)
    switch_conductances_before = copy(state.switch_conductances)
    retriangularization_count_before = state.retriangularization_count
    preview = over16_switch_admittance_update(
        state.base_admittance,
        from_nodes,
        to_nodes,
        closed_mask;
        closed_conductances = closed_conductances,
        open_conductances = open_conductances,
        closed_conductance = closed_conductance,
        open_conductance = open_conductance,
        request_retriangularization = request_retriangularization,
    )

    state.admittance .= preview.admittance
    resize!(state.switch_conductances, length(preview.switch_conductances))
    state.switch_conductances .= preview.switch_conductances
    if preview.should_retriangularize
        state.retriangularization_count += 1
    end

    admittance_mutated = state.admittance != admittance_before
    switch_conductances_mutated =
        state.switch_conductances != switch_conductances_before
    retriangularization_count_mutated =
        state.retriangularization_count != retriangularization_count_before
    switch_admittance_state_mutated =
        admittance_mutated || switch_conductances_mutated ||
        retriangularization_count_mutated
    return merge(
        preview,
        (
            base_admittance = copy(state.base_admittance),
            admittance = copy(state.admittance),
            switch_conductances = copy(state.switch_conductances),
            retriangularization_count = state.retriangularization_count,
            admittance_mutated = admittance_mutated,
            switch_conductances_mutated = switch_conductances_mutated,
            retriangularization_count_mutated = retriangularization_count_mutated,
            switch_admittance_state_mutated = switch_admittance_state_mutated,
            topology_mutated = false,
            sparse_factor_mutated = false,
            retriangularized = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function switch_admittance_update_lean!(
    state::OVER16SwitchAdmittanceState,
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    closed_mask::AbstractVector{Bool};
    closed_conductances::AbstractVector{<:Real}=Float64[],
    open_conductances::AbstractVector{<:Real}=Float64[],
    closed_conductance::Real=1.0e9,
    open_conductance::Real=0.0,
    request_retriangularization::Bool=true,
)
    base_admittance = state.base_admittance
    size(base_admittance, 1) == size(base_admittance, 2) ||
        throw(ArgumentError("base_admittance must be square"))
    node_count = size(base_admittance, 1)
    node_count > 0 || throw(ArgumentError("base_admittance must not be empty"))
    size(state.admittance) == size(base_admittance) ||
        throw(ArgumentError("admittance and base_admittance sizes must match"))
    switch_count = length(closed_mask)
    _over16_check_length("from_nodes", from_nodes, switch_count)
    _over16_check_length("to_nodes", to_nodes, switch_count)
    _over16_check_optional_length(
        "closed_conductances",
        closed_conductances,
        switch_count,
    )
    _over16_check_optional_length(
        "open_conductances",
        open_conductances,
        switch_count,
    )
    closed_default = Float64(closed_conductance)
    open_default = Float64(open_conductance)
    isfinite(closed_default) && closed_default >= 0.0 ||
        throw(ArgumentError("closed_conductances default must be finite and nonnegative"))
    isfinite(open_default) && open_default >= 0.0 ||
        throw(ArgumentError("open_conductances default must be finite and nonnegative"))

    for value in base_admittance
        isfinite(value) ||
            throw(ArgumentError("base_admittance entries must be finite"))
    end
    @inbounds for row in 1:switch_count
        from_node = from_nodes[row]
        to_node = to_nodes[row]
        0 <= from_node <= node_count ||
            throw(ArgumentError("from_nodes entries must be between 0 and node count"))
        0 <= to_node <= node_count ||
            throw(ArgumentError("to_nodes entries must be between 0 and node count"))
        from_node != to_node ||
            throw(ArgumentError("switch admittance endpoints must be distinct"))
        from_node != 0 || to_node != 0 ||
            throw(ArgumentError("switch admittance cannot connect ground to ground"))
        closed_value = isempty(closed_conductances) ?
            closed_default : Float64(closed_conductances[row])
        open_value = isempty(open_conductances) ?
            open_default : Float64(open_conductances[row])
        isfinite(closed_value) && closed_value >= 0.0 ||
            throw(ArgumentError("closed_conductances entries must be finite and nonnegative"))
        isfinite(open_value) && open_value >= 0.0 ||
            throw(ArgumentError("open_conductances entries must be finite and nonnegative"))
    end

    old_conductance_count = length(state.switch_conductances)
    switch_conductances_mutated = old_conductance_count != switch_count
    resize!(state.switch_conductances, switch_count)
    @inbounds for row in 1:switch_count
        conductance = closed_mask[row] ?
            (isempty(closed_conductances) ? closed_default : Float64(closed_conductances[row])) :
            (isempty(open_conductances) ? open_default : Float64(open_conductances[row]))
        if row <= old_conductance_count
            switch_conductances_mutated |=
                state.switch_conductances[row] != conductance
        end
        state.switch_conductances[row] = conductance
    end

    if size(state.admittance_workspace) != size(base_admittance)
        state.admittance_workspace = similar(base_admittance)
    end
    copyto!(state.admittance_workspace, base_admittance)
    @inbounds for row in 1:switch_count
        conductance = state.switch_conductances[row]
        conductance == 0.0 && continue
        stamp_conductance!(
            state.admittance_workspace,
            from_nodes[row],
            to_nodes[row],
            conductance,
        )
    end
    admittance_mutated = state.admittance != state.admittance_workspace
    copyto!(state.admittance, state.admittance_workspace)

    if request_retriangularization
        state.retriangularization_count += 1
    end
    state_mutated =
        admittance_mutated || switch_conductances_mutated ||
        request_retriangularization
    return SwitchAdmittanceStepResult(
        request_retriangularization,
        admittance_mutated,
        state_mutated,
    )
end

function over16_switch_topology_admittance_update!(
    topology_state::OVER16SwitchTopologyState,
    admittance_state::OVER16SwitchAdmittanceState,
    modswt::AbstractVector{Int},
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    partition_boundary::Int;
    strict_consistency::Bool=true,
    node_count::Int=size(admittance_state.base_admittance, 1),
    reference_node::Int=1,
    closed_conductances::AbstractVector{<:Real}=Float64[],
    open_conductances::AbstractVector{<:Real}=Float64[],
    closed_conductance::Real=1.0e9,
    open_conductance::Real=0.0,
    request_retriangularization::Bool=!isempty(modswt),
)
    status_result = over16_switch_status_update!(
        topology_state,
        modswt;
        strict_consistency = strict_consistency,
    )
    should_rebuild_order =
        status_result.requires_order_rebuild ||
        (topology_state.closed_switch_count > 0 && all(entry -> entry == 0, topology_state.nextsw))
    order_result = nothing
    if should_rebuild_order
        order_result = over16_switch_simple_ordering!(
            topology_state,
            from_nodes,
            to_nodes,
            partition_boundary;
            node_count = node_count,
            reference_node = reference_node,
        )
    end
    admittance_result = over16_switch_admittance_update!(
        admittance_state,
        from_nodes,
        to_nodes,
        topology_state.closed_mask;
        closed_conductances = closed_conductances,
        open_conductances = open_conductances,
        closed_conductance = closed_conductance,
        open_conductance = open_conductance,
        request_retriangularization = request_retriangularization,
    )

    switch_topology_state_mutated =
        status_result.switch_status_state_mutated ||
        (order_result !== nothing && order_result.switch_order_state_mutated)
    switch_admittance_state_mutated =
        admittance_result.switch_admittance_state_mutated
    return (
        source = :over16_switch_topology_admittance_update,
        outcome = :integration_boundary,
        fortran_files = (:OVER16_FOR,),
        fortran_routines = (:SWITCH, :LAST14),
        fortran_labels = (312, 317, 601, 1010, 3111),
        status_result = status_result,
        order_result = order_result,
        admittance_result = admittance_result,
        closed_mask = copy(topology_state.closed_mask),
        closed_switch_count = topology_state.closed_switch_count,
        first_group_head = topology_state.first_group_head,
        nextsw = copy(topology_state.nextsw),
        kode = copy(topology_state.kode),
        admittance = copy(admittance_state.admittance),
        switch_conductances = copy(admittance_state.switch_conductances),
        retriangularization_count = admittance_state.retriangularization_count,
        switch_topology_state_mutated = switch_topology_state_mutated,
        switch_admittance_state_mutated = switch_admittance_state_mutated,
        switch_topology_admittance_state_mutated =
            switch_topology_state_mutated || switch_admittance_state_mutated,
        topology_mutated = switch_topology_state_mutated,
        admittance_mutated = admittance_result.admittance_mutated,
        sparse_factor_mutated = false,
        should_retriangularize = admittance_result.should_retriangularize,
        retriangularized = false,
        deferred_calls = admittance_result.deferred_calls,
        tacs_executed = false,
        solvum_executed = false,
    )
end

function switch_topology_admittance_update_lean!(
    topology_state::OVER16SwitchTopologyState,
    admittance_state::OVER16SwitchAdmittanceState,
    modswt::AbstractVector{Int},
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    partition_boundary::Int;
    strict_consistency::Bool=true,
    node_count::Int=size(admittance_state.base_admittance, 1),
    reference_node::Int=1,
    closed_conductances::AbstractVector{<:Real}=Float64[],
    open_conductances::AbstractVector{<:Real}=Float64[],
    closed_conductance::Real=1.0e9,
    open_conductance::Real=0.0,
    request_retriangularization::Bool=!isempty(modswt),
)
    status_result = switch_status_update_lean!(
        topology_state,
        modswt;
        strict_consistency = strict_consistency,
    )
    should_rebuild_order =
        status_result.requires_order_rebuild ||
        (topology_state.closed_switch_count > 0 &&
         all(iszero, topology_state.nextsw))
    order_result = should_rebuild_order ?
        switch_simple_ordering_lean!(
            topology_state,
            from_nodes,
            to_nodes,
            partition_boundary;
            node_count = node_count,
            reference_node = reference_node,
        ) : nothing
    admittance_result = switch_admittance_update_lean!(
        admittance_state,
        from_nodes,
        to_nodes,
        topology_state.closed_mask;
        closed_conductances = closed_conductances,
        open_conductances = open_conductances,
        closed_conductance = closed_conductance,
        open_conductance = open_conductance,
        request_retriangularization = request_retriangularization,
    )
    topology_mutated =
        status_result.switch_status_state_mutated ||
        (order_result !== nothing && order_result.switch_order_state_mutated)
    admittance_mutated = admittance_result.switch_admittance_state_mutated
    return SwitchTopologyAdmittanceStepResult(
        status_result,
        order_result,
        admittance_result,
        topology_mutated,
        admittance_mutated,
        topology_mutated || admittance_mutated,
        topology_mutated,
        admittance_result.admittance_mutated,
        admittance_result.should_retriangularize,
    )
end

function over16_switch_retriangularization_update(
    admittance::AbstractMatrix{<:Real};
    pivot_tolerance::Real=0.0,
)
    size(admittance, 1) == size(admittance, 2) ||
        throw(ArgumentError("admittance must be square"))
    node_count = size(admittance, 1)
    node_count > 0 || throw(ArgumentError("admittance must not be empty"))
    tolerance = Float64(pivot_tolerance)
    isfinite(tolerance) && tolerance >= 0.0 ||
        throw(ArgumentError("pivot_tolerance must be finite and nonnegative"))

    factor = Float64.(admittance)
    for value in factor
        isfinite(value) ||
            throw(ArgumentError("admittance entries must be finite"))
    end
    pivot_values = zeros(Float64, node_count)

    for k in 1:node_count
        pivot = factor[k, k]
        pivot_values[k] = pivot
        abs(pivot) > tolerance ||
            throw(ArgumentError("switch retriangularization pivot is zero or below tolerance"))
        if k < node_count
            for i in (k + 1):node_count
                multiplier = factor[i, k] / pivot
                factor[i, k] = multiplier
                for j in (k + 1):node_count
                    factor[i, j] -= multiplier * factor[k, j]
                end
            end
        end
    end

    return (
        source = :over16_switch_retriangularization_update,
        outcome = :integration_boundary,
        fortran_files = (:OVER16_FOR,),
        fortran_routines = (:LAST14, :SUBTS1),
        fortran_labels = (929, 930, 934, 935, 936, 940, 946, 963, 996, 1069, 1071, 1072, 1076, 1077, 1082, 1096),
        node_count = node_count,
        factor = factor,
        pivot_values = pivot_values,
        pivot_tolerance = tolerance,
        should_retriangularize = true,
        dense_factor_mutated = false,
        sparse_factor_mutated = false,
        fortran_sparse_factor_mutated = false,
        retriangularized = true,
        fortran_sparse_retriangularized = false,
        deferred_calls = [:last14, :fortran_sparse_factor_workspace, :nonlinear_inverse_columns],
        topology_mutated = false,
        admittance_mutated = false,
        tacs_executed = false,
        solvum_executed = false,
    )
end

function over16_switch_retriangularization_update!(
    state::OVER16SwitchRetriangularizationState,
    admittance::AbstractMatrix{<:Real};
    pivot_tolerance::Real=0.0,
)
    factor_before = copy(state.factor)
    pivots_before = copy(state.pivot_values)
    factorization_count_before = state.factorization_count
    preview = over16_switch_retriangularization_update(
        admittance;
        pivot_tolerance = pivot_tolerance,
    )

    if size(state.factor) == size(preview.factor)
        state.factor .= preview.factor
    else
        state.factor = copy(preview.factor)
    end
    resize!(state.pivot_values, length(preview.pivot_values))
    state.pivot_values .= preview.pivot_values
    state.factorization_count += 1

    dense_factor_mutated = state.factor != factor_before
    pivot_values_mutated = state.pivot_values != pivots_before
    factorization_count_mutated =
        state.factorization_count != factorization_count_before
    switch_retriangularization_state_mutated =
        dense_factor_mutated || pivot_values_mutated || factorization_count_mutated
    return merge(
        preview,
        (
            factor = copy(state.factor),
            pivot_values = copy(state.pivot_values),
            factorization_count = state.factorization_count,
            dense_factor_mutated = dense_factor_mutated,
            pivot_values_mutated = pivot_values_mutated,
            factorization_count_mutated = factorization_count_mutated,
            switch_retriangularization_state_mutated =
                switch_retriangularization_state_mutated,
            sparse_factor_mutated = false,
            fortran_sparse_factor_mutated = false,
            retriangularized = true,
            fortran_sparse_retriangularized = false,
            topology_mutated = false,
            admittance_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function over16_switch_retriangularization_update!(
    state::OVER16SwitchRetriangularizationState,
    admittance_state::OVER16SwitchAdmittanceState;
    kwargs...,
)
    return over16_switch_retriangularization_update!(
        state,
        admittance_state.admittance;
        kwargs...,
    )
end

function over16_switch_retriangularization_solve(
    factor::AbstractMatrix{<:Real},
    rhs::AbstractVector{<:Real};
    pivot_tolerance::Real=0.0,
)
    size(factor, 1) == size(factor, 2) ||
        throw(ArgumentError("factor must be square"))
    node_count = size(factor, 1)
    length(rhs) == node_count ||
        throw(ArgumentError("rhs length must match factor size"))
    tolerance = Float64(pivot_tolerance)
    isfinite(tolerance) && tolerance >= 0.0 ||
        throw(ArgumentError("pivot_tolerance must be finite and nonnegative"))

    dense_factor = Float64.(factor)
    for value in dense_factor
        isfinite(value) ||
            throw(ArgumentError("factor entries must be finite"))
    end
    solution = Float64.(rhs)

    for i in 1:node_count
        for j in 1:(i - 1)
            solution[i] -= dense_factor[i, j] * solution[j]
        end
    end
    for i in node_count:-1:1
        total = solution[i]
        for j in (i + 1):node_count
            total -= dense_factor[i, j] * solution[j]
        end
        pivot = dense_factor[i, i]
        abs(pivot) > tolerance ||
            throw(ArgumentError("switch retriangularization factor pivot is zero or below tolerance"))
        solution[i] = total / pivot
    end
    return solution
end

function over16_switch_retriangularization_solve(
    state::OVER16SwitchRetriangularizationState,
    rhs::AbstractVector{<:Real};
    kwargs...,
)
    return over16_switch_retriangularization_solve(state.factor, rhs; kwargs...)
end

function over16_switch_sparse_factor_matrix(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int},
)
    length(km) == length(ykm) ||
        throw(ArgumentError("km and ykm lengths must match"))
    node_count = length(kk)
    node_count > 0 || throw(ArgumentError("kk must not be empty"))
    dense = zeros(Float64, node_count, node_count)
    previous_end = 0
    for row in 1:node_count
        row_end = kk[row]
        previous_end <= row_end <= length(km) ||
            throw(ArgumentError("kk entries must be nondecreasing row ends within km"))
        row_start = previous_end + 1
        row_start <= row_end ||
            throw(ArgumentError("each sparse factor row must include a diagonal marker"))
        diagonal_seen = false
        for index in row_start:row_end
            marker = km[index]
            marker != 0 ||
                throw(ArgumentError("km entries must not be zero"))
            column = marker < 0 ? -marker : marker
            1 <= column <= node_count ||
                throw(ArgumentError("km column entries must be within node count"))
            value = Float64(ykm[index])
            isfinite(value) ||
                throw(ArgumentError("ykm entries must be finite"))
            if marker < 0
                column == row ||
                    throw(ArgumentError("negative km entries must mark the row diagonal"))
                !diagonal_seen ||
                    throw(ArgumentError("sparse factor row has duplicate diagonal markers"))
                diagonal_seen = true
            end
            dense[row, column] = value
        end
        diagonal_seen ||
            throw(ArgumentError("sparse factor row is missing a diagonal marker"))
        previous_end = row_end
    end
    previous_end == length(km) ||
        throw(ArgumentError("kk must consume every km/ykm entry"))
    return dense
end

function over16_switch_sparse_factor_matrix(
    state::OVER16SwitchSparseFactorWorkspaceState,
)
    return over16_switch_sparse_factor_matrix(state.km, state.ykm, state.kk)
end

function over16_switch_sparse_factor_update(
    factor::AbstractMatrix{<:Real};
    zero_tolerance::Real=0.0,
)
    size(factor, 1) == size(factor, 2) ||
        throw(ArgumentError("factor must be square"))
    node_count = size(factor, 1)
    node_count > 0 || throw(ArgumentError("factor must not be empty"))
    tolerance = Float64(zero_tolerance)
    isfinite(tolerance) && tolerance >= 0.0 ||
        throw(ArgumentError("zero_tolerance must be finite and nonnegative"))

    km = Int[]
    ykm = Float64[]
    kk = zeros(Int, node_count)
    for row in 1:node_count
        diagonal = Float64(factor[row, row])
        isfinite(diagonal) ||
            throw(ArgumentError("factor entries must be finite"))
        push!(km, -row)
        push!(ykm, diagonal)
        for column in 1:node_count
            column == row && continue
            value = Float64(factor[row, column])
            isfinite(value) ||
                throw(ArgumentError("factor entries must be finite"))
            if abs(value) > tolerance
                push!(km, column)
                push!(ykm, value)
            end
        end
        kk[row] = length(km)
    end

    return (
        source = :over16_switch_sparse_factor_update,
        outcome = :state_mutation,
        fortran_files = (:OVER16_FOR,),
        fortran_routines = (:LAST14, :SUBTS1),
        fortran_labels = (929, 930, 934, 935, 936, 940, 946, 963, 996, 1069, 1071, 1072, 1076, 1077, 1082, 1096),
        node_count = node_count,
        km = km,
        ykm = ykm,
        kk = kk,
        zero_tolerance = tolerance,
        sparse_factor_workspace_built = true,
        sparse_factor_workspace_mutated = false,
        switch_sparse_factor_workspace_state_mutated = false,
        sparse_factor_mutated = false,
        fortran_sparse_factor_mutated = false,
        deferred_calls = [:fortran_sparse_factor_ordering, :nonlinear_inverse_columns],
        tacs_executed = false,
        solvum_executed = false,
    )
end

function over16_switch_sparse_factor_update(
    state::OVER16SwitchRetriangularizationState;
    kwargs...,
)
    return over16_switch_sparse_factor_update(state.factor; kwargs...)
end

function over16_switch_sparse_factor_update!(
    state::OVER16SwitchSparseFactorWorkspaceState,
    factor::AbstractMatrix{<:Real};
    zero_tolerance::Real=0.0,
)
    km_before = copy(state.km)
    ykm_before = copy(state.ykm)
    kk_before = copy(state.kk)
    count_before = state.workspace_update_count
    preview = over16_switch_sparse_factor_update(
        factor;
        zero_tolerance = zero_tolerance,
    )

    resize!(state.km, length(preview.km))
    state.km .= preview.km
    resize!(state.ykm, length(preview.ykm))
    state.ykm .= preview.ykm
    resize!(state.kk, length(preview.kk))
    state.kk .= preview.kk
    state.workspace_update_count += 1

    km_mutated = state.km != km_before
    ykm_mutated = state.ykm != ykm_before
    kk_mutated = state.kk != kk_before
    workspace_update_count_mutated =
        state.workspace_update_count != count_before
    sparse_factor_workspace_mutated =
        km_mutated || ykm_mutated || kk_mutated || workspace_update_count_mutated
    return merge(
        preview,
        (
            km = copy(state.km),
            ykm = copy(state.ykm),
            kk = copy(state.kk),
            workspace_update_count = state.workspace_update_count,
            km_mutated = km_mutated,
            ykm_mutated = ykm_mutated,
            kk_mutated = kk_mutated,
            workspace_update_count_mutated = workspace_update_count_mutated,
            sparse_factor_workspace_mutated = sparse_factor_workspace_mutated,
            switch_sparse_factor_workspace_state_mutated =
                sparse_factor_workspace_mutated,
            sparse_factor_mutated = sparse_factor_workspace_mutated,
            fortran_sparse_factor_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function over16_switch_sparse_factor_update!(
    state::OVER16SwitchSparseFactorWorkspaceState,
    retriangularization_state::OVER16SwitchRetriangularizationState;
    kwargs...,
)
    return over16_switch_sparse_factor_update!(
        state,
        retriangularization_state.factor;
        kwargs...,
    )
end

function over16_switch_sparse_factor_solve(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int},
    rhs::AbstractVector{<:Real};
    kwargs...,
)
    factor = over16_switch_sparse_factor_matrix(km, ykm, kk)
    return over16_switch_retriangularization_solve(factor, rhs; kwargs...)
end

function over16_switch_sparse_factor_solve(
    state::OVER16SwitchSparseFactorWorkspaceState,
    rhs::AbstractVector{<:Real};
    kwargs...,
)
    return over16_switch_sparse_factor_solve(state.km, state.ykm, state.kk, rhs; kwargs...)
end

const OVER16_SPARSE_FACTOR_LABELS = (
    2205, 2220, 2227, 2237, 2229, 2230, 2231, 2232, 2233, 2240,
    2260, 2265, 2270, 2278, 2272, 2273, 2275, 2276, 2277, 2279,
    2280, 4312, 2283, 2288, 2285, 2290,
)

function sparse_network_admittance_matrix(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    row_end_pointers::AbstractVector{Int};
    node_count::Int=length(row_end_pointers),
)
    length(km) == length(ykm) ||
        throw(ArgumentError("km and ykm lengths must match"))
    node_count == length(row_end_pointers) ||
        throw(ArgumentError("node_count must match row_end_pointers length"))
    node_count > 1 ||
        throw(ArgumentError("row_end_pointers must include reference and at least one sparse row"))
    previous_end = 0
    admittance = zeros(Float64, node_count, node_count)
    for row in 1:node_count
        one_past_row_end = row_end_pointers[row]
        1 <= one_past_row_end <= length(km) + 1 ||
            throw(ArgumentError("row_end_pointers entries must point one past a sparse row"))
        row_end = one_past_row_end - 1
        previous_end <= row_end ||
            throw(ArgumentError("row_end_pointers entries must be nondecreasing"))
        for index in (previous_end + 1):row_end
            column = abs(km[index])
            1 <= column <= node_count ||
                throw(ArgumentError("km column entries must be within node count"))
            value = Float64(ykm[index])
            isfinite(value) || throw(ArgumentError("ykm entries must be finite"))
            admittance[row, column] += value
        end
        previous_end = row_end
    end
    previous_end == length(km) ||
        throw(ArgumentError("row_end_pointers must consume every km/ykm entry"))
    return admittance
end

function _check_fortran_sparse_workspace(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int},
    first_factor_row::Int,
    partition_boundary::Int,
)
    length(km) == length(ykm) ||
        throw(ArgumentError("km and ykm lengths must match"))
    node_count = length(kk)
    node_count > 1 || throw(ArgumentError("kk must cover at least the reference and one factor row"))
    2 <= first_factor_row <= partition_boundary <= node_count ||
        throw(ArgumentError("invalid first_factor_row or partition_boundary"))
    all(index -> index == 0, kk[1:(first_factor_row - 1)]) ||
        throw(ArgumentError("kk entries before first_factor_row must be zero"))
    previous_end = 0
    for row in first_factor_row:partition_boundary
        row_end = kk[row]
        previous_end < row_end <= length(km) ||
            throw(ArgumentError("kk entries must provide nonempty increasing factor rows"))
        row_start = previous_end + 1
        km[row_start] == -row ||
            throw(ArgumentError("each Fortran sparse factor row must start with its negative diagonal marker"))
        for index in row_start:row_end
            marker = km[index]
            marker != 0 || throw(ArgumentError("km entries must not be zero"))
            column = abs(marker)
            1 <= column <= node_count ||
                throw(ArgumentError("km column entries must be within node count"))
            isfinite(Float64(ykm[index])) ||
                throw(ArgumentError("ykm entries must be finite"))
            if index > row_start
                marker > row ||
                    throw(ArgumentError("off-diagonal Fortran sparse factor entries must point to later columns"))
            end
        end
        previous_end = row_end
    end
    previous_end == length(km) ||
        throw(ArgumentError("kk must consume every km/ykm entry"))
    return nothing
end

function _node_group_representatives(
    node_group_successors::AbstractVector{<:Integer},
    node_count::Int,
)
    length(node_group_successors) >= node_count ||
        throw(ArgumentError("node_group_successors must cover every network node"))
    representatives = collect(1:node_count)
    for start in 1:node_count
        successor = Int(node_group_successors[start])
        successor == 0 && continue
        1 <= successor <= node_count ||
            throw(ArgumentError("node_group_successors entries must be valid node indices"))
        path = Int[]
        seen = Dict{Int,Int}()
        node = start
        while node != 0 && !haskey(seen, node)
            1 <= node <= node_count ||
                throw(ArgumentError("node group chain left the valid node range"))
            seen[node] = length(path) + 1
            push!(path, node)
            node = Int(node_group_successors[node])
        end
        node == 0 && continue
        group = path[seen[node]:end]
        representative = maximum(group)
        for grouped_node in group
            representatives[grouped_node] = representative
        end
    end
    return representatives
end

function _node_group_contracted_admittance(
    admittance::AbstractMatrix{<:Real},
    node_group_successors::AbstractVector{<:Integer},
)
    size(admittance, 1) == size(admittance, 2) ||
        throw(ArgumentError("admittance must be square"))
    node_count = size(admittance, 1)
    dense = Float64.(admittance)
    for value in dense
        isfinite(value) ||
            throw(ArgumentError("admittance entries must be finite"))
    end
    representatives = _node_group_representatives(node_group_successors, node_count)
    contracted = zeros(Float64, node_count, node_count)
    for column in 1:node_count
        representative_column = representatives[column]
        for row in 1:node_count
            contracted[representatives[row], representative_column] += dense[row, column]
        end
    end
    return contracted, representatives
end

function _grouped_sparse_network_factorization_update(
    admittance::AbstractMatrix{<:Real},
    node_group_successors::AbstractVector{<:Integer};
    partition_boundary::Int=size(admittance, 1),
    first_factor_row::Int=2,
    pivot_tolerance::Real=0.0,
    zero_tolerance::Real=0.0,
)
    size(admittance, 1) == size(admittance, 2) ||
        throw(ArgumentError("admittance must be square"))
    node_count = size(admittance, 1)
    node_count > 1 || throw(ArgumentError("admittance must include a reference node and factor rows"))
    2 <= first_factor_row <= partition_boundary <= node_count ||
        throw(ArgumentError("invalid first_factor_row or partition_boundary"))
    pivot_tol = Float64(pivot_tolerance)
    zero_tol = Float64(zero_tolerance)
    isfinite(pivot_tol) && pivot_tol >= 0.0 ||
        throw(ArgumentError("pivot_tolerance must be finite and nonnegative"))
    isfinite(zero_tol) && zero_tol >= 0.0 ||
        throw(ArgumentError("zero_tolerance must be finite and nonnegative"))

    dense, representatives =
        _node_group_contracted_admittance(admittance, node_group_successors)
    active_rows = Int[
        row for row in first_factor_row:partition_boundary
        if representatives[row] == row
    ]
    km = Int[]
    ykm = Float64[]
    kk = zeros(Int, node_count)
    row_starts = zeros(Int, node_count)
    pivot_values = Float64[]
    inverse_diagonal_values = Float64[]

    for row in active_rows
        f = zeros(Float64, node_count)
        f[1] = dense[row, row]
        for column in first_factor_row:node_count
            column == row && continue
            representatives[column] == column || continue
            value = dense[row, column]
            if abs(value) > zero_tol
                f[column] = value
            end
        end

        for previous_row in active_rows
            previous_row < row || break
            a = f[previous_row]
            abs(a) > zero_tol || continue
            previous_start = row_starts[previous_row]
            previous_end = kk[previous_row]
            previous_start <= previous_end ||
                throw(ArgumentError("invalid previously built grouped sparse factor row"))
            for index in (previous_start + 1):previous_end
                column = km[index]
                if column == row
                    f[1] -= a * ykm[index]
                else
                    f[column] -= a * ykm[index]
                end
            end
        end

        pivot = f[1]
        abs(pivot) > pivot_tol ||
            throw(ArgumentError("grouped sparse factor pivot is zero or below tolerance"))
        inverse_pivot = 1.0 / pivot
        push!(pivot_values, pivot)
        push!(inverse_diagonal_values, inverse_pivot)
        row_starts[row] = length(km) + 1
        push!(km, -row)
        push!(ykm, inverse_pivot)
        for column in (row + 1):node_count
            representatives[column] == column || continue
            value = f[column]
            if abs(value) > zero_tol
                push!(km, column)
                push!(ykm, value * inverse_pivot)
            end
        end
        kk[row] = length(km)
    end

    return (
        source = :grouped_sparse_network_factorization_update,
        outcome = :equation_oracle,
        fortran_files = (:OVER16_FOR,),
        fortran_routines = (:SUBTS1, :SUBTS3),
        fortran_labels = OVER16_SPARSE_FACTOR_LABELS,
        node_count = node_count,
        partition_boundary = partition_boundary,
        first_factor_row = first_factor_row,
        km = km,
        ykm = ykm,
        kk = kk,
        pivot_values = pivot_values,
        inverse_diagonal_values = inverse_diagonal_values,
        pivot_tolerance = pivot_tol,
        zero_tolerance = zero_tol,
        node_group_successors = Int.(node_group_successors[1:node_count]),
        node_group_representatives = representatives,
        sparse_factor_workspace_built = true,
        fortran_sparse_factor_mutated = false,
        deferred_calls = [:nonlinear_inverse_columns, :full_timestep_oracle],
        tacs_executed = false,
        solvum_executed = false,
    )
end

function over16_fortran_sparse_factor_update(
    admittance::AbstractMatrix{<:Real};
    partition_boundary::Int=size(admittance, 1),
    first_factor_row::Int=2,
    pivot_tolerance::Real=0.0,
    zero_tolerance::Real=0.0,
    node_group_successors::Union{Nothing,AbstractVector{<:Integer}}=nothing,
)
    if node_group_successors !== nothing
        return _grouped_sparse_network_factorization_update(
            admittance,
            node_group_successors;
            partition_boundary = partition_boundary,
            first_factor_row = first_factor_row,
            pivot_tolerance = pivot_tolerance,
            zero_tolerance = zero_tolerance,
        )
    end
    size(admittance, 1) == size(admittance, 2) ||
        throw(ArgumentError("admittance must be square"))
    node_count = size(admittance, 1)
    node_count > 1 || throw(ArgumentError("admittance must include a reference node and factor rows"))
    2 <= first_factor_row <= partition_boundary <= node_count ||
        throw(ArgumentError("invalid first_factor_row or partition_boundary"))
    pivot_tol = Float64(pivot_tolerance)
    zero_tol = Float64(zero_tolerance)
    isfinite(pivot_tol) && pivot_tol >= 0.0 ||
        throw(ArgumentError("pivot_tolerance must be finite and nonnegative"))
    isfinite(zero_tol) && zero_tol >= 0.0 ||
        throw(ArgumentError("zero_tolerance must be finite and nonnegative"))

    dense = Float64.(admittance)
    for value in dense
        isfinite(value) ||
            throw(ArgumentError("admittance entries must be finite"))
    end

    km = Int[]
    ykm = Float64[]
    kk = zeros(Int, node_count)
    pivot_values = Float64[]
    inverse_diagonal_values = Float64[]

    for row in first_factor_row:partition_boundary
        f = zeros(Float64, node_count)
        f[1] = dense[row, row]
        for column in first_factor_row:node_count
            column == row && continue
            value = dense[row, column]
            if abs(value) > zero_tol
                f[column] = value
            end
        end

        for previous_row in first_factor_row:(row - 1)
            a = f[previous_row]
            abs(a) > zero_tol || continue
            previous_start = previous_row == first_factor_row ? 1 : kk[previous_row - 1] + 1
            previous_end = kk[previous_row]
            previous_start <= previous_end ||
                throw(ArgumentError("invalid previously built Fortran sparse factor row"))
            for index in (previous_start + 1):previous_end
                column = km[index]
                if column == row
                    f[1] -= a * ykm[index]
                else
                    f[column] -= a * ykm[index]
                end
            end
        end

        pivot = f[1]
        abs(pivot) > pivot_tol ||
            throw(ArgumentError("Fortran sparse factor pivot is zero or below tolerance"))
        inverse_pivot = 1.0 / pivot
        push!(pivot_values, pivot)
        push!(inverse_diagonal_values, inverse_pivot)
        push!(km, -row)
        push!(ykm, inverse_pivot)
        for column in (row + 1):node_count
            value = f[column]
            if abs(value) > zero_tol
                push!(km, column)
                push!(ykm, value * inverse_pivot)
            end
        end
        kk[row] = length(km)
    end

    return (
        source = :over16_fortran_sparse_factor_update,
        outcome = :equation_oracle,
        fortran_files = (:OVER16_FOR,),
        fortran_routines = (:SUBTS1,),
        fortran_labels = OVER16_SPARSE_FACTOR_LABELS,
        node_count = node_count,
        partition_boundary = partition_boundary,
        first_factor_row = first_factor_row,
        km = km,
        ykm = ykm,
        kk = kk,
        pivot_values = pivot_values,
        inverse_diagonal_values = inverse_diagonal_values,
        pivot_tolerance = pivot_tol,
        zero_tolerance = zero_tol,
        sparse_factor_workspace_built = true,
        fortran_sparse_factor_mutated = false,
        deferred_calls = [:kode_group_ordering, :nonlinear_inverse_columns, :full_timestep_oracle],
        tacs_executed = false,
        solvum_executed = false,
    )
end

function over16_fortran_sparse_factor_update!(
    state::OVER16FortranSparseFactorWorkspaceState,
    admittance::AbstractMatrix{<:Real};
    kwargs...,
)
    km_before = copy(state.km)
    ykm_before = copy(state.ykm)
    kk_before = copy(state.kk)
    count_before = state.workspace_update_count
    preview = over16_fortran_sparse_factor_update(admittance; kwargs...)

    resize!(state.km, length(preview.km))
    state.km .= preview.km
    resize!(state.ykm, length(preview.ykm))
    state.ykm .= preview.ykm
    resize!(state.kk, length(preview.kk))
    state.kk .= preview.kk
    state.workspace_update_count += 1

    workspace_mutated =
        state.km != km_before ||
        state.ykm != ykm_before ||
        state.kk != kk_before ||
        state.workspace_update_count != count_before
    return merge(
        preview,
        (
            km = copy(state.km),
            ykm = copy(state.ykm),
            kk = copy(state.kk),
            workspace_update_count = state.workspace_update_count,
            fortran_sparse_factor_workspace_state_mutated = workspace_mutated,
            fortran_sparse_factor_mutated = workspace_mutated,
        ),
    )
end

function over16_fortran_sparse_factor_solve(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int},
    rhs::AbstractVector{<:Real};
    first_factor_row::Int=2,
    partition_boundary::Int=length(kk),
)
    _check_fortran_sparse_workspace(km, ykm, kk, first_factor_row, partition_boundary)
    node_count = length(kk)
    length(rhs) == node_count ||
        throw(ArgumentError("rhs length must match kk node count"))
    solution = Float64.(rhs)
    for value in solution
        isfinite(value) || throw(ArgumentError("rhs entries must be finite"))
    end
    solution[1] = 0.0

    for row in first_factor_row:partition_boundary
        row_start = row == first_factor_row ? 1 : kk[row - 1] + 1
        row_end = kk[row]
        a = solution[row]
        solution[row] = a * Float64(ykm[row_start])
        for index in (row_start + 1):row_end
            column = km[index]
            if column <= partition_boundary
                solution[column] -= a * Float64(ykm[index])
            end
        end
    end

    for row in partition_boundary:-1:first_factor_row
        row_start = row == first_factor_row ? 1 : kk[row - 1] + 1
        row_end = kk[row]
        correction = 0.0
        for index in row_end:-1:(row_start + 1)
            column = km[index]
            correction -= solution[column] * Float64(ykm[index])
        end
        solution[row] += correction
    end

    return solution
end

function over16_fortran_sparse_factor_solve(
    state::OVER16FortranSparseFactorWorkspaceState,
    rhs::AbstractVector{<:Real};
    kwargs...,
)
    return over16_fortran_sparse_factor_solve(state.km, state.ykm, state.kk, rhs; kwargs...)
end

function _check_node_group_successors(
    node_group_successors::AbstractVector{<:Integer},
    node_count::Int,
)
    length(node_group_successors) >= node_count ||
        throw(ArgumentError("node_group_successors must cover every network node"))
    for index in 1:node_count
        successor = node_group_successors[index]
        0 <= successor <= node_count ||
            throw(ArgumentError("node_group_successors entries must be zero or a valid node index"))
    end
    return nothing
end

function _node_group_representatives!(
    workspace::GroupedSparseNetworkWorkspace,
    node_group_successors::AbstractVector{<:Integer},
)
    node_count = length(workspace.node_group_representatives)
    _check_node_group_successors(node_group_successors, node_count)
    representatives = workspace.node_group_representatives
    for node in 1:node_count
        representatives[node] = node
    end

    marks = workspace.visit_marks
    positions = workspace.visit_positions
    path = workspace.visit_path
    for start in 1:node_count
        generation = start
        path_length = 0
        node = start
        while node != 0 && marks[node] != generation
            path_length += 1
            path[path_length] = node
            marks[node] = generation
            positions[node] = path_length
            node = Int(node_group_successors[node])
        end
        node == 0 && continue
        cycle_start = positions[node]
        representative = path[cycle_start]
        for path_index in (cycle_start + 1):path_length
            representative = max(representative, path[path_index])
        end
        for path_index in cycle_start:path_length
            representatives[path[path_index]] = representative
        end
    end
    return representatives
end

function grouped_sparse_network_factorization_update!(
    workspace::GroupedSparseNetworkWorkspace,
    admittance::AbstractMatrix{<:Real},
    node_group_successors::AbstractVector{<:Integer};
    partition_boundary::Int=size(admittance, 1),
    first_factor_row::Int=2,
    pivot_tolerance::Real=0.0,
    zero_tolerance::Real=0.0,
)
    size(admittance, 1) == size(admittance, 2) ||
        throw(ArgumentError("admittance must be square"))
    node_count = size(admittance, 1)
    node_count == length(workspace.kk) || throw(ArgumentError(
        "grouped sparse workspace size must match admittance",
    ))
    2 <= first_factor_row <= partition_boundary <= node_count ||
        throw(ArgumentError("invalid first_factor_row or partition_boundary"))
    pivot_tol = Float64(pivot_tolerance)
    zero_tol = Float64(zero_tolerance)
    isfinite(pivot_tol) && pivot_tol >= 0.0 ||
        throw(ArgumentError("pivot_tolerance must be finite and nonnegative"))
    isfinite(zero_tol) && zero_tol >= 0.0 ||
        throw(ArgumentError("zero_tolerance must be finite and nonnegative"))

    representatives =
        _node_group_representatives!(workspace, node_group_successors)
    contracted = workspace.contracted_admittance
    fill!(contracted, 0.0)
    for column in 1:node_count
        representative_column = representatives[column]
        for row in 1:node_count
            value = Float64(admittance[row, column])
            isfinite(value) ||
                throw(ArgumentError("admittance entries must be finite"))
            contracted[representatives[row], representative_column] += value
        end
    end

    active_rows = workspace.active_rows
    resize!(active_rows, 0)
    for row in first_factor_row:partition_boundary
        representatives[row] == row && push!(active_rows, row)
    end
    km = workspace.km
    ykm = workspace.ykm
    kk = workspace.kk
    row_starts = workspace.row_starts
    pivot_values = workspace.pivot_values
    inverse_diagonal_values = workspace.inverse_diagonal_values
    factor_row = workspace.factor_row
    resize!(km, 0)
    resize!(ykm, 0)
    resize!(pivot_values, 0)
    resize!(inverse_diagonal_values, 0)
    fill!(kk, 0)
    fill!(row_starts, 0)

    for row in active_rows
        fill!(factor_row, 0.0)
        factor_row[1] = contracted[row, row]
        for column in first_factor_row:node_count
            column == row && continue
            representatives[column] == column || continue
            value = contracted[row, column]
            abs(value) > zero_tol && (factor_row[column] = value)
        end

        for previous_row in active_rows
            previous_row < row || break
            a = factor_row[previous_row]
            abs(a) > zero_tol || continue
            previous_start = row_starts[previous_row]
            previous_end = kk[previous_row]
            previous_start <= previous_end || throw(ArgumentError(
                "invalid previously built grouped sparse factor row",
            ))
            for index in (previous_start + 1):previous_end
                column = km[index]
                if column == row
                    factor_row[1] -= a * ykm[index]
                else
                    factor_row[column] -= a * ykm[index]
                end
            end
        end

        pivot = factor_row[1]
        abs(pivot) > pivot_tol || throw(ArgumentError(
            "grouped sparse factor pivot is zero or below tolerance",
        ))
        inverse_pivot = 1.0 / pivot
        push!(pivot_values, pivot)
        push!(inverse_diagonal_values, inverse_pivot)
        row_starts[row] = length(km) + 1
        push!(km, -row)
        push!(ykm, inverse_pivot)
        for column in (row + 1):node_count
            representatives[column] == column || continue
            value = factor_row[column]
            if abs(value) > zero_tol
                push!(km, column)
                push!(ykm, value * inverse_pivot)
            end
        end
        kk[row] = length(km)
    end
    workspace.factorization_count += 1
    return workspace
end

function grouped_sparse_network_solution!(
    workspace::GroupedSparseNetworkWorkspace,
    rhs::AbstractVector{<:Real},
    node_group_successors::AbstractVector{<:Integer},
    initial_solution::AbstractVector{<:Real};
    first_factor_row::Int=2,
    partition_boundary::Int=length(workspace.kk),
)
    first_factor_row >= 1 ||
        throw(ArgumentError("first_factor_row must be positive"))
    node_count = length(workspace.solution)
    length(initial_solution) == node_count || throw(ArgumentError(
        "initial solution length must match grouped sparse workspace",
    ))
    length(rhs) >= partition_boundary ||
        throw(ArgumentError("rhs must cover the active solved partition"))
    1 <= partition_boundary <= node_count ||
        throw(ArgumentError("partition_boundary must be within network nodes"))
    _check_node_group_successors(node_group_successors, node_count)

    solution = workspace.solution
    for index in 1:node_count
        value = Float64(initial_solution[index])
        isfinite(value) || throw(ArgumentError("solution entries must be finite"))
        solution[index] = value
    end
    for index in 1:partition_boundary
        value = Float64(rhs[index])
        isfinite(value) || throw(ArgumentError("rhs entries must be finite"))
        solution[index] = value
    end
    alias_successors = workspace.alias_successors
    for index in 1:node_count
        alias_successors[index] = Int(node_group_successors[index])
    end

    active_node = 1
    visited = 0
    while true
        solution[active_node] = 0.0
        alias_successors[active_node] <= active_node && break
        active_node = alias_successors[active_node]
        visited += 1
        visited <= node_count || throw(ArgumentError(
            "node group reference-chain search did not terminate",
        ))
    end

    for node in 2:partition_boundary
        successor = alias_successors[node]
        successor == 0 && continue
        successor > partition_boundary && continue
        successor > node && (solution[successor] += solution[node])
    end

    km = workspace.km
    ykm = workspace.ykm
    kk = workspace.kk
    entry_index = 1
    entry_count = length(km)
    while entry_index <= entry_count
        row = abs(km[entry_index])
        scale = solution[row]
        solution[row] = scale * ykm[entry_index]
        row_end = kk[row]
        while true
            entry_index += 1
            entry_index > row_end && break
            column = km[entry_index]
            column > partition_boundary && continue
            solution[column] -= scale * ykm[entry_index]
        end
    end

    node = node_count
    while true
        successor = alias_successors[node]
        if successor != 0 && successor <= node
            copied_from = node
            copied_to = successor
            visited = 0
            while true
                solution[copied_to] = solution[copied_from]
                copied_from = copied_to
                copied_to = alias_successors[copied_from]
                copied_to == node && break
                visited += 1
                visited <= node_count || throw(ArgumentError(
                    "node group copy cycle did not terminate",
                ))
            end
        end
        if node <= partition_boundary
            entry_index == 1 && break
            correction = 0.0
            while true
                entry_index -= 1
                column = km[entry_index]
                if column < 0
                    node = abs(column)
                    solution[node] += correction
                    break
                end
                correction -= solution[column] * ykm[entry_index]
            end
        else
            node -= 1
        end
    end
    return solution
end

function _grouped_sparse_network_solution(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int},
    rhs::AbstractVector{<:Real},
    node_group_successors::AbstractVector{<:Integer},
    initial_solution::AbstractVector{<:Real};
    first_factor_row::Int=2,
    partition_boundary::Int=length(kk),
)
    first_factor_row >= 1 ||
        throw(ArgumentError("first_factor_row must be positive"))
    length(km) == length(ykm) ||
        throw(ArgumentError("km and ykm lengths must match"))
    node_count = length(initial_solution)
    length(kk) >= node_count ||
        throw(ArgumentError("kk must cover every network node"))
    length(rhs) >= partition_boundary ||
        throw(ArgumentError("rhs must cover the active solved partition"))
    1 <= partition_boundary <= node_count ||
        throw(ArgumentError("partition_boundary must be within network nodes"))
    _check_node_group_successors(node_group_successors, node_count)
    for entry in eachindex(km, ykm)
        node = abs(km[entry])
        1 <= node <= node_count ||
            throw(ArgumentError("km entries must point to a valid network node"))
        isfinite(Float64(ykm[entry])) ||
            throw(ArgumentError("ykm entries must be finite"))
    end

    solution = Float64.(initial_solution)
    solution[1:partition_boundary] .= Float64.(rhs[1:partition_boundary])
    for value in solution
        isfinite(value) || throw(ArgumentError("solution entries must be finite"))
    end

    alias_successors = Int.(node_group_successors)
    active_node = 1
    visited = 0
    while true
        solution[active_node] = 0.0
        alias_successors[active_node] <= active_node && break
        active_node = alias_successors[active_node]
        visited += 1
        visited <= node_count ||
            throw(ArgumentError("node group reference-chain search did not terminate"))
    end

    for node in 2:partition_boundary
        successor = alias_successors[node]
        successor == 0 && continue
        successor > partition_boundary && continue
        successor > node && (solution[successor] += solution[node])
    end

    entry_index = 1
    entry_count = length(km)
    while entry_index <= entry_count
        row = abs(km[entry_index])
        scale = solution[row]
        solution[row] = scale * Float64(ykm[entry_index])
        row_end = kk[row]
        while true
            entry_index += 1
            entry_index > row_end && break
            column = km[entry_index]
            column > partition_boundary && continue
            solution[column] -= scale * Float64(ykm[entry_index])
        end
    end

    node = node_count
    while true
        successor = alias_successors[node]
        if successor != 0 && successor <= node
            copied_from = node
            copied_to = successor
            visited = 0
            while true
                solution[copied_to] = solution[copied_from]
                copied_from = copied_to
                copied_to = alias_successors[copied_from]
                copied_to == node && break
                visited += 1
                visited <= node_count ||
                    throw(ArgumentError("node group copy cycle did not terminate"))
            end
        end
        if node <= partition_boundary
            entry_index == 1 && break
            correction = 0.0
            while true
                entry_index -= 1
                column = km[entry_index]
                if column < 0
                    node = abs(column)
                    solution[node] += correction
                    break
                end
                correction -= solution[column] * Float64(ykm[entry_index])
            end
        else
            node -= 1
        end
    end

    return solution
end

function over16_network_solution_update!(
    state::OVER16SwitchCurrentState,
    factor::OVER16FortranSparseFactorWorkspaceState;
    kwargs...,
)
    before = copy(state.network_solution)
    solution = sparse_network_solution(factor, state.rhs; kwargs...)
    resize!(state.network_solution, length(solution))
    state.network_solution .= solution
    network_solution_mutated = state.network_solution != before
    return (
        rhs = copy(state.rhs),
        network_solution = copy(state.network_solution),
        network_solution_mutated = network_solution_mutated,
        switch_current_state_mutated = network_solution_mutated,
    )
end

function sparse_network_factorization_update(admittance::AbstractMatrix{<:Real}; kwargs...)
    return over16_fortran_sparse_factor_update(admittance; kwargs...)
end

function sparse_network_factorization_update!(
    state::OVER16FortranSparseFactorWorkspaceState,
    admittance::AbstractMatrix{<:Real};
    kwargs...,
)
    return over16_fortran_sparse_factor_update!(state, admittance; kwargs...)
end

function sparse_network_solution(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int},
    rhs::AbstractVector{<:Real};
    node_group_successors::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    initial_solution::Union{Nothing,AbstractVector{<:Real}}=nothing,
    kwargs...,
)
    if node_group_successors !== nothing
        seed = initial_solution === nothing ? zeros(Float64, length(kk)) : initial_solution
        return _grouped_sparse_network_solution(
            km,
            ykm,
            kk,
            rhs,
            node_group_successors,
            seed;
            kwargs...,
        )
    end
    return over16_fortran_sparse_factor_solve(km, ykm, kk, rhs; kwargs...)
end

function sparse_network_solution(
    state::OVER16FortranSparseFactorWorkspaceState,
    rhs::AbstractVector{<:Real};
    kwargs...,
)
    return sparse_network_solution(state.km, state.ykm, state.kk, rhs; kwargs...)
end

function sparse_network_solution_update!(
    state::OVER16SwitchCurrentState,
    factor::OVER16FortranSparseFactorWorkspaceState;
    kwargs...,
)
    return over16_network_solution_update!(state, factor; kwargs...)
end

function _over16_sparse_row_balance(
    node::Int,
    rhs::AbstractVector{<:Real},
    kks::AbstractVector{Int},
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    voltages::AbstractVector{<:Real},
)
    node_count = length(rhs)
    length(kks) >= node_count ||
        throw(ArgumentError("kks length must cover rhs nodes"))
    length(voltages) >= node_count ||
        throw(ArgumentError("voltages length must cover rhs nodes"))
    length(km) == length(ykm) ||
        throw(ArgumentError("km and ykm lengths must match"))

    pointer = kks[node]
    2 <= pointer <= length(km) + 1 ||
        throw(ArgumentError("kks row pointer must reference one past a stored row"))

    balance_current = -Float64(rhs[node])
    while pointer > 1
        pointer -= 1
        row_node = km[pointer]
        voltage_node = abs(row_node)
        1 <= voltage_node <= length(voltages) ||
            throw(ArgumentError("km row node is outside the voltage vector"))
        balance_current += Float64(ykm[pointer]) * Float64(voltages[voltage_node])
        if row_node < 0
            return balance_current
        end
    end
    throw(ArgumentError("km row is missing the negative row-start marker"))
end

function over16_switch_current_reconstruction(
    from_node::Int,
    to_node::Int,
    nextsw_entry::Int,
    rhs::AbstractVector{<:Real},
    kks::AbstractVector{Int},
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    voltages::AbstractVector{<:Real},
    previous_switch_current::Real=0.0,
)
    node_count = length(rhs)
    1 <= from_node <= node_count ||
        throw(ArgumentError("from_node must be within rhs"))
    1 <= to_node <= node_count ||
        throw(ArgumentError("to_node must be within rhs"))
    from_node != to_node ||
        throw(ArgumentError("switch endpoints must be distinct"))
    nextsw_entry != 0 ||
        throw(ArgumentError("nextsw_entry must be nonzero for a closed switch"))

    kcl_node = nextsw_entry < 0 ? to_node : from_node
    opposite_node = kcl_node == from_node ? to_node : from_node
    row_balance_current =
        _over16_sparse_row_balance(kcl_node, rhs, kks, km, ykm, voltages)

    updated_rhs = Float64.(rhs)
    updated_rhs[opposite_node] -= row_balance_current
    switch_current = nextsw_entry > 0 ? -row_balance_current : row_balance_current
    current_product = switch_current * Float64(previous_switch_current)
    if current_product == 0.0 && previous_switch_current != 0.0
        current_product = -1.0
    end

    return (
        kcl_node = kcl_node,
        opposite_node = opposite_node,
        selected_endpoint_index = kcl_node == from_node ? 1 : 2,
        row_balance_current = row_balance_current,
        switch_current = switch_current,
        previous_current_product = current_product,
        updated_rhs = updated_rhs,
        topology_mutated = false,
        admittance_mutated = false,
    )
end

function over16_switch_current_reconstruction_table(
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    nextsw::AbstractVector{Int},
    first_group_head::Int,
    rhs::AbstractVector{<:Real},
    kks::AbstractVector{Int},
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    voltages::AbstractVector{<:Real};
    previous_switch_currents::AbstractVector{<:Real}=Float64[],
)
    switch_count = length(nextsw)
    _over16_check_length("from_nodes", from_nodes, switch_count)
    _over16_check_length("to_nodes", to_nodes, switch_count)
    _over16_check_optional_length("previous_switch_currents", previous_switch_currents, switch_count)
    0 <= first_group_head <= switch_count ||
        throw(ArgumentError("first_group_head must be between 0 and switch count"))

    updated_rhs = Float64.(rhs)
    switch_currents = zeros(Float64, switch_count)
    row_balance_currents = zeros(Float64, switch_count)
    current_products = zeros(Float64, switch_count)
    kcl_endpoint_indices = zeros(Int, switch_count)
    kcl_nodes = zeros(Int, switch_count)
    opposite_nodes = zeros(Int, switch_count)
    ordered_rows = Int[]
    visited = falses(switch_count)

    if first_group_head == 0
        any(!=(0), nextsw) &&
            throw(ArgumentError("first_group_head is zero but nextsw contains closed rows"))
        return (
            switch_currents = switch_currents,
            row_balance_currents = row_balance_currents,
            current_products = current_products,
            updated_rhs = updated_rhs,
            ordered_rows = ordered_rows,
            kcl_endpoint_indices = kcl_endpoint_indices,
            kcl_nodes = kcl_nodes,
            opposite_nodes = opposite_nodes,
            topology_mutated = false,
            admittance_mutated = false,
        )
    end
    nextsw[first_group_head] != 0 ||
        throw(ArgumentError("first_group_head must reference a closed nextsw row"))

    row = first_group_head
    for _ in 1:switch_count
        nextsw[row] != 0 ||
            throw(ArgumentError("NEXTSW ring reached an open row"))
        !visited[row] ||
            throw(ArgumentError("NEXTSW ring repeats before returning to first_group_head"))
        previous_current =
            isempty(previous_switch_currents) ? 0.0 : previous_switch_currents[row]
        reconstruction = over16_switch_current_reconstruction(
            from_nodes[row],
            to_nodes[row],
            nextsw[row],
            updated_rhs,
            kks,
            km,
            ykm,
            voltages,
            previous_current,
        )

        visited[row] = true
        push!(ordered_rows, row)
        updated_rhs = reconstruction.updated_rhs
        switch_currents[row] = reconstruction.switch_current
        row_balance_currents[row] = reconstruction.row_balance_current
        current_products[row] = reconstruction.previous_current_product
        kcl_endpoint_indices[row] = reconstruction.selected_endpoint_index
        kcl_nodes[row] = reconstruction.kcl_node
        opposite_nodes[row] = reconstruction.opposite_node

        next_row = abs(nextsw[row])
        1 <= next_row <= switch_count ||
            throw(ArgumentError("NEXTSW successor row is out of range"))
        if next_row == first_group_head
            for i in 1:switch_count
                if nextsw[i] != 0 && !visited[i]
                    throw(ArgumentError("NEXTSW contains a closed row outside the circular order"))
                end
            end
            return (
                switch_currents = switch_currents,
                row_balance_currents = row_balance_currents,
                current_products = current_products,
                updated_rhs = updated_rhs,
                ordered_rows = ordered_rows,
                kcl_endpoint_indices = kcl_endpoint_indices,
                kcl_nodes = kcl_nodes,
                opposite_nodes = opposite_nodes,
                topology_mutated = false,
                admittance_mutated = false,
            )
        end
        row = next_row
    end
    throw(ArgumentError("NEXTSW ring did not return to first_group_head"))
end

function over16_switch_current_reconstruction_table!(
    state::OVER16SwitchCurrentState,
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    nextsw::AbstractVector{Int},
    first_group_head::Int,
    kks::AbstractVector{Int},
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    voltages::AbstractVector{<:Real},
)
    rhs_before = copy(state.rhs)
    currents_before = copy(state.switch_currents)
    products_before = copy(state.current_products)
    preview = over16_switch_current_reconstruction_table(
        from_nodes,
        to_nodes,
        nextsw,
        first_group_head,
        state.rhs,
        kks,
        km,
        ykm,
        voltages;
        previous_switch_currents = state.switch_currents,
    )

    resize!(state.rhs, length(preview.updated_rhs))
    state.rhs .= preview.updated_rhs
    resize!(state.switch_currents, length(preview.switch_currents))
    state.switch_currents .= preview.switch_currents
    resize!(state.current_products, length(preview.current_products))
    state.current_products .= preview.current_products

    rhs_mutated = state.rhs != rhs_before
    switch_currents_mutated = state.switch_currents != currents_before
    current_products_mutated = state.current_products != products_before
    switch_current_state_mutated =
        rhs_mutated || switch_currents_mutated || current_products_mutated
    return merge(
        preview,
        (
            updated_rhs = copy(state.rhs),
            switch_currents = copy(state.switch_currents),
            current_products = copy(state.current_products),
            rhs_mutated = rhs_mutated,
            switch_currents_mutated = switch_currents_mutated,
            current_products_mutated = current_products_mutated,
            switch_current_state_mutated = switch_current_state_mutated,
            topology_mutated = false,
            admittance_mutated = false,
            sparse_factor_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function over16_switch_post_current_transition(
    position::Int,
    switch_current::Real,
    current_product::Real,
    energy::Real,
    source_voltage_difference::Real,
    source_index::Int,
    open_time::Real,
    critical_current::Real,
    delay_time::Real,
    t::Real,
    dt::Real,
)
    dt > 0.0 || throw(ArgumentError("dt must be positive"))

    state = abs(position)
    current = Float64(switch_current)
    product_signal = state == 3 ? -current : Float64(current_product)
    vsl_value = position <= 0 ? current : 0.0
    updated_energy = Float64(energy)
    source_power_increment = 0.0
    blocked_by_topen = false
    blocked_by_delay = false
    open_reason = :none

    if !(1 <= state <= 3)
        return (
            position = position,
            switch_current = current,
            energy = updated_energy,
            opened = false,
            open_reason = :skipped_state,
            current_product_signal = product_signal,
            source_power_increment = source_power_increment,
            vsl_value = vsl_value,
            tstop_reset = false,
            energy_report_value = 0.0,
            blocked_by_topen = false,
            blocked_by_delay = false,
            modswt_sign = 0,
            altered = false,
        )
    end

    if source_index > 0
        source_power_increment = Float64(source_voltage_difference) * current
        updated_energy += source_power_increment
        if source_power_increment < 0.0
            open_reason = :source_energy_reversal
        end
    end

    if open_reason == :none
        if state > 1 && Float64(t) < Float64(open_time)
            blocked_by_topen = true
        else
            critical_requested = abs(current) < Float64(critical_current)
            if critical_requested
                product_signal = -1.0
            end
            if Float64(t) < Float64(delay_time)
                product_signal = 1.0
                blocked_by_delay = true
            end

            if product_signal < 0.0
                if critical_requested
                    open_reason = :critical_current
                elseif state == 3
                    open_reason = :state3_current
                else
                    open_reason = :current_reversal
                end
            end
        end
    end

    if open_reason == :none
        return (
            position = position,
            switch_current = current,
            energy = updated_energy,
            opened = false,
            open_reason = blocked_by_topen ? :blocked_by_topen :
                          blocked_by_delay ? :blocked_by_delay : :none,
            current_product_signal = product_signal,
            source_power_increment = source_power_increment,
            vsl_value = vsl_value,
            tstop_reset = false,
            energy_report_value = 0.0,
            blocked_by_topen = blocked_by_topen,
            blocked_by_delay = blocked_by_delay,
            modswt_sign = 0,
            altered = false,
        )
    end

    new_state = state + 1
    if new_state == 2
        new_state = 10
    elseif new_state == 3
        new_state = 5
    end
    updated_position = position < 0 ? -new_state : new_state
    return (
        position = updated_position,
        switch_current = 0.0,
        energy = 0.0,
        opened = true,
        open_reason = open_reason,
        current_product_signal = product_signal,
        source_power_increment = source_power_increment,
        vsl_value = vsl_value,
        tstop_reset = source_index != 0 && state == 1,
        energy_report_value = state == 1 ? updated_energy * Float64(dt) : 0.0,
        blocked_by_topen = blocked_by_topen,
        blocked_by_delay = blocked_by_delay,
        modswt_sign = -1,
        altered = true,
    )
end

function over16_switch_post_current_transition_table(
    positions::AbstractVector{Int},
    nextsw::AbstractVector{Int},
    switch_currents::AbstractVector{<:Real},
    current_products::AbstractVector{<:Real},
    energies::AbstractVector{<:Real},
    source_voltage_differences::AbstractVector{<:Real},
    source_indices::AbstractVector{Int},
    open_times::AbstractVector{<:Real},
    critical_currents::AbstractVector{<:Real},
    delay_times::AbstractVector{<:Real},
    t::Real,
    dt::Real;
    scan_rows::AbstractVector{Int}=Int[],
)
    switch_count = length(positions)
    _over16_check_length("nextsw", nextsw, switch_count)
    _over16_check_length("switch_currents", switch_currents, switch_count)
    _over16_check_length("current_products", current_products, switch_count)
    _over16_check_length("energies", energies, switch_count)
    _over16_check_length(
        "source_voltage_differences",
        source_voltage_differences,
        switch_count,
    )
    _over16_check_length("source_indices", source_indices, switch_count)
    _over16_check_length("open_times", open_times, switch_count)
    _over16_check_length("critical_currents", critical_currents, switch_count)
    _over16_check_length("delay_times", delay_times, switch_count)
    dt > 0.0 || throw(ArgumentError("dt must be positive"))

    ordered_rows = isempty(scan_rows) ? findall(!=(0), nextsw) : collect(scan_rows)
    seen_rows = falses(switch_count)
    for row in ordered_rows
        1 <= row <= switch_count ||
            throw(ArgumentError("scan_rows entries must be within switch table"))
        !seen_rows[row] ||
            throw(ArgumentError("scan_rows entries must not repeat"))
        seen_rows[row] = true
        nextsw[row] != 0 ||
            throw(ArgumentError("scan_rows entries must reference closed NEXTSW rows"))
    end

    updated_positions = collect(positions)
    updated_switch_currents = Float64.(switch_currents)
    updated_energies = Float64.(energies)
    source_power_increments = zeros(Float64, switch_count)
    vsl_values = zeros(Float64, switch_count)
    energy_report_values = zeros(Float64, switch_count)
    open_reasons = fill(:none, switch_count)
    blocked_by_topen = falses(switch_count)
    blocked_by_delay = falses(switch_count)
    opened_rows = Int[]
    tstop_reset_rows = Int[]
    modswt = Int[]

    for row in ordered_rows
        transition = over16_switch_post_current_transition(
            updated_positions[row],
            updated_switch_currents[row],
            current_products[row],
            updated_energies[row],
            source_voltage_differences[row],
            source_indices[row],
            open_times[row],
            critical_currents[row],
            delay_times[row],
            t,
            dt,
        )
        updated_positions[row] = transition.position
        updated_switch_currents[row] = transition.switch_current
        updated_energies[row] = transition.energy
        source_power_increments[row] = transition.source_power_increment
        vsl_values[row] = transition.vsl_value
        energy_report_values[row] = transition.energy_report_value
        open_reasons[row] = transition.open_reason
        blocked_by_topen[row] = transition.blocked_by_topen
        blocked_by_delay[row] = transition.blocked_by_delay
        if transition.opened
            push!(opened_rows, row)
            push!(modswt, -row)
        end
        if transition.tstop_reset
            push!(tstop_reset_rows, row)
        end
    end

    return (
        positions = updated_positions,
        switch_currents = updated_switch_currents,
        energies = updated_energies,
        source_power_increments = source_power_increments,
        vsl_values = vsl_values,
        energy_report_values = energy_report_values,
        open_reasons = open_reasons,
        blocked_by_topen = blocked_by_topen,
        blocked_by_delay = blocked_by_delay,
        opened_rows = opened_rows,
        tstop_reset_rows = tstop_reset_rows,
        modswt = modswt,
        scan_rows = ordered_rows,
        altered = !isempty(modswt),
        topology_mutated = false,
        admittance_mutated = false,
        output_mutated = false,
        tacs_executed = false,
    )
end

function over16_switch_post_current_transition_table!(
    state::OVER16SwitchPostCurrentState,
    nextsw::AbstractVector{Int},
    current_products::AbstractVector{<:Real},
    source_voltage_differences::AbstractVector{<:Real},
    source_indices::AbstractVector{Int},
    open_times::AbstractVector{<:Real},
    critical_currents::AbstractVector{<:Real},
    delay_times::AbstractVector{<:Real},
    t::Real,
    dt::Real;
    scan_rows::AbstractVector{Int}=Int[],
)
    positions_before = copy(state.positions)
    currents_before = copy(state.switch_currents)
    energies_before = copy(state.energies)
    modswt_before = copy(state.modswt)
    preview = over16_switch_post_current_transition_table(
        state.positions,
        nextsw,
        state.switch_currents,
        current_products,
        state.energies,
        source_voltage_differences,
        source_indices,
        open_times,
        critical_currents,
        delay_times,
        t,
        dt;
        scan_rows = scan_rows,
    )

    state.positions .= preview.positions
    state.switch_currents .= preview.switch_currents
    state.energies .= preview.energies
    empty!(state.modswt)
    append!(state.modswt, preview.modswt)

    positions_mutated = state.positions != positions_before
    switch_currents_mutated = state.switch_currents != currents_before
    energies_mutated = state.energies != energies_before
    modswt_mutated = state.modswt != modswt_before
    switch_post_current_state_mutated =
        positions_mutated || switch_currents_mutated ||
        energies_mutated || modswt_mutated
    return merge(
        preview,
        (
            positions = copy(state.positions),
            switch_currents = copy(state.switch_currents),
            energies = copy(state.energies),
            modswt = copy(state.modswt),
            positions_mutated = positions_mutated,
            switch_currents_mutated = switch_currents_mutated,
            energies_mutated = energies_mutated,
            modswt_mutated = modswt_mutated,
            switch_post_current_state_mutated = switch_post_current_state_mutated,
            switch_graph_state_mutated = positions_mutated || modswt_mutated,
            topology_mutated = false,
            admittance_mutated = false,
            sparse_factor_mutated = false,
            output_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function over16_switch_alteration_rebuild_intent(
    ialter::Int,
    ktrlsw1::Int,
    ktrlsw2::Int,
    ktrlsw3::Int,
    ktrlsw4::Int=0,
    ktrlsw5::Int=0,
    ktrlsw6::Int=1;
    m4plot::Bool=false,
    yserlc_altered::Bool=false,
    kanal::Int=0,
)
    ialter in (0, 1) ||
        throw(ArgumentError("ialter must be 0 or 1"))
    ktrlsw1 >= 0 ||
        throw(ArgumentError("ktrlsw1 operation count must be nonnegative"))
    ktrlsw2 >= 0 ||
        throw(ArgumentError("ktrlsw2 closed-switch count must be nonnegative"))
    ktrlsw3 >= 0 ||
        throw(ArgumentError("ktrlsw3 triangularization count must be nonnegative"))
    ktrlsw4 >= 0 ||
        throw(ArgumentError("ktrlsw4 first-group head must be nonnegative"))
    ktrlsw5 >= 0 ||
        throw(ArgumentError("ktrlsw5 total operation count must be nonnegative"))
    ktrlsw6 in (0, 1) ||
        throw(ArgumentError("ktrlsw6 switch-logic mode must be 0 or 1"))
    (!yserlc_altered || m4plot) ||
        throw(ArgumentError("yserlc_altered requires m4plot=true"))

    effective_ialter = (ialter == 1 || (m4plot && yserlc_altered)) ? 1 : 0
    should_call_switch = effective_ialter != 0 && ktrlsw1 > 0
    switch_path = should_call_switch ? (ktrlsw6 == 0 ? :sophisticated : :simple) : :none
    updated_triangularization_count =
        effective_ialter != 0 ? ktrlsw3 + 1 : ktrlsw3
    updated_total_operation_count =
        should_call_switch ? ktrlsw5 + ktrlsw1 : ktrlsw5

    deferred_calls = Symbol[]
    if should_call_switch
        push!(deferred_calls, :switch)
    end
    if effective_ialter != 0 && kanal == 2
        push!(deferred_calls, :last14)
    end
    if effective_ialter != 0
        push!(deferred_calls, :retriangularization)
    end

    return (
        should_call_yserlc = m4plot,
        series_rlc_parameter_mutated = m4plot && yserlc_altered,
        ialter = ialter,
        effective_ialter = effective_ialter,
        ialter_after_last14 = effective_ialter != 0 && kanal == 2 ? 1 : effective_ialter,
        altered = effective_ialter != 0,
        operation_count = ktrlsw1,
        closed_switch_count = ktrlsw2,
        triangularization_count = updated_triangularization_count,
        first_group_head = ktrlsw4,
        total_operation_count = updated_total_operation_count,
        ktrlsw6 = ktrlsw6,
        should_call_switch = should_call_switch,
        switch_path = switch_path,
        lastsw_required = should_call_switch && ktrlsw6 == 0,
        should_call_last14 = effective_ialter != 0 && kanal == 2,
        should_retriangularize = effective_ialter != 0,
        should_clear_factor_workspace = effective_ialter != 0,
        should_continue_to_subts1_exit = true,
        normal_next_nchain = 17,
        kill_next_nchain = 51,
        deferred_calls = deferred_calls,
        topology_mutated = false,
        admittance_mutated = m4plot && yserlc_altered,
        sparse_factor_mutated = false,
        tacs_executed = false,
        solvum_executed = false,
    )
end

function over16_switch_alteration_rebuild_update!(
    state::OVER16SwitchAlterationState;
    m4plot::Bool=false,
    yserlc_altered::Bool=false,
    kanal::Int=0,
)
    ialter_before = state.ialter
    triangularization_count_before = state.triangularization_count
    total_operation_count_before = state.total_operation_count
    preview = over16_switch_alteration_rebuild_intent(
        state.ialter,
        state.operation_count,
        state.closed_switch_count,
        state.triangularization_count,
        state.first_group_head,
        state.total_operation_count,
        state.ktrlsw6;
        m4plot = m4plot,
        yserlc_altered = yserlc_altered,
        kanal = kanal,
    )

    state.ialter = preview.ialter_after_last14
    state.triangularization_count = preview.triangularization_count
    state.total_operation_count = preview.total_operation_count

    ialter_mutated = state.ialter != ialter_before
    triangularization_count_mutated =
        state.triangularization_count != triangularization_count_before
    total_operation_count_mutated =
        state.total_operation_count != total_operation_count_before
    switch_alteration_state_mutated =
        ialter_mutated || triangularization_count_mutated ||
        total_operation_count_mutated
    return merge(
        preview,
        (
            ialter = state.ialter,
            triangularization_count = state.triangularization_count,
            total_operation_count = state.total_operation_count,
            ialter_mutated = ialter_mutated,
            triangularization_count_mutated = triangularization_count_mutated,
            total_operation_count_mutated = total_operation_count_mutated,
            switch_alteration_state_mutated = switch_alteration_state_mutated,
            topology_mutated = false,
            admittance_mutated = false,
            sparse_factor_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function over16_switch_bvalue_export(
    kpos::AbstractVector{Int},
    nextsw::AbstractVector{Int},
    switch_currents::AbstractVector{<:Real},
    base_count::Int=0,
)
    switch_count = length(kpos)
    _over16_check_length("nextsw", nextsw, switch_count)
    _over16_check_length("switch_currents", switch_currents, switch_count)
    base_count >= 0 || throw(ArgumentError("base_count must be nonnegative"))

    export_rows = Int[]
    bvalue_indices = Int[]
    values = Float64[]
    count = base_count
    for row in 1:switch_count
        if kpos[row] >= 0
            continue
        end
        count += 1
        push!(export_rows, row)
        push!(bvalue_indices, count)
        push!(values, nextsw[row] == 0 ? 0.0 : Float64(switch_currents[row]))
    end

    return (
        export_rows = export_rows,
        bvalue_indices = bvalue_indices,
        values = values,
        final_count = count,
        output_mutated = false,
        tacs_executed = false,
    )
end

function over16_switch_bvalue_export!(
    state::OVER16SwitchBValueExportState,
    kpos::AbstractVector{Int},
    nextsw::AbstractVector{Int},
    switch_currents::AbstractVector{<:Real},
)
    bvalue_before = copy(state.bvalue)
    output_count_before = state.output_count
    preview = over16_switch_bvalue_export(
        kpos,
        nextsw,
        switch_currents,
        state.output_count,
    )

    if length(state.bvalue) < preview.final_count
        resize!(state.bvalue, preview.final_count)
    end
    for (index, value) in zip(preview.bvalue_indices, preview.values)
        state.bvalue[index] = value
    end
    state.output_count = preview.final_count

    bvalue_mutated = state.bvalue != bvalue_before
    output_count_mutated = state.output_count != output_count_before
    switch_bvalue_state_mutated = bvalue_mutated || output_count_mutated
    return merge(
        preview,
        (
            bvalue = copy(state.bvalue),
            output_count = state.output_count,
            bvalue_mutated = bvalue_mutated,
            output_count_mutated = output_count_mutated,
            switch_bvalue_state_mutated = switch_bvalue_state_mutated,
            output_mutated = switch_bvalue_state_mutated,
            topology_mutated = false,
            admittance_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function stamp!(y::AbstractMatrix{Float64}, _rhs::AbstractVector{Float64},
                s::Union{IdealSwitch,TimeSwitch,CurrentZeroSwitch}, t::Float64, _dt::Float64)
    stamp_conductance!(y, s.a, s.b, switch_conductance(s, t))
    return nothing
end

backward_euler_companion_supported(
    ::Union{IdealSwitch,TimeSwitch,CurrentZeroSwitch},
) = true

update!(::Union{IdealSwitch,TimeSwitch}, _voltages::AbstractVector{Float64},
        _dt::Float64) = nothing

update!(::CurrentZeroSwitch, _voltages::AbstractVector{Float64}, _dt::Float64) = nothing

end

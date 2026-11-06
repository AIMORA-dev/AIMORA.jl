module Branches

using ..Companion
using ..StudyCore: ContinuousEvaluationTime, continuous_time_seconds
using LinearAlgebra: I, Symmetric, cholesky, isposdef, mul!

export ConductanceBranch,
       CompanionIntegrationMethod,
       TrapezoidalCompanion,
       BackwardEulerCompanion,
       IdealTransformerVoltageConstraint,
       SeriesRLBranch,
       SeriesRLCBranch,
       CoupledInductiveBranch,
       CoupledSeriesRLBranch,
       CapacitorBranch,
       GeneratorEquivalentModalBranch,
       BreqivHistoryInjection,
       BranchCompanionSnapshot,
       TheveninSource,
       TwoTerminalTheveninSource,
       CurrentInjection,
       branch_companion_snapshot,
       branch_current_value,
       single_phase_breqiv_history_injection,
       three_phase_breqiv_history_injection,
       generator_equivalent_history_injection,
       seed_breqiv_frequency_histories!,
       trace_output_channel_count,
       trace_output_channel_names!,
       trace_output_is_public,
       trace_output_values!,
       advance_breqiv_history_current!,
       stamp!,
       update!,
       refresh_accepted_dependents!

abstract type EMTElement end

"""Whether a stamp depends on step/history/topology, but not on absolute time.

Unknown owners must implement precise timed stamping or explicitly establish
this property; the nonlinear solver never drops a nonzero time correction.
"""
continuous_time_invariant(::EMTElement) = false

"""Refresh cached dependent fields after all coupled nonlinear devices accept a step."""
refresh_accepted_dependents!(::EMTElement, _voltages, _dt_s::Float64) = nothing

@enum CompanionIntegrationMethod begin
    TrapezoidalCompanion
    BackwardEulerCompanion
end

struct ConductanceBranch <: EMTElement
    a::Int
    b::Int
    g::Float64
end

"""
    IdealTransformerVoltageConstraint(primary_positive, secondary_positive,
                                      secondary_negative, primary_negative,
                                      constraint_node, turns_ratio)

Modified-nodal voltage constraint used by a four-terminal ideal transformer.
The auxiliary `constraint_node` is the modified-nodal Lagrange multiplier for
transformer current, while the terminal equation is

`Vpp + Vsp / ratio - Vsn / ratio - Vpn = -Vsource`.

The sign follows the BPA type-18 source topology: the source value is assigned
to the auxiliary right-hand-side entry and the terminal coefficients are
stamped with the opposite sign. The turns ratio is dimensionless.
"""
struct IdealTransformerVoltageConstraint <: EMTElement
    terminal_nodes::NTuple{4,Int}
    constraint_node::Int
    turns_ratio::Float64
    terminal_coefficients::NTuple{4,Float64}

    function IdealTransformerVoltageConstraint(
        primary_positive::Integer,
        secondary_positive::Integer,
        secondary_negative::Integer,
        primary_negative::Integer,
        constraint_node::Integer,
        turns_ratio::Real,
    )
        nodes = (
            Int(primary_positive),
            Int(secondary_positive),
            Int(secondary_negative),
            Int(primary_negative),
        )
        auxiliary = Int(constraint_node)
        all(>(0), nodes) || throw(ArgumentError(
            "ideal-transformer voltage-constraint terminals must be non-ground nodes",
        ))
        auxiliary > 0 || throw(ArgumentError(
            "ideal-transformer voltage-constraint auxiliary node must be positive",
        ))
        auxiliary in nodes && throw(ArgumentError(
            "ideal-transformer voltage-constraint auxiliary node must be distinct",
        ))
        ratio = Float64(turns_ratio)
        isfinite(ratio) && ratio != 0.0 || throw(ArgumentError(
            "ideal-transformer turns ratio must be finite and nonzero",
        ))
        inverse_ratio = inv(ratio)
        return new(
            nodes,
            auxiliary,
            ratio,
            (1.0, inverse_ratio, -inverse_ratio, -1.0),
        )
    end
end

struct BranchCompanionSnapshot
    kind::Symbol
    a::Int
    b::Int
    conductance::Float64
    history_current::Float64
    branch_voltage::Float64
    branch_current::Float64
    previous_current::Float64
    previous_voltage::Float64
end

mutable struct SeriesRLBranch <: EMTElement
    a::Int
    b::Int
    r::Float64
    l::Float64
    # The compensated path retains the physical current before Float64 history
    # rounding. Keep the same seven-field layout for existing raw checkpoints;
    # ordinary starts and compatibility updates still use Float64.
    i_prev::Union{Float64,BigFloat}
    v_prev::Float64
    i_last::Union{Float64,BigFloat}

    function SeriesRLBranch(
        a::Int,
        b::Int,
        r::Float64,
        l::Float64,
        i_prev::Float64,
        v_prev::Float64,
        i_last::Float64,
    )
        isfinite(r) ||
            throw(ArgumentError("series R-L resistance must be finite"))
        isfinite(l) && l >= 0.0 ||
            throw(ArgumentError("series R-L inductance must be finite and nonnegative"))
        l > 0.0 || r > 0.0 ||
            throw(ArgumentError("series resistor-only branch resistance must be positive"))
        all(isfinite, (i_prev, v_prev, i_last)) ||
            throw(ArgumentError("series R-L accepted state must be finite"))
        return new(a, b, r, l, i_prev, v_prev, i_last)
    end
end

SeriesRLBranch(a::Int, b::Int, r::Float64, l::Float64) = SeriesRLBranch(a, b, r, l, 0.0, 0.0, 0.0)

mutable struct SeriesRLCBranch <: EMTElement
    a::Int
    b::Int
    r::Float64
    l::Float64
    c::Float64
    i_prev::Float64
    inductor_voltage_prev::Float64
    capacitor_voltage_prev::Float64
    v_prev::Float64
    i_last::Float64
end

function SeriesRLCBranch(a::Int, b::Int, r::Float64, l::Float64, c::Float64)
    r >= 0.0 || throw(ArgumentError("series RLC resistance must be nonnegative"))
    l >= 0.0 || throw(ArgumentError("series RLC inductance must be nonnegative"))
    c > 0.0 || throw(ArgumentError("series RLC capacitance must be positive"))
    return SeriesRLCBranch(a, b, r, l, c, 0.0, 0.0, 0.0, 0.0, 0.0)
end

mutable struct CoupledInductiveBranch <: EMTElement
    a::Vector{Int}
    b::Vector{Int}
    susceptance::Matrix{Float64}
    angular_frequency::Float64
    series_resistance::Float64
    resistance_reference_port::Int
    previous_current::Vector{Float64}
    previous_voltage::Vector{Float64}
    last_current::Vector{Float64}
    conductance_workspace::Matrix{Float64}
    history_current_workspace::Vector{Float64}
    port_voltage_workspace::Vector{Float64}
    current_workspace::Vector{Float64}
end

"""
    CoupledSeriesRLBranch(a, b, resistance_matrix, inductance_matrix)

An explicitly coupled, passive series R-L branch whose port equation is
`v = R*i + L*di/dt`. The trapezoidal companion is
`G = inv(R + 2L/dt)` with history current
`G * (v_previous + (2L/dt - R) * i_previous)`.
"""
mutable struct CoupledSeriesRLBranch <: EMTElement
    a::Vector{Int}
    b::Vector{Int}
    resistance_matrix::Matrix{Float64}
    inductance_matrix::Matrix{Float64}
    previous_current::Vector{Float64}
    previous_voltage::Vector{Float64}
    last_current::Vector{Float64}
    conductance_workspace::Matrix{Float64}
    history_current_workspace::Vector{Float64}
    port_voltage_workspace::Vector{Float64}
    current_workspace::Vector{Float64}
    cached_dt_s::Float64
end

function _symmetric_coupled_matrix(
    values::AbstractMatrix{<:Real},
    port_count::Int,
    field::AbstractString,
)
    size(values) == (port_count, port_count) ||
        throw(ArgumentError("$field matrix size must match coupled R-L ports"))
    matrix = Matrix{Float64}(values)
    all(isfinite, matrix) ||
        throw(ArgumentError("$field matrix entries must be finite"))
    tolerance = 64.0 * eps(Float64) * max(maximum(abs, matrix; init=0.0), 1.0)
    maximum(abs, matrix - transpose(matrix); init=0.0) <= tolerance ||
        throw(ArgumentError("$field matrix must be symmetric"))
    return 0.5 .* (matrix .+ transpose(matrix))
end

function CoupledSeriesRLBranch(
    a::AbstractVector{<:Integer},
    b::AbstractVector{<:Integer},
    resistance_matrix::AbstractMatrix{<:Real},
    inductance_matrix::AbstractMatrix{<:Real},
)
    port_count = length(a)
    port_count > 0 ||
        throw(ArgumentError("coupled series R-L branch requires at least one port"))
    length(b) == port_count ||
        throw(ArgumentError("coupled series R-L terminal counts must match"))
    from_nodes = Int.(a)
    to_nodes = Int.(b)
    all(>=(0), from_nodes) && all(>=(0), to_nodes) ||
        throw(ArgumentError("coupled series R-L nodes must be nonnegative"))
    resistance =
        _symmetric_coupled_matrix(resistance_matrix, port_count, "resistance")
    inductance =
        _symmetric_coupled_matrix(inductance_matrix, port_count, "inductance")
    resistance_tolerance =
        64.0 * eps(Float64) * max(maximum(abs, resistance; init=0.0), 1.0)
    isposdef(Symmetric(resistance + resistance_tolerance * I)) ||
        throw(ArgumentError("coupled series R-L resistance matrix must be positive semidefinite"))
    isposdef(Symmetric(inductance)) ||
        throw(ArgumentError("coupled series R-L inductance matrix must be positive definite"))
    return CoupledSeriesRLBranch(
        from_nodes,
        to_nodes,
        resistance,
        inductance,
        zeros(Float64, port_count),
        zeros(Float64, port_count),
        zeros(Float64, port_count),
        zeros(Float64, port_count, port_count),
        zeros(Float64, port_count),
        zeros(Float64, port_count),
        zeros(Float64, port_count),
        NaN,
    )
end

function CoupledInductiveBranch(
    a::AbstractVector{<:Integer},
    b::AbstractVector{<:Integer},
    susceptance::AbstractMatrix{<:Real},
    angular_frequency::Real,
    ;
    series_resistance::Real=0.0,
    resistance_reference_port::Integer=1,
)
    port_count = length(a)
    port_count > 0 ||
        throw(ArgumentError("coupled inductive branch requires at least one port"))
    length(b) == port_count ||
        throw(ArgumentError("coupled inductive branch terminal counts must match"))
    size(susceptance, 1) == port_count && size(susceptance, 2) == port_count ||
        throw(ArgumentError("coupled inductive branch susceptance matrix size must match ports"))
    omega = Float64(angular_frequency)
    isfinite(omega) && omega > 0.0 ||
        throw(ArgumentError("coupled inductive branch angular_frequency must be positive"))
    matrix = Float64.(susceptance)
    all(isfinite, matrix) ||
        throw(ArgumentError("coupled inductive branch susceptance entries must be finite"))
    resistance = Float64(series_resistance)
    isfinite(resistance) && resistance >= 0.0 ||
        throw(ArgumentError("coupled inductive branch series resistance must be finite and nonnegative"))
    reference_port = Int(resistance_reference_port)
    1 <= reference_port <= port_count ||
        throw(ArgumentError("coupled inductive branch resistance reference port is outside its ports"))
    resistance > 0.0 && matrix[reference_port, reference_port] >= 0.0 &&
        throw(ArgumentError("resistive coupled inductive branch reference susceptance must be negative"))
    return CoupledInductiveBranch(
        Int.(a),
        Int.(b),
        matrix,
        omega,
        resistance,
        reference_port,
        zeros(Float64, port_count),
        zeros(Float64, port_count),
        zeros(Float64, port_count),
        zeros(Float64, port_count, port_count),
        zeros(Float64, port_count),
        zeros(Float64, port_count),
        zeros(Float64, port_count),
    )
end

mutable struct CapacitorBranch <: EMTElement
    a::Int
    b::Int
    c::Float64
    i_prev::Float64
    v_prev::Float64
    i_last::Float64
end

CapacitorBranch(a::Int, b::Int, c::Float64) = CapacitorBranch(a, b, c, 0.0, 0.0, 0.0)

mutable struct BreqivHistoryInjection <: EMTElement
    a::Vector{Int}
    b::Vector{Int}
    state::BreqivBranchSetState
    initial_phase_voltage::Vector{Float64}
    phase_voltage::Vector{Float64}
    phase_current::Vector{Float64}
    phase_admittance::Matrix{Float64}
    history_current_scale::Float64
    history_voltage_scale::Float64
    frequency_history_voltage_scale::Float64
    history_current_consumed_for_step::Bool
    initialized::Bool
end

struct GeneratorEquivalentModalBranch
    resistance_ohm::Float64
    inductance_h::Float64
    capacitance_f::Float64
    damping_resistance_ohm::Float64

    function GeneratorEquivalentModalBranch(
        resistance_ohm::Real,
        inductance_h::Real,
        capacitance_f::Real,
        damping_resistance_ohm::Real,
    )
        values = Float64[
            value for value in (
                resistance_ohm,
                inductance_h,
                capacitance_f,
                damping_resistance_ohm,
            )
        ]
        all(isfinite, values) ||
            throw(ArgumentError("generator-equivalent modal branch values must be finite"))
        all(>=(0.0), values) ||
            throw(ArgumentError("generator-equivalent modal branch values must be nonnegative"))
        any(>(0.0), values) ||
            throw(ArgumentError("generator-equivalent modal branch must contain a positive component"))
        return new(values...)
    end
end

function BreqivHistoryInjection(
    a::AbstractVector{<:Integer},
    b::AbstractVector{<:Integer},
    state::BreqivBranchSetState,
    initial_phase_voltage::AbstractVector{Float64},
    history_current_scale::Real = 1.0,
    ;
    history_voltage_scale::Real = 1.0,
    frequency_history_voltage_scale::Real = 1.0,
)
    nph = length(initial_phase_voltage)
    nph > 0 || throw(ArgumentError("initial_phase_voltage must not be empty"))
    length(a) == nph || throw(ArgumentError("from-node count must match initial_phase_voltage"))
    length(b) == nph || throw(ArgumentError("to-node count must match initial_phase_voltage"))
    length(state.modal) == nph || throw(ArgumentError("BREQIV state phase count mismatch"))
    scale = Float64(history_current_scale)
    isfinite(scale) || throw(ArgumentError("history_current_scale must be finite"))
    voltage_scale = Float64(history_voltage_scale)
    isfinite(voltage_scale) || throw(ArgumentError("history_voltage_scale must be finite"))
    frequency_voltage_scale = Float64(frequency_history_voltage_scale)
    isfinite(frequency_voltage_scale) ||
        throw(ArgumentError("frequency_history_voltage_scale must be finite"))
    return BreqivHistoryInjection(
        Int.(a),
        Int.(b),
        state,
        collect(initial_phase_voltage),
        zeros(Float64, nph),
        zeros(Float64, nph),
        zeros(Float64, nph, nph),
        scale,
        voltage_scale,
        frequency_voltage_scale,
        false,
        false,
    )
end

struct TheveninSource{F} <: EMTElement
    node::Int
    g::Float64
    value::F
end

"""A finite Thevenin source between two nonidentical terminals, with positive source voltage and delivered current oriented from `b` to `a`."""
struct TwoTerminalTheveninSource{F} <: EMTElement
    a::Int
    b::Int
    g::Float64
    value::F

    function TwoTerminalTheveninSource(
        a::Integer,
        b::Integer,
        conductance::Real,
        value::F,
    ) where {F}
        from_node = Int(a)
        to_node = Int(b)
        from_node >= 0 && to_node >= 0 && from_node != to_node || throw(ArgumentError(
            "two-terminal Thevenin source nodes must be nonnegative and distinct",
        ))
        g = Float64(conductance)
        isfinite(g) && g > 0.0 || throw(ArgumentError(
            "two-terminal Thevenin source conductance must be finite and positive",
        ))
        return new{F}(from_node, to_node, g, value)
    end
end

struct CurrentInjection{F} <: EMTElement
    node::Int
    value::F
end

backward_euler_companion_supported(::EMTElement) = false

backward_euler_companion_supported(
    ::Union{
        ConductanceBranch,
        IdealTransformerVoltageConstraint,
        SeriesRLBranch,
        SeriesRLCBranch,
        CapacitorBranch,
        TheveninSource,
        TwoTerminalTheveninSource,
        CurrentInjection,
    },
) = true

function stamp_conductance!(y::AbstractMatrix{<:Real}, a::Int, b::Int, g::Float64)
    if a != 0
        y[a, a] += g
    end
    if b != 0
        y[b, b] += g
    end
    if a != 0 && b != 0
        y[a, b] -= g
        y[b, a] -= g
    end
end

function stamp_history_current!(rhs::AbstractVector{<:Real}, a::Int, b::Int, ih::Float64)
    # ih is positive from a to b. RHS is current injected into the node.
    a != 0 && (rhs[a] -= ih)
    b != 0 && (rhs[b] += ih)
end

function stamp_history_current!(rhs::AbstractVector{<:Real}, a::Int, b::Int, ih::BigFloat)
    # Preserve the continuous source contribution until RHS storage rounds it.
    a != 0 && (rhs[a] -= ih)
    b != 0 && (rhs[b] += ih)
end

branch_voltage(v::AbstractVector{Float64}, a::Int, b::Int) = (a == 0 ? 0.0 : v[a]) - (b == 0 ? 0.0 : v[b])

function branch_companion_snapshot(b::ConductanceBranch, voltage::AbstractVector{Float64}, dt::Float64)
    vb = branch_voltage(voltage, b.a, b.b)
    return BranchCompanionSnapshot(
        :conductance,
        b.a,
        b.b,
        b.g,
        0.0,
        vb,
        b.g * vb,
        0.0,
        0.0,
    )
end

function branch_companion_snapshot(b::SeriesRLBranch, voltage::AbstractVector{Float64}, dt::Float64)
    g, ih = companion(b, dt)
    vb = branch_voltage(voltage, b.a, b.b)
    isfinite(vb) ||
        throw(ArgumentError("series R-L branch voltage must be finite"))
    current = g * vb + ih
    isfinite(current) ||
        throw(ArgumentError("series R-L branch current must be finite"))
    return BranchCompanionSnapshot(
        :series_rl,
        b.a,
        b.b,
        g,
        ih,
        vb,
        current,
        b.i_prev,
        b.v_prev,
    )
end

function branch_companion_snapshot(
    b::SeriesRLCBranch,
    voltage::AbstractVector{Float64},
    dt::Float64,
)
    g, ih = companion(b, dt)
    vb = branch_voltage(voltage, b.a, b.b)
    return BranchCompanionSnapshot(
        :series_rlc,
        b.a,
        b.b,
        g,
        ih,
        vb,
        g * vb + ih,
        b.i_prev,
        b.v_prev,
    )
end

function branch_companion_snapshot(
    b::CoupledInductiveBranch,
    voltage::AbstractVector{Float64},
    dt::Float64,
)
    conductance = _coupled_inductive_conductance_matrix!(
        b.conductance_workspace,
        b,
        dt,
    )
    history_current = _coupled_inductive_history_current!(
        b.history_current_workspace,
        b,
        conductance,
        dt,
    )
    branch_voltages = _coupled_branch_port_voltages!(
        b.port_voltage_workspace,
        voltage,
        b.a,
        b.b,
    )
    _coupled_matrix_vector_mul!(b.current_workspace, conductance, branch_voltages)
    @inbounds for index in eachindex(b.current_workspace, history_current)
        b.current_workspace[index] += history_current[index]
    end
    return BranchCompanionSnapshot(
        :coupled_inductive,
        b.a[1],
        b.b[1],
        conductance[1, 1],
        history_current[1],
        branch_voltages[1],
        b.current_workspace[1],
        b.previous_current[1],
        b.previous_voltage[1],
    )
end

function branch_companion_snapshot(
    b::CoupledSeriesRLBranch,
    voltage::AbstractVector{Float64},
    dt::Float64,
)
    conductance = _coupled_series_rl_conductance_matrix!(
        b.conductance_workspace,
        b,
        dt,
    )
    history_current = _coupled_series_rl_history_current!(
        b.history_current_workspace,
        b,
        conductance,
        dt,
    )
    branch_voltages = _coupled_branch_port_voltages!(
        b.port_voltage_workspace,
        voltage,
        b.a,
        b.b,
    )
    _coupled_matrix_vector_mul!(b.current_workspace, conductance, branch_voltages)
    b.current_workspace .+= history_current
    return BranchCompanionSnapshot(
        :coupled_series_rl,
        b.a[1],
        b.b[1],
        conductance[1, 1],
        history_current[1],
        branch_voltages[1],
        b.current_workspace[1],
        b.previous_current[1],
        b.previous_voltage[1],
    )
end

function branch_companion_snapshot(b::CapacitorBranch, voltage::AbstractVector{Float64}, dt::Float64)
    g, ih = companion(b, dt)
    vb = branch_voltage(voltage, b.a, b.b)
    return BranchCompanionSnapshot(
        :capacitor,
        b.a,
        b.b,
        g,
        ih,
        vb,
        g * vb + ih,
        b.i_prev,
        b.v_prev,
    )
end

branch_companion_snapshot(element::EMTElement, voltage::AbstractVector{Float64}, dt::Float64) = nothing

branch_current_value(
    b::ConductanceBranch,
    voltage::AbstractVector{Float64},
    ::Float64,
) = b.g * branch_voltage(voltage, b.a, b.b)

branch_current_value(
    b::SeriesRLBranch,
    ::AbstractVector{Float64},
    ::Float64,
) = isfinite(b.i_last) ?
    b.i_last :
    throw(ArgumentError("series R-L accepted current must be finite"))

branch_current_value(
    b::SeriesRLCBranch,
    ::AbstractVector{Float64},
    ::Float64,
) = b.i_last

branch_current_value(
    b::CoupledInductiveBranch,
    ::AbstractVector{Float64},
    ::Float64,
) = b.last_current[1]

branch_current_value(
    b::CoupledSeriesRLBranch,
    ::AbstractVector{Float64},
    ::Float64,
) = b.last_current[1]

branch_current_value(
    b::CapacitorBranch,
    ::AbstractVector{Float64},
    ::Float64,
) = b.i_last

function branch_current_value(
    element::EMTElement,
    voltage::AbstractVector{Float64},
    dt::Float64,
)
    snapshot = branch_companion_snapshot(element, voltage, dt)
    snapshot === nothing &&
        throw(ArgumentError("branch current output requires a current-owning element"))
    return snapshot.branch_current
end

trace_output_channel_count(::EMTElement) = 0

trace_output_is_public(::EMTElement) = false

trace_output_channel_count(::IdealTransformerVoltageConstraint) = 1

trace_output_channel_names!(names::Vector{Symbol}, element_name::Symbol, element::EMTElement) = names

function trace_output_channel_names!(
    names::Vector{Symbol},
    element_name::Symbol,
    ::IdealTransformerVoltageConstraint,
)
    push!(names, Symbol(String(element_name), "_constraint_voltage"))
    return names
end

function trace_output_values!(
    output::AbstractMatrix{Float64},
    first_channel::Int,
    sample::Int,
    element::EMTElement,
    voltage::AbstractVector{Float64},
)
    return first_channel
end

function trace_output_values!(
    output::AbstractMatrix{Float64},
    first_channel::Int,
    sample::Int,
    constraint::IdealTransformerVoltageConstraint,
    voltage::AbstractVector{Float64},
)
    value = 0.0
    for (node, coefficient) in
        zip(constraint.terminal_nodes, constraint.terminal_coefficients)
        value += coefficient * voltage[node]
    end
    output[first_channel, sample] = value
    return first_channel + 1
end

function single_phase_breqiv_history_injection(
    a::Int,
    b::Int,
    r::Real,
    l::Real,
    c::Real,
    rl::Real,
    previous_current::Real,
    capacitor_voltage::Real,
    initial_voltage::Real,
)
    rmfd = [Float64(r), Float64(l), Float64(c), Float64(rl), 0.0]
    cikfd = [0.0, Float64(previous_current), Float64(capacitor_voltage)]
    state = breqiv_branch_state_from_arrays(rmfd, cikfd, [1, 0], 1, 0, 0, 0)
    return BreqivHistoryInjection([a], [b], state, [Float64(initial_voltage)])
end

function three_phase_breqiv_history_injection(
    a1::Int,
    b1::Int,
    a2::Int,
    b2::Int,
    a3::Int,
    b3::Int,
    zero_r::Real,
    zero_l::Real,
    zero_c::Real,
    zero_rl::Real,
    positive_r::Real,
    positive_l::Real,
    positive_c::Real,
    positive_rl::Real,
    initial_v1::Real,
    initial_v2::Real,
    initial_v3::Real,
    ;
    zero_previous_current::Real=0.0,
    zero_capacitor_voltage::Real=0.0,
    positive_previous_current_1::Real=0.0,
    positive_capacitor_voltage_1::Real=0.0,
    positive_previous_current_2::Real=0.0,
    positive_capacitor_voltage_2::Real=0.0,
    history_current_scale::Real=1.0,
    history_voltage_scale::Real=1.0,
)
    rmfd = [
        Float64(zero_r),
        Float64(zero_l),
        Float64(zero_c),
        Float64(zero_rl),
        0.0,
        Float64(positive_r),
        Float64(positive_l),
        Float64(positive_c),
        Float64(positive_rl),
        0.0,
    ]
    cikfd = [
        0.0,
        Float64(zero_previous_current),
        Float64(zero_capacitor_voltage),
        0.0,
        Float64(positive_previous_current_1),
        Float64(positive_capacitor_voltage_1),
        0.0,
        Float64(positive_previous_current_2),
        Float64(positive_capacitor_voltage_2),
    ]
    state = breqiv_branch_state_from_arrays(rmfd, cikfd, [1, 1], 3, 0, 0, 0)
    return BreqivHistoryInjection(
        [a1, a2, a3],
        [b1, b2, b3],
        state,
        [Float64(initial_v1), Float64(initial_v2), Float64(initial_v3)],
        history_current_scale;
        history_voltage_scale = history_voltage_scale,
    )
end

function generator_equivalent_history_injection(
    from_nodes::AbstractVector{<:Integer},
    to_nodes::AbstractVector{<:Integer},
    zero_branches::AbstractVector{GeneratorEquivalentModalBranch},
    positive_branches::AbstractVector{GeneratorEquivalentModalBranch};
    initial_phase_voltage::AbstractVector{<:Real} = zeros(Float64, length(from_nodes)),
)
    nph = length(from_nodes)
    nph >= 2 || throw(ArgumentError("generator equivalent requires at least two phases"))
    length(to_nodes) == nph ||
        throw(ArgumentError("generator-equivalent terminal counts must match"))
    length(initial_phase_voltage) == nph ||
        throw(ArgumentError("generator-equivalent initial voltage count must match phases"))
    isempty(zero_branches) &&
        throw(ArgumentError("generator equivalent requires at least one zero-mode branch"))
    isempty(positive_branches) &&
        throw(ArgumentError("generator equivalent requires at least one positive-mode branch"))

    modal_branches = [zero_branches; positive_branches]
    rmfd = Float64[]
    sizehint!(rmfd, 5 * length(modal_branches))
    for branch in modal_branches
        append!(
            rmfd,
            (
                branch.resistance_ohm,
                branch.inductance_h,
                branch.capacitance_f,
                branch.damping_resistance_ohm,
                0.0,
            ),
        )
    end
    history_length =
        3 * (length(zero_branches) + (nph - 1) * length(positive_branches))
    state = breqiv_branch_state_from_arrays(
        rmfd,
        zeros(Float64, history_length),
        [length(zero_branches), length(positive_branches)],
        nph,
        0,
        0,
        0,
    )
    return BreqivHistoryInjection(
        from_nodes,
        to_nodes,
        state,
        Float64.(initial_phase_voltage),
        1.0;
        history_voltage_scale = 1.0,
        frequency_history_voltage_scale = -1.0,
    )
end

function _breqiv_frequency_domain_admittance(
    record::BreqivBranchRecord,
    angular_frequency::Float64,
)
    angular_frequency > 0.0 || throw(ArgumentError("angular frequency must be positive"))
    inductive_reactance = record.l * angular_frequency
    damping_resistance = record.rl
    real_impedance = record.r
    imaginary_impedance = inductive_reactance
    if damping_resistance != 0.0 && inductive_reactance != 0.0
        denominator = inv(damping_resistance^2 + inductive_reactance^2)
        real_impedance += damping_resistance * inductive_reactance^2 * denominator
        imaginary_impedance = inductive_reactance * damping_resistance^2 * denominator
    end
    capacitive_reactance =
        record.c > 0.0 ? inv(record.c * angular_frequency) : 0.0
    imaginary_impedance -= capacitive_reactance
    admittance = inv(complex(real_impedance, imaginary_impedance))
    return admittance, capacitive_reactance
end

function _breqiv_complex_modal_transform(phase_voltage::AbstractVector{<:Complex})
    nph = length(phase_voltage)
    nph > 0 || throw(ArgumentError("phase voltage phasors must not be empty"))
    modal = Vector{ComplexF64}(undef, nph)
    modal[1] = sum(phase_voltage) / nph
    for phase in 2:nph
        modal[phase] = (phase_voltage[1] - phase_voltage[phase]) / nph
    end
    return modal
end

function seed_breqiv_frequency_histories!(
    element::BreqivHistoryInjection,
    terminal_voltage_phasors::AbstractVector{<:Complex},
    angular_frequency::Real,
)
    nph = length(element.a)
    length(terminal_voltage_phasors) == nph ||
        throw(ArgumentError("terminal voltage phasor count must match BREQIV phases"))
    omega = Float64(angular_frequency)
    omega > 0.0 && isfinite(omega) ||
        throw(ArgumentError("angular frequency must be positive and finite"))
    oriented_phasors =
        element.frequency_history_voltage_scale .* ComplexF64.(terminal_voltage_phasors)
    element.initial_phase_voltage .= real.(oriented_phasors)
    modal_phasors = _breqiv_complex_modal_transform(oriented_phasors)

    fill!(element.state.zero_history, 0.0)
    fill!(element.state.positive_history, 0.0)
    for branch in eachindex(element.state.zero)
        admittance, capacitive_reactance =
            _breqiv_frequency_domain_admittance(element.state.zero[branch], omega)
        current_phasor = modal_phasors[1] * admittance
        element.state.zero_history[2, branch] = real(current_phasor)
        element.state.zero_history[3, branch] =
            -imag(current_phasor) * capacitive_reactance
    end
    for branch in eachindex(element.state.positive)
        admittance, capacitive_reactance =
            _breqiv_frequency_domain_admittance(element.state.positive[branch], omega)
        for mode in axes(element.state.positive_history, 3)
            current_phasor = modal_phasors[mode + 1] * admittance
            element.state.positive_history[2, branch, mode] = real(current_phasor)
            element.state.positive_history[3, branch, mode] =
                -imag(current_phasor) * capacitive_reactance
        end
    end
    element.initialized = false
    return element
end

function breqiv_history_current_to_phase!(phase_current::AbstractVector{Float64}, state::BreqivBranchSetState)
    fill!(state.modal_current, 0.0)
    for branch in eachindex(state.zero)
        state.modal_current[1] += state.zero_history[1, branch]
    end
    for mode in 1:size(state.positive_history, 3)
        mode_current = 0.0
        for branch in axes(state.positive_history, 2)
            mode_current += state.positive_history[1, branch, mode]
        end
        state.modal_current[mode + 1] = mode_current
    end
    return modal_current_to_phase!(phase_current, state.modal_current)
end

trace_output_channel_count(b::BreqivHistoryInjection) = length(b.phase_current)

function trace_output_channel_names!(names::Vector{Symbol}, element_name::Symbol, b::BreqivHistoryInjection)
    if length(b.phase_current) == 1
        push!(names, Symbol(String(element_name), "_history_i_pu"))
    else
        for phase in eachindex(b.phase_current)
            push!(names, Symbol(String(element_name), "_phase", phase, "_history_i_pu"))
        end
    end
    return names
end

function trace_output_values!(
    output::AbstractMatrix{Float64},
    first_channel::Int,
    sample::Int,
    b::BreqivHistoryInjection,
    voltage::AbstractVector{Float64},
)
    for phase in eachindex(b.phase_current)
        output[first_channel + phase - 1, sample] = b.phase_current[phase]
    end
    return first_channel + length(b.phase_current)
end

function initialize_breqiv_history_injection!(b::BreqivHistoryInjection, dt::Float64)
    b.initialized && return b
    breqiv_branch_state!(b.state, b.initial_phase_voltage, dt / 2.0)
    breqiv_history_current_to_phase!(b.phase_current, b.state)
    update_breqiv_phase_admittance!(b.phase_admittance, b.state)
    b.initialized = true
    return b
end

function advance_breqiv_history_current!(
    b::BreqivHistoryInjection,
    voltage::AbstractVector{<:Real},
    dt::Float64;
    consumed_for_step::Bool = true,
)
    initialize_breqiv_history_injection!(b, dt)
    for phase in eachindex(b.phase_voltage)
        b.phase_voltage[phase] =
            b.history_voltage_scale * branch_voltage(voltage, b.a[phase], b.b[phase])
    end
    fdcinj_branch_state!(b.phase_current, b.state, b.phase_voltage)
    b.history_current_consumed_for_step = consumed_for_step
    return b
end

function update_breqiv_phase_admittance!(phase_admittance::AbstractMatrix{Float64}, state::BreqivBranchSetState)
    nph = length(state.modal)
    size(phase_admittance, 1) == nph && size(phase_admittance, 2) == nph ||
        throw(ArgumentError("phase_admittance size does not match BREQIV state"))

    zero_admittance = 0.0
    for record in state.zero
        zero_admittance += record.azr
    end
    positive_admittance = 0.0
    for record in state.positive
        positive_admittance += record.azr
    end

    off_diagonal = (zero_admittance - positive_admittance) / nph
    diagonal = (zero_admittance + (nph - 1) * positive_admittance) / nph
    for col in 1:nph
        for row in 1:nph
            phase_admittance[row, col] = row == col ? diagonal : off_diagonal
        end
    end
    return phase_admittance
end

function stamp_admittance_entry!(y::AbstractMatrix{<:Real}, a_row::Int, b_row::Int, a_col::Int, b_col::Int, value::Float64)
    if a_row != 0 && a_col != 0
        y[a_row, a_col] += value
    end
    if a_row != 0 && b_col != 0
        y[a_row, b_col] -= value
    end
    if b_row != 0 && a_col != 0
        y[b_row, a_col] -= value
    end
    if b_row != 0 && b_col != 0
        y[b_row, b_col] += value
    end
    return nothing
end

function _coupled_inductive_conductance_matrix(
    b::CoupledInductiveBranch,
    dt::Float64,
)
    dt > 0.0 || throw(ArgumentError("dt must be positive"))
    inductive_conductance = (-0.5 * dt * b.angular_frequency) .* b.susceptance
    b.series_resistance == 0.0 && return inductive_conductance
    reference_conductance =
        inductive_conductance[b.resistance_reference_port, b.resistance_reference_port]
    reference_conductance > 0.0 ||
        throw(ArgumentError("coupled inductive resistance reference conductance must be positive"))
    return inductive_conductance ./
           (1.0 + b.series_resistance * reference_conductance)
end

function _coupled_inductive_conductance_matrix!(
    destination::Matrix{Float64},
    b::CoupledInductiveBranch,
    dt::Float64,
)
    dt > 0.0 || throw(ArgumentError("dt must be positive"))
    size(destination) == size(b.susceptance) ||
        throw(ArgumentError("coupled inductive conductance workspace size must match susceptance"))
    scale = -0.5 * dt * b.angular_frequency
    @inbounds @simd for index in eachindex(destination, b.susceptance)
        destination[index] = scale * b.susceptance[index]
    end
    if b.series_resistance != 0.0
        reference_conductance =
            destination[b.resistance_reference_port, b.resistance_reference_port]
        reference_conductance > 0.0 || throw(ArgumentError(
            "coupled inductive resistance reference conductance must be positive",
        ))
        denominator = 1.0 + b.series_resistance * reference_conductance
        @inbounds @simd for index in eachindex(destination)
            destination[index] /= denominator
        end
    end
    return destination
end

function _coupled_inductive_history_current(
    b::CoupledInductiveBranch,
    conductance::AbstractMatrix{Float64},
    dt::Float64,
)
    if b.series_resistance == 0.0
        return b.previous_current .+ conductance * b.previous_voltage
    end
    inductive_reference_conductance =
        -0.5 * dt * b.angular_frequency *
        b.susceptance[b.resistance_reference_port, b.resistance_reference_port]
    resistance_ratio = b.series_resistance * inductive_reference_conductance
    history_scale = (1.0 - resistance_ratio) / (1.0 + resistance_ratio)
    return history_scale .* b.previous_current .+ conductance * b.previous_voltage
end

function _coupled_inductive_history_current!(
    destination::Vector{Float64},
    b::CoupledInductiveBranch,
    conductance::AbstractMatrix{Float64},
    dt::Float64,
)
    length(destination) == length(b.previous_current) || throw(ArgumentError(
        "coupled inductive history-current workspace size must match ports",
    ))
    _coupled_matrix_vector_mul!(destination, conductance, b.previous_voltage)
    if b.series_resistance == 0.0
        @inbounds @simd for index in eachindex(destination, b.previous_current)
            destination[index] = b.previous_current[index] + destination[index]
        end
        return destination
    end
    inductive_reference_conductance =
        -0.5 * dt * b.angular_frequency *
        b.susceptance[b.resistance_reference_port, b.resistance_reference_port]
    resistance_ratio = b.series_resistance * inductive_reference_conductance
    history_scale = (1.0 - resistance_ratio) / (1.0 + resistance_ratio)
    @inbounds @simd for index in eachindex(destination, b.previous_current)
        destination[index] =
            history_scale * b.previous_current[index] + destination[index]
    end
    return destination
end

function _coupled_series_rl_conductance_matrix!(
    destination::Matrix{Float64},
    b::CoupledSeriesRLBranch,
    dt::Float64,
)
    dt > 0.0 || throw(ArgumentError("dt must be positive"))
    size(destination) == size(b.resistance_matrix) ||
        throw(ArgumentError("coupled series R-L conductance workspace size must match ports"))
    b.cached_dt_s == dt && return destination
    companion_impedance =
        b.resistance_matrix .+ (2.0 / dt) .* b.inductance_matrix
    factor = cholesky(Symmetric(companion_impedance))
    destination .= factor \ I
    b.cached_dt_s = dt
    return destination
end

function _coupled_series_rl_history_current!(
    destination::Vector{Float64},
    b::CoupledSeriesRLBranch,
    conductance::AbstractMatrix{Float64},
    dt::Float64,
)
    dt > 0.0 || throw(ArgumentError("dt must be positive"))
    length(destination) == length(b.previous_current) ||
        throw(ArgumentError("coupled series R-L history workspace size must match ports"))
    companion_history =
        b.previous_voltage .+
        ((2.0 / dt) .* b.inductance_matrix .- b.resistance_matrix) *
        b.previous_current
    _coupled_matrix_vector_mul!(destination, conductance, companion_history)
    return destination
end

function _coupled_branch_port_voltages(
    voltage::AbstractVector{Float64},
    a::AbstractVector{Int},
    b::AbstractVector{Int},
)
    return Float64[branch_voltage(voltage, a[index], b[index]) for index in eachindex(a)]
end

function _coupled_branch_port_voltages!(
    destination::Vector{Float64},
    voltage::AbstractVector{Float64},
    a::AbstractVector{Int},
    b::AbstractVector{Int},
)
    length(destination) == length(a) == length(b) || throw(ArgumentError(
        "coupled inductive port-voltage workspace size must match terminals",
    ))
    @inbounds @simd for index in eachindex(destination, a, b)
        destination[index] = branch_voltage(voltage, a[index], b[index])
    end
    return destination
end

function _coupled_matrix_vector_mul!(
    destination::Vector{Float64},
    matrix::Matrix{Float64},
    source::Vector{Float64},
)
    size(matrix, 1) == length(destination) || throw(DimensionMismatch(
        "coupled matrix row count must match the destination length",
    ))
    size(matrix, 2) == length(source) || throw(DimensionMismatch(
        "coupled matrix column count must match the source length",
    ))
    @inbounds for row in axes(matrix, 1)
        value = 0.0
        for column in axes(matrix, 2)
            value += matrix[row, column] * source[column]
        end
        destination[row] = value
    end
    return destination
end

function stamp_history_current!(
    rhs::AbstractVector{<:Real},
    b::CoupledInductiveBranch,
    dt::Float64,
)
    conductance = _coupled_inductive_conductance_matrix!(
        b.conductance_workspace,
        b,
        dt,
    )
    history_current = _coupled_inductive_history_current!(
        b.history_current_workspace,
        b,
        conductance,
        dt,
    )
    for index in eachindex(history_current)
        stamp_history_current!(rhs, b.a[index], b.b[index], history_current[index])
    end
    return rhs
end

function stamp_history_current!(
    rhs::AbstractVector{<:Real},
    b::CoupledSeriesRLBranch,
    dt::Float64,
)
    conductance = _coupled_series_rl_conductance_matrix!(
        b.conductance_workspace,
        b,
        dt,
    )
    history_current = _coupled_series_rl_history_current!(
        b.history_current_workspace,
        b,
        conductance,
        dt,
    )
    @inbounds for index in eachindex(history_current)
        stamp_history_current!(rhs, b.a[index], b.b[index], history_current[index])
    end
    return rhs
end

function stamp_breqiv_phase_admittance!(y::AbstractMatrix{<:Real}, b::BreqivHistoryInjection)
    for col in eachindex(b.a)
        for row in eachindex(b.a)
            stamp_admittance_entry!(
                y,
                b.a[row],
                b.b[row],
                b.a[col],
                b.b[col],
                b.phase_admittance[row, col],
            )
        end
    end
    return nothing
end

function companion(b::SeriesRLBranch, dt::Float64)
    isfinite(dt) && dt > 0.0 ||
        throw(ArgumentError("series R-L timestep must be finite and positive"))
    isfinite(b.r) ||
        throw(ArgumentError("series R-L resistance must be finite"))
    isfinite(b.l) && b.l >= 0.0 ||
        throw(ArgumentError("series R-L inductance must be finite and nonnegative"))
    b.l > 0.0 || b.r > 0.0 ||
        throw(ArgumentError("series resistor-only branch resistance must be positive"))
    (isfinite(b.i_prev) && isfinite(b.v_prev) && isfinite(b.i_last)) ||
        throw(ArgumentError("series R-L accepted state must be finite"))
    denominator = b.r + 2.0 * b.l / dt
    isfinite(denominator) && denominator != 0.0 ||
        throw(ArgumentError("series R-L companion denominator must be finite and nonzero"))
    conductance, history_current =
        series_rl_companion(b.r, b.l, Float64(b.i_prev), b.v_prev, dt)
    all(isfinite, (conductance, history_current)) ||
        throw(ArgumentError("series R-L companion values must be finite"))
    return conductance, history_current
end

function companion(b::SeriesRLCBranch, dt::Float64)
    return series_rlc_companion(
        b.r,
        b.l,
        b.c,
        b.i_prev,
        b.inductor_voltage_prev,
        b.capacitor_voltage_prev,
        dt,
    )
end

function companion(b::CapacitorBranch, dt::Float64)
    return capacitor_companion(b.c, b.i_prev, b.v_prev, dt)
end

companion(b::SeriesRLBranch, dt::Float64, ::Val{TrapezoidalCompanion}) =
    companion(b, dt)

function companion(b::SeriesRLBranch, dt::Float64, ::Val{BackwardEulerCompanion})
    isfinite(dt) && dt > 0.0 || throw(ArgumentError(
        "series R-L backward-Euler timestep must be finite and positive",
    ))
    isfinite(b.r) || throw(ArgumentError("series R-L resistance must be finite"))
    isfinite(b.l) && b.l >= 0.0 || throw(ArgumentError(
        "series R-L inductance must be finite and nonnegative",
    ))
    b.l > 0.0 || b.r > 0.0 || throw(ArgumentError(
        "series resistor-only branch resistance must be positive",
    ))
    all(isfinite, (b.i_prev, b.v_prev, b.i_last)) || throw(ArgumentError(
        "series R-L accepted state must be finite",
    ))
    conductance, history_current = series_rl_backward_euler_companion(
        b.r,
        b.l,
        Float64(b.i_prev),
        dt,
    )
    all(isfinite, (conductance, history_current)) || throw(ArgumentError(
        "series R-L backward-Euler companion values must be finite",
    ))
    return conductance, history_current
end

companion(b::SeriesRLCBranch, dt::Float64, ::Val{TrapezoidalCompanion}) =
    companion(b, dt)

function companion(b::SeriesRLCBranch, dt::Float64, ::Val{BackwardEulerCompanion})
    isfinite(dt) && dt > 0.0 || throw(ArgumentError(
        "series RLC backward-Euler timestep must be finite and positive",
    ))
    conductance, history_current = series_rlc_backward_euler_companion(
        b.r,
        b.l,
        b.c,
        b.i_prev,
        b.capacitor_voltage_prev,
        dt,
    )
    all(isfinite, (conductance, history_current)) || throw(ArgumentError(
        "series RLC backward-Euler companion values must be finite",
    ))
    return conductance, history_current
end

companion(b::CapacitorBranch, dt::Float64, ::Val{TrapezoidalCompanion}) =
    companion(b, dt)

function companion(b::CapacitorBranch, dt::Float64, ::Val{BackwardEulerCompanion})
    isfinite(dt) && dt > 0.0 || throw(ArgumentError(
        "capacitor backward-Euler timestep must be finite and positive",
    ))
    isfinite(b.c) && b.c > 0.0 || throw(ArgumentError(
        "capacitor capacitance must be finite and positive",
    ))
    all(isfinite, (b.i_prev, b.v_prev, b.i_last)) || throw(ArgumentError(
        "capacitor accepted state must be finite",
    ))
    conductance, history_current = capacitor_backward_euler_companion(
        b.c,
        b.v_prev,
        dt,
    )
    all(isfinite, (conductance, history_current)) || throw(ArgumentError(
        "capacitor backward-Euler companion values must be finite",
    ))
    return conductance, history_current
end

function stamp!(y, rhs, b::ConductanceBranch, t::Float64, dt::Float64)
    stamp_conductance!(y, b.a, b.b, b.g)
end

function stamp!(
    y,
    _rhs,
    constraint::IdealTransformerVoltageConstraint,
    _t::Float64,
    _dt::Float64,
)
    auxiliary = constraint.constraint_node
    for (node, coefficient) in
        zip(constraint.terminal_nodes, constraint.terminal_coefficients)
        value = -coefficient
        y[node, auxiliary] += value
        y[auxiliary, node] += value
    end
    return nothing
end

function stamp!(y, rhs, b::SeriesRLBranch, t::Float64, dt::Float64)
    g, ih = companion(b, dt)
    stamp_conductance!(y, b.a, b.b, g)
    stamp_history_current!(rhs, b.a, b.b, ih)
end

function stamp!(y, rhs, b::SeriesRLCBranch, t::Float64, dt::Float64)
    g, ih = companion(b, dt)
    stamp_conductance!(y, b.a, b.b, g)
    stamp_history_current!(rhs, b.a, b.b, ih)
end

function stamp!(y, rhs, b::CoupledInductiveBranch, t::Float64, dt::Float64)
    conductance = _coupled_inductive_conductance_matrix!(
        b.conductance_workspace,
        b,
        dt,
    )
    for col in eachindex(b.a)
        for row in eachindex(b.a)
            stamp_admittance_entry!(
                y,
                b.a[row],
                b.b[row],
                b.a[col],
                b.b[col],
                conductance[row, col],
            )
        end
    end
    history_current = _coupled_inductive_history_current!(
        b.history_current_workspace,
        b,
        conductance,
        dt,
    )
    @inbounds for index in eachindex(history_current)
        stamp_history_current!(rhs, b.a[index], b.b[index], history_current[index])
    end
    return nothing
end

function stamp!(y, rhs, b::CoupledSeriesRLBranch, t::Float64, dt::Float64)
    conductance = _coupled_series_rl_conductance_matrix!(
        b.conductance_workspace,
        b,
        dt,
    )
    for col in eachindex(b.a), row in eachindex(b.a)
        stamp_admittance_entry!(
            y,
            b.a[row],
            b.b[row],
            b.a[col],
            b.b[col],
            conductance[row, col],
        )
    end
    history_current = _coupled_series_rl_history_current!(
        b.history_current_workspace,
        b,
        conductance,
        dt,
    )
    for index in eachindex(history_current)
        stamp_history_current!(rhs, b.a[index], b.b[index], history_current[index])
    end
    return nothing
end

function stamp!(y, rhs, b::CapacitorBranch, t::Float64, dt::Float64)
    g, ih = companion(b, dt)
    stamp_conductance!(y, b.a, b.b, g)
    stamp_history_current!(rhs, b.a, b.b, ih)
end

function stamp!(
    y,
    rhs,
    element::EMTElement,
    t::Float64,
    dt::Float64,
    ::Val{TrapezoidalCompanion},
)
    return stamp!(y, rhs, element, t, dt)
end

function stamp!(
    y,
    rhs,
    element::EMTElement,
    t::Float64,
    dt::Float64,
    ::Val{BackwardEulerCompanion},
)
    backward_euler_companion_supported(element) || throw(ArgumentError(
        "backward-Euler companion treatment is not declared for $(typeof(element))",
    ))
    return stamp!(y, rhs, element, t, dt)
end

function _stamp_scalar_companion!(
    y,
    rhs,
    branch::Union{SeriesRLBranch,SeriesRLCBranch,CapacitorBranch},
    t::Float64,
    dt::Float64,
    method,
)
    conductance, history_current = companion(branch, dt, method)
    stamp_conductance!(y, branch.a, branch.b, conductance)
    stamp_history_current!(rhs, branch.a, branch.b, history_current)
    return nothing
end

stamp!(
    y,
    rhs,
    branch::Union{SeriesRLBranch,SeriesRLCBranch,CapacitorBranch},
    t::Float64,
    dt::Float64,
    method::Val{TrapezoidalCompanion},
) = _stamp_scalar_companion!(y, rhs, branch, t, dt, method)

stamp!(
    y,
    rhs,
    branch::Union{SeriesRLBranch,SeriesRLCBranch,CapacitorBranch},
    t::Float64,
    dt::Float64,
    method::Val{BackwardEulerCompanion},
) = _stamp_scalar_companion!(y, rhs, branch, t, dt, method)

function stamp!(y, rhs, b::BreqivHistoryInjection, t::Float64, dt::Float64)
    initialize_breqiv_history_injection!(b, dt)
    stamp_breqiv_phase_admittance!(y, b)
    for phase in eachindex(b.phase_current)
        stamp_history_current!(
            rhs,
            b.a[phase],
            b.b[phase],
            b.history_current_scale * b.phase_current[phase],
        )
    end
    return nothing
end

function stamp!(y, rhs, s::TheveninSource, t::Float64, dt::Float64)
    stamp_conductance!(y, s.node, 0, s.g)
    rhs[s.node] += s.g * s.value(t)
end

function stamp!(y, rhs, s::TwoTerminalTheveninSource, t::Float64, dt::Float64)
    source_voltage = Float64(s.value(t))
    isfinite(source_voltage) || throw(ArgumentError(
        "two-terminal Thevenin source voltage must be finite",
    ))
    stamp_conductance!(y, s.a, s.b, s.g)
    stamp_history_current!(rhs, s.a, s.b, -s.g * source_voltage)
    return nothing
end

function stamp!(y, rhs, s::CurrentInjection, t::Float64, dt::Float64)
    rhs[s.node] += s.value(t)
end

continuous_time_invariant(::Union{ConductanceBranch,IdealTransformerVoltageConstraint,
    SeriesRLBranch,SeriesRLCBranch,CoupledInductiveBranch,CoupledSeriesRLBranch,CapacitorBranch}) = true

function _continuous_source_value(signal, time::ContinuousEvaluationTime)
    # An untyped numeric closure accepts metadata by dispatch but cannot do
    # arithmetic with it. Dedicated metadata methods stay authoritative;
    # otherwise evaluate numeric sources at the unrounded physical seconds.
    metadata_method = applicable(signal,time) ? which(signal,(ContinuousEvaluationTime,)) : nothing
    dedicated_metadata = metadata_method !== nothing &&
        Base.unwrap_unionall(metadata_method.sig).parameters[2] !== Any
    if dedicated_metadata
        value = signal(time)
    elseif applicable(signal,BigFloat(time.high_s))
        value = setprecision(BigFloat,106) do
            signal(continuous_time_seconds(time))
        end
    elseif iszero(time.low_s)
        value = signal(time.high_s)
    else
        throw(ArgumentError("source signal does not implement continuous-time evaluation"))
    end
    isfinite(value) || throw(ArgumentError("continuous source value must be finite"))
    value
end

function stamp!(y,rhs,source::TheveninSource,time::ContinuousEvaluationTime,step::Float64,method)
    stamp_conductance!(y,source.node,0,source.g)
    rhs[source.node] += source.g * _continuous_source_value(source.value,time)
    nothing
end

function stamp!(y,rhs,source::TwoTerminalTheveninSource,time::ContinuousEvaluationTime,step::Float64,method)
    stamp_conductance!(y,source.a,source.b,source.g)
    stamp_history_current!(rhs,source.a,source.b,-source.g * _continuous_source_value(source.value,time))
    nothing
end

function stamp!(y,rhs,source::CurrentInjection,time::ContinuousEvaluationTime,step::Float64,method)
    rhs[source.node] += _continuous_source_value(source.value,time)
    nothing
end

update!(b::ConductanceBranch, v, dt::Float64) = nothing
update!(::IdealTransformerVoltageConstraint, _v, _dt::Float64) = nothing
update!(s::TheveninSource, v, dt::Float64) = nothing
update!(s::TwoTerminalTheveninSource, v, dt::Float64) = nothing
update!(s::CurrentInjection, v, dt::Float64) = nothing

function update!(b::BreqivHistoryInjection, v, dt::Float64)
    if b.history_current_consumed_for_step
        b.history_current_consumed_for_step = false
        return nothing
    end
    initialize_breqiv_history_injection!(b, dt)
    for phase in eachindex(b.phase_voltage)
        b.phase_voltage[phase] =
            b.history_voltage_scale * branch_voltage(v, b.a[phase], b.b[phase])
    end
    fdcinj_branch_state!(b.phase_current, b.state, b.phase_voltage)
    return nothing
end

function update!(b::SeriesRLBranch, v, dt::Float64)
    g, ih = companion(b, dt)
    vb = branch_voltage(v, b.a, b.b)
    isfinite(vb) ||
        throw(ArgumentError("series R-L branch voltage must be finite"))
    i = g * vb + ih
    isfinite(i) ||
        throw(ArgumentError("series R-L accepted current must be finite"))
    b.v_prev = vb
    b.i_prev = i
    b.i_last = i
    return nothing
end

function update!(b::SeriesRLCBranch, v, dt::Float64)
    previous_current = b.i_prev
    previous_capacitor_voltage = b.capacitor_voltage_prev
    g, ih = companion(b, dt)
    branch_value = branch_voltage(v, b.a, b.b)
    current = g * branch_value + ih
    capacitor_voltage =
        previous_capacitor_voltage + dt / (2.0 * b.c) * (current + previous_current)
    b.inductor_voltage_prev = branch_value - b.r * current - capacitor_voltage
    b.capacitor_voltage_prev = capacitor_voltage
    b.v_prev = branch_value
    b.i_prev = current
    b.i_last = current
    return nothing
end

function update!(b::CoupledInductiveBranch, v, dt::Float64)
    conductance = _coupled_inductive_conductance_matrix!(
        b.conductance_workspace,
        b,
        dt,
    )
    history_current = _coupled_inductive_history_current!(
        b.history_current_workspace,
        b,
        conductance,
        dt,
    )
    voltage = _coupled_branch_port_voltages!(
        b.port_voltage_workspace,
        v,
        b.a,
        b.b,
    )
    _coupled_matrix_vector_mul!(b.last_current, conductance, voltage)
    @inbounds @simd for index in eachindex(b.last_current, history_current)
        b.last_current[index] += history_current[index]
    end
    copyto!(b.previous_voltage, voltage)
    copyto!(b.previous_current, b.last_current)
    return nothing
end

function update!(b::CoupledSeriesRLBranch, v, dt::Float64)
    conductance = _coupled_series_rl_conductance_matrix!(
        b.conductance_workspace,
        b,
        dt,
    )
    history_current = _coupled_series_rl_history_current!(
        b.history_current_workspace,
        b,
        conductance,
        dt,
    )
    voltage = _coupled_branch_port_voltages!(
        b.port_voltage_workspace,
        v,
        b.a,
        b.b,
    )
    _coupled_matrix_vector_mul!(b.last_current, conductance, voltage)
    b.last_current .+= history_current
    copyto!(b.previous_voltage, voltage)
    copyto!(b.previous_current, b.last_current)
    return nothing
end

function update!(b::CapacitorBranch, v, dt::Float64)
    g, ih = companion(b, dt)
    vb = branch_voltage(v, b.a, b.b)
    i = g * vb + ih
    b.v_prev = vb
    b.i_prev = i
    b.i_last = i
    return nothing
end

"""Recover terminal-subtraction bits only for a nodal-bound history high part."""
function _companion_history_high_voltage(branch::Union{CapacitorBranch,SeriesRLBranch}, incoming_voltage)
    # Explicit independently initialized history remains authoritative. When
    # its high part is the accepted nodal branch voltage, retain the terminal
    # subtraction's lost bits rather than treating the rounded difference as
    # an additional physical voltage change.
    incoming_voltage === nothing && return BigFloat(branch.v_prev)
    branch.v_prev == branch_voltage(incoming_voltage,branch.a,branch.b) ||
        return BigFloat(branch.v_prev)
    high(node) = node == 0 ? BigFloat(0) : BigFloat(incoming_voltage[node])
    high(branch.a)-high(branch.b)
end

_capacitor_history_high_voltage(branch::CapacitorBranch, incoming_voltage) =
    _companion_history_high_voltage(branch,incoming_voltage)

"""Evaluate R-L current as a change from its authoritative incoming current.

Use the stamped conductance, but do not round the Norton history before
combining it with the endpoint voltage. For an inductive branch this is
`i₀ + g*(v-r*i₀)` (BE) or `i₀ + g*(v+v₀-2r*i₀)` (trapezoidal).
The resistor-only limit has no incoming history.
"""
function _series_rl_companion_current(branch::SeriesRLBranch, conductance::Real,
    endpoint_voltage::Real, incoming_voltage::Real, ::Val{M},
) where {M}
    M in (TrapezoidalCompanion,BackwardEulerCompanion) || throw(ArgumentError(
        "unsupported compensated series R-L method"))
    g,voltage = BigFloat(conductance),BigFloat(endpoint_voltage)
    branch.l <= 0.0 && return g*voltage
    current,resistance = BigFloat(branch.i_prev),BigFloat(branch.r)
    drive = M === TrapezoidalCompanion ?
        voltage+BigFloat(incoming_voltage)-2resistance*current :
        voltage-resistance*current
    return current+g*drive
end

"""Accept the same R-L current used by compensated nodal KCL, before history rounding."""
function _series_rl_accepted_candidate(branch::SeriesRLBranch, voltage, endpoint_correction,
    incoming_correction, dt::Float64, method::Val; incoming_voltage=nothing,
)
    conductance,_ = companion(branch,dt,method)
    branch_voltage_v = branch_voltage(voltage,branch.a,branch.b)
    physical = setprecision(BigFloat,106) do
        high(node) = node == 0 ? BigFloat(0) : BigFloat(voltage[node])
        low(values,node) = node == 0 ? BigFloat(0) : BigFloat(values[node])
        endpoint = high(branch.a)+low(endpoint_correction,branch.a)-
            high(branch.b)-low(endpoint_correction,branch.b)
        previous = _companion_history_high_voltage(branch,incoming_voltage)+
            low(incoming_correction,branch.a)-low(incoming_correction,branch.b)
        current=_series_rl_companion_current(branch,conductance,endpoint,previous,method)
        (;current_a=current,endpoint_voltage_v=endpoint,incoming_voltage_v=previous)
    end
    current=physical.current_a
    all(isfinite,(branch_voltage_v,current)) || throw(ArgumentError(
        "accepted compensated series R-L voltage and current must be finite"))
    return (; voltage_v=branch_voltage_v,current_a=current,physical.endpoint_voltage_v,physical.incoming_voltage_v)
end
function update!(branch::SeriesRLBranch, voltage, endpoint_correction,
    incoming_correction, dt::Float64, method::Val; incoming_voltage=nothing,
)
    candidate = _series_rl_accepted_candidate(branch,voltage,endpoint_correction,
        incoming_correction,dt,method;incoming_voltage)
    branch.v_prev = candidate.voltage_v
    branch.i_prev = branch.i_last = candidate.current_a
    return nothing
end

# Exact-sum STORAGE only: never evaluates a branch law above kernel106.
function _three_wire_rl_current_storage(first::Real,second::Real)
    exact = -Rational{BigInt}(first)-Rational{BigInt}(second)
    bits = max(106,iszero(exact) ? 1 : ndigits(abs(numerator(exact));base=2))
    setprecision(BigFloat,bits) do
        value=BigFloat(numerator(exact))/BigFloat(denominator(exact))
        Rational{BigInt}(value)==exact || error("dependent R-L current storage lost bits")
        value
    end
end

function _three_wire_rl_joint_candidate(branches,candidates,current_bounds_a,step_s::Float64,method::Val;relative_flux_bound::Float64,relative_energy_bound::Float64)
    length(branches)==length(candidates)==length(current_bounds_a)==3 ||
        throw(DimensionMismatch("three-wire R-L joint acceptance requires three phases"))
    currents=[candidate.current_a for candidate in candidates]
    constrained=copy(currents)
    constrained[3]=_three_wire_rl_current_storage(currents[1],currents[2])
    errors=[abs(Rational{BigInt}(constrained[i])-Rational{BigInt}(currents[i])) for i in 1:3]
    all(i->errors[i]<=Rational{BigInt}(current_bounds_a[i]),1:3) ||
        throw(DomainError(errors,"joint R-L acceptance exceeded original companion/current bounds"))
    sum(Rational{BigInt}.(constrained))==0 || error("joint R-L acceptance violated exact KCL")
    # Original public phase reports, resistive power and stored energy must be
    # bitwise unchanged; physical bounds remain checked by ordinary integration.
    Float64.(constrained)==Float64.(currents) ||
        throw(DomainError(constrained,"joint R-L acceptance changed phase reporting"))
    stored(values)=sum(i->0.5*branches[i].l*Float64(values[i])^2,1:3)
    loss(values)=sum(i->branches[i].r*Float64(values[i])^2,1:3)
    stored(constrained)==stored(currents) && loss(constrained)==loss(currents) ||
        throw(DomainError(constrained,"joint R-L acceptance changed original energy/power report"))
    # Exact represented-value conservation checks are guards, not a wider
    # integration kernel or persistent Rational state. The staged candidate
    # above is always106; only dependent current STORAGE may exceed106.
    exact=Rational{BigInt}
    dt=exact(step_s)
    old=exact.([branch.i_prev for branch in branches]);accepted=exact.(constrained)
    r=exact.([branch.r for branch in branches]);l=exact.([branch.l for branch in branches])
    endpoint=exact.([candidate.endpoint_voltage_v for candidate in candidates])
    incoming=exact.([candidate.incoming_voltage_v for candidate in candidates])
    trap=method isa Val{TrapezoidalCompanion}
    method isa Union{Val{TrapezoidalCompanion},Val{BackwardEulerCompanion}} || throw(ArgumentError("unknown joint R-L method"))
    mean_current=trap ? (accepted.+old)./2 : accepted
    drive=trap ? (endpoint.+incoming.-r.*(accepted.+old))./2 : endpoint.-r.*accepted
    flux_change=l.*(accepted.-old);flux_integral=dt.*drive
    flux_scales=max.(abs.(l.*accepted),abs.(l.*old),abs.(flux_integral),exact(eps(Float64)))
    flux_errors=abs.(flux_change.-flux_integral)./flux_scales
    all(x->x<=exact(relative_flux_bound),flux_errors) || throw(DomainError(flux_errors,"joint R-L acceptance exceeded original three flux/companion bounds"))
    energy_change=sum(l.*(accepted.^2 .- old.^2))/2
    mean_voltage=trap ? (endpoint.+incoming)./2 : endpoint
    work=dt*sum(mean_voltage.*mean_current)
    loss=dt*sum(r.*mean_current.^2)
    numerical_loss=trap ? zero(exact) : sum(l.*(accepted.-old).^2)/2
    energy_defect=energy_change-work+loss+numerical_loss
    energy_scale=max(abs(work),abs(loss+numerical_loss),abs(energy_change),exact(eps(Float64)))
    abs(energy_defect)<=exact(relative_energy_bound)*energy_scale ||
        throw(DomainError(energy_defect,"joint R-L acceptance exceeded original energy bound"))
    return (;currents=constrained,companion_errors_a=errors,flux_errors,energy_defect,energy_scale)
end


"""Accept the capacitor current of a compensated nodal solve.

The incoming high/low pair belongs to the preceding accepted voltage. Use the
same rounded conductance as the stamp, but form its physical voltage change
before multiplying. Supplying the incoming high vector also retains terminal
subtraction bits; independently initialized component history stays authoritative.
Without it, the compatibility API uses the declared component history high part.
"""
function update!(branch::CapacitorBranch, voltage, endpoint_correction,
    incoming_correction, dt::Float64, method::Val{M},
    ; incoming_voltage=nothing,
) where {M}
    M in (TrapezoidalCompanion, BackwardEulerCompanion) || throw(ArgumentError(
        "unsupported compensated capacitor acceptance method"))
    conductance, _ = companion(branch, dt, method)
    branch_voltage_v = branch_voltage(voltage, branch.a, branch.b)
    current = setprecision(BigFloat, 106) do
        high(node) = node == 0 ? BigFloat(0) : BigFloat(voltage[node])
        low(values, node) = node == 0 ? BigFloat(0) : BigFloat(values[node])
        endpoint = high(branch.a) + low(endpoint_correction, branch.a) -
            high(branch.b) - low(endpoint_correction, branch.b)
        previous = _capacitor_history_high_voltage(branch,incoming_voltage) + low(incoming_correction, branch.a) -
            low(incoming_correction, branch.b)
        Float64(BigFloat(conductance) * (endpoint - previous) -
            (M === TrapezoidalCompanion ? BigFloat(branch.i_prev) : BigFloat(0)))
    end
    all(isfinite, (branch_voltage_v, current)) || throw(ArgumentError(
        "accepted compensated capacitor voltage and current must be finite"))
    branch.v_prev = branch_voltage_v
    branch.i_prev = current
    branch.i_last = current
    return nothing
end

function update!(
    element::EMTElement,
    voltage,
    dt::Float64,
    ::Val{TrapezoidalCompanion},
)
    return update!(element, voltage, dt)
end

function update!(
    element::EMTElement,
    voltage,
    dt::Float64,
    ::Val{BackwardEulerCompanion},
)
    backward_euler_companion_supported(element) || throw(ArgumentError(
        "backward-Euler companion acceptance is not declared for $(typeof(element))",
    ))
    return update!(element, voltage, dt)
end

function _update_scalar_companion!(
    branch::Union{SeriesRLBranch,CapacitorBranch},
    voltage,
    dt::Float64,
    method,
)
    conductance, history_current = companion(branch, dt, method)
    branch_voltage_v = branch_voltage(voltage, branch.a, branch.b)
    branch_current_a = conductance * branch_voltage_v + history_current
    all(isfinite, (branch_voltage_v, branch_current_a)) || throw(ArgumentError(
        "accepted backward-Euler branch voltage and current must be finite",
    ))
    branch.v_prev = branch_voltage_v
    branch.i_prev = branch_current_a
    branch.i_last = branch_current_a
    return nothing
end


update!(
    branch::Union{SeriesRLBranch,CapacitorBranch},
    voltage,
    dt::Float64,
    method::Val{BackwardEulerCompanion},
) = _update_scalar_companion!(branch, voltage, dt, method)

update!(
    branch::Union{SeriesRLBranch,CapacitorBranch},
    voltage,
    dt::Float64,
    ::Val{TrapezoidalCompanion},
) = update!(branch, voltage, dt)

function update!(
    branch::SeriesRLCBranch,
    voltage,
    dt::Float64,
    ::Val{BackwardEulerCompanion},
)
    conductance, history_current = companion(
        branch,
        dt,
        Val(BackwardEulerCompanion),
    )
    branch_voltage_v = branch_voltage(voltage, branch.a, branch.b)
    branch_current_a = conductance * branch_voltage_v + history_current
    capacitor_voltage_v = branch.capacitor_voltage_prev +
        dt * branch_current_a / branch.c
    inductor_voltage_v = branch_voltage_v -
        branch.r * branch_current_a - capacitor_voltage_v
    all(isfinite, (
        branch_voltage_v,
        branch_current_a,
        capacitor_voltage_v,
        inductor_voltage_v,
    )) || throw(ArgumentError(
        "accepted backward-Euler series RLC state must be finite",
    ))
    branch.inductor_voltage_prev = inductor_voltage_v
    branch.capacitor_voltage_prev = capacitor_voltage_v
    branch.v_prev = branch_voltage_v
    branch.i_prev = branch_current_a
    branch.i_last = branch_current_a
    return nothing
end

update!(
    branch::SeriesRLCBranch,
    voltage,
    dt::Float64,
    ::Val{TrapezoidalCompanion},
) = update!(branch, voltage, dt)

end

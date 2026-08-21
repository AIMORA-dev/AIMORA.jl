module EMTPartitioning

using SHA
using ..EMTTaskPlatform: EMTLogicalTime

export CausalLinearReconstruction,
       CausalZeroOrderReconstruction,
       DirectCoupledExchange,
       EMTDeckPartitionCheckpoint,
       EMTDeckPartitionResult,
       EMTDeckRegion,
       EMTInterfacePort,
       EMTInterfacePortKind,
       EMTPartitionCheckpoint,
       EMTPartitionExchangeMethod,
       EMTPartitionExchangePolicy,
       EMTPartitionFailure,
       EMTPartitionInterpolation,
       EMTPartitionPlan,
       EMTPartitionRegion,
       EMTPartitionResult,
       EMTPartitionTolerancePolicy,
       IteratedWaveformExchange,
       NortonInterfacePort,
       PartitionedDeckEMTStudy,
       PassiveTwoRegionRLCStudy,
       ScatteringInterfacePort,
       TheveninInterfacePort,
       TravelingWaveInterfacePort,
       VoltageCurrentInterfacePort,
       emt_partition_plan,
       partition_checkpoint_signature_sha256,
       partition_plan_signature_sha256,
       partition_result,
       partition_study_signature_sha256

const _PARTITION_IDENTITY = r"^[A-Za-z][A-Za-z0-9_.:-]*$"
const _MAXIMUM_REGIONS = 64
const _MAXIMUM_LOCAL_RATES = 16
const _MAXIMUM_RATE_RATIO = 10_000

@enum EMTPartitionExchangeMethod::UInt8 begin
    DirectCoupledExchange = 0x01
    IteratedWaveformExchange = 0x02
end

@enum EMTPartitionInterpolation::UInt8 begin
    CausalZeroOrderReconstruction = 0x01
    CausalLinearReconstruction = 0x02
end

@enum EMTInterfacePortKind::UInt8 begin
    VoltageCurrentInterfacePort = 0x01
    NortonInterfacePort = 0x02
    TheveninInterfacePort = 0x03
    ScatteringInterfacePort = 0x04
    TravelingWaveInterfacePort = 0x05
end

function _partition_identity(value::AbstractString, owner::AbstractString)
    identity = String(value)
    occursin(_PARTITION_IDENTITY, identity) || throw(ArgumentError(
        "$owner must be portable nonempty text beginning with a letter",
    ))
    return identity
end

function _positive_finite(value::Real, owner::AbstractString)
    normalized = Float64(value)
    isfinite(normalized) && normalized > 0.0 || throw(ArgumentError(
        "$owner must be finite and positive",
    ))
    return normalized
end

function _nonnegative_finite(value::Real, owner::AbstractString)
    normalized = Float64(value)
    isfinite(normalized) && normalized >= 0.0 || throw(ArgumentError(
        "$owner must be finite and nonnegative",
    ))
    return normalized
end

function _logical_ratio(numerator::EMTLogicalTime, denominator::EMTLogicalTime)
    denominator > zero(denominator) || throw(ArgumentError(
        "partition logical-time denominator must be positive",
    ))
    wide_numerator = BigInt(numerator.numerator) * denominator.denominator
    wide_denominator = BigInt(numerator.denominator) * denominator.numerator
    wide_denominator > 0 || throw(ArgumentError(
        "partition logical-time ratio must be positive",
    ))
    iszero(rem(wide_numerator, wide_denominator)) || return nothing
    quotient = div(wide_numerator, wide_denominator)
    typemin(Int) <= quotient <= typemax(Int) || throw(OverflowError(
        "partition logical-time ratio exceeds Int",
    ))
    return Int(quotient)
end

"""One exact fixed-step regional owner; each physical model identity may occur in only one region."""
struct EMTPartitionRegion
    identity::String
    model_identities::Tuple
    local_step::EMTLogicalTime

    function EMTPartitionRegion(
        identity::AbstractString,
        model_identities,
        local_step::EMTLogicalTime,
    )
        name = _partition_identity(identity, "partition region identity")
        local_step > zero(local_step) || throw(ArgumentError(
            "partition region local step must be positive",
        ))
        models = Tuple(sort!(String[
            _partition_identity(String(model), "partition model identity")
            for model in model_identities
        ]))
        isempty(models) && throw(ArgumentError(
            "partition region must own at least one model",
        ))
        length(models) == length(unique(models)) || throw(ArgumentError(
            "partition region repeats a model identity",
        ))
        return new(name, models, local_step)
    end
end

"""One explicitly oriented interface; both regional currents are positive outward toward the interface."""
struct EMTInterfacePort
    identity::String
    kind::EMTInterfacePortKind
    positive_region::String
    negative_region::String
    positive_terminal::String
    negative_terminal::String
    voltage_unit::String
    current_unit::String
    voltage_base_v::Float64
    current_base_a::Float64
    reference_impedance_ohm::Float64

    function EMTInterfacePort(
        identity::AbstractString,
        kind::EMTInterfacePortKind,
        positive_region::AbstractString,
        negative_region::AbstractString,
        positive_terminal::AbstractString,
        negative_terminal::AbstractString;
        voltage_unit::AbstractString = "V",
        current_unit::AbstractString = "A",
        voltage_base_v::Real = 1.0,
        current_base_a::Real = 1.0,
        reference_impedance_ohm::Real = 1.0,
    )
        positive = _partition_identity(positive_region, "positive port region")
        negative = _partition_identity(negative_region, "negative port region")
        positive != negative || throw(ArgumentError(
            "partition interface must join two distinct regions",
        ))
        String(voltage_unit) == "V" || throw(ArgumentError(
            "partition voltage/current ports currently require SI volts",
        ))
        String(current_unit) == "A" || throw(ArgumentError(
            "partition voltage/current ports currently require SI amperes",
        ))
        return new(
            _partition_identity(identity, "partition interface identity"),
            kind,
            positive,
            negative,
            _partition_identity(positive_terminal, "positive port terminal"),
            _partition_identity(negative_terminal, "negative port terminal"),
            String(voltage_unit),
            String(current_unit),
            _positive_finite(voltage_base_v, "partition voltage base"),
            _positive_finite(current_base_a, "partition current base"),
            _positive_finite(
                reference_impedance_ohm,
                "partition reference impedance",
            ),
        )
    end
end

"""Independent physical and scaled acceptance budgets for one interface window."""
struct EMTPartitionTolerancePolicy
    voltage_absolute_v::Float64
    voltage_relative::Float64
    current_absolute_a::Float64
    current_relative::Float64
    kcl_absolute_a::Float64
    kcl_relative::Float64
    interface_energy_absolute_j::Float64
    interface_energy_relative::Float64

    function EMTPartitionTolerancePolicy(;
        voltage_absolute_v::Real = 1.0e-8,
        voltage_relative::Real = 1.0e-7,
        current_absolute_a::Real = 1.0e-9,
        current_relative::Real = 1.0e-7,
        kcl_absolute_a::Real = 1.0e-9,
        kcl_relative::Real = 1.0e-8,
        interface_energy_absolute_j::Real = 1.0e-9,
        interface_energy_relative::Real = 1.0e-6,
    )
        values = Float64.(tuple(
            voltage_absolute_v,
            voltage_relative,
            current_absolute_a,
            current_relative,
            kcl_absolute_a,
            kcl_relative,
            interface_energy_absolute_j,
            interface_energy_relative,
        ))
        all(isfinite, values) && all(>=(0.0), values) || throw(ArgumentError(
            "partition tolerances must be finite and nonnegative",
        ))
        values[1] > 0.0 && values[3] > 0.0 && values[5] > 0.0 &&
            values[7] > 0.0 || throw(ArgumentError(
                "partition absolute tolerances must be positive",
            ))
        return new(values...)
    end
end

"""A declared deterministic exchange algorithm and its bounded convergence policy."""
struct EMTPartitionExchangePolicy
    method::EMTPartitionExchangeMethod
    interpolation::EMTPartitionInterpolation
    maximum_iterations::Int
    relaxation::Float64
    tolerances::EMTPartitionTolerancePolicy

    function EMTPartitionExchangePolicy(
        method::EMTPartitionExchangeMethod = IteratedWaveformExchange;
        interpolation::EMTPartitionInterpolation = CausalLinearReconstruction,
        maximum_iterations::Integer = 32,
        relaxation::Real = 1.0,
        tolerances::EMTPartitionTolerancePolicy = EMTPartitionTolerancePolicy(),
    )
        iteration_count = Int(maximum_iterations)
        1 <= iteration_count <= 10_000 || throw(ArgumentError(
            "partition maximum iterations must be from 1 through 10,000",
        ))
        damping = Float64(relaxation)
        isfinite(damping) && 0.0 < damping <= 1.0 || throw(ArgumentError(
            "partition relaxation must be in (0, 1]",
        ))
        return new(method, interpolation, iteration_count, damping, tolerances)
    end
end

struct EMTPartitionPlan{R<:Tuple,P<:Tuple}
    start::EMTLogicalTime
    stop::EMTLogicalTime
    communication_step::EMTLogicalTime
    regions::R
    ports::P
    exchange::EMTPartitionExchangePolicy
    rate_ratios::Tuple
    communication_window_count::Int
    signature_sha256::String
end

function _partition_signature_lines(
    start,
    stop,
    communication_step,
    regions,
    ports,
    exchange,
    ratios,
)
    lines = String[
        "schema=aimora.emt.partition_plan.v1",
        "start=$(start.numerator)/$(start.denominator)",
        "stop=$(stop.numerator)/$(stop.denominator)",
        "communication_step=$(communication_step.numerator)/$(communication_step.denominator)",
        "method=$(UInt8(exchange.method))",
        "interpolation=$(UInt8(exchange.interpolation))",
        "maximum_iterations=$(exchange.maximum_iterations)",
        "relaxation=$(repr(exchange.relaxation))",
    ]
    tolerance = exchange.tolerances
    for field in fieldnames(EMTPartitionTolerancePolicy)
        push!(lines, "tolerance.$field=$(repr(getfield(tolerance, field)))")
    end
    for (region, ratio) in zip(regions, ratios)
        push!(lines, "region=$(region.identity)|$(join(region.model_identities, ','))|$(region.local_step.numerator)/$(region.local_step.denominator)|$ratio")
    end
    for port in ports
        push!(lines, "port=$(port.identity)|$(UInt8(port.kind))|$(port.positive_region)|$(port.negative_region)|$(port.positive_terminal)|$(port.negative_terminal)|$(port.voltage_unit)|$(port.current_unit)|$(repr(port.voltage_base_v))|$(repr(port.current_base_a))|$(repr(port.reference_impedance_ohm))")
    end
    return lines
end

"""Solver-free declaration of one canonical EMT deck owned by one partition region."""
struct EMTDeckRegion
    identity::String
    deck_lines::Tuple
    source_identity::String
    initial_voltage_source::Symbol
    saturated_transformer_runtime::Bool
    coupled_lumped_history_runtime::Bool
    distributed_line_runtime::Bool
    signature_sha256::String

    function EMTDeckRegion(
        identity::AbstractString,
        deck_lines;
        source_identity::AbstractString = identity,
        initial_voltage_source::Symbol = :none,
        saturated_transformer_runtime::Bool = false,
        coupled_lumped_history_runtime::Bool = false,
        distributed_line_runtime::Bool = true,
    )
        region_identity = _partition_identity(identity, "deck region identity")
        source = _partition_identity(
            source_identity,
            "deck region source identity",
        )
        initial_voltage_source in (:none, :steady_state) || throw(ArgumentError(
            "deck region initial voltage source must be :none or :steady_state",
        ))
        lines = Tuple(String(line) for line in deck_lines)
        isempty(lines) && throw(ArgumentError(
            "deck region must contain at least one deck line",
        ))
        all(line -> !occursin('\0', line), lines) || throw(ArgumentError(
            "deck region lines must not contain NUL bytes",
        ))
        signature_lines = String[
            "schema=aimora.emt.deck_region.v1",
            "identity=$region_identity",
            "source_identity=$source",
            "initial_voltage_source=$initial_voltage_source",
            "saturated_transformer_runtime=$saturated_transformer_runtime",
            "coupled_lumped_history_runtime=$coupled_lumped_history_runtime",
            "distributed_line_runtime=$distributed_line_runtime",
        ]
        for (index, line) in enumerate(lines)
            push!(
                signature_lines,
                "line=$index|$(ncodeunits(line))|$(bytes2hex(sha256(line)))",
            )
        end
        return new(
            region_identity,
            lines,
            source,
            initial_voltage_source,
            saturated_transformer_runtime,
            coupled_lumped_history_runtime,
            distributed_line_runtime,
            bytes2hex(sha256(join(signature_lines, '\n'))),
        )
    end
end

"""A partition plan bound to complete canonical regional decks and initial port currents."""
struct PartitionedDeckEMTStudy{P<:EMTPartitionPlan,R<:Tuple}
    plan::P
    regions::R
    initial_interface_current_a::Tuple
    signature_sha256::String

    function PartitionedDeckEMTStudy(
        plan::P,
        regions;
        initial_interface_current_a = zeros(length(plan.ports)),
    ) where {P<:EMTPartitionPlan}
        declarations = Tuple(regions)
        all(region -> region isa EMTDeckRegion, declarations) || throw(
            ArgumentError("partitioned deck regions must use EMTDeckRegion"),
        )
        declared_names = getfield.(declarations, :identity)
        length(declared_names) == length(unique(declared_names)) || throw(
            ArgumentError("partitioned deck study repeats a region identity"),
        )
        plan_names = getfield.(plan.regions, :identity)
        Set(declared_names) == Set(plan_names) || throw(ArgumentError(
            "partitioned deck declarations must match every plan region exactly",
        ))
        ordered = Tuple(
            declarations[only(findall(==(name), declared_names))]
            for name in plan_names
        )
        currents = Tuple(Float64.(initial_interface_current_a))
        length(currents) == length(plan.ports) || throw(DimensionMismatch(
            "partitioned deck initial currents must contain one value per interface",
        ))
        all(isfinite, currents) || throw(ArgumentError(
            "partitioned deck initial interface currents must be finite",
        ))
        signature_lines = String[
            "schema=aimora.emt.partitioned_deck_study.v1",
            "plan=$(plan.signature_sha256)",
        ]
        for region in ordered
            push!(signature_lines, "region=$(region.identity)|$(region.signature_sha256)")
        end
        for (port, current) in zip(plan.ports, currents)
            push!(signature_lines, "initial_current=$(port.identity)|$(repr(current))")
        end
        return new{P,typeof(ordered)}(
            plan,
            ordered,
            currents,
            bytes2hex(sha256(join(signature_lines, '\n'))),
        )
    end
end

partition_study_signature_sha256(study::PartitionedDeckEMTStudy) =
    study.signature_sha256

"""Accepted synchronization samples and conservative interface diagnostics."""
struct EMTDeckPartitionResult
    plan_signature_sha256::String
    study_signature_sha256::String
    region_identities::Tuple
    port_identities::Tuple
    time_s::Vector{Float64}
    positive_terminal_voltage_v::Matrix{Float64}
    negative_terminal_voltage_v::Matrix{Float64}
    interface_current_a::Matrix{Float64}
    voltage_residual_v::Matrix{Float64}
    kcl_residual_a::Matrix{Float64}
    interface_energy_defect_j::Matrix{Float64}
    fixed_point_iterations::Vector{Int}
    regional_local_step_counts::Tuple
    accepted_window_count::Int
    rejected_window_count::Int
    accepted::Bool
    deterministic_signature_sha256::String
end

function emt_partition_plan(
    regions,
    ports;
    start::EMTLogicalTime,
    stop::EMTLogicalTime,
    communication_step::EMTLogicalTime,
    exchange::EMTPartitionExchangePolicy = EMTPartitionExchangePolicy(),
)
    start < stop || throw(ArgumentError(
        "partition stop must be later than start",
    ))
    communication_step > zero(communication_step) || throw(ArgumentError(
        "partition communication step must be positive",
    ))
    typed_regions = Tuple(regions)
    typed_ports = Tuple(ports)
    all(region -> region isa EMTPartitionRegion, typed_regions) || throw(ArgumentError(
        "partition regions must use EMTPartitionRegion",
    ))
    all(port -> port isa EMTInterfacePort, typed_ports) || throw(ArgumentError(
        "partition ports must use EMTInterfacePort",
    ))
    1 <= length(typed_regions) <= _MAXIMUM_REGIONS || throw(ArgumentError(
        "partition plan requires from 1 through 64 regions",
    ))
    region_names = getfield.(typed_regions, :identity)
    length(region_names) == length(unique(region_names)) || throw(ArgumentError(
        "partition plan repeats a region identity",
    ))
    model_names = reduce(vcat, (collect(region.model_identities) for region in typed_regions))
    length(model_names) == length(unique(model_names)) || throw(ArgumentError(
        "partition plan assigns one physical model to multiple regions",
    ))
    port_names = getfield.(typed_ports, :identity)
    length(port_names) == length(unique(port_names)) || throw(ArgumentError(
        "partition plan repeats an interface identity",
    ))
    for port in typed_ports
        port.positive_region in region_names && port.negative_region in region_names ||
            throw(ArgumentError("partition interface names an absent region"))
    end
    ratios = Tuple(map(typed_regions) do region
        ratio = _logical_ratio(communication_step, region.local_step)
        ratio === nothing && throw(ArgumentError(
            "partition local steps must divide the communication step exactly",
        ))
        1 <= ratio <= _MAXIMUM_RATE_RATIO || throw(ArgumentError(
            "partition local-to-communication rate ratio must be from 1 through 10,000",
        ))
        ratio
    end)
    length(unique(ratios)) <= _MAXIMUM_LOCAL_RATES || throw(ArgumentError(
        "partition plan exceeds 16 distinct local rates",
    ))
    window_count = _logical_ratio(stop - start, communication_step)
    window_count === nothing && throw(ArgumentError(
        "partition horizon must contain an exact integer number of communication windows",
    ))
    window_count > 0 || throw(ArgumentError(
        "partition horizon must contain at least one communication window",
    ))
    if exchange.method == DirectCoupledExchange
        length(unique(getfield.(typed_regions, :local_step))) == 1 || throw(ArgumentError(
            "direct coupled exchange requires one equal regional step",
        ))
    end
    signature = bytes2hex(sha256(join(
        _partition_signature_lines(
            start,
            stop,
            communication_step,
            typed_regions,
            typed_ports,
            exchange,
            ratios,
        ),
        '\n',
    )))
    return EMTPartitionPlan(
        start,
        stop,
        communication_step,
        typed_regions,
        typed_ports,
        exchange,
        ratios,
        window_count,
        signature,
    )
end

partition_plan_signature_sha256(plan::EMTPartitionPlan) = plan.signature_sha256

"""A synthetic passive source-RL and load-RC network split by one oriented voltage/current interface."""
struct PassiveTwoRegionRLCStudy{P<:EMTPartitionPlan}
    plan::P
    source_voltage_v::Float64
    source_resistance_ohm::Float64
    source_inductance_h::Float64
    load_resistance_ohm::Float64
    load_capacitance_f::Float64
    initial_source_current_a::Float64
    initial_interface_voltage_v::Float64
    signature_sha256::String

    function PassiveTwoRegionRLCStudy(
        plan::P;
        source_voltage_v::Real,
        source_resistance_ohm::Real,
        source_inductance_h::Real,
        load_resistance_ohm::Real,
        load_capacitance_f::Real,
        initial_source_current_a::Real = 0.0,
        initial_interface_voltage_v::Real = 0.0,
    ) where {P<:EMTPartitionPlan}
        length(plan.regions) == 2 || throw(ArgumentError(
            "passive two-region study requires exactly two regions",
        ))
        length(plan.ports) == 1 || throw(ArgumentError(
            "passive two-region study requires exactly one interface",
        ))
        only(plan.ports).kind == VoltageCurrentInterfacePort || throw(ArgumentError(
            "passive two-region study requires a voltage/current interface",
        ))
        source_voltage = Float64(source_voltage_v)
        source_resistance = _nonnegative_finite(
            source_resistance_ohm,
            "source resistance",
        )
        source_inductance = _positive_finite(source_inductance_h, "source inductance")
        load_resistance = _positive_finite(load_resistance_ohm, "load resistance")
        load_capacitance = _positive_finite(load_capacitance_f, "load capacitance")
        initial_current = Float64(initial_source_current_a)
        initial_voltage = Float64(initial_interface_voltage_v)
        all(isfinite, (source_voltage, initial_current, initial_voltage)) || throw(
            ArgumentError("passive two-region sources and initial state must be finite"),
        )
        lines = (
            "schema=aimora.emt.passive_two_region_rlc.v1",
            "plan=$(plan.signature_sha256)",
            "source_voltage_v=$(repr(source_voltage))",
            "source_resistance_ohm=$(repr(source_resistance))",
            "source_inductance_h=$(repr(source_inductance))",
            "load_resistance_ohm=$(repr(load_resistance))",
            "load_capacitance_f=$(repr(load_capacitance))",
            "initial_source_current_a=$(repr(initial_current))",
            "initial_interface_voltage_v=$(repr(initial_voltage))",
        )
        signature = bytes2hex(sha256(join(lines, '\n')))
        return new{P}(
            plan,
            source_voltage,
            source_resistance,
            source_inductance,
            load_resistance,
            load_capacitance,
            initial_current,
            initial_voltage,
            signature,
        )
    end
end

partition_study_signature_sha256(study::PassiveTwoRegionRLCStudy) =
    study.signature_sha256

struct EMTPartitionFailure <: Exception
    code::Symbol
    last_accepted_time_s::Float64
    message::String
end

function Base.showerror(io::IO, failure::EMTPartitionFailure)
    print(
        io,
        String(failure.code),
        ": ",
        failure.message,
        " [last_accepted_time_s=",
        failure.last_accepted_time_s,
        ']'
    )
end

struct EMTPartitionCheckpoint
    schema::Symbol
    plan_signature_sha256::String
    study_signature_sha256::String
    accepted_window_count::Int
    rejected_window_count::Int
    source_local_step_count::Int
    load_local_step_count::Int
    time_s::Float64
    source_current_a::Float64
    interface_voltage_v::Float64
    load_capacitor_current_a::Float64
    previous_voltage_slope_v_per_s::Float64
    time_trace_s::Tuple
    current_trace_a::Tuple
    voltage_trace_v::Tuple
    voltage_residual_trace_v::Tuple
    kcl_residual_trace_a::Tuple
    interface_energy_defect_trace_j::Tuple
    fixed_point_iteration_trace::Tuple
    signature_sha256::String
end

function partition_checkpoint_signature_sha256(checkpoint::EMTPartitionCheckpoint)
    return checkpoint.signature_sha256
end

"""One accepted multi-region synchronization point with portable regional state."""
struct EMTDeckPartitionCheckpoint
    schema::Symbol
    plan_signature_sha256::String
    study_signature_sha256::String
    region_identities::Tuple
    port_identities::Tuple
    accepted_window_count::Int
    rejected_window_count::Int
    regional_local_step_counts::Tuple
    time_s::Float64
    interface_current_a::Tuple
    positive_terminal_voltage_v::Tuple
    negative_terminal_voltage_v::Tuple
    time_trace_s::Tuple
    positive_terminal_voltage_trace_v::Tuple
    negative_terminal_voltage_trace_v::Tuple
    interface_current_trace_a::Tuple
    voltage_residual_trace_v::Tuple
    kcl_residual_trace_a::Tuple
    interface_energy_defect_trace_j::Tuple
    fixed_point_iteration_trace::Tuple
    regional_snapshots::Tuple
    last_failure::Union{Nothing,String}
    signature_sha256::String
end

function partition_checkpoint_signature_sha256(
    checkpoint::EMTDeckPartitionCheckpoint,
)
    return checkpoint.signature_sha256
end

struct EMTPartitionResult
    plan_signature_sha256::String
    study_signature_sha256::String
    time_s::Vector{Float64}
    source_current_a::Vector{Float64}
    interface_voltage_v::Vector{Float64}
    voltage_residual_v::Vector{Float64}
    kcl_residual_a::Vector{Float64}
    interface_energy_defect_j::Vector{Float64}
    fixed_point_iterations::Vector{Int}
    accepted_window_count::Int
    rejected_window_count::Int
    source_local_step_count::Int
    load_local_step_count::Int
    accepted::Bool
    deterministic_signature_sha256::String
end

function partition_result(
    plan_signature::AbstractString,
    study_signature::AbstractString,
    time_s,
    source_current_a,
    interface_voltage_v,
    voltage_residual_v,
    kcl_residual_a,
    interface_energy_defect_j,
    fixed_point_iterations;
    accepted_window_count::Integer,
    rejected_window_count::Integer,
    source_local_step_count::Integer,
    load_local_step_count::Integer,
    accepted::Bool,
)
    times = Float64[time_s...]
    currents = Float64[source_current_a...]
    voltages = Float64[interface_voltage_v...]
    voltage_residuals = Float64[voltage_residual_v...]
    kcl_residuals = Float64[kcl_residual_a...]
    energy_defects = Float64[interface_energy_defect_j...]
    iterations = Int[fixed_point_iterations...]
    sample_count = length(times)
    length(currents) == sample_count && length(voltages) == sample_count || throw(
        DimensionMismatch("partition state traces must have equal lengths"),
    )
    diagnostic_count = max(sample_count - 1, 0)
    all(length(values) == diagnostic_count for values in (
        voltage_residuals,
        kcl_residuals,
        energy_defects,
        iterations,
    )) || throw(DimensionMismatch(
        "partition window diagnostics must contain one row per accepted interval",
    ))
    all(isfinite, vcat(
        times,
        currents,
        voltages,
        voltage_residuals,
        kcl_residuals,
        energy_defects,
    )) || throw(ArgumentError("partition result traces must be finite"))
    lines = String[
        "schema=aimora.emt.partition_result.v1",
        "plan=$(String(plan_signature))",
        "study=$(String(study_signature))",
        "accepted_window_count=$(Int(accepted_window_count))",
        "rejected_window_count=$(Int(rejected_window_count))",
        "source_local_step_count=$(Int(source_local_step_count))",
        "load_local_step_count=$(Int(load_local_step_count))",
        "accepted=$accepted",
    ]
    for index in eachindex(times)
        push!(lines, "state=$index|$(repr(times[index]))|$(repr(currents[index]))|$(repr(voltages[index]))")
    end
    for index in eachindex(iterations)
        push!(lines, "window=$index|$(repr(voltage_residuals[index]))|$(repr(kcl_residuals[index]))|$(repr(energy_defects[index]))|$(iterations[index])")
    end
    signature = bytes2hex(sha256(join(lines, '\n')))
    return EMTPartitionResult(
        String(plan_signature),
        String(study_signature),
        times,
        currents,
        voltages,
        voltage_residuals,
        kcl_residuals,
        energy_defects,
        iterations,
        Int(accepted_window_count),
        Int(rejected_window_count),
        Int(source_local_step_count),
        Int(load_local_step_count),
        accepted,
        signature,
    )
end

end

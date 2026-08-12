module EMTTaskPlatform

using SHA

export AbstractEMTTask,
       CarrierEMTTask,
       ConverterControlEMTTask,
       EMTLogicalTime,
       EMTTaskCheckpoint,
       EMTTaskComputeStage,
       EMTTaskEffect,
       EMTTaskEnqueueStage,
       EMTTaskFamily,
       EMTTaskHoldStage,
       EMTTaskOccurrence,
       EMTTaskPlan,
       EMTTaskPlanEntry,
       EMTTaskPlatformFailure,
       EMTTaskReadStage,
       EMTTaskResult,
       EMTTaskSpec,
       EMTTaskStage,
       EMTTaskWriteStage,
       InterfaceEMTTask,
       InvalidateEMTInterface,
       InvalidateEMTOutput,
       InvalidateEMTPowerHistory,
       InvalidateEMTTopology,
       MechanicalEMTTask,
       ProtectionEMTTask,
       SourceEMTTask,
       ThermalEMTTask,
       UserDefinedEMTTask,
       compute_emt_task,
       emt_logical_time,
       emt_task_checkpoint,
       emt_task_name,
       emt_task_plan,
       emt_task_spec,
       emt_task_result_signature,
       read_emt_task,
       restore_emt_task_checkpoint!,
       write_emt_task!

const _PORTABLE_TASK_CODE = r"^[a-z][a-z0-9_]*$"
const _MAXIMUM_TASKS = 1_024
const _MAXIMUM_EDGES = 4_096
const _MAXIMUM_TASK_ACTIVATIONS = 1_000_000
const _MAXIMUM_HORIZON_TICKS = 1_000_000_000

@enum EMTTaskFamily::UInt8 begin
    ProtectionEMTTask = 0x01
    CarrierEMTTask = 0x02
    ConverterControlEMTTask = 0x03
    MechanicalEMTTask = 0x04
    SourceEMTTask = 0x05
    ThermalEMTTask = 0x06
    InterfaceEMTTask = 0x07
    UserDefinedEMTTask = 0x08
end

@enum EMTTaskEffect::UInt8 begin
    InvalidateEMTPowerHistory = 0x01
    InvalidateEMTTopology = 0x02
    InvalidateEMTInterface = 0x03
    InvalidateEMTOutput = 0x04
end

@enum EMTTaskStage::UInt8 begin
    EMTTaskReadStage = 0x01
    EMTTaskComputeStage = 0x02
    EMTTaskEnqueueStage = 0x03
    EMTTaskWriteStage = 0x04
    EMTTaskHoldStage = 0x05
end

"""A normalized rational SI second used without binary floating-point calendar arithmetic."""
struct EMTLogicalTime
    numerator::Int128
    denominator::Int128

    function EMTLogicalTime(numerator::Integer, denominator::Integer = 1)
        wide_numerator = BigInt(numerator)
        wide_denominator = BigInt(denominator)
        iszero(wide_denominator) && throw(ArgumentError(
            "EMT logical-time denominator must not be zero",
        ))
        if wide_denominator < 0
            wide_numerator = -wide_numerator
            wide_denominator = -wide_denominator
        end
        divisor = gcd(abs(wide_numerator), wide_denominator)
        wide_numerator = div(wide_numerator, divisor)
        wide_denominator = div(wide_denominator, divisor)
        typemin(Int128) <= wide_numerator <= typemax(Int128) ||
            throw(OverflowError("EMT logical-time numerator exceeds Int128"))
        wide_denominator <= typemax(Int128) ||
            throw(OverflowError("EMT logical-time denominator exceeds Int128"))
        return new(Int128(wide_numerator), Int128(wide_denominator))
    end
end

emt_logical_time(numerator::Integer, denominator::Integer = 1) =
    EMTLogicalTime(numerator, denominator)

function emt_logical_time(value::Rational)
    return EMTLogicalTime(numerator(value), denominator(value))
end

emt_logical_time(::AbstractFloat) = throw(ArgumentError(
    "EMT logical time must be supplied as an exact integer or rational value",
))

Base.zero(::Type{EMTLogicalTime}) = EMTLogicalTime(0)
Base.zero(::EMTLogicalTime) = EMTLogicalTime(0)
Base.iszero(value::EMTLogicalTime) = iszero(value.numerator)
Base.:(==)(left::EMTLogicalTime, right::EMTLogicalTime) =
    left.numerator == right.numerator && left.denominator == right.denominator
Base.isless(left::EMTLogicalTime, right::EMTLogicalTime) =
    BigInt(left.numerator) * right.denominator < BigInt(right.numerator) * left.denominator
Base.:<(left::EMTLogicalTime, right::EMTLogicalTime) = isless(left, right)
Base.:<=(left::EMTLogicalTime, right::EMTLogicalTime) = !(right < left)
Base.hash(value::EMTLogicalTime, seed::UInt) =
    hash((value.numerator, value.denominator), seed)
Base.Float64(value::EMTLogicalTime) = Float64(value.numerator) / Float64(value.denominator)
Base.show(io::IO, value::EMTLogicalTime) =
    print(io, value.numerator, '/', value.denominator, " s")

Base.:-(value::EMTLogicalTime) = EMTLogicalTime(-BigInt(value.numerator), value.denominator)
Base.:+(left::EMTLogicalTime, right::EMTLogicalTime) = EMTLogicalTime(
    BigInt(left.numerator) * right.denominator + BigInt(right.numerator) * left.denominator,
    BigInt(left.denominator) * right.denominator,
)
Base.:-(left::EMTLogicalTime, right::EMTLogicalTime) = left + (-right)
Base.:*(left::EMTLogicalTime, right::EMTLogicalTime) = EMTLogicalTime(
    BigInt(left.numerator) * right.numerator,
    BigInt(left.denominator) * right.denominator,
)
Base.:*(factor::Integer, value::EMTLogicalTime) = EMTLogicalTime(
    BigInt(factor) * value.numerator,
    value.denominator,
)
Base.:*(value::EMTLogicalTime, factor::Integer) = factor * value

"""A stable typed failure that retains the exact last accepted scheduler boundary."""
struct EMTTaskPlatformFailure <: Exception
    code::Symbol
    task::Union{Nothing,String}
    family::Union{Nothing,EMTTaskFamily}
    instant::Union{Nothing,EMTLogicalTime}
    stage::Union{Nothing,EMTTaskStage}
    last_accepted_instant::Union{Nothing,EMTLogicalTime}
    message::String

    function EMTTaskPlatformFailure(
        code::Symbol,
        message::AbstractString;
        task::Union{Nothing,AbstractString} = nothing,
        family::Union{Nothing,EMTTaskFamily} = nothing,
        instant::Union{Nothing,EMTLogicalTime} = nothing,
        stage::Union{Nothing,EMTTaskStage} = nothing,
        last_accepted_instant::Union{Nothing,EMTLogicalTime} = nothing,
    )
        occursin(_PORTABLE_TASK_CODE, String(code)) ||
            throw(ArgumentError("EMT task failure code is not portable"))
        normalized_message = String(message)
        isempty(normalized_message) &&
            throw(ArgumentError("EMT task failure message must not be empty"))
        return new(
            code,
            isnothing(task) ? nothing : String(task),
            family,
            instant,
            stage,
            last_accepted_instant,
            normalized_message,
        )
    end
end

function Base.showerror(io::IO, failure::EMTTaskPlatformFailure)
    print(io, String(failure.code), ": ", failure.message)
    isnothing(failure.task) || print(io, " [task=", failure.task, ']')
    isnothing(failure.instant) || print(io, " [instant=", failure.instant, ']')
end

_task_fail(code::Symbol, message::AbstractString; kwargs...) =
    throw(EMTTaskPlatformFailure(code, message; kwargs...))

function _owned_names(values, label::AbstractString)
    names = sort!(String[String(value) for value in values])
    any(name -> isempty(strip(name)) || occursin('\0', name), names) &&
        throw(ArgumentError("EMT task $label contains an empty or NUL identity"))
    length(names) == length(unique(names)) ||
        throw(ArgumentError("EMT task repeats a $label identity"))
    return names
end

"""One inert public task contract; the named owner retains all physical equations and state."""
struct EMTTaskSpec
    name::String
    family::EMTTaskFamily
    epoch::EMTLogicalTime
    period::EMTLogicalTime
    phase::EMTLogicalTime
    computational_delay::EMTLogicalTime
    priority::Int
    read_resources::Tuple
    write_resources::Tuple
    predecessors::Tuple
    effects::Tuple

    function EMTTaskSpec(
        name::AbstractString,
        family::EMTTaskFamily,
        epoch::EMTLogicalTime,
        period::EMTLogicalTime,
        phase::EMTLogicalTime,
        computational_delay::EMTLogicalTime;
        priority::Integer = 0,
        read_resources = String[],
        write_resources = String[],
        predecessors = String[],
        effects::AbstractVector{EMTTaskEffect} = EMTTaskEffect[],
    )
        identity = String(name)
        isempty(strip(identity)) || occursin('\0', identity) ?
            throw(ArgumentError("EMT task identity must be nonempty portable text")) : nothing
        typed_priority = try
            Int(priority)
        catch
            _task_fail(:invalid_task_priority, "EMT task priority exceeds Int"; task = identity)
        end
        minimum_period = EMTLogicalTime(1, 1_000_000_000)
        maximum_period = EMTLogicalTime(1_000)
        minimum_period <= period <= maximum_period || _task_fail(
            :task_period_out_of_domain,
            "EMT task period must be from 1 ns through 1,000 s";
            task = identity,
            family,
        )
        zero_time = EMTLogicalTime(0)
        zero_time <= phase < period || _task_fail(
            :invalid_task_phase,
            "EMT task phase must be in [0, period)";
            task = identity,
            family,
        )
        zero_time <= computational_delay <= 100 * period || _task_fail(
            :invalid_task_delay,
            "EMT task computational delay must be from zero through 100 periods";
            task = identity,
            family,
        )
        dependencies = _owned_names(predecessors, "predecessor")
        identity in dependencies && _task_fail(
            :task_self_dependency,
            "EMT task cannot declare itself as a predecessor";
            task = identity,
            family,
        )
        typed_effects = sort!(collect(effects); by = UInt8)
        length(typed_effects) == length(unique(typed_effects)) || _task_fail(
            :duplicate_task_effect,
            "EMT task repeats an invalidation effect";
            task = identity,
            family,
        )
        return new(
            identity,
            family,
            epoch,
            period,
            phase,
            computational_delay,
            typed_priority,
            Tuple(_owned_names(read_resources, "read resource")),
            Tuple(_owned_names(write_resources, "write resource")),
            Tuple(dependencies),
            Tuple(typed_effects),
        )
    end
end

"""One normalized immutable-by-contract task row consumed by the private dispatcher."""
struct EMTTaskPlanEntry
    spec::EMTTaskSpec
    registration_index::Int
    first_activation_tick::Int64
    period_ticks::Int64
    delay_ticks::Int64
    execution_rank::Int
end

"""A solver-free, exact, dependency-checked task plan for one bounded EMT horizon."""
struct EMTTaskPlan
    start::EMTLogicalTime
    stop::EMTLogicalTime
    quantum::EMTLogicalTime
    horizon_ticks::Int64
    entries::Vector{EMTTaskPlanEntry}
    execution_order::Vector{Int}
    signature_sha256::String
end

function _big_rational_values(values::Vector{EMTLogicalTime})
    denominator = foldl(
        lcm,
        (BigInt(value.denominator) for value in values);
        init = BigInt(1),
    )
    integers = BigInt[
        BigInt(value.numerator) * div(denominator, value.denominator) for value in values
    ]
    return denominator, integers
end

function _logical_quantum(values::Vector{EMTLogicalTime})
    denominator, integers = _big_rational_values(values)
    divisor = foldl(gcd, (abs(value) for value in integers); init = BigInt(0))
    iszero(divisor) && _task_fail(
        :empty_logical_calendar,
        "EMT task calendar has no positive logical-time quantity",
    )
    try
        return EMTLogicalTime(divisor, denominator)
    catch error
        error isa OverflowError || rethrow()
        _task_fail(:logical_calendar_overflow, sprint(showerror, error))
    end
end

function _logical_ticks(value::EMTLogicalTime, quantum::EMTLogicalTime, label::String)
    numerator = BigInt(value.numerator) * quantum.denominator
    denominator = BigInt(value.denominator) * quantum.numerator
    iszero(rem(numerator, denominator)) && denominator > 0 || _task_fail(
        :nonintegral_logical_boundary,
        "$label is not an integer multiple of the normalized logical quantum",
    )
    ticks = div(numerator, denominator)
    typemin(Int64) <= ticks <= typemax(Int64) || _task_fail(
        :logical_tick_overflow,
        "$label exceeds the exact Int64 tick domain",
    )
    return Int64(ticks)
end

function _task_graph(specs::Vector{EMTTaskSpec})
    names = getfield.(specs, :name)
    length(names) == length(unique(names)) || _task_fail(
        :duplicate_task_identity,
        "EMT task plan repeats a task identity",
    )
    indices = Dict(name => index for (index, name) in pairs(names))
    successors = [Int[] for _ in specs]
    indegree = zeros(Int, length(specs))
    edge_count = 0
    for (target, spec) in pairs(specs)
        for predecessor in spec.predecessors
            haskey(indices, predecessor) || _task_fail(
                :unknown_task_predecessor,
                "EMT task predecessor is absent from the plan";
                task = spec.name,
                family = spec.family,
            )
            source = indices[predecessor]
            push!(successors[source], target)
            indegree[target] += 1
            edge_count += 1
        end
    end
    edge_count <= _MAXIMUM_EDGES || _task_fail(
        :task_edge_limit_exceeded,
        "EMT task plan exceeds 4,096 explicit predecessor edges",
    )
    order = Int[]
    while length(order) < length(specs)
        selected = 0
        selected_key = nothing
        for index in eachindex(specs)
            indegree[index] == 0 || continue
            index in order && continue
            spec = specs[index]
            key = (spec.priority, UInt8(spec.family), spec.name, index)
            if selected == 0 || key < selected_key
                selected = index
                selected_key = key
            end
        end
        selected == 0 && _task_fail(
            :cyclic_task_dependencies,
            "EMT task predecessor graph contains a cycle",
        )
        push!(order, selected)
        for target in successors[selected]
            indegree[target] -= 1
        end
    end
    reachable = falses(length(specs), length(specs))
    for source in reverse(order)
        for target in successors[source]
            reachable[source, target] = true
            for downstream in eachindex(specs)
                reachable[target, downstream] && (reachable[source, downstream] = true)
            end
        end
    end
    return order, reachable
end

function _calendar_collides(
    left_start::Int64,
    left_period::Int64,
    right_start::Int64,
    right_period::Int64,
    horizon_ticks::Int64,
)
    left_start > horizon_ticks && return false
    right_start > horizon_ticks && return false
    common = gcd(left_period, right_period)
    delta = BigInt(right_start) - left_start
    iszero(rem(delta, common)) || return false
    reduced_right = div(BigInt(right_period), common)
    multiplier = if reduced_right == 1
        BigInt(0)
    else
        mod(
            div(delta, common) * invmod(mod(div(BigInt(left_period), common), reduced_right), reduced_right),
            reduced_right,
        )
    end
    collision = BigInt(left_start) + BigInt(left_period) * multiplier
    cycle = lcm(BigInt(left_period), BigInt(right_period))
    earliest = max(BigInt(left_start), BigInt(right_start))
    collision < earliest && (collision += cld(earliest - collision, cycle) * cycle)
    return collision <= horizon_ticks
end

function _resource_conflict(left::EMTTaskSpec, right::EMTTaskSpec)
    left_reads = Set(left.read_resources)
    left_writes = Set(left.write_resources)
    right_reads = Set(right.read_resources)
    right_writes = Set(right.write_resources)
    return !isempty(intersect(left_writes, union(right_reads, right_writes))) ||
        !isempty(intersect(right_writes, left_reads))
end

function _task_plan_signature(
    start::EMTLogicalTime,
    stop::EMTLogicalTime,
    quantum::EMTLogicalTime,
    entries,
)
    io = IOBuffer()
    println(io, "aimora-emt-task-plan-v1")
    println(io, start.numerator, '/', start.denominator)
    println(io, stop.numerator, '/', stop.denominator)
    println(io, quantum.numerator, '/', quantum.denominator)
    for entry in entries
        spec = entry.spec
        println(
            io,
            join((
                spec.name,
                UInt8(spec.family),
                spec.epoch.numerator,
                spec.epoch.denominator,
                spec.period.numerator,
                spec.period.denominator,
                spec.phase.numerator,
                spec.phase.denominator,
                spec.computational_delay.numerator,
                spec.computational_delay.denominator,
                spec.priority,
                join(spec.read_resources, ','),
                join(spec.write_resources, ','),
                join(spec.predecessors, ','),
                join(UInt8.(spec.effects), ','),
                entry.registration_index,
                entry.first_activation_tick,
                entry.period_ticks,
                entry.delay_ticks,
                entry.execution_rank,
            ), '\t'),
        )
    end
    return bytes2hex(sha256(take!(io)))
end

function emt_task_plan(
    declarations::AbstractVector{EMTTaskSpec};
    start::EMTLogicalTime,
    stop::EMTLogicalTime,
)
    specs = collect(declarations)
    isempty(specs) && _task_fail(:empty_task_plan, "EMT task plan requires declarations")
    length(specs) <= _MAXIMUM_TASKS || _task_fail(
        :task_limit_exceeded,
        "EMT task plan exceeds 1,024 task declarations",
    )
    start <= stop || _task_fail(
        :invalid_task_horizon,
        "EMT task plan stop must not precede its start",
    )
    quantum_values = EMTLogicalTime[stop - start]
    for spec in specs
        append!(
            quantum_values,
            (spec.period, spec.phase, spec.computational_delay, spec.epoch - start),
        )
    end
    quantum = _logical_quantum(quantum_values)
    horizon_ticks = _logical_ticks(stop - start, quantum, "EMT task horizon")
    0 <= horizon_ticks <= _MAXIMUM_HORIZON_TICKS || _task_fail(
        :task_horizon_limit_exceeded,
        "EMT task horizon exceeds 10^9 normalized logical ticks",
    )
    order, reachable = _task_graph(specs)
    rank = zeros(Int, length(specs))
    for (execution_rank, index) in pairs(order)
        rank[index] = execution_rank
    end
    entries = EMTTaskPlanEntry[]
    for (index, spec) in pairs(specs)
        period_ticks = _logical_ticks(spec.period, quantum, "task period")
        delay_ticks = _logical_ticks(spec.computational_delay, quantum, "task delay")
        initial_tick = _logical_ticks(spec.epoch + spec.phase - start, quantum, "task first activation")
        if initial_tick < 0
            initial_tick += cld(-initial_tick, period_ticks) * period_ticks
        end
        activation_count = initial_tick > horizon_ticks ? 0 :
            div(horizon_ticks - initial_tick, period_ticks) + 1
        activation_count <= _MAXIMUM_TASK_ACTIVATIONS || _task_fail(
            :task_activation_limit_exceeded,
            "EMT task plan exceeds 1,000,000 activations for one task";
            task = spec.name,
            family = spec.family,
        )
        push!(entries, EMTTaskPlanEntry(
            spec,
            index,
            initial_tick,
            period_ticks,
            delay_ticks,
            rank[index],
        ))
    end
    for left_index in eachindex(entries)
        left = entries[left_index]
        for right_index in (left_index + 1):length(entries)
            right = entries[right_index]
            _resource_conflict(left.spec, right.spec) || continue
            _calendar_collides(
                left.first_activation_tick,
                left.period_ticks,
                right.first_activation_tick,
                right.period_ticks,
                horizon_ticks,
            ) || continue
            (reachable[left_index, right_index] || reachable[right_index, left_index]) ||
                _task_fail(
                    :unordered_task_resource_conflict,
                    "same-instant EMT tasks with conflicting access require a predecessor path";
                    task = right.spec.name,
                    family = right.spec.family,
                )
        end
    end
    signature = _task_plan_signature(start, stop, quantum, entries)
    return EMTTaskPlan(
        start,
        stop,
        quantum,
        horizon_ticks,
        entries,
        order,
        signature,
    )
end

"""One accepted read/compute/enqueue/write/hold occurrence in exact and display time."""
struct EMTTaskOccurrence
    task::String
    family::EMTTaskFamily
    exact_instant::EMTLogicalTime
    instant_s::Float64
    stage::EMTTaskStage
    priority::Int
    activation_index::Int
    sample_index::Int
    release_index::Int
    execution_index::Int
end

function Base.:(==)(left::EMTTaskOccurrence, right::EMTTaskOccurrence)
    return left.task == right.task &&
        left.family == right.family &&
        left.exact_instant == right.exact_instant &&
        left.instant_s == right.instant_s &&
        left.stage == right.stage &&
        left.priority == right.priority &&
        left.activation_index == right.activation_index &&
        left.sample_index == right.sample_index &&
        left.release_index == right.release_index &&
        left.execution_index == right.execution_index
end

Base.isequal(left::EMTTaskOccurrence, right::EMTTaskOccurrence) =
    isequal(
        (
            left.task,
            left.family,
            left.exact_instant,
            left.instant_s,
            left.stage,
            left.priority,
            left.activation_index,
            left.sample_index,
            left.release_index,
            left.execution_index,
        ),
        (
            right.task,
            right.family,
            right.exact_instant,
            right.instant_s,
            right.stage,
            right.priority,
            right.activation_index,
            right.sample_index,
            right.release_index,
            right.execution_index,
        ),
    )

Base.hash(occurrence::EMTTaskOccurrence, seed::UInt) = hash(
    (
        occurrence.task,
        occurrence.family,
        occurrence.exact_instant,
        occurrence.instant_s,
        occurrence.stage,
        occurrence.priority,
        occurrence.activation_index,
        occurrence.sample_index,
        occurrence.release_index,
        occurrence.execution_index,
    ),
    seed,
)

"""A public checkpoint envelope whose concrete private state remains type-parameterized."""
struct EMTTaskCheckpoint{S}
    schema_version::VersionNumber
    plan_signature_sha256::String
    last_accepted_instant::EMTLogicalTime
    state::S
    state_signature_sha256::String

    function EMTTaskCheckpoint(
        plan_signature_sha256::AbstractString,
        last_accepted_instant::EMTLogicalTime,
        state::S,
        state_signature_sha256::AbstractString;
        schema_version::VersionNumber = v"1.0.0",
    ) where {S}
        plan_signature = String(plan_signature_sha256)
        state_signature = String(state_signature_sha256)
        all(signature -> occursin(r"^[0-9a-f]{64}$", signature), (plan_signature, state_signature)) ||
            throw(ArgumentError("EMT task checkpoint signatures must be lowercase SHA-256"))
        schema_version.major == 1 ||
            throw(ArgumentError("EMT task checkpoint major version is unsupported"))
        return new{S}(
            schema_version,
            plan_signature,
            last_accepted_instant,
            state,
            state_signature,
        )
    end
end

"""A typed accepted-or-failed public task-platform result."""
struct EMTTaskResult{C}
    accepted::Bool
    plan_signature_sha256::String
    occurrences::Vector{EMTTaskOccurrence}
    checkpoint::Union{Nothing,C}
    failure::Union{Nothing,EMTTaskPlatformFailure}
    task_counts::Vector{Pair{String,Int}}
    maximum_pending_depth::Int
    effects::Vector{EMTTaskEffect}
    deterministic_signature_sha256::String

    function EMTTaskResult(
        accepted::Bool,
        plan_signature_sha256::AbstractString,
        occurrences::AbstractVector{EMTTaskOccurrence};
        checkpoint = nothing,
        failure::Union{Nothing,EMTTaskPlatformFailure} = nothing,
        task_counts::AbstractVector{<:Pair{String,<:Integer}} = Pair{String,Int}[],
        maximum_pending_depth::Integer = 0,
        effects::AbstractVector{EMTTaskEffect} = EMTTaskEffect[],
        deterministic_signature_sha256::AbstractString,
    )
        accepted == isnothing(failure) || throw(ArgumentError(
            "accepted EMT task result and failure state disagree",
        ))
        accepted && isnothing(checkpoint) && throw(ArgumentError(
            "accepted EMT task result requires a final checkpoint",
        ))
        plan_signature = String(plan_signature_sha256)
        deterministic_signature = String(deterministic_signature_sha256)
        all(signature -> occursin(r"^[0-9a-f]{64}$", signature), (plan_signature, deterministic_signature)) ||
            throw(ArgumentError("EMT task result signatures must be lowercase SHA-256"))
        pending_depth = Int(maximum_pending_depth)
        pending_depth >= 0 || throw(ArgumentError(
            "EMT task maximum pending depth must be nonnegative",
        ))
        counts = Pair{String,Int}[
            String(pair.first) => Int(pair.second) for pair in task_counts
        ]
        all(pair -> pair.second >= 0, counts) || throw(ArgumentError(
            "EMT task result counts must be nonnegative",
        ))
        sort!(counts; by = first)
        length(first.(counts)) == length(unique(first.(counts))) || throw(ArgumentError(
            "EMT task result repeats a task count",
        ))
        typed_effects = sort!(unique(collect(effects)); by = UInt8)
        return new{typeof(checkpoint)}(
            accepted,
            plan_signature,
            collect(occurrences),
            checkpoint,
            failure,
            counts,
            pending_depth,
            typed_effects,
            deterministic_signature,
        )
    end
end

function emt_task_result_signature(plan::EMTTaskPlan, occurrences, checkpoint_signature::AbstractString)
    context = SHA.SHA2_256_CTX()
    function update_line(values...)
        line = join(values, '\t') * '\n'
        SHA.update!(context, codeunits(line))
        return nothing
    end
    update_line("aimora-emt-task-result-v1")
    update_line(plan.signature_sha256)
    update_line(checkpoint_signature)
    for occurrence in occurrences
        update_line(
            occurrence.task,
            UInt8(occurrence.family),
            occurrence.exact_instant.numerator,
            occurrence.exact_instant.denominator,
            UInt8(occurrence.stage),
            occurrence.priority,
            occurrence.activation_index,
            occurrence.sample_index,
            occurrence.release_index,
            occurrence.execution_index,
        )
    end
    return bytes2hex(SHA.digest!(context))
end

abstract type AbstractEMTTask end

emt_task_spec(project, declaration) = throw(MethodError(emt_task_spec, (project, declaration)))
emt_task_name(task::AbstractEMTTask) = throw(MethodError(emt_task_name, (task,)))
emt_task_checkpoint(task::AbstractEMTTask) = throw(MethodError(emt_task_checkpoint, (task,)))
restore_emt_task_checkpoint!(task::AbstractEMTTask, checkpoint) =
    throw(MethodError(restore_emt_task_checkpoint!, (task, checkpoint)))
read_emt_task(task::AbstractEMTTask, owner, instant::EMTLogicalTime, activation_index::Int) =
    throw(MethodError(read_emt_task, (task, owner, instant, activation_index)))
compute_emt_task(task::AbstractEMTTask, owner, input, instant::EMTLogicalTime, activation_index::Int) =
    throw(MethodError(compute_emt_task, (task, owner, input, instant, activation_index)))
write_emt_task!(task::AbstractEMTTask, owner, value, instant::EMTLogicalTime, activation_index::Int) =
    throw(MethodError(write_emt_task!, (task, owner, value, instant, activation_index)))

end

using Random
using Statistics

export EMTSwitchTimingRule,
       EMTSwitchEventSchedule,
       EMTDistributionFit,
       EMTEnsembleRunResult,
       EMTEnsembleChannelStatistics,
       EMTEnsembleResult,
       emt_switch_schedules,
       run_deck_emt_ensemble

struct EMTSwitchTimingRule
    switch_name::Symbol
    event::Symbol
    mode::Symbol
    center_time_s::Union{Nothing,Float64}
    spread_s::Float64
    dependency_switch::Union{Nothing,Symbol}
    dependency_event::Symbol
    dependency_delay_s::Float64

    function EMTSwitchTimingRule(
        switch_name,
        event::Symbol;
        mode::Symbol = :fixed,
        center_time_s::Union{Nothing,Real} = nothing,
        spread_s::Real = 0.0,
        dependency_switch = nothing,
        dependency_event::Symbol = :close,
        dependency_delay_s::Real = 0.0,
    )
        name = Symbol(switch_name)
        isempty(String(name)) &&
            throw(ArgumentError("ensemble switch name must not be empty"))
        event in (:close, :open) ||
            throw(ArgumentError("ensemble switch event must be close or open"))
        mode in (:fixed, :normal, :uniform, :systematic, :dependent) ||
            throw(ArgumentError("unsupported ensemble switch timing mode $mode"))
        center = center_time_s === nothing ?
            nothing : Float64(center_time_s)
        center === nothing ||
            (isfinite(center) && center >= 0.0) ||
            throw(ArgumentError("ensemble switch center time must be finite and nonnegative"))
        spread = Float64(spread_s)
        isfinite(spread) && spread >= 0.0 ||
            throw(ArgumentError("ensemble switch spread must be finite and nonnegative"))
        dependency = dependency_switch === nothing ?
            nothing : Symbol(dependency_switch)
        dependency_event in (:close, :open) ||
            throw(ArgumentError("ensemble dependency event must be close or open"))
        delay = Float64(dependency_delay_s)
        isfinite(delay) ||
            throw(ArgumentError("ensemble dependency delay must be finite"))
        if mode == :dependent
            dependency === nothing &&
                throw(ArgumentError("dependent switch timing requires dependency_switch"))
        elseif dependency !== nothing
            throw(ArgumentError("dependency_switch is valid only for dependent timing"))
        end
        return new(
            name,
            event,
            mode,
            center,
            spread,
            dependency,
            dependency_event,
            delay,
        )
    end
end

struct EMTSwitchEventSchedule
    energization_index::Int
    switch_names::Vector{Symbol}
    events::Vector{Symbol}
    event_times_s::Vector{Float64}
    dependency_invariants_passed::Bool
end

struct EMTDistributionFit
    distribution::Symbol
    location::Float64
    scale::Float64
    kolmogorov_smirnov_statistic::Float64
    fit_valid::Bool
end

struct EMTEnsembleRunResult
    schedule::EMTSwitchEventSchedule
    node_maximum_values::Vector{Float64}
    node_maximum_times_s::Vector{Float64}
    node_minimum_values::Vector{Float64}
    node_minimum_times_s::Vector{Float64}
    output_maximum_values::Vector{Float64}
    output_maximum_times_s::Vector{Float64}
    output_minimum_values::Vector{Float64}
    output_minimum_times_s::Vector{Float64}
    terminal_state::EMTTerminalState
    extrema_recomputed::Bool
end

struct EMTEnsembleChannelStatistics
    name::Symbol
    quantity_kind::Symbol
    maximum_samples::Vector{Float64}
    minimum_samples::Vector{Float64}
    maximum_mean::Float64
    maximum_standard_deviation::Float64
    minimum_mean::Float64
    minimum_standard_deviation::Float64
    absolute_extreme::Float64
    sorted_maximum_samples::Vector{Float64}
    maximum_exceedance_probabilities::Vector{Float64}
    normal_fit::EMTDistributionFit
    gumbel_fit::EMTDistributionFit
    preferred_distribution::Symbol
end

struct EMTEnsembleResult
    source::String
    seed::UInt64
    energization_count::Int
    timing_rules::Vector{EMTSwitchTimingRule}
    schedules::Vector{EMTSwitchEventSchedule}
    runs::Vector{EMTEnsembleRunResult}
    channel_statistics::Vector{EMTEnsembleChannelStatistics}
    deterministic_replay_signature::UInt64
    dependency_invariants_passed::Bool
    extrema_recomputed::Bool
    complete_emt_path::Bool
    physical_checks_passed::Bool
end

function _ensemble_switch_rows_by_name(parsed::DeckParser.DeckParseResult)
    rows = Dict{Symbol,Any}()
    for row in DeckParser.deck_over5_switch_rows(parsed)
        haskey(rows, row.name) &&
            throw(ArgumentError("ensemble switch names must be unique"))
        rows[row.name] = row
    end
    return rows
end

function _ensemble_base_event_time(row, event::Symbol)
    return event == :close ?
        Float64(row.close_time_s) :
        Float64(row.open_time_s)
end

function _ensemble_default_timing_rules(parsed::DeckParser.DeckParseResult)
    rules = EMTSwitchTimingRule[]
    for row in DeckParser.deck_over5_switch_rows(parsed)
        standard_deviation =
            Float64(row.random_opening_standard_deviation_s)
        standard_deviation > 0.0 || continue
        isfinite(row.open_time_s) ||
            throw(ArgumentError(
                "random-opening switch $(row.name) requires a finite base open time",
            ))
        push!(
            rules,
            EMTSwitchTimingRule(
                row.name,
                :open;
                mode = :normal,
                center_time_s = row.open_time_s,
                spread_s = standard_deviation,
            ),
        )
    end
    return rules
end

function _ensemble_checked_rule_order(
    parsed::DeckParser.DeckParseResult,
    timing_rules::AbstractVector{EMTSwitchTimingRule},
)
    rows = _ensemble_switch_rows_by_name(parsed)
    rules = collect(timing_rules)
    keys = Tuple{Symbol,Symbol}[]
    for rule in rules
        haskey(rows, rule.switch_name) ||
            throw(ArgumentError("ensemble rule references unknown switch $(rule.switch_name)"))
        key = (rule.switch_name, rule.event)
        key in keys &&
            throw(ArgumentError("ensemble switch event $key has duplicate timing rules"))
        push!(keys, key)
    end
    ordered = EMTSwitchTimingRule[]
    pending = copy(rules)
    while !isempty(pending)
        progress = false
        for index in reverse(eachindex(pending))
            rule = pending[index]
            if rule.mode != :dependent ||
               (rule.dependency_switch, rule.dependency_event) in [
                   (prior.switch_name, prior.event) for prior in ordered
               ]
                push!(ordered, rule)
                deleteat!(pending, index)
                progress = true
            end
        end
        progress ||
            throw(ArgumentError("ensemble switch timing dependencies contain a cycle or missing rule"))
    end
    return ordered, rows
end

function _ensemble_rule_time(
    rule::EMTSwitchTimingRule,
    run_index::Int,
    run_count::Int,
    rng::AbstractRNG,
    base_time_s::Float64,
    resolved_times::Dict{Tuple{Symbol,Symbol},Float64},
)
    center = rule.center_time_s === nothing ? base_time_s : rule.center_time_s
    if rule.mode != :dependent
        isfinite(center) ||
            throw(ArgumentError(
                "ensemble $(rule.event) event for $(rule.switch_name) requires a finite center time",
            ))
    end
    time = if rule.mode == :fixed
        center
    elseif rule.mode == :normal
        center + rule.spread_s * randn(rng)
    elseif rule.mode == :uniform
        center + rule.spread_s * (2.0 * rand(rng) - 1.0)
    elseif rule.mode == :systematic
        fraction = run_count == 1 ?
            0.0 : -1.0 + 2.0 * (run_index - 1) / (run_count - 1)
        center + rule.spread_s * fraction
    else
        dependency_key =
            (something(rule.dependency_switch), rule.dependency_event)
        haskey(resolved_times, dependency_key) ||
            throw(ArgumentError("ensemble dependent event is missing its upstream event"))
        resolved_times[dependency_key] + rule.dependency_delay_s
    end
    isfinite(time) ||
        throw(ArgumentError("ensemble switch event time must be finite"))
    return max(time, 0.0)
end

function emt_switch_schedules(
    parsed::DeckParser.DeckParseResult;
    energization_count::Union{Nothing,Integer} = nothing,
    seed::Integer = 0,
    timing_rules::AbstractVector{EMTSwitchTimingRule} =
        EMTSwitchTimingRule[],
)
    DeckParser.assert_deck_valid!(parsed)
    seed >= 0 ||
        throw(ArgumentError("EMT ensemble seed must be nonnegative"))
    deck_count = DeckParser.deck_output_schedule_options(parsed).energization_count
    run_count = energization_count === nothing ? deck_count : Int(energization_count)
    run_count > 0 ||
        throw(ArgumentError("EMT ensemble energization_count must be positive"))
    rules = isempty(timing_rules) ?
        _ensemble_default_timing_rules(parsed) :
        collect(timing_rules)
    ordered_rules, rows = _ensemble_checked_rule_order(parsed, rules)
    rng = MersenneTwister(seed)
    schedules = EMTSwitchEventSchedule[]
    for run_index in 1:run_count
        resolved_times = Dict{Tuple{Symbol,Symbol},Float64}()
        names = Symbol[]
        events = Symbol[]
        times = Float64[]
        dependency_checks = true
        for rule in ordered_rules
            row = rows[rule.switch_name]
            base_time = _ensemble_base_event_time(row, rule.event)
            time = _ensemble_rule_time(
                rule,
                run_index,
                run_count,
                rng,
                base_time,
                resolved_times,
            )
            key = (rule.switch_name, rule.event)
            resolved_times[key] = time
            push!(names, rule.switch_name)
            push!(events, rule.event)
            push!(times, time)
            if rule.mode == :dependent
                upstream =
                    resolved_times[
                        (something(rule.dependency_switch), rule.dependency_event)
                    ]
                dependency_checks &=
                    isapprox(
                        time,
                        max(upstream + rule.dependency_delay_s, 0.0);
                        atol = 16.0 * eps(Float64),
                        rtol = 0.0,
                    )
            end
        end
        push!(
            schedules,
            EMTSwitchEventSchedule(
                run_index,
                names,
                events,
                times,
                dependency_checks,
            ),
        )
    end
    return schedules
end

function _apply_ensemble_schedule!(
    parsed::DeckParser.DeckParseResult,
    schedule::EMTSwitchEventSchedule,
)
    length(schedule.switch_names) ==
        length(schedule.events) ==
        length(schedule.event_times_s) ||
        throw(ArgumentError("ensemble schedule vectors must have equal length"))
    for index in eachindex(schedule.switch_names)
        name = schedule.switch_names[index]
        element_index = findfirst(==(name), parsed.element_names)
        element_index === nothing &&
            throw(ArgumentError("ensemble schedule references unknown element $name"))
        element = parsed.elements[element_index]
        element isa TimeSwitch ||
            throw(ArgumentError("ensemble schedule element $name is not a time switch"))
        event = schedule.events[index]
        time = schedule.event_times_s[index]
        event == :close ?
            (element.close_time_s = time) :
            event == :open ?
            (element.open_time_s = time) :
            throw(ArgumentError("ensemble schedule event must be close or open"))
    end
    for element in parsed.elements
        element isa TimeSwitch || continue
        if !element.initially_closed &&
           isfinite(element.close_time_s) &&
           isfinite(element.open_time_s) &&
           element.open_time_s < element.close_time_s
            throw(ArgumentError("ensemble switch opens before it closes"))
        end
    end
    return parsed
end

function _ensemble_extrema_recomputed(trace::DeckEMTTrace)
    node_maximum = vec(maximum(trace.voltage_pu; dims = 2))
    node_minimum = vec(minimum(trace.voltage_pu; dims = 2))
    output_maximum = isempty(trace.output_pu) ?
        Float64[] : vec(maximum(trace.output_pu; dims = 2))
    output_minimum = isempty(trace.output_pu) ?
        Float64[] : vec(minimum(trace.output_pu; dims = 2))
    return (
        node_maximum == trace.node_maximum_values &&
        node_minimum == trace.node_minimum_values &&
        output_maximum == trace.output_maximum_values &&
        output_minimum == trace.output_minimum_values
    )
end

function _normal_cdf(value::Float64)
    x = abs(value)
    t = inv(1.0 + 0.2316419 * x)
    density = 0.3989422804014327 * exp(-0.5 * x * x)
    polynomial =
        t * (
            0.319381530 +
            t * (
                -0.356563782 +
                t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))
            )
        )
    upper = 1.0 - density * polynomial
    return value >= 0.0 ? upper : 1.0 - upper
end

function _ensemble_ks_statistic(
    samples::Vector{Float64},
    cdf::F,
) where {F}
    sorted = sort(samples)
    count = length(sorted)
    statistic = 0.0
    for (index, value) in enumerate(sorted)
        probability = clamp(cdf(value), 0.0, 1.0)
        statistic = max(
            statistic,
            abs(probability - (index - 1) / count),
            abs(index / count - probability),
        )
    end
    return statistic
end

function _ensemble_distribution_fits(samples::Vector{Float64})
    count = length(samples)
    count > 0 || throw(ArgumentError("ensemble distribution samples must not be empty"))
    location = mean(samples)
    sample_scale = count > 1 ? std(samples; corrected = true) : 0.0
    scale_floor = max(maximum(abs, samples; init = 0.0), 1.0) * eps(Float64)
    normal_scale = max(sample_scale, scale_floor)
    normal_ks = _ensemble_ks_statistic(
        samples,
        value -> _normal_cdf((value - location) / normal_scale),
    )
    normal = EMTDistributionFit(
        :normal,
        location,
        normal_scale,
        normal_ks,
        count >= 2 && sample_scale > 0.0,
    )
    gumbel_scale = max(sample_scale * sqrt(6.0) / pi, scale_floor)
    gumbel_location = location - 0.5772156649015329 * gumbel_scale
    gumbel_ks = _ensemble_ks_statistic(
        samples,
        value -> exp(-exp(-(value - gumbel_location) / gumbel_scale)),
    )
    gumbel = EMTDistributionFit(
        :gumbel,
        gumbel_location,
        gumbel_scale,
        gumbel_ks,
        count >= 2 && sample_scale > 0.0,
    )
    return normal, gumbel
end

function _ensemble_channel_statistics(
    name::Symbol,
    quantity_kind::Symbol,
    maximum_samples::Vector{Float64},
    minimum_samples::Vector{Float64},
)
    length(maximum_samples) == length(minimum_samples) > 0 ||
        throw(ArgumentError("ensemble channel extrema samples must be aligned and nonempty"))
    normal, gumbel = _ensemble_distribution_fits(maximum_samples)
    sorted_maximum = sort(maximum_samples)
    count = length(sorted_maximum)
    exceedance = Float64[
        1.0 - (index - 0.5) / count for index in 1:count
    ]
    preferred = normal.kolmogorov_smirnov_statistic <=
        gumbel.kolmogorov_smirnov_statistic ? :normal : :gumbel
    return EMTEnsembleChannelStatistics(
        name,
        quantity_kind,
        maximum_samples,
        minimum_samples,
        mean(maximum_samples),
        count > 1 ? std(maximum_samples; corrected = true) : 0.0,
        mean(minimum_samples),
        count > 1 ? std(minimum_samples; corrected = true) : 0.0,
        max(maximum(abs, maximum_samples), maximum(abs, minimum_samples)),
        sorted_maximum,
        exceedance,
        normal,
        gumbel,
        preferred,
    )
end

function _ensemble_replay_signature(
    schedules::Vector{EMTSwitchEventSchedule},
)
    hash_value = UInt64(0xcbf29ce484222325)
    for schedule in schedules
        for index in eachindex(schedule.switch_names)
            for byte in codeunits(String(schedule.switch_names[index]))
                hash_value = (hash_value ⊻ UInt64(byte)) * UInt64(0x100000001b3)
            end
            for byte in codeunits(String(schedule.events[index]))
                hash_value = (hash_value ⊻ UInt64(byte)) * UInt64(0x100000001b3)
            end
            hash_value =
                (hash_value ⊻ reinterpret(UInt64, schedule.event_times_s[index])) *
                UInt64(0x100000001b3)
        end
    end
    return hash_value
end

function run_deck_emt_ensemble(
    parsed::DeckParser.DeckParseResult;
    energization_count::Union{Nothing,Integer} = nothing,
    seed::Integer = 0,
    timing_rules::AbstractVector{EMTSwitchTimingRule} =
        EMTSwitchTimingRule[],
    schedules::Union{Nothing,AbstractVector{EMTSwitchEventSchedule}} = nothing,
    execution_kwargs::NamedTuple = (time_horizon = :deck,),
)
    DeckParser.assert_deck_valid!(parsed)
    seed >= 0 ||
        throw(ArgumentError("EMT ensemble seed must be nonnegative"))
    resolved_rules = isempty(timing_rules) ?
        _ensemble_default_timing_rules(parsed) :
        collect(timing_rules)
    resolved_schedules = schedules === nothing ?
        emt_switch_schedules(
            parsed;
            energization_count,
            seed,
            timing_rules = resolved_rules,
        ) :
        collect(schedules)
    !isempty(resolved_schedules) ||
        throw(ArgumentError("EMT ensemble schedules must not be empty"))
    for (index, schedule) in enumerate(resolved_schedules)
        schedule.energization_index == index ||
            throw(ArgumentError("EMT ensemble schedule indices must be consecutive"))
    end
    runs = EMTEnsembleRunResult[]
    node_names = Symbol[]
    output_names = Symbol[]
    for schedule in resolved_schedules
        run_deck = deepcopy(parsed)
        _apply_ensemble_schedule!(run_deck, schedule)
        execution = run_deck_emt_execution(run_deck; execution_kwargs...)
        trace = execution.trace
        if isempty(node_names)
            node_names = copy(trace.node_names)
            output_names = copy(trace.output_channel_names)
        else
            node_names == trace.node_names &&
                output_names == trace.output_channel_names ||
                throw(ArgumentError("EMT ensemble runs changed their output schema"))
        end
        push!(
            runs,
            EMTEnsembleRunResult(
                schedule,
                copy(trace.node_maximum_values),
                copy(trace.node_maximum_times_s),
                copy(trace.node_minimum_values),
                copy(trace.node_minimum_times_s),
                copy(trace.output_maximum_values),
                copy(trace.output_maximum_times_s),
                copy(trace.output_minimum_values),
                copy(trace.output_minimum_times_s),
                execution.terminal_state,
                _ensemble_extrema_recomputed(trace),
            ),
        )
    end
    statistics = EMTEnsembleChannelStatistics[]
    for (channel, name) in enumerate(node_names)
        push!(
            statistics,
            _ensemble_channel_statistics(
                name,
                :node_voltage,
                Float64[run.node_maximum_values[channel] for run in runs],
                Float64[run.node_minimum_values[channel] for run in runs],
            ),
        )
    end
    for (channel, name) in enumerate(output_names)
        push!(
            statistics,
            _ensemble_channel_statistics(
                name,
                :output,
                Float64[run.output_maximum_values[channel] for run in runs],
                Float64[run.output_minimum_values[channel] for run in runs],
            ),
        )
    end
    dependency_checks = all(
        schedule -> schedule.dependency_invariants_passed,
        resolved_schedules,
    )
    extrema_checks = all(run -> run.extrema_recomputed, runs)
    complete_path = all(
        run -> run.terminal_state.physical_checks_passed,
        runs,
    )
    return EMTEnsembleResult(
        parsed.source,
        UInt64(seed),
        length(resolved_schedules),
        resolved_rules,
        resolved_schedules,
        runs,
        statistics,
        _ensemble_replay_signature(resolved_schedules),
        dependency_checks,
        extrema_checks,
        complete_path,
        dependency_checks && extrema_checks && complete_path,
    )
end

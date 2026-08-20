export MeasurementAnalogCandidate,
       DelayedMeasurementSample,
       MeasurementSample,
       MeasurementChainRuntime,
       MeasurementChainSnapshot,
       initialize_measurement_chain_at_tick!,
       prepare_measurement_analog_step!,
       discard_measurement_analog_step!,
       accept_measurement_analog_step!,
       measurement_chain_snapshot,
       restore_measurement_chain_snapshot!,
       measurement_chain_result_signature,
       measurement_due_at_tick

mutable struct MeasurementAnalogCandidate
    state::Matrix{Float64}
    input::Vector{Float64}
    output::Vector{Float64}
    timestep_s::Float64
end

struct DelayedMeasurementSample
    source_tick::Int
    release_tick::Int
    values::Vector{Float64}
    codes::Union{Nothing,Vector{Int64}}
    clipped::BitVector
end

struct MeasurementSample
    source_tick::Int
    release_tick::Int
    source_time_s::Float64
    release_time_s::Float64
    instantaneous::Vector{Float64}
    codes::Union{Nothing,Vector{Int64}}
    clipped::BitVector
    sliding_rms::Union{Nothing,Vector{Float64}}
    fundamental_rms_phasors::Union{Nothing,Vector{ComplexF64}}
    sequence_phasors::Union{Nothing,NamedTuple{(:zero, :positive, :negative),NTuple{3,ComplexF64}}}
    frequency_hz::Union{Nothing,Float64}
    quality::Symbol
    deterministic_signature_sha256::String
end

mutable struct MeasurementChainRuntime{S<:MeasurementChainSpecification}
    specification::S
    specification_signature_sha256::String
    analog_state::Matrix{Float64}
    previous_input::Vector{Float64}
    analog_output::Vector{Float64}
    candidate::MeasurementAnalogCandidate
    candidate_active::Bool
    delayed_samples::Vector{DelayedMeasurementSample}
    delayed_sample_head::Int
    held_values::Vector{Float64}
    held_codes::Union{Nothing,Vector{Int64}}
    window_values::Matrix{Float64}
    window_ticks::Vector{Int}
    window_next_index::Int
    retained_window_count::Int
    squared_sums::Vector{Float64}
    positive_sequence_history::Vector{ComplexF64}
    positive_sequence_time_history_s::Vector{Float64}
    samples::Vector{MeasurementSample}
    last_accepted_tick::Int
    accepted_analog_step_count::Int
    accepted_sample_count::Int
    released_sample_count::Int
end

function MeasurementChainRuntime(
    specification::MeasurementChainSpecification;
    initial_input::AbstractVector{<:Real}=zeros(length(specification.channel_names)),
    initial_analog_state::AbstractMatrix{<:Real}=zeros(
        size(specification.conditioning.state_matrix_per_s, 1),
        length(specification.channel_names),
    ),
)
    channel_count = length(specification.channel_names)
    state_count = size(specification.conditioning.state_matrix_per_s, 1)
    input = Float64.(initial_input)
    state = Matrix{Float64}(initial_analog_state)
    length(input) == channel_count || throw(DimensionMismatch(
        "initial measurement input count must match the declared channels",
    ))
    size(state) == (state_count, channel_count) || throw(DimensionMismatch(
        "initial analog measurement state must match state and channel counts",
    ))
    all(isfinite, input) && all(isfinite, state) || throw(ArgumentError(
        "initial measurement input and analog state must be finite",
    ))
    all(value -> specification.minimum_input <= value <= specification.maximum_input, input) ||
        throw(ArgumentError("initial measurement input is outside the declared domain"))
    model = specification.conditioning
    output = Float64[
        dot(model.output_vector, view(state, :, channel)) + model.direct_gain * input[channel]
        for channel in 1:channel_count
    ]
    window_count = length(specification.acquisition.window_weights_newest_first)
    return MeasurementChainRuntime(
        specification,
        measurement_chain_signature(specification),
        state,
        input,
        output,
        MeasurementAnalogCandidate(
            similar(state),
            similar(input),
            similar(output),
            NaN,
        ),
        false,
        DelayedMeasurementSample[],
        1,
        copy(output),
        nothing,
        zeros(window_count, channel_count),
        fill(-1, window_count),
        1,
        0,
        zeros(channel_count),
        ComplexF64[],
        Float64[],
        MeasurementSample[],
        -1,
        0,
        0,
        0,
    )
end

measurement_due_at_tick(settings::MeasurementAcquisitionSettings, tick::Integer) =
    tick >= settings.first_sample_tick &&
    (tick - settings.first_sample_tick) % settings.sample_period_ticks == 0

function _first_measurement_due_tick_after(
    settings::MeasurementAcquisitionSettings,
    tick::Int,
)
    tick < settings.first_sample_tick && return settings.first_sample_tick
    elapsed_periods = fld(
        tick - settings.first_sample_tick,
        settings.sample_period_ticks,
    )
    return settings.first_sample_tick +
        (elapsed_periods + 1) * settings.sample_period_ticks
end

function initialize_measurement_chain_at_tick!(
    runtime::MeasurementChainRuntime,
    accepted_tick::Integer,
)
    !runtime.candidate_active || throw(ArgumentError(
        "measurement initialization is unavailable during an active analog trial",
    ))
    runtime.last_accepted_tick == -1 && runtime.accepted_analog_step_count == 0 &&
        runtime.accepted_sample_count == 0 && runtime.released_sample_count == 0 ||
        throw(ArgumentError("measurement runtime has already been initialized or advanced"))
    tick = Int(accepted_tick)
    tick >= 0 || throw(ArgumentError("measurement initialization tick must be nonnegative"))
    settings = runtime.specification.acquisition
    if settings.first_sample_tick < tick
        throw(MeasurementChainRefusal(
            :missed_measurement_acquisition,
            :initialize,
            runtime.specification.id,
            runtime.specification.family,
            "measurement initialization begins after its first required acquisition",
            (expected_tick=settings.first_sample_tick, observed_tick=tick),
        ))
    end
    runtime.last_accepted_tick = tick
    if measurement_due_at_tick(settings, tick)
        values, codes, clipped = _quantize_measurement(
            runtime.specification.quantizer,
            runtime.analog_output,
        )
        push!(
            runtime.delayed_samples,
            DelayedMeasurementSample(
                tick,
                tick + settings.delay_ticks,
                values,
                codes,
                clipped,
            ),
        )
        runtime.accepted_sample_count += 1
    end
    return _release_due_measurement_samples!(MeasurementSample[], runtime, tick)
end

function _prepare_measurement_analog_step!(
    runtime::MeasurementChainRuntime,
    timestep_s::Real,
)
    specification = runtime.specification
    input = runtime.candidate.input
    all(isfinite, input) || throw(ArgumentError("measurement input must be finite"))
    all(value -> specification.minimum_input <= value <= specification.maximum_input, input) ||
        throw(MeasurementChainRefusal(
            :input_outside_validity_domain,
            :prepare_analog_step,
            specification.id,
            specification.family,
            "accepted-time input is outside the declared measurement domain",
            (minimum=minimum(input), maximum=maximum(input)),
        ))
    timestep = Float64(timestep_s)
    isfinite(timestep) && 0.0 < timestep <= specification.maximum_timestep_s ||
        throw(MeasurementChainRefusal(
            :timestep_outside_validity_domain,
            :prepare_analog_step,
            specification.id,
            specification.family,
            "measurement timestep is outside the declared fixed-step domain",
            (timestep_s=timestep, maximum_timestep_s=specification.maximum_timestep_s),
        ))
    model = specification.conditioning
    state_count = size(model.state_matrix_per_s, 1)
    candidate_state = if state_count == 0
        runtime.candidate.state
    elseif state_count == 1
        half_timestep = 0.5 * timestep
        state_coefficient = model.state_matrix_per_s[1, 1]
        input_coefficient = model.input_vector_per_s[1]
        left_coefficient = 1.0 - half_timestep * state_coefficient
        right_coefficient = 1.0 + half_timestep * state_coefficient
        next_state = runtime.candidate.state
        for channel in eachindex(input)
            forcing = half_timestep * input_coefficient *
                (runtime.previous_input[channel] + input[channel])
            next_state[1, channel] = (
                right_coefficient * runtime.analog_state[1, channel] + forcing
            ) / left_coefficient
        end
        next_state
    else
        identity_matrix = Matrix{Float64}(I, state_count, state_count)
        left = identity_matrix .- 0.5 * timestep .* model.state_matrix_per_s
        right = (identity_matrix .+ 0.5 * timestep .* model.state_matrix_per_s) *
            runtime.analog_state
        forcing = 0.5 * timestep .* model.input_vector_per_s *
            transpose(runtime.previous_input .+ input)
        runtime.candidate.state .= left \ (right .+ forcing)
    end
    candidate_output = runtime.candidate.output
    for channel in eachindex(input)
        state_output = 0.0
        for state_index in eachindex(model.output_vector)
            state_output += model.output_vector[state_index] *
                candidate_state[state_index, channel]
        end
        candidate_output[channel] = state_output + model.direct_gain * input[channel]
    end
    all(isfinite, candidate_state) && all(isfinite, candidate_output) ||
        throw(MeasurementChainRefusal(
            :nonfinite_analog_candidate,
            :prepare_analog_step,
            specification.id,
            specification.family,
            "analog measurement candidate contains a nonfinite value",
            (timestep_s=timestep,),
        ))
    runtime.candidate.timestep_s = timestep
    runtime.candidate_active = true
    return runtime.candidate
end

function prepare_measurement_analog_step!(
    runtime::MeasurementChainRuntime,
    next_input::AbstractVector{<:Real},
    timestep_s::Real,
)
    !runtime.candidate_active || throw(ArgumentError(
        "measurement analog trial is already active",
    ))
    length(next_input) == length(runtime.specification.channel_names) ||
        throw(DimensionMismatch(
            "measurement input count must match the declared channels",
        ))
    for channel in eachindex(runtime.candidate.input)
        runtime.candidate.input[channel] = Float64(next_input[channel])
    end
    return _prepare_measurement_analog_step!(runtime, timestep_s)
end

function prepare_measurement_analog_step!(
    runtime::MeasurementChainRuntime,
    next_input::Real,
    timestep_s::Real,
)
    !runtime.candidate_active || throw(ArgumentError(
        "measurement analog trial is already active",
    ))
    length(runtime.specification.channel_names) == 1 || throw(DimensionMismatch(
        "scalar measurement input requires exactly one declared channel",
    ))
    runtime.candidate.input[1] = Float64(next_input)
    return _prepare_measurement_analog_step!(runtime, timestep_s)
end

function discard_measurement_analog_step!(runtime::MeasurementChainRuntime)
    runtime.candidate_active = false
    return runtime
end

function _round_measurement_code(value::Float64, rule::MeasurementQuantizerTieRule)
    if rule === MeasurementTiesToEven
        return round(Int64, value, RoundNearest)
    elseif rule === MeasurementTiesAwayFromZero
        return round(Int64, value, RoundNearestTiesAway)
    end
    nearest = round(value)
    if abs(value - nearest) == 0.5
        return trunc(Int64, value)
    end
    return round(Int64, value, RoundNearest)
end

function _quantize_measurement(
    quantizer::UnquantizedMeasurement,
    values::Vector{Float64},
)
    quantized = similar(values)
    clipped = falses(length(values))
    for channel in eachindex(values)
        value = values[channel]
        quantized[channel] = clamp(
            value,
            quantizer.lower_limit,
            quantizer.upper_limit,
        )
        clipped[channel] = value != quantized[channel]
    end
    return quantized, nothing, clipped
end

function _quantize_measurement(
    quantizer::UniformMeasurementQuantizer,
    values::Vector{Float64},
)
    quantized = similar(values)
    codes = Vector{Int64}(undef, length(values))
    clipped = falses(length(values))
    for channel in eachindex(values)
        value = values[channel]
        clipped_value = clamp(value, quantizer.lower_limit, quantizer.upper_limit)
        clipped[channel] = value != clipped_value
        code = clamp(
            _round_measurement_code(
                (clipped_value - quantizer.engineering_offset) /
                    quantizer.engineering_step,
                quantizer.tie_rule,
            ),
            quantizer.minimum_code,
            quantizer.maximum_code,
        )
        codes[channel] = code
        quantized[channel] = quantizer.engineering_offset +
            quantizer.engineering_step * code
    end
    return quantized, codes, clipped
end

function _sample_signature(
    runtime::MeasurementChainRuntime,
    source_tick,
    release_tick,
    values,
    codes,
    rms,
    phasors,
    sequence,
    frequency,
    quality,
)
    io = IOBuffer()
    for value in (
        runtime.specification_signature_sha256,
        source_tick,
        release_tick,
        values,
        codes,
        rms,
        phasors,
        sequence,
        frequency,
        quality,
    )
        _write_measurement_signature_value(io, value)
        print(io, '\n')
    end
    return bytes2hex(sha256(take!(io)))
end

function _measurement_window_estimates!(
    runtime::MeasurementChainRuntime,
    delayed::DelayedMeasurementSample,
)
    window_count = size(runtime.window_values, 1)
    insertion_index = runtime.window_next_index
    if runtime.retained_window_count == window_count
        for channel in axes(runtime.window_values, 2)
            runtime.squared_sums[channel] -=
                abs2(runtime.window_values[insertion_index, channel])
        end
    else
        runtime.retained_window_count += 1
    end
    runtime.window_values[insertion_index, :] .= delayed.values
    runtime.window_ticks[insertion_index] = delayed.source_tick
    for channel in eachindex(delayed.values)
        runtime.squared_sums[channel] += abs2(delayed.values[channel])
    end
    runtime.window_next_index = mod1(insertion_index + 1, window_count)
    runtime.retained_window_count == window_count || return nothing, nothing, nothing, nothing

    rms = sqrt.(max.(runtime.squared_sums, 0.0) ./ window_count)
    settings = runtime.specification.acquisition
    coherent_gain = sum(settings.window_weights_newest_first)
    phasors = zeros(ComplexF64, size(runtime.window_values, 2))
    for age in 0:(window_count - 1)
        index = mod1(insertion_index - age, window_count)
        tick = runtime.window_ticks[index]
        time_s = tick * settings.tick_s
        kernel = settings.window_weights_newest_first[age + 1] *
            cis(-2.0 * pi * settings.nominal_frequency_hz * time_s)
        for channel in eachindex(phasors)
            phasors[channel] += runtime.window_values[index, channel] * kernel
        end
    end
    phasors .*= sqrt(2.0) / coherent_gain

    sequence = nothing
    frequency = nothing
    if runtime.specification.phase_order == [:a, :b, :c]
        rotation = cis(2.0 * pi / 3.0)
        phase_a, phase_b, phase_c = phasors
        zero = (phase_a + phase_b + phase_c) / 3.0
        positive = (phase_a + rotation * phase_b + rotation^2 * phase_c) / 3.0
        negative = (phase_a + rotation^2 * phase_b + rotation * phase_c) / 3.0
        sequence = (zero=zero, positive=positive, negative=negative)
        if abs(positive) >= settings.positive_sequence_threshold
            push!(runtime.positive_sequence_history, positive)
            push!(runtime.positive_sequence_time_history_s, delayed.source_tick * settings.tick_s)
            required_count = settings.frequency_update_separation + 1
            while length(runtime.positive_sequence_history) > required_count
                popfirst!(runtime.positive_sequence_history)
                popfirst!(runtime.positive_sequence_time_history_s)
            end
            if length(runtime.positive_sequence_history) == required_count
                phase_increment = angle(
                    runtime.positive_sequence_history[end] *
                    conj(runtime.positive_sequence_history[1]),
                )
                elapsed = runtime.positive_sequence_time_history_s[end] -
                    runtime.positive_sequence_time_history_s[1]
                elapsed > 0.0 || throw(ArgumentError(
                    "measurement frequency estimator timestamps must increase",
                ))
                frequency = settings.nominal_frequency_hz +
                    phase_increment / (2.0 * pi * elapsed)
            end
        else
            empty!(runtime.positive_sequence_history)
            empty!(runtime.positive_sequence_time_history_s)
        end
    end
    return rms, phasors, sequence, frequency
end

function _release_measurement_sample!(
    runtime::MeasurementChainRuntime,
    delayed::DelayedMeasurementSample,
)
    rms, phasors, sequence, frequency = _measurement_window_estimates!(runtime, delayed)
    runtime.held_values .= delayed.values
    runtime.held_codes = delayed.codes === nothing ? nothing : copy(delayed.codes)
    quality = rms === nothing ? :window_incomplete :
        (sequence !== nothing && frequency === nothing ? :frequency_unavailable : :valid)
    signature = _sample_signature(
        runtime,
        delayed.source_tick,
        delayed.release_tick,
        delayed.values,
        delayed.codes,
        rms,
        phasors,
        sequence,
        frequency,
        quality,
    )
    settings = runtime.specification.acquisition
    sample = MeasurementSample(
        delayed.source_tick,
        delayed.release_tick,
        delayed.source_tick * settings.tick_s,
        delayed.release_tick * settings.tick_s,
        delayed.values,
        delayed.codes,
        delayed.clipped,
        rms,
        phasors,
        sequence,
        frequency,
        quality,
        signature,
    )
    if length(runtime.samples) ==
            runtime.specification.acquisition.maximum_retained_samples
        popfirst!(runtime.samples)
    end
    push!(runtime.samples, sample)
    runtime.released_sample_count += 1
    return sample
end

function _release_due_measurement_samples!(
    released::Vector{MeasurementSample},
    runtime::MeasurementChainRuntime,
    tick::Int,
)
    empty!(released)
    while runtime.delayed_sample_head <= length(runtime.delayed_samples)
        delayed = runtime.delayed_samples[runtime.delayed_sample_head]
        delayed.release_tick > tick && break
        delayed.release_tick == tick || throw(MeasurementChainRefusal(
            :missed_measurement_release,
            :accept_analog_step,
            runtime.specification.id,
            runtime.specification.family,
            "a delayed measurement sample release tick was skipped",
            (expected_tick=delayed.release_tick, observed_tick=tick),
        ))
        runtime.delayed_sample_head += 1
        push!(released, _release_measurement_sample!(runtime, delayed))
    end
    consumed_count = runtime.delayed_sample_head - 1
    if consumed_count == length(runtime.delayed_samples)
        empty!(runtime.delayed_samples)
        runtime.delayed_sample_head = 1
    elseif consumed_count >= 1024 && consumed_count >= length(runtime.delayed_samples) ÷ 2
        deleteat!(runtime.delayed_samples, 1:consumed_count)
        runtime.delayed_sample_head = 1
    end
    return released
end

function _validated_measurement_acceptance_tick(
    runtime::MeasurementChainRuntime,
    accepted_tick::Integer,
)
    runtime.candidate_active || throw(ArgumentError(
        "measurement analog step must be prepared before acceptance",
    ))
    tick = Int(accepted_tick)
    tick > runtime.last_accepted_tick || throw(ArgumentError(
        "accepted measurement ticks must increase strictly",
    ))
    missed_sample_tick = _first_measurement_due_tick_after(
        runtime.specification.acquisition,
        runtime.last_accepted_tick,
    )
    missed_sample_tick < tick && throw(MeasurementChainRefusal(
        :missed_measurement_acquisition,
        :accept_analog_step,
        runtime.specification.id,
        runtime.specification.family,
        "an exact measurement acquisition tick was skipped",
        (expected_tick=missed_sample_tick, observed_tick=tick),
    ))
    if runtime.delayed_sample_head <= length(runtime.delayed_samples)
        expected_release_tick =
            runtime.delayed_samples[runtime.delayed_sample_head].release_tick
        expected_release_tick < tick && throw(MeasurementChainRefusal(
            :missed_measurement_release,
            :accept_analog_step,
            runtime.specification.id,
            runtime.specification.family,
            "a delayed measurement sample release tick was skipped",
            (expected_tick=expected_release_tick, observed_tick=tick),
        ))
    end
    return tick
end

function accept_measurement_analog_step!(
    released::Vector{MeasurementSample},
    runtime::MeasurementChainRuntime,
    accepted_tick::Integer,
)
    tick = _validated_measurement_acceptance_tick(runtime, accepted_tick)
    candidate = runtime.candidate
    runtime.analog_state .= candidate.state
    runtime.previous_input .= candidate.input
    runtime.analog_output .= candidate.output
    runtime.candidate_active = false
    runtime.last_accepted_tick = tick
    runtime.accepted_analog_step_count += 1
    settings = runtime.specification.acquisition
    if measurement_due_at_tick(settings, tick)
        values, codes, clipped = _quantize_measurement(
            runtime.specification.quantizer,
            runtime.analog_output,
        )
        push!(
            runtime.delayed_samples,
            DelayedMeasurementSample(
                tick,
                tick + settings.delay_ticks,
                values,
                codes,
                clipped,
            ),
        )
        runtime.accepted_sample_count += 1
    end
    return _release_due_measurement_samples!(released, runtime, tick)
end

function accept_measurement_analog_step!(
    runtime::MeasurementChainRuntime,
    accepted_tick::Integer,
)
    return accept_measurement_analog_step!(
        MeasurementSample[],
        runtime,
        accepted_tick,
    )
end

struct MeasurementChainSnapshot
    schema_version::Int
    specification_signature_sha256::String
    state::NamedTuple
    deterministic_signature_sha256::String
end

function _measurement_runtime_state(runtime::MeasurementChainRuntime)
    return (
        analog_state=copy(runtime.analog_state),
        previous_input=copy(runtime.previous_input),
        analog_output=copy(runtime.analog_output),
        delayed_samples=deepcopy(
            runtime.delayed_samples[runtime.delayed_sample_head:end],
        ),
        held_values=copy(runtime.held_values),
        held_codes=runtime.held_codes === nothing ? nothing : copy(runtime.held_codes),
        window_values=copy(runtime.window_values),
        window_ticks=copy(runtime.window_ticks),
        window_next_index=runtime.window_next_index,
        retained_window_count=runtime.retained_window_count,
        squared_sums=copy(runtime.squared_sums),
        positive_sequence_history=copy(runtime.positive_sequence_history),
        positive_sequence_time_history_s=copy(runtime.positive_sequence_time_history_s),
        samples=deepcopy(runtime.samples),
        last_accepted_tick=runtime.last_accepted_tick,
        accepted_analog_step_count=runtime.accepted_analog_step_count,
        accepted_sample_count=runtime.accepted_sample_count,
        released_sample_count=runtime.released_sample_count,
    )
end

function _measurement_snapshot_signature(specification_signature, state)
    io = IOBuffer()
    _write_measurement_signature_value(io, specification_signature)
    _write_measurement_signature_value(io, state)
    return bytes2hex(sha256(take!(io)))
end

function measurement_chain_snapshot(runtime::MeasurementChainRuntime)
    !runtime.candidate_active || throw(ArgumentError(
        "measurement snapshot cannot be taken during an active analog trial",
    ))
    state = _measurement_runtime_state(runtime)
    return MeasurementChainSnapshot(
        1,
        runtime.specification_signature_sha256,
        state,
        _measurement_snapshot_signature(runtime.specification_signature_sha256, state),
    )
end

const _MEASUREMENT_SNAPSHOT_STATE_FIELDS = (
    :analog_state,
    :previous_input,
    :analog_output,
    :delayed_samples,
    :held_values,
    :held_codes,
    :window_values,
    :window_ticks,
    :window_next_index,
    :retained_window_count,
    :squared_sums,
    :positive_sequence_history,
    :positive_sequence_time_history_s,
    :samples,
    :last_accepted_tick,
    :accepted_analog_step_count,
    :accepted_sample_count,
    :released_sample_count,
)

function _validated_measurement_snapshot_state(
    runtime::MeasurementChainRuntime,
    state::NamedTuple,
)
    all(field -> hasproperty(state, field), _MEASUREMENT_SNAPSHOT_STATE_FIELDS) ||
        throw(ArgumentError("measurement snapshot state inventory is incomplete"))
    channel_count = length(runtime.specification.channel_names)
    analog_state = Matrix{Float64}(state.analog_state)
    previous_input = Vector{Float64}(state.previous_input)
    analog_output = Vector{Float64}(state.analog_output)
    delayed_samples = Vector{DelayedMeasurementSample}(deepcopy(state.delayed_samples))
    held_values = Vector{Float64}(state.held_values)
    held_codes = state.held_codes === nothing ?
        nothing : Vector{Int64}(state.held_codes)
    window_values = Matrix{Float64}(state.window_values)
    window_ticks = Vector{Int}(state.window_ticks)
    window_next_index = Int(state.window_next_index)
    retained_window_count = Int(state.retained_window_count)
    squared_sums = Vector{Float64}(state.squared_sums)
    positive_sequence_history = Vector{ComplexF64}(state.positive_sequence_history)
    positive_sequence_time_history_s =
        Vector{Float64}(state.positive_sequence_time_history_s)
    samples = Vector{MeasurementSample}(deepcopy(state.samples))
    last_accepted_tick = Int(state.last_accepted_tick)
    accepted_analog_step_count = Int(state.accepted_analog_step_count)
    accepted_sample_count = Int(state.accepted_sample_count)
    released_sample_count = Int(state.released_sample_count)

    size(analog_state) == size(runtime.analog_state) || throw(DimensionMismatch(
        "measurement snapshot analog state size is incompatible",
    ))
    length(previous_input) == channel_count || throw(DimensionMismatch(
        "measurement snapshot previous-input size is incompatible",
    ))
    length(analog_output) == channel_count || throw(DimensionMismatch(
        "measurement snapshot analog-output size is incompatible",
    ))
    length(held_values) == channel_count || throw(DimensionMismatch(
        "measurement snapshot held-value size is incompatible",
    ))
    held_codes === nothing || length(held_codes) == channel_count ||
        throw(DimensionMismatch("measurement snapshot held-code size is incompatible"))
    size(window_values) == size(runtime.window_values) || throw(DimensionMismatch(
        "measurement snapshot estimator window size is incompatible",
    ))
    length(window_ticks) == size(window_values, 1) || throw(DimensionMismatch(
        "measurement snapshot estimator tick-window size is incompatible",
    ))
    length(squared_sums) == channel_count || throw(DimensionMismatch(
        "measurement snapshot squared-sum size is incompatible",
    ))
    length(positive_sequence_history) ==
        length(positive_sequence_time_history_s) || throw(DimensionMismatch(
        "measurement snapshot frequency history sizes are incompatible",
    ))
    length(positive_sequence_history) <=
        runtime.specification.acquisition.frequency_update_separation + 1 ||
        throw(ArgumentError("measurement snapshot frequency history is over-retained"))
    1 <= window_next_index <= size(window_values, 1) || throw(ArgumentError(
        "measurement snapshot estimator window index is outside its domain",
    ))
    0 <= retained_window_count <= size(window_values, 1) || throw(ArgumentError(
        "measurement snapshot retained estimator count is outside its domain",
    ))
    last_accepted_tick >= -1 || throw(ArgumentError(
        "measurement snapshot accepted tick is outside its domain",
    ))
    all(>=(0), (
        accepted_analog_step_count,
        accepted_sample_count,
        released_sample_count,
    )) || throw(ArgumentError("measurement snapshot counters must be nonnegative"))
    released_sample_count <= accepted_sample_count || throw(ArgumentError(
        "measurement snapshot released-sample count exceeds acquisitions",
    ))
    length(samples) <=
        runtime.specification.acquisition.maximum_retained_samples || throw(ArgumentError(
        "measurement snapshot retained output exceeds its declared limit",
    ))
    all(isfinite, analog_state) &&
        all(isfinite, previous_input) &&
        all(isfinite, analog_output) &&
        all(isfinite, held_values) &&
        all(isfinite, window_values) &&
        all(isfinite, squared_sums) &&
        all(isfinite, positive_sequence_history) &&
        all(isfinite, positive_sequence_time_history_s) || throw(ArgumentError(
        "measurement snapshot numerical state must be finite",
    ))
    issorted(positive_sequence_time_history_s) || throw(ArgumentError(
        "measurement snapshot frequency-history times must be monotone",
    ))
    for delayed in delayed_samples
        delayed.source_tick >= 0 && delayed.release_tick >= delayed.source_tick ||
            throw(ArgumentError("measurement snapshot delayed-sample ticks are invalid"))
        length(delayed.values) == channel_count &&
            length(delayed.clipped) == channel_count &&
            (delayed.codes === nothing || length(delayed.codes) == channel_count) ||
            throw(DimensionMismatch(
                "measurement snapshot delayed-sample channel count is incompatible",
            ))
        all(isfinite, delayed.values) || throw(ArgumentError(
            "measurement snapshot delayed-sample values must be finite",
        ))
    end
    issorted(getfield.(delayed_samples, :release_tick)) || throw(ArgumentError(
        "measurement snapshot delayed-sample releases must be monotone",
    ))
    for sample in samples
        sample.source_tick >= 0 && sample.release_tick >= sample.source_tick ||
            throw(ArgumentError("measurement snapshot retained-sample ticks are invalid"))
        length(sample.instantaneous) == channel_count &&
            length(sample.clipped) == channel_count &&
            (sample.codes === nothing || length(sample.codes) == channel_count) &&
            (sample.sliding_rms === nothing ||
                length(sample.sliding_rms) == channel_count) &&
            (sample.fundamental_rms_phasors === nothing ||
                length(sample.fundamental_rms_phasors) == channel_count) ||
            throw(DimensionMismatch(
                "measurement snapshot retained-sample channel count is incompatible",
            ))
        all(isfinite, sample.instantaneous) &&
            (sample.sliding_rms === nothing || all(isfinite, sample.sliding_rms)) &&
            (sample.fundamental_rms_phasors === nothing ||
                all(isfinite, sample.fundamental_rms_phasors)) &&
            (sample.frequency_hz === nothing || isfinite(sample.frequency_hz)) ||
            throw(ArgumentError(
                "measurement snapshot retained-sample values must be finite",
            ))
    end
    return (
        analog_state=analog_state,
        previous_input=previous_input,
        analog_output=analog_output,
        delayed_samples=delayed_samples,
        held_values=held_values,
        held_codes=held_codes,
        window_values=window_values,
        window_ticks=window_ticks,
        window_next_index=window_next_index,
        retained_window_count=retained_window_count,
        squared_sums=squared_sums,
        positive_sequence_history=positive_sequence_history,
        positive_sequence_time_history_s=positive_sequence_time_history_s,
        samples=samples,
        last_accepted_tick=last_accepted_tick,
        accepted_analog_step_count=accepted_analog_step_count,
        accepted_sample_count=accepted_sample_count,
        released_sample_count=released_sample_count,
    )
end

function restore_measurement_chain_snapshot!(
    runtime::MeasurementChainRuntime,
    snapshot::MeasurementChainSnapshot,
)
    !runtime.candidate_active || throw(ArgumentError(
        "measurement snapshot cannot be restored during an active analog trial",
    ))
    snapshot.schema_version == 1 || throw(ArgumentError(
        "measurement snapshot schema version is unsupported",
    ))
    snapshot.specification_signature_sha256 == runtime.specification_signature_sha256 ||
        throw(ArgumentError("measurement snapshot specification identity is stale"))
    expected_signature = _measurement_snapshot_signature(
        snapshot.specification_signature_sha256,
        snapshot.state,
    )
    snapshot.deterministic_signature_sha256 == expected_signature || throw(ArgumentError(
        "measurement snapshot integrity signature does not match",
    ))
    state = _validated_measurement_snapshot_state(runtime, snapshot.state)
    runtime.analog_state .= state.analog_state
    runtime.previous_input .= state.previous_input
    runtime.analog_output .= state.analog_output
    runtime.delayed_samples = deepcopy(state.delayed_samples)
    runtime.delayed_sample_head = 1
    runtime.held_values .= state.held_values
    runtime.held_codes = state.held_codes === nothing ? nothing : copy(state.held_codes)
    runtime.window_values .= state.window_values
    runtime.window_ticks .= state.window_ticks
    runtime.window_next_index = state.window_next_index
    runtime.retained_window_count = state.retained_window_count
    runtime.squared_sums .= state.squared_sums
    runtime.positive_sequence_history = copy(state.positive_sequence_history)
    runtime.positive_sequence_time_history_s = copy(
        state.positive_sequence_time_history_s,
    )
    runtime.samples = deepcopy(state.samples)
    runtime.last_accepted_tick = state.last_accepted_tick
    runtime.accepted_analog_step_count = state.accepted_analog_step_count
    runtime.accepted_sample_count = state.accepted_sample_count
    runtime.released_sample_count = state.released_sample_count
    return runtime
end

function measurement_chain_result_signature(runtime::MeasurementChainRuntime)
    !runtime.candidate_active || throw(ArgumentError(
        "measurement result signature is unavailable during an active analog trial",
    ))
    state = _measurement_runtime_state(runtime)
    return _measurement_snapshot_signature(runtime.specification_signature_sha256, state)
end

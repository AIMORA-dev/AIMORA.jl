export ElectronicSensorCouplingMode,
       NonloadingElectronicObservation,
       VoltageShuntElectronicLoading,
       CurrentSeriesElectronicInsertion,
       ElectronicSensorDefinition,
       ElectronicSensorRuntime,
       ElectronicSensorOutput,
       ElectronicSensorSnapshot,
       electronic_sensor_loading_branches,
       electronic_sensor_runtime,
       prepare_electronic_sensor_step!,
       discard_electronic_sensor_step!,
       accept_electronic_sensor_step!,
       electronic_sensor_output,
       electronic_sensor_snapshot,
       restore_electronic_sensor_snapshot!,
       electronic_sensor_signature

@enum ElectronicSensorCouplingMode begin
    NonloadingElectronicObservation
    VoltageShuntElectronicLoading
    CurrentSeriesElectronicInsertion
end

struct ElectronicSensorDefinition
    id::Symbol
    family::MeasurementProductFamily
    coupling::ElectronicSensorCouplingMode
    transducer::AnalogMeasurementStateSpace
    input_scale::Float64
    input_offset::Float64
    minimum_observed_input::Float64
    maximum_observed_input::Float64
    loading_resistance_ohm::Float64
    loading_inductance_h::Float64
    loading_capacitance_f::Float64
    provenance::ParameterProvenance
    deterministic_signature_sha256::String
end

function ElectronicSensorDefinition(
    id::Symbol,
    family::MeasurementProductFamily;
    coupling::ElectronicSensorCouplingMode,
    transducer::AnalogMeasurementStateSpace,
    input_scale::Real,
    input_offset::Real=0.0,
    minimum_observed_input::Real,
    maximum_observed_input::Real,
    loading_resistance_ohm::Real=0.0,
    loading_inductance_h::Real=0.0,
    loading_capacitance_f::Real=0.0,
    provenance::ParameterProvenance,
)
    family in (
        ElectronicCurrentSensorMeasurement,
        ElectronicVoltageSensorMeasurement,
    ) || throw(ArgumentError(
        "electronic sensor definition requires the electronic current or voltage family",
    ))
    isempty(String(id)) && throw(ArgumentError("electronic sensor id must not be empty"))
    family === ElectronicCurrentSensorMeasurement &&
        coupling === VoltageShuntElectronicLoading && throw(ArgumentError(
            "electronic current sensors cannot silently use voltage-shunt loading",
        ))
    family === ElectronicVoltageSensorMeasurement &&
        coupling === CurrentSeriesElectronicInsertion && throw(ArgumentError(
            "electronic voltage sensors cannot silently use current-series insertion",
        ))
    scale = Float64(input_scale)
    offset = Float64(input_offset)
    limits = Float64.((minimum_observed_input, maximum_observed_input))
    loading = Float64.((
        loading_resistance_ohm,
        loading_inductance_h,
        loading_capacitance_f,
    ))
    isfinite(scale) && scale != 0.0 && isfinite(offset) || throw(ArgumentError(
        "electronic sensor scale must be finite and nonzero and offset must be finite",
    ))
    all(isfinite, limits) && limits[1] < limits[2] || throw(ArgumentError(
        "electronic sensor observed-input limits must be finite and ordered",
    ))
    all(value -> isfinite(value) && value >= 0.0, loading) || throw(ArgumentError(
        "electronic sensor loading RLC parameters must be finite and nonnegative",
    ))
    if coupling === NonloadingElectronicObservation
        all(iszero, loading) || throw(ArgumentError(
            "nonloading electronic observation cannot retain hidden insertion or shunt data",
        ))
    elseif coupling === VoltageShuntElectronicLoading
        loading[1] > 0.0 || loading[3] > 0.0 || throw(ArgumentError(
            "voltage-shunt electronic loading requires resistance or capacitance",
        ))
        loading[2] == 0.0 || throw(ArgumentError(
            "voltage-shunt electronic loading does not infer a hidden series inductance",
        ))
    else
        loading[1] > 0.0 || loading[2] > 0.0 || throw(ArgumentError(
            "current-series electronic insertion requires resistance or inductance",
        ))
        loading[3] == 0.0 || throw(ArgumentError(
            "current-series electronic insertion does not infer a hidden shunt capacitance",
        ))
    end
    io = IOBuffer()
    println(io, "aimora.electronic_sensor.v1")
    println(io, id)
    println(io, Int(family), ',', Int(coupling))
    for value in (scale, offset, limits..., loading...)
        println(io, bitstring(value))
    end
    for value in (
        transducer.state_matrix_per_s,
        transducer.input_vector_per_s,
        transducer.output_vector,
        transducer.direct_gain,
        transducer.stability_margin_per_s,
        provenance,
    )
        _write_measurement_signature_value(io, value)
        print(io, '\n')
    end
    signature = bytes2hex(sha256(take!(io)))
    return ElectronicSensorDefinition(
        id,
        family,
        coupling,
        transducer,
        scale,
        offset,
        limits...,
        loading...,
        provenance,
        signature,
    )
end

function electronic_sensor_loading_branches(
    definition::ElectronicSensorDefinition,
    positive_node::Integer,
    negative_node::Integer,
)
    definition.coupling === NonloadingElectronicObservation &&
        return MeasurementBurdenBranchDefinition[]
    branches = MeasurementBurdenBranchDefinition[]
    if definition.coupling === CurrentSeriesElectronicInsertion
        if definition.loading_inductance_h > 0.0
            push!(
                branches,
                MeasurementBurdenBranchDefinition(
                    :series_rl,
                    positive_node,
                    negative_node,
                    definition.loading_resistance_ohm,
                    definition.loading_inductance_h,
                    0.0,
                ),
            )
        else
            push!(
                branches,
                MeasurementBurdenBranchDefinition(
                    :resistance,
                    positive_node,
                    negative_node,
                    definition.loading_resistance_ohm,
                    0.0,
                    0.0,
                ),
            )
        end
    else
        definition.loading_resistance_ohm > 0.0 && push!(
            branches,
            MeasurementBurdenBranchDefinition(
                :resistance,
                positive_node,
                negative_node,
                definition.loading_resistance_ohm,
                0.0,
                0.0,
            ),
        )
        definition.loading_capacitance_f > 0.0 && push!(
            branches,
            MeasurementBurdenBranchDefinition(
                :shunt_capacitance,
                positive_node,
                negative_node,
                0.0,
                0.0,
                definition.loading_capacitance_f,
            ),
        )
    end
    return branches
end

function _same_analog_state_space(left, right)
    return left.state_matrix_per_s == right.state_matrix_per_s &&
        left.input_vector_per_s == right.input_vector_per_s &&
        left.output_vector == right.output_vector &&
        left.direct_gain == right.direct_gain &&
        left.stability_margin_per_s == right.stability_margin_per_s
end

mutable struct ElectronicSensorRuntime{C<:MeasurementChainRuntime}
    definition::ElectronicSensorDefinition
    measurement_runtime::C
    accepted_observed_input::Float64
    candidate_observed_input::Union{Nothing,Float64}
    released_samples::Vector{MeasurementSample}
end

function _electronic_sensor_conditioned_input(definition, observed_input)
    observed = Float64(observed_input)
    isfinite(observed) || throw(ArgumentError(
        "electronic sensor observed input must be finite",
    ))
    definition.minimum_observed_input <= observed <= definition.maximum_observed_input ||
        throw(MeasurementChainRefusal(
            :electronic_sensor_input_outside_domain,
            :prepare_transducer_step,
            definition.id,
            definition.family,
            "electronic sensor observed input is outside its declared domain",
            (observed_input=observed,),
        ))
    return muladd(definition.input_scale, observed, definition.input_offset)
end

function electronic_sensor_runtime(
    definition::ElectronicSensorDefinition,
    measurement_specification::MeasurementChainSpecification;
    initial_observed_input::Real=0.0,
    initial_analog_state=nothing,
    initial_tick::Integer=0,
)
    measurement_specification.family === definition.family || throw(ArgumentError(
        "electronic sensor definition and measurement chain families must match",
    ))
    length(measurement_specification.channel_names) == 1 || throw(ArgumentError(
        "one electronic sensor runtime requires exactly one sampled output channel",
    ))
    expected_quantity = definition.family === ElectronicCurrentSensorMeasurement ?
        :current : :voltage
    measurement_specification.quantity === expected_quantity || throw(ArgumentError(
        "electronic sensor measurement quantity does not match its family",
    ))
    _same_analog_state_space(
        definition.transducer,
        measurement_specification.conditioning,
    ) || throw(ArgumentError(
        "electronic sensor transducer state must be the canonical sampled-chain conditioning state",
    ))
    observed = Float64(initial_observed_input)
    transducer_input = _electronic_sensor_conditioned_input(definition, observed)
    state_count = size(definition.transducer.state_matrix_per_s, 1)
    analog_state = initial_analog_state === nothing ?
        zeros(state_count, 1) : Matrix{Float64}(initial_analog_state)
    measurement = MeasurementChainRuntime(
        measurement_specification;
        initial_input=[transducer_input],
        initial_analog_state=analog_state,
    )
    released = initialize_measurement_chain_at_tick!(measurement, initial_tick)
    return ElectronicSensorRuntime(definition, measurement, observed, nothing, released)
end

function prepare_electronic_sensor_step!(
    runtime::ElectronicSensorRuntime,
    observed_input::Real,
    timestep_s::Real,
)
    runtime.candidate_observed_input === nothing || throw(ArgumentError(
        "electronic sensor trial is already active",
    ))
    observed = Float64(observed_input)
    transducer_input = _electronic_sensor_conditioned_input(runtime.definition, observed)
    prepare_measurement_analog_step!(
        runtime.measurement_runtime,
        transducer_input,
        timestep_s,
    )
    runtime.candidate_observed_input = observed
    return runtime
end

function discard_electronic_sensor_step!(runtime::ElectronicSensorRuntime)
    runtime.candidate_observed_input === nothing && throw(ArgumentError(
        "electronic sensor has no active trial to discard",
    ))
    discard_measurement_analog_step!(runtime.measurement_runtime)
    runtime.candidate_observed_input = nothing
    return runtime
end

function accept_electronic_sensor_step!(
    released::Vector{MeasurementSample},
    runtime::ElectronicSensorRuntime,
    accepted_tick::Integer,
)
    observed = runtime.candidate_observed_input
    observed === nothing && throw(ArgumentError(
        "electronic sensor step must be prepared before acceptance",
    ))
    released = accept_measurement_analog_step!(
        released,
        runtime.measurement_runtime,
        accepted_tick,
    )
    runtime.accepted_observed_input = observed
    runtime.candidate_observed_input = nothing
    runtime.released_samples = released
    return released
end

function accept_electronic_sensor_step!(
    runtime::ElectronicSensorRuntime,
    accepted_tick::Integer,
)
    return accept_electronic_sensor_step!(
        MeasurementSample[],
        runtime,
        accepted_tick,
    )
end

struct ElectronicSensorOutput
    id::Symbol
    family::MeasurementProductFamily
    coupling::ElectronicSensorCouplingMode
    observed_input::Float64
    transducer_input::Float64
    conditioned_output::Float64
    held_measurement::Float64
    latest_sample::Union{Nothing,MeasurementSample}
    deterministic_signature_sha256::String
end

function electronic_sensor_output(runtime::ElectronicSensorRuntime)
    runtime.candidate_observed_input === nothing || throw(ArgumentError(
        "electronic sensor output is unavailable during an active trial",
    ))
    definition = runtime.definition
    transducer_input = _electronic_sensor_conditioned_input(
        definition,
        runtime.accepted_observed_input,
    )
    latest = isempty(runtime.measurement_runtime.samples) ?
        nothing : last(runtime.measurement_runtime.samples)
    io = IOBuffer()
    println(io, definition.deterministic_signature_sha256)
    println(io, measurement_chain_result_signature(runtime.measurement_runtime))
    println(io, bitstring(runtime.accepted_observed_input))
    signature = bytes2hex(sha256(take!(io)))
    return ElectronicSensorOutput(
        definition.id,
        definition.family,
        definition.coupling,
        runtime.accepted_observed_input,
        transducer_input,
        only(runtime.measurement_runtime.analog_output),
        only(runtime.measurement_runtime.held_values),
        latest,
        signature,
    )
end

struct ElectronicSensorSnapshot
    schema_version::Int
    definition_signature_sha256::String
    accepted_observed_input::Float64
    measurement_snapshot::MeasurementChainSnapshot
    deterministic_signature_sha256::String
end

function electronic_sensor_snapshot(runtime::ElectronicSensorRuntime)
    runtime.candidate_observed_input === nothing || throw(ArgumentError(
        "electronic sensor snapshot is unavailable during an active trial",
    ))
    measurement_snapshot = measurement_chain_snapshot(runtime.measurement_runtime)
    signature = bytes2hex(sha256(
        runtime.definition.deterministic_signature_sha256 *
        bitstring(runtime.accepted_observed_input) *
        measurement_snapshot.deterministic_signature_sha256,
    ))
    return ElectronicSensorSnapshot(
        1,
        runtime.definition.deterministic_signature_sha256,
        runtime.accepted_observed_input,
        measurement_snapshot,
        signature,
    )
end

function restore_electronic_sensor_snapshot!(
    runtime::ElectronicSensorRuntime,
    snapshot::ElectronicSensorSnapshot,
)
    runtime.candidate_observed_input === nothing || throw(ArgumentError(
        "electronic sensor snapshot cannot be restored during an active trial",
    ))
    snapshot.schema_version == 1 || throw(ArgumentError(
        "electronic sensor snapshot schema is unsupported",
    ))
    snapshot.definition_signature_sha256 ==
        runtime.definition.deterministic_signature_sha256 || throw(ArgumentError(
        "electronic sensor snapshot definition is stale",
    ))
    expected_signature = bytes2hex(sha256(
        snapshot.definition_signature_sha256 *
        bitstring(snapshot.accepted_observed_input) *
        snapshot.measurement_snapshot.deterministic_signature_sha256,
    ))
    snapshot.deterministic_signature_sha256 == expected_signature || throw(
        ArgumentError("electronic sensor snapshot integrity failed"),
    )
    measurement_probe = deepcopy(runtime.measurement_runtime)
    restore_measurement_chain_snapshot!(
        measurement_probe,
        snapshot.measurement_snapshot,
    )
    restore_measurement_chain_snapshot!(
        runtime.measurement_runtime,
        snapshot.measurement_snapshot,
    )
    runtime.accepted_observed_input = snapshot.accepted_observed_input
    empty!(runtime.released_samples)
    return runtime
end

electronic_sensor_signature(runtime::ElectronicSensorRuntime) =
    electronic_sensor_snapshot(runtime).deterministic_signature_sha256

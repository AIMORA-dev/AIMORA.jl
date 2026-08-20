export MeasurementBurden,
       MeasurementBurdenBranchDefinition,
       measurement_burden_impedance,
       measurement_burden_branches,
       InstrumentTransformerMeasurementDefinition,
       InstrumentTransformerMeasurementReadiness,
       InstrumentTransformerMeasurementRuntime,
       InstrumentTransformerMeasurementOutput,
       InstrumentTransformerMeasurementSnapshot,
       instrument_transformer_measurement_readiness,
       instrument_transformer_measurement_runtime,
       instrument_transformer_measurement_output,
       instrument_transformer_measurement_snapshot,
       restore_instrument_transformer_measurement_snapshot!,
       instrument_transformer_measurement_signature

struct MeasurementBurden
    connected::Bool
    series_resistance_ohm::Float64
    series_inductance_h::Float64
    shunt_capacitance_f::Float64
    cable_resistance_ohm::Float64
    cable_inductance_h::Float64
    cable_capacitance_f::Float64
    provenance::ParameterProvenance

    function MeasurementBurden(;
        connected::Bool=true,
        series_resistance_ohm::Real=0.0,
        series_inductance_h::Real=0.0,
        shunt_capacitance_f::Real=0.0,
        cable_resistance_ohm::Real=0.0,
        cable_inductance_h::Real=0.0,
        cable_capacitance_f::Real=0.0,
        provenance::ParameterProvenance,
    )
        values = Float64.((
            series_resistance_ohm,
            series_inductance_h,
            shunt_capacitance_f,
            cable_resistance_ohm,
            cable_inductance_h,
            cable_capacitance_f,
        ))
        all(value -> isfinite(value) && value >= 0.0, values) || throw(ArgumentError(
            "measurement burden RLC and cable parameters must be finite and nonnegative",
        ))
        if connected
            any(>(0.0), values) || throw(ArgumentError(
                "a connected measurement burden must contain a finite passive R, L, or C owner",
            ))
        else
            all(iszero, values) || throw(ArgumentError(
                "a disconnected measurement burden cannot retain hidden connected parameters",
            ))
        end
        return new(connected, values..., provenance)
    end
end

struct MeasurementBurdenBranchDefinition
    kind::Symbol
    positive_node::Int
    negative_node::Int
    resistance_ohm::Float64
    inductance_h::Float64
    capacitance_f::Float64

    function MeasurementBurdenBranchDefinition(
        kind::Symbol,
        positive_node::Integer,
        negative_node::Integer,
        resistance_ohm::Real,
        inductance_h::Real,
        capacitance_f::Real,
    )
        kind in (:series_rl, :resistance, :shunt_capacitance) || throw(
            ArgumentError("unknown measurement-burden branch kind"),
        )
        positive = Int(positive_node)
        negative = Int(negative_node)
        positive >= 0 && negative >= 0 && positive != negative || throw(ArgumentError(
            "measurement-burden nodes must be distinct and nonnegative",
        ))
        values = Float64.((resistance_ohm, inductance_h, capacitance_f))
        all(value -> isfinite(value) && value >= 0.0, values) || throw(ArgumentError(
            "measurement-burden branch values must be finite and nonnegative",
        ))
        if kind === :series_rl
            values[2] > 0.0 && values[3] == 0.0 || throw(ArgumentError(
                "series-RL burden branch requires positive inductance and zero capacitance",
            ))
        elseif kind === :resistance
            values[1] > 0.0 && values[2] == 0.0 && values[3] == 0.0 || throw(
                ArgumentError("resistive burden branch requires resistance only"),
            )
        else
            values[1] == 0.0 && values[2] == 0.0 && values[3] > 0.0 || throw(
                ArgumentError("shunt-capacitance burden branch requires capacitance only"),
            )
        end
        return new(kind, positive, negative, values...)
    end
end

function measurement_burden_impedance(burden::MeasurementBurden, frequency_hz::Real)
    frequency = Float64(frequency_hz)
    isfinite(frequency) && frequency >= 0.0 || throw(ArgumentError(
        "measurement burden frequency must be finite and nonnegative",
    ))
    burden.connected || return ComplexF64(Inf, 0.0)
    resistance = burden.series_resistance_ohm + burden.cable_resistance_ohm
    inductance = burden.series_inductance_h + burden.cable_inductance_h
    capacitance = burden.shunt_capacitance_f + burden.cable_capacitance_f
    angular_frequency = 2.0 * pi * frequency
    series_impedance = ComplexF64(resistance, angular_frequency * inductance)
    if capacitance == 0.0
        return series_impedance
    elseif iszero(series_impedance)
        return zero(ComplexF64)
    end
    return inv(inv(series_impedance) + im * angular_frequency * capacitance)
end

function measurement_burden_branches(
    burden::MeasurementBurden,
    positive_node::Integer,
    negative_node::Integer,
)
    burden.connected || return MeasurementBurdenBranchDefinition[]
    resistance = burden.series_resistance_ohm + burden.cable_resistance_ohm
    inductance = burden.series_inductance_h + burden.cable_inductance_h
    capacitance = burden.shunt_capacitance_f + burden.cable_capacitance_f
    branches = MeasurementBurdenBranchDefinition[]
    if inductance > 0.0
        push!(
            branches,
            MeasurementBurdenBranchDefinition(
                :series_rl,
                Int(positive_node),
                Int(negative_node),
                resistance,
                inductance,
                0.0,
            ),
        )
    elseif resistance > 0.0
        push!(
            branches,
            MeasurementBurdenBranchDefinition(
                :resistance,
                Int(positive_node),
                Int(negative_node),
                resistance,
                0.0,
                0.0,
            ),
        )
    end
    capacitance > 0.0 && push!(
        branches,
        MeasurementBurdenBranchDefinition(
            :shunt_capacitance,
            Int(positive_node),
            Int(negative_node),
            0.0,
            0.0,
            capacitance,
        ),
    )
    return branches
end

struct InstrumentTransformerMeasurementDefinition{S<:TransformerApparatusSpecification}
    id::Symbol
    family::MeasurementProductFamily
    apparatus::S
    primary_coil_index::Int
    secondary_coil_index::Int
    primary_turns::Float64
    secondary_turns::Float64
    burden::MeasurementBurden
    maximum_secondary_impedance_ohm::Float64
    secondary_output_sign::Float64
    deterministic_signature_sha256::String
end

function _instrument_transformer_magnetic_graph(model)
    model isa HybridTransformerModel && return model.magnetic_graph
    model isa MagneticEquivalentCircuitModel && return model.magnetic_graph
    return nothing
end

function InstrumentTransformerMeasurementDefinition(
    id::Symbol,
    family::MeasurementProductFamily,
    apparatus::S;
    primary_coil_index::Integer,
    secondary_coil_index::Integer,
    primary_turns::Real,
    secondary_turns::Real,
    burden::MeasurementBurden,
    maximum_secondary_impedance_ohm::Real,
    secondary_output_sign::Real,
) where {S<:TransformerApparatusSpecification}
    family in (
        LinearCurrentTransformerMeasurement,
        MagneticCurrentTransformerMeasurement,
        InductiveVoltageTransformerMeasurement,
        CouplingCapacitorVoltageTransformerMeasurement,
    ) || throw(ArgumentError(
        "instrument-transformer definition requires a CT, inductive VT, or CVT electromagnetic-unit family",
    ))
    apparatus.tier in (
        LowFrequencyTerminalTier,
        BCTRANTerminalTier,
        HybridTransformerTier,
        MagneticEquivalentCircuitTier,
    ) || throw(ArgumentError(
        "instrument-transformer measurement requires a tier with explicit coil current and voltage",
    ))
    isempty(String(id)) && throw(ArgumentError(
        "instrument-transformer measurement id must not be empty",
    ))
    coil_count = length(apparatus.connection.coil_order)
    primary_index = Int(primary_coil_index)
    secondary_index = Int(secondary_coil_index)
    1 <= primary_index <= coil_count && 1 <= secondary_index <= coil_count ||
        throw(BoundsError(apparatus.connection.coil_order, (primary_index, secondary_index)))
    primary_index != secondary_index || throw(ArgumentError(
        "instrument-transformer primary and secondary coils must be distinct",
    ))
    turns = Float64.((primary_turns, secondary_turns))
    all(value -> isfinite(value) && value > 0.0, turns) || throw(ArgumentError(
        "instrument-transformer primary and secondary turns must be finite and positive",
    ))
    maximum_impedance = Float64(maximum_secondary_impedance_ohm)
    isfinite(maximum_impedance) && maximum_impedance > 0.0 || throw(ArgumentError(
        "instrument-transformer maximum secondary impedance must be finite and positive",
    ))
    output_sign = Float64(secondary_output_sign)
    output_sign in (-1.0, 1.0) || throw(ArgumentError(
        "instrument-transformer secondary output sign must be exactly minus or plus one",
    ))
    graph = _instrument_transformer_magnetic_graph(apparatus.model)
    if family === MagneticCurrentTransformerMeasurement
        graph !== nothing || throw(ArgumentError(
            "saturating/remanent CT requires the canonical transformer magnetic graph",
        ))
        any(material -> material isa TellinenTransformerMagneticMaterial, graph.materials) ||
            throw(ArgumentError(
                "saturating/remanent CT requires an explicit Tellinen hysteresis material",
            ))
    elseif family === LinearCurrentTransformerMeasurement && graph !== nothing
        all(material -> material isa LinearTransformerMagneticMaterial, graph.materials) ||
            throw(ArgumentError(
                "linear CT cannot silently contain nonlinear or hysteretic magnetic material",
            ))
    end
    if graph !== nothing
        graph_turns = graph.winding_turns
        primary_graph_turns = maximum(abs, view(graph_turns, :, primary_index))
        secondary_graph_turns = maximum(abs, view(graph_turns, :, secondary_index))
        tolerance = 64.0 * eps(Float64) * max(turns..., 1.0)
        abs(primary_graph_turns - turns[1]) <= tolerance || throw(ArgumentError(
            "instrument-transformer primary turns disagree with its magnetic graph",
        ))
        abs(secondary_graph_turns - turns[2]) <= tolerance || throw(ArgumentError(
            "instrument-transformer secondary turns disagree with its magnetic graph",
        ))
    end
    io = IOBuffer()
    println(io, "aimora.instrument_transformer_measurement.v1")
    println(io, id)
    println(io, Int(family))
    println(io, apparatus.deterministic_signature_sha256)
    println(io, primary_index, ',', secondary_index)
    println(io, bitstring(turns[1]), ',', bitstring(turns[2]))
    println(io, burden.connected)
    for value in (
        burden.series_resistance_ohm,
        burden.series_inductance_h,
        burden.shunt_capacitance_f,
        burden.cable_resistance_ohm,
        burden.cable_inductance_h,
        burden.cable_capacitance_f,
        maximum_impedance,
        output_sign,
    )
        println(io, bitstring(value))
    end
    for field in fieldnames(ParameterProvenance)
        println(io, getfield(burden.provenance, field))
    end
    signature = bytes2hex(sha256(take!(io)))
    return InstrumentTransformerMeasurementDefinition{S}(
        id,
        family,
        apparatus,
        primary_index,
        secondary_index,
        turns...,
        burden,
        maximum_impedance,
        output_sign,
        signature,
    )
end

struct InstrumentTransformerMeasurementReadiness
    ready::Bool
    reasons::Tuple
    burden_impedance_ohm::ComplexF64
    deterministic_signature_sha256::String
end

function instrument_transformer_measurement_readiness(
    definition::InstrumentTransformerMeasurementDefinition,
)
    burden_impedance = measurement_burden_impedance(
        definition.burden,
        definition.apparatus.rated_frequency_hz,
    )
    reasons = Symbol[]
    if definition.family in (
        LinearCurrentTransformerMeasurement,
        MagneticCurrentTransformerMeasurement,
    )
        definition.burden.connected || push!(reasons, :current_transformer_secondary_open)
        abs(burden_impedance) <= definition.maximum_secondary_impedance_ohm ||
            push!(reasons, :current_transformer_burden_outside_domain)
    end
    return InstrumentTransformerMeasurementReadiness(
        isempty(reasons),
        Tuple(reasons),
        burden_impedance,
        definition.deterministic_signature_sha256,
    )
end

mutable struct InstrumentTransformerMeasurementRuntime{D,A,C} <:
               AbstractNonlinearCurrentDevice
    definition::D
    apparatus_runtime::A
    measurement_runtime::C
    released_samples::Vector{MeasurementSample}
end

function _instrument_measurement_input(definition, coil_current, coil_voltage)
    if definition.family in (
        LinearCurrentTransformerMeasurement,
        MagneticCurrentTransformerMeasurement,
    )
        return definition.secondary_output_sign *
            coil_current[definition.secondary_coil_index]
    end
    return definition.secondary_output_sign *
        coil_voltage[definition.secondary_coil_index]
end

function _instrument_measurement_tick(definition, settings, time_s)
    raw_tick = time_s / settings.tick_s
    isfinite(raw_tick) && 0.0 <= raw_tick <= typemax(Int) || throw(
        MeasurementChainRefusal(
            :unrepresentable_measurement_time,
            :accept_physical_state,
            definition.id,
            definition.family,
            "physical accepted time is outside the representable measurement tick domain",
            (time_s=time_s, tick_s=settings.tick_s),
        ),
    )
    tick = round(Int, raw_tick)
    tolerance = 64.0 * eps(Float64) * max(abs(time_s), settings.tick_s, 1.0)
    abs(time_s - tick * settings.tick_s) <= tolerance || throw(ArgumentError(
        "instrument-transformer accepted time is not exactly representable on its measurement clock",
    ))
    return tick
end

function instrument_transformer_measurement_runtime(
    definition::InstrumentTransformerMeasurementDefinition,
    preparation::TransformerApparatusPreparation,
    terminal_nodes,
    measurement_specification::MeasurementChainSpecification,
)
    preparation.specification.deterministic_signature_sha256 ==
        definition.apparatus.deterministic_signature_sha256 || throw(ArgumentError(
        "instrument-transformer preparation does not match its measurement definition",
    ))
    measurement_specification.family === definition.family || throw(ArgumentError(
        "instrument-transformer physical and sampled-chain families must match",
    ))
    length(measurement_specification.channel_names) == 1 || throw(ArgumentError(
        "one instrument-transformer runtime requires exactly one sampled output channel",
    ))
    expected_quantity = definition.family in (
        InductiveVoltageTransformerMeasurement,
        CouplingCapacitorVoltageTransformerMeasurement,
    ) ? :voltage : :current
    measurement_specification.quantity === expected_quantity || throw(ArgumentError(
        "instrument-transformer sampled quantity does not match its physical family",
    ))
    settings = measurement_specification.acquisition
    physical_step = definition.apparatus.settings.timestep_s
    tolerance = 64.0 * eps(Float64) * max(physical_step, settings.tick_s, 1.0)
    abs(physical_step - settings.tick_s) <= tolerance || throw(ArgumentError(
        "instrument-transformer coupling currently requires one exact measurement tick per physical step",
    ))
    readiness = instrument_transformer_measurement_readiness(definition)
    readiness.ready || throw(MeasurementChainRefusal(
        first(readiness.reasons),
        :prepare_physical_measurement,
        definition.id,
        definition.family,
        "instrument-transformer physical measurement is outside its burden or apparatus domain",
        (reasons=readiness.reasons, burden_impedance_ohm=readiness.burden_impedance_ohm),
    ))
    apparatus = transformer_apparatus_runtime(preparation, terminal_nodes)
    state = apparatus.accepted_state
    initial_input = _instrument_measurement_input(
        definition,
        state.coil_current_a,
        state.coil_voltage_v,
    )
    measurement = MeasurementChainRuntime(
        measurement_specification;
        initial_input=[initial_input],
    )
    initial_tick = _instrument_measurement_tick(
        definition,
        settings,
        preparation.initial_time_s,
    )
    released = initialize_measurement_chain_at_tick!(measurement, initial_tick)
    return InstrumentTransformerMeasurementRuntime(
        definition,
        apparatus,
        measurement,
        released,
    )
end

nonlinear_terminal_nodes(runtime::InstrumentTransformerMeasurementRuntime) =
    nonlinear_terminal_nodes(runtime.apparatus_runtime)

nonlinear_device_formulation(::InstrumentTransformerMeasurementRuntime) =
    PhysicalConstitutiveCurrent

nonlinear_device_provenance(runtime::InstrumentTransformerMeasurementRuntime) =
    nonlinear_device_provenance(runtime.apparatus_runtime)

function prepare_nonlinear_device_step!(
    runtime::InstrumentTransformerMeasurementRuntime,
    time_s::Float64,
    step_s::Float64,
    companion_method::Symbol,
)
    !runtime.measurement_runtime.candidate_active || throw(ArgumentError(
        "instrument-transformer measurement retained an unfinished analog trial",
    ))
    return prepare_nonlinear_device_step!(
        runtime.apparatus_runtime,
        time_s,
        step_s,
        companion_method,
    )
end

function nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    runtime::InstrumentTransformerMeasurementRuntime,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::Float64,
)
    return nonlinear_current_jacobian!(
        terminal_current_a,
        terminal_jacobian_s,
        runtime.apparatus_runtime,
        terminal_voltage_v,
        time_s,
    )
end

function _instrument_candidate_input(runtime::InstrumentTransformerMeasurementRuntime)
    candidate = runtime.apparatus_runtime.candidate
    hasproperty(candidate, :coil_current_a) && hasproperty(candidate, :coil_voltage_v) ||
        throw(ArgumentError(
            "selected transformer apparatus tier does not expose physical winding quantities",
        ))
    return _instrument_measurement_input(
        runtime.definition,
        candidate.coil_current_a,
        candidate.coil_voltage_v,
    )
end

function accept_nonlinear_device_state!(
    runtime::InstrumentTransformerMeasurementRuntime,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    time_s::Float64,
)
    tick = _instrument_measurement_tick(
        runtime.definition,
        runtime.measurement_runtime.specification.acquisition,
        time_s,
    )
    prepare_measurement_analog_step!(
        runtime.measurement_runtime,
        _instrument_candidate_input(runtime),
        runtime.apparatus_runtime.candidate_step_s,
    )
    try
        _validated_measurement_acceptance_tick(runtime.measurement_runtime, tick)
        accept_nonlinear_device_state!(
            runtime.apparatus_runtime,
            terminal_voltage_v,
            terminal_current_a,
            terminal_jacobian_s,
            time_s,
        )
        released = accept_measurement_analog_step!(runtime.measurement_runtime, tick)
        runtime.released_samples = released
    catch
        if runtime.measurement_runtime.candidate_active
            discard_measurement_analog_step!(runtime.measurement_runtime)
        end
        rethrow()
    end
    return nothing
end

function accept_nonlinear_device_state!(
    runtime::InstrumentTransformerMeasurementRuntime,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    time_s::Float64,
)
    candidate = runtime.apparatus_runtime.candidate
    hasproperty(candidate, :terminal_jacobian_s) || throw(ArgumentError(
        "instrument-transformer candidate does not expose its converged terminal Jacobian",
    ))
    return accept_nonlinear_device_state!(
        runtime,
        terminal_voltage_v,
        terminal_current_a,
        candidate.terminal_jacobian_s,
        time_s,
    )
end

function finish_nonlinear_device_step!(runtime::InstrumentTransformerMeasurementRuntime)
    finish_nonlinear_device_step!(runtime.apparatus_runtime)
    return runtime
end

function nonlinear_device_event_surfaces(runtime::InstrumentTransformerMeasurementRuntime)
    return Tuple(
        NonlinearDeviceEventSurface(
            surface.name,
            (wrapper, time_s) -> nonlinear_device_event_value(
                surface,
                wrapper.apparatus_runtime,
                time_s,
            ),
            (wrapper, time_s) -> apply_nonlinear_device_event!(
                surface,
                wrapper.apparatus_runtime,
                time_s,
            );
            direction=surface.direction,
            priority=surface.priority,
            topology_invalidating=surface.topology_invalidating,
            candidate_time=wrapper -> nonlinear_device_event_candidate_time(
                surface,
                wrapper.apparatus_runtime,
            ),
        ) for surface in nonlinear_device_event_surfaces(runtime.apparatus_runtime)
    )
end

struct InstrumentTransformerMeasurementOutput
    id::Symbol
    family::MeasurementProductFamily
    accepted_time_s::Float64
    primary_current_a::Float64
    primary_voltage_v::Float64
    secondary_current_a::Float64
    secondary_voltage_v::Float64
    primary_referred_measurement::Float64
    ampere_turn_residual_at::Float64
    terminal_power_w::Float64
    held_measurement::Float64
    latest_sample::Union{Nothing,MeasurementSample}
    deterministic_signature_sha256::String
end

function instrument_transformer_measurement_output(
    runtime::InstrumentTransformerMeasurementRuntime,
)
    runtime.apparatus_runtime.prepared && throw(ArgumentError(
        "instrument-transformer output is unavailable during an active physical trial",
    ))
    state = runtime.apparatus_runtime.accepted_state
    definition = runtime.definition
    primary_current = state.coil_current_a[definition.primary_coil_index]
    secondary_current = state.coil_current_a[definition.secondary_coil_index]
    primary_voltage = state.coil_voltage_v[definition.primary_coil_index]
    secondary_voltage = state.coil_voltage_v[definition.secondary_coil_index]
    referred = if definition.family in (
        LinearCurrentTransformerMeasurement,
        MagneticCurrentTransformerMeasurement,
    )
        -definition.secondary_turns / definition.primary_turns * secondary_current
    else
        definition.primary_turns / definition.secondary_turns * secondary_voltage
    end
    ampere_turn_residual = definition.primary_turns * primary_current +
        definition.secondary_turns * secondary_current
    latest_sample = isempty(runtime.measurement_runtime.samples) ?
        nothing : last(runtime.measurement_runtime.samples)
    signature_io = IOBuffer()
    println(signature_io, "aimora.instrument_transformer_measurement_output.v1")
    println(signature_io, definition.deterministic_signature_sha256)
    println(
        signature_io,
        measurement_chain_result_signature(runtime.measurement_runtime),
    )
    for value in (
        state.accepted_time_s,
        primary_current,
        primary_voltage,
        secondary_current,
        secondary_voltage,
        referred,
        ampere_turn_residual,
        state.terminal_power_w,
        only(runtime.measurement_runtime.held_values),
    )
        println(signature_io, bitstring(value))
    end
    signature = bytes2hex(sha256(take!(signature_io)))
    return InstrumentTransformerMeasurementOutput(
        definition.id,
        definition.family,
        state.accepted_time_s,
        primary_current,
        primary_voltage,
        secondary_current,
        secondary_voltage,
        referred,
        ampere_turn_residual,
        state.terminal_power_w,
        only(runtime.measurement_runtime.held_values),
        latest_sample,
        signature,
    )
end

struct InstrumentTransformerMeasurementSnapshot
    schema_version::Int
    definition_signature_sha256::String
    apparatus_snapshot::TransformerApparatusRuntimeSnapshot
    measurement_snapshot::MeasurementChainSnapshot
    deterministic_signature_sha256::String
end

function instrument_transformer_measurement_snapshot(
    runtime::InstrumentTransformerMeasurementRuntime,
)
    apparatus_snapshot = transformer_apparatus_runtime_snapshot(
        runtime.apparatus_runtime,
    )
    measurement_snapshot = measurement_chain_snapshot(runtime.measurement_runtime)
    signature = bytes2hex(sha256(
        runtime.definition.deterministic_signature_sha256 *
        apparatus_snapshot.deterministic_signature_sha256 *
        measurement_snapshot.deterministic_signature_sha256,
    ))
    return InstrumentTransformerMeasurementSnapshot(
        1,
        runtime.definition.deterministic_signature_sha256,
        apparatus_snapshot,
        measurement_snapshot,
        signature,
    )
end

function restore_instrument_transformer_measurement_snapshot!(
    runtime::InstrumentTransformerMeasurementRuntime,
    snapshot::InstrumentTransformerMeasurementSnapshot,
)
    snapshot.schema_version == 1 || throw(ArgumentError(
        "instrument-transformer measurement snapshot schema is unsupported",
    ))
    snapshot.definition_signature_sha256 ==
        runtime.definition.deterministic_signature_sha256 || throw(ArgumentError(
        "instrument-transformer measurement snapshot definition is stale",
    ))
    expected_signature = bytes2hex(sha256(
        snapshot.definition_signature_sha256 *
        snapshot.apparatus_snapshot.deterministic_signature_sha256 *
        snapshot.measurement_snapshot.deterministic_signature_sha256,
    ))
    snapshot.deterministic_signature_sha256 == expected_signature || throw(
        ArgumentError("instrument-transformer measurement snapshot integrity failed"),
    )
    apparatus_probe = deepcopy(runtime.apparatus_runtime)
    measurement_probe = deepcopy(runtime.measurement_runtime)
    restore_transformer_apparatus_runtime_snapshot!(
        apparatus_probe,
        snapshot.apparatus_snapshot,
    )
    restore_measurement_chain_snapshot!(
        measurement_probe,
        snapshot.measurement_snapshot,
    )
    restore_transformer_apparatus_runtime_snapshot!(
        runtime.apparatus_runtime,
        snapshot.apparatus_snapshot,
    )
    restore_measurement_chain_snapshot!(
        runtime.measurement_runtime,
        snapshot.measurement_snapshot,
    )
    empty!(runtime.released_samples)
    return runtime
end

function instrument_transformer_measurement_signature(
    runtime::InstrumentTransformerMeasurementRuntime,
)
    snapshot = instrument_transformer_measurement_snapshot(runtime)
    return snapshot.deterministic_signature_sha256
end

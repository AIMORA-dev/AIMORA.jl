export CvtNetworkBranchDefinition,
       CouplingCapacitorVoltageTransformerDefinition,
       CouplingCapacitorVoltageTransformerReadiness,
       CouplingCapacitorVoltageTransformerOutput,
       cvt_divider_ratio,
       cvt_series_resonance_hz,
       cvt_network_branches,
       cvt_measurement_readiness,
       cvt_measurement_runtime,
       cvt_measurement_output,
       cvt_stored_energy_j

struct CvtNetworkBranchDefinition
    id::Symbol
    kind::Symbol
    positive_node::Int
    negative_node::Int
    resistance_ohm::Float64
    inductance_h::Float64
    capacitance_f::Float64

    function CvtNetworkBranchDefinition(
        id::Symbol,
        kind::Symbol,
        positive_node::Integer,
        negative_node::Integer;
        resistance_ohm::Real=0.0,
        inductance_h::Real=0.0,
        capacitance_f::Real=0.0,
    )
        isempty(String(id)) && throw(ArgumentError("CVT network branch id must not be empty"))
        kind in (:series_rl, :resistance, :capacitance) || throw(ArgumentError(
            "CVT network branch kind must be series-RL, resistance, or capacitance",
        ))
        positive = Int(positive_node)
        negative = Int(negative_node)
        positive >= 0 && negative >= 0 && positive != negative || throw(ArgumentError(
            "CVT network branch nodes must be distinct and nonnegative",
        ))
        values = Float64.((resistance_ohm, inductance_h, capacitance_f))
        all(value -> isfinite(value) && value >= 0.0, values) || throw(ArgumentError(
            "CVT network branch RLC values must be finite and nonnegative",
        ))
        if kind === :series_rl
            values[2] > 0.0 && values[3] == 0.0 || throw(ArgumentError(
                "CVT series-RL branch requires positive inductance and zero capacitance",
            ))
        elseif kind === :resistance
            values[1] > 0.0 && values[2] == 0.0 && values[3] == 0.0 || throw(
                ArgumentError("CVT resistance branch requires resistance only"),
            )
        else
            values[1] == 0.0 && values[2] == 0.0 && values[3] > 0.0 || throw(
                ArgumentError("CVT capacitance branch requires capacitance only"),
            )
        end
        return new(id, kind, positive, negative, values...)
    end
end

struct CouplingCapacitorVoltageTransformerDefinition{
    I<:InstrumentTransformerMeasurementDefinition,
}
    id::Symbol
    electromagnetic_unit::I
    high_voltage_capacitance_f::Float64
    intermediate_voltage_capacitance_f::Float64
    compensation_resistance_ohm::Float64
    compensation_inductance_h::Float64
    suppression_resistance_ohm::Float64
    suppression_capacitance_f::Float64
    maximum_spectral_frequency_hz::Float64
    maximum_timestep_s::Float64
    provenance::ParameterProvenance
    deterministic_signature_sha256::String
end

function cvt_divider_ratio(
    definition::CouplingCapacitorVoltageTransformerDefinition,
)
    return definition.high_voltage_capacitance_f / (
        definition.high_voltage_capacitance_f +
        definition.intermediate_voltage_capacitance_f
    )
end

function cvt_series_resonance_hz(
    definition::CouplingCapacitorVoltageTransformerDefinition,
)
    equivalent_capacitance = inv(
        inv(definition.high_voltage_capacitance_f) +
        inv(definition.intermediate_voltage_capacitance_f),
    )
    return inv(2.0 * pi * sqrt(
        definition.compensation_inductance_h * equivalent_capacitance,
    ))
end

function CouplingCapacitorVoltageTransformerDefinition(
    id::Symbol,
    electromagnetic_unit::I;
    high_voltage_capacitance_f::Real,
    intermediate_voltage_capacitance_f::Real,
    compensation_resistance_ohm::Real,
    compensation_inductance_h::Real,
    suppression_resistance_ohm::Real,
    suppression_capacitance_f::Real,
    maximum_spectral_frequency_hz::Real,
    maximum_timestep_s::Real,
    provenance::ParameterProvenance,
) where {I<:InstrumentTransformerMeasurementDefinition}
    electromagnetic_unit.family === CouplingCapacitorVoltageTransformerMeasurement ||
        throw(ArgumentError(
            "CVT electromagnetic unit must carry the CVT measurement family identity",
        ))
    isempty(String(id)) && throw(ArgumentError("CVT id must not be empty"))
    capacitances = Float64.((
        high_voltage_capacitance_f,
        intermediate_voltage_capacitance_f,
        suppression_capacitance_f,
    ))
    resistances = Float64.((
        compensation_resistance_ohm,
        suppression_resistance_ohm,
    ))
    compensation_inductance = Float64(compensation_inductance_h)
    all(value -> isfinite(value) && value > 0.0, capacitances[1:2]) || throw(
        ArgumentError("CVT divider capacitances must be finite and positive"),
    )
    isfinite(capacitances[3]) && capacitances[3] >= 0.0 || throw(ArgumentError(
        "CVT suppression capacitance must be finite and nonnegative",
    ))
    all(value -> isfinite(value) && value >= 0.0, resistances) || throw(
        ArgumentError("CVT compensation and suppression resistance must be finite and nonnegative"),
    )
    compensation_inductance > 0.0 && isfinite(compensation_inductance) || throw(
        ArgumentError("CVT compensation inductance must be finite and positive"),
    )
    resistances[2] > 0.0 || capacitances[3] > 0.0 || throw(ArgumentError(
        "CVT must declare a passive ferroresonance-suppression branch",
    ))
    maximum_frequency = Float64(maximum_spectral_frequency_hz)
    maximum_timestep = Float64(maximum_timestep_s)
    isfinite(maximum_frequency) && maximum_frequency > 0.0 || throw(ArgumentError(
        "CVT maximum spectral frequency must be finite and positive",
    ))
    isfinite(maximum_timestep) && maximum_timestep > 0.0 || throw(ArgumentError(
        "CVT maximum timestep must be finite and positive",
    ))
    electromagnetic_unit.apparatus.settings.timestep_s <= maximum_timestep || throw(
        ArgumentError("CVT electromagnetic-unit timestep exceeds the CVT domain"),
    )
    provisional = CouplingCapacitorVoltageTransformerDefinition{I}(
        id,
        electromagnetic_unit,
        capacitances[1],
        capacitances[2],
        resistances[1],
        compensation_inductance,
        resistances[2],
        capacitances[3],
        maximum_frequency,
        maximum_timestep,
        provenance,
        "",
    )
    resonance = cvt_series_resonance_hz(provisional)
    resonance <= maximum_frequency || throw(ArgumentError(
        "CVT internal series resonance lies outside its registered spectral band",
    ))
    maximum_timestep * maximum_frequency <= 0.05 || throw(ArgumentError(
        "CVT timestep does not provide at least twenty steps per highest registered cycle",
    ))
    io = IOBuffer()
    println(io, "aimora.cvt_measurement.v1")
    println(io, id)
    println(io, electromagnetic_unit.deterministic_signature_sha256)
    for value in (
        capacitances...,
        resistances...,
        compensation_inductance,
        maximum_frequency,
        maximum_timestep,
    )
        println(io, bitstring(value))
    end
    _write_measurement_signature_value(io, provenance)
    signature = bytes2hex(sha256(take!(io)))
    return CouplingCapacitorVoltageTransformerDefinition{I}(
        id,
        electromagnetic_unit,
        capacitances[1],
        capacitances[2],
        resistances[1],
        compensation_inductance,
        resistances[2],
        capacitances[3],
        maximum_frequency,
        maximum_timestep,
        provenance,
        signature,
    )
end

function cvt_network_branches(
    definition::CouplingCapacitorVoltageTransformerDefinition;
    line_node::Integer,
    divider_node::Integer,
    electromagnetic_primary_node::Integer,
    secondary_node::Integer,
    reference_node::Integer=0,
)
    branches = CvtNetworkBranchDefinition[
        CvtNetworkBranchDefinition(
            :high_voltage_coupling_capacitor,
            :capacitance,
            line_node,
            divider_node;
            capacitance_f=definition.high_voltage_capacitance_f,
        ),
        CvtNetworkBranchDefinition(
            :intermediate_voltage_divider_capacitor,
            :capacitance,
            divider_node,
            reference_node;
            capacitance_f=definition.intermediate_voltage_capacitance_f,
        ),
        CvtNetworkBranchDefinition(
            :compensation_reactor,
            :series_rl,
            divider_node,
            electromagnetic_primary_node;
            resistance_ohm=definition.compensation_resistance_ohm,
            inductance_h=definition.compensation_inductance_h,
        ),
    ]
    definition.suppression_resistance_ohm > 0.0 && push!(
        branches,
        CvtNetworkBranchDefinition(
            :ferroresonance_suppression_resistance,
            :resistance,
            electromagnetic_primary_node,
            reference_node;
            resistance_ohm=definition.suppression_resistance_ohm,
        ),
    )
    definition.suppression_capacitance_f > 0.0 && push!(
        branches,
        CvtNetworkBranchDefinition(
            :ferroresonance_suppression_capacitance,
            :capacitance,
            electromagnetic_primary_node,
            reference_node;
            capacitance_f=definition.suppression_capacitance_f,
        ),
    )
    for (index, burden) in enumerate(measurement_burden_branches(
        definition.electromagnetic_unit.burden,
        secondary_node,
        reference_node,
    ))
        kind = burden.kind === :shunt_capacitance ? :capacitance : burden.kind
        push!(
            branches,
            CvtNetworkBranchDefinition(
                Symbol(:secondary_burden_, index),
                kind,
                burden.positive_node,
                burden.negative_node;
                resistance_ohm=burden.resistance_ohm,
                inductance_h=burden.inductance_h,
                capacitance_f=burden.capacitance_f,
            ),
        )
    end
    return branches
end

struct CouplingCapacitorVoltageTransformerReadiness
    ready::Bool
    reasons::Tuple
    divider_ratio::Float64
    series_resonance_hz::Float64
    electromagnetic_unit::InstrumentTransformerMeasurementReadiness
    deterministic_signature_sha256::String
end

function cvt_measurement_readiness(
    definition::CouplingCapacitorVoltageTransformerDefinition,
)
    electromagnetic_readiness = instrument_transformer_measurement_readiness(
        definition.electromagnetic_unit,
    )
    reasons = Symbol[]
    electromagnetic_readiness.ready || push!(reasons, :electromagnetic_unit_not_ready)
    resonance = cvt_series_resonance_hz(definition)
    resonance <= definition.maximum_spectral_frequency_hz ||
        push!(reasons, :cvt_resonance_outside_registered_band)
    return CouplingCapacitorVoltageTransformerReadiness(
        isempty(reasons),
        Tuple(reasons),
        cvt_divider_ratio(definition),
        resonance,
        electromagnetic_readiness,
        definition.deterministic_signature_sha256,
    )
end

function cvt_measurement_runtime(
    definition::CouplingCapacitorVoltageTransformerDefinition,
    preparation::TransformerApparatusPreparation,
    terminal_nodes,
    measurement_specification::MeasurementChainSpecification,
)
    readiness = cvt_measurement_readiness(definition)
    readiness.ready || throw(MeasurementChainRefusal(
        first(readiness.reasons),
        :prepare_cvt_measurement,
        definition.id,
        CouplingCapacitorVoltageTransformerMeasurement,
        "CVT definition is outside its declared electromagnetic or resonance domain",
        (reasons=readiness.reasons,),
    ))
    return instrument_transformer_measurement_runtime(
        definition.electromagnetic_unit,
        preparation,
        terminal_nodes,
        measurement_specification,
    )
end

function cvt_stored_energy_j(
    definition::CouplingCapacitorVoltageTransformerDefinition;
    line_voltage_v::Real,
    divider_voltage_v::Real,
    electromagnetic_primary_voltage_v::Real,
    compensation_current_a::Real,
)
    values = Float64.((
        line_voltage_v,
        divider_voltage_v,
        electromagnetic_primary_voltage_v,
        compensation_current_a,
    ))
    all(isfinite, values) || throw(ArgumentError(
        "CVT energy inputs must be finite",
    ))
    coupling_voltage = values[1] - values[2]
    return 0.5 * definition.high_voltage_capacitance_f * coupling_voltage^2 +
        0.5 * definition.intermediate_voltage_capacitance_f * values[2]^2 +
        0.5 * definition.compensation_inductance_h * values[4]^2 +
        0.5 * definition.suppression_capacitance_f * values[3]^2
end

struct CouplingCapacitorVoltageTransformerOutput
    measurement::InstrumentTransformerMeasurementOutput
    line_voltage_v::Float64
    divider_voltage_v::Float64
    electromagnetic_primary_voltage_v::Float64
    compensation_current_a::Float64
    stored_energy_j::Float64
    divider_ratio::Float64
    series_resonance_hz::Float64
    deterministic_signature_sha256::String
end

function cvt_measurement_output(
    definition::CouplingCapacitorVoltageTransformerDefinition,
    runtime::InstrumentTransformerMeasurementRuntime;
    line_voltage_v::Real,
    divider_voltage_v::Real,
    electromagnetic_primary_voltage_v::Real,
    compensation_current_a::Real,
)
    runtime.definition.deterministic_signature_sha256 ==
        definition.electromagnetic_unit.deterministic_signature_sha256 || throw(
        ArgumentError("CVT runtime does not belong to the requested definition"),
    )
    values = Float64.((
        line_voltage_v,
        divider_voltage_v,
        electromagnetic_primary_voltage_v,
        compensation_current_a,
    ))
    energy = cvt_stored_energy_j(
        definition;
        line_voltage_v=values[1],
        divider_voltage_v=values[2],
        electromagnetic_primary_voltage_v=values[3],
        compensation_current_a=values[4],
    )
    measurement = instrument_transformer_measurement_output(runtime)
    io = IOBuffer()
    println(io, definition.deterministic_signature_sha256)
    println(io, measurement.deterministic_signature_sha256)
    for value in (values..., energy)
        println(io, bitstring(value))
    end
    signature = bytes2hex(sha256(take!(io)))
    return CouplingCapacitorVoltageTransformerOutput(
        measurement,
        values...,
        energy,
        cvt_divider_ratio(definition),
        cvt_series_resonance_hz(definition),
        signature,
    )
end

export MeasurementProductFamily,
       LinearCurrentTransformerMeasurement,
       MagneticCurrentTransformerMeasurement,
       InductiveVoltageTransformerMeasurement,
       CouplingCapacitorVoltageTransformerMeasurement,
       ElectronicCurrentSensorMeasurement,
       ElectronicVoltageSensorMeasurement,
       ThreePhaseSampledMeasurement,
       MeasurementQuantizerTieRule,
       MeasurementTiesToEven,
       MeasurementTiesAwayFromZero,
       MeasurementTiesTowardZero,
       AbstractMeasurementQuantizer,
       UnquantizedMeasurement,
       UniformMeasurementQuantizer,
       AnalogMeasurementStateSpace,
       MeasurementAcquisitionSettings,
       MeasurementChainSpecification,
       MeasurementChainRefusal,
       measurement_chain_contract,
       measurement_chain_signature,
       measurement_sample_rate_hz

@enum MeasurementProductFamily begin
    LinearCurrentTransformerMeasurement
    MagneticCurrentTransformerMeasurement
    InductiveVoltageTransformerMeasurement
    CouplingCapacitorVoltageTransformerMeasurement
    ElectronicCurrentSensorMeasurement
    ElectronicVoltageSensorMeasurement
    ThreePhaseSampledMeasurement
end

@enum MeasurementQuantizerTieRule begin
    MeasurementTiesToEven
    MeasurementTiesAwayFromZero
    MeasurementTiesTowardZero
end

const _MEASUREMENT_PRODUCT_IDS = Dict(
    LinearCurrentTransformerMeasurement => :linear_current_transformer,
    MagneticCurrentTransformerMeasurement => :magnetic_current_transformer,
    InductiveVoltageTransformerMeasurement => :inductive_voltage_transformer,
    CouplingCapacitorVoltageTransformerMeasurement => :coupling_capacitor_voltage_transformer,
    ElectronicCurrentSensorMeasurement => :electronic_current_sensor,
    ElectronicVoltageSensorMeasurement => :electronic_voltage_sensor,
    ThreePhaseSampledMeasurement => :three_phase_sampled_measurement,
)

struct MeasurementChainRefusal <: Exception
    code::Symbol
    operation::Symbol
    measurement::Symbol
    family::MeasurementProductFamily
    message::String
    diagnostics::NamedTuple
end

Base.showerror(io::IO, refusal::MeasurementChainRefusal) = print(
    io,
    String(refusal.code),
    " during ",
    String(refusal.operation),
    " for measurement chain ",
    String(refusal.measurement),
    " (",
    String(_MEASUREMENT_PRODUCT_IDS[refusal.family]),
    "): ",
    refusal.message,
)

abstract type AbstractMeasurementQuantizer end

struct UnquantizedMeasurement <: AbstractMeasurementQuantizer
    lower_limit::Float64
    upper_limit::Float64

    function UnquantizedMeasurement(
        lower_limit::Real=-Inf,
        upper_limit::Real=Inf,
    )
        lower = Float64(lower_limit)
        upper = Float64(upper_limit)
        isnan(lower) && throw(ArgumentError("measurement lower limit must not be NaN"))
        isnan(upper) && throw(ArgumentError("measurement upper limit must not be NaN"))
        lower < upper || throw(ArgumentError(
            "measurement lower limit must be less than its upper limit",
        ))
        return new(lower, upper)
    end
end

struct UniformMeasurementQuantizer <: AbstractMeasurementQuantizer
    lower_limit::Float64
    upper_limit::Float64
    engineering_offset::Float64
    engineering_step::Float64
    minimum_code::Int64
    maximum_code::Int64
    tie_rule::MeasurementQuantizerTieRule

    function UniformMeasurementQuantizer(;
        lower_limit::Real,
        upper_limit::Real,
        engineering_offset::Real=0.0,
        engineering_step::Real,
        minimum_code::Integer,
        maximum_code::Integer,
        tie_rule::MeasurementQuantizerTieRule=MeasurementTiesToEven,
    )
        lower = Float64(lower_limit)
        upper = Float64(upper_limit)
        offset = Float64(engineering_offset)
        step = Float64(engineering_step)
        codes = (Int64(minimum_code), Int64(maximum_code))
        all(isfinite, (lower, upper, offset, step)) || throw(ArgumentError(
            "uniform measurement quantizer values must be finite",
        ))
        lower < upper || throw(ArgumentError(
            "uniform measurement quantizer limits must be strictly ordered",
        ))
        step > 0.0 || throw(ArgumentError(
            "uniform measurement quantizer step must be positive",
        ))
        codes[1] < codes[2] || throw(ArgumentError(
            "uniform measurement quantizer codes must be strictly ordered",
        ))
        representable_lower = offset + step * codes[1]
        representable_upper = offset + step * codes[2]
        tolerance = 32.0 * eps(Float64) * max(
            1.0,
            abs(lower),
            abs(upper),
            abs(representable_lower),
            abs(representable_upper),
        )
        lower >= representable_lower - tolerance || throw(ArgumentError(
            "measurement lower clipping limit is outside the quantizer code range",
        ))
        upper <= representable_upper + tolerance || throw(ArgumentError(
            "measurement upper clipping limit is outside the quantizer code range",
        ))
        return new(lower, upper, offset, step, codes..., tie_rule)
    end
end

"""One stable real continuous-time SISO model applied independently to every channel."""
struct AnalogMeasurementStateSpace
    state_matrix_per_s::Matrix{Float64}
    input_vector_per_s::Vector{Float64}
    output_vector::Vector{Float64}
    direct_gain::Float64
    stability_margin_per_s::Float64

    function AnalogMeasurementStateSpace(
        state_matrix_per_s::AbstractMatrix{<:Real},
        input_vector_per_s::AbstractVector{<:Real},
        output_vector::AbstractVector{<:Real},
        direct_gain::Real;
        stability_margin_per_s::Real=0.0,
    )
        state_matrix = Matrix{Float64}(state_matrix_per_s)
        input_vector = Float64.(input_vector_per_s)
        output = Float64.(output_vector)
        direct = Float64(direct_gain)
        margin = Float64(stability_margin_per_s)
        state_count = size(state_matrix, 1)
        size(state_matrix, 2) == state_count || throw(DimensionMismatch(
            "analog measurement state matrix must be square",
        ))
        length(input_vector) == state_count || throw(DimensionMismatch(
            "analog measurement input vector must match the state count",
        ))
        length(output) == state_count || throw(DimensionMismatch(
            "analog measurement output vector must match the state count",
        ))
        all(isfinite, state_matrix) && all(isfinite, input_vector) &&
            all(isfinite, output) && isfinite(direct) && isfinite(margin) ||
            throw(ArgumentError("analog measurement state-space data must be finite"))
        margin >= 0.0 || throw(ArgumentError(
            "analog measurement stability margin must be nonnegative",
        ))
        if state_count > 0
            largest_real_pole = maximum(real, eigvals(state_matrix))
            largest_real_pole < -margin || throw(ArgumentError(
                "analog measurement state space must be strictly stable at the declared margin",
            ))
        end
        return new(state_matrix, input_vector, output, direct, margin)
    end
end

AnalogMeasurementStateSpace(direct_gain::Real=1.0) = AnalogMeasurementStateSpace(
    zeros(0, 0),
    Float64[],
    Float64[],
    direct_gain,
)

struct MeasurementAcquisitionSettings
    tick_s::Float64
    sample_period_ticks::Int
    first_sample_tick::Int
    delay_ticks::Int
    window_weights_newest_first::Vector{Float64}
    nominal_frequency_hz::Float64
    positive_sequence_threshold::Float64
    frequency_update_separation::Int
    maximum_retained_samples::Int

    function MeasurementAcquisitionSettings(;
        tick_s::Real,
        sample_period_ticks::Integer,
        first_sample_tick::Integer=0,
        delay_ticks::Integer=0,
        window_weights_newest_first::AbstractVector{<:Real},
        nominal_frequency_hz::Real,
        positive_sequence_threshold::Real=0.0,
        frequency_update_separation::Integer=1,
        maximum_retained_samples::Integer=1,
    )
        tick = Float64(tick_s)
        period = Int(sample_period_ticks)
        first = Int(first_sample_tick)
        delay = Int(delay_ticks)
        weights = Float64.(window_weights_newest_first)
        nominal_frequency = Float64(nominal_frequency_hz)
        threshold = Float64(positive_sequence_threshold)
        separation = Int(frequency_update_separation)
        retained_samples = Int(maximum_retained_samples)
        isfinite(tick) && tick > 0.0 || throw(ArgumentError(
            "measurement scheduler tick must be finite and positive",
        ))
        period > 0 || throw(ArgumentError(
            "measurement sample period must be a positive tick count",
        ))
        first >= 0 || throw(ArgumentError(
            "measurement first sample tick must be nonnegative",
        ))
        delay >= 0 || throw(ArgumentError(
            "measurement delay must be a nonnegative tick count",
        ))
        4 <= length(weights) <= 4096 || throw(ArgumentError(
            "measurement estimator window must contain 4 through 4096 samples",
        ))
        all(isfinite, weights) || throw(ArgumentError(
            "measurement estimator weights must be finite",
        ))
        sum(weights) > 0.0 || throw(ArgumentError(
            "measurement estimator coherent gain must be positive",
        ))
        isfinite(nominal_frequency) && nominal_frequency > 0.0 ||
            throw(ArgumentError("measurement nominal frequency must be finite and positive"))
        isfinite(threshold) && threshold >= 0.0 || throw(ArgumentError(
            "measurement positive-sequence threshold must be finite and nonnegative",
        ))
        separation > 0 || throw(ArgumentError(
            "measurement frequency update separation must be positive",
        ))
        retained_samples > 0 || throw(ArgumentError(
            "measurement retained-sample limit must be positive",
        ))
        sample_rate = inv(tick * period)
        nominal_frequency < 0.4 * sample_rate || throw(ArgumentError(
            "measurement nominal frequency must remain below 40 percent of sample rate",
        ))
        return new(
            tick,
            period,
            first,
            delay,
            weights,
            nominal_frequency,
            threshold,
            separation,
            retained_samples,
        )
    end
end

measurement_sample_rate_hz(settings::MeasurementAcquisitionSettings) =
    inv(settings.tick_s * settings.sample_period_ticks)

struct MeasurementChainSpecification{Q<:AbstractMeasurementQuantizer}
    id::Symbol
    family::MeasurementProductFamily
    channel_names::Vector{Symbol}
    quantity::Symbol
    unit::String
    orientation::String
    phase_order::Vector{Symbol}
    conditioning::AnalogMeasurementStateSpace
    quantizer::Q
    acquisition::MeasurementAcquisitionSettings
    minimum_input::Float64
    maximum_input::Float64
    maximum_spectral_frequency_hz::Float64
    maximum_timestep_s::Float64
    provenance::ParameterProvenance

    function MeasurementChainSpecification(
        id::Symbol,
        family::MeasurementProductFamily;
        channel_names::AbstractVector{Symbol},
        quantity::Symbol,
        unit::AbstractString,
        orientation::AbstractString,
        phase_order::AbstractVector{Symbol}=Symbol[],
        conditioning::AnalogMeasurementStateSpace=AnalogMeasurementStateSpace(),
        quantizer::Q=UnquantizedMeasurement(),
        acquisition::MeasurementAcquisitionSettings,
        minimum_input::Real,
        maximum_input::Real,
        maximum_spectral_frequency_hz::Real,
        maximum_timestep_s::Real,
        provenance::ParameterProvenance,
    ) where {Q<:AbstractMeasurementQuantizer}
        channels = Symbol.(channel_names)
        phases = Symbol.(phase_order)
        unit_string = String(unit)
        orientation_string = String(orientation)
        input_limits = (Float64(minimum_input), Float64(maximum_input))
        maximum_frequency = Float64(maximum_spectral_frequency_hz)
        maximum_timestep = Float64(maximum_timestep_s)
        isempty(String(id)) && throw(ArgumentError("measurement-chain id must not be empty"))
        isempty(channels) && throw(ArgumentError(
            "measurement chain must declare at least one channel",
        ))
        length(channels) == length(unique(channels)) || throw(ArgumentError(
            "measurement channel names must be unique",
        ))
        isempty(String(quantity)) && throw(ArgumentError(
            "measurement quantity must not be empty",
        ))
        isempty(strip(unit_string)) && throw(ArgumentError(
            "measurement unit must not be empty",
        ))
        isempty(strip(orientation_string)) && throw(ArgumentError(
            "measurement orientation must not be empty",
        ))
        if isempty(phases)
            length(channels) == 3 && throw(ArgumentError(
                "a three-channel measurement chain must declare abc phase order",
            ))
        else
            phases == [:a, :b, :c] || throw(ArgumentError(
                "three-phase measurement chains currently require explicit abc phase order",
            ))
            length(channels) == 3 || throw(ArgumentError(
                "phase order is valid only for exactly three measurement channels",
            ))
        end
        if family === ThreePhaseSampledMeasurement
            length(channels) == 3 && phases == [:a, :b, :c] || throw(ArgumentError(
                "the complete three-phase sampled product requires three explicit abc channels",
            ))
        end
        all(isfinite, input_limits) && input_limits[1] < input_limits[2] ||
            throw(ArgumentError("measurement input limits must be finite and ordered"))
        isfinite(maximum_frequency) && maximum_frequency >= 0.0 ||
            throw(ArgumentError("measurement maximum spectral frequency must be finite and nonnegative"))
        maximum_frequency <= 0.4 * measurement_sample_rate_hz(acquisition) ||
            throw(ArgumentError("measurement spectral frequency exceeds the acquisition domain"))
        isfinite(maximum_timestep) && maximum_timestep > 0.0 || throw(ArgumentError(
            "measurement maximum timestep must be finite and positive",
        ))
        return new{Q}(
            id,
            family,
            channels,
            quantity,
            unit_string,
            orientation_string,
            phases,
            conditioning,
            quantizer,
            acquisition,
            input_limits...,
            maximum_frequency,
            maximum_timestep,
            provenance,
        )
    end
end

function measurement_chain_contract(specification::MeasurementChainSpecification)
    family_id = _MEASUREMENT_PRODUCT_IDS[specification.family]
    return ScientificModelContract(
        Symbol(family_id, :_measurement_chain),
        :instantaneous_emt_physical_to_sampled_measurement;
        owner="AIMORA.MeasurementChains",
        maturity=:implemented,
        fidelity=FieldCoupledDetailed,
        validity_domain=ModelValidityDomain(
            Symbol(family_id, :_fixed_step_measurement_domain);
            description="Explicit physical or electronic measurement input, stable analog conditioning, exact accepted-time acquisition, finite causal estimators, typed outputs, and no automatic family or fidelity fallback.",
            bounds=(
                NumericDomainBound(
                    :channel_count;
                    unit="count",
                    lower=1,
                    upper=8192,
                ),
                NumericDomainBound(
                    :fundamental_frequency_hz;
                    unit="Hz",
                    lower=45.0,
                    upper=65.0,
                ),
                NumericDomainBound(
                    :sample_rate_hz;
                    unit="sample/s",
                    lower=1.0e3,
                    upper=1.0e5,
                ),
                NumericDomainBound(
                    :timestep_s;
                    unit="s",
                    lower=0.0,
                    upper=specification.maximum_timestep_s,
                    lower_inclusive=false,
                ),
            ),
            unsupported_phenomena=(
                :protected_standard_conformance,
                :accuracy_class_or_calibration,
                :vendor_or_field_equivalence,
                :protection_logic_or_selectivity,
                :arbitrary_comtrade_extension,
                :atp_or_pscad_equivalence,
                :hil_or_realtime_qualification,
                :certification,
            ),
        ),
        state_inventory=DynamicStateInventory(
            differential=(:analog_conditioning_state,),
            algebraic=(:oriented_input, :conditioned_output),
            discrete=(
                :clip_state,
                :quantizer_code,
                :window_validity,
                :channel_quality,
            ),
            delayed_history=(
                :delayed_sample_queue,
                :held_sample,
                :rms_window,
                :phasor_window,
                :frequency_phase_history,
            ),
            scheduler=(
                :next_sample_tick,
                :sample_count,
                :release_count,
            ),
            random=(),
        ),
        inputs=(
            ContractQuantity(
                :accepted_oriented_input;
                unit=specification.unit,
                orientation=specification.orientation,
            ),
        ),
        outputs=(
            ContractQuantity(
                :conditioned_value;
                unit=specification.unit,
                orientation=specification.orientation,
            ),
            ContractQuantity(
                :sliding_rms;
                unit=specification.unit,
                orientation="nonnegative causal full window",
            ),
            ContractQuantity(
                :fundamental_rms_phasor;
                unit=specification.unit,
                orientation="exp(j*omega*t) at declared reference time",
            ),
            ContractQuantity(
                :frequency;
                unit="Hz",
                orientation="positive-sequence phase increment",
            ),
        ),
        assumptions=(
            "The caller supplies an accepted-time physical or electronic input with explicit units, orientation, source, uncertainty, and validity.",
            "The analog state-space model is real, causal, and strictly stable; acquisition clocks and delays are exact integer ticks.",
            "Three-phase sequence and frequency outputs require explicit abc channels from one aligned quantity, window, and time base.",
            "The generic estimators make no protected-standard, relay, vendor, field, ATP/PSCAD, HIL, or certification claim.",
        ),
        mutation_order=(
            :validate_identity_and_domain,
            :prepare_pure_analog_candidate,
            :accept_physical_state_once,
            :sample_exact_due_tick,
            :release_exact_delayed_sample,
            :update_causal_estimators,
            :publish_typed_output_once,
        ),
    )
end

function _write_measurement_signature_value(io::IO, value)
    if value === nothing
        write(io, UInt8('n'))
    elseif value isa Bool
        write(io, UInt8('b'), value ? UInt8(1) : UInt8(0))
    elseif value isa Float64
        write(io, UInt8('f'))
        write(io, hton(reinterpret(UInt64, value)))
    elseif value isa ComplexF64
        write(io, UInt8('c'))
        write(io, hton(reinterpret(UInt64, real(value))))
        write(io, hton(reinterpret(UInt64, imag(value))))
    elseif value isa Int
        write(io, UInt8('i'))
        write(io, hton(reinterpret(UInt64, Int64(value))))
    elseif value isa Symbol
        bytes = codeunits(String(value))
        write(io, UInt8('y'), hton(UInt64(length(bytes))), bytes)
    elseif value isa AbstractString
        bytes = codeunits(value)
        write(io, UInt8('s'), hton(UInt64(length(bytes))), bytes)
    elseif value isa AbstractArray
        print(io, '[', size(value), ':')
        for item in value
            _write_measurement_signature_value(io, item)
            print(io, ',')
        end
        print(io, ']')
    elseif value isa ParameterProvenance
        for field in fieldnames(ParameterProvenance)
            _write_measurement_signature_value(io, getfield(value, field))
            print(io, '|')
        end
    elseif value isa NamedTuple
        print(io, '{')
        for field in keys(value)
            show(io, field)
            print(io, '=')
            _write_measurement_signature_value(io, getfield(value, field))
            print(io, ',')
        end
        print(io, '}')
    elseif value isa Tuple
        print(io, '(')
        for item in value
            _write_measurement_signature_value(io, item)
            print(io, ',')
        end
        print(io, ')')
    else
        show(io, value)
    end
    return io
end

function measurement_chain_signature(specification::MeasurementChainSpecification)
    io = IOBuffer()
    for field in fieldnames(typeof(specification))
        _write_measurement_signature_value(io, getfield(specification, field))
        print(io, '\n')
    end
    return bytes2hex(sha256(take!(io)))
end

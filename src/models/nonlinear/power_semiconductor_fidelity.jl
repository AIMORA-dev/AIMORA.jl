abstract type AbstractPowerSemiconductorFidelityComponent end

function _generic_semiconductor_provenance(
    units::AbstractString,
    validity_domain::AbstractString,
)
    return ParameterProvenance(
        "AIMORA generic caller-supplied semiconductor parameters",
        String(units),
        "converted explicitly to Float64 SI quantities without inferred device scaling",
        "unknown unless narrowed by the caller's parameter record",
        String(validity_domain),
        PhysicalModelParameter,
    )
end

function _validate_semiconductor_provenance(
    provenance::NonlinearParameterProvenance,
    owner::AbstractString,
)
    provenance.nature === PhysicalModelParameter || throw(ArgumentError(
        "$owner provenance must describe physical model parameters",
    ))
    isempty(strip(provenance.source)) && throw(ArgumentError(
        "$owner provenance source must not be empty",
    ))
    isempty(strip(provenance.units)) && throw(ArgumentError(
        "$owner provenance units must not be empty",
    ))
    isempty(strip(provenance.transformation)) && throw(ArgumentError(
        "$owner provenance transformation must not be empty",
    ))
    isempty(strip(provenance.uncertainty)) && throw(ArgumentError(
        "$owner provenance uncertainty must not be empty",
    ))
    isempty(strip(provenance.validity_domain)) && throw(ArgumentError(
        "$owner provenance validity domain must not be empty",
    ))
    return provenance
end

"""First-order stored-charge state for a diode or an explicitly selected antiparallel path."""
mutable struct RecoveredChargeFidelity <: AbstractPowerSemiconductorFidelityComponent
    lifetime_s::Float64
    stored_charge_c::Float64
    previous_stored_charge_c::Float64
    recovery_active::Bool
    last_recovery_current_a::Float64
    peak_reverse_current_a::Float64
    cumulative_recovered_charge_c::Float64
    recovery_start_time_s::Float64
    last_recovery_duration_s::Float64
    recovery_zero_event_count::Int
    last_recovery_zero_time_s::Float64
    provenance::NonlinearParameterProvenance
end

function RecoveredChargeFidelity(
    lifetime_s::Real;
    initial_charge_c::Real=0.0,
    provenance::NonlinearParameterProvenance=_generic_semiconductor_provenance(
        "second and coulomb",
        "first-order recovered charge with positive lifetime through 10 ms and charge from zero through 10 C",
    ),
)
    lifetime = Float64(lifetime_s)
    initial_charge = Float64(initial_charge_c)
    isfinite(lifetime) && 0.0 < lifetime <= 1.0e-2 || throw(ArgumentError(
        "recovered-charge lifetime must be finite, positive, and no greater than 10 ms",
    ))
    isfinite(initial_charge) && 0.0 <= initial_charge <= 10.0 || throw(ArgumentError(
        "initial recovered charge must be finite and between zero and 10 C",
    ))
    _validate_semiconductor_provenance(provenance, "recovered-charge")
    return RecoveredChargeFidelity(
        lifetime,
        initial_charge,
        initial_charge,
        false,
        0.0,
        0.0,
        0.0,
        Inf,
        0.0,
        0,
        -Inf,
        provenance,
    )
end

"""Positive continuous voltage-dependent junction charge with an accepted voltage history."""
mutable struct NonlinearJunctionChargeFidelity <: AbstractPowerSemiconductorFidelityComponent
    zero_bias_capacitance_f::Float64
    junction_voltage_v::Float64
    grading_exponent::Float64
    minimum_voltage_v::Float64
    maximum_voltage_v::Float64
    previous_voltage_v::Float64
    previous_charge_c::Float64
    last_capacitance_f::Float64
    last_charge_c::Float64
    last_displacement_current_a::Float64
    provenance::NonlinearParameterProvenance
end

function power_semiconductor_junction_capacitance(
    zero_bias_capacitance_f::Real,
    junction_voltage_v::Real,
    grading_exponent::Real,
    voltage_v::Real,
)
    capacitance = Float64(zero_bias_capacitance_f)
    junction_voltage = Float64(junction_voltage_v)
    exponent = Float64(grading_exponent)
    voltage = Float64(voltage_v)
    isfinite(capacitance) && capacitance > 0.0 || throw(ArgumentError(
        "zero-bias junction capacitance must be finite and positive",
    ))
    isfinite(junction_voltage) && junction_voltage > 0.0 || throw(ArgumentError(
        "junction voltage must be finite and positive",
    ))
    isfinite(exponent) && 0.0 <= exponent < 1.0 || throw(ArgumentError(
        "junction grading exponent must be finite and in [0, 1)",
    ))
    isfinite(voltage) || throw(ArgumentError("junction voltage evaluation must be finite"))
    return voltage >= 0.0 ? capacitance :
        capacitance * (1.0 - voltage / junction_voltage)^(-exponent)
end

function power_semiconductor_junction_charge(
    zero_bias_capacitance_f::Real,
    junction_voltage_v::Real,
    grading_exponent::Real,
    voltage_v::Real,
)
    capacitance = Float64(zero_bias_capacitance_f)
    junction_voltage = Float64(junction_voltage_v)
    exponent = Float64(grading_exponent)
    voltage = Float64(voltage_v)
    power_semiconductor_junction_capacitance(
        capacitance,
        junction_voltage,
        exponent,
        voltage,
    )
    voltage >= 0.0 && return capacitance * voltage
    normalized = 1.0 - voltage / junction_voltage
    return -capacitance * junction_voltage *
        (normalized^(1.0 - exponent) - 1.0) / (1.0 - exponent)
end

function _power_semiconductor_junction_charge_capacitance(
    fidelity::NonlinearJunctionChargeFidelity,
    voltage_v::Float64,
)
    fidelity.minimum_voltage_v <= voltage_v <= fidelity.maximum_voltage_v ||
        throw(DomainError(
            voltage_v,
            "junction voltage lies outside its declared domain",
        ))
    voltage_v == fidelity.previous_voltage_v && return (
        fidelity.previous_charge_c,
        fidelity.last_capacitance_f,
    )
    capacitance = fidelity.zero_bias_capacitance_f
    voltage_v >= 0.0 && return capacitance * voltage_v, capacitance
    normalized = 1.0 - voltage_v / fidelity.junction_voltage_v
    normalized_negative_power = normalized^(-fidelity.grading_exponent)
    charge = -capacitance * fidelity.junction_voltage_v * (
        normalized * normalized_negative_power - 1.0
    ) / (1.0 - fidelity.grading_exponent)
    return charge, capacitance * normalized_negative_power
end

function power_semiconductor_junction_stored_energy(
    zero_bias_capacitance_f::Real,
    junction_voltage_v::Real,
    grading_exponent::Real,
    voltage_v::Real,
)
    capacitance = Float64(zero_bias_capacitance_f)
    junction_voltage = Float64(junction_voltage_v)
    exponent = Float64(grading_exponent)
    voltage = Float64(voltage_v)
    power_semiconductor_junction_capacitance(
        capacitance,
        junction_voltage,
        exponent,
        voltage,
    )
    voltage >= 0.0 && return 0.5 * capacitance * voltage^2
    normalized = 1.0 - voltage / junction_voltage
    return capacitance * junction_voltage^2 * (
        (normalized^(2.0 - exponent) - 1.0) / (2.0 - exponent) -
        (normalized^(1.0 - exponent) - 1.0) / (1.0 - exponent)
    )
end

function NonlinearJunctionChargeFidelity(
    zero_bias_capacitance_f::Real,
    junction_voltage_v::Real,
    grading_exponent::Real;
    voltage_domain_v::Tuple{<:Real,<:Real}=(-20_000.0, 20_000.0),
    initial_voltage_v::Real=0.0,
    provenance::NonlinearParameterProvenance=_generic_semiconductor_provenance(
        "farad, volt, and dimensionless grading exponent",
        "continuous positive junction charge over one explicitly bounded terminal-voltage domain",
    ),
)
    capacitance = Float64(zero_bias_capacitance_f)
    junction_voltage = Float64(junction_voltage_v)
    exponent = Float64(grading_exponent)
    minimum_voltage, maximum_voltage = Float64.(voltage_domain_v)
    initial_voltage = Float64(initial_voltage_v)
    1.0e-12 <= capacitance <= 1.0e-2 || throw(ArgumentError(
        "zero-bias junction capacitance must be between 1 pF and 10 mF",
    ))
    isfinite(junction_voltage) && junction_voltage > 0.0 || throw(ArgumentError(
        "junction voltage must be finite and positive",
    ))
    isfinite(exponent) && 0.0 <= exponent < 1.0 || throw(ArgumentError(
        "junction grading exponent must be finite and in [0, 1)",
    ))
    isfinite(minimum_voltage) && isfinite(maximum_voltage) &&
        minimum_voltage < maximum_voltage || throw(ArgumentError(
        "junction voltage domain must be finite and increasing",
    ))
    minimum_voltage <= initial_voltage <= maximum_voltage || throw(ArgumentError(
        "initial junction voltage lies outside its declared domain",
    ))
    _validate_semiconductor_provenance(provenance, "nonlinear-junction-charge")
    initial_charge = power_semiconductor_junction_charge(
        capacitance,
        junction_voltage,
        exponent,
        initial_voltage,
    )
    initial_capacitance = power_semiconductor_junction_capacitance(
        capacitance,
        junction_voltage,
        exponent,
        initial_voltage,
    )
    return NonlinearJunctionChargeFidelity(
        capacitance,
        junction_voltage,
        exponent,
        minimum_voltage,
        maximum_voltage,
        initial_voltage,
        initial_charge,
        initial_capacitance,
        initial_charge,
        0.0,
        provenance,
    )
end

function power_semiconductor_junction_capacitance(
    fidelity::NonlinearJunctionChargeFidelity,
    voltage_v::Real,
)
    voltage = Float64(voltage_v)
    fidelity.minimum_voltage_v <= voltage <= fidelity.maximum_voltage_v || throw(DomainError(
        voltage,
        "junction voltage lies outside its declared domain",
    ))
    return power_semiconductor_junction_capacitance(
        fidelity.zero_bias_capacitance_f,
        fidelity.junction_voltage_v,
        fidelity.grading_exponent,
        voltage,
    )
end

function power_semiconductor_junction_charge(
    fidelity::NonlinearJunctionChargeFidelity,
    voltage_v::Real,
)
    voltage = Float64(voltage_v)
    fidelity.minimum_voltage_v <= voltage <= fidelity.maximum_voltage_v || throw(DomainError(
        voltage,
        "junction voltage lies outside its declared domain",
    ))
    return power_semiconductor_junction_charge(
        fidelity.zero_bias_capacitance_f,
        fidelity.junction_voltage_v,
        fidelity.grading_exponent,
        voltage,
    )
end

function power_semiconductor_junction_stored_energy(
    fidelity::NonlinearJunctionChargeFidelity,
    voltage_v::Real,
)
    voltage = Float64(voltage_v)
    fidelity.minimum_voltage_v <= voltage <= fidelity.maximum_voltage_v || throw(DomainError(
        voltage,
        "junction voltage lies outside its declared domain",
    ))
    return power_semiconductor_junction_stored_energy(
        fidelity.zero_bias_capacitance_f,
        fidelity.junction_voltage_v,
        fidelity.grading_exponent,
        voltage,
    )
end

"""Generic exponential forward-current tail active only after an accepted turn-off."""
mutable struct TurnOffTailFidelity <: AbstractPowerSemiconductorFidelityComponent
    decay_time_s::Float64
    cutoff_current_a::Float64
    active::Bool
    current_a::Float64
    initial_current_a::Float64
    turn_off_time_s::Float64
    last_duration_s::Float64
    cutoff_event_count::Int
    last_cutoff_time_s::Float64
    provenance::NonlinearParameterProvenance
end

function TurnOffTailFidelity(
    decay_time_s::Real;
    cutoff_current_a::Real=0.0,
    provenance::NonlinearParameterProvenance=_generic_semiconductor_provenance(
        "second and ampere",
        "generic nonnegative exponential forward-current tail through 10 ms",
    ),
)
    decay_time = Float64(decay_time_s)
    cutoff = Float64(cutoff_current_a)
    isfinite(decay_time) && 0.0 < decay_time <= 1.0e-2 || throw(ArgumentError(
        "turn-off-tail time constant must be finite, positive, and no greater than 10 ms",
    ))
    isfinite(cutoff) && cutoff >= 0.0 || throw(ArgumentError(
        "turn-off-tail cutoff current must be finite and nonnegative",
    ))
    _validate_semiconductor_provenance(provenance, "turn-off-tail")
    return TurnOffTailFidelity(
        decay_time,
        cutoff,
        false,
        0.0,
        0.0,
        Inf,
        0.0,
        0,
        -Inf,
        provenance,
    )
end

function _validate_energy_axis(axis::AbstractVector{<:Real}, name::AbstractString)
    values = Float64.(axis)
    2 <= length(values) <= 64 || throw(ArgumentError(
        "$name switching-energy axis must contain between two and 64 points",
    ))
    all(isfinite, values) || throw(ArgumentError("$name switching-energy axis must be finite"))
    all(diff(values) .> 0.0) || throw(ArgumentError(
        "$name switching-energy axis must be strictly increasing",
    ))
    return values
end

function _validate_energy_values(
    values::AbstractArray{<:Real,3},
    dimensions::NTuple{3,Int},
    name::AbstractString,
)
    size(values) == dimensions || throw(DimensionMismatch(
        "$name switching-energy table dimensions must match current, voltage, and temperature axes",
    ))
    energies = Float64.(values)
    all(value -> isfinite(value) && value >= 0.0, energies) || throw(ArgumentError(
        "$name switching-energy values must be finite and nonnegative",
    ))
    return energies
end

"""Deterministic trilinear turn-on, turn-off, and reverse-recovery event energy."""
mutable struct SwitchingEnergyTable <: AbstractPowerSemiconductorFidelityComponent
    current_axis_a::Vector{Float64}
    blocking_voltage_axis_v::Vector{Float64}
    junction_temperature_axis_k::Vector{Float64}
    turn_on_energy_j::Array{Float64,3}
    turn_off_energy_j::Array{Float64,3}
    reverse_recovery_energy_j::Array{Float64,3}
    cumulative_turn_on_energy_j::Float64
    cumulative_turn_off_energy_j::Float64
    cumulative_reverse_recovery_energy_j::Float64
    last_event_kind::Symbol
    last_event_energy_j::Float64
    last_event_transition_count::Int
    last_reverse_recovery_start_time_s::Float64
    provenance::NonlinearParameterProvenance
end

function SwitchingEnergyTable(
    current_axis_a::AbstractVector{<:Real},
    blocking_voltage_axis_v::AbstractVector{<:Real},
    junction_temperature_axis_k::AbstractVector{<:Real};
    turn_on_energy_j::AbstractArray{<:Real,3},
    turn_off_energy_j::AbstractArray{<:Real,3},
    reverse_recovery_energy_j::AbstractArray{<:Real,3},
    provenance::NonlinearParameterProvenance=_generic_semiconductor_provenance(
        "ampere, volt, kelvin, and joule",
        "deterministic nonnegative trilinear event energy with no extrapolation",
    ),
)
    current = _validate_energy_axis(current_axis_a, "current")
    voltage = _validate_energy_axis(blocking_voltage_axis_v, "blocking-voltage")
    temperature = _validate_energy_axis(junction_temperature_axis_k, "temperature")
    first(current) >= 0.0 || throw(ArgumentError(
        "switching-energy current axis must be nonnegative",
    ))
    first(voltage) >= 0.0 || throw(ArgumentError(
        "switching-energy blocking-voltage axis must be nonnegative",
    ))
    first(temperature) >= 200.0 && last(temperature) <= 600.0 || throw(ArgumentError(
        "switching-energy temperature axis must lie within 200 K through 600 K",
    ))
    dimensions = (length(current), length(voltage), length(temperature))
    on_energy = _validate_energy_values(turn_on_energy_j, dimensions, "turn-on")
    off_energy = _validate_energy_values(turn_off_energy_j, dimensions, "turn-off")
    recovery_energy = _validate_energy_values(
        reverse_recovery_energy_j,
        dimensions,
        "reverse-recovery",
    )
    _validate_semiconductor_provenance(provenance, "switching-energy")
    return SwitchingEnergyTable(
        current,
        voltage,
        temperature,
        on_energy,
        off_energy,
        recovery_energy,
        0.0,
        0.0,
        0.0,
        :none,
        0.0,
        0,
        Inf,
        provenance,
    )
end

function _energy_axis_cell(axis::Vector{Float64}, coordinate::Float64, name::AbstractString)
    first(axis) <= coordinate <= last(axis) || throw(DomainError(
        coordinate,
        "$name lies outside the switching-energy table domain; extrapolation is unavailable",
    ))
    coordinate == last(axis) && return length(axis) - 1, 1.0
    lower = searchsortedlast(axis, coordinate)
    lower = clamp(lower, 1, length(axis) - 1)
    fraction = (coordinate - axis[lower]) / (axis[lower + 1] - axis[lower])
    return lower, fraction
end

function _trilinear_energy(
    table::Array{Float64,3},
    current_axis::Vector{Float64},
    voltage_axis::Vector{Float64},
    temperature_axis::Vector{Float64},
    current_a::Float64,
    voltage_v::Float64,
    temperature_k::Float64,
)
    current_index, current_fraction = _energy_axis_cell(current_axis, current_a, "current")
    voltage_index, voltage_fraction = _energy_axis_cell(voltage_axis, voltage_v, "blocking voltage")
    temperature_index, temperature_fraction = _energy_axis_cell(
        temperature_axis,
        temperature_k,
        "junction temperature",
    )
    energy = 0.0
    for current_offset in 0:1, voltage_offset in 0:1, temperature_offset in 0:1
        weight = (current_offset == 0 ? 1.0 - current_fraction : current_fraction) *
            (voltage_offset == 0 ? 1.0 - voltage_fraction : voltage_fraction) *
            (temperature_offset == 0 ? 1.0 - temperature_fraction : temperature_fraction)
        energy += weight * table[
            current_index + current_offset,
            voltage_index + voltage_offset,
            temperature_index + temperature_offset,
        ]
    end
    return energy
end

function power_semiconductor_switching_energy(
    table::SwitchingEnergyTable,
    event_kind::Symbol,
    current_a::Real,
    blocking_voltage_v::Real,
    junction_temperature_k::Real,
)
    current = Float64(current_a)
    voltage = Float64(blocking_voltage_v)
    temperature = Float64(junction_temperature_k)
    all(isfinite, (current, voltage, temperature)) || throw(ArgumentError(
        "switching-energy coordinates must be finite",
    ))
    current >= 0.0 && voltage >= 0.0 || throw(ArgumentError(
        "switching-energy current and blocking-voltage magnitudes must be nonnegative",
    ))
    values = if event_kind === :turn_on
        table.turn_on_energy_j
    elseif event_kind === :turn_off
        table.turn_off_energy_j
    elseif event_kind === :reverse_recovery
        table.reverse_recovery_energy_j
    else
        throw(ArgumentError("unknown switching-energy event kind $event_kind"))
    end
    return _trilinear_energy(
        values,
        table.current_axis_a,
        table.blocking_voltage_axis_v,
        table.junction_temperature_axis_k,
        current,
        voltage,
        temperature,
    )
end

"""Passive one-through-eight-stage Cauer thermal state with ambient rejection."""
mutable struct CauerThermalFidelity <: AbstractPowerSemiconductorFidelityComponent
    capacitance_j_per_k::Vector{Float64}
    resistance_k_per_w::Vector{Float64}
    ambient_temperature_k::Float64
    node_temperature_k::Vector{Float64}
    minimum_temperature_k::Float64
    maximum_temperature_k::Float64
    last_loss_power_w::Float64
    last_ambient_heat_flow_w::Float64
    last_stored_energy_j::Float64
    cumulative_input_energy_j::Float64
    cumulative_ambient_energy_j::Float64
    provenance::NonlinearParameterProvenance
    trial_lower_conductance_w_per_k::Vector{Float64}
    trial_diagonal_conductance_w_per_k::Vector{Float64}
    trial_upper_conductance_w_per_k::Vector{Float64}
    trial_right_hand_side_w::Vector{Float64}
    trial_temperature_rise_k::Vector{Float64}
    trial_temperature_k::Vector{Float64}
    trial_storage_conductance_w_per_k::Vector{Float64}
    factorized_step_s::Float64
end

function CauerThermalFidelity(
    capacitance_j_per_k::AbstractVector{<:Real},
    resistance_k_per_w::AbstractVector{<:Real};
    ambient_temperature_k::Real=293.15,
    initial_temperature_k::AbstractVector{<:Real}=fill(
        Float64(ambient_temperature_k),
        length(capacitance_j_per_k),
    ),
    temperature_domain_k::Tuple{<:Real,<:Real}=(200.0, 600.0),
    provenance::NonlinearParameterProvenance=_generic_semiconductor_provenance(
        "joule per kelvin, kelvin per watt, kelvin, watt, and joule",
        "one-through-eight-stage passive lumped Cauer network over 200 K through 600 K",
    ),
)
    capacitance = Float64.(capacitance_j_per_k)
    resistance = Float64.(resistance_k_per_w)
    1 <= length(capacitance) <= 8 || throw(ArgumentError(
        "Cauer thermal network must contain between one and eight stages",
    ))
    length(resistance) == length(capacitance) || throw(DimensionMismatch(
        "Cauer network requires one resistance per thermal stage, including ambient rejection",
    ))
    all(value -> isfinite(value) && value > 0.0, capacitance) || throw(ArgumentError(
        "Cauer thermal capacitances must be finite and positive",
    ))
    all(value -> isfinite(value) && value > 0.0, resistance) || throw(ArgumentError(
        "Cauer thermal resistances must be finite and positive",
    ))
    ambient = Float64(ambient_temperature_k)
    temperatures = Float64.(initial_temperature_k)
    length(temperatures) == length(capacitance) || throw(DimensionMismatch(
        "initial thermal temperatures must cover every Cauer stage",
    ))
    minimum_temperature, maximum_temperature = Float64.(temperature_domain_k)
    200.0 <= minimum_temperature < maximum_temperature <= 600.0 || throw(ArgumentError(
        "thermal temperature domain must be increasing within 200 K through 600 K",
    ))
    minimum_temperature <= ambient <= maximum_temperature || throw(ArgumentError(
        "ambient temperature lies outside the declared thermal domain",
    ))
    all(temperature -> isfinite(temperature) &&
        minimum_temperature <= temperature <= maximum_temperature, temperatures) ||
        throw(ArgumentError("initial thermal-node temperature lies outside its domain"))
    _validate_semiconductor_provenance(provenance, "Cauer-thermal")
    rises = temperatures .- ambient
    stored_energy = sum(capacitance .* rises)
    ambient_flow = rises[end] / resistance[end]
    return CauerThermalFidelity(
        capacitance,
        resistance,
        ambient,
        temperatures,
        minimum_temperature,
        maximum_temperature,
        0.0,
        ambient_flow,
        stored_energy,
        0.0,
        0.0,
        provenance,
        zeros(Float64, max(0, length(capacitance) - 1)),
        zeros(Float64, length(capacitance)),
        zeros(Float64, max(0, length(capacitance) - 1)),
        zeros(Float64, length(capacitance)),
        zeros(Float64, length(capacitance)),
        copy(temperatures),
        zeros(Float64, length(capacitance)),
        NaN,
    )
end

function power_semiconductor_thermal_stored_energy(thermal::CauerThermalFidelity)
    stored_energy_j = 0.0
    for stage in eachindex(
        thermal.capacitance_j_per_k,
        thermal.node_temperature_k,
    )
        stored_energy_j += thermal.capacitance_j_per_k[stage] *
            (thermal.node_temperature_k[stage] - thermal.ambient_temperature_k)
    end
    return stored_energy_j
end

function _prepare_cauer_thermal_factorization!(
    thermal::CauerThermalFidelity,
    step_s::Float64,
)
    step_s == thermal.factorized_step_s && return thermal
    stage_count = length(thermal.capacitance_j_per_k)
    lower = thermal.trial_lower_conductance_w_per_k
    diagonal = thermal.trial_diagonal_conductance_w_per_k
    upper = thermal.trial_upper_conductance_w_per_k
    storage = thermal.trial_storage_conductance_w_per_k
    for stage in 1:stage_count
        storage[stage] = thermal.capacitance_j_per_k[stage] / step_s
        diagonal[stage] = storage[stage]
        if stage > 1
            coupling = inv(thermal.resistance_k_per_w[stage - 1])
            diagonal[stage] += coupling
            lower[stage - 1] = -coupling
        end
        if stage < stage_count
            coupling = inv(thermal.resistance_k_per_w[stage])
            diagonal[stage] += coupling
            upper[stage] = -coupling
        else
            diagonal[stage] += inv(thermal.resistance_k_per_w[stage])
        end
    end
    for stage in 2:stage_count
        elimination_factor = lower[stage - 1] / diagonal[stage - 1]
        lower[stage - 1] = elimination_factor
        diagonal[stage] -= elimination_factor * upper[stage - 1]
    end
    thermal.factorized_step_s = step_s
    return thermal
end

function _cauer_thermal_trial(
    thermal::CauerThermalFidelity,
    loss_power_w::Float64,
    step_s::Float64,
)
    isfinite(loss_power_w) && loss_power_w >= 0.0 || throw(ArgumentError(
        "accepted semiconductor thermal input power must be finite and nonnegative",
    ))
    isfinite(step_s) && step_s > 0.0 || throw(ArgumentError(
        "Cauer thermal timestep must be finite and positive",
    ))
    stage_count = length(thermal.capacitance_j_per_k)
    lower = thermal.trial_lower_conductance_w_per_k
    diagonal = thermal.trial_diagonal_conductance_w_per_k
    upper = thermal.trial_upper_conductance_w_per_k
    right_hand_side = thermal.trial_right_hand_side_w
    rise = thermal.trial_temperature_rise_k
    temperature = thermal.trial_temperature_k
    storage = thermal.trial_storage_conductance_w_per_k
    _prepare_cauer_thermal_factorization!(thermal, step_s)
    if stage_count == 4
        ambient = thermal.ambient_temperature_k
        @inbounds begin
            right_hand_side[1] = storage[1] *
                (thermal.node_temperature_k[1] - ambient) + loss_power_w
            right_hand_side[2] = storage[2] *
                (thermal.node_temperature_k[2] - ambient) -
                lower[1] * right_hand_side[1]
            right_hand_side[3] = storage[3] *
                (thermal.node_temperature_k[3] - ambient) -
                lower[2] * right_hand_side[2]
            right_hand_side[4] = storage[4] *
                (thermal.node_temperature_k[4] - ambient) -
                lower[3] * right_hand_side[3]
            rise[4] = right_hand_side[4] / diagonal[4]
            rise[3] = (right_hand_side[3] - upper[3] * rise[4]) / diagonal[3]
            rise[2] = (right_hand_side[2] - upper[2] * rise[3]) / diagonal[2]
            rise[1] = (right_hand_side[1] - upper[1] * rise[2]) / diagonal[1]
            temperature[1] = rise[1] + ambient
            temperature[2] = rise[2] + ambient
            temperature[3] = rise[3] + ambient
            temperature[4] = rise[4] + ambient
        end
        minimum_temperature = thermal.minimum_temperature_k
        maximum_temperature = thermal.maximum_temperature_k
        @inbounds valid_temperature =
            isfinite(temperature[1]) && minimum_temperature <= temperature[1] <= maximum_temperature &&
            isfinite(temperature[2]) && minimum_temperature <= temperature[2] <= maximum_temperature &&
            isfinite(temperature[3]) && minimum_temperature <= temperature[3] <= maximum_temperature &&
            isfinite(temperature[4]) && minimum_temperature <= temperature[4] <= maximum_temperature
        valid_temperature || throw(DomainError(
            copy(temperature),
            "Cauer thermal state left its declared domain",
        ))
        @inbounds stored_energy =
            thermal.capacitance_j_per_k[1] * rise[1] +
            thermal.capacitance_j_per_k[2] * rise[2] +
            thermal.capacitance_j_per_k[3] * rise[3] +
            thermal.capacitance_j_per_k[4] * rise[4]
        @inbounds ambient_heat_flow = rise[4] / thermal.resistance_k_per_w[4]
        return ambient_heat_flow, stored_energy
    end
    @inbounds for stage in 1:stage_count
        right_hand_side[stage] = storage[stage] *
            (thermal.node_temperature_k[stage] - thermal.ambient_temperature_k)
        if stage == 1
            right_hand_side[stage] += loss_power_w
        end
    end
    @inbounds for stage in 2:stage_count
        right_hand_side[stage] -= lower[stage - 1] * right_hand_side[stage - 1]
    end
    rise[stage_count] = right_hand_side[stage_count] / diagonal[stage_count]
    @inbounds for stage in (stage_count - 1):-1:1
        rise[stage] = (
            right_hand_side[stage] - upper[stage] * rise[stage + 1]
        ) / diagonal[stage]
    end
    stored_energy = 0.0
    @inbounds for stage in 1:stage_count
        stage_temperature = rise[stage] + thermal.ambient_temperature_k
        isfinite(stage_temperature) &&
            thermal.minimum_temperature_k <= stage_temperature <=
                thermal.maximum_temperature_k || throw(DomainError(
            stage_temperature,
            "Cauer thermal state left its declared domain",
        ))
        temperature[stage] = stage_temperature
        stored_energy += thermal.capacitance_j_per_k[stage] * rise[stage]
    end
    ambient_heat_flow = rise[end] / thermal.resistance_k_per_w[end]
    return ambient_heat_flow, stored_energy
end

function _accept_cauer_thermal_step!(
    thermal::CauerThermalFidelity,
    loss_power_w::Float64,
    step_s::Float64,
)
    ambient_heat_flow, stored_energy = _cauer_thermal_trial(
        thermal,
        loss_power_w,
        step_s,
    )
    copyto!(thermal.node_temperature_k, thermal.trial_temperature_k)
    thermal.last_loss_power_w = loss_power_w
    thermal.last_ambient_heat_flow_w = ambient_heat_flow
    thermal.last_stored_energy_j = stored_energy
    thermal.cumulative_input_energy_j += step_s * loss_power_w
    thermal.cumulative_ambient_energy_j += step_s * ambient_heat_flow
    return thermal
end

function accept_power_semiconductor_thermal_step!(
    thermal::CauerThermalFidelity,
    loss_power_w::Real,
    step_s::Real,
)
    return _accept_cauer_thermal_step!(
        thermal,
        Float64(loss_power_w),
        Float64(step_s),
    )
end

function _deposit_cauer_event_energy!(
    thermal::CauerThermalFidelity,
    energy_j::Float64,
)
    isfinite(energy_j) && energy_j >= 0.0 || throw(ArgumentError(
        "thermal event energy must be finite and nonnegative",
    ))
    temperature = thermal.node_temperature_k[1] +
        energy_j / thermal.capacitance_j_per_k[1]
    thermal.minimum_temperature_k <= temperature <= thermal.maximum_temperature_k ||
        throw(DomainError(temperature, "event energy leaves the thermal domain"))
    thermal.node_temperature_k[1] = temperature
    thermal.cumulative_input_energy_j += energy_j
    thermal.last_stored_energy_j = power_semiconductor_thermal_stored_energy(thermal)
    return thermal
end

"""Hash- and rights-bound identity for one supported AIMORA semiconductor equation schema."""
struct DeclaredSemiconductorModelIdentity
    schema::Symbol
    version::VersionNumber
    source_identity::String
    content_sha256::String
    licence::String
    redistribution::Symbol
    provenance::NonlinearParameterProvenance

    function DeclaredSemiconductorModelIdentity(
        schema::Symbol,
        version::VersionNumber,
        source_identity::AbstractString,
        content_sha256::AbstractString,
        licence::AbstractString,
        redistribution::Symbol,
        provenance::NonlinearParameterProvenance,
    )
        schema in (:aimora_generic_semiconductor_v1,) || throw(ArgumentError(
            "unsupported declared semiconductor equation schema $schema",
        ))
        isempty(strip(source_identity)) && throw(ArgumentError(
            "declared semiconductor source identity must not be empty",
        ))
        occursin(r"^[0-9a-f]{64}$", lowercase(content_sha256)) || throw(ArgumentError(
            "declared semiconductor content hash must be 64 hexadecimal characters",
        ))
        isempty(strip(licence)) && throw(ArgumentError(
            "declared semiconductor licence must not be empty",
        ))
        redistribution in (:public, :private, :prohibited) || throw(ArgumentError(
            "declared semiconductor redistribution must be :public, :private, or :prohibited",
        ))
        _validate_semiconductor_provenance(provenance, "declared-semiconductor-model")
        return new(
            schema,
            version,
            String(source_identity),
            lowercase(String(content_sha256)),
            String(licence),
            redistribution,
            provenance,
        )
    end
end

"""Compositional extended fidelity attached to one accepted semiconductor owner."""
mutable struct PowerSemiconductorExtendedFidelity
    recovered_charge::Union{Nothing,RecoveredChargeFidelity}
    junction_charge::Union{Nothing,NonlinearJunctionChargeFidelity}
    turn_off_tail::Union{Nothing,TurnOffTailFidelity}
    switching_energy::Union{Nothing,SwitchingEnergyTable}
    thermal::Union{Nothing,CauerThermalFidelity}
    declared_model::Union{Nothing,DeclaredSemiconductorModelIdentity}
    provenance::NonlinearParameterProvenance
    candidate_time_s::Float64
    candidate_step_s::Float64
    candidate_method::Symbol
    candidate_prepared::Bool
    previous_terminal_voltage_v::Float64
    previous_terminal_current_a::Float64
    accepted_topology_transition_count::Int
    pending_event_current_a::Float64
    pending_event_blocking_voltage_v::Float64
    candidate_recovery_charge_c::Float64
end

function PowerSemiconductorExtendedFidelity(;
    recovered_charge::Union{Nothing,RecoveredChargeFidelity}=nothing,
    junction_charge::Union{Nothing,NonlinearJunctionChargeFidelity}=nothing,
    turn_off_tail::Union{Nothing,TurnOffTailFidelity}=nothing,
    switching_energy::Union{Nothing,SwitchingEnergyTable}=nothing,
    thermal::Union{Nothing,CauerThermalFidelity}=nothing,
    declared_model::Union{Nothing,DeclaredSemiconductorModelIdentity}=nothing,
    provenance::NonlinearParameterProvenance=_generic_semiconductor_provenance(
        "component-specific SI units",
        "composition of explicitly selected generic semiconductor fidelity components",
    ),
)
    any(!isnothing, (
        recovered_charge,
        junction_charge,
        turn_off_tail,
        switching_energy,
        thermal,
        declared_model,
    )) || throw(ArgumentError(
        "extended semiconductor fidelity must select at least one component",
    ))
    _validate_semiconductor_provenance(provenance, "extended-semiconductor")
    return PowerSemiconductorExtendedFidelity(
        recovered_charge,
        junction_charge,
        turn_off_tail,
        switching_energy,
        thermal,
        declared_model,
        provenance,
        0.0,
        0.0,
        :unprepared,
        false,
        junction_charge === nothing ? 0.0 : junction_charge.previous_voltage_v,
        0.0,
        0,
        0.0,
        0.0,
        recovered_charge === nothing ? 0.0 : recovered_charge.stored_charge_c,
    )
end

"""Versioned extended physical state that leaves the accepted baseline terminal state unchanged."""
struct PowerSemiconductorExtendedState
    schema_version::Int
    fidelity_components::Tuple{Vararg{Symbol}}
    conduction_direction::Int8
    gate_turn_off_disposition::Symbol
    recovery_active::Bool
    stored_recovery_charge_c::Float64
    recovery_current_a::Float64
    peak_reverse_current_a::Float64
    cumulative_recovered_charge_c::Float64
    recovery_zero_event_count::Int
    last_recovery_zero_time_s::Float64
    junction_capacitance_f::Float64
    junction_charge_c::Float64
    displacement_current_a::Float64
    junction_stored_energy_j::Float64
    tail_active::Bool
    tail_current_a::Float64
    tail_cutoff_event_count::Int
    last_tail_cutoff_time_s::Float64
    last_event_kind::Symbol
    last_event_energy_j::Float64
    cumulative_turn_on_energy_j::Float64
    cumulative_turn_off_energy_j::Float64
    cumulative_reverse_recovery_energy_j::Float64
    junction_temperature_k::Float64
    thermal_node_temperature_k::Vector{Float64}
    ambient_heat_flow_w::Float64
    thermal_stored_energy_j::Float64
end

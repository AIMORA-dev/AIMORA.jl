module SwitchDetailedVSC

using ..StudyCore: ContractQuantity,
                   DynamicStateInventory,
                   ModelValidityDomain,
                   NumericDomainBound,
                   ScientificModelContract,
                   SwitchingDetailed,
                   assess_validity,
                   assert_validity

export PulseWidthModulationKind,
       SinusoidalPulseWidthModulation,
       ZeroSequenceInjectedPulseWidthModulation,
       StationaryReferenceFrame,
       SynchronousReferenceFrame,
       ThreePhaseTwoLevelVSCParameters,
       ThreePhaseGridFollowingControllerState,
       ThreePhaseVSCMeasurement,
       ThreePhaseVSCModulationCommand,
       switch_detailed_vsc_contract,
       validate_three_phase_vsc_parameters,
       clarke_transform,
       inverse_clarke_transform,
       park_transform,
       inverse_park_transform,
       instantaneous_three_phase_power,
       modulation_duties,
       compute_grid_following_current_control!

"""Carrier modulation families: phase sinusoidal PWM or minimum-maximum zero-sequence injection with centered-space-vector-equivalent line voltages in the linear region."""
@enum PulseWidthModulationKind begin
    SinusoidalPulseWidthModulation
    ZeroSequenceInjectedPulseWidthModulation
end

"""Amplitude-invariant Clarke coordinates. With a zero-sequence-free balanced sinusoid, the stationary-vector magnitude equals the phase crest."""
struct StationaryReferenceFrame
    alpha::Float64
    beta::Float64
    zero::Float64
end

"""Synchronous d-q-zero coordinates using a positive-sequence angle measured from the stationary alpha axis."""
struct SynchronousReferenceFrame
    direct::Float64
    quadrature::Float64
    zero::Float64
end

"""Generic, redistributable parameters for one grounded-wye, three-wire, switch-detailed two-level grid-following VSC slice."""
Base.@kwdef struct ThreePhaseTwoLevelVSCParameters
    frequency_hz::Float64 = 50.0
    rated_power_va::Float64 = 50.0e3
    grid_line_line_rms_v::Float64 = 400.0
    transformer_grid_to_converter_ratio::Float64 = 1.0
    filter_resistance_ohm::Float64 = 0.12
    filter_inductance_h::Float64 = 2.0e-3
    transformer_leakage_resistance_ohm::Float64 = 0.03
    transformer_leakage_inductance_h::Float64 = 0.30e-3
    grid_source_resistance_ohm::Float64 = 0.08
    dc_source_voltage_v::Float64 = 800.0
    dc_source_resistance_ohm::Float64 = 0.60
    dc_link_capacitance_f::Float64 = 10.0e-3
    active_power_reference_w::Float64 = 20.0e3
    reactive_power_reference_var::Float64 = 0.0
    current_controller_proportional_gain_v_per_a::Float64 = 4.0
    current_controller_integral_gain_v_per_as::Float64 = 600.0
    current_limit_a::Float64 = 180.0
    modulation::PulseWidthModulationKind = ZeroSequenceInjectedPulseWidthModulation
    minimum_duty::Float64 = 0.02
    maximum_duty::Float64 = 0.98
    carrier_frequency_hz::Float64 = 10.0e3
    control_period_s::Float64 = 100.0e-6
    control_delay_s::Float64 = 2.0e-6
    scheduler_tick_s::Float64 = 1.0e-6
    timestep_s::Float64 = 1.0e-6
    end_time_s::Float64 = 60.0e-3
    gate_turn_on_delay_s::Float64 = 0.0
    gate_turn_off_delay_s::Float64 = 0.0
    commutation_dead_time_s::Float64 = 2.0e-6
    minimum_gate_pulse_width_s::Float64 = 1.0e-6
    semiconductor_on_conductance_s::Float64 = 200.0
    semiconductor_off_conductance_s::Float64 = 1.0e-6
    semiconductor_forward_voltage_v::Float64 = 1.2
    diode_on_conductance_s::Float64 = 200.0
    diode_forward_voltage_v::Float64 = 0.9
    sag_start_s::Float64 = 10.0e-3
    sag_end_s::Float64 = 15.0e-3
    sag_voltage_factor::Float64 = 0.70
    fault_start_s::Float64 = 25.0e-3
    fault_end_s::Float64 = 27.0e-3
    faulted_phase::Int = 1
    fault_voltage_factor::Float64 = 0.08
    block_time_s::Float64 = 25.5e-3
    restart_time_s::Float64 = 28.0e-3
    harmonic_window_start_s::Float64 = 40.0e-3
    harmonic_window_end_s::Float64 = 60.0e-3
    maximum_harmonic_order::Int = 50
end

"""Explicit continuous and discrete state of the synchronous-reference-frame current controller."""
mutable struct ThreePhaseGridFollowingControllerState
    direct_integral_as::Float64
    quadrature_integral_as::Float64
    held_duties::NTuple{3,Float64}
    held_phase_voltage_reference_v::NTuple{3,Float64}
    sample_count::Int
    write_count::Int
    saturation_count::Int
    overcurrent_count::Int
end

ThreePhaseGridFollowingControllerState() = ThreePhaseGridFollowingControllerState(
    0.0,
    0.0,
    (0.5, 0.5, 0.5),
    (0.0, 0.0, 0.0),
    0,
    0,
    0,
    0,
)

"""One exact sampled-controller input with converter-side grid terminal voltage, filter current, DC-link voltage, and synchronous angle."""
struct ThreePhaseVSCMeasurement
    phase_voltage_v::NTuple{3,Float64}
    phase_current_a::NTuple{3,Float64}
    dc_link_voltage_v::Float64
    synchronous_angle_rad::Float64
end

"""One held three-phase modulation result produced by the sampled current controller."""
struct ThreePhaseVSCModulationCommand
    duties::NTuple{3,Float64}
    phase_voltage_reference_v::NTuple{3,Float64}
    direct_current_reference_a::Float64
    quadrature_current_reference_a::Float64
    direct_current_a::Float64
    quadrature_current_a::Float64
    saturated::Bool
end

const SWITCH_DETAILED_VSC_CONTRACT = ScientificModelContract(
    :three_phase_two_level_switch_detailed_vsc,
    :synchronous_reference_frame_grid_following_current_control;
    owner = "AIMORA.SwitchDetailedVSC and AIMORA.EMTStudy",
    maturity = :implemented,
    fidelity = SwitchingDetailed,
    validity_domain = ModelValidityDomain(
        :grounded_wye_three_wire_grid_following_vsc;
        description = "Three-phase, three-wire, two-level IGBT VSC with ideal zero-recovery antiparallel diodes, a dynamic DC capacitor supplied through a finite Thevenin resistance, exact sampled synchronous-reference-frame current control, exact trailing-edge carrier PWM, a converter-side series L filter, and a grounded-wye ideal-transformer/grid interface represented on the converter side.",
        bounds = (
            NumericDomainBound(:frequency_hz; unit = "Hz", lower = 45.0, upper = 65.0),
            NumericDomainBound(:rated_power_va; unit = "VA", lower = 0.0, lower_inclusive = false),
            NumericDomainBound(:grid_line_line_rms_v; unit = "V", lower = 320.0, upper = 440.0),
            NumericDomainBound(:dc_source_voltage_v; unit = "V", lower = 650.0, upper = 900.0),
            NumericDomainBound(:dc_link_capacitance_f; unit = "F", lower = 0.0, lower_inclusive = false),
            NumericDomainBound(:filter_inductance_h; unit = "H", lower = 0.0, lower_inclusive = false),
            NumericDomainBound(:carrier_frequency_hz; unit = "Hz", lower = 0.0, lower_inclusive = false),
            NumericDomainBound(:control_period_s; unit = "s", lower = 0.0, lower_inclusive = false),
            NumericDomainBound(:current_limit_a; unit = "A", lower = 0.0, lower_inclusive = false),
            NumericDomainBound(:timestep_s; unit = "s", lower = 0.0, upper = 1.0e-6, lower_inclusive = false),
            NumericDomainBound(:end_time_s; unit = "s", lower = 0.0, lower_inclusive = false),
        ),
        unsupported_phenomena = (
            :phase_locked_loop_dynamics,
            :grid_forming_control,
            :four_wire_zero_sequence_current,
            :lcl_filter_resonance,
            :transformer_magnetizing_saturation,
            :semiconductor_reverse_recovery,
            :nonlinear_device_capacitance,
            :switching_energy_maps,
            :electrothermal_state,
            :manufacturer_parameter_prediction,
            :adaptive_global_timestep,
            :standard_conformance_or_certification,
        ),
    ),
    state_inventory = DynamicStateInventory(
        differential = (
            :dc_link_capacitor_voltage,
            :phase_filter_currents,
            :direct_current_integral,
            :quadrature_current_integral,
        ),
        algebraic = (
            :nodal_voltages,
            :device_terminal_currents,
            :instantaneous_active_power,
            :instantaneous_reactive_power,
        ),
        discrete = (
            :six_gate_commands,
            :six_conduction_states,
            :six_antiparallel_diode_states,
            :blocked_state,
            :protection_counters,
        ),
        delayed_history = (
            :filter_companion_histories,
            :dc_capacitor_companion_history,
            :sampled_control_pending_queue,
        ),
        scheduler = (
            :control_sample_tick,
            :control_release_tick,
            :carrier_boundary_ticks,
            :pwm_edge_ticks,
            :disturbance_event_ticks,
        ),
    ),
    inputs = (
        ContractQuantity(:grid_phase_voltage_v; unit = "V", orientation = "grid_neutral_to_phase"),
        ContractQuantity(:dc_source_voltage_v; unit = "V", orientation = "dc_negative_to_positive"),
        ContractQuantity(:active_power_reference_w; unit = "W", orientation = "positive_converter_to_grid"),
        ContractQuantity(:reactive_power_reference_var; unit = "var", orientation = "positive_converter_to_grid"),
    ),
    outputs = (
        ContractQuantity(:pole_voltage_v; unit = "V", orientation = "dc_negative_to_pole"),
        ContractQuantity(:line_voltage_v; unit = "V", orientation = "ab_bc_ca"),
        ContractQuantity(:filter_current_a; unit = "A", orientation = "converter_to_grid"),
        ContractQuantity(:dc_link_voltage_v; unit = "V", orientation = "dc_negative_to_positive"),
        ContractQuantity(:active_power_w; unit = "W", orientation = "positive_converter_to_grid"),
        ContractQuantity(:reactive_power_var; unit = "var", orientation = "positive_converter_to_grid"),
        ContractQuantity(:device_state; unit = "state"),
        ContractQuantity(:terminal_energy_j; unit = "J"),
        ContractQuantity(:dc_source_energy_j; unit = "J", orientation = "positive_source_to_converter"),
        ContractQuantity(:ac_terminal_energy_j; unit = "J", orientation = "positive_converter_to_grid"),
        ContractQuantity(:dc_ac_energy_residual_j; unit = "J", orientation = "dc_input_minus_ac_output_loss_and_storage"),
        ContractQuantity(:harmonic_distortion; unit = "ratio"),
    ),
    assumptions = (
        "The synchronous angle is supplied by the known ideal grid-source phase; this release does not contain a phase-locked loop.",
        "The minimum-maximum zero-sequence-injected carrier formulation is line-voltage-equivalent to centered space-vector PWM in the linear modulation region; no separate vector dwell-time sequencer is claimed.",
        "The grounded-wye transformer is an ideal ratio with leakage resistance and inductance referred to the converter side; magnetizing and saturation branches are excluded.",
        "The canonical fault is a prescribed source-side phase-to-ground voltage collapse behind the declared grid Thevenin resistance.",
        "Semiconductor and diode behavior retains the generic piecewise-linear, ideal zero-recovery limitations of the accepted device and bridge owners.",
    ),
    mutation_order = (
        :capture_reversible_state,
        :apply_exact_disturbance_and_protection_boundaries,
        :sample_grid_voltage_and_filter_current,
        :compute_synchronous_reference_frame_control,
        :release_delayed_modulation_command,
        :apply_exact_carrier_edges_through_bridge_interlocks,
        :stamp_devices_dc_link_filter_and_grid,
        :solve_three_phase_nodal_network,
        :update_companion_device_and_energy_state,
        :record_typed_output,
        :commit_or_restore,
    ),
)

switch_detailed_vsc_contract() = SWITCH_DETAILED_VSC_CONTRACT

function _finite_positive(value::Real, name::AbstractString)
    isfinite(value) && value > 0.0 || throw(ArgumentError("$name must be finite and positive"))
    return Float64(value)
end

function _finite_nonnegative(value::Real, name::AbstractString)
    isfinite(value) && value >= 0.0 || throw(ArgumentError("$name must be finite and nonnegative"))
    return Float64(value)
end

function _is_integer_multiple(value::Float64, quantum::Float64)
    represented = round(value / quantum) * quantum
    tolerance = 64.0 * eps(Float64) * max(abs(value), quantum)
    return abs(represented - value) <= tolerance
end

"""Validate every released converter domain, timing, modulation, disturbance, and analysis-window bound."""
function validate_three_phase_vsc_parameters(parameters::ThreePhaseTwoLevelVSCParameters)
    p = parameters
    assessment = assess_validity(
        switch_detailed_vsc_contract(),
        (
            frequency_hz = p.frequency_hz,
            rated_power_va = p.rated_power_va,
            grid_line_line_rms_v = p.grid_line_line_rms_v,
            dc_source_voltage_v = p.dc_source_voltage_v,
            dc_link_capacitance_f = p.dc_link_capacitance_f,
            filter_inductance_h = p.filter_inductance_h,
            carrier_frequency_hz = p.carrier_frequency_hz,
            control_period_s = p.control_period_s,
            current_limit_a = p.current_limit_a,
            timestep_s = p.timestep_s,
            end_time_s = p.end_time_s,
        );
        requested_fidelity = SwitchingDetailed,
    )
    assert_validity(assessment)
    for (value, name) in (
        (p.transformer_grid_to_converter_ratio, "transformer ratio"),
        (p.grid_source_resistance_ohm, "grid-source resistance"),
        (p.dc_source_resistance_ohm, "DC-source resistance"),
        (p.current_limit_a, "current limit"),
        (p.scheduler_tick_s, "scheduler tick"),
        (p.semiconductor_on_conductance_s, "semiconductor on conductance"),
        (p.diode_on_conductance_s, "diode on conductance"),
    )
        _finite_positive(value, name)
    end
    for (value, name) in (
        (p.filter_resistance_ohm, "filter resistance"),
        (p.transformer_leakage_resistance_ohm, "transformer leakage resistance"),
        (p.transformer_leakage_inductance_h, "transformer leakage inductance"),
        (p.current_controller_proportional_gain_v_per_a, "controller proportional gain"),
        (p.current_controller_integral_gain_v_per_as, "controller integral gain"),
        (p.control_delay_s, "control delay"),
        (p.gate_turn_on_delay_s, "gate turn-on delay"),
        (p.gate_turn_off_delay_s, "gate turn-off delay"),
        (p.commutation_dead_time_s, "commutation dead time"),
        (p.minimum_gate_pulse_width_s, "minimum gate pulse width"),
        (p.semiconductor_off_conductance_s, "semiconductor off conductance"),
        (p.semiconductor_forward_voltage_v, "semiconductor forward voltage"),
        (p.diode_forward_voltage_v, "diode forward voltage"),
    )
        _finite_nonnegative(value, name)
    end
    0.0 <= p.minimum_duty < p.maximum_duty <= 1.0 || throw(ArgumentError(
        "VSC duty bounds must satisfy 0 <= minimum < maximum <= 1",
    ))
    0.0 < p.sag_voltage_factor <= 1.0 || throw(ArgumentError(
        "sag voltage factor must be within (0, 1]",
    ))
    0.0 <= p.fault_voltage_factor <= 1.0 || throw(ArgumentError(
        "fault voltage factor must be within [0, 1]",
    ))
    1 <= p.faulted_phase <= 3 || throw(ArgumentError(
        "faulted phase must be 1, 2, or 3",
    ))
    0.0 <= p.sag_start_s < p.sag_end_s <= p.end_time_s || throw(ArgumentError(
        "sag interval must be ordered inside the simulation horizon",
    ))
    0.0 <= p.fault_start_s < p.fault_end_s <= p.end_time_s || throw(ArgumentError(
        "fault interval must be ordered inside the simulation horizon",
    ))
    p.fault_start_s <= p.block_time_s < p.fault_end_s || throw(ArgumentError(
        "block time must lie inside the fault interval",
    ))
    p.fault_end_s <= p.restart_time_s <= p.end_time_s || throw(ArgumentError(
        "restart time must follow fault clearance inside the horizon",
    ))
    0.0 <= p.harmonic_window_start_s < p.harmonic_window_end_s <= p.end_time_s ||
        throw(ArgumentError("harmonic window must be ordered inside the horizon"))
    window_cycles = (p.harmonic_window_end_s - p.harmonic_window_start_s) * p.frequency_hz
    abs(window_cycles - round(window_cycles)) <= 1.0e-10 || throw(ArgumentError(
        "harmonic window must span an integer number of fundamental cycles",
    ))
    p.maximum_harmonic_order >= 2 || throw(ArgumentError(
        "maximum harmonic order must be at least two",
    ))
    carrier_period_s = inv(p.carrier_frequency_hz)
    for (value, name) in (
        (p.timestep_s, "timestep"),
        (p.control_period_s, "control period"),
        (p.control_delay_s, "control delay"),
        (carrier_period_s, "carrier period"),
        (p.gate_turn_on_delay_s, "gate turn-on delay"),
        (p.gate_turn_off_delay_s, "gate turn-off delay"),
        (p.commutation_dead_time_s, "commutation dead time"),
        (p.minimum_gate_pulse_width_s, "minimum gate pulse width"),
        (p.sag_start_s, "sag start"),
        (p.sag_end_s, "sag end"),
        (p.fault_start_s, "fault start"),
        (p.fault_end_s, "fault end"),
        (p.block_time_s, "block time"),
        (p.restart_time_s, "restart time"),
    )
        _is_integer_multiple(value, p.scheduler_tick_s) || throw(ArgumentError(
            "$name must be an integer multiple of scheduler_tick_s",
        ))
    end
    _is_integer_multiple(p.end_time_s, p.timestep_s) || throw(ArgumentError(
        "end time must be an integer multiple of timestep_s",
    ))
    _is_integer_multiple(p.timestep_s, p.scheduler_tick_s) || throw(ArgumentError(
        "timestep_s must be an integer multiple of scheduler_tick_s",
    ))
    return p
end

function clarke_transform(a::Real, b::Real, c::Real)
    phase_a = Float64(a)
    phase_b = Float64(b)
    phase_c = Float64(c)
    all(isfinite, (phase_a, phase_b, phase_c)) || throw(ArgumentError(
        "Clarke-transform inputs must be finite",
    ))
    return StationaryReferenceFrame(
        (2.0 / 3.0) * (phase_a - 0.5 * phase_b - 0.5 * phase_c),
        (sqrt(3.0) / 3.0) * (phase_b - phase_c),
        (phase_a + phase_b + phase_c) / 3.0,
    )
end

clarke_transform(values::NTuple{3,<:Real}) = clarke_transform(values...)

function inverse_clarke_transform(frame::StationaryReferenceFrame)
    half_sqrt_three_beta = 0.5 * sqrt(3.0) * frame.beta
    return (
        frame.alpha + frame.zero,
        -0.5 * frame.alpha + half_sqrt_three_beta + frame.zero,
        -0.5 * frame.alpha - half_sqrt_three_beta + frame.zero,
    )
end

function park_transform(frame::StationaryReferenceFrame, angle_rad::Real)
    angle = Float64(angle_rad)
    isfinite(angle) || throw(ArgumentError("Park-transform angle must be finite"))
    cosine = cos(angle)
    sine = sin(angle)
    return SynchronousReferenceFrame(
        cosine * frame.alpha + sine * frame.beta,
        -sine * frame.alpha + cosine * frame.beta,
        frame.zero,
    )
end

function inverse_park_transform(frame::SynchronousReferenceFrame, angle_rad::Real)
    angle = Float64(angle_rad)
    isfinite(angle) || throw(ArgumentError("inverse-Park angle must be finite"))
    cosine = cos(angle)
    sine = sin(angle)
    return StationaryReferenceFrame(
        cosine * frame.direct - sine * frame.quadrature,
        sine * frame.direct + cosine * frame.quadrature,
        frame.zero,
    )
end

function instantaneous_three_phase_power(
    phase_voltage_v::NTuple{3,<:Real},
    phase_current_a::NTuple{3,<:Real},
    angle_rad::Real,
)
    voltage = park_transform(clarke_transform(phase_voltage_v), angle_rad)
    current = park_transform(clarke_transform(phase_current_a), angle_rad)
    active = 1.5 * (
        voltage.direct * current.direct +
        voltage.quadrature * current.quadrature
    ) + 3.0 * voltage.zero * current.zero
    reactive = 1.5 * (
        voltage.quadrature * current.direct -
        voltage.direct * current.quadrature
    )
    return (active_w = active, reactive_var = reactive)
end

function modulation_duties(
    phase_voltage_reference_v::NTuple{3,<:Real},
    dc_link_voltage_v::Real,
    modulation::PulseWidthModulationKind;
    minimum_duty::Real = 0.0,
    maximum_duty::Real = 1.0,
)
    dc_voltage = _finite_positive(dc_link_voltage_v, "DC-link voltage")
    lower = Float64(minimum_duty)
    upper = Float64(maximum_duty)
    0.0 <= lower < upper <= 1.0 || throw(ArgumentError(
        "modulation duty bounds must satisfy 0 <= minimum < maximum <= 1",
    ))
    references = Float64.(phase_voltage_reference_v)
    all(isfinite, references) || throw(ArgumentError(
        "phase-voltage modulation references must be finite",
    ))
    common_mode = modulation == SinusoidalPulseWidthModulation ? 0.0 :
        -0.5 * (maximum(references) + minimum(references))
    duties = ntuple(3) do phase
        clamp(0.5 + (references[phase] + common_mode) / dc_voltage, lower, upper)
    end
    return duties
end

"""Apply one sampled synchronous-reference-frame PI current-control update and return the held phase-voltage/duty command."""
function compute_grid_following_current_control!(
    state::ThreePhaseGridFollowingControllerState,
    measurement::ThreePhaseVSCMeasurement,
    parameters::ThreePhaseTwoLevelVSCParameters,
)
    p = validate_three_phase_vsc_parameters(parameters)
    voltage_dq = park_transform(
        clarke_transform(measurement.phase_voltage_v),
        measurement.synchronous_angle_rad,
    )
    current_dq = park_transform(
        clarke_transform(measurement.phase_current_a),
        measurement.synchronous_angle_rad,
    )
    direct_voltage = max(abs(voltage_dq.direct), 0.05 * p.grid_line_line_rms_v)
    direct_reference = (2.0 / 3.0) * p.active_power_reference_w / direct_voltage
    quadrature_reference = -(2.0 / 3.0) * p.reactive_power_reference_var / direct_voltage
    reference_magnitude = hypot(direct_reference, quadrature_reference)
    if reference_magnitude > p.current_limit_a
        scale = p.current_limit_a / reference_magnitude
        direct_reference *= scale
        quadrature_reference *= scale
        state.overcurrent_count += 1
    end
    direct_error = direct_reference - current_dq.direct
    quadrature_error = quadrature_reference - current_dq.quadrature
    direct_integral_candidate = state.direct_integral_as + p.control_period_s * direct_error
    quadrature_integral_candidate =
        state.quadrature_integral_as + p.control_period_s * quadrature_error
    omega = 2.0 * pi * p.frequency_hz
    series_resistance = p.filter_resistance_ohm + p.transformer_leakage_resistance_ohm
    series_inductance = p.filter_inductance_h + p.transformer_leakage_inductance_h
    direct_reference_voltage = voltage_dq.direct +
        series_resistance * current_dq.direct -
        omega * series_inductance * current_dq.quadrature +
        p.current_controller_proportional_gain_v_per_a * direct_error +
        p.current_controller_integral_gain_v_per_as * direct_integral_candidate
    quadrature_reference_voltage = voltage_dq.quadrature +
        series_resistance * current_dq.quadrature +
        omega * series_inductance * current_dq.direct +
        p.current_controller_proportional_gain_v_per_a * quadrature_error +
        p.current_controller_integral_gain_v_per_as * quadrature_integral_candidate
    reference_magnitude_v = hypot(direct_reference_voltage, quadrature_reference_voltage)
    voltage_limit_v = 0.95 * measurement.dc_link_voltage_v / sqrt(3.0)
    saturated = reference_magnitude_v > voltage_limit_v
    if saturated && reference_magnitude_v > 0.0
        scale = voltage_limit_v / reference_magnitude_v
        direct_reference_voltage *= scale
        quadrature_reference_voltage *= scale
        state.saturation_count += 1
    else
        state.direct_integral_as = direct_integral_candidate
        state.quadrature_integral_as = quadrature_integral_candidate
    end
    phase_reference = inverse_clarke_transform(inverse_park_transform(
        SynchronousReferenceFrame(
            direct_reference_voltage,
            quadrature_reference_voltage,
            0.0,
        ),
        measurement.synchronous_angle_rad,
    ))
    duties = modulation_duties(
        phase_reference,
        measurement.dc_link_voltage_v,
        p.modulation;
        minimum_duty = p.minimum_duty,
        maximum_duty = p.maximum_duty,
    )
    state.sample_count += 1
    return ThreePhaseVSCModulationCommand(
        duties,
        phase_reference,
        direct_reference,
        quadrature_reference,
        current_dq.direct,
        current_dq.quadrature,
        saturated,
    )
end

end

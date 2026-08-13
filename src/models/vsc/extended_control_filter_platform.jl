export ExtendedVSCControllerFamily,
       SynchronousPLLGridFollowing,
       StationaryResonantGridFollowing,
       PowerDroopGridForming,
       VirtualSynchronousGridForming,
       ExtendedVSCFilterFamily,
       SeriesLFilter,
       ShuntLCFilter,
       LCLFilter,
       ExtendedVSCWireForm,
       ThreeWireForm,
       FourWireForm,
       VSCCurrentPriority,
       ActiveCurrentPriority,
       ReactiveCurrentPriority,
       VectorMagnitudePriority,
       VSCPLLLossPolicy,
       FreezePLLOnVoltageLoss,
       BlockVSCOnVoltageLoss,
       ExtendedVSCOperatingMode,
       VSCNormalOperation,
       VSCCurrentLimitedOperation,
       VSCBlockedOperation,
       VSCPlantRequestDisposition,
       VSCPlantRequestApplied,
       VSCPlantRequestLimited,
       VSCPlantRequestStale,
       VSCPlantRequestRefused,
       ExtendedVSCFilterParameters,
       ExtendedVSCControlParameters,
       ExtendedVSCProtectionParameters,
       ExtendedVSCScenarioParameters,
       ExtendedVSCParameters,
       ExtendedVSCPlantRequest,
       ExtendedVSCMeasurement,
       VSCRotatingSequenceState,
       ExtendedVSCControlState,
       ExtendedVSCControlCommand,
       extended_vsc_contract,
       validate_extended_vsc_parameters,
       supported_extended_vsc_combination,
       advance_vsc_sequence_extractor!,
       project_vsc_current_reference,
       extended_vsc_modulation_duties,
       compute_extended_vsc_control!

@enum ExtendedVSCControllerFamily begin
    SynchronousPLLGridFollowing
    StationaryResonantGridFollowing
    PowerDroopGridForming
    VirtualSynchronousGridForming
end

@enum ExtendedVSCFilterFamily begin
    SeriesLFilter
    ShuntLCFilter
    LCLFilter
end

@enum ExtendedVSCWireForm begin
    ThreeWireForm
    FourWireForm
end

@enum VSCCurrentPriority begin
    ActiveCurrentPriority
    ReactiveCurrentPriority
    VectorMagnitudePriority
end

@enum VSCPLLLossPolicy begin
    FreezePLLOnVoltageLoss
    BlockVSCOnVoltageLoss
end

@enum ExtendedVSCOperatingMode begin
    VSCNormalOperation
    VSCCurrentLimitedOperation
    VSCBlockedOperation
end

@enum VSCPlantRequestDisposition begin
    VSCPlantRequestApplied
    VSCPlantRequestLimited
    VSCPlantRequestStale
    VSCPlantRequestRefused
end

const _EXTENDED_VSC_HARMONIC_ORDERS = (1, 5, 7, 11, 13)

function _generic_vsc_provenance(owner::AbstractString, units::AbstractString)
    return ParameterProvenance(
        "AIMORA-authored synthetic generic $owner parameters",
        units,
        "direct SI values; no hidden per-unit or vendor conversion",
        "deterministic generic input with explicit unknown physical uncertainty",
        "REQ-EMT-VSC-002 frozen generic validity domain",
        PhysicalModelParameter,
    )
end

function _validate_vsc_physical_provenance(
    provenance::ParameterProvenance,
    owner::AbstractString,
)
    provenance.nature === PhysicalModelParameter || throw(ArgumentError(
        "$owner provenance must describe physical model parameters",
    ))
    return provenance
end

Base.@kwdef struct ExtendedVSCFilterParameters
    family::ExtendedVSCFilterFamily = SeriesLFilter
    converter_resistance_ohm::Float64 = 0.12
    converter_inductance_h::Float64 = 2.0e-3
    shunt_capacitance_f::Float64 = 20.0e-6
    shunt_damping_conductance_s::Float64 = 0.02
    grid_resistance_ohm::Float64 = 0.05
    grid_inductance_h::Float64 = 0.50e-3
    neutral_resistance_ohm::Float64 = 0.08
    neutral_inductance_h::Float64 = 0.50e-3
    provenance::ParameterProvenance = _generic_vsc_provenance(
        "VSC filter and neutral",
        "ohm, henry, farad, siemens",
    )
end

Base.@kwdef struct ExtendedVSCControlParameters
    family::ExtendedVSCControllerFamily = SynchronousPLLGridFollowing
    current_priority::VSCCurrentPriority = ActiveCurrentPriority
    pll_loss_policy::VSCPLLLossPolicy = BlockVSCOnVoltageLoss
    pll_proportional_gain_rad_per_s::Float64 = 110.0
    pll_integral_gain_rad_per_s2::Float64 = 4.0e3
    pll_voltage_floor_v::Float64 = 20.0
    pll_minimum_frequency_hz::Float64 = 45.0
    pll_maximum_frequency_hz::Float64 = 65.0
    current_proportional_gain_v_per_a::Float64 = 4.0
    current_integral_gain_v_per_as::Float64 = 600.0
    current_antiwindup_gain_per_s::Float64 = 120.0
    voltage_proportional_gain_a_per_v::Float64 = 0.20
    voltage_integral_gain_a_per_vs::Float64 = 20.0
    resonant_gain_v_per_a::Float64 = 80.0
    resonant_bandwidth_rad_per_s::Float64 = 12.0
    selected_harmonic_orders::Tuple{Vararg{Int}} = (1, 5, 7)
    power_filter_time_constant_s::Float64 = 2.0e-3
    active_power_frequency_droop_rad_per_ws::Float64 = 3.0e-5
    reactive_power_voltage_droop_v_per_var::Float64 = 2.0e-3
    virtual_inertia_w_s2_per_rad::Float64 = 4.0
    virtual_damping_w_s_per_rad::Float64 = 120.0
    virtual_resistance_ohm::Float64 = 0.08
    virtual_inductance_h::Float64 = 0.40e-3
    virtual_current_filter_time_constant_s::Float64 = 0.50e-3
    control_period_s::Float64 = 50.0e-6
    control_delay_s::Float64 = 2.0e-6
    provenance::ParameterProvenance = _generic_vsc_provenance(
        "VSC controller",
        "field-specific SI units",
    )
end

Base.@kwdef struct ExtendedVSCProtectionParameters
    ac_overcurrent_a::Float64 = 190.0
    dc_undervoltage_v::Float64 = 600.0
    dc_overvoltage_v::Float64 = 950.0
    minimum_frequency_hz::Float64 = 44.0
    maximum_frequency_hz::Float64 = 66.0
    minimum_phase_rms_voltage_v::Float64 = 120.0
    maximum_phase_rms_voltage_v::Float64 = 500.0
    trip_delay_s::Float64 = 100.0e-6
    restart_delay_s::Float64 = 500.0e-6
    provenance::ParameterProvenance = _generic_vsc_provenance(
        "generic VSC protection",
        "ampere, volt, hertz, second",
    )
end

Base.@kwdef struct ExtendedVSCScenarioParameters
    negative_sequence_voltage_ratio::Float64 = 0.05
    zero_sequence_voltage_ratio::Float64 = 0.0
    harmonic_orders::Tuple{Vararg{Int}} = (5, 7)
    harmonic_voltage_ratios::Tuple{Vararg{Float64}} = (0.02, 0.015)
    fault_start_s::Float64 = 12.0e-3
    fault_end_s::Float64 = 15.0e-3
    fault_kind::Symbol = :phase_to_ground
    faulted_phases::Tuple{Vararg{Int}} = (1,)
    fault_voltage_factor::Float64 = 0.15
    dc_sag_start_s::Float64 = 20.0e-3
    dc_sag_end_s::Float64 = 23.0e-3
    dc_sag_factor::Float64 = 0.80
    island_start_s::Float64 = 28.0e-3
    reconnect_time_s::Float64 = 34.0e-3
    block_time_s::Float64 = 22.0e-3
    restart_time_s::Float64 = 25.0e-3
    end_time_s::Float64 = 40.0e-3
    provenance::ParameterProvenance = _generic_vsc_provenance(
        "VSC disturbance scenario",
        "ratio and second",
    )
end

Base.@kwdef struct ExtendedVSCParameters
    controller::ExtendedVSCControlParameters = ExtendedVSCControlParameters()
    filter::ExtendedVSCFilterParameters = ExtendedVSCFilterParameters()
    wire_form::ExtendedVSCWireForm = ThreeWireForm
    protection::ExtendedVSCProtectionParameters = ExtendedVSCProtectionParameters()
    scenario::ExtendedVSCScenarioParameters = ExtendedVSCScenarioParameters()
    frequency_hz::Float64 = 50.0
    rated_power_va::Float64 = 50.0e3
    grid_line_line_rms_v::Float64 = 400.0
    grid_source_resistance_ohm::Float64 = 0.08
    grid_load_resistance_ohm::Float64 = 24.0
    neutral_grounding_resistance_ohm::Float64 = 1.0e-3
    dc_source_voltage_v::Float64 = 800.0
    dc_source_resistance_ohm::Float64 = 0.60
    dc_link_capacitance_f::Float64 = 10.0e-3
    current_limit_a::Float64 = 180.0
    carrier_frequency_hz::Float64 = 10.0e3
    scheduler_tick_s::Float64 = 1.0e-6
    timestep_s::Float64 = 1.0e-6
    minimum_duty::Float64 = 0.02
    maximum_duty::Float64 = 0.98
    semiconductor_on_conductance_s::Float64 = 100.0
    semiconductor_off_conductance_s::Float64 = 1.0e-6
    semiconductor_forward_voltage_v::Float64 = 1.2
    diode_on_conductance_s::Float64 = 100.0
    diode_forward_voltage_v::Float64 = 0.9
    snubber_resistance_ohm::Float64 = 20.0
    snubber_capacitance_f::Float64 = 50.0e-9
    commutation_dead_time_s::Float64 = 2.0e-6
    minimum_gate_pulse_width_s::Float64 = 1.0e-6
end

struct ExtendedVSCPlantRequest
    timestamp_s::Float64
    valid_until_s::Float64
    available_active_power_w::Float64
    active_power_reference_w::Float64
    reactive_power_reference_var::Float64
    voltage_reference_v::Float64
    frequency_reference_hz::Float64
    active_power_ramp_w_per_s::Float64
    reactive_power_ramp_var_per_s::Float64
    current_priority::VSCCurrentPriority
    provenance::ParameterProvenance
end

function ExtendedVSCPlantRequest(;
    timestamp_s::Real=0.0,
    valid_until_s::Real=Inf,
    available_active_power_w::Real=50.0e3,
    active_power_reference_w::Real=20.0e3,
    reactive_power_reference_var::Real=0.0,
    voltage_reference_v::Real=sqrt(2.0 / 3.0) * 400.0,
    frequency_reference_hz::Real=50.0,
    active_power_ramp_w_per_s::Real=1.0e9,
    reactive_power_ramp_var_per_s::Real=1.0e9,
    current_priority::VSCCurrentPriority=ActiveCurrentPriority,
    provenance::ParameterProvenance=_generic_vsc_provenance(
        "VSC plant request",
        "watt, var, volt, hertz, second",
    ),
)
    values = Float64[
        timestamp_s,
        available_active_power_w,
        active_power_reference_w,
        reactive_power_reference_var,
        voltage_reference_v,
        frequency_reference_hz,
        active_power_ramp_w_per_s,
        reactive_power_ramp_var_per_s,
    ]
    all(isfinite, values) || throw(ArgumentError(
        "VSC plant request finite fields must be finite",
    ))
    valid_until = Float64(valid_until_s)
    (isfinite(valid_until) || valid_until == Inf) && valid_until >= timestamp_s ||
        throw(ArgumentError("VSC plant request validity must follow its timestamp"))
    available_active_power_w >= 0.0 || throw(ArgumentError(
        "VSC available active power must be nonnegative",
    ))
    voltage_reference_v > 0.0 || throw(ArgumentError(
        "VSC voltage reference must be positive",
    ))
    frequency_reference_hz > 0.0 || throw(ArgumentError(
        "VSC frequency reference must be positive",
    ))
    active_power_ramp_w_per_s >= 0.0 && reactive_power_ramp_var_per_s >= 0.0 ||
        throw(ArgumentError("VSC plant-request ramps must be nonnegative"))
    return ExtendedVSCPlantRequest(
        Float64(timestamp_s),
        valid_until,
        Float64(available_active_power_w),
        Float64(active_power_reference_w),
        Float64(reactive_power_reference_var),
        Float64(voltage_reference_v),
        Float64(frequency_reference_hz),
        Float64(active_power_ramp_w_per_s),
        Float64(reactive_power_ramp_var_per_s),
        current_priority,
        provenance,
    )
end

struct ExtendedVSCMeasurement
    phase_voltage_v::NTuple{3,Float64}
    converter_current_a::NTuple{3,Float64}
    grid_current_a::NTuple{3,Float64}
    filter_voltage_v::NTuple{3,Float64}
    neutral_current_a::Float64
    dc_link_voltage_v::Float64
    time_s::Float64

    function ExtendedVSCMeasurement(
        phase_voltage_v::NTuple{3,<:Real},
        converter_current_a::NTuple{3,<:Real},
        grid_current_a::NTuple{3,<:Real},
        filter_voltage_v::NTuple{3,<:Real},
        neutral_current_a::Real,
        dc_link_voltage_v::Real,
        time_s::Real,
    )
        values = Float64[
            phase_voltage_v...,
            converter_current_a...,
            grid_current_a...,
            filter_voltage_v...,
            neutral_current_a,
            dc_link_voltage_v,
            time_s,
        ]
        all(isfinite, values) || throw(ArgumentError(
            "extended VSC measurements must be finite",
        ))
        dc_link_voltage_v > 0.0 || throw(ArgumentError(
            "extended VSC DC-link measurement must be positive",
        ))
        return new(
            Tuple(Float64.(phase_voltage_v)),
            Tuple(Float64.(converter_current_a)),
            Tuple(Float64.(grid_current_a)),
            Tuple(Float64.(filter_voltage_v)),
            Float64(neutral_current_a),
            Float64(dc_link_voltage_v),
            Float64(time_s),
        )
    end
end

mutable struct VSCResonantAxisState
    derivative_state::Float64
    integral_state::Float64
end

struct VSCRotatingSequenceSample
    positive_direct_v::Float64
    positive_quadrature_v::Float64
    negative_direct_v::Float64
    negative_quadrature_v::Float64
    zero_direct_v::Float64
    zero_quadrature_v::Float64
end

const _ZERO_VSC_ROTATING_SEQUENCE_SAMPLE = VSCRotatingSequenceSample(
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
)

"""Causal one-cycle rotating-frame sequence state retained in every controller checkpoint."""
mutable struct VSCRotatingSequenceState
    samples::Vector{VSCRotatingSequenceSample}
    next_sample_index::Int
    retained_sample_count::Int
    positive_direct_sum_v::Float64
    positive_quadrature_sum_v::Float64
    negative_direct_sum_v::Float64
    negative_quadrature_sum_v::Float64
    zero_direct_sum_v::Float64
    zero_quadrature_sum_v::Float64
end

function VSCRotatingSequenceState(parameters::ExtendedVSCParameters)
    sample_count = max(
        1,
        round(Int, inv(parameters.frequency_hz * parameters.controller.control_period_s)),
    )
    return VSCRotatingSequenceState(
        fill(_ZERO_VSC_ROTATING_SEQUENCE_SAMPLE, sample_count),
        1,
        0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    )
end

mutable struct ExtendedVSCControlState
    angle_rad::Float64
    frequency_rad_per_s::Float64
    pll_integral_rad_per_s::Float64
    pll_locked::Bool
    direct_current_integral_v::Float64
    quadrature_current_integral_v::Float64
    zero_current_integral_v::Float64
    direct_voltage_integral_a::Float64
    quadrature_voltage_integral_a::Float64
    sequence::VSCRotatingSequenceState
    resonant_alpha::Vector{VSCResonantAxisState}
    resonant_beta::Vector{VSCResonantAxisState}
    resonant_zero::Vector{VSCResonantAxisState}
    filtered_active_power_w::Float64
    filtered_reactive_power_var::Float64
    previous_current_alpha_a::Float64
    previous_current_beta_a::Float64
    applied_active_power_reference_w::Float64
    applied_reactive_power_reference_var::Float64
    held_duties::NTuple{4,Float64}
    held_pole_voltage_reference_v::NTuple{4,Float64}
    mode::ExtendedVSCOperatingMode
    request_disposition::VSCPlantRequestDisposition
    sample_count::Int
    saturation_count::Int
    pll_loss_count::Int
    refusal_count::Int
end

function ExtendedVSCControlState(parameters::ExtendedVSCParameters)
    p = validate_extended_vsc_parameters(parameters)
    omega = 2.0 * pi * p.frequency_hz
    resonators() = [VSCResonantAxisState(0.0, 0.0) for _ in _EXTENDED_VSC_HARMONIC_ORDERS]
    return ExtendedVSCControlState(
        0.0,
        omega,
        0.0,
        false,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        VSCRotatingSequenceState(p),
        resonators(),
        resonators(),
        resonators(),
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        (0.5, 0.5, 0.5, 0.5),
        (0.0, 0.0, 0.0, 0.0),
        VSCNormalOperation,
        VSCPlantRequestApplied,
        0,
        0,
        0,
        0,
    )
end

struct ExtendedVSCControlCommand
    controller_family::ExtendedVSCControllerFamily
    wire_form::ExtendedVSCWireForm
    duties::NTuple{4,Float64}
    pole_voltage_reference_v::NTuple{4,Float64}
    phase_current_reference_a::NTuple{3,Float64}
    angle_rad::Float64
    frequency_hz::Float64
    active_power_w::Float64
    reactive_power_var::Float64
    positive_sequence_voltage_v::Float64
    negative_sequence_voltage_v::Float64
    zero_sequence_voltage_v::Float64
    sequence_extractor_settled::Bool
    limited::Bool
    mode::ExtendedVSCOperatingMode
    request_disposition::VSCPlantRequestDisposition
end

const EXTENDED_VSC_CONTRACT = ScientificModelContract(
    :extended_switch_detailed_vsc_platform,
    :typed_controller_filter_wire_matrix;
    owner="AIMORA.SwitchDetailedVSC and AIMORA.EMTStudy",
    maturity=:implemented,
    fidelity=SwitchingDetailed,
    validity_domain=ModelValidityDomain(
        :extended_two_level_vsc_platform;
        description="Four exact GFL/GFM control families with L/LC/LCL and explicit three-/four-wire switch-detailed network forms.",
        bounds=(
            NumericDomainBound(:frequency_hz; unit="Hz", lower=45.0, upper=65.0),
            NumericDomainBound(:rated_power_va; unit="VA", lower=10.0e3, upper=5.0e6),
            NumericDomainBound(:grid_line_line_rms_v; unit="V", lower=208.0, upper=690.0),
            NumericDomainBound(:dc_source_voltage_v; unit="V", lower=500.0, upper=1500.0),
            NumericDomainBound(:carrier_frequency_hz; unit="Hz", lower=2.0e3, upper=50.0e3),
            NumericDomainBound(:timestep_s; unit="s", lower=0.1e-6, upper=5.0e-6),
        ),
        unsupported_phenomena=(
            :three_level_or_multipulse_complete_converter,
            :facts_hvdc_or_mmc_system,
            :renewable_plant_physics,
            :transformer_saturation,
            :adaptive_or_dassl_execution,
            :realtime_hil,
            :vendor_or_grid_code_equivalence,
            :atp_pscad_equivalence,
            :standard_conformance_or_certification,
        ),
    ),
    state_inventory=DynamicStateInventory(
        differential=(
            :dc_link_and_filter_energy,
            :pll_and_current_control,
            :resonant_control,
            :filtered_power_and_grid_forming_oscillator,
            :virtual_impedance_and_voltage_control,
            :selected_device_state,
        ),
        algebraic=(
            :terminal_filter_neutral_voltage_current,
            :frame_sequence_and_power_quantities,
            :limited_references_and_modulation,
            :kcl_charge_energy_and_domain_residuals,
        ),
        discrete=(
            :controller_filter_wire_and_schema_identity,
            :pll_lock_and_current_limit_mode,
            :protection_island_reconnect_and_request_disposition,
            :bridge_gate_and_conduction_state,
        ),
        delayed_history=(
            :measurement_and_control_release,
            :plant_request_and_protection_delay,
            :pwm_gate_and_device_delay,
            :filter_controller_and_sequence_history,
        ),
        scheduler=(
            :measurement_control_dispatch_protection_ticks,
            :carrier_and_event_ticks,
        ),
    ),
    inputs=(
        ContractQuantity(:plant_request; unit="typed SI request"),
        ContractQuantity(:grid_phase_neutral_voltage_v; unit="V"),
        ContractQuantity(:dc_source_voltage_v; unit="V"),
        ContractQuantity(:fault_island_reconnect_schedule; unit="s"),
    ),
    outputs=(
        ContractQuantity(:controller_filter_wire_state; unit="typed state"),
        ContractQuantity(:terminal_filter_neutral_quantities; unit="SI"),
        ContractQuantity(:sequence_harmonic_quantities; unit="SI"),
        ContractQuantity(:plant_and_protection_disposition; unit="typed state"),
        ContractQuantity(:power_loss_and_energy; unit="W and J"),
        ContractQuantity(:adequacy_domain_and_deterministic_signature; unit="typed result"),
    ),
    assumptions=(
        "Every controller/filter/wire combination is explicit and no family inherits another family's evidence.",
        "All filters and neutral paths are finite physical branches in the accepted global nodal solve.",
        "Protection settings are synthetic generic thresholds and make no grid-code or certification claim.",
    ),
    mutation_order=(
        :capture_complete_transaction,
        :apply_protection_fault_island_and_reconnect,
        :acquire_due_plant_request,
        :sample_physical_measurements,
        :update_sequence_pll_power_and_controller_state,
        :project_current_and_modulation_limits,
        :release_pwm_and_apply_bridge_device_events,
        :solve_coupled_nonlinear_network,
        :accept_filter_controller_device_energy_and_output_or_restore,
    ),
)

extended_vsc_contract() = EXTENDED_VSC_CONTRACT

supported_extended_vsc_combination(
    ::ExtendedVSCControllerFamily,
    ::ExtendedVSCFilterFamily,
    ::ExtendedVSCWireForm,
) = true

function _require_finite_positive(value::Real, field::AbstractString)
    isfinite(value) && value > 0.0 || throw(ArgumentError(
        "$field must be finite and positive",
    ))
    return Float64(value)
end

function _require_finite_nonnegative(value::Real, field::AbstractString)
    isfinite(value) && value >= 0.0 || throw(ArgumentError(
        "$field must be finite and nonnegative",
    ))
    return Float64(value)
end

function validate_extended_vsc_parameters(parameters::ExtendedVSCParameters)
    p = parameters
    assert_validity(assess_validity(
        extended_vsc_contract(),
        (
            frequency_hz=p.frequency_hz,
            rated_power_va=p.rated_power_va,
            grid_line_line_rms_v=p.grid_line_line_rms_v,
            dc_source_voltage_v=p.dc_source_voltage_v,
            carrier_frequency_hz=p.carrier_frequency_hz,
            timestep_s=p.timestep_s,
        );
        requested_fidelity=SwitchingDetailed,
    ))
    supported_extended_vsc_combination(
        p.controller.family,
        p.filter.family,
        p.wire_form,
    ) || throw(ArgumentError("unsupported extended VSC controller/filter/wire combination"))
    _validate_vsc_physical_provenance(p.filter.provenance, "VSC filter")
    _validate_vsc_physical_provenance(p.controller.provenance, "VSC controller")
    _validate_vsc_physical_provenance(p.protection.provenance, "VSC protection")
    _validate_vsc_physical_provenance(p.scenario.provenance, "VSC scenario")
    for (value, field) in (
        (p.grid_source_resistance_ohm, "grid-source resistance"),
        (p.grid_load_resistance_ohm, "grid-load resistance"),
        (p.neutral_grounding_resistance_ohm, "neutral-grounding resistance"),
        (p.dc_source_resistance_ohm, "DC-source resistance"),
        (p.dc_link_capacitance_f, "DC-link capacitance"),
        (p.current_limit_a, "current limit"),
        (p.snubber_resistance_ohm, "snubber resistance"),
        (p.snubber_capacitance_f, "snubber capacitance"),
        (p.scheduler_tick_s, "scheduler tick"),
        (p.filter.converter_inductance_h, "converter-side inductance"),
        (p.filter.converter_resistance_ohm, "converter-side resistance"),
        (p.controller.control_period_s, "controller period"),
        (p.controller.power_filter_time_constant_s, "power-filter time constant"),
    )
        _require_finite_positive(value, field)
    end
    if p.filter.family !== SeriesLFilter
        _require_finite_positive(p.filter.shunt_capacitance_f, "shunt capacitance")
        _require_finite_positive(
            p.filter.shunt_damping_conductance_s,
            "shunt damping conductance",
        )
    end
    if p.filter.family === LCLFilter
        _require_finite_positive(p.filter.grid_resistance_ohm, "grid-side resistance")
        _require_finite_positive(p.filter.grid_inductance_h, "grid-side inductance")
        resonance_rad_per_s = sqrt(
            (p.filter.converter_inductance_h + p.filter.grid_inductance_h) /
            (p.filter.converter_inductance_h * p.filter.grid_inductance_h *
             p.filter.shunt_capacitance_f),
        )
        resonance_rad_per_s < pi / p.controller.control_period_s || throw(ArgumentError(
            "LCL resonance must remain below the sampled-controller Nyquist frequency",
        ))
        resonance_rad_per_s < pi / p.timestep_s || throw(ArgumentError(
            "LCL resonance must remain below the electrical-step Nyquist frequency",
        ))
    end
    if p.wire_form === FourWireForm
        _require_finite_positive(p.filter.neutral_resistance_ohm, "neutral resistance")
        _require_finite_nonnegative(p.filter.neutral_inductance_h, "neutral inductance")
    elseif p.scenario.zero_sequence_voltage_ratio != 0.0
        throw(ArgumentError("three-wire VSC refuses a nonzero zero-sequence source"))
    end
    for (value, field) in (
        (p.controller.pll_proportional_gain_rad_per_s, "PLL proportional gain"),
        (p.controller.pll_integral_gain_rad_per_s2, "PLL integral gain"),
        (p.controller.pll_voltage_floor_v, "PLL voltage floor"),
        (p.controller.current_proportional_gain_v_per_a, "current proportional gain"),
        (p.controller.current_integral_gain_v_per_as, "current integral gain"),
        (p.controller.resonant_bandwidth_rad_per_s, "resonant bandwidth"),
        (p.controller.virtual_inertia_w_s2_per_rad, "virtual inertia"),
        (p.controller.virtual_current_filter_time_constant_s, "virtual-current filter time constant"),
        (p.protection.ac_overcurrent_a, "protection AC overcurrent threshold"),
        (p.protection.dc_undervoltage_v, "protection DC undervoltage threshold"),
        (p.protection.dc_overvoltage_v, "protection DC overvoltage threshold"),
        (p.protection.minimum_phase_rms_voltage_v, "protection minimum phase RMS voltage"),
        (p.protection.maximum_phase_rms_voltage_v, "protection maximum phase RMS voltage"),
        (p.protection.trip_delay_s, "protection trip delay"),
        (p.protection.restart_delay_s, "protection restart delay"),
    )
        _require_finite_positive(value, field)
    end
    p.controller.pll_minimum_frequency_hz < p.controller.pll_maximum_frequency_hz ||
        throw(ArgumentError("PLL frequency bounds must be strictly increasing"))
    p.protection.minimum_frequency_hz < p.protection.maximum_frequency_hz ||
        throw(ArgumentError("protection frequency bounds must be strictly increasing"))
    p.protection.minimum_phase_rms_voltage_v <
        p.protection.maximum_phase_rms_voltage_v || throw(ArgumentError(
        "protection voltage bounds must be strictly increasing",
    ))
    p.protection.dc_undervoltage_v < p.protection.dc_overvoltage_v ||
        throw(ArgumentError("protection DC voltage bounds must be strictly increasing"))
    for (value, field) in (
        (p.controller.current_antiwindup_gain_per_s, "current anti-windup gain"),
        (p.controller.resonant_gain_v_per_a, "resonant gain"),
        (p.controller.active_power_frequency_droop_rad_per_ws, "active-power droop"),
        (p.controller.reactive_power_voltage_droop_v_per_var, "reactive-power droop"),
        (p.controller.virtual_damping_w_s_per_rad, "virtual damping"),
        (p.controller.virtual_resistance_ohm, "virtual resistance"),
        (p.controller.virtual_inductance_h, "virtual inductance"),
        (p.controller.control_delay_s, "control delay"),
    )
        _require_finite_nonnegative(value, field)
    end
    all(order -> order in _EXTENDED_VSC_HARMONIC_ORDERS,
        p.controller.selected_harmonic_orders) || throw(ArgumentError(
        "VSC PR harmonic orders must be selected from 1, 5, 7, 11, and 13",
    ))
    length(unique(p.controller.selected_harmonic_orders)) ==
        length(p.controller.selected_harmonic_orders) || throw(ArgumentError(
        "VSC PR harmonic orders must be unique",
    ))
    length(p.scenario.harmonic_orders) == length(p.scenario.harmonic_voltage_ratios) ||
        throw(DimensionMismatch("VSC source harmonic orders and ratios must match"))
    all(order -> 2 <= order <= 25, p.scenario.harmonic_orders) || throw(ArgumentError(
        "VSC source harmonic orders must be between 2 and 25",
    ))
    0.0 <= p.scenario.negative_sequence_voltage_ratio <= 0.2 || throw(ArgumentError(
        "negative-sequence voltage ratio must be within [0, 0.2]",
    ))
    0.0 <= p.scenario.zero_sequence_voltage_ratio <= 0.2 || throw(ArgumentError(
        "zero-sequence voltage ratio must be within [0, 0.2]",
    ))
    0.0 <= p.scenario.fault_voltage_factor <= 1.0 || throw(ArgumentError(
        "VSC fault voltage factor must lie within [0, 1]",
    ))
    0.0 < p.scenario.dc_sag_factor <= 1.0 || throw(ArgumentError(
        "VSC DC-sag factor must lie within (0, 1]",
    ))
    all(ratio -> isfinite(ratio) && 0.0 <= ratio <= 0.1,
        p.scenario.harmonic_voltage_ratios) || throw(ArgumentError(
        "source harmonic voltage ratios must be finite and within [0, 0.1]",
    ))
    0.0 < p.minimum_duty < p.maximum_duty < 1.0 || throw(ArgumentError(
        "extended VSC duty bounds must lie strictly inside [0, 1]",
    ))
    carrier_period_s = inv(p.carrier_frequency_hz)
    carrier_period_s / p.timestep_s >= 40.0 || throw(ArgumentError(
        "extended VSC timestep must provide at least 40 points per carrier period",
    ))
    for (value, field) in (
        (p.timestep_s, "timestep"),
        (p.controller.control_period_s, "control period"),
        (p.controller.control_delay_s, "control delay"),
        (carrier_period_s, "carrier period"),
        (p.commutation_dead_time_s, "commutation dead time"),
        (p.minimum_gate_pulse_width_s, "minimum gate pulse width"),
    )
        _is_integer_multiple(value, p.scheduler_tick_s) || throw(ArgumentError(
            "$field must be an integer multiple of scheduler_tick_s",
        ))
    end
    scenario = p.scenario
    0.0 <= scenario.fault_start_s < scenario.fault_end_s <= scenario.end_time_s ||
        throw(ArgumentError("VSC fault interval must be ordered inside the horizon"))
    0.0 <= scenario.dc_sag_start_s < scenario.dc_sag_end_s <= scenario.end_time_s ||
        throw(ArgumentError("VSC DC-sag interval must be ordered inside the horizon"))
    0.0 <= scenario.block_time_s < scenario.restart_time_s <= scenario.end_time_s ||
        throw(ArgumentError("VSC block/restart interval must be ordered inside the horizon"))
    0.0 <= scenario.island_start_s < scenario.reconnect_time_s <= scenario.end_time_s ||
        throw(ArgumentError("VSC island/reconnect interval must be ordered inside the horizon"))
    scenario.fault_kind in (:phase_to_ground, :phase_to_phase, :three_phase) ||
        throw(ArgumentError("unsupported VSC fault kind"))
    all(phase -> 1 <= phase <= 3, scenario.faulted_phases) || throw(ArgumentError(
        "VSC faulted phase indices must lie in 1:3",
    ))
    return p
end

function _vsc_rotating_sequence_sample(
    voltage::StationaryReferenceFrame,
    angle_rad::Float64,
)
    positive = park_transform(voltage, angle_rad)
    negative = park_transform(voltage, -angle_rad)
    cosine = cos(angle_rad)
    sine = sin(angle_rad)
    return VSCRotatingSequenceSample(
        positive.direct,
        positive.quadrature,
        negative.direct,
        negative.quadrature,
        2.0 * voltage.zero * cosine,
        -2.0 * voltage.zero * sine,
    )
end

"""
    advance_vsc_sequence_extractor!(state, phase_voltage_v, angle_rad)

Advance a causal rotating-frame, one-cycle sequence extractor. The positive and
negative components are demodulated in counter-rotating frames; the scalar zero
component is synchronously demodulated into an orthogonal pair. The retained
window rejects the opposite fundamental sequence and integer harmonics after one
nominal cycle. All delay history and running sums live in `state`, so transaction
rollback and restart can restore the extractor without hidden signal history.
"""
function advance_vsc_sequence_extractor!(
    state::VSCRotatingSequenceState,
    phase_voltage_v::NTuple{3,<:Real},
    angle_rad::Real,
)
    angle = Float64(angle_rad)
    isfinite(angle) || throw(ArgumentError("VSC sequence angle must be finite"))
    sample = _vsc_rotating_sequence_sample(clarke_transform(phase_voltage_v), angle)
    index = state.next_sample_index
    previous = state.samples[index]
    window_full = state.retained_sample_count == length(state.samples)
    if window_full
        state.positive_direct_sum_v -= previous.positive_direct_v
        state.positive_quadrature_sum_v -= previous.positive_quadrature_v
        state.negative_direct_sum_v -= previous.negative_direct_v
        state.negative_quadrature_sum_v -= previous.negative_quadrature_v
        state.zero_direct_sum_v -= previous.zero_direct_v
        state.zero_quadrature_sum_v -= previous.zero_quadrature_v
    else
        state.retained_sample_count += 1
    end
    state.samples[index] = sample
    state.positive_direct_sum_v += sample.positive_direct_v
    state.positive_quadrature_sum_v += sample.positive_quadrature_v
    state.negative_direct_sum_v += sample.negative_direct_v
    state.negative_quadrature_sum_v += sample.negative_quadrature_v
    state.zero_direct_sum_v += sample.zero_direct_v
    state.zero_quadrature_sum_v += sample.zero_quadrature_v
    state.next_sample_index = index == length(state.samples) ? 1 : index + 1
    divisor = Float64(state.retained_sample_count)
    positive_direct = state.positive_direct_sum_v / divisor
    positive_quadrature = state.positive_quadrature_sum_v / divisor
    negative_direct = state.negative_direct_sum_v / divisor
    negative_quadrature = state.negative_quadrature_sum_v / divisor
    zero_direct = state.zero_direct_sum_v / divisor
    zero_quadrature = state.zero_quadrature_sum_v / divisor
    return (
        positive=SynchronousReferenceFrame(positive_direct, positive_quadrature, 0.0),
        negative=SynchronousReferenceFrame(negative_direct, negative_quadrature, 0.0),
        zero=SynchronousReferenceFrame(zero_direct, zero_quadrature, 0.0),
        positive_magnitude_v=hypot(positive_direct, positive_quadrature),
        negative_magnitude_v=hypot(negative_direct, negative_quadrature),
        zero_magnitude_v=hypot(zero_direct, zero_quadrature),
        settled=state.retained_sample_count == length(state.samples),
    )
end

function project_vsc_current_reference(
    direct_a::Real,
    quadrature_a::Real,
    limit_a::Real,
    priority::VSCCurrentPriority,
)
    direct = Float64(direct_a)
    quadrature = Float64(quadrature_a)
    limit = _require_finite_positive(limit_a, "VSC current projection limit")
    all(isfinite, (direct, quadrature)) || throw(ArgumentError(
        "VSC current reference must be finite",
    ))
    magnitude = hypot(direct, quadrature)
    magnitude <= limit && return (direct=direct, quadrature=quadrature, limited=false)
    if priority === ActiveCurrentPriority
        direct_limited = clamp(direct, -limit, limit)
        quadrature_limited = clamp(
            quadrature,
            -sqrt(max(0.0, limit^2 - direct_limited^2)),
            sqrt(max(0.0, limit^2 - direct_limited^2)),
        )
    elseif priority === ReactiveCurrentPriority
        quadrature_limited = clamp(quadrature, -limit, limit)
        direct_limited = clamp(
            direct,
            -sqrt(max(0.0, limit^2 - quadrature_limited^2)),
            sqrt(max(0.0, limit^2 - quadrature_limited^2)),
        )
    else
        scale = limit / magnitude
        direct_limited = scale * direct
        quadrature_limited = scale * quadrature
    end
    return (direct=direct_limited, quadrature=quadrature_limited, limited=true)
end

function extended_vsc_modulation_duties(
    phase_reference_v::NTuple{3,<:Real},
    neutral_reference_v::Real,
    dc_link_voltage_v::Real,
    wire_form::ExtendedVSCWireForm;
    minimum_duty::Real=0.02,
    maximum_duty::Real=0.98,
)
    dc = _require_finite_positive(dc_link_voltage_v, "extended VSC DC-link voltage")
    minimum = Float64(minimum_duty)
    maximum = Float64(maximum_duty)
    0.0 <= minimum < maximum <= 1.0 || throw(ArgumentError(
        "extended VSC duty bounds must satisfy 0 <= minimum < maximum <= 1",
    ))
    phase = Float64.(phase_reference_v)
    all(isfinite, phase) || throw(ArgumentError("phase references must be finite"))
    if wire_form === ThreeWireForm
        duties = modulation_duties(
            Tuple(phase),
            dc,
            ZeroSequenceInjectedPulseWidthModulation;
            minimum_duty=minimum,
            maximum_duty=maximum,
        )
        return (duties..., 0.5)
    end
    neutral = Float64(neutral_reference_v)
    isfinite(neutral) || throw(ArgumentError("neutral reference must be finite"))
    pole = (phase[1], phase[2], phase[3], neutral)
    minimum_pole, maximum_pole = extrema(pole)
    common = -0.5 * (maximum_pole + minimum_pole)
    return ntuple(4) do index
        clamp(0.5 + (pole[index] + common) / dc, minimum, maximum)
    end
end

function _advance_pll!(
    state::ExtendedVSCControlState,
    positive_sequence::SynchronousReferenceFrame,
    sequence_extractor_settled::Bool,
    parameters::ExtendedVSCParameters,
)
    control = parameters.controller
    dt = control.control_period_s
    if !sequence_extractor_settled
        state.angle_rad = rem2pi(
            state.angle_rad + dt * state.frequency_rad_per_s,
            RoundNearest,
        )
        state.pll_locked = false
        return state.angle_rad, state.frequency_rad_per_s
    end
    magnitude = hypot(positive_sequence.direct, positive_sequence.quadrature)
    if magnitude < control.pll_voltage_floor_v
        state.pll_locked = false
        state.pll_loss_count += 1
        control.pll_loss_policy === BlockVSCOnVoltageLoss &&
            (state.mode = VSCBlockedOperation)
        return state.angle_rad, state.frequency_rad_per_s
    end
    error = positive_sequence.quadrature /
        max(magnitude, control.pll_voltage_floor_v)
    state.pll_integral_rad_per_s += control.pll_integral_gain_rad_per_s2 * dt * error
    minimum_omega = 2.0 * pi * control.pll_minimum_frequency_hz
    maximum_omega = 2.0 * pi * control.pll_maximum_frequency_hz
    state.frequency_rad_per_s = clamp(
        2.0 * pi * parameters.frequency_hz +
        control.pll_proportional_gain_rad_per_s * error +
        state.pll_integral_rad_per_s,
        minimum_omega,
        maximum_omega,
    )
    state.angle_rad = rem2pi(
        state.angle_rad + dt * state.frequency_rad_per_s,
        RoundNearest,
    )
    state.pll_locked = abs(error) <= 0.05
    return state.angle_rad, state.frequency_rad_per_s
end

function _filtered_power_update(previous, measured, time_constant, step)
    decay = exp(-step / time_constant)
    return decay * previous + (1.0 - decay) * measured
end

function _resonant_derivative(state, error, frequency_rad_per_s, order, control)
    return (
        -2.0 * control.resonant_bandwidth_rad_per_s * state.derivative_state -
        (order * frequency_rad_per_s)^2 * state.integral_state +
        2.0 * control.resonant_gain_v_per_a *
        control.resonant_bandwidth_rad_per_s * error,
        state.derivative_state,
    )
end

function _advance_resonator!(state, error, frequency_rad_per_s, order, control)
    dt = control.control_period_s
    x1 = state.derivative_state
    x2 = state.integral_state
    derivative(a, b) = _resonant_derivative(
        VSCResonantAxisState(a, b),
        error,
        frequency_rad_per_s,
        order,
        control,
    )
    k1 = derivative(x1, x2)
    k2 = derivative(x1 + 0.5 * dt * k1[1], x2 + 0.5 * dt * k1[2])
    k3 = derivative(x1 + 0.5 * dt * k2[1], x2 + 0.5 * dt * k2[2])
    k4 = derivative(x1 + dt * k3[1], x2 + dt * k3[2])
    state.derivative_state = x1 + (dt / 6.0) *
        (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
    state.integral_state = x2 + (dt / 6.0) *
        (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2])
    return state.derivative_state
end

function _validated_request(
    state::ExtendedVSCControlState,
    request::ExtendedVSCPlantRequest,
    time_s::Float64,
    control_period_s::Float64,
)
    _validate_vsc_physical_provenance(request.provenance, "VSC plant request")
    if time_s < request.timestamp_s || time_s > request.valid_until_s
        state.request_disposition = VSCPlantRequestStale
        state.refusal_count += 1
        return nothing
    end
    active_target = clamp(
        request.active_power_reference_w,
        -request.available_active_power_w,
        request.available_active_power_w,
    )
    active_delta = clamp(
        active_target - state.applied_active_power_reference_w,
        -request.active_power_ramp_w_per_s * control_period_s,
        request.active_power_ramp_w_per_s * control_period_s,
    )
    reactive_delta = clamp(
        request.reactive_power_reference_var - state.applied_reactive_power_reference_var,
        -request.reactive_power_ramp_var_per_s * control_period_s,
        request.reactive_power_ramp_var_per_s * control_period_s,
    )
    state.applied_active_power_reference_w += active_delta
    state.applied_reactive_power_reference_var += reactive_delta
    limited = active_target != request.active_power_reference_w ||
        active_delta != active_target - (state.applied_active_power_reference_w - active_delta) ||
        reactive_delta != request.reactive_power_reference_var -
            (state.applied_reactive_power_reference_var - reactive_delta)
    state.request_disposition = limited ? VSCPlantRequestLimited : VSCPlantRequestApplied
    return (
        active_power_w=state.applied_active_power_reference_w,
        reactive_power_var=state.applied_reactive_power_reference_var,
        voltage_v=request.voltage_reference_v,
        frequency_hz=request.frequency_reference_hz,
        priority=request.current_priority,
    )
end

function _current_reference_from_power(active_power, reactive_power, voltage_dq, limit, priority)
    voltage_squared = max(voltage_dq.direct^2 + voltage_dq.quadrature^2, 1.0)
    direct = (2.0 / 3.0) * (
        voltage_dq.direct * active_power + voltage_dq.quadrature * reactive_power
    ) / voltage_squared
    quadrature = (2.0 / 3.0) * (
        voltage_dq.quadrature * active_power - voltage_dq.direct * reactive_power
    ) / voltage_squared
    return project_vsc_current_reference(direct, quadrature, limit, priority)
end

function _apply_voltage_limit(
    direct,
    quadrature,
    zero,
    dc_voltage,
    state::ExtendedVSCControlState,
    control::ExtendedVSCControlParameters,
)
    voltage_limit = 0.95 * dc_voltage / sqrt(3.0)
    magnitude = sqrt(direct^2 + quadrature^2 + zero^2)
    if magnitude <= voltage_limit || magnitude == 0.0
        return direct, quadrature, zero, false
    end
    scale = voltage_limit / magnitude
    state.saturation_count += 1
    return scale * direct, scale * quadrature, scale * zero, true
end

function compute_extended_vsc_control!(
    state::ExtendedVSCControlState,
    measurement::ExtendedVSCMeasurement,
    request::ExtendedVSCPlantRequest,
    parameters::ExtendedVSCParameters,
)
    p = validate_extended_vsc_parameters(parameters)
    control = p.controller
    dt = control.control_period_s
    validated_request = _validated_request(state, request, measurement.time_s, dt)
    if validated_request === nothing
        state.mode = VSCBlockedOperation
        state.sample_count += 1
        return ExtendedVSCControlCommand(
            control.family,
            p.wire_form,
            (0.5, 0.5, 0.5, 0.5),
            (0.0, 0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
            state.angle_rad,
            state.frequency_rad_per_s / (2.0 * pi),
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            false,
            false,
            state.mode,
            state.request_disposition,
        )
    end
    voltage_stationary = clarke_transform(measurement.phase_voltage_v)
    current_stationary = clarke_transform(measurement.grid_current_a)
    sequence = advance_vsc_sequence_extractor!(
        state.sequence,
        measurement.phase_voltage_v,
        state.angle_rad,
    )
    power = instantaneous_three_phase_power(
        measurement.phase_voltage_v,
        measurement.grid_current_a,
        state.angle_rad,
    )
    state.filtered_active_power_w = _filtered_power_update(
        state.filtered_active_power_w,
        power.active_w,
        control.power_filter_time_constant_s,
        dt,
    )
    state.filtered_reactive_power_var = _filtered_power_update(
        state.filtered_reactive_power_var,
        power.reactive_var,
        control.power_filter_time_constant_s,
        dt,
    )
    if control.family in (SynchronousPLLGridFollowing, StationaryResonantGridFollowing)
        _advance_pll!(state, sequence.positive, sequence.settled, p)
    elseif control.family === PowerDroopGridForming
        state.frequency_rad_per_s = 2.0 * pi * validated_request.frequency_hz -
            control.active_power_frequency_droop_rad_per_ws *
            (state.filtered_active_power_w - validated_request.active_power_w)
        state.angle_rad = rem2pi(
            state.angle_rad + dt * state.frequency_rad_per_s,
            RoundNearest,
        )
        state.pll_locked = false
    else
        acceleration = (
            validated_request.active_power_w - state.filtered_active_power_w -
            control.virtual_damping_w_s_per_rad *
            (state.frequency_rad_per_s - 2.0 * pi * validated_request.frequency_hz)
        ) / control.virtual_inertia_w_s2_per_rad
        state.frequency_rad_per_s += dt * acceleration
        state.angle_rad = rem2pi(
            state.angle_rad + dt * state.frequency_rad_per_s,
            RoundNearest,
        )
        state.pll_locked = false
    end
    voltage_dq = if control.family in (
        SynchronousPLLGridFollowing,
        StationaryResonantGridFollowing,
    ) && sequence.settled
        sequence.positive
    else
        park_transform(voltage_stationary, state.angle_rad)
    end
    current_dq = park_transform(current_stationary, state.angle_rad)
    virtual_current_decay = exp(-dt / control.virtual_current_filter_time_constant_s)
    filtered_current_alpha = virtual_current_decay * state.previous_current_alpha_a +
        (1.0 - virtual_current_decay) * current_stationary.alpha
    filtered_current_beta = virtual_current_decay * state.previous_current_beta_a +
        (1.0 - virtual_current_decay) * current_stationary.beta
    filtered_current_derivative = StationaryReferenceFrame(
        (filtered_current_alpha - state.previous_current_alpha_a) / dt,
        (filtered_current_beta - state.previous_current_beta_a) / dt,
        0.0,
    )
    state.previous_current_alpha_a = filtered_current_alpha
    state.previous_current_beta_a = filtered_current_beta
    current_reference = if control.family in (
        SynchronousPLLGridFollowing,
        StationaryResonantGridFollowing,
    )
        _current_reference_from_power(
            validated_request.active_power_w,
            validated_request.reactive_power_var,
            voltage_dq,
            p.current_limit_a,
            validated_request.priority,
        )
    else
        voltage_magnitude_reference = validated_request.voltage_v -
            control.reactive_power_voltage_droop_v_per_var *
            (state.filtered_reactive_power_var - validated_request.reactive_power_var)
        virtual_current = park_transform(
            StationaryReferenceFrame(filtered_current_alpha, filtered_current_beta, 0.0),
            state.angle_rad,
        )
        virtual_current_rate = park_transform(
            filtered_current_derivative,
            state.angle_rad,
        )
        virtual_drop_direct = control.virtual_resistance_ohm * virtual_current.direct +
            control.virtual_inductance_h * virtual_current_rate.direct
        virtual_drop_quadrature = control.virtual_resistance_ohm *
            virtual_current.quadrature + control.virtual_inductance_h *
            virtual_current_rate.quadrature
        voltage_error_direct = voltage_magnitude_reference - virtual_drop_direct -
            voltage_dq.direct
        voltage_error_quadrature = -virtual_drop_quadrature - voltage_dq.quadrature
        state.direct_voltage_integral_a +=
            dt * control.voltage_integral_gain_a_per_vs * voltage_error_direct
        state.quadrature_voltage_integral_a +=
            dt * control.voltage_integral_gain_a_per_vs * voltage_error_quadrature
        project_vsc_current_reference(
            control.voltage_proportional_gain_a_per_v * voltage_error_direct +
                state.direct_voltage_integral_a,
            control.voltage_proportional_gain_a_per_v * voltage_error_quadrature +
                state.quadrature_voltage_integral_a,
            p.current_limit_a,
            validated_request.priority,
        )
    end
    current_reference.limited && begin
        state.mode = VSCCurrentLimitedOperation
        state.saturation_count += 1
    end
    direct_error = current_reference.direct - current_dq.direct
    quadrature_error = current_reference.quadrature - current_dq.quadrature
    zero_reference = p.wire_form === FourWireForm ?
        clamp(-voltage_stationary.zero * control.voltage_proportional_gain_a_per_v,
            -0.25 * p.current_limit_a, 0.25 * p.current_limit_a) : 0.0
    zero_error = zero_reference - current_stationary.zero
    omega = state.frequency_rad_per_s
    resistance = p.filter.converter_resistance_ohm
    inductance = p.filter.converter_inductance_h
    direct_voltage = 0.0
    quadrature_voltage = 0.0
    zero_voltage = 0.0
    if control.family === StationaryResonantGridFollowing
        current_reference_stationary = inverse_park_transform(
            SynchronousReferenceFrame(
                current_reference.direct,
                current_reference.quadrature,
                zero_reference,
            ),
            state.angle_rad,
        )
        alpha_error = current_reference_stationary.alpha - current_stationary.alpha
        beta_error = current_reference_stationary.beta - current_stationary.beta
        alpha_resonant = 0.0
        beta_resonant = 0.0
        zero_resonant = 0.0
        for (index, order) in enumerate(_EXTENDED_VSC_HARMONIC_ORDERS)
            order in control.selected_harmonic_orders || continue
            alpha_resonant += _advance_resonator!(
                state.resonant_alpha[index], alpha_error, omega, order, control,
            )
            beta_resonant += _advance_resonator!(
                state.resonant_beta[index], beta_error, omega, order, control,
            )
            p.wire_form === FourWireForm && (zero_resonant += _advance_resonator!(
                state.resonant_zero[index], zero_error, omega, order, control,
            ))
        end
        stationary_voltage = StationaryReferenceFrame(
            voltage_stationary.alpha + resistance * current_stationary.alpha +
                control.current_proportional_gain_v_per_a * alpha_error + alpha_resonant,
            voltage_stationary.beta + resistance * current_stationary.beta +
                control.current_proportional_gain_v_per_a * beta_error + beta_resonant,
            p.wire_form === FourWireForm ?
                voltage_stationary.zero + resistance * current_stationary.zero +
                control.current_proportional_gain_v_per_a * zero_error + zero_resonant : 0.0,
        )
        synchronous_voltage = park_transform(stationary_voltage, state.angle_rad)
        direct_voltage = synchronous_voltage.direct
        quadrature_voltage = synchronous_voltage.quadrature
        zero_voltage = synchronous_voltage.zero
    else
        direct_integral_candidate = state.direct_current_integral_v +
            dt * control.current_integral_gain_v_per_as * direct_error
        quadrature_integral_candidate = state.quadrature_current_integral_v +
            dt * control.current_integral_gain_v_per_as * quadrature_error
        zero_integral_candidate = state.zero_current_integral_v +
            dt * control.current_integral_gain_v_per_as * zero_error
        direct_voltage = voltage_dq.direct + resistance * current_dq.direct -
            omega * inductance * current_dq.quadrature +
            control.current_proportional_gain_v_per_a * direct_error +
            direct_integral_candidate
        quadrature_voltage = voltage_dq.quadrature + resistance * current_dq.quadrature +
            omega * inductance * current_dq.direct +
            control.current_proportional_gain_v_per_a * quadrature_error +
            quadrature_integral_candidate
        zero_voltage = p.wire_form === FourWireForm ?
            voltage_stationary.zero + resistance * current_stationary.zero +
            control.current_proportional_gain_v_per_a * zero_error +
            zero_integral_candidate : 0.0
        limited_direct, limited_quadrature, limited_zero, voltage_limited =
            _apply_voltage_limit(
                direct_voltage,
                quadrature_voltage,
                zero_voltage,
                measurement.dc_link_voltage_v,
                state,
                control,
            )
        if voltage_limited
            state.direct_current_integral_v = direct_integral_candidate + dt *
                control.current_antiwindup_gain_per_s *
                (limited_direct - direct_voltage)
            state.quadrature_current_integral_v = quadrature_integral_candidate + dt *
                control.current_antiwindup_gain_per_s *
                (limited_quadrature - quadrature_voltage)
            state.zero_current_integral_v = zero_integral_candidate + dt *
                control.current_antiwindup_gain_per_s *
                (limited_zero - zero_voltage)
            direct_voltage = limited_direct
            quadrature_voltage = limited_quadrature
            zero_voltage = limited_zero
            state.mode = VSCCurrentLimitedOperation
        else
            state.direct_current_integral_v = direct_integral_candidate
            state.quadrature_current_integral_v = quadrature_integral_candidate
            state.zero_current_integral_v = zero_integral_candidate
        end
    end
    phase_voltage_reference = inverse_clarke_transform(inverse_park_transform(
        SynchronousReferenceFrame(direct_voltage, quadrature_voltage, zero_voltage),
        state.angle_rad,
    ))
    neutral_reference = p.wire_form === FourWireForm ? -zero_voltage : 0.0
    duties = extended_vsc_modulation_duties(
        phase_voltage_reference,
        neutral_reference,
        measurement.dc_link_voltage_v,
        p.wire_form;
        minimum_duty=p.minimum_duty,
        maximum_duty=p.maximum_duty,
    )
    pole_reference = (
        phase_voltage_reference[1],
        phase_voltage_reference[2],
        phase_voltage_reference[3],
        neutral_reference,
    )
    state.held_duties = duties
    state.held_pole_voltage_reference_v = pole_reference
    state.sample_count += 1
    if state.mode !== VSCBlockedOperation && !current_reference.limited
        state.mode = VSCNormalOperation
    end
    phase_current_reference = inverse_clarke_transform(inverse_park_transform(
        SynchronousReferenceFrame(
            current_reference.direct,
            current_reference.quadrature,
            zero_reference,
        ),
        state.angle_rad,
    ))
    return ExtendedVSCControlCommand(
        control.family,
        p.wire_form,
        duties,
        pole_reference,
        phase_current_reference,
        state.angle_rad,
        state.frequency_rad_per_s / (2.0 * pi),
        power.active_w,
        power.reactive_var,
        sequence.positive_magnitude_v,
        sequence.negative_magnitude_v,
        sequence.zero_magnitude_v,
        sequence.settled,
        state.mode === VSCCurrentLimitedOperation,
        state.mode,
        state.request_disposition,
    )
end

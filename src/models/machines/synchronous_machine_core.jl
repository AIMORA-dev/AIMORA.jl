
using LinearAlgebra
using Printf

using ..Companion: PASTMachineHistoryState,
                   past_machine_history_state,
                   past_machine_history_update!

export SynchronousMachineTACSInterfaceState,
       SynchronousMachineRotorMassTACSState,
       SynchronousMachineEquationState,
       SynchronousMachineExciterPortState,
       SynchronousMachineDynamicState,
       SynchronousMachineRotorMass,
       SynchronousMachineParameters,
       SynchronousMachineInitialization,
       UniversalMachineRotorCurrentState,
       UniversalMachineStatorExcitationCurrentState,
       UniversalMachineDirectCouplingState,
       UniversalMachineHistoryCurrentState,
       UniversalMachineMechanicalIterationState,
       UniversalMachineFluxSaturationState,
       UniversalMachineFluxCurrentJacobianState,
       UniversalMachineSharedTorquePredictionState,
       UniversalMachinePostsolveState,
       UniversalMachineReportScheduleState,
       InductionMachineParameters,
       InductionMachineState,
       MachineTerminalPredictionState,
       MachineConvergenceError,
       CoupledDQMachineParameters,
       CoupledDQMachineState,
       DirectCurrentMachineParameters,
       DirectCurrentMachineState,
       UniversalMachineType4Parameters,
       UniversalMachineType4State,
       InductionMachineSteadyStateEquivalent,
       MachineTerminalCurrentState,
       MachineNetworkCouplingState,
       machine_terminal_current_state_preview,
       machine_terminal_current_state_update!,
       machine_network_coupling_state_preview,
       machine_network_coupling_state_update!,
       synchronous_machine_equation_step_preview,
       synchronous_machine_equation_step!,
       synchronous_machine_exciter_port_update!,
       synchronous_machine_companion_update!,
       synchronous_machine_dynamic_step!,
       synchronous_machine_initial_state,
       universal_machine_rotor_current_solution_preview,
       universal_machine_rotor_current_solution!,
       universal_machine_stator_excitation_current_solution_preview,
       universal_machine_stator_excitation_current_solution!,
       universal_machine_direct_coupling_preview,
       universal_machine_direct_coupling!,
       universal_machine_history_current_update_preview,
       universal_machine_history_current_update!,
       universal_machine_mechanical_iteration_preview,
       universal_machine_mechanical_iteration!,
       universal_machine_flux_saturation_preview,
       universal_machine_flux_saturation!,
       universal_machine_flux_current_jacobian_preview,
       universal_machine_flux_current_jacobian!,
       universal_machine_shared_torque_prediction_preview,
       universal_machine_shared_torque_prediction!,
       universal_machine_postsolve_update_preview,
       universal_machine_postsolve_update!,
       universal_machine_report_schedule_preview,
       universal_machine_report_schedule_update!,
       induction_machine_axis_fluxes,
       induction_machine_step!,
       coupled_dq_machine_step!,
       coupled_dq_power_transform,
       predict_machine_terminal_currents!,
       universal_machine_type4_step!,
       induction_machine_steady_state_equivalent,
       induction_machine_initial_state,
       single_phase_induction_steady_state_initialization,
       synchronous_machine_tacs_transfer_preview,
       synchronous_machine_tacs_transfer_update!,
       synchronous_machine_tacs_transfer_from_past_preview,
       synchronous_machine_tacs_transfer_from_past_update!,
       synchronous_machine_rotor_mass_tacs_preview,
       synchronous_machine_rotor_mass_tacs_update!

Base.@kwdef mutable struct SynchronousMachineTACSInterfaceState
    etac_values::Vector{Float64} = Float64[]
    lmset::Int = 0
    transfer_count::Int = 0
    transfer_mutated::Bool = false
end

Base.@kwdef mutable struct SynchronousMachineRotorMassTACSState
    machine_output_table::Vector{Float64} = Float64[]
    tacs_values::Vector{Float64} = Float64[]
    transfer_pass_marker::Int = 0
    mass_history_initialized::Bool = false
    transfer_mutated::Bool = false
end

Base.@kwdef mutable struct SynchronousMachineEquationState
    cu_values::Vector{Float64} = Float64[]
    histq_values::Vector{Float64} = Float64[]
    shp_values::Vector{Float64} = Float64[]
    iteration_count::Int = 0
    converged::Bool = false
    equation_mutated::Bool = false
end

Base.@kwdef mutable struct SynchronousMachineExciterPortState
    field_winding_current_a::Float64 = 0.0
    current_reduction_factor::Float64 = 1.0
    terminal_current_a::Float64 = 0.0
    voltage_input_pu::Float64 = 0.0
    previous_voltage_input_pu::Float64 = 0.0
    sensor_closed::Bool = true
    update_count::Int = 0
    voltage_update_count::Int = 0
end

Base.@kwdef mutable struct SynchronousMachineDynamicState
    equation_state::SynchronousMachineEquationState
    current_history::Vector{Float64}
    electrical_coefficients::Vector{Float64}
    applied_torque_distribution::Vector{Float64} = Float64[]
    exciter_port_state::SynchronousMachineExciterPortState =
        SynchronousMachineExciterPortState()
    call_count::Int = 0
    companion_history_mutated::Bool = false
    saturation_enabled::Bool = false
    d_axis_saturation_region::Int = 0
    q_axis_saturation_region::Int = 0
    saturation_refactor_count::Int = 0
end

struct SynchronousMachineRotorMass
    torque_fraction::Float64
    inertia::Float64
    speed_deviation_damping::Float64
    mutual_damping::Float64
    shaft_stiffness::Float64
    absolute_speed_damping::Float64
end

struct SynchronousMachineParameters
    generator_mass_index::Int
    exciter_mass_index::Int
    pole_count::Int
    rated_power_mva::Float64
    rated_voltage_kv::Float64
    saturation_base_current::Float64
    d_axis_saturation_current_1::Float64
    d_axis_saturation_current_2::Float64
    saturation_voltage_ratios::NTuple{7,Float64}
    armature_resistance_pu::Float64
    leakage_reactance_pu::Float64
    d_axis_synchronous_reactance_pu::Float64
    q_axis_synchronous_reactance_pu::Float64
    d_axis_transient_reactance_pu::Float64
    q_axis_transient_reactance_pu::Float64
    d_axis_subtransient_reactance_pu::Float64
    q_axis_subtransient_reactance_pu::Float64
    d_axis_open_circuit_time_constant_s::Float64
    q_axis_open_circuit_time_constant_s::Float64
    d_axis_subtransient_time_constant_s::Float64
    q_axis_subtransient_time_constant_s::Float64
    zero_sequence_reactance_pu::Float64
    neutral_resistance_pu::Float64
    neutral_reactance_pu::Float64
    parameter_fitting_value::Float64
    delta_connected::Bool
    rotor_masses::Vector{SynchronousMachineRotorMass}
end

struct SynchronousMachineInitialization
    state::SynchronousMachineDynamicState
    phase_current_phasors::Vector{ComplexF64}
    numask::Int
    nlocg::Int
    nloce::Int
    angle_half_step_inverse::Float64
    speed_tolerance::Float64
    omega_tolerance::Float64
    speed_floor::Float64
    max_iterations::Int
    damping_ratio::Float64
    rotor_angle_extrapolation_interval::Float64
    speed_voltage_factor::Float64
    electrical_speed_rad_s::Float64
    electrical_angle_increment::Float64
end

Base.@kwdef mutable struct UniversalMachineRotorCurrentState
    current_values::Vector{Float64} = Float64[]
    solution_matrix::Matrix{Float64} = zeros(Float64, 0, 0)
    right_hand_side::Vector{Float64} = Float64[]
    rotor_current_solution_mutated::Bool = false
end

Base.@kwdef mutable struct UniversalMachineStatorExcitationCurrentState
    current_values::Vector{Float64} = Float64[]
    solution_matrix::Matrix{Float64} = zeros(Float64, 0, 0)
    right_hand_side::Vector{Float64} = Float64[]
    axis_kinds::Vector{Symbol} = Symbol[]
    stator_excitation_solution_mutated::Bool = false
end

Base.@kwdef mutable struct UniversalMachineDirectCouplingState
    current_values::Vector{Float64} = Float64[]
    coupling_fraction::Float64 = 0.0
    voltage_coupling_fraction::Float64 = 0.0
    resistance_coupling_fraction::Float64 = 0.0
    direct_coupling_mutated::Bool = false
end

Base.@kwdef mutable struct UniversalMachineHistoryCurrentState
    history_currents::Vector{Float64} = Float64[]
    rotor_history_matrix::Matrix{Float64} = zeros(Float64, 0, 0)
    stator_thevenin_drop::Vector{Float64} = Float64[]
    axis_kinds::Vector{Symbol} = Symbol[]
    history_current_update_mutated::Bool = false
end

Base.@kwdef mutable struct UniversalMachineMechanicalIterationState
    mechanical_speed_rad_s::Float64 = 0.0
    previous_mechanical_speed_rad_s::Float64 = 0.0
    mechanical_angle_rad::Float64 = 0.0
    predicted_speed_rad_s::Float64 = 0.0
    solved_speed_rad_s::Float64 = 0.0
    generated_torque::Float64 = 0.0
    torque_increment::Float64 = 0.0
    iteration_count::Int = 0
    converged::Bool = false
    mechanical_iteration_mutated::Bool = false
end

Base.@kwdef mutable struct UniversalMachineFluxSaturationState
    d_axis_flux::Float64 = 0.0
    q_axis_flux::Float64 = 0.0
    d_axis_current::Float64 = 0.0
    q_axis_current::Float64 = 0.0
    d_axis_saturation_offset::Float64 = 0.0
    q_axis_saturation_offset::Float64 = 0.0
    flux_saturation_mutated::Bool = false
end

Base.@kwdef mutable struct UniversalMachineFluxCurrentJacobianState
    current_values::Vector{Float64} = Float64[]
    d_axis_flux::Float64 = 0.0
    q_axis_flux::Float64 = 0.0
    d_axis_current::Float64 = 0.0
    q_axis_current::Float64 = 0.0
    base_d_axis_current::Float64 = 0.0
    base_q_axis_current::Float64 = 0.0
    d_axis_d_flux_sensitivity::Float64 = 0.0
    q_axis_d_flux_sensitivity::Float64 = 0.0
    d_axis_q_flux_sensitivity::Float64 = 0.0
    q_axis_q_flux_sensitivity::Float64 = 0.0
    d_axis_saturated::Bool = false
    q_axis_saturated::Bool = false
    iteration_count::Int = 0
    flux_current_jacobian_mutated::Bool = false
end

Base.@kwdef mutable struct UniversalMachineSharedTorquePredictionState
    predicted_torques::Vector{Float64} = Float64[]
    torque_history::Vector{Float64} = Float64[]
    current_substitution_values::Vector{Float64} = Float64[]
    active_machine_indices::Vector{Int} = Int[]
    written_substitution_indices::Vector{Int} = Int[]
    shared_torque_prediction_mutated::Bool = false
end

Base.@kwdef mutable struct UniversalMachinePostsolveState
    coil_parameters::Vector{Float64} = Float64[]
    source_crests::Vector{Float64} = Float64[]
    current_values::Vector{Float64} = Float64[]
    prediction_values::Vector{Float64} = Float64[]
    history_values::Vector{Float64} = Float64[]
    prediction_report_text_lines::Vector{String} = String[]
    d_axis_flux::Float64 = 0.0
    q_axis_flux::Float64 = 0.0
    theta_electric_rad::Float64 = 0.0
    machine_type::Int = 0
    input_mode::Int = 0
    prediction_loop_marker::Int = 0
    postsolve_update_mutated::Bool = false
    prediction_report_text_mutated::Bool = false
end

Base.@kwdef mutable struct UniversalMachineReportScheduleState
    last_report_step::Int = 0
    report_text_lines::Vector{String} = String[]
    report_schedule_mutated::Bool = false
    report_text_mutated::Bool = false
end

struct InductionMachineParameters
    machine_type::Int
    coil_conductances::Vector{Float64}
    coil_reactances::Vector{Float64}
    coil_predictor_factors::Vector{Float64}
    d_axis_coil_count::Int
    q_axis_coil_count::Int
    d_axis_unsaturated_inductance::Float64
    d_axis_saturated_inductance::Float64
    q_axis_unsaturated_inductance::Float64
    q_axis_saturated_inductance::Float64
    d_axis_saturation_flux::Float64
    q_axis_saturation_flux::Float64
    d_axis_remanent_flux::Float64
    q_axis_remanent_flux::Float64
    d_axis_saturation_mode::Int
    q_axis_saturation_mode::Int
    pole_pair_count::Int
    synchronous_electrical_speed_rad_s::Float64
    speed_tolerance_rad_s::Float64
    half_step_interval_s::Float64
    series_path_leakage_inductance_h::Float64
    effective_armature_leakage_inductance_h::Float64
    effective_compound_field_leakage_inductance_h::Float64
end

mutable struct InductionMachineState
    current_values::Vector{Float64}
    history_currents::Vector{Float64}
    mechanical_speed_rad_s::Float64
    previous_mechanical_speed_rad_s::Float64
    mechanical_angle_rad::Float64
    d_axis_flux::Float64
    q_axis_flux::Float64
    generated_torque::Float64
    output_values::Vector{Float64}
    call_count::Int
end

mutable struct MachineTerminalPredictionState
    previous_d_axis_flux_Wb::Float64
    previous_q_axis_flux_Wb::Float64
    previous_d_axis_internal_flux_Wb::Float64
    previous_q_axis_internal_flux_Wb::Float64
    phase_current_injections_A::Vector{Float64}
    update_count::Int
    initialized::Bool
end

MachineTerminalPredictionState() = MachineTerminalPredictionState(
    0.0,
    0.0,
    0.0,
    0.0,
    zeros(Float64, 3),
    0,
    false,
)

struct MachineConvergenceError <: Exception
    machine_type::Int
    call_index::Int
    iteration_count::Int
    speed_residual_rad_s::Float64
    speed_tolerance_rad_s::Float64
end

function Base.showerror(io::IO, error::MachineConvergenceError)
    print(
        io,
        "type-",
        error.machine_type,
        " coupled d/q machine failed to converge on call ",
        error.call_index,
        " after ",
        error.iteration_count,
        " iterations (speed residual ",
        error.speed_residual_rad_s,
        " rad/s, tolerance ",
        error.speed_tolerance_rad_s,
        " rad/s)",
    )
end

const UniversalMachineType4Parameters = InductionMachineParameters
const UniversalMachineType4State = InductionMachineState
const CoupledDQMachineParameters = InductionMachineParameters
const CoupledDQMachineState = InductionMachineState
const DirectCurrentMachineParameters = InductionMachineParameters
const DirectCurrentMachineState = InductionMachineState

struct InductionMachineSteadyStateEquivalent
    terminal_admittance::ComplexF64
    terminal_impedance::ComplexF64
    magnetizing_impedance::ComplexF64
    slip_referred_rotor_impedance::ComplexF64
    power_coil_resistance::Float64
end

function induction_machine_steady_state_equivalent(;
    power_coil_resistance::Real,
    rotor_resistance::Real,
    rotor_leakage_inductance::Real,
    main_inductance::Real,
    slip::Real,
    frequency_hz::Real,
)
    resistance = Float64(power_coil_resistance)
    rotor_resistance_value = Float64(rotor_resistance)
    rotor_leakage = Float64(rotor_leakage_inductance)
    main = Float64(main_inductance)
    slip_value = Float64(slip)
    frequency = Float64(frequency_hz)
    resistance >= 0.0 || throw(ArgumentError("power_coil_resistance must be nonnegative"))
    rotor_resistance_value > 0.0 || throw(ArgumentError("rotor_resistance must be positive"))
    rotor_leakage >= 0.0 || throw(ArgumentError("rotor_leakage_inductance must be nonnegative"))
    main > 0.0 || throw(ArgumentError("main_inductance must be positive"))
    abs(slip_value) > eps(Float64) || throw(ArgumentError("slip must be nonzero"))
    frequency > 0.0 || throw(ArgumentError("frequency_hz must be positive"))
    angular_frequency = 2.0 * pi * frequency
    magnetizing_impedance = complex(0.0, angular_frequency * main)
    rotor_impedance = complex(
        rotor_resistance_value / slip_value,
        angular_frequency * rotor_leakage,
    )
    parallel_impedance = inv(inv(magnetizing_impedance) + inv(rotor_impedance))
    terminal_impedance = resistance + parallel_impedance
    return InductionMachineSteadyStateEquivalent(
        inv(terminal_impedance),
        terminal_impedance,
        magnetizing_impedance,
        rotor_impedance,
        resistance,
    )
end

function induction_machine_initial_state(
    terminal_voltage_phasors::AbstractVector{<:Complex},
    equivalent::InductionMachineSteadyStateEquivalent;
    mechanical_speed_rad_s::Real,
    mechanical_angle_rad::Real,
    pole_pair_count::Integer,
    machine_type::Integer=4,
    d_axis_coil_count::Integer=1,
    q_axis_coil_count::Integer=1,
)
    length(terminal_voltage_phasors) == 3 ||
        throw(ArgumentError("induction-machine initialization requires three terminal phasors"))
    pole_pair_count > 0 || throw(ArgumentError("pole_pair_count must be positive"))
    machine_type in (3, 4, 5, 6, 7) ||
        throw(ArgumentError("induction-machine type must be between 3 and 7"))
    terminal_voltages = ComplexF64.(terminal_voltage_phasors)
    input_current_phasors = equivalent.terminal_admittance .* terminal_voltages
    parallel_voltage_phasors =
        terminal_voltages .- equivalent.power_coil_resistance .* input_current_phasors
    rotor_current_phasors =
        parallel_voltage_phasors ./ equivalent.slip_referred_rotor_impedance
    transform = _induction_machine_power_transform(
        Int(pole_pair_count) * Float64(mechanical_angle_rad),
        Int(machine_type),
    )
    power_currents = transform * (-real.(input_current_phasors))
    rotor_coordinates = transform * real.(rotor_current_phasors)
    d_count = Int(d_axis_coil_count)
    q_count = Int(q_axis_coil_count)
    d_count >= 0 && q_count >= 0 ||
        throw(ArgumentError("induction-machine axis coil counts must be nonnegative"))
    excitation_currents = vcat(
        fill(rotor_coordinates[2], d_count),
        fill(rotor_coordinates[3], q_count),
        machine_type == 4 ? [rotor_coordinates[1]] : Float64[],
    )
    current_values = vcat(power_currents, excitation_currents)
    coil_count = 3 + d_count + q_count + (machine_type == 4 ? 1 : 0)
    length(current_values) == coil_count ||
        throw(ArgumentError("induction-machine initial current count mismatch"))
    history_currents = vcat(current_values[1:3], -current_values[4:coil_count])
    return InductionMachineState(
        current_values,
        history_currents;
        mechanical_speed_rad_s = mechanical_speed_rad_s,
        previous_mechanical_speed_rad_s = mechanical_speed_rad_s,
        mechanical_angle_rad = mechanical_angle_rad,
    )
end

function single_phase_induction_steady_state_initialization(;
    power_thevenin_voltage_phasor::Complex,
    power_thevenin_impedance::Complex,
    excitation_voltage_phasor::Complex,
    excitation_thevenin_impedance::Complex,
    q_axis_excitation_voltage_phasor::Union{Nothing,Complex}=nothing,
    q_axis_excitation_thevenin_impedance::Complex=0.0 + 0.0im,
    coil_resistances::AbstractVector{<:Real},
    coil_inductances::AbstractVector{<:Real},
    d_axis_main_inductance::Real,
    q_axis_main_inductance::Real,
    slip::Real,
    frequency_hz::Real,
    mechanical_speed_rad_s::Real,
    mechanical_angle_rad::Real,
    pole_pair_count::Integer,
)
    excitation_count = q_axis_excitation_voltage_phasor === nothing ? 1 : 2
    coil_count = 3 + excitation_count
    length(coil_resistances) == coil_count && length(coil_inductances) == coil_count ||
        throw(ArgumentError("single-phase induction initialization requires $coil_count coils"))
    resistances = Float64.(coil_resistances)
    inductances = Float64.(coil_inductances)
    resistances[3] > 0.0 && all(>(0.0), resistances[4:end]) ||
        throw(ArgumentError("single-phase power and excitation resistances must be positive"))
    all(>=(0.0), inductances) ||
        throw(ArgumentError("single-phase coil inductances must be nonnegative"))
    d_main = Float64(d_axis_main_inductance)
    q_main = Float64(q_axis_main_inductance)
    d_main > 0.0 && q_main > 0.0 ||
        throw(ArgumentError("single-phase main inductances must be positive"))
    slip_value = Float64(slip)
    abs(slip_value) > eps(Float64) ||
        throw(ArgumentError("single-phase induction initialization requires nonzero slip"))
    frequency = Float64(frequency_hz)
    frequency > 0.0 || throw(ArgumentError("frequency_hz must be positive"))
    pole_pair_count > 0 || throw(ArgumentError("pole_pair_count must be positive"))

    power_voltage = ComplexF64(power_thevenin_voltage_phasor)
    power_impedance = ComplexF64(power_thevenin_impedance)
    excitation_voltages = ComplexF64[excitation_voltage_phasor]
    excitation_impedances = ComplexF64[excitation_thevenin_impedance]
    if q_axis_excitation_voltage_phasor !== nothing
        push!(excitation_voltages, q_axis_excitation_voltage_phasor)
        push!(excitation_impedances, q_axis_excitation_thevenin_impedance)
    end
    all(
        isfinite,
        (
            real(power_voltage),
            imag(power_voltage),
            real(power_impedance),
            imag(power_impedance),
            [real(value) for value in excitation_voltages]...,
            [imag(value) for value in excitation_voltages]...,
            [real(value) for value in excitation_impedances]...,
            [imag(value) for value in excitation_impedances]...,
        ),
    ) || throw(ArgumentError("single-phase steady-state phasors must be finite"))

    angular_frequency = 2.0 * pi * frequency
    power_d_reactance = angular_frequency * inductances[2]
    power_q_reactance = angular_frequency * inductances[3]
    positive_sequence_reactance =
        imag(power_impedance) + power_q_reactance +
        (power_d_reactance - power_q_reactance) / 2.0
    positive_sequence_resistance = real(power_impedance) + resistances[3]
    positive_sequence_denominator = 2.0 * (
        positive_sequence_reactance^2 + positive_sequence_resistance^2
    )
    positive_sequence_denominator > 0.0 ||
        throw(ArgumentError("single-phase positive-sequence impedance is singular"))

    function currents_from_flux(flux_values::AbstractVector{<:Real})
        d_flux_real, d_flux_imag, q_flux_real, q_flux_imag = Float64.(flux_values)
        induced_power_real = angular_frequency * (d_flux_real + q_flux_imag)
        induced_power_imag = angular_frequency * (d_flux_imag - q_flux_real)
        power_rhs_real = -sqrt(2.0) * real(power_voltage) + induced_power_real
        power_rhs_imag = -sqrt(2.0) * imag(power_voltage) + induced_power_imag
        d_power_real = (
            power_rhs_imag * positive_sequence_resistance -
            power_rhs_real * positive_sequence_reactance
        ) / positive_sequence_denominator
        q_power_real = (
            power_rhs_real * positive_sequence_resistance +
            power_rhs_imag * positive_sequence_reactance
        ) / positive_sequence_denominator
        d_power_imag = -q_power_real
        q_power_imag = d_power_real

        excitation_real_currents = Float64[]
        excitation_imaginary_currents = Float64[]
        axis_fluxes = [(d_flux_real, d_flux_imag), (q_flux_real, q_flux_imag)]
        for excitation_index in 1:excitation_count
            coil_index = 3 + excitation_index
            excitation_voltage = excitation_voltages[excitation_index]
            excitation_impedance = excitation_impedances[excitation_index]
            axis_flux_real, axis_flux_imag = axis_fluxes[excitation_index]
            excitation_conductance = inv(resistances[coil_index])
            slip_reactance_ratio = slip_value * angular_frequency *
                                    inductances[coil_index] * excitation_conductance
            excitation_factor = inv(1.0 + slip_reactance_ratio^2)
            excitation_admittance = excitation_conductance * excitation_factor
            flux_excitation_admittance = slip_value * angular_frequency *
                                         excitation_conductance * excitation_factor
            excitation_real = -excitation_admittance * (
                real(excitation_voltage) +
                slip_reactance_ratio * imag(excitation_voltage)
            ) + flux_excitation_admittance * (
                axis_flux_imag - slip_reactance_ratio * axis_flux_real
            )
            excitation_imag = -excitation_admittance * (
                imag(excitation_voltage) -
                slip_reactance_ratio * real(excitation_voltage)
            ) - flux_excitation_admittance * (
                slip_reactance_ratio * axis_flux_imag + axis_flux_real
            )

            excitation_network_real =
                excitation_admittance * real(excitation_impedance) +
                slip_reactance_ratio * excitation_admittance *
                imag(excitation_impedance)
            excitation_network_imag =
                excitation_admittance * imag(excitation_impedance) -
                slip_reactance_ratio * excitation_admittance *
                real(excitation_impedance)
            excitation_matrix = [
                1.0 + excitation_network_real -excitation_network_imag
                excitation_network_imag 1.0 + excitation_network_real
            ]
            excitation_real, excitation_imag =
                excitation_matrix \ [excitation_real, excitation_imag]
            push!(excitation_real_currents, excitation_real)
            push!(excitation_imaginary_currents, excitation_imag)
        end

        real_currents = vcat(
            0.0,
            d_power_real,
            q_power_real,
            excitation_real_currents,
        )
        imaginary_currents = vcat(
            0.0,
            d_power_imag,
            q_power_imag,
            excitation_imaginary_currents,
        )
        axis_currents = [
            d_power_real + excitation_real_currents[1],
            d_power_imag + excitation_imaginary_currents[1],
            q_power_real + get(excitation_real_currents, 2, 0.0),
            q_power_imag + get(excitation_imaginary_currents, 2, 0.0),
        ]
        return axis_currents, real_currents, imaginary_currents
    end

    base_axis_currents = first(currents_from_flux(zeros(4)))
    flux_sensitivity = hcat([
        first(currents_from_flux([index == probe ? 1.0 : 0.0 for index in 1:4])) .-
        base_axis_currents for probe in 1:4
    ]...)
    main_inductances = Diagonal([d_main, d_main, q_main, q_main])
    flux_phasor_values = (
        Matrix{Float64}(I, 4, 4) - main_inductances * flux_sensitivity
    ) \ (main_inductances * base_axis_currents)
    _, real_currents, imaginary_currents = currents_from_flux(flux_phasor_values)

    state = InductionMachineState(
        real_currents,
        vcat(real_currents[1:3], -real_currents[4:end]);
        mechanical_speed_rad_s = mechanical_speed_rad_s,
        previous_mechanical_speed_rad_s = mechanical_speed_rad_s,
        mechanical_angle_rad = mechanical_angle_rad,
    )
    state.d_axis_flux = flux_phasor_values[1]
    state.q_axis_flux = flux_phasor_values[3]
    terminal_current_phasor = complex(
        -imaginary_currents[2] + real_currents[3],
        real_currents[2] + imaginary_currents[3],
    )
    terminal_voltage_phasor =
        power_voltage + power_impedance * terminal_current_phasor
    return (
        source = :single_phase_induction_steady_state_initialization,
        state = state,
        power_terminal_voltage_phasor = terminal_voltage_phasor,
        power_terminal_current_phasor = terminal_current_phasor,
        excitation_current_phasor = complex(real_currents[4], imaginary_currents[4]),
        d_axis_excitation_current_phasor =
            complex(real_currents[4], imaginary_currents[4]),
        q_axis_excitation_current_phasor = excitation_count == 2 ?
            complex(real_currents[5], imaginary_currents[5]) : 0.0 + 0.0im,
        d_axis_flux_phasor = complex(flux_phasor_values[1], flux_phasor_values[2]),
        q_axis_flux_phasor = complex(flux_phasor_values[3], flux_phasor_values[4]),
        real_current_values = real_currents,
        imaginary_current_values = imaginary_currents,
        initialized = true,
    )
end

function InductionMachineParameters(
    coil_resistances::AbstractVector{<:Real},
    coil_inductances::AbstractVector{<:Real};
    time_step_s::Real,
    d_axis_unsaturated_inductance::Real,
    q_axis_unsaturated_inductance::Real,
    d_axis_saturated_inductance::Real=0.0,
    q_axis_saturated_inductance::Real=0.0,
    d_axis_saturation_flux::Real=0.0,
    q_axis_saturation_flux::Real=0.0,
    d_axis_remanent_flux::Real=0.0,
    q_axis_remanent_flux::Real=0.0,
    d_axis_saturation_mode::Integer=0,
    q_axis_saturation_mode::Integer=0,
    pole_pair_count::Integer,
    d_axis_coil_count::Integer=1,
    q_axis_coil_count::Integer=1,
    synchronous_electrical_speed_rad_s::Real=0.0,
    speed_tolerance_rad_s::Real,
    machine_type::Integer=length(coil_resistances) == 5 ? 3 : 4,
    power_leakage_owned_by_network::Bool=true,
)
    machine_type in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) ||
        throw(ArgumentError("coupled d/q machine type must be between 1 and 12"))
    d_count = Int(d_axis_coil_count)
    q_count = Int(q_axis_coil_count)
    d_count >= 0 && q_count >= 0 ||
        throw(ArgumentError("coupled d/q axis coil counts must be nonnegative"))
    machine_type == 6 && (d_count != 1 || q_count != 0) &&
        throw(ArgumentError("type-6 machine requires one d-axis and zero q-axis field coils"))
    machine_type == 7 && (d_count != 1 || q_count != 1) &&
        throw(ArgumentError("type-7 machine requires one d-axis and one q-axis rotor coil"))
    machine_type == 8 && (d_count != 1 || q_count != 0) &&
        throw(ArgumentError("type-8 separately excited DC machine requires one d-axis field coil"))
    machine_type == 9 && (d_count != 2 || q_count != 0) &&
        throw(ArgumentError("type-9 series-compound DC machine requires two d-axis field coils"))
    machine_type == 10 && (d_count != 2 || q_count != 0) &&
        throw(ArgumentError("type-10 series-field DC machine requires one inactive and one series d-axis field slot"))
    machine_type == 11 && (d_count != 2 || q_count != 0) &&
        throw(ArgumentError("type-11 parallel-compound DC machine requires compound and series d-axis field coils"))
    machine_type == 12 && (d_count != 2 || q_count != 0) &&
        throw(ArgumentError("type-12 self-excited shunt DC machine requires shunt and inactive-series d-axis field slots"))
    coil_count = 3 + d_count + q_count + (machine_type == 4 ? 1 : 0)
    length(coil_resistances) == coil_count && length(coil_inductances) == coil_count ||
        throw(ArgumentError("type-$machine_type coupled d/q machine parameters require $coil_count coils"))
    resistances = _machine_real_vector("coil_resistances", coil_resistances, coil_count)
    reactances = _machine_real_vector("coil_inductances", coil_inductances, coil_count)
    if power_leakage_owned_by_network &&
       machine_type ∉ (6, 7, 8, 9, 10, 11, 12)
        leakage_shift = min(reactances[2], reactances[3])
        reactances[2] -= leakage_shift
        reactances[3] -= leakage_shift
    end
    series_leakage = machine_type in (9, 10, 11, 12) ? reactances[5] : 0.0
    if machine_type in (9, 10, 11, 12)
        reactances[3] += series_leakage
        reactances[4] += series_leakage
    end
    conductances = Float64[value == 0.0 ? 0.0 : inv(value) for value in resistances]
    step = Float64(time_step_s)
    step > 0.0 || throw(ArgumentError("time_step_s must be positive"))
    half_step = step / 2.0
    predictors = Float64[
        inv(1.0 + conductances[index] * reactances[index] / half_step)
        for index in eachindex(conductances)
    ]
    pole_pair_count > 0 || throw(ArgumentError("pole_pair_count must be positive"))
    synchronous_speed = Float64(synchronous_electrical_speed_rad_s)
    synchronous_speed >= 0.0 ||
        throw(ArgumentError("synchronous_electrical_speed_rad_s must be nonnegative"))
    machine_type in (1, 2) && synchronous_speed <= 0.0 &&
        throw(ArgumentError("type-$machine_type machine synchronous electrical speed must be positive"))
    tolerance = Float64(speed_tolerance_rad_s)
    tolerance >= 0.0 || throw(ArgumentError("speed_tolerance_rad_s must be nonnegative"))
    return InductionMachineParameters(
        Int(machine_type),
        conductances,
        reactances,
        predictors,
        d_count,
        q_count,
        Float64(d_axis_unsaturated_inductance),
        Float64(d_axis_saturated_inductance),
        Float64(q_axis_unsaturated_inductance),
        Float64(q_axis_saturated_inductance),
        Float64(d_axis_saturation_flux),
        Float64(q_axis_saturation_flux),
        Float64(d_axis_remanent_flux),
        Float64(q_axis_remanent_flux),
        Int(d_axis_saturation_mode),
        Int(q_axis_saturation_mode),
        Int(pole_pair_count),
        synchronous_speed,
        tolerance,
        half_step,
        series_leakage,
        reactances[3],
        reactances[4],
    )
end

function InductionMachineState(
    current_values::AbstractVector{<:Real},
    history_currents::AbstractVector{<:Real};
    mechanical_speed_rad_s::Real,
    previous_mechanical_speed_rad_s::Real=mechanical_speed_rad_s,
    mechanical_angle_rad::Real,
)
    coil_count = length(current_values)
    coil_count in (4, 5, 6) ||
        throw(ArgumentError("coupled d/q machine state requires four, five, or six currents"))
    currents = _machine_real_vector("current_values", current_values, coil_count)
    histories = _machine_real_vector("history_currents", history_currents, coil_count)
    return InductionMachineState(
        currents,
        histories,
        Float64(mechanical_speed_rad_s),
        Float64(previous_mechanical_speed_rad_s),
        Float64(mechanical_angle_rad),
        0.0,
        0.0,
        0.0,
        zeros(Float64, coil_count + 3),
        0,
    )
end

Base.@kwdef mutable struct MachineTerminalCurrentState
    call_indices::Vector{Int} = Int[]
    row_indices::Vector{Int} = Int[]
    times_s::Vector{Float64} = Float64[]
    phase_nodes::Vector{Int} = Int[]
    terminal_nodes::Vector{Int} = Int[]
    output_slots::Vector{Int} = Int[]
    coupling_flags::Vector{Int} = Int[]
    terminal_currents::Vector{Float64} = Float64[]
    output_values::Vector{Float64} = Float64[]
    row_count::Int = 0
    call_count::Int = 0
    rows_per_call::Int = 0
    terminal_output_slot_count::Int = 0
    first_time_s::Float64 = 0.0
    last_time_s::Float64 = 0.0
    terminal_rows_mutated::Bool = false
end

Base.@kwdef mutable struct MachineNetworkCouplingState
    header_call_indices::Vector{Int} = Int[]
    machine_counts::Vector{Int} = Int[]
    coupling_row_counts::Vector{Int} = Int[]
    output_counts::Vector{Int} = Int[]
    source_constant_starts::Vector{Int} = Int[]
    branch_counts::Vector{Int} = Int[]
    vector_call_indices::Vector{Int} = Int[]
    vector_names::Vector{Symbol} = Symbol[]
    vector_indices::Vector{Int} = Int[]
    vector_values::Vector{Float64} = Float64[]
    source_call_indices::Vector{Int} = Int[]
    source_stage_names::Vector{Symbol} = Symbol[]
    source_kinds::Vector{Symbol} = Symbol[]
    source_indices::Vector{Int} = Int[]
    source_has_nodes::Vector{Bool} = Bool[]
    source_nodes::Vector{Int} = Int[]
    source_constant_indices::Vector{Int} = Int[]
    source_frequencies_hz::Vector{Float64} = Float64[]
    source_has_crests::Vector{Bool} = Bool[]
    source_crests::Vector{Float64} = Float64[]
    source_has_time1::Vector{Bool} = Bool[]
    source_time1_s::Vector{Float64} = Float64[]
    source_has_start_times::Vector{Bool} = Bool[]
    source_start_times_s::Vector{Float64} = Float64[]
    source_has_stop_times::Vector{Bool} = Bool[]
    source_stop_times_s::Vector{Float64} = Float64[]
    exit_call_indices::Vector{Int} = Int[]
    exit_start_indices::Vector{Int} = Int[]
    exit_source_constant_starts::Vector{Int} = Int[]
    exit_branch_counts::Vector{Int} = Int[]
    exit_reference_frequencies_rad_s::Vector{Float64} = Float64[]
    exit_twopi_values::Vector{Float64} = Float64[]
    header_count::Int = 0
    vector_row_count::Int = 0
    created_source_count::Int = 0
    changed_source_count::Int = 0
    slack_source_count::Int = 0
    exit_count::Int = 0
    machine_network_rows_mutated::Bool = false
end

function SynchronousMachineEquationState(
    cu_values::AbstractVector{<:Real},
    histq_values::AbstractVector{<:Real},
    shp_values::AbstractVector{<:Real},
)
    return SynchronousMachineEquationState(
        cu_values = _machine_float_vector(cu_values),
        histq_values = _machine_float_vector(histq_values),
        shp_values = _machine_float_vector(shp_values),
        iteration_count = 0,
        converged = false,
        equation_mutated = false,
    )
end

function SynchronousMachineDynamicState(
    equation_state::SynchronousMachineEquationState,
    current_history::AbstractVector{<:Real},
    electrical_coefficients::AbstractVector{<:Real},
    ;
    saturation_enabled::Union{Nothing,Bool}=nothing,
    d_axis_saturation_region::Union{Nothing,Int}=nothing,
    q_axis_saturation_region::Union{Nothing,Int}=nothing,
    saturation_refactor_count::Int=0,
    applied_torque_distribution::AbstractVector{<:Real}=Float64[],
    exciter_port_state::SynchronousMachineExciterPortState=
        SynchronousMachineExciterPortState(),
)
    length(current_history) == 3 ||
        throw(ArgumentError("current_history must contain three phase values"))
    coefficients = _machine_float_vector(electrical_coefficients)
    torque_distribution = _machine_float_vector(applied_torque_distribution)
    all(isfinite, torque_distribution) || throw(ArgumentError(
        "applied_torque_distribution entries must be finite",
    ))
    enabled = if saturation_enabled === nothing
        length(coefficients) >= 30 && coefficients[22] > 0.0 &&
        coefficients[23] != 0.0 && coefficients[24] > 0.0 && coefficients[25] != 0.0
    else
        saturation_enabled
    end
    flux_magnitude = 0.0
    if enabled
        for (factor_index, threshold_index, slope_index) in ((29, 22, 23), (30, 24, 25))
            factor = coefficients[factor_index]
            factor < 1.0 && (flux_magnitude = max(
                flux_magnitude,
                coefficients[threshold_index] +
                (inv(factor) - 1.0) / coefficients[slope_index],
            ))
        end
    end
    d_region = d_axis_saturation_region === nothing ?
               (enabled ? _synchronous_machine_saturation_region(
                   flux_magnitude,
                   coefficients[22],
               ) : 0) : d_axis_saturation_region
    q_region = q_axis_saturation_region === nothing ?
               (enabled ? _synchronous_machine_saturation_region(
                   flux_magnitude,
                   coefficients[24],
               ) : 0) : q_axis_saturation_region
    return SynchronousMachineDynamicState(
        equation_state = equation_state,
        current_history = _machine_float_vector(current_history),
        electrical_coefficients = coefficients,
        applied_torque_distribution = torque_distribution,
        exciter_port_state = exciter_port_state,
        saturation_enabled = enabled,
        d_axis_saturation_region = d_region,
        q_axis_saturation_region = q_region,
        saturation_refactor_count = saturation_refactor_count,
    )
end

function SynchronousMachineExciterPortState(current_reduction_factor::Real)
    factor = Float64(current_reduction_factor)
    isfinite(factor) && factor > 0.0 || throw(ArgumentError(
        "current_reduction_factor must be finite and positive",
    ))
    return SynchronousMachineExciterPortState(current_reduction_factor = factor)
end

function synchronous_machine_exciter_port_update!(
    state::SynchronousMachineExciterPortState,
    field_winding_current_a::Real,
    ;
    external_field_voltage_input_pu::Union{Nothing,Real}=nothing,
)
    state.sensor_closed || throw(ArgumentError(
        "the detailed-machine exciter-current sensor must remain closed",
    ))
    current = Float64(field_winding_current_a)
    isfinite(current) || throw(ArgumentError(
        "field_winding_current_a must be finite",
    ))
    voltage_input = external_field_voltage_input_pu === nothing ? nothing :
        Float64(external_field_voltage_input_pu)
    voltage_input === nothing || isfinite(voltage_input) || throw(ArgumentError(
        "external_field_voltage_input_pu must be finite",
    ))
    state.field_winding_current_a = current
    state.terminal_current_a = current * state.current_reduction_factor
    state.update_count += 1
    if voltage_input !== nothing
        state.previous_voltage_input_pu = state.voltage_input_pu
        state.voltage_input_pu = voltage_input
        state.voltage_update_count += 1
    end
    return (
        field_winding_current_a = state.field_winding_current_a,
        current_reduction_factor = state.current_reduction_factor,
        terminal_current_a = state.terminal_current_a,
        voltage_input_pu = state.voltage_input_pu,
        previous_voltage_input_pu = state.previous_voltage_input_pu,
        sensor_closed = state.sensor_closed,
        update_count = state.update_count,
        voltage_update_count = state.voltage_update_count,
        current_mutated = true,
    )
end

function _synchronous_machine_fit_axis(
    synchronous_reactance::Float64,
    transient_reactance::Float64,
    subtransient_reactance::Float64,
    subtransient_time_rad::Float64,
    transient_time_rad::Float64,
    leakage_reactance::Float64,
    improve_resistances::Bool,
    coincident_reactance_factor::Float64,
)
    d1 = synchronous_reactance
    h3 = leakage_reactance
    h2 = transient_reactance
    d1 == h2 && (h2 *= coincident_reactance_factor)
    a = subtransient_reactance - h3
    c = d1 - h3
    u2 = c * c
    if d1 == h2
        f1 = 1.0
        f2 = u2 / (d1 - subtransient_reactance)
        f3 = 0.0
        f4 = f2 / subtransient_time_rad
    else
        b = u2 / (d1 - h2) - c
        f1 = b + c
        f3 = f1 / transient_time_rad
        d = b * c
        a = -a / (a * f1 / d - 1.0)
        f2 = a + c
        f4 = (a + d / f1) / subtransient_time_rad
        if improve_resistances
            u = (subtransient_time_rad + transient_time_rad) / 2.0
            discriminant =
                u * u - subtransient_time_rad * transient_time_rad /
                        (1.0 - u2 / (f1 * f2))
            if discriminant >= 0.0
                d = u - sqrt(discriminant)
                f4 = f2 / d
                f3 = f1 / (2.0 * u - d)
            end
        end
    end
    return (mutual_reactance = c, rotor_self_1 = f1, rotor_self_2 = f2,
            rotor_resistance_1 = f3, rotor_resistance_2 = f4)
end

function _synchronous_machine_saturation_constants!(
    elp::Vector{Float64},
    parameters::SynchronousMachineParameters,
    leakage_reactance::Float64,
)
    base_current = abs(parameters.saturation_base_current)
    base_current > 0.0 || return false
    d_current_1 = parameters.d_axis_saturation_current_1
    d_current_2 = parameters.d_axis_saturation_current_2
    d_current_1 >= base_current && d_current_2 > d_current_1 ||
        throw(ArgumentError("inconsistent d-axis saturation currents"))
    ad1, ad2, aq1, aq2, qbase, qcurrent1, qcurrent2 =
        parameters.saturation_voltage_ratios
    ad1 <= 0.0 && (ad1 = 1.0)
    ad2 <= 0.0 && (ad2 = Float64(Float32(1.2)))
    ratio1 = d_current_1 / (base_current * ad1)
    ratio2 = d_current_2 / (base_current * ad2)
    slope = (ratio2 - ratio1) / (d_current_2 - d_current_1)
    threshold = d_current_2 + (1.0 - ratio2) / slope
    elp[22] = threshold
    elp[23] = slope
    if qbase <= 0.0
        axis_ratio = (elp[9] - leakage_reactance) / (elp[1] - leakage_reactance)
        elp[24] = threshold * axis_ratio
        elp[25] = slope
    else
        qcurrent1 >= qbase && qcurrent2 > qcurrent1 ||
            throw(ArgumentError("inconsistent q-axis saturation currents"))
        aq1 <= 0.0 && (aq1 = 1.0)
        aq2 <= 0.0 && (aq2 = Float64(Float32(1.2)))
        ratio1 = qcurrent1 / (qbase * aq1)
        ratio2 = qcurrent2 / (qbase * aq2)
        slope = (ratio2 - ratio1) / (qcurrent2 - qcurrent1)
        elp[24] = qcurrent2 + (1.0 - ratio2) / slope
        elp[25] = slope
    end
    return true
end

function _synchronous_machine_physical_storage(
    parameters::SynchronousMachineParameters,
    electrical_speed_rad_s::Float64,
)
    electrical_speed_rad_s > 0.0 ||
        throw(ArgumentError("electrical_speed_rad_s must be positive"))
    parameters.rated_power_mva > 0.0 ||
        throw(ArgumentError("rated_power_mva must be positive"))
    parameters.rated_voltage_kv > 0.0 ||
        throw(ArgumentError("rated_voltage_kv must be positive"))
    parameters.pole_count > 0 && iseven(parameters.pole_count) ||
        throw(ArgumentError("pole_count must be a positive even integer"))
    numask = length(parameters.rotor_masses)
    numask > 0 || throw(ArgumentError("at least one rotor mass is required"))
    1 <= parameters.generator_mass_index <= numask ||
        throw(ArgumentError("generator_mass_index is outside rotor masses"))
    0 <= parameters.exciter_mass_index <= numask ||
        throw(ArgumentError("exciter_mass_index is outside rotor masses"))

    base_impedance = parameters.rated_voltage_kv^2 / parameters.rated_power_mva
    emf_scale = parameters.rated_voltage_kv /
                (1.0e-3 * abs(parameters.saturation_base_current))
    isfinite(emf_scale) || throw(ArgumentError("saturation_base_current must be nonzero"))
    leakage = parameters.leakage_reactance_pu * base_impedance
    d_sync = parameters.d_axis_synchronous_reactance_pu * base_impedance
    q_sync = parameters.q_axis_synchronous_reactance_pu * base_impedance
    armature_resistance = parameters.armature_resistance_pu * base_impedance
    fit_value = parameters.parameter_fitting_value
    improve_resistances = fit_value <= 1.0
    coincident_factor = fit_value > 1.0 ? 1.0 : fit_value

    d_axis = _synchronous_machine_fit_axis(
        d_sync,
        parameters.d_axis_transient_reactance_pu * base_impedance,
        parameters.d_axis_subtransient_reactance_pu * base_impedance,
        parameters.d_axis_subtransient_time_constant_s * electrical_speed_rad_s,
        parameters.d_axis_open_circuit_time_constant_s * electrical_speed_rad_s,
        leakage,
        improve_resistances,
        coincident_factor,
    )
    transform_ratio = emf_scale / (d_axis.mutual_reactance / sqrt(3.0 / 2.0))
    rotor_scale = 2.0 * transform_ratio^2 / 3.0
    mutual_scale = emf_scale * transform_ratio / sqrt(3.0 / 2.0)
    armature_to_field_ratio = emf_scale / d_axis.mutual_reactance

    q_reduced =
        parameters.q_axis_synchronous_reactance_pu == parameters.q_axis_transient_reactance_pu ==
        parameters.q_axis_subtransient_reactance_pu
    q_axis =
        q_reduced ? nothing :
        _synchronous_machine_fit_axis(
            q_sync,
            parameters.q_axis_transient_reactance_pu * base_impedance,
            parameters.q_axis_subtransient_reactance_pu * base_impedance,
            parameters.q_axis_subtransient_time_constant_s * electrical_speed_rad_s,
            parameters.q_axis_open_circuit_time_constant_s * electrical_speed_rad_s,
            leakage,
            improve_resistances,
            coincident_factor,
        )

    elp = zeros(Float64, 101)
    elp[1] = d_sync
    elp[2] = emf_scale
    elp[3] = d_axis.rotor_self_1 * rotor_scale
    elp[4] = emf_scale
    elp[5] = mutual_scale
    elp[6] = d_axis.rotor_self_2 * rotor_scale
    elp[7] = d_axis.rotor_resistance_1 * rotor_scale
    elp[8] = d_axis.rotor_resistance_2 * rotor_scale
    elp[9] = q_sync
    if q_axis === nothing
        elp[11] = rotor_scale
        elp[14] = rotor_scale
    else
        q_mutual = q_axis.mutual_reactance * transform_ratio / sqrt(3.0 / 2.0)
        elp[10] = q_mutual
        elp[11] = q_axis.rotor_self_1 * rotor_scale
        elp[12] = q_mutual
        elp[13] = q_axis.mutual_reactance * rotor_scale
        elp[14] = q_axis.rotor_self_2 * rotor_scale
        elp[15] = q_axis.rotor_resistance_1 * rotor_scale
        elp[16] = q_axis.rotor_resistance_2 * rotor_scale
    end
    elp[17] = parameters.zero_sequence_reactance_pu * base_impedance +
              3.0 * parameters.neutral_reactance_pu
    elp[18] = armature_resistance + 3.0 * parameters.neutral_resistance_pu
    elp[19] = leakage
    elp[20] = armature_resistance
    elp[21] = armature_to_field_ratio
    elp[26] = parameters.pole_count / 2.0
    saturation_enabled = parameters.saturation_base_current < 0.0
    saturation_enabled &&
        _synchronous_machine_saturation_constants!(elp, parameters, leakage)

    shp = zeros(Float64, 24 * numask)
    num2 = 2 * numask
    inertia_scale = 0.0421409745
    mechanical_scale = 1.356306e6
    total_fraction = sum(mass.torque_fraction for mass in parameters.rotor_masses)
    numask == 1 && (total_fraction = 1.0)
    abs(total_fraction) > eps(Float64) ||
        throw(ArgumentError("rotor mass torque fractions must have a nonzero sum"))
    for (index, mass) in enumerate(parameters.rotor_masses)
        shp[num2 + index] = mass.torque_fraction / total_fraction
        shp[num2 + numask + index] = mass.inertia * inertia_scale
        shp[num2 + 2 * numask + index] = mass.mutual_damping * mechanical_scale
        shp[num2 + 3 * numask + index] = mass.shaft_stiffness * mechanical_scale / 1.0e6
        shp[num2 + 4 * numask + index] = mass.absolute_speed_damping * mechanical_scale
        shp[num2 + 5 * numask + index] = mass.speed_deviation_damping * mechanical_scale
    end
    return (; elp, shp, saturation_enabled)
end

function _synchronous_machine_reduce_matrix!(values::Vector{Float64}, order::Int, retained::Int)
    j = order
    ik = order^2
    nk = ik - order
    while j > retained
        reciprocal = inv(values[ik])
        row = copy(values[(nk + 1):(nk + order)])
        k = 1
        while k <= order
            if k == j
                k += 1
                continue
            end
            base = (k - 1) * order
            multiplier = -values[base + j] * reciprocal
            for i in 1:order
                values[base + i] += multiplier * row[i]
            end
            values[base + j] = multiplier
            k += 1
        end
        for k in 1:order
            values[nk + k] = row[k] * reciprocal
        end
        values[ik] = reciprocal
        j -= 1
        ik -= order + 1
        nk -= order
    end
    return values
end

function _synchronous_machine_saturation_region(
    flux_magnitude::Float64,
    threshold::Float64,
)
    threshold > 0.0 || throw(ArgumentError("saturation threshold must be positive"))
    flux_magnitude <= threshold && return 0
    # The machine characteristic is defined on the source model's promoted REAL grid.
    segment_origin = Float64(Float32(0.9))
    raw_segment = 10.0 * (flux_magnitude / threshold - segment_origin)
    coarse_segment = if typemin(Int32) <= raw_segment <= typemax(Int32)
        Int(trunc(Int32, raw_segment))
    else
        Int(typemin(Int32))
    end
    numerator = coarse_segment + 1
    return numerator >= 0 ? div(numerator, 2) : -div(-numerator, 2)
end

function _synchronous_machine_saturation_factor(
    flux_magnitude::Float64,
    threshold::Float64,
    slope::Float64,
)
    flux_magnitude <= threshold && return 1.0
    return inv(1.0 + slope * (flux_magnitude - threshold))
end

function _synchronous_machine_incremental_saturation_factor(
    threshold::Float64,
    slope::Float64,
    region::Int,
)
    region == 0 && return 1.0
    segment = Float64(region) * Float64(Float32(0.1))
    lower_flux = threshold * (Float64(Float32(0.9)) + segment)
    upper_flux = threshold * (Float64(Float32(1.1)) + segment)
    lower_current = lower_flux / (1.0 + slope * (lower_flux - threshold))
    upper_current = upper_flux / (1.0 + slope * (upper_flux - threshold))
    return (upper_current - lower_current) / (upper_flux - lower_flux)
end

function _synchronous_machine_refactor_electrical_companion!(
    elp::Vector{Float64},
    d_axis_region::Int,
    q_axis_region::Int,
    damping_ratio::Float64,
)
    length(elp) >= 101 ||
        throw(ArgumentError("electrical_coefficients must contain at least 101 values"))
    d_increment = _synchronous_machine_incremental_saturation_factor(
        elp[22],
        elp[23],
        d_axis_region,
    )
    q_increment = _synchronous_machine_incremental_saturation_factor(
        elp[24],
        elp[25],
        q_axis_region,
    )
    companion_factor = elp[66] / elp[2]
    z = zeros(Float64, 55)
    z[50] = (elp[1] * d_increment + elp[19]) * companion_factor
    z[51] = (elp[9] * q_increment + elp[19]) * companion_factor
    z[52] = elp[3] + elp[62] * d_increment
    z[53] = elp[6] + elp[63] * d_increment
    z[54] = elp[11] + elp[64] * q_increment
    z[55] = elp[14] + elp[65] * q_increment
    z[1] = z[50] + elp[20]
    z[8] = z[51] + elp[20]
    z[15] = z[52] + elp[7]
    z[22] = z[53] + elp[8]
    z[29] = z[54] + elp[15]
    z[36] = z[55] + elp[16]
    z[11] = elp[68] * q_increment
    z[26] = z[11]
    z[12] = elp[69] * q_increment
    z[32] = z[12]
    z[30] = elp[13] * q_increment
    z[35] = z[30]
    z[3] = elp[66] * d_increment
    z[13] = z[3]
    z[4] = elp[67] * d_increment
    z[19] = z[4]
    z[16] = elp[5] * d_increment
    z[21] = z[16]

    reduced = copy(z[1:36])
    _synchronous_machine_reduce_matrix!(reduced, 6, 2)
    d_axis_reduced = reduced[1]
    q_axis_reduced = reduced[8]
    reduced[8] -= d_axis_reduced
    phase_positive_inverse = inv(d_axis_reduced)
    zero_inverse = elp[17]
    elp[27] = (zero_inverse + 2.0 * phase_positive_inverse) / 3.0
    elp[28] = (zero_inverse - phase_positive_inverse) / 3.0
    elp[80] = phase_positive_inverse / sqrt(3.0 / 2.0)
    elp[81] = zero_inverse / sqrt(3.0)
    for (destination, source) in zip(32:43, (3, 4, 11, 12, 15, 21, 16, 22, 29, 35, 30, 36))
        elp[destination] = reduced[source]
    end
    elp[70] = reduced[13]
    elp[71] = reduced[19]
    elp[72] = reduced[26]
    elp[73] = reduced[32]
    elp[74] = phase_positive_inverse
    elp[75] = elp[19] * elp[66]
    elp[76] = q_axis_reduced
    elp[79] = reduced[8]
    elp[44] = z[50] - elp[20] * damping_ratio
    elp[45] = z[13]
    elp[46] = z[19]
    elp[47] = z[51] - elp[20] * damping_ratio
    elp[48] = z[26]
    elp[49] = z[32]
    elp[50] = z[3]
    elp[51] = z[52] - elp[7] * damping_ratio
    elp[52] = z[21]
    elp[53] = z[4]
    elp[54] = z[16]
    elp[55] = z[53] - elp[8] * damping_ratio
    elp[56] = z[11]
    elp[57] = z[54] - elp[15] * damping_ratio
    elp[58] = z[35]
    elp[59] = z[12]
    elp[60] = z[30]
    elp[61] = z[55] - elp[16] * damping_ratio
    return (; d_increment, q_increment, phase_positive_inverse, zero_inverse)
end

function _synchronous_machine_saturation_update!(
    elp::Vector{Float64};
    d_axis_current::Float64,
    q_axis_current::Float64,
    d_axis_rotor_current_1::Float64,
    d_axis_rotor_current_2::Float64,
    q_axis_rotor_current_1::Float64,
    q_axis_rotor_current_2::Float64,
    d_axis_region::Int,
    q_axis_region::Int,
    damping_ratio::Float64,
)
    d_flux = d_axis_current * elp[21] +
             d_axis_rotor_current_1 + d_axis_rotor_current_2
    q_flux = (q_axis_current * elp[21] +
              q_axis_rotor_current_1 + q_axis_rotor_current_2) * elp[31]
    flux_magnitude = hypot(d_flux, q_flux)
    next_d_region = _synchronous_machine_saturation_region(flux_magnitude, elp[22])
    next_q_region = _synchronous_machine_saturation_region(flux_magnitude, elp[24])
    elp[29] = _synchronous_machine_saturation_factor(flux_magnitude, elp[22], elp[23])
    elp[30] = _synchronous_machine_saturation_factor(flux_magnitude, elp[24], elp[25])
    refactorized = next_d_region != d_axis_region || next_q_region != q_axis_region
    refactorized && _synchronous_machine_refactor_electrical_companion!(
        elp,
        next_d_region,
        next_q_region,
        damping_ratio,
    )
    return (;
        flux_magnitude,
        d_axis_region = next_d_region,
        q_axis_region = next_q_region,
        d_axis_factor = elp[29],
        q_axis_factor = elp[30],
        refactorized,
    )
end

function _synchronous_machine_electrical_companion!(
    elp::Vector{Float64},
    cu::Vector{Float64},
    electrical_speed_rad_s::Float64,
    half_step_s::Float64,
    damping_control::Float64,
    saturation_enabled::Bool=false,
)
    factor = inv(electrical_speed_rad_s * half_step_s)
    damping_base = damping_control < 2.0 ? 100.0 : damping_control
    damping_inverse = inv(damping_base)
    damping_normalizer = inv(1.0 + damping_inverse)
    damping_ratio = 1.0 - 2.0 * damping_normalizer * damping_inverse
    factor *= damping_normalizer
    leakage = elp[19]
    elp[1] -= leakage
    elp[9] -= leakage
    z = zeros(Float64, 68)
    z[65] = elp[2] * factor
    z[66] = elp[4] * factor
    elp[5] *= factor
    z[67] = elp[10] * factor
    z[68] = elp[12] * factor
    elp[13] *= factor
    ratio = elp[21]
    z[61] = z[65] * ratio
    z[62] = z[66] * ratio
    z[63] = z[67] * ratio
    z[64] = z[68] * ratio
    elp[3] = elp[3] * factor - z[61]
    elp[6] = elp[6] * factor - z[62]
    elp[11] = elp[11] * factor - z[63]
    elp[14] = elp[14] * factor - z[64]
    zero_reactance = elp[17] * factor
    elp[17] = zero_reactance + elp[18]
    elp[18] = zero_reactance - elp[18] * damping_ratio
    elp[21] = inv(ratio)
    d_to_q_ratio = elp[9] / elp[1]
    z[50] = (elp[1] + leakage) * factor
    z[51] = (elp[9] + leakage) * factor
    z[52] = elp[3] + z[61]
    z[53] = elp[6] + z[62]
    z[54] = elp[11] + z[63]
    z[55] = elp[14] + z[64]
    z[1] = z[50] + elp[20]
    z[8] = z[51] + elp[20]
    z[15] = z[52] + elp[7]
    z[22] = z[53] + elp[8]
    z[29] = z[54] + elp[15]
    z[36] = z[55] + elp[16]
    z[11] = z[67]
    z[26] = z[11]
    z[12] = z[68]
    z[32] = z[12]
    z[30] = elp[13]
    z[35] = z[30]
    z[3] = z[65]
    z[13] = z[3]
    z[4] = z[66]
    z[19] = z[4]
    z[16] = elp[5]
    z[21] = z[16]
    for index in (1, 2, 4, 9, 10, 12)
        elp[index] /= electrical_speed_rad_s
    end
    elp[19] = leakage / electrical_speed_rad_s

    reduced = copy(z[1:36])
    _synchronous_machine_reduce_matrix!(reduced, 6, 2)
    d_axis_reduced = reduced[1]
    q_axis_reduced = reduced[8]
    reduced[8] -= d_axis_reduced
    zero_sequence = elp[17]
    elp[17] = inv(zero_sequence)
    branch_mutual = (zero_sequence - d_axis_reduced) / 3.0
    branch_direct = (zero_sequence + 2.0 * d_axis_reduced) / 3.0
    phase_positive_inverse = inv(d_axis_reduced)
    zero_inverse = inv(zero_sequence)
    elp[27] = (zero_inverse + 2.0 * phase_positive_inverse) / 3.0
    elp[28] = (zero_inverse - phase_positive_inverse) / 3.0
    elp[29] = 1.0
    elp[30] = 1.0
    elp[80] = phase_positive_inverse / sqrt(3.0 / 2.0)
    elp[81] = zero_inverse / sqrt(3.0)
    elp[31] = d_to_q_ratio
    for (destination, source) in zip(32:43, (3, 4, 11, 12, 15, 21, 16, 22, 29, 35, 30, 36))
        elp[destination] = reduced[source]
    end
    elp[70] = reduced[13]
    elp[71] = reduced[19]
    elp[72] = reduced[26]
    elp[73] = reduced[32]
    elp[74] = inv(d_axis_reduced)
    elp[75] = (leakage / electrical_speed_rad_s) / d_axis_reduced
    elp[76] = q_axis_reduced
    elp[79] = reduced[8]
    elp[44] = z[50] - elp[20] * damping_ratio
    elp[45] = z[13]
    elp[46] = z[19]
    elp[47] = z[51] - elp[20] * damping_ratio
    elp[48] = z[26]
    elp[49] = z[32]
    elp[50] = z[3]
    elp[51] = z[52] - elp[7] * damping_ratio
    elp[52] = z[21]
    elp[53] = z[4]
    elp[54] = z[16]
    elp[55] = z[53] - elp[8] * damping_ratio
    elp[56] = z[11]
    elp[57] = z[54] - elp[15] * damping_ratio
    elp[58] = z[35]
    elp[59] = z[12]
    elp[60] = z[30]
    elp[61] = z[55] - elp[16] * damping_ratio
    for offset in 35:42
        elp[27 + offset] = z[offset + 26]
    end
    saturation = if saturation_enabled
        _synchronous_machine_saturation_update!(
            elp;
            d_axis_current = cu[1],
            q_axis_current = cu[2],
            d_axis_rotor_current_1 = cu[4],
            d_axis_rotor_current_2 = cu[5],
            q_axis_rotor_current_1 = cu[6],
            q_axis_rotor_current_2 = cu[7],
            d_axis_region = 0,
            q_axis_region = 0,
            damping_ratio,
        )
    else
        (
            flux_magnitude = 0.0,
            d_axis_region = 0,
            q_axis_region = 0,
            d_axis_factor = 1.0,
            q_axis_factor = 1.0,
            refactorized = false,
        )
    end
    phase_positive_inverse = elp[74]
    saturation_enabled && (elp[75] = elp[19] * phase_positive_inverse)
    d_axis_reduced = inv(phase_positive_inverse)
    zero_sequence = inv(elp[17])
    branch_mutual = (zero_sequence - d_axis_reduced) / 3.0
    branch_direct = (zero_sequence + 2.0 * d_axis_reduced) / 3.0
    d_flux = (elp[1] * cu[1] + elp[2] * cu[4] + elp[4] * cu[5]) * elp[29] +
             elp[19] * cu[1]
    q_flux = (elp[9] * cu[2] + elp[10] * cu[6] + elp[12] * cu[7]) * elp[30] +
             elp[19] * cu[2]
    cu[15] = d_flux
    cu[16] = d_flux
    cu[17] = q_flux
    cu[18] = q_flux
    return (; damping_ratio, branch_direct, branch_mutual,
            d_flux, q_flux, phase_positive_inverse, zero_inverse,
            saturation_enabled,
            d_axis_saturation_region = saturation.d_axis_region,
            q_axis_saturation_region = saturation.q_axis_region,
            saturation_refactorized = saturation.refactorized)
end

function _synchronous_machine_mechanical_companion!(
    histq::Vector{Float64},
    shp::Vector{Float64},
    numask::Int,
    electrical_speed_rad_s::Float64,
    pole_pair_count::Float64,
    half_step_s::Float64,
)
    num2 = 2 * numask
    num4 = 4 * numask
    num8 = 8 * numask
    fill!(@view(shp[1:num2]), 0.0)
    fill!(@view(shp[(num8 + 1):(num8 + num4)]), 0.0)
    l2 = 0
    l4 = num8
    matrix4_start = l4
    n21 = num2 + numask
    n22 = n21 + num2
    n24 = n21 + num4
    if numask > 1
        n4 = n22 - 1
        for k in 2:numask
            l4 += 4
            l2 += 2
            n3 = n4 + k
            coupling = -(shp[n3] * half_step_s + shp[n3 - numask])
            shp[l2] = coupling
            shp[l2 + 1] = -coupling
            shp[l4 - 2] = -2.0 * shp[n3]
            shp[l4] = coupling
            shp[l4 + 1] = -shp[l4 - 2]
            shp[l4 + 3] = -coupling
        end
        l2 = 0
        l4 = matrix4_start
        for _ in 1:(numask - 1)
            l2 += 2
            l4 += 4
            shp[l2 - 1] -= shp[l2]
            shp[l4 - 3] -= shp[l4 - 2]
            shp[l4 - 1] -= shp[l4]
        end
    end
    l2 = 0
    l4 = matrix4_start
    for k in 1:numask
        n24 += 1
        local_n22 = n21 + k
        n23 = n24 - numask
        l2 += 2
        l4 += 4
        shp[l2 - 1] += shp[n23] + shp[n24] + shp[local_n22] / half_step_s
        shp[l4 - 1] += shp[n23] + shp[n24] - shp[local_n22] / half_step_s
    end
    scale = electrical_speed_rad_s / pole_pair_count
    for k in 1:numask
        histq[numask + num4 + k] *= scale
        shp[num2 + k] *= 2.0 * scale
    end
    return shp
end

function _synchronous_machine_sequence_components(values::Vector{ComplexF64})
    length(values) == 3 || throw(ArgumentError("three phase phasors are required"))
    ar, ai = reim(values[1])
    br, bi = reim(values[2])
    cr, ci = reim(values[3])
    half_sum_real = (br + cr) / 2.0
    quadrature_imag = (ci - bi) * sqrt(3.0) / 2.0
    half_sum_imag = (ci + bi) / 2.0
    quadrature_real = (br - cr) * sqrt(3.0) / 2.0
    positive = complex(
        (ar - half_sum_real + quadrature_imag) / 3.0,
        (ai - half_sum_imag + quadrature_real) / 3.0,
    )
    zero = sum(values) / 3.0
    negative = complex(
        (ar - half_sum_real - quadrature_imag) / 3.0,
        (ai - half_sum_imag - quadrature_real) / 3.0,
    )
    return (; positive, negative, zero)
end

function _synchronous_machine_steady_state_initialization(
    parameters::SynchronousMachineParameters,
    physical_elp::Vector{Float64},
    raw_shp::Vector{Float64},
    terminal_voltage_phasors::Vector{ComplexF64},
    phase_current_phasors::Vector{ComplexF64},
    electrical_speed_rad_s::Float64,
    initialization_tolerance::Float64,
    max_iterations::Int,
)
    numask = length(parameters.rotor_masses)
    pole_pairs = physical_elp[26]
    current_sequence = _synchronous_machine_sequence_components(phase_current_phasors)
    voltage_phasors = copy(terminal_voltage_phasors)
    if parameters.delta_connected
        voltage_phasors = ComplexF64[
            terminal_voltage_phasors[1] - terminal_voltage_phasors[2],
            terminal_voltage_phasors[2] - terminal_voltage_phasors[3],
            terminal_voltage_phasors[3] - terminal_voltage_phasors[1],
        ]
    end
    winding_current_phasors = parameters.delta_connected ? ComplexF64[
        (phase_current_phasors[1] - phase_current_phasors[2]) / 3.0,
        (phase_current_phasors[2] - phase_current_phasors[3]) / 3.0,
        (phase_current_phasors[3] - phase_current_phasors[1]) / 3.0,
    ] : phase_current_phasors
    voltage_sequence = _synchronous_machine_sequence_components(voltage_phasors)
    current_gain = parameters.delta_connected ? inv(sqrt(3.0)) : 1.0
    current_angle_shift = parameters.delta_connected ? pi / 6.0 : 0.0
    positive_current = current_sequence.positive * current_gain * cis(current_angle_shift)
    current_magnitude = abs(positive_current)
    current_angle = angle(positive_current)
    terminal_voltage = voltage_sequence.positive
    voltage_magnitude = abs(terminal_voltage)
    voltage_angle = angle(terminal_voltage)
    scaled_voltage_magnitude = voltage_magnitude * sqrt(3.0 / 2.0)

    es = physical_elp[1]
    ds = physical_elp[2]
    cs = physical_elp[9]
    leakage = physical_elp[19]
    resistance = physical_elp[20]
    d_main = es - leakage
    q_main = cs - leakage
    q_to_d_ratio = q_main / d_main
    rotor_angle = angle(terminal_voltage +
                        complex(resistance, cs) * positive_current)
    scaled_current_magnitude = current_magnitude * sqrt(3.0 / 2.0)
    relative_current_angle = current_angle - rotor_angle
    d_current = scaled_current_magnitude * sin(relative_current_angle)
    q_current = scaled_current_magnitude * cos(relative_current_angle)
    relative_voltage_angle = voltage_angle - rotor_angle
    q_voltage = scaled_voltage_magnitude * cos(relative_voltage_angle)
    field_current = (q_voltage + resistance * q_current - es * d_current) / ds
    negative_ratio = abs(current_sequence.negative) /
                     max(current_magnitude, eps(Float64))
    negative_ratio <= 10.0 * initialization_tolerance ||
        throw(ArgumentError("unbalanced synchronous-machine initialization requires negative-sequence rotor correction"))

    d_saturation_factor = 1.0
    q_saturation_factor = 1.0
    effective_d_synchronous = es
    effective_q_synchronous = cs
    effective_d_mutual = ds
    if parameters.saturation_base_current < 0.0
        d_threshold = physical_elp[22]
        d_slope = physical_elp[23]
        q_threshold = physical_elp[24]
        q_slope = physical_elp[25]
        armature_to_field_ratio = physical_elp[21]
        max_iterations > 0 || throw(ArgumentError(
            "synchronous-machine saturation initialization requires a positive iteration limit",
        ))

        function operating_point_from_saturation_factors(
            d_factor::Float64,
            q_factor::Float64,
        )
            isfinite(d_factor) && d_factor > 0.0 || throw(ArgumentError(
                "synchronous-machine d-axis saturation factor must be finite and positive",
            ))
            isfinite(q_factor) && q_factor > 0.0 || throw(ArgumentError(
                "synchronous-machine q-axis saturation factor must be finite and positive",
            ))
            d_mutual = ds * d_factor
            d_synchronous = d_main * d_factor + leakage
            q_synchronous = q_main * q_factor + leakage
            d_mutual != 0.0 || throw(ArgumentError(
                "synchronous-machine saturated d-axis mutual reactance is zero",
            ))
            angle_value = angle(
                terminal_voltage +
                complex(resistance, q_synchronous) * positive_current,
            )
            relative_current = current_angle - angle_value
            d_current_value =
                scaled_current_magnitude * sin(relative_current)
            q_current_value =
                scaled_current_magnitude * cos(relative_current)
            relative_voltage = voltage_angle - angle_value
            q_voltage_value =
                scaled_voltage_magnitude * cos(relative_voltage)
            field_current_value = (
                q_voltage_value + resistance * q_current_value -
                d_synchronous * d_current_value
            ) / d_mutual
            flux_value = hypot(
                d_current_value / armature_to_field_ratio +
                field_current_value,
                (q_current_value / armature_to_field_ratio) * q_to_d_ratio,
            )
            all(isfinite, (
                angle_value,
                d_current_value,
                q_current_value,
                field_current_value,
                flux_value,
            )) || throw(ArgumentError(
                "synchronous-machine saturated operating point is nonfinite",
            ))
            return (;
                d_saturation_factor=d_factor,
                q_saturation_factor=q_factor,
                effective_d_mutual=d_mutual,
                effective_d_synchronous=d_synchronous,
                effective_q_synchronous=q_synchronous,
                rotor_angle=angle_value,
                d_current=d_current_value,
                q_current=q_current_value,
                field_current=field_current_value,
                flux_magnitude=flux_value,
            )
        end

        function saturated_operating_point(flux_value::Float64)
            isfinite(flux_value) && flux_value >= 0.0 || throw(ArgumentError(
                "synchronous-machine saturation flux must be finite and nonnegative",
            ))
            return operating_point_from_saturation_factors(
                _synchronous_machine_saturation_factor(
                    flux_value,
                    d_threshold,
                    d_slope,
                ),
                _synchronous_machine_saturation_factor(
                    flux_value,
                    q_threshold,
                    q_slope,
                ),
            )
        end

        initial_flux_magnitude = hypot(
            d_current / armature_to_field_ratio + field_current,
            (q_current / armature_to_field_ratio) * q_to_d_ratio,
        )
        predictor = operating_point_from_saturation_factors(
            initial_flux_magnitude <= d_threshold ? 1.0 :
                (1.0 - d_slope * initial_flux_magnitude) /
                (1.0 - d_slope * d_threshold),
            initial_flux_magnitude <= q_threshold ? 1.0 :
                (1.0 - q_slope * initial_flux_magnitude) /
                (1.0 - q_slope * q_threshold),
        )
        trial_flux_magnitude = predictor.flux_magnitude
        converged_point = nothing
        last_residual = Inf
        for _ in 1:max_iterations
            point = saturated_operating_point(trial_flux_magnitude)
            residual = point.flux_magnitude - trial_flux_magnitude
            flux_scale = max(
                abs(point.flux_magnitude),
                abs(trial_flux_magnitude),
                d_threshold,
                q_threshold,
                1.0,
            )
            last_residual = residual
            if abs(residual) <= initialization_tolerance * flux_scale
                converged_point = point
                break
            end

            derivative_interval = sqrt(eps(Float64)) * flux_scale
            perturbed_flux = trial_flux_magnitude + derivative_interval
            perturbed_point = saturated_operating_point(perturbed_flux)
            perturbed_residual = perturbed_point.flux_magnitude - perturbed_flux
            residual_derivative =
                (perturbed_residual - residual) / derivative_interval
            newton_flux = if isfinite(residual_derivative) &&
                             abs(residual_derivative) > sqrt(eps(Float64))
                trial_flux_magnitude - residual / residual_derivative
            else
                point.flux_magnitude
            end
            if !isfinite(newton_flux) || newton_flux < 0.0
                newton_flux = point.flux_magnitude
            end
            newton_point = saturated_operating_point(newton_flux)
            newton_residual = newton_point.flux_magnitude - newton_flux
            newton_flux_scale = max(
                abs(newton_point.flux_magnitude),
                abs(newton_flux),
                d_threshold,
                q_threshold,
                1.0,
            )
            if abs(newton_residual) <=
               initialization_tolerance * newton_flux_scale
                last_residual = newton_residual
                converged_point = newton_point
                break
            end
            trial_flux_magnitude =
                abs(newton_residual) < abs(residual) ?
                newton_flux : point.flux_magnitude
        end
        converged_point === nothing && throw(ArgumentError(
            "synchronous-machine saturation initialization did not converge; " *
            "flux residual $(last_residual)",
        ))
        d_saturation_factor = converged_point.d_saturation_factor
        q_saturation_factor = converged_point.q_saturation_factor
        effective_d_mutual = converged_point.effective_d_mutual
        effective_d_synchronous = converged_point.effective_d_synchronous
        effective_q_synchronous = converged_point.effective_q_synchronous
        rotor_angle = converged_point.rotor_angle
        d_current = converged_point.d_current
        q_current = converged_point.q_current
        field_current = converged_point.field_current
    end

    zero_current = real(current_sequence.zero) * sqrt(3.0)
    field_voltage = -physical_elp[7] * field_current
    cu = zeros(Float64, 24)
    cu[1:7] .= (d_current, q_current, zero_current, field_current, 0.0, 0.0, 0.0)
    cu[8:10] .= (d_current, q_current, zero_current)
    cu[11] = field_voltage
    phase_currents = real.(winding_current_phasors)
    cu[12:14] .= phase_currents

    mechanical_angle = (rotor_angle + pi / 2.0) / pole_pairs
    mechanical_speed = electrical_speed_rad_s / pole_pairs
    d_flux_voltage = (
        effective_d_synchronous * d_current + effective_d_mutual * field_current
    ) / electrical_speed_rad_s
    q_flux_voltage = effective_q_synchronous * q_current / electrical_speed_rad_s
    generator_torque =
        (d_flux_voltage * q_current - q_flux_voltage * d_current) * pole_pairs * 1.0e-6
    exciter_torque =
        -(field_voltage * pole_pairs * 1.0e-6 * field_current) / mechanical_speed
    dc_generator_torque =
        ((effective_d_synchronous - effective_q_synchronous) * d_current +
         effective_d_mutual * field_current) * q_current *
        1.0e-6 / mechanical_speed
    total_external_torque = dc_generator_torque +
                            (parameters.exciter_mass_index == 0 ? 0.0 : exciter_torque)

    torque_fractions = raw_shp[(2 * numask + 1):(3 * numask)]
    applied_torques = torque_fractions .* total_external_torque
    applied_torques[parameters.generator_mass_index] -= dc_generator_torque
    parameters.exciter_mass_index != 0 &&
        (applied_torques[parameters.exciter_mass_index] -= exciter_torque)
    shaft_torques = cumsum(applied_torques)[1:max(numask - 1, 0)]
    angles = fill(mechanical_angle, numask)
    for index in (parameters.generator_mass_index - 1):-1:1
        stiffness = raw_shp[2 * numask + 3 * numask + index]
        stiffness != 0.0 || throw(ArgumentError("zero shaft stiffness left of generator mass"))
        angles[index] = angles[index + 1] + shaft_torques[index] / stiffness
    end
    for index in parameters.generator_mass_index:(numask - 1)
        stiffness = raw_shp[2 * numask + 3 * numask + index]
        stiffness != 0.0 || throw(ArgumentError("zero shaft stiffness right of generator mass"))
        angles[index + 1] = angles[index] - shaft_torques[index] / stiffness
    end
    histq = zeros(Float64, 24 * numask)
    histq[1:numask] .= angles
    histq[(numask + 1):(2 * numask)] .= mechanical_speed
    histq[(4 * numask + 1):(4 * numask + length(shaft_torques))] .= shaft_torques
    histq[(5 * numask + 1):(6 * numask)] .= torque_fractions .* total_external_torque
    return (; cu, histq, d_current, q_current, zero_current, field_current,
            field_voltage, d_flux_voltage, q_flux_voltage, generator_torque,
            exciter_torque, dc_generator_torque, phase_currents,
            phase_voltage_phasors = voltage_phasors)
end

function synchronous_machine_initial_state(
    parameters::SynchronousMachineParameters,
    terminal_voltage_phasors::AbstractVector{<:Complex},
    phase_current_phasors::AbstractVector{<:Complex};
    time_step_s::Real,
    frequency_hz::Real,
    damping_control::Real=100.0,
    initialization_tolerance::Real=1.0e-9,
    omega_tolerance::Real=1.0e-4,
    speed_tolerance::Real=1.0e-5,
    max_iterations::Int=10,
    speed_floor::Real=1.0e-12,
    terminal_network_admittance=nothing,
    initial_mechanical_angle_rad::Union{Nothing,Real}=nothing,
)
    length(terminal_voltage_phasors) == 3 && length(phase_current_phasors) == 3 ||
        throw(ArgumentError("synchronous-machine initialization requires three voltage and current phasors"))
    dt_s = Float64(time_step_s)
    dt_s > 0.0 || throw(ArgumentError("time_step_s must be positive"))
    frequency = Float64(frequency_hz)
    frequency > 0.0 || throw(ArgumentError("frequency_hz must be positive"))
    omega = 2.0 * pi * frequency
    half_step = dt_s / 2.0
    physical = _synchronous_machine_physical_storage(parameters, omega)
    initial_mechanical_angle = initial_mechanical_angle_rad === nothing ? nothing :
        Float64(initial_mechanical_angle_rad)
    initial_mechanical_angle === nothing || isfinite(initial_mechanical_angle) ||
        throw(ArgumentError("initial_mechanical_angle_rad must be finite"))
    steady = _synchronous_machine_steady_state_initialization(
        parameters,
        physical.elp,
        physical.shp,
        ComplexF64.(terminal_voltage_phasors),
        ComplexF64.(phase_current_phasors),
        omega,
        Float64(initialization_tolerance),
        max_iterations,
    )
    elp = copy(physical.elp)
    shp = copy(physical.shp)
    histq = copy(steady.histq)
    cu = copy(steady.cu)
    initial_mechanical_angle === nothing ||
        (histq[1:length(parameters.rotor_masses)] .= initial_mechanical_angle)
    companion = _synchronous_machine_electrical_companion!(
        elp,
        cu,
        omega,
        half_step,
        Float64(damping_control),
        physical.saturation_enabled,
    )
    if terminal_network_admittance !== nothing
        size(terminal_network_admittance) == (3, 3) ||
            throw(ArgumentError("terminal_network_admittance must be a 3 x 3 matrix"))
        for col in 1:3, row in 1:3
            elp[81 + (col - 1) * 3 + row] =
                Float64(terminal_network_admittance[row, col])
        end
    end
    pole_pairs = elp[26]
    cu[19] = histq[parameters.generator_mass_index] * pole_pairs
    cu[20] = cos(cu[19])
    cu[21] = sin(cu[19])
    electrical_increment = omega * dt_s
    cu[22] = cu[19] - electrical_increment
    cu[23] = cu[22] - electrical_increment
    cu[24] = omega
    _synchronous_machine_mechanical_companion!(
        histq,
        shp,
        length(parameters.rotor_masses),
        omega,
        pole_pairs,
        half_step,
    )

    past_state = past_machine_history_state(histq = histq, shp = shp, elp = elp)
    past = past_machine_history_update!(
        past_state;
        stator = (
            voltages = real.(steady.phase_voltage_phasors),
            current_terms = steady.phase_currents,
            a1 = elp[27],
            a2 = elp[28],
        ),
        rotor = (
            d_stator_current = steady.d_current,
            q_stator_current = steady.q_current,
            d_rotor_currents = (steady.field_current, 0.0),
            q_rotor_currents = (0.0, 0.0),
            d_gains = Tuple(elp[50:55]),
            q_gains = Tuple(elp[56:61]),
            damping_gains = (elp[7], elp[8], elp[15], elp[16]),
            damrat = companion.damping_ratio,
        ),
        internal_stator = (
            d_current = steady.d_current,
            q_current = steady.q_current,
            zero_sequence_history = steady.zero_current,
            zero_sequence_voltage = sum(real.(steady.phase_voltage_phasors)) / sqrt(3.0),
            zero_sequence_gain = elp[18],
            d_axis_history = cu[9],
            q_axis_history = cu[8],
            d_axis_gain = elp[66],
            q_axis_gain = elp[68],
            field_d = cu[15],
            field_q = cu[17],
            omega = omega,
        ),
        mechanical_backoff = (
            numask = length(parameters.rotor_masses),
            nlocg = parameters.generator_mass_index,
            nloce = parameters.exciter_mass_index,
            a1_current = steady.d_current,
            a2_current = steady.q_current,
            a3_current = steady.field_current,
            afd = cu[15] * omega,
            afq = cu[17] * omega,
            coupling_gain = pole_pairs,
            omega = omega,
            delta2 = half_step,
            exciter_input = steady.field_voltage,
        ),
        mechanical_matrix = (numask = length(parameters.rotor_masses),),
    )
    cu[1:3] .= past.internal_stator_history[1:3]
    cu[4:7] .= past.rotor_history
    cu[12] = 1.0
    cu[13] = steady.d_current
    cu[14] = steady.q_current
    equation_state = SynchronousMachineEquationState(cu, past.histq, past.shp)
    current_reduction_factor = past.elp[21]
    state = SynchronousMachineDynamicState(
        equation_state,
        past.stator_history,
        past.elp,
        applied_torque_distribution =
            Float64[mass.torque_fraction for mass in parameters.rotor_masses],
        exciter_port_state =
            SynchronousMachineExciterPortState(current_reduction_factor),
        saturation_enabled = physical.saturation_enabled,
        d_axis_saturation_region = companion.d_axis_saturation_region,
        q_axis_saturation_region = companion.q_axis_saturation_region,
    )
    return SynchronousMachineInitialization(
        state,
        ComplexF64.(phase_current_phasors),
        length(parameters.rotor_masses),
        parameters.generator_mass_index,
        parameters.exciter_mass_index,
        inv(sqrt(3.0 / 2.0)),
        Float64(speed_tolerance),
        Float64(omega_tolerance),
        Float64(speed_floor),
        max_iterations,
        companion.damping_ratio,
        6.0 * dt_s,
        omega / 2.0,
        omega,
        electrical_increment,
    )
end

function SynchronousMachineTACSInterfaceState(etac_values::AbstractVector{<:Real}; lmset::Int=0)
    lmset >= 0 || throw(ArgumentError("lmset must be nonnegative"))
    return SynchronousMachineTACSInterfaceState(
        etac_values = Float64[Float64(value) for value in etac_values],
        lmset = lmset,
        transfer_count = 0,
        transfer_mutated = false,
    )
end

function SynchronousMachineRotorMassTACSState(
    machine_output_table::AbstractVector{<:Real},
    tacs_values::AbstractVector{<:Real};
    transfer_pass_marker::Int=0,
)
    return SynchronousMachineRotorMassTACSState(
        machine_output_table = _machine_float_vector(machine_output_table),
        tacs_values = _machine_float_vector(tacs_values),
        transfer_pass_marker = transfer_pass_marker,
        mass_history_initialized = false,
        transfer_mutated = false,
    )
end

function UniversalMachineRotorCurrentState(current_values::AbstractVector{<:Real})
    currents = _machine_real_vector("current_values", current_values, 3)
    return UniversalMachineRotorCurrentState(
        current_values = currents,
        solution_matrix = zeros(Float64, 3, 3),
        right_hand_side = zeros(Float64, 3),
        rotor_current_solution_mutated = false,
    )
end

function UniversalMachineStatorExcitationCurrentState(current_values::AbstractVector{<:Real})
    currents = _machine_float_vector(current_values)
    all(isfinite, currents) || throw(ArgumentError("current_values entries must be finite"))
    return UniversalMachineStatorExcitationCurrentState(
        current_values = currents,
        solution_matrix = zeros(Float64, length(currents), length(currents)),
        right_hand_side = zeros(Float64, length(currents)),
        axis_kinds = Symbol[],
        stator_excitation_solution_mutated = false,
    )
end

function UniversalMachineDirectCouplingState(current_values::AbstractVector{<:Real})
    currents = _machine_float_vector(current_values)
    all(isfinite, currents) || throw(ArgumentError("current_values entries must be finite"))
    return UniversalMachineDirectCouplingState(
        current_values = currents,
        coupling_fraction = 0.0,
        voltage_coupling_fraction = 0.0,
        resistance_coupling_fraction = 0.0,
        direct_coupling_mutated = false,
    )
end

function UniversalMachineHistoryCurrentState(history_currents::AbstractVector{<:Real})
    currents = _machine_float_vector(history_currents)
    all(isfinite, currents) || throw(ArgumentError("history_currents entries must be finite"))
    return UniversalMachineHistoryCurrentState(
        history_currents = currents,
        rotor_history_matrix = zeros(Float64, 3, 3),
        stator_thevenin_drop = zeros(Float64, max(length(currents) - 3, 0)),
        axis_kinds = Symbol[],
        history_current_update_mutated = false,
    )
end

function UniversalMachineMechanicalIterationState(
    mechanical_speed_rad_s::Real,
    previous_mechanical_speed_rad_s::Real,
    mechanical_angle_rad::Real,
)
    speed = Float64(mechanical_speed_rad_s)
    previous_speed = Float64(previous_mechanical_speed_rad_s)
    angle = Float64(mechanical_angle_rad)
    all(isfinite, (speed, previous_speed, angle)) ||
        throw(ArgumentError("mechanical iteration state values must be finite"))
    return UniversalMachineMechanicalIterationState(
        mechanical_speed_rad_s = speed,
        previous_mechanical_speed_rad_s = previous_speed,
        mechanical_angle_rad = angle,
        predicted_speed_rad_s = 0.0,
        solved_speed_rad_s = speed,
        generated_torque = 0.0,
        torque_increment = 0.0,
        iteration_count = 0,
        converged = false,
        mechanical_iteration_mutated = false,
    )
end

function UniversalMachineFluxSaturationState(d_axis_flux::Real, q_axis_flux::Real)
    d_flux = Float64(d_axis_flux)
    q_flux = Float64(q_axis_flux)
    all(isfinite, (d_flux, q_flux)) || throw(ArgumentError("flux values must be finite"))
    return UniversalMachineFluxSaturationState(
        d_axis_flux = d_flux,
        q_axis_flux = q_flux,
        d_axis_current = 0.0,
        q_axis_current = 0.0,
        d_axis_saturation_offset = 0.0,
        q_axis_saturation_offset = 0.0,
        flux_saturation_mutated = false,
    )
end

function UniversalMachinePostsolveState(
    coil_parameters::AbstractVector{<:Real},
    source_crests::AbstractVector{<:Real},
    current_values::AbstractVector{<:Real},
    prediction_values::AbstractVector{<:Real},
    history_values::AbstractVector{<:Real};
    d_axis_flux::Real=0.0,
    q_axis_flux::Real=0.0,
    theta_electric_rad::Real=0.0,
    machine_type::Int=0,
    input_mode::Int=0,
    prediction_loop_marker::Int=0,
)
    return UniversalMachinePostsolveState(
        coil_parameters = _machine_float_vector(coil_parameters),
        source_crests = _machine_float_vector(source_crests),
        current_values = _machine_float_vector(current_values),
        prediction_values = _machine_float_vector(prediction_values),
        history_values = _machine_float_vector(history_values),
        prediction_report_text_lines = String[],
        d_axis_flux = Float64(d_axis_flux),
        q_axis_flux = Float64(q_axis_flux),
        theta_electric_rad = Float64(theta_electric_rad),
        machine_type = machine_type,
        input_mode = input_mode,
        prediction_loop_marker = prediction_loop_marker,
        postsolve_update_mutated = false,
        prediction_report_text_mutated = false,
    )
end

function UniversalMachineReportScheduleState(last_report_step::Int)
    return UniversalMachineReportScheduleState(
        last_report_step = last_report_step,
        report_text_lines = String[],
        report_schedule_mutated = false,
        report_text_mutated = false,
    )
end

function _machine_float_vector(values)
    return Float64[Float64(value) for value in values]
end

function _machine_check_index(index::Int, values, name::String)
    1 <= index <= length(values) || throw(ArgumentError("$name index $index is outside supplied storage"))
    return index
end

function _universal_machine_axis_current_totals(
    currents::AbstractVector{Float64},
    d_axis_stator_current_indices,
    q_axis_stator_current_indices,
    d_axis_subtract_current_index::Int,
)
    length(currents) >= 3 ||
        throw(ArgumentError("universal-machine axis current totals require at least three currents"))
    d_current = currents[2]
    q_current = currents[3]
    for index in d_axis_stator_current_indices
        d_current += currents[_machine_check_index(Int(index), currents, "d-axis stator current")]
    end
    for index in q_axis_stator_current_indices
        q_current += currents[_machine_check_index(Int(index), currents, "q-axis stator current")]
    end
    if d_axis_subtract_current_index != 0
        d_current -= currents[_machine_check_index(d_axis_subtract_current_index, currents, "d-axis subtract current")]
    end
    return (; d_axis_current = d_current, q_axis_current = q_current)
end

function _machine_bansol_segment!(
    ab::AbstractVector{Float64},
    ab_offset::Int,
    x::AbstractVector{Float64},
    x_offset::Int,
    n::Int,
)
    n > 0 || throw(ArgumentError("BANSOL dimension must be positive"))
    ab_offset >= 0 || throw(ArgumentError("BANSOL matrix offset must be nonnegative"))
    x_offset >= 0 || throw(ArgumentError("BANSOL vector offset must be nonnegative"))
    length(ab) >= ab_offset + 2 * n ||
        throw(ArgumentError("BANSOL matrix storage is shorter than 2*numask entries"))
    length(x) >= x_offset + n ||
        throw(ArgumentError("BANSOL vector storage is shorter than numask entries"))

    i2 = 2
    i1 = 1
    d = x[x_offset + 1]
    while true
        x[x_offset + i1] = d * ab[ab_offset + i2 - 1]
        i1 == n && break
        i1 += 1
        d = x[x_offset + i1] - d * ab[ab_offset + i2]
        i2 += 2
    end
    while i1 != 1
        i2 -= 2
        i1 -= 1
        x[x_offset + i1] = x[x_offset + i1] - x[x_offset + i1 + 1] * ab[ab_offset + i2]
    end
    return x
end

function _synchronous_machine_equation_check_storage(
    cu::AbstractVector{Float64},
    histq::AbstractVector{Float64},
    shp::AbstractVector{Float64},
    elp::AbstractVector{Float64},
    numask::Int,
    nlocg::Int,
    nloce::Int,
)
    numask > 0 || throw(ArgumentError("numask must be positive"))
    nlocg >= 0 || throw(ArgumentError("nlocg must be nonnegative"))
    nloce >= 0 || throw(ArgumentError("nloce must be nonnegative"))
    num2 = 2 * numask
    num4 = 2 * num2
    ikv = numask
    ikw = ikv + 1
    ikp = ikv + numask
    n26 = num4
    n27 = n26 + numask
    n22 = num2
    ksg = nlocg + numask
    ksex = ikv + nloce
    length(cu) >= 21 || throw(ArgumentError("cu_values must contain at least 21 entries"))
    length(histq) >= max(n27 + numask, ikp + num2, n22 + 6 * numask, ksg, ksex) ||
        throw(ArgumentError("histq_values is too short for the synchronous-machine equation step"))
    length(shp) >= n22 + 6 * numask ||
        throw(ArgumentError("shp_values is too short for the synchronous-machine equation step"))
    length(elp) >= 80 || throw(ArgumentError("electrical_coefficients must contain at least 80 ELP entries"))
    return (; num2, num4, ikv, ikw, ikp, n26, n27, n22, ksg, ksex)
end

function _machine_int_vector(name::String, values, expected_count::Int)
    length(values) == expected_count ||
        throw(ArgumentError("$name must contain $expected_count entries"))
    return Int[Int(value) for value in values]
end

function _machine_real_vector(name::String, values, expected_count::Int)
    length(values) == expected_count ||
        throw(ArgumentError("$name must contain $expected_count entries"))
    result = Float64[Float64(value) for value in values]
    all(isfinite, result) || throw(ArgumentError("$name entries must be finite"))
    return result
end

function _machine_real_matrix(name::String, values, expected_rows::Int, expected_cols::Int)
    size(values) == (expected_rows, expected_cols) ||
        throw(ArgumentError("$name must be a $expected_rows x $expected_cols matrix"))
    result = Matrix{Float64}(undef, expected_rows, expected_cols)
    for col in 1:expected_cols, row in 1:expected_rows
        value = Float64(values[row, col])
        isfinite(value) || throw(ArgumentError("$name entries must be finite"))
        result[row, col] = value
    end
    return result
end

function _universal_machine_three_phase_solution(
    matrix::Matrix{Float64},
    right_hand_side::Vector{Float64},
)
    augmented = zeros(Float64, 3, 4)
    for row in 1:3
        for col in 1:3
            augmented[row, col] = matrix[row, col]
        end
        augmented[row, 4] = right_hand_side[row]
    end

    for pivot in 1:3
        pivot_value = augmented[pivot, pivot]
        pivot_value != 0.0 ||
            throw(ArgumentError("rotor-current equation matrix has a zero pivot"))
        pivot == 3 && continue
        for row in (pivot + 1):3
            numerator = augmented[row, pivot]
            for col in 4:-1:(pivot + 1)
                augmented[row, col] -= numerator * augmented[pivot, col] / pivot_value
            end
            augmented[row, pivot] = 0.0
        end
    end

    solution = zeros(Float64, 3)
    for offset in 1:3
        row = 4 - offset
        pivot_value = augmented[row, row]
        pivot_value != 0.0 ||
            throw(ArgumentError("rotor-current equation matrix has a zero pivot"))
        solution[row] = augmented[row, 4] / pivot_value
        for col in (row + 1):3
            solution[row] -= augmented[row, col] * solution[col] / pivot_value
        end
    end
    return solution
end

function _universal_machine_active_coil_solution(
    matrix::Matrix{Float64},
    right_hand_side::Vector{Float64},
)
    count = length(right_hand_side)
    size(matrix) == (count, count) || throw(ArgumentError("active-coil matrix dimension mismatch"))
    count > 0 || throw(ArgumentError("at least one stator/excitation coil is required"))
    if count == 1
        matrix[1, 1] != 0.0 ||
            throw(ArgumentError("stator/excitation current equation matrix has a zero pivot"))
        return [right_hand_side[1] / matrix[1, 1]]
    elseif count == 2
        determinant = matrix[1, 1] * matrix[2, 2] - matrix[1, 2] * matrix[2, 1]
        determinant != 0.0 ||
            throw(ArgumentError("stator/excitation current equation matrix is singular"))
        return [
            (matrix[2, 2] * right_hand_side[1] - matrix[1, 2] * right_hand_side[2]) / determinant,
            (matrix[1, 1] * right_hand_side[2] - matrix[2, 1] * right_hand_side[1]) / determinant,
        ]
    elseif count == 3
        return _universal_machine_three_phase_solution(matrix, right_hand_side)
    end
    throw(ArgumentError("stator/excitation current solve supports one to three active coils"))
end

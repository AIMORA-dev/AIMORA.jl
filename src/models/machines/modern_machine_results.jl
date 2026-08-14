export ModernMachineResult,
       ModernMachineTrace,
       simulate_modern_machine,
       modern_machine_result,
       modern_machine_result_quantity,
       write_modern_machine_trace_csv,
       modern_machine_report_text

struct ModernMachineResult
    schema::Symbol
    machine::Symbol
    family::Symbol
    operating_mode::Symbol
    deterministic_signature_sha256::String
    status::Symbol
    available_quantities::Dict{Symbol,Any}
    unavailable_quantities::Dict{Symbol,String}
    diagnostics::NamedTuple
    limitations::Vector{String}
end

struct ModernMachineTrace
    time_s::Vector{Float64}
    phase_voltage_v::Matrix{Float64}
    phase_current_a::Matrix{Float64}
    zero_d_q_voltage_v::Matrix{Float64}
    zero_d_q_current_a::Matrix{Float64}
    electromagnetic_torque_nm::Vector{Float64}
    mechanical_speed_rad_s::Matrix{Float64}
    mechanical_angle_rad::Matrix{Float64}
    shaft_coupling_torque_nm::Matrix{Float64}
    field_voltage_v::Vector{Float64}
    rotor_port_power_w::Vector{Float64}
    terminal_power_w::Vector{Float64}
    magnetic_coenergy_j::Vector{Float64}
    kinetic_energy_j::Vector{Float64}
    elastic_energy_j::Vector{Float64}
    dissipated_energy_j::Vector{Float64}
    energy_residual_j::Vector{Float64}
    energy_quadrature_defect_j::Vector{Float64}
    event_count::Vector{Int}
    control_sample_count::Vector{Int}
    result::ModernMachineResult
end

function modern_machine_result_quantity(result::ModernMachineResult, key::Symbol)
    haskey(result.available_quantities, key) && return result.available_quantities[key]
    haskey(result.unavailable_quantities, key) || throw(KeyError(key))
    return (
        available=false,
        quantity=key,
        reason=result.unavailable_quantities[key],
    )
end

function modern_machine_result(runtime::ModernMachineRuntime)
    preparation = runtime.preparation
    specification = preparation.specification
    accepted = runtime.accepted_state
    layout = preparation.layout
    available = Dict{Symbol,Any}(
        :phase_terminal_voltage_v => copy(accepted.terminal_voltage_v),
        :phase_terminal_current_a => copy(accepted.terminal_current_a),
        :electrical_flux_wb => copy(accepted.flux_wb),
        :winding_current_a => copy(accepted.winding_current_a),
        :electromagnetic_torque_nm => accepted.electromagnetic_torque_nm,
        :shaft_angle_rad => copy(runtime.shaft_state.angle_rad),
        :shaft_speed_rad_s => copy(runtime.shaft_state.speed_rad_s),
        :shaft_coupling_torque_nm => copy(runtime.shaft_state.coupling_torque_nm),
        :terminal_power_w => accepted.terminal_power_w,
        :field_power_w => accepted.field_power_w,
        :rotor_port_power_w => accepted.rotor_port_power_w,
        :magnetic_coenergy_j => accepted.magnetic_coenergy_j,
        :kinetic_energy_j => runtime.shaft_state.kinetic_energy_j,
        :elastic_energy_j => runtime.shaft_state.elastic_energy_j,
        :dissipated_energy_j => accepted.dissipated_energy_j,
        :maximum_energy_residual_j => accepted.maximum_energy_residual_j,
        :maximum_energy_quadrature_defect_j =>
            accepted.maximum_energy_quadrature_defect_j,
        :field_voltage_v => runtime.control_state.field_voltage_v,
        :mechanical_torque_nm => runtime.control_state.mechanical_torque_nm,
        :control_sample_count => runtime.control_state.sample_count,
    )
    unavailable = Dict{Symbol,String}()
    if layout.field_index === nothing
        unavailable[:field_current_a] = "selected family has no field winding"
    else
        available[:field_current_a] = accepted.winding_current_a[layout.field_index]
    end
    if isempty(layout.rotor_d_indices)
        unavailable[:rotor_current_dq_a] = "selected family has no induction rotor circuit"
    else
        available[:rotor_current_dq_a] = vcat(
            permutedims(accepted.winding_current_a[layout.rotor_d_indices]),
            permutedims(accepted.winding_current_a[layout.rotor_q_indices]),
        )
    end
    if specification.family !== PermanentMagnetSynchronousMachine
        unavailable[:permanent_magnet_flux_wb] =
            "selected family has no permanent-magnet flux"
    else
        available[:permanent_magnet_flux_wb] =
            specification.electrical.permanent_magnet_flux_wb
    end
    return ModernMachineResult(
        :aimora_modern_machine_result_v1,
        specification.id,
        _MODERN_MACHINE_FAMILY_IDS[specification.family],
        _MODERN_MACHINE_MODE_IDS[specification.operating_mode],
        preparation.deterministic_signature_sha256,
        :accepted,
        available,
        unavailable,
        modern_machine_runtime_diagnostics(runtime),
        String[
            "generic fixed-step phase-domain machine; not a vendor or standard model",
            "thermal, hysteresis, spatial-field, arbitrary internal-fault, protection, and certification behavior are outside this result",
        ],
    )
end

function _machine_trace_sample!(
    runtime::ModernMachineRuntime,
    sample::Int,
    time_s::Vector{Float64},
    phase_voltage_v::Matrix{Float64},
    phase_current_a::Matrix{Float64},
    zero_d_q_voltage_v::Matrix{Float64},
    zero_d_q_current_a::Matrix{Float64},
    electromagnetic_torque_nm::Vector{Float64},
    mechanical_speed_rad_s::Matrix{Float64},
    mechanical_angle_rad::Matrix{Float64},
    shaft_coupling_torque_nm::Matrix{Float64},
    field_voltage_v::Vector{Float64},
    rotor_port_power_w::Vector{Float64},
    terminal_power_w::Vector{Float64},
    magnetic_coenergy_j::Vector{Float64},
    kinetic_energy_j::Vector{Float64},
    elastic_energy_j::Vector{Float64},
    dissipated_energy_j::Vector{Float64},
    energy_residual_j::Vector{Float64},
    energy_quadrature_defect_j::Vector{Float64},
    event_count::Vector{Int},
    control_sample_count::Vector{Int},
)
    accepted = runtime.accepted_state
    specification = runtime.preparation.specification
    mass_index = _machine_electromagnetic_mass_index(specification)
    electrical_angle = specification.pole_pairs * runtime.shaft_state.angle_rad[mass_index]
    phase_terminal = runtime.preparation.terminal_voltage_map * accepted.terminal_voltage_v
    phase_current = accepted.terminal_current_a[1:3]
    time_s[sample] = accepted.accepted_time_s
    phase_voltage_v[:, sample] .= phase_terminal
    phase_current_a[:, sample] .= phase_current
    zero_d_q_voltage_v[:, sample] .= machine_phase_to_rotor(
        phase_terminal,
        electrical_angle,
    )
    zero_d_q_current_a[:, sample] .= machine_phase_to_rotor(
        phase_current,
        electrical_angle,
    )
    electromagnetic_torque_nm[sample] = accepted.electromagnetic_torque_nm
    mechanical_speed_rad_s[:, sample] .= runtime.shaft_state.speed_rad_s
    mechanical_angle_rad[:, sample] .= runtime.shaft_state.angle_rad
    isempty(runtime.shaft_state.coupling_torque_nm) ||
        (shaft_coupling_torque_nm[:, sample] .= runtime.shaft_state.coupling_torque_nm)
    field_voltage_v[sample] = accepted.field_voltage_v
    rotor_port_power_w[sample] = accepted.rotor_port_power_w
    terminal_power_w[sample] = accepted.terminal_power_w
    magnetic_coenergy_j[sample] = accepted.magnetic_coenergy_j
    kinetic_energy_j[sample] = runtime.shaft_state.kinetic_energy_j
    elastic_energy_j[sample] = runtime.shaft_state.elastic_energy_j
    dissipated_energy_j[sample] = accepted.dissipated_energy_j
    energy_residual_j[sample] = accepted.maximum_energy_residual_j
    energy_quadrature_defect_j[sample] =
        accepted.maximum_energy_quadrature_defect_j
    event_count[sample] = accepted.event_count
    control_sample_count[sample] = runtime.control_state.sample_count
    return nothing
end

function simulate_modern_machine(
    preparation::ModernMachinePreparation,
    terminal_voltage;
    duration_s::Real,
    events=ModernMachineEvent[],
    terminal_nodes=(1, 2, 3, 0),
)
    duration = Float64(duration_s)
    isfinite(duration) && duration >= 0.0 || throw(ArgumentError(
        "machine simulation duration must be finite and nonnegative",
    ))
    step = preparation.specification.settings.timestep_s
    step_count_float = duration / step
    step_count = round(Int, step_count_float)
    abs(step_count_float - step_count) <= 256.0 * eps(Float64) *
        max(step_count_float, 1.0) || throw(ArgumentError(
            "machine simulation duration must be an integer number of timesteps",
        ))
    runtime = modern_machine_runtime(preparation, terminal_nodes; events=events)
    sample_count = step_count + 1
    mass_count = length(preparation.specification.shaft_masses)
    coupling_count = length(preparation.specification.shaft_couplings)
    time_s = zeros(sample_count)
    phase_voltage_v = zeros(3, sample_count)
    phase_current_a = zeros(3, sample_count)
    zero_d_q_voltage_v = zeros(3, sample_count)
    zero_d_q_current_a = zeros(3, sample_count)
    electromagnetic_torque_nm = zeros(sample_count)
    mechanical_speed_rad_s = zeros(mass_count, sample_count)
    mechanical_angle_rad = zeros(mass_count, sample_count)
    shaft_coupling_torque_nm = zeros(coupling_count, sample_count)
    field_voltage_v = zeros(sample_count)
    rotor_port_power_w = zeros(sample_count)
    terminal_power_w = zeros(sample_count)
    magnetic_coenergy_j = zeros(sample_count)
    kinetic_energy_j = zeros(sample_count)
    elastic_energy_j = zeros(sample_count)
    dissipated_energy_j = zeros(sample_count)
    energy_residual_j = zeros(sample_count)
    energy_quadrature_defect_j = zeros(sample_count)
    event_count = zeros(Int, sample_count)
    control_sample_count = zeros(Int, sample_count)
    _machine_trace_sample!(
        runtime,
        1,
        time_s,
        phase_voltage_v,
        phase_current_a,
        zero_d_q_voltage_v,
        zero_d_q_current_a,
        electromagnetic_torque_nm,
        mechanical_speed_rad_s,
        mechanical_angle_rad,
        shaft_coupling_torque_nm,
        field_voltage_v,
        rotor_port_power_w,
        terminal_power_w,
        magnetic_coenergy_j,
        kinetic_energy_j,
        elastic_energy_j,
        dissipated_energy_j,
        energy_residual_j,
        energy_quadrature_defect_j,
        event_count,
        control_sample_count,
    )
    for step_index in 1:step_count
        time = step_index * step
        voltage = applicable(terminal_voltage, time, runtime) ?
            terminal_voltage(time, runtime) : terminal_voltage(time)
        advance_modern_machine!(runtime, voltage; time_s=time)
        _machine_trace_sample!(
            runtime,
            step_index + 1,
            time_s,
            phase_voltage_v,
            phase_current_a,
            zero_d_q_voltage_v,
            zero_d_q_current_a,
            electromagnetic_torque_nm,
            mechanical_speed_rad_s,
            mechanical_angle_rad,
            shaft_coupling_torque_nm,
            field_voltage_v,
            rotor_port_power_w,
            terminal_power_w,
            magnetic_coenergy_j,
            kinetic_energy_j,
            elastic_energy_j,
            dissipated_energy_j,
            energy_residual_j,
            energy_quadrature_defect_j,
            event_count,
            control_sample_count,
        )
    end
    return ModernMachineTrace(
        time_s,
        phase_voltage_v,
        phase_current_a,
        zero_d_q_voltage_v,
        zero_d_q_current_a,
        electromagnetic_torque_nm,
        mechanical_speed_rad_s,
        mechanical_angle_rad,
        shaft_coupling_torque_nm,
        field_voltage_v,
        rotor_port_power_w,
        terminal_power_w,
        magnetic_coenergy_j,
        kinetic_energy_j,
        elastic_energy_j,
        dissipated_energy_j,
        energy_residual_j,
        energy_quadrature_defect_j,
        event_count,
        control_sample_count,
        modern_machine_result(runtime),
    )
end

function write_modern_machine_trace_csv(path::AbstractString, trace::ModernMachineTrace)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        mass_count = size(trace.mechanical_speed_rad_s, 1)
        coupling_count = size(trace.shaft_coupling_torque_nm, 1)
        headers = String[
            "time_s",
            "phase_a_voltage_v",
            "phase_b_voltage_v",
            "phase_c_voltage_v",
            "phase_a_current_a",
            "phase_b_current_a",
            "phase_c_current_a",
            "zero_current_a",
            "d_axis_current_a",
            "q_axis_current_a",
            "electromagnetic_torque_nm",
        ]
        append!(headers, ["shaft_mass_$(index)_speed_rad_s" for index in 1:mass_count])
        append!(headers, ["shaft_coupling_$(index)_torque_nm" for index in 1:coupling_count])
        append!(headers, [
            "terminal_power_w",
            "rotor_port_power_w",
            "magnetic_coenergy_j",
            "kinetic_energy_j",
            "elastic_energy_j",
            "dissipated_energy_j",
            "maximum_energy_residual_j",
            "maximum_energy_quadrature_defect_j",
            "event_count",
            "control_sample_count",
        ])
        println(io, join(headers, ','))
        for sample in eachindex(trace.time_s)
            values = Any[
                trace.time_s[sample],
                trace.phase_voltage_v[:, sample]...,
                trace.phase_current_a[:, sample]...,
                trace.zero_d_q_current_a[:, sample]...,
                trace.electromagnetic_torque_nm[sample],
                trace.mechanical_speed_rad_s[:, sample]...,
                trace.shaft_coupling_torque_nm[:, sample]...,
                trace.terminal_power_w[sample],
                trace.rotor_port_power_w[sample],
                trace.magnetic_coenergy_j[sample],
                trace.kinetic_energy_j[sample],
                trace.elastic_energy_j[sample],
                trace.dissipated_energy_j[sample],
                trace.energy_residual_j[sample],
                trace.energy_quadrature_defect_j[sample],
                trace.event_count[sample],
                trace.control_sample_count[sample],
            ]
            println(io, join(values, ','))
        end
    end
    return String(path)
end

function modern_machine_report_text(trace::ModernMachineTrace)
    result = trace.result
    diagnostics = result.diagnostics
    return join((
        "AIMORA modern EMT machine report",
        "machine=$(result.machine)",
        "family=$(result.family)",
        "operating_mode=$(result.operating_mode)",
        "signature_sha256=$(result.deterministic_signature_sha256)",
        "accepted_steps=$(diagnostics.accepted_step_count)",
        "events=$(diagnostics.event_count)",
        "control_samples=$(diagnostics.control_sample_count)",
        "maximum_flux_residual_wb=$(diagnostics.maximum_flux_residual_wb)",
        "maximum_kcl_residual_a=$(diagnostics.maximum_kcl_residual_a)",
        "maximum_energy_residual_j=$(diagnostics.maximum_energy_residual_j)",
        "maximum_energy_quadrature_defect_j=$(diagnostics.maximum_energy_quadrature_defect_j)",
        "maximum_angular_momentum_residual_nms=$(diagnostics.maximum_angular_momentum_residual_nms)",
        "limitation=generic fixed-step phase-domain family; no vendor, standard, protection, thermal, field, HIL, or certification equivalence",
    ), '\n') * "\n"
end

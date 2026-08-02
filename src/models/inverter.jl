module Inverter

using Printf
using ..StudyCore: AverageValue,
                   ContractQuantity,
                   DynamicStateInventory,
                   ModelFidelity,
                   ModelValidityDomain,
                   NumericDomainBound,
                   ScientificModelContract,
                   assess_validity,
                   assert_validity

export InverterParams,
       InverterState,
       inverter_contract,
       initial_inverter_state,
       step_inverter,
       inverter_row,
       simulate_inverter,
       inverter_summary,
       write_inverter_csv,
       write_inverter_summary

Base.@kwdef struct InverterParams
    f_hz::Float64 = 60.0
    s_base_va::Float64 = 1.0e6
    v_ll_rms_v::Float64 = 4160.0
    r_filter_ohm::Float64 = 0.05
    l_filter_h::Float64 = 1.5e-3
    v_dc_v::Float64 = 7000.0
    kp_i::Float64 = 7.5
    ki_i::Float64 = 1800.0
    p_ref_0_pu::Float64 = 0.45
    p_ref_1_pu::Float64 = 0.85
    q_ref_pu::Float64 = 0.10
    step_time_s::Float64 = 0.050
end

const AVERAGE_VALUE_INVERTER_CONTRACT = ScientificModelContract(
    :average_value_grid_following_inverter,
    :dq_average_current_control;
    owner = "AIMORA.Inverter",
    maturity = :prototype,
    fidelity = AverageValue,
    validity_domain = ModelValidityDomain(
        :positive_finite_fixed_step;
        description = "Balanced positive-sequence dq average-value current controller with a stiff sinusoidal grid, constant DC voltage, lumped series R-L filter, fixed-step RK4 integration, and no switching or network-coupled event state.",
        bounds = (
            NumericDomainBound(:frequency_hz; unit = "Hz", lower = 0.0, lower_inclusive = false),
            NumericDomainBound(:base_power_va; unit = "VA", lower = 0.0, lower_inclusive = false),
            NumericDomainBound(:line_voltage_rms_v; unit = "V", lower = 0.0, lower_inclusive = false),
            NumericDomainBound(:filter_resistance_ohm; unit = "ohm", lower = 0.0),
            NumericDomainBound(:filter_inductance_h; unit = "H", lower = 0.0, lower_inclusive = false),
            NumericDomainBound(:dc_voltage_v; unit = "V", lower = 0.0, lower_inclusive = false),
            NumericDomainBound(:proportional_gain; unit = "V/A", lower = 0.0),
            NumericDomainBound(:integral_gain; unit = "V/(A*s)", lower = 0.0),
            NumericDomainBound(:initial_active_power_pu; unit = "pu"),
            NumericDomainBound(:final_active_power_pu; unit = "pu"),
            NumericDomainBound(:reactive_power_pu; unit = "pu"),
            NumericDomainBound(:reference_step_time_s; unit = "s"),
            NumericDomainBound(:timestep_s; unit = "s", lower = 0.0, lower_inclusive = false),
            NumericDomainBound(:end_time_s; unit = "s", lower = 0.0),
        ),
        unsupported_phenomena = (
            :individual_semiconductor_state,
            :switching_ripple,
            :dc_link_energy_dynamics,
            :sampled_control_delay,
            :network_coupled_events,
            :within_step_rollback,
            :unbalanced_sequence_dynamics,
            :fault_blocking_and_restart,
            :device_loss_or_thermal_state,
        ),
    ),
    state_inventory = DynamicStateInventory(
        differential = (:id_a, :iq_a, :xid, :xiq),
    ),
    inputs = (
        ContractQuantity(:time_s; unit = "s"),
        ContractQuantity(:active_power_reference_pu; unit = "pu", base = "s_base_va", orientation = "positive_grid_injection"),
        ContractQuantity(:reactive_power_reference_pu; unit = "pu", base = "s_base_va", orientation = "positive_grid_injection"),
        ContractQuantity(:grid_d_axis_voltage_v; unit = "V", orientation = "grid_to_converter_reference"),
    ),
    outputs = (
        ContractQuantity(:d_axis_current_pu; unit = "pu", base = "peak_current_base", orientation = "positive_grid_injection"),
        ContractQuantity(:q_axis_current_pu; unit = "pu", base = "peak_current_base", orientation = "positive_grid_injection"),
        ContractQuantity(:active_power_pu; unit = "pu", base = "s_base_va", orientation = "positive_grid_injection"),
        ContractQuantity(:reactive_power_pu; unit = "pu", base = "s_base_va", orientation = "positive_grid_injection"),
    ),
    assumptions = (
        "The grid voltage is a stiff balanced positive-sequence sinusoid aligned with the d axis.",
        "The DC voltage is constant and only limits the commanded dq voltage magnitude.",
        "The requested end time is represented by the nearest integer number of fixed timesteps.",
        "This prototype is not a switching-detailed converter and cannot satisfy a switching-detailed request.",
    ),
    mutation_order = (
        :sample_outputs,
        :read_power_references,
        :evaluate_rk4_stages,
        :apply_voltage_limit,
        :commit_state,
    ),
)

inverter_contract() = AVERAGE_VALUE_INVERTER_CONTRACT

struct InverterState
    id_a::Float64
    iq_a::Float64
    xid::Float64
    xiq::Float64
end

function base_values(p::InverterParams)
    omega = 2.0 * pi * p.f_hz
    v_phase_peak = sqrt(2.0) * p.v_ll_rms_v / sqrt(3.0)
    i_base_peak = sqrt(2.0) * p.s_base_va / (sqrt(3.0) * p.v_ll_rms_v)
    return omega, v_phase_peak, i_base_peak
end

function references(t::Float64, p::InverterParams)
    _, vd, _ = base_values(p)
    p_ref = (t < p.step_time_s ? p.p_ref_0_pu : p.p_ref_1_pu) * p.s_base_va
    q_ref = p.q_ref_pu * p.s_base_va
    id_ref = 2.0 * p_ref / (3.0 * vd)
    iq_ref = -2.0 * q_ref / (3.0 * vd)
    return p_ref, q_ref, id_ref, iq_ref
end

function limit_voltage(vd::Float64, vq::Float64, p::InverterParams)
    limit = p.v_dc_v / sqrt(3.0)
    mag = hypot(vd, vq)
    if mag <= limit || mag == 0.0
        return vd, vq
    end
    scale = limit / mag
    return vd * scale, vq * scale
end

function derivative(x::InverterState, t::Float64, p::InverterParams)
    omega, vd_grid, _ = base_values(p)
    _, _, id_ref, iq_ref = references(t, p)

    ed = id_ref - x.id_a
    eq = iq_ref - x.iq_a

    vd_cmd = vd_grid + p.r_filter_ohm * x.id_a - omega * p.l_filter_h * x.iq_a + p.kp_i * ed + p.ki_i * x.xid
    vq_cmd = p.r_filter_ohm * x.iq_a + omega * p.l_filter_h * x.id_a + p.kp_i * eq + p.ki_i * x.xiq
    vd_inv, vq_inv = limit_voltage(vd_cmd, vq_cmd, p)

    did = (vd_inv - vd_grid - p.r_filter_ohm * x.id_a + omega * p.l_filter_h * x.iq_a) / p.l_filter_h
    diq = (vq_inv - p.r_filter_ohm * x.iq_a - omega * p.l_filter_h * x.id_a) / p.l_filter_h
    return InverterState(did, diq, ed, eq)
end

function add_state(a::InverterState, b::InverterState, scale::Float64)
    return InverterState(
        a.id_a + scale * b.id_a,
        a.iq_a + scale * b.iq_a,
        a.xid + scale * b.xid,
        a.xiq + scale * b.xiq,
    )
end

function rk4_step(x::InverterState, t::Float64, dt::Float64, p::InverterParams)
    k1 = derivative(x, t, p)
    k2 = derivative(add_state(x, k1, 0.5 * dt), t + 0.5 * dt, p)
    k3 = derivative(add_state(x, k2, 0.5 * dt), t + 0.5 * dt, p)
    k4 = derivative(add_state(x, k3, dt), t + dt, p)
    return InverterState(
        x.id_a + dt * (k1.id_a + 2.0 * k2.id_a + 2.0 * k3.id_a + k4.id_a) / 6.0,
        x.iq_a + dt * (k1.iq_a + 2.0 * k2.iq_a + 2.0 * k3.iq_a + k4.iq_a) / 6.0,
        x.xid + dt * (k1.xid + 2.0 * k2.xid + 2.0 * k3.xid + k4.xid) / 6.0,
        x.xiq + dt * (k1.xiq + 2.0 * k2.xiq + 2.0 * k3.xiq + k4.xiq) / 6.0,
    )
end

function initial_inverter_state(p::InverterParams = InverterParams())
    _, _, id0, iq0 = references(0.0, p)
    return InverterState(id0, iq0, 0.0, 0.0)
end

step_inverter(x::InverterState, t::Float64, dt::Float64, p::InverterParams) = rk4_step(x, t, dt, p)

function inverter_row(x::InverterState, t::Float64, p::InverterParams)
    _, vd, i_base = base_values(p)
    p_ref, q_ref, id_ref, iq_ref = references(t, p)
    p_out = 1.5 * vd * x.id_a
    q_out = -1.5 * vd * x.iq_a
    return (
        t,
        x.id_a / i_base,
        x.iq_a / i_base,
        id_ref / i_base,
        iq_ref / i_base,
        p_out / p.s_base_va,
        q_out / p.s_base_va,
        p_ref / p.s_base_va,
        q_ref / p.s_base_va,
    )
end

function simulate_inverter(;
    t_end::Float64 = 0.150,
    dt::Float64 = 20e-6,
    p::InverterParams = InverterParams(),
    fidelity::ModelFidelity = AverageValue,
)
    contract = inverter_contract()
    assessment = assess_validity(
        contract,
        (
            frequency_hz = p.f_hz,
            base_power_va = p.s_base_va,
            line_voltage_rms_v = p.v_ll_rms_v,
            filter_resistance_ohm = p.r_filter_ohm,
            filter_inductance_h = p.l_filter_h,
            dc_voltage_v = p.v_dc_v,
            proportional_gain = p.kp_i,
            integral_gain = p.ki_i,
            initial_active_power_pu = p.p_ref_0_pu,
            final_active_power_pu = p.p_ref_1_pu,
            reactive_power_pu = p.q_ref_pu,
            reference_step_time_s = p.step_time_s,
            timestep_s = dt,
            end_time_s = t_end,
        );
        requested_fidelity = fidelity,
    )
    assert_validity(assessment)
    steps = Int(round(t_end / dt))
    x = initial_inverter_state(p)

    rows = Vector{NTuple{9, Float64}}(undef, steps + 1)
    for n in 0:steps
        t = n * dt
        rows[n + 1] = inverter_row(x, t, p)
        n < steps && (x = step_inverter(x, t, dt, p))
    end
    return rows
end

function inverter_summary(rows; elapsed_s::Float64 = 0.0)
    last = rows[end]
    return (
        samples = length(rows),
        t_end_s = last[1],
        elapsed_s = elapsed_s,
        final_p_pu = last[6],
        final_q_pu = last[7],
        final_id_pu = last[2],
        final_iq_pu = last[3],
        max_p_error_pu = maximum(abs(r[6] - r[8]) for r in rows),
        max_q_error_pu = maximum(abs(r[7] - r[9]) for r in rows),
    )
end

function ensure_dir(path::AbstractString)
    isdir(path) || mkpath(path)
end

function write_inverter_csv(path::AbstractString, rows)
    ensure_dir(dirname(path))
    open(path, "w") do io
        println(io, "time_s,id_pu,iq_pu,id_ref_pu,iq_ref_pu,p_pu,q_pu,p_ref_pu,q_ref_pu")
        for r in rows
            @printf(io, "%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n", r...)
        end
    end
end

function write_inverter_summary(path::AbstractString, summary)
    ensure_dir(dirname(path))
    open(path, "w") do io
        println(io, "{")
        @printf(io, "  \"samples\": %d,\n", summary.samples)
        @printf(io, "  \"t_end_s\": %.9f,\n", summary.t_end_s)
        @printf(io, "  \"elapsed_s\": %.9f,\n", summary.elapsed_s)
        @printf(io, "  \"final_p_pu\": %.9f,\n", summary.final_p_pu)
        @printf(io, "  \"final_q_pu\": %.9f,\n", summary.final_q_pu)
        @printf(io, "  \"final_id_pu\": %.9f,\n", summary.final_id_pu)
        @printf(io, "  \"final_iq_pu\": %.9f,\n", summary.final_iq_pu)
        @printf(io, "  \"max_p_error_pu\": %.9f,\n", summary.max_p_error_pu)
        @printf(io, "  \"max_q_error_pu\": %.9f\n", summary.max_q_error_pu)
        println(io, "}")
    end
end

end

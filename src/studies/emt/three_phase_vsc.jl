import ..Branches
import ..BridgeTopologies
import ..Nodal
import ..Nonlinear
import ..OVER16TimestepIntegration
import ..SwitchDetailedVSC

export ThreePhaseVSCBoundaryOccurrence,
       ThreePhaseVSCCommutationOccurrence,
       ThreePhaseVSCHarmonicMetadata,
       ThreePhaseVSCMetrics,
       ThreePhaseVSCTrace,
       ThreePhaseVSCIntegrator,
       configure_three_phase_vsc,
       advance_three_phase_vsc!,
       run_three_phase_vsc!,
       three_phase_vsc_trace,
       simulate_three_phase_vsc,
       write_three_phase_vsc_checkpoint,
       read_three_phase_vsc_checkpoint

const _THREE_PHASE_VSC_CHECKPOINT_MAGIC = collect(codeunits("AIMORA-THREE-PHASE-VSC"))
const _THREE_PHASE_VSC_CHECKPOINT_SCHEMA = UInt32(1)
const _VSC_DC_POSITIVE_NODE = 1
const _VSC_DC_NEGATIVE_NODE = 2
const ThreePhaseVSCBoundaryPlan = NTuple{6,Tuple{Symbol,Float64,Bool}}

_three_phase_vsc_pole_node(phase::Int) = phase + 2
_three_phase_vsc_grid_terminal_node(phase::Int) = phase + 5

"""Exact converter disturbance or protection boundary recorded on the scheduler tick calendar."""
struct ThreePhaseVSCBoundaryOccurrence
    name::Symbol
    time_s::Float64
    tick::Int
    topology_invalidating::Bool
end

"""One conduction-state transition resolved against the coupled converter network at an exact accepted scheduler tick."""
struct ThreePhaseVSCCommutationOccurrence
    name::Symbol
    time_s::Float64
    tick::Int
    phase::Int
    position::Symbol
    transition::Symbol
    residual::Float64
    topology_iteration::Int
end

"""Synchronous-window Fourier metadata for the three phase currents and converter line voltages."""
struct ThreePhaseVSCHarmonicMetadata
    window_start_s::Float64
    window_end_s::Float64
    sample_count::Int
    fundamental_frequency_hz::Float64
    maximum_harmonic_order::Int
    phase_current_fundamental_rms_a::NTuple{3,Float64}
    phase_current_thd::NTuple{3,Float64}
    line_voltage_fundamental_rms_v::NTuple{3,Float64}
    line_voltage_thd::NTuple{3,Float64}
end

"""Accepted physical, numerical, protection, timing, distortion, and energy summary of one converter horizon."""
struct ThreePhaseVSCMetrics
    mean_active_power_w::Float64
    mean_reactive_power_var::Float64
    dc_link_mean_voltage_v::Float64
    dc_link_peak_to_peak_ripple_v::Float64
    maximum_nodal_kcl_residual_a::Float64
    dc_source_energy_j::Float64
    ac_terminal_energy_j::Float64
    converter_dissipated_energy_j::Float64
    dc_ac_energy_residual_j::Float64
    relative_dc_ac_energy_residual::Float64
    external_energy_j::Float64
    dissipated_energy_j::Float64
    stored_energy_change_j::Float64
    energy_residual_j::Float64
    relative_energy_residual::Float64
    control_sample_count::Int
    control_write_count::Int
    pwm_cycle_counts::NTuple{3,Int}
    pwm_edge_counts::NTuple{3,Int}
    bridge_topology_transition_counts::NTuple{3,Int}
    shoot_through_rejection_count::Int
    all_off_dead_time_sample_count::Int
    antiparallel_diode_sample_count::Int
    block_count::Int
    restart_count::Int
    boundary_count::Int
    commutation_count::Int
    maximum_topology_iterations::Int
    exact_boundary_alignment::Bool
    finite_output::Bool
    harmonic::ThreePhaseVSCHarmonicMetadata
end

"""Typed device, terminal, control, power, energy, event, and harmonic result for a switch-detailed three-phase VSC horizon."""
struct ThreePhaseVSCTrace
    parameters::SwitchDetailedVSC.ThreePhaseTwoLevelVSCParameters
    time_s::Vector{Float64}
    pole_voltage_v::Matrix{Float64}
    line_voltage_v::Matrix{Float64}
    grid_terminal_voltage_v::Matrix{Float64}
    filter_current_a::Matrix{Float64}
    dc_link_voltage_v::Vector{Float64}
    dc_source_current_a::Vector{Float64}
    active_power_w::Vector{Float64}
    reactive_power_var::Vector{Float64}
    duty::Matrix{Float64}
    upper_gate_applied::Matrix{Float64}
    lower_gate_applied::Matrix{Float64}
    upper_forward_conducting::Matrix{Float64}
    lower_forward_conducting::Matrix{Float64}
    upper_antiparallel_diode_conducting::Matrix{Float64}
    lower_antiparallel_diode_conducting::Matrix{Float64}
    semiconductor_loss_w::Vector{Float64}
    resistive_loss_w::Vector{Float64}
    stored_energy_j::Vector{Float64}
    external_source_power_w::Vector{Float64}
    nodal_kcl_residual_a::Vector{Float64}
    blocked::Vector{Float64}
    boundary_occurrences::Vector{ThreePhaseVSCBoundaryOccurrence}
    commutation_occurrences::Vector{ThreePhaseVSCCommutationOccurrence}
    sampled_task_occurrences::Vector{OVER16TimestepIntegration.SampledTaskOccurrence}
    pwm_edge_ticks::NTuple{3,Vector{Int}}
    metrics::ThreePhaseVSCMetrics
end

struct ThreePhaseGridVoltageSignal
    parameters::SwitchDetailedVSC.ThreePhaseTwoLevelVSCParameters
    phase::Int
end

function (signal::ThreePhaseGridVoltageSignal)(time_s::Real)
    p = signal.parameters
    time = Float64(time_s)
    angle = 2.0 * pi * p.frequency_hz * time -
        (signal.phase - 1) * 2.0 * pi / 3.0
    phase_crest_v = sqrt(2.0 / 3.0) * p.grid_line_line_rms_v *
        p.transformer_grid_to_converter_ratio
    factor = p.sag_start_s <= time < p.sag_end_s ? p.sag_voltage_factor : 1.0
    if signal.phase == p.faulted_phase && p.fault_start_s <= time < p.fault_end_s
        factor *= p.fault_voltage_factor
    end
    return factor * phase_crest_v * cos(angle)
end

struct ConstantDCVoltageSignal
    voltage_v::Float64
end

(signal::ConstantDCVoltageSignal)(_time_s::Real) = signal.voltage_v

mutable struct ThreePhaseVSCRuntime{S,B,C,F}
    parameters::SwitchDetailedVSC.ThreePhaseTwoLevelVSCParameters
    controller::SwitchDetailedVSC.ThreePhaseGridFollowingControllerState
    system::S
    bridges::B
    dc_link_capacitor::C
    phase_filters::F
    boundary_plan::ThreePhaseVSCBoundaryPlan
    boundary_occurrences::Vector{ThreePhaseVSCBoundaryOccurrence}
    commutation_occurrences::Vector{ThreePhaseVSCCommutationOccurrence}
    next_boundary_index::Int
    block_count::Int
    restart_count::Int
    protection_trip_count::Int
    maximum_topology_iterations::Int
end

struct ThreePhaseVSCControlReader end
struct ThreePhaseVSCControlComputer end
struct ThreePhaseVSCControlWriter end

struct ThreePhaseVSCDutyReader
    phase::Int
end

struct ThreePhaseVSCBridgeGateWriter
    phase::Int
end

function _three_phase_vsc_measurement(runtime::ThreePhaseVSCRuntime, time_s::Float64)
    system = runtime.system
    phase_voltage = ntuple(
        phase -> system.v[_three_phase_vsc_grid_terminal_node(phase)],
        3,
    )
    phase_current = ntuple(phase -> runtime.phase_filters[phase].i_last, 3)
    return SwitchDetailedVSC.ThreePhaseVSCMeasurement(
        phase_voltage,
        phase_current,
        system.v[_VSC_DC_POSITIVE_NODE] - system.v[_VSC_DC_NEGATIVE_NODE],
        2.0 * pi * runtime.parameters.frequency_hz * time_s,
    )
end

function (::ThreePhaseVSCControlReader)(
    runtime::ThreePhaseVSCRuntime,
    time_s::Real,
    _sample_index::Int,
)
    return _three_phase_vsc_measurement(runtime, Float64(time_s))
end

function (::ThreePhaseVSCControlComputer)(
    runtime::ThreePhaseVSCRuntime,
    measurement::SwitchDetailedVSC.ThreePhaseVSCMeasurement,
    _time_s::Real,
    _sample_index::Int,
)
    maximum(abs, measurement.phase_current_a) > runtime.parameters.current_limit_a &&
        (runtime.controller.overcurrent_count += 1)
    return SwitchDetailedVSC.compute_grid_following_current_control!(
        runtime.controller,
        measurement,
        runtime.parameters,
    )
end

function (::ThreePhaseVSCControlWriter)(
    runtime::ThreePhaseVSCRuntime,
    command::SwitchDetailedVSC.ThreePhaseVSCModulationCommand,
    _time_s::Real,
    _sample_index::Int,
)
    runtime.controller.held_duties = command.duties
    runtime.controller.held_phase_voltage_reference_v = command.phase_voltage_reference_v
    runtime.controller.write_count += 1
    return nothing
end

function (reader::ThreePhaseVSCDutyReader)(
    runtime::ThreePhaseVSCRuntime,
    _time_s::Real,
    _cycle_index::Int,
)
    return runtime.controller.held_duties[reader.phase]
end

function (writer::ThreePhaseVSCBridgeGateWriter)(
    runtime::ThreePhaseVSCRuntime,
    high::Bool,
    time_s::Real,
    _edge_index::Int,
)
    Nonlinear.request_power_semiconductor_bridge_pole!(
        runtime.bridges[writer.phase],
        high,
        time_s,
    )
    return nothing
end

mutable struct ThreePhaseVSCRecorder
    time_s::Vector{Float64}
    pole_voltage_v::Matrix{Float64}
    line_voltage_v::Matrix{Float64}
    grid_terminal_voltage_v::Matrix{Float64}
    filter_current_a::Matrix{Float64}
    dc_link_voltage_v::Vector{Float64}
    dc_source_current_a::Vector{Float64}
    active_power_w::Vector{Float64}
    reactive_power_var::Vector{Float64}
    duty::Matrix{Float64}
    upper_gate_applied::Matrix{Float64}
    lower_gate_applied::Matrix{Float64}
    upper_forward_conducting::Matrix{Float64}
    lower_forward_conducting::Matrix{Float64}
    upper_antiparallel_diode_conducting::Matrix{Float64}
    lower_antiparallel_diode_conducting::Matrix{Float64}
    semiconductor_loss_w::Vector{Float64}
    resistive_loss_w::Vector{Float64}
    stored_energy_j::Vector{Float64}
    external_source_power_w::Vector{Float64}
    nodal_kcl_residual_a::Vector{Float64}
    blocked::Vector{Float64}
    write_index::Int
end

function ThreePhaseVSCRecorder(sample_count::Int)
    sample_count > 1 || throw(ArgumentError("VSC trace requires at least two samples"))
    vector() = Vector{Float64}(undef, sample_count)
    matrix() = Matrix{Float64}(undef, 3, sample_count)
    return ThreePhaseVSCRecorder(
        vector(), matrix(), matrix(), matrix(), matrix(), vector(), vector(),
        vector(), vector(), matrix(), matrix(), matrix(), matrix(), matrix(),
        matrix(), matrix(), vector(), vector(), vector(), vector(), vector(),
        vector(), 0,
    )
end

mutable struct ThreePhaseVSCIntegrator{R,S,T,C,P}
    runtime::R
    scheduler::S
    transaction::T
    controller_task::C
    pwm_tasks::P
    recorder::ThreePhaseVSCRecorder
    accepted_step_index::Int
    completed::Bool
    failed::Bool
    last_failure::Union{Nothing,String}
end

# Keep the historical integrator and its serialized layout intact. The
# explicit causal route retains its selection across checkpoints.
struct CausalThreePhaseVSCIntegrator{I<:ThreePhaseVSCIntegrator}
    state::I
end

_three_phase_vsc_state(integrator::ThreePhaseVSCIntegrator) = integrator
_three_phase_vsc_state(integrator::CausalThreePhaseVSCIntegrator) = integrator.state

function _three_phase_vsc_bridge(
    parameters::SwitchDetailedVSC.ThreePhaseTwoLevelVSCParameters,
    pole_node::Int,
)
    upper_driver = Nonlinear.PowerSemiconductorGateDriver(
        turn_on_delay_s = parameters.gate_turn_on_delay_s,
        turn_off_delay_s = parameters.gate_turn_off_delay_s,
        minimum_pulse_width_s = parameters.minimum_gate_pulse_width_s,
    )
    lower_driver = Nonlinear.PowerSemiconductorGateDriver(
        turn_on_delay_s = parameters.gate_turn_on_delay_s,
        turn_off_delay_s = parameters.gate_turn_off_delay_s,
        minimum_pulse_width_s = parameters.minimum_gate_pulse_width_s,
    )
    diode = Nonlinear.AntiparallelDiodeParameters(
        forward_voltage_v = parameters.diode_forward_voltage_v,
        on_conductance_s = parameters.diode_on_conductance_s,
    )
    upper = Nonlinear.IGBTSwitch(
        _VSC_DC_POSITIVE_NODE,
        pole_node;
        gate_driver = upper_driver,
        forward_voltage_drop_v = parameters.semiconductor_forward_voltage_v,
        on_conductance = parameters.semiconductor_on_conductance_s,
        off_conductance = parameters.semiconductor_off_conductance_s,
        antiparallel_diode = diode,
    )
    lower = Nonlinear.IGBTSwitch(
        pole_node,
        _VSC_DC_NEGATIVE_NODE;
        gate_driver = lower_driver,
        forward_voltage_drop_v = parameters.semiconductor_forward_voltage_v,
        on_conductance = parameters.semiconductor_on_conductance_s,
        off_conductance = parameters.semiconductor_off_conductance_s,
        antiparallel_diode = diode,
    )
    return Nonlinear.PowerSemiconductorBridgeLeg(
        upper,
        lower;
        commutation_dead_time_s = parameters.commutation_dead_time_s,
    )
end

function _three_phase_vsc_runtime(
    parameters::SwitchDetailedVSC.ThreePhaseTwoLevelVSCParameters,
)
    p = SwitchDetailedVSC.validate_three_phase_vsc_parameters(parameters)
    bridges = ntuple(
        phase -> _three_phase_vsc_bridge(p, _three_phase_vsc_pole_node(phase)),
        3,
    )
    foreach(Nonlinear.power_semiconductor_event_localization!, bridges)
    dc_link_capacitor = Branches.CapacitorBranch(
        _VSC_DC_POSITIVE_NODE,
        _VSC_DC_NEGATIVE_NODE,
        p.dc_link_capacitance_f,
        0.0,
        p.dc_source_voltage_v,
        0.0,
    )
    series_resistance = p.filter_resistance_ohm +
        p.transformer_leakage_resistance_ohm
    series_inductance = p.filter_inductance_h +
        p.transformer_leakage_inductance_h
    phase_filters = ntuple(
        phase -> Branches.SeriesRLBranch(
            _three_phase_vsc_pole_node(phase),
            _three_phase_vsc_grid_terminal_node(phase),
            series_resistance,
            series_inductance,
        ),
        3,
    )
    dc_source = Branches.TwoTerminalTheveninSource(
        _VSC_DC_POSITIVE_NODE,
        _VSC_DC_NEGATIVE_NODE,
        inv(p.dc_source_resistance_ohm),
        ConstantDCVoltageSignal(p.dc_source_voltage_v),
    )
    grid_sources = ntuple(
        phase -> Branches.TheveninSource(
            _three_phase_vsc_grid_terminal_node(phase),
            inv(p.grid_source_resistance_ohm),
            ThreePhaseGridVoltageSignal(p, phase),
        ),
        3,
    )
    elements = Any[
        bridges...,
        dc_link_capacitor,
        phase_filters...,
        dc_source,
        grid_sources...,
    ]
    system = Nodal.NodalSystem(8, elements)
    system.v[_VSC_DC_POSITIVE_NODE] = 0.5 * p.dc_source_voltage_v
    system.v[_VSC_DC_NEGATIVE_NODE] = -0.5 * p.dc_source_voltage_v
    for phase in 1:3
        grid_voltage = grid_sources[phase].value(0.0)
        system.v[_three_phase_vsc_pole_node(phase)] = grid_voltage
        system.v[_three_phase_vsc_grid_terminal_node(phase)] = grid_voltage
    end
    initial_phase_reference = ntuple(phase -> grid_sources[phase].value(0.0), 3)
    initial_duties = SwitchDetailedVSC.modulation_duties(
        initial_phase_reference,
        p.dc_source_voltage_v,
        p.modulation;
        minimum_duty = p.minimum_duty,
        maximum_duty = p.maximum_duty,
    )
    controller = SwitchDetailedVSC.ThreePhaseGridFollowingControllerState()
    controller.held_duties = initial_duties
    controller.held_phase_voltage_reference_v = initial_phase_reference
    return ThreePhaseVSCRuntime(
        p,
        controller,
        system,
        bridges,
        dc_link_capacitor,
        phase_filters,
        _three_phase_vsc_boundary_plan(p),
        ThreePhaseVSCBoundaryOccurrence[],
        ThreePhaseVSCCommutationOccurrence[],
        1,
        0,
        0,
        0,
        0,
    )
end

function _three_phase_vsc_boundary_plan(
    parameters::SwitchDetailedVSC.ThreePhaseTwoLevelVSCParameters,
)
    return Tuple(sort!(
        [
            (:grid_sag_begins, parameters.sag_start_s, false),
            (:grid_sag_ends, parameters.sag_end_s, false),
            (:phase_to_ground_fault_begins, parameters.fault_start_s, false),
            (:converter_protection_blocks, parameters.block_time_s, true),
            (:phase_to_ground_fault_clears, parameters.fault_end_s, false),
            (:converter_protection_restarts, parameters.restart_time_s, true),
        ];
        by = row -> (row[2], String(row[1])),
    ))
end

function _three_phase_vsc_scheduler(runtime::ThreePhaseVSCRuntime)
    p = runtime.parameters
    controller = ExactSampledControlTask(
        :grid_following_current_controller,
        p.control_period_s,
        ThreePhaseVSCControlReader(),
        ThreePhaseVSCControlComputer(),
        ThreePhaseVSCControlWriter();
        tick_s = p.scheduler_tick_s,
        computational_delay_s = p.control_delay_s,
        initial_output = SwitchDetailedVSC.ThreePhaseVSCModulationCommand(
            runtime.controller.held_duties,
            runtime.controller.held_phase_voltage_reference_v,
            0.0,
            0.0,
            0.0,
            0.0,
            false,
        ),
        priority = -10,
        power_history_invalidating = true,
    )
    carrier_period_s = inv(p.carrier_frequency_hz)
    pwm_tasks = ntuple(3) do phase
        ExactPWMTask(
            Symbol(:phase_, Char(Int('a') + phase - 1), :_carrier),
            carrier_period_s,
            ThreePhaseVSCDutyReader(phase),
            ThreePhaseVSCBridgeGateWriter(phase);
            tick_s = p.scheduler_tick_s,
            priority = 0,
            power_history_invalidating = true,
        )
    end
    scheduler = ExactSampledTaskScheduler(
        p.scheduler_tick_s;
        tasks = [controller, pwm_tasks...],
    )
    return scheduler, controller, pwm_tasks
end

function _three_phase_vsc_task_states(integrator::ThreePhaseVSCIntegrator)
    return [sampled_task_checkpoint(task) for task in integrator.scheduler.tasks]
end

function _restore_three_phase_vsc_task_states!(integrator, states, occurrence_count)
    for (task, state) in zip(integrator.scheduler.tasks, states)
        restore_sampled_task_checkpoint!(task, state)
    end
    resize!(integrator.scheduler.occurrences, occurrence_count)
    return integrator
end

function _record_three_phase_vsc_boundary!(
    runtime::ThreePhaseVSCRuntime,
    name::Symbol,
    time_s::Float64,
    topology_invalidating::Bool,
)
    tick = round(Int, time_s / runtime.parameters.scheduler_tick_s)
    push!(
        runtime.boundary_occurrences,
        ThreePhaseVSCBoundaryOccurrence(name, time_s, tick, topology_invalidating),
    )
    return runtime
end

function _resynchronize_three_phase_vsc_gates!(integrator, time_s::Float64)
    for phase in 1:3
        request = integrator.pwm_tasks[phase].gate_high
        Nonlinear.request_power_semiconductor_bridge_pole!(
            integrator.runtime.bridges[phase],
            request,
            time_s,
        )
    end
    return integrator
end

function _apply_three_phase_vsc_boundaries!(
    integrator::ThreePhaseVSCIntegrator,
    time_s::Float64,
)
    runtime = integrator.runtime
    plan = runtime.boundary_plan
    tolerance = 64.0 * eps(Float64) * max(1.0, abs(time_s))
    while runtime.next_boundary_index <= length(plan)
        name, boundary_time_s, topology_invalidating = plan[runtime.next_boundary_index]
        boundary_time_s > time_s + tolerance && break
        boundary_time_s < time_s - tolerance && throw(ArgumentError(
            "converter boundary $name was missed before $time_s s",
        ))
        if name === :converter_protection_blocks
            changed = false
            for bridge in runtime.bridges
                changed |= Nonlinear.block_power_semiconductor_bridge!(bridge, time_s)
            end
            changed && (runtime.block_count += 1)
        elseif name === :converter_protection_restarts
            changed = false
            for bridge in runtime.bridges
                changed |= Nonlinear.restart_power_semiconductor_bridge!(bridge, time_s)
            end
            if changed
                runtime.restart_count += 1
                _resynchronize_three_phase_vsc_gates!(integrator, time_s)
            end
        end
        _record_three_phase_vsc_boundary!(runtime, name, time_s, topology_invalidating)
        runtime.next_boundary_index += 1
    end
    return runtime
end

function _next_three_phase_vsc_boundary_time(
    runtime::ThreePhaseVSCRuntime,
    current_time_s::Float64,
    endpoint_time_s::Float64,
)
    plan = runtime.boundary_plan
    runtime.next_boundary_index <= length(plan) || return Inf
    boundary_time_s = plan[runtime.next_boundary_index][2]
    tolerance = 64.0 * eps(Float64) * max(
        1.0,
        abs(current_time_s),
        abs(endpoint_time_s),
    )
    boundary_time_s < current_time_s - tolerance && throw(ArgumentError(
        "converter boundary was missed before $current_time_s s",
    ))
    return boundary_time_s <= endpoint_time_s + tolerance ? boundary_time_s : Inf
end

function _next_three_phase_vsc_task_time(
    integrator::ThreePhaseVSCIntegrator,
    current_time_s::Float64,
    endpoint_time_s::Float64,
    tolerance_s::Float64,
)
    scheduler = integrator.scheduler
    next_tick = sampled_task_next_tick(integrator.controller_task)
    for task in integrator.pwm_tasks
        next_tick = min(next_tick, sampled_task_next_tick(task))
    end
    minimum_tick = ceil(
        Int,
        (current_time_s - scheduler.origin_s - tolerance_s) / scheduler.tick_s,
    )
    next_tick < minimum_tick && throw(ArgumentError(
        "converter sampled task missed tick $next_tick before requested time $current_time_s s",
    ))
    next_tick == typemax(Int) && return Inf
    task_time_s = scheduler.origin_s + next_tick * scheduler.tick_s
    return task_time_s <= endpoint_time_s + tolerance_s ? task_time_s : Inf
end

function _three_phase_vsc_applied_gates(runtime::ThreePhaseVSCRuntime)
    return map(runtime.bridges) do bridge
        (bridge.upper_switch.gate_driver.applied_on, bridge.lower_switch.gate_driver.applied_on)
    end
end

function _three_phase_vsc_source_edge(runtime::ThreePhaseVSCRuntime, time_s::Float64)
    p = runtime.parameters
    tolerance = 64 * eps(max(1.0, abs(time_s)))
    for edge_s in (p.sag_start_s, p.sag_end_s, p.fault_start_s, p.fault_end_s)
        abs(time_s - edge_s) <= tolerance && return edge_s
    end
    return nothing
end

struct ThreePhaseVSCIntervalVoltage
    start_time_s::Float64
    end_time_s::Float64
    start_filter_voltage_v::NTuple{3,Float64}
    end_filter_voltage_v::NTuple{3,Float64}
end

function _advance_three_phase_vsc_interval!(
    integrator::ThreePhaseVSCIntegrator,
    start_time_s::Float64,
    endpoint_time_s::Float64,
    intervals::Union{Nothing,Vector{ThreePhaseVSCIntervalVoltage}},
)
    runtime = integrator.runtime
    dt_s = endpoint_time_s - start_time_s
    source_edge = _three_phase_vsc_source_edge(runtime, endpoint_time_s)
    left_evaluation_s = source_edge === nothing ? endpoint_time_s : prevfloat(source_edge)

    _advance_three_phase_vsc_causal_storage!(
        runtime, start_time_s, endpoint_time_s, left_evaluation_s, intervals,
    )

    gates_before = _three_phase_vsc_applied_gates(runtime)
    _apply_three_phase_vsc_boundaries!(integrator, endpoint_time_s)
    if source_edge !== nothing
        # A measurement coincident with a source edge sees the new grid
        # voltage and the same accepted inductor current.
        p = runtime.parameters
        for (phase, filter) in enumerate(runtime.phase_filters)
            source_v = ThreePhaseGridVoltageSignal(p, phase)(source_edge)
            runtime.system.v[filter.b] = source_v + p.grid_source_resistance_ohm * filter.i_prev
        end
    end
    run_due_sampled_tasks!(integrator.scheduler, runtime, endpoint_time_s)
    for bridge in runtime.bridges
        Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(bridge, endpoint_time_s)
    end
    if _three_phase_vsc_applied_gates(runtime) != gates_before || source_edge !== nothing
        right_evaluation_s = source_edge === nothing ? endpoint_time_s : source_edge
        _resolve_three_phase_vsc_boundary!(
            runtime, endpoint_time_s; source_time_s = right_evaluation_s,
        )
    end
    return runtime
end

function _advance_historical_three_phase_vsc_window!(
    integrator::ThreePhaseVSCIntegrator,
    start_time_s::Float64,
    endpoint_time_s::Float64,
)
    runtime = integrator.runtime
    scheduler = integrator.scheduler
    current_time_s = start_time_s
    tolerance = runtime.parameters.scheduler_tick_s * 1.0e-9
    if endpoint_time_s - start_time_s <= scheduler.tick_s + tolerance
        _apply_three_phase_vsc_boundaries!(integrator, endpoint_time_s)
        run_due_sampled_tasks!(scheduler, runtime, endpoint_time_s)
        _stabilize_historical_three_phase_vsc_topology!(
            runtime,
            endpoint_time_s,
            endpoint_time_s - start_time_s,
        )
        Nodal.accept_algebraic_state!(
            runtime.system,
            endpoint_time_s - start_time_s,
        )
        return runtime
    end
    while current_time_s < endpoint_time_s - tolerance
        task_time_s = _next_three_phase_vsc_task_time(
            integrator,
            current_time_s,
            endpoint_time_s,
            tolerance,
        )
        boundary_time_s = _next_three_phase_vsc_boundary_time(
            runtime,
            current_time_s,
            endpoint_time_s,
        )
        segment_endpoint_s = min(endpoint_time_s, task_time_s, boundary_time_s)
        segment_timestep_s = segment_endpoint_s - current_time_s
        segment_timestep_s > tolerance || throw(ArgumentError(
            "converter scheduler produced a nonadvancing network interval at $current_time_s s",
        ))
        _apply_three_phase_vsc_boundaries!(integrator, segment_endpoint_s)
        run_due_sampled_tasks!(scheduler, runtime, segment_endpoint_s)
        _stabilize_historical_three_phase_vsc_topology!(
            runtime,
            segment_endpoint_s,
            segment_timestep_s,
        )
        Nodal.accept_algebraic_state!(runtime.system, segment_timestep_s)
        current_time_s = segment_endpoint_s
    end
    abs(current_time_s - endpoint_time_s) <= tolerance || throw(ArgumentError(
        "converter network window did not reach its requested endpoint",
    ))
    return runtime
end

function _advance_three_phase_vsc_network_window!(
    integrator::ThreePhaseVSCIntegrator,
    start_time_s::Float64,
    endpoint_time_s::Float64,
    intervals::Union{Nothing,Vector{ThreePhaseVSCIntervalVoltage}},
)
    runtime = integrator.runtime
    scheduler = integrator.scheduler
    tolerance = runtime.parameters.scheduler_tick_s * 1.0e-9
    if endpoint_time_s - start_time_s <= scheduler.tick_s + tolerance
        return _advance_three_phase_vsc_interval!(integrator, start_time_s, endpoint_time_s, intervals)
    end
    current_time_s = start_time_s
    while current_time_s < endpoint_time_s - tolerance
        task_time_s = _next_three_phase_vsc_task_time(
            integrator, current_time_s, endpoint_time_s, tolerance,
        )
        boundary_time_s = _next_three_phase_vsc_boundary_time(
            runtime, current_time_s, endpoint_time_s,
        )
        gate_time_s = minimum(runtime.bridges) do bridge
            something(Nonlinear.power_semiconductor_bridge_gate_transition_time(bridge), Inf)
        end
        segment_endpoint_s = min(endpoint_time_s, task_time_s, boundary_time_s, gate_time_s)
        segment_endpoint_s - current_time_s > tolerance || throw(ArgumentError(
            "converter scheduler produced a nonadvancing network interval at $current_time_s s",
        ))
        _advance_three_phase_vsc_interval!(integrator, current_time_s, segment_endpoint_s, intervals)
        current_time_s = segment_endpoint_s
    end
    abs(current_time_s - endpoint_time_s) <= tolerance || throw(ArgumentError(
        "converter network window did not reach its requested endpoint",
    ))
    return runtime
end

function _three_phase_vsc_device_current(device, voltages)
    terminal_voltage = Branches.branch_voltage(voltages, device.a, device.b)
    terminal_current = device.last_conductance * terminal_voltage +
        device.last_history_current_a
    return terminal_voltage, terminal_current
end

function _three_phase_vsc_topology_action(device, voltages)
    terminal_voltage, terminal_current = _three_phase_vsc_device_current(
        device,
        voltages,
    )
    current_tolerance_a = 1.0e-9
    voltage_tolerance_v = 1.0e-9
    if device.closed
        residual = terminal_current - device.holding_current
        residual < -current_tolerance_a && return (
            transition = :forward_extinction,
            residual,
        )
        return nothing
    end
    diode = device.antiparallel_diode
    if device.reverse_diode_conducting
        residual = -terminal_current - diode.holding_current_a
        residual < -current_tolerance_a && return (
            transition = :reverse_diode_extinction,
            residual,
        )
        return nothing
    end
    driver = device.gate_driver
    if driver === nothing
        forward_threshold = max(
            device.threshold_v,
            device.forward_voltage_drop_v + device.holding_current / device.on_conductance,
        )
        residual = terminal_voltage - forward_threshold
        residual > voltage_tolerance_v && return (
            transition = :forward_turn_on,
            residual,
        )
        return nothing
    end
    if driver !== nothing && driver.applied_on
        residual = terminal_voltage -
            Nonlinear.power_semiconductor_forward_turn_on_voltage(device)
        residual > voltage_tolerance_v && return (
            transition = :forward_turn_on,
            residual,
        )
    end
    if diode !== nothing
        residual = -terminal_voltage -
            Nonlinear.power_semiconductor_reverse_turn_on_voltage(device)
        residual > voltage_tolerance_v && return (
            transition = :reverse_diode_turn_on,
            residual,
        )
    end
    return nothing
end

function _apply_three_phase_vsc_topology_action!(
    device,
    transition::Symbol,
    time_s::Float64,
)
    transition === :forward_extinction &&
        return Nonlinear.apply_power_semiconductor_forward_extinction!(device, time_s)
    transition === :reverse_diode_extinction &&
        return Nonlinear.apply_power_semiconductor_reverse_extinction!(device, time_s)
    transition === :forward_turn_on &&
        return Nonlinear.apply_power_semiconductor_forward_turn_on!(device, time_s)
    transition === :reverse_diode_turn_on &&
        return Nonlinear.apply_power_semiconductor_reverse_turn_on!(device, time_s)
    throw(ArgumentError("unsupported converter topology transition $transition"))
end

function _three_phase_vsc_guard(device, voltage, transition)
    v, i = _three_phase_vsc_device_current(device, voltage)
    transition === :forward_extinction && return device.holding_current - i
    transition === :reverse_diode_extinction &&
        return device.antiparallel_diode.holding_current_a + i
    transition === :forward_turn_on &&
        return v - Nonlinear.power_semiconductor_forward_turn_on_voltage(device)
    return -v - Nonlinear.power_semiconductor_reverse_turn_on_voltage(device)
end

function _three_phase_vsc_incoming_guards(runtime::ThreePhaseVSCRuntime)
    return [
        (; phase, position, device, transition,
           left = _three_phase_vsc_guard(device, runtime.system.v, transition))
        for (phase, bridge) in enumerate(runtime.bridges), position in (:upper, :lower)
        for device in (Nonlinear.power_semiconductor_bridge_switch(bridge, position),)
        for transition in (device.closed ? (:forward_extinction,) :
            device.reverse_diode_conducting ? (:reverse_diode_extinction,) :
            device.gate_driver.applied_on ? (:forward_turn_on, :reverse_diode_turn_on) :
            (:reverse_diode_turn_on,))
    ]
end

function _solve_three_phase_vsc_causal_segment!(runtime, time_s, dt_s, matrix, rhs, solution)
    system = runtime.system
    capacitor = runtime.dc_link_capacitor
    fill!(matrix, 0.0)
    fill!(rhs, 0.0)
    y = @view matrix[1:system.node_count, 1:system.node_count]
    nodal_rhs = @view rhs[1:system.node_count]
    for element in system.elements
        element === capacitor || Branches.stamp!(y, nodal_rhs, element, time_s, dt_s)
    end
    # Keep capacitor current as an unknown. The equivalent voltage equation
    # avoids the very large 2C/dt Norton coefficient in short event probes.
    slot = length(rhs)
    matrix[capacitor.a, slot] = matrix[slot, capacitor.a] = 1.0
    matrix[capacitor.b, slot] = matrix[slot, capacitor.b] = -1.0
    scale = dt_s / (2.0 * capacitor.c)
    matrix[slot, slot] = -scale
    rhs[slot] = muladd(scale, capacitor.i_prev, capacitor.v_prev)
    copyto!(system.y, y)
    copyto!(system.rhs, nodal_rhs)
    Nodal.solve_dense!(solution, matrix, rhs)
    copyto!(system.v, 1, solution, 1, system.node_count)
    Branches.stamp_history_current!(system.rhs, capacitor.a, capacitor.b, solution[slot])
    return solution[slot]
end

function _accept_three_phase_vsc_causal_segment!(runtime, start_s, end_s, capacitor_i, intervals)
    start_v = intervals === nothing ? nothing : map(filter -> filter.v_prev, runtime.phase_filters)
    dt_s = end_s - start_s
    system = runtime.system
    capacitor = runtime.dc_link_capacitor
    for element in system.elements
        element === capacitor || Branches.update!(element, system.v, dt_s)
    end
    capacitor.v_prev = muladd(dt_s / (2.0 * capacitor.c),
        capacitor.i_prev + capacitor_i, capacitor.v_prev)
    capacitor.i_prev = capacitor.i_last = capacitor_i
    if start_v !== nothing
        push!(intervals, ThreePhaseVSCIntervalVoltage(
            start_s, end_s, start_v, map(filter -> filter.v_prev, runtime.phase_filters),
        ))
    end
    return nothing
end

function _three_phase_vsc_guard_bounds(root, bracket_start_s)
    upper_s = root.time_s
    # Width and endpoint are Float64. Round the reconstructed lower
    # endpoint outwards so cancellation cannot overstate separation.
    lower_s = max(bracket_start_s, prevfloat(upper_s - nextfloat(root.bracket_width_s)))
    return (; lower_s, upper_s)
end

function _unique_three_phase_vsc_guard(candidates, tolerance_s)
    isempty(candidates) && return nothing
    winner = argmin(candidate -> candidate.upper_s, candidates)
    for rival in candidates
        rival.index == winner.index && continue
        # The upward-rounded threshold gives a conservative strict comparison.
        rival.lower_s > nextfloat(winner.upper_s + tolerance_s) || throw(ArgumentError(
            "converter intrinsic guard ordering is unresolved at t=$(winner.upper_s) s",
        ))
    end
    return winner
end

function _advance_three_phase_vsc_causal_storage!(
    runtime, start_s::Float64, end_s::Float64, source_end_s::Float64, intervals,
)
    hybrid = OVER16TimestepIntegration
    policy = hybrid.HybridEventPolicy()
    cursor_s = start_s
    event_count = 0
    size = runtime.system.node_count + 1
    matrix = zeros(size, size)
    rhs = zeros(size)
    solution = zeros(size)
    while cursor_s < end_s
        guards = _three_phase_vsc_incoming_guards(runtime)
        solve(t) = _solve_three_phase_vsc_causal_segment!(
            runtime, t == end_s ? source_end_s : t, t - cursor_s, matrix, rhs, solution,
        )
        solve(end_s)
        # Capture all endpoint values before any root probe changes the trial arrays.
        right = map(guards) do guard
            _three_phase_vsc_guard(guard.device, runtime.system.v, guard.transition)
        end
        candidates = NamedTuple{(:index, :lower_s, :upper_s),Tuple{Int,Float64,Float64}}[]
        for (index, (guard, value)) in enumerate(zip(guards, right))
            value > 1.0e-9 || continue
            bounds = if guard.left >= 0.0
                guard.left <= 1.0e-9 || throw(ArgumentError(
                    "infeasible converter guard at the incoming boundary for $(guard.transition)",
                ))
                (; lower_s = cursor_s, upper_s = cursor_s)
            else
                bracket = hybrid.HybridEventBracket(cursor_s, guard.left, end_s, value)
                root = hybrid.localize_hybrid_event!(
                    hybrid.HYBRID_EVENT_RISING, bracket, policy,
                ) do probe_s
                    solve(probe_s)
                    _three_phase_vsc_guard(guard.device, runtime.system.v, guard.transition)
                end
                root.converged || throw(ArgumentError("converter guard root did not converge"))
                _three_phase_vsc_guard_bounds(root, cursor_s)
            end
            push!(candidates, (; index, bounds...))
        end
        winner = _unique_three_phase_vsc_guard(candidates, policy.simultaneity_tolerance_s)
        if winner === nothing
            capacitor_i = solve(end_s)
            _accept_three_phase_vsc_causal_segment!(runtime, cursor_s, end_s, capacitor_i, intervals)
            return runtime
        end
        guard = guards[winner.index]
        earliest_s = winner.upper_s
        if earliest_s > cursor_s
            earliest_s - cursor_s >= policy.minimum_progress_s || throw(ArgumentError(
                "converter intrinsic guard spacing is below the minimum progression interval",
            ))
            capacitor_i = solve(earliest_s)
            _accept_three_phase_vsc_causal_segment!(runtime, cursor_s, earliest_s, capacitor_i, intervals)
        else
            _solve_three_phase_vsc_boundary!(runtime, cursor_s)
        end
        residual = _three_phase_vsc_guard(guard.device, runtime.system.v, guard.transition)
        _apply_three_phase_vsc_topology_action!(guard.device, guard.transition, earliest_s)
        tick = round(Int, earliest_s / runtime.parameters.scheduler_tick_s)
        phase_name = Char(Int('a') + guard.phase - 1)
        push!(runtime.commutation_occurrences, ThreePhaseVSCCommutationOccurrence(
            Symbol(:phase_, phase_name, :_, guard.position, :_, guard.transition),
            earliest_s, tick, guard.phase, guard.position, guard.transition, residual, 1,
        ))
        evaluation_s = earliest_s == end_s ? source_end_s : earliest_s
        _stage_three_phase_vsc_boundary!(runtime, earliest_s,
            _solve_three_phase_vsc_boundary!(runtime, evaluation_s))
        cursor_s = earliest_s
        event_count += 1
        event_count <= policy.maximum_events_per_step || throw(ArgumentError(
            "converter intrinsic events did not re-arm with finite progress",
        ))
    end
    return runtime
end

function _stabilize_three_phase_vsc_topology!(
    solve_state!,
    runtime::ThreePhaseVSCRuntime,
    time_s::Float64;
    maximum_iterations::Int = 32,
    record_actions::Bool = true,
)
    maximum_iterations > 0 || throw(ArgumentError(
        "converter topology maximum iterations must be positive",
    ))
    for iteration in 1:maximum_iterations
        solve_state!()
        selected = nothing
        for phase in 1:3, position in (:upper, :lower)
            device = Nonlinear.power_semiconductor_bridge_switch(
                runtime.bridges[phase],
                position,
            )
            action = _three_phase_vsc_topology_action(device, runtime.system.v)
            action === nothing && continue
            selected = (phase, position, device, action)
            break
        end
        if selected === nothing
            runtime.maximum_topology_iterations = max(
                runtime.maximum_topology_iterations,
                iteration,
            )
            return iteration
        end
        phase, position, device, action = selected
        _apply_three_phase_vsc_topology_action!(
            device,
            action.transition,
            time_s,
        )
        tick = round(Int, time_s / runtime.parameters.scheduler_tick_s)
        record_actions && push!(
            runtime.commutation_occurrences,
            ThreePhaseVSCCommutationOccurrence(
                Symbol(
                    :phase_,
                    Char(Int('a') + phase - 1),
                    :_,
                    position,
                    :_,
                    action.transition,
                ),
                time_s,
                tick,
                phase,
                position,
                action.transition,
                action.residual,
                iteration,
            ),
        )
    end
    throw(ArgumentError(
        "converter semiconductor topology did not stabilize within $maximum_iterations algebraic solves at t=$time_s s",
    ))
end


function _stabilize_three_phase_vsc_topology!(
    runtime::ThreePhaseVSCRuntime, time_s::Float64, dt_s::Float64;
    source_time_s::Float64 = time_s,
)
    return _stabilize_three_phase_vsc_topology!(runtime, time_s) do
        Nodal.solve_algebraic_state!(runtime.system, source_time_s, dt_s)
    end
end

function _stabilize_historical_three_phase_vsc_topology!(
    runtime::ThreePhaseVSCRuntime, time_s::Float64, dt_s::Float64,
)
    for bridge in runtime.bridges
        Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(bridge, time_s)
    end
    return _stabilize_three_phase_vsc_topology!(runtime, time_s, dt_s)
end

function _solve_three_phase_vsc_boundary!(runtime::ThreePhaseVSCRuntime, time_s::Float64)
    system = runtime.system
    p = runtime.parameters
    positive = _VSC_DC_POSITIVE_NODE
    negative = _VSC_DC_NEGATIVE_NODE
    y, rhs, v = system.y, system.rhs, system.v
    fill!(y, 0.0)
    fill!(rhs, 0.0)
    for bridge in runtime.bridges
        Branches.stamp!(y, rhs, bridge, time_s, 0.0)
    end
    dc_g = inv(p.dc_source_resistance_ohm)
    Branches.stamp_conductance!(y, positive, negative, dc_g)
    Branches.stamp_history_current!(rhs, positive, negative, -dc_g * p.dc_source_voltage_v)

    # Use the negative rail as a temporary voltage reference. With storage
    # held fixed, the pole KCL equations are scalar and the DC link voltage
    # is known. Recover common mode from d(sum(i_phase))/dt = 0.
    v[negative] = 0.0
    v[positive] = runtime.dc_link_capacitor.v_prev
    common_mode_sum_v = 0.0
    for (phase, filter) in enumerate(runtime.phase_filters)
        pole, terminal = filter.a, filter.b
        current_a = filter.i_prev
        source_v = ThreePhaseGridVoltageSignal(p, phase)(time_s)
        grid_g = inv(p.grid_source_resistance_ohm)
        Branches.stamp_conductance!(y, terminal, 0, grid_g)
        Branches.stamp_history_current!(rhs, terminal, 0, -grid_g * source_v)
        Branches.stamp_history_current!(rhs, pole, terminal, current_a)
        v[terminal] = source_v + p.grid_source_resistance_ohm * current_a
        if y[pole, pole] > 0.0
            v[pole] = (rhs[pole] - y[pole, positive] * v[positive]) / y[pole, pole]
        else
            # With ideal off switches a disconnected pole needs a nonzero
            # filter current to select its diode. At zero current keep its
            # previous voltage within the two nonconducting thresholds.
            diode_drop = p.diode_forward_voltage_v
            previous_pole_v = runtime.bridges[phase].last_ac_terminal_voltage_v -
                runtime.bridges[phase].last_dc_negative_voltage_v
            v[pole] = current_a > 0 ? -diode_drop - 1.0 :
                current_a < 0 ? v[positive] + diode_drop + 1.0 :
                clamp(previous_pole_v, -diode_drop, v[positive] + diode_drop)
        end
        common_mode_sum_v += v[terminal] + filter.r * current_a - v[pole]
    end
    common_mode_v = common_mode_sum_v / length(runtime.phase_filters)
    for node in (positive, negative, _three_phase_vsc_pole_node.(1:3)...)
        v[node] += common_mode_v
    end
    # The fixed-voltage capacitor supplies the current required by rail KCL.
    return rhs[positive] - sum(y[positive, node] * v[node] for node in eachindex(v))
end

function _resolve_three_phase_vsc_boundary!(
    runtime::ThreePhaseVSCRuntime, time_s::Float64;
    source_time_s::Float64 = time_s,
)
    capacitor_current_a = 0.0
    _stabilize_three_phase_vsc_topology!(runtime, time_s; record_actions = false) do
        capacitor_current_a = _solve_three_phase_vsc_boundary!(runtime, source_time_s)
    end
    return _stage_three_phase_vsc_boundary!(runtime, time_s, capacitor_current_a)
end

function _stage_three_phase_vsc_boundary!(runtime, time_s, capacitor_current_a)
    system = runtime.system
    capacitor = runtime.dc_link_capacitor
    Branches.stamp_history_current!(system.rhs, capacitor.a, capacitor.b, capacitor_current_a)
    capacitor.i_prev = capacitor_current_a
    capacitor.i_last = capacitor_current_a
    for filter in runtime.phase_filters
        filter.v_prev = Branches.branch_voltage(system.v, filter.a, filter.b)
    end
    for bridge in runtime.bridges
        Branches.update!(bridge, system.v, 0.0)
    end
    return runtime
end

function _three_phase_vsc_source_quantities(runtime::ThreePhaseVSCRuntime, time_s::Float64)
    p = runtime.parameters
    system = runtime.system
    dc_link_voltage = system.v[_VSC_DC_POSITIVE_NODE] -
        system.v[_VSC_DC_NEGATIVE_NODE]
    dc_source_current = (p.dc_source_voltage_v - dc_link_voltage) /
        p.dc_source_resistance_ohm
    dc_external_power = p.dc_source_voltage_v * dc_source_current
    dc_resistor_loss = p.dc_source_resistance_ohm * dc_source_current^2
    grid_external_power = 0.0
    grid_resistor_loss = 0.0
    for phase in 1:3
        source_voltage = ThreePhaseGridVoltageSignal(p, phase)(time_s)
        source_current = (
            source_voltage -
            system.v[_three_phase_vsc_grid_terminal_node(phase)]
        ) /
            p.grid_source_resistance_ohm
        grid_external_power += source_voltage * source_current
        grid_resistor_loss += p.grid_source_resistance_ohm * source_current^2
    end
    return (
        dc_source_current_a = dc_source_current,
        external_power_w = dc_external_power + grid_external_power,
        source_resistor_loss_w = dc_resistor_loss + grid_resistor_loss,
    )
end

function _three_phase_vsc_stored_energy(runtime::ThreePhaseVSCRuntime)
    p = runtime.parameters
    dc_link_voltage = runtime.system.v[_VSC_DC_POSITIVE_NODE] -
        runtime.system.v[_VSC_DC_NEGATIVE_NODE]
    dc_energy = 0.5 * p.dc_link_capacitance_f * dc_link_voltage^2
    filter_inductance = p.filter_inductance_h + p.transformer_leakage_inductance_h
    filter_energy = 0.5 * filter_inductance * sum(
        filter.i_last^2 for filter in runtime.phase_filters
    )
    bridge_energy = sum(
        Nonlinear.power_semiconductor_bridge_terminal_state(bridge).stored_energy_j
        for bridge in runtime.bridges
    )
    return dc_energy + filter_energy + bridge_energy
end

function _record_three_phase_vsc_sample!(
    recorder::ThreePhaseVSCRecorder,
    runtime::ThreePhaseVSCRuntime,
    time_s::Float64,
)
    sample = recorder.write_index + 1
    sample <= length(recorder.time_s) || throw(ArgumentError(
        "converter trace storage is exhausted",
    ))
    system = runtime.system
    p = runtime.parameters
    dc_negative_voltage = system.v[_VSC_DC_NEGATIVE_NODE]
    pole = ntuple(
        phase -> system.v[_three_phase_vsc_pole_node(phase)] - dc_negative_voltage,
        3,
    )
    grid = ntuple(
        phase -> system.v[_three_phase_vsc_grid_terminal_node(phase)],
        3,
    )
    current = ntuple(phase -> runtime.phase_filters[phase].i_last, 3)
    line = (pole[1] - pole[2], pole[2] - pole[3], pole[3] - pole[1])
    angle = 2.0 * pi * p.frequency_hz * time_s
    power = SwitchDetailedVSC.instantaneous_three_phase_power(grid, current, angle)
    source = _three_phase_vsc_source_quantities(runtime, time_s)
    bridge_states = ntuple(
        phase -> Nonlinear.power_semiconductor_bridge_terminal_state(
            runtime.bridges[phase],
        ),
        3,
    )
    semiconductor_loss = sum(
        state.semiconductor_loss_w + state.snubber_resistor_loss_w
        for state in bridge_states
    )
    series_resistance = p.filter_resistance_ohm +
        p.transformer_leakage_resistance_ohm
    filter_resistor_loss = series_resistance * sum(abs2, current)
    kcl_residual = system.y * system.v - system.rhs
    recorder.time_s[sample] = time_s
    for phase in 1:3
        bridge = bridge_states[phase]
        upper = bridge.upper_switch
        lower = bridge.lower_switch
        recorder.pole_voltage_v[phase, sample] = pole[phase]
        recorder.line_voltage_v[phase, sample] = line[phase]
        recorder.grid_terminal_voltage_v[phase, sample] = grid[phase]
        recorder.filter_current_a[phase, sample] = current[phase]
        recorder.duty[phase, sample] = runtime.controller.held_duties[phase]
        recorder.upper_gate_applied[phase, sample] = upper.gate_applied_on ? 1.0 : 0.0
        recorder.lower_gate_applied[phase, sample] = lower.gate_applied_on ? 1.0 : 0.0
        recorder.upper_forward_conducting[phase, sample] = upper.forward_conducting ? 1.0 : 0.0
        recorder.lower_forward_conducting[phase, sample] = lower.forward_conducting ? 1.0 : 0.0
        recorder.upper_antiparallel_diode_conducting[phase, sample] =
            upper.reverse_diode_conducting ? 1.0 : 0.0
        recorder.lower_antiparallel_diode_conducting[phase, sample] =
            lower.reverse_diode_conducting ? 1.0 : 0.0
    end
    recorder.dc_link_voltage_v[sample] =
        system.v[_VSC_DC_POSITIVE_NODE] - dc_negative_voltage
    recorder.dc_source_current_a[sample] = source.dc_source_current_a
    recorder.active_power_w[sample] = power.active_w
    recorder.reactive_power_var[sample] = power.reactive_var
    recorder.semiconductor_loss_w[sample] = semiconductor_loss
    recorder.resistive_loss_w[sample] =
        source.source_resistor_loss_w + filter_resistor_loss
    recorder.stored_energy_j[sample] = _three_phase_vsc_stored_energy(runtime)
    recorder.external_source_power_w[sample] = source.external_power_w
    recorder.nodal_kcl_residual_a[sample] = maximum(abs, kcl_residual; init = 0.0)
    recorder.blocked[sample] = all(bridge.blocked for bridge in runtime.bridges) ? 1.0 : 0.0
    recorder.write_index = sample
    return recorder
end

"""
Construct the converter integrator. The default preserves the historical
sampled topology execution. Set `causal_events=true` to localize directed
intrinsic guards between scheduled boundaries and preserve that selection
across checkpoint/restart.
"""
function configure_three_phase_vsc(
    parameters::SwitchDetailedVSC.ThreePhaseTwoLevelVSCParameters =
        SwitchDetailedVSC.ThreePhaseTwoLevelVSCParameters();
    causal_events::Bool = false,
)
    runtime = _three_phase_vsc_runtime(parameters)
    scheduler, controller, pwm_tasks = _three_phase_vsc_scheduler(runtime)
    sample_count = round(Int, parameters.end_time_s / parameters.timestep_s) + 1
    integrator = ThreePhaseVSCIntegrator(
        runtime,
        scheduler,
        TimestepTransaction(runtime),
        controller,
        pwm_tasks,
        ThreePhaseVSCRecorder(sample_count),
        0,
        false,
        false,
        nothing,
    )
    task_states = _three_phase_vsc_task_states(integrator)
    occurrence_count = length(scheduler.occurrences)
    begin_timestep_transaction!(integrator.transaction)
    try
        if causal_events
            _apply_three_phase_vsc_boundaries!(integrator, 0.0)
            run_due_sampled_tasks!(scheduler, runtime, 0.0)
            for bridge in runtime.bridges
                Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(bridge, 0.0)
            end
            _resolve_three_phase_vsc_boundary!(runtime, 0.0)
        else
            run_due_sampled_tasks!(scheduler, runtime, 0.0)
            _stabilize_historical_three_phase_vsc_topology!(runtime, 0.0, parameters.timestep_s)
        end
        commit_timestep_transaction!(integrator.transaction)
    catch error
        restore_timestep_transaction!(integrator.transaction)
        _restore_three_phase_vsc_task_states!(integrator, task_states, occurrence_count)
        rethrow(error)
    end
    _record_three_phase_vsc_sample!(integrator.recorder, runtime, 0.0)
    return causal_events ? CausalThreePhaseVSCIntegrator(integrator) : integrator
end

"""Advance exactly one accepted converter timestep with atomic scheduler, protection, network, companion, and device mutation."""
advance_three_phase_vsc!(integrator::ThreePhaseVSCIntegrator) =
    _advance_three_phase_vsc!(integrator, false)

function advance_three_phase_vsc!(
    integrator::CausalThreePhaseVSCIntegrator; accepted_intervals = nothing,
)
    return _advance_three_phase_vsc!(integrator.state, true; accepted_intervals)
end

function _advance_three_phase_vsc!(
    integrator::ThreePhaseVSCIntegrator, causal_events::Bool;
    accepted_intervals::Union{Nothing,Vector{ThreePhaseVSCIntervalVoltage}} = nothing,
)
    integrator.failed && throw(ArgumentError(
        "three-phase VSC integrator is terminally failed: $(integrator.last_failure)",
    ))
    integrator.completed && return false
    p = integrator.runtime.parameters
    final_step = round(Int, p.end_time_s / p.timestep_s)
    integrator.accepted_step_index >= final_step && begin
        integrator.completed = true
        return false
    end
    next_step = integrator.accepted_step_index + 1
    time_s = next_step * p.timestep_s
    task_states = _three_phase_vsc_task_states(integrator)
    occurrence_count = length(integrator.scheduler.occurrences)
    recorder_write_index = integrator.recorder.write_index
    intervals = accepted_intervals === nothing ? nothing : ThreePhaseVSCIntervalVoltage[]
    begin_timestep_transaction!(integrator.transaction)
    try
        start_s = (next_step - 1) * p.timestep_s
        if causal_events
            _advance_three_phase_vsc_network_window!(integrator, start_s, time_s, intervals)
        else
            _advance_historical_three_phase_vsc_window!(integrator, start_s, time_s)
        end
        _record_three_phase_vsc_sample!(
            integrator.recorder,
            integrator.runtime,
            time_s,
        )
        commit_timestep_transaction!(integrator.transaction)
    catch error
        timestep_transaction_active(integrator.transaction) &&
            restore_timestep_transaction!(integrator.transaction)
        _restore_three_phase_vsc_task_states!(integrator, task_states, occurrence_count)
        integrator.recorder.write_index = recorder_write_index
        integrator.failed = true
        integrator.last_failure = sprint(showerror, error)
        rethrow(error)
    end
    integrator.accepted_step_index = next_step
    integrator.completed = next_step == final_step
    accepted_intervals === nothing || append!(accepted_intervals, intervals)
    return true
end

function run_three_phase_vsc!(
    integrator::Union{ThreePhaseVSCIntegrator,CausalThreePhaseVSCIntegrator};
    stop_time_s::Union{Nothing,Real} = nothing,
)
    state = _three_phase_vsc_state(integrator)
    stop_step = stop_time_s === nothing ? typemax(Int) : begin
        stop_time = Float64(stop_time_s)
        isfinite(stop_time) && stop_time >= 0.0 || throw(ArgumentError(
            "converter stop time must be finite and nonnegative",
        ))
        round(Int, stop_time / state.runtime.parameters.timestep_s)
    end
    while !state.completed && state.accepted_step_index < stop_step
        advance_three_phase_vsc!(integrator)
    end
    return integrator
end

function _integrate_trapezoidal(values::AbstractVector{Float64}, dt_s::Float64)
    length(values) >= 2 || return 0.0
    total = 0.0
    @inbounds for index in 2:length(values)
        total += 0.5 * dt_s * (values[index - 1] + values[index])
    end
    return total
end

function _three_phase_vsc_converter_dissipated_energy(
    recorder::ThreePhaseVSCRecorder,
    parameters::SwitchDetailedVSC.ThreePhaseTwoLevelVSCParameters,
)
    length(recorder.time_s) >= 2 || return 0.0
    series_resistance = parameters.filter_resistance_ohm +
        parameters.transformer_leakage_resistance_ohm
    function converter_loss_power(sample::Int)
        filter_loss_w = series_resistance * sum(
            abs2,
            @view(recorder.filter_current_a[:, sample]),
        )
        dc_source_loss_w = parameters.dc_source_resistance_ohm *
            recorder.dc_source_current_a[sample]^2
        return recorder.semiconductor_loss_w[sample] + filter_loss_w +
            dc_source_loss_w
    end
    energy_j = 0.0
    previous_loss_w = converter_loss_power(1)
    @inbounds for sample in 2:length(recorder.time_s)
        current_loss_w = converter_loss_power(sample)
        energy_j += 0.5 * parameters.timestep_s *
            (previous_loss_w + current_loss_w)
        previous_loss_w = current_loss_w
    end
    return energy_j
end

function _synchronous_harmonics(
    values::AbstractVector{Float64},
    time_s::AbstractVector{Float64},
    frequency_hz::Float64,
    maximum_order::Int,
)
    sample_count = length(values)
    sample_count == length(time_s) || throw(DimensionMismatch(
        "harmonic values and timestamps must have equal lengths",
    ))
    sample_count >= 4 || throw(ArgumentError("harmonic window has too few samples"))
    rms = Vector{Float64}(undef, maximum_order)
    for order in 1:maximum_order
        cosine_sum = 0.0
        sine_sum = 0.0
        angular_frequency = 2.0 * pi * order * frequency_hz
        @inbounds for index in eachindex(values, time_s)
            angle = angular_frequency * time_s[index]
            cosine_sum += values[index] * cos(angle)
            sine_sum += values[index] * sin(angle)
        end
        scale = 2.0 / sample_count
        rms[order] = hypot(scale * cosine_sum, scale * sine_sum) / sqrt(2.0)
    end
    fundamental = rms[1]
    thd = fundamental > eps(Float64) ?
        sqrt(sum(abs2, @view rms[2:end])) / fundamental : Inf
    return fundamental, thd
end

function _three_phase_vsc_harmonic_metadata(
    recorder::ThreePhaseVSCRecorder,
    parameters::SwitchDetailedVSC.ThreePhaseTwoLevelVSCParameters,
)
    indices = findall(time ->
        parameters.harmonic_window_start_s <= time < parameters.harmonic_window_end_s,
        recorder.time_s,
    )
    times = recorder.time_s[indices]
    maximum_order = min(parameters.maximum_harmonic_order, div(length(indices), 2) - 1)
    maximum_order >= 2 || throw(ArgumentError("harmonic window cannot resolve two harmonics"))
    current = ntuple(3) do phase
        _synchronous_harmonics(
            recorder.filter_current_a[phase, indices],
            times,
            parameters.frequency_hz,
            maximum_order,
        )
    end
    voltage = ntuple(3) do line
        _synchronous_harmonics(
            recorder.line_voltage_v[line, indices],
            times,
            parameters.frequency_hz,
            maximum_order,
        )
    end
    return ThreePhaseVSCHarmonicMetadata(
        parameters.harmonic_window_start_s,
        parameters.harmonic_window_end_s,
        length(indices),
        parameters.frequency_hz,
        maximum_order,
        ntuple(phase -> current[phase][1], 3),
        ntuple(phase -> current[phase][2], 3),
        ntuple(line -> voltage[line][1], 3),
        ntuple(line -> voltage[line][2], 3),
    )
end

function _three_phase_vsc_metrics(integrator::ThreePhaseVSCIntegrator)
    integrator.completed || throw(ArgumentError(
        "converter metrics require a completed horizon",
    ))
    runtime = integrator.runtime
    recorder = integrator.recorder
    p = runtime.parameters
    harmonic = _three_phase_vsc_harmonic_metadata(recorder, p)
    window_indices = findall(time ->
        p.harmonic_window_start_s <= time < p.harmonic_window_end_s,
        recorder.time_s,
    )
    external_energy = _integrate_trapezoidal(
        recorder.external_source_power_w,
        p.timestep_s,
    )
    dc_source_energy = p.dc_source_voltage_v * _integrate_trapezoidal(
        recorder.dc_source_current_a,
        p.timestep_s,
    )
    ac_terminal_energy = _integrate_trapezoidal(
        recorder.active_power_w,
        p.timestep_s,
    )
    converter_dissipated_energy = _three_phase_vsc_converter_dissipated_energy(
        recorder,
        p,
    )
    dissipated_power = recorder.semiconductor_loss_w .+ recorder.resistive_loss_w
    dissipated_energy = _integrate_trapezoidal(dissipated_power, p.timestep_s)
    stored_change = recorder.stored_energy_j[end] - recorder.stored_energy_j[1]
    dc_ac_energy_residual = dc_source_energy - ac_terminal_energy -
        converter_dissipated_energy - stored_change
    dc_ac_energy_scale = max(
        abs(dc_source_energy),
        abs(ac_terminal_energy) + converter_dissipated_energy + abs(stored_change),
        eps(Float64),
    )
    energy_residual = external_energy - dissipated_energy - stored_change
    energy_scale = max(abs(external_energy), abs(dissipated_energy) + abs(stored_change), eps(Float64))
    all_off_count = count(1:length(recorder.time_s)) do sample
        any(phase ->
            recorder.upper_gate_applied[phase, sample] == 0.0 &&
            recorder.lower_gate_applied[phase, sample] == 0.0,
            1:3,
        )
    end
    diode_count = count(1:length(recorder.time_s)) do sample
        any(phase ->
            recorder.upper_antiparallel_diode_conducting[phase, sample] > 0.5 ||
            recorder.lower_antiparallel_diode_conducting[phase, sample] > 0.5,
            1:3,
        )
    end
    finite_output = all(isfinite, recorder.pole_voltage_v) &&
        all(isfinite, recorder.line_voltage_v) &&
        all(isfinite, recorder.grid_terminal_voltage_v) &&
        all(isfinite, recorder.filter_current_a) &&
        all(isfinite, recorder.dc_link_voltage_v) &&
        all(isfinite, recorder.active_power_w) &&
        all(isfinite, recorder.reactive_power_var)
    exact_alignment = all(runtime.boundary_occurrences) do occurrence
        occurrence.time_s == occurrence.tick * p.scheduler_tick_s
    end && all(runtime.commutation_occurrences) do occurrence
        occurrence.time_s == occurrence.tick * p.scheduler_tick_s
    end
    return ThreePhaseVSCMetrics(
        sum(recorder.active_power_w[window_indices]) / length(window_indices),
        sum(recorder.reactive_power_var[window_indices]) / length(window_indices),
        sum(recorder.dc_link_voltage_v[window_indices]) / length(window_indices),
        maximum(recorder.dc_link_voltage_v[window_indices]) -
            minimum(recorder.dc_link_voltage_v[window_indices]),
        maximum(recorder.nodal_kcl_residual_a),
        dc_source_energy,
        ac_terminal_energy,
        converter_dissipated_energy,
        dc_ac_energy_residual,
        abs(dc_ac_energy_residual) / dc_ac_energy_scale,
        external_energy,
        dissipated_energy,
        stored_change,
        energy_residual,
        abs(energy_residual) / energy_scale,
        integrator.controller_task.sample_count,
        integrator.controller_task.write_count,
        ntuple(phase -> integrator.pwm_tasks[phase].cycle_count, 3),
        ntuple(phase -> integrator.pwm_tasks[phase].edge_count, 3),
        ntuple(
            phase -> Nonlinear.power_semiconductor_bridge_topology_transition_count(
                runtime.bridges[phase],
            ),
            3,
        ),
        sum(bridge.shoot_through_rejection_count for bridge in runtime.bridges),
        all_off_count,
        diode_count,
        runtime.block_count,
        runtime.restart_count,
        length(runtime.boundary_occurrences),
        length(runtime.commutation_occurrences),
        runtime.maximum_topology_iterations,
        exact_alignment,
        finite_output,
        harmonic,
    )
end

three_phase_vsc_trace(integrator::CausalThreePhaseVSCIntegrator) =
    three_phase_vsc_trace(integrator.state)

function three_phase_vsc_trace(integrator::ThreePhaseVSCIntegrator)
    integrator.completed || throw(ArgumentError(
        "converter trace requires a completed horizon",
    ))
    recorder = integrator.recorder
    recorder.write_index == length(recorder.time_s) || throw(ArgumentError(
        "converter trace storage is incomplete",
    ))
    metrics = _three_phase_vsc_metrics(integrator)
    return ThreePhaseVSCTrace(
        integrator.runtime.parameters,
        copy(recorder.time_s),
        copy(recorder.pole_voltage_v),
        copy(recorder.line_voltage_v),
        copy(recorder.grid_terminal_voltage_v),
        copy(recorder.filter_current_a),
        copy(recorder.dc_link_voltage_v),
        copy(recorder.dc_source_current_a),
        copy(recorder.active_power_w),
        copy(recorder.reactive_power_var),
        copy(recorder.duty),
        copy(recorder.upper_gate_applied),
        copy(recorder.lower_gate_applied),
        copy(recorder.upper_forward_conducting),
        copy(recorder.lower_forward_conducting),
        copy(recorder.upper_antiparallel_diode_conducting),
        copy(recorder.lower_antiparallel_diode_conducting),
        copy(recorder.semiconductor_loss_w),
        copy(recorder.resistive_loss_w),
        copy(recorder.stored_energy_j),
        copy(recorder.external_source_power_w),
        copy(recorder.nodal_kcl_residual_a),
        copy(recorder.blocked),
        copy(integrator.runtime.boundary_occurrences),
        copy(integrator.runtime.commutation_occurrences),
        copy(integrator.scheduler.occurrences),
        ntuple(phase -> Int[getfield(edge, :tick) for edge in integrator.pwm_tasks[phase].occurrences], 3),
        metrics,
    )
end

function simulate_three_phase_vsc(
    parameters::SwitchDetailedVSC.ThreePhaseTwoLevelVSCParameters =
        SwitchDetailedVSC.ThreePhaseTwoLevelVSCParameters();
    causal_events::Bool = false,
)
    integrator = configure_three_phase_vsc(parameters; causal_events)
    run_three_phase_vsc!(integrator)
    return three_phase_vsc_trace(integrator)
end

struct ThreePhaseVSCCheckpointPayload{I}
    schema::UInt32
    julia_version::VersionNumber
    integrator::I
end

"""Write an integrity-checked, trusted Julia checkpoint for a non-active converter integrator."""
function write_three_phase_vsc_checkpoint(
    path::AbstractString,
    integrator::Union{ThreePhaseVSCIntegrator,CausalThreePhaseVSCIntegrator},
)
    timestep_transaction_active(_three_phase_vsc_state(integrator).transaction) && throw(ArgumentError(
        "converter checkpoint cannot be written during an active timestep transaction",
    ))
    payload_io = IOBuffer()
    serialize(
        payload_io,
        ThreePhaseVSCCheckpointPayload(
            _THREE_PHASE_VSC_CHECKPOINT_SCHEMA,
            VERSION,
            integrator,
        ),
    )
    payload = take!(payload_io)
    digest = sha256(payload)
    output_path = abspath(String(path))
    mkpath(dirname(output_path))
    temporary_path, io = mktemp(dirname(output_path))
    try
        write(io, _THREE_PHASE_VSC_CHECKPOINT_MAGIC)
        _write_checkpoint_integer(io, UInt64(length(payload)))
        write(io, digest)
        write(io, payload)
        close(io)
        mv(temporary_path, output_path; force = true)
    catch
        isopen(io) && close(io)
        isfile(temporary_path) && rm(temporary_path; force = true)
        rethrow()
    end
    return output_path
end

"""Read and integrity-check a trusted Julia converter checkpoint in a compatible Julia process."""
function read_three_phase_vsc_checkpoint(
    path::AbstractString;
    max_payload_bytes::Integer = 2_000_000_000,
)
    input_path = abspath(String(path))
    isfile(input_path) || throw(ArgumentError(
        "converter checkpoint does not exist: $input_path",
    ))
    max_payload_bytes > 0 || throw(ArgumentError(
        "converter checkpoint payload limit must be positive",
    ))
    payload = open(input_path, "r") do io
        magic = read(io, length(_THREE_PHASE_VSC_CHECKPOINT_MAGIC))
        magic == _THREE_PHASE_VSC_CHECKPOINT_MAGIC || throw(ArgumentError(
            "file is not an AIMORA three-phase VSC checkpoint",
        ))
        length_value = _read_checkpoint_integer(io, UInt64)
        length_value <= UInt64(max_payload_bytes) || throw(ArgumentError(
            "converter checkpoint exceeds the configured payload limit",
        ))
        expected_digest = read(io, 32)
        bytes = read(io, Int(length_value))
        eof(io) || throw(ArgumentError("converter checkpoint contains trailing bytes"))
        sha256(bytes) == expected_digest || throw(ArgumentError(
            "converter checkpoint integrity digest does not match",
        ))
        return bytes
    end
    payload = deserialize(IOBuffer(payload))
    payload isa ThreePhaseVSCCheckpointPayload || throw(ArgumentError(
        "converter checkpoint payload has an unexpected root type",
    ))
    payload.schema == _THREE_PHASE_VSC_CHECKPOINT_SCHEMA || throw(ArgumentError(
        "converter checkpoint schema is unsupported",
    ))
    payload.julia_version == VERSION || throw(ArgumentError(
        "converter checkpoint Julia version $(payload.julia_version) does not match $VERSION",
    ))
    integrator = payload.integrator
    integrator isa Union{ThreePhaseVSCIntegrator,CausalThreePhaseVSCIntegrator} || throw(ArgumentError(
        "converter checkpoint does not contain a converter integrator",
    ))
    state = _three_phase_vsc_state(integrator)
    timestep_transaction_active(state.transaction) && throw(ArgumentError(
        "converter checkpoint contains an active timestep transaction",
    ))
    SwitchDetailedVSC.validate_three_phase_vsc_parameters(
        state.runtime.parameters,
    )
    return integrator
end

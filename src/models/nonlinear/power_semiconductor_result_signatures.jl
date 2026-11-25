using ..PortableSnapshots: PortableSnapshotRecord, _write_portable_value

const POWER_SEMICONDUCTOR_RESULT_SIGNATURE_SCHEMA = "aimora.semiconductor.result.v1"

"""The public identity of a device result for one committed network advance."""
struct PowerSemiconductorResultSignature
    schema_id::String
    device_index::Int
    interval_start_s::Float64
    interval_end_s::Float64
    sha256::String
end

"""Project the approved generic source and uncertainty declarations for public results."""
function power_semiconductor_public_provenance(provenance::Union{Nothing,ParameterProvenance})
    provenance === nothing && return (source="unavailable", uncertainty="unavailable")
    source = provenance.source in (
        "AIMORA-authored synthetic exact case input",
        "AIMORA generic caller-supplied semiconductor parameters",
    ) ? provenance.source : "undisclosed"
    uncertainty = provenance.uncertainty in (
        "unknown",
        "unknown; no manufacturer or measurement inference",
        "unknown unless narrowed by the caller's parameter record",
    ) ? provenance.uncertainty : "undisclosed"
    return (; source, uncertainty)
end

const _SEMICONDUCTOR_TERMINAL_RESULT_FIELDS = (
    :device_kind, :forward_conducting, :reverse_diode_conducting, :gate_commanded_on,
    :gate_applied_on, :terminal_voltage_v, :terminal_current_a, :forward_current_a,
    :reverse_diode_current_a, :snubber_current_a, :semiconductor_loss_w,
    :snubber_resistor_loss_w, :semiconductor_dissipated_energy_j, :snubber_dissipated_energy_j,
    :snubber_capacitor_energy_j,
)
const _SEMICONDUCTOR_EXTENDED_RESULT_FIELDS = (
    :schema_version, :fidelity_components, :conduction_direction, :gate_turn_off_disposition,
    :recovery_active, :stored_recovery_charge_c, :recovery_current_a, :peak_reverse_current_a,
    :cumulative_recovered_charge_c, :recovery_zero_event_count, :last_recovery_zero_time_s,
    :junction_capacitance_f, :junction_charge_c, :displacement_current_a, :junction_stored_energy_j,
    :companion_energy_residual_j, :tail_active, :tail_current_a, :tail_cutoff_event_count,
    :last_tail_cutoff_time_s, :last_event_kind, :last_event_energy_j, :cumulative_turn_on_energy_j,
    :cumulative_turn_off_energy_j, :cumulative_reverse_recovery_energy_j, :junction_temperature_k,
    :thermal_node_temperature_k, :ambient_heat_flow_w, :thermal_stored_energy_j,
)
const _SEMICONDUCTOR_EVENT_RESULT_FIELDS = (
    :device_index, :event_kind, :event_time_s, :accepted_side, :interval_start_s, :interval_end_s,
    :order, :current_magnitude_a, :blocking_voltage_magnitude_v, :pre_impulse_temperature_k,
    :energy_j,
)
const _SEMICONDUCTOR_RESOLUTION_RESULT_FIELDS = (
    :device_index, :component, :interval_start_s, :interval_end_s, :step_s, :interval_method,
    :publication_start_s, :publication_end_s, :law_evaluated, :participating, :participation,
    :scale_source, :scale_s, :scale_over_step, :meets_twenty_points,
)

_semiconductor_public_scalar(value) = value
_semiconductor_public_scalar(value::Symbol) = String(value)
_semiconductor_public_scalar(values::Tuple{Vararg{Symbol}}) = map(String, values)

function _semiconductor_public_record(schema, value, fields)
    return PortableSnapshotRecord(schema, [
        String(name) => _semiconductor_public_scalar(getfield(value, name)) for name in fields
    ])
end

function _semiconductor_extended_public_value(state, field)
    value = getfield(state, field)
    field === :last_recovery_zero_time_s && value == -Inf &&
        state.recovery_zero_event_count == 0 && return nothing
    field === :last_tail_cutoff_time_s && value == -Inf &&
        state.tail_cutoff_event_count == 0 && return nothing
    return _semiconductor_public_scalar(value)
end

function _semiconductor_provenance_record(provenance)
    provenance === nothing && return nothing
    public_values = power_semiconductor_public_provenance(provenance)
    return PortableSnapshotRecord("semiconductor.provenance.v1", [
        "source" => public_values.source, "uncertainty" => public_values.uncertainty,
    ])
end

function _semiconductor_event_margin_record(row)
    schema = "semiconductor.axis-margin.v1"
    fields = (:domain_lower, :domain_upper, :cell_lower, :cell_upper)
    return PortableSnapshotRecord("semiconductor.event-margin.v1", [
        "lookup_order" => row.lookup.order,
        "current_a" => _semiconductor_public_record(schema, row.current_a, fields),
        "blocking_voltage_v" => _semiconductor_public_record(schema, row.blocking_voltage_v, fields),
        "temperature_k" => _semiconductor_public_record(schema, row.temperature_k, fields),
    ])
end

"""
Canonical bytes of the public result for one device and interval.
Pass empty sequences for recorded outputs with no rows and nothing for unavailable history.
"""
function power_semiconductor_public_result_bytes(
    device_index::Int, terminals::NTuple{2,Int}, interval_start_s::Float64, interval_end_s::Float64,
    step_s::Float64, terminal::PowerSemiconductorTerminalState, extended::PowerSemiconductorExtendedResult;
    terminal_current_margins=nothing, events=nothing, event_margins=nothing, resolution=nothing,
)
    terminal_state = _semiconductor_public_record(
        "semiconductor.terminal.v1", terminal, _SEMICONDUCTOR_TERMINAL_RESULT_FIELDS,
    )
    state = extended.state
    endpoint = PortableSnapshotRecord("semiconductor.extended-state.v1", [
        String(field) => _semiconductor_extended_public_value(state, field)
        for field in _SEMICONDUCTOR_EXTENDED_RESULT_FIELDS
    ])
    components = extended.component_provenance
    component_provenance = components === nothing ? nothing :
        PortableSnapshotRecord("semiconductor.component-provenance.v1", [
            String(key) => _semiconductor_provenance_record(components[key]) for key in (
                :recovered_charge, :junction_charge, :turn_off_tail, :switching_energy,
                :thermal, :declared_model,
            )
        ])
    current = terminal_current_margins === nothing ? nothing : _semiconductor_public_record(
        "semiconductor.current-margin.v1", terminal_current_margins, (:current_a, :lower_a, :upper_a),
    )
    event_records = events === nothing ? nothing : [
        _semiconductor_public_record("semiconductor.energy-lookup.v1", row, _SEMICONDUCTOR_EVENT_RESULT_FIELDS)
        for row in events
    ]
    margin_records = event_margins === nothing ? nothing : map(_semiconductor_event_margin_record, event_margins)
    resolution_records = resolution === nothing ? nothing : [
        _semiconductor_public_record("semiconductor.timestep-resolution.v1", row, _SEMICONDUCTOR_RESOLUTION_RESULT_FIELDS)
        for row in resolution
    ]
    record = PortableSnapshotRecord(POWER_SEMICONDUCTOR_RESULT_SIGNATURE_SCHEMA, [
        "device_index" => device_index,
        "terminal_a" => terminals[1],
        "terminal_b" => terminals[2],
        "interval_start_s" => interval_start_s,
        "interval_end_s" => interval_end_s,
        "step_s" => step_s,
        "terminal_state" => terminal_state,
        "extended_result_version" => extended.schema_version,
        "extended_state" => endpoint,
        "device_provenance" => _semiconductor_provenance_record(extended.device_provenance),
        "component_provenance" => component_provenance,
        "terminal_current_margins" => current,
        "events" => event_records,
        "event_margins" => margin_records,
        "resolution" => resolution_records,
    ])
    io = IOBuffer()
    _write_portable_value(io, record)
    return take!(io)
end

function _power_semiconductor_public_result_signature(
    device_index, interval_start_s, interval_end_s, step_s, device, current, events, margins, resolution,
)
    bytes = power_semiconductor_public_result_bytes(
        device_index, (device.a, device.b), interval_start_s, interval_end_s, step_s,
        power_semiconductor_terminal_state(device), power_semiconductor_extended_result(device);
        terminal_current_margins=current, events, event_margins=margins, resolution,
    )
    return PowerSemiconductorResultSignature(
        POWER_SEMICONDUCTOR_RESULT_SIGNATURE_SCHEMA, device_index, interval_start_s, interval_end_s,
        bytes2hex(SHA.sha256(bytes)),
    )
end

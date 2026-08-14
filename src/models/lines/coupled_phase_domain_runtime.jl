using ..CoupledLineRuntime
using ..CoupledLineFitting: CoupledLineFitResult
using TOML

import ..Branches:
    trace_output_channel_count,
    trace_output_channel_names!,
    trace_output_is_public,
    trace_output_values!

export CoupledFrequencyDependentLine,
       coupled_frequency_dependent_line,
       initialize_coupled_frequency_dependent_line_deenergized!,
       initialize_coupled_frequency_dependent_line_sinusoidal!,
       CoupledFrequencyDependentLineSnapshot,
       coupled_frequency_dependent_line_snapshot,
       restore_coupled_frequency_dependent_line_snapshot!,
       write_coupled_frequency_dependent_line_snapshot,
       read_coupled_frequency_dependent_line_snapshot

mutable struct CoupledFrequencyDependentLine <: EMTElement
    from_nodes::Vector{Int}
    to_nodes::Vector{Int}
    port_nodes::Vector{Int}
    runtime_state::CoupledLineRuntimeState
    voltage_workspace_v::Vector{Float64}
end

function _checked_coupled_line_runtime_nodes(
    preparation::CoupledLineRuntimePreparation,
    from_nodes,
    to_nodes,
)
    phase_count = length(preparation.phase_order)
    length(from_nodes) == phase_count && length(to_nodes) == phase_count ||
        throw(ArgumentError(
            "coupled frequency-dependent line node counts must match its phases",
        ))
    sending = Int.(from_nodes)
    receiving = Int.(to_nodes)
    all(>=(0), sending) && all(>=(0), receiving) ||
        throw(ArgumentError(
            "coupled frequency-dependent line nodes must be nonnegative",
        ))
    all(index -> sending[index] != 0 || receiving[index] != 0, eachindex(sending)) ||
        throw(ArgumentError(
            "coupled frequency-dependent line cannot ground both ends of one phase",
        ))
    active_nodes = filter(!iszero, vcat(sending, receiving))
    length(unique(active_nodes)) == length(active_nodes) ||
        throw(ArgumentError(
            "coupled frequency-dependent line active terminal nodes must be unique",
        ))
    return sending, receiving, vcat(sending, receiving)
end

function coupled_frequency_dependent_line(
    preparation::CoupledLineRuntimePreparation,
    from_nodes,
    to_nodes,
)
    sending, receiving, ports = _checked_coupled_line_runtime_nodes(
        preparation,
        from_nodes,
        to_nodes,
    )
    state = coupled_line_runtime_state(preparation)
    return CoupledFrequencyDependentLine(
        sending,
        receiving,
        ports,
        state,
        zeros(length(ports)),
    )
end

function coupled_frequency_dependent_line(
    fit::CoupledLineFitResult,
    settings::CoupledLineRuntimeSettings,
    from_nodes,
    to_nodes,
)
    return coupled_frequency_dependent_line(
        prepare_coupled_line_runtime(fit, settings),
        from_nodes,
        to_nodes,
    )
end

function initialize_coupled_frequency_dependent_line_deenergized!(
    line::CoupledFrequencyDependentLine,
)
    initialize_coupled_line_runtime_deenergized!(line.runtime_state)
    fill!(line.voltage_workspace_v, 0.0)
    return line
end

function initialize_coupled_frequency_dependent_line_sinusoidal!(
    line::CoupledFrequencyDependentLine,
    from_voltage_phasor_v::AbstractVector{<:Complex},
    to_voltage_phasor_v::AbstractVector{<:Complex},
    frequency_hz::Real,
)
    phase_count = length(line.from_nodes)
    length(from_voltage_phasor_v) == phase_count &&
        length(to_voltage_phasor_v) == phase_count ||
        throw(ArgumentError(
            "coupled frequency-dependent line initialization phasors must match phases",
        ))
    initialize_coupled_line_runtime_sinusoidal!(
        line.runtime_state,
        vcat(from_voltage_phasor_v, to_voltage_phasor_v),
        frequency_hz,
    )
    line.voltage_workspace_v .= line.runtime_state.terminal_voltage_v
    return line
end

function _coupled_line_runtime_timestep_matches(
    line::CoupledFrequencyDependentLine,
    timestep_s::Float64,
)
    expected = line.runtime_state.preparation.settings.timestep_s
    abs(timestep_s - expected) <=
        64.0 * eps(Float64) * max(abs(timestep_s), abs(expected), 1.0) ||
        throw(ArgumentError(
            "coupled frequency-dependent line timestep changed after preparation",
        ))
    return nothing
end

function stamp!(
    admittance::AbstractMatrix{Float64},
    right_hand_side::AbstractVector{Float64},
    line::CoupledFrequencyDependentLine,
    _time_s::Float64,
    timestep_s::Float64,
)
    _coupled_line_runtime_timestep_matches(line, timestep_s)
    preparation = line.runtime_state.preparation
    for column in eachindex(line.port_nodes), row in eachindex(line.port_nodes)
        stamp_admittance_entry!(
            admittance,
            line.port_nodes[row],
            0,
            line.port_nodes[column],
            0,
            preparation.companion_admittance_s[row, column],
        )
    end
    for port in eachindex(line.port_nodes)
        stamp_history_current!(
            right_hand_side,
            line.port_nodes[port],
            0,
            line.runtime_state.history_current_a[port],
        )
    end
    return nothing
end

function update!(
    line::CoupledFrequencyDependentLine,
    voltage::AbstractVector{Float64},
    timestep_s::Float64,
)
    _coupled_line_runtime_timestep_matches(line, timestep_s)
    @inbounds for port in eachindex(line.port_nodes)
        node = line.port_nodes[port]
        line.voltage_workspace_v[port] = node == 0 ? 0.0 : voltage[node]
    end
    accept_coupled_line_runtime_step!(
        line.runtime_state,
        line.voltage_workspace_v,
    )
    return nothing
end

line_terminal_voltages(line::CoupledFrequencyDependentLine) = (
    from=copy(@view(line.runtime_state.terminal_voltage_v[1:length(line.from_nodes)])),
    to=copy(@view(line.runtime_state.terminal_voltage_v[(length(line.from_nodes) + 1):end])),
)

line_terminal_currents(line::CoupledFrequencyDependentLine) = (
    from=copy(@view(line.runtime_state.terminal_current_a[1:length(line.from_nodes)])),
    to=copy(@view(line.runtime_state.terminal_current_a[(length(line.from_nodes) + 1):end])),
)

line_history_currents(line::CoupledFrequencyDependentLine) = (
    from=copy(@view(line.runtime_state.history_current_a[1:length(line.from_nodes)])),
    to=copy(@view(line.runtime_state.history_current_a[(length(line.from_nodes) + 1):end])),
)

struct CoupledFrequencyDependentLineSnapshot
    schema_version::Int
    from_nodes::Vector{Int}
    to_nodes::Vector{Int}
    runtime_snapshot::CoupledLineRuntimeSnapshot
end

function coupled_frequency_dependent_line_snapshot(
    line::CoupledFrequencyDependentLine,
)
    return CoupledFrequencyDependentLineSnapshot(
        1,
        copy(line.from_nodes),
        copy(line.to_nodes),
        coupled_line_runtime_snapshot(line.runtime_state),
    )
end

function restore_coupled_frequency_dependent_line_snapshot!(
    line::CoupledFrequencyDependentLine,
    snapshot::CoupledFrequencyDependentLineSnapshot,
)
    snapshot.schema_version == 1 || throw(ArgumentError(
        "coupled frequency-dependent line snapshot schema is unsupported",
    ))
    snapshot.from_nodes == line.from_nodes && snapshot.to_nodes == line.to_nodes ||
        throw(ArgumentError(
            "coupled frequency-dependent line snapshot terminal nodes are stale",
        ))
    restore_coupled_line_runtime_snapshot!(line.runtime_state, snapshot.runtime_snapshot)
    line.voltage_workspace_v .= line.runtime_state.terminal_voltage_v
    return line
end

function write_coupled_frequency_dependent_line_snapshot(
    path::AbstractString,
    snapshot::CoupledFrequencyDependentLineSnapshot,
)
    output_path = abspath(path)
    data = Dict{String,Any}(
        "schema" => "aimora.coupled_frequency_dependent_line_snapshot.v1",
        "schema_version" => snapshot.schema_version,
        "from_nodes" => snapshot.from_nodes,
        "to_nodes" => snapshot.to_nodes,
        "runtime" => CoupledLineRuntime._coupled_line_runtime_snapshot_dictionary(
            snapshot.runtime_snapshot,
        ),
    )
    mkpath(dirname(output_path))
    mktemp(dirname(output_path)) do temporary_path, io
        TOML.print(io, data; sorted=true)
        close(io)
        mv(temporary_path, output_path; force=true)
    end
    return output_path
end

function read_coupled_frequency_dependent_line_snapshot(path::AbstractString)
    isfile(path) || throw(ArgumentError(
        "coupled frequency-dependent line snapshot file does not exist",
    ))
    data = try
        TOML.parsefile(path)
    catch
        throw(ArgumentError(
            "coupled frequency-dependent line snapshot is not valid TOML",
        ))
    end
    get(data, "schema", "") ==
        "aimora.coupled_frequency_dependent_line_snapshot.v1" ||
        throw(ArgumentError(
            "coupled frequency-dependent line snapshot schema is unsupported",
        ))
    try
        return CoupledFrequencyDependentLineSnapshot(
            Int(data["schema_version"]),
            Int.(data["from_nodes"]),
            Int.(data["to_nodes"]),
            CoupledLineRuntime._coupled_line_runtime_snapshot_from_dictionary(
                data["runtime"],
            ),
        )
    catch error
        error isa ArgumentError && rethrow()
        throw(ArgumentError(
            "coupled frequency-dependent line snapshot is malformed",
        ))
    end
end

trace_output_channel_count(line::CoupledFrequencyDependentLine) =
    8 * length(line.from_nodes) + 3

trace_output_is_public(::CoupledFrequencyDependentLine) = true

function trace_output_channel_names!(
    names::Vector{Symbol},
    element_name::Symbol,
    line::CoupledFrequencyDependentLine,
)
    phase_count = length(line.from_nodes)
    for phase in 1:phase_count
        append!(names, (
            Symbol(element_name, :_from_current_, phase, :_a),
            Symbol(element_name, :_to_current_, phase, :_a),
            Symbol(element_name, :_from_history_, phase, :_a),
            Symbol(element_name, :_to_history_, phase, :_a),
            Symbol(element_name, :_from_incident_wave_, phase),
            Symbol(element_name, :_to_incident_wave_, phase),
            Symbol(element_name, :_from_outgoing_wave_, phase),
            Symbol(element_name, :_to_outgoing_wave_, phase),
        ))
    end
    append!(names, (
        Symbol(element_name, :_terminal_power_w),
        Symbol(element_name, :_cumulative_supplied_energy_j),
        Symbol(element_name, :_maximum_kcl_residual_a),
    ))
    return names
end

function trace_output_values!(
    output::AbstractMatrix{Float64},
    first_channel::Int,
    sample::Int,
    line::CoupledFrequencyDependentLine,
    _voltage::AbstractVector{Float64},
)
    state = line.runtime_state
    phase_count = length(line.from_nodes)
    channel = first_channel
    for phase in 1:phase_count
        receiving = phase_count + phase
        output[channel, sample] = state.terminal_current_a[phase]
        output[channel + 1, sample] = state.terminal_current_a[receiving]
        output[channel + 2, sample] = state.history_current_a[phase]
        output[channel + 3, sample] = state.history_current_a[receiving]
        output[channel + 4, sample] = state.incident_wave[phase]
        output[channel + 5, sample] = state.incident_wave[receiving]
        output[channel + 6, sample] = state.outgoing_wave[phase]
        output[channel + 7, sample] = state.outgoing_wave[receiving]
        channel += 8
    end
    output[channel, sample] = state.terminal_power_w
    output[channel + 1, sample] = state.cumulative_supplied_energy_j
    output[channel + 2, sample] = state.maximum_kcl_residual_a
    return channel + 3
end

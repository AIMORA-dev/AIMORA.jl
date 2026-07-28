using SHA
using ..DeckParser
using ..EMTStudy:
    DeckEMTTrace,
    DeckEMTExecution,
    EMTTerminalState,
    EMTStudyWorkspace,
    DeckOutputChannelMetadata,
    deck_report_output_trace,
    deck_output_channel_metadata,
    deck_output_step_indices,
    evaluate_emt_study!,
    prepare_emt_study,
    read_emt_checkpoint,
    restart_emt_checkpoint,
    run_deck_emt,
    run_deck_emt_execution,
    write_emt_checkpoint

export DeckEMTOutputBundle,
       DeckBinaryPlotChannelMetadata,
       DeckEMTResumableOutputArchive,
       deck_binary_plot_channel_metadata,
       deck_emt_output_bundle,
       run_deck_emt_output_bundle,
       run_deck_emt_resumable_output_archive

struct DeckBinaryPlotChannelMetadata
    name::Symbol
    physical_quantity::Symbol
    physical_unit::String
    storage_class::Symbol
    upper_name::String
    lower_name::String
end

struct DeckEMTOutputBundle
    source::String
    trace::DeckEMTTrace
    channel_metadata::Vector{DeckOutputChannelMetadata}
    binary_channel_metadata::Vector{DeckBinaryPlotChannelMetadata}
    report_table::ElectromagneticReportOutputTable
    binary_plot_table::ElectromagneticBinaryPlotTable
    column_width::Int
    print_step_indices::Vector{Int}
    plot_step_indices::Vector{Int}
    plot_file_retention_mode::Int
    restart_snapshot_requested::Bool
    physical_checks_passed::Bool
    deferred_effects::Tuple{Vararg{Symbol}}
end

struct DeckEMTResumableOutputArchive{I,C,R}
    source::String
    initial_output::I
    continued_output::C
    restart_run::R
    initial_checkpoint_path::String
    continued_checkpoint_path::String
    physical_checks_passed::Bool
end

const _SYNCHRONOUS_MACHINE_BINARY_LABELS = Dict(
    1 => "ID",
    2 => "IQ",
    3 => "I0",
    4 => "IF",
    5 => "IKD",
    6 => "IG",
    7 => "IKQ",
    8 => "IA",
    9 => "IB",
    10 => "IC",
    11 => "EFD",
    12 => "MFORCE",
    13 => "MANG",
    14 => "TQ GEN",
    15 => "TQ EXC",
)

function _binary_plot_storage_name(name::AbstractString)
    text = String(strip(String(name)))
    isempty(text) && return ""
    if all(character -> Int(character) < 0x80, text) && ncodeunits(text) <= 8
        return text
    end
    digest = bytes2hex(SHA.sha1(codeunits(text)))
    prefix = String(first(text, min(length(text), 3)))
    ascii_prefix = String(character for character in prefix if Int(character) < 0x80)
    return string(rpad(ascii_prefix, 3, '_')[1:3], uppercase(digest[1:5]))
end

function _synchronous_machine_binary_labels(
    parsed::DeckParser.DeckParseResult,
)
    labels = Tuple{Int,String}[]
    for row in DeckParser.deck_synchronous_machine_output_request_rows(parsed)
        if row.group_index == 1
            for code in row.output_codes
                label = get(_SYNCHRONOUS_MACHINE_BINARY_LABELS, code, nothing)
                label === nothing && continue
                push!(labels, (row.machine_index, label))
            end
        elseif row.group_index in (2, 3, 4)
            prefix = row.group_index == 2 ? "ANG" :
                row.group_index == 3 ? "VEL" : "TOR"
            for index in row.output_codes
                push!(labels, (row.machine_index, string(prefix, " ", index)))
            end
        end
    end
    return labels
end

function deck_binary_plot_channel_metadata(
    parsed::DeckParser.DeckParseResult,
    metadata::AbstractVector{DeckOutputChannelMetadata},
)
    machine_labels = _synchronous_machine_binary_labels(parsed)
    machine_label_index = 0
    storage = DeckBinaryPlotChannelMetadata[]
    sizehint!(storage, length(metadata))
    for row in metadata
        if row.quantity in (:node_voltage, :branch_voltage, :branch_current)
            storage_class = row.quantity
            upper_name = _binary_plot_storage_name(row.upper_name)
            lower_name = _binary_plot_storage_name(row.lower_name)
        elseif row.quantity == :tacs
            storage_class = :branch_current
            upper_name = "TACS"
            lower_name = _binary_plot_storage_name(row.lower_name)
        elseif row.quantity == :synchronous_machine
            machine_label_index += 1
            machine_label_index <= length(machine_labels) || throw(ArgumentError(
                "synchronous-machine trace has more channels than its deck requests",
            ))
            machine_index, label = machine_labels[machine_label_index]
            storage_class = :branch_current
            upper_name = _binary_plot_storage_name("MACH $machine_index")
            lower_name = _binary_plot_storage_name(label)
        else
            throw(ArgumentError(
                "output quantity $(row.quantity) has no binary plot storage mapping",
            ))
        end
        push!(
            storage,
            DeckBinaryPlotChannelMetadata(
                row.name,
                row.quantity,
                row.unit,
                storage_class,
                upper_name,
                lower_name,
            ),
        )
    end
    machine_label_index == length(machine_labels) || throw(ArgumentError(
        "synchronous-machine deck requests do not match the emitted trace channels",
    ))
    ranks = Dict(:node_voltage => 1, :branch_voltage => 2, :branch_current => 3)
    issorted([ranks[row.storage_class] for row in storage]) || throw(ArgumentError(
        "binary plot channels must retain node-voltage, branch-voltage, branch-current order",
    ))
    return storage
end

function _deck_output_column_width(parsed::DeckParser.DeckParseResult)
    requests = DeckParser.deck_output_width_request_rows(parsed)
    width = isempty(requests) ? 132 : last(requests).column_width
    width in (80, 132) ||
        throw(ArgumentError("deck output width must resolve to 80 or 132 columns"))
    return width
end

function _deck_plot_file_retention_mode(
    parsed::DeckParser.DeckParseResult,
)
    options = DeckParser.deck_output_schedule_options(parsed)
    requests = DeckParser.deck_plot_file_request_rows(parsed)
    mode = isempty(requests) ?
        options.plot_file_retention_mode :
        last(requests).plot_file_mode
    mode in (0, 1, 2) || throw(ArgumentError(
        "deck plot-file retention mode must resolve to 0, 1, or 2",
    ))
    return mode
end

function _trace_output_columns(
    trace::DeckEMTTrace,
    step_indices::AbstractVector{<:Integer},
)
    columns = Dict{Int,Int}()
    for column in eachindex(trace.time_s)
        step = Int(round(trace.time_s[column] / trace.dt_s))
        haskey(columns, step) &&
            throw(ArgumentError("trace contains duplicate recorded output step $step"))
        columns[step] = column
    end
    selected = Int[]
    for step in step_indices
        column = get(columns, Int(step), 0)
        column > 0 || throw(ArgumentError("trace does not contain requested output step $step"))
        push!(selected, column)
    end
    return selected
end

function _deck_output_extrema(trace::DeckEMTTrace)
    channel_count = length(trace.output_channel_names)
    kinds = Symbol[]
    indices = Int[]
    values = Float64[]
    for (kind, source) in (
        (:maxima, trace.output_maximum_values),
        (:times_of_maxima, trace.output_maximum_times_s),
        (:minima, trace.output_minimum_values),
        (:times_of_minima, trace.output_minimum_times_s),
    )
        length(source) == channel_count ||
            throw(ArgumentError("trace extrema must cover every output channel"))
        append!(kinds, fill(kind, channel_count))
        append!(indices, 1:channel_count)
        append!(values, Float64.(source))
    end
    return kinds, indices, values
end

function _deck_output_label(row::DeckOutputChannelMetadata)
    isempty(row.lower_name) && return row.upper_name
    return string(row.upper_name, ".", row.lower_name)
end

function _deck_emt_report_components(
    parsed::DeckParser.DeckParseResult,
    trace::DeckEMTTrace;
    title::AbstractString = "AIMORA electromagnetic output",
)
    DeckParser.assert_deck_valid!(parsed)
    trace.source == parsed.source ||
        throw(ArgumentError("trace source must match the parsed deck source"))
    metadata = deck_output_channel_metadata(parsed, trace)
    print_steps = deck_output_step_indices(
        parsed,
        trace.dt_s,
        trace.t_end_s;
        schedule = :print,
    )
    plot_steps = deck_output_step_indices(
        parsed,
        trace.dt_s,
        trace.t_end_s;
        schedule = :plot,
    )
    print_columns = _trace_output_columns(trace, print_steps)
    plot_columns = _trace_output_columns(trace, plot_steps)
    extrema_kinds, extrema_indices, extrema_values = _deck_output_extrema(trace)
    report = electromagnetic_report_output_table(
        title = String(title),
        source = parsed.source,
        channel_classes = [row.quantity for row in metadata],
        upper_names = [row.upper_name for row in metadata],
        lower_names = [row.lower_name for row in metadata],
        labels = [_deck_output_label(row) for row in metadata],
        sample_steps = print_steps,
        sample_times_s = trace.time_s[print_columns],
        sample_values = permutedims(trace.output_pu[:, print_columns]),
        extrema_kinds = extrema_kinds,
        extrema_indices = extrema_indices,
        extrema_values = extrema_values,
        full_text_report_compatibility = true,
        deferred_effects = (),
    )
    checks = issorted(print_steps) && issorted(plot_steps) &&
        all(isfinite, report.sample_values) &&
        all(isfinite, report.extrema_values) &&
        size(report.sample_values, 2) == length(metadata)
    return (
        metadata = metadata,
        report = report,
        print_steps = print_steps,
        plot_steps = plot_steps,
        plot_columns = plot_columns,
        column_width = _deck_output_column_width(parsed),
        physical_checks_passed = checks,
    )
end

function deck_emt_output_bundle(
    parsed::DeckParser.DeckParseResult,
    trace::DeckEMTTrace;
    title::AbstractString = "AIMORA electromagnetic output",
)
    report_trace = deck_report_output_trace(parsed, trace)
    components = _deck_emt_report_components(parsed, report_trace; title)
    metadata = components.metadata
    binary_metadata = deck_binary_plot_channel_metadata(parsed, metadata)
    plot = electromagnetic_binary_plot_table(
        title = String(title),
        source = parsed.source,
        channel_classes = [row.storage_class for row in binary_metadata],
        upper_names = [row.upper_name for row in binary_metadata],
        lower_names = [row.lower_name for row in binary_metadata],
        sample_times_s = report_trace.time_s[components.plot_columns],
        sample_values =
            permutedims(report_trace.output_pu[:, components.plot_columns]),
        binary_plot_format_compatibility = true,
        deferred_effects = (),
    )
    schedule_options = DeckParser.deck_output_schedule_options(parsed)
    checks = components.physical_checks_passed &&
        all(isfinite, plot.sample_values) &&
        size(plot.sample_values, 2) == length(metadata) &&
        length(binary_metadata) == length(metadata) &&
        plot.channel_classes == [row.storage_class for row in binary_metadata]
    return DeckEMTOutputBundle(
        parsed.source,
        report_trace,
        metadata,
        binary_metadata,
        components.report,
        plot,
        components.column_width,
        components.print_steps,
        components.plot_steps,
        _deck_plot_file_retention_mode(parsed),
        schedule_options.restart_snapshot_enabled,
        checks,
        (),
    )
end

function _write_deck_emt_output_bundle(
    parsed::DeckParser.DeckParseResult,
    trace::DeckEMTTrace,
    output_dir::AbstractString;
    basename::AbstractString = "deck_output",
    title::AbstractString = "AIMORA electromagnetic output",
    replay_commands = nothing,
)
    bundle = deck_emt_output_bundle(parsed, trace; title = title)
    bundle.physical_checks_passed ||
        throw(ErrorException("deck electromagnetic output bundle failed physical checks"))
    mkpath(output_dir)
    report_path = joinpath(output_dir, string(basename, ".out"))
    plot_path = joinpath(output_dir, string(basename, ".PL4"))
    write_electromagnetic_report_output_fixed_width_text(
        report_path,
        bundle.report_table;
        column_width = bundle.column_width,
    )
    write_electromagnetic_binary_plot_file(plot_path, bundle.binary_plot_table)
    roundtrip = read_electromagnetic_binary_plot_file(plot_path)
    roundtrip.sample_times_s == Float32.(bundle.binary_plot_table.sample_times_s) ||
        throw(ErrorException("generated binary plot time roundtrip mismatch"))
    roundtrip.sample_values == Float32.(bundle.binary_plot_table.sample_values) ||
        throw(ErrorException("generated binary plot value roundtrip mismatch"))
    commands = replay_commands === nothing ? (
        "SPY",
        "LUNIT4",
        "OPEN",
        plot_path,
        "OLD",
        "TOP",
        "NEXT",
        "BACK",
        "TIME",
        "0.0",
        "BOT",
        "CLOSE",
        "N",
    ) : replay_commands
    replay = electromagnetic_plot_replay_session(plot_path; commands)
    return (
        bundle = bundle,
        report_path = report_path,
        plot_path = plot_path,
        binary_plot_roundtrip = roundtrip,
        replay = replay,
    )
end

function run_deck_emt_output_bundle(
    parsed::DeckParser.DeckParseResult,
    output_dir::AbstractString;
    basename::AbstractString = "deck_output",
    title::AbstractString = "AIMORA electromagnetic output",
    replay_commands = nothing,
    runtime_kwargs...,
)
    trace = run_deck_emt(
        parsed;
        time_horizon = :deck,
        output_schedule = :print_and_plot,
        runtime_kwargs...,
    )
    return _write_deck_emt_output_bundle(
        parsed,
        trace,
        output_dir;
        basename,
        title,
        replay_commands,
    )
end

function run_deck_emt_resumable_output_archive(
    parsed::DeckParser.DeckParseResult,
    request::DeckParser.DeckRestartRequest,
    output_dir::AbstractString;
    additional_time_s::Real,
    basename::AbstractString = "deck_output",
    title::AbstractString = "AIMORA electromagnetic output",
    replay_commands = nothing,
    runtime_kwargs...,
)
    prepared = prepare_emt_study(
        parsed;
        time_horizon = :deck,
        output_schedule = :print_and_plot,
        runtime_kwargs...,
    )
    workspace = EMTStudyWorkspace(prepared)
    initial_trace = evaluate_emt_study!(workspace)
    initial_output = _write_deck_emt_output_bundle(
        parsed,
        initial_trace,
        output_dir;
        basename = string(basename, "_initial"),
        title,
        replay_commands,
    )

    initial_checkpoint_path = joinpath(
        output_dir,
        string(basename, "_initial.aimora-checkpoint"),
    )
    write_emt_checkpoint(initial_checkpoint_path, workspace)
    final_time_s = initial_trace.t_end_s + Float64(additional_time_s)
    checkpoint_step = Int(round(initial_trace.t_end_s / initial_trace.dt_s))
    future_steps = filter(
        step -> step > checkpoint_step,
        deck_output_step_indices(
            parsed,
            initial_trace.dt_s,
            final_time_s;
            schedule = :print_and_plot,
        ),
    )
    isempty(future_steps) && throw(ArgumentError(
        "resumable output archive has no scheduled continuation samples",
    ))
    continued_checkpoint_path = joinpath(
        output_dir,
        string(basename, "_continued.aimora-checkpoint"),
    )
    restart_report_path = joinpath(
        output_dir,
        string(basename, "_restart.json"),
    )
    continuation = restart_emt_checkpoint(
        initial_checkpoint_path,
        request;
        additional_time_s,
        recorded_step_indices = future_steps,
        report_path = restart_report_path,
        continued_checkpoint_path,
    )
    continued_output = _write_deck_emt_output_bundle(
        parsed,
        continuation.run.trace,
        output_dir;
        basename = string(basename, "_continued"),
        title,
        replay_commands,
    )
    read_emt_checkpoint(continued_checkpoint_path)
    physical_checks = initial_output.bundle.physical_checks_passed &&
        continued_output.bundle.physical_checks_passed &&
        continuation.run.checkpoint_state_error <= 1.0e-9 &&
        isfinite(continuation.run.final_kcl_error) &&
        continuation.run.final_time_s > continuation.run.checkpoint_time_s &&
        first(continuation.run.trace.time_s) == first(initial_trace.time_s) &&
        last(continuation.run.trace.time_s) == continuation.run.final_time_s
    physical_checks || throw(ErrorException(
        "resumable electromagnetic output archive failed physical checks",
    ))
    return DeckEMTResumableOutputArchive(
        parsed.source,
        initial_output,
        continued_output,
        continuation.run,
        initial_checkpoint_path,
        continued_checkpoint_path,
        physical_checks,
    )
end

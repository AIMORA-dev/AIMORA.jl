export ElectromagneticTerminalPunchRow,
       ElectromagneticChannelStatistics,
       ElectromagneticExpressionReportRow,
       ElectromagneticTextOutputSuite,
       electromagnetic_text_output_suite,
       write_electromagnetic_terminal_punch,
       write_electromagnetic_statistics_text,
       write_electromagnetic_spy_transcript,
       write_electromagnetic_expression_report,
       run_deck_emt_text_output_suite

struct ElectromagneticTerminalPunchRow
    component::Symbol
    name::Symbol
    from_node::Symbol
    to_node::Symbol
    state::Symbol
    unit::String
    value::Float64
end

struct ElectromagneticChannelStatistics
    name::Symbol
    quantity::Symbol
    unit::String
    sample_count::Int
    mean_value::Float64
    root_mean_square_value::Float64
    variance_value::Float64
    standard_deviation_value::Float64
    maximum_value::Float64
    maximum_time_s::Float64
    minimum_value::Float64
    minimum_time_s::Float64
end

struct ElectromagneticExpressionReportRow
    name::Symbol
    expression::String
    unit::String
    sample_count::Int
    final_value::Float64
    maximum_value::Float64
    maximum_time_s::Float64
    minimum_value::Float64
    minimum_time_s::Float64
end

struct ElectromagneticTextOutputSuite
    source::String
    report_table::ElectromagneticReportOutputTable
    channel_metadata::Vector{DeckOutputChannelMetadata}
    terminal_state::EMTTerminalState
    terminal_punch_rows::Vector{ElectromagneticTerminalPunchRow}
    statistics_rows::Vector{ElectromagneticChannelStatistics}
    expression_rows::Vector{ElectromagneticExpressionReportRow}
    report_path::String
    terminal_punch_path::String
    statistics_path::String
    emtspy_path::String
    expression_report_path::String
    terminal_punch_requested::Bool
    physical_checks_passed::Bool
end

function _text_output_extrema(
    values::AbstractVector{<:Real},
    times_s::AbstractVector{<:Real},
)
    length(values) == length(times_s) ||
        throw(ArgumentError("text-output values and times must have equal lengths"))
    isempty(values) && throw(ArgumentError("text-output rows require at least one sample"))
    maximum_index = argmax(values)
    minimum_index = argmin(values)
    return (
        maximum_value = Float64(values[maximum_index]),
        maximum_time_s = Float64(times_s[maximum_index]),
        minimum_value = Float64(values[minimum_index]),
        minimum_time_s = Float64(times_s[minimum_index]),
    )
end

function _text_output_statistics_row(
    metadata::DeckOutputChannelMetadata,
    values::AbstractVector{<:Real},
    times_s::AbstractVector{<:Real},
)
    extrema = _text_output_extrema(values, times_s)
    float_values = Float64.(values)
    mean_value = sum(float_values) / length(float_values)
    variance_value = length(float_values) == 1 ? 0.0 :
        sum(value -> abs2(value - mean_value), float_values) /
        (length(float_values) - 1)
    return ElectromagneticChannelStatistics(
        metadata.name,
        metadata.quantity,
        metadata.unit,
        length(float_values),
        mean_value,
        sqrt(sum(abs2, float_values) / length(float_values)),
        variance_value,
        sqrt(variance_value),
        extrema.maximum_value,
        extrema.maximum_time_s,
        extrema.minimum_value,
        extrema.minimum_time_s,
    )
end

function _terminal_punch_row(
    component::Symbol,
    name::Symbol,
    from_node::Symbol,
    to_node::Symbol,
    state::Symbol,
    unit::AbstractString,
    value::Real,
)
    return ElectromagneticTerminalPunchRow(
        component,
        name,
        from_node,
        to_node,
        state,
        String(unit),
        Float64(value),
    )
end

function _terminal_punch_rows(state::EMTTerminalState)
    rows = ElectromagneticTerminalPunchRow[]
    for node in state.nodes
        push!(rows, _terminal_punch_row(
            :node, node.name, node.name, :ground, :voltage, "V", node.voltage_v,
        ))
    end
    branch_fields = (
        (:conductance, "S", row -> row.conductance_s),
        (:history_current, "A", row -> row.history_current_a),
        (:voltage, "V", row -> row.voltage_v),
        (:current, "A", row -> row.current_a),
        (:power, "W", row -> row.voltage_v * row.current_a),
        (:previous_current, "A", row -> row.previous_current_a),
        (:previous_voltage, "V", row -> row.previous_voltage_v),
    )
    for branch in state.branches, (field, unit, value) in branch_fields
        push!(rows, _terminal_punch_row(
            :branch,
            branch.name,
            branch.from_node,
            branch.to_node,
            field,
            unit,
            value(branch),
        ))
    end
    for branch in state.branches
        isfinite(branch.energy_j) || continue
        push!(rows, _terminal_punch_row(
            :branch,
            branch.name,
            branch.from_node,
            branch.to_node,
            :energy,
            "J",
            branch.energy_j,
        ))
    end
    nonlinear_fields = (
        (:voltage, "V", row -> row.voltage_v),
        (:current, "A", row -> row.current_a),
        (:power, "W", row -> row.voltage_v * row.current_a),
        (:active_segment, "count", row -> row.active_segment),
        (:energy, "J", row -> row.energy_j),
    )
    for nonlinear in state.nonlinear_elements
        for (field, unit, value) in nonlinear_fields
            push!(rows, _terminal_punch_row(
                :nonlinear,
                nonlinear.name,
                nonlinear.from_node,
                nonlinear.to_node,
                field,
                unit,
                value(nonlinear),
            ))
        end
        isfinite(nonlinear.flux_wb) && push!(rows, _terminal_punch_row(
            :nonlinear,
            nonlinear.name,
            nonlinear.from_node,
            nonlinear.to_node,
            :flux,
            "Wb",
            nonlinear.flux_wb,
        ))
    end
    switch_fields = (
        (:closed, "bool", row -> row.closed ? 1.0 : 0.0),
        (:conductance, "S", row -> row.conductance_s),
        (:voltage, "V", row -> row.voltage_v),
        (:current, "A", row -> row.current_a),
        (:power, "W", row -> row.power_w),
        (:close_time, "s", row -> row.close_time_s),
        (:open_time, "s", row -> row.open_time_s),
    )
    for switch in state.switches, (field, unit, value) in switch_fields
        push!(rows, _terminal_punch_row(
            :switch,
            switch.name,
            switch.from_node,
            switch.to_node,
            field,
            unit,
            value(switch),
        ))
    end
    for switch in state.switches
        isfinite(switch.energy_j) || continue
        push!(rows, _terminal_punch_row(
            :switch,
            switch.name,
            switch.from_node,
            switch.to_node,
            :energy,
            "J",
            switch.energy_j,
        ))
    end
    return rows
end

function _expression_report_rows(
    parsed::DeckParser.DeckParseResult,
    trace::DeckEMTTrace,
    metadata::Vector{DeckOutputChannelMetadata},
)
    channel_index = Dict(name => index for (index, name) in enumerate(trace.output_channel_names))
    metadata_by_name = Dict(row.name => row for row in metadata)
    rows = ElectromagneticExpressionReportRow[]
    for expression in DeckParser.deck_control_system_expression_rows(parsed)
        index = get(channel_index, expression.name, 0)
        index > 0 || continue
        values = view(trace.output_pu, index, :)
        extrema = _text_output_extrema(values, trace.time_s)
        channel = get(metadata_by_name, expression.name, nothing)
        push!(
            rows,
            ElectromagneticExpressionReportRow(
                expression.name,
                expression.expression,
                channel === nothing ? "1" : channel.unit,
                length(values),
                Float64(last(values)),
                extrema.maximum_value,
                extrema.maximum_time_s,
                extrema.minimum_value,
                extrema.minimum_time_s,
            ),
        )
    end
    return rows
end

function _text_output_punch_requested(parsed::DeckParser.DeckParseResult)
    options = DeckParser.deck_output_schedule_options(parsed)
    return options.terminal_conditions_punch_enabled
end

function _text_output_physical_checks(
    trace::DeckEMTTrace,
    report_table::ElectromagneticReportOutputTable,
    metadata,
    terminal_state,
    punch_rows,
    statistics_rows,
)
    length(metadata) == size(trace.output_pu, 1) || return false
    all(row -> !isempty(row.unit), metadata) || return false
    all(isfinite, trace.output_pu) || return false
    all(
        row -> row.sample_count == length(report_table.sample_times_s),
        statistics_rows,
    ) || return false
    all(row -> isfinite(row.value) || row.value == Inf, punch_rows) || return false
    terminal_state.time_s == last(trace.time_s) || return false
    for (index, row) in enumerate(statistics_rows)
        values = view(report_table.sample_values, :, index)
        row.maximum_value == maximum(values) || return false
        row.minimum_value == minimum(values) || return false
        isapprox(
            row.root_mean_square_value,
            sqrt(sum(abs2, values) / length(values));
            atol = 0.0,
            rtol = 4eps(Float64),
        ) || return false
        mean_value = sum(values) / length(values)
        variance_value = length(values) == 1 ? 0.0 :
            sum(value -> abs2(value - mean_value), values) /
            (length(values) - 1)
        isapprox(row.mean_value, mean_value; atol = 0.0, rtol = 4eps(Float64)) ||
            return false
        isapprox(
            row.variance_value,
            variance_value;
            atol = 0.0,
            rtol = 8eps(Float64),
        ) || return false
        isapprox(
            row.standard_deviation_value,
            sqrt(variance_value);
            atol = 0.0,
            rtol = 4eps(Float64),
        ) || return false
    end
    return true
end

function electromagnetic_text_output_suite(
    parsed::DeckParser.DeckParseResult,
    execution::DeckEMTExecution;
    title::AbstractString = "AIMORA electromagnetic output",
)
    DeckParser.assert_deck_valid!(parsed)
    trace = execution.trace
    components = _deck_emt_report_components(parsed, trace; title)
    metadata = components.metadata
    statistics_rows = ElectromagneticChannelStatistics[
        _text_output_statistics_row(
            metadata[index],
            view(components.report.sample_values, :, index),
            components.report.sample_times_s,
        )
        for index in eachindex(metadata)
    ]
    punch_rows = _terminal_punch_rows(execution.terminal_state)
    expression_rows = _expression_report_rows(parsed, trace, metadata)
    punch_requested = _text_output_punch_requested(parsed)
    physical_checks =
        components.physical_checks_passed &&
        _text_output_physical_checks(
            trace,
            components.report,
            metadata,
            execution.terminal_state,
            punch_rows,
            statistics_rows,
        ) &&
        (!punch_requested || execution.terminal_state.physical_checks_passed)
    return ElectromagneticTextOutputSuite(
        parsed.source,
        components.report,
        copy(metadata),
        execution.terminal_state,
        punch_rows,
        statistics_rows,
        expression_rows,
        "",
        "",
        "",
        "",
        "",
        punch_requested,
        physical_checks,
    )
end

function write_electromagnetic_terminal_punch(
    path::AbstractString,
    rows::AbstractVector{ElectromagneticTerminalPunchRow},
)
    _report_ensure_dir(path)
    open(path, "w") do io
        println(io, "AIMORA ELECTROMAGNETIC TERMINAL CONDITIONS")
        println(io, "COMPONENT    NAME                 FROM       TO         STATE                  UNIT              VALUE")
        for row in rows
            @printf(
                io,
                "%-12s %-20s %-10s %-10s %-22s %-8s %16.8E\n",
                uppercase(String(row.component)),
                uppercase(String(row.name)),
                uppercase(String(row.from_node)),
                uppercase(String(row.to_node)),
                uppercase(String(row.state)),
                row.unit,
                row.value,
            )
        end
    end
    return abspath(String(path))
end

function write_electromagnetic_statistics_text(
    path::AbstractString,
    rows::AbstractVector{ElectromagneticChannelStatistics},
)
    _report_ensure_dir(path)
    open(path, "w") do io
        println(io, "AIMORA ELECTROMAGNETIC CHANNEL SUMMARY STATISTICS")
        println(io, "NAME             QUANTITY       UNIT      N         MEAN          RMS          MAX       T(MAX)          MIN       T(MIN)")
        for row in rows
            @printf(
                io,
                "%-16s %-14s %-5s %5d %12.4E %12.4E %12.4E %11.4E %12.4E %11.4E\n",
                uppercase(String(row.name)),
                uppercase(String(row.quantity)),
                row.unit,
                row.sample_count,
                row.mean_value,
                row.root_mean_square_value,
                row.maximum_value,
                row.maximum_time_s,
                row.minimum_value,
                row.minimum_time_s,
            )
            @printf(
                io,
                "  DISTRIBUTION PARAMETERS: VARIANCE=%16.8E STANDARD_DEVIATION=%16.8E\n",
                row.variance_value,
                row.standard_deviation_value,
            )
        end
    end
    return abspath(String(path))
end

function write_electromagnetic_spy_transcript(
    path::AbstractString,
    suite::ElectromagneticTextOutputSuite,
)
    _report_ensure_dir(path)
    final_values = Dict(
        suite.channel_metadata[index].name =>
            suite.report_table.sample_values[end, index]
        for index in eachindex(suite.channel_metadata)
    )
    open(path, "w") do io
        println(io, "EMTP BEGINS.  SEND (SPY, HELP, MODULE, STOP) :")
        println(io, "SPY:")
        println(io, "AIMORA TEXT OUTPUT STATUS")
        println(io, "SOURCE: ", suite.source)
        println(io, "CHANNELS: ", length(suite.channel_metadata))
        println(io, "PRINTED SAMPLES: ", length(suite.report_table.sample_steps))
        println(io, "TERMINAL CONDITIONS: ", suite.terminal_punch_requested ? "SAVED" : "NOT REQUESTED")
        for row in suite.statistics_rows
            @printf(
                io,
                "CHANNEL %s %s [%s] FINAL=%12.4E\n",
                uppercase(String(row.name)),
                uppercase(String(row.quantity)),
                row.unit,
                final_values[row.name],
            )
            @printf(
                io,
                "  MAX=%12.4E AT=%10.3E S MIN=%12.4E AT=%10.3E S\n",
                row.maximum_value,
                row.maximum_time_s,
                row.minimum_value,
                row.minimum_time_s,
            )
        end
        println(io, "SPY RETURNS TO EMTP.")
    end
    return abspath(String(path))
end

function _write_expression_source(io::IO, expression::AbstractString)
    characters = collect(String(expression))
    isempty(characters) && return println(io, "  EXPRESSION:")
    for first_index in 1:112:length(characters)
        last_index = min(first_index + 111, length(characters))
        prefix = first_index == 1 ? "  EXPRESSION: " : "              "
        println(io, prefix, String(characters[first_index:last_index]))
    end
    return nothing
end

function write_electromagnetic_expression_report(
    path::AbstractString,
    rows::AbstractVector{ElectromagneticExpressionReportRow},
)
    _report_ensure_dir(path)
    open(path, "w") do io
        println(io, "AIMORA XPR CONTROL-EXPRESSION REPORT")
        println(io, "NAME             UNIT       N        FINAL      MAXIMUM       T(MAX)      MINIMUM       T(MIN)")
        for row in rows
            @printf(
                io,
                "%-16s %-6s %5d %12.4E %12.4E %12.4E %12.4E %12.4E\n",
                uppercase(String(row.name)),
                row.unit,
                row.sample_count,
                row.final_value,
                row.maximum_value,
                row.maximum_time_s,
                row.minimum_value,
                row.minimum_time_s,
            )
            _write_expression_source(io, row.expression)
        end
    end
    return abspath(String(path))
end

function run_deck_emt_text_output_suite(
    parsed::DeckParser.DeckParseResult,
    output_dir::AbstractString;
    basename::AbstractString = "deck_text_output",
    title::AbstractString = "AIMORA electromagnetic output",
    runtime_kwargs...,
)
    execution = run_deck_emt_execution(
        parsed;
        time_horizon = :deck,
        output_schedule = :print_and_plot,
        runtime_kwargs...,
    )
    components = _deck_emt_report_components(parsed, execution.trace; title)
    suite = electromagnetic_text_output_suite(parsed, execution; title)
    suite.physical_checks_passed ||
        throw(ErrorException("electromagnetic text-output suite failed physical checks"))
    output_path = abspath(String(output_dir))
    mkpath(output_path)
    report_path = write_electromagnetic_report_output_fixed_width_text(
        joinpath(output_path, string(basename, ".out")),
        suite.report_table;
        column_width = components.column_width,
    )
    punch_path = suite.terminal_punch_requested ?
        write_electromagnetic_terminal_punch(
            joinpath(output_path, string(basename, ".pch")),
            suite.terminal_punch_rows,
        ) : ""
    statistics_path = write_electromagnetic_statistics_text(
        joinpath(output_path, string(basename, ".stats")),
        suite.statistics_rows,
    )
    expression_path = write_electromagnetic_expression_report(
        joinpath(output_path, string(basename, ".xpr")),
        suite.expression_rows,
    )
    completed = ElectromagneticTextOutputSuite(
        suite.source,
        suite.report_table,
        suite.channel_metadata,
        suite.terminal_state,
        suite.terminal_punch_rows,
        suite.statistics_rows,
        suite.expression_rows,
        report_path,
        punch_path,
        statistics_path,
        "",
        expression_path,
        suite.terminal_punch_requested,
        suite.physical_checks_passed,
    )
    spy_path = write_electromagnetic_spy_transcript(
        joinpath(output_path, string(basename, ".spy")),
        completed,
    )
    return ElectromagneticTextOutputSuite(
        completed.source,
        completed.report_table,
        completed.channel_metadata,
        completed.terminal_state,
        completed.terminal_punch_rows,
        completed.statistics_rows,
        completed.expression_rows,
        completed.report_path,
        completed.terminal_punch_path,
        completed.statistics_path,
        spy_path,
        completed.expression_report_path,
        completed.terminal_punch_requested,
        completed.physical_checks_passed,
    )
end

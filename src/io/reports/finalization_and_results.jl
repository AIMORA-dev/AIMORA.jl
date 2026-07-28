function emt_trace_report_artifact(
    trace::DeckEMTTrace;
    title::AbstractString="AIMORA EMT trace report",
)
    headers = _trace_report_headers(trace)
    table = _report_table(trace)
    sample_indices = collect(0:(size(table, 1) - 1))
    final_voltages = Dict{Symbol,Float64}(
        name => final_voltage_pu(trace, name) for name in trace.node_names
    )
    final_outputs = Dict{Symbol,Float64}(
        channel => final_output_pu(trace, channel) for channel in trace.output_channel_names
    )
    return EMTTraceReportArtifact(
        String(title),
        trace.source,
        copy(trace.node_names),
        copy(trace.output_channel_names),
        headers,
        table,
        sample_indices,
        copy(trace.time_s),
        final_voltages,
        final_outputs,
        (
            :over16_output_writer_labels_1643_1654,
            :over16_fdcinj_history_current_trace_channels,
        ),
        false,
        false,
        (
            :full_over29_overlay,
            :full_over31_plot_setup,
            :dekplt_interactive_plotting,
            :dekspy_interactive_storage,
            :binary_plot_file_compatibility,
            :full_bpa_text_report_compatibility,
            :xpr_expression_report_formatting,
        ),
    )
end
function emt_trace_report_result(
    artifact::EMTTraceReportArtifact;
    csv_path::AbstractString="",
    text_path::AbstractString="",
    svg_path::AbstractString="",
    manifest_path::AbstractString="",
)
    quantities = ResultQuantity[
        result_quantity(:samples, size(artifact.table, 1); unit = "count", description = "Recorded fixed-step EMT samples."),
        result_quantity(:report_columns, size(artifact.table, 2); unit = "count", description = "CSV/text report columns including time."),
        result_quantity(:trace_sample_index_count, length(artifact.sample_indices); unit = "count", description = "Explicit trace sample index count."),
        result_quantity(:trace_time_value_count, length(artifact.time_s_values); unit = "count", description = "Explicit trace time value count."),
        result_quantity(:trace_node_column_count, length(artifact.node_names); unit = "count", description = "Node-voltage trace column count."),
        result_quantity(:trace_output_column_count, length(artifact.output_channel_names); unit = "count", description = "Accepted non-voltage trace output column count."),
        result_quantity(:trace_table_value_count, length(artifact.table); unit = "count", description = "Total numeric values in the trace report table."),
        result_quantity(:node_count, length(artifact.node_names); unit = "count", description = "Reported node count."),
        result_quantity(:output_channel_count, length(artifact.output_channel_names); unit = "count", description = "Reported accepted non-voltage trace channel count."),
    ]
    for node in artifact.node_names
        push!(
            quantities,
            result_quantity(
                Symbol("final_", String(node), "_v_pu"),
                artifact.final_voltages_pu[node];
                unit = "pu",
                description = "Final reported node voltage for $(String(node)).",
            ),
        )
    end
    for channel in artifact.output_channel_names
        push!(
            quantities,
            result_quantity(
                Symbol("final_", String(channel)),
                artifact.final_outputs_pu[channel];
                unit = "pu",
                description = "Final reported accepted output channel $(String(channel)).",
            ),
        )
    end

    return study_result(
        :emt;
        status = :warning,
        quantities = quantities,
        assumptions = [
            study_assumption(:artifact_boundary, "DeckEMTTrace"; description = "Report artifacts consume the accepted Julia fixed-step trace."),
            study_assumption(:legacy_fortran_in_loop, artifact.legacy_fortran_in_loop; description = "The legacy Fortran executable is not called by the artifact writer."),
            study_assumption(:full_bpa_report_compatibility, artifact.full_bpa_report_compatibility; description = "The text report is validation-oriented and not full BPA report compatibility."),
        ],
        warnings = [
            study_warning(:prototype_report_artifact, "CSV/SVG/text artifacts cover only the accepted Julia trace subset."),
        ],
        metadata = Dict{Symbol,Any}(
            :engine => "Julia EMT report artifact boundary",
            :source => artifact.source,
            :title => artifact.title,
            :node_names => copy(artifact.node_names),
            :output_channel_names => copy(artifact.output_channel_names),
            :headers => copy(artifact.headers),
            :sample_indices => copy(artifact.sample_indices),
            :time_s_values => copy(artifact.time_s_values),
            :trace_table_row_count => size(artifact.table, 1),
            :trace_table_column_count => size(artifact.table, 2),
            :trace_node_column_count => length(artifact.node_names),
            :trace_output_column_count => length(artifact.output_channel_names),
            :final_voltages_pu => copy(artifact.final_voltages_pu),
            :final_outputs_pu => copy(artifact.final_outputs_pu),
            :legacy_fortran_in_loop => artifact.legacy_fortran_in_loop,
            :full_bpa_report_compatibility => artifact.full_bpa_report_compatibility,
            :csv_path => String(csv_path),
            :text_path => String(text_path),
            :svg_path => String(svg_path),
            :manifest_path => String(manifest_path),
            :fortran_scope => collect(artifact.fortran_scope),
            :deferred_effects => collect(artifact.deferred_effects),
        ),
    )
end

function emt_finalization_report_artifact(
    trace::DeckEMTTrace;
    title::AbstractString="AIMORA OVER20 finalization report",
    terminal_conditions_requested::Bool=true,
    binary_plot_record_requested::Bool=true,
    interactive_plot_requested::Bool=false,
    catalog_requested::Bool=false,
    table_save_requested::Bool=false,
)
    isempty(trace.time_s) && throw(ArgumentError("OVER20 finalization artifact requires at least one trace sample"))
    isempty(trace.node_names) && throw(ArgumentError("OVER20 finalization artifact requires at least one node voltage"))
    headers = _trace_report_headers(trace)
    table = _report_table(trace)
    sample_indices = collect(0:(size(table, 1) - 1))
    variable_names = Symbol.(headers[2:end])

    node_count = length(trace.node_names)
    output_count = length(trace.output_channel_names)
    length(trace.node_maximum_values) == node_count &&
        length(trace.node_maximum_times_s) == node_count &&
        length(trace.node_minimum_values) == node_count &&
        length(trace.node_minimum_times_s) == node_count ||
        throw(ArgumentError("trace node extrema must cover every node"))
    length(trace.output_maximum_values) == output_count &&
        length(trace.output_maximum_times_s) == output_count &&
        length(trace.output_minimum_values) == output_count &&
        length(trace.output_minimum_times_s) == output_count ||
        throw(ArgumentError("trace output extrema must cover every output channel"))

    node_peak_magnitudes = max.(
        abs.(trace.node_maximum_values),
        abs.(trace.node_minimum_values),
    )
    peak_node_index = argmax(node_peak_magnitudes)
    peak_node_name = trace.node_names[peak_node_index]
    maximum_is_peak =
        abs(trace.node_maximum_values[peak_node_index]) >=
        abs(trace.node_minimum_values[peak_node_index])
    peak_node_voltage =
        maximum_is_peak ?
        trace.node_maximum_values[peak_node_index] :
        trace.node_minimum_values[peak_node_index]
    peak_time =
        maximum_is_peak ?
        trace.node_maximum_times_s[peak_node_index] :
        trace.node_minimum_times_s[peak_node_index]

    variable_count = length(variable_names)
    variable_maxima = vcat(
        copy(trace.node_maximum_values),
        copy(trace.output_maximum_values),
    )
    variable_maxima_time = vcat(
        copy(trace.node_maximum_times_s),
        copy(trace.output_maximum_times_s),
    )
    variable_minima = vcat(
        copy(trace.node_minimum_values),
        copy(trace.output_minimum_values),
    )
    variable_minima_time = vcat(
        copy(trace.node_minimum_times_s),
        copy(trace.output_minimum_times_s),
    )
    length(variable_maxima) == variable_count ||
        throw(ArgumentError("trace extrema count must match report variable count"))

    terminal_node_voltages = Dict{Symbol,Float64}(
        node => final_voltage_pu(trace, node) for node in trace.node_names
    )

    return EMTFinalizationReportArtifact(
        String(title),
        trace.source,
        copy(trace.node_names),
        copy(trace.output_channel_names),
        variable_names,
        headers,
        table,
        sample_indices,
        copy(trace.time_s),
        trace.time_s[end],
        terminal_node_voltages,
        peak_node_name,
        peak_node_voltage,
        peak_time,
        variable_maxima,
        variable_maxima_time,
        variable_minima,
        variable_minima_time,
        -9999.0,
        terminal_conditions_requested,
        binary_plot_record_requested,
        interactive_plot_requested,
        catalog_requested,
        table_save_requested,
        false,
        false,
        false,
        false,
        false,
        (
            :over20_peak_node_report_labels_5011_5019,
            :over20_plot_end_sentinel_labels_5019_5022_8005,
            :over20_extrema_report_labels_8002_8004_18005,
            :over20_terminal_node_voltage_labels_5759_7012,
            :over20_final_route_labels_9800_9850,
        ),
        (
            :full_terminal_branch_state_punch,
            :full_terminal_nonlinear_switch_state_punch,
            :binary_plot_file_compatibility,
            :pltfil_interactive_plotting,
            :katalg_tables_save_restore,
            :emtspy_spying_interactive_storage,
            :over20_statistics_route,
            :full_bpa_text_report_compatibility,
            :full_over29_overlay,
            :full_over31_plot_setup,
            :dekplt_interactive_plotting,
            :dekspy_interactive_storage,
            :xpr_expression_report_formatting,
        ),
    )
end

function emt_finalization_report_result(
    artifact::EMTFinalizationReportArtifact;
    csv_path::AbstractString="",
    text_path::AbstractString="",
    manifest_path::AbstractString="",
)
    quantities = ResultQuantity[
        result_quantity(:samples, size(artifact.table, 1); unit = "count", description = "Recorded fixed-step EMT samples consumed by the finalization report."),
        result_quantity(:report_variable_count, length(artifact.variable_names); unit = "count", description = "OVER20-style extrema report variable count."),
        result_quantity(:terminal_node_count, length(artifact.node_names); unit = "count", description = "Terminal node voltage count."),
        result_quantity(:over20_sample_index_count, length(artifact.sample_indices); unit = "count", description = "Explicit OVER20 sample index count."),
        result_quantity(:over20_time_value_count, length(artifact.time_s_values); unit = "count", description = "Explicit OVER20 time value count."),
        result_quantity(:over20_terminal_time_s, artifact.terminal_time_s; unit = "s", description = "Terminal trace time represented by the OVER20 finalization artifact."),
        result_quantity(:over20_variable_extrema_row_count, 2 * length(artifact.variable_names); unit = "count", description = "OVER20 variable extrema rows represented by the artifact."),
        result_quantity(
            :over20_final_route_request_count,
            count(
                identity,
                (
                    artifact.terminal_conditions_requested,
                    artifact.binary_plot_record_requested,
                    artifact.interactive_plot_requested,
                    artifact.catalog_requested,
                    artifact.table_save_requested,
                ),
            );
            unit = "count",
            description = "Requested OVER20 final-route options represented by the artifact.",
        ),
        result_quantity(:peak_node_voltage_pu, artifact.peak_node_voltage_pu; unit = "pu", description = "Overall peak node voltage value from accepted trace samples."),
        result_quantity(:peak_node_time_s, artifact.peak_node_time_s; unit = "s", description = "Trace time of the overall peak node voltage."),
        result_quantity(:extrema_row_count, 2 * length(artifact.variable_names); unit = "count", description = "Maximum plus minimum extrema rows."),
        result_quantity(:plot_end_sentinel_present, artifact.plot_end_sentinel == -9999.0 ? 1 : 0; unit = "bool", description = "Whether the OVER20 plot end sentinel is represented."),
        result_quantity(:terminal_conditions_requested, artifact.terminal_conditions_requested ? 1 : 0; unit = "bool", description = "Whether terminal node-voltage conditions are represented."),
        result_quantity(:binary_plot_record_requested, artifact.binary_plot_record_requested ? 1 : 0; unit = "bool", description = "Whether the legacy binary plot-end record path is requested."),
        result_quantity(:interactive_plot_requested, artifact.interactive_plot_requested ? 1 : 0; unit = "bool", description = "Whether the interactive PLTFIL plot path is requested."),
        result_quantity(:catalog_requested, artifact.catalog_requested ? 1 : 0; unit = "bool", description = "Whether the KATALG catalog path is requested."),
        result_quantity(:table_save_requested, artifact.table_save_requested ? 1 : 0; unit = "bool", description = "Whether TABLES save/restore is requested."),
        result_quantity(:full_bpa_report_compatibility, artifact.full_bpa_report_compatibility ? 1 : 0; unit = "bool", description = "Whether this artifact is full BPA text report compatibility."),
        result_quantity(:binary_plot_format_compatibility, artifact.binary_plot_format_compatibility ? 1 : 0; unit = "bool", description = "Whether this artifact writes legacy binary plot format."),
        result_quantity(:emtspy_compatibility, artifact.emtspy_compatibility ? 1 : 0; unit = "bool", description = "Whether this artifact executes EMTSPY/SPYING compatibility."),
        result_quantity(:statistics_route_compatibility, artifact.statistics_route_compatibility ? 1 : 0; unit = "bool", description = "Whether this artifact executes the NENERG statistics route."),
        result_quantity(:legacy_fortran_in_loop, artifact.legacy_fortran_in_loop ? 1 : 0; unit = "bool", description = "Whether legacy Fortran is called by this artifact."),
    ]
    for node in artifact.node_names
        push!(
            quantities,
            result_quantity(
                Symbol("terminal_", String(node), "_v_pu"),
                artifact.terminal_node_voltages_pu[node];
                unit = "pu",
                description = "Terminal node voltage for $(String(node)).",
            ),
        )
    end

    return study_result(
        :emt;
        status = :warning,
        quantities = quantities,
        assumptions = [
            study_assumption(:artifact_boundary, "DeckEMTTrace"; description = "OVER20 finalization artifacts consume the accepted Julia fixed-step trace."),
            study_assumption(:legacy_fortran_in_loop, artifact.legacy_fortran_in_loop; description = "The legacy Fortran executable is not called by the artifact writer."),
            study_assumption(:full_bpa_report_compatibility, artifact.full_bpa_report_compatibility; description = "The text report is validation-oriented and not full BPA report compatibility."),
        ],
        warnings = [
            study_warning(:bounded_over20_finalization_artifact, "OVER20 finalization artifacts cover peak node, plot sentinel, extrema, final route, and terminal node-voltage rows only."),
        ],
        metadata = Dict{Symbol,Any}(
            :engine => "Julia OVER20 finalization report artifact boundary",
            :source => artifact.source,
            :title => artifact.title,
            :node_names => copy(artifact.node_names),
            :output_channel_names => copy(artifact.output_channel_names),
            :variable_names => copy(artifact.variable_names),
            :headers => copy(artifact.headers),
            :sample_indices => copy(artifact.sample_indices),
            :time_s_values => copy(artifact.time_s_values),
            :terminal_time_s => artifact.terminal_time_s,
            :variable_extrema_row_count => 2 * length(artifact.variable_names),
            :peak_node_name => artifact.peak_node_name,
            :peak_node_voltage_pu => artifact.peak_node_voltage_pu,
            :peak_node_time_s => artifact.peak_node_time_s,
            :terminal_node_voltages_pu => copy(artifact.terminal_node_voltages_pu),
            :variable_maxima_pu => copy(artifact.variable_maxima_pu),
            :variable_maxima_time_s => copy(artifact.variable_maxima_time_s),
            :variable_minima_pu => copy(artifact.variable_minima_pu),
            :variable_minima_time_s => copy(artifact.variable_minima_time_s),
            :plot_end_sentinel => artifact.plot_end_sentinel,
            :terminal_conditions_requested => artifact.terminal_conditions_requested,
            :binary_plot_record_requested => artifact.binary_plot_record_requested,
            :interactive_plot_requested => artifact.interactive_plot_requested,
            :catalog_requested => artifact.catalog_requested,
            :table_save_requested => artifact.table_save_requested,
            :final_route_requested_flags => Dict{Symbol,Bool}(
                :terminal_conditions => artifact.terminal_conditions_requested,
                :binary_plot_record => artifact.binary_plot_record_requested,
                :interactive_plot => artifact.interactive_plot_requested,
                :catalog => artifact.catalog_requested,
                :table_save => artifact.table_save_requested,
            ),
            :full_bpa_report_compatibility => artifact.full_bpa_report_compatibility,
            :binary_plot_format_compatibility => artifact.binary_plot_format_compatibility,
            :emtspy_compatibility => artifact.emtspy_compatibility,
            :statistics_route_compatibility => artifact.statistics_route_compatibility,
            :csv_path => String(csv_path),
            :text_path => String(text_path),
            :manifest_path => String(manifest_path),
            :fortran_scope => collect(artifact.fortran_scope),
            :deferred_effects => collect(artifact.deferred_effects),
        ),
    )
end

function write_deck_trace_csv(path::AbstractString, artifact::EMTTraceReportArtifact)
    _report_ensure_dir(path)
    open(path, "w") do io
        println(io, join(artifact.headers, ","))
        for row in axes(artifact.table, 1)
            @printf(io, "%.9f", artifact.table[row, 1])
            for col in 2:size(artifact.table, 2)
                @printf(io, ",%.9f", artifact.table[row, col])
            end
            println(io)
        end
    end
    return path
end

function write_deck_trace_csv(path::AbstractString, trace::DeckEMTTrace; title::AbstractString="AIMORA EMT trace report")
    return write_deck_trace_csv(path, emt_trace_report_artifact(trace; title = title))
end

function write_deck_trace_emtp_text(path::AbstractString, artifact::EMTTraceReportArtifact)
    _report_ensure_dir(path)
    open(path, "w") do io
        println(io, "AIMORA EMTP-COMPATIBLE TRACE REPORT")
        println(io, "SOURCE: $(artifact.source)")
        println(io, "TITLE: $(artifact.title)")
        println(io, "LEGACY FORTRAN IN LOOP: NO")
        println(io, "FULL BPA REPORT COMPATIBILITY: NO")
        println(io, "FORTRAN SCOPE:")
        for scope in artifact.fortran_scope
            println(io, "  - $(scope)")
        end
        println(io, "TRACE OWNER SUMMARY:")
        println(io, "  SAMPLE INDICES: $(_report_json_int_array(artifact.sample_indices))")
        println(io, "  TIME S VALUES: $(_report_json_float_array(artifact.time_s_values))")
        println(io, "  TRACE TABLE ROWS: $(size(artifact.table, 1))")
        println(io, "  TRACE TABLE COLUMNS: $(size(artifact.table, 2))")
        println(io, "  TRACE NODE COLUMNS: $(length(artifact.node_names))")
        println(io, "  TRACE OUTPUT COLUMNS: $(length(artifact.output_channel_names))")
        println(io, "OUTPUT COLUMNS:")
        println(io, join(artifact.headers, ","))
        print(io, @sprintf("%8s %13s", "STEP", "TIME_S"))
        for header in artifact.headers[2:end]
            print(io, @sprintf(" %13s", uppercase(header)))
        end
        println(io)
        for row in axes(artifact.table, 1)
            print(io, @sprintf("%8d %13.6E", row - 1, artifact.table[row, 1]))
            for col in 2:size(artifact.table, 2)
                print(io, @sprintf(" %13.6E", artifact.table[row, col]))
            end
            println(io)
        end
        println(io, "FINAL VOLTAGES PU:")
        for name in artifact.node_names
            @printf(io, "  %-24s %13.6E\n", uppercase(String(name)), artifact.final_voltages_pu[name])
        end
        println(io, "FINAL OUTPUTS PU:")
        for channel in artifact.output_channel_names
            @printf(io, "  %-24s %13.6E\n", uppercase(String(channel)), artifact.final_outputs_pu[channel])
        end
        println(io, "DEFERRED EFFECTS:")
        for effect in artifact.deferred_effects
            println(io, "  - $(effect)")
        end
    end
    return path
end

function write_deck_trace_emtp_text(path::AbstractString, trace::DeckEMTTrace; title::AbstractString="AIMORA EMT trace report")
    return write_deck_trace_emtp_text(path, emt_trace_report_artifact(trace; title = title))
end

function write_deck_trace_svg(path::AbstractString, artifact::EMTTraceReportArtifact)
    _report_ensure_dir(path)
    width = 1200.0
    height = 680.0
    left = 84.0
    right = 34.0
    top = 76.0
    bottom = 86.0
    plot_w = width - left - right
    plot_h = height - top - bottom
    xs = artifact.table[:, 1] .* 1000.0
    ys = artifact.table[:, 2:end]
    xmin = minimum(xs)
    xmax = maximum(xs)
    ymin, ymax = _nice_limits(vec(ys))
    xscale(x) = left + (x - xmin) / max(xmax - xmin, 1.0e-9) * plot_w
    yscale(y) = top + (ymax - y) / max(ymax - ymin, 1.0e-9) * plot_h
    colors = ["#176b87", "#9a4d9e", "#5f7f52", "#b86b2b", "#4d648d", "#7b4f7f"]

    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$(Int(width))" height="$(Int(height))" viewBox="0 0 $(Int(width)) $(Int(height))">""")
        println(io, """<rect width="100%" height="100%" fill="#f7f7f4"/>""")
        println(io, """<rect x="$left" y="$top" width="$plot_w" height="$plot_h" fill="#ffffff" stroke="#c8c8c0"/>""")
        for i in 0:5
            y = top + i * plot_h / 5
            value = ymax - i * (ymax - ymin) / 5
            println(io, """<line x1="$left" y1="$y" x2="$(left + plot_w)" y2="$y" stroke="#ddddd5" stroke-width="1"/>""")
            @printf(io, "<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"end\" font-family=\"Arial\" font-size=\"13\" fill=\"#444\">%.3f</text>\n", left - 10, y + 4, value)
        end
        for i in 0:6
            x = left + i * plot_w / 6
            value = xmin + i * (xmax - xmin) / 6
            println(io, """<line x1="$x" y1="$top" x2="$x" y2="$(top + plot_h)" stroke="#eeeeea" stroke-width="1"/>""")
            @printf(io, "<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"13\" fill=\"#444\">%.3f</text>\n", x, top + plot_h + 28, value)
        end
        legend_x = left + plot_w - 230
        legend_y = top + 26
        channel_labels = artifact.headers[2:end]
        for (index, name) in enumerate(channel_labels)
            color = colors[mod1(index, length(colors))]
            y = legend_y + (index - 1) * 24
            println(io, """<line x1="$legend_x" y1="$y" x2="$(legend_x + 34)" y2="$y" stroke="$color" stroke-width="3"/>""")
            println(io, """<text x="$(legend_x + 44)" y="$(y + 5)" font-family="Arial" font-size="14" fill="#222">$(_report_escape(name))</text>""")
        end
        for (index, name) in enumerate(channel_labels)
            color = colors[mod1(index, length(colors))]
            parts = String[]
            for row in axes(artifact.table, 1)
                push!(parts, @sprintf("%.2f,%.2f", xscale(xs[row]), yscale(artifact.table[row, index + 1])))
            end
            println(io, """<polyline fill="none" stroke="$color" stroke-width="2.5" points="$(join(parts, " "))"/>""")
        end
        println(io, """<text x="$(width / 2)" y="36" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700" fill="#1f1f1d">$(_report_escape(artifact.title))</text>""")
        println(io, """<text x="$(width / 2)" y="$(height - 24)" text-anchor="middle" font-family="Arial" font-size="16" fill="#333">Time (ms)</text>""")
        println(io, """<text transform="translate(24 $(height / 2)) rotate(-90)" text-anchor="middle" font-family="Arial" font-size="16" fill="#333">Trace value (pu)</text>""")
        println(io, "</svg>")
    end
    return path
end

function write_deck_trace_svg(path::AbstractString, trace::DeckEMTTrace; title::AbstractString="AIMORA EMT trace report")
    return write_deck_trace_svg(path, emt_trace_report_artifact(trace; title = title))
end

function write_deck_trace_report_manifest(
    path::AbstractString,
    artifact::EMTTraceReportArtifact;
    csv_path::AbstractString="",
    text_path::AbstractString="",
    svg_path::AbstractString="",
    finalization_csv_path::AbstractString="",
    finalization_text_path::AbstractString="",
    finalization_manifest_path::AbstractString="",
)
    _report_ensure_dir(path)
    open(path, "w") do io
        println(io, "{")
        @printf(io, "  \"schema\": \"aimora.emt_trace_report_artifact.v1\",\n")
        @printf(io, "  \"title\": \"%s\",\n", _report_json_escape(artifact.title))
        @printf(io, "  \"source\": \"%s\",\n", _report_json_escape(artifact.source))
        @printf(io, "  \"samples\": %d,\n", size(artifact.table, 1))
        @printf(io, "  \"report_columns\": %d,\n", size(artifact.table, 2))
        @printf(io, "  \"sample_indices\": %s,\n", _report_json_int_array(artifact.sample_indices))
        @printf(io, "  \"time_s_values\": %s,\n", _report_json_float_array(artifact.time_s_values))
        @printf(io, "  \"trace_table_row_count\": %d,\n", size(artifact.table, 1))
        @printf(io, "  \"trace_table_column_count\": %d,\n", size(artifact.table, 2))
        @printf(io, "  \"trace_node_column_count\": %d,\n", length(artifact.node_names))
        @printf(io, "  \"trace_output_column_count\": %d,\n", length(artifact.output_channel_names))
        @printf(io, "  \"node_names\": %s,\n", _report_json_string_array(artifact.node_names))
        @printf(io, "  \"output_channel_names\": %s,\n", _report_json_string_array(artifact.output_channel_names))
        @printf(io, "  \"final_voltage_names\": %s,\n", _report_json_string_array(artifact.node_names))
        @printf(io, "  \"final_output_names\": %s,\n", _report_json_string_array(artifact.output_channel_names))
        @printf(io, "  \"headers\": %s,\n", _report_json_string_array(artifact.headers))
        @printf(io, "  \"paths\": {\n")
        @printf(io, "    \"csv_path\": \"%s\",\n", _report_json_escape(csv_path))
        @printf(io, "    \"text_path\": \"%s\",\n", _report_json_escape(text_path))
        @printf(io, "    \"svg_path\": \"%s\",\n", _report_json_escape(svg_path))
        @printf(io, "    \"manifest_path\": \"%s\"\n", _report_json_escape(path))
        @printf(io, "  },\n")
        @printf(io, "  \"linked_finalization\": {\n")
        @printf(io, "    \"csv_path\": \"%s\",\n", _report_json_escape(finalization_csv_path))
        @printf(io, "    \"text_path\": \"%s\",\n", _report_json_escape(finalization_text_path))
        @printf(io, "    \"manifest_path\": \"%s\"\n", _report_json_escape(finalization_manifest_path))
        @printf(io, "  },\n")
        @printf(io, "  \"final_voltages_pu\": %s,\n", _report_json_symbol_float_dict(artifact.final_voltages_pu))
        @printf(io, "  \"final_outputs_pu\": %s,\n", _report_json_symbol_float_dict(artifact.final_outputs_pu))
        @printf(io, "  \"legacy_fortran_in_loop\": %s,\n", _report_json_bool(artifact.legacy_fortran_in_loop))
        @printf(io, "  \"full_bpa_report_compatibility\": %s,\n", _report_json_bool(artifact.full_bpa_report_compatibility))
        @printf(io, "  \"fortran_scope\": %s,\n", _report_json_string_array(artifact.fortran_scope))
        @printf(io, "  \"deferred_effects\": %s\n", _report_json_string_array(artifact.deferred_effects))
        println(io, "}")
    end
    return path
end

function write_deck_trace_report_manifest(path::AbstractString, trace::DeckEMTTrace; title::AbstractString="AIMORA EMT trace report")
    return write_deck_trace_report_manifest(path, emt_trace_report_artifact(trace; title = title))
end

function write_deck_trace_report_artifacts(
    output_dir::AbstractString,
    trace::DeckEMTTrace;
    basename::AbstractString="deck_trace",
    title::AbstractString="AIMORA EMT trace report",
)
    mkpath(output_dir)
    artifact = emt_trace_report_artifact(trace; title = title)
    csv_path = write_deck_trace_csv(joinpath(output_dir, string(basename, ".csv")), artifact)
    text_path = write_deck_trace_emtp_text(joinpath(output_dir, string(basename, "_emtp_report.txt")), artifact)
    svg_path = write_deck_trace_svg(joinpath(output_dir, string(basename, ".svg")), artifact)
    finalization = write_emt_finalization_report_artifacts(
        output_dir,
        trace;
        basename = string(basename, "_over20_finalization"),
        title = string(title, " OVER20 Finalization"),
    )
    manifest_path = write_deck_trace_report_manifest(
        joinpath(output_dir, string(basename, "_report_manifest.json")),
        artifact;
        csv_path = csv_path,
        text_path = text_path,
        svg_path = svg_path,
        finalization_csv_path = finalization.csv_path,
        finalization_text_path = finalization.text_path,
        finalization_manifest_path = finalization.manifest_path,
    )
    return (
        artifact = artifact,
        csv_path = csv_path,
        text_path = text_path,
        svg_path = svg_path,
        manifest_path = manifest_path,
        finalization_artifact = finalization.artifact,
        finalization_csv_path = finalization.csv_path,
        finalization_text_path = finalization.text_path,
        finalization_manifest_path = finalization.manifest_path,
        finalization_result = finalization.result,
        result = emt_trace_report_result(
            artifact;
            csv_path = csv_path,
            text_path = text_path,
            svg_path = svg_path,
            manifest_path = manifest_path,
        ),
    )
end

function write_emt_finalization_report_csv(path::AbstractString, artifact::EMTFinalizationReportArtifact)
    _report_ensure_dir(path)
    terminal_time = artifact.terminal_time_s
    open(path, "w") do io
        println(io, "kind,name,value,time_s")
        @printf(
            io,
            "peak_node_voltage,%s,%.9f,%.9f\n",
            String(artifact.peak_node_name),
            artifact.peak_node_voltage_pu,
            artifact.peak_node_time_s,
        )
        @printf(io, "plot_end_sentinel,VOLTI(1),%.9f,%.9f\n", artifact.plot_end_sentinel, terminal_time)
        for index in eachindex(artifact.variable_names)
            @printf(
                io,
                "variable_maximum,%s,%.9f,%.9f\n",
                String(artifact.variable_names[index]),
                artifact.variable_maxima_pu[index],
                artifact.variable_maxima_time_s[index],
            )
        end
        for index in eachindex(artifact.variable_names)
            @printf(
                io,
                "variable_minimum,%s,%.9f,%.9f\n",
                String(artifact.variable_names[index]),
                artifact.variable_minima_pu[index],
                artifact.variable_minima_time_s[index],
            )
        end
        for node in artifact.node_names
            @printf(
                io,
                "terminal_voltage,%s,%.9f,%.9f\n",
                String(node),
                artifact.terminal_node_voltages_pu[node],
                terminal_time,
            )
        end
        @printf(io, "final_route,terminal_conditions_requested,%d,%.9f\n", artifact.terminal_conditions_requested ? 1 : 0, terminal_time)
        @printf(io, "final_route,binary_plot_record_requested,%d,%.9f\n", artifact.binary_plot_record_requested ? 1 : 0, terminal_time)
        @printf(io, "final_route,interactive_plot_requested,%d,%.9f\n", artifact.interactive_plot_requested ? 1 : 0, terminal_time)
        @printf(io, "final_route,catalog_requested,%d,%.9f\n", artifact.catalog_requested ? 1 : 0, terminal_time)
        @printf(io, "final_route,table_save_requested,%d,%.9f\n", artifact.table_save_requested ? 1 : 0, terminal_time)
    end
    return path
end

function write_emt_finalization_report_csv(path::AbstractString, trace::DeckEMTTrace; title::AbstractString="AIMORA OVER20 finalization report")
    return write_emt_finalization_report_csv(path, emt_finalization_report_artifact(trace; title = title))
end

function write_emt_finalization_report_text(path::AbstractString, artifact::EMTFinalizationReportArtifact)
    _report_ensure_dir(path)
    open(path, "w") do io
        println(io, "AIMORA OVER20 FINALIZATION REPORT")
        println(io, "SOURCE: $(artifact.source)")
        println(io, "TITLE: $(artifact.title)")
        println(io, "LEGACY FORTRAN IN LOOP: NO")
        println(io, "FULL BPA REPORT COMPATIBILITY: NO")
        println(io, "BINARY PLOT FORMAT COMPATIBILITY: NO")
        println(io, "EMTSPY COMPATIBILITY: NO")
        println(io, "STATISTICS ROUTE COMPATIBILITY: NO")
        println(io, "FORTRAN SCOPE:")
        for scope in artifact.fortran_scope
            println(io, "  - $(scope)")
        end
        println(io, "OVER20 OWNER SUMMARY:")
        println(io, "  SAMPLE INDICES: $(_report_json_int_array(artifact.sample_indices))")
        println(io, "  TIME S VALUES: $(_report_json_float_array(artifact.time_s_values))")
        @printf(io, "  TERMINAL TIME S: %.6E\n", artifact.terminal_time_s)
        println(io, "  VARIABLE EXTREMA ROWS: $(2 * length(artifact.variable_names))")
        println(io, "  FINAL ROUTE REQUESTS:")
        println(io, "    TERMINAL CONDITIONS: $(_report_json_bool(artifact.terminal_conditions_requested))")
        println(io, "    BINARY PLOT RECORD: $(_report_json_bool(artifact.binary_plot_record_requested))")
        println(io, "    INTERACTIVE PLOT: $(_report_json_bool(artifact.interactive_plot_requested))")
        println(io, "    CATALOG: $(_report_json_bool(artifact.catalog_requested))")
        println(io, "    TABLE SAVE: $(_report_json_bool(artifact.table_save_requested))")
        println(io, "OVERALL SIMULATION PEAK NODE VOLTAGE")
        @printf(
            io,
            "  %-24s %16.6E %16.6E\n",
            uppercase(String(artifact.peak_node_name)),
            artifact.peak_node_voltage_pu,
            artifact.peak_node_time_s,
        )
        @printf(io, "PLOT SENTINEL VOLTI(1): %.0f\n", artifact.plot_end_sentinel)
        println(io, "VARIABLE MAXIMA")
        for index in eachindex(artifact.variable_names)
            @printf(
                io,
                "  %-32s %16.6E\n",
                uppercase(String(artifact.variable_names[index])),
                artifact.variable_maxima_pu[index],
            )
        end
        println(io, "TIMES OF MAXIMA")
        for index in eachindex(artifact.variable_names)
            @printf(
                io,
                "  %-32s %16.6E\n",
                uppercase(String(artifact.variable_names[index])),
                artifact.variable_maxima_time_s[index],
            )
        end
        println(io, "VARIABLE MINIMA")
        for index in eachindex(artifact.variable_names)
            @printf(
                io,
                "  %-32s %16.6E\n",
                uppercase(String(artifact.variable_names[index])),
                artifact.variable_minima_pu[index],
            )
        end
        println(io, "TIMES OF MINIMA")
        for index in eachindex(artifact.variable_names)
            @printf(
                io,
                "  %-32s %16.6E\n",
                uppercase(String(artifact.variable_names[index])),
                artifact.variable_minima_time_s[index],
            )
        end
        println(io, "TERMINAL NODE VOLTAGES")
        for node in artifact.node_names
            @printf(io, "  %-24s %16.6E\n", uppercase(String(node)), artifact.terminal_node_voltages_pu[node])
        end
        println(io, "DEFERRED EFFECTS:")
        for effect in artifact.deferred_effects
            println(io, "  - $(effect)")
        end
    end
    return path
end

function write_emt_finalization_report_text(path::AbstractString, trace::DeckEMTTrace; title::AbstractString="AIMORA OVER20 finalization report")
    return write_emt_finalization_report_text(path, emt_finalization_report_artifact(trace; title = title))
end

function write_emt_finalization_report_manifest(
    path::AbstractString,
    artifact::EMTFinalizationReportArtifact;
    csv_path::AbstractString="",
    text_path::AbstractString="",
)
    _report_ensure_dir(path)
    open(path, "w") do io
        println(io, "{")
        @printf(io, "  \"schema\": \"aimora.over20_finalization_report_artifact.v1\",\n")
        @printf(io, "  \"title\": \"%s\",\n", _report_json_escape(artifact.title))
        @printf(io, "  \"source\": \"%s\",\n", _report_json_escape(artifact.source))
        @printf(io, "  \"samples\": %d,\n", size(artifact.table, 1))
        @printf(io, "  \"report_variable_count\": %d,\n", length(artifact.variable_names))
        @printf(io, "  \"terminal_node_count\": %d,\n", length(artifact.node_names))
        @printf(io, "  \"sample_indices\": %s,\n", _report_json_int_array(artifact.sample_indices))
        @printf(io, "  \"time_s_values\": %s,\n", _report_json_float_array(artifact.time_s_values))
        @printf(io, "  \"terminal_time_s\": %.17g,\n", artifact.terminal_time_s)
        @printf(io, "  \"variable_extrema_row_count\": %d,\n", 2 * length(artifact.variable_names))
        @printf(io, "  \"node_names\": %s,\n", _report_json_string_array(artifact.node_names))
        @printf(io, "  \"output_channel_names\": %s,\n", _report_json_string_array(artifact.output_channel_names))
        @printf(io, "  \"variable_names\": %s,\n", _report_json_string_array(artifact.variable_names))
        @printf(io, "  \"headers\": %s,\n", _report_json_string_array(artifact.headers))
        @printf(io, "  \"paths\": {\n")
        @printf(io, "    \"csv_path\": \"%s\",\n", _report_json_escape(csv_path))
        @printf(io, "    \"text_path\": \"%s\",\n", _report_json_escape(text_path))
        @printf(io, "    \"manifest_path\": \"%s\"\n", _report_json_escape(path))
        @printf(io, "  },\n")
        @printf(io, "  \"peak_node\": {\n")
        @printf(io, "    \"name\": \"%s\",\n", _report_json_escape(artifact.peak_node_name))
        @printf(io, "    \"voltage_pu\": %.17g,\n", artifact.peak_node_voltage_pu)
        @printf(io, "    \"time_s\": %.17g\n", artifact.peak_node_time_s)
        @printf(io, "  },\n")
        @printf(io, "  \"plot_end_sentinel\": %.17g,\n", artifact.plot_end_sentinel)
        @printf(io, "  \"terminal_node_voltages_pu\": %s,\n", _report_json_symbol_float_dict(artifact.terminal_node_voltages_pu))
        @printf(io, "  \"variable_extrema\": [\n")
        for index in eachindex(artifact.variable_names)
            @printf(
                io,
                "    {\"name\": \"%s\", \"maximum_pu\": %.17g, \"maximum_time_s\": %.17g, \"minimum_pu\": %.17g, \"minimum_time_s\": %.17g}%s\n",
                _report_json_escape(artifact.variable_names[index]),
                artifact.variable_maxima_pu[index],
                artifact.variable_maxima_time_s[index],
                artifact.variable_minima_pu[index],
                artifact.variable_minima_time_s[index],
                index == length(artifact.variable_names) ? "" : ",",
            )
        end
        @printf(io, "  ],\n")
        @printf(io, "  \"requests\": {\n")
        @printf(io, "    \"terminal_conditions\": %s,\n", _report_json_bool(artifact.terminal_conditions_requested))
        @printf(io, "    \"binary_plot_record\": %s,\n", _report_json_bool(artifact.binary_plot_record_requested))
        @printf(io, "    \"interactive_plot\": %s,\n", _report_json_bool(artifact.interactive_plot_requested))
        @printf(io, "    \"catalog\": %s,\n", _report_json_bool(artifact.catalog_requested))
        @printf(io, "    \"table_save\": %s\n", _report_json_bool(artifact.table_save_requested))
        @printf(io, "  },\n")
        @printf(io, "  \"final_route_requested_flags\": {\n")
        @printf(io, "    \"terminal_conditions\": %s,\n", _report_json_bool(artifact.terminal_conditions_requested))
        @printf(io, "    \"binary_plot_record\": %s,\n", _report_json_bool(artifact.binary_plot_record_requested))
        @printf(io, "    \"interactive_plot\": %s,\n", _report_json_bool(artifact.interactive_plot_requested))
        @printf(io, "    \"catalog\": %s,\n", _report_json_bool(artifact.catalog_requested))
        @printf(io, "    \"table_save\": %s\n", _report_json_bool(artifact.table_save_requested))
        @printf(io, "  },\n")
        @printf(io, "  \"compatibility\": {\n")
        @printf(io, "    \"legacy_fortran_in_loop\": %s,\n", _report_json_bool(artifact.legacy_fortran_in_loop))
        @printf(io, "    \"full_bpa_report\": %s,\n", _report_json_bool(artifact.full_bpa_report_compatibility))
        @printf(io, "    \"binary_plot_format\": %s,\n", _report_json_bool(artifact.binary_plot_format_compatibility))
        @printf(io, "    \"emtspy\": %s,\n", _report_json_bool(artifact.emtspy_compatibility))
        @printf(io, "    \"statistics_route\": %s\n", _report_json_bool(artifact.statistics_route_compatibility))
        @printf(io, "  },\n")
        @printf(io, "  \"fortran_scope\": %s,\n", _report_json_string_array(artifact.fortran_scope))
        @printf(io, "  \"deferred_effects\": %s\n", _report_json_string_array(artifact.deferred_effects))
        println(io, "}")
    end
    return path
end

function write_emt_finalization_report_manifest(path::AbstractString, trace::DeckEMTTrace; title::AbstractString="AIMORA OVER20 finalization report")
    return write_emt_finalization_report_manifest(path, emt_finalization_report_artifact(trace; title = title))
end

function write_emt_finalization_report_artifacts(
    output_dir::AbstractString,
    trace::DeckEMTTrace;
    basename::AbstractString="emt_finalization",
    title::AbstractString="AIMORA OVER20 finalization report",
)
    mkpath(output_dir)
    artifact = emt_finalization_report_artifact(trace; title = title)
    csv_path = write_emt_finalization_report_csv(joinpath(output_dir, string(basename, ".csv")), artifact)
    text_path = write_emt_finalization_report_text(joinpath(output_dir, string(basename, "_report.txt")), artifact)
    manifest_path = write_emt_finalization_report_manifest(
        joinpath(output_dir, string(basename, "_manifest.json")),
        artifact;
        csv_path = csv_path,
        text_path = text_path,
    )
    return (
        artifact = artifact,
        csv_path = csv_path,
        text_path = text_path,
        manifest_path = manifest_path,
        result = emt_finalization_report_result(
            artifact;
            csv_path = csv_path,
            text_path = text_path,
            manifest_path = manifest_path,
        ),
    )
end

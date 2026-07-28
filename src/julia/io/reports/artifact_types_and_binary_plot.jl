
using Printf

using ..EMTStudy:
    DeckEMTTrace,
    final_output_pu,
    final_voltage_pu
using ..StudyCore

export EMTTraceReportArtifact,
       EMTFinalizationReportArtifact,
       ElectromagneticCaseStorageReportRow,
       ElectromagneticControlSystemStorageReportRow,
       ElectromagneticCaseTimingReportRow,
       ElectromagneticCaseStatisticsReport,
       ElectromagneticReportOutputTable,
       ElectromagneticMixedOutputReportTable,
       ElectromagneticBinaryPlotTable,
       ElectromagneticPlotReplaySearch,
       ElectromagneticPlotReplaySession,
       emt_trace_report_artifact,
       emt_trace_report_result,
       emt_finalization_report_artifact,
       emt_finalization_report_result,
       electromagnetic_case_statistics_report,
       electromagnetic_report_output_table,
       electromagnetic_report_channel_class_counts,
       electromagnetic_mixed_output_report_table,
       write_electromagnetic_case_statistics_fixed_width_text,
       write_electromagnetic_report_output_text,
       read_electromagnetic_report_output_text,
       write_electromagnetic_report_output_fixed_width_text,
       electromagnetic_binary_plot_table,
       write_electromagnetic_binary_plot_file,
       read_electromagnetic_binary_plot_file,
       electromagnetic_plot_replay_search,
       electromagnetic_plot_replay_session,
       write_deck_trace_csv,
       write_deck_trace_emtp_text,
       write_deck_trace_svg,
       write_deck_trace_report_manifest,
       write_deck_trace_report_artifacts,
       write_emt_finalization_report_csv,
       write_emt_finalization_report_text,
       write_emt_finalization_report_manifest,
       write_emt_finalization_report_artifacts

struct EMTTraceReportArtifact
    title::String
    source::String
    node_names::Vector{Symbol}
    output_channel_names::Vector{Symbol}
    headers::Vector{String}
    table::Matrix{Float64}
    sample_indices::Vector{Int}
    time_s_values::Vector{Float64}
    final_voltages_pu::Dict{Symbol,Float64}
    final_outputs_pu::Dict{Symbol,Float64}
    fortran_scope::Tuple{Vararg{Symbol}}
    legacy_fortran_in_loop::Bool
    full_bpa_report_compatibility::Bool
    deferred_effects::Tuple{Vararg{Symbol}}
end

struct EMTFinalizationReportArtifact
    title::String
    source::String
    node_names::Vector{Symbol}
    output_channel_names::Vector{Symbol}
    variable_names::Vector{Symbol}
    headers::Vector{String}
    table::Matrix{Float64}
    sample_indices::Vector{Int}
    time_s_values::Vector{Float64}
    terminal_time_s::Float64
    terminal_node_voltages_pu::Dict{Symbol,Float64}
    peak_node_name::Symbol
    peak_node_voltage_pu::Float64
    peak_node_time_s::Float64
    variable_maxima_pu::Vector{Float64}
    variable_maxima_time_s::Vector{Float64}
    variable_minima_pu::Vector{Float64}
    variable_minima_time_s::Vector{Float64}
    plot_end_sentinel::Float64
    terminal_conditions_requested::Bool
    binary_plot_record_requested::Bool
    interactive_plot_requested::Bool
    catalog_requested::Bool
    table_save_requested::Bool
    legacy_fortran_in_loop::Bool
    full_bpa_report_compatibility::Bool
    binary_plot_format_compatibility::Bool
    emtspy_compatibility::Bool
    statistics_route_compatibility::Bool
    fortran_scope::Tuple{Vararg{Symbol}}
    deferred_effects::Tuple{Vararg{Symbol}}
end

struct ElectromagneticCaseStorageReportRow
    index::Int
    description::String
    present_figure::Int
    program_limit::Int
    limit_name::String
end

struct ElectromagneticControlSystemStorageReportRow
    table_number::Int
    present_figure::Int
    program_limit::Int
end

struct ElectromagneticCaseTimingReportRow
    label::String
    cpu_s::Float64
    input_output_s::Float64
    total_s::Float64
end

struct ElectromagneticCaseStatisticsReport
    title::String
    source::String
    storage_rows::Vector{ElectromagneticCaseStorageReportRow}
    control_system_storage_rows::Vector{ElectromagneticControlSystemStorageReportRow}
    timing_rows::Vector{ElectromagneticCaseTimingReportRow}
    timing_totals::NamedTuple{(:cpu_s, :input_output_s, :total_s),Tuple{Float64,Float64,Float64}}
    case_statistics_report_compatibility::Bool
    deferred_effects::Tuple{Vararg{Symbol}}
end

struct ElectromagneticReportOutputTable
    title::String
    source::String
    channel_classes::Vector{Symbol}
    upper_names::Vector{String}
    lower_names::Vector{String}
    labels::Vector{String}
    sample_steps::Vector{Int}
    sample_times_s::Vector{Float64}
    sample_values::Matrix{Float64}
    extrema_kinds::Vector{Symbol}
    extrema_indices::Vector{Int}
    extrema_values::Vector{Float64}
    full_text_report_compatibility::Bool
    deferred_effects::Tuple{Vararg{Symbol}}
end

struct ElectromagneticMixedOutputReportTable
    report::ElectromagneticReportOutputTable
    channel_class_counts::Dict{Symbol,Int}
    sample_count::Int
    extrema_count::Int
    mixed_output_report_compatibility::Bool
    deferred_effects::Tuple{Vararg{Symbol}}
end

function _case_storage_report_row(row::ElectromagneticCaseStorageReportRow)
    return row
end

function _case_storage_report_row(row)
    return ElectromagneticCaseStorageReportRow(
        Int(row.index),
        String(row.description),
        Int(row.present_figure),
        Int(row.program_limit),
        String(row.limit_name),
    )
end

function _control_system_storage_report_row(row::ElectromagneticControlSystemStorageReportRow)
    return row
end

function _control_system_storage_report_row(row)
    return ElectromagneticControlSystemStorageReportRow(
        Int(row.table_number),
        Int(row.present_figure),
        Int(row.program_limit),
    )
end

function _case_timing_report_row(row::ElectromagneticCaseTimingReportRow)
    return row
end

function _case_timing_report_row(row)
    return ElectromagneticCaseTimingReportRow(
        String(row.label),
        Float64(row.cpu_s),
        Float64(row.input_output_s),
        Float64(row.total_s),
    )
end

function _case_timing_totals(totals)
    return (
        cpu_s = Float64(totals.cpu_s),
        input_output_s = Float64(totals.input_output_s),
        total_s = Float64(totals.total_s),
    )
end

function electromagnetic_case_statistics_report(;
    title::AbstractString,
    source::AbstractString,
    storage_rows,
    control_system_storage_rows=ElectromagneticControlSystemStorageReportRow[],
    timing_rows,
    timing_totals,
    case_statistics_report_compatibility::Bool=false,
    deferred_effects::Tuple{Vararg{Symbol}}=(
        :terminal_branch_state_punch,
        :terminal_nonlinear_switch_state_punch,
        :binary_plot_file_compatibility,
        :interactive_plot_replay,
        :broad_report_suite_qualification,
    ),
)
    typed_storage_rows = [_case_storage_report_row(row) for row in storage_rows]
    typed_control_rows =
        [_control_system_storage_report_row(row) for row in control_system_storage_rows]
    typed_timing_rows = [_case_timing_report_row(row) for row in timing_rows]
    isempty(typed_storage_rows) &&
        throw(ArgumentError("case statistics report requires at least one storage row"))
    isempty(typed_timing_rows) &&
        throw(ArgumentError("case statistics report requires at least one timing row"))
    length(unique(row.index for row in typed_storage_rows)) == length(typed_storage_rows) ||
        throw(ArgumentError("case statistics report storage row indices must be unique"))
    length(unique(row.table_number for row in typed_control_rows)) == length(typed_control_rows) ||
        throw(ArgumentError("control-system storage table numbers must be unique"))
    for row in typed_storage_rows
        row.index >= 1 ||
            throw(ArgumentError("case statistics storage row index must be positive"))
        isempty(row.description) &&
            throw(ArgumentError("case statistics storage row description must be nonempty"))
        isempty(row.limit_name) &&
            throw(ArgumentError("case statistics storage row limit name must be nonempty"))
    end
    for row in typed_control_rows
        row.table_number >= 1 ||
            throw(ArgumentError("control-system storage table number must be positive"))
    end
    return ElectromagneticCaseStatisticsReport(
        String(title),
        String(source),
        typed_storage_rows,
        typed_control_rows,
        typed_timing_rows,
        _case_timing_totals(timing_totals),
        case_statistics_report_compatibility,
        deferred_effects,
    )
end

function electromagnetic_report_output_table(;
    title::AbstractString,
    source::AbstractString,
    channel_classes::AbstractVector{Symbol},
    upper_names::AbstractVector{<:AbstractString},
    lower_names::AbstractVector{<:AbstractString},
    labels::AbstractVector{<:AbstractString},
    sample_steps::AbstractVector{<:Integer},
    sample_times_s::AbstractVector{<:Real},
    sample_values::AbstractMatrix{<:Real},
    extrema_kinds::AbstractVector{Symbol}=Symbol[],
    extrema_indices::AbstractVector{<:Integer}=Int[],
    extrema_values::AbstractVector{<:Real}=Float64[],
    full_text_report_compatibility::Bool=false,
    deferred_effects::Tuple{Vararg{Symbol}}=(
        :full_text_report_format_compatibility,
        :binary_plot_format_compatibility,
        :interactive_plot_replay,
    ),
)
    channel_count = length(channel_classes)
    length(upper_names) == channel_count ||
        throw(ArgumentError("upper_names length must match channel_classes"))
    length(lower_names) == channel_count ||
        throw(ArgumentError("lower_names length must match channel_classes"))
    length(labels) == channel_count ||
        throw(ArgumentError("labels length must match channel_classes"))
    sample_count = length(sample_steps)
    length(sample_times_s) == sample_count ||
        throw(ArgumentError("sample_times_s length must match sample_steps"))
    size(sample_values, 1) == sample_count ||
        throw(ArgumentError("sample_values row count must match sample_steps"))
    size(sample_values, 2) == channel_count ||
        throw(ArgumentError("sample_values column count must match channel count"))
    length(extrema_kinds) == length(extrema_indices) == length(extrema_values) ||
        throw(ArgumentError("extrema vectors must have equal length"))
    values = Float64.(sample_values)
    for value in values
        isfinite(value) || throw(ArgumentError("sample_values entries must be finite"))
    end
    times = Float64.(sample_times_s)
    for value in times
        isfinite(value) || throw(ArgumentError("sample_times_s entries must be finite"))
    end
    extrema = Float64.(extrema_values)
    for value in extrema
        isfinite(value) || throw(ArgumentError("extrema_values entries must be finite"))
    end
    return ElectromagneticReportOutputTable(
        String(title),
        String(source),
        collect(channel_classes),
        String.(upper_names),
        String.(lower_names),
        String.(labels),
        Int.(sample_steps),
        times,
        Matrix{Float64}(values),
        collect(extrema_kinds),
        Int.(extrema_indices),
        extrema,
        full_text_report_compatibility,
        deferred_effects,
    )
end

function electromagnetic_report_channel_class_counts(table::ElectromagneticReportOutputTable)
    counts = Dict{Symbol,Int}()
    for kind in table.channel_classes
        counts[kind] = get(counts, kind, 0) + 1
    end
    return counts
end

function electromagnetic_mixed_output_report_table(
    report::ElectromagneticReportOutputTable;
    required_channel_class_counts::Dict{Symbol,<:Integer}=Dict{Symbol,Int}(),
    mixed_output_report_compatibility::Bool=false,
    deferred_effects::Tuple{Vararg{Symbol}}=(
        :full_text_report_format_compatibility,
        :full_control_system_execution,
        :full_machine_solution,
        :julia_waveform_generation,
    ),
)
    counts = electromagnetic_report_channel_class_counts(report)
    length(counts) >= 2 ||
        throw(ArgumentError("mixed-output report requires at least two channel classes"))
    for (kind, expected) in required_channel_class_counts
        get(counts, kind, 0) == Int(expected) ||
            throw(ArgumentError("mixed-output report channel count for $kind does not match"))
    end
    haskey(counts, :tacs) || haskey(counts, :universal_machine) ||
        throw(ArgumentError("mixed-output report requires TACS or universal-machine channels"))
    return ElectromagneticMixedOutputReportTable(
        report,
        counts,
        length(report.sample_steps),
        length(report.extrema_kinds),
        mixed_output_report_compatibility,
        deferred_effects,
    )
end

const _CASE_STORAGE_PRESENT_COLUMN_PREFIX_WIDTH = 103
const _CASE_TIMING_LABEL_WIDTH = 105

function _case_statistics_storage_prefix(index::Integer)
    suffix = index < 10 ? "   " : "  "
    return "     SIZE LIST $(Int(index)).$(suffix)"
end

function _case_statistics_checked_rpad(text::AbstractString, width::Integer, field::Symbol)
    value = String(text)
    length(value) <= width ||
        throw(ArgumentError("case statistics $field exceeds fixed-width report column"))
    return rpad(value, width)
end

function _write_case_statistics_storage_row(io::IO, row::ElectromagneticCaseStorageReportRow)
    prefix = _case_statistics_storage_prefix(row.index)
    description_width = _CASE_STORAGE_PRESENT_COLUMN_PREFIX_WIDTH - length(prefix)
    description =
        _case_statistics_checked_rpad(row.description, description_width, :storage_description)
    @printf(
        io,
        "%s%s%7d%10d (%s)\n",
        prefix,
        description,
        row.present_figure,
        row.program_limit,
        row.limit_name,
    )
end

function _write_case_statistics_integer_row(
    io::IO,
    label::AbstractString,
    values::AbstractVector{<:Integer},
)
    print(io, rpad(String(label), 21))
    for value in values
        @printf(io, "%10d", Int(value))
    end
    println(io)
end

function _write_case_statistics_control_system_rows(
    io::IO,
    rows::AbstractVector{ElectromagneticControlSystemStorageReportRow},
)
    isempty(rows) && return
    _write_case_statistics_integer_row(io, "       TACS TABLE NO.", [row.table_number for row in rows])
    _write_case_statistics_integer_row(io, "       PRESENT FIGURE", [row.present_figure for row in rows])
    _write_case_statistics_integer_row(io, "       PROGRAM LIMIT", [row.program_limit for row in rows])
end

function _write_case_statistics_timing_row(io::IO, row::ElectromagneticCaseTimingReportRow)
    label = _case_statistics_checked_rpad(row.label, _CASE_TIMING_LABEL_WIDTH, :timing_label)
    @printf(io, "%s%5.3f%10.3f%10.3f\n", label, row.cpu_s, row.input_output_s, row.total_s)
end

function _write_case_statistics_timing_totals(io::IO, totals)
    label = rpad(lpad("TOTALS", 99), _CASE_TIMING_LABEL_WIDTH)
    @printf(io, "%s%5.3f%10.3f%10.3f\n", label, totals.cpu_s, totals.input_output_s, totals.total_s)
end

function write_electromagnetic_case_statistics_fixed_width_text(
    path::AbstractString,
    report::ElectromagneticCaseStatisticsReport,
)
    _report_ensure_dir(path)
    open(path, "w") do io
        println(
            io,
            " CORE STORAGE FIGURES FOR PRECEDING DATA CASE NOW COMPLETED.  ---------------------------------------  PRESENT   PROGRAM",
        )
        println(
            io,
            " A VALUE OF  -9999 INDICATES DEFAULT, WITH NO FIGURE AVAILABLE.                                         FIGURE     LIMIT (NAME)",
        )
        for row in report.storage_rows
            _write_case_statistics_storage_row(io, row)
            if row.index == 19
                _write_case_statistics_control_system_rows(io, report.control_system_storage_rows)
            end
        end
        println(
            io,
            " TIMING FIGURES (DECIMAL) CHARACTERIZING CASE SOLUTION SPEED.  -------------------------------------    CP SEC   I/O SEC   SUM SEC",
        )
        for row in report.timing_rows
            _write_case_statistics_timing_row(io, row)
        end
        println(io, "                                                                                                     -----------------------------")
        _write_case_statistics_timing_totals(io, report.timing_totals)
    end
    return path
end

function write_electromagnetic_report_output_text(
    path::AbstractString,
    table::ElectromagneticReportOutputTable,
)
    _report_ensure_dir(path)
    open(path, "w") do io
        println(io, "AIMORA ELECTROMAGNETIC OUTPUT REPORT")
        println(io, "TITLE\t$(table.title)")
        println(io, "SOURCE\t$(table.source)")
        println(io, "OUTPUT_VARIABLES\t$(length(table.channel_classes))")
        println(io, "PRINTED_SAMPLES\t$(length(table.sample_steps))")
        println(io, "EXTREMA_ROWS\t$(length(table.extrema_kinds))")
        println(io, "FULL_TEXT_REPORT_COMPATIBILITY\t$(_report_json_bool(table.full_text_report_compatibility))")
        println(io, "CHANNELS")
        println(io, "index\tclass\tupper_name\tlower_name\tlabel")
        for index in eachindex(table.channel_classes)
            println(
                io,
                join((
                    string(index),
                    String(table.channel_classes[index]),
                    table.upper_names[index],
                    table.lower_names[index],
                    table.labels[index],
                ), '\t'),
            )
        end
        println(io, "SAMPLES")
        print(io, "step\ttime_s")
        for index in eachindex(table.channel_classes)
            @printf(io, "\tvalue_%02d", index)
        end
        println(io)
        for row in eachindex(table.sample_steps)
            @printf(io, "%d\t%.12e", table.sample_steps[row], table.sample_times_s[row])
            for column in axes(table.sample_values, 2)
                @printf(io, "\t%.12e", table.sample_values[row, column])
            end
            println(io)
        end
        println(io, "EXTREMA")
        println(io, "kind\tindex\tvalue")
        for index in eachindex(table.extrema_kinds)
            @printf(
                io,
                "%s\t%d\t%.12e\n",
                String(table.extrema_kinds[index]),
                table.extrema_indices[index],
                table.extrema_values[index],
            )
        end
        println(io, "DEFERRED_EFFECTS")
        for effect in table.deferred_effects
            println(io, String(effect))
        end
    end
    return path
end

function _report_fixed_width_exponential(value::Real)
    number = Float64(value)
    isfinite(number) ||
        throw(ArgumentError("fixed-width report values must be finite"))
    if number == 0.0
        return @sprintf("%13s", "0.000000E+00")
    end
    exponent = floor(Int, log10(abs(number))) + 1
    mantissa = number / 10.0^exponent
    rounded = round(abs(mantissa); digits = 6)
    if rounded >= 1.0
        exponent += 1
        mantissa /= 10.0
    end
    body = @sprintf("%.6fE%+03d", abs(mantissa), exponent)
    text = number < 0.0 ? string("-", body) : body
    return @sprintf("%13s", text)
end

function _write_fixed_width_value_rows(
    io::IO,
    values::AbstractVector{<:Real};
    first_prefix::AbstractString="",
    continuation_prefix::AbstractString=repeat(" ", 15),
    values_per_line::Integer=9,
)
    isempty(values) && return nothing
    start = firstindex(values)
    while start <= lastindex(values)
        stop = min(start + Int(values_per_line) - 1, lastindex(values))
        print(io, start == firstindex(values) ? first_prefix : continuation_prefix)
        for index in start:stop
            print(io, _report_fixed_width_exponential(values[index]))
        end
        println(io)
        start = stop + 1
    end
    return nothing
end

function _report_fixed_width_name(name::AbstractString)
    text = first(String(name), min(ncodeunits(String(name)), 8))
    return @sprintf("%13s", text)
end

function _write_fixed_width_heading_names(
    io::IO,
    table::ElectromagneticReportOutputTable;
    values_per_line::Integer,
)
    start = firstindex(table.upper_names)
    while start <= lastindex(table.upper_names)
        stop = min(start + Int(values_per_line) - 1, lastindex(table.upper_names))
        prefix = start == firstindex(table.upper_names) ? "  STEP     TIME" : repeat(" ", 15)
        println(
            io,
            prefix,
            join((_report_fixed_width_name(table.upper_names[index]) for index in start:stop), ""),
        )
        lower = join(
            (_report_fixed_width_name(table.lower_names[index]) for index in start:stop),
            "",
        )
        isempty(strip(lower)) || println(io, repeat(" ", 15), lower)
        start = stop + 1
    end
    println(io)
    return nothing
end

function _report_extrema_values(table::ElectromagneticReportOutputTable, kind::Symbol)
    values = fill(NaN, length(table.channel_classes))
    for index in eachindex(table.extrema_kinds)
        table.extrema_kinds[index] == kind || continue
        channel_index = table.extrema_indices[index]
        1 <= channel_index <= length(values) ||
            throw(ArgumentError("extrema index out of range"))
        values[channel_index] = table.extrema_values[index]
    end
    all(isfinite, values) ||
        throw(ArgumentError("fixed-width report extrema are incomplete for $kind"))
    return values
end

function write_electromagnetic_report_output_fixed_width_text(
    path::AbstractString,
    table::ElectromagneticReportOutputTable,
    ;
    column_width::Integer = 132,
)
    _report_ensure_dir(path)
    width = Int(column_width)
    width in (80, 132) ||
        throw(ArgumentError("fixed-width report column_width must be 80 or 132"))
    nominal_values_per_line = width == 132 ? 9 : 5
    sample_prefix_width = maximum(
        ncodeunits(@sprintf("%6d %8.6f", step, time_s)) for
        (step, time_s) in zip(table.sample_steps, table.sample_times_s)
    )
    values_per_line = min(
        nominal_values_per_line,
        max(fld(width - sample_prefix_width, 13), 1),
    )
    counts = electromagnetic_report_channel_class_counts(table)
    output_count = length(table.channel_classes)
    open(path, "w") do io
        println(io)
        if width == 132
            @printf(
                io,
                " COLUMN HEADINGS FOR THE%4d  EMTP OUTPUT VARIABLES FOLLOW.   THESE ARE ORDERED ACCORDING TO THE FIVE\n",
                output_count,
            )
            println(io, " POSSIBLE EMTP OUTPUT-VARIABLE CLASSES, AS FOLLOWS ....")
            @printf(io, "   FIRST%4d  OUTPUT VARIABLES ARE ELECTRIC-NETWORK NODE VOLTAGES (WITH RESPECT TO LOCAL GROUND)|\n", get(counts, :node_voltage, 0))
            @printf(io, "    NEXT%4d  OUTPUT VARIABLES ARE BRANCH VOLTAGES (VOLTAGE OF UPPER NODE MINUS VOLTAGE OF LOWER NODE)|\n", get(counts, :branch_voltage, 0))
            @printf(io, "    NEXT%4d  OUTPUT VARIABLES ARE BRANCH CURRENTS (FLOWING FROM THE UPPER EMTP NODE TO THE LOWER)|\n", get(counts, :branch_current, 0))
            @printf(io, "    NEXT%4d  OUTPUT VARIABLES PERTAIN TO DYNAMIC SYNCHRONOUS MACHINES, WITH NAMES GENERATED INTERNALLY|\n", get(counts, :synchronous_machine, 0))
            @printf(io, "   FINAL%4d  OUTPUT VARIABLES BELONG TO  'TACS'  (NOTE INTERNALLY-ADDED UPPER NAME OF PAIR).\n", get(counts, :tacs, 0))
            println(io, " BRANCH POWER  CONSUMPTION (POWER  FLOW, IF A SWITCH) IS TREATED LIKE A BRANCH VOLTAGE FOR THIS GROUPING|")
            println(io, " BRANCH ENERGY CONSUMPTION (ENERGY FLOW, IF A SWITCH) IS TREATED LIKE A BRANCH CURRENT FOR THIS GROUPING.")
        else
            println(io, " TIME-STEP LOOP OUTPUT COUNTS (NODE V, BRANCH V, BRANCH I, S.M., TACS):")
            @printf(
                io,
                " %6d%6d%6d%6d%6d   TOTAL%6d\n",
                get(counts, :node_voltage, 0),
                get(counts, :branch_voltage, 0),
                get(counts, :branch_current, 0),
                get(counts, :synchronous_machine, 0),
                get(counts, :tacs, 0),
                output_count,
            )
        end
        println(io)
        _write_fixed_width_heading_names(io, table; values_per_line = values_per_line)
        for row in eachindex(table.sample_steps)
            _write_fixed_width_value_rows(
                io,
                vec(table.sample_values[row, :]);
                first_prefix = @sprintf("%6d %8.6f", table.sample_steps[row], table.sample_times_s[row]),
                values_per_line = values_per_line,
            )
        end
        println(io)
        if width == 132
            println(io, " MAXIMA AND MINIMA WHICH OCCURRED DURING THE SIMULATION FOLLOW.   THE ORDER AND COLUMN POSITIONING ARE THE")
            println(io, " SAME AS FOR THE REGULAR PRINTED OUTPUT VS. TIME.")
        else
            println(io, " MAXIMA AND MINIMA FOLLOW IN REGULAR OUTPUT COLUMN ORDER.")
        end
        for (title, kind) in (
            ("VARIABLE MAXIMA :", :maxima),
            ("TIMES OF MAXIMA :", :times_of_maxima),
            ("VARIABLE MINIMA :", :minima),
            ("TIMES OF MINIMA :", :times_of_minima),
        )
            println(io, " ", title)
            _write_fixed_width_value_rows(
                io,
                _report_extrema_values(table, kind);
                first_prefix = repeat(" ", 15),
                values_per_line = values_per_line,
            )
        end
    end
    all(ncodeunits(line) <= width for line in eachline(path)) ||
        throw(ErrorException("fixed-width report emitted a line wider than $width columns"))
    return path
end

function _report_output_section(lines::Vector{String}, name::AbstractString)
    start = findfirst(==(String(name)), lines)
    start === nothing && return String[]
    stop = findnext(
        line -> line in ("CHANNELS", "SAMPLES", "EXTREMA", "DEFERRED_EFFECTS"),
        lines,
        start + 1,
    )
    last_index = stop === nothing ? length(lines) : stop - 1
    return lines[(start + 1):last_index]
end

function read_electromagnetic_report_output_text(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && throw(ArgumentError("report output text is empty"))
    lines[1] == "AIMORA ELECTROMAGNETIC OUTPUT REPORT" ||
        throw(ArgumentError("unsupported report output text header"))
    properties = Dict{String,String}()
    for line in lines[2:end]
        line in ("CHANNELS", "SAMPLES", "EXTREMA", "DEFERRED_EFFECTS") && break
        fields = split(line, '\t'; limit = 2)
        length(fields) == 2 && (properties[fields[1]] = fields[2])
    end

    channel_rows = _report_output_section(lines, "CHANNELS")
    sample_rows = _report_output_section(lines, "SAMPLES")
    extrema_rows = _report_output_section(lines, "EXTREMA")
    effect_rows = _report_output_section(lines, "DEFERRED_EFFECTS")
    isempty(channel_rows) && throw(ArgumentError("missing channel rows"))
    isempty(sample_rows) && throw(ArgumentError("missing sample rows"))
    isempty(extrema_rows) && throw(ArgumentError("missing extrema rows"))

    channel_classes = Symbol[]
    upper_names = String[]
    lower_names = String[]
    labels = String[]
    for row in channel_rows[2:end]
        isempty(row) && continue
        fields = split(row, '\t'; keepempty = true)
        length(fields) == 5 || throw(ArgumentError("invalid channel row"))
        push!(channel_classes, Symbol(fields[2]))
        push!(upper_names, fields[3])
        push!(lower_names, fields[4])
        push!(labels, fields[5])
    end

    sample_steps = Int[]
    sample_times = Float64[]
    sample_values = zeros(Float64, 0, length(channel_classes))
    for row in sample_rows[2:end]
        isempty(row) && continue
        fields = split(row, '\t'; keepempty = true)
        length(fields) == length(channel_classes) + 2 ||
            throw(ArgumentError("invalid sample row"))
        push!(sample_steps, parse(Int, fields[1]))
        push!(sample_times, parse(Float64, fields[2]))
        values = [parse(Float64, fields[index + 2]) for index in eachindex(channel_classes)]
        sample_values = vcat(sample_values, reshape(values, 1, length(values)))
    end

    extrema_kinds = Symbol[]
    extrema_indices = Int[]
    extrema_values = Float64[]
    for row in extrema_rows[2:end]
        isempty(row) && continue
        fields = split(row, '\t'; keepempty = true)
        length(fields) == 3 || throw(ArgumentError("invalid extrema row"))
        push!(extrema_kinds, Symbol(fields[1]))
        push!(extrema_indices, parse(Int, fields[2]))
        push!(extrema_values, parse(Float64, fields[3]))
    end
    return electromagnetic_report_output_table(
        title = get(properties, "TITLE", ""),
        source = get(properties, "SOURCE", ""),
        channel_classes = channel_classes,
        upper_names = upper_names,
        lower_names = lower_names,
        labels = labels,
        sample_steps = sample_steps,
        sample_times_s = sample_times,
        sample_values = sample_values,
        extrema_kinds = extrema_kinds,
        extrema_indices = extrema_indices,
        extrema_values = extrema_values,
        full_text_report_compatibility =
            lowercase(get(properties, "FULL_TEXT_REPORT_COMPATIBILITY", "false")) ==
            "true",
        deferred_effects = Tuple(Symbol(strip(row)) for row in effect_rows if !isempty(strip(row))),
    )
end

struct ElectromagneticBinaryPlotTable
    title::String
    source::String
    channel_classes::Vector{Symbol}
    upper_names::Vector{String}
    lower_names::Vector{String}
    sample_times_s::Vector{Float64}
    sample_values::Matrix{Float64}
    sentinel_value::Float64
    sentinel_sample_values::Matrix{Float64}
    sentinel_record_count::Int
    header_record_prefix::String
    name_table::Vector{String}
    fixed_name_table::Vector{String}
    node_name_indices::Vector{Int}
    branch_upper_name_indices::Vector{Int}
    branch_lower_name_indices::Vector{Int}
    binary_plot_format_compatibility::Bool
    deferred_effects::Tuple{Vararg{Symbol}}
end

function _binary_plot_channel_rank(kind::Symbol)
    kind == :node_voltage && return 1
    kind == :branch_voltage && return 2
    kind == :branch_current && return 3
    throw(ArgumentError("unsupported binary plot channel class $kind"))
end

function _binary_plot_name(name::AbstractString)
    text = String(strip(String(name)))
    all(character -> Int(character) < 0x80, text) ||
        throw(ArgumentError("binary plot channel names must be ASCII"))
    ncodeunits(text) <= 8 ||
        throw(ArgumentError("binary plot channel names must be at most 8 bytes"))
    occursin('\0', text) && throw(ArgumentError("binary plot channel names cannot contain NUL"))
    return text
end

function _binary_plot_fixed_name_entry(name::AbstractString)
    text = String(name)
    all(character -> Int(character) < 0x80, text) ||
        throw(ArgumentError("binary plot fixed channel names must be ASCII"))
    ncodeunits(text) == 8 ||
        throw(ArgumentError("binary plot fixed channel names must be exactly 8 bytes"))
    occursin('\0', text) &&
        throw(ArgumentError("binary plot fixed channel names cannot contain NUL"))
    return text
end

function _binary_plot_record_prefix(prefix::AbstractString)
    text = String(prefix)
    all(character -> Int(character) < 0x80, text) ||
        throw(ArgumentError("binary plot header prefix must be ASCII"))
    ncodeunits(text) <= 32 ||
        throw(ArgumentError("binary plot header prefix must be at most 32 bytes"))
    occursin('\0', text) &&
        throw(ArgumentError("binary plot header prefix cannot contain NUL"))
    return rpad(text, 32)[1:32]
end

function _binary_plot_indices(
    values::AbstractVector{<:Integer},
    expected_count::Int,
    table_count::Int,
    field_name::AbstractString,
)
    indices = Int.(values)
    length(indices) == expected_count ||
        throw(ArgumentError("$field_name length must be $expected_count"))
    for index in indices
        1 <= index <= table_count ||
            throw(ArgumentError("$field_name contains an index outside the binary plot name table"))
    end
    return indices
end

function _binary_plot_default_name_table(
    channel_classes::AbstractVector{Symbol},
    upper_names::AbstractVector{String},
    lower_names::AbstractVector{String},
)
    names = String[""]
    indices = Dict{String,Int}("" => 1)
    function index_for(name::String)
        if !haskey(indices, name)
            push!(names, name)
            indices[name] = length(names)
        end
        return indices[name]
    end
    upper_indices = Int[]
    lower_indices = Int[]
    node_indices = Int[]
    for index in eachindex(channel_classes)
        upper_index = index_for(upper_names[index])
        lower_index = index_for(lower_names[index])
        if channel_classes[index] == :node_voltage
            push!(node_indices, upper_index)
        else
            push!(upper_indices, upper_index)
            push!(lower_indices, lower_index)
        end
    end
    return (
        names = names,
        fixed_names = [_binary_plot_fixed_name(name) for name in names],
        node_indices = node_indices,
        upper_indices = upper_indices,
        lower_indices = lower_indices,
    )
end

function _binary_plot_layout(
    channel_classes::AbstractVector{Symbol},
    upper_names::AbstractVector{String},
    lower_names::AbstractVector{String};
    name_table,
    fixed_name_table,
    node_name_indices,
    branch_upper_name_indices,
    branch_lower_name_indices,
)
    node_count = count(==(:node_voltage), channel_classes)
    branch_pair_count = length(channel_classes) - node_count
    all_layout_fields_missing =
        name_table === nothing &&
        node_name_indices === nothing &&
        branch_upper_name_indices === nothing &&
        branch_lower_name_indices === nothing
    all_layout_fields_present =
        name_table !== nothing &&
        node_name_indices !== nothing &&
        branch_upper_name_indices !== nothing &&
        branch_lower_name_indices !== nothing
    if all_layout_fields_missing
        fixed_name_table === nothing ||
            throw(ArgumentError("fixed_name_table requires a binary plot name table"))
        return _binary_plot_default_name_table(channel_classes, upper_names, lower_names)
    elseif !all_layout_fields_present
        throw(ArgumentError("binary plot layout fields must be provided together"))
    end
    names = [_binary_plot_name(name) for name in name_table]
    !isempty(names) || throw(ArgumentError("binary plot name table cannot be empty"))
    fixed_names = if fixed_name_table === nothing
        [_binary_plot_fixed_name(name) for name in names]
    else
        fixed = [_binary_plot_fixed_name_entry(name) for name in fixed_name_table]
        length(fixed) == length(names) ||
            throw(ArgumentError("fixed_name_table length must match name_table"))
        for index in eachindex(names)
            strip(fixed[index]) == names[index] ||
                throw(ArgumentError("fixed_name_table entries must match stripped name_table entries"))
        end
        fixed
    end
    node_indices = _binary_plot_indices(
        node_name_indices,
        node_count,
        length(names),
        "node_name_indices",
    )
    upper_indices = _binary_plot_indices(
        branch_upper_name_indices,
        branch_pair_count,
        length(names),
        "branch_upper_name_indices",
    )
    lower_indices = _binary_plot_indices(
        branch_lower_name_indices,
        branch_pair_count,
        length(names),
        "branch_lower_name_indices",
    )

    node_index = 1
    branch_index = 1
    for channel_index in eachindex(channel_classes)
        if channel_classes[channel_index] == :node_voltage
            names[node_indices[node_index]] == upper_names[channel_index] ||
                throw(ArgumentError("node_name_indices do not match node-voltage channel names"))
            node_index += 1
        else
            names[upper_indices[branch_index]] == upper_names[channel_index] ||
                throw(ArgumentError("branch_upper_name_indices do not match branch channel names"))
            names[lower_indices[branch_index]] == lower_names[channel_index] ||
                throw(ArgumentError("branch_lower_name_indices do not match branch channel names"))
            branch_index += 1
        end
    end
    return (
        names = names,
        fixed_names = fixed_names,
        node_indices = node_indices,
        upper_indices = upper_indices,
        lower_indices = lower_indices,
    )
end

function electromagnetic_binary_plot_table(;
    title::AbstractString,
    source::AbstractString,
    channel_classes::AbstractVector{Symbol},
    upper_names::AbstractVector{<:AbstractString},
    lower_names::AbstractVector{<:AbstractString},
    sample_times_s::AbstractVector{<:Real},
    sample_values::AbstractMatrix{<:Real},
    sentinel_value::Real=-9999.0,
    sentinel_sample_values=nothing,
    sentinel_record_count::Integer=1,
    header_record_prefix::AbstractString="AIMORA BINARY PLOT",
    name_table=nothing,
    fixed_name_table=nothing,
    node_name_indices=nothing,
    branch_upper_name_indices=nothing,
    branch_lower_name_indices=nothing,
    binary_plot_format_compatibility::Bool=true,
    deferred_effects::Tuple{Vararg{Symbol}}=(
        :byte_for_byte_plot_file_parity,
        :interactive_plot_replay,
        :mixed_output_report_equivalence,
    ),
)
    channel_count = length(channel_classes)
    channel_count > 0 || throw(ArgumentError("binary plot table requires at least one channel"))
    length(upper_names) == channel_count ||
        throw(ArgumentError("upper_names length must match channel_classes"))
    length(lower_names) == channel_count ||
        throw(ArgumentError("lower_names length must match channel_classes"))
    size(sample_values, 2) == channel_count ||
        throw(ArgumentError("sample_values column count must match channel count"))
    sample_count = length(sample_times_s)
    sample_count > 0 || throw(ArgumentError("binary plot table requires at least one sample"))
    size(sample_values, 1) == sample_count ||
        throw(ArgumentError("sample_values row count must match sample_times_s"))
    sentinel_count = Int(sentinel_record_count)
    sentinel_count >= 1 || throw(ArgumentError("binary plot table requires a sentinel record"))
    ranks = [_binary_plot_channel_rank(kind) for kind in channel_classes]
    issorted(ranks) ||
        throw(ArgumentError("binary plot channels must be ordered as node voltage, branch voltage, branch current"))
    names_upper = [_binary_plot_name(name) for name in upper_names]
    names_lower = [_binary_plot_name(name) for name in lower_names]
    for (index, kind) in enumerate(channel_classes)
        !isempty(names_upper[index]) ||
            throw(ArgumentError("binary plot upper channel names cannot be empty"))
        kind == :node_voltage && !isempty(names_lower[index]) &&
            throw(ArgumentError("node-voltage binary plot channels cannot have a lower name"))
    end
    times = Float64.(sample_times_s)
    for value in times
        isfinite(value) || throw(ArgumentError("sample_times_s entries must be finite"))
        value != Float64(sentinel_value) ||
            throw(ArgumentError("sample time cannot equal the binary plot sentinel"))
    end
    values = Float64.(sample_values)
    for value in values
        isfinite(value) || throw(ArgumentError("sample_values entries must be finite"))
    end
    sentinel_values = if sentinel_sample_values === nothing
        repeat(reshape(values[end, :], 1, channel_count), sentinel_count, 1)
    else
        Matrix{Float64}(sentinel_sample_values)
    end
    size(sentinel_values) == (sentinel_count, channel_count) ||
        throw(ArgumentError("sentinel_sample_values must have sentinel_record_count rows and channel_count columns"))
    for value in sentinel_values
        isfinite(value) || throw(ArgumentError("sentinel_sample_values entries must be finite"))
    end
    sentinel = Float64(sentinel_value)
    isfinite(sentinel) || throw(ArgumentError("sentinel_value must be finite"))
    layout = _binary_plot_layout(
        channel_classes,
        names_upper,
        names_lower;
        name_table = name_table,
        fixed_name_table = fixed_name_table,
        node_name_indices = node_name_indices,
        branch_upper_name_indices = branch_upper_name_indices,
        branch_lower_name_indices = branch_lower_name_indices,
    )
    return ElectromagneticBinaryPlotTable(
        String(title),
        String(source),
        collect(channel_classes),
        names_upper,
        names_lower,
        times,
        Matrix{Float64}(values),
        sentinel,
        sentinel_values,
        sentinel_count,
        _binary_plot_record_prefix(header_record_prefix),
        layout.names,
        layout.fixed_names,
        layout.node_indices,
        layout.upper_indices,
        layout.lower_indices,
        binary_plot_format_compatibility,
        deferred_effects,
    )
end

function electromagnetic_binary_plot_table(
    report::ElectromagneticReportOutputTable;
    kwargs...,
)
    return electromagnetic_binary_plot_table(;
        title = report.title,
        source = report.source,
        channel_classes = report.channel_classes,
        upper_names = report.upper_names,
        lower_names = report.lower_names,
        sample_times_s = report.sample_times_s,
        sample_values = report.sample_values,
        binary_plot_format_compatibility = report.full_text_report_compatibility,
        kwargs...,
    )
end

function _binary_plot_int32(bytes::Vector{UInt8}, offset::Int)
    offset + 3 <= length(bytes) || throw(ArgumentError("binary plot Int32 is truncated"))
    return Int(reinterpret(Int32, bytes[offset:(offset + 3)])[1])
end

function _binary_plot_records(path::AbstractString)
    bytes = read(path)
    records = Vector{Vector{UInt8}}()
    offset = 1
    while offset <= length(bytes)
        record_length = _binary_plot_int32(bytes, offset)
        record_length >= 0 || throw(ArgumentError("binary plot record length is negative"))
        record_start = offset + 4
        record_stop = record_start + record_length - 1
        trailing_offset = record_stop + 1
        trailing_offset + 3 <= length(bytes) ||
            throw(ArgumentError("binary plot record is truncated before trailing marker"))
        trailing_length = _binary_plot_int32(bytes, trailing_offset)
        trailing_length == record_length ||
            throw(ArgumentError("binary plot record markers do not match"))
        push!(records, bytes[record_start:record_stop])
        offset = trailing_offset + 4
    end
    return records
end

function _binary_plot_bytes(values::AbstractVector{Int32})
    buffer = IOBuffer()
    for value in values
        write(buffer, value)
    end
    return take!(buffer)
end

function _binary_plot_bytes(values::AbstractVector{Float32})
    buffer = IOBuffer()
    for value in values
        write(buffer, value)
    end
    return take!(buffer)
end

function _write_binary_plot_record(io::IO, bytes::Vector{UInt8})
    length(bytes) <= typemax(Int32) ||
        throw(ArgumentError("binary plot record is too large"))
    marker = Int32(length(bytes))
    write(io, marker)
    write(io, bytes)
    write(io, marker)
    return io
end

function _binary_plot_fixed_name(name::AbstractString)
    text = _binary_plot_name(name)
    return rpad(text, 8)[1:8]
end

function _binary_plot_name_table(table::ElectromagneticBinaryPlotTable)
    return (
        names = table.name_table,
        fixed_names = table.fixed_name_table,
        node_indices = table.node_name_indices,
        upper_indices = table.branch_upper_name_indices,
        lower_indices = table.branch_lower_name_indices,
    )
end

function write_electromagnetic_binary_plot_file(
    path::AbstractString,
    table::ElectromagneticBinaryPlotTable,
)
    _report_ensure_dir(path)
    name_table = _binary_plot_name_table(table)
    node_count = count(==(:node_voltage), table.channel_classes)
    branch_voltage_count = count(==(:branch_voltage), table.channel_classes)
    branch_current_count = count(==(:branch_current), table.channel_classes)
    branch_pair_count = branch_voltage_count + branch_current_count
    header = IOBuffer()
    write(header, codeunits(table.header_record_prefix))
    write(header, Int32(length(name_table.names)))
    write(header, Int32(node_count))
    write(header, Int32(branch_current_count))
    write(header, Int32(branch_pair_count))
    for name in name_table.fixed_names
        write(header, codeunits(name))
    end
    open(path, "w") do io
        _write_binary_plot_record(io, take!(header))
        if node_count > 0
            _write_binary_plot_record(io, _binary_plot_bytes(Int32.(name_table.node_indices)))
        end
        if branch_pair_count > 0
            branch_indices = vcat(name_table.upper_indices, name_table.lower_indices)
            _write_binary_plot_record(io, _binary_plot_bytes(Int32.(branch_indices)))
        end
        for row in eachindex(table.sample_times_s)
            values = Vector{Float32}(undef, size(table.sample_values, 2) + 1)
            values[1] = Float32(table.sample_times_s[row])
            for column in axes(table.sample_values, 2)
                values[column + 1] = Float32(table.sample_values[row, column])
            end
            _write_binary_plot_record(io, _binary_plot_bytes(values))
        end
        for sentinel_row in axes(table.sentinel_sample_values, 1)
            sentinel_record = Vector{Float32}(undef, size(table.sample_values, 2) + 1)
            sentinel_record[1] = Float32(table.sentinel_value)
            for column in axes(table.sentinel_sample_values, 2)
                sentinel_record[column + 1] =
                    Float32(table.sentinel_sample_values[sentinel_row, column])
            end
            _write_binary_plot_record(io, _binary_plot_bytes(sentinel_record))
        end
    end
    return path
end

function read_electromagnetic_binary_plot_file(path::AbstractString)
    records = _binary_plot_records(path)
    length(records) >= 3 || throw(ArgumentError("binary plot file has too few records"))
    header = records[1]
    length(header) >= 48 || throw(ArgumentError("binary plot header is too short"))
    header_prefix = String(header[1:32])
    name_count = _binary_plot_int32(header, 33)
    node_count = _binary_plot_int32(header, 37)
    branch_current_count = _binary_plot_int32(header, 41)
    branch_pair_count = _binary_plot_int32(header, 45)
    branch_voltage_count = branch_pair_count - branch_current_count
    branch_voltage_count >= 0 ||
        throw(ArgumentError("binary plot branch-current count exceeds branch-pair count"))
    length(header) >= 48 + 8 * name_count ||
        throw(ArgumentError("binary plot header is too short for the name table"))
    fixed_names = [
        String(header[(49 + 8 * (index - 1)):(56 + 8 * (index - 1))])
        for index in 1:name_count
    ]
    names = [strip(name) for name in fixed_names]
    record_index = 2
    node_indices = Int[]
    if node_count > 0
        length(records) >= record_index ||
            throw(ArgumentError("binary plot node-index record is missing"))
        length(records[record_index]) == 4 * node_count ||
            throw(ArgumentError("binary plot node-index record has the wrong size"))
        node_indices = Int.(reinterpret(Int32, records[record_index]))
        record_index += 1
    end
    if branch_pair_count > 0 &&
            record_index <= length(records) &&
            isempty(records[record_index])
        record_index += 1
    end
    branch_indices = Int[]
    if branch_pair_count > 0
        length(records) >= record_index ||
            throw(ArgumentError("binary plot branch-index record is missing"))
        length(records[record_index]) == 8 * branch_pair_count ||
            throw(ArgumentError("binary plot branch-index record has the wrong size"))
        branch_indices = Int.(reinterpret(Int32, records[record_index]))
        record_index += 1
    end
    upper_branch_indices = branch_pair_count == 0 ? Int[] : branch_indices[1:branch_pair_count]
    lower_branch_indices = branch_pair_count == 0 ? Int[] : branch_indices[(branch_pair_count + 1):end]
    name_at(index::Int) = 1 <= index <= length(names) ? names[index] :
        throw(ArgumentError("binary plot name index is outside the name table"))

    channel_classes = Symbol[]
    upper_names = String[]
    lower_names = String[]
    for index in node_indices
        push!(channel_classes, :node_voltage)
        push!(upper_names, name_at(index))
        push!(lower_names, "")
    end
    for pair_index in 1:branch_pair_count
        push!(
            channel_classes,
            pair_index <= branch_voltage_count ? :branch_voltage : :branch_current,
        )
        push!(upper_names, name_at(upper_branch_indices[pair_index]))
        push!(lower_names, name_at(lower_branch_indices[pair_index]))
    end

    channel_count = length(channel_classes)
    expected_data_bytes = 4 * (channel_count + 1)
    times = Float64[]
    value_rows = Vector{Vector{Float64}}()
    sentinel_rows = Vector{Vector{Float64}}()
    sentinel_count = 0
    sentinel = -9999.0
    for record in records[record_index:end]
        length(record) == expected_data_bytes ||
            throw(ArgumentError("binary plot data record has the wrong size"))
        values = reinterpret(Float32, record)
        time_s = Float64(values[1])
        if time_s == sentinel
            sentinel_count += 1
            push!(sentinel_rows, Float64.(values[2:end]))
            continue
        end
        push!(times, time_s)
        push!(value_rows, Float64.(values[2:end]))
    end
    sample_values = zeros(Float64, length(value_rows), channel_count)
    for (row_index, row) in enumerate(value_rows)
        sample_values[row_index, :] .= row
    end
    sentinel_values = zeros(Float64, length(sentinel_rows), channel_count)
    for (row_index, row) in enumerate(sentinel_rows)
        sentinel_values[row_index, :] .= row
    end
    return electromagnetic_binary_plot_table(
        title = basename(path),
        source = path,
        channel_classes = channel_classes,
        upper_names = upper_names,
        lower_names = lower_names,
        sample_times_s = times,
        sample_values = sample_values,
        sentinel_value = sentinel,
        sentinel_sample_values = sentinel_values,
        sentinel_record_count = sentinel_count,
        header_record_prefix = header_prefix,
        name_table = names,
        fixed_name_table = fixed_names,
        node_name_indices = node_indices,
        branch_upper_name_indices = upper_branch_indices,
        branch_lower_name_indices = lower_branch_indices,
    )
end

struct ElectromagneticPlotReplaySearch
    target_time_s::Float64
    selected_sample_index::Int
    selected_time_s::Float64
    selected_values::Vector{Float64}
    clean_return::Bool
    reached_end::Bool
    deferred_effects::Tuple{Vararg{Symbol}}
end

function electromagnetic_plot_replay_search(
    table::ElectromagneticBinaryPlotTable;
    target_time_s::Real,
    deferred_effects::Tuple{Vararg{Symbol}}=(
        :interactive_command_loop,
        :terminal_plot_user_interface,
        :xpr_compatibility,
    ),
)
    target = Float64(target_time_s)
    isfinite(target) || throw(ArgumentError("target_time_s must be finite"))
    isempty(table.sample_times_s) && throw(ArgumentError("binary plot table has no samples"))
    selected = findfirst(time_s -> time_s >= target, table.sample_times_s)
    reached_end = selected === nothing
    selected_index = reached_end ? length(table.sample_times_s) : selected
    return ElectromagneticPlotReplaySearch(
        target,
        selected_index,
        table.sample_times_s[selected_index],
        vec(table.sample_values[selected_index, :]),
        true,
        reached_end,
        deferred_effects,
    )
end

struct ElectromagneticPlotReplaySession
    commands::Vector{Symbol}
    opened_plot_file::String
    target_time_s::Float64
    raw_record_index::Int
    raw_record_time_s::Float64
    selected_sample_index::Int
    selected_time_s::Float64
    selected_values::Vector{Float64}
    record_read_count::Int
    operation_prompt_count::Int
    terminal_fragments::Vector{String}
    clean_return::Bool
    expected_eof_after_replay::Bool
    reached_end::Bool
    blockers::Tuple{Vararg{Symbol}}
    deferred_effects::Tuple{Vararg{Symbol}}
    save_prompt_count::Int
    end_message_count::Int
    full_interactive_command_handling::Bool
end

const _PLOT_REPLAY_START_PROMPT =
    "EMTP BEGINS.  SEND (SPY, \$ATTACH, DEBUG, HELP, MODULE, JUNK, STOP) :"
const _PLOT_REPLAY_OPERATION_PROMPT =
    "OPERATION (OPEN, CLOSE, TOP, BOT, NEXT, BACK, TIME) :"
const _PLOT_REPLAY_FILE_PROMPT = "DESIRED DISK FILE NAME :"
const _PLOT_REPLAY_STATUS_PROMPT = "DESIRED STATUS (NEW, OLD) :"
const _PLOT_REPLAY_TIME_PROMPT = "SEND DESIRED TIME"
const _PLOT_REPLAY_RECORD_READ_PREFIX = "Ok, record of  t ="
const _PLOT_REPLAY_RECORD_READ_SUFFIX = "just read on try number"
const _PLOT_REPLAY_END_PREFIX = "Ok, at end LUNIT4."
const _PLOT_REPLAY_SAVE_PROMPT = "SAVE PERMANENTLY? (Y OR N) :"

function _plot_replay_command_lines(commands)
    if commands isa AbstractString
        return [strip(line) for line in split(String(commands), '\n') if !isempty(strip(line))]
    end
    return [strip(String(command)) for command in commands if !isempty(strip(String(command)))]
end

function _plot_replay_command_symbol(command::AbstractString)
    token = uppercase(strip(String(command)))
    startswith(token, "SPY") && return :spy
    startswith(token, "LUNIT4") && return :lunit4
    startswith(token, "OPEN") && return :open
    startswith(token, "CLOSE") && return :close
    startswith(token, "TOP") && return :top
    startswith(token, "BOT") && return :bot
    startswith(token, "NEXT") && return :next
    startswith(token, "BACK") && return :back
    startswith(token, "TIME") && return :time
    return :unsupported
end

function _plot_replay_save_response_symbol(command::AbstractString)
    token = uppercase(strip(String(command)))
    startswith(token, "Y") && return :save_yes
    startswith(token, "N") && return :save_no
    return :unsupported
end

function _plot_replay_record_first_time(record::Vector{UInt8})
    length(record) >= 4 || return NaN
    return Float64(reinterpret(Float32, record[1:4])[1])
end

_plot_replay_record_value_fragment(time_s::Real) = @sprintf("%.4E", Float64(time_s))

function _plot_replay_raw_time_search(
    records::Vector{Vector{UInt8}},
    start_index::Int,
    target_time_s::Float64,
)
    for record_index in max(start_index, 1):length(records)
        time_s = _plot_replay_record_first_time(records[record_index])
        if isfinite(time_s) && time_s >= target_time_s
            return (
                record_index = record_index,
                time_s = time_s,
                try_number = record_index - max(start_index, 1) + 1,
                next_record_index = record_index + 1,
                reached_end = false,
            )
        end
    end
    return (
        record_index = 0,
        time_s = NaN,
        try_number = max(length(records) - max(start_index, 1) + 1, 0),
        next_record_index = length(records) + 1,
        reached_end = true,
    )
end

function _plot_replay_first_data_record_index(
    records::Vector{Vector{UInt8}},
    table::ElectromagneticBinaryPlotTable,
)
    data_record_count = length(table.sample_times_s) + table.sentinel_record_count
    return max(length(records) - data_record_count + 1, 1)
end

function electromagnetic_plot_replay_session(
    path::AbstractString;
    commands,
    deferred_effects::Tuple{Vararg{Symbol}}=(
        :terminal_plot_user_interface,
        :xpr_compatibility,
    ),
)
    table = read_electromagnetic_binary_plot_file(path)
    records = _binary_plot_records(path)
    command_lines = _plot_replay_command_lines(commands)
    command_symbols = Symbol[]
    blockers = Symbol[]
    index = 1
    opened = false
    closed = false
    opened_plot_file = ""
    record_pointer = 1
    operation_prompt_count = 0
    record_read_count = 0
    raw_record_index = 0
    raw_record_time_s = NaN
    raw_reached_end = false
    target_time_s = NaN
    selected_sample_index = 0
    selected_time_s = NaN
    selected_values = Float64[]
    save_prompt_count = 0
    end_message_count = 0
    terminal_fragments = String[_PLOT_REPLAY_START_PROMPT]

    if index <= length(command_lines) &&
            _plot_replay_command_symbol(command_lines[index]) == :spy
        push!(command_symbols, :spy)
        push!(terminal_fragments, "SPY:")
        index += 1
    end
    if index <= length(command_lines) &&
            _plot_replay_command_symbol(command_lines[index]) == :lunit4
        push!(command_symbols, :lunit4)
        push!(terminal_fragments, _PLOT_REPLAY_OPERATION_PROMPT)
        operation_prompt_count += 1
        index += 1
    else
        push!(blockers, :lunit4_command_missing)
    end

    while index <= length(command_lines)
        command = _plot_replay_command_symbol(command_lines[index])
        push!(command_symbols, command)
        index += 1
        if command == :open
            if index + 1 > length(command_lines)
                push!(blockers, :open_arguments_missing)
                break
            end
            push!(terminal_fragments, _PLOT_REPLAY_FILE_PROMPT)
            opened_plot_file = command_lines[index]
            push!(terminal_fragments, _PLOT_REPLAY_STATUS_PROMPT)
            status = uppercase(strip(command_lines[index + 1]))
            status == "OLD" || push!(blockers, :open_status_not_old)
            basename(opened_plot_file) == basename(path) ||
                push!(blockers, :open_plot_file_mismatch)
            opened = true
            closed = false
            record_pointer = 1
            index += 2
            push!(terminal_fragments, _PLOT_REPLAY_OPERATION_PROMPT)
            operation_prompt_count += 1
        elseif command == :close
            opened = false
            closed = true
            push!(terminal_fragments, _PLOT_REPLAY_SAVE_PROMPT)
            save_prompt_count += 1
            if index <= length(command_lines)
                save_response = _plot_replay_save_response_symbol(command_lines[index])
                if save_response == :unsupported
                    push!(blockers, :close_save_response_invalid)
                else
                    push!(command_symbols, save_response)
                    index += 1
                end
            else
                push!(blockers, :close_save_response_missing)
            end
            push!(terminal_fragments, _PLOT_REPLAY_OPERATION_PROMPT)
            operation_prompt_count += 1
        elseif command == :top
            opened || push!(blockers, :top_without_open_plot_file)
            record_pointer = _plot_replay_first_data_record_index(records, table)
            push!(terminal_fragments, _PLOT_REPLAY_OPERATION_PROMPT)
            operation_prompt_count += 1
        elseif command == :bot
            opened || push!(blockers, :bot_without_open_plot_file)
            record_pointer = length(records)
            raw_reached_end = true
            push!(terminal_fragments, _PLOT_REPLAY_END_PREFIX)
            end_message_count += 1
            push!(terminal_fragments, _PLOT_REPLAY_OPERATION_PROMPT)
            operation_prompt_count += 1
        elseif command == :next
            opened || push!(blockers, :next_without_open_plot_file)
            raw = _plot_replay_raw_time_search(records, record_pointer, -Inf)
            record_pointer = raw.next_record_index
            raw_record_index = raw.record_index
            raw_record_time_s = raw.time_s
            raw_reached_end = raw.reached_end
            record_read_count += raw.reached_end ? 0 : 1
            raw.reached_end ||
                push!(terminal_fragments, _plot_replay_record_value_fragment(raw.time_s))
            push!(terminal_fragments, _PLOT_REPLAY_OPERATION_PROMPT)
            operation_prompt_count += 1
        elseif command == :back
            opened || push!(blockers, :back_without_open_plot_file)
            record_pointer = max(
                record_pointer - 1,
                _plot_replay_first_data_record_index(records, table),
            )
            raw = _plot_replay_raw_time_search(records, record_pointer, -Inf)
            record_pointer = raw.next_record_index
            raw_record_index = raw.record_index
            raw_record_time_s = raw.time_s
            raw_reached_end = raw.reached_end
            record_read_count += raw.reached_end ? 0 : 1
            raw.reached_end ||
                push!(terminal_fragments, _plot_replay_record_value_fragment(raw.time_s))
            push!(terminal_fragments, _PLOT_REPLAY_OPERATION_PROMPT)
            operation_prompt_count += 1
        elseif command == :time
            opened || push!(blockers, :time_without_open_plot_file)
            if index > length(command_lines)
                push!(blockers, :time_argument_missing)
                break
            end
            push!(terminal_fragments, _PLOT_REPLAY_TIME_PROMPT)
            target_time_s = try
                parse(Float64, replace(command_lines[index], 'D' => 'E', 'd' => 'e'))
            catch
                push!(blockers, :time_argument_invalid)
                NaN
            end
            index += 1
            if isfinite(target_time_s)
                raw = _plot_replay_raw_time_search(records, record_pointer, target_time_s)
                record_pointer = raw.next_record_index
                raw_record_index = raw.record_index
                raw_record_time_s = raw.time_s
                raw_reached_end = raw.reached_end
                record_read_count += raw.reached_end ? 0 : 1
                raw.reached_end || append!(
                    terminal_fragments,
                    [_PLOT_REPLAY_RECORD_READ_PREFIX, _PLOT_REPLAY_RECORD_READ_SUFFIX],
                )
                search = electromagnetic_plot_replay_search(
                    table;
                    target_time_s = target_time_s,
                    deferred_effects = (),
                )
                selected_sample_index = search.selected_sample_index
                selected_time_s = search.selected_time_s
                selected_values = search.selected_values
            end
            push!(terminal_fragments, _PLOT_REPLAY_OPERATION_PROMPT)
            operation_prompt_count += 1
        else
            push!(blockers, :unsupported_plot_replay_command)
            push!(terminal_fragments, _PLOT_REPLAY_OPERATION_PROMPT)
            operation_prompt_count += 1
        end
    end
    selected_sample_index > 0 || push!(blockers, :time_search_not_executed)
    expected_eof_after_replay = opened && !closed
    required_commands = (:open, :close, :top, :bot, :next, :back, :time)
    full_interactive_command_handling =
        all(command -> command in command_symbols, required_commands) &&
        any(command -> command in (:save_yes, :save_no), command_symbols) &&
        save_prompt_count > 0 &&
        end_message_count > 0 &&
        isempty(blockers)
    final_deferred_effects = full_interactive_command_handling ?
        Tuple(effect for effect in deferred_effects if effect != :terminal_plot_user_interface) :
        deferred_effects
    return ElectromagneticPlotReplaySession(
        command_symbols,
        opened_plot_file,
        target_time_s,
        raw_record_index,
        raw_record_time_s,
        selected_sample_index,
        selected_time_s,
        selected_values,
        record_read_count,
        operation_prompt_count,
        terminal_fragments,
        isempty(blockers),
        expected_eof_after_replay,
        raw_reached_end,
        Tuple(unique(blockers)),
        final_deferred_effects,
        save_prompt_count,
        end_message_count,
        full_interactive_command_handling,
    )
end

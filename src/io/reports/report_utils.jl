function _report_ensure_dir(path::AbstractString)
    isdir(dirname(path)) || mkpath(dirname(path))
    return path
end

_report_escape(text) = replace(
    string(text),
    "&" => "&amp;",
    "<" => "&lt;",
    ">" => "&gt;",
    "\"" => "&quot;",
)

function _report_json_escape(value)
    return replace(
        string(value),
        "\\" => "\\\\",
        "\"" => "\\\"",
        "\n" => "\\n",
    )
end

_report_json_bool(value::Bool) = value ? "true" : "false"

function _report_json_float(value::Real)
    return isfinite(value) ? @sprintf("%.12g", Float64(value)) : "null"
end

_report_json_string_array(values) =
    "[" * join(("\"" * _report_json_escape(value) * "\"" for value in values), ", ") * "]"

_report_json_int_array(values) =
    "[" * join(string.(values), ", ") * "]"

_report_json_bool_array(values) =
    "[" * join(_report_json_bool.(values), ", ") * "]"

_report_json_float_array(values) =
    "[" * join(_report_json_float.(values), ", ") * "]"

_report_json_nested_float_array(values) =
    "[" * join(_report_json_float_array.(values), ", ") * "]"

function _report_json_symbol_float_dict(values::Dict{Symbol,Float64})
    entries = [
        "\"" * _report_json_escape(key) * "\": " * _report_json_float(values[key])
        for key in sort!(collect(keys(values)); by = string)
    ]
    return "{" * join(entries, ", ") * "}"
end

function _nice_limits(values; pad::Float64 = 0.05)
    lo = minimum(values)
    hi = maximum(values)
    if lo == hi
        span = max(abs(lo), 1.0)
        return lo - pad * span, hi + pad * span
    end
    span = hi - lo
    return lo - pad * span, hi + pad * span
end

function _report_table(trace::DeckEMTTrace)
    samples = length(trace.time_s)
    table = Matrix{Float64}(
        undef,
        samples,
        1 + length(trace.node_names) + length(trace.output_channel_names),
    )
    table[:, 1] = trace.time_s
    node_count = length(trace.node_names)
    node_count > 0 && (table[:, 2:(1 + node_count)] = transpose(trace.voltage_pu))
    !isempty(trace.output_channel_names) &&
        (table[:, (2 + node_count):end] = transpose(trace.output_pu))
    return table
end

function _trace_report_headers(trace::DeckEMTTrace)
    return vcat(
        ["time_s"],
        ["$(name)_v_pu" for name in trace.node_names],
        String.(trace.output_channel_names),
    )
end

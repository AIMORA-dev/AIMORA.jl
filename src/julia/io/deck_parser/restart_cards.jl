"""A time-switch mutation accepted after a `START AGAIN` request."""
struct DeckRestartSwitchMutation
    line_no::Int
    selector_kind::Symbol
    switch_index::Union{Missing,Int}
    from_node::Union{Missing,Symbol}
    to_node::Union{Missing,Symbol}
    close_time_s::Union{Missing,Float64}
    open_time_s::Union{Missing,Float64}
    open_delay_s::Union{Missing,Float64}
    critical_current::Union{Missing,Float64}
    overwrite_zero_fields::Bool
    raw_text::String
end

"""A control-source mutation accepted after a `START AGAIN` request."""
struct DeckRestartControlSourceMutation
    line_no::Int
    name::Symbol
    source_type::Union{Missing,Int}
    amplitude::Union{Missing,Float64}
    frequency_or_delay::Union{Missing,Float64}
    phase_or_width::Union{Missing,Float64}
    activation_start_time_s::Union{Missing,Float64}
    activation_stop_time_s::Union{Missing,Float64}
    overwrite_zero_fields::Bool
    raw_text::String
end

"""
Typed continuation request for a preserved Julia EMT state.

The request consumes the fixed-column switch and control-source modification
cards used by `START AGAIN`, but the checkpoint itself remains a typed Julia
workspace rather than a binary table dump.
"""
struct DeckRestartRequest
    line_no::Int
    checkpoint_path::Union{Nothing,String}
    switch_mutations::Vector{DeckRestartSwitchMutation}
    control_source_mutations::Vector{DeckRestartControlSourceMutation}
    terminator_line_no::Int
end

function _restart_checkpoint_path(request_line::AbstractString)
    text = strip(String(request_line))
    match_result = match(r"(?i)^START\s+AGAIN(?:\s*,\s*|\s+)(.*)$", text)
    match_result === nothing && return nothing
    path = strip(only(match_result.captures))
    startswith(uppercase(path), "CHECKPOINT=") &&
        (path = strip(path[length("CHECKPOINT=") + 1:end]))
    if length(path) >= 2 &&
       ((first(path) == '"' && last(path) == '"') ||
        (first(path) == '\'' && last(path) == '\''))
        path = path[2:end - 1]
    end
    isempty(path) && throw(ArgumentError(
        "START AGAIN checkpoint path must not be empty",
    ))
    occursin(',', path) && throw(ArgumentError(
        "START AGAIN accepts one AIMORA checkpoint path; plot-file routing is separate",
    ))
    return path
end

function _restart_fixed_image(line::AbstractString)
    value = String(line)
    width = ncodeunits(value)
    return width >= 80 ? first(value, 80) : rpad(value, 80)
end

function _restart_fixed_field(image::AbstractString, first_col::Int, last_col::Int)
    return strip(String(SubString(image, first_col, last_col)))
end

function _restart_integer_field(
    image::AbstractString,
    first_col::Int,
    last_col::Int,
    line_no::Int,
    field_name::AbstractString;
    default::Int = 0,
)
    text = _restart_fixed_field(image, first_col, last_col)
    isempty(text) && return default
    value = tryparse(Int, text)
    value === nothing && throw(ArgumentError(
        "restart line $line_no has invalid $field_name '$text'",
    ))
    return value
end

function _restart_float_field(
    image::AbstractString,
    first_col::Int,
    last_col::Int,
    line_no::Int,
    field_name::AbstractString,
)
    text = _restart_fixed_field(image, first_col, last_col)
    isempty(text) && return 0.0
    value = tryparse(Float64, replace(replace(text, 'D' => 'E'), 'd' => 'e'))
    value === nothing && throw(ArgumentError(
        "restart line $line_no has invalid $field_name '$text'",
    ))
    return value
end

_restart_optional(value::Float64, overwrite::Bool) =
    overwrite || value != 0.0 ? value : missing

function _parse_restart_control_source_mutation(
    image::AbstractString,
    line::AbstractString,
    line_no::Int,
    control_code::Int,
)
    name_text = _restart_fixed_field(image, 3, 8)
    isempty(name_text) && throw(ArgumentError(
        "restart line $line_no is missing the control-source name in columns 3-8",
    ))
    overwrite = control_code < 0
    source_type = _restart_integer_field(
        image, 1, 2, line_no, "control-source type";
        default = 0,
    )
    amplitude = _restart_float_field(image, 11, 20, line_no, "source amplitude")
    frequency_or_delay =
        _restart_float_field(image, 21, 30, line_no, "source frequency or delay")
    phase_or_width =
        _restart_float_field(image, 31, 40, line_no, "source phase or width")
    start_time = _restart_float_field(image, 61, 70, line_no, "source start time")
    stop_time = _restart_float_field(image, 71, 80, line_no, "source stop time")
    return DeckRestartControlSourceMutation(
        line_no,
        Symbol(name_text),
        source_type == 0 ? missing : source_type,
        _restart_optional(amplitude, overwrite),
        _restart_optional(frequency_or_delay, overwrite),
        _restart_optional(phase_or_width, overwrite),
        _restart_optional(start_time, overwrite),
        stop_time == 0.0 ? missing : stop_time,
        overwrite,
        String(line),
    )
end

function _parse_restart_named_switch_mutation(
    image::AbstractString,
    line::AbstractString,
    line_no::Int,
    control_code::Int,
)
    switch_type = _restart_integer_field(
        image, 1, 2, line_no, "switch type";
        default = 0,
    )
    switch_type == 0 || throw(ArgumentError(
        "restart line $line_no requests unsupported topology-changing switch type $switch_type",
    ))
    from_text = _restart_fixed_field(image, 3, 8)
    to_text = _restart_fixed_field(image, 9, 14)
    isempty(from_text) && throw(ArgumentError(
        "restart line $line_no is missing the switch from-node in columns 3-8",
    ))
    isempty(to_text) && throw(ArgumentError(
        "restart line $line_no is missing the switch to-node in columns 9-14",
    ))
    overwrite = control_code < 0
    close_time = _restart_float_field(image, 15, 24, line_no, "switch close time")
    second_time = _restart_float_field(image, 25, 34, line_no, "switch open time")
    critical_current =
        _restart_float_field(image, 35, 44, line_no, "switch critical current")
    alternate_open_time =
        _restart_float_field(image, 45, 54, line_no, "switch alternate open time")
    delayed_open = alternate_open_time != 0.0
    return DeckRestartSwitchMutation(
        line_no,
        :endpoints,
        missing,
        Symbol(from_text),
        Symbol(to_text),
        _restart_optional(delayed_open && close_time < 0.0 ? 0.0 : close_time, overwrite),
        _restart_optional(delayed_open ? abs(alternate_open_time) : second_time, overwrite),
        delayed_open ? abs(second_time) : missing,
        _restart_optional(critical_current, overwrite),
        overwrite,
        String(line),
    )
end

function _parse_restart_index_switch_mutation(
    image::AbstractString,
    line::AbstractString,
    line_no::Int,
)
    switch_index = _restart_integer_field(
        image, 1, 8, line_no, "switch index";
        default = 0,
    )
    switch_index > 0 || throw(ArgumentError(
        "restart line $line_no requires a positive switch index",
    ))
    close_time = _restart_float_field(image, 9, 24, line_no, "switch close time")
    open_time = _restart_float_field(image, 25, 40, line_no, "switch open time")
    return DeckRestartSwitchMutation(
        line_no,
        :index,
        switch_index,
        missing,
        missing,
        close_time,
        open_time > 0.0 ? open_time : missing,
        missing,
        missing,
        true,
        String(line),
    )
end

"""Parse a fixed-card `START AGAIN` mutation request."""
function parse_emt_restart_request(lines; source::AbstractString = "restart request")
    records = [(line_no, String(line)) for (line_no, line) in enumerate(lines) if !isempty(strip(String(line)))]
    isempty(records) && throw(ArgumentError("$source is empty"))
    request_line_no, request_line = first(records)
    startswith(uppercase(strip(request_line)), "START AGAIN") || throw(ArgumentError(
        "$source must begin with START AGAIN",
    ))
    switch_mutations = DeckRestartSwitchMutation[]
    source_mutations = DeckRestartControlSourceMutation[]
    terminator_line_no = 0
    for (line_no, line) in records[2:end]
        image = _restart_fixed_image(line)
        control_code = _restart_integer_field(
            image, 55, 60, line_no, "restart control code";
            default = 0,
        )
        if control_code == 0
            leading_integer = tryparse(Int, _restart_fixed_field(image, 1, 8))
            if leading_integer == 9999
                terminator_line_no = line_no
                break
            end
            push!(
                switch_mutations,
                _parse_restart_index_switch_mutation(image, line, line_no),
            )
        elseif abs(control_code) == 1111
            push!(
                source_mutations,
                _parse_restart_control_source_mutation(
                    image,
                    line,
                    line_no,
                    control_code,
                ),
            )
        else
            push!(
                switch_mutations,
                _parse_restart_named_switch_mutation(
                    image,
                    line,
                    line_no,
                    control_code,
                ),
            )
        end
    end
    terminator_line_no > 0 || throw(ArgumentError(
        "$source is missing the 9999 restart terminator",
    ))
    return DeckRestartRequest(
        request_line_no,
        _restart_checkpoint_path(request_line),
        switch_mutations,
        source_mutations,
        terminator_line_no,
    )
end

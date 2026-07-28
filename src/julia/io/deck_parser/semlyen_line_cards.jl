mutable struct SemlyenModeCardBuilder
    header_line_no::Int
    detail_line_numbers::Vector{Int}
    from_node::Symbol
    to_node::Symbol
    from_node_index::Int
    to_node_index::Int
    mode_index::Int
    characteristic_admittance_s::Float64
    travel_time_s::Float64
    propagation_triple_count::Int
    admittance_triple_count::Int
    phasor_values::Vector{Float64}
    propagation_values::Vector{Float64}
    admittance_values::Vector{Float64}
end

mutable struct SemlyenLineCardParseState
    phase_count::Int
    imaginary_transform_parts_supplied::Bool
    modes::Vector{DeckSemlyenModeRow}
    current_mode::Union{Nothing,SemlyenModeCardBuilder}
    stage::Symbol
    voltage_transform_values::Vector{Float64}
    current_transform_values::Vector{Float64}
    transform_line_numbers::Vector{Int}
end

function semlyen_transform_imaginary_parts_marker(tokens)
    values = lowercase.(deck_token_value.(tokens))
    length(values) >= 3 || return nothing
    joined = join(values, " ")
    occursin("no imaginary part", joined) && return false
    occursin("yes imaginary part", joined) && return true
    return nothing
end

function _semlyen_field(image::AbstractString, first::Int, last::Int)
    padded = rpad(String(image), last)
    return strip(padded[first:last])
end

function _semlyen_header_integer(image::AbstractString, first::Int, last::Int)
    value = _semlyen_field(image, first, last)
    isempty(value) && return 0
    return tryparse(Int, value)
end

function _semlyen_header_float(image::AbstractString, first::Int, last::Int)
    value = replace(_semlyen_field(image, first, last), 'D' => 'E', 'd' => 'e')
    isempty(value) && return nothing
    return tryparse(Float64, value)
end

function semlyen_line_header_card(line::AbstractString)
    image = rpad(String(line), 78)
    _semlyen_field(image, 1, 2) == "-1" || return false
    phase_count = _semlyen_header_integer(image, 76, 78)
    mode_row = _semlyen_header_integer(image, 63, 65)
    mode_column = _semlyen_header_integer(image, 66, 68)
    propagation_count = _semlyen_header_integer(image, 69, 71)
    admittance_count = _semlyen_header_integer(image, 72, 74)
    return phase_count !== nothing && phase_count > 0 &&
        mode_row !== nothing && mode_row > 0 &&
        mode_column !== nothing && mode_column > 0 &&
        propagation_count !== nothing && propagation_count >= 0 &&
        admittance_count !== nothing && admittance_count >= 0 &&
        _semlyen_header_float(image, 27, 38) !== nothing &&
        _semlyen_header_float(image, 39, 50) !== nothing
end

function _semlyen_node!(result::DeckParseResult, raw::AbstractString)
    value = strip(String(raw))
    isempty(value) && return (:ground, 0)
    name = Symbol(value)
    return name, node_id!(result, value)
end

function _semlyen_mode_builder!(result::DeckParseResult, line, line_no::Int, phase_count::Int)
    image = rpad(String(line), 78)
    semlyen_line_header_card(image) || begin
        add_issue!(result.validation, invalid_value("line $line_no", "expected a Semlyen mode header"))
        return nothing
    end
    declared_phase_count = something(_semlyen_header_integer(image, 76, 78), 0)
    declared_phase_count == phase_count || begin
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "Semlyen mode phase count $declared_phase_count does not match group count $phase_count",
            ),
        )
        return nothing
    end
    from_node, from_index = _semlyen_node!(result, _semlyen_field(image, 3, 8))
    to_node, to_index = _semlyen_node!(result, _semlyen_field(image, 9, 14))
    mode_row = something(_semlyen_header_integer(image, 63, 65), 0)
    mode_column = something(_semlyen_header_integer(image, 66, 68), 0)
    mode_row == mode_column || begin
        add_issue!(result.validation, invalid_value("line $line_no", "Semlyen mode row and column must match"))
        return nothing
    end
    return SemlyenModeCardBuilder(
        line_no,
        Int[],
        from_node,
        to_node,
        from_index,
        to_index,
        mode_row,
        something(_semlyen_header_float(image, 27, 38), NaN),
        something(_semlyen_header_float(image, 39, 50), NaN),
        something(_semlyen_header_integer(image, 69, 71), 0),
        something(_semlyen_header_integer(image, 72, 74), 0),
        Float64[],
        Float64[],
        Float64[],
    )
end

function _semlyen_numeric_values(line::AbstractString)
    data = first(split(first(split(String(line), '{')), '!'))
    values = Float64[]
    for matched in eachmatch(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?", data)
        push!(values, parse(Float64, replace(matched.match, 'D' => 'E', 'd' => 'e')))
    end
    return values
end

function start_semlyen_line_parse!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    imaginary_transform_parts_supplied::Bool,
)
    image = rpad(String(line), 78)
    phase_count = something(_semlyen_header_integer(image, 76, 78), 0)
    builder = _semlyen_mode_builder!(result, image, line_no, phase_count)
    builder === nothing && return nothing
    record_card!(result, :semlyen_frequency_dependent_line_header)
    return SemlyenLineCardParseState(
        phase_count,
        imaginary_transform_parts_supplied,
        DeckSemlyenModeRow[],
        builder,
        :phasor,
        Float64[],
        Float64[],
        Int[],
    )
end

function _semlyen_triples(values::Vector{Float64})
    length(values) % 3 == 0 || throw(ArgumentError("Semlyen exponential data must contain triples"))
    return [(values[index], values[index + 1], values[index + 2]) for index in 1:3:length(values)]
end

function _finish_semlyen_mode!(result::DeckParseResult, state::SemlyenLineCardParseState)
    builder = something(state.current_mode)
    length(builder.phasor_values) == 5 ||
        throw(ArgumentError("Semlyen phasor row must contain five values"))
    propagation_terms = semlyen_rational_terms(_semlyen_triples(builder.propagation_values))
    admittance_terms = semlyen_rational_terms(_semlyen_triples(builder.admittance_values))
    phasor = builder.phasor_values
    parameters = SemlyenModeParameters(
        builder.characteristic_admittance_s,
        builder.travel_time_s,
        complex(phasor[1], phasor[2]),
        complex(phasor[3], phasor[4]),
        phasor[5],
        propagation_terms,
        admittance_terms,
    )
    push!(
        state.modes,
        DeckSemlyenModeRow(
            builder.header_line_no,
            copy(builder.detail_line_numbers),
            builder.from_node,
            builder.to_node,
            builder.from_node_index,
            builder.to_node_index,
            builder.mode_index,
            parameters,
        ),
    )
    record_card!(result, :semlyen_frequency_dependent_line_mode)
    state.current_mode = nothing
    state.stage = length(state.modes) == state.phase_count ? :voltage_transform : :mode_header
    return state
end

function _semlyen_transform_matrix(values::Vector{Float64}, phase_count::Int, complex_parts::Bool)
    entry_count = phase_count * phase_count
    if complex_parts
        length(values) == 2 * entry_count ||
            throw(ArgumentError("Semlyen complex transform has the wrong value count"))
        transform_values = complex.(values[1:2:end], values[2:2:end])
    else
        length(values) == entry_count ||
            throw(ArgumentError("Semlyen real transform has the wrong value count"))
        transform_values = complex.(values)
    end
    return Matrix(transpose(reshape(transform_values, phase_count, phase_count)))
end

function _finish_semlyen_group!(result::DeckParseResult, state::SemlyenLineCardParseState)
    voltage_transform = _semlyen_transform_matrix(
        state.voltage_transform_values,
        state.phase_count,
        state.imaginary_transform_parts_supplied,
    )
    current_transform = _semlyen_transform_matrix(
        state.current_transform_values,
        state.phase_count,
        state.imaginary_transform_parts_supplied,
    )
    semlyen_line_physical_checks(
        getfield.(state.modes, :parameters),
        voltage_transform,
        current_transform,
        inactive_phase_indices = [
            index for index in eachindex(state.modes)
            if state.modes[index].from_node_index == 0 && state.modes[index].to_node_index == 0
        ],
    )
    indices = sort(getfield.(state.modes, :mode_index))
    indices == collect(1:state.phase_count) ||
        throw(ArgumentError("Semlyen mode indices must cover 1:$(state.phase_count)"))
    ordered_modes = sort(state.modes; by = row -> row.mode_index)
    push!(
        result.semlyen_line_groups,
        DeckSemlyenLineGroupRow(
            Symbol("semlyen_line_", length(result.semlyen_line_groups) + 1),
            state.phase_count,
            ordered_modes,
            voltage_transform,
            current_transform,
            copy(state.transform_line_numbers),
            state.imaginary_transform_parts_supplied,
        ),
    )
    record_card!(result, :semlyen_frequency_dependent_line_group)
    return nothing
end

function parse_semlyen_line_card!(
    result::DeckParseResult,
    state::SemlyenLineCardParseState,
    line::AbstractString,
    line_no::Int,
)
    try
        marker = semlyen_transform_imaginary_parts_marker(deck_tokens(line))
        if marker !== nothing
            state.stage in (:voltage_transform, :current_transform) ||
                throw(ArgumentError("Semlyen transform representation marker appeared inside mode data"))
            state.imaginary_transform_parts_supplied = marker
            record_card!(
                result,
                marker ?
                    :semlyen_complex_modal_transform_enabled :
                    :semlyen_real_modal_transform_enabled,
            )
            return state
        end
        if state.stage == :mode_header
            state.current_mode = _semlyen_mode_builder!(result, line, line_no, state.phase_count)
            state.current_mode === nothing && return nothing
            state.stage = :phasor
            record_card!(result, :semlyen_frequency_dependent_line_header)
            return state
        end
        values = _semlyen_numeric_values(line)
        isempty(values) && throw(ArgumentError("expected numeric Semlyen continuation data"))
        if state.stage == :phasor
            length(values) >= 5 || throw(ArgumentError("Semlyen phasor row must contain five values"))
            append!(something(state.current_mode).phasor_values, values[1:5])
            push!(something(state.current_mode).detail_line_numbers, line_no)
            builder = something(state.current_mode)
            state.stage = builder.propagation_triple_count > 0 ? :propagation :
                (builder.admittance_triple_count > 0 ? :admittance : :finish_mode)
        elseif state.stage == :propagation
            builder = something(state.current_mode)
            append!(builder.propagation_values, values)
            push!(builder.detail_line_numbers, line_no)
            length(builder.propagation_values) <= 3 * builder.propagation_triple_count ||
                throw(ArgumentError("too many Semlyen propagation values"))
            length(builder.propagation_values) == 3 * builder.propagation_triple_count &&
                (state.stage = builder.admittance_triple_count > 0 ? :admittance : :finish_mode)
        elseif state.stage == :admittance
            builder = something(state.current_mode)
            append!(builder.admittance_values, values)
            push!(builder.detail_line_numbers, line_no)
            length(builder.admittance_values) <= 3 * builder.admittance_triple_count ||
                throw(ArgumentError("too many Semlyen admittance values"))
            length(builder.admittance_values) == 3 * builder.admittance_triple_count &&
                (state.stage = :finish_mode)
        elseif state.stage == :voltage_transform
            append!(state.voltage_transform_values, values)
            push!(state.transform_line_numbers, line_no)
            expected = state.phase_count^2 * (state.imaginary_transform_parts_supplied ? 2 : 1)
            length(state.voltage_transform_values) <= expected ||
                throw(ArgumentError("too many Semlyen voltage-transform values"))
            length(state.voltage_transform_values) == expected && (state.stage = :current_transform)
        elseif state.stage == :current_transform
            append!(state.current_transform_values, values)
            push!(state.transform_line_numbers, line_no)
            expected = state.phase_count^2 * (state.imaginary_transform_parts_supplied ? 2 : 1)
            length(state.current_transform_values) <= expected ||
                throw(ArgumentError("too many Semlyen current-transform values"))
            length(state.current_transform_values) == expected && return _finish_semlyen_group!(result, state)
        end
        state.stage == :finish_mode && _finish_semlyen_mode!(result, state)
        return state
    catch err
        add_issue!(result.validation, invalid_value("line $line_no", sprint(showerror, err)))
        return nothing
    end
end

function deck_semlyen_line_elements(result::DeckParseResult, timestep_s::Real)
    return SemlyenFrequencyDependentLine[
        semlyen_frequency_dependent_line(
            getfield.(group.modes, :from_node_index),
            getfield.(group.modes, :to_node_index),
            getfield.(group.modes, :parameters),
            group.voltage_modal_to_phase,
            group.current_modal_to_phase,
            timestep_s,
        )
        for group in result.semlyen_line_groups
    ]
end

deck_semlyen_line_element_names(result::DeckParseResult) =
    getfield.(result.semlyen_line_groups, :name)

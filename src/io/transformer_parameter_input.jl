module TransformerParameterInput

using ..TransformerParameters:
    AbstractTransformerParameterCase,
    MultiphaseTransformerParameterCase,
    MultiphaseTransformerWinding,
    SaturableTransformerParameterCase,
    TransformerSequenceShortCircuitTest,
    TransformerParameterWinding,
    TransformerShortCircuitCase

export OVER41TransformerParameterDeck,
       parse_over41_transformer_parameter_file,
       parse_over41_transformer_parameter_lines

struct OVER41TransformerParameterDeck
    source::String
    cases::Vector{AbstractTransformerParameterCase}
end

function _fixed_image(line::AbstractString)
    text = String(line)
    ncodeunits(text) <= 80 || (text = first(text, 80))
    return rpad(text, 80)
end

_field(image::AbstractString, first_column::Int, last_column::Int) =
    strip(String(SubString(image, first_column, last_column)))

function _float_field(
    image::AbstractString,
    first_column::Int,
    last_column::Int,
)
    raw = _field(image, first_column, last_column)
    isempty(raw) && return nothing
    return tryparse(Float64, replace(raw, 'D' => 'E', 'd' => 'e'))
end

function _int_field(
    image::AbstractString,
    first_column::Int,
    last_column::Int,
)
    raw = _field(image, first_column, last_column)
    isempty(raw) && return nothing
    value = tryparse(Int, raw)
    value !== nothing && return value
    numeric = tryparse(Float64, replace(raw, 'D' => 'E', 'd' => 'e'))
    return numeric === nothing ? nothing : trunc(Int, numeric)
end

function _numeric_tokens(line::AbstractString)
    return Float64[
        value
        for token in split(strip(String(line)))
        for value in (tryparse(Float64, replace(token, 'D' => 'E', 'd' => 'e')),)
        if value !== nothing
    ]
end

function _required_float(
    image::AbstractString,
    first_column::Int,
    last_column::Int,
    label::AbstractString,
    line_number::Int,
)
    value = _float_field(image, first_column, last_column)
    value === nothing &&
        throw(ArgumentError("line $line_number: expected $label in columns $first_column-$last_column"))
    return value
end

function _xformer_header(line::AbstractString, line_number::Int)
    image = _fixed_image(line)
    startswith(uppercase(strip(String(line))), "XFORMER") ||
        throw(ArgumentError("line $line_number: expected XFORMER request"))
    frequency = _float_field(image, 33, 40)
    inductance_frequency = _float_field(image, 41, 48)
    if frequency === nothing
        values = _numeric_tokens(line)
        frequency = isempty(values) ? 0.0 : values[1]
        inductance_frequency =
            length(values) >= 2 ? values[2] : 0.0
    end
    return Float64(frequency), something(inductance_frequency, 0.0)
end

function _required_int(
    image::AbstractString,
    first_column::Int,
    last_column::Int,
    label::AbstractString,
    line_number::Int,
)
    value = _int_field(image, first_column, last_column)
    value === nothing &&
        throw(ArgumentError("line $line_number: expected $label in columns $first_column-$last_column"))
    return value
end

function _node_pairs(line::AbstractString)
    image = _fixed_image(line)
    uppercase(_field(image, 1, 6)) == "BRANCH" || return nothing
    names = Symbol[]
    for first_column in 9:6:80
        name = _field(image, first_column, min(first_column + 5, 80))
        isempty(name) || push!(names, Symbol(name))
    end
    length(names) >= 2 || return NTuple{2,Symbol}[]
    return NTuple{2,Symbol}[
        (names[index], names[index + 1])
        for index in 1:2:(length(names) - 1)
    ]
end

function _short_circuit_case(
    lines,
    start_index::Int,
    source_line::Int,
    inherited_pairs::Vector{NTuple{2,Symbol}},
)
    config_image = _fixed_image(lines[start_index])
    winding_count = _int_field(config_image, 1, 1)
    if winding_count === nothing
        values = _numeric_tokens(lines[start_index])
        isempty(values) && return nothing, start_index + 1
        winding_count = trunc(Int, values[1])
    end
    winding_count in (2, 3) || return nothing, start_index
    magnetizing_current =
        something(_float_field(config_image, 2, 10), 0.0)
    magnetizing_rating =
        something(_float_field(config_image, 11, 20), 0.0)
    if magnetizing_rating <= 0.0
        values = _numeric_tokens(lines[start_index])
        length(values) >= 3 ||
            throw(ArgumentError("line $source_line: XFORMER case requires winding count, magnetizing current, and rating"))
        magnetizing_current = values[2]
        magnetizing_rating = values[3]
    end

    row_count = winding_count == 2 ? 1 : 3
    last_index = start_index + row_count
    last_index <= length(lines) ||
        throw(ArgumentError("line $source_line: incomplete XFORMER short-circuit rows"))
    voltages = Float64[]
    losses = Float64[]
    impedances = Float64[]
    ratings = Float64[]
    for index in 1:row_count
        line_number = start_index + index
        image = _fixed_image(lines[line_number])
        if winding_count == 2
            push!(voltages, _required_float(image, 1, 10, "primary_voltage_kv", line_number))
            push!(voltages, _required_float(image, 11, 20, "secondary_voltage_kv", line_number))
            push!(losses, _required_float(image, 21, 30, "load_loss_kw", line_number))
            push!(impedances, _required_float(image, 31, 40, "impedance_percent", line_number))
            push!(ratings, _required_float(image, 41, 50, "rating_mva", line_number))
        else
            push!(voltages, _required_float(image, 1, 10, "winding_voltage_kv", line_number))
            push!(losses, _required_float(image, 11, 20, "load_loss_kw", line_number))
            push!(impedances, _required_float(image, 21, 30, "impedance_percent", line_number))
            push!(ratings, _required_float(image, 31, 40, "rating_mva", line_number))
        end
    end
    pairs = isempty(inherited_pairs) ?
        fill((Symbol(""), Symbol("")), winding_count) :
        inherited_pairs[1:min(length(inherited_pairs), winding_count)]
    case = TransformerShortCircuitCase(
        source_line,
        voltages,
        losses,
        impedances,
        ratings,
        magnetizing_current,
        magnetizing_rating,
        pairs,
    )
    return case, last_index + 1
end

function _saturable_header(
    line::AbstractString,
    line_number::Int,
)
    image = _fixed_image(line)
    uppercase(_field(image, 3, 8)) == "TRANSF" ||
        throw(ArgumentError("line $line_number: expected saturable TRANSF card"))
    current = _required_float(image, 27, 32, "magnetizing_current", line_number)
    flux = _required_float(image, 33, 38, "magnetizing_flux", line_number)
    resistance = something(_float_field(image, 45, 50), 0.0)
    return current, flux, resistance <= 0.0 ? Inf : resistance
end

function _saturable_case(
    lines,
    start_index::Int,
    frequency_hz::Float64,
    inductance_frequency_hz::Float64,
)
    current, flux, magnetizing_resistance =
        _saturable_header(lines[start_index], start_index)
    windings = TransformerParameterWinding[]
    index = start_index + 1
    while index <= length(lines)
        isempty(strip(String(lines[index]))) && break
        image = _fixed_image(lines[index])
        winding_number = _int_field(image, 1, 2)
        winding_number == length(windings) + 1 || break
        from_node = _field(image, 3, 8)
        to_node = _field(image, 9, 14)
        isempty(from_node) &&
            throw(ArgumentError("line $index: transformer winding requires a from-node"))
        isempty(to_node) && (to_node = "0")
        push!(
            windings,
            TransformerParameterWinding(
                winding_number,
                from_node,
                to_node,
                _required_float(image, 27, 32, "winding_resistance_ohm", index),
                _required_float(image, 33, 38, "winding_inductance_or_reactance", index),
                _required_float(image, 39, 44, "winding_voltage_kv", index),
            ),
        )
        index += 1
    end
    case = SaturableTransformerParameterCase(
        start_index,
        frequency_hz,
        inductance_frequency_hz,
        current,
        flux,
        magnetizing_resistance,
        windings,
    )
    return case, index
end

function _bctran_header(
    line::AbstractString,
    line_number::Int,
)
    image = _fixed_image(line)
    winding_count = _int_field(image, 1, 2)
    winding_count === nothing && return nothing
    2 <= winding_count <= 10 || return nothing
    excitation = ntuple(Val(7)) do field
        first_column = 3 + (field - 1) * 10
        return something(
            _float_field(image, first_column, first_column + 9),
            0.0,
        )
    end
    flags = ntuple(Val(4)) do field
        first_column = 73 + (field - 1) * 2
        return something(_int_field(image, first_column, first_column + 1), 0)
    end
    return (
        winding_count = winding_count,
        frequency_hz = excitation[1],
        positive_exciting_current_percent = excitation[2],
        positive_excitation_rating_mva = excitation[3],
        positive_excitation_loss_kw = excitation[4],
        zero_exciting_current_percent = excitation[5],
        zero_excitation_rating_mva = excitation[6],
        zero_excitation_loss_kw = excitation[7],
        phase_count = flags[1] == 1 ? 1 : 3,
        excitation_test_winding = flags[2],
        magnetizing_output_winding = flags[3],
        output_representation =
            flags[4] > 0 ? :reactance : :inverse_inductance,
    )
end

function _bctran_winding(
    line::AbstractString,
    line_number::Int,
)
    image = _fixed_image(line)
    winding_number =
        _required_int(image, 1, 3, "winding_number", line_number)
    voltage =
        _required_float(image, 4, 13, "winding_voltage_kv", line_number)
    resistance =
        _required_float(image, 14, 23, "winding_resistance_ohm", line_number)
    names = ntuple(Val(6)) do field
        first_column = 25 + (field - 1) * 6
        return Symbol(_field(image, first_column, first_column + 5))
    end
    node_pairs = (
        (names[1], names[2]),
        (names[3], names[4]),
        (names[5], names[6]),
    )
    return MultiphaseTransformerWinding(
        winding_number,
        voltage,
        resistance,
        node_pairs,
    )
end

function _bctran_test(
    line::AbstractString,
    line_number::Int,
    phase_count::Int,
)
    image = _fixed_image(line)
    from_winding =
        _required_int(image, 1, 2, "from_winding", line_number)
    to_winding =
        _required_int(image, 3, 4, "to_winding", line_number)
    values = ntuple(Val(5)) do field
        first_column = 5 + (field - 1) * 10
        return something(
            _float_field(image, first_column, first_column + 9),
            0.0,
        )
    end
    zero_impedance =
        phase_count == 1 ? values[2] : values[4]
    zero_rating =
        phase_count == 1 ? values[3] : values[5]
    closed_delta =
        something(_int_field(image, 55, 56), 0)
    calculate_resistance =
        something(_int_field(image, 57, 58), 0) != 0
    return (
        test = TransformerSequenceShortCircuitTest(
            from_winding,
            to_winding,
            values[1],
            values[2],
            values[3],
            zero_impedance,
            zero_rating,
            closed_delta,
        ),
        calculate_resistance = calculate_resistance,
    )
end

function _next_bctran_data_index(
    lines,
    start_index::Int,
)
    index = start_index
    while index <= length(lines)
        stripped = strip(String(lines[index]))
        if isempty(stripped) ||
           startswith(stripped, "C") ||
           startswith(stripped, "\$")
            index += 1
            continue
        end
        return index
    end
    return nothing
end

function _bctran_case(
    lines,
    start_index::Int,
)
    header = _bctran_header(lines[start_index], start_index)
    header === nothing &&
        throw(ArgumentError("line $start_index: expected a complete BCTRAN case header"))
    windings = MultiphaseTransformerWinding[]
    index = start_index + 1
    for _ in 1:header.winding_count
        data_index = _next_bctran_data_index(lines, index)
        data_index === nothing &&
            throw(ArgumentError("line $start_index: incomplete BCTRAN winding block"))
        push!(windings, _bctran_winding(lines[data_index], data_index))
        index = data_index + 1
    end
    test_count =
        header.winding_count * (header.winding_count - 1) ÷ 2
    tests = TransformerSequenceShortCircuitTest[]
    calculate_resistance = false
    for test_number in 1:test_count
        data_index = _next_bctran_data_index(lines, index)
        data_index === nothing &&
            throw(ArgumentError("line $start_index: incomplete BCTRAN short-circuit block"))
        parsed = _bctran_test(
            lines[data_index],
            data_index,
            header.phase_count,
        )
        push!(tests, parsed.test)
        test_number == 1 &&
            (calculate_resistance = parsed.calculate_resistance)
        index = data_index + 1
    end
    parameter_case = MultiphaseTransformerParameterCase(
        start_index,
        header.frequency_hz,
        header.positive_exciting_current_percent,
        header.positive_excitation_rating_mva,
        header.positive_excitation_loss_kw,
        header.zero_exciting_current_percent,
        header.zero_excitation_rating_mva,
        header.zero_excitation_loss_kw,
        header.phase_count,
        header.excitation_test_winding,
        header.magnetizing_output_winding,
        header.output_representation,
        calculate_resistance,
        windings,
        tests,
    )
    return parameter_case, index
end

function _bctran_cases(
    lines,
    start_index::Int,
)
    cases = MultiphaseTransformerParameterCase[]
    index = start_index
    while index <= length(lines)
        stripped = strip(String(lines[index]))
        uppercase_line = uppercase(stripped)
        if isempty(stripped) ||
           startswith(stripped, "C") ||
           startswith(stripped, "\$") ||
           startswith(uppercase_line, "BLANK")
            index += 1
            continue
        end
        startswith(uppercase_line, "BEGIN NEW DATA CASE") && break
        _bctran_header(lines[index], index) === nothing && break
        parameter_case, index = _bctran_case(lines, index)
        push!(cases, parameter_case)
    end
    isempty(cases) &&
        throw(ArgumentError("BCTRAN request contains no complete transformer cases"))
    return cases, index
end

function _ordinary_xformer_cases(
    lines,
    start_index::Int,
    frequency_hz::Float64,
    inductance_frequency_hz::Float64,
)
    cases = AbstractTransformerParameterCase[]
    index = start_index
    node_pairs = NTuple{2,Symbol}[]
    while index <= length(lines)
        stripped = strip(String(lines[index]))
        uppercase_line = uppercase(stripped)
        startswith(uppercase_line, "BEGIN NEW DATA CASE") && break
        startswith(uppercase_line, "BLANK") && break
        if isempty(stripped)
            isempty(cases) || break
            index += 1
            continue
        end
        if startswith(stripped, "C") ||
           startswith(stripped, "/") ||
           startswith(stripped, "\$")
            index += 1
            continue
        end
        if frequency_hz > 0.0
            image = _fixed_image(lines[index])
            uppercase(_field(image, 3, 8)) == "TRANSF" || break
            parameter_case, index = _saturable_case(
                lines,
                index,
                frequency_hz,
                inductance_frequency_hz,
            )
            push!(cases, parameter_case)
            continue
        end
        parsed_pairs = _node_pairs(lines[index])
        if parsed_pairs !== nothing
            node_pairs = parsed_pairs
            index += 1
            continue
        end
        parameter_case, next_index = _short_circuit_case(
            lines,
            index,
            index,
            node_pairs,
        )
        parameter_case === nothing && break
        push!(cases, parameter_case)
        index = next_index
    end
    isempty(cases) &&
        throw(ArgumentError("XFORMER request contains no complete transformer cases"))
    return cases, index
end

function parse_over41_transformer_parameter_lines(
    lines;
    source::AbstractString = "deck",
)
    records = String.(lines)
    cases = AbstractTransformerParameterCase[]
    index = 1
    while index <= length(records)
        line = records[index]
        stripped = strip(line)
        uppercase_line = uppercase(stripped)
        if isempty(stripped) || startswith(stripped, "/") || startswith(stripped, "C ")
            index += 1
            continue
        end
        if startswith(uppercase_line, "ACCESS MODULE BCTRAN")
            parsed_cases, index = _bctran_cases(records, index + 1)
            append!(cases, parsed_cases)
            continue
        end
        startswith(uppercase_line, "XFORMER") || begin
            index += 1
            continue
        end
        frequency, inductance_frequency = _xformer_header(line, index)
        index += 1
        if frequency == 44.0
            parsed_cases, index = _bctran_cases(records, index)
            append!(cases, parsed_cases)
        else
            parsed_cases, index = _ordinary_xformer_cases(
                records,
                index,
                frequency,
                inductance_frequency,
            )
            append!(cases, parsed_cases)
        end
    end
    isempty(cases) &&
        throw(ArgumentError("deck contains no accepted OVER41 XFORMER parameter cases"))
    return OVER41TransformerParameterDeck(String(source), cases)
end

parse_over41_transformer_parameter_file(path::AbstractString) =
    parse_over41_transformer_parameter_lines(
        readlines(path);
        source = String(path),
    )

end

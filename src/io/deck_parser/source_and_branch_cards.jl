
function parse_arrester_constant_row!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    nonlinear_row_index::Int,
    first_constant_index::Int,
)::Int
    initial_issues = length(result.validation.issues)
    values = arrester_constant_values(line)
    remaining = 19 - first_constant_index
    if values === nothing || isempty(values) || length(values) > remaining
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "expected arrester constants within the 18-value A-constant table",
            ),
        )
        return 0
    end
    push!(
        result.arrester_constant_rows,
        DeckArresterConstantRow(
            nonlinear_row_index,
            first_constant_index,
            line_no,
            Float64.(values),
            String(line),
        ),
    )
    record_fixed_card!(result, :arrester_nonlinear, :arrester_constants, initial_issues)
    return length(values)
end

function fixed_card_coupled_line_descriptor(line_type::Int)
    if 51 <= line_type <= 90
        return (kind = :mutual_source_equivalent, phase_index = line_type - 50)
    elseif -3 <= line_type <= -1
        return (kind = :distributed_transmission_line, phase_index = abs(line_type))
    end
    return nothing
end

function fixed_card_coupled_lumped_pair_values!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
)
    resistance_fields = Union{Missing,Float64}[]
    inductance_fields = Union{Missing,Float64}[]
    for slot in 1:3
        first_column = 27 + 18 * (slot - 1)
        resistance = fixed_card_optional_float_field!(
            result,
            image,
            line_no,
            first_column,
            first_column + 5,
            "coupled_lumped_resistance_$slot",
        )
        inductance = fixed_card_optional_float_field!(
            result,
            image,
            line_no,
            first_column + 6,
            first_column + 17,
            "coupled_lumped_inductance_$slot",
        )
        push!(
            resistance_fields,
            ismissing(resistance) ? missing : Float64(resistance),
        )
        push!(
            inductance_fields,
            ismissing(inductance) ? missing : Float64(inductance),
        )
    end
    last_value_slot = findlast(eachindex(resistance_fields)) do slot
        !ismissing(resistance_fields[slot]) ||
            !ismissing(inductance_fields[slot])
    end
    last_value_slot === nothing &&
        return (resistance_values = Float64[], inductance_values = Float64[])
    resistance_values = Float64[
        coalesce(resistance_fields[slot], 0.0)
        for slot in 1:last_value_slot
    ]
    inductance_values = Float64[
        coalesce(inductance_fields[slot], 0.0)
        for slot in 1:last_value_slot
    ]
    return (; resistance_values, inductance_values)
end

function coupled_lumped_numeric_continuation_expected(result::DeckParseResult)
    isempty(result.coupled_line_rows) && return false
    row = last(result.coupled_line_rows)
    return row.line_kind == :mutual_source_equivalent &&
           ismissing(row.reference_from_node_value) &&
           ismissing(row.reference_to_node_value) &&
           length(row.triangular_resistance_values) < row.phase_index
end

function coupled_lumped_numeric_continuation_row(
    result::DeckParseResult,
    image::AbstractString,
)
    coupled_lumped_numeric_continuation_expected(result) || return false
    isempty(fixed_field(image, 1, 26)) || return false
    return any(
        !isempty(fixed_field(image, first_column, first_column + 17))
        for first_column in (27, 45, 63)
    )
end

function append_coupled_lumped_continuation_row!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    initial_issues::Int,
)
    coupled_lumped_numeric_continuation_expected(result) || begin
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "coupled R-L continuation row has no incomplete owner row",
            ),
        )
        return true
    end
    row = last(result.coupled_line_rows)
    values = fixed_card_coupled_lumped_pair_values!(result, image, line_no)
    isempty(values.resistance_values) && begin
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "coupled R-L continuation row requires numeric resistance/inductance pairs",
            ),
        )
        return true
    end
    remaining = row.phase_index - length(row.triangular_resistance_values)
    length(values.resistance_values) <= remaining || begin
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "coupled R-L continuation supplies more than $remaining remaining triangular pairs",
            ),
        )
        return true
    end
    append!(row.triangular_resistance_values, values.resistance_values)
    append!(row.triangular_inductance_values, values.inductance_values)
    record_fixed_card!(
        result,
        :bpa_fixed_branch,
        :fixed_card_coupled_lumped_continuation,
        initial_issues,
    )
    return true
end

function cascaded_pi_header_card(line::AbstractString)::Bool
    image = fixed_image(line)
    return uppercase(fixed_field(image, 3, 8)) == "CASCAD" &&
           uppercase(fixed_field(image, 9, 14)) == "ED PI"
end

function parse_cascaded_pi_header!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
)
    image = fixed_image(line)
    phase_count = fixed_int_field!(
        result,
        image,
        line_no,
        27,
        32,
        "cascaded_pi_phase_count",
    )
    frequency_hz = fixed_float_field!(
        result,
        image,
        line_no,
        33,
        38,
        "cascaded_pi_frequency_hz",
    )
    if phase_count !== nothing && !(1 <= phase_count <= 30)
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "cascaded PI phase count must be between 1 and 30",
        ))
    end
    if frequency_hz !== nothing && !(isfinite(frequency_hz) && frequency_hz > 0.0)
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "cascaded PI frequency must be finite and positive",
        ))
    end
    record_card!(result, :fixed_field)
    record_card!(result, :cascaded_pi_header)
    return (
        header_line_no = line_no,
        phase_count = something(phase_count, 0),
        frequency_hz = something(frequency_hz, 0.0),
        first_section_row_index = length(result.coupled_phase_pi_section_rows) + 1,
        configuration = nothing,
        blocks = DeckCascadedPiBlock[],
    )
end

function cascaded_pi_source_section_complete(result::DeckParseResult, pending)
    first_index = pending.first_section_row_index
    last_index = first_index + pending.phase_count - 1
    last_index <= length(result.coupled_phase_pi_section_rows) || return false
    rows = result.coupled_phase_pi_section_rows[first_index:last_index]
    return all(rows) do row
        row.reference_kind != :none ||
            (
                length(row.raw_resistance_values) >= row.phase_index &&
                length(row.raw_inductance_values) >= row.phase_index &&
                length(row.raw_capacitance_values) >= row.phase_index
            )
    end
end

function _cascaded_pi_detail_stage(
    series_modifier_count::Int,
    shunt_modifier_count::Int,
    explicit_section_count::Int,
)
    series_modifier_count > 0 && return :series
    shunt_modifier_count > 0 && return :shunt
    explicit_section_count > 0 && return :explicit_section
    return :complete
end

function _cascaded_pi_next_detail_stage(configuration, completed::Symbol)
    completed == :series && configuration.shunt_modifier_count > 0 && return :shunt
    completed in (:series, :shunt) && configuration.explicit_section_count > 0 &&
        return :explicit_section
    return :complete
end

function _cascaded_pi_validate_phase_map!(
    result::DeckParseResult,
    phase_map::AbstractVector{Int},
    phase_count::Int,
    line_no::Int,
)
    sort(phase_map) == collect(1:phase_count) ||
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "cascaded PI phase map must be a permutation of its phase indices",
        ))
    return phase_map
end

function _cascaded_pi_configuration_start_card(pending, image::AbstractString)
    phase_slots = min(pending.phase_count, 14)
    phase_map = Union{Nothing,Int}[
        fixed_int_value(image, 25 + 4 * (phase - 1), 28 + 4 * (phase - 1))
        for phase in 1:phase_slots
    ]
    all(value -> value !== nothing, phase_map) || return false
    values = Int[value for value in phase_map]
    all(phase -> 1 <= phase <= pending.phase_count, values) || return false
    length(unique(values)) == length(values) || return false
    return fixed_float_value(image, 3, 8) !== nothing &&
        fixed_int_value(image, 9, 12) !== nothing
end

function parse_cascaded_pi_configuration!(
    result::DeckParseResult,
    pending,
    line::AbstractString,
    line_no::Int,
)
    image = fixed_image(line)
    raw_scale = fixed_float_value(image, 3, 8)
    raw_scale === nothing && add_issue!(result.validation, missing_data(
        "line $line_no",
        "expected cascaded PI section scale in columns 3-8",
    ))
    scale = raw_scale === nothing || raw_scale == 0.0 ? 1.0 : Float64(raw_scale)
    section_count = something(fixed_int_value(image, 9, 12), 1)
    series_modifier_count = something(fixed_int_value(image, 13, 16), 0)
    shunt_modifier_count = something(fixed_int_value(image, 17, 20), 0)
    explicit_section_count = something(fixed_int_value(image, 21, 24), 0)
    phase_map = Int[]
    phase_slots = min(pending.phase_count, 14)
    for phase in 1:phase_slots
        value = fixed_int_value(image, 25 + 4 * (phase - 1), 28 + 4 * (phase - 1))
        value === nothing && add_issue!(result.validation, missing_data(
            "line $line_no",
            "expected cascaded PI phase map value $phase",
        ))
        push!(phase_map, something(value, 0))
    end
    isfinite(scale) && scale > 0.0 || add_issue!(result.validation, invalid_value(
        "line $line_no",
        "cascaded PI section scale must be finite and positive",
    ))
    section_count > 0 || add_issue!(result.validation, invalid_value(
        "line $line_no",
        "cascaded PI multiplicity must be positive",
    ))
    pending.phase_count <= 14 &&
        _cascaded_pi_validate_phase_map!(result, phase_map, pending.phase_count, line_no)
    all(>=(0), (series_modifier_count, shunt_modifier_count, explicit_section_count)) ||
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "cascaded PI modifier flags must be nonnegative",
        ))
    if isempty(pending.blocks) &&
       (series_modifier_count > 0 || shunt_modifier_count > 0 || explicit_section_count > 0)
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "the first cascaded PI block must use its admitted source section without preceding modifiers",
        ))
    end
    record_card!(result, :fixed_field)
    record_card!(result, :cascaded_pi_configuration)
    return merge(
        pending,
        (
            configuration = (
                line_no = line_no,
                section_scale = scale,
                multiplicity = section_count,
                phase_map = phase_map,
                series_modifier_count = series_modifier_count,
                shunt_modifier_count = shunt_modifier_count,
                explicit_section_count = explicit_section_count,
                series_impedances = DeckCascadedPiSeriesImpedanceRow[],
                shunt_impedances = DeckCascadedPiShuntImpedanceRow[],
                explicit_values = NTuple{3,Float64}[],
                detail_line_numbers = Int[],
                detail_stage = pending.phase_count > 14 ? :phase_map :
                    _cascaded_pi_detail_stage(
                        series_modifier_count,
                        shunt_modifier_count,
                        explicit_section_count,
                    ),
            ),
        ),
    )
end

function _parse_cascaded_pi_phase_map_continuation!(
    result::DeckParseResult,
    pending,
    line::AbstractString,
    line_no::Int,
)
    configuration = pending.configuration
    image = fixed_image(line)
    remaining = pending.phase_count - length(configuration.phase_map)
    remaining > 0 || throw(ArgumentError("cascade phase map is already complete"))
    value_count = min(remaining, 14)
    phase_map = copy(configuration.phase_map)
    for offset in 0:(value_count - 1)
        value = fixed_int_value(image, 25 + 4 * offset, 28 + 4 * offset)
        value === nothing && add_issue!(result.validation, missing_data(
            "line $line_no",
            "expected cascaded PI phase map continuation value $(length(phase_map) + 1)",
        ))
        push!(phase_map, something(value, 0))
    end
    complete = length(phase_map) == pending.phase_count
    complete && _cascaded_pi_validate_phase_map!(
        result,
        phase_map,
        pending.phase_count,
        line_no,
    )
    detail_stage = complete ?
        _cascaded_pi_detail_stage(
            configuration.series_modifier_count,
            configuration.shunt_modifier_count,
            configuration.explicit_section_count,
        ) : :phase_map
    push!(configuration.detail_line_numbers, line_no)
    record_card!(result, :cascaded_pi_phase_map_continuation)
    return merge(pending, (configuration = merge(configuration, (
        phase_map = phase_map,
        detail_stage = detail_stage,
    )),))
end

function _cascaded_pi_impedance_values!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    label::AbstractString,
)
    resistance = fixed_float_field!(
        result,
        image,
        line_no,
        27,
        32,
        "$(label)_resistance_ohm",
    )
    inductance = fixed_float_field!(
        result,
        image,
        line_no,
        33,
        38,
        "$(label)_inductance_value",
    )
    capacitance = fixed_float_field!(
        result,
        image,
        line_no,
        39,
        44,
        "$(label)_capacitance_value",
    )
    return (
        something(resistance, 0.0),
        something(inductance, 0.0),
        something(capacitance, 0.0),
    )
end

function _parse_cascaded_pi_series_impedance!(
    result::DeckParseResult,
    pending,
    line::AbstractString,
    line_no::Int,
)
    configuration = pending.configuration
    image = fixed_image(line)
    phase_index = something(fixed_int_value(image, 1, 2), 0)
    if phase_index == 0
        record_card!(result, :cascaded_pi_series_terminator)
        updated = merge(pending, (configuration = merge(configuration, (
            detail_stage = _cascaded_pi_next_detail_stage(configuration, :series),
        )),))
        strip(line) == "0" && return updated
        return parse_cascaded_pi_detail!(result, updated, line, line_no)
    end
    1 <= phase_index <= pending.phase_count || add_issue!(
        result.validation,
        invalid_value(
            "line $line_no",
            "cascaded PI series modifier phase must be between one and $(pending.phase_count)",
        ),
    )
    any(row -> row.phase_index == phase_index, configuration.series_impedances) &&
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "cascaded PI series modifier phase $phase_index is duplicated",
        ))
    resistance, inductance, capacitance = _cascaded_pi_impedance_values!(
        result,
        image,
        line_no,
        "cascaded_pi_series",
    )
    open_circuit = resistance == 999999.0
    node_continuity =
        !open_circuit && resistance == 0.0 && inductance == 0.0 && capacitance == 0.0
    node_continuity || push!(
        configuration.series_impedances,
        DeckCascadedPiSeriesImpedanceRow(
            line_no,
            phase_index,
            open_circuit ? 0.0 : resistance,
            inductance,
            capacitance,
            open_circuit,
        ),
    )
    push!(configuration.detail_line_numbers, line_no)
    record_card!(result, :cascaded_pi_series_impedance)
    return pending
end

function _parse_cascaded_pi_shunt_impedance!(
    result::DeckParseResult,
    pending,
    line::AbstractString,
    line_no::Int,
)
    configuration = pending.configuration
    image = fixed_image(line)
    raw_from_terminal = fixed_int_value(image, 3, 8)
    raw_to_terminal = fixed_int_value(image, 9, 14)
    next_configuration = _cascaded_pi_configuration_start_card(pending, image)
    if (raw_from_terminal === nothing && raw_to_terminal === nothing) ||
       next_configuration
        record_card!(result, :cascaded_pi_shunt_terminator)
        updated = merge(pending, (configuration = merge(configuration, (
            detail_stage = _cascaded_pi_next_detail_stage(configuration, :shunt),
        )),))
        return parse_cascaded_pi_detail!(result, updated, line, line_no)
    end
    from_terminal = fixed_int_field!(
        result, image, line_no, 3, 8, "cascaded_pi_shunt_from_terminal",
    )
    to_terminal = fixed_int_field!(
        result, image, line_no, 9, 14, "cascaded_pi_shunt_to_terminal",
    )
    from_index = something(from_terminal, 0)
    to_index = something(to_terminal, 0)
    abs(from_index) <= pending.phase_count && abs(to_index) <= pending.phase_count ||
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "cascaded PI shunt terminal magnitude must not exceed $(pending.phase_count)",
        ))
    from_index != to_index || add_issue!(result.validation, invalid_value(
        "line $line_no",
        "cascaded PI shunt terminals must be distinct",
    ))
    resistance, inductance, capacitance = _cascaded_pi_impedance_values!(
        result,
        image,
        line_no,
        "cascaded_pi_shunt",
    )
    resistance > 0.0 || inductance > 0.0 || capacitance > 0.0 ||
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "cascaded PI shunt impedance cannot be zero",
        ))
    push!(
        configuration.shunt_impedances,
        DeckCascadedPiShuntImpedanceRow(
            line_no,
            from_index,
            to_index,
            resistance,
            inductance,
            capacitance,
        ),
    )
    push!(configuration.detail_line_numbers, line_no)
    record_card!(result, :cascaded_pi_shunt_impedance)
    return pending
end

function _parse_cascaded_pi_explicit_section!(
    result::DeckParseResult,
    pending,
    line::AbstractString,
    line_no::Int,
)
    configuration = pending.configuration
    image = fixed_image(line)
    values = Float64[]
    for column in 27:6:80
        field = fixed_field(image, column, min(column + 5, 80))
        isempty(field) && continue
        value = tryparse_deck_float(field)
        value === nothing && add_issue!(result.validation, invalid_value(
            "line $line_no",
            "invalid cascaded PI explicit-section matrix value in columns $column-$(min(column + 5, 80))",
        ))
        value === nothing || push!(values, Float64(value))
    end
    length(values) % 3 == 0 || add_issue!(result.validation, invalid_value(
        "line $line_no",
        "cascaded PI explicit-section values must be resistance/inductance/capacitance triples",
    ))
    for index in 1:3:(length(values) - 2)
        push!(configuration.explicit_values, Tuple(values[index:(index + 2)]))
    end
    push!(configuration.detail_line_numbers, line_no)
    expected = pending.phase_count * (pending.phase_count + 1) ÷ 2
    length(configuration.explicit_values) <= expected || add_issue!(
        result.validation,
        invalid_value(
            "line $line_no",
            "cascaded PI explicit section contains more than $expected matrix triples",
        ),
    )
    if length(configuration.explicit_values) >= expected
        configuration = merge(configuration, (detail_stage = :complete,))
    end
    record_card!(result, :cascaded_pi_explicit_section)
    return merge(pending, (configuration = configuration,))
end

function _complete_cascaded_pi_block!(result::DeckParseResult, pending)
    configuration = pending.configuration
    configuration.detail_stage == :complete ||
        throw(ArgumentError("cascaded PI block details are incomplete"))
    resistance = nothing
    inductance = nothing
    capacitance = nothing
    if configuration.explicit_section_count > 0
        expected = pending.phase_count * (pending.phase_count + 1) ÷ 2
        length(configuration.explicit_values) == expected || add_issue!(
            result.validation,
            missing_data(
                "line $(configuration.line_no)",
                "cascaded PI explicit section requires $expected matrix triples",
            ),
        )
        resistance = zeros(pending.phase_count, pending.phase_count)
        inductance = zeros(pending.phase_count, pending.phase_count)
        capacitance = zeros(pending.phase_count, pending.phase_count)
        value_index = 1
        for row in 1:pending.phase_count, column in 1:row
            r, x, c = configuration.explicit_values[value_index]
            resistance[row, column] = resistance[column, row] = r
            inductance[row, column] = inductance[column, row] = x
            capacitance[row, column] = capacitance[column, row] = c
            value_index += 1
        end
    end
    push!(
        pending.blocks,
        DeckCascadedPiBlock(
            configuration.line_no,
            configuration.section_scale,
            configuration.multiplicity,
            copy(configuration.phase_map),
            copy(configuration.series_impedances),
            copy(configuration.shunt_impedances),
            resistance,
            inductance,
            capacitance,
            copy(configuration.detail_line_numbers),
        ),
    )
    return merge(pending, (configuration = nothing,))
end

function parse_cascaded_pi_detail!(
    result::DeckParseResult,
    pending,
    line::AbstractString,
    line_no::Int,
)
    stage = pending.configuration.detail_stage
    stage == :phase_map &&
        return _parse_cascaded_pi_phase_map_continuation!(result, pending, line, line_no)
    stage == :series &&
        return _parse_cascaded_pi_series_impedance!(result, pending, line, line_no)
    stage == :shunt &&
        return _parse_cascaded_pi_shunt_impedance!(result, pending, line, line_no)
    stage == :explicit_section &&
        return _parse_cascaded_pi_explicit_section!(result, pending, line, line_no)
    stage == :complete || throw(ArgumentError("unknown cascaded PI detail stage $stage"))
    completed = _complete_cascaded_pi_block!(result, pending)
    if cascaded_pi_terminator_card(line)
        finish_cascaded_pi_request!(result, completed, line_no)
        return nothing
    end
    return parse_cascaded_pi_configuration!(result, completed, line, line_no)
end

function cascaded_pi_terminator_card(line::AbstractString)::Bool
    image = fixed_image(line)
    return uppercase(fixed_field(image, 3, 14)) == "STOP CASCADE"
end

function finish_cascaded_pi_request!(
    result::DeckParseResult,
    pending,
    line_no::Int,
)
    isempty(pending.blocks) && throw(ArgumentError("cascaded PI configuration is missing"))
    first_configuration = first(pending.blocks)
    last_row_index = length(result.coupled_phase_pi_section_rows)
    row_indices = collect(pending.first_section_row_index:last_row_index)
    length(row_indices) == pending.phase_count || add_issue!(result.validation, invalid_value(
        "line $line_no",
        "cascaded PI source section row count must match its phase count",
    ))
    if length(row_indices) == pending.phase_count
        rows = result.coupled_phase_pi_section_rows[row_indices]
        [row.phase_index for row in rows] == collect(1:pending.phase_count) ||
            add_issue!(result.validation, invalid_value(
                "line $line_no",
                "cascaded PI source section phases must be contiguous from one",
            ))
    end
    push!(
        result.cascaded_pi_request_rows,
        DeckCascadedPiRequestRow(
            Symbol("cascaded_phase_pi_", length(result.cascaded_pi_request_rows) + 1),
            pending.header_line_no,
            first_configuration.configuration_line_no,
            line_no,
            pending.phase_count,
            pending.frequency_hz,
            first_configuration.section_scale,
            sum(block.multiplicity for block in pending.blocks),
            copy(first_configuration.phase_map),
            row_indices,
            copy(pending.blocks),
        ),
    )
    record_card!(result, :fixed_field)
    record_card!(result, :cascaded_pi_terminator)
    record_card!(result, :cascaded_pi_request)
    return result
end

function start_sampled_frequency_line_parse(
    result::DeckParseResult,
    branch_row_index::Int,
)
    1 <= branch_row_index <= length(result.coupled_line_rows) ||
        throw(ArgumentError("sampled frequency line branch index is outside the coupled-line table"))
    row = result.coupled_line_rows[branch_row_index]
    return (
        branch_row_index = branch_row_index,
        branch_line_no = row.line_no,
        configuration = nothing,
        propagation_time_s = Float64[],
        propagation_amplitude = Float64[],
        admittance_time_s = Float64[],
        admittance_amplitude = Float64[],
        sample_line_numbers = Int[],
    )
end

function parse_sampled_frequency_line_configuration!(
    result::DeckParseResult,
    pending,
    line::AbstractString,
    line_no::Int,
)
    image = fixed_image(line)
    propagation_count = fixed_int_field!(
        result, image, line_no, 1, 8, "sampled_line_propagation_point_count",
    )
    admittance_count = fixed_int_field!(
        result, image, line_no, 9, 16, "sampled_line_admittance_point_count",
    )
    propagation_peak_index = fixed_int_field!(
        result, image, line_no, 17, 24, "sampled_line_propagation_peak_index",
    )
    admittance_rise_index = fixed_int_field!(
        result, image, line_no, 25, 32, "sampled_line_admittance_rise_index",
    )
    characteristic_impedance = fixed_float_field!(
        result, image, line_no, 33, 40, "sampled_line_characteristic_impedance_ohm",
    )
    propagation_cutoff = something(fixed_float_value(image, 41, 48), 0.01)
    admittance_cutoff = something(fixed_float_value(image, 49, 56), 0.1)
    total_resistance = something(fixed_float_value(image, 65, 72), 0.0)
    maximum_iterations = something(fixed_int_value(image, 73, 80), 100)
    if propagation_count !== nothing && propagation_count < 2
        add_issue!(result.validation, invalid_value(
            "line $line_no", "sampled propagation weighting requires at least two points",
        ))
    end
    if admittance_count !== nothing && admittance_count < 2
        add_issue!(result.validation, invalid_value(
            "line $line_no", "sampled admittance weighting requires at least two points",
        ))
    end
    if propagation_count !== nothing && propagation_peak_index !== nothing &&
       !(1 <= propagation_peak_index <= propagation_count)
        add_issue!(result.validation, invalid_value(
            "line $line_no", "sampled propagation peak index is outside its point table",
        ))
    end
    if admittance_count !== nothing && admittance_rise_index !== nothing &&
       !(1 <= admittance_rise_index <= admittance_count)
        add_issue!(result.validation, invalid_value(
            "line $line_no", "sampled admittance rise index is outside its point table",
        ))
    end
    if characteristic_impedance !== nothing &&
       !(isfinite(characteristic_impedance) && characteristic_impedance > 0.0)
        add_issue!(result.validation, invalid_value(
            "line $line_no", "sampled line characteristic impedance must be positive",
        ))
    end
    propagation_cutoff > 0.0 || (propagation_cutoff = 0.01)
    admittance_cutoff > 0.0 || (admittance_cutoff = 0.1)
    maximum_iterations > 0 || (maximum_iterations = 100)
    record_card!(result, :sampled_frequency_line_configuration)
    return merge(pending, (configuration = (
        line_no = line_no,
        propagation_count = something(propagation_count, 0),
        admittance_count = something(admittance_count, 0),
        propagation_peak_index = something(propagation_peak_index, 0),
        admittance_rise_index = something(admittance_rise_index, 0),
        characteristic_impedance_ohm = something(characteristic_impedance, 0.0),
        propagation_cutoff_fraction = Float64(propagation_cutoff),
        admittance_cutoff_fraction = Float64(admittance_cutoff),
        total_resistance_ohm = Float64(total_resistance),
        maximum_tail_iterations = Int(maximum_iterations),
    ),))
end

function _sampled_frequency_line_point_pairs(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
)
    image = fixed_image(line)
    pairs = Tuple{Float64,Float64}[]
    for pair_index in 1:4
        first_col = 20 * (pair_index - 1) + 1
        time_field = fixed_field(image, first_col, first_col + 9)
        amplitude_field = fixed_field(image, first_col + 10, first_col + 19)
        isempty(time_field) && isempty(amplitude_field) && continue
        if isempty(time_field) || isempty(amplitude_field)
            add_issue!(result.validation, missing_data(
                "line $line_no", "sampled line weighting points require time/amplitude pairs",
            ))
            continue
        end
        time_us = tryparse_deck_float(time_field)
        amplitude = tryparse_deck_float(amplitude_field)
        if time_us === nothing || amplitude === nothing
            add_issue!(result.validation, invalid_value(
                "line $line_no", "invalid sampled line weighting time/amplitude pair",
            ))
            continue
        end
        push!(pairs, (Float64(time_us) / 1.0e6, Float64(amplitude)))
    end
    isempty(pairs) && add_issue!(result.validation, missing_data(
        "line $line_no", "expected at least one sampled line weighting point",
    ))
    return pairs
end

function parse_sampled_frequency_line_points!(
    result::DeckParseResult,
    pending,
    line::AbstractString,
    line_no::Int,
)
    config = pending.configuration
    config === nothing && throw(ArgumentError("sampled line configuration must precede points"))
    pairs = _sampled_frequency_line_point_pairs(result, line, line_no)
    push!(pending.sample_line_numbers, line_no)
    for (time_s, amplitude) in pairs
        if length(pending.propagation_time_s) < config.propagation_count
            push!(pending.propagation_time_s, time_s)
            push!(pending.propagation_amplitude, amplitude)
        elseif length(pending.admittance_time_s) < config.admittance_count
            push!(pending.admittance_time_s, time_s)
            push!(pending.admittance_amplitude, amplitude)
        else
            add_issue!(result.validation, invalid_value(
                "line $line_no", "sampled line point cards contain more values than declared",
            ))
        end
    end
    complete = length(pending.propagation_time_s) == config.propagation_count &&
               length(pending.admittance_time_s) == config.admittance_count
    complete || return pending

    branch = result.coupled_line_rows[pending.branch_row_index]
    if !all(diff(pending.propagation_time_s) .> 0.0) ||
       !all(diff(pending.admittance_time_s) .> 0.0)
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "sampled frequency-dependent line times must increase strictly within each table",
        ))
    end
    if any(value -> value < 0.0, pending.propagation_amplitude) ||
       any(value -> value < 0.0, pending.admittance_amplitude)
        add_issue!(result.validation, invalid_value(
            "line $line_no",
            "sampled frequency-dependent line amplitudes must be nonnegative",
        ))
    end
    push!(
        result.sampled_frequency_line_rows,
        DeckSampledFrequencyLineRow(
            Symbol("sampled_frequency_line_", length(result.sampled_frequency_line_rows) + 1),
            pending.branch_row_index,
            pending.branch_line_no,
            config.line_no,
            copy(pending.sample_line_numbers),
            branch.from_node,
            branch.to_node,
            branch.from_node_value,
            branch.to_node_value,
            config.propagation_peak_index,
            config.admittance_rise_index,
            config.characteristic_impedance_ohm,
            config.propagation_cutoff_fraction,
            config.admittance_cutoff_fraction,
            config.total_resistance_ohm,
            config.maximum_tail_iterations,
            copy(pending.propagation_time_s),
            copy(pending.propagation_amplitude),
            copy(pending.admittance_time_s),
            copy(pending.admittance_amplitude),
        ),
    )
    record_card!(result, :sampled_frequency_line_points)
    record_card!(result, :sampled_frequency_line_request)
    return nothing
end

function fixed_card_branch_vintage_mode(tokens)
    length(tokens) >= 2 || return nothing
    first_token = replace(compact_deck_keyword(tokens[1]), "\$" => "")
    first_token == "vintage" || return nothing
    return tryparse(Int, deck_token_value(tokens[2]))
end

function fixed_card_numeric_token_values(tokens)
    isempty(tokens) && return nothing
    values = Float64[]
    for token in tokens
        value = tryparse_deck_float(token)
        value === nothing && return nothing
        push!(values, Float64(value))
    end
    return values
end

function fixed_card_kc_lee_transform_row_candidate(tokens)::Bool
    values = fixed_card_numeric_token_values(tokens)
    return values !== nothing && 1 <= length(values) <= 30
end

function parse_fixed_card_kc_lee_transform_row!(
    result::DeckParseResult,
    tokens,
    line_no::Int,
)::Bool
    values = fixed_card_numeric_token_values(tokens)
    values === nothing && return false
    record_card!(result, :fixed_field)
    record_card!(result, :fixed_card_kc_lee_modal_transform_row)
    for _ in values
        record_card!(result, :fixed_card_kc_lee_modal_transform_value)
    end
    push!(
        result.line_modal_transform_rows,
        DeckLineModalTransformRow(line_no, Float64.(values), join(tokens, " ")),
    )
    return true
end

function compact_leading_sign_field(raw::AbstractString)::String
    text = strip(String(raw))
    isempty(text) && return text
    first_char = first(text)
    if first_char == '+' || first_char == '-'
        rest = lstrip(text[nextind(text, firstindex(text)):lastindex(text)])
        return string(first_char, rest)
    end
    return text
end

function fixed_card_kc_lee_optional_float_field(
    image::AbstractString,
    first_col::Int,
    last_col::Int,
)
    raw = fixed_field(image, first_col, last_col)
    isempty(raw) && return missing
    value = tryparse_deck_float(compact_leading_sign_field(raw))
    return value === nothing ? missing : Float64(value)
end

function fixed_card_optional_float_field!(result::DeckParseResult,
                                          image::AbstractString,
                                          line_no::Int,
                                          first_col::Int,
                                          last_col::Int,
                                          field::AbstractString)
    isempty(fixed_field(image, first_col, last_col)) && return missing
    value = fixed_float_field!(result, image, line_no, first_col, last_col, field)
    return value === nothing ? missing : Float64(value)
end

function push_fixed_card_deferred_coupled_line_row!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    line_type::Int,
    line_kind::Symbol,
    phase_index::Int,
    from_node::AbstractString,
    to_node::AbstractString,
    reference_from::AbstractString,
    reference_to::AbstractString,
    raw_resistance,
    raw_inductance,
    raw_capacitance,
    line_length,
)::Bool
    from_index = node_id!(result, from_node)
    to_index = node_id!(result, to_node)
    reference_from_name, reference_from_index =
        fixed_card_optional_node!(result, reference_from)
    reference_to_name, reference_to_index =
        fixed_card_optional_node!(result, reference_to)
    push!(
        result.coupled_line_rows,
        DeckCoupledLineRow(
            Symbol("coupled_line_fixed_", length(result.coupled_line_rows) + 1),
            Symbol(String(from_node)),
            Symbol(String(to_node)),
            from_index,
            to_index,
            reference_from_name,
            reference_to_name,
            reference_from_index,
            reference_to_index,
            line_no,
            line_type,
            line_kind,
            phase_index,
            missing,
            missing,
            Float64[],
            Float64[],
            raw_resistance,
            raw_inductance,
            raw_capacitance,
            line_length,
            false,
        ),
    )
    return true
end

function parse_fixed_card_terminal_surge_impedance_row!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    branch_type::Int,
    from_node::AbstractString,
    to_node::AbstractString,
    initial_issues::Int,
)::Bool
    values = branch_numeric_tail_values(image)
    if isempty(values)
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "expected terminal surge-impedance numeric values in columns 27-74",
            ),
        )
        return true
    end
    push_fixed_card_deferred_coupled_line_row!(
        result,
        image,
        line_no,
        branch_type,
        :terminal_surge_impedance,
        branch_type,
        from_node,
        to_node,
        "",
        "",
        Float64(values[1]),
        length(values) >= 2 ? Float64(values[2]) : missing,
        length(values) >= 3 ? Float64(values[3]) : missing,
        missing,
    )
    record_fixed_card!(
        result,
        :bpa_fixed_branch,
        :fixed_card_terminal_surge_impedance_row,
        initial_issues,
    )
    return true
end

function parse_fixed_card_kc_lee_untransposed_line_row!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    line_type::Int,
    from_node::AbstractString,
    to_node::AbstractString,
    reference_from::AbstractString,
    reference_to::AbstractString,
    initial_issues::Int,
)::Bool
    phase_index = abs(line_type)
    push_fixed_card_deferred_coupled_line_row!(
        result,
        image,
        line_no,
        line_type,
        :kc_lee_untransposed_line,
        phase_index,
        from_node,
        to_node,
        reference_from,
        reference_to,
        fixed_card_kc_lee_optional_float_field(image, 27, 38),
        fixed_card_kc_lee_optional_float_field(image, 39, 50),
        fixed_card_kc_lee_optional_float_field(image, 51, 62),
        fixed_card_kc_lee_optional_float_field(image, 63, 74),
    )
    record_fixed_card!(
        result,
        :bpa_fixed_branch,
        :fixed_card_kc_lee_untransposed_line_row,
        initial_issues,
    )
    record_card!(result, Symbol("fixed_card_kc_lee_untransposed_line_phase_", phase_index))
    if isempty(reference_from) && isempty(reference_to)
        record_card!(result, :fixed_card_kc_lee_untransposed_line_modal_parameters)
    else
        record_card!(result, :fixed_card_kc_lee_untransposed_line_copy_reference)
    end
    return true
end

function fixed_card_optional_node!(
    result::DeckParseResult,
    node::AbstractString,
)
    isempty(node) && return missing, missing
    node_index = node_id!(result, node)
    return Symbol(String(node)), node_index
end

function parse_fixed_card_coupled_line_row!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    line_type::Int,
    from_node::AbstractString,
    to_node::AbstractString,
    reference_from::AbstractString,
    reference_to::AbstractString,
)::Bool
    descriptor = fixed_card_coupled_line_descriptor(line_type)
    descriptor === nothing && return false
    from_index = node_id!(result, from_node)
    to_index = node_id!(result, to_node)
    reference_from_name, reference_from_index =
        fixed_card_optional_node!(result, reference_from)
    reference_to_name, reference_to_index =
        fixed_card_optional_node!(result, reference_to)
    if descriptor.kind == :mutual_source_equivalent
        values = fixed_card_coupled_lumped_pair_values!(result, image, line_no)
        triangular_resistance_values = values.resistance_values
        triangular_inductance_values = values.inductance_values
        sequence_resistance =
            isempty(triangular_resistance_values) ?
            missing :
            first(triangular_resistance_values)
        sequence_inductance =
            isempty(triangular_inductance_values) ?
            missing :
            first(triangular_inductance_values)
        raw_resistance = sequence_resistance
        raw_inductance = sequence_inductance
        raw_capacitance = missing
        line_length = missing
    else
        sequence_resistance = missing
        sequence_inductance = missing
        triangular_resistance_values = Float64[]
        triangular_inductance_values = Float64[]
        raw_resistance = fixed_card_optional_float_field!(
            result,
            image,
            line_no,
            27,
            32,
            "coupled_line_resistance",
        )
        raw_inductance = fixed_card_optional_float_field!(
            result,
            image,
            line_no,
            33,
            38,
            "coupled_line_inductance",
        )
        raw_capacitance = fixed_card_optional_float_field!(
            result,
            image,
            line_no,
            39,
            44,
            "coupled_line_capacitance",
        )
        line_length = fixed_card_optional_float_field!(
            result,
            image,
            line_no,
            45,
            50,
            "coupled_line_length",
        )
    end
    sampled_frequency_data_requested = descriptor.kind == :distributed_transmission_line &&
        something(fixed_int_value(image, 53, 54), 0) < 0
    name = Symbol("coupled_line_fixed_", length(result.coupled_line_rows) + 1)
    push!(
        result.coupled_line_rows,
        DeckCoupledLineRow(
            name,
            Symbol(String(from_node)),
            Symbol(String(to_node)),
            from_index,
            to_index,
            reference_from_name,
            reference_to_name,
            reference_from_index,
            reference_to_index,
            line_no,
            line_type,
            descriptor.kind,
            descriptor.phase_index,
            sequence_resistance,
            sequence_inductance,
            triangular_resistance_values,
            triangular_inductance_values,
            raw_resistance,
            raw_inductance,
            raw_capacitance,
            line_length,
            sampled_frequency_data_requested,
        ),
    )
    record_card!(result, :fixed_card_coupled_line)
    record_card!(result, Symbol("fixed_card_coupled_line_phase_", descriptor.phase_index))
    record_card!(result, Symbol("fixed_card_coupled_line_", descriptor.kind))
    if !ismissing(reference_from_index) || !ismissing(reference_to_index)
        if ismissing(reference_from_index) || ismissing(reference_to_index)
            add_issue!(
                result.validation,
                missing_data(
                    "line $line_no",
                    "coupled-line COPY requires both reference terminals",
                ),
            )
        else
            record_card!(result, :fixed_card_coupled_line_copy_reference)
            record_card!(
                result,
                descriptor.kind == :mutual_source_equivalent ?
                :fixed_card_coupled_lumped_copy_reference :
                :fixed_card_distributed_line_copy_reference,
            )
        end
    end
    return true
end

function push_over2_branch_row!(
    result::DeckParseResult,
    element_index::Int,
    line_no::Int,
    branch_type::Int,
    layout_kind::Symbol,
    source_kind::Symbol,
    reference_kind::Symbol,
    reference_name::Symbol,
    reference_line_no::Int,
    raw_resistance::Float64,
    raw_inductance::Float64,
    raw_capacitance::Float64,
    output_code::Int,
)
    element = result.elements[element_index]
    node_names = _deck_node_names_for_indices(
        result,
        [_deck_branch_from_node(element), _deck_branch_to_node(element)],
    )
    push!(
        result.over2_branch_rows,
        DeckOVER2BranchRow(
            result.element_names[element_index],
            node_names[1],
            node_names[2],
            _deck_branch_from_node(element),
            _deck_branch_to_node(element),
            line_no,
            branch_type,
            deck_element_kind(element),
            layout_kind,
            source_kind,
            reference_kind,
            reference_name,
            reference_line_no,
            raw_resistance,
            raw_inductance,
            raw_capacitance,
            _deck_branch_conductance_value(element),
            _deck_branch_resistance_value(element),
            _deck_branch_inductance_value(element),
            _deck_branch_capacitance_value(element),
            output_code,
        ),
    )
    return result
end

function push_bergeron_line_row!(
    result::DeckParseResult,
    element_index::Int,
    line_no::Int,
)
    line = result.elements[element_index]
    node_names = _deck_node_names_for_indices(result, [line.a, line.b])
    history = line_history_currents(line)
    voltages = line_terminal_voltages(line)
    currents = line_terminal_currents(line)
    waves = line_traveling_waves(line)
    push!(
        result.bergeron_line_rows,
        DeckBergeronLineRow(
            result.element_names[element_index],
            node_names[1],
            node_names[2],
            line.a,
            line.b,
            line_no,
            line.zc,
            line_surge_admittance(line),
            line.travel_time_s,
            line.dt_s,
            line.attenuation,
            line.delay_steps,
            line.write_index,
            history.from,
            history.to,
            voltages.from,
            voltages.to,
            currents.from,
            currents.to,
            waves.from,
            waves.to,
        ),
    )
    return result
end

function parse_bpa_fixed_branch_copy_reference!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    from_node::AbstractString,
    to_node::AbstractString,
    branch_type::Int,
    reference_index::Int,
    reference_kind::Symbol,
    initial_issues::Int,
    ;
    output_code::Union{Nothing,Int}=nothing,
    branch_layout_kind::Symbol=:fixed_copy_reference,
)::Bool
    name = bpa_fixed_owner_name(result, "branch")
    from_index = node_id!(result, from_node)
    to_index = node_id!(result, to_node)
    element, copied_kind = bpa_fixed_branch_copy_element(
        result.elements[reference_index],
        from_index,
        to_index,
    )
    if element === nothing
        record_fixed_blocker!(result, :bpa_fixed_branch_blocked,
                              :bpa_fixed_branch_blocked_copy_reference)
        add_issue!(
            result.validation,
            unknown_field(
                "line $line_no",
                "Unsupported OVER2 fixed-field branch COPY reference target: only accepted scalar R/L/C owners are translated",
            ),
        )
        return true
    end
    push_element!(result, ["branch", name, from_node, to_node], element, line_no)
    record_fixed_card!(result, :bpa_fixed_branch, :bpa_fixed_branch_copy_reference, initial_issues)
    if length(result.validation.issues) == initial_issues
        record_card!(result, reference_kind)
        record_card!(result, copied_kind)
        resolved_output_code =
            output_code === nothing ? bpa_fixed_branch_output_code_from_image(image) : output_code
        if output_code === nothing
            record_bpa_fixed_branch_output_requests!(result, image, line_no, name)
        else
            record_bpa_branch_output_requests_by_code!(result, line_no, name, output_code)
        end
        if length(result.validation.issues) == initial_issues
            reference_name = result.element_names[reference_index]
            reference_line_no = result.element_line_numbers[reference_index]
            reference_row_index = findlast(result.over2_branch_rows) do row
                row.name == reference_name && row.line_no == reference_line_no
            end
            reference_row_index === nothing && throw(ArgumentError(
                "accepted OVER2 COPY reference is missing its typed branch row",
            ))
            reference_row = result.over2_branch_rows[reference_row_index]
            push_over2_branch_row!(
                result,
                length(result.elements),
                line_no,
                branch_type,
                branch_layout_kind,
                :copy_reference,
                reference_kind,
                result.element_names[reference_index],
                result.element_line_numbers[reference_index],
                reference_row.raw_resistance,
                reference_row.raw_inductance,
                reference_row.raw_capacitance,
                resolved_output_code,
            )
        end
    end
    return true
end

function enqueue_fixed_owner_name!(result::DeckParseResult, owner_kind::Symbol,
                                   name::AbstractString)
    queue = get!(result.pending_fixed_owner_names, owner_kind, String[])
    push!(queue, String(name))
    return result
end

function bpa_fixed_owner_name(result::DeckParseResult, prefix::AbstractString;
                              explicit_name::Union{Nothing,AbstractString}=nothing)::String
    count_key = Symbol("bpa_fixed_", prefix, "_owner")
    next_index = get(result.card_counts, count_key, 0) + 1
    record_card!(result, count_key)
    owner_kind = Symbol(prefix)
    queue = get!(result.pending_fixed_owner_names, owner_kind, String[])
    if explicit_name !== nothing
        isempty(queue) || popfirst!(queue)
        return String(explicit_name)
    end
    if !isempty(queue)
        return popfirst!(queue)
    end
    return string(prefix, "_fixed_", next_index)
end

function bpa_fixed_nonlinear_owner_name(
    result::DeckParseResult,
    fallback_name::AbstractString;
    explicit_name::Union{Nothing,AbstractString}=nothing,
)::Symbol
    queue = get!(result.pending_fixed_owner_names, :nonlinear, String[])
    if explicit_name !== nothing
        isempty(queue) || popfirst!(queue)
        return Symbol(String(explicit_name))
    end
    isempty(queue) || return Symbol(popfirst!(queue))
    return Symbol(String(fallback_name))
end

function record_fixed_card!(result::DeckParseResult, kind::Symbol, subtype::Symbol,
                            initial_issue_count::Int)
    length(result.validation.issues) == initial_issue_count || return result
    record_card!(result, :fixed_field)
    record_card!(result, kind)
    return record_card!(result, subtype)
end

function record_fixed_blocker!(result::DeckParseResult, family::Symbol, blocker::Symbol)
    record_card!(result, :fixed_field)
    record_card!(result, family)
    return record_card!(result, blocker)
end

function push_bpa_fixed_branch_scalar_row!(
    result::DeckParseResult,
    from_node::AbstractString,
    to_node::AbstractString,
    branch_type::Int,
    resistance::Float64,
    inductance::Float64,
    capacitance::Float64,
    line_no::Int,
    initial_issues::Int;
    branch_layout_count::Union{Nothing,Symbol}=nothing,
    branch_layout_kind::Union{Nothing,Symbol}=nothing,
    inline_name::Union{Nothing,AbstractString}=nothing,
    output_image::Union{Nothing,AbstractString}=nothing,
    output_code::Union{Nothing,Int}=nothing,
)::Bool
    name = bpa_fixed_owner_name(result, "branch"; explicit_name=inline_name)
    element_count_before = length(result.elements)
    timestep_inductance = fixed_card_branch_timestep_inductance(result, inductance)
    timestep_capacitance = fixed_card_branch_timestep_capacitance(
        result,
        capacitance;
        legacy_microfarad_units = branch_layout_kind != :free_field,
    )
    if !all(
        isfinite,
        (
            resistance,
            inductance,
            capacitance,
            timestep_inductance,
            timestep_capacitance,
        ),
    )
        add_issue!(result.validation,
                   invalid_value("line $line_no",
                                 "OVER2 scalar branch R/L/C values must be finite"))
        return true
    end
    if resistance > 0.0 && inductance == 0.0 && capacitance == 0.0
        parse_resistor!(result, ["resistor", name, from_node, to_node, string(resistance)], line_no)
        record_fixed_card!(result, :bpa_fixed_branch, :bpa_fixed_branch_resistor, initial_issues)
    elseif resistance == 0.0 && inductance > 0.0 && capacitance == 0.0
        parse_inductor!(result, ["inductor", name, from_node, to_node, string(timestep_inductance)], line_no)
        record_fixed_card!(result, :bpa_fixed_branch, :bpa_fixed_branch_inductor, initial_issues)
    elseif resistance == 0.0 && inductance == 0.0 && capacitance > 0.0
        parse_capacitor!(
            result,
            ["capacitor", name, from_node, to_node, string(timestep_capacitance)],
            line_no,
        )
        record_fixed_card!(result, :bpa_fixed_branch, :bpa_fixed_branch_capacitor, initial_issues)
    elseif resistance != 0.0 && inductance > 0.0 && capacitance == 0.0
        active_denominator_valid = true
        if resistance < 0.0
            timestep_s = Float64(deck_fixed_time_horizon_options(result).dt_s)
            denominator =
                timestep_s > 0.0 ?
                resistance + 2.0 * timestep_inductance / timestep_s :
                NaN
            active_denominator_valid =
                isfinite(timestep_s) &&
                timestep_s > 0.0 &&
                isfinite(denominator) &&
                denominator != 0.0
            active_denominator_valid || add_issue!(
                result.validation,
                invalid_value(
                    "line $line_no",
                    "active series R-L requires a positive timestep and a finite nonzero companion denominator",
                ),
            )
        end
        if active_denominator_valid
            parse_series_rl!(
                result,
                [
                    "rl",
                    name,
                    from_node,
                    to_node,
                    string(resistance),
                    string(timestep_inductance),
                ],
                line_no,
            )
            record_fixed_card!(
                result,
                :bpa_fixed_branch,
                :bpa_fixed_branch_rl,
                initial_issues,
            )
            if resistance < 0.0 &&
               length(result.validation.issues) == initial_issues
                record_card!(result, :bpa_fixed_branch_active_rl)
            end
        end
    elseif resistance >= 0.0 && inductance >= 0.0 && capacitance > 0.0
        from_index = node_id!(result, from_node)
        to_index = node_id!(result, to_node)
        push_element!(
            result,
            ["rlc", name, from_node, to_node],
            SeriesRLCBranch(
                from_index,
                to_index,
                resistance,
                timestep_inductance,
                timestep_capacitance,
            ),
            line_no,
        )
        record_fixed_card!(result, :bpa_fixed_branch, :bpa_fixed_branch_rlc, initial_issues)
    else
        add_issue!(result.validation,
                   invalid_value("line $line_no",
                                 "unsupported accepted-owner OVER2 branch R/L/C combination"))
    end
    if length(result.validation.issues) == initial_issues
        branch_layout_count === nothing || record_card!(result, branch_layout_count)
        inline_name === nothing || record_card!(result, :bpa_fixed_branch_inline_name)
        resolved_output_code =
            output_code === nothing && output_image !== nothing ?
            bpa_fixed_branch_output_code_from_image(output_image) :
            (output_code === nothing ? 0 : output_code)
        if output_code !== nothing
            record_bpa_branch_output_requests_by_code!(result, line_no, name, output_code)
        elseif output_image !== nothing
            record_bpa_fixed_branch_output_requests!(result, output_image, line_no, name)
        end
        if length(result.validation.issues) == initial_issues &&
           length(result.elements) == element_count_before + 1
            push_over2_branch_row!(
                result,
                length(result.elements),
                line_no,
                branch_type,
                branch_layout_kind === nothing ?
                bpa_fixed_branch_layout_kind(branch_layout_count) :
                branch_layout_kind,
                :scalar,
                :none,
                :none,
                0,
                resistance,
                inductance,
                capacitance,
                resolved_output_code,
            )
        end
    end
    return true
end

function bpa_fixed_source_node_reference!(result::DeckParseResult, raw_node::AbstractString,
                                          line_no::Int)
    parsed_index = tryparse(Int, raw_node)
    if parsed_index !== nothing
        parsed_index != 0 || begin
            add_issue!(result.validation,
                       invalid_value("line $line_no",
                                     "OVER5A fixed-field source node index must be nonzero"))
            return Symbol(raw_node), 0
        end
        name = Symbol("N", abs(parsed_index))
        get!(result.node_map, name, length(result.node_map) + 1)
        return name, parsed_index
    end
    node_index = node_id!(result, raw_node)
    return Symbol(raw_node), node_index
end

function bpa_fixed_source_effective_iform(result::DeckParseResult, node_index::Int,
                                          source_type::Int)
    if node_index > 0 && source_type in BPA_FIXED_INCREMENTAL_IFORM_SOURCE_TYPES &&
       any(row -> row.node_value == node_index, result.over5a_source_rows)
        return -source_type
    end
    return source_type
end

function parse_bpa_fixed_type18_ungrounded_source_card!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    bus1::AbstractString,
    initial_issues::Int,
)::Bool
    ignored_crest = fixed_float_field!(result, image, line_no, 11, 20, "type18_ignored_crest")
    ignored_crest === nothing && return true

    bus2 = fixed_field(image, 21, 26)
    bus3 = fixed_field(image, 27, 32)
    if isempty(bus2) || isempty(bus3)
        record_fixed_blocker!(result, :bpa_fixed_source_blocked,
                              :bpa_fixed_source_blocked_type18_ungrounded_terminal)
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "OVER5A type-18 ungrounded source card requires BUS2 and BUS3 terminal names in columns 21-32",
            ),
        )
        return true
    end
    if isempty(result.over5a_source_rows)
        record_fixed_blocker!(result, :bpa_fixed_source_blocked,
                              :bpa_fixed_source_blocked_type18_ungrounded_terminal)
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "OVER5A type-18 ungrounded source card requires a previous source row to mirror",
            ),
        )
        return true
    end

    node_name, node_index = bpa_fixed_source_node_reference!(result, bus1, line_no)
    _, bus2_index = bpa_fixed_source_node_reference!(result, bus2, line_no)
    _, bus3_index = bpa_fixed_source_node_reference!(result, bus3, line_no)
    (node_index == 0 || bus2_index == 0 || bus3_index == 0) && return true

    previous = last(result.over5a_source_rows)
    name = bpa_fixed_owner_name(result, "source")
    push!(
        result.over5a_source_rows,
        DeckOVER5ASourceRow(
            Symbol(name),
            node_name,
            abs(node_index),
            line_no,
            :fixed_field,
            previous.iform,
            -previous.crest,
            previous.time1,
            Float64(abs(bus2_index)),
            previous.sfreq,
            previous.tstart,
            previous.tstop,
        ),
    )
    record_fixed_card!(result, :bpa_fixed_source, :bpa_fixed_source_type18, initial_issues)
    if length(result.validation.issues) == initial_issues
        record_card!(result, :bpa_fixed_source_type18_ungrounded_terminal)
        record_card!(result, :bpa_fixed_source_ungrounded_terminal_state)
        record_card!(result, :bpa_fixed_source_over16_state)
        if previous.tstart != 0.0 || isfinite(previous.tstop)
            record_card!(result, :bpa_fixed_source_window)
        end
    end
    return true
end

function remove_analytic_source_element!(
    result::DeckParseResult,
    source_name::Symbol,
)
    index = findlast(==(source_name), result.element_names)
    index === nothing && return false
    element = result.elements[index]
    element isa Union{TheveninSource,CurrentInjection} || return false
    deleteat!(result.elements, index)
    deleteat!(result.element_line_numbers, index)
    deleteat!(result.element_names, index)
    return true
end

function parse_bpa_fixed_type18_ideal_transformer_card!(
    result::DeckParseResult,
    image::AbstractString,
    line_no::Int,
    bus1::AbstractString,
    initial_issues::Int,
)::Bool
    ratio = fixed_float_field!(result, image, line_no, 11, 20, "type18_transformer_ratio")
    ratio === nothing && return true
    bus2 = fixed_field(image, 21, 26)
    bus3 = fixed_field(image, 27, 32)
    bus4 = fixed_field(image, 33, 38)
    if isempty(bus2) || isempty(bus3) || isempty(bus4)
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "OVER5A type-18 ideal-transformer source card requires BUS2, BUS3, and BUS4 terminal names in columns 21-38",
            ),
        )
        return true
    end
    isfinite(ratio) && ratio != 0.0 || begin
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "OVER5A type-18 ideal-transformer turns ratio must be finite and nonzero",
            ),
        )
        return true
    end
    isempty(result.over5a_source_rows) && begin
        add_issue!(
            result.validation,
            missing_data(
                "line $line_no",
                "OVER5A type-18 ideal-transformer source card requires a preceding source row",
            ),
        )
        return true
    end
    internal_name = Symbol(bus4)
    haskey(result.node_map, internal_name) && begin
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "OVER5A type-18 auxiliary node $bus4 must be a new node name",
            ),
        )
        return true
    end

    bus1_name, bus1_index = bpa_fixed_source_node_reference!(result, bus1, line_no)
    _, bus2_index = bpa_fixed_source_node_reference!(result, bus2, line_no)
    _, bus3_index = bpa_fixed_source_node_reference!(result, bus3, line_no)
    previous = last(result.over5a_source_rows)
    previous_node_index = abs(previous.node_value)
    previous_node_index > 0 || begin
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "OVER5A type-18 preceding source must own a non-ground node",
            ),
        )
        return true
    end
    internal_node_index = node_id!(result, bus4)
    name = Symbol(bpa_fixed_owner_name(result, "source"))
    constraint = IdealTransformerVoltageConstraint(
        bus1_index,
        bus2_index,
        bus3_index,
        previous_node_index,
        internal_node_index,
        ratio,
    )

    remove_analytic_source_element!(result, previous.name)
    push_element!(result, ["ideal_transformer_voltage_constraint", String(name)], constraint, line_no)
    push!(
        result.over5a_source_rows,
        DeckOVER5ASourceRow(
            name,
            bus1_name,
            bus1_index,
            line_no,
            :fixed_field,
            18,
            Float64(ratio),
            Float64(internal_node_index),
            Float64(bus2_index),
            Float64(bus3_index),
            0.0,
            Inf,
        ),
    )
    record_fixed_card!(result, :bpa_fixed_source, :bpa_fixed_source_type18, initial_issues)
    if length(result.validation.issues) == initial_issues
        record_card!(result, :bpa_fixed_source_type18_ideal_transformer)
        record_card!(result, :bpa_fixed_source_ideal_transformer_topology)
        record_card!(result, :bpa_fixed_source_over16_state)
    end
    return true
end

function bpa_fixed_source_card_image_candidate(image::AbstractString)::Bool
    raw = fixed_field(image, 1, 8)
    tryparse_deck_float(raw) === nothing && return false
    occursin(r"[.dDeE]", raw) && return true
    return isempty(fixed_field(image, 1, 2))
end

function bpa_fixed_source_card_e8_values!(result::DeckParseResult,
                                          image::AbstractString,
                                          line_no::Int)
    values = Float64[]
    provided_value_count = 0
    for first_col in 1:8:73
        raw = fixed_field(image, first_col, first_col + 7)
        if isempty(raw)
            push!(values, 0.0)
            continue
        end
        value = tryparse_deck_float(raw)
        if value === nothing
            add_issue!(
                result.validation,
                invalid_value(
                    "line $line_no",
                    "source_card_value_$(length(values) + 1)=$(raw): expected OVER16 fixed-field Float64 in FORMAT(10E8.0)",
                ),
            )
            return nothing
        end
        push!(values, Float64(value))
        provided_value_count += 1
    end
    return values, provided_value_count
end

function parse_bpa_fixed_source_card_image!(result::DeckParseResult,
                                            image::AbstractString,
                                            line_no::Int)::Bool
    value_row = bpa_fixed_source_card_e8_values!(result, image, line_no)
    value_row === nothing && return true
    values, provided_value_count = value_row
    push!(
        result.over16_source_card_rows,
        DeckOVER16SourceCardRow(:fixed_field, values, provided_value_count, line_no),
    )
    record_card!(result, :fixed_field)
    record_card!(result, :over16_source_card)
    record_card!(result, :over16_source_card_fixed_field)
    record_card!(result, :over16_source_card_fixed_field_image)
    provided_value_count == 10 || record_card!(result, :over16_source_card_zero_filled)
    return true
end

function bpa_fixed_source_card_free_field_image_candidate(line::AbstractString)::Bool
    occursin(',', String(line)) || return false
    fields = bpa_fixed_source_card_free_field_image_fields(line)
    return !isempty(fields)
end

function bpa_fixed_source_card_free_field_image_fields(line::AbstractString)
    raw_fields = split(String(line), ','; keepempty = true)
    last_nonempty = findlast(field -> !isempty(strip(field)), raw_fields)
    last_nonempty === nothing && return String[]
    return [replace(strip(field), r"\s+" => "") for field in raw_fields[1:last_nonempty]]
end

function parse_bpa_fixed_source_card_free_field_image!(result::DeckParseResult,
                                                       line::AbstractString,
                                                       line_no::Int)::Bool
    fields = bpa_fixed_source_card_free_field_image_fields(line)
    if length(fields) > 10
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "OVER16 free-field source-card rows accept at most ten numeric values"),
        )
        return true
    end
    values = Float64[]
    provided_value_count = 0
    for (index, field) in enumerate(fields)
        if isempty(field)
            push!(values, 0.0)
            continue
        end
        value = parse_float!(result, field, line_no, "source_card_value_$index")
        value === nothing && return true
        push!(values, Float64(value))
        provided_value_count += 1
    end
    while length(values) < 10
        push!(values, 0.0)
    end
    push!(
        result.over16_source_card_rows,
        DeckOVER16SourceCardRow(:free_field, values, provided_value_count, line_no),
    )
    record_card!(result, :fixed_field)
    record_card!(result, :over16_source_card)
    record_card!(result, :over16_source_card_free_field)
    record_card!(result, :over16_source_card_free_field_image)
    provided_value_count == 10 || record_card!(result, :over16_source_card_zero_filled)
    return true
end

function bpa_fixed_source_free_field_fields(line::AbstractString)
    raw_fields = split(String(line), ','; keepempty = true)
    last_nonempty = findlast(field -> !isempty(strip(field)), raw_fields)
    last_nonempty === nothing && return String[]
    return [String(strip(String(field))) for field in raw_fields[1:last_nonempty]]
end

bpa_fixed_branch_free_field_fields(line::AbstractString) =
    bpa_fixed_source_free_field_fields(line)
bpa_fixed_switch_free_field_fields(line::AbstractString) =
    bpa_fixed_source_free_field_fields(line)

function bpa_fixed_source_free_field_row_candidate(line::AbstractString)::Bool
    occursin(',', String(line)) || return false
    fields = bpa_fixed_source_free_field_fields(line)
    length(fields) >= 2 || return false
    source_type = tryparse(Int, fields[1])
    source_type === nothing && return false
    bpa_fixed_source_type_accepted(source_type) || return false
    !isempty(fields[2]) || return false
    tryparse_deck_float(fields[2]) === nothing && return true
    length(fields) >= 4 || return false
    secondary_control = strip(fields[3])
    isempty(secondary_control) && return true
    parsed_secondary_control = tryparse(Int, secondary_control)
    return parsed_secondary_control === 0
end

function free_field_int_or_default!(result::DeckParseResult, fields::Vector{String},
                                    index::Int, line_no::Int,
                                    field::AbstractString, default::Int)
    index > length(fields) && return default
    raw = strip(fields[index])
    isempty(raw) && return default
    value = tryparse(Int, raw)
    if value === nothing
        add_issue!(result.validation,
                   invalid_value("line $line_no",
                                 "$field=$(raw): expected free-field Int"))
        return nothing
    end
    return value
end

function free_field_float_or_default!(result::DeckParseResult, fields::Vector{String},
                                      index::Int, line_no::Int,
                                      field::AbstractString, default::Real)
    index > length(fields) && return Float64(default)
    raw = strip(fields[index])
    isempty(raw) && return Float64(default)
    value = tryparse_deck_float(raw)
    if value === nothing
        add_issue!(result.validation,
                   invalid_value("line $line_no",
                                 "$field=$(raw): expected free-field Float64"))
        return nothing
    end
    return Float64(value)
end

function bpa_fixed_source_stop_time(source_type::Int, tstop::Float64)
    if tstop == 0.0 && source_type in BPA_FIXED_INFINITE_STOP_SOURCE_TYPES
        return Inf
    end
    return tstop
end

function push_fixed_source_table_row!(
    result::DeckParseResult,
    source_name::Symbol,
    node_name::Symbol,
    node_index::Int,
    line_no::Int,
    source_layout_kind::Symbol,
    source_iform::Int,
    amplitude::Float64,
    source_time1::Float64,
    source_time2::Float64,
    source_sfreq::Float64,
    source_tstart::Float64,
    source_tstop::Float64,
)
    push!(
        result.over5a_source_rows,
        DeckOVER5ASourceRow(
            source_name,
            node_name,
            node_index,
            line_no,
            source_layout_kind,
            source_iform,
            Float64(amplitude),
            Float64(source_time1),
            Float64(source_time2),
            Float64(source_sfreq),
            Float64(source_tstart),
            Float64(source_tstop),
        ),
    )
    return result
end

struct DCSimulatorPrimaryParseState
    line_no::Int
    layout_kind::Symbol
    node::String
    control_mode::Int
    control_gain::Float64
    initial_control_value::Float64
    first_time_coefficient_s::Float64
    numerator_time_constant_s::Float64
    second_time_coefficient_s::Float64
    switch_close_delay_s::Float64
    valid::Bool
end

function dc_simulator_primary_card_candidate(line::AbstractString)::Bool
    if occursin(',', String(line))
        fields = bpa_fixed_source_free_field_fields(line)
        length(fields) >= 3 || return false
        source_type = tryparse(Int, fields[1])
        control_mode = tryparse(Int, fields[3])
        return source_type == 16 && control_mode !== nothing && control_mode in 1:3
    end
    image = fixed_image(line)
    return fixed_int_value(image, 1, 2) == 16 &&
           something(fixed_int_value(image, 9, 10), 0) in 1:3
end

function parse_dc_simulator_primary_card(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
)
    initial_issues = length(result.validation.issues)
    free_field = occursin(',', String(line))
    layout_kind = free_field ? :free_field : :fixed_field
    if free_field
        fields = bpa_fixed_source_free_field_fields(line)
        if length(fields) > 10
            add_issue!(
                result.validation,
                invalid_value(
                    "line $line_no",
                    "free-field DC-simulator primary cards accept at most ten fields",
                ),
            )
        end
        source_type = free_field_int_or_default!(
            result, fields, 1, line_no, "dc_simulator_source_type", 0,
        )
        node = length(fields) >= 2 ? strip(fields[2]) : ""
        control_mode = free_field_int_or_default!(
            result, fields, 3, line_no, "dc_simulator_control_mode", 0,
        )
        control_gain = free_field_float_or_default!(
            result, fields, 4, line_no, "dc_simulator_control_gain", 0.0,
        )
        initial_control_value = free_field_float_or_default!(
            result, fields, 5, line_no, "dc_simulator_initial_control_value", 0.0,
        )
        first_time_coefficient_s = free_field_float_or_default!(
            result, fields, 6, line_no, "dc_simulator_first_time_coefficient_s", 0.0,
        )
        numerator_time_constant_s = free_field_float_or_default!(
            result, fields, 7, line_no, "dc_simulator_numerator_time_constant_s", 0.0,
        )
        second_time_coefficient_s = free_field_float_or_default!(
            result, fields, 8, line_no, "dc_simulator_second_time_coefficient_s", 0.0,
        )
        switch_close_delay_s = free_field_float_or_default!(
            result, fields, 9, line_no, "dc_simulator_switch_close_delay_s", 0.0,
        )
        # OVER5A reads the tenth field into GUS4, then replaces it with the
        # discretized gain before any state mutation. Decode it to preserve
        # FREFLD validation even though it has no result-affecting role.
        trailing_value = free_field_float_or_default!(
            result, fields, 10, line_no, "dc_simulator_replaced_terminal_value", 0.0,
        )
        source_type == 16 || add_issue!(
            result.validation,
            invalid_value("line $line_no", "free-field DC-simulator source type must be 16"),
        )
    else
        image = fixed_image(line)
        node = fixed_field(image, 3, 8)
        control_mode = fixed_int_field!(
            result, image, line_no, 9, 10, "dc_simulator_control_mode",
        )
        control_gain = fixed_float_field!(
            result, image, line_no, 11, 20, "dc_simulator_control_gain",
        )
        initial_control_value = fixed_float_or_default!(
            result, image, line_no, 21, 30, "dc_simulator_initial_control_value", 0.0,
        )
        first_time_coefficient_s = fixed_float_or_default!(
            result, image, line_no, 31, 40, "dc_simulator_first_time_coefficient_s", 0.0,
        )
        numerator_time_constant_s = fixed_float_or_default!(
            result, image, line_no, 41, 50, "dc_simulator_numerator_time_constant_s", 0.0,
        )
        second_time_coefficient_s = fixed_float_or_default!(
            result, image, line_no, 51, 60, "dc_simulator_second_time_coefficient_s", 0.0,
        )
        switch_close_delay_s = fixed_float_or_default!(
            result, image, line_no, 61, 70, "dc_simulator_switch_close_delay_s", 0.0,
        )
        trailing_value = fixed_float_or_default!(
            result, image, line_no, 71, 80, "dc_simulator_replaced_terminal_value", 0.0,
        )
    end
    isempty(node) && add_issue!(
        result.validation,
        missing_data("line $line_no", "expected DC-simulator controlled node in columns 3-8"),
    )
    values = (
        control_gain,
        initial_control_value,
        first_time_coefficient_s,
        numerator_time_constant_s,
        second_time_coefficient_s,
        switch_close_delay_s,
        trailing_value,
    )
    valid = control_mode !== nothing && control_mode in 1:3 &&
            all(value -> value !== nothing && isfinite(value), values) &&
            length(result.validation.issues) == initial_issues
    control_mode !== nothing && !(control_mode in 1:3) && add_issue!(
        result.validation,
        invalid_value("line $line_no", "DC-simulator control mode must be 1, 2, or 3"),
    )
    return DCSimulatorPrimaryParseState(
        line_no,
        layout_kind,
        node,
        something(control_mode, 0),
        something(control_gain, NaN),
        something(initial_control_value, NaN),
        something(first_time_coefficient_s, NaN),
        something(numerator_time_constant_s, NaN),
        something(second_time_coefficient_s, NaN),
        something(switch_close_delay_s, NaN),
        valid,
    )
end

function push_dc_simulator_resistor!(
    result::DeckParseResult,
    name::String,
    from_node::String,
    to_node::String,
    resistance_ohm::Float64,
    line_no::Int;
    output_code::Int=0,
)
    from_index = node_id!(result, from_node)
    to_index = node_id!(result, to_node)
    push_element!(
        result,
        ["resistor", name, from_node, to_node],
        ConductanceBranch(from_index, to_index, inv(resistance_ohm)),
        line_no,
    )
    push_over2_branch_row!(
        result,
        length(result.elements),
        line_no,
        0,
        :dc_simulator_generated,
        :generated_topology,
        :none,
        :none,
        0,
        resistance_ohm,
        0.0,
        0.0,
        output_code,
    )
    output_code == 0 ||
        record_bpa_branch_output_requests_by_code!(result, line_no, name, output_code)
    return length(result.elements)
end

function push_dc_simulator_voltage_probe!(
    result::DeckParseResult,
    name::String,
    positive_node::String,
    negative_node::String,
    line_no::Int,
)
    positive_index = node_id!(result, positive_node)
    negative_index = node_id!(result, negative_node)
    push_element!(
        result,
        ["voltage_probe", name, positive_node, negative_node],
        ConductanceBranch(positive_index, negative_index, 0.0),
        line_no,
    )
    push_over2_branch_row!(
        result,
        length(result.elements),
        line_no,
        0,
        :dc_simulator_voltage_probe,
        :generated_measurement,
        :none,
        :none,
        0,
        Inf,
        0.0,
        0.0,
        2,
    )
    record_bpa_branch_output_requests_by_code!(result, line_no, name, 2)
    record_card!(result, :dc_simulator_terminal_voltage_probe)
    return length(result.elements)
end

function parse_dc_simulator_secondary_card!(
    result::DeckParseResult,
    primary::DCSimulatorPrimaryParseState,
    line::AbstractString,
    line_no::Int,
)::Bool
    initial_issues = length(result.validation.issues)
    free_field = occursin(',', String(line))
    secondary_layout_kind = free_field ? :free_field : :fixed_field
    if free_field
        fields = bpa_fixed_source_free_field_fields(line)
        if length(fields) > 10
            add_issue!(
                result.validation,
                invalid_value(
                    "line $line_no",
                    "free-field DC-simulator secondary cards accept at most ten fields",
                ),
            )
        end
        # OVER5A/FREFLD decodes field 1 into N3 and deliberately replaces it
        # with field 3. Keep the decode and validation order explicit.
        _replaced_record_value = free_field_int_or_default!(
            result, fields, 1, line_no, "dc_simulator_replaced_record_value", 0,
        )
        driven_node = length(fields) >= 2 ? strip(fields[2]) : ""
        _output_selector = free_field_int_or_default!(
            result, fields, 3, line_no, "dc_simulator_output_selector", 0,
        )
        scaling_value = free_field_float_or_default!(
            result, fields, 4, line_no, "dc_simulator_scaling_value", 0.0,
        )
        scaling_denominator = free_field_float_or_default!(
            result, fields, 5, line_no, "dc_simulator_scaling_denominator", 0.0,
        )
        lower_limit = free_field_float_or_default!(
            result, fields, 6, line_no, "dc_simulator_lower_limit", 0.0,
        )
        upper_limit = free_field_float_or_default!(
            result, fields, 7, line_no, "dc_simulator_upper_limit", 0.0,
        )
        isolation_resistance_ohm = free_field_float_or_default!(
            result, fields, 8, line_no, "dc_simulator_isolation_resistance_ohm", 0.0,
        )
        balancing_frequency_hz = free_field_float_or_default!(
            result, fields, 9, line_no, "dc_simulator_balancing_frequency_hz", 0.0,
        )
        branch_output_code = free_field_int_or_default!(
            result, fields, 10, line_no, "dc_simulator_branch_output_code", 0,
        )
    else
        image = fixed_image(line)
        driven_node = fixed_field(image, 3, 8)
        _output_selector = fixed_int_or_default!(
            result, image, line_no, 9, 10, "dc_simulator_output_selector", 0,
        )
        scaling_value = fixed_float_field!(
            result, image, line_no, 11, 20, "dc_simulator_scaling_value",
        )
        scaling_denominator = fixed_float_field!(
            result, image, line_no, 21, 30, "dc_simulator_scaling_denominator",
        )
        lower_limit = fixed_float_or_default!(
            result, image, line_no, 31, 40, "dc_simulator_lower_limit", 0.0,
        )
        upper_limit = fixed_float_or_default!(
            result, image, line_no, 41, 50, "dc_simulator_upper_limit", 0.0,
        )
        isolation_resistance_ohm = fixed_float_or_default!(
            result, image, line_no, 51, 60, "dc_simulator_isolation_resistance_ohm", 0.0,
        )
        balancing_frequency_hz = fixed_float_or_default!(
            result, image, line_no, 61, 70, "dc_simulator_balancing_frequency_hz", 0.0,
        )
        branch_output_code = fixed_int_or_default!(
            result, image, line_no, 80, 80, "dc_simulator_branch_output_code", 0,
        )
    end
    isempty(driven_node) && add_issue!(
        result.validation,
        missing_data("line $line_no", "expected DC-simulator driven node in columns 3-8"),
    )
    branch_output_code !== nothing && !(branch_output_code in 0:3) && add_issue!(
        result.validation,
        invalid_value(
            "line $line_no",
            "DC-simulator branch output code must be 0, 1, 2, or 3",
        ),
    )
    values = (
        scaling_value,
        scaling_denominator,
        lower_limit,
        upper_limit,
        isolation_resistance_ohm,
        balancing_frequency_hz,
    )
    if !primary.valid || !all(value -> value !== nothing && isfinite(value), values) ||
       length(result.validation.issues) != initial_issues
        return true
    end

    horizon = deck_fixed_time_horizon_options(result)
    dt_s = horizon.dt_s
    dt_s > 0.0 || begin
        add_issue!(
            result.validation,
            missing_data(
                "line $(primary.line_no)",
                "DC-simulator topology requires a positive deck time step",
            ),
        )
        return true
    end
    delta2_s = dt_s / 2.0
    denominator = 1.0 +
        (primary.first_time_coefficient_s + primary.second_time_coefficient_s) / delta2_s +
        primary.first_time_coefficient_s * primary.second_time_coefficient_s / delta2_s^2
    discrete_gain =
        (1.0 + primary.numerator_time_constant_s / delta2_s) * primary.control_gain
    scaling_denominator != 0.0 || begin
        add_issue!(result.validation, invalid_value(
            "line $line_no", "DC-simulator scaling denominator must be nonzero",
        ))
        return true
    end
    denominator != 0.0 && discrete_gain != 0.0 || begin
        add_issue!(result.validation, invalid_value(
            "line $(primary.line_no)", "DC-simulator discrete coefficients are singular",
        ))
        return true
    end
    equivalent_resistance_ohm = discrete_gain * scaling_denominator / denominator
    equivalent_resistance_ohm > 0.0 && isfinite(equivalent_resistance_ohm) || begin
        add_issue!(result.validation, invalid_value(
            "line $line_no", "DC-simulator equivalent resistance must be finite and positive",
        ))
        return true
    end
    isolation_resistance_ohm >= 0.0 || begin
        add_issue!(result.validation, invalid_value(
            "line $line_no", "DC-simulator isolation resistance must be nonnegative",
        ))
        return true
    end
    epsilon = horizon.epsilon > 0.0 ? horizon.epsilon : eps(Float64)
    effective_isolation_resistance = isolation_resistance_ohm == 0.0 ?
        sqrt(epsilon) : isolation_resistance_ohm

    primary_time1 = denominator
    primary_tstart = 4.0 * primary.first_time_coefficient_s *
        primary.second_time_coefficient_s / dt_s / denominator
    selected_control_value = primary.control_mode == 1 ? primary.initial_control_value :
        primary.control_mode == 2 ? primary_tstart : primary_time1
    if primary.control_mode == 1 &&
       !(primary_time1 <= selected_control_value <= primary_tstart)
        add_issue!(result.validation, invalid_value(
            "line $(primary.line_no)",
            "DC-simulator mode-1 initial control value must lie between its discrete limits",
        ))
        return true
    end

    simulator_index = get(result.card_counts, :dc_simulator_topology, 0) + 1
    bridge_node = "dc_simulator_$(simulator_index)_bridge"
    reference_node = "dc_simulator_$(simulator_index)_reference"
    haskey(result.node_map, Symbol(bridge_node)) ||
        node_id!(result, bridge_node)
    haskey(result.node_map, Symbol(reference_node)) ||
        node_id!(result, reference_node)
    primary_name, primary_index = bpa_fixed_source_node_reference!(
        result, primary.node, primary.line_no,
    )
    driven_name, _ = bpa_fixed_source_node_reference!(result, driven_node, line_no)
    bridge_index = result.node_map[Symbol(bridge_node)]
    reference_index = result.node_map[Symbol(reference_node)]

    switch_name = "dc_simulator_$(simulator_index)_startup_switch"
    push_element!(
        result,
        ["time_switch", switch_name, reference_node, bridge_node],
        TimeSwitch(
            reference_index,
            bridge_index;
            close_time_s = -max(primary.switch_close_delay_s, 0.0),
            open_time_s = 1.5 * dt_s,
            initially_closed = true,
        ),
        line_no,
    )
    driven_branch_name = "dc_simulator_$(simulator_index)_isolation"
    requested_output_code = something(branch_output_code, 0)
    push_dc_simulator_resistor!(
        result,
        driven_branch_name,
        String(driven_name),
        bridge_node,
        effective_isolation_resistance,
        line_no;
        output_code = requested_output_code in (1, 3) ? 1 : 0,
    )
    push_dc_simulator_resistor!(
        result,
        "dc_simulator_$(simulator_index)_equivalent",
        String(primary_name),
        reference_node,
        equivalent_resistance_ohm,
        line_no,
    )
    requested_output_code >= 2 && push_dc_simulator_voltage_probe!(
        result,
        "dc_simulator_$(simulator_index)_terminal_voltage",
        String(driven_name),
        String(primary_name),
        line_no,
    )

    balancing_frequency = balancing_frequency_hz == 0.0 ? 1.0e-3 : balancing_frequency_hz
    balancing_crest = selected_control_value / equivalent_resistance_ohm +
        max(primary.switch_close_delay_s, 0.0)
    source_rows = (
        DeckOVER5ASourceRow(
            Symbol("dc_simulator_$(simulator_index)_controller"),
            primary_name,
            -abs(primary_index),
            primary.line_no,
            :dc_simulator_primary,
            16,
            scaling_value * denominator / scaling_denominator / discrete_gain,
            primary_time1,
            primary.second_time_coefficient_s,
            primary.control_gain,
            primary_tstart,
            inv(equivalent_resistance_ohm),
        ),
        DeckOVER5ASourceRow(
            Symbol("dc_simulator_$(simulator_index)_controlled"),
            Symbol(reference_node),
            -reference_index,
            line_no,
            :dc_simulator_successor,
            primary.control_mode,
            discrete_gain,
            lower_limit,
            0.0,
            0.0,
            upper_limit,
            primary.initial_control_value,
        ),
        DeckOVER5ASourceRow(
            Symbol("dc_simulator_$(simulator_index)_balancing_input"),
            primary_name,
            -abs(primary_index),
            line_no,
            :dc_simulator_balancing_source,
            14,
            balancing_crest,
            0.0,
            0.0,
            2.0 * pi * balancing_frequency,
            0.0,
            delta2_s,
        ),
        DeckOVER5ASourceRow(
            Symbol("dc_simulator_$(simulator_index)_balancing_reference"),
            Symbol(reference_node),
            -reference_index,
            line_no,
            :dc_simulator_balancing_source,
            14,
            -balancing_crest,
            0.0,
            0.0,
            2.0 * pi * balancing_frequency,
            0.0,
            delta2_s,
        ),
    )
    append!(result.over5a_source_rows, source_rows)
    record_card!(result, primary.layout_kind)
    record_card!(result, secondary_layout_kind)
    record_card!(result, :bpa_fixed_source)
    record_card!(result, :bpa_fixed_source_type16)
    record_card!(result, :dc_simulator_primary_card)
    record_card!(result, :dc_simulator_secondary_card)
    record_card!(result, :dc_simulator_topology)
    record_card!(result, :dc_simulator_timed_switch)
    if primary.layout_kind == :free_field || secondary_layout_kind == :free_field
        record_card!(result, :bpa_fixed_source_free_field)
    end
    primary.layout_kind == :free_field &&
        record_card!(result, :dc_simulator_primary_free_field)
    secondary_layout_kind == :free_field &&
        record_card!(result, :dc_simulator_secondary_free_field)
    primary.layout_kind == secondary_layout_kind ||
        record_card!(result, :dc_simulator_mixed_field_layout)
    for _ in 1:2
        record_card!(result, :dc_simulator_generated_resistor)
        record_card!(result, :bpa_fixed_source_controlled_state)
    end
    for _ in source_rows
        record_card!(result, :bpa_fixed_source_over16_state)
        record_card!(result, :dc_simulator_source_row)
    end
    return true
end

function push_bpa_fixed_source_row!(
    result::DeckParseResult,
    source_type::Int,
    node::AbstractString,
    secondary_control::Int,
    amplitude::Float64,
    frequency::Float64,
    time1_or_phase::Float64,
    h1::Float64,
    h2::Float64,
    tstart::Float64,
    tstop::Float64,
    line_no::Int,
    initial_issues::Int;
    source_layout_count::Union{Nothing,Symbol}=nothing,
    source_layout_kind::Symbol=:fixed_field,
)::Bool
    source_tstop = bpa_fixed_source_stop_time(source_type, tstop)
    source_tstart = source_type == 14 ? tstart : max(tstart, 0.0)
    source_time1 = time1_or_phase
    source_time2 = h2
    source_sfreq = frequency
    if source_type == 13
        if h2 == time1_or_phase
            add_issue!(result.validation,
                       invalid_value("line $line_no",
                                     "OVER5A type-13 fixed-field source requires distinct H2 and TIME1 fields"))
            return true
        end
        source_sfreq = (h1 - amplitude) / (h2 - time1_or_phase)
    elseif source_type == 14
        source_time1 = time1_or_phase * (pi / 180.0)
        source_sfreq = 2.0 * pi * frequency
    end

    name = bpa_fixed_owner_name(result, "source")
    node_name, node_index = bpa_fixed_source_node_reference!(result, node, line_no)
    node_index == 0 && return true
    source_iform = bpa_fixed_source_effective_iform(result, node_index, source_type)
    source_name = Symbol(name)

    signed_routing_source =
        source_type in 1:15 || source_type == 17 || source_type >= 60
    if secondary_control != 0 && !signed_routing_source
        push_fixed_source_table_row!(
            result,
            source_name,
            node_name,
            node_index,
            line_no,
            source_layout_kind,
            source_iform,
            Float64(amplitude),
            Float64(source_time1),
            Float64(source_time2),
            Float64(source_sfreq),
            Float64(source_tstart),
            Float64(source_tstop),
        )
        record_fixed_blocker!(result, :bpa_fixed_source_blocked,
                              :bpa_fixed_source_blocked_secondary_control)
        add_issue!(result.validation,
                   unknown_field("line $line_no",
                                 "Unsupported OVER5A fixed-field source secondary/control field $secondary_control: secondary control requires coupled source-state ownership"))
        return true
    end
    current_injection_routing = secondary_control < 0

    controlled_successor =
        !isempty(result.over5a_source_rows) && last(result.over5a_source_rows).iform == 16
    controlled_model_source =
        source_type in BPA_FIXED_CONTROLLED_SOURCE_TYPES || controlled_successor
    voltbc_model_source = source_type in BPA_FIXED_VOLTBC_SOURCE_TYPES
    tacs_model_source = source_type >= 60 && !controlled_model_source
    analytic_model_source =
        source_type in BPA_FIXED_ANALYTIC_SOURCE_TYPES && !controlled_model_source &&
        !voltbc_model_source && !tacs_model_source
    source_node_value =
        current_injection_routing ? -abs(node_index) :
        secondary_control > 0 ? abs(node_index) : node_index
    if analytic_model_source
        source_element =
            source_node_value < 0 ?
            analytic_current_injection_source(
                abs(source_node_value),
                source_type,
                amplitude,
                source_time1,
                source_sfreq,
                source_tstart,
                source_tstop,
            ) :
            analytic_thevenin_source(
                abs(source_node_value),
                BPA_FIXED_SOURCE_CONDUCTANCE,
                source_type,
                amplitude,
                source_time1,
                source_sfreq,
                source_tstart,
                source_tstop,
            )
        push_element!(
            result,
            ["source", name],
            source_element,
            line_no,
        )
    end
    push_fixed_source_table_row!(
        result,
        source_name,
        node_name,
        source_node_value,
        line_no,
        source_layout_kind,
        source_iform,
        Float64(amplitude),
        Float64(source_time1),
        Float64(source_time2),
        Float64(source_sfreq),
        Float64(source_tstart),
        Float64(source_tstop),
    )
    record_fixed_card!(result, :bpa_fixed_source, Symbol("bpa_fixed_source_type", string(source_type)), initial_issues)
    if length(result.validation.issues) == initial_issues
        source_layout_count === nothing || record_card!(result, source_layout_count)
        if analytic_model_source
            record_card!(result, :source)
            record_card!(result, :bpa_fixed_source_analytic)
        elseif voltbc_model_source
            record_card!(result, :bpa_fixed_source_voltbc_state)
        elseif tacs_model_source
            record_card!(result, :bpa_fixed_source_tacs_state)
        else
            record_card!(result, :bpa_fixed_source_controlled_state)
        end
        record_card!(result, :bpa_fixed_source_over16_state)
        if source_tstart != 0.0 || isfinite(source_tstop)
            record_card!(result, :bpa_fixed_source_window)
        end
        if tstop == 0.0 && isinf(source_tstop)
            record_card!(result, :bpa_fixed_source_zero_tstop_infinite)
        end
        current_injection_routing &&
            record_card!(result, :bpa_fixed_source_negative_secondary_control)
        source_iform == source_type ||
            record_card!(result, :bpa_fixed_source_duplicate_node_incremental_iform)
    end
    return true
end

function parse_bpa_fixed_source_free_field_card!(result::DeckParseResult,
                                                 line::AbstractString,
                                                 line_no::Int)::Bool
    fields = bpa_fixed_source_free_field_fields(line)
    initial_issues = length(result.validation.issues)
    if length(fields) > 10
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "OVER5A free-field source rows accept at most ten fields"),
        )
        return true
    end
    source_type = free_field_int_or_default!(result, fields, 1, line_no, "source_type", 0)
    source_type === nothing && return true
    if source_type == 0
        add_issue!(result.validation,
                   invalid_value("line $line_no",
                                 "OVER5A free-field source type must be nonzero"))
        return true
    end
    node = length(fields) >= 2 ? strip(fields[2]) : ""
    if isempty(node)
        add_issue!(result.validation,
                   missing_data("line $line_no",
                                "expected OVER5A free-field source node in field 2"))
        return true
    end
    secondary_control = free_field_int_or_default!(
        result,
        fields,
        3,
        line_no,
        "source_secondary_control",
        0,
    )
    secondary_control === nothing && return true
    if !bpa_fixed_source_type_accepted(source_type)
        blocker, message = bpa_fixed_source_blocker(source_type)
        record_fixed_blocker!(result, :bpa_fixed_source_blocked, blocker)
        add_issue!(result.validation,
                   unknown_field("line $line_no", message))
        return true
    end
    amplitude = free_field_float_or_default!(result, fields, 4, line_no, "source_crest", 0.0)
    frequency = free_field_float_or_default!(result, fields, 5, line_no, "source_frequency", 0.0)
    time1_or_phase = free_field_float_or_default!(result, fields, 6, line_no, "source_time1_or_phase", 0.0)
    h1 = free_field_float_or_default!(result, fields, 7, line_no, "source_h1", 0.0)
    h2 = free_field_float_or_default!(result, fields, 8, line_no, "source_h2", 0.0)
    tstart = free_field_float_or_default!(result, fields, 9, line_no, "source_tstart", 0.0)
    tstop = free_field_float_or_default!(result, fields, 10, line_no, "source_tstop", 0.0)
    if amplitude === nothing || frequency === nothing || time1_or_phase === nothing ||
       h1 === nothing || h2 === nothing || tstart === nothing || tstop === nothing
        return true
    end

    return push_bpa_fixed_source_row!(
        result,
        source_type,
        node,
        secondary_control,
        amplitude,
        frequency,
        time1_or_phase,
        h1,
        h2,
        tstart,
        tstop,
        line_no,
        initial_issues;
        source_layout_count = :bpa_fixed_source_free_field,
        source_layout_kind = :free_field,
    )
end

function bpa_fixed_branch_free_field_row_candidate(line::AbstractString)::Bool
    occursin(',', String(line)) || return false
    fields = bpa_fixed_branch_free_field_fields(line)
    isempty(fields) && return false
    uppercase(fields[1]) == "NAME" && return true
    return tryparse(Int, fields[1]) !== nothing
end

function parse_bpa_fixed_branch_free_field_card!(result::DeckParseResult,
                                                 line::AbstractString,
                                                 line_no::Int)::Bool
    fields = bpa_fixed_branch_free_field_fields(line)
    initial_issues = length(result.validation.issues)
    if isempty(fields)
        return true
    elseif uppercase(fields[1]) == "NAME"
        if length(fields) < 2 || isempty(strip(fields[2]))
            add_issue!(result.validation,
                       missing_data("line $line_no",
                                    "expected OVER2 free-field branch moniker in field 2"))
            return true
        end
        enqueue_fixed_owner_name!(result, :branch, fields[2])
        record_card!(result, :fixed_field)
        record_card!(result, :bpa_fixed_branch_name)
        record_card!(result, :bpa_fixed_branch_free_field)
        return true
    elseif length(fields) > 7
        add_issue!(
            result.validation,
            invalid_value("line $line_no", "OVER2 free-field branch rows accept at most seven fields"),
        )
        return true
    end

    itype = free_field_int_or_default!(result, fields, 1, line_no, "branch_type", 0)
    itype === nothing && return true
    from_node = length(fields) >= 2 ? strip(fields[2]) : ""
    to_node = length(fields) >= 3 ? strip(fields[3]) : ""
    if isempty(from_node) || isempty(to_node)
        add_issue!(result.validation,
                   missing_data("line $line_no", "expected OVER2 free-field branch bus names in fields 2 and 3"))
        return true
    end
    if !(itype in (0, 1))
        record_fixed_blocker!(result, :bpa_fixed_branch_blocked,
                              :bpa_fixed_branch_blocked_type)
        add_issue!(result.validation,
                   unknown_field("line $line_no",
                                 "Unsupported OVER2 branch type $itype"))
        return true
    end
    branch_field_4 = length(fields) >= 4 ? strip(fields[4]) : ""
    branch_field_5 = length(fields) >= 5 ? strip(fields[5]) : ""
    if uppercase(branch_field_4) == "COPY"
        if isempty(branch_field_5)
            add_issue!(result.validation,
                       missing_data("line $line_no",
                                    "expected OVER2 free-field COPY reference name in field 5"))
            return true
        end
        branch_field_6 = length(fields) >= 6 ? strip(fields[6]) : ""
        name_reference_index = findfirst(==(Symbol(branch_field_5)), result.element_names)
        node_pair_reference =
            !isempty(branch_field_6) &&
            (length(fields) >= 7 ||
             (name_reference_index === nothing &&
              (tryparse(Int, branch_field_6) === nothing ||
               (bpa_fixed_existing_node_index(result, branch_field_5) !== nothing &&
                bpa_fixed_existing_node_index(result, branch_field_6) !== nothing))))
        if node_pair_reference
            output_code = free_field_int_or_default!(result, fields, 7, line_no, "branch_output_code", 0)
            output_code === nothing && return true
            copy_reference_index = bpa_fixed_branch_reference_by_node_pair!(
                result,
                branch_field_5,
                branch_field_6,
                line_no,
            )
            if copy_reference_index === nothing
                record_fixed_blocker!(result, :bpa_fixed_branch_blocked,
                                      :bpa_fixed_branch_missing_copy_reference)
                return true
            end
            accepted = parse_bpa_fixed_branch_copy_reference!(
                result,
                "",
                line_no,
                from_node,
                to_node,
                itype,
                copy_reference_index,
                :bpa_fixed_branch_copy_node_pair_reference,
                initial_issues;
                output_code = output_code,
                branch_layout_kind = :free_field_copy_reference,
            )
            if length(result.validation.issues) == initial_issues
                record_card!(result, :bpa_fixed_branch_free_field)
            end
            return accepted
        end
        output_code = free_field_int_or_default!(result, fields, 6, line_no, "branch_output_code", 0)
        output_code === nothing && return true
        copy_reference_index = bpa_fixed_branch_reference_by_name!(result, branch_field_5, line_no)
        if copy_reference_index === nothing
            record_fixed_blocker!(result, :bpa_fixed_branch_blocked,
                                  :bpa_fixed_branch_missing_copy_reference)
            return true
        end
        accepted = parse_bpa_fixed_branch_copy_reference!(
            result,
            "",
            line_no,
            from_node,
            to_node,
            itype,
            copy_reference_index,
            :bpa_fixed_branch_copy_name_reference,
            initial_issues;
            output_code = output_code,
            branch_layout_kind = :free_field_copy_reference,
        )
        if length(result.validation.issues) == initial_issues
            record_card!(result, :bpa_fixed_branch_free_field)
        end
        return accepted
    elseif !isempty(branch_field_4) && tryparse_deck_float(branch_field_4) === nothing
        if isempty(branch_field_5)
            add_issue!(result.validation,
                       missing_data("line $line_no",
                                    "expected OVER2 free-field branch reference node pair in fields 4 and 5"))
            return true
        end
        output_code = free_field_int_or_default!(result, fields, 6, line_no, "branch_output_code", 0)
        output_code === nothing && return true
        copy_reference_index = bpa_fixed_branch_reference_by_node_pair!(
            result,
            branch_field_4,
            branch_field_5,
            line_no,
        )
        if copy_reference_index === nothing
            record_fixed_blocker!(result, :bpa_fixed_branch_blocked,
                                  :bpa_fixed_branch_missing_copy_reference)
            return true
        end
        accepted = parse_bpa_fixed_branch_copy_reference!(
            result,
            "",
            line_no,
            from_node,
            to_node,
            itype,
            copy_reference_index,
            :bpa_fixed_branch_copy_node_pair_reference,
            initial_issues;
            output_code = output_code,
            branch_layout_kind = :free_field_copy_reference,
        )
        if length(result.validation.issues) == initial_issues
            record_card!(result, :bpa_fixed_branch_free_field)
        end
        return accepted
    end
    resistance = free_field_float_or_default!(result, fields, 4, line_no, "branch_resistance", 0.0)
    inductance = free_field_float_or_default!(result, fields, 5, line_no, "branch_inductance", 0.0)
    capacitance = free_field_float_or_default!(result, fields, 6, line_no, "branch_capacitance", 0.0)
    output_code = free_field_int_or_default!(result, fields, 7, line_no, "branch_output_code", 0)
    if resistance === nothing || inductance === nothing || capacitance === nothing
        return true
    end
    output_code === nothing && return true
    accepted = push_bpa_fixed_branch_scalar_row!(
        result,
        from_node,
        to_node,
        itype,
        resistance,
        inductance,
        capacitance,
        line_no,
        initial_issues;
        branch_layout_kind = :free_field,
        output_code = output_code,
    )
    if length(result.validation.issues) == initial_issues
        record_card!(result, :bpa_fixed_branch_free_field)
    end
    return accepted
end

function parse_bpa_fixed_branch_card!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int;
    branch_vintage_mode::Int = 0,
)::Bool
    image = fixed_image(line)
    initial_issues = length(result.validation.issues)
    if coupled_lumped_numeric_continuation_row(result, image)
        return append_coupled_lumped_continuation_row!(
            result,
            image,
            line_no,
            initial_issues,
        )
    end
    if bpa_fixed_branch_free_field_row_candidate(line)
        return parse_bpa_fixed_branch_free_field_card!(result, line, line_no)
    end
    from_node, to_node = fixed_card_branch_terminal_pair!(result, image, line_no)
    aux_1 = fixed_field(image, 15, 20)
    aux_2 = fixed_field(image, 21, 26)
    if coupled_phase_pi_numeric_continuation_row(image)
        return append_coupled_phase_pi_continuation_row!(
            result,
            image,
            line_no,
            initial_issues,
        )
    end
    if uppercase(from_node) == "BRANCH" && uppercase(to_node) == "NAME:"
        if isempty(aux_1)
            add_issue!(result.validation,
                       missing_data("line $line_no",
                                    "expected OVER2 BRANCH NAME fixed-field owner name in columns 15-20"))
            return true
        end
        enqueue_fixed_owner_name!(result, :branch, aux_1)
        record_card!(result, :fixed_field)
        record_card!(result, :bpa_fixed_branch_name)
        return true
    elseif uppercase(from_node) == "NONLIN" && uppercase(to_node) == "NAME:"
        if isempty(aux_1)
            add_issue!(
                result.validation,
                missing_data(
                    "line $line_no",
                    "expected OVER2 NONLIN NAME fixed-field owner name in columns 15-20",
                ),
            )
            return true
        end
        enqueue_fixed_owner_name!(result, :nonlinear, aux_1)
        record_card!(result, :fixed_field)
        record_card!(result, :bpa_fixed_nonlinear_name)
        return true
    end
    inline_name = nothing
    copy_reference_kind = nothing
    copy_reference_index = nothing
    copy_reference_name = nothing
    copy_reference_from_node = nothing
    copy_reference_to_node = nothing
    if single_terminal_capacitance_row(image, from_node, to_node, aux_1, aux_2)
        return parse_bpa_fixed_single_terminal_capacitance_row!(
            result,
            image,
            line_no,
            from_node,
            to_node,
            initial_issues,
        )
    elseif single_terminal_capacitance_continuation_row(image, from_node, to_node, aux_1, aux_2)
        return parse_fixed_single_terminal_capacitance_continuation_row!(
            result,
            image,
            line_no,
            from_node,
            to_node,
            aux_1,
            initial_issues,
        )
    end
    if uppercase(aux_1) == "BRANCH"
        if isempty(aux_2)
            add_issue!(result.validation,
                       missing_data("line $line_no",
                                    "expected OVER2 fixed-field inline BRANCH name in columns 21-26"))
            return true
        end
        inline_name = aux_2
    elseif uppercase(aux_1) == "COPY"
        if isempty(aux_2)
            add_issue!(result.validation,
                       missing_data("line $line_no",
                                    "expected OVER2 fixed-field COPY reference name in columns 21-26"))
            return true
        end
        copy_reference_kind = :bpa_fixed_branch_copy_name_reference
        copy_reference_name = String(aux_2)
    elseif !isempty(aux_1) || !isempty(aux_2)
        if isempty(aux_1) || isempty(aux_2)
            add_issue!(result.validation,
                       missing_data("line $line_no",
                                    "expected OVER2 fixed-field branch reference node pair in columns 15-26"))
            return true
        end
        copy_reference_kind = :bpa_fixed_branch_copy_node_pair_reference
        copy_reference_from_node = String(aux_1)
        copy_reference_to_node = String(aux_2)
    end

    itype = fixed_card_branch_type!(result, image, line_no, from_node, to_node)
    if itype === nothing
        return true
    end
    if isempty(from_node) || isempty(to_node)
        add_issue!(result.validation,
                   missing_data("line $line_no", "expected OVER2 fixed-field branch bus names in columns 3-14"))
        return true
    end
    if branch_vintage_mode != 0 && itype < 0
        return parse_fixed_card_kc_lee_untransposed_line_row!(
            result,
            image,
            line_no,
            itype,
            from_node,
            to_node,
            aux_1,
            aux_2,
            initial_issues,
        )
    end
    if fixed_card_coupled_line_descriptor(itype) !== nothing
        return parse_fixed_card_coupled_line_row!(
            result,
            image,
            line_no,
            itype,
            from_node,
            to_node,
            aux_1,
            aux_2,
        )
    end
    coupled_phase_slices = fixed_card_compact_matrix_triplet_slices(image)
    coupled_phase_numeric =
        copy_reference_kind === nothing &&
        !isempty(coupled_phase_slices.resistance) &&
        (
            itype > 1 ||
            (
                !isempty(coupled_phase_slices.capacitance) &&
                coupled_phase_slices.capacitance[1] != 0.0 &&
                (
                    coupled_phase_slices.resistance[1] != 0.0 ||
                    coupled_phase_slices.inductance[1] != 0.0
                )
            )
        )
    coupled_phase_copy =
        copy_reference_kind == :bpa_fixed_branch_copy_node_pair_reference &&
        coupled_phase_pi_reference_by_node_pair(
            result,
            copy_reference_from_node,
            copy_reference_to_node,
        ) !== nothing
    coupled_phase_copy_continuation =
        isempty(coupled_phase_slices.resistance) &&
        copy_reference_kind === nothing &&
        coupled_phase_pi_copy_continuation_expected(result, itype)
    if 1 <= itype <= 50 &&
       (coupled_phase_numeric || coupled_phase_copy || coupled_phase_copy_continuation)
        return push_bpa_fixed_coupled_phase_pi_row!(
            result,
            image,
            line_no,
            from_node,
            to_node,
            itype,
            initial_issues;
            reference_kind =
                coupled_phase_copy ? :bpa_fixed_coupled_phase_pi_copy_reference : :none,
            reference_from_node =
                coupled_phase_copy ? copy_reference_from_node : nothing,
            reference_to_node =
                coupled_phase_copy ? copy_reference_to_node : nothing,
        )
    end
    if itype == 2 && copy_reference_kind === nothing
        return parse_fixed_card_terminal_surge_impedance_row!(
            result,
            image,
            line_no,
            itype,
            from_node,
            to_node,
            initial_issues,
        )
    end
    if itype == 91
        if copy_reference_kind !== nothing
            reference_index = bpa_fixed_nonlinear_reference_index!(
                result,
                result.nonlinear_resistance_rows,
                "time-varying resistance",
                copy_reference_kind,
                copy_reference_name,
                copy_reference_from_node,
                copy_reference_to_node,
                line_no,
            )
            reference_index === nothing && return true
            return push_nonlinear_resistance_row!(
                result,
                image,
                line_no,
                from_node,
                to_node,
                itype,
                initial_issues;
                inline_name = inline_name,
                reference_kind = copy_reference_kind,
                reference_index = reference_index,
            )
        end
        return push_nonlinear_resistance_row!(
            result,
            image,
            line_no,
            from_node,
            to_node,
            itype,
            initial_issues,
            inline_name = inline_name,
        )
    end
    if itype == 92
        wide_marker = fixed_float_value(image, 59, 74)
        table_marker =
            wide_marker === nothing ? fixed_float_value(image, 39, 44) : wide_marker
        if table_marker == 5555.0 && copy_reference_kind !== nothing
            reference_index = bpa_fixed_nonlinear_reference_index!(
                result,
                result.zinc_oxide_nonlinear_rows,
                "zinc-oxide",
                copy_reference_kind,
                copy_reference_name,
                copy_reference_from_node,
                copy_reference_to_node,
                line_no,
            )
            reference_index === nothing && return true
            return push_zinc_oxide_nonlinear_row!(
                result,
                image,
                line_no,
                from_node,
                to_node,
                itype,
                initial_issues;
                inline_name = inline_name,
                reference_kind = copy_reference_kind,
                reference_index = reference_index,
            )
        end
        if table_marker == 4444.0
            if copy_reference_kind !== nothing
                reference_index = bpa_fixed_nonlinear_reference_index!(
                    result,
                    result.nonlinear_resistance_rows,
                    "piecewise resistance",
                    copy_reference_kind,
                    copy_reference_name,
                    copy_reference_from_node,
                    copy_reference_to_node,
                    line_no,
                )
                reference_index === nothing && return true
                return push_nonlinear_resistance_row!(
                    result,
                    image,
                    line_no,
                    from_node,
                    to_node,
                    itype,
                    initial_issues;
                    inline_name = inline_name,
                    reference_kind = copy_reference_kind,
                    reference_index = reference_index,
                )
            end
            return push_nonlinear_resistance_row!(
                result,
                image,
                line_no,
                from_node,
                to_node,
                itype,
                initial_issues,
                inline_name = inline_name,
            )
        end
        return push_zinc_oxide_nonlinear_row!(
            result,
            image,
            line_no,
            from_node,
            to_node,
            itype,
            initial_issues,
            inline_name = inline_name,
        )
    end
    if itype in (93, 98)
        if copy_reference_kind !== nothing
            reference_index = bpa_fixed_nonlinear_reference_index!(
                result,
                result.piecewise_nonlinear_inductor_rows,
                itype == 93 ?
                    "piecewise nonlinear inductor" :
                    "pseudo-nonlinear inductor",
                copy_reference_kind,
                copy_reference_name,
                copy_reference_from_node,
                copy_reference_to_node,
                line_no,
            )
            reference_index === nothing && return true
            return push_piecewise_nonlinear_inductor_row!(
                result,
                image,
                line_no,
                from_node,
                to_node,
                itype,
                initial_issues;
                inline_name = inline_name,
                reference_kind = copy_reference_kind,
                reference_index = reference_index,
            )
        end
        return push_piecewise_nonlinear_inductor_row!(
            result,
            image,
            line_no,
            from_node,
            to_node,
            itype,
            initial_issues;
            inline_name = inline_name,
        )
    end
    if itype == 94
        if copy_reference_kind !== nothing
            reference_index = bpa_fixed_nonlinear_reference_index!(
                result,
                result.arrester_nonlinear_rows,
                "arrester",
                copy_reference_kind,
                copy_reference_name,
                copy_reference_from_node,
                copy_reference_to_node,
                line_no,
            )
            reference_index === nothing && return true
            return push_arrester_nonlinear_row!(
                result,
                image,
                line_no,
                from_node,
                to_node,
                itype,
                initial_issues;
                inline_name = inline_name,
                reference_kind = copy_reference_kind,
                reference_index = reference_index,
            )
        end
        return push_arrester_nonlinear_row!(
            result,
            image,
            line_no,
            from_node,
            to_node,
            itype,
            initial_issues,
            inline_name = inline_name,
        )
    end
    if itype == 96
        if copy_reference_kind !== nothing
            reference_index = bpa_fixed_nonlinear_reference_index!(
                result,
                result.hysteretic_inductor_rows,
                "hysteretic inductor",
                copy_reference_kind,
                copy_reference_name,
                copy_reference_from_node,
                copy_reference_to_node,
                line_no,
            )
            reference_index === nothing && return true
            return push_hysteretic_inductor_row!(
                result,
                image,
                line_no,
                from_node,
                to_node,
                initial_issues;
                inline_name = inline_name,
                reference_kind = copy_reference_kind,
                reference_index = reference_index,
            )
        end
        return push_hysteretic_inductor_row!(
            result,
            image,
            line_no,
            from_node,
            to_node,
            initial_issues;
            inline_name = inline_name,
        )
    end
    if itype == 97
        if copy_reference_kind !== nothing
            reference_index = bpa_fixed_nonlinear_reference_index!(
                result,
                result.triggered_timed_resistance_rows,
                "triggered timed resistance",
                copy_reference_kind,
                copy_reference_name,
                copy_reference_from_node,
                copy_reference_to_node,
                line_no,
            )
            reference_index === nothing && return true
            return push_triggered_timed_resistance_row!(
                result,
                image,
                line_no,
                from_node,
                to_node,
                initial_issues;
                inline_name = inline_name,
                reference_kind = copy_reference_kind,
                reference_index = reference_index,
            )
        end
        return push_triggered_timed_resistance_row!(
            result,
            image,
            line_no,
            from_node,
            to_node,
            initial_issues;
            inline_name = inline_name,
        )
    end
    if itype == 99
        if copy_reference_kind !== nothing
            reference_index = bpa_fixed_nonlinear_reference_index!(
                result,
                result.switching_nonlinear_resistor_rows,
                "switching nonlinear resistor",
                copy_reference_kind,
                copy_reference_name,
                copy_reference_from_node,
                copy_reference_to_node,
                line_no,
            )
            reference_index === nothing && return true
            return push_switching_nonlinear_resistor_row!(
                result,
                image,
                line_no,
                from_node,
                to_node,
                initial_issues;
                inline_name = inline_name,
                reference_kind = copy_reference_kind,
                reference_index = reference_index,
            )
        end
        return push_switching_nonlinear_resistor_row!(
            result,
            image,
            line_no,
            from_node,
            to_node,
            initial_issues;
            inline_name = inline_name,
        )
    end
    if itype == 95
        record_fixed_blocker!(
            result,
            :bpa_fixed_branch_blocked,
            :bpa_fixed_branch_undefined_nonlinear_type,
        )
        add_issue!(
            result.validation,
            unknown_field(
                "line $line_no",
                "OVER2 nonlinear type 95 has no defined physical model",
            ),
        )
        return true
    end
    if !(itype in (0, 1))
        record_fixed_blocker!(result, :bpa_fixed_branch_blocked,
                              :bpa_fixed_branch_blocked_type)
        add_issue!(result.validation,
                   unknown_field("line $line_no",
                                 "Unsupported OVER2 branch type $itype"))
        return true
    end
    if copy_reference_kind !== nothing
        copy_reference_index =
            copy_reference_kind == :bpa_fixed_branch_copy_name_reference ?
            bpa_fixed_branch_reference_by_name!(result, copy_reference_name, line_no) :
            bpa_fixed_branch_reference_by_node_pair!(
                result,
                copy_reference_from_node,
                copy_reference_to_node,
                line_no,
            )
        if copy_reference_index === nothing
            record_fixed_blocker!(result, :bpa_fixed_branch_blocked,
                                  :bpa_fixed_branch_missing_copy_reference)
            return true
        end
        return parse_bpa_fixed_branch_copy_reference!(
            result,
            image,
            line_no,
            from_node,
            to_node,
            itype,
            copy_reference_index,
            copy_reference_kind,
            initial_issues,
        )
    end

    resistance, inductance, capacitance, layout_count =
        bpa_fixed_branch_triplet!(result, image, line_no)
    if resistance === nothing || inductance === nothing || capacitance === nothing
        return true
    end

    return push_bpa_fixed_branch_scalar_row!(
        result,
        from_node,
        to_node,
        itype,
        resistance,
        inductance,
        capacitance,
        line_no,
        initial_issues;
        branch_layout_count = layout_count,
        inline_name = inline_name,
        output_image = image,
    )
end

function parse_bpa_fixed_source_card!(result::DeckParseResult, line::AbstractString,
                                      line_no::Int)::Bool
    image = fixed_image(line)
    if bpa_fixed_source_card_image_candidate(image)
        return parse_bpa_fixed_source_card_image!(result, image, line_no)
    elseif bpa_fixed_source_free_field_row_candidate(line)
        return parse_bpa_fixed_source_free_field_card!(result, line, line_no)
    elseif bpa_fixed_source_card_free_field_image_candidate(line)
        return parse_bpa_fixed_source_card_free_field_image!(result, line, line_no)
    end
    initial_issues = length(result.validation.issues)
    source_type = fixed_int_field!(result, image, line_no, 1, 2, "source_type")
    node = fixed_field(image, 3, 8)
    secondary_control = fixed_int_or_default!(result, image, line_no, 9, 10,
                                              "source_secondary_control", 0)
    if source_type === nothing || secondary_control === nothing
        return true
    end
    if isempty(node)
        add_issue!(result.validation,
                   missing_data("line $line_no", "expected OVER5A fixed-field source node in columns 3-8"))
        return true
    end
    if !bpa_fixed_source_type_accepted(source_type)
        blocker, message = bpa_fixed_source_blocker(source_type)
        record_fixed_blocker!(result, :bpa_fixed_source_blocked, blocker)
        add_issue!(result.validation,
                   unknown_field("line $line_no", message))
        return true
    end
    if source_type == 18 && secondary_control == -1
        return parse_bpa_fixed_type18_ungrounded_source_card!(
            result,
            image,
            line_no,
            node,
            initial_issues,
        )
    elseif source_type == 18 && (secondary_control != 0 || fixed_float_value(image, 21, 30) === nothing)
        return parse_bpa_fixed_type18_ideal_transformer_card!(
            result,
            image,
            line_no,
            node,
            initial_issues,
        )
    end
    amplitude = fixed_float_field!(result, image, line_no, 11, 20, "source_crest")
    frequency = fixed_float_or_default!(result, image, line_no, 21, 30, "source_frequency", 0.0)
    time1_or_phase = fixed_float_or_default!(result, image, line_no, 31, 40, "source_time1_or_phase", 0.0)
    h1 = fixed_float_or_default!(result, image, line_no, 41, 50, "source_h1", 0.0)
    h2 = fixed_float_or_default!(result, image, line_no, 51, 60, "source_h2", 0.0)
    tstart = fixed_float_or_default!(result, image, line_no, 61, 70, "source_tstart", 0.0)
    tstop = fixed_float_or_default!(result, image, line_no, 71, 80, "source_tstop", Inf)
    if amplitude === nothing || frequency === nothing || time1_or_phase === nothing ||
       h1 === nothing || h2 === nothing || tstart === nothing || tstop === nothing
        return true
    end

    return push_bpa_fixed_source_row!(
        result,
        source_type,
        node,
        secondary_control,
        Float64(amplitude),
        Float64(frequency),
        Float64(time1_or_phase),
        Float64(h1),
        Float64(h2),
        Float64(tstart),
        Float64(tstop),
        line_no,
        initial_issues,
    )
end

function validate_over5a_controlled_source_rows!(result::DeckParseResult)
    rows = result.over5a_source_rows
    for (index, row) in enumerate(rows)
        if row.iform == 16 && index == length(rows)
            add_issue!(
                result.validation,
                invalid_value(
                    "line $(row.line_no)",
                    "OVER5A type-16 source row requires a successor source row",
                ),
            )
        elseif row.iform == 17 && (!isinteger(row.sfreq) || row.sfreq <= 0.0)
            add_issue!(
                result.validation,
                invalid_value(
                    "line $(row.line_no)",
                    "OVER5A type-17 source row requires positive integer XTCS index in SFREQ",
                ),
            )
        elseif row.iform == 18 && (!isinteger(row.time1) || row.time1 <= 0.0)
            add_issue!(
                result.validation,
                invalid_value(
                    "line $(row.line_no)",
                    "OVER5A type-18 source row requires positive integer F index in TIME1",
                ),
            )
        elseif abs(row.iform) >= 60 && (!isinteger(row.sfreq) || row.sfreq <= 0.0)
            add_issue!(
                result.validation,
                invalid_value(
                    "line $(row.line_no)",
                    "OVER5A type-60-plus source row requires positive integer XTCS index in SFREQ",
                ),
            )
        end
    end
    return result
end

function validate_coupled_line_rows!(result::DeckParseResult)
    rows = result.coupled_line_rows
    index = 1
    accepted_terminal_surge_index = 0
    accepted_single_phase_distributed_index = 0
    accepted_group_index = 0
    accepted_distributed_group_index = 0
    accepted_coupled_impedances = CoupledLumpedSequenceImpedance[]
    accepted_distributed_constants = DistributedTransposedLineConstants[]
    while index <= length(rows)
        row = rows[index]
        if row.sampled_frequency_data_requested
            record_card!(result, :fixed_card_sampled_frequency_line_intake)
            if row.phase_index == 1 && index + 2 <= length(rows)
                candidate = rows[index:(index + 2)]
                complete_group = all(
                    candidate_row -> candidate_row.sampled_frequency_data_requested,
                    candidate,
                ) && [candidate_row.phase_index for candidate_row in candidate] == [1, 2, 3]
                if complete_group
                    record_card!(result, :fixed_card_sampled_frequency_line_intake)
                    record_card!(result, :fixed_card_sampled_frequency_line_intake)
                    index += 3
                    continue
                end
            end
            next_is_continuation = index < length(rows) &&
                rows[index + 1].sampled_frequency_data_requested &&
                rows[index + 1].phase_index != 1
            if row.phase_index != 1 || next_is_continuation
                record_fixed_blocker!(
                    result,
                    :bpa_fixed_branch_blocked,
                    :sampled_frequency_line_multiphase_group_incomplete,
                )
                add_issue!(
                    result.validation,
                    invalid_value(
                        "line $(row.line_no)",
                        "sampled frequency-dependent line phases must be a direct phase-1 row or contiguous rows [-1, -2, -3]",
                    ),
                )
            end
            index += 1
            continue
        elseif row.line_kind == :terminal_surge_impedance
            accepted_terminal_surge_index += 1
            try
                accept_terminal_surge_impedance_row!(
                    result,
                    row,
                    accepted_terminal_surge_index,
                )
            catch err
                record_fixed_blocker!(
                    result,
                    :bpa_fixed_branch_blocked,
                    :fixed_card_terminal_surge_impedance_runtime_deferred,
                )
                add_issue!(
                    result.validation,
                    unknown_field(
                        "line $(row.line_no)",
                        "Deferred terminal surge impedance equations: $(sprint(showerror, err))",
                    ),
                )
            end
            index += 1
            continue
        elseif row.line_kind == :kc_lee_untransposed_line
            group = _kc_lee_group_at(rows, index)
            if isempty(group)
                record_fixed_blocker!(result, :bpa_fixed_branch_blocked,
                                      :fixed_card_kc_lee_untransposed_line_incomplete)
                index += 1
            else
                index += length(group)
            end
            continue
        end
        if row.line_kind == :distributed_transmission_line
            if index + 2 > length(rows)
                if row.phase_index == 1
                    accepted_single_phase_distributed_index += 1
                    try
                        accept_single_phase_distributed_line_row!(
                            result,
                            row,
                            accepted_single_phase_distributed_index,
                        )
                    catch err
                        record_fixed_blocker!(
                            result,
                            :bpa_fixed_branch_blocked,
                            :fixed_card_single_phase_distributed_line_runtime_deferred,
                        )
                        add_issue!(
                            result.validation,
                            unknown_field(
                                "line $(row.line_no)",
                                "Deferred single-phase distributed line equations: $(sprint(showerror, err))",
                            ),
                        )
                    end
                    index += 1
                    continue
                end
                record_fixed_blocker!(result, :bpa_fixed_branch_blocked,
                                      :fixed_card_distributed_transposed_line_incomplete)
                add_issue!(
                    result.validation,
                    missing_data(
                        "line $(row.line_no)",
                        "expected fixed-card distributed transposed line rows -1, -2, and -3",
                    ),
                )
                index += 1
                continue
            end
            candidate = rows[index:(index + 2)]
            if !all(candidate_row -> candidate_row.line_kind == :distributed_transmission_line,
                    candidate) ||
               [candidate_row.phase_index for candidate_row in candidate] != [1, 2, 3]
                if row.phase_index == 1 &&
                   (
                       candidate[2].line_kind != :distributed_transmission_line ||
                       candidate[2].phase_index == 1
                   )
                    accepted_single_phase_distributed_index += 1
                    try
                        accept_single_phase_distributed_line_row!(
                            result,
                            row,
                            accepted_single_phase_distributed_index,
                        )
                    catch err
                        record_fixed_blocker!(
                            result,
                            :bpa_fixed_branch_blocked,
                            :fixed_card_single_phase_distributed_line_runtime_deferred,
                        )
                        add_issue!(
                            result.validation,
                            unknown_field(
                                "line $(row.line_no)",
                                "Deferred single-phase distributed line equations: $(sprint(showerror, err))",
                            ),
                        )
                    end
                    index += 1
                    continue
                end
                record_fixed_blocker!(result, :bpa_fixed_branch_blocked,
                                      :fixed_card_distributed_transposed_line_incomplete)
                add_issue!(
                    result.validation,
                    invalid_value(
                        "line $(row.line_no)",
                        "expected contiguous fixed-card distributed transposed line phase rows [-1, -2, -3]",
                    ),
                )
                index += 1
                continue
            end
            accepted_distributed_group_index += 1
            try
                constants = _distributed_transposed_line_constants_from_rows(
                    result,
                    candidate,
                    accepted_distributed_group_index,
                    accepted_distributed_constants,
                )
                push!(accepted_distributed_constants, constants)
                modal_state = distributed_transposed_line_modal_branch_state(
                    constants;
                    name = Symbol(
                        "distributed_transposed_line_modal_branch_state_",
                        accepted_distributed_group_index,
                    ),
                )
                distributed_transposed_line_steady_state_pi_equivalent(
                    constants;
                    steady_state_frequency_hz = _deck_nominal_steady_state_frequency_hz(result),
                    storage_start_index = 1,
                    name = Symbol(
                        "distributed_transposed_line_steady_state_pi_equivalent_",
                        accepted_distributed_group_index,
                    ),
                )
                distributed_transposed_line_history_state(
                    modal_state;
                    timestep_s = deck_fixed_time_horizon_options(result).dt_s,
                    steady_state_frequency_hz = _deck_nominal_steady_state_frequency_hz(result),
                    history_storage_start_index = 1,
                    initialized_from_steady_state = false,
                    name = Symbol(
                        "distributed_transposed_line_history_state_",
                        accepted_distributed_group_index,
                    ),
                )
                admittance = distributed_transposed_line_companion_admittance(
                    modal_state;
                    name = Symbol(
                        "distributed_transposed_line_companion_admittance_",
                        accepted_distributed_group_index,
                    ),
                )
                push!(result.elements, admittance)
                push!(result.element_line_numbers, first(admittance.line_numbers))
                push!(result.element_names, admittance.name)
                record_card!(result, :fixed_card_distributed_transposed_line_constants)
                record_card!(result, :fixed_card_distributed_transposed_modal_branch_state)
                record_card!(result,
                             :fixed_card_distributed_transposed_steady_state_pi_equivalent)
                record_card!(result, :fixed_card_distributed_transposed_history_state)
                record_card!(result,
                             :fixed_card_distributed_transposed_modal_timestep_update)
                record_card!(result,
                             :fixed_card_distributed_transposed_phase_current_injection)
                record_card!(result,
                             :fixed_card_distributed_transposed_companion_admittance)
            catch err
                record_fixed_blocker!(result, :bpa_fixed_branch_blocked,
                                      :fixed_card_distributed_line_constants_deferred)
                add_issue!(
                    result.validation,
                    unknown_field(
                        "line $(row.line_no)",
                        "Deferred fixed-card distributed line constants: $(sprint(showerror, err))",
                    ),
                )
            end
            index += 3
            continue
        elseif row.line_kind != :mutual_source_equivalent
            record_fixed_blocker!(result, :bpa_fixed_branch_blocked,
                                  :fixed_card_coupled_line_equation_deferred)
            add_issue!(
                result.validation,
                unknown_field(
                    "line $(row.line_no)",
                    "Deferred fixed-card coupled-line equations for branch type $(row.line_type)",
                ),
            )
            index += 1
            continue
        end

        candidate = _coupled_lumped_sequence_group_at(rows, index)
        if isempty(candidate)
            record_fixed_blocker!(result, :bpa_fixed_branch_blocked,
                                  :fixed_card_coupled_lumped_sequence_incomplete)
            add_issue!(
                result.validation,
                invalid_value(
                    "line $(row.line_no)",
                    "expected contiguous fixed-card coupled R-L phase rows beginning with type 51",
                ),
            )
            index += 1
            continue
        end

        accepted_group_index += 1
        try
            impedance = _coupled_lumped_sequence_impedance_from_rows(
                candidate,
                accepted_group_index,
                accepted_coupled_impedances,
            )
            push!(accepted_coupled_impedances, impedance)
            record_card!(result, :fixed_card_coupled_lumped_sequence_impedance)
        catch err
            record_fixed_blocker!(result, :bpa_fixed_branch_blocked,
                                  :fixed_card_coupled_lumped_sequence_explicit_matrix_deferred)
            add_issue!(
                result.validation,
                unknown_field(
                    "line $(row.line_no)",
                    "Deferred fixed-card coupled lumped R-L explicit matrix input: $(sprint(showerror, err))",
                ),
            )
        end
        index += length(candidate)
    end
    try
        accept_kc_lee_modal_line_groups!(result)
    catch err
        record_fixed_blocker!(result, :bpa_fixed_branch_blocked,
                              :fixed_card_kc_lee_untransposed_line_runtime_deferred)
        add_issue!(
            result.validation,
            unknown_field(
                "fixed-card modal line",
                "Deferred fixed-card modal line equations: $(sprint(showerror, err))",
            ),
        )
    end
    return result
end

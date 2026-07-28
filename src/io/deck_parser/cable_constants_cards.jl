mutable struct CableConstantsParseState
    request_line_no::Int
    phase::Symbol
    case_line_no::Int
    cable_kind_code::Int
    surface_position_code::Int
    phase_count::Int
    earth_model_code::Int
    modal_output_flag::Int
    impedance_output_flag::Int
    admittance_output_flag::Int
    pipe_count::Int
    grounding_selector::Int
    layer_counts::Vector{Int}
    boundary_radii_m::Matrix{Float64}
    resistivity_ohm_m::Matrix{Float64}
    conductor_relative_permeability::Matrix{Float64}
    insulation_relative_permeability::Matrix{Float64}
    insulation_relative_permittivity::Matrix{Float64}
    depths_m::Vector{Float64}
    horizontal_positions_m::Vector{Float64}
    frequency_cards::Vector{DeckCableConstantsFrequencyCard}
    phase_index::Int
    position_index::Int
end

function CableConstantsParseState(request_line_no::Int)
    return CableConstantsParseState(
        request_line_no,
        :miscellaneous,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        Int[],
        zeros(0, 0),
        zeros(0, 0),
        zeros(0, 0),
        zeros(0, 0),
        zeros(0, 0),
        Float64[],
        Float64[],
        DeckCableConstantsFrequencyCard[],
        0,
        0,
    )
end

function _cable_constants_marker(tokens, required_words)::Bool
    isempty(tokens) && return false
    normalized_deck_token(tokens[1]) == "blank" || return false
    words = Set(normalized_deck_token(token) for token in tokens[2:end])
    return all(word -> word in words, required_words)
end

_cable_frequency_terminator(tokens) =
    _cable_constants_marker(tokens, ("frequency", "cards"))

_cable_request_terminator(tokens) =
    _cable_constants_marker(tokens, ("cable", "constants", "cases"))

function _cable_required_ints!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    count::Int,
    label::AbstractString,
)
    image = fixed_image(line)
    values = Union{Nothing,Int}[
        fixed_int_field!(
            result,
            image,
            line_no,
            5 * (index - 1) + 1,
            5 * index,
            "$(label)_$index",
        ) for index in 1:count
    ]
    any(isnothing, values) && return nothing
    return Int[value for value in values]
end

function _cable_required_floats!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    count::Int,
    label::AbstractString,
)
    image = fixed_image(line)
    values = Union{Nothing,Float64}[
        fixed_float_field!(
            result,
            image,
            line_no,
            10 * (index - 1) + 1,
            10 * index,
            "$(label)_$index",
        ) for index in 1:count
    ]
    any(isnothing, values) && return nothing
    return Float64[value for value in values]
end

function _cable_positive_values!(
    result::DeckParseResult,
    values,
    line_no::Int,
    label::AbstractString,
)
    for (index, value) in enumerate(values)
        if !isfinite(value) || value <= 0.0
            add_issue!(
                result.validation,
                invalid_value(
                    "line $line_no",
                    "$(label)_$index=$value: expected a finite positive cable value",
                ),
            )
        end
    end
    return values
end

function _reset_cable_case!(state::CableConstantsParseState)
    state.phase = :miscellaneous
    state.case_line_no = 0
    empty!(state.layer_counts)
    state.boundary_radii_m = zeros(0, 0)
    state.resistivity_ohm_m = zeros(0, 0)
    state.conductor_relative_permeability = zeros(0, 0)
    state.insulation_relative_permeability = zeros(0, 0)
    state.insulation_relative_permittivity = zeros(0, 0)
    empty!(state.depths_m)
    empty!(state.horizontal_positions_m)
    empty!(state.frequency_cards)
    state.phase_index = 0
    state.position_index = 0
    return state
end

function _parse_cable_miscellaneous!(
    result::DeckParseResult,
    state::CableConstantsParseState,
    line::AbstractString,
    line_no::Int,
)
    values = _cable_required_ints!(result, line, line_no, 9, "cable_miscellaneous")
    values === nothing && return
    state.case_line_no = line_no
    state.cable_kind_code = values[1]
    state.surface_position_code = values[2]
    state.phase_count = values[3]
    state.earth_model_code = values[4]
    state.modal_output_flag = values[5]
    state.impedance_output_flag = values[6]
    state.admittance_output_flag = values[7]
    state.pipe_count = values[8]
    state.grounding_selector = values[9]
    if state.cable_kind_code != 2
        add_issue!(
            result.validation,
            unknown_field(
                "line $line_no",
                "CABLE CONSTANTS type $(state.cable_kind_code) is not accepted; this owner currently requires type 2 concentric cables",
            ),
        )
        state.phase = :unsupported
        return
    end
    if !(1 <= state.phase_count <= 20)
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "cable phase_count=$(state.phase_count): expected 1-20",
            ),
        )
        state.phase = :unsupported
        return
    end
    state.phase = :layer_counts
    record_card!(result, :fixed_field)
    record_card!(result, :cable_constants_miscellaneous)
    return
end

function _parse_cable_layer_counts!(
    result::DeckParseResult,
    state::CableConstantsParseState,
    line::AbstractString,
    line_no::Int,
)
    values = _cable_required_ints!(
        result,
        line,
        line_no,
        state.phase_count,
        "cable_layer_count",
    )
    values === nothing && return
    if any(!=(3), values)
        add_issue!(
            result.validation,
            unknown_field(
                "line $line_no",
                "accepted type-2 CABLE CONSTANTS intake currently requires three concentric conductors per cable",
            ),
        )
    end
    state.layer_counts = values
    state.boundary_radii_m = fill(NaN, state.phase_count, 7)
    state.resistivity_ohm_m = fill(NaN, state.phase_count, 3)
    state.conductor_relative_permeability = fill(NaN, state.phase_count, 3)
    state.insulation_relative_permeability = fill(NaN, state.phase_count, 3)
    state.insulation_relative_permittivity = fill(NaN, state.phase_count, 3)
    state.depths_m = fill(NaN, state.phase_count)
    state.horizontal_positions_m = fill(NaN, state.phase_count)
    state.phase_index = 1
    state.phase = :boundary_radii
    record_card!(result, :fixed_field)
    record_card!(result, :cable_constants_layer_counts)
    return
end

function _parse_cable_boundary_radii!(
    result::DeckParseResult,
    state::CableConstantsParseState,
    line::AbstractString,
    line_no::Int,
)
    values = _cable_required_floats!(result, line, line_no, 7, "cable_boundary_radius")
    values === nothing && return
    _cable_positive_values!(result, values, line_no, "cable_boundary_radius")
    if !issorted(values) || any(diff(values) .<= 0.0)
        add_issue!(
            result.validation,
            invalid_value(
                "line $line_no",
                "cable boundary radii must be strictly increasing from core to outer insulation",
            ),
        )
    end
    state.boundary_radii_m[state.phase_index, :] .= values
    state.phase = :primary_materials
    record_card!(result, :fixed_field)
    record_card!(result, :cable_constants_boundary_radii)
    return
end

function _parse_cable_primary_materials!(
    result::DeckParseResult,
    state::CableConstantsParseState,
    line::AbstractString,
    line_no::Int,
)
    values = _cable_required_floats!(result, line, line_no, 8, "cable_material")
    values === nothing && return
    _cable_positive_values!(result, values, line_no, "cable_material")
    phase = state.phase_index
    state.resistivity_ohm_m[phase, 1] = values[1]
    state.conductor_relative_permeability[phase, 1] = values[2]
    state.insulation_relative_permeability[phase, 1] = values[3]
    state.insulation_relative_permittivity[phase, 1] = values[4]
    state.resistivity_ohm_m[phase, 2] = values[5]
    state.conductor_relative_permeability[phase, 2] = values[6]
    state.insulation_relative_permeability[phase, 2] = values[7]
    state.insulation_relative_permittivity[phase, 2] = values[8]
    state.phase = :armor_materials
    record_card!(result, :fixed_field)
    record_card!(result, :cable_constants_materials)
    return
end

function _parse_cable_armor_materials!(
    result::DeckParseResult,
    state::CableConstantsParseState,
    line::AbstractString,
    line_no::Int,
)
    values = _cable_required_floats!(result, line, line_no, 4, "cable_armor_material")
    values === nothing && return
    _cable_positive_values!(result, values, line_no, "cable_armor_material")
    phase = state.phase_index
    state.resistivity_ohm_m[phase, 3] = values[1]
    state.conductor_relative_permeability[phase, 3] = values[2]
    state.insulation_relative_permeability[phase, 3] = values[3]
    state.insulation_relative_permittivity[phase, 3] = values[4]
    if phase < state.phase_count
        state.phase_index += 1
        state.phase = :boundary_radii
    else
        state.position_index = 1
        state.phase = :positions
    end
    record_card!(result, :fixed_field)
    record_card!(result, :cable_constants_materials)
    return
end

function _parse_cable_positions!(
    result::DeckParseResult,
    state::CableConstantsParseState,
    line::AbstractString,
    line_no::Int,
)
    count = min(4, state.phase_count - state.position_index + 1)
    values = _cable_required_floats!(result, line, line_no, 2 * count, "cable_position")
    values === nothing && return
    for local_index in 1:count
        phase = state.position_index + local_index - 1
        depth = values[2 * local_index - 1]
        horizontal = values[2 * local_index]
        if !isfinite(depth) || depth <= state.boundary_radii_m[phase, 7] / 2.0
            add_issue!(
                result.validation,
                invalid_value(
                    "line $line_no",
                    "cable depth for phase $phase must be finite and exceed half its outer diameter",
                ),
            )
        end
        isfinite(horizontal) || add_issue!(
            result.validation,
            invalid_value("line $line_no", "cable horizontal position for phase $phase must be finite"),
        )
        state.depths_m[phase] = depth
        state.horizontal_positions_m[phase] = horizontal
    end
    state.position_index += count
    state.phase = state.position_index > state.phase_count ? :frequency : :positions
    record_card!(result, :fixed_field)
    record_card!(result, :cable_constants_positions)
    return
end

function _parse_cable_frequency!(
    result::DeckParseResult,
    state::CableConstantsParseState,
    line::AbstractString,
    line_no::Int,
)
    image = fixed_image(line)
    earth = fixed_float_field!(result, image, line_no, 1, 15, "cable_earth_resistivity")
    frequency = fixed_float_field!(result, image, line_no, 16, 30, "cable_start_frequency")
    decades = fixed_int_or_default!(result, image, line_no, 31, 35, "cable_frequency_decades", 0)
    points = fixed_int_or_default!(result, image, line_no, 36, 40, "cable_frequency_points", 1)
    distance = fixed_float_or_default!(result, image, line_no, 41, 48, "cable_distance", 0.0)
    output_flag = fixed_int_or_default!(result, image, line_no, 49, 58, "cable_card_output_flag", 0)
    transform_flag = fixed_int_or_default!(result, image, line_no, 59, 60, "cable_transform_flag", 0)
    any(isnothing, (earth, frequency, decades, points, distance, output_flag, transform_flag)) && return
    if !isfinite(earth) || earth <= 0.0
        add_issue!(result.validation, invalid_value("line $line_no", "cable earth resistivity must be finite and positive"))
    end
    if !isfinite(frequency) || frequency <= 0.0
        add_issue!(result.validation, invalid_value("line $line_no", "cable start frequency must be finite and positive"))
    end
    decades >= 0 || add_issue!(result.validation, invalid_value("line $line_no", "cable frequency decade count must be nonnegative"))
    points > 0 || add_issue!(result.validation, invalid_value("line $line_no", "cable frequency points per decade must be positive"))
    push!(
        state.frequency_cards,
        DeckCableConstantsFrequencyCard(
            line_no,
            earth,
            frequency,
            decades,
            points,
            distance,
            output_flag,
            transform_flag,
        ),
    )
    record_card!(result, :fixed_field)
    record_card!(result, :cable_constants_frequency_card)
    return
end

function _finish_cable_case!(result::DeckParseResult, state::CableConstantsParseState, line_no::Int)
    isempty(state.frequency_cards) && add_issue!(
        result.validation,
        missing_data("line $line_no", "expected at least one CABLE CONSTANTS frequency card"),
    )
    push!(
        result.cable_constants_cases,
        DeckCableConstantsCase(
            state.case_line_no,
            state.cable_kind_code,
            state.surface_position_code,
            state.phase_count,
            state.earth_model_code,
            state.modal_output_flag,
            state.impedance_output_flag,
            state.admittance_output_flag,
            state.pipe_count,
            state.grounding_selector,
            copy(state.layer_counts),
            copy(state.boundary_radii_m),
            copy(state.resistivity_ohm_m),
            copy(state.conductor_relative_permeability),
            copy(state.insulation_relative_permeability),
            copy(state.insulation_relative_permittivity),
            copy(state.depths_m),
            copy(state.horizontal_positions_m),
            copy(state.frequency_cards),
        ),
    )
    record_card!(result, :cable_constants_case)
    state.phase = :between_cases
    return state
end

function parse_cable_constants_card!(
    result::DeckParseResult,
    state::CableConstantsParseState,
    line::AbstractString,
    tokens,
    line_no::Int,
)::Bool
    if state.phase == :unsupported
        return _cable_request_terminator(tokens)
    elseif state.phase == :between_cases
        if _cable_request_terminator(tokens)
            record_control_card!(result, :cable_constants_section_end, tokens, line_no)
            return true
        end
        _reset_cable_case!(state)
    end

    if _cable_request_terminator(tokens)
        add_issue!(
            result.validation,
            missing_data("line $line_no", "CABLE CONSTANTS case ended before its frequency-card terminator"),
        )
        return true
    elseif _cable_frequency_terminator(tokens)
        if state.phase != :frequency
            add_issue!(
                result.validation,
                missing_data("line $line_no", "CABLE CONSTANTS frequency terminator arrived during $(state.phase) input"),
            )
        end
        _finish_cable_case!(result, state, line_no)
        record_control_card!(result, :cable_constants_frequency_section_end, tokens, line_no)
        return false
    end

    if state.phase == :miscellaneous
        _parse_cable_miscellaneous!(result, state, line, line_no)
    elseif state.phase == :layer_counts
        _parse_cable_layer_counts!(result, state, line, line_no)
    elseif state.phase == :boundary_radii
        _parse_cable_boundary_radii!(result, state, line, line_no)
    elseif state.phase == :primary_materials
        _parse_cable_primary_materials!(result, state, line, line_no)
    elseif state.phase == :armor_materials
        _parse_cable_armor_materials!(result, state, line, line_no)
    elseif state.phase == :positions
        _parse_cable_positions!(result, state, line, line_no)
    elseif state.phase == :frequency
        _parse_cable_frequency!(result, state, line, line_no)
    else
        add_issue!(result.validation, unknown_field("line $line_no", "unexpected CABLE CONSTANTS parser phase $(state.phase)"))
    end
    return false
end

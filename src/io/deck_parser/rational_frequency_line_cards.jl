mutable struct RationalLineModeCardBuilder
    header_line_no::Int
    detail_line_numbers::Vector{Int}
    from_node::Symbol
    to_node::Symbol
    from_node_index::Int
    to_node_index::Int
    mode_index::Int
    pole_reduction_threshold::Float64
    characteristic_impedance_order::Int
    characteristic_impedance_infinity_ohm::Float64
    characteristic_impedance_residues::Vector{Float64}
    characteristic_impedance_poles::Vector{Float64}
    propagation_order::Int
    travel_time_s::Float64
    propagation_residues::Vector{Float64}
    propagation_poles::Vector{Float64}
end

mutable struct RationalLineCardParseState
    phase_count::Int
    modes::Vector{DeckRationalLineModeRow}
    current_mode::Union{Nothing,RationalLineModeCardBuilder}
    stage::Symbol
    transform_row::Int
    transform_row_real_values::Vector{Float64}
    transform_row_imaginary_values::Vector{Float64}
    transform_real_values::Vector{Float64}
    transform_imaginary_values::Vector{Float64}
    transform_line_numbers::Vector{Int}
end

function _rational_line_fixed_integer(image::AbstractString, first::Int, last::Int)
    value = strip(rpad(String(image), last)[first:last])
    isempty(value) && return 0
    return tryparse(Int, value)
end

function rational_frequency_line_header_card(
    line::AbstractString;
    require_phase_count::Bool=true,
)
    image = rpad(String(line), 80)
    mode_type = _rational_line_fixed_integer(image, 1, 2)
    model_kind = _rational_line_fixed_integer(image, 53, 54)
    phase_count = _rational_line_fixed_integer(image, 55, 56)
    mode_type !== nothing && mode_type < 0 || return false
    model_kind == -2 || return false
    return !require_phase_count ||
        (mode_type == -1 && phase_count !== nothing && phase_count > 0)
end

function _rational_line_mode_builder!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
    phase_count::Int,
    expected_mode::Int,
)
    rational_frequency_line_header_card(line; require_phase_count = false) ||
        throw(ArgumentError("expected a rational frequency-dependent line mode header"))
    image = rpad(String(line), 80)
    mode_index = -something(_rational_line_fixed_integer(image, 1, 2), 0)
    mode_index == expected_mode || throw(ArgumentError(
        "rational frequency-dependent line modes must be ordered 1:$phase_count",
    ))
    declared_phase_count = something(_rational_line_fixed_integer(image, 55, 56), 0)
    if expected_mode == 1
        declared_phase_count == phase_count || throw(ArgumentError(
            "first rational line mode must declare the group phase count",
        ))
    elseif declared_phase_count != 0 && declared_phase_count != phase_count
        throw(ArgumentError("rational line mode phase counts do not match"))
    end
    from_node, from_index = _semlyen_node!(result, _semlyen_field(image, 3, 8))
    to_node, to_index = _semlyen_node!(result, _semlyen_field(image, 9, 14))
    pole_reduction_threshold = something(_semlyen_header_float(image, 33, 38), NaN)
    isfinite(pole_reduction_threshold) && pole_reduction_threshold >= 0.0 ||
        throw(ArgumentError("rational line pole-reduction threshold must be nonnegative"))
    return RationalLineModeCardBuilder(
        line_no,
        Int[],
        from_node,
        to_node,
        from_index,
        to_index,
        mode_index,
        pole_reduction_threshold,
        0,
        NaN,
        Float64[],
        Float64[],
        0,
        NaN,
        Float64[],
        Float64[],
    )
end

function start_rational_frequency_line_parse!(
    result::DeckParseResult,
    line::AbstractString,
    line_no::Int,
)
    try
        image = rpad(String(line), 80)
        phase_count = something(_rational_line_fixed_integer(image, 55, 56), 0)
        1 <= phase_count <= 18 || throw(ArgumentError(
            "rational frequency-dependent line phase count must be in 1:18",
        ))
        builder = _rational_line_mode_builder!(result, image, line_no, phase_count, 1)
        record_card!(result, :marti_frequency_dependent_line_header)
        return RationalLineCardParseState(
            phase_count,
            DeckRationalLineModeRow[],
            builder,
            :characteristic_impedance_header,
            1,
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Int[],
        )
    catch err
        add_issue!(result.validation, invalid_value("line $line_no", sprint(showerror, err)))
        return nothing
    end
end

function _rational_line_order_and_value(line::AbstractString, label::AbstractString)
    image = rpad(String(line), 40)
    order_text = strip(image[1:8])
    value_text = replace(strip(image[9:40]), 'D' => 'E', 'd' => 'e')
    order = tryparse(Int, order_text)
    value = tryparse(Float64, value_text)
    order !== nothing && order >= 0 ||
        throw(ArgumentError("$label order must be a nonnegative integer"))
    value !== nothing && isfinite(value) ||
        throw(ArgumentError("$label value must be finite"))
    return order, value
end

function _append_exact_rational_values!(
    destination::Vector{Float64},
    line::AbstractString,
    expected_count::Int,
    label::AbstractString,
)
    values = _semlyen_numeric_values(line)
    isempty(values) && throw(ArgumentError("expected numeric $label values"))
    append!(destination, values)
    length(destination) <= expected_count ||
        throw(ArgumentError("too many $label values"))
    return length(destination) == expected_count
end

function _finish_rational_line_mode!(
    result::DeckParseResult,
    state::RationalLineCardParseState,
)
    builder = something(state.current_mode)
    impedance_cards = PoleResidueTransfer(
        builder.characteristic_impedance_infinity_ohm,
        builder.characteristic_impedance_residues,
        builder.characteristic_impedance_poles,
    )
    propagation_cards = PoleResidueTransfer(
        0.0,
        builder.propagation_residues,
        builder.propagation_poles,
    )
    timestep_s = deck_fixed_time_horizon_options(result).dt_s
    impedance_reduction = pole_residue_transfer_for_timestep(
        impedance_cards,
        timestep_s,
        builder.pole_reduction_threshold,
    )
    propagation_reduction = pole_residue_transfer_for_timestep(
        propagation_cards,
        timestep_s,
        builder.pole_reduction_threshold,
    )
    conversion = rational_frequency_dependent_mode_parameters(
        impedance_reduction.response,
        propagation_reduction.response,
        builder.travel_time_s,
        _deck_nominal_steady_state_frequency_hz(result),
    )
    push!(
        state.modes,
        DeckRationalLineModeRow(
            builder.header_line_no,
            copy(builder.detail_line_numbers),
            builder.from_node,
            builder.to_node,
            builder.from_node_index,
            builder.to_node_index,
            builder.mode_index,
            conversion,
            impedance_reduction,
            propagation_reduction,
        ),
    )
    record_card!(result, :marti_frequency_dependent_line_mode)
    state.current_mode = nothing
    state.stage = length(state.modes) == state.phase_count ? :transform_real : :mode_header
    return state
end

function _rational_line_transform_matrix(state::RationalLineCardParseState)
    expected = state.phase_count^2
    length(state.transform_real_values) == expected &&
        length(state.transform_imaginary_values) == expected ||
        throw(ArgumentError("rational line complex modal transform has the wrong value count"))
    values = complex.(state.transform_real_values, state.transform_imaginary_values)
    return Matrix(transpose(reshape(values, state.phase_count, state.phase_count)))
end

function _finish_rational_line_group!(
    result::DeckParseResult,
    state::RationalLineCardParseState,
)
    voltage_transform = _rational_line_transform_matrix(state)
    current_transform = inv(transpose(voltage_transform))
    parameters = [mode.conversion.parameters for mode in state.modes]
    semlyen_line_physical_checks(
        parameters,
        voltage_transform,
        current_transform;
        inactive_phase_indices = [
            index for index in eachindex(state.modes)
            if state.modes[index].from_node_index == 0 &&
                state.modes[index].to_node_index == 0
        ],
    )
    push!(
        result.rational_frequency_line_groups,
        DeckRationalLineGroupRow(
            Symbol(
                "rational_frequency_line_",
                length(result.rational_frequency_line_groups) + 1,
            ),
            state.phase_count,
            copy(state.modes),
            voltage_transform,
            current_transform,
            copy(state.transform_line_numbers),
            :marti_pole_residue_cards,
        ),
    )
    record_card!(result, :marti_frequency_dependent_line_group)
    return nothing
end

function _parse_rational_line_transform!(
    result::DeckParseResult,
    state::RationalLineCardParseState,
    line::AbstractString,
    line_no::Int,
)
    values = _semlyen_numeric_values(line)
    isempty(values) && throw(ArgumentError("expected rational line modal-transform values"))
    destination = state.stage == :transform_real ?
        state.transform_row_real_values : state.transform_row_imaginary_values
    append!(destination, values)
    length(destination) <= state.phase_count ||
        throw(ArgumentError("too many rational line modal-transform row values"))
    push!(state.transform_line_numbers, line_no)
    length(destination) < state.phase_count && return state
    if state.stage == :transform_real
        append!(state.transform_real_values, destination)
        empty!(destination)
        state.stage = :transform_imaginary
        return state
    end
    append!(state.transform_imaginary_values, destination)
    empty!(destination)
    if state.transform_row == state.phase_count
        return _finish_rational_line_group!(result, state)
    end
    state.transform_row += 1
    state.stage = :transform_real
    return state
end

function parse_rational_frequency_line_card!(
    result::DeckParseResult,
    state::RationalLineCardParseState,
    line::AbstractString,
    line_no::Int,
)
    try
        if state.stage == :mode_header
            next_mode = length(state.modes) + 1
            state.current_mode = _rational_line_mode_builder!(
                result,
                line,
                line_no,
                state.phase_count,
                next_mode,
            )
            state.stage = :characteristic_impedance_header
            record_card!(result, :marti_frequency_dependent_line_header)
            return state
        elseif state.stage in (:transform_real, :transform_imaginary)
            return _parse_rational_line_transform!(result, state, line, line_no)
        end

        builder = something(state.current_mode)
        push!(builder.detail_line_numbers, line_no)
        if state.stage == :characteristic_impedance_header
            order, infinity_value = _rational_line_order_and_value(
                line,
                "characteristic impedance",
            )
            builder.characteristic_impedance_order = order
            builder.characteristic_impedance_infinity_ohm = infinity_value
            state.stage = order == 0 ? :propagation_header : :characteristic_impedance_residues
        elseif state.stage == :characteristic_impedance_residues
            _append_exact_rational_values!(
                builder.characteristic_impedance_residues,
                line,
                builder.characteristic_impedance_order,
                "characteristic-impedance residue",
            ) && (state.stage = :characteristic_impedance_poles)
        elseif state.stage == :characteristic_impedance_poles
            _append_exact_rational_values!(
                builder.characteristic_impedance_poles,
                line,
                builder.characteristic_impedance_order,
                "characteristic-impedance pole",
            ) && (state.stage = :propagation_header)
        elseif state.stage == :propagation_header
            order, delay = _rational_line_order_and_value(line, "propagation")
            builder.propagation_order = order
            builder.travel_time_s = delay
            state.stage = order == 0 ? :finish_mode : :propagation_residues
        elseif state.stage == :propagation_residues
            _append_exact_rational_values!(
                builder.propagation_residues,
                line,
                builder.propagation_order,
                "propagation residue",
            ) && (state.stage = :propagation_poles)
        elseif state.stage == :propagation_poles
            _append_exact_rational_values!(
                builder.propagation_poles,
                line,
                builder.propagation_order,
                "propagation pole",
            ) && (state.stage = :finish_mode)
        end
        state.stage == :finish_mode && _finish_rational_line_mode!(result, state)
        return state
    catch err
        add_issue!(result.validation, invalid_value("line $line_no", sprint(showerror, err)))
        return nothing
    end
end

function deck_rational_frequency_line_elements(result::DeckParseResult, timestep_s::Real)
    return SemlyenFrequencyDependentLine[
        semlyen_frequency_dependent_line(
            getfield.(group.modes, :from_node_index),
            getfield.(group.modes, :to_node_index),
            [mode.conversion.parameters for mode in group.modes],
            group.voltage_modal_to_phase,
            group.current_modal_to_phase,
            timestep_s,
        )
        for group in result.rational_frequency_line_groups
    ]
end

deck_rational_frequency_line_element_names(result::DeckParseResult) =
    getfield.(result.rational_frequency_line_groups, :name)

function _deck_control_system_sinusoidal_signal(
    row::DeckParser.DeckControlSystemSourceRow,
)
    row.source_type == 14 ||
        throw(ArgumentError("control source $(row.name) is not sinusoidal type 14"))
    length(row.numeric_values) >= 3 || throw(ArgumentError(
        "type-14 control source $(row.name) requires amplitude, frequency, and phase",
    ))
    start_time_s = length(row.numeric_values) >= 4 ? row.numeric_values[4] : 0.0
    stop_time_s = length(row.numeric_values) >= 5 && row.numeric_values[5] != 0.0 ?
        row.numeric_values[5] : Inf
    return SinusoidalControlSignal(
        row.name,
        row.numeric_values[1],
        row.numeric_values[2],
        deg2rad(row.numeric_values[3]);
        start_time_s = start_time_s,
        stop_time_s = stop_time_s,
    )
end

function _deck_control_system_sinusoidal_sources(
    parsed::DeckParser.DeckParseResult,
)
    return SinusoidalControlSignal[
        _deck_control_system_sinusoidal_signal(row)
        for row in DeckParser.deck_control_system_source_rows(parsed)
        if row.source_type == 14
    ]
end

function _control_system_frequency_weighted_input(
    phasors::AbstractVector{ComplexF64},
    row::OVER16CSUPDeviceRow,
)
    row.device_type in (60, 61, 63, 67) && return 0.0 + 0.0im
    return sum(
        term.xtcs_index == 0 ? 0.0 + 0.0im :
        phasors[term.xtcs_index] * term.scale
        for term in row.input_terms;
        init = 0.0 + 0.0im,
    )
end

function _control_system_frequency_device_value(
    phasors::AbstractVector{ComplexF64},
    state::OVER16CSUPState,
    row::OVER16CSUPDeviceRow,
    angular_frequency::Float64,
)
    input = _control_system_frequency_weighted_input(phasors, row)
    p1 = state.parsup_values[row.parsup_index]
    p2 = row.device_type == 66 ? 0.0 : state.parsup_values[row.parsup_index + 1]
    p3 = row.device_type in (60, 66) ? 0.0 :
        state.parsup_values[row.parsup_index + 2]
    if row.device_type in (51, 52)
        return p3 <= -1.0 ? (p1 == 0.0 ? input : p1 * input) : 0.0 + 0.0im
    elseif row.device_type == 53
        state.parsup_values[row.parsup_index + 3] == 0.0 && return 0.0 + 0.0im
        phase_shift = angular_frequency * p3 * state.parsup_values[row.parsup_index + 1]
        return input * cis(-phase_shift)
    elseif row.device_type == 59
        return p3 == 1.0 ? im * angular_frequency * input : 0.0 + 0.0im
    elseif row.device_type == 60
        selector = trunc(Int, p2)
        1 <= selector <= 3 || return 0.0 + 0.0im
        term = row.input_terms[4 - selector]
        return term.xtcs_index == 0 ? 0.0 + 0.0im :
            phasors[term.xtcs_index] * term.scale
    elseif row.device_type == 61
        selector = trunc(Int, p1)
        if 1 <= selector <= 5
            index = length(row.input_terms) - selector + 1
            1 <= index <= length(row.input_terms) || return 0.0 + 0.0im
            term = row.input_terms[index]
            return term.xtcs_index == 0 ? 0.0 + 0.0im :
                phasors[term.xtcs_index] * term.scale
        elseif selector == 6 && row.control_index > 0
            return phasors[row.control_index]
        end
        return 0.0 + 0.0im
    elseif row.device_type == 62
        return p2 >= 1.5 ? input : 0.0 + 0.0im
    end
    return 0.0 + 0.0im
end

function _control_system_frequency_initializations(
    state::ControlSystemExecutionState,
    supplemental::Union{Nothing,ControlSystemSupplementalDeviceRuntime},
    sources::Vector{SinusoidalControlSignal},
)
    frequencies = sort!(unique(
        source.frequency_hz for source in sources if source.start_time_s < 0.0
    ))
    ordinary_signal_names =
        supplemental === nothing ?
        sort!(collect(keys(state.values)); by = String) :
        supplemental.ordinary_signal_names
    ordinary_slots = Dict(
        name => index for (index, name) in enumerate(ordinary_signal_names)
    )
    initializations = ControlSystemFrequencyInitialization[]
    for frequency_hz in frequencies
        phasor_count =
            supplemental === nothing ?
            length(ordinary_signal_names) :
            length(supplemental.state.xtcs_values)
        phasors = zeros(ComplexF64, phasor_count)
        source_names = Symbol[]
        for source in sources
            source.start_time_s < 0.0 && source.frequency_hz == frequency_hz || continue
            slot = get(ordinary_slots, source.name, 0)
            slot > 0 || throw(ArgumentError(
                "sinusoidal control source $(source.name) has no runtime signal slot",
            ))
            phasors[slot] += sinusoidal_control_signal_phasor(source)
            push!(source_names, source.name)
        end
        angular_frequency = 2.0 * pi * frequency_hz
        function_count = length(state.functions)
        if function_count > 0
            matrix = Matrix{ComplexF64}(I, function_count, function_count)
            right_hand_side = zeros(ComplexF64, function_count)
            output_indices = Dict(
                row.output_name => index
                for (index, row) in enumerate(state.functions)
            )
            for (index, function_row) in enumerate(state.functions)
                frequency_variable = im * angular_frequency
                numerator = sum(
                    function_row.numerator_coefficients[order + 1] *
                    frequency_variable^order
                    for order in 0:function_row.order;
                    init = 0.0 + 0.0im,
                )
                denominator = sum(
                    function_row.denominator_coefficients[order + 1] *
                    frequency_variable^order
                    for order in 0:function_row.order;
                    init = 0.0 + 0.0im,
                )
                abs(denominator) > eps(Float64) || throw(ArgumentError(
                    "control function $(function_row.output_name) has a singular " *
                    "response at $frequency_hz Hz",
                ))
                response = function_row.gain * numerator / denominator
                for term in function_row.input_terms
                    coupled_index = get(output_indices, term.name, 0)
                    if coupled_index > 0
                        matrix[index, coupled_index] -= response * term.polarity
                    else
                        input_slot = get(ordinary_slots, term.name, 0)
                        input_slot > 0 || throw(ArgumentError(
                            "control function $(function_row.output_name) has no " *
                            "phasor input slot for $(term.name)",
                        ))
                        right_hand_side[index] +=
                            response * term.polarity * phasors[input_slot]
                    end
                end
            end
            outputs = matrix \ right_hand_side
            all(isfinite, outputs) || throw(ArgumentError(
                "control function frequency solution is nonfinite at $frequency_hz Hz",
            ))
            for index in eachindex(state.functions)
                slot = get(ordinary_slots, state.functions[index].output_name, 0)
                slot > 0 || continue
                phasors[slot] = outputs[index]
            end
        end
        if supplemental !== nothing
            for (index, row) in enumerate(supplemental.rows)
                phasors[supplemental.device_output_slots[index]] =
                    _control_system_frequency_device_value(
                        phasors,
                        supplemental.state,
                        row,
                        angular_frequency,
                    )
            end
        end
        history_mutation_count = supplemental === nothing ? 0 : sum(
            index -> supplemental.initialize_transport_delay_from_input[index] ?
                trunc(Int, supplemental.state.parsup_values[
                    supplemental.rows[index].parsup_index + 2
                ]) : 0,
            eachindex(supplemental.rows);
            init = 0,
        )
        push!(
            initializations,
            ControlSystemFrequencyInitialization(
                frequency_hz,
                supplemental === nothing ?
                    copy(ordinary_signal_names) :
                    vcat(
                        copy(supplemental.ordinary_signal_names),
                        copy(supplemental.device_output_names),
                    ),
                copy(phasors),
                source_names,
                supplemental === nothing ?
                    Symbol[] : copy(supplemental.device_output_names),
                history_mutation_count,
            ),
        )
    end
    return initializations
end

function _apply_control_system_frequency_initializations!(
    runtime::ControlSystemNetworkRuntime,
)
    supplemental = runtime.supplemental_devices
    for initialization in runtime.frequency_initializations
        angular_frequency = 2.0 * pi * initialization.frequency_hz
        ordinary_names =
            supplemental === nothing ?
            initialization.signal_names :
            supplemental.ordinary_signal_names
        ordinary_count = length(ordinary_names)
        for index in 1:ordinary_count
            name = ordinary_names[index]
            runtime.state.values[name] = get(runtime.state.values, name, 0.0) +
                                         real(initialization.signal_phasors[index])
        end
        phasor_slots = Dict(
            name => index for (index, name) in enumerate(initialization.signal_names)
        )
        for (index, function_row) in enumerate(runtime.state.functions)
            function_slot = get(phasor_slots, function_row.output_name, 0)
            function_slot > 0 || continue
            input_phasor = 0.0 + 0.0im
            for term in function_row.input_terms
                input_slot = get(phasor_slots, term.name, 0)
                input_slot > 0 || throw(ArgumentError(
                    "control function $(function_row.output_name) has no " *
                    "frequency-initialization slot for $(term.name)",
                ))
                input_phasor +=
                    term.polarity * initialization.signal_phasors[input_slot]
            end
            output_phasor = initialization.signal_phasors[function_slot]
            function_state = runtime.state.function_states[index]
            order = length(function_state.history_terms)
            order == 0 && continue
            advance = cis(angular_frequency * runtime.state.deltat_s)
            history_phasors = zeros(ComplexF64, order)
            history_phasors[end] = (
                function_state.feedforward_coefficients[end] * input_phasor -
                function_state.feedback_coefficients[end] * output_phasor
            ) / advance
            abs(advance + 1.0) > eps(Float64) || throw(ArgumentError(
                "control function $(function_row.output_name) cannot initialize " *
                "at the bilinear Nyquist frequency",
            ))
            for history_index in (order - 1):-1:1
                history_phasors[history_index] = (
                    function_state.feedforward_coefficients[history_index + 1] *
                    input_phasor -
                    function_state.feedback_coefficients[history_index + 1] *
                    output_phasor +
                    history_phasors[history_index + 1]
                ) / (advance + 1.0)
            end
            for history_index in eachindex(history_phasors)
                function_state.history_terms[history_index] +=
                    real(history_phasors[history_index])
            end
        end
        supplemental === nothing && continue
        for (index, row) in enumerate(supplemental.rows)
            input = _control_system_frequency_weighted_input(
                initialization.signal_phasors,
                row,
            )
            if supplemental.initialize_transport_delay_from_input[index]
                history_start = trunc(Int, supplemental.state.parsup_values[row.parsup_index])
                history_count = trunc(Int, supplemental.state.parsup_values[row.parsup_index + 2])
                for lag in 1:history_count
                    history_index = history_start + history_count - lag
                    supplemental.state.parsup_values[history_index] +=
                        real(input * cis(-angular_frequency * runtime.state.deltat_s * lag))
                end
                pointer = supplemental.transport_delay_pointers[index]
                supplemental.state.parsup_values[pointer] = real(input)
                pointer += 1
                pointer == history_start + history_count && (pointer = history_start)
                supplemental.transport_delay_pointers[index] = pointer
                supplemental.rows[index] = over16_csup_device_row(
                    row.supplemental_index,
                    row.device_type,
                    row.parsup_index;
                    input_terms = row.input_terms,
                    control_index = row.control_index,
                    reference_index = pointer,
                )
            elseif row.device_type == 59
                supplemental.state.parsup_values[row.parsup_index + 1] += real(input)
            end
            output_slot = supplemental.device_output_slots[index]
            supplemental.state.xtcs_values[output_slot] +=
                real(initialization.signal_phasors[output_slot])
            runtime.state.values[supplemental.device_output_names[index]] =
                supplemental.state.xtcs_values[output_slot]
        end
    end
    return runtime
end

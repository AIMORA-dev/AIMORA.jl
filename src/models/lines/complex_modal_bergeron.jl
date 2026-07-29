export ComplexModalBergeronLine,
       complex_modal_bergeron_phase_admittance,
       complex_modal_bergeron_power_invariance_error,
       complex_modal_bergeron_steady_state_terminal_admittance,
       initialize_complex_modal_bergeron_steady_state!,
       complex_modal_bergeron_terminal_currents,
       complex_modal_bergeron_history_currents

"""
    ComplexModalBergeronLine(from_nodes, to_nodes, transform,
                             characteristic_admittance, travel_time_s, dt_s)

A multiphase Bergeron element for general complex K.C. Lee modal transforms.
Phase voltages are transformed with `transform.phase_to_modal`; modal currents
return through its transpose, which is the power-invariant current transform.
Each mode owns an independent complex traveling-wave history. The nodal stamp
is the real phase-domain projection, while the complex modal histories retain
the quadrature state needed by nonsymmetric cable and untransposed-line modes.
"""
mutable struct ComplexModalBergeronLine <: EMTElement
    from_nodes::Vector{Int}
    to_nodes::Vector{Int}
    transform::LineModalTransform
    characteristic_admittance::Vector{ComplexF64}
    travel_time_s::Vector{Float64}
    attenuation::Vector{Float64}
    delay_steps::Vector{Int}
    delay_interpolation_factors::Vector{Float64}
    from_wave_history::Vector{Vector{ComplexF64}}
    to_wave_history::Vector{Vector{ComplexF64}}
    write_indices::Vector{Int}
    phase_admittance_complex::Matrix{ComplexF64}
    phase_admittance::Matrix{Float64}
    modal_voltage_from::Vector{ComplexF64}
    modal_voltage_to::Vector{ComplexF64}
    modal_current_from::Vector{ComplexF64}
    modal_current_to::Vector{ComplexF64}
    modal_history_from::Vector{ComplexF64}
    modal_history_to::Vector{ComplexF64}
    phase_current_from::Vector{Float64}
    phase_current_to::Vector{Float64}
    phase_history_from::Vector{Float64}
    phase_history_to::Vector{Float64}
    power_invariance_error::Float64
    imaginary_stamp_residual::Float64
    dt_s::Float64
end

function _complex_modal_checked_nodes(
    values::AbstractVector{<:Integer},
    phase_count::Int,
    label::AbstractString,
)
    length(values) == phase_count ||
        throw(ArgumentError("$label must have one node per phase"))
    nodes = Int.(values)
    all(>=(0), nodes) || throw(ArgumentError("$label entries must be nonnegative"))
    return nodes
end

function _complex_modal_checked_vector(
    values::AbstractVector,
    mode_count::Int,
    label::AbstractString,
)
    length(values) == mode_count ||
        throw(ArgumentError("$label must have one entry per mode"))
    result = ComplexF64.(values)
    all(value -> isfinite(real(value)) && isfinite(imag(value)), result) ||
        throw(ArgumentError("$label entries must be finite"))
    return result
end

function _complex_modal_positive_vector(
    values::AbstractVector{<:Real},
    mode_count::Int,
    label::AbstractString,
)
    length(values) == mode_count ||
        throw(ArgumentError("$label must have one entry per mode"))
    result = Float64.(values)
    all(value -> isfinite(value) && value > 0.0, result) ||
        throw(ArgumentError("$label entries must be finite and positive"))
    return result
end

function ComplexModalBergeronLine(
    from_nodes::AbstractVector{<:Integer},
    to_nodes::AbstractVector{<:Integer},
    transform::LineModalTransform,
    characteristic_admittance::AbstractVector,
    travel_time_s::AbstractVector{<:Real},
    dt_s::Real;
    attenuation::AbstractVector{<:Real}=ones(length(characteristic_admittance)),
    initial_from_wave::AbstractVector=zeros(ComplexF64, length(characteristic_admittance)),
    initial_to_wave::AbstractVector=zeros(ComplexF64, length(characteristic_admittance)),
)
    phase_to_modal = transform.phase_to_modal
    modal_to_phase = transform.modal_to_phase
    mode_count = size(phase_to_modal, 1)
    size(phase_to_modal) == (mode_count, mode_count) &&
        size(modal_to_phase) == (mode_count, mode_count) ||
        throw(ArgumentError("complex modal transform matrices must be equally square"))
    from = _complex_modal_checked_nodes(from_nodes, mode_count, "from_nodes")
    to = _complex_modal_checked_nodes(to_nodes, mode_count, "to_nodes")
    admittance = _complex_modal_checked_vector(
        characteristic_admittance,
        mode_count,
        "characteristic_admittance",
    )
    all(value -> real(value) > 0.0, admittance) ||
        throw(ArgumentError("modal characteristic admittance must have positive real parts"))
    travel = _complex_modal_positive_vector(travel_time_s, mode_count, "travel_time_s")
    dt = Float64(dt_s)
    isfinite(dt) && dt > 0.0 || throw(ArgumentError("dt_s must be finite and positive"))
    atten = Float64.(attenuation)
    length(atten) == mode_count &&
        all(value -> isfinite(value) && 0.0 <= value <= 1.0, atten) ||
        throw(ArgumentError("attenuation entries must be finite and between zero and one"))
    initial_from = _complex_modal_checked_vector(
        initial_from_wave,
        mode_count,
        "initial_from_wave",
    )
    initial_to = _complex_modal_checked_vector(
        initial_to_wave,
        mode_count,
        "initial_to_wave",
    )
    inverse_error = maximum(abs.(phase_to_modal * modal_to_phase - I))
    inverse_error <= 1.0e-9 ||
        throw(ArgumentError("complex modal voltage transforms must be mutual inverses"))
    current_to_phase = transpose(phase_to_modal)
    power_error = maximum(abs.(transpose(modal_to_phase) * current_to_phase - I))
    power_error <= 1.0e-9 ||
        throw(ArgumentError("complex modal current transform is not power invariant"))
    phase_admittance_complex =
        current_to_phase * Diagonal(admittance) * phase_to_modal
    imaginary_residual = maximum(abs, imag.(phase_admittance_complex); init=0.0)
    real_scale = maximum(abs, real.(phase_admittance_complex); init=1.0)
    imaginary_residual <= 1.0e-9 * real_scale || throw(ArgumentError(
        "complex modal characteristic admittance does not form a real phase-domain nodal stamp",
    ))
    phase_admittance = Matrix{Float64}(real.(phase_admittance_complex))
    symmetry_error = maximum(abs, phase_admittance - transpose(phase_admittance); init=0.0)
    symmetry_error <= 1.0e-9 * max(maximum(abs, phase_admittance; init=0.0), 1.0) ||
        throw(ArgumentError("complex modal phase-domain nodal stamp must be reciprocal"))
    delay_ratios = travel ./ dt
    all(>=(1.0), delay_ratios) ||
        throw(ArgumentError("travel_time_s must be at least one timestep"))
    interpolation_factors = delay_ratios .- floor.(delay_ratios)
    delays = Int.(floor.(delay_ratios)) .+ 1
    return ComplexModalBergeronLine(
        from,
        to,
        transform,
        admittance,
        travel,
        atten,
        delays,
        interpolation_factors,
        [fill(initial_from[index], delays[index]) for index in 1:mode_count],
        [fill(initial_to[index], delays[index]) for index in 1:mode_count],
        ones(Int, mode_count),
        phase_admittance_complex,
        phase_admittance,
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(ComplexF64, mode_count),
        zeros(Float64, mode_count),
        zeros(Float64, mode_count),
        zeros(Float64, mode_count),
        zeros(Float64, mode_count),
        power_error,
        imaginary_residual,
        dt,
    )
end

complex_modal_bergeron_phase_admittance(line::ComplexModalBergeronLine) =
    copy(line.phase_admittance)

complex_modal_bergeron_power_invariance_error(line::ComplexModalBergeronLine) =
    line.power_invariance_error

function complex_modal_bergeron_steady_state_terminal_admittance(
    line::ComplexModalBergeronLine,
    frequency_hz::Real,
)
    frequency = Float64(frequency_hz)
    isfinite(frequency) && frequency > 0.0 ||
        throw(ArgumentError("frequency_hz must be finite and positive"))
    mode_count = length(line.characteristic_admittance)
    modal_self = zeros(ComplexF64, mode_count)
    modal_transfer = zeros(ComplexF64, mode_count)
    for mode in 1:mode_count
        propagation = line.attenuation[mode] *
            cis(-2.0 * pi * frequency * line.travel_time_s[mode])
        denominator = 1.0 - propagation^2
        abs(denominator) > 64.0 * eps(Float64) || throw(ArgumentError(
            "complex modal Bergeron steady-state terminal admittance is resonant",
        ))
        modal_self[mode] = line.characteristic_admittance[mode] *
            (1.0 + propagation^2) / denominator
        modal_transfer[mode] = -2.0 * line.characteristic_admittance[mode] *
            propagation / denominator
    end
    phase_to_modal = line.transform.phase_to_modal
    current_to_phase = transpose(phase_to_modal)
    self = current_to_phase * Diagonal(modal_self) * phase_to_modal
    transfer = current_to_phase * Diagonal(modal_transfer) * phase_to_modal
    return [self transfer; transfer self]
end

"""
    initialize_complex_modal_bergeron_steady_state!(
        line, from_voltage_phasors, to_voltage_phasors, frequency_hz)

Seed every modal delay slot from a phase-domain sinusoidal steady state. The
phase terminal currents come from the same frequency-domain terminal
admittance used by the initial network solve, while the stored modal voltage
and current samples retain the complex K.C. Lee transform at each historical
instant.
"""
function initialize_complex_modal_bergeron_steady_state!(
    line::ComplexModalBergeronLine,
    from_voltage_phasors::AbstractVector{<:Complex},
    to_voltage_phasors::AbstractVector{<:Complex},
    frequency_hz::Real,
)
    mode_count = length(line.characteristic_admittance)
    length(from_voltage_phasors) == mode_count ||
        throw(ArgumentError("from_voltage_phasors must have one value per phase"))
    length(to_voltage_phasors) == mode_count ||
        throw(ArgumentError("to_voltage_phasors must have one value per phase"))
    frequency = Float64(frequency_hz)
    isfinite(frequency) && frequency > 0.0 ||
        throw(ArgumentError("frequency_hz must be finite and positive"))
    from_phasors = ComplexF64.(from_voltage_phasors)
    to_phasors = ComplexF64.(to_voltage_phasors)
    all(value -> isfinite(real(value)) && isfinite(imag(value)), from_phasors) &&
        all(value -> isfinite(real(value)) && isfinite(imag(value)), to_phasors) ||
        throw(ArgumentError("terminal voltage phasors must be finite"))
    terminal_admittance =
        complex_modal_bergeron_steady_state_terminal_admittance(line, frequency)
    terminal_voltages = vcat(from_phasors, to_phasors)
    terminal_currents = terminal_admittance * terminal_voltages
    from_current_phasors = view(terminal_currents, 1:mode_count)
    to_current_phasors = view(terminal_currents, (mode_count + 1):(2 * mode_count))
    phase_to_modal = line.transform.phase_to_modal
    phase_current_to_modal = transpose(line.transform.modal_to_phase)
    angular_frequency = 2.0 * pi * frequency
    for mode in 1:mode_count
        delay = line.delay_steps[mode]
        for slot in 1:delay
            historical_time_s = -(delay - slot) * line.dt_s
            rotation = cis(angular_frequency * historical_time_s)
            from_phase_voltage =
                complex.(real.(from_phasors .* rotation), 0.0)
            to_phase_voltage =
                complex.(real.(to_phasors .* rotation), 0.0)
            from_phase_current =
                complex.(real.(from_current_phasors .* rotation), 0.0)
            to_phase_current =
                complex.(real.(to_current_phasors .* rotation), 0.0)
            from_modal_voltage = phase_to_modal * from_phase_voltage
            to_modal_voltage = phase_to_modal * to_phase_voltage
            from_modal_current = phase_current_to_modal * from_phase_current
            to_modal_current = phase_current_to_modal * to_phase_current
            line.from_wave_history[mode][slot] =
                line.characteristic_admittance[mode] * from_modal_voltage[mode] +
                from_modal_current[mode]
            line.to_wave_history[mode][slot] =
                line.characteristic_admittance[mode] * to_modal_voltage[mode] +
                to_modal_current[mode]
        end
    end
    line.modal_voltage_from .= phase_to_modal * complex.(real.(from_phasors), 0.0)
    line.modal_voltage_to .= phase_to_modal * complex.(real.(to_phasors), 0.0)
    line.modal_current_from .=
        phase_current_to_modal * complex.(real.(from_current_phasors), 0.0)
    line.modal_current_to .=
        phase_current_to_modal * complex.(real.(to_current_phasors), 0.0)
    line.phase_current_from .= real.(from_current_phasors)
    line.phase_current_to .= real.(to_current_phasors)
    fill!(line.write_indices, 1)
    _complex_modal_history_phase_currents!(line)
    return line
end

complex_modal_bergeron_terminal_currents(line::ComplexModalBergeronLine) = (
    from = copy(line.phase_current_from),
    to = copy(line.phase_current_to),
)

complex_modal_bergeron_history_currents(line::ComplexModalBergeronLine) = (
    from = copy(line.phase_history_from),
    to = copy(line.phase_history_to),
)

function _complex_modal_phase_values!(
    destination::Vector{ComplexF64},
    nodes::Vector{Int},
    voltage::AbstractVector{Float64},
)
    @inbounds for index in eachindex(nodes)
        node = nodes[index]
        destination[index] = complex(node == 0 ? 0.0 : voltage[node], 0.0)
    end
    return destination
end

function _complex_modal_history_phase_currents!(line::ComplexModalBergeronLine)
    @inbounds for mode in eachindex(line.write_indices)
        cursor = line.write_indices[mode]
        next_cursor = cursor == line.delay_steps[mode] ? 1 : cursor + 1
        factor = line.delay_interpolation_factors[mode]
        delayed_to_wave =
            factor * line.to_wave_history[mode][cursor] +
            (1.0 - factor) * line.to_wave_history[mode][next_cursor]
        delayed_from_wave =
            factor * line.from_wave_history[mode][cursor] +
            (1.0 - factor) * line.from_wave_history[mode][next_cursor]
        line.modal_history_from[mode] =
            -line.attenuation[mode] * delayed_to_wave
        line.modal_history_to[mode] =
            -line.attenuation[mode] * delayed_from_wave
    end
    current_to_phase = transpose(line.transform.phase_to_modal)
    line.phase_history_from .= real.(current_to_phase * line.modal_history_from)
    line.phase_history_to .= real.(current_to_phase * line.modal_history_to)
    return line
end

function _stamp_complex_modal_phase_admittance!(
    y::AbstractMatrix{Float64},
    nodes::Vector{Int},
    admittance::Matrix{Float64},
)
    for column in eachindex(nodes), row in eachindex(nodes)
        row_node = nodes[row]
        column_node = nodes[column]
        row_node == 0 && continue
        column_node == 0 && continue
        y[row_node, column_node] += admittance[row, column]
    end
    return y
end

function stamp!(
    y::AbstractMatrix{Float64},
    rhs::AbstractVector{Float64},
    line::ComplexModalBergeronLine,
    _t::Float64,
    dt::Float64,
)
    dt == line.dt_s ||
        throw(ArgumentError("complex modal Bergeron runtime timestep must match construction"))
    _complex_modal_history_phase_currents!(line)
    _stamp_complex_modal_phase_admittance!(y, line.from_nodes, line.phase_admittance)
    _stamp_complex_modal_phase_admittance!(y, line.to_nodes, line.phase_admittance)
    @inbounds for phase in eachindex(line.from_nodes)
        stamp_history_current!(
            rhs,
            line.from_nodes[phase],
            0,
            line.phase_history_from[phase],
        )
        stamp_history_current!(
            rhs,
            line.to_nodes[phase],
            0,
            line.phase_history_to[phase],
        )
    end
    return nothing
end

function update!(
    line::ComplexModalBergeronLine,
    voltage::AbstractVector{Float64},
    dt::Float64,
)
    dt == line.dt_s ||
        throw(ArgumentError("complex modal Bergeron runtime timestep must match construction"))
    phase_from = _complex_modal_phase_values!(
        line.modal_current_from,
        line.from_nodes,
        voltage,
    )
    phase_to = _complex_modal_phase_values!(
        line.modal_current_to,
        line.to_nodes,
        voltage,
    )
    mul!(line.modal_voltage_from, line.transform.phase_to_modal, phase_from)
    mul!(line.modal_voltage_to, line.transform.phase_to_modal, phase_to)
    @inbounds for mode in eachindex(line.characteristic_admittance)
        line.modal_current_from[mode] =
            line.characteristic_admittance[mode] * line.modal_voltage_from[mode] +
            line.modal_history_from[mode]
        line.modal_current_to[mode] =
            line.characteristic_admittance[mode] * line.modal_voltage_to[mode] +
            line.modal_history_to[mode]
    end
    current_to_phase = transpose(line.transform.phase_to_modal)
    line.phase_current_from .= real.(current_to_phase * line.modal_current_from)
    line.phase_current_to .= real.(current_to_phase * line.modal_current_to)
    @inbounds for mode in eachindex(line.write_indices)
        cursor = line.write_indices[mode]
        line.from_wave_history[mode][cursor] =
            line.characteristic_admittance[mode] * line.modal_voltage_from[mode] +
            line.modal_current_from[mode]
        line.to_wave_history[mode][cursor] =
            line.characteristic_admittance[mode] * line.modal_voltage_to[mode] +
            line.modal_current_to[mode]
        line.write_indices[mode] =
            cursor == line.delay_steps[mode] ? 1 : cursor + 1
    end
    return nothing
end

export CascadedPiSeriesImpedance,
       CascadedPiShuntImpedance,
       CascadedPhasePiBlock,
       CascadedPhasePiEquivalent,
       cascaded_pi_series_impedance,
       cascaded_pi_shunt_impedance,
       cascaded_phase_pi_block,
       cascaded_phase_pi_equivalent

struct CascadedPiSeriesImpedance
    phase_index::Int
    resistance_ohm::Float64
    inductance_h::Float64
    capacitance_f::Float64
    open_circuit::Bool
end

struct CascadedPiShuntImpedance
    from_terminal::Int
    to_terminal::Int
    resistance_ohm::Float64
    inductance_h::Float64
    capacitance_f::Float64
end

struct CascadedPhasePiBlock
    section::CoupledLumpedPhasePiSection
    multiplicity::Int
    section_scale::Float64
    series_impedances::Vector{CascadedPiSeriesImpedance}
    shunt_impedances::Vector{CascadedPiShuntImpedance}
end

function _cascade_impedance_values(
    resistance_ohm::Real,
    inductance_h::Real,
    capacitance_f::Real,
)
    resistance = Float64(resistance_ohm)
    inductance = Float64(inductance_h)
    capacitance = Float64(capacitance_f)
    all(isfinite, (resistance, inductance, capacitance)) ||
        throw(ArgumentError("cascade impedance values must be finite"))
    resistance >= 0.0 && inductance >= 0.0 && capacitance >= 0.0 ||
        throw(ArgumentError("cascade impedance values must be nonnegative"))
    return resistance, inductance, capacitance
end

function cascaded_pi_series_impedance(
    phase_index::Integer,
    resistance_ohm::Real,
    inductance_h::Real,
    capacitance_f::Real;
    open_circuit::Bool=false,
)
    phase = Int(phase_index)
    phase > 0 || throw(ArgumentError("cascade series phase index must be positive"))
    resistance, inductance, capacitance = _cascade_impedance_values(
        resistance_ohm,
        inductance_h,
        capacitance_f,
    )
    open_circuit || resistance > 0.0 || inductance > 0.0 || capacitance > 0.0 ||
        throw(ArgumentError("a zero cascade series impedance is represented by node continuity"))
    return CascadedPiSeriesImpedance(
        phase,
        resistance,
        inductance,
        capacitance,
        open_circuit,
    )
end

function cascaded_pi_shunt_impedance(
    from_terminal::Integer,
    to_terminal::Integer,
    resistance_ohm::Real,
    inductance_h::Real,
    capacitance_f::Real,
)
    from_index = Int(from_terminal)
    to_index = Int(to_terminal)
    from_index != to_index ||
        throw(ArgumentError("cascade shunt terminals must be distinct"))
    resistance, inductance, capacitance = _cascade_impedance_values(
        resistance_ohm,
        inductance_h,
        capacitance_f,
    )
    resistance > 0.0 || inductance > 0.0 || capacitance > 0.0 ||
        throw(ArgumentError("cascade shunt impedance cannot be zero"))
    return CascadedPiShuntImpedance(
        from_index,
        to_index,
        resistance,
        inductance,
        capacitance,
    )
end

function cascaded_phase_pi_block(
    section::CoupledLumpedPhasePiSection,
    multiplicity::Integer,
    section_scale::Real;
    series_impedances::AbstractVector{CascadedPiSeriesImpedance}=CascadedPiSeriesImpedance[],
    shunt_impedances::AbstractVector{CascadedPiShuntImpedance}=CascadedPiShuntImpedance[],
)
    count = Int(multiplicity)
    count > 0 || throw(ArgumentError("cascade block multiplicity must be positive"))
    scale = Float64(section_scale)
    isfinite(scale) && scale > 0.0 ||
        throw(ArgumentError("cascade block section scale must be finite and positive"))
    phases = Int[modifier.phase_index for modifier in series_impedances]
    all(phase -> 1 <= phase <= section.phase_count, phases) ||
        throw(ArgumentError("cascade series phase index exceeds the section phase count"))
    length(unique(phases)) == length(phases) ||
        throw(ArgumentError("cascade block has duplicate series phase impedances"))
    for branch in shunt_impedances
        abs(branch.from_terminal) <= section.phase_count &&
            abs(branch.to_terminal) <= section.phase_count ||
            throw(ArgumentError("cascade shunt terminal exceeds the section phase count"))
    end
    return CascadedPhasePiBlock(
        section,
        count,
        scale,
        collect(series_impedances),
        collect(shunt_impedances),
    )
end

struct CascadedPhasePiEquivalent
    name::Symbol
    phase_count::Int
    section_count::Int
    section_scale::Float64
    frequency_hz::Float64
    from_nodes::Vector{Symbol}
    to_nodes::Vector{Symbol}
    from_node_indices::Vector{Int}
    to_node_indices::Vector{Int}
    terminal_admittance_s::Matrix{ComplexF64}
    internal_node_count::Int
    symmetry_max_abs_error::Float64
    minimum_real_power_eigenvalue_s::Float64
    blocks::Vector{CascadedPhasePiBlock}
    series_modifier_count::Int
    shunt_modifier_count::Int
end

function _cascade_series_admittance(
    resistance::Float64,
    inductance::Float64,
    capacitance::Float64,
    angular_frequency::Float64,
)
    impedance = ComplexF64(resistance, angular_frequency * inductance)
    capacitance > 0.0 && (impedance += inv(im * angular_frequency * capacitance))
    abs(impedance) > 0.0 || throw(ArgumentError("cascade branch impedance is singular"))
    return inv(impedance)
end

function _push_cascade_admittance!(
    stamps::Vector{Tuple{Int,Int,Int,Int,ComplexF64}},
    row_from::Int,
    row_to::Int,
    column_from::Int,
    column_to::Int,
    value::ComplexF64,
)
    push!(stamps, (row_from, row_to, column_from, column_to, value))
    return stamps
end

function _push_cascade_section!(
    stamps::Vector{Tuple{Int,Int,Int,Int,ComplexF64}},
    from_nodes::Vector{Int},
    to_nodes::Vector{Int},
    block::CascadedPhasePiBlock,
    angular_frequency::Float64,
)
    section = block.section
    series_impedance = block.section_scale .* ComplexF64.(
        section.phase_resistance_matrix,
        angular_frequency .* section.phase_inductance_matrix,
    )
    series_admittance = inv(series_impedance)
    shunt_admittance =
        (0.5im * angular_frequency * block.section_scale) .* section.phase_capacitance_matrix
    for column in eachindex(from_nodes), row in eachindex(from_nodes)
        _push_cascade_admittance!(
            stamps,
            from_nodes[row],
            to_nodes[row],
            from_nodes[column],
            to_nodes[column],
            series_admittance[row, column],
        )
        _push_cascade_admittance!(
            stamps,
            from_nodes[row],
            0,
            from_nodes[column],
            0,
            shunt_admittance[row, column],
        )
        _push_cascade_admittance!(
            stamps,
            to_nodes[row],
            0,
            to_nodes[column],
            0,
            shunt_admittance[row, column],
        )
    end
    return stamps
end

function _stamp_complex_admittance_entry!(
    admittance::Matrix{ComplexF64},
    row_from::Int,
    row_to::Int,
    column_from::Int,
    column_to::Int,
    value::ComplexF64,
)
    row_from != 0 && column_from != 0 && (admittance[row_from, column_from] += value)
    row_from != 0 && column_to != 0 && (admittance[row_from, column_to] -= value)
    row_to != 0 && column_from != 0 && (admittance[row_to, column_from] -= value)
    row_to != 0 && column_to != 0 && (admittance[row_to, column_to] += value)
    return admittance
end

function _cascade_shunt_node(
    terminal::Int,
    station::Vector{Int},
    auxiliary_nodes::Dict{Int,Int},
)
    terminal == 0 && return 0
    terminal > 0 && return station[terminal]
    return auxiliary_nodes[-terminal]
end

function cascaded_phase_pi_equivalent(
    name::Symbol,
    blocks::AbstractVector{CascadedPhasePiBlock},
    frequency_hz::Real,
)
    isempty(blocks) && throw(ArgumentError("cascaded PI requires at least one block"))
    frequency = Float64(frequency_hz)
    isfinite(frequency) && frequency > 0.0 ||
        throw(ArgumentError("cascaded PI frequency_hz must be finite and positive"))
    phase_count = first(blocks).section.phase_count
    all(block -> block.section.phase_count == phase_count, blocks) ||
        throw(ArgumentError("all cascade blocks must have the same phase count"))
    isempty(first(blocks).series_impedances) && isempty(first(blocks).shunt_impedances) ||
        throw(ArgumentError("the first cascade block cannot precede itself with modifiers"))

    angular_frequency = 2.0 * pi * frequency
    stamps = Tuple{Int,Int,Int,Int,ComplexF64}[]
    node_count = phase_count
    initial_station = collect(1:phase_count)
    current_station = copy(initial_station)
    first_section = true
    for block in blocks
        for _ in 1:block.multiplicity
            section_from = current_station
            if !first_section
                series_by_phase = Dict(
                    modifier.phase_index => modifier for modifier in block.series_impedances
                )
                if !isempty(series_by_phase)
                    section_from = copy(current_station)
                    for phase in 1:phase_count
                        haskey(series_by_phase, phase) || continue
                        modifier = series_by_phase[phase]
                        node_count += 1
                        section_from[phase] = node_count
                        modifier.open_circuit && continue
                        branch_admittance = _cascade_series_admittance(
                            modifier.resistance_ohm,
                            modifier.inductance_h,
                            modifier.capacitance_f,
                            angular_frequency,
                        )
                        _push_cascade_admittance!(
                            stamps,
                            current_station[phase],
                            section_from[phase],
                            current_station[phase],
                            section_from[phase],
                            branch_admittance,
                        )
                    end
                end
                auxiliary_nodes = Dict{Int,Int}()
                for branch in block.shunt_impedances
                    for terminal in (branch.from_terminal, branch.to_terminal)
                        terminal < 0 || continue
                        get!(auxiliary_nodes, -terminal) do
                            node_count += 1
                            node_count
                        end
                    end
                    from_node = _cascade_shunt_node(
                        branch.from_terminal,
                        section_from,
                        auxiliary_nodes,
                    )
                    to_node = _cascade_shunt_node(
                        branch.to_terminal,
                        section_from,
                        auxiliary_nodes,
                    )
                    branch_admittance = _cascade_series_admittance(
                        branch.resistance_ohm,
                        branch.inductance_h,
                        branch.capacitance_f,
                        angular_frequency,
                    )
                    _push_cascade_admittance!(
                        stamps,
                        from_node,
                        to_node,
                        from_node,
                        to_node,
                        branch_admittance,
                    )
                end
            end
            next_station = collect((node_count + 1):(node_count + phase_count))
            node_count += phase_count
            _push_cascade_section!(
                stamps,
                section_from,
                next_station,
                block,
                angular_frequency,
            )
            current_station = next_station
            first_section = false
        end
    end

    assembled = zeros(ComplexF64, node_count, node_count)
    for stamp in stamps
        _stamp_complex_admittance_entry!(assembled, stamp...)
    end
    terminal_indices = vcat(initial_station, current_station)
    terminal_set = Set(terminal_indices)
    internal_indices = Int[index for index in 1:node_count if !(index in terminal_set)]
    terminal_admittance = if isempty(internal_indices)
        copy(assembled[terminal_indices, terminal_indices])
    else
        @views assembled[terminal_indices, terminal_indices] -
            assembled[terminal_indices, internal_indices] *
            (assembled[internal_indices, internal_indices] \
             assembled[internal_indices, terminal_indices])
    end
    symmetry_error = maximum(
        abs,
        terminal_admittance - transpose(terminal_admittance);
        init = 0.0,
    )
    symmetry_error <= 1.0e-11 ||
        throw(ArgumentError("cascaded PI terminal admittance must be complex symmetric"))
    real_power_matrix = Hermitian(0.5 .* (terminal_admittance + adjoint(terminal_admittance)))
    minimum_real_power_eigenvalue = minimum(eigvals(real_power_matrix))
    minimum_real_power_eigenvalue >= -1.0e-11 ||
        throw(ArgumentError("cascaded PI terminal admittance must be passive"))

    first_section_owner = first(blocks).section
    return CascadedPhasePiEquivalent(
        name,
        phase_count,
        sum(block.multiplicity for block in blocks),
        first(blocks).section_scale,
        frequency,
        copy(first_section_owner.from_nodes),
        copy(first_section_owner.to_nodes),
        copy(first_section_owner.from_node_indices),
        copy(first_section_owner.to_node_indices),
        terminal_admittance,
        length(internal_indices),
        symmetry_error,
        minimum_real_power_eigenvalue,
        collect(blocks),
        sum(length(block.series_impedances) for block in blocks),
        sum(length(block.shunt_impedances) for block in blocks),
    )
end

function cascaded_phase_pi_equivalent(
    name::Symbol,
    section::CoupledLumpedPhasePiSection,
    section_count::Integer,
    section_scale::Real,
    frequency_hz::Real,
)
    block = cascaded_phase_pi_block(section, section_count, section_scale)
    return cascaded_phase_pi_equivalent(
        name,
        CascadedPhasePiBlock[block],
        frequency_hz,
    )
end

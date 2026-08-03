
using LinearAlgebra

using ..Branches
using ..Companion

import ..Branches: EMTElement,
                   stamp!,
                   stamp_conductance!,
                   stamp_history_current!,
                   trace_output_channel_count,
                   trace_output_channel_names!,
                   trace_output_is_public,
                   trace_output_values!,
                   update!

export PowerSemiconductorSwitch,
       DiodeValveSwitch,
       ThyristorValveSwitch,
       IGBTSwitch,
       MOSFETSwitch,
       PowerSemiconductorGateDriver,
       AntiparallelDiodeParameters,
       SeriesRCSnubber,
       PowerSemiconductorTerminalState,
       request_power_semiconductor_gate!,
       apply_power_semiconductor_gate_transition!,
       apply_power_semiconductor_forward_turn_on!,
       apply_power_semiconductor_forward_extinction!,
       apply_power_semiconductor_reverse_turn_on!,
       apply_power_semiconductor_reverse_extinction!,
       power_semiconductor_gate_transition_time,
       power_semiconductor_forward_turn_on_residual,
       power_semiconductor_forward_extinction_residual,
       power_semiconductor_reverse_turn_on_residual,
       power_semiconductor_reverse_extinction_residual,
       power_semiconductor_event_localization!,
       power_semiconductor_terminal_state,
       OVER16NonlinearInverseColumnState,
       SaturableInductorBranch,
       SaturatedTransformerCharacteristicPoint,
       SaturatedTransformerNonlinearArrays,
       SaturatedTransformerNonlinearIntake,
       SaturatedTransformerNonlinearSeed,
       SaturatedTransformerNonlinearSlopeBranch,
       SaturatedTransformerSegmentUpdate,
       SaturatedTransformerWindingSeed,
       diode_next_closed,
       diode_conductance,
       saturated_transformer_nonlinear_arrays,
       saturated_transformer_current_compensation_table,
       saturated_transformer_nonlinear_current_config,
       saturated_transformer_nonlinear_intake,
       saturated_transformer_nonlinear_slope_branches,
       set_saturated_transformer_nonlinear_slope!,
       saturated_transformer_sparse_admittance_update,
       saturated_transformer_segment_update,
       saturated_transformer_steady_state_branch_admittance,
       saturated_transformer_winding_branch_assembly,
       saturated_transformer_winding_branch_parameters,
       saturated_transformer_winding_current_config,
       saturated_transformer_winding_terminal_nodes,
       hysteretic_inductor_current_update,
       over16_nonlinear_source_column_assembly,
       over16_nonlinear_source_column_assembly!,
       over16_nonlinear_inverse_column_solution,
       over16_nonlinear_inverse_column_solution!,
       over16_simultaneous_zno_solution,
       over16_simultaneous_zno_solution!,
       over16_nonlinear_current_compensation_update,
       over16_nonlinear_current_compensation_update!,
       effective_inductance,
       saturable_inductor_companion

const OVER16_NONLINEAR_SOURCE_COLUMN_LABELS = (
    2300, 2320, 2321, 2306, 4377, 4319, 2315, 4323, 2322,
)
const OVER16_NONLINEAR_INVERSE_COLUMN_LABELS = (
    2410, 2413, 2420, 2422, 2423, 2450, 2500, 2510, 2513, 2520,
    2526, 2528, 2531, 2550, 2557, 2561, 2600, 2616,
)
const OVER16_NONLINEAR_CURRENT_COMPENSATION_LABELS = (
    1550, 1551, 3100, 3200, 3300, 3400, 3450, 3475, 3500,
    3534, 3800, 8259, 1555, 1559, 1570,
)
const OVER16_SATURATED_TRANSFORMER_SEGMENT_LABELS = (
    7642, 73910, 73908, 73912, 73913, 73915, 73925, 73926,
    73930, 73931, 73932, 73935, 73937, 3973, 3975, 73978,
)
const OVER16_SATURATED_TRANSFORMER_ADMITTANCE_LABELS = (
    3950, 4113, 4114, 3963, 3964, 3965, 3967, 3970, 3971, 3972, 3974,
)
const OVER16_SIMULTANEOUS_ZNO_LABELS = (
    3425, 3431, 3452, 3621, 3634, 3645, 3662, 3674, 3474,
    41, 55, 65, 66, 68, 3499, 3501, 80, 85, 92, 95, 100,
    3522, 139, 155, 165, 166, 180, 195, 196, 200,
)
const HYSTERETIC_INDUCTOR_REFERENCE_LABELS = (
    4317, 1100, 1116, 1118, 1110, 1119, 1120, 1123, 1125,
    1127, 1130, 1140, 1141, 1147, 1148, 1149, 1150, 1195,
    1212, 1312, 1315, 1319, 1322, 4372, 3950, 3975,
)

const SATURATED_TRANSFORMER_NONLINEAR_TYPE = -98
const HYSTERETIC_INDUCTOR_NONLINEAR_TYPE = -96
const EMPTY_SATURATED_TRANSFORMER_REFERENCE = Symbol("")

struct SaturatedTransformerCharacteristicPoint
    transformer_name::Symbol
    point_index::Int
    current::Float64
    flux::Float64
end

struct SaturatedTransformerWindingSeed
    transformer_name::Symbol
    winding_number::Int
    from_node::Symbol
    to_node::Symbol
    resistance::Union{Missing,Float64}
    inductance::Union{Missing,Float64}
    turns::Union{Missing,Float64}
    inherited_parameters::Bool
end

struct SaturatedTransformerNonlinearSeed
    transformer_name::Symbol
    reference_name::Symbol
    nonlinear_type::Int
    table_start_index::Int
    table_end_index::Int
    steady_state_current::Float64
    steady_state_flux::Float64
    magnetizing_resistance::Union{Missing,Float64}
    current_seed::Float64
    inherited_characteristic::Bool
end

struct SaturatedTransformerNonlinearIntake
    source::String
    seeds::Vector{SaturatedTransformerNonlinearSeed}
    characteristic_points::Vector{SaturatedTransformerCharacteristicPoint}
    windings::Vector{SaturatedTransformerWindingSeed}
    mutation_order::Tuple{Vararg{Symbol}}
    deferred_calls::Vector{Symbol}
    replacement_ready::Bool
end

struct SaturatedTransformerNonlinearArrays
    source::String
    transformer_names::Vector{Symbol}
    reference_names::Vector{Symbol}
    nonlinear_types::Vector{Int}
    table_start_indices::Vector{Int}
    table_end_indices::Vector{Int}
    steady_state_currents::Vector{Float64}
    steady_state_fluxes::Vector{Float64}
    magnetizing_resistances::Vector{Union{Missing,Float64}}
    current_seeds::Vector{Float64}
    inherited_characteristic_flags::Vector{Bool}
    characteristic_transformer_names::Vector{Symbol}
    characteristic_indices::Vector{Int}
    characteristic_currents::Vector{Float64}
    characteristic_fluxes::Vector{Float64}
    winding_transformer_names::Vector{Symbol}
    winding_numbers::Vector{Int}
    winding_from_nodes::Vector{Symbol}
    winding_to_nodes::Vector{Symbol}
    winding_resistances::Vector{Union{Missing,Float64}}
    winding_inductances::Vector{Union{Missing,Float64}}
    winding_turns::Vector{Union{Missing,Float64}}
    winding_inherited_flags::Vector{Bool}
end

struct SaturatedTransformerSegmentUpdate
    source::Symbol
    outcome::Symbol
    fortran_files::Tuple{Symbol}
    fortran_routines::Tuple{Symbol}
    fortran_labels::Tuple{Vararg{Int}}
    previous_current_segment::Float64
    current_segment::Float64
    previous_companion_current::Float64
    companion_current::Float64
    previous_stored_voltage::Float64
    stored_voltage::Float64
    branch_voltage::Float64
    voltage_increment::Float64
    effective_voltage::Float64
    companion_current_delta::Float64
    equivalent_current::Float64
    admittance_delta::Float64
    previous_table_index::Int
    current_table_index::Int
    segment_changed::Bool
    polarity_mismatch_warning::Bool
    polarity_reversed::Bool
    finitial_from_delta::Float64
    finitial_to_delta::Float64
    mutation_order::Tuple{Vararg{Symbol}}
    deferred_calls::Vector{Symbol}
    replacement_ready::Bool
end

mutable struct SaturatedTransformerNonlinearSlopeBranch <: EMTElement
    from_node::Int
    to_node::Int
    conductance::Float64
end

function set_saturated_transformer_nonlinear_slope!(
    branch::SaturatedTransformerNonlinearSlopeBranch,
    conductance::Real,
)
    g = Float64(conductance)
    isfinite(g) ||
        throw(ArgumentError("saturated transformer nonlinear slope conductance must be finite"))
    changed = branch.conductance != g
    branch.conductance = g
    return changed
end

function SaturatedTransformerNonlinearSlopeBranch(
    from_node::Integer,
    to_node::Integer,
    conductance::Real,
)
    source_node = Int(from_node)
    sink_node = Int(to_node)
    source_node > 0 ||
        throw(ArgumentError("saturated transformer nonlinear slope from_node must be positive"))
    sink_node >= 0 ||
        throw(ArgumentError("saturated transformer nonlinear slope to_node must be nonnegative"))
    g = Float64(conductance)
    isfinite(g) || throw(ArgumentError("saturated transformer nonlinear slope conductance must be finite"))
    return SaturatedTransformerNonlinearSlopeBranch(source_node, sink_node, g)
end

mutable struct OVER16NonlinearInverseColumnState
    ntot::Int
    ncomp::Int
    znonl::Vector{Float64}
    anonl::Vector{Float64}
    voltbc::Vector{Float64}
    vzero::Vector{Float64}
    vnonl::Vector{Float64}
    ilast::Vector{Int}
    curr::Vector{Float64}
    cursub::Vector{Float64}
    cchar::Vector{Float64}
    vchar::Vector{Float64}
    gslope::Vector{Float64}
    source_column_update_count::Int
    update_count::Int
    current_update_count::Int
end

function OVER16NonlinearInverseColumnState(
    ntot::Int,
    ncomp::Int;
    nonlinear_count::Int=0,
)
    _check_over16_nonlinear_dimensions(ntot, ncomp)
    nonlinear_count >= 0 || throw(ArgumentError("nonlinear_count must be nonnegative"))
    return OVER16NonlinearInverseColumnState(
        ntot,
        ncomp,
        zeros(Float64, ntot * ncomp),
        zeros(Float64, nonlinear_count),
        zeros(Float64, ncomp),
        zeros(Float64, nonlinear_count),
        Float64[],
        zeros(Int, nonlinear_count),
        zeros(Float64, nonlinear_count),
        zeros(Float64, nonlinear_count),
        Float64[],
        Float64[],
        Float64[],
        0,
        0,
        0,
    )
end

function OVER16NonlinearInverseColumnState(
    znonl::AbstractVector{<:Real},
    ntot::Int,
    ncomp::Int;
    anonl::AbstractVector{<:Real}=Float64[],
    voltbc::AbstractVector{<:Real}=Float64[],
    vzero::AbstractVector{<:Real}=Float64[],
    vnonl::AbstractVector{<:Real}=Float64[],
    ilast::AbstractVector{Int}=Int[],
    curr::AbstractVector{<:Real}=Float64[],
    cursub::AbstractVector{<:Real}=Float64[],
    cchar::AbstractVector{<:Real}=Float64[],
    vchar::AbstractVector{<:Real}=Float64[],
    gslope::AbstractVector{<:Real}=Float64[],
    source_column_update_count::Int=0,
    update_count::Int=0,
    current_update_count::Int=0,
)
    _check_over16_nonlinear_dimensions(ntot, ncomp)
    _check_over16_znonl_layout(znonl, ntot, ncomp)
    source_column_update_count >= 0 ||
        throw(ArgumentError("source_column_update_count must be nonnegative"))
    update_count >= 0 || throw(ArgumentError("update_count must be nonnegative"))
    current_update_count >= 0 ||
        throw(ArgumentError("current_update_count must be nonnegative"))
    nonlinear_values = Float64.(anonl)
    scratch = isempty(voltbc) ? zeros(Float64, ncomp) : Float64.(voltbc)
    vzero_values = isempty(vzero) ? zeros(Float64, length(nonlinear_values)) : Float64.(vzero)
    vnonl_values = Float64.(vnonl)
    ilast_values = isempty(ilast) ? zeros(Int, length(nonlinear_values)) : collect(ilast)
    current_values = isempty(curr) ? zeros(Float64, length(nonlinear_values)) : Float64.(curr)
    cursub_values = isempty(cursub) ? zeros(Float64, length(nonlinear_values)) : Float64.(cursub)
    cchar_values = Float64.(cchar)
    vchar_values = Float64.(vchar)
    gslope_values = Float64.(gslope)
    length(scratch) == ncomp || throw(ArgumentError("voltbc length must match ncomp"))
    length(vzero_values) == length(nonlinear_values) ||
        throw(ArgumentError("vzero length must match anonl length"))
    isempty(vnonl_values) || length(vnonl_values) == length(nonlinear_values) ||
        throw(ArgumentError("vnonl length must match anonl length"))
    length(ilast_values) == length(nonlinear_values) ||
        throw(ArgumentError("ilast length must match anonl length"))
    length(current_values) == length(nonlinear_values) ||
        throw(ArgumentError("curr length must match anonl length"))
    _check_over16_finite_vector("anonl", nonlinear_values)
    _check_over16_finite_vector("voltbc", scratch)
    _check_over16_vzero_vector(vzero_values)
    _check_over16_vnonl_vector(vnonl_values)
    _check_over16_finite_vector("curr", current_values)
    _check_over16_finite_vector("cursub", cursub_values)
    _check_over16_finite_vector("cchar", cchar_values)
    _check_over16_finite_vector("vchar", vchar_values)
    _check_over16_finite_vector("gslope", gslope_values)
    return OVER16NonlinearInverseColumnState(
        ntot,
        ncomp,
        Float64.(znonl),
        nonlinear_values,
        scratch,
        vzero_values,
        vnonl_values,
        ilast_values,
        current_values,
        cursub_values,
        cchar_values,
        vchar_values,
        gslope_values,
        source_column_update_count,
        update_count,
        current_update_count,
    )
end

function stamp!(
    y::AbstractMatrix{Float64},
    rhs::AbstractVector{Float64},
    branch::SaturatedTransformerNonlinearSlopeBranch,
    t::Float64,
    dt::Float64,
)
    stamp_conductance!(y, branch.from_node, branch.to_node, branch.conductance)
    return nothing
end

update!(
    branch::SaturatedTransformerNonlinearSlopeBranch,
    voltage,
    dt::Float64,
) = nothing

function _required_saturated_transformer_float(value, label::AbstractString)
    ismissing(value) && throw(ArgumentError("$label is required"))
    return Float64(value)
end

function _saturated_transformer_rows_by_name(rows)
    by_name = Dict{Symbol,Any}()
    for row in rows
        name = getproperty(row, :name)
        haskey(by_name, name) && throw(ArgumentError("duplicate saturated transformer name $name"))
        by_name[name] = row
    end
    return by_name
end

function _saturated_transformer_reference_seed(
    seeds::Dict{Symbol,SaturatedTransformerNonlinearSeed},
    reference_name::Symbol,
)
    haskey(seeds, reference_name) ||
        throw(ArgumentError("missing saturated transformer copy reference $reference_name"))
    return seeds[reference_name]
end

function saturated_transformer_nonlinear_intake(parsed_intake)
    source = String(getproperty(parsed_intake, :source))
    transformer_rows = collect(getproperty(parsed_intake, :transformers))
    breakpoint_rows = collect(getproperty(parsed_intake, :breakpoints))
    winding_rows = collect(getproperty(parsed_intake, :windings))

    by_name = _saturated_transformer_rows_by_name(transformer_rows)
    breakpoint_groups = Dict{Symbol,Vector{Any}}()
    for row in breakpoint_rows
        name = getproperty(row, :transformer_name)
        push!(get!(breakpoint_groups, name, Any[]), row)
    end

    characteristic_points = SaturatedTransformerCharacteristicPoint[]
    windings = SaturatedTransformerWindingSeed[]
    seeds = SaturatedTransformerNonlinearSeed[]
    seed_by_name = Dict{Symbol,SaturatedTransformerNonlinearSeed}()
    next_table_index = 1

    for row in transformer_rows
        name = getproperty(row, :name)
        reference_name = getproperty(row, :reference_name)
        inherited = reference_name != EMPTY_SATURATED_TRANSFORMER_REFERENCE
        local table_start::Int
        local table_end::Int
        local steady_state_current::Float64
        local steady_state_flux::Float64
        local magnetizing_resistance::Union{Missing,Float64}
        if inherited
            haskey(by_name, reference_name) ||
                throw(ArgumentError("missing saturated transformer reference $reference_name"))
            reference_seed = _saturated_transformer_reference_seed(seed_by_name, reference_name)
            table_start = reference_seed.table_start_index
            table_end = reference_seed.table_end_index
            steady_state_current = reference_seed.steady_state_current
            steady_state_flux = reference_seed.steady_state_flux
            magnetizing_resistance = reference_seed.magnetizing_resistance
        else
            points = get(breakpoint_groups, name, Any[])
            isempty(points) &&
                throw(ArgumentError("saturated transformer $name has no characteristic points"))
            table_start = next_table_index
            for (point_offset, point) in enumerate(points)
                push!(
                    characteristic_points,
                    SaturatedTransformerCharacteristicPoint(
                        name,
                        table_start + point_offset - 1,
                        Float64(getproperty(point, :current)),
                        Float64(getproperty(point, :flux)),
                    ),
                )
            end
            table_end = table_start + length(points) - 1
            next_table_index = table_end + 1
            steady_state_current = _required_saturated_transformer_float(
                getproperty(row, :initial_current),
                "saturated transformer $name steady-state current",
            )
            steady_state_flux = _required_saturated_transformer_float(
                getproperty(row, :initial_flux),
                "saturated transformer $name steady-state flux",
            )
            magnetizing_resistance = getproperty(row, :magnetizing_resistance)
        end

        seed = SaturatedTransformerNonlinearSeed(
            name,
            reference_name,
            SATURATED_TRANSFORMER_NONLINEAR_TYPE,
            table_start,
            table_end,
            steady_state_current,
            steady_state_flux,
            magnetizing_resistance,
            1.0,
            inherited,
        )
        push!(seeds, seed)
        seed_by_name[name] = seed
    end

    for row in winding_rows
        push!(
            windings,
            SaturatedTransformerWindingSeed(
                getproperty(row, :transformer_name),
                Int(getproperty(row, :winding_number)),
                getproperty(row, :from_node),
                getproperty(row, :to_node),
                getproperty(row, :resistance),
                getproperty(row, :inductance),
                getproperty(row, :turns),
                Bool(getproperty(row, :inherited_parameters)),
            ),
        )
    end

    return SaturatedTransformerNonlinearIntake(
        source,
        seeds,
        characteristic_points,
        windings,
        (
            :saturated_transformer_header,
            :characteristic_points,
            :winding_rows,
            :copy_rows,
        ),
        [
            :over2_branch_table_mutation,
            :last14_initial_flux_solution,
            :over16_current_compensation,
        ],
        false,
    )
end

function saturated_transformer_nonlinear_arrays(intake::SaturatedTransformerNonlinearIntake)
    return SaturatedTransformerNonlinearArrays(
        intake.source,
        [seed.transformer_name for seed in intake.seeds],
        [seed.reference_name for seed in intake.seeds],
        [seed.nonlinear_type for seed in intake.seeds],
        [seed.table_start_index for seed in intake.seeds],
        [seed.table_end_index for seed in intake.seeds],
        [seed.steady_state_current for seed in intake.seeds],
        [seed.steady_state_flux for seed in intake.seeds],
        [seed.magnetizing_resistance for seed in intake.seeds],
        [seed.current_seed for seed in intake.seeds],
        [seed.inherited_characteristic for seed in intake.seeds],
        [point.transformer_name for point in intake.characteristic_points],
        [point.point_index for point in intake.characteristic_points],
        [point.current for point in intake.characteristic_points],
        [point.flux for point in intake.characteristic_points],
        [winding.transformer_name for winding in intake.windings],
        [winding.winding_number for winding in intake.windings],
        [winding.from_node for winding in intake.windings],
        [winding.to_node for winding in intake.windings],
        [winding.resistance for winding in intake.windings],
        [winding.inductance for winding in intake.windings],
        [winding.turns for winding in intake.windings],
        [winding.inherited_parameters for winding in intake.windings],
    )
end

function saturated_transformer_nonlinear_arrays(parsed_intake)
    return saturated_transformer_nonlinear_arrays(
        saturated_transformer_nonlinear_intake(parsed_intake),
    )
end

function _saturated_transformer_positive_delta(delta2::Real)
    delta = Float64(delta2)
    isfinite(delta) && delta > 0.0 ||
        throw(ArgumentError("delta2 must be finite and positive"))
    return delta
end

function _saturated_transformer_indexed_float_values(
    indices::AbstractVector{Int},
    values::AbstractVector{<:Real},
    label::AbstractString,
)
    length(indices) == length(values) ||
        throw(ArgumentError("$label index and value lengths must match"))
    isempty(indices) && return Float64[]
    point_count = maximum(indices)
    point_count > 0 || throw(ArgumentError("$label indices must be positive"))
    result = fill(NaN, point_count)
    seen = falses(point_count)
    for (position, index) in enumerate(indices)
        1 <= index <= point_count ||
            throw(ArgumentError("$label indices must be positive"))
        value = Float64(values[position])
        isfinite(value) || throw(ArgumentError("$label values must be finite"))
        if seen[index]
            result[index] == value ||
                throw(ArgumentError("$label index $index has inconsistent duplicate values"))
        else
            result[index] = value
            seen[index] = true
        end
    end
    all(seen) || throw(ArgumentError("$label indices must be contiguous"))
    return result
end

function saturated_transformer_current_compensation_table(
    arrays::SaturatedTransformerNonlinearArrays;
    delta2::Real=1.0,
)
    delta = _saturated_transformer_positive_delta(delta2)
    characteristic_indices = Int.(arrays.characteristic_indices)
    currents = _saturated_transformer_indexed_float_values(
        characteristic_indices,
        arrays.characteristic_currents,
        "saturated transformer current characteristic",
    )
    fluxes = _saturated_transformer_indexed_float_values(
        characteristic_indices,
        arrays.characteristic_fluxes,
        "saturated transformer flux characteristic",
    )
    length(currents) == length(fluxes) ||
        throw(ArgumentError("saturated transformer current and flux tables must match"))

    cchar = zeros(Float64, length(currents))
    gslope = zeros(Float64, length(currents))
    processed_ranges = Set{Tuple{Int,Int}}()
    for (table_start, table_end) in zip(
        arrays.table_start_indices,
        arrays.table_end_indices,
    )
        start_index = Int(table_start)
        end_index = Int(table_end)
        1 <= start_index <= end_index <= length(currents) ||
            throw(ArgumentError("saturated transformer table range is outside characteristics"))
        range_key = (start_index, end_index)
        range_key in processed_ranges && continue
        push!(processed_ranges, range_key)
        previous_current = 0.0
        previous_flux = 0.0
        for table_index in start_index:end_index
            current = currents[table_index]
            flux = fluxes[table_index]
            current > previous_current ||
                throw(ArgumentError("saturated transformer currents must strictly increase"))
            flux > previous_flux ||
                throw(ArgumentError("saturated transformer fluxes must strictly increase"))
            current_delta = current - previous_current
            flux_delta = flux - previous_flux
            flux_per_current = flux_delta / current_delta
            gslope[table_index] = delta * current_delta / flux_delta
            cchar[table_index] = previous_flux - flux_per_current * previous_current
            previous_current = current
            previous_flux = flux
        end
    end
    return (
        source = :saturated_transformer_current_compensation_table,
        outcome = :nonlinear_current_characteristic_arrays,
        cchar = cchar,
        gslope = gslope,
        vchar = fluxes,
        delta2 = delta,
        fortran_files = (:OVER2_FOR,),
        fortran_labels = (73420, 73424, 6472, 73430, 73436),
        mutation_order = (:breakpoint_read, :segment_slope, :segment_intercept),
        replacement_ready = false,
    )
end

function saturated_transformer_current_compensation_table(parsed_intake; kwargs...)
    return saturated_transformer_current_compensation_table(
        saturated_transformer_nonlinear_arrays(parsed_intake);
        kwargs...,
    )
end

function _saturated_transformer_int_vector(
    values::AbstractVector{<:Integer},
    label::AbstractString,
    count::Int,
)
    length(values) == count ||
        throw(ArgumentError("$label length must match saturated transformer count"))
    return Int.(values)
end

function _saturated_transformer_winding_by_key(arrays::SaturatedTransformerNonlinearArrays)
    by_key = Dict{Tuple{Symbol,Int},Int}()
    for (position, (name, winding_number)) in enumerate(zip(
        arrays.winding_transformer_names,
        arrays.winding_numbers,
    ))
        key = (name, Int(winding_number))
        haskey(by_key, key) &&
            throw(ArgumentError("duplicate saturated transformer winding $(key[1])/$(key[2])"))
        by_key[key] = position
    end
    return by_key
end

function _saturated_transformer_resolve_node(
    node::Symbol,
    node_map::AbstractDict{Symbol,<:Integer},
    label::AbstractString;
    reference_node_index::Integer=0,
)
    if node == EMPTY_SATURATED_TRANSFORMER_REFERENCE
        return Int(reference_node_index)
    end
    haskey(node_map, node) ||
        throw(ArgumentError("missing saturated transformer $label node $(String(node))"))
    node_index = Int(node_map[node])
    node_index >= 0 ||
        throw(ArgumentError("saturated transformer $label node index must be nonnegative"))
    return node_index
end

function saturated_transformer_winding_terminal_nodes(
    arrays::SaturatedTransformerNonlinearArrays,
    node_map::AbstractDict{Symbol,<:Integer};
    winding_number::Integer=1,
    reference_node_index::Integer=0,
)
    requested_winding = Int(winding_number)
    requested_winding > 0 ||
        throw(ArgumentError("saturated transformer winding_number must be positive"))
    winding_by_key = _saturated_transformer_winding_by_key(arrays)
    terminal_names = Symbol[]
    terminal_indices = Int[]
    matched_windings = Int[]
    reference_terminal_count = 0
    for transformer_name in arrays.transformer_names
        key = (transformer_name, requested_winding)
        haskey(winding_by_key, key) ||
            throw(ArgumentError("missing saturated transformer winding $(transformer_name)/$(requested_winding)"))
        row_index = winding_by_key[key]
        terminal_name = arrays.winding_to_nodes[row_index]
        terminal_index = _saturated_transformer_resolve_node(
            terminal_name,
            node_map,
            "terminal";
            reference_node_index = reference_node_index,
        )
        terminal_index == Int(reference_node_index) && (reference_terminal_count += 1)
        push!(terminal_names, terminal_name)
        push!(terminal_indices, terminal_index)
        push!(matched_windings, requested_winding)
    end
    return (
        source = :saturated_transformer_winding_terminal_nodes,
        outcome = :terminal_binding,
        transformer_names = copy(arrays.transformer_names),
        winding_numbers = matched_windings,
        terminal_node_names = terminal_names,
        terminal_node_indices = terminal_indices,
        reference_terminal_count = reference_terminal_count,
        fortran_files = (:OVER2_FOR,),
        fortran_labels = (4056,),
        mutation_order = (:winding_row_scan, :terminal_node_lookup),
        replacement_ready = false,
    )
end

function saturated_transformer_winding_branch_parameters(
    arrays::SaturatedTransformerNonlinearArrays,
    node_map::AbstractDict{Symbol,<:Integer};
    reference_node_index::Integer=0,
)
    winding_count = length(arrays.winding_transformer_names)
    length(arrays.winding_numbers) == winding_count ||
        throw(ArgumentError("saturated transformer winding numbers length mismatch"))
    length(arrays.winding_from_nodes) == winding_count ||
        throw(ArgumentError("saturated transformer winding from-node length mismatch"))
    length(arrays.winding_to_nodes) == winding_count ||
        throw(ArgumentError("saturated transformer winding to-node length mismatch"))
    from_indices = Int[]
    to_indices = Int[]
    reference_count = 0
    for (from_name, to_name) in zip(arrays.winding_from_nodes, arrays.winding_to_nodes)
        from_index = _saturated_transformer_resolve_node(
            from_name,
            node_map,
            "winding from";
            reference_node_index = reference_node_index,
        )
        to_index = _saturated_transformer_resolve_node(
            to_name,
            node_map,
            "winding to";
            reference_node_index = reference_node_index,
        )
        from_index == Int(reference_node_index) && (reference_count += 1)
        to_index == Int(reference_node_index) && (reference_count += 1)
        push!(from_indices, from_index)
        push!(to_indices, to_index)
    end
    return (
        source = :saturated_transformer_winding_branch_parameters,
        outcome = :winding_branch_parameters,
        transformer_names = copy(arrays.winding_transformer_names),
        winding_numbers = Int.(arrays.winding_numbers),
        from_node_names = copy(arrays.winding_from_nodes),
        to_node_names = copy(arrays.winding_to_nodes),
        from_node_indices = from_indices,
        to_node_indices = to_indices,
        resistances = copy(arrays.winding_resistances),
        inductances = copy(arrays.winding_inductances),
        turns = copy(arrays.winding_turns),
        inherited_parameter_flags = Bool.(arrays.winding_inherited_flags),
        winding_count = winding_count,
        parameterized_winding_count = Base.count(!, arrays.winding_inherited_flags),
        reference_terminal_count = reference_count,
        fortran_files = (:OVER2_FOR,),
        fortran_labels = (4010, 74003, 74046, 54047, 181, 54118, 4056),
        mutation_order = (
            :winding_row_scan,
            :node_lookup,
            :parameter_binding,
        ),
        replacement_ready = false,
    )
end

function saturated_transformer_winding_branch_parameters(parsed_intake, node_map; kwargs...)
    return saturated_transformer_winding_branch_parameters(
        saturated_transformer_nonlinear_arrays(parsed_intake),
        node_map;
        kwargs...,
    )
end

function _saturated_transformer_reference_by_name(arrays::SaturatedTransformerNonlinearArrays)
    reference_by_name = Dict{Symbol,Symbol}()
    for (name, reference_name) in zip(arrays.transformer_names, arrays.reference_names)
        haskey(reference_by_name, name) &&
            throw(ArgumentError("duplicate saturated transformer seed $name"))
        reference_by_name[name] = reference_name
    end
    return reference_by_name
end

function _saturated_transformer_required_winding_value(
    arrays::SaturatedTransformerNonlinearArrays,
    winding_by_key::AbstractDict{Tuple{Symbol,Int},Int},
    reference_by_name::AbstractDict{Symbol,Symbol},
    row_index::Int,
    values::AbstractVector,
    field::AbstractString,
)
    winding_name = arrays.winding_transformer_names[row_index]
    winding_number = Int(arrays.winding_numbers[row_index])
    value = values[row_index]
    if !arrays.winding_inherited_flags[row_index]
        ismissing(value) && throw(ArgumentError(
            "saturated transformer $(winding_name)/$(winding_number) missing $field",
        ))
        return Float64(value)
    end

    haskey(reference_by_name, winding_name) ||
        throw(ArgumentError("missing saturated transformer reference for $winding_name"))
    reference_name = reference_by_name[winding_name]
    reference_name != EMPTY_SATURATED_TRANSFORMER_REFERENCE ||
        throw(ArgumentError("saturated transformer $winding_name marks inherited winding without a reference"))
    seen = Set{Tuple{Symbol,Int}}()
    while true
        key = (reference_name, winding_number)
        key in seen && throw(ArgumentError(
            "cyclic saturated transformer winding reference $(winding_name)/$(winding_number)",
        ))
        push!(seen, key)
        haskey(winding_by_key, key) ||
            throw(ArgumentError("missing saturated transformer reference winding $(key[1])/$(key[2])"))
        reference_index = winding_by_key[key]
        reference_value = values[reference_index]
        if !arrays.winding_inherited_flags[reference_index]
            ismissing(reference_value) && throw(ArgumentError(
                "saturated transformer reference winding $(key[1])/$(key[2]) missing $field",
            ))
            return Float64(reference_value)
        end
        haskey(reference_by_name, key[1]) ||
            throw(ArgumentError("missing saturated transformer nested reference for $(key[1])"))
        reference_name = reference_by_name[key[1]]
        reference_name != EMPTY_SATURATED_TRANSFORMER_REFERENCE ||
            throw(ArgumentError("saturated transformer $(key[1]) inherited winding has no reference"))
    end
end

function _saturated_transformer_internal_node_indices(
    names::AbstractVector{Symbol},
    node_map::AbstractDict{Symbol,<:Integer},
    internal_node_map::AbstractDict{Symbol,<:Integer},
    reference_node_index::Integer,
)
    assigned_node_indices = Int[Int(value) for value in values(node_map)]
    append!(assigned_node_indices, Int(value) for value in values(internal_node_map))
    next_index = isempty(assigned_node_indices) ?
        Int(reference_node_index) + 1 :
        maximum(assigned_node_indices) + 1
    indices = Int[]
    for name in names
        if haskey(internal_node_map, name)
            node_index = Int(internal_node_map[name])
        elseif haskey(node_map, name)
            node_index = Int(node_map[name])
        else
            while next_index in assigned_node_indices
                next_index += 1
            end
            node_index = next_index
            push!(assigned_node_indices, node_index)
            next_index += 1
        end
        node_index > Int(reference_node_index) ||
            throw(ArgumentError("saturated transformer internal top node must be above the reference node"))
        push!(indices, node_index)
    end
    return indices
end

function saturated_transformer_winding_branch_assembly(
    arrays::SaturatedTransformerNonlinearArrays,
    node_map::AbstractDict{Symbol,<:Integer};
    internal_node_map::AbstractDict{Symbol,<:Integer}=Dict{Symbol,Int}(),
    nonlinear_winding_number::Integer=1,
    reference_node_index::Integer=0,
)
    winding_count = length(arrays.winding_transformer_names)
    length(arrays.winding_numbers) == winding_count ||
        throw(ArgumentError("saturated transformer winding numbers length mismatch"))
    length(arrays.winding_from_nodes) == winding_count ||
        throw(ArgumentError("saturated transformer winding from-node length mismatch"))
    length(arrays.winding_to_nodes) == winding_count ||
        throw(ArgumentError("saturated transformer winding to-node length mismatch"))
    requested_winding = Int(nonlinear_winding_number)
    requested_winding > 0 ||
        throw(ArgumentError("saturated transformer nonlinear_winding_number must be positive"))

    winding_by_key = _saturated_transformer_winding_by_key(arrays)
    reference_by_name = _saturated_transformer_reference_by_name(arrays)
    internal_indices = _saturated_transformer_internal_node_indices(
        arrays.transformer_names,
        node_map,
        internal_node_map,
        reference_node_index,
    )
    internal_by_name = Dict(
        name => internal_indices[index]
        for (index, name) in enumerate(arrays.transformer_names)
    )

    winding_from_indices = Int[]
    winding_terminal_indices = Int[]
    winding_internal_indices = Int[]
    resolved_resistances = Float64[]
    resolved_inductances = Float64[]
    resolved_turns = Float64[]
    for row_index in eachindex(arrays.winding_transformer_names)
        transformer_name = arrays.winding_transformer_names[row_index]
        haskey(internal_by_name, transformer_name) ||
            throw(ArgumentError("missing saturated transformer internal node for $transformer_name"))
        push!(winding_internal_indices, internal_by_name[transformer_name])
        push!(
            winding_from_indices,
            _saturated_transformer_resolve_node(
                arrays.winding_from_nodes[row_index],
                node_map,
                "winding from";
                reference_node_index = reference_node_index,
            ),
        )
        push!(
            winding_terminal_indices,
            _saturated_transformer_resolve_node(
                arrays.winding_to_nodes[row_index],
                node_map,
                "winding to";
                reference_node_index = reference_node_index,
            ),
        )
        push!(
            resolved_resistances,
            _saturated_transformer_required_winding_value(
                arrays,
                winding_by_key,
                reference_by_name,
                row_index,
                arrays.winding_resistances,
                "resistance",
            ),
        )
        push!(
            resolved_inductances,
            _saturated_transformer_required_winding_value(
                arrays,
                winding_by_key,
                reference_by_name,
                row_index,
                arrays.winding_inductances,
                "inductance",
            ),
        )
        push!(
            resolved_turns,
            _saturated_transformer_required_winding_value(
                arrays,
                winding_by_key,
                reference_by_name,
                row_index,
                arrays.winding_turns,
                "turns",
            ),
        )
    end

    nonlinear_from_indices = Int[]
    nonlinear_to_indices = Int[]
    magnetizing_from_indices = Int[]
    magnetizing_to_indices = Int[]
    magnetizing_resistances = Float64[]
    for (transformer_index, transformer_name) in enumerate(arrays.transformer_names)
        key = (transformer_name, requested_winding)
        haskey(winding_by_key, key) ||
            throw(ArgumentError("missing saturated transformer nonlinear winding $(key[1])/$(key[2])"))
        winding_index = winding_by_key[key]
        top_index = internal_indices[transformer_index]
        terminal_index = winding_terminal_indices[winding_index]
        push!(nonlinear_from_indices, top_index)
        push!(nonlinear_to_indices, terminal_index)
        magnetizing = arrays.magnetizing_resistances[transformer_index]
        if !ismissing(magnetizing) && Float64(magnetizing) != 0.0
            push!(magnetizing_from_indices, top_index)
            push!(magnetizing_to_indices, terminal_index)
            push!(magnetizing_resistances, Float64(magnetizing))
        end
    end

    return (
        source = :saturated_transformer_winding_branch_assembly,
        outcome = :winding_branch_assembly,
        transformer_names = copy(arrays.transformer_names),
        reference_names = copy(arrays.reference_names),
        internal_top_node_names = copy(arrays.transformer_names),
        internal_top_node_indices = internal_indices,
        winding_transformer_names = copy(arrays.winding_transformer_names),
        winding_numbers = Int.(arrays.winding_numbers),
        winding_from_node_names = copy(arrays.winding_from_nodes),
        winding_to_node_names = copy(arrays.winding_to_nodes),
        winding_from_node_indices = winding_from_indices,
        winding_terminal_node_indices = winding_terminal_indices,
        winding_internal_top_node_indices = winding_internal_indices,
        series_branch_from_node_indices = winding_from_indices,
        series_branch_to_node_indices = winding_internal_indices,
        series_branch_resistances = resolved_resistances,
        series_branch_inductances = resolved_inductances,
        series_branch_turns = resolved_turns,
        series_branch_storage_capacitances = zeros(Float64, winding_count),
        inherited_parameter_flags = Bool.(arrays.winding_inherited_flags),
        magnetizing_branch_from_node_indices = magnetizing_from_indices,
        magnetizing_branch_to_node_indices = magnetizing_to_indices,
        magnetizing_branch_resistances = magnetizing_resistances,
        magnetizing_branch_count = length(magnetizing_resistances),
        nonlinear_winding_number = requested_winding,
        nonlinear_from_node_indices = nonlinear_from_indices,
        nonlinear_to_node_indices = nonlinear_to_indices,
        nonlinear_binding_count = length(nonlinear_from_indices),
        fortran_files = (:OVER2_FOR,),
        fortran_labels = (4044, 4047, 4049, 4048, 4052, 4054, 4056),
        mutation_order = (
            :resolve_inherited_winding_parameters,
            :allocate_internal_top_nodes,
            :series_branch_storage,
            :magnetizing_branch_storage,
            :nonlinear_terminal_binding,
        ),
        replacement_ready = false,
    )
end

function saturated_transformer_winding_branch_assembly(parsed_intake, node_map; kwargs...)
    return saturated_transformer_winding_branch_assembly(
        saturated_transformer_nonlinear_arrays(parsed_intake),
        node_map;
        kwargs...,
    )
end

function _saturated_transformer_positive_reactance(value::Real, label::AbstractString)
    reactance = Float64(value)
    isfinite(reactance) && reactance > 0.0 ||
        throw(ArgumentError("$label must be finite and positive"))
    return reactance
end

function _saturated_transformer_branch_admittance(
    resistance::Real,
    reactance::Real,
    label::AbstractString,
)
    r = Float64(resistance)
    x = _saturated_transformer_positive_reactance(reactance, label)
    return inv(complex(r, x))
end

function saturated_transformer_steady_state_branch_admittance(
    arrays::SaturatedTransformerNonlinearArrays,
    node_map::AbstractDict{Symbol,<:Integer};
    internal_node_map::AbstractDict{Symbol,<:Integer}=Dict{Symbol,Int}(),
    nonlinear_winding_number::Integer=1,
    reference_node_index::Integer=0,
    primary_winding_number::Integer=1,
    ideal_winding_number::Integer=2,
    reactance_units::Real=2.0 * pi * 60.0,
)
    xunits = Float64(reactance_units)
    isfinite(xunits) && xunits > 0.0 ||
        throw(ArgumentError("reactance_units must be finite and positive"))
    primary_winding = Int(primary_winding_number)
    ideal_winding = Int(ideal_winding_number)
    primary_winding > 0 || throw(ArgumentError("primary_winding_number must be positive"))
    ideal_winding > 0 || throw(ArgumentError("ideal_winding_number must be positive"))
    primary_winding != ideal_winding ||
        throw(ArgumentError("primary and ideal winding numbers must differ"))

    branch_assembly = saturated_transformer_winding_branch_assembly(
        arrays,
        node_map;
        internal_node_map = internal_node_map,
        nonlinear_winding_number = nonlinear_winding_number,
        reference_node_index = reference_node_index,
    )
    winding_by_key = _saturated_transformer_winding_by_key(arrays)

    scalar_from_nodes = Int[]
    scalar_to_nodes = Int[]
    scalar_admittances = ComplexF64[]
    scalar_winding_numbers = Int[]
    scalar_reactances = Float64[]
    ideal_low_from_nodes = Int[]
    ideal_low_to_nodes = Int[]
    ideal_internal_from_nodes = Int[]
    ideal_internal_to_nodes = Int[]
    ideal_low_admittances = ComplexF64[]
    ideal_mutual_admittances = ComplexF64[]
    ideal_internal_admittances = ComplexF64[]
    ideal_primary_storage_coefficients = Float64[]
    ideal_mutual_storage_coefficients = Float64[]
    ideal_secondary_storage_coefficients = Float64[]
    ideal_resistive_storage_coefficients = Float64[]
    ideal_low_resistances = Float64[]
    ideal_low_reactances = Float64[]
    ideal_internal_reactances = Float64[]
    ideal_primary_turns = Float64[]
    ideal_turns_ratios = Float64[]
    ideal_transformer_names = Symbol[]
    linearized_nonlinear_from_nodes = Int[]
    linearized_nonlinear_to_nodes = Int[]
    linearized_nonlinear_admittances = ComplexF64[]
    linearized_nonlinear_reactances = Float64[]

    handled_windings = falses(length(arrays.winding_transformer_names))
    for (transformer_index, transformer_name) in enumerate(arrays.transformer_names)
        primary_key = (transformer_name, primary_winding)
        haskey(winding_by_key, primary_key) || continue
        primary_index = winding_by_key[primary_key]
        push!(scalar_from_nodes, branch_assembly.series_branch_from_node_indices[primary_index])
        push!(scalar_to_nodes, branch_assembly.series_branch_to_node_indices[primary_index])
        primary_reactance = _saturated_transformer_positive_reactance(
            branch_assembly.series_branch_inductances[primary_index],
            "saturated transformer primary leakage reactance",
        )
        push!(
            scalar_admittances,
            _saturated_transformer_branch_admittance(
                branch_assembly.series_branch_resistances[primary_index],
                primary_reactance,
                "saturated transformer primary leakage reactance",
            ),
        )
        push!(scalar_winding_numbers, primary_winding)
        push!(scalar_reactances, primary_reactance)
        handled_windings[primary_index] = true

        steady_current = _saturated_transformer_positive_reactance(
            arrays.steady_state_currents[transformer_index],
            "saturated transformer steady-state current",
        )
        steady_flux = _saturated_transformer_positive_reactance(
            arrays.steady_state_fluxes[transformer_index],
            "saturated transformer steady-state flux",
        )
        nonlinear_reactance = xunits * steady_flux / steady_current
        push!(linearized_nonlinear_from_nodes, branch_assembly.nonlinear_from_node_indices[transformer_index])
        push!(linearized_nonlinear_to_nodes, branch_assembly.nonlinear_to_node_indices[transformer_index])
        push!(linearized_nonlinear_admittances, inv(complex(0.0, nonlinear_reactance)))
        push!(linearized_nonlinear_reactances, nonlinear_reactance)

        ideal_key = (transformer_name, ideal_winding)
        haskey(winding_by_key, ideal_key) || continue
        ideal_index = winding_by_key[ideal_key]
        low_resistance = Float64(branch_assembly.series_branch_resistances[ideal_index])
        isfinite(low_resistance) && low_resistance >= 0.0 ||
            throw(ArgumentError("saturated transformer ideal winding resistance must be finite and nonnegative"))
        low_reactance = _saturated_transformer_positive_reactance(
            branch_assembly.series_branch_inductances[ideal_index],
            "saturated transformer ideal winding reactance",
        )
        primary_turns = _saturated_transformer_positive_reactance(
            branch_assembly.series_branch_turns[primary_index],
            "saturated transformer primary turns",
        )
        ideal_turns = _saturated_transformer_positive_reactance(
            branch_assembly.series_branch_turns[ideal_index],
            "saturated transformer ideal winding turns",
        )
        turns_ratio = primary_turns / ideal_turns
        internal_reactance = low_reactance * turns_ratio^2
        low_admittance = inv(complex(low_resistance, low_reactance))
        push!(ideal_low_from_nodes, branch_assembly.winding_from_node_indices[ideal_index])
        push!(ideal_low_to_nodes, branch_assembly.winding_terminal_node_indices[ideal_index])
        push!(ideal_internal_from_nodes, branch_assembly.winding_internal_top_node_indices[ideal_index])
        push!(ideal_internal_to_nodes, branch_assembly.winding_terminal_node_indices[primary_index])
        push!(ideal_low_admittances, low_admittance)
        push!(ideal_mutual_admittances, -low_admittance / turns_ratio)
        push!(ideal_internal_admittances, low_admittance / turns_ratio^2)
        push!(ideal_primary_storage_coefficients, xunits / low_reactance)
        push!(ideal_mutual_storage_coefficients, -xunits / (turns_ratio * low_reactance))
        push!(ideal_secondary_storage_coefficients, xunits / internal_reactance)
        push!(ideal_resistive_storage_coefficients,
              -low_resistance * xunits / low_reactance)
        push!(ideal_low_resistances, low_resistance)
        push!(ideal_low_reactances, low_reactance)
        push!(ideal_internal_reactances, internal_reactance)
        push!(ideal_primary_turns, primary_turns)
        push!(ideal_turns_ratios, turns_ratio)
        push!(ideal_transformer_names, transformer_name)
        handled_windings[ideal_index] = true
    end

    for row_index in eachindex(arrays.winding_transformer_names)
        handled_windings[row_index] && continue
        reactance = _saturated_transformer_positive_reactance(
            branch_assembly.series_branch_inductances[row_index],
            "saturated transformer leakage reactance",
        )
        push!(scalar_from_nodes, branch_assembly.series_branch_from_node_indices[row_index])
        push!(scalar_to_nodes, branch_assembly.series_branch_to_node_indices[row_index])
        push!(
            scalar_admittances,
            _saturated_transformer_branch_admittance(
                branch_assembly.series_branch_resistances[row_index],
                reactance,
                "saturated transformer leakage reactance",
            ),
        )
        push!(scalar_winding_numbers, Int(arrays.winding_numbers[row_index]))
        push!(scalar_reactances, reactance)
    end

    magnetizing_admittances = ComplexF64[]
    for resistance in branch_assembly.magnetizing_branch_resistances
        r = Float64(resistance)
        r != 0.0 ||
            throw(ArgumentError("saturated transformer magnetizing resistance must be nonzero"))
        push!(magnetizing_admittances, complex(inv(r), 0.0))
    end

    return (
        source = :saturated_transformer_steady_state_branch_admittance,
        outcome = :steady_state_admittance,
        transformer_names = copy(arrays.transformer_names),
        internal_top_node_names = copy(branch_assembly.internal_top_node_names),
        internal_top_node_indices = copy(branch_assembly.internal_top_node_indices),
        scalar_branch_from_node_indices = scalar_from_nodes,
        scalar_branch_to_node_indices = scalar_to_nodes,
        scalar_branch_admittances = scalar_admittances,
        scalar_branch_winding_numbers = scalar_winding_numbers,
        scalar_branch_reactances = scalar_reactances,
        ideal_transformer_names = ideal_transformer_names,
        ideal_low_branch_from_node_indices = ideal_low_from_nodes,
        ideal_low_branch_to_node_indices = ideal_low_to_nodes,
        ideal_internal_branch_from_node_indices = ideal_internal_from_nodes,
        ideal_internal_branch_to_node_indices = ideal_internal_to_nodes,
        ideal_low_branch_admittances = ideal_low_admittances,
        ideal_mutual_branch_admittances = ideal_mutual_admittances,
        ideal_internal_branch_admittances = ideal_internal_admittances,
        ideal_primary_storage_coefficients = ideal_primary_storage_coefficients,
        ideal_mutual_storage_coefficients = ideal_mutual_storage_coefficients,
        ideal_secondary_storage_coefficients = ideal_secondary_storage_coefficients,
        ideal_resistive_storage_coefficients = ideal_resistive_storage_coefficients,
        ideal_low_branch_resistances = ideal_low_resistances,
        ideal_low_branch_reactances = ideal_low_reactances,
        ideal_internal_branch_reactances = ideal_internal_reactances,
        ideal_primary_turns = ideal_primary_turns,
        ideal_turns_ratios = ideal_turns_ratios,
        linearized_nonlinear_branch_from_node_indices = linearized_nonlinear_from_nodes,
        linearized_nonlinear_branch_to_node_indices = linearized_nonlinear_to_nodes,
        linearized_nonlinear_branch_admittances = linearized_nonlinear_admittances,
        linearized_nonlinear_branch_reactances = linearized_nonlinear_reactances,
        magnetizing_branch_from_node_indices =
            copy(branch_assembly.magnetizing_branch_from_node_indices),
        magnetizing_branch_to_node_indices =
            copy(branch_assembly.magnetizing_branch_to_node_indices),
        magnetizing_branch_admittances = magnetizing_admittances,
        scalar_branch_count = length(scalar_admittances),
        ideal_branch_pair_count = length(ideal_low_admittances),
        ideal_branch_count = 2 * length(ideal_low_admittances),
        linearized_nonlinear_branch_count = length(linearized_nonlinear_admittances),
        magnetizing_branch_count = length(magnetizing_admittances),
        branch_count =
            length(scalar_admittances) +
            2 * length(ideal_low_admittances) +
            length(linearized_nonlinear_admittances) +
            length(magnetizing_admittances),
        reactance_units = xunits,
        fortran_files = (:OVER2_FOR, :OVER10_FOR),
        fortran_labels = (4011, 4056, 4031, 4033, 4034, 4047, 4049, 74035, 3300, 3504),
        mutation_order = (
            :primary_leakage_branch_storage,
            :ideal_winding_storage,
            :reference_copy_branch_storage,
            :linearized_nonlinear_magnetizing_inductance,
            :positive_kodebr_steady_state_admittance,
            :magnetizing_branch_storage,
        ),
        replacement_ready = false,
    )
end

function saturated_transformer_steady_state_branch_admittance(parsed_intake, node_map; kwargs...)
    return saturated_transformer_steady_state_branch_admittance(
        saturated_transformer_nonlinear_arrays(parsed_intake),
        node_map;
        kwargs...,
    )
end

function saturated_transformer_nonlinear_current_config(
    arrays::SaturatedTransformerNonlinearArrays,
    nonlinear_from_nodes::AbstractVector{<:Integer},
    nonlinear_to_nodes::AbstractVector{<:Integer};
    delta2::Real=1.0,
    nonlinear_admittance_nodes::AbstractVector{<:Integer}=arrays.table_start_indices,
    nonlinear_table_end_indices::AbstractVector{<:Integer}=arrays.table_end_indices,
    nonlinear_subsystem_indices::AbstractVector{<:Integer}=
        ones(Int, length(arrays.nonlinear_types)),
    subsystem_begin_indices::AbstractVector{<:Integer}=Int[1],
    subsystem_owner_rows::AbstractVector{<:Integer}=Int[0],
    subsystem_simultaneous_flags::AbstractVector{<:Integer}=Int[0],
    saturated_transformer_sparse_config::Union{Nothing,NamedTuple}=nothing,
)
    count = length(arrays.nonlinear_types)
    count > 0 || throw(ArgumentError("saturated transformer config requires at least one seed"))
    table = saturated_transformer_current_compensation_table(arrays; delta2 = delta2)
    length(subsystem_owner_rows) == length(subsystem_simultaneous_flags) ||
        throw(ArgumentError("subsystem owner rows and flags lengths must match"))
    base_config = (
        nonlinear_types = Int.(arrays.nonlinear_types),
        nonlinear_from_nodes = _saturated_transformer_int_vector(
            nonlinear_from_nodes,
            "nonlinear_from_nodes",
            count,
        ),
        nonlinear_to_nodes = _saturated_transformer_int_vector(
            nonlinear_to_nodes,
            "nonlinear_to_nodes",
            count,
        ),
        nonlinear_admittance_nodes = _saturated_transformer_int_vector(
            nonlinear_admittance_nodes,
            "nonlinear_admittance_nodes",
            count,
        ),
        nonlinear_table_end_indices = _saturated_transformer_int_vector(
            nonlinear_table_end_indices,
            "nonlinear_table_end_indices",
            count,
        ),
        nonlinear_subsystem_indices = _saturated_transformer_int_vector(
            nonlinear_subsystem_indices,
            "nonlinear_subsystem_indices",
            count,
        ),
        subsystem_begin_indices = Int.(subsystem_begin_indices),
        subsystem_owner_rows = Int.(subsystem_owner_rows),
        subsystem_simultaneous_flags = Int.(subsystem_simultaneous_flags),
        cchar = table.cchar,
        vchar = table.vchar,
        gslope = table.gslope,
        nonlinear_current_segments = Float64.(arrays.current_seeds),
        nonlinear_steady_state_current_values = Float64.(arrays.steady_state_currents),
        delta2 = table.delta2,
    )
    return saturated_transformer_sparse_config === nothing ?
        base_config :
        merge(
            base_config,
            (saturated_transformer_sparse_config = saturated_transformer_sparse_config,),
        )
end

function saturated_transformer_winding_current_config(
    arrays::SaturatedTransformerNonlinearArrays,
    internal_top_nodes::AbstractVector{<:Integer},
    node_map::AbstractDict{Symbol,<:Integer};
    winding_number::Integer=1,
    reference_node_index::Integer=0,
    kwargs...,
)
    count = length(arrays.nonlinear_types)
    top_nodes = _saturated_transformer_int_vector(
        internal_top_nodes,
        "internal_top_nodes",
        count,
    )
    any(node -> node < 1, top_nodes) &&
        throw(ArgumentError("saturated transformer internal_top_nodes must be positive"))
    terminal_binding = saturated_transformer_winding_terminal_nodes(
        arrays,
        node_map;
        winding_number = winding_number,
        reference_node_index = reference_node_index,
    )
    base_config = saturated_transformer_nonlinear_current_config(
        arrays,
        top_nodes,
        terminal_binding.terminal_node_indices;
        kwargs...,
    )
    return merge(
        base_config,
        (
            source = :saturated_transformer_winding_current_config,
            outcome = :nonlinear_current_config,
            nonlinear_top_nodes = top_nodes,
            nonlinear_terminal_node_names = terminal_binding.terminal_node_names,
            nonlinear_terminal_node_indices = terminal_binding.terminal_node_indices,
            nonlinear_terminal_reference_count = terminal_binding.reference_terminal_count,
            nonlinear_winding_numbers = terminal_binding.winding_numbers,
            fortran_files = (:OVER2_FOR,),
            fortran_labels = (4056,),
            mutation_order = (
                :terminal_binding,
                :current_characteristic_arrays,
                :accepted_current_config,
            ),
            replacement_ready = false,
        ),
    )
end

function saturated_transformer_winding_current_config(
    parsed_intake,
    internal_top_nodes::AbstractVector{<:Integer},
    node_map::AbstractDict{Symbol,<:Integer};
    kwargs...,
)
    return saturated_transformer_winding_current_config(
        saturated_transformer_nonlinear_arrays(parsed_intake),
        internal_top_nodes,
        node_map;
        kwargs...,
    )
end

function saturated_transformer_nonlinear_current_config(
    parsed_intake,
    nonlinear_from_nodes::AbstractVector{<:Integer},
    nonlinear_to_nodes::AbstractVector{<:Integer};
    kwargs...,
)
    return saturated_transformer_nonlinear_current_config(
        saturated_transformer_nonlinear_arrays(parsed_intake),
        nonlinear_from_nodes,
        nonlinear_to_nodes;
        kwargs...,
    )
end

function saturated_transformer_nonlinear_slope_branches(
    nonlinear_current_config::NamedTuple,
)
    nonlinear_types = Int.(get(nonlinear_current_config, :nonlinear_types, Int[]))
    count = length(nonlinear_types)
    if count == 0
        return (
            source = :saturated_transformer_nonlinear_slope_branches,
            outcome = :base_nonlinear_slope_admittance,
            elements = Any[],
            element_names = Symbol[],
            branch_count = 0,
            from_nodes = Int[],
            to_nodes = Int[],
            table_indices = Int[],
            conductances = Float64[],
            mutation_order = (:segment_slope_lookup, :base_admittance_stamp),
            replacement_ready = false,
        )
    end
    from_nodes = _saturated_transformer_int_vector(
        get(nonlinear_current_config, :nonlinear_from_nodes, Int[]),
        "nonlinear_from_nodes",
        count,
    )
    to_nodes = _saturated_transformer_int_vector(
        get(nonlinear_current_config, :nonlinear_to_nodes, Int[]),
        "nonlinear_to_nodes",
        count,
    )
    table_start_indices = _saturated_transformer_int_vector(
        get(nonlinear_current_config, :nonlinear_admittance_nodes, Int[]),
        "nonlinear_admittance_nodes",
        count,
    )
    table_end_indices = _saturated_transformer_int_vector(
        get(nonlinear_current_config, :nonlinear_table_end_indices, Int[]),
        "nonlinear_table_end_indices",
        count,
    )
    current_segments = Float64.(
        get(
            nonlinear_current_config,
            :nonlinear_current_segments,
            ones(Float64, count),
        ),
    )
    length(current_segments) == count ||
        throw(ArgumentError("nonlinear_current_segments length must match nonlinear_types"))
    gslope = Float64.(get(nonlinear_current_config, :gslope, Float64[]))

    elements = Any[]
    element_names = Symbol[]
    stamped_from_nodes = Int[]
    stamped_to_nodes = Int[]
    stamped_table_indices = Int[]
    conductances = Float64[]
    for index in eachindex(nonlinear_types)
        _is_pseudo_nonlinear_inductor_type(nonlinear_types[index]) || continue
        segment = Int(abs(current_segments[index]))
        segment > 0 ||
            throw(ArgumentError("pseudo-nonlinear inductor current segment must be nonzero"))
        table_index = table_start_indices[index] + segment - 1
        table_start_indices[index] <= table_index <= table_end_indices[index] ||
            throw(ArgumentError("pseudo-nonlinear inductor current segment must address gslope table"))
        table_index <= length(gslope) ||
            throw(ArgumentError("gslope must cover pseudo-nonlinear inductor slope branches"))
        conductance = gslope[table_index]
        conductance == 0.0 && continue
        to_node = abs(to_nodes[index])
        push!(
            elements,
            SaturatedTransformerNonlinearSlopeBranch(
                from_nodes[index],
                to_node,
                conductance,
            ),
        )
        push!(element_names, Symbol("nonlinear_magnetizing_slope_", string(index)))
        push!(stamped_from_nodes, from_nodes[index])
        push!(stamped_to_nodes, to_node)
        push!(stamped_table_indices, table_index)
        push!(conductances, conductance)
    end
    return (
        source = :saturated_transformer_nonlinear_slope_branches,
        outcome = :base_nonlinear_slope_admittance,
        elements = elements,
        element_names = element_names,
        branch_count = length(elements),
        from_nodes = stamped_from_nodes,
        to_nodes = stamped_to_nodes,
        table_indices = stamped_table_indices,
        conductances = conductances,
        fortran_files = (:OVER2_FOR, :OVER16_FOR),
        fortran_routines = (:SUBTS3,),
        fortran_labels = (73430, 73436, 3950, 4113, 4114),
        mutation_order = (:segment_slope_lookup, :base_admittance_stamp),
        replacement_ready = false,
    )
end

function saturated_transformer_segment_update(
    current_segment::Real,
    companion_current::Real,
    stored_voltage::Real,
    branch_voltage::Real,
    table_start_index::Int,
    table_end_index::Int,
    cchar::AbstractVector{<:Real},
    gslope::AbstractVector{<:Real},
    vchar::AbstractVector{<:Real};
    delta2::Real=1.0,
    voltage_tolerance::Real=0.0,
)
    table_start_index >= 1 ||
        throw(ArgumentError("table_start_index must be positive"))
    table_start_index <= table_end_index ||
        throw(ArgumentError("table_start_index must not exceed table_end_index"))
    table_end_index <= length(cchar) ||
        throw(ArgumentError("saturated transformer table range must address cchar"))
    table_end_index <= length(gslope) ||
        throw(ArgumentError("saturated transformer table range must address gslope"))
    table_end_index <= length(vchar) ||
        throw(ArgumentError("saturated transformer table range must address vchar"))
    _check_over16_finite_vector("cchar", cchar)
    _check_over16_finite_vector("gslope", gslope)
    _check_over16_finite_vector("vchar", vchar)

    previous_current = Float64(current_segment)
    previous_companion = Float64(companion_current)
    previous_stored = Float64(stored_voltage)
    voltage = Float64(branch_voltage)
    delta = Float64(delta2)
    zero = Float64(voltage_tolerance)
    isfinite(previous_current) || throw(ArgumentError("current_segment must be finite"))
    isfinite(previous_companion) || throw(ArgumentError("companion_current must be finite"))
    isfinite(previous_stored) || throw(ArgumentError("stored_voltage must be finite"))
    isfinite(voltage) || throw(ArgumentError("branch_voltage must be finite"))
    isfinite(delta) && delta > 0.0 || throw(ArgumentError("delta2 must be finite and positive"))
    isfinite(zero) && zero >= 0.0 ||
        throw(ArgumentError("voltage_tolerance must be finite and nonnegative"))

    absolute_segment = Int(abs(previous_current))
    absolute_segment > 0 ||
        throw(ArgumentError("saturated transformer current segment must be nonzero"))
    table_index = table_start_index + absolute_segment - 1
    table_start_index <= table_index <= table_end_index ||
        throw(ArgumentError("saturated transformer current segment must lie inside the table"))

    current = previous_current
    voltage_increment = delta * voltage
    equivalent_current = voltage * Float64(gslope[table_index]) + previous_companion
    effective_voltage = previous_stored + voltage_increment
    admittance_delta = 0.0
    segment_changed = false
    polarity_warning = false
    polarity_reversed = false

    polarity_reversal =
        zero == 0.0 ?
        current * effective_voltage <= 0.0 :
        current * effective_voltage < -zero
    if polarity_reversal
        if absolute_segment == 1
            current = -current
            polarity_reversed = true
        else
            polarity_warning = true
            next_index = table_index - 1
            if current > 0.0
                current -= 1.0
            else
                current += 1.0
            end
            admittance_delta = Float64(gslope[next_index]) - Float64(gslope[table_index])
            signed_intercept = Float64(cchar[next_index])
            effective_voltage < 0.0 && (signed_intercept = -signed_intercept)
            target_companion =
                (effective_voltage + voltage_increment - signed_intercept) *
                Float64(gslope[next_index]) / delta
            equivalent_current = (previous_companion + target_companion) / 2.0
            table_index = next_index
            segment_changed = true
        end
    elseif table_index < table_end_index &&
            abs(effective_voltage) > Float64(vchar[table_index])
        next_index = table_index + 1
        if current < 0.0
            current -= 1.0
        else
            current += 1.0
        end
        admittance_delta = Float64(gslope[next_index]) - Float64(gslope[table_index])
        signed_intercept = Float64(cchar[next_index])
        effective_voltage < 0.0 && (signed_intercept = -signed_intercept)
        target_companion =
            (effective_voltage + voltage_increment - signed_intercept) *
            Float64(gslope[next_index]) / delta
        equivalent_current = (previous_companion + target_companion) / 2.0
        table_index = next_index
        segment_changed = true
    elseif absolute_segment > 1
        previous_index = table_index - 1
        if abs(effective_voltage) < Float64(vchar[previous_index])
            if current < 0.0
                current += 1.0
            else
                current -= 1.0
            end
            admittance_delta = Float64(gslope[previous_index]) - Float64(gslope[table_index])
            signed_intercept = Float64(cchar[previous_index])
            effective_voltage < 0.0 && (signed_intercept = -signed_intercept)
            target_companion =
                (effective_voltage + voltage_increment - signed_intercept) *
                Float64(gslope[previous_index]) / delta
            equivalent_current = (previous_companion + target_companion) / 2.0
            table_index = previous_index
            segment_changed = true
        end
    end

    companion_delta = 2.0 * (equivalent_current - previous_companion)
    new_stored_voltage = effective_voltage + voltage_increment
    new_companion = previous_companion + companion_delta
    return SaturatedTransformerSegmentUpdate(
        :saturated_transformer_segment_update,
        :state_mutation,
        (:OVER16_FOR,),
        (:SUBTS3,),
        OVER16_SATURATED_TRANSFORMER_SEGMENT_LABELS,
        previous_current,
        current,
        previous_companion,
        new_companion,
        previous_stored,
        new_stored_voltage,
        voltage,
        voltage_increment,
        effective_voltage,
        companion_delta,
        equivalent_current,
        admittance_delta,
        table_start_index + absolute_segment - 1,
        table_index,
        segment_changed,
        polarity_warning,
        polarity_reversed,
        -companion_delta,
        companion_delta,
        (
            :branch_voltage,
            :segment_threshold_test,
            :companion_current_update,
            :stored_voltage_update,
            :finit_current_injection,
        ),
        segment_changed ?
            [:sparse_admittance_restamp, :full_last14_card_execution, :bulk_last14_oracle] :
            [:full_last14_card_execution, :bulk_last14_oracle],
        false,
    )
end

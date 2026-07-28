abstract type AbstractSourceSignalProvider end

struct IdentitySourceSignalProvider <: AbstractSourceSignalProvider end

"""
    TabulatedSourceSignalProvider(times_s, values; extrapolation=:hold)
    TabulatedSourceSignalProvider(path; extrapolation=:hold)

Own ten ordered source-signal slots as a strictly increasing time table. Values
between samples are linearly interpolated. The file form accepts a comma- or
tab-delimited header `time_s,source_1,...,source_10`; relative paths are resolved
by the caller before construction.
"""
struct TabulatedSourceSignalProvider <: AbstractSourceSignalProvider
    times_s::Vector{Float64}
    values::Matrix{Float64}
    extrapolation::Symbol
    source_path::Union{Nothing,String}

    function TabulatedSourceSignalProvider(
        times_s::AbstractVector{<:Real},
        values::AbstractMatrix{<:Real};
        extrapolation::Symbol = :hold,
        source_path::Union{Nothing,AbstractString} = nothing,
    )
        times = Float64.(times_s)
        samples = Matrix{Float64}(values)
        size(samples, 1) == 10 ||
            throw(ArgumentError("source signal values must have ten rows"))
        size(samples, 2) == length(times) ||
            throw(ArgumentError("source signal times and sample columns must have the same length"))
        isempty(times) && throw(ArgumentError("source signal table must contain at least one sample"))
        all(diff(times) .> 0.0) ||
            throw(ArgumentError("source signal times must be strictly increasing"))
        extrapolation in (:hold, :zero, :error) ||
            throw(ArgumentError("source signal extrapolation must be :hold, :zero, or :error"))
        all(isfinite, times) || throw(ArgumentError("source signal times must be finite"))
        all(isfinite, samples) || throw(ArgumentError("source signal values must be finite"))
        provenance = source_path === nothing ? nothing : abspath(String(source_path))
        return new(times, samples, extrapolation, provenance)
    end
end

function TabulatedSourceSignalProvider(
    times_s::AbstractVector{<:Real},
    values::AbstractVector{<:AbstractVector{<:Real}};
    extrapolation::Symbol = :hold,
)
    all(row -> length(row) == 10, values) ||
        throw(ArgumentError("each source signal sample must contain ten values"))
    samples = isempty(values) ?
        zeros(Float64, 10, 0) :
        hcat((Float64.(row) for row in values)...)
    return TabulatedSourceSignalProvider(times_s, samples; extrapolation = extrapolation)
end

function _source_signal_file_records(path::AbstractString)
    isfile(path) || throw(ArgumentError("source signal file does not exist: $path"))
    return [
        strip(line) for line in readlines(path)
        if !isempty(strip(line)) && !startswith(strip(line), '#')
    ]
end

function _source_signal_file_fields(line::AbstractString, delimiter::Char)
    return strip.(split(String(line), delimiter; keepempty = true))
end

function TabulatedSourceSignalProvider(
    path::AbstractString;
    extrapolation::Symbol = :hold,
)
    records = _source_signal_file_records(path)
    length(records) >= 2 ||
        throw(ArgumentError("source signal file must contain a header and at least one sample"))
    delimiter = occursin('\t', records[1]) ? '\t' : ','
    header = _source_signal_file_fields(records[1], delimiter)
    expected_header = ["time_s"; ["source_$index" for index in 1:10]]
    header == expected_header || throw(ArgumentError(
        "source signal header must be time_s followed by source_1 through source_10",
    ))
    times = Float64[]
    samples = Vector{Vector{Float64}}()
    for (record_index, record) in enumerate(records[2:end])
        fields = _source_signal_file_fields(record, delimiter)
        length(fields) == 11 || throw(ArgumentError(
            "source signal file record $(record_index + 1) must contain eleven fields",
        ))
        values = try
            parse.(Float64, fields)
        catch error
            throw(ArgumentError(
                "source signal file record $(record_index + 1) contains a nonnumeric field: " *
                sprint(showerror, error),
            ))
        end
        push!(times, values[1])
        push!(samples, values[2:11])
    end
    return TabulatedSourceSignalProvider(
        times,
        hcat(samples...);
        extrapolation = extrapolation,
        source_path = path,
    )
end

"""
    AnalyticSourceSlot(slot, signal; assignment=:replace, unit=:V)

Bind one typed analytic source equation to an ordered source-card slot.
`assignment=:if_zero` preserves a nonzero card, interpolation, or TACS value.
The unit is explicit report metadata and must be `:V` or `:A`.
"""
struct AnalyticSourceSlot
    slot::Int
    signal::AnalyticSourceSignal
    assignment::Symbol
    unit::Symbol

    function AnalyticSourceSlot(
        slot::Integer,
        signal::AnalyticSourceSignal;
        assignment::Symbol = :replace,
        unit::Symbol = :V,
    )
        source_slot = Int(slot)
        1 <= source_slot <= 10 ||
            throw(ArgumentError("analytic source slot must be between 1 and 10"))
        assignment in (:replace, :if_zero) ||
            throw(ArgumentError("analytic source assignment must be :replace or :if_zero"))
        unit in (:V, :A) ||
            throw(ArgumentError("analytic source slot unit must be :V or :A"))
        return new(source_slot, signal, assignment, unit)
    end
end

"""
    SourceSignalProgram(interpolation, analytic_slots)

Compose source-table interpolation before TACS feedback and typed analytic slot
assignments after TACS feedback.
"""
struct SourceSignalProgram{P<:AbstractSourceSignalProvider} <:
       AbstractSourceSignalProvider
    interpolation::P
    analytic_slots::Vector{AnalyticSourceSlot}

    function SourceSignalProgram(
        interpolation::P,
        analytic_slots::AbstractVector{AnalyticSourceSlot},
    ) where {P<:AbstractSourceSignalProvider}
        slots = AnalyticSourceSlot[slot for slot in analytic_slots]
        indices = [slot.slot for slot in slots]
        allunique(indices) ||
            throw(ArgumentError("analytic source program contains duplicate slots"))
        return new{P}(interpolation, slots)
    end
end

SourceSignalProgram(analytic_slots::AbstractVector{AnalyticSourceSlot}) =
    SourceSignalProgram(IdentitySourceSignalProvider(), analytic_slots)

function source_signal_values(
    ::IdentitySourceSignalProvider,
    input_values::AbstractVector{<:Real},
    ::Real,
)
    length(input_values) == 10 ||
        throw(ArgumentError("source signal input must contain ten values"))
    return Float64.(input_values)
end

function source_signal_values(
    provider::TabulatedSourceSignalProvider,
    input_values::AbstractVector{<:Real},
    time_s::Real,
)
    length(input_values) == 10 ||
        throw(ArgumentError("source signal input must contain ten values"))
    time = Float64(time_s)
    isfinite(time) || throw(ArgumentError("source signal time must be finite"))
    times = provider.times_s
    if time < times[1]
        provider.extrapolation == :error &&
            throw(ArgumentError("source signal time precedes the first sample"))
        return provider.extrapolation == :zero ? zeros(Float64, 10) : copy(provider.values[:, 1])
    elseif time > times[end]
        provider.extrapolation == :error &&
            throw(ArgumentError("source signal time exceeds the last sample"))
        return provider.extrapolation == :zero ? zeros(Float64, 10) : copy(provider.values[:, end])
    end
    upper = searchsortedfirst(times, time)
    (upper == 1 || times[upper] == time) && return copy(provider.values[:, upper])
    lower = upper - 1
    weight = (time - times[lower]) / (times[upper] - times[lower])
    return (1.0 - weight) .* provider.values[:, lower] .+
           weight .* provider.values[:, upper]
end

source_signal_values(
    program::SourceSignalProgram,
    input_values::AbstractVector{<:Real},
    time_s::Real,
) = source_signal_values(program.interpolation, input_values, time_s)

source_signal_interpolation_active(::IdentitySourceSignalProvider) = false
source_signal_interpolation_active(::TabulatedSourceSignalProvider) = true
source_signal_interpolation_active(program::SourceSignalProgram) =
    source_signal_interpolation_active(program.interpolation)

source_signal_analytic_active(::AbstractSourceSignalProvider) = false
source_signal_analytic_active(program::SourceSignalProgram) =
    !isempty(program.analytic_slots)

source_signal_provenance(::IdentitySourceSignalProvider) = nothing
source_signal_provenance(provider::TabulatedSourceSignalProvider) =
    provider.source_path
source_signal_provenance(program::SourceSignalProgram) =
    source_signal_provenance(program.interpolation)

function _source_signal_owner_units(
    node_values::AbstractVector{Int},
    source_types::AbstractVector{Int},
)
    length(node_values) == length(source_types) ||
        throw(ArgumentError("source node/type vectors must have equal length"))
    units = Dict{Int,Symbol}()
    previous_type = 0
    for index in eachindex(node_values, source_types)
        source_type = abs(source_types[index])
        controlled_successor = previous_type == 16
        previous_type = source_types[index]
        (1 <= source_type <= 10 && !controlled_successor) || continue
        unit = node_values[index] < 0 ? :A : :V
        haskey(units, source_type) && units[source_type] != unit &&
            throw(ArgumentError(
                "source signal slot $source_type cannot own both voltage and current units",
            ))
        units[source_type] = unit
    end
    return units
end

function validate_source_signal_program_units(
    provider::AbstractSourceSignalProvider,
    node_values::AbstractVector{Int},
    source_types::AbstractVector{Int},
)
    source_signal_interpolation_active(provider) ||
        return nothing
    _source_signal_owner_units(node_values, source_types)
    return nothing
end

function validate_source_signal_program_units(
    program::SourceSignalProgram,
    node_values::AbstractVector{Int},
    source_types::AbstractVector{Int},
)
    units = _source_signal_owner_units(node_values, source_types)
    for assignment in program.analytic_slots
        haskey(units, assignment.slot) || throw(ArgumentError(
            "analytic source slot $(assignment.slot) has no executable source owner",
        ))
    end
    for assignment in program.analytic_slots
        assignment.unit == units[assignment.slot] || throw(ArgumentError(
            "analytic source slot $(assignment.slot) declares $(assignment.unit) " *
            "but its deck owner requires $(units[assignment.slot])",
        ))
    end
    return nothing
end

function source_signal_analytic_values(
    ::AbstractSourceSignalProvider,
    input_values::AbstractVector{<:Real},
    ::Real,
)
    length(input_values) == 10 ||
        throw(ArgumentError("source signal input must contain ten values"))
    return (
        values = Float64.(input_values),
        assignment_indices = Int[],
        assignment_units = Symbol[],
    )
end

function source_signal_analytic_values(
    program::SourceSignalProgram,
    input_values::AbstractVector{<:Real},
    time_s::Real,
)
    length(input_values) == 10 ||
        throw(ArgumentError("source signal input must contain ten values"))
    time = Float64(time_s)
    isfinite(time) || throw(ArgumentError("source signal time must be finite"))
    values = Float64.(input_values)
    assignment_indices = Int[]
    assignment_units = Symbol[]
    for assignment in program.analytic_slots
        assignment.assignment == :if_zero &&
            !iszero(values[assignment.slot]) &&
            continue
        value = assignment.signal(time)
        isfinite(value) ||
            throw(ArgumentError("analytic source slot $(assignment.slot) became nonfinite"))
        values[assignment.slot] = value
        push!(assignment_indices, assignment.slot)
        push!(assignment_units, assignment.unit)
    end
    return (
        values = values,
        assignment_indices = assignment_indices,
        assignment_units = assignment_units,
    )
end

struct SourceSignalStageSample
    step_index::Int
    time_s::Float64
    card_input_values::Vector{Float64}
    interpolation_output_values::Vector{Float64}
    tacs_output_values::Vector{Float64}
    analytic_output_values::Vector{Float64}
    accepted_values::Vector{Float64}
    card_read::Bool
    interpolation_applied::Bool
    tacs_override_count::Int
    analytic_assignment_indices::Vector{Int}
end

function _source_signal_slot_units(runtime)
    units = fill(:unused, 10)
    previous_type = 0
    for index in eachindex(runtime.plan.source_iform_values)
        source_type = abs(runtime.plan.source_iform_values[index])
        controlled_successor = previous_type == 16
        previous_type = runtime.plan.source_iform_values[index]
        (1 <= source_type <= 10 && !controlled_successor) || continue
        units[source_type] =
            runtime.plan.source_node_values[index] < 0 ? :A : :V
    end
    return units
end

function _source_signal_csv_text(value::AbstractString)
    return "\"" * replace(String(value), "\"" => "\"\"") * "\""
end

"""
    write_source_signal_stage_report(path, context)

Write the accepted reader/interpolation/TACS/analytic source stages in physical
slot order. The report is sourced only from Julia-owned accepted runtime state.
"""
function write_source_signal_stage_report(
    path::AbstractString,
    context,
)
    runtime = context.source_function_runtime
    runtime === nothing &&
        throw(ArgumentError("EMT context has no source-function runtime"))
    mkpath(dirname(abspath(String(path))))
    units = _source_signal_slot_units(runtime)
    provenance = source_signal_provenance(runtime.signal_provider)
    provenance_text = provenance === nothing ? "" : String(provenance)
    open(path, "w") do io
        println(
            io,
            "step_index,time_s,slot,unit,card_input,interpolated,post_tacs," *
            "post_analytic,accepted,card_read,interpolation_applied," *
            "tacs_override_count,analytic_assigned,source_path",
        )
        for sample in runtime.stage_samples
            assignments = Set(sample.analytic_assignment_indices)
            for slot in 1:10
                @printf(
                    io,
                    "%d,%.17g,%d,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%s,%s,%d,%s,%s\n",
                    sample.step_index,
                    sample.time_s,
                    slot,
                    String(units[slot]),
                    sample.card_input_values[slot],
                    sample.interpolation_output_values[slot],
                    sample.tacs_output_values[slot],
                    sample.analytic_output_values[slot],
                    sample.accepted_values[slot],
                    sample.card_read,
                    sample.interpolation_applied,
                    sample.tacs_override_count,
                    slot in assignments,
                    _source_signal_csv_text(provenance_text),
                )
            end
        end
    end
    return String(path)
end

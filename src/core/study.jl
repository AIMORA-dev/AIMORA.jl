module StudyCore

export STUDY_STATUSES,
       RESULT_STATUSES,
       StudyDescriptor,
       StudyRunRequest,
       ResultQuantity,
       StudyAssumption,
       StudyWarning,
       StudyResult,
       result_quantity,
       study_assumption,
       study_warning,
       study_result,
       add_quantity!,
       add_assumption!,
       add_warning!,
       result_ok,
       has_warnings,
       study_not_implemented

const STUDY_STATUSES = (:implemented, :prototype, :planned, :legacy_reference)
const RESULT_STATUSES = (:ok, :warning, :failed, :not_implemented, :invalid_input)

struct StudyDescriptor
    id::Symbol
    name::String
    domain::Symbol
    status::Symbol
    source_path::String

    function StudyDescriptor(id::Symbol, name::String, domain::Symbol, status::Symbol, source_path::String)
        status in STUDY_STATUSES || error("Unsupported study status $(status). Allowed statuses: $(join(STUDY_STATUSES, ", ")).")
        return new(id, name, domain, status, source_path)
    end
end

Base.@kwdef struct StudyRunRequest
    study::Symbol
    case_path::Union{Nothing,String} = nothing
    output_dir::Union{Nothing,String} = nothing
end

struct ResultQuantity
    key::Symbol
    value::Any
    unit::String
    base::Union{Nothing,String}
    description::String
end

function ResultQuantity(
    key::Symbol,
    value;
    unit::AbstractString,
    base::Union{Nothing,AbstractString} = nothing,
    description::AbstractString = "",
)
    return ResultQuantity(
        key,
        value,
        String(unit),
        base === nothing ? nothing : String(base),
        String(description),
    )
end

struct StudyAssumption
    key::Symbol
    value::Any
    description::String
end

function StudyAssumption(key::Symbol, value; description::AbstractString = "")
    return StudyAssumption(key, value, String(description))
end

struct StudyWarning
    code::Symbol
    severity::Symbol
    message::String
end

function StudyWarning(code::Symbol, message::AbstractString; severity::Symbol = :warning)
    severity in (:info, :warning, :error) || error("Unsupported warning severity $(severity).")
    return StudyWarning(code, severity, String(message))
end

Base.@kwdef mutable struct StudyResult
    study::Symbol
    status::Symbol
    quantities::Dict{Symbol,ResultQuantity} = Dict{Symbol,ResultQuantity}()
    assumptions::Vector{StudyAssumption} = StudyAssumption[]
    warnings::Vector{StudyWarning} = StudyWarning[]
    metadata::Dict{Symbol,Any} = Dict{Symbol,Any}()

    function StudyResult(
        study::Symbol,
        status::Symbol,
        quantities::Dict{Symbol,ResultQuantity},
        assumptions::Vector{StudyAssumption},
        warnings::Vector{StudyWarning},
        metadata::Dict{Symbol,Any},
    )
        status in RESULT_STATUSES || error("Unsupported result status $(status). Allowed statuses: $(join(RESULT_STATUSES, ", ")).")
        if status == :ok && any(w -> w.severity != :info, warnings)
            status = :warning
        end
        return new(study, status, quantities, assumptions, warnings, metadata)
    end
end

result_quantity(args...; kwargs...) = ResultQuantity(args...; kwargs...)
study_assumption(args...; kwargs...) = StudyAssumption(args...; kwargs...)
study_warning(args...; kwargs...) = StudyWarning(args...; kwargs...)

function study_result(
    study::Symbol;
    status::Symbol = :ok,
    quantities = ResultQuantity[],
    assumptions = StudyAssumption[],
    warnings = StudyWarning[],
    metadata = Dict{Symbol,Any}(),
)
    quantity_map = Dict{Symbol,ResultQuantity}()
    for quantity in quantities
        quantity_map[quantity.key] = quantity
    end
    return StudyResult(
        study = study,
        status = status,
        quantities = quantity_map,
        assumptions = collect(assumptions),
        warnings = collect(warnings),
        metadata = Dict{Symbol,Any}(metadata),
    )
end

function add_quantity!(result::StudyResult, quantity::ResultQuantity)
    result.quantities[quantity.key] = quantity
    return result
end

function add_assumption!(result::StudyResult, assumption::StudyAssumption)
    push!(result.assumptions, assumption)
    return result
end

function add_warning!(result::StudyResult, warning::StudyWarning)
    push!(result.warnings, warning)
    if result.status == :ok && warning.severity != :info
        result.status = :warning
    end
    return result
end

result_ok(result::StudyResult) = result.status in (:ok, :warning)
has_warnings(result::StudyResult) = any(w -> w.severity != :info, result.warnings)

function study_not_implemented(desc::StudyDescriptor)
    error("Study $(desc.id) ($(desc.name)) has status $(desc.status) and is not implemented as a runnable study yet.")
end

end

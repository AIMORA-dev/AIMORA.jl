module ValidationCore

export VALIDATION_SEVERITIES,
       VALIDATION_KINDS,
       ValidationIssue,
       ValidationResult,
       ValidationError,
       validation_issue,
       validation_result,
       add_issue!,
       has_errors,
       has_warnings,
       is_valid,
       issue_message,
       issue_summary,
       assert_valid!,
       missing_data,
       invalid_type,
       invalid_value,
       invalid_unit,
       unknown_field,
       unsupported_assumption,
       not_implemented_issue

const VALIDATION_SEVERITIES = (:info, :warning, :error)
const VALIDATION_KINDS = (
    :missing_data,
    :invalid_type,
    :invalid_value,
    :invalid_unit,
    :unknown_field,
    :unsupported_assumption,
    :not_implemented,
)

struct ValidationIssue
    code::Symbol
    severity::Symbol
    kind::Symbol
    subject::String
    message::String
    context::Dict{Symbol,Any}

    function ValidationIssue(
        code::Symbol,
        severity::Symbol,
        kind::Symbol,
        subject::AbstractString,
        message::AbstractString,
        context::Dict{Symbol,Any},
    )
        severity in VALIDATION_SEVERITIES || error("Unsupported validation severity $(severity).")
        kind in VALIDATION_KINDS || error("Unsupported validation kind $(kind).")
        return new(code, severity, kind, String(subject), String(message), context)
    end
end

function ValidationIssue(
    code::Symbol;
    severity::Symbol = :error,
    kind::Symbol,
    subject::AbstractString = "",
    message::AbstractString,
    context = Dict{Symbol,Any}(),
)
    return ValidationIssue(
        code,
        severity,
        kind,
        subject,
        message,
        Dict{Symbol,Any}(Symbol(key) => value for (key, value) in context),
    )
end

Base.@kwdef mutable struct ValidationResult
    source::String = ""
    issues::Vector{ValidationIssue} = ValidationIssue[]
end

struct ValidationError <: Exception
    result::ValidationResult
end

validation_issue(args...; kwargs...) = ValidationIssue(args...; kwargs...)
validation_result(; source::AbstractString = "") = ValidationResult(source = String(source))

function add_issue!(result::ValidationResult, issue::ValidationIssue)
    push!(result.issues, issue)
    return result
end

has_errors(result::ValidationResult) = any(issue -> issue.severity == :error, result.issues)
has_warnings(result::ValidationResult) = any(issue -> issue.severity == :warning, result.issues)
is_valid(result::ValidationResult) = !has_errors(result)

function issue_message(issue::ValidationIssue)
    isempty(issue.subject) && return issue.message
    return "$(issue.subject): $(issue.message)"
end

function issue_summary(result::ValidationResult)
    isempty(result.issues) && return "No validation issues."
    lines = String[]
    prefix = isempty(result.source) ? "Validation failed" : "Validation failed for $(result.source)"
    push!(lines, prefix * ":")
    for issue in result.issues
        push!(lines, "[$(issue.severity)/$(issue.kind)/$(issue.code)] $(issue_message(issue))")
    end
    return join(lines, "\n")
end

function assert_valid!(result::ValidationResult)
    has_errors(result) && throw(ValidationError(result))
    return true
end

Base.showerror(io::IO, err::ValidationError) = print(io, issue_summary(err.result))

function missing_data(subject::AbstractString, message::AbstractString; code::Symbol = :missing_data, context = Dict{Symbol,Any}())
    return ValidationIssue(code; kind = :missing_data, subject = subject, message = message, context = context)
end

function invalid_type(subject::AbstractString, message::AbstractString; code::Symbol = :invalid_type, context = Dict{Symbol,Any}())
    return ValidationIssue(code; kind = :invalid_type, subject = subject, message = message, context = context)
end

function invalid_value(subject::AbstractString, message::AbstractString; code::Symbol = :invalid_value, context = Dict{Symbol,Any}())
    return ValidationIssue(code; kind = :invalid_value, subject = subject, message = message, context = context)
end

function invalid_unit(subject::AbstractString, message::AbstractString; code::Symbol = :invalid_unit, context = Dict{Symbol,Any}())
    return ValidationIssue(code; kind = :invalid_unit, subject = subject, message = message, context = context)
end

function unknown_field(subject::AbstractString, message::AbstractString; code::Symbol = :unknown_field, context = Dict{Symbol,Any}())
    return ValidationIssue(code; kind = :unknown_field, subject = subject, message = message, context = context)
end

function unsupported_assumption(subject::AbstractString, message::AbstractString; code::Symbol = :unsupported_assumption, context = Dict{Symbol,Any}())
    return ValidationIssue(code; kind = :unsupported_assumption, subject = subject, message = message, context = context)
end

function not_implemented_issue(subject::AbstractString, message::AbstractString; code::Symbol = :not_implemented, context = Dict{Symbol,Any}())
    return ValidationIssue(code; kind = :not_implemented, subject = subject, message = message, context = context)
end

end

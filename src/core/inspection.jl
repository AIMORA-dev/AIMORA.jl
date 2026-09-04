module InspectionCore

export INSPECTION_FIELD_KINDS,
       INSPECTION_SECTION_IDS,
       INSPECTION_SCHEMA_VERSION,
       CanonicalInspectionEdit,
       InspectionChoice,
       InspectionCommitResult,
       InspectionField,
       InspectionIdentity,
       InspectionIssue,
       InspectionRange,
       InspectionSchema,
       InspectionSection,
       InspectionSession,
       InspectionTransaction,
       InspectionUnit,
       InspectionValue,
       commit_inspection!,
       commit_inspection_selection!,
       inspection_document,
       inspection_result_document,
       merge_inspection_documents,
       redo_inspection!,
       redo_inspection_selection!,
       undo_inspection!,
       undo_inspection_selection!,
       validate_inspection_schema

const INSPECTION_SCHEMA_VERSION = "1.0.0"
const INSPECTION_SECTION_IDS = (
    "general",
    "connections_terminals",
    "ratings",
    "equipment_construction",
    "study_facets",
    "curves_tables",
    "controls_protection",
    "scenarios_events",
    "results_freshness",
    "drawing_view",
    "validation_readiness",
    "provenance_audit",
)
const INSPECTION_FIELD_KINDS = (
    :boolean,
    :integer,
    :number,
    :text,
    :choice,
    :reference,
    :table,
    :curve,
)
const INSPECTION_COMMIT_STATUSES = (:accepted, :rejected, :conflict, :unavailable)
const INSPECTION_ISSUE_SEVERITIES = (:info, :warning, :error)
const MAXIMUM_SECTIONS = 64
const MAXIMUM_FIELDS = 4096
const MAXIMUM_CHOICES = 4096
const MAXIMUM_EDITS = 4096
const MAXIMUM_TABLE_ROWS = 100_000
const MAXIMUM_TEXT_BYTES = 1024 * 1024

Base.@kwdef struct InspectionIdentity
    project_id::String
    asset_id::String
    projection_id::String
    view_id::String
    equipment_class::String
    result_bindings::Vector{String} = String[]
end

Base.@kwdef struct InspectionUnit
    symbol::String
    scale_to_canonical::Float64 = 1.0
    offset_to_canonical::Float64 = 0.0
end

Base.@kwdef struct InspectionRange
    lower::Union{Nothing,Float64} = nothing
    upper::Union{Nothing,Float64} = nothing
    lower_inclusive::Bool = true
    upper_inclusive::Bool = true
end

Base.@kwdef struct InspectionChoice
    value::Any
    label::String
end

Base.@kwdef struct InspectionIssue
    path::String = ""
    code::String
    severity::Symbol = :error
    message::String
    owner_id::String = ""
end

Base.@kwdef struct InspectionField
    path::String
    label::String
    kind::Symbol
    dimension::String = ""
    canonical_unit::String = ""
    display_units::Vector{InspectionUnit} = InspectionUnit[]
    choices::Vector{InspectionChoice} = InspectionChoice[]
    range::Union{Nothing,InspectionRange} = nothing
    dependencies::Vector{String} = String[]
    read_only::Bool = false
    required::Bool = false
    uncertainty::Union{Nothing,Float64} = nothing
    provenance::String = ""
    affects_model_paths::Vector{String} = String[]
    affects_view_ids::Vector{String} = String[]
    invalidates_result_ids::Vector{String} = String[]
end

Base.@kwdef struct InspectionSection
    id::String
    title::String
    fields::Vector{InspectionField} = InspectionField[]
    available::Bool = true
    unavailable_reason::String = ""
end

Base.@kwdef struct InspectionSchema
    version::String = INSPECTION_SCHEMA_VERSION
    equipment_class::String
    sections::Vector{InspectionSection}
end

Base.@kwdef struct InspectionValue
    value::Any
    canonical_unit::String = ""
    provenance::String = ""
    issues::Vector{InspectionIssue} = InspectionIssue[]
end

Base.@kwdef struct CanonicalInspectionEdit
    path::String
    value::Any
    display_unit::String = ""
end

Base.@kwdef struct InspectionTransaction
    base_revision::UInt64
    edits::Vector{CanonicalInspectionEdit}
end

Base.@kwdef struct InspectionCommitResult
    status::Symbol
    base_revision::UInt64
    revision::UInt64
    issues::Vector{InspectionIssue} = InspectionIssue[]
    affected_model_paths::Vector{String} = String[]
    affected_view_ids::Vector{String} = String[]
    invalidated_result_ids::Vector{String} = String[]
end

mutable struct InspectionSession
    identity::InspectionIdentity
    schema::InspectionSchema
    values::Dict{String,InspectionValue}
    revision::UInt64
    history::Vector{Dict{String,InspectionValue}}
    future::Vector{Dict{String,InspectionValue}}
end

function _valid_identifier(value::AbstractString)
    return !isempty(value) && ncodeunits(value) <= 512 &&
           occursin(r"^[A-Za-z0-9][A-Za-z0-9._:/-]*$", value)
end

function _validate_identity(identity::InspectionIdentity)
    for (name, value) in (
        ("project", identity.project_id),
        ("asset", identity.asset_id),
        ("projection", identity.projection_id),
        ("view", identity.view_id),
        ("equipment class", identity.equipment_class),
    )
        _valid_identifier(value) || throw(ArgumentError("invalid $(name) inspection identity"))
    end
    length(identity.result_bindings) <= MAXIMUM_FIELDS ||
        throw(ArgumentError("too many inspection result bindings"))
    all(_valid_identifier, identity.result_bindings) ||
        throw(ArgumentError("invalid inspection result binding"))
    return identity
end

function _validate_unit(unit::InspectionUnit)
    _valid_identifier(unit.symbol) || throw(ArgumentError("invalid inspection unit"))
    isfinite(unit.scale_to_canonical) && unit.scale_to_canonical != 0.0 ||
        throw(ArgumentError("inspection unit scale must be finite and nonzero"))
    isfinite(unit.offset_to_canonical) ||
        throw(ArgumentError("inspection unit offset must be finite"))
    return unit
end

function _validate_range(range::InspectionRange)
    range.lower === nothing || isfinite(range.lower) ||
        throw(ArgumentError("inspection range lower bound must be finite"))
    range.upper === nothing || isfinite(range.upper) ||
        throw(ArgumentError("inspection range upper bound must be finite"))
    if range.lower !== nothing && range.upper !== nothing
        range.lower <= range.upper || throw(ArgumentError("inspection range is reversed"))
    end
    return range
end

function validate_inspection_schema(schema::InspectionSchema)
    schema.version == INSPECTION_SCHEMA_VERSION ||
        throw(ArgumentError("unsupported inspection schema version"))
    _valid_identifier(schema.equipment_class) ||
        throw(ArgumentError("invalid inspection equipment class"))
    0 < length(schema.sections) <= MAXIMUM_SECTIONS ||
        throw(ArgumentError("inspection schema section count is invalid"))
    section_ids = String[]
    field_paths = String[]
    for section in schema.sections
        section.id in INSPECTION_SECTION_IDS ||
            throw(ArgumentError("unknown inspection section $(section.id)"))
        push!(section_ids, section.id)
        !isempty(strip(section.title)) || throw(ArgumentError("inspection section title is empty"))
        if section.available
            isempty(section.unavailable_reason) ||
                throw(ArgumentError("available inspection section has an unavailable reason"))
        else
            isempty(section.fields) ||
                throw(ArgumentError("unavailable inspection section contains fields"))
            !isempty(strip(section.unavailable_reason)) ||
                throw(ArgumentError("unavailable inspection section has no reason"))
        end
        for field in section.fields
            _valid_identifier(field.path) || throw(ArgumentError("invalid inspection field path"))
            !isempty(strip(field.label)) || throw(ArgumentError("inspection field label is empty"))
            field.kind in INSPECTION_FIELD_KINDS ||
                throw(ArgumentError("unknown inspection field kind $(field.kind)"))
            length(field.display_units) <= MAXIMUM_CHOICES ||
                throw(ArgumentError("too many inspection display units"))
            foreach(_validate_unit, field.display_units)
            unit_symbols = getfield.(field.display_units, :symbol)
            length(unique(unit_symbols)) == length(unit_symbols) ||
                throw(ArgumentError("duplicate inspection display unit"))
            isempty(field.display_units) == isempty(field.canonical_unit) ||
                throw(ArgumentError("inspection canonical/display units are inconsistent"))
            isempty(field.display_units) || field.canonical_unit in unit_symbols ||
                throw(ArgumentError("canonical inspection unit is unavailable"))
            length(field.choices) <= MAXIMUM_CHOICES ||
                throw(ArgumentError("too many inspection choices"))
            all(choice -> !isempty(strip(choice.label)), field.choices) ||
                throw(ArgumentError("inspection choice label is empty"))
            field.kind == :choice || isempty(field.choices) ||
                throw(ArgumentError("non-choice inspection field defines choices"))
            field.kind != :choice || !isempty(field.choices) ||
                throw(ArgumentError("choice inspection field has no choices"))
            field.range === nothing || _validate_range(field.range)
            field.uncertainty === nothing ||
                (isfinite(field.uncertainty) && field.uncertainty >= 0.0) ||
                throw(ArgumentError("inspection uncertainty is invalid"))
            push!(field_paths, field.path)
        end
    end
    length(unique(section_ids)) == length(section_ids) ||
        throw(ArgumentError("duplicate inspection section"))
    length(field_paths) <= MAXIMUM_FIELDS || throw(ArgumentError("too many inspection fields"))
    length(unique(field_paths)) == length(field_paths) ||
        throw(ArgumentError("duplicate inspection field path"))
    path_set = Set(field_paths)
    for section in schema.sections, field in section.fields
        all(path -> path in path_set, field.dependencies) ||
            throw(ArgumentError("unknown inspection field dependency"))
    end
    return schema
end

function _field_map(schema::InspectionSchema)
    return Dict(field.path => field for section in schema.sections for field in section.fields)
end

function _issue(path, code, message; owner_id = "")
    return InspectionIssue(
        path = String(path),
        code = String(code),
        severity = :error,
        message = String(message),
        owner_id = String(owner_id),
    )
end

function _normalize_unit(field::InspectionField, raw_value, display_unit::AbstractString)
    isempty(field.display_units) && return raw_value, "", InspectionIssue[]
    raw_value isa Real && !(raw_value isa Bool) ||
        return raw_value,
        field.canonical_unit,
        [_issue(field.path, "INVALID_TYPE", "unit-bearing field requires a numeric value")]
    unit_name = isempty(display_unit) ? field.canonical_unit : String(display_unit)
    unit_index = findfirst(unit -> unit.symbol == unit_name, field.display_units)
    unit_index === nothing &&
        return raw_value,
        field.canonical_unit,
        [_issue(field.path, "INVALID_UNIT", "display unit is not allowed for this field")]
    unit = field.display_units[unit_index]
    value = Float64(raw_value)
    isfinite(value) ||
        return value,
        field.canonical_unit,
        [_issue(field.path, "NONFINITE_VALUE", "numeric field must be finite")]
    return value * unit.scale_to_canonical + unit.offset_to_canonical,
    field.canonical_unit,
    InspectionIssue[]
end

function _validate_table_rows(field::InspectionField, value::AbstractVector)
    length(value) <= MAXIMUM_TABLE_ROWS ||
        return [_issue(field.path, "TOO_MANY_ROWS", "table or curve exceeds its row limit")]
    row_ids = String[]
    for row in value
        row isa AbstractDict ||
            return [_issue(field.path, "INVALID_ROW", "table or curve row must be an object")]
        row_id = get(row, "row_id", get(row, :row_id, nothing))
        row_id isa AbstractString && _valid_identifier(row_id) ||
            return [_issue(field.path, "INVALID_ROW_ID", "row requires a stable identifier")]
        push!(row_ids, String(row_id))
    end
    length(unique(row_ids)) == length(row_ids) ||
        return [_issue(field.path, "DUPLICATE_ROW_ID", "row identifiers must be unique")]
    return InspectionIssue[]
end

function _normalize_edit(field::InspectionField, edit::CanonicalInspectionEdit)
    field.read_only &&
        return nothing, [_issue(field.path, "READ_ONLY", "field is read-only")]
    value, canonical_unit, unit_issues =
        _normalize_unit(field, edit.value, edit.display_unit)
    isempty(unit_issues) || return nothing, unit_issues
    issues = InspectionIssue[]
    if field.kind == :boolean
        value isa Bool || push!(issues, _issue(field.path, "INVALID_TYPE", "expected boolean"))
    elseif field.kind == :integer
        value isa Integer && !(value isa Bool) ||
            push!(issues, _issue(field.path, "INVALID_TYPE", "expected integer"))
    elseif field.kind == :number
        value isa Real && !(value isa Bool) && isfinite(value) ||
            push!(issues, _issue(field.path, "INVALID_TYPE", "expected finite number"))
        value isa Real && !(value isa Bool) && (value = Float64(value))
    elseif field.kind in (:text, :reference)
        value isa AbstractString ||
            push!(issues, _issue(field.path, "INVALID_TYPE", "expected text"))
        if value isa AbstractString
            ncodeunits(value) <= MAXIMUM_TEXT_BYTES ||
                push!(issues, _issue(field.path, "VALUE_TOO_LARGE", "text value is too large"))
            value = String(value)
        end
    elseif field.kind == :choice
        any(choice -> isequal(choice.value, value), field.choices) ||
            push!(issues, _issue(field.path, "INVALID_CHOICE", "value is not an allowed choice"))
    elseif field.kind in (:table, :curve)
        value isa AbstractVector ||
            push!(issues, _issue(field.path, "INVALID_TYPE", "expected row array"))
        value isa AbstractVector && append!(issues, _validate_table_rows(field, value))
        value isa AbstractVector && (value = deepcopy(collect(value)))
    end
    if isempty(issues) && field.range !== nothing && value isa Real
        range = field.range
        lower_ok = range.lower === nothing ||
                   (range.lower_inclusive ? value >= range.lower : value > range.lower)
        upper_ok = range.upper === nothing ||
                   (range.upper_inclusive ? value <= range.upper : value < range.upper)
        lower_ok && upper_ok ||
            push!(issues, _issue(field.path, "OUT_OF_RANGE", "numeric value is outside its range"))
    end
    isempty(issues) || return nothing, issues
    return InspectionValue(
        value = value,
        canonical_unit = canonical_unit,
        provenance = field.provenance,
    ), issues
end

function InspectionSession(
    identity::InspectionIdentity,
    schema::InspectionSchema,
    values::AbstractDict;
    revision::Integer = 1,
)
    _validate_identity(identity)
    validate_inspection_schema(schema)
    revision > 0 || throw(ArgumentError("inspection revision must be positive"))
    revision <= typemax(UInt64) || throw(ArgumentError("inspection revision exceeds UInt64"))
    fields = _field_map(schema)
    normalized = Dict{String,InspectionValue}()
    for (raw_path, raw_value) in pairs(values)
        path = String(raw_path)
        field = get(fields, path, nothing)
        field === nothing && throw(ArgumentError("unknown initial inspection field $(path)"))
        edit = raw_value isa InspectionValue ?
               CanonicalInspectionEdit(
            path = path,
            value = raw_value.value,
            display_unit = raw_value.canonical_unit,
        ) : CanonicalInspectionEdit(path = path, value = raw_value)
        value, issues = _normalize_edit(
            InspectionField(
                path = field.path,
                label = field.label,
                kind = field.kind,
                dimension = field.dimension,
                canonical_unit = field.canonical_unit,
                display_units = field.display_units,
                choices = field.choices,
                range = field.range,
                dependencies = field.dependencies,
                read_only = false,
                required = field.required,
                uncertainty = field.uncertainty,
                provenance = field.provenance,
                affects_model_paths = field.affects_model_paths,
                affects_view_ids = field.affects_view_ids,
                invalidates_result_ids = field.invalidates_result_ids,
            ),
            edit,
        )
        isempty(issues) || throw(ArgumentError("invalid initial inspection value for $(path)"))
        normalized[path] = value
    end
    required_paths = [
        field.path for section in schema.sections for field in section.fields if field.required
    ]
    all(path -> haskey(normalized, path), required_paths) ||
        throw(ArgumentError("required inspection value is missing"))
    return InspectionSession(
        identity,
        schema,
        normalized,
        UInt64(revision),
        Dict{String,InspectionValue}[],
        Dict{String,InspectionValue}[],
    )
end

function _result(
    status::Symbol,
    base_revision::UInt64,
    revision::UInt64;
    issues = InspectionIssue[],
    affected_model_paths = String[],
    affected_view_ids = String[],
    invalidated_result_ids = String[],
)
    status in INSPECTION_COMMIT_STATUSES || error("invalid inspection commit status")
    return InspectionCommitResult(
        status = status,
        base_revision = base_revision,
        revision = revision,
        issues = collect(issues),
        affected_model_paths = sort!(unique!(collect(affected_model_paths))),
        affected_view_ids = sort!(unique!(collect(affected_view_ids))),
        invalidated_result_ids = sort!(unique!(collect(invalidated_result_ids))),
    )
end

function commit_inspection!(
    session::InspectionSession,
    transaction::InspectionTransaction;
    validator::Function = (_, _) -> InspectionIssue[],
)
    transaction.base_revision == session.revision ||
        return _result(:conflict, transaction.base_revision, session.revision;
                       issues = [_issue("", "REVISION_CONFLICT", "base revision is stale")])
    0 < length(transaction.edits) <= MAXIMUM_EDITS ||
        return _result(:rejected, transaction.base_revision, session.revision;
                       issues = [_issue("", "EDIT_COUNT_INVALID", "edit count is invalid")])
    paths = getfield.(transaction.edits, :path)
    length(unique(paths)) == length(paths) ||
        return _result(:rejected, transaction.base_revision, session.revision;
                       issues = [_issue("", "DUPLICATE_EDIT", "field is edited more than once")])
    fields = _field_map(session.schema)
    next_values = deepcopy(session.values)
    issues = InspectionIssue[]
    affected_model_paths = String[]
    affected_view_ids = String[]
    invalidated_result_ids = String[]
    for edit in transaction.edits
        field = get(fields, edit.path, nothing)
        if field === nothing
            push!(issues, _issue(edit.path, "UNKNOWN_FIELD", "field path is not in the schema"))
            continue
        end
        value, field_issues = _normalize_edit(field, edit)
        append!(issues, field_issues)
        value === nothing && continue
        next_values[edit.path] = value
        append!(affected_model_paths, field.affects_model_paths)
        append!(affected_view_ids, field.affects_view_ids)
        append!(invalidated_result_ids, field.invalidates_result_ids)
    end
    isempty(issues) && append!(issues, validator(session.schema, deepcopy(next_values)))
    all(issue -> issue isa InspectionIssue, issues) ||
        throw(ArgumentError("inspection validator must return InspectionIssue values"))
    all(issue -> issue.severity in INSPECTION_ISSUE_SEVERITIES, issues) ||
        throw(ArgumentError("inspection validator returned an invalid severity"))
    any(issue -> issue.severity == :error, issues) &&
        return _result(:rejected, transaction.base_revision, session.revision; issues = issues)
    push!(session.history, deepcopy(session.values))
    empty!(session.future)
    session.values = next_values
    session.revision += 1
    return _result(
        :accepted,
        transaction.base_revision,
        session.revision;
        issues = issues,
        affected_model_paths = affected_model_paths,
        affected_view_ids = affected_view_ids,
        invalidated_result_ids = invalidated_result_ids,
    )
end

function undo_inspection!(session::InspectionSession, base_revision::Integer)
    requested_revision = try
        UInt64(base_revision)
    catch
        return _result(:conflict, session.revision, session.revision;
                       issues = [_issue("", "REVISION_INVALID", "base revision is invalid")])
    end
    requested_revision == session.revision ||
        return _result(:conflict, requested_revision, session.revision;
                       issues = [_issue("", "REVISION_CONFLICT", "base revision is stale")])
    isempty(session.history) &&
        return _result(:unavailable, session.revision, session.revision;
                       issues = [_issue("", "UNDO_UNAVAILABLE", "no inspection edit can be undone")])
    push!(session.future, deepcopy(session.values))
    session.values = pop!(session.history)
    previous = session.revision
    session.revision += 1
    return _result(:accepted, previous, session.revision)
end

function redo_inspection!(session::InspectionSession, base_revision::Integer)
    requested_revision = try
        UInt64(base_revision)
    catch
        return _result(:conflict, session.revision, session.revision;
                       issues = [_issue("", "REVISION_INVALID", "base revision is invalid")])
    end
    requested_revision == session.revision ||
        return _result(:conflict, requested_revision, session.revision;
                       issues = [_issue("", "REVISION_CONFLICT", "base revision is stale")])
    isempty(session.future) &&
        return _result(:unavailable, session.revision, session.revision;
                       issues = [_issue("", "REDO_UNAVAILABLE", "no inspection edit can be redone")])
    push!(session.history, deepcopy(session.values))
    session.values = pop!(session.future)
    previous = session.revision
    session.revision += 1
    return _result(:accepted, previous, session.revision)
end

function _selection_result(results::AbstractVector{InspectionCommitResult})
    first_result = first(results)
    status = any(result -> result.status == :conflict, results) ? :conflict :
             all(result -> result.status == :unavailable, results) ? :unavailable :
             any(result -> result.status != :accepted, results) ? :rejected : :accepted
    return _result(
        status,
        first_result.base_revision,
        status == :accepted ? first_result.revision : maximum(getfield.(results, :revision));
        issues = reduce(vcat, getfield.(results, :issues); init = InspectionIssue[]),
        affected_model_paths = reduce(
            vcat,
            getfield.(results, :affected_model_paths);
            init = String[],
        ),
        affected_view_ids = reduce(
            vcat,
            getfield.(results, :affected_view_ids);
            init = String[],
        ),
        invalidated_result_ids = reduce(
            vcat,
            getfield.(results, :invalidated_result_ids);
            init = String[],
        ),
    )
end

function _apply_selection_trials!(sessions, trials, results)
    result = _selection_result(results)
    result.status == :accepted || return result
    for (session, trial) in zip(sessions, trials)
        session.values = trial.values
        session.revision = trial.revision
        session.history = trial.history
        session.future = trial.future
    end
    return result
end

function commit_inspection_selection!(
    sessions::AbstractVector{<:InspectionSession},
    transaction::InspectionTransaction;
    validator::Function = (_, _) -> InspectionIssue[],
)
    isempty(sessions) && throw(ArgumentError("inspection selection is empty"))
    revisions = getfield.(sessions, :revision)
    all(==(first(revisions)), revisions) ||
        return _result(:conflict, transaction.base_revision, maximum(revisions);
                       issues = [_issue("", "REVISION_CONFLICT", "selection revisions differ")])
    trials = deepcopy.(sessions)
    results = [
        commit_inspection!(trial, transaction; validator = validator) for trial in trials
    ]
    return _apply_selection_trials!(sessions, trials, results)
end

function _history_inspection_selection!(sessions, base_revision::Integer, operation::Function)
    isempty(sessions) && throw(ArgumentError("inspection selection is empty"))
    trials = deepcopy.(sessions)
    results = [operation(trial, base_revision) for trial in trials]
    return _apply_selection_trials!(sessions, trials, results)
end

undo_inspection_selection!(sessions, base_revision::Integer) =
    _history_inspection_selection!(sessions, base_revision, undo_inspection!)

redo_inspection_selection!(sessions, base_revision::Integer) =
    _history_inspection_selection!(sessions, base_revision, redo_inspection!)

function _transport_value(value)
    if value isa Symbol
        return String(value)
    elseif value isa AbstractDict
        return Dict{String,Any}(String(key) => _transport_value(item) for (key, item) in pairs(value))
    elseif value isa Tuple || value isa AbstractVector
        return Any[_transport_value(item) for item in value]
    elseif value === nothing || value isa Bool || value isa Number || value isa AbstractString
        return value
    end
    throw(ArgumentError("inspection value is not transport-safe: $(typeof(value))"))
end

function _issue_document(issue::InspectionIssue)
    return Dict{String,Any}(
        "path" => issue.path,
        "code" => issue.code,
        "severity" => String(issue.severity),
        "message" => issue.message,
        "owner_id" => issue.owner_id,
    )
end

function _field_document(field::InspectionField)
    return Dict{String,Any}(
        "path" => field.path,
        "label" => field.label,
        "kind" => String(field.kind),
        "dimension" => field.dimension,
        "canonical_unit" => field.canonical_unit,
        "display_units" => Any[
            Dict{String,Any}(
                "symbol" => unit.symbol,
                "scale_to_canonical" => unit.scale_to_canonical,
                "offset_to_canonical" => unit.offset_to_canonical,
            ) for unit in field.display_units
        ],
        "choices" => Any[
            Dict{String,Any}("value" => _transport_value(choice.value), "label" => choice.label)
            for choice in field.choices
        ],
        "range" => field.range === nothing ? nothing : Dict{String,Any}(
            "lower" => field.range.lower,
            "upper" => field.range.upper,
            "lower_inclusive" => field.range.lower_inclusive,
            "upper_inclusive" => field.range.upper_inclusive,
        ),
        "dependencies" => copy(field.dependencies),
        "read_only" => field.read_only,
        "required" => field.required,
        "uncertainty" => field.uncertainty,
        "provenance" => field.provenance,
        "affects_model_paths" => copy(field.affects_model_paths),
        "affects_view_ids" => copy(field.affects_view_ids),
        "invalidates_result_ids" => copy(field.invalidates_result_ids),
    )
end

function inspection_document(session::InspectionSession)
    return Dict{String,Any}(
        "schema_version" => session.schema.version,
        "revision" => string(session.revision),
        "identity" => Dict{String,Any}(
            "project_id" => session.identity.project_id,
            "asset_id" => session.identity.asset_id,
            "projection_id" => session.identity.projection_id,
            "view_id" => session.identity.view_id,
            "equipment_class" => session.identity.equipment_class,
            "result_bindings" => copy(session.identity.result_bindings),
        ),
        "sections" => Any[
            Dict{String,Any}(
                "id" => section.id,
                "title" => section.title,
                "available" => section.available,
                "unavailable_reason" => section.unavailable_reason,
                "fields" => Any[_field_document(field) for field in section.fields],
            ) for section in session.schema.sections
        ],
        "values" => Dict{String,Any}(
            path => Dict{String,Any}(
                "value" => _transport_value(value.value),
                "canonical_unit" => value.canonical_unit,
                "provenance" => value.provenance,
                "issues" => Any[_issue_document(issue) for issue in value.issues],
            ) for (path, value) in session.values
        ),
        "undo_available" => !isempty(session.history),
        "redo_available" => !isempty(session.future),
    )
end

function inspection_result_document(result::InspectionCommitResult)
    return Dict{String,Any}(
        "status" => String(result.status),
        "base_revision" => string(result.base_revision),
        "revision" => string(result.revision),
        "issues" => Any[_issue_document(issue) for issue in result.issues],
        "affected_model_paths" => copy(result.affected_model_paths),
        "affected_view_ids" => copy(result.affected_view_ids),
        "invalidated_result_ids" => copy(result.invalidated_result_ids),
    )
end

function merge_inspection_documents(sessions::AbstractVector{<:InspectionSession})
    isempty(sessions) && throw(ArgumentError("inspection selection is empty"))
    first_session = first(sessions)
    all(session -> session.schema.version == first_session.schema.version, sessions) ||
        throw(ArgumentError("inspection schema versions differ"))
    field_maps = _field_map.(getfield.(sessions, :schema))
    common_paths = reduce(intersect, Set.(keys.(field_maps)))
    fields = Dict{String,Any}()
    for path in sort!(collect(common_paths))
        descriptors = [field_map[path] for field_map in field_maps]
        first_field = first(descriptors)
        all(field -> field.kind == first_field.kind &&
                     field.canonical_unit == first_field.canonical_unit, descriptors) || continue
        values = [session.values[path].value for session in sessions if haskey(session.values, path)]
        length(values) == length(sessions) || continue
        fields[path] = Dict{String,Any}(
            "kind" => String(first_field.kind),
            "canonical_unit" => first_field.canonical_unit,
            "mixed" => !all(value -> isequal(value, first(values)), values),
            "value" => all(value -> isequal(value, first(values)), values) ?
                       _transport_value(first(values)) : nothing,
        )
    end
    document = inspection_document(first_session)
    asset_ids = sort!(getfield.(getfield.(sessions, :identity), :asset_id))
    document["selection_count"] = length(sessions)
    document["asset_ids"] = asset_ids
    document["identity"]["asset_ids"] = asset_ids
    document["values"] = Dict{String,Any}(
        path => Dict{String,Any}(
            "value" => field["value"],
            "canonical_unit" => field["canonical_unit"],
            "mixed" => field["mixed"],
            "provenance" => "",
            "issues" => Any[],
        ) for (path, field) in fields
    )
    for section in document["sections"]
        section["fields"] = Any[
            field for field in section["fields"] if haskey(fields, field["path"])
        ]
    end
    document["undo_available"] = all(!isempty(session.history) for session in sessions)
    document["redo_available"] = all(!isempty(session.future) for session in sessions)
    return document
end

end

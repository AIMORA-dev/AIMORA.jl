module AssetTables

using ..ValidationCore

export ColumnSpec,
       TableSchema,
       normalize_row,
       apply_defaults,
       missing_required_columns,
       unknown_columns,
       validate_row_issues,
       validate_table_issues,
       validate_row,
       validate_table,
       table_columns,
       table_matrix,
       update_rows!,
       format_table

struct ColumnSpec
    key::Symbol
    label::String
    unit::Union{Nothing,String}
    required::Bool
    default::Any
    value_type::Any
    min_value::Any
    max_value::Any
    allowed_values::Union{Nothing,Vector{Any}}
    description::String
end

function ColumnSpec(
    key::Symbol;
    label::AbstractString = String(key),
    unit::Union{Nothing,AbstractString} = nothing,
    required::Bool = false,
    default = missing,
    value_type = nothing,
    min_value = nothing,
    max_value = nothing,
    allowed_values = nothing,
    description::AbstractString = "",
)
    return ColumnSpec(
        key,
        String(label),
        unit === nothing ? nothing : String(unit),
        required,
        default,
        value_type,
        min_value,
        max_value,
        allowed_values === nothing ? nothing : Any[allowed_values...],
        String(description),
    )
end

struct TableSchema
    table::Symbol
    columns::Vector{ColumnSpec}
    notes::String
end

function TableSchema(table::Symbol, columns::Vector{ColumnSpec}; notes::AbstractString = "")
    return TableSchema(table, columns, String(notes))
end

function normalize_row(row)
    if row isa AbstractDict
        return Dict{Symbol,Any}(Symbol(k) => v for (k, v) in row)
    elseif row isa NamedTuple
        return Dict{Symbol,Any}(Symbol(k) => getfield(row, k) for k in keys(row))
    else
        return Dict{Symbol,Any}(Symbol(k) => getproperty(row, k) for k in propertynames(row))
    end
end

function apply_defaults(schema::TableSchema, row::AbstractDict)
    normalized = normalize_row(row)
    for column in schema.columns
        if !haskey(normalized, column.key) && column.default !== missing
            normalized[column.key] = column.default
        end
    end
    return normalized
end

function missing_required_columns(schema::TableSchema, row::AbstractDict)
    normalized = normalize_row(row)
    return [column for column in schema.columns if column.required && !haskey(normalized, column.key)]
end

function unknown_columns(schema::TableSchema, row::AbstractDict)
    normalized = normalize_row(row)
    allowed = Set(column.key for column in schema.columns)
    return [key for key in keys(normalized) if !(key in allowed)]
end

function validate_column_value!(result::ValidationResult, schema::TableSchema, column::ColumnSpec, value)
    subject = "$(schema.table).$(column.key)"
    if ismissing(value)
        if column.required
            add_issue!(
                result,
                missing_data(subject, "Missing value for required column $(column.key)."; context = Dict(:table => schema.table, :column => column.key)),
            )
        end
        return result
    end

    if column.value_type !== nothing && !(value isa column.value_type)
        add_issue!(
            result,
            invalid_type(
                subject,
                "Expected $(column.value_type), got $(typeof(value)).";
                context = Dict(:table => schema.table, :column => column.key, :expected => column.value_type, :actual => typeof(value)),
            ),
        )
    end

    if column.min_value !== nothing && value < column.min_value
        add_issue!(
            result,
            invalid_value(
                subject,
                "$(value) is below minimum $(column.min_value).";
                context = Dict(:table => schema.table, :column => column.key, :value => value, :minimum => column.min_value),
            ),
        )
    end

    if column.max_value !== nothing && value > column.max_value
        add_issue!(
            result,
            invalid_value(
                subject,
                "$(value) is above maximum $(column.max_value).";
                context = Dict(:table => schema.table, :column => column.key, :value => value, :maximum => column.max_value),
            ),
        )
    end

    if column.allowed_values !== nothing && !(value in column.allowed_values)
        add_issue!(
            result,
            invalid_value(
                subject,
                "$(value) is not one of $(join(column.allowed_values, ", ")).";
                context = Dict(:table => schema.table, :column => column.key, :value => value, :allowed_values => column.allowed_values),
            ),
        )
    end

    return result
end

function validate_row_issues(schema::TableSchema, row::AbstractDict; strict::Bool = false)
    result = validation_result(source = "table $(schema.table)")
    normalized = normalize_row(row)
    missing_columns = missing_required_columns(schema, row)
    if !isempty(missing_columns)
        for column in missing_columns
            add_issue!(
                result,
                missing_data(
                    "$(schema.table).$(column.key)",
                    "Missing required column $(column.key).";
                    context = Dict(:table => schema.table, :column => column.key),
                ),
            )
        end
    end

    extra_columns = unknown_columns(schema, normalized)
    if strict
        for key in extra_columns
            add_issue!(
                result,
                unknown_field(
                    "$(schema.table).$(key)",
                    "Unknown column $(key).";
                    context = Dict(:table => schema.table, :column => key),
                ),
            )
        end
    end

    for column in schema.columns
        haskey(normalized, column.key) && validate_column_value!(result, schema, column, normalized[column.key])
    end

    return result
end

function validate_table_issues(schema::TableSchema, rows; strict::Bool = false)
    result = validation_result(source = "table $(schema.table)")
    for (index, row) in enumerate(rows)
        row_result = validate_row_issues(schema, normalize_row(row); strict = strict)
        for issue in row_result.issues
            context = copy(issue.context)
            context[:row] = index
            add_issue!(
                result,
                validation_issue(
                    issue.code;
                    severity = issue.severity,
                    kind = issue.kind,
                    subject = isempty(issue.subject) ? "row $(index)" : "row $(index) $(issue.subject)",
                    message = issue.message,
                    context = context,
                ),
            )
        end
    end
    return result
end

function validate_row(schema::TableSchema, row::AbstractDict; strict::Bool = false)
    assert_valid!(validate_row_issues(schema, row; strict = strict))
    return true
end

function validate_table(schema::TableSchema, rows; strict::Bool = false)
    assert_valid!(validate_table_issues(schema, rows; strict = strict))
    return true
end

function table_columns(rows::Vector{Dict{Symbol,Any}})
    cols = Symbol[]
    seen = Set{Symbol}()
    for row in rows
        for key in keys(row)
            if !(key in seen)
                push!(cols, key)
                push!(seen, key)
            end
        end
    end
    return cols
end

function table_matrix(rows::Vector{Dict{Symbol,Any}}, columns::Vector{Symbol} = table_columns(rows))
    return [[get(row, col, missing) for col in columns] for row in rows]
end

function update_rows!(rows::Vector{Dict{Symbol,Any}}, predicate::Function, changes::AbstractDict)
    normalized_changes = Dict{Symbol,Any}(Symbol(k) => v for (k, v) in changes)
    count = 0
    for row in rows
        if predicate(row)
            merge!(row, normalized_changes)
            count += 1
        end
    end
    return count
end

function format_table(rows::Vector{Dict{Symbol,Any}}, columns::Vector{Symbol} = table_columns(rows))
    isempty(columns) && return ""

    body = [[string(get(row, col, "")) for col in columns] for row in rows]
    header = [String(col) for col in columns]
    widths = [length(header[i]) for i in eachindex(header)]

    for row in body
        for i in eachindex(row)
            widths[i] = max(widths[i], length(row[i]))
        end
    end

    lines = String[]
    push!(lines, join([rpad(header[i], widths[i]) for i in eachindex(header)], "  "))
    push!(lines, join([repeat("-", widths[i]) for i in eachindex(header)], "  "))
    for row in body
        push!(lines, join([rpad(row[i], widths[i]) for i in eachindex(row)], "  "))
    end
    return join(lines, "\n")
end

end

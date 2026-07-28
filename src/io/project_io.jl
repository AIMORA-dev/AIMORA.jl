module ProjectIO

using Dates
using TOML

using ..AssetTables
using ..ProjectData

export write_table_csv,
       read_table_csv,
       write_project_file,
       read_project_file

const PROJECT_FORMAT = "aimora-project-v1"

function csv_escape(value)
    ismissing(value) && return ""
    text = string(value)
    if occursin(',', text) || occursin('"', text) || occursin('\n', text) || occursin('\r', text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end
    return text
end

function parse_csv_line(line::AbstractString)
    fields = String[]
    buffer = IOBuffer()
    in_quotes = false
    i = firstindex(line)
    while i <= lastindex(line)
        ch = line[i]
        if in_quotes
            if ch == '"'
                next_i = nextind(line, i)
                if next_i <= lastindex(line) && line[next_i] == '"'
                    write(buffer, '"')
                    i = next_i
                else
                    in_quotes = false
                end
            else
                write(buffer, ch)
            end
        elseif ch == '"'
            in_quotes = true
        elseif ch == ','
            push!(fields, String(take!(buffer)))
        else
            write(buffer, ch)
        end
        i = nextind(line, i)
    end
    push!(fields, String(take!(buffer)))
    in_quotes && error("Unclosed quote in CSV line: $(line)")
    return fields
end

function write_table_csv(path::AbstractString, rows::Vector{Dict{Symbol,Any}}, columns::Vector{Symbol} = AssetTables.table_columns(rows))
    dir = dirname(path)
    isempty(dir) || isdir(dir) || mkpath(dir)
    open(path, "w") do io
        println(io, join(string.(columns), ","))
        for row in rows
            println(io, join((csv_escape(get(row, column, missing)) for column in columns), ","))
        end
    end
    return path
end

function column_spec(schema::AssetTables.TableSchema, key::Symbol)
    for column in schema.columns
        column.key == key && return column
    end
    return nothing
end

function parse_bool(text::AbstractString)
    lower = lowercase(strip(text))
    lower == "true" && return true
    lower == "false" && return false
    error("Cannot parse Bool value: $(text)")
end

function parse_untyped_value(text::AbstractString)
    stripped = strip(text)
    isempty(stripped) && return missing
    lower = lowercase(stripped)
    lower == "true" && return true
    lower == "false" && return false
    if startswith(stripped, ":") && length(stripped) > 1
        return Symbol(stripped[2:end])
    end
    int_value = tryparse(Int, stripped)
    int_value !== nothing && return int_value
    float_value = tryparse(Float64, stripped)
    float_value !== nothing && return float_value
    return stripped
end

function parse_typed_value(text::AbstractString, value_type)
    stripped = strip(text)
    isempty(stripped) && return missing
    value_type === nothing && return parse_untyped_value(stripped)
    value_type === String && return String(stripped)
    value_type === Symbol && return Symbol(startswith(stripped, ":") ? stripped[2:end] : stripped)
    value_type === Bool && return parse_bool(stripped)
    value_type <: Integer && return parse(Int, stripped)
    value_type <: Real && return parse(Float64, stripped)
    return parse_untyped_value(stripped)
end

function read_table_csv(path::AbstractString; schema::Union{Nothing,AssetTables.TableSchema} = nothing, strict::Bool = false)
    lines = readlines(path)
    isempty(lines) && return Vector{Dict{Symbol,Any}}()

    columns = Symbol.(parse_csv_line(lines[1]))
    rows = Vector{Dict{Symbol,Any}}()
    for line in lines[2:end]
        isempty(strip(line)) && continue
        fields = parse_csv_line(line)
        length(fields) == length(columns) || error("CSV row has $(length(fields)) fields, expected $(length(columns)): $(line)")
        row = Dict{Symbol,Any}()
        for (key, field) in zip(columns, fields)
            spec = schema === nothing ? nothing : column_spec(schema, key)
            row[key] = parse_typed_value(field, spec === nothing ? nothing : spec.value_type)
        end
        if schema !== nothing
            row = AssetTables.apply_defaults(schema, row)
            AssetTables.validate_row(schema, row; strict = strict)
        end
        push!(rows, row)
    end
    return rows
end

toml_key(key::Symbol) = String(key)

function toml_value(value)
    if value isa Symbol
        return String(value)
    elseif value isa DateTime
        return Dates.format(value, dateformat"yyyy-mm-ddTHH:MM:SS.s")
    elseif value isa AbstractString || value isa Number || value isa Bool
        return value
    elseif value isa AbstractVector
        return [toml_value(item) for item in value]
    elseif value isa AbstractDict
        return Dict(string(key) => toml_value(val) for (key, val) in value)
    elseif ismissing(value) || value === nothing
        return ""
    else
        return string(value)
    end
end

function revision_to_dict(revision::ProjectData.Revision)
    data = Dict{String,Any}(
        "id" => String(revision.id),
        "created_at" => Dates.format(revision.created_at, dateformat"yyyy-mm-ddTHH:MM:SS.s"),
        "author" => revision.author,
        "description" => revision.description,
    )
    revision.parent_id !== nothing && (data["parent_id"] = String(revision.parent_id))
    return data
end

function scenario_to_dict(scenario::ProjectData.Scenario)
    return Dict{String,Any}(
        "id" => String(scenario.id),
        "name" => scenario.name,
        "asset_tables" => sort!(String[String(table) for table in keys(scenario.asset_tables)]),
        "study_settings" => [
            Dict{String,Any}(
                "study" => String(settings.study),
                "parameters" => toml_value(settings.parameters),
            )
            for settings in values(scenario.study_settings)
        ],
    )
end

function case_to_dict(case::ProjectData.Case)
    return Dict{String,Any}(
        "id" => String(case.id),
        "name" => case.name,
        "revisions" => [revision_to_dict(revision) for revision in case.revisions],
        "scenarios" => [scenario_to_dict(scenario) for scenario in values(case.scenarios)],
    )
end

function project_to_dict(project::ProjectData.Project)
    return Dict{String,Any}(
        "format" => PROJECT_FORMAT,
        "project" => Dict{String,Any}(
            "id" => String(project.id),
            "name" => project.name,
            "metadata" => toml_value(project.metadata),
        ),
        "cases" => [case_to_dict(case) for case in values(project.cases)],
    )
end

function write_project_file(path::AbstractString, project::ProjectData.Project)
    dir = dirname(path)
    isempty(dir) || isdir(dir) || mkpath(dir)
    open(path, "w") do io
        TOML.print(io, project_to_dict(project); sorted = true)
    end
    return path
end

function symbol_dict(data)
    return Dict{Symbol,Any}(Symbol(key) => value for (key, value) in data)
end

function parse_datetime(text::AbstractString)
    return DateTime(String(text), dateformat"yyyy-mm-ddTHH:MM:SS.s")
end

function revision_from_dict(data)
    return Revision(
        id = Symbol(data["id"]),
        parent_id = haskey(data, "parent_id") ? Symbol(data["parent_id"]) : nothing,
        created_at = parse_datetime(data["created_at"]),
        author = get(data, "author", ""),
        description = get(data, "description", ""),
    )
end

function scenario_from_dict(data)
    scenario = Scenario(id = Symbol(data["id"]), name = data["name"])
    for settings_data in get(data, "study_settings", Any[])
        parameters = haskey(settings_data, "parameters") ? symbol_dict(settings_data["parameters"]) : Dict{Symbol,Any}()
        set_study_settings!(
            scenario,
            StudySettings(study = Symbol(settings_data["study"]), parameters = parameters),
        )
    end
    return scenario
end

function case_from_dict(data)
    case = Case(id = Symbol(data["id"]), name = data["name"])
    for revision_data in get(data, "revisions", Any[])
        add_revision!(case, revision_from_dict(revision_data))
    end
    for scenario_data in get(data, "scenarios", Any[])
        add_scenario!(case, scenario_from_dict(scenario_data))
    end
    return case
end

function read_project_file(path::AbstractString)
    data = TOML.parsefile(path)
    get(data, "format", "") == PROJECT_FORMAT || error("Unsupported AIMORA project format: $(get(data, "format", ""))")
    project_data = data["project"]
    project = Project(
        id = Symbol(project_data["id"]),
        name = project_data["name"],
        metadata = haskey(project_data, "metadata") ? symbol_dict(project_data["metadata"]) : Dict{Symbol,Any}(),
    )
    for case_data in get(data, "cases", Any[])
        add_case!(project, case_from_dict(case_data))
    end
    return project
end

end

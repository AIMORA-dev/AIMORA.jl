module ProjectData

using Dates

using ..AssetTables

export StudySettings,
       Revision,
       Scenario,
       Case,
       Project,
       ColumnSpec,
       TableSchema,
       add_case!,
       add_scenario!,
       add_revision!,
       set_study_settings!,
       study_settings,
       set_asset_table_schema!,
       asset_table_schema,
       set_asset_table!,
       asset_table,
       add_asset!,
       update_assets!,
       table_columns,
       table_matrix,
       format_table

Base.@kwdef mutable struct StudySettings
    study::Symbol
    parameters::Dict{Symbol,Any} = Dict{Symbol,Any}()
end

Base.@kwdef struct Revision
    id::Symbol
    parent_id::Union{Nothing,Symbol} = nothing
    created_at::DateTime = now(UTC)
    author::String = ""
    description::String = ""
end

Base.@kwdef mutable struct Scenario
    id::Symbol
    name::String
    asset_tables::Dict{Symbol,Vector{Dict{Symbol,Any}}} = Dict{Symbol,Vector{Dict{Symbol,Any}}}()
    table_schemas::Dict{Symbol,TableSchema} = Dict{Symbol,TableSchema}()
    study_settings::Dict{Symbol,StudySettings} = Dict{Symbol,StudySettings}()
end

Base.@kwdef mutable struct Case
    id::Symbol
    name::String
    scenarios::Dict{Symbol,Scenario} = Dict{Symbol,Scenario}()
    revisions::Vector{Revision} = Revision[]
end

Base.@kwdef mutable struct Project
    id::Symbol
    name::String
    cases::Dict{Symbol,Case} = Dict{Symbol,Case}()
    metadata::Dict{Symbol,Any} = Dict{Symbol,Any}()
end

function add_case!(project::Project, case::Case)
    project.cases[case.id] = case
    return case
end

function add_scenario!(case::Case, scenario::Scenario)
    case.scenarios[scenario.id] = scenario
    return scenario
end

function add_revision!(case::Case, revision::Revision)
    push!(case.revisions, revision)
    return revision
end

function set_study_settings!(scenario::Scenario, settings::StudySettings)
    scenario.study_settings[settings.study] = settings
    return settings
end

function study_settings(scenario::Scenario, study::Symbol)
    return get(scenario.study_settings, study, nothing)
end

function set_asset_table_schema!(scenario::Scenario, schema::TableSchema)
    scenario.table_schemas[schema.table] = schema
    return schema
end

function asset_table_schema(scenario::Scenario, table::Symbol)
    return get(scenario.table_schemas, table, nothing)
end

function set_asset_table!(scenario::Scenario, table::Symbol, rows)
    schema = asset_table_schema(scenario, table)
    scenario.asset_tables[table] = if schema === nothing
        [AssetTables.normalize_row(row) for row in rows]
    else
        [AssetTables.apply_defaults(schema, row) for row in rows]
    end
    schema !== nothing && AssetTables.validate_table(schema, scenario.asset_tables[table]; strict = true)
    return scenario.asset_tables[table]
end

function asset_table(scenario::Scenario, table::Symbol)
    return get!(scenario.asset_tables, table, Vector{Dict{Symbol,Any}}())
end

function add_asset!(scenario::Scenario, table::Symbol, row)
    rows = asset_table(scenario, table)
    schema = asset_table_schema(scenario, table)
    normalized = schema === nothing ? AssetTables.normalize_row(row) : AssetTables.apply_defaults(schema, row)
    schema !== nothing && AssetTables.validate_row(schema, normalized; strict = true)
    push!(rows, normalized)
    return rows[end]
end

function update_assets!(scenario::Scenario, table::Symbol, predicate::Function, changes::AbstractDict)
    rows = asset_table(scenario, table)
    count = AssetTables.update_rows!(rows, predicate, changes)
    schema = asset_table_schema(scenario, table)
    schema !== nothing && AssetTables.validate_table(schema, rows; strict = true)
    return count
end

table_columns(rows::Vector{Dict{Symbol,Any}}) = AssetTables.table_columns(rows)
table_matrix(rows::Vector{Dict{Symbol,Any}}, columns::Vector{Symbol} = table_columns(rows)) =
    AssetTables.table_matrix(rows, columns)
format_table(rows::Vector{Dict{Symbol,Any}}, columns::Vector{Symbol} = table_columns(rows)) =
    AssetTables.format_table(rows, columns)

end

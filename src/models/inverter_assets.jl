module InverterAssets

using ..AssetTables

export inverter_row,
       inverter_schema,
       set_inverter_table!,
       inverter_table,
       update_all_inverters!,
       update_inverters!,
       inverter_table_text

const INVERTER_SCHEMA = TableSchema(
    :inverters,
    [
        ColumnSpec(:id; label = "ID", required = true, value_type = Symbol, description = "Unique inverter identifier."),
        ColumnSpec(:bus; label = "Bus", required = true, value_type = Symbol, description = "Connected bus."),
        ColumnSpec(:rated_kva; label = "Rated kVA", unit = "kVA", required = true, value_type = Real, min_value = 0.0),
        ColumnSpec(:v_ll_rms_v; label = "Line-line RMS voltage", unit = "V", required = true, value_type = Real, min_value = 0.0),
        ColumnSpec(:p_ref_pu; label = "P reference", unit = "pu", default = 1.0, value_type = Real),
        ColumnSpec(:q_ref_pu; label = "Q reference", unit = "pu", default = 0.0, value_type = Real),
        ColumnSpec(:enabled; label = "Enabled", default = true, value_type = Bool),
        ColumnSpec(:model; label = "Model", default = :grid_following, value_type = Symbol, allowed_values = [:grid_following, :grid_forming]),
    ];
    notes = "Scenario-level inverter asset table for bulk editing and study setup.",
)

inverter_schema() = INVERTER_SCHEMA

function inverter_row(;
    id::Symbol,
    bus::Symbol,
    rated_kva::Float64,
    v_ll_rms_v::Float64,
    p_ref_pu::Float64 = 1.0,
    q_ref_pu::Float64 = 0.0,
    enabled::Bool = true,
    model::Symbol = :grid_following,
)
    return AssetTables.apply_defaults(INVERTER_SCHEMA, Dict{Symbol,Any}(
        :id => id,
        :bus => bus,
        :rated_kva => rated_kva,
        :v_ll_rms_v => v_ll_rms_v,
        :p_ref_pu => p_ref_pu,
        :q_ref_pu => q_ref_pu,
        :enabled => enabled,
        :model => model,
    ))
end

function set_inverter_table!(scenario, rows)
    scenario.asset_tables[:inverters] = [AssetTables.apply_defaults(INVERTER_SCHEMA, row) for row in rows]
    AssetTables.validate_table(INVERTER_SCHEMA, scenario.asset_tables[:inverters]; strict = true)
    return scenario.asset_tables[:inverters]
end

inverter_table(scenario) =
    get!(scenario.asset_tables, :inverters, Vector{Dict{Symbol,Any}}())

function update_all_inverters!(scenario; kwargs...)
    changes = Dict{Symbol,Any}(Symbol(k) => v for (k, v) in kwargs)
    return update_inverters!(scenario, row -> true; changes...)
end

function update_inverters!(scenario, predicate::Function; kwargs...)
    changes = Dict{Symbol,Any}(Symbol(k) => v for (k, v) in kwargs)
    count = AssetTables.update_rows!(inverter_table(scenario), predicate, changes)
    AssetTables.validate_table(INVERTER_SCHEMA, inverter_table(scenario); strict = true)
    return count
end

function inverter_table_text(scenario)
    columns = [:id, :bus, :rated_kva, :v_ll_rms_v, :p_ref_pu, :q_ref_pu, :enabled, :model]
    return AssetTables.format_table(inverter_table(scenario), columns)
end

end

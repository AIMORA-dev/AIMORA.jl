using Test

const EquipmentCatalogCore = AIMORA.EquipmentCatalog

function catalog_equipment_fixture(category::Symbol, index::Int)
    return Dict{String,Any}(
        "id" => "aimora://catalog/system/$(category)@1.0.0",
        "scope" => "system",
        "category" => String(category),
        "label" => titlecase(replace(String(category), '_' => ' ')),
        "description" => "Synthetic $(category) catalogue contract fixture.",
        "equipment_class" => String(category),
        "symbol_id" => "fixture.$(category)",
        "designator_prefix" => "X$(index)",
        "keywords" => Any[String(category), "fixture"],
        "terminal" => Any[
            Dict{String,Any}(
                "id" => "terminal",
                "role" => category == :annotation ? "annotation" : "electrical",
                "direction" => category == :annotation ? "none" : "bidirectional",
            ),
        ],
        "part" => Any[
            Dict{String,Any}(
                "number" => "PART-$(index)",
                "description" => "Fixture part $(index)",
                "unit" => "item",
                "quantity" => 1,
            ),
        ],
    )
end

function equipment_catalog_fixture()
    equipment = Any[
        catalog_equipment_fixture(category, index)
        for (index, category) in enumerate(EquipmentCatalogCore.EQUIPMENT_CATEGORIES)
    ]
    return Dict{String,Any}(
        "schema" => EquipmentCatalogCore.EQUIPMENT_CATALOG_SCHEMA,
        "version" => string(EquipmentCatalogCore.EQUIPMENT_CATALOG_VERSION),
        "licence" => "PolyForm-Noncommercial-1.0.0",
        "collection" => Any[
            Dict{String,Any}("scope" => "system", "mutable" => false),
            Dict{String,Any}("scope" => "user", "mutable" => true),
            Dict{String,Any}("scope" => "project", "mutable" => true),
        ],
        "equipment" => equipment,
        "assembly" => Any[
            Dict{String,Any}(
                "id" => "aimora://catalog/system/assembly.fixture_bay@1.0.0",
                "scope" => "system",
                "category" => "switching",
                "label" => "Fixture bay",
                "description" => "Reusable test bay.",
                "keywords" => Any["bay", "fixture"],
                "member" => Any[
                    Dict{String,Any}(
                        "id" => "primary",
                        "equipment_id" => equipment[1]["id"],
                    ),
                    Dict{String,Any}(
                        "id" => "secondary",
                        "equipment_id" => equipment[4]["id"],
                        "parent" => "primary",
                        "cross_reference" => Any[
                            Dict{String,Any}("field" => "upstream", "target" => "primary"),
                        ],
                    ),
                ],
            ),
        ],
    )
end

@testset "Julia-owned equipment catalog and assembly contracts" begin
    document = equipment_catalog_fixture()
    catalog = EquipmentCatalogCore.parse_equipment_catalog(document)
    @test catalog.version == v"1.0.0"
    @test getfield.(catalog.collections, :scope) == [:system, :user, :project]
    @test Set(getfield.(catalog.equipment, :category)) ==
          Set(EquipmentCatalogCore.EQUIPMENT_CATEGORIES)
    @test length(catalog.assemblies) == 1

    results = EquipmentCatalogCore.search_equipment_catalog(catalog, "fixture bay")
    @test length(results) == 1
    @test only(results).kind == :assembly
    @test only(results).member_count == 2

    transformer = EquipmentCatalogCore.search_equipment_catalog(
        catalog,
        "transformer";
        scope = :system,
        category = :transformer,
        include_assemblies = false,
    )
    @test length(transformer) == 1
    @test only(transformer).equipment_class == "transformer"

    native = EquipmentCatalogCore.native_equipment_catalog_document(catalog)
    @test native["schema"] == "aimora-equipment-library-v1"
    @test length(native["entries"]) == length(catalog.equipment) + length(catalog.assemblies)
    @test Dict(row["scope"] => row["count"] for row in native["collections"]) ==
          Dict("system" => 14, "user" => 0, "project" => 0)

    duplicate = deepcopy(document)
    push!(duplicate["equipment"], deepcopy(first(duplicate["equipment"])))
    @test_throws ArgumentError EquipmentCatalogCore.parse_equipment_catalog(duplicate)

    missing_member = deepcopy(document)
    missing_member["assembly"][1]["member"][2]["equipment_id"] =
        "aimora://catalog/system/missing@1.0.0"
    @test_throws ArgumentError EquipmentCatalogCore.parse_equipment_catalog(missing_member)

    cyclic = deepcopy(document)
    cyclic["assembly"][1]["member"][1]["parent"] = "secondary"
    @test_throws ArgumentError EquipmentCatalogCore.parse_equipment_catalog(cyclic)
end

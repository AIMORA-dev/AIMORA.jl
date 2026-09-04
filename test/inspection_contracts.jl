using Test
using AIMORA

const Inspector = AIMORA.InspectionCore

function inspection_fixture(asset_id::String = "asset.breaker.1")
    units = [
        Inspector.InspectionUnit(symbol = "V", scale_to_canonical = 1.0),
        Inspector.InspectionUnit(symbol = "kV", scale_to_canonical = 1000.0),
    ]
    sections = Inspector.InspectionSection[]
    for (index, section_id) in enumerate(Inspector.INSPECTION_SECTION_IDS)
        if section_id == "study_facets"
            push!(
                sections,
                Inspector.InspectionSection(
                    id = section_id,
                    title = "Study facets",
                    available = false,
                    unavailable_reason = "No qualified harmonic model is bound to this asset.",
                ),
            )
            continue
        elseif section_id == "ratings"
            field = Inspector.InspectionField(
                path = "ratings.field",
                label = "Rated voltage",
                kind = :number,
                canonical_unit = "V",
                display_units = units,
                range = Inspector.InspectionRange(lower = 0.0),
                affects_model_paths = [asset_id],
                affects_view_ids = ["view.sld.1"],
                invalidates_result_ids = ["result.short-circuit.1"],
                provenance = "fixture.original",
            )
        elseif section_id == "curves_tables"
            field = Inspector.InspectionField(
                path = "curves.trip",
                label = "Trip curve",
                kind = :curve,
                provenance = "catalog.breaker",
            )
        elseif section_id == "results_freshness"
            field = Inspector.InspectionField(
                path = "results.current",
                label = "Current",
                kind = :number,
                canonical_unit = "V",
                display_units = units,
                read_only = true,
                provenance = "result.short-circuit.1",
            )
        else
            field = Inspector.InspectionField(
                path = "$(section_id).field",
                label = replace(titlecase(replace(section_id, "_" => " ")), " " => " "),
                kind = index == 1 ? :text : :number,
                range = index == 1 ? nothing : Inspector.InspectionRange(lower = 0.0),
                required = index == 1,
                affects_model_paths = [asset_id],
                affects_view_ids = ["view.sld.1"],
                invalidates_result_ids = ["result.short-circuit.1"],
                provenance = "fixture.original",
            )
        end
        push!(
            sections,
            Inspector.InspectionSection(
                id = section_id,
                title = titlecase(replace(section_id, "_" => " ")),
                fields = [field],
            ),
        )
    end
    schema = Inspector.InspectionSchema(
        equipment_class = "switching.breaker",
        sections = sections,
    )
    identity = Inspector.InspectionIdentity(
        project_id = "project.demo",
        asset_id = asset_id,
        projection_id = "projection.$(asset_id)",
        view_id = "view.sld.1",
        equipment_class = "switching.breaker",
        result_bindings = ["result.short-circuit.1"],
    )
    values = Dict{String,Any}(
        "general.field" => "Q1",
        "connections_terminals.field" => 2.0,
        "ratings.field" => 12_000.0,
        "equipment_construction.field" => 1.0,
        "curves.trip" => [Dict("row_id" => "point.1", "x" => 0.0, "y" => 1.0)],
        "controls_protection.field" => 1.0,
        "scenarios_events.field" => 1.0,
        "results.current" => Inspector.InspectionValue(value = 11_500.0, canonical_unit = "V"),
        "drawing_view.field" => 1.0,
        "validation_readiness.field" => 1.0,
        "provenance_audit.field" => 1.0,
    )
    return Inspector.InspectionSession(identity, schema, values; revision = 7)
end

@testset "Julia-owned schema-driven inspection" begin
    session = inspection_fixture()
    @test Inspector.validate_inspection_schema(session.schema) === session.schema
    @test length(session.schema.sections) == 12
    @test Set(Inspector.INSPECTION_SECTION_IDS) ⊆ Set(getfield.(session.schema.sections, :id))

    document = Inspector.inspection_document(session)
    @test document["schema_version"] == "1.0.0"
    @test document["revision"] == "7"
    @test document["identity"]["asset_id"] == "asset.breaker.1"
    @test document["identity"]["projection_id"] == "projection.asset.breaker.1"
    @test document["identity"]["view_id"] == "view.sld.1"
    @test document["identity"]["result_bindings"] == ["result.short-circuit.1"]
    @test !document["undo_available"]
    @test !document["redo_available"]
    @test any(section -> !section["available"], document["sections"])

    accepted = Inspector.commit_inspection!(
        session,
        Inspector.InspectionTransaction(
            base_revision = 7,
            edits = [
                Inspector.CanonicalInspectionEdit(
                    path = "ratings.field",
                    value = 13.8,
                    display_unit = "kV",
                ),
            ],
        ),
    )
    @test accepted.status == :accepted
    @test accepted.revision == 8
    @test session.values["ratings.field"].value == 13_800.0
    @test session.values["ratings.field"].canonical_unit == "V"
    @test accepted.affected_model_paths == ["asset.breaker.1"]
    @test accepted.affected_view_ids == ["view.sld.1"]
    @test accepted.invalidated_result_ids == ["result.short-circuit.1"]

    stale = Inspector.commit_inspection!(
        session,
        Inspector.InspectionTransaction(
            base_revision = 7,
            edits = [Inspector.CanonicalInspectionEdit(path = "ratings.field", value = 1.0)],
        ),
    )
    @test stale.status == :conflict
    @test stale.revision == 8
    @test only(stale.issues).code == "REVISION_CONFLICT"

    invalid_unit = Inspector.commit_inspection!(
        session,
        Inspector.InspectionTransaction(
            base_revision = 8,
            edits = [
                Inspector.CanonicalInspectionEdit(
                    path = "ratings.field",
                    value = 1.0,
                    display_unit = "A",
                ),
            ],
        ),
    )
    @test invalid_unit.status == :rejected
    @test only(invalid_unit.issues).code == "INVALID_UNIT"
    @test session.revision == 8

    read_only = Inspector.commit_inspection!(
        session,
        Inspector.InspectionTransaction(
            base_revision = 8,
            edits = [
                Inspector.CanonicalInspectionEdit(path = "results.current", value = 3.0),
            ],
        ),
    )
    @test read_only.status == :rejected
    @test only(read_only.issues).code == "READ_ONLY"

    validator_rejection = Inspector.commit_inspection!(
        session,
        Inspector.InspectionTransaction(
            base_revision = 8,
            edits = [Inspector.CanonicalInspectionEdit(path = "ratings.field", value = 15.0)],
        );
        validator = (_, _) -> [
            Inspector.InspectionIssue(
                path = "ratings.field",
                code = "OWNER_REJECTED",
                message = "The canonical asset owner rejected this rating.",
                owner_id = "asset.breaker.1",
            ),
        ],
    )
    @test validator_rejection.status == :rejected
    @test session.revision == 8

    undo = Inspector.undo_inspection!(session, 8)
    @test undo.status == :accepted
    @test undo.revision == 9
    @test session.values["ratings.field"].value == 12_000.0
    redo = Inspector.redo_inspection!(session, 9)
    @test redo.status == :accepted
    @test redo.revision == 10
    @test session.values["ratings.field"].value == 13_800.0

    second = inspection_fixture("asset.breaker.2")
    second.values["ratings.field"] = Inspector.InspectionValue(
        value = 11_000.0,
        canonical_unit = "V",
    )
    merged = Inspector.merge_inspection_documents([session, second])
    @test merged["selection_count"] == 2
    @test merged["asset_ids"] == ["asset.breaker.1", "asset.breaker.2"]
    @test merged["identity"]["asset_ids"] == merged["asset_ids"]
    @test merged["values"]["ratings.field"]["mixed"]
    @test merged["values"]["general.field"]["value"] == "Q1"
    @test length(merged["sections"]) == 12

    selection = [inspection_fixture("asset.breaker.1"), inspection_fixture("asset.breaker.2")]
    selection_result = Inspector.commit_inspection_selection!(
        selection,
        Inspector.InspectionTransaction(
            base_revision = 7,
            edits = [Inspector.CanonicalInspectionEdit(path = "general.field", value = "Q2")],
        ),
    )
    @test selection_result.status == :accepted
    @test all(item -> item.revision == 8, selection)
    @test all(item -> item.values["general.field"].value == "Q2", selection)
    @test Inspector.undo_inspection_selection!(selection, 8).status == :accepted
    @test all(item -> item.values["general.field"].value == "Q1", selection)
    @test Inspector.redo_inspection_selection!(selection, 9).status == :accepted
    @test all(item -> item.values["general.field"].value == "Q2", selection)

    diverged = [inspection_fixture("asset.breaker.1"), inspection_fixture("asset.breaker.2")]
    diverged[2].revision = 8
    rejected_selection = Inspector.commit_inspection_selection!(
        diverged,
        Inspector.InspectionTransaction(
            base_revision = 7,
            edits = [Inspector.CanonicalInspectionEdit(path = "general.field", value = "Q3")],
        ),
    )
    @test rejected_selection.status == :conflict
    @test all(item -> item.values["general.field"].value == "Q1", diverged)
    @test Inspector.undo_inspection!(session, -1).status == :conflict

    result_document = Inspector.inspection_result_document(accepted)
    @test result_document["status"] == "accepted"
    @test result_document["base_revision"] == "7"
    @test result_document["revision"] == "8"
end

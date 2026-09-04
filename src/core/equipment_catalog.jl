module EquipmentCatalog

export EQUIPMENT_CATALOG_SCHEMA,
       EQUIPMENT_CATALOG_VERSION,
       EQUIPMENT_COLLECTION_SCOPES,
       EQUIPMENT_CATEGORIES,
       CatalogCollection,
       CatalogTerminalDefinition,
       CatalogPartDefinition,
       EquipmentDefinition,
       CatalogCrossReferenceDefinition,
       CatalogAssemblyMemberDefinition,
       ReusableAssemblyDefinition,
       EquipmentCatalogDocument,
       CatalogSearchResult,
       parse_equipment_catalog,
       catalog_entries,
       search_equipment_catalog,
       native_equipment_catalog_document

const EQUIPMENT_CATALOG_SCHEMA = "aimora-equipment-library-v1"
const EQUIPMENT_CATALOG_VERSION = v"1.0.0"
const EQUIPMENT_COLLECTION_SCOPES = (:system, :user, :project)
const EQUIPMENT_CATEGORIES = (
    :transformer,
    :current_transformer,
    :voltage_transformer,
    :switching,
    :cable,
    :bus,
    :load,
    :grounding,
    :machine,
    :converter,
    :storage,
    :source,
    :annotation,
)
const _CATALOG_ID = r"^aimora://catalog/(system|user|project)/[a-z0-9._-]+@[1-9][0-9]*\.[0-9]+\.[0-9]+$"
const _SYMBOL_ID = r"^[a-z][a-z0-9_]*(?:\.[a-z][a-z0-9_]*)+$"
const _DESIGNATOR_PREFIX = r"^[A-Z][A-Z0-9]{0,7}$"
const _MAXIMUM_EQUIPMENT = 10_000
const _MAXIMUM_ASSEMBLIES = 2_000
const _MAXIMUM_MEMBERS = 1_000

struct CatalogCollection
    scope::Symbol
    mutable::Bool
end

struct CatalogTerminalDefinition
    id::String
    role::Symbol
    direction::Symbol
end

struct CatalogPartDefinition
    number::String
    description::String
    unit::String
    quantity::Int
end

struct EquipmentDefinition
    id::String
    scope::Symbol
    category::Symbol
    label::String
    description::String
    equipment_class::String
    symbol_id::String
    designator_prefix::String
    keywords::Vector{String}
    terminals::Vector{CatalogTerminalDefinition}
    parts::Vector{CatalogPartDefinition}
end

struct CatalogCrossReferenceDefinition
    field::String
    target::String
end

struct CatalogAssemblyMemberDefinition
    local_id::String
    equipment_id::String
    parent::Union{Nothing,String}
    cross_references::Vector{CatalogCrossReferenceDefinition}
end

struct ReusableAssemblyDefinition
    id::String
    scope::Symbol
    category::Symbol
    label::String
    description::String
    keywords::Vector{String}
    members::Vector{CatalogAssemblyMemberDefinition}
end

struct EquipmentCatalogDocument
    version::VersionNumber
    licence::String
    collections::Vector{CatalogCollection}
    equipment::Vector{EquipmentDefinition}
    assemblies::Vector{ReusableAssemblyDefinition}
end

struct CatalogSearchResult
    kind::Symbol
    id::String
    scope::Symbol
    category::Symbol
    label::String
    description::String
    equipment_class::String
    symbol_id::String
    designator_prefix::String
    keywords::Vector{String}
    parts::Vector{CatalogPartDefinition}
    member_count::Int
end

function _required(table::AbstractDict, key::AbstractString, context::AbstractString)
    haskey(table, key) || throw(ArgumentError("$context is missing '$key'"))
    return table[key]
end

function _required_string(table::AbstractDict, key::AbstractString, context::AbstractString)
    value = _required(table, key, context)
    value isa AbstractString || throw(ArgumentError("$context '$key' must be text"))
    normalized = String(value)
    isempty(strip(normalized)) && throw(ArgumentError("$context '$key' must not be empty"))
    ncodeunits(normalized) <= 4096 || throw(ArgumentError("$context '$key' is too large"))
    return normalized
end

function _string_array(value, context::AbstractString)
    value isa AbstractVector || throw(ArgumentError("$context must be an array"))
    length(value) <= 256 || throw(ArgumentError("$context is too large"))
    result = String[]
    for item in value
        item isa AbstractString || throw(ArgumentError("$context must contain text"))
        normalized = strip(String(item))
        isempty(normalized) && throw(ArgumentError("$context contains an empty value"))
        push!(result, normalized)
    end
    length(result) == length(unique(result)) ||
        throw(ArgumentError("$context contains duplicate values"))
    return result
end

function _catalog_scope(value, context)
    value isa AbstractString || throw(ArgumentError("$context scope must be text"))
    scope = Symbol(value)
    scope in EQUIPMENT_COLLECTION_SCOPES ||
        throw(ArgumentError("$context has an unsupported collection scope"))
    return scope
end

function _catalog_category(value, context)
    value isa AbstractString || throw(ArgumentError("$context category must be text"))
    category = Symbol(value)
    category in EQUIPMENT_CATEGORIES ||
        throw(ArgumentError("$context has an unsupported equipment category"))
    return category
end

function _parse_collection(table::AbstractDict)
    scope = _catalog_scope(_required(table, "scope", "catalog collection"), "catalog collection")
    mutable = _required(table, "mutable", "catalog collection")
    mutable isa Bool || throw(ArgumentError("catalog collection mutable flag must be boolean"))
    mutable == (scope != :system) ||
        throw(ArgumentError("system collection must be read-only and overlays must be mutable"))
    return CatalogCollection(scope, mutable)
end

function _parse_terminal(table::AbstractDict, equipment_id::String)
    context = "equipment '$equipment_id' terminal"
    id = _required_string(table, "id", context)
    occursin(r"^[a-z][a-z0-9_]*$", id) ||
        throw(ArgumentError("$context has an invalid ID"))
    role = Symbol(_required_string(table, "role", context))
    role in (:electrical, :control, :mechanical, :thermal, :annotation) ||
        throw(ArgumentError("$context has an unsupported role"))
    direction = Symbol(_required_string(table, "direction", context))
    direction in (:input, :output, :bidirectional, :none) ||
        throw(ArgumentError("$context has an unsupported direction"))
    return CatalogTerminalDefinition(id, role, direction)
end

function _parse_part(table::AbstractDict, context::String)
    number = _required_string(table, "number", context)
    description = _required_string(table, "description", context)
    unit = _required_string(table, "unit", context)
    quantity = _required(table, "quantity", context)
    quantity isa Integer && !(quantity isa Bool) && quantity > 0 ||
        throw(ArgumentError("$context quantity must be a positive integer"))
    quantity <= 1_000_000 || throw(ArgumentError("$context quantity is too large"))
    return CatalogPartDefinition(number, description, unit, Int(quantity))
end

function _parse_equipment(table::AbstractDict)
    id = _required_string(table, "id", "equipment")
    occursin(_CATALOG_ID, id) || throw(ArgumentError("equipment has an invalid stable catalog ID"))
    context = "equipment '$id'"
    scope = _catalog_scope(_required(table, "scope", context), context)
    category = _catalog_category(_required(table, "category", context), context)
    symbol_id = _required_string(table, "symbol_id", context)
    occursin(_SYMBOL_ID, symbol_id) || throw(ArgumentError("$context has an invalid symbol ID"))
    prefix = _required_string(table, "designator_prefix", context)
    occursin(_DESIGNATOR_PREFIX, prefix) ||
        throw(ArgumentError("$context has an invalid designator prefix"))
    keywords = _string_array(get(table, "keywords", Any[]), "$context keywords")
    terminals = CatalogTerminalDefinition[
        _parse_terminal(item, id) for item in get(table, "terminal", Any[])
    ]
    terminal_ids = getfield.(terminals, :id)
    length(terminal_ids) == length(unique(terminal_ids)) ||
        throw(ArgumentError("$context repeats a terminal ID"))
    parts = CatalogPartDefinition[
        _parse_part(item, "$context part") for item in get(table, "part", Any[])
    ]
    isempty(parts) && throw(ArgumentError("$context must define parts-list data"))
    return EquipmentDefinition(
        id,
        scope,
        category,
        _required_string(table, "label", context),
        _required_string(table, "description", context),
        _required_string(table, "equipment_class", context),
        symbol_id,
        prefix,
        keywords,
        terminals,
        parts,
    )
end

function _parse_member(table::AbstractDict, assembly_id::String)
    context = "assembly '$assembly_id' member"
    local_id = _required_string(table, "id", context)
    occursin(r"^[a-z][a-z0-9_]*(?:\.[a-z][a-z0-9_]*)*$", local_id) ||
        throw(ArgumentError("$context has an invalid local ID"))
    parent_value = get(table, "parent", nothing)
    parent = if parent_value === nothing
        nothing
    elseif parent_value isa AbstractString && !isempty(strip(parent_value))
        String(parent_value)
    else
        throw(ArgumentError("$context parent must be a non-empty local ID"))
    end
    references = CatalogCrossReferenceDefinition[]
    for item in get(table, "cross_reference", Any[])
        push!(
            references,
            CatalogCrossReferenceDefinition(
                _required_string(item, "field", "$context cross-reference"),
                _required_string(item, "target", "$context cross-reference"),
            ),
        )
    end
    fields = getfield.(references, :field)
    length(fields) == length(unique(fields)) ||
        throw(ArgumentError("$context repeats a cross-reference field"))
    return CatalogAssemblyMemberDefinition(
        local_id,
        _required_string(table, "equipment_id", context),
        parent,
        references,
    )
end

function _parse_assembly(table::AbstractDict)
    id = _required_string(table, "id", "assembly")
    occursin(_CATALOG_ID, id) || throw(ArgumentError("assembly has an invalid stable catalog ID"))
    context = "assembly '$id'"
    members = CatalogAssemblyMemberDefinition[
        _parse_member(item, id) for item in get(table, "member", Any[])
    ]
    0 < length(members) <= _MAXIMUM_MEMBERS ||
        throw(ArgumentError("$context has an invalid member count"))
    local_ids = getfield.(members, :local_id)
    length(local_ids) == length(unique(local_ids)) ||
        throw(ArgumentError("$context repeats a member local ID"))
    return ReusableAssemblyDefinition(
        id,
        _catalog_scope(_required(table, "scope", context), context),
        _catalog_category(_required(table, "category", context), context),
        _required_string(table, "label", context),
        _required_string(table, "description", context),
        _string_array(get(table, "keywords", Any[]), "$context keywords"),
        members,
    )
end

function _validate_member_graph(assembly::ReusableAssemblyDefinition, equipment_ids::Set{String})
    members = Dict(member.local_id => member for member in assembly.members)
    for member in assembly.members
        member.equipment_id in equipment_ids ||
            throw(ArgumentError("assembly '$(assembly.id)' references unknown equipment"))
        member.parent === nothing || haskey(members, member.parent) ||
            throw(ArgumentError("assembly '$(assembly.id)' references an unknown parent"))
        all(reference -> haskey(members, reference.target), member.cross_references) ||
            throw(ArgumentError("assembly '$(assembly.id)' has an unknown cross-reference"))
        visited = Set{String}([member.local_id])
        parent = member.parent
        while parent !== nothing
            parent in visited &&
                throw(ArgumentError("assembly '$(assembly.id)' contains a parent cycle"))
            push!(visited, parent)
            parent = members[parent].parent
        end
    end
end

function parse_equipment_catalog(document::AbstractDict)
    _required_string(document, "schema", "equipment catalog") == EQUIPMENT_CATALOG_SCHEMA ||
        throw(ArgumentError("equipment catalog schema is unsupported"))
    version = VersionNumber(_required_string(document, "version", "equipment catalog"))
    version == EQUIPMENT_CATALOG_VERSION ||
        throw(ArgumentError("equipment catalog version is unsupported"))
    collections = CatalogCollection[
        _parse_collection(item) for item in get(document, "collection", Any[])
    ]
    scopes = getfield.(collections, :scope)
    Set(scopes) == Set(EQUIPMENT_COLLECTION_SCOPES) && length(scopes) == 3 ||
        throw(ArgumentError("equipment catalog must declare system, user, and project collections"))
    equipment = EquipmentDefinition[
        _parse_equipment(item) for item in get(document, "equipment", Any[])
    ]
    0 < length(equipment) <= _MAXIMUM_EQUIPMENT ||
        throw(ArgumentError("equipment catalog has an invalid equipment count"))
    equipment_ids = getfield.(equipment, :id)
    length(equipment_ids) == length(unique(equipment_ids)) ||
        throw(ArgumentError("equipment catalog repeats an equipment ID"))
    Set(getfield.(equipment, :category)) == Set(EQUIPMENT_CATEGORIES) ||
        throw(ArgumentError("equipment catalog does not cover every required category"))
    assemblies = ReusableAssemblyDefinition[
        _parse_assembly(item) for item in get(document, "assembly", Any[])
    ]
    length(assemblies) <= _MAXIMUM_ASSEMBLIES ||
        throw(ArgumentError("equipment catalog has too many assemblies"))
    assembly_ids = getfield.(assemblies, :id)
    length(assembly_ids) == length(unique(assembly_ids)) ||
        throw(ArgumentError("equipment catalog repeats an assembly ID"))
    isempty(intersect(Set(equipment_ids), Set(assembly_ids))) ||
        throw(ArgumentError("equipment and assembly IDs overlap"))
    known_scopes = Set(scopes)
    all(item -> item.scope in known_scopes, equipment) ||
        throw(ArgumentError("equipment uses an undeclared collection scope"))
    all(item -> item.scope in known_scopes, assemblies) ||
        throw(ArgumentError("assembly uses an undeclared collection scope"))
    known_equipment = Set(equipment_ids)
    foreach(assembly -> _validate_member_graph(assembly, known_equipment), assemblies)
    return EquipmentCatalogDocument(
        version,
        _required_string(document, "licence", "equipment catalog"),
        sort!(collections; by = item -> findfirst(==(item.scope), EQUIPMENT_COLLECTION_SCOPES)),
        sort!(equipment; by = item -> item.id),
        sort!(assemblies; by = item -> item.id),
    )
end

function _assembly_parts(
    assembly::ReusableAssemblyDefinition,
    equipment_by_id::Dict{String,EquipmentDefinition},
)
    quantities = Dict{Tuple{String,String,String},Int}()
    for member in assembly.members, part in equipment_by_id[member.equipment_id].parts
        key = (part.number, part.description, part.unit)
        quantities[key] = get(quantities, key, 0) + part.quantity
    end
    return CatalogPartDefinition[
        CatalogPartDefinition(number, description, unit, quantity)
        for ((number, description, unit), quantity) in sort!(collect(quantities); by = first)
    ]
end

function catalog_entries(catalog::EquipmentCatalogDocument)
    equipment_rows = CatalogSearchResult[
        CatalogSearchResult(
            :equipment,
            item.id,
            item.scope,
            item.category,
            item.label,
            item.description,
            item.equipment_class,
            item.symbol_id,
            item.designator_prefix,
            copy(item.keywords),
            copy(item.parts),
            0,
        ) for item in catalog.equipment
    ]
    equipment_by_id = Dict(item.id => item for item in catalog.equipment)
    assembly_rows = CatalogSearchResult[
        CatalogSearchResult(
            :assembly,
            item.id,
            item.scope,
            item.category,
            item.label,
            item.description,
            "",
            "",
            "",
            copy(item.keywords),
            _assembly_parts(item, equipment_by_id),
            length(item.members),
        ) for item in catalog.assemblies
    ]
    return sort!(vcat(equipment_rows, assembly_rows); by = item -> (item.label, item.id))
end

function search_equipment_catalog(
    catalog::EquipmentCatalogDocument,
    query::AbstractString = "";
    scope::Union{Nothing,Symbol} = nothing,
    category::Union{Nothing,Symbol} = nothing,
    include_assemblies::Bool = true,
)
    scope === nothing || scope in EQUIPMENT_COLLECTION_SCOPES ||
        throw(ArgumentError("unsupported equipment collection scope"))
    category === nothing || category in EQUIPMENT_CATEGORIES ||
        throw(ArgumentError("unsupported equipment category"))
    terms = split(lowercase(strip(String(query))))
    return filter(catalog_entries(catalog)) do item
        scope !== nothing && item.scope != scope && return false
        category !== nothing && item.category != category && return false
        !include_assemblies && item.kind == :assembly && return false
        haystack = lowercase(join(
            (
                item.id,
                item.label,
                item.description,
                item.equipment_class,
                item.symbol_id,
                join(item.keywords, ' '),
            ),
            ' ',
        ))
        return all(term -> occursin(term, haystack), terms)
    end
end

_part_document(part::CatalogPartDefinition) = Dict{String,Any}(
    "number" => part.number,
    "description" => part.description,
    "unit" => part.unit,
    "quantity" => part.quantity,
)

function native_equipment_catalog_document(catalog::EquipmentCatalogDocument)
    entries = Any[
        Dict{String,Any}(
            "kind" => String(item.kind),
            "id" => item.id,
            "scope" => String(item.scope),
            "category" => String(item.category),
            "label" => item.label,
            "description" => item.description,
            "equipment_class" => item.equipment_class,
            "symbol_id" => item.symbol_id,
            "designator_prefix" => item.designator_prefix,
            "keywords" => copy(item.keywords),
            "parts" => Any[_part_document(part) for part in item.parts],
            "member_count" => item.member_count,
        ) for item in catalog_entries(catalog)
    ]
    return Dict{String,Any}(
        "schema" => EQUIPMENT_CATALOG_SCHEMA,
        "version" => string(catalog.version),
        "source_owner" => "AIMORAResources/AIMORACatalogs",
        "collections" => Any[
            Dict{String,Any}(
                "scope" => String(collection.scope),
                "mutable" => collection.mutable,
                "count" => count(
                    item -> item["scope"] == String(collection.scope),
                    entries,
                ),
            ) for collection in catalog.collections
        ],
        "entries" => entries,
    )
end

end

module NativeExtensions

using SHA
using UUIDs

using ..NonlinearNetwork
using ..StudyCore: ParameterProvenance,
                   PhysicalModelParameter

import ..NonlinearNetwork: accept_nonlinear_device_state!,
                           nonlinear_current_jacobian!,
                           nonlinear_device_formulation,
                           nonlinear_device_provenance,
                           nonlinear_terminal_nodes

export AbstractNativeExtension,
       AbstractExtensionControlBlock,
       AbstractExtensionElectricalDevice,
       AbstractExtensionNonlinearDevice,
       AbstractExtensionSource,
       DirectedExtensionEvent,
       ExtensionComponentCheckpoint,
       ExtensionContract,
       ExtensionExecutionResult,
       ExtensionFailure,
       ExtensionIdentity,
       ExtensionOutputValue,
       ExtensionRegistration,
       ExtensionRegistry,
       ExtensionStateFamily,
       ExtensionStateInventory,
       NativeCubicCurrentBranch,
       NativeSeriesRLCompanion,
       SampledSaturatingLag,
       accept_extension_event!,
       accept_extension_state!,
       construct_extension,
       extension_checkpoint,
       extension_companion,
       extension_contract,
       extension_event_crossed,
       extension_event_value,
       extension_identity,
       extension_outputs,
       validate_extension_result,
       extension_source_value,
       extension_state_signature,
       extension_terminal_nodes,
       migrate_extension_state,
       register_extension!,
       register_extension_migration!,
       registered_extension_identities,
       release_extension_task_output!,
       resolve_extension,
       restore_extension_checkpoint!,
       sample_extension_task!

const EXTENSION_API_VERSION = v"1.0.0"
const AIMORA_PACKAGE_UUID = UUID("bb4e3f50-0ac0-4ed3-9c38-7355ab8ea495")
const _STATE_FAMILIES = (
    :continuous,
    :algebraic,
    :discrete,
    :delayed,
    :scheduler,
    :random,
    :history,
    :output,
    :checkpoint,
)
const _SERVICES = Set((
    :initialize,
    :nonlinear_current,
    :jacobian,
    :companion_stamp,
    :state_acceptance,
    :event,
    :sampled_task,
    :source,
    :output,
    :checkpoint,
    :reusable_definition,
))

struct ExtensionFailure <: Exception
    code::Symbol
    operation::Symbol
    identity::Union{Nothing,Any}
    message::String

    function ExtensionFailure(
        code::Symbol,
        operation::Symbol,
        identity,
        message::AbstractString,
    )
        occursin(r"^[a-z][a-z0-9_]*$", String(code)) || throw(ArgumentError(
            "extension failure code is not portable",
        ))
        isempty(strip(message)) && throw(ArgumentError(
            "extension failure message must not be empty",
        ))
        return new(code, operation, identity, String(message))
    end
end

Base.showerror(io::IO, failure::ExtensionFailure) =
    print(io, String(failure.code), " during ", String(failure.operation), ": ", failure.message)

_extension_fail(code, operation, identity, message) =
    throw(ExtensionFailure(code, operation, identity, message))

struct ExtensionIdentity
    package_uuid::UUID
    namespace::String
    semantic_type::String
    semantic_version::VersionNumber
    implementation_symbol::Symbol
    api_version::VersionNumber
    content_sha256::String

    function ExtensionIdentity(
        package_uuid::UUID,
        namespace::AbstractString,
        semantic_type::AbstractString,
        semantic_version::VersionNumber,
        implementation_symbol::Symbol,
        api_version::VersionNumber,
        content_sha256::AbstractString,
    )
        package_uuid == UUID(UInt128(0)) && throw(ArgumentError(
            "extension package UUID must not be nil",
        ))
        namespace_text = String(namespace)
        occursin(r"^[a-z][a-z0-9_-]*(?:\.[a-z][a-z0-9_-]*)*$", namespace_text) ||
            throw(ArgumentError("extension namespace is not portable"))
        type_text = String(semantic_type)
        occursin(r"^[a-z][a-z0-9_]*(?:\.[a-z][a-z0-9_]*)*$", type_text) ||
            throw(ArgumentError("extension semantic type is not portable"))
        semantic_version.major > 0 || throw(ArgumentError(
            "extension semantic version major must be positive",
        ))
        isempty(String(implementation_symbol)) && throw(ArgumentError(
            "extension implementation symbol must not be empty",
        ))
        api_version.major == EXTENSION_API_VERSION.major || throw(ArgumentError(
            "extension API major version is unsupported",
        ))
        digest = String(content_sha256)
        occursin(r"^[0-9a-f]{64}$", digest) || throw(ArgumentError(
            "extension content identity must be lowercase SHA-256 hexadecimal",
        ))
        return new(
            package_uuid,
            namespace_text,
            type_text,
            semantic_version,
            implementation_symbol,
            api_version,
            digest,
        )
    end
end

Base.:(==)(left::ExtensionIdentity, right::ExtensionIdentity) =
    left.package_uuid == right.package_uuid && left.namespace == right.namespace &&
    left.semantic_type == right.semantic_type &&
    left.semantic_version == right.semantic_version &&
    left.implementation_symbol == right.implementation_symbol &&
    left.api_version == right.api_version && left.content_sha256 == right.content_sha256
Base.hash(identity::ExtensionIdentity, seed::UInt) = hash((
    identity.package_uuid,
    identity.namespace,
    identity.semantic_type,
    identity.semantic_version,
    identity.implementation_symbol,
    identity.api_version,
    identity.content_sha256,
), seed)

struct ExtensionStateFamily
    family::Symbol
    names::Tuple{Vararg{Symbol}}
    not_applicable_reason::Union{Nothing,Symbol}

    function ExtensionStateFamily(
        family::Symbol,
        names = ();
        not_applicable_reason::Union{Nothing,Symbol} = nothing,
    )
        family in _STATE_FAMILIES || throw(ArgumentError(
            "extension state family is unknown",
        ))
        owned_names = Tuple(Symbol.(names))
        length(owned_names) == length(unique(owned_names)) || throw(ArgumentError(
            "extension state family repeats a state name",
        ))
        if isempty(owned_names)
            isnothing(not_applicable_reason) && throw(ArgumentError(
                "empty extension state family requires a not-applicable reason",
            ))
        elseif !isnothing(not_applicable_reason)
            throw(ArgumentError(
                "extension state family cannot own state and be not applicable",
            ))
        end
        return new(family, owned_names, not_applicable_reason)
    end
end

Base.:(==)(left::ExtensionStateFamily, right::ExtensionStateFamily) =
    left.family === right.family && left.names == right.names &&
    left.not_applicable_reason == right.not_applicable_reason

struct ExtensionStateInventory
    families::NTuple{9,ExtensionStateFamily}

    function ExtensionStateInventory(families)
        owned = sort!(collect(families); by = family -> findfirst(==(family.family), _STATE_FAMILIES))
        getfield.(owned, :family) == collect(_STATE_FAMILIES) || throw(ArgumentError(
            "extension state inventory must dispose every state family exactly once",
        ))
        return new(Tuple(owned))
    end
end

Base.:(==)(left::ExtensionStateInventory, right::ExtensionStateInventory) =
    left.families == right.families

struct ExtensionContract
    identity::ExtensionIdentity
    category::Symbol
    representation::Symbol
    fidelity::Symbol
    terminal_count::Int
    services::Tuple{Vararg{Symbol}}
    state::ExtensionStateInventory
    validity_domain::String
    unsupported::Tuple{Vararg{Symbol}}

    function ExtensionContract(
        identity::ExtensionIdentity,
        category::Symbol,
        representation::Symbol,
        fidelity::Symbol,
        terminal_count::Integer,
        services,
        state::ExtensionStateInventory,
        validity_domain::AbstractString;
        unsupported = (),
    )
        category in (:control, :electrical, :nonlinear_electrical, :source, :composed) ||
            throw(ArgumentError("extension category is unsupported"))
        representation === :instantaneous_emt || throw(ArgumentError(
            "native extension contract currently admits only instantaneous EMT",
        ))
        fidelity === :switching_detailed || throw(ArgumentError(
            "native extension contract currently admits only switching-detailed fidelity",
        ))
        terminals = Int(terminal_count)
        terminals >= 1 || throw(ArgumentError("extension terminal count must be positive"))
        service_tuple = Tuple(Symbol.(services))
        isempty(service_tuple) && throw(ArgumentError("extension contract requires services"))
        all(in(_SERVICES), service_tuple) || throw(ArgumentError(
            "extension contract contains an unknown service",
        ))
        length(service_tuple) == length(unique(service_tuple)) || throw(ArgumentError(
            "extension contract repeats a service",
        ))
        :initialize in service_tuple || throw(ArgumentError(
            "extension contract requires initialization",
        ))
        :checkpoint in service_tuple || throw(ArgumentError(
            "extension contract requires checkpoint state",
        ))
        :nonlinear_current in service_tuple && !(:jacobian in service_tuple) &&
            throw(ArgumentError("nonlinear extension requires an analytic Jacobian"))
        validity = String(validity_domain)
        isempty(strip(validity)) && throw(ArgumentError(
            "extension validity domain must not be empty",
        ))
        unsupported_tuple = Tuple(Symbol.(unsupported))
        return new(
            identity,
            category,
            representation,
            fidelity,
            terminals,
            service_tuple,
            state,
            validity,
            unsupported_tuple,
        )
    end
end

Base.:(==)(left::ExtensionContract, right::ExtensionContract) =
    left.identity == right.identity && left.category === right.category &&
    left.representation === right.representation && left.fidelity === right.fidelity &&
    left.terminal_count == right.terminal_count && left.services == right.services &&
    left.state == right.state && left.validity_domain == right.validity_domain &&
    left.unsupported == right.unsupported

abstract type AbstractNativeExtension end
abstract type AbstractExtensionControlBlock <: AbstractNativeExtension end
abstract type AbstractExtensionElectricalDevice <: AbstractNativeExtension end
abstract type AbstractExtensionNonlinearDevice <: NonlinearNetwork.AbstractNonlinearCurrentDevice end
abstract type AbstractExtensionSource <: AbstractNativeExtension end

extension_identity(::Type{T}) where {T<:Union{AbstractNativeExtension,AbstractExtensionNonlinearDevice}} =
    throw(MethodError(extension_identity, (T,)))
extension_contract(::Type{T}) where {T<:Union{AbstractNativeExtension,AbstractExtensionNonlinearDevice}} =
    throw(MethodError(extension_contract, (T,)))
extension_contract(component::Union{AbstractNativeExtension,AbstractExtensionNonlinearDevice}) =
    extension_contract(typeof(component))
extension_identity(declaration) = throw(MethodError(extension_identity, (declaration,)))
extension_source_value(component::AbstractNativeExtension, time_s::Real) =
    throw(MethodError(extension_source_value, (component, time_s)))

struct ExtensionComponentCheckpoint{S}
    identity::ExtensionIdentity
    state::S
    state_sha256::String
end

function _require_extension_method(
    type::DataType,
    operation::Function,
    signature::Type{<:Tuple},
    identity::ExtensionIdentity,
)
    hasmethod(operation, signature) || _extension_fail(
        :missing_extension_method,
        :register,
        identity,
        "$(nameof(type)) does not implement $(nameof(operation)) for its declared service",
    )
    return nothing
end

function _validate_extension_methods!(type::DataType, contract::ExtensionContract)
    identity = contract.identity
    _require_extension_method(type, extension_checkpoint, Tuple{type}, identity)
    _require_extension_method(
        type,
        restore_extension_checkpoint!,
        Tuple{type,ExtensionComponentCheckpoint},
        identity,
    )
    if contract.category === :control
        _require_extension_method(
            type,
            sample_extension_task!,
            Tuple{type,Float64,Int,Float64},
            identity,
        )
        :output in contract.services && _require_extension_method(
            type,
            extension_outputs,
            Tuple{type,Float64},
            identity,
        )
    elseif contract.category === :electrical
        _require_extension_method(type, extension_terminal_nodes, Tuple{type}, identity)
        _require_extension_method(
            type,
            extension_companion,
            Tuple{type,Float64},
            identity,
        )
        _require_extension_method(
            type,
            accept_extension_state!,
            Tuple{type,Vector{Float64},Float64},
            identity,
        )
        :output in contract.services && _require_extension_method(
            type,
            extension_outputs,
            Tuple{type,Float64},
            identity,
        )
    elseif contract.category === :nonlinear_electrical
        _require_extension_method(
            type,
            NonlinearNetwork.nonlinear_terminal_nodes,
            Tuple{type},
            identity,
        )
        _require_extension_method(
            type,
            NonlinearNetwork.nonlinear_current_jacobian!,
            Tuple{Vector{Float64},Matrix{Float64},type,Vector{Float64},Float64},
            identity,
        )
        :output in contract.services && _require_extension_method(
            type,
            extension_outputs,
            Tuple{type,Vector{Float64},Float64},
            identity,
        )
    end
    :source in contract.services && _require_extension_method(
        type,
        extension_source_value,
        Tuple{type,Float64},
        identity,
    )
    return nothing
end

struct ExtensionRegistration{T}
    identity::ExtensionIdentity
    implementation_type::Type{T}
    contract::ExtensionContract
end

struct ExtensionMigration{F}
    from::ExtensionIdentity
    to::ExtensionIdentity
    operation::F
end

mutable struct ExtensionRegistry
    registrations::Dict{ExtensionIdentity,Any}
    migrations::Dict{Tuple{ExtensionIdentity,ExtensionIdentity},Any}
end

ExtensionRegistry() = ExtensionRegistry(
    Dict{ExtensionIdentity,Any}(),
    Dict{Tuple{ExtensionIdentity,ExtensionIdentity},Any}(),
)

function _implementation_package_uuid(type::DataType)
    root = Base.moduleroot(parentmodule(type))
    return Base.PkgId(root).uuid
end

function register_extension!(registry::ExtensionRegistry, type::Type{T}) where {
    T<:Union{AbstractNativeExtension,AbstractExtensionNonlinearDevice},
}
    identity = extension_identity(type)
    contract = extension_contract(type)
    contract.identity == identity || _extension_fail(
        :extension_contract_identity_mismatch,
        :register,
        identity,
        "extension contract does not use its registered identity",
    )
    _implementation_package_uuid(type) == identity.package_uuid || _extension_fail(
        :extension_package_uuid_mismatch,
        :register,
        identity,
        "loaded implementation package UUID does not match the declared identity",
    )
    _validate_extension_methods!(type, contract)
    haskey(registry.registrations, identity) && begin
        existing = registry.registrations[identity]
        existing.implementation_type === type && return existing
        _extension_fail(
            :duplicate_extension_registration,
            :register,
            identity,
            "extension identity is already bound to another loaded type",
        )
    end
    registration = ExtensionRegistration(identity, type, contract)
    registry.registrations[identity] = registration
    return registration
end

function resolve_extension(registry::ExtensionRegistry, identity::ExtensionIdentity)
    haskey(registry.registrations, identity) || _extension_fail(
        :unregistered_extension,
        :resolve,
        identity,
        "load and explicitly register the exact extension implementation first",
    )
    return registry.registrations[identity]
end

function registered_extension_identities(registry::ExtensionRegistry)
    return sort!(collect(keys(registry.registrations)); by = identity -> (
        identity.namespace,
        identity.semantic_type,
        identity.semantic_version,
        string(identity.package_uuid),
        identity.content_sha256,
    ))
end

function construct_extension(registry::ExtensionRegistry, identity::ExtensionIdentity, args...; kwargs...)
    registration = resolve_extension(registry, identity)
    component = registration.implementation_type(args...; kwargs...)
    extension_contract(component).identity == identity || _extension_fail(
        :constructed_extension_identity_mismatch,
        :construct,
        identity,
        "constructed extension does not retain the registered identity",
    )
    return component
end

function register_extension_migration!(
    registry::ExtensionRegistry,
    from::ExtensionIdentity,
    to::ExtensionIdentity,
    operation::F,
) where {F}
    from.package_uuid == to.package_uuid || _extension_fail(
        :extension_migration_package_mismatch,
        :register_migration,
        from,
        "extension migration cannot change package UUID",
    )
    from.namespace == to.namespace && from.semantic_type == to.semantic_type ||
        _extension_fail(
            :extension_migration_type_mismatch,
            :register_migration,
            from,
            "extension migration cannot change semantic type",
        )
    from.semantic_version < to.semantic_version || _extension_fail(
        :extension_migration_not_forward,
        :register_migration,
        from,
        "extension migration target must be newer than its source",
    )
    key = (from, to)
    haskey(registry.migrations, key) && _extension_fail(
        :duplicate_extension_migration,
        :register_migration,
        from,
        "extension migration is already registered",
    )
    registry.migrations[key] = ExtensionMigration(from, to, operation)
    return registry.migrations[key]
end

function migrate_extension_state(
    registry::ExtensionRegistry,
    from::ExtensionIdentity,
    to::ExtensionIdentity,
    source,
)
    key = (from, to)
    haskey(registry.migrations, key) || _extension_fail(
        :missing_extension_migration,
        :migrate,
        from,
        "no exact one-step migration is registered",
    )
    migration = registry.migrations[key]
    before = deepcopy(source)
    first_result = migration.operation(deepcopy(source))
    isequal(source, before) || _extension_fail(
        :extension_migration_mutated_source,
        :migrate,
        from,
        "extension migration mutated its source value",
    )
    second_result = migration.operation(deepcopy(source))
    isequal(first_result, second_result) || _extension_fail(
        :nondeterministic_extension_migration,
        :migrate,
        from,
        "extension migration does not replay deterministically",
    )
    return first_result
end

struct ExtensionOutputValue
    name::Symbol
    value::Float64
    unit::String
    time_s::Float64
    validity::Symbol
    identity::ExtensionIdentity

    function ExtensionOutputValue(
        name::Symbol,
        value::Real,
        unit::AbstractString,
        time_s::Real,
        validity::Symbol,
        identity::ExtensionIdentity,
    )
        numeric_value = Float64(value)
        time = Float64(time_s)
        all(isfinite, (numeric_value, time)) || throw(ArgumentError(
            "extension output value and time must be finite",
        ))
        time >= 0.0 || throw(ArgumentError("extension output time must be nonnegative"))
        isempty(strip(unit)) && throw(ArgumentError("extension output unit must not be empty"))
        validity in (:valid, :warning, :invalid) || throw(ArgumentError(
            "extension output validity is unknown",
        ))
        return new(name, numeric_value, String(unit), time, validity, identity)
    end
end

"""Typed terminal record for one accepted or failed extension execution boundary."""
struct ExtensionExecutionResult{O<:Tuple,D<:NamedTuple}
    accepted::Bool
    identity::ExtensionIdentity
    representation::Symbol
    fidelity::Symbol
    terminal_nodes::Tuple{Vararg{Int}}
    state_inventory::ExtensionStateInventory
    project_sha256::String
    outputs::O
    event_cursor::Int
    task_cursor::Int
    output_cursor::Int
    checkpoint_sha256::String
    diagnostics::D
    units_and_bases::Tuple{Vararg{String}}
    warnings::Tuple{Vararg{Symbol}}
    failure::Union{Nothing,ExtensionFailure}
    deterministic_sha256::String

    function ExtensionExecutionResult(
        accepted::Bool,
        contract::ExtensionContract,
        terminal_nodes,
        project_sha256::AbstractString,
        outputs::O,
        event_cursor::Integer,
        task_cursor::Integer,
        checkpoint_sha256::AbstractString,
        diagnostics::D,
        units_and_bases,
        warnings = ();
        failure::Union{Nothing,ExtensionFailure} = nothing,
    ) where {O<:Tuple,D<:NamedTuple}
        accepted == isnothing(failure) || throw(ArgumentError(
            "accepted extension result and failure state disagree",
        ))
        terminals = Tuple(Int.(terminal_nodes))
        length(terminals) == contract.terminal_count || throw(ArgumentError(
            "extension result terminal count differs from its contract",
        ))
        all(>=(0), terminals) && length(unique(terminals)) == length(terminals) ||
            throw(ArgumentError("extension result terminals must be distinct nonnegative nodes"))
        project_digest = String(project_sha256)
        checkpoint_digest = String(checkpoint_sha256)
        all(digest -> occursin(r"^[0-9a-f]{64}$", digest), (project_digest, checkpoint_digest)) ||
            throw(ArgumentError("extension result signatures must be lowercase SHA-256 hexadecimal"))
        all(output -> output isa ExtensionOutputValue, outputs) || throw(ArgumentError(
            "extension result outputs must use ExtensionOutputValue",
        ))
        all(output -> output.identity == contract.identity, outputs) || throw(ArgumentError(
            "extension result output identity differs from its contract",
        ))
        event_index = Int(event_cursor)
        task_index = Int(task_cursor)
        event_index >= 0 && task_index >= 0 || throw(ArgumentError(
            "extension result event and task cursors must be nonnegative",
        ))
        unit_contract = Tuple(String.(units_and_bases))
        all(unit -> !isempty(strip(unit)), unit_contract) || throw(ArgumentError(
            "extension result unit/base declarations must not be empty",
        ))
        warning_codes = Tuple(Symbol.(warnings))
        signature_payload = (
            accepted = accepted,
            identity = contract.identity,
            representation = contract.representation,
            fidelity = contract.fidelity,
            terminal_nodes = terminals,
            state_inventory = contract.state,
            project_sha256 = project_digest,
            outputs = outputs,
            event_cursor = event_index,
            task_cursor = task_index,
            output_cursor = length(outputs),
            checkpoint_sha256 = checkpoint_digest,
            diagnostics = diagnostics,
            units_and_bases = unit_contract,
            warnings = warning_codes,
            failure_code = isnothing(failure) ? nothing : failure.code,
        )
        deterministic_digest = extension_state_signature(signature_payload)
        return new{O,D}(
            accepted,
            contract.identity,
            contract.representation,
            contract.fidelity,
            terminals,
            contract.state,
            project_digest,
            outputs,
            event_index,
            task_index,
            length(outputs),
            checkpoint_digest,
            diagnostics,
            unit_contract,
            warning_codes,
            failure,
            deterministic_digest,
        )
    end
end

"""Verify one typed result against the expected extension, Project, checkpoint, and output-unit identities."""
function validate_extension_result(
    result::ExtensionExecutionResult,
    contract::ExtensionContract,
    project_sha256::AbstractString,
    checkpoint_sha256::AbstractString,
    output_units,
)
    result.identity == contract.identity || _extension_fail(
        :stale_extension_identity,
        :validate_result,
        contract.identity,
        "extension result identity differs from the expected contract",
    )
    result.representation === contract.representation &&
        result.fidelity === contract.fidelity &&
        result.state_inventory == contract.state || _extension_fail(
            :stale_extension_contract,
            :validate_result,
            contract.identity,
            "extension result representation, fidelity, or state inventory is stale",
        )
    expected_project = String(project_sha256)
    expected_checkpoint = String(checkpoint_sha256)
    result.project_sha256 == expected_project || _extension_fail(
        :stale_extension_project_signature,
        :validate_result,
        contract.identity,
        "extension result Project signature is stale",
    )
    result.checkpoint_sha256 == expected_checkpoint || _extension_fail(
        :stale_extension_checkpoint_signature,
        :validate_result,
        contract.identity,
        "extension result checkpoint signature is stale",
    )
    expected_units = Tuple(String.(output_units))
    actual_units = getfield.(result.outputs, :unit)
    actual_units == expected_units || _extension_fail(
        :extension_output_unit_mismatch,
        :validate_result,
        contract.identity,
        "extension result output units differ from the declared result contract",
    )
    return result
end

function _state_text(value)
    if value isa NamedTuple
        return join((String(name) * "=" * _state_text(getfield(value, name)) for name in keys(value)), ";")
    elseif value isa Tuple || value isa AbstractVector
        return "[" * join((_state_text(item) for item in value), ",") * "]"
    elseif value isa Union{Nothing,Bool,Integer,AbstractFloat,Symbol,AbstractString,VersionNumber,UUID}
        return repr(value)
    end
    names = fieldnames(typeof(value))
    return string(nameof(typeof(value)), "{", join((String(name) * "=" * _state_text(getfield(value, name)) for name in names), ";"), "}")
end

extension_state_signature(value) = bytes2hex(sha256("aimora-extension-state-v1\n" * _state_text(value)))

function _component_checkpoint(component, state)
    return ExtensionComponentCheckpoint(
        extension_identity(typeof(component)),
        state,
        extension_state_signature(state),
    )
end

function _validate_checkpoint(component, checkpoint::ExtensionComponentCheckpoint)
    identity = extension_identity(typeof(component))
    checkpoint.identity == identity || _extension_fail(
        :incompatible_extension_checkpoint,
        :restore_checkpoint,
        identity,
        "checkpoint identity does not match the receiving extension",
    )
    extension_state_signature(checkpoint.state) == checkpoint.state_sha256 || _extension_fail(
        :corrupt_extension_checkpoint,
        :restore_checkpoint,
        identity,
        "checkpoint state digest does not match its payload",
    )
    return checkpoint.state
end

function _example_identity(type_name::String, symbol::Symbol, equation::String)
    return ExtensionIdentity(
        AIMORA_PACKAGE_UUID,
        "aimora.examples",
        type_name,
        v"1.0.0",
        symbol,
        EXTENSION_API_VERSION,
        bytes2hex(sha256("aimora-native-extension-example-v1\n" * equation)),
    )
end

function _state_inventory(owned::AbstractDict{Symbol,<:Tuple})
    return ExtensionStateInventory([
        haskey(owned, family) ?
        ExtensionStateFamily(family, owned[family]) :
        ExtensionStateFamily(family; not_applicable_reason = Symbol("not_applicable_", family))
        for family in _STATE_FAMILIES
    ])
end

const _SAMPLED_LAG_IDENTITY = _example_identity(
    "sampled_saturating_lag",
    :SampledSaturatingLag,
    "a=1-exp(-h/tau);x_next=x+a*(u-x);y=clamp(gain*x_next,y_min,y_max)",
)

mutable struct SampledSaturatingLag <: AbstractExtensionControlBlock
    gain::Float64
    time_constant_s::Float64
    minimum_output::Float64
    maximum_output::Float64
    period_ticks::Int
    delay_ticks::Int
    state::Float64
    held_input::Float64
    held_output::Float64
    pending_output::Float64
    next_sample_tick::Int
    pending_release_tick::Int
    sample_count::Int
    write_count::Int
    provenance::ParameterProvenance

    function SampledSaturatingLag(
        gain::Real,
        time_constant_s::Real,
        minimum_output::Real,
        maximum_output::Real,
        period_ticks::Integer;
        delay_ticks::Integer = 0,
        initial_state::Real = 0.0,
        provenance::ParameterProvenance = ParameterProvenance(
            "caller-supplied sampled lag parameters",
            "declared control unit and seconds",
            "converted to finite Float64 values and exact integer ticks",
            "not supplied; caller retains uncertainty ownership",
            "positive time constant/period and bounded output",
            PhysicalModelParameter,
        ),
    )
        values = Float64.((gain, time_constant_s, minimum_output, maximum_output, initial_state))
        all(isfinite, values) || throw(ArgumentError("sampled lag parameters must be finite"))
        gain_value, time_constant, minimum, maximum, state = values
        time_constant > 0.0 || throw(ArgumentError("sampled lag time constant must be positive"))
        minimum <= maximum || throw(ArgumentError("sampled lag output bounds are reversed"))
        period = Int(period_ticks)
        delay = Int(delay_ticks)
        period > 0 || throw(ArgumentError("sampled lag period must be positive"))
        0 <= delay < period || throw(ArgumentError(
            "sampled lag delay must be nonnegative and shorter than its period",
        ))
        provenance.nature === PhysicalModelParameter || throw(ArgumentError(
            "sampled lag provenance must describe physical model parameters",
        ))
        output = clamp(gain_value * state, minimum, maximum)
        return new(
            gain_value,
            time_constant,
            minimum,
            maximum,
            period,
            delay,
            state,
            state,
            output,
            output,
            0,
            -1,
            0,
            0,
            provenance,
        )
    end
end

extension_identity(::Type{SampledSaturatingLag}) = _SAMPLED_LAG_IDENTITY
extension_contract(::Type{SampledSaturatingLag}) = ExtensionContract(
    _SAMPLED_LAG_IDENTITY,
    :control,
    :instantaneous_emt,
    :switching_detailed,
    1,
    (:initialize, :sampled_task, :event, :source, :output, :checkpoint, :reusable_definition),
    _state_inventory(Dict(
        :continuous => (:state,),
        :discrete => (:saturation_mode,),
        :delayed => (:pending_output, :pending_release_tick),
        :scheduler => (:next_sample_tick, :sample_count, :write_count),
        :history => (:held_input, :held_output),
        :output => (:held_output,),
        :checkpoint => (:complete_component_state,),
    )),
    "finite scalar control input; positive time constant and exact period; delay shorter than period",
    unsupported = (:continuous_time_callback, :random_state, :implicit_algebraic_loop),
)

function sample_extension_task!(
    control::SampledSaturatingLag,
    input::Real,
    tick::Integer,
    tick_s::Real,
)
    input_value = Float64(input)
    tick_duration = Float64(tick_s)
    current_tick = Int(tick)
    isfinite(input_value) || throw(ArgumentError("sampled lag input must be finite"))
    isfinite(tick_duration) && tick_duration > 0.0 || throw(ArgumentError(
        "sampled lag scheduler tick duration must be finite and positive",
    ))
    current_tick == control.next_sample_tick || _extension_fail(
        :missed_extension_task,
        :sample_task,
        _SAMPLED_LAG_IDENTITY,
        "sampled lag must execute on its exact next tick",
    )
    step_s = control.period_ticks * tick_duration
    blend = -expm1(-step_s / control.time_constant_s)
    control.held_input = input_value
    control.state += blend * (input_value - control.state)
    control.pending_output = clamp(
        control.gain * control.state,
        control.minimum_output,
        control.maximum_output,
    )
    control.pending_release_tick = current_tick + control.delay_ticks
    control.next_sample_tick += control.period_ticks
    control.sample_count += 1
    control.delay_ticks == 0 && release_extension_task_output!(control, current_tick)
    return control.pending_output
end

function release_extension_task_output!(control::SampledSaturatingLag, tick::Integer)
    current_tick = Int(tick)
    current_tick == control.pending_release_tick || _extension_fail(
        :extension_task_release_mismatch,
        :release_task_output,
        _SAMPLED_LAG_IDENTITY,
        "sampled lag output must release on its exact pending tick",
    )
    control.held_output = control.pending_output
    control.pending_release_tick = -1
    control.write_count += 1
    return control.held_output
end

extension_checkpoint(control::SampledSaturatingLag) = _component_checkpoint(control, (
    state = control.state,
    held_input = control.held_input,
    held_output = control.held_output,
    pending_output = control.pending_output,
    next_sample_tick = control.next_sample_tick,
    pending_release_tick = control.pending_release_tick,
    sample_count = control.sample_count,
    write_count = control.write_count,
))

function restore_extension_checkpoint!(
    control::SampledSaturatingLag,
    checkpoint::ExtensionComponentCheckpoint,
)
    state = _validate_checkpoint(control, checkpoint)
    control.state = state.state
    control.held_input = state.held_input
    control.held_output = state.held_output
    control.pending_output = state.pending_output
    control.next_sample_tick = state.next_sample_tick
    control.pending_release_tick = state.pending_release_tick
    control.sample_count = state.sample_count
    control.write_count = state.write_count
    return control
end

extension_outputs(control::SampledSaturatingLag, time_s::Real) = (
    ExtensionOutputValue(:state, control.state, "declared_control_unit", time_s, :valid, _SAMPLED_LAG_IDENTITY),
    ExtensionOutputValue(:output, control.held_output, "declared_control_unit", time_s, :valid, _SAMPLED_LAG_IDENTITY),
)

function extension_source_value(control::SampledSaturatingLag, time_s::Real)
    time = Float64(time_s)
    isfinite(time) && time >= 0.0 || throw(ArgumentError(
        "sampled lag source evaluation time must be finite and nonnegative",
    ))
    isfinite(control.held_output) || throw(ArgumentError(
        "sampled lag held source output must be finite",
    ))
    return control.held_output
end

const _CUBIC_CURRENT_IDENTITY = _example_identity(
    "passive_cubic_current_branch",
    :NativeCubicCurrentBranch,
    "v=vp-vn;i=g*v+k*v^3;d=g+3*k*v^2",
)

struct NativeCubicCurrentBranch <: AbstractExtensionNonlinearDevice
    positive_node::Int
    negative_node::Int
    linear_conductance_s::Float64
    cubic_coefficient_a_per_v3::Float64
    provenance::ParameterProvenance

    function NativeCubicCurrentBranch(
        positive_node::Integer,
        negative_node::Integer,
        linear_conductance_s::Real,
        cubic_coefficient_a_per_v3::Real;
        provenance::ParameterProvenance = ParameterProvenance(
            "caller-supplied passive cubic branch parameters",
            "siemens and ampere per volt cubed",
            "converted to finite nonnegative Float64 SI values",
            "not supplied; caller retains uncertainty ownership",
            "two distinct nonnegative node indices and at least one positive coefficient",
            PhysicalModelParameter,
        ),
    )
        positive = Int(positive_node)
        negative = Int(negative_node)
        positive >= 0 && negative >= 0 && positive != negative || throw(ArgumentError(
            "cubic extension terminals must be distinct nonnegative nodes",
        ))
        conductance = Float64(linear_conductance_s)
        cubic = Float64(cubic_coefficient_a_per_v3)
        isfinite(conductance) && conductance >= 0.0 || throw(ArgumentError(
            "cubic extension conductance must be finite and nonnegative",
        ))
        isfinite(cubic) && cubic >= 0.0 || throw(ArgumentError(
            "cubic extension coefficient must be finite and nonnegative",
        ))
        conductance > 0.0 || cubic > 0.0 || throw(ArgumentError(
            "cubic extension requires a positive constitutive coefficient",
        ))
        provenance.nature === PhysicalModelParameter || throw(ArgumentError(
            "cubic extension provenance must describe physical model parameters",
        ))
        return new(positive, negative, conductance, cubic, provenance)
    end
end

extension_identity(::Type{NativeCubicCurrentBranch}) = _CUBIC_CURRENT_IDENTITY
extension_contract(::Type{NativeCubicCurrentBranch}) = ExtensionContract(
    _CUBIC_CURRENT_IDENTITY,
    :nonlinear_electrical,
    :instantaneous_emt,
    :switching_detailed,
    2,
    (:initialize, :nonlinear_current, :jacobian, :output, :checkpoint),
    _state_inventory(Dict(
        :algebraic => (:terminal_voltage, :terminal_current),
        :output => (:terminal_current, :absorbed_power),
        :checkpoint => (:identity_only,),
    )),
    "passive two-terminal branch with finite g>=0 and k>=0, not both zero",
    unsupported = (:active_power_generation, :hidden_state, :topology_change),
)

nonlinear_terminal_nodes(device::NativeCubicCurrentBranch) =
    (device.positive_node, device.negative_node)
nonlinear_device_formulation(::NativeCubicCurrentBranch) =
    NonlinearNetwork.PhysicalConstitutiveCurrent
nonlinear_device_provenance(device::NativeCubicCurrentBranch) = device.provenance

function nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    device::NativeCubicCurrentBranch,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::Float64,
)
    length(terminal_current_a) >= 2 && length(terminal_voltage_v) >= 2 ||
        throw(DimensionMismatch("cubic extension requires two terminal values"))
    size(terminal_jacobian_s, 1) >= 2 && size(terminal_jacobian_s, 2) >= 2 ||
        throw(DimensionMismatch("cubic extension requires a 2x2 Jacobian"))
    isfinite(time_s) || throw(ArgumentError("cubic extension evaluation time must be finite"))
    voltage = terminal_voltage_v[1] - terminal_voltage_v[2]
    current = device.linear_conductance_s * voltage +
        device.cubic_coefficient_a_per_v3 * voltage^3
    derivative = device.linear_conductance_s +
        3.0 * device.cubic_coefficient_a_per_v3 * voltage^2
    terminal_current_a[1] = current
    terminal_current_a[2] = -current
    terminal_jacobian_s[1, 1] = derivative
    terminal_jacobian_s[1, 2] = -derivative
    terminal_jacobian_s[2, 1] = -derivative
    terminal_jacobian_s[2, 2] = derivative
    return nothing
end

extension_checkpoint(device::NativeCubicCurrentBranch) =
    _component_checkpoint(device, (stateless = true,))
function restore_extension_checkpoint!(
    device::NativeCubicCurrentBranch,
    checkpoint::ExtensionComponentCheckpoint,
)
    state = _validate_checkpoint(device, checkpoint)
    state.stateless === true || _extension_fail(
        :incompatible_extension_checkpoint,
        :restore_checkpoint,
        _CUBIC_CURRENT_IDENTITY,
        "cubic extension checkpoint is not stateless",
    )
    return device
end

function extension_outputs(
    device::NativeCubicCurrentBranch,
    terminal_voltage_v::AbstractVector{<:Real},
    time_s::Real,
)
    length(terminal_voltage_v) >= 2 || throw(DimensionMismatch(
        "cubic extension output requires two terminal voltages",
    ))
    voltage = Float64(terminal_voltage_v[1] - terminal_voltage_v[2])
    current = device.linear_conductance_s * voltage +
        device.cubic_coefficient_a_per_v3 * voltage^3
    return (
        ExtensionOutputValue(:terminal_current, current, "A", time_s, :valid, _CUBIC_CURRENT_IDENTITY),
        ExtensionOutputValue(:absorbed_power, voltage * current, "W", time_s, :valid, _CUBIC_CURRENT_IDENTITY),
    )
end

const _SERIES_RL_IDENTITY = _example_identity(
    "series_rl_trapezoidal_companion",
    :NativeSeriesRLCompanion,
    "G=1/(R+2*L/h);Ihist=G*(vprev+(2*L/h-R)*iprev);i=G*v+Ihist",
)

mutable struct NativeSeriesRLCompanion <: AbstractExtensionElectricalDevice
    positive_node::Int
    negative_node::Int
    resistance_ohm::Float64
    inductance_h::Float64
    previous_current_a::Float64
    previous_voltage_v::Float64
    last_current_a::Float64
    provenance::ParameterProvenance

    function NativeSeriesRLCompanion(
        positive_node::Integer,
        negative_node::Integer,
        resistance_ohm::Real,
        inductance_h::Real;
        initial_current_a::Real = 0.0,
        initial_voltage_v::Real = 0.0,
        provenance::ParameterProvenance = ParameterProvenance(
            "caller-supplied series R-L parameters",
            "ohm, henry, ampere, and volt",
            "converted to finite Float64 SI values",
            "not supplied; caller retains uncertainty ownership",
            "R>=0, L>0, two distinct nonnegative nodes, and positive timestep",
            PhysicalModelParameter,
        ),
    )
        positive = Int(positive_node)
        negative = Int(negative_node)
        positive >= 0 && negative >= 0 && positive != negative || throw(ArgumentError(
            "series R-L extension terminals must be distinct nonnegative nodes",
        ))
        resistance, inductance, current, voltage = Float64.((
            resistance_ohm,
            inductance_h,
            initial_current_a,
            initial_voltage_v,
        ))
        all(isfinite, (resistance, inductance, current, voltage)) || throw(ArgumentError(
            "series R-L extension parameters and initial state must be finite",
        ))
        resistance >= 0.0 || throw(ArgumentError(
            "series R-L extension resistance must be nonnegative",
        ))
        inductance > 0.0 || throw(ArgumentError(
            "series R-L extension inductance must be positive",
        ))
        provenance.nature === PhysicalModelParameter || throw(ArgumentError(
            "series R-L extension provenance must describe physical model parameters",
        ))
        return new(
            positive,
            negative,
            resistance,
            inductance,
            current,
            voltage,
            current,
            provenance,
        )
    end
end

extension_identity(::Type{NativeSeriesRLCompanion}) = _SERIES_RL_IDENTITY
extension_contract(::Type{NativeSeriesRLCompanion}) = ExtensionContract(
    _SERIES_RL_IDENTITY,
    :electrical,
    :instantaneous_emt,
    :switching_detailed,
    2,
    (:initialize, :companion_stamp, :state_acceptance, :output, :checkpoint, :reusable_definition),
    _state_inventory(Dict(
        :continuous => (:inductor_current,),
        :algebraic => (:terminal_voltage, :terminal_current, :companion_conductance),
        :history => (:previous_voltage, :previous_current),
        :output => (:terminal_current, :stored_energy, :dissipated_power),
        :checkpoint => (:complete_component_state,),
    )),
    "two-terminal series R-L branch with R>=0, L>0, finite accepted state, and positive timestep",
    unsupported = (:zero_inductance, :negative_resistance, :hidden_acceptance),
)

extension_terminal_nodes(component::NativeSeriesRLCompanion) =
    (component.positive_node, component.negative_node)

function extension_companion(component::NativeSeriesRLCompanion, step_s::Real)
    step = Float64(step_s)
    isfinite(step) && step > 0.0 || throw(ArgumentError(
        "series R-L extension timestep must be finite and positive",
    ))
    inductive_resistance = 2.0 * component.inductance_h / step
    conductance = inv(component.resistance_ohm + inductive_resistance)
    history_current = conductance * (
        component.previous_voltage_v +
        (inductive_resistance - component.resistance_ohm) * component.previous_current_a
    )
    return (conductance_s = conductance, history_current_a = history_current)
end

function accept_extension_state!(
    component::NativeSeriesRLCompanion,
    terminal_voltage_v,
    step_s::Real,
)
    length(terminal_voltage_v) >= 2 || throw(DimensionMismatch(
        "series R-L extension acceptance requires two terminal voltages",
    ))
    voltage = Float64(terminal_voltage_v[1] - terminal_voltage_v[2])
    isfinite(voltage) || throw(ArgumentError(
        "series R-L extension accepted voltage must be finite",
    ))
    companion = extension_companion(component, step_s)
    current = companion.conductance_s * voltage + companion.history_current_a
    isfinite(current) || throw(ArgumentError(
        "series R-L extension accepted current must be finite",
    ))
    component.previous_voltage_v = voltage
    component.previous_current_a = current
    component.last_current_a = current
    return current
end

extension_checkpoint(component::NativeSeriesRLCompanion) = _component_checkpoint(component, (
    previous_current_a = component.previous_current_a,
    previous_voltage_v = component.previous_voltage_v,
    last_current_a = component.last_current_a,
))

function restore_extension_checkpoint!(
    component::NativeSeriesRLCompanion,
    checkpoint::ExtensionComponentCheckpoint,
)
    state = _validate_checkpoint(component, checkpoint)
    component.previous_current_a = state.previous_current_a
    component.previous_voltage_v = state.previous_voltage_v
    component.last_current_a = state.last_current_a
    return component
end

function extension_outputs(component::NativeSeriesRLCompanion, time_s::Real)
    stored_energy = 0.5 * component.inductance_h * component.last_current_a^2
    dissipated_power = component.resistance_ohm * component.last_current_a^2
    return (
        ExtensionOutputValue(:terminal_current, component.last_current_a, "A", time_s, :valid, _SERIES_RL_IDENTITY),
        ExtensionOutputValue(:stored_energy, stored_energy, "J", time_s, :valid, _SERIES_RL_IDENTITY),
        ExtensionOutputValue(:dissipated_power, dissipated_power, "W", time_s, :valid, _SERIES_RL_IDENTITY),
    )
end

mutable struct DirectedExtensionEvent
    quantity::Symbol
    threshold::Float64
    direction::Symbol
    tolerance::Float64
    priority::Int
    maximum_occurrences::Int
    occurrence_count::Int

    function DirectedExtensionEvent(
        quantity::Symbol,
        threshold::Real;
        direction::Symbol = :rising,
        tolerance::Real = 1.0e-12,
        priority::Integer = 0,
        maximum_occurrences::Integer = 1,
    )
        threshold_value = Float64(threshold)
        tolerance_value = Float64(tolerance)
        isfinite(threshold_value) || throw(ArgumentError("event threshold must be finite"))
        isfinite(tolerance_value) && tolerance_value > 0.0 || throw(ArgumentError(
            "event tolerance must be finite and positive",
        ))
        direction in (:rising, :falling, :either) || throw(ArgumentError(
            "event direction must be rising, falling, or either",
        ))
        maximum = Int(maximum_occurrences)
        maximum > 0 || throw(ArgumentError("event occurrence bound must be positive"))
        return new(
            quantity,
            threshold_value,
            direction,
            tolerance_value,
            Int(priority),
            maximum,
            0,
        )
    end
end

extension_event_value(event::DirectedExtensionEvent, quantity_value::Real) =
    Float64(quantity_value) - event.threshold

function extension_event_crossed(
    event::DirectedExtensionEvent,
    left_value::Real,
    right_value::Real,
)
    left = extension_event_value(event, left_value)
    right = extension_event_value(event, right_value)
    all(isfinite, (left, right)) || throw(ArgumentError("event surface values must be finite"))
    rising = left < -event.tolerance && right >= -event.tolerance
    falling = left > event.tolerance && right <= event.tolerance
    return event.direction === :rising ? rising :
           event.direction === :falling ? falling : rising || falling
end

function accept_extension_event!(event::DirectedExtensionEvent)
    event.occurrence_count < event.maximum_occurrences || _extension_fail(
        :extension_event_occurrence_limit,
        :accept_event,
        nothing,
        "extension event exceeded its declared occurrence bound",
    )
    event.occurrence_count += 1
    return event.occurrence_count
end

end

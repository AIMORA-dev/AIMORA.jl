module StudyCore

export STUDY_STATUSES,
       RESULT_STATUSES,
       ModelFidelity,
       LegacyDetailed,
       SwitchingDetailed,
       SwitchingStateEquivalent,
       AverageValue,
       FieldCoupledDetailed,
       NumericDomainBound,
       ModelValidityDomain,
       DynamicStateInventory,
       ContractQuantity,
       ScientificModelContract,
       ValidityViolation,
       ValidityAssessment,
       ValidityDomainError,
       assess_validity,
       assert_validity,
       StudyDescriptor,
       StudyRunRequest,
       ResultQuantity,
       StudyAssumption,
       StudyWarning,
       StudyResult,
       result_quantity,
       study_assumption,
       study_warning,
       study_result,
       add_quantity!,
       add_assumption!,
       add_warning!,
       result_ok,
       has_warnings,
       study_not_implemented

const STUDY_STATUSES = (:implemented, :prototype, :planned, :legacy_reference)
const RESULT_STATUSES = (:ok, :warning, :failed, :not_implemented, :invalid_input)

@enum ModelFidelity begin
    LegacyDetailed
    SwitchingDetailed
    SwitchingStateEquivalent
    AverageValue
    FieldCoupledDetailed
end

struct NumericDomainBound
    quantity::Symbol
    unit::String
    lower::Union{Nothing,Float64}
    upper::Union{Nothing,Float64}
    lower_inclusive::Bool
    upper_inclusive::Bool

    function NumericDomainBound(
        quantity::Symbol,
        unit::AbstractString,
        lower::Union{Nothing,Real},
        upper::Union{Nothing,Real},
        lower_inclusive::Bool,
        upper_inclusive::Bool,
    )
        isempty(String(quantity)) && throw(ArgumentError("domain-bound quantity cannot be empty"))
        isempty(strip(unit)) && throw(ArgumentError("domain-bound unit cannot be empty"))
        lower_value = lower === nothing ? nothing : Float64(lower)
        upper_value = upper === nothing ? nothing : Float64(upper)
        lower_value === nothing || isfinite(lower_value) ||
            throw(ArgumentError("domain-bound lower limit must be finite"))
        upper_value === nothing || isfinite(upper_value) ||
            throw(ArgumentError("domain-bound upper limit must be finite"))
        lower_value === nothing || upper_value === nothing || lower_value <= upper_value ||
            throw(ArgumentError("domain-bound lower limit must not exceed its upper limit"))
        return new(
            quantity,
            String(unit),
            lower_value,
            upper_value,
            lower_inclusive,
            upper_inclusive,
        )
    end
end

function NumericDomainBound(
    quantity::Symbol;
    unit::AbstractString,
    lower::Union{Nothing,Real} = nothing,
    upper::Union{Nothing,Real} = nothing,
    lower_inclusive::Bool = true,
    upper_inclusive::Bool = true,
)
    return NumericDomainBound(
        quantity,
        unit,
        lower,
        upper,
        lower_inclusive,
        upper_inclusive,
    )
end

struct ModelValidityDomain{B<:Tuple,U<:Tuple}
    id::Symbol
    description::String
    bounds::B
    unsupported_phenomena::U
end

function ModelValidityDomain(
    id::Symbol;
    description::AbstractString,
    bounds = (),
    unsupported_phenomena = (),
)
    isempty(String(id)) && throw(ArgumentError("validity-domain id cannot be empty"))
    isempty(strip(description)) && throw(ArgumentError("validity-domain description cannot be empty"))
    bound_tuple = Tuple(bounds)
    all(bound -> bound isa NumericDomainBound, bound_tuple) ||
        throw(ArgumentError("validity-domain bounds must be NumericDomainBound values"))
    bound_names = getfield.(bound_tuple, :quantity)
    length(unique(bound_names)) == length(bound_names) ||
        throw(ArgumentError("validity-domain bounds must have unique quantities"))
    unsupported = Tuple(Symbol.(unsupported_phenomena))
    return ModelValidityDomain(id, String(description), bound_tuple, unsupported)
end

struct DynamicStateInventory{D<:Tuple,A<:Tuple,X<:Tuple,H<:Tuple,S<:Tuple,R<:Tuple}
    differential::D
    algebraic::A
    discrete::X
    delayed_history::H
    scheduler::S
    random::R
end

function DynamicStateInventory(;
    differential = (),
    algebraic = (),
    discrete = (),
    delayed_history = (),
    scheduler = (),
    random = (),
)
    inventories = (
        Tuple(Symbol.(differential)),
        Tuple(Symbol.(algebraic)),
        Tuple(Symbol.(discrete)),
        Tuple(Symbol.(delayed_history)),
        Tuple(Symbol.(scheduler)),
        Tuple(Symbol.(random)),
    )
    all_names = Symbol[name for inventory in inventories for name in inventory]
    length(unique(all_names)) == length(all_names) ||
        throw(ArgumentError("dynamic-state names must be unique across state families"))
    return DynamicStateInventory(inventories...)
end

struct ContractQuantity
    key::Symbol
    unit::String
    base::String
    orientation::String

    function ContractQuantity(
        key::Symbol;
        unit::AbstractString,
        base::AbstractString = "absolute",
        orientation::AbstractString = "not_applicable",
    )
        isempty(String(key)) && throw(ArgumentError("contract quantity key cannot be empty"))
        isempty(strip(unit)) && throw(ArgumentError("contract quantity unit cannot be empty"))
        isempty(strip(base)) && throw(ArgumentError("contract quantity base cannot be empty"))
        isempty(strip(orientation)) && throw(ArgumentError("contract quantity orientation cannot be empty"))
        return new(key, String(unit), String(base), String(orientation))
    end
end

struct ScientificModelContract{D<:ModelValidityDomain,S<:DynamicStateInventory,I<:Tuple,O<:Tuple,A<:Tuple,M<:Tuple}
    id::Symbol
    capability::Symbol
    owner::String
    maturity::Symbol
    fidelity::ModelFidelity
    validity_domain::D
    state_inventory::S
    inputs::I
    outputs::O
    assumptions::A
    mutation_order::M
end

function ScientificModelContract(
    id::Symbol,
    capability::Symbol;
    owner::AbstractString,
    maturity::Symbol,
    fidelity::ModelFidelity,
    validity_domain::ModelValidityDomain,
    state_inventory::DynamicStateInventory,
    inputs = (),
    outputs = (),
    assumptions = (),
    mutation_order = (),
)
    maturity in STUDY_STATUSES ||
        throw(ArgumentError("scientific-contract maturity must be one of $(STUDY_STATUSES)"))
    isempty(String(id)) && throw(ArgumentError("scientific-contract id cannot be empty"))
    isempty(String(capability)) && throw(ArgumentError("scientific-contract capability cannot be empty"))
    isempty(strip(owner)) && throw(ArgumentError("scientific-contract owner cannot be empty"))
    input_tuple = Tuple(inputs)
    output_tuple = Tuple(outputs)
    all(quantity -> quantity isa ContractQuantity, (input_tuple..., output_tuple...)) ||
        throw(ArgumentError("scientific-contract inputs and outputs must be ContractQuantity values"))
    quantity_keys = getfield.((input_tuple..., output_tuple...), :key)
    length(unique(quantity_keys)) == length(quantity_keys) ||
        throw(ArgumentError("scientific-contract quantity keys must be unique"))
    assumption_tuple = Tuple(String.(assumptions))
    all(assumption -> !isempty(strip(assumption)), assumption_tuple) ||
        throw(ArgumentError("scientific-contract assumptions cannot be empty"))
    mutation_tuple = Tuple(Symbol.(mutation_order))
    isempty(mutation_tuple) && throw(ArgumentError("scientific-contract mutation order cannot be empty"))
    return ScientificModelContract(
        id,
        capability,
        String(owner),
        maturity,
        fidelity,
        validity_domain,
        state_inventory,
        input_tuple,
        output_tuple,
        assumption_tuple,
        mutation_tuple,
    )
end

struct ValidityViolation
    code::Symbol
    quantity::Symbol
    observed::Union{Nothing,Float64}
    message::String
end

struct ValidityAssessment{V<:Tuple}
    contract_id::Symbol
    requested_fidelity::ModelFidelity
    accepted_fidelity::ModelFidelity
    violations::V
end

struct ValidityDomainError <: Exception
    assessment::ValidityAssessment
end

Base.showerror(io::IO, error::ValidityDomainError) = print(
    io,
    "scientific contract $(error.assessment.contract_id) rejected the request: ",
    join((violation.message for violation in error.assessment.violations), "; "),
)

function _contract_value(values::NamedTuple, quantity::Symbol)
    hasproperty(values, quantity) || return false, nothing
    return true, getproperty(values, quantity)
end

function _contract_value(values::AbstractDict, quantity::Symbol)
    haskey(values, quantity) || return false, nothing
    return true, values[quantity]
end

function assess_validity(
    contract::ScientificModelContract,
    values::Union{NamedTuple,AbstractDict};
    requested_fidelity::ModelFidelity = contract.fidelity,
)
    violations = ValidityViolation[]
    requested_fidelity == contract.fidelity || push!(
        violations,
        ValidityViolation(
            :fidelity_mismatch,
            :fidelity,
            nothing,
            "requested fidelity $(requested_fidelity) is not provided; this owner is $(contract.fidelity)",
        ),
    )
    for bound in contract.validity_domain.bounds
        found, raw_value = _contract_value(values, bound.quantity)
        if !found
            push!(
                violations,
                ValidityViolation(
                    :missing_domain_quantity,
                    bound.quantity,
                    nothing,
                    "missing validity-domain quantity $(bound.quantity)",
                ),
            )
            continue
        end
        if !(raw_value isa Real)
            push!(
                violations,
                ValidityViolation(
                    :invalid_domain_quantity_type,
                    bound.quantity,
                    nothing,
                    "validity-domain quantity $(bound.quantity) must be real",
                ),
            )
            continue
        end
        value = Float64(raw_value)
        if !isfinite(value)
            push!(
                violations,
                ValidityViolation(
                    :nonfinite_domain_quantity,
                    bound.quantity,
                    value,
                    "validity-domain quantity $(bound.quantity) must be finite",
                ),
            )
            continue
        end
        lower_ok = bound.lower === nothing ||
                   (bound.lower_inclusive ? value >= bound.lower : value > bound.lower)
        upper_ok = bound.upper === nothing ||
                   (bound.upper_inclusive ? value <= bound.upper : value < bound.upper)
        lower_ok && upper_ok && continue
        push!(
            violations,
            ValidityViolation(
                :outside_validity_domain,
                bound.quantity,
                value,
                "validity-domain quantity $(bound.quantity)=$(value) $(bound.unit) is outside its declared bounds",
            ),
        )
    end
    return ValidityAssessment(
        contract.id,
        requested_fidelity,
        contract.fidelity,
        Tuple(violations),
    )
end

function assert_validity(assessment::ValidityAssessment)
    isempty(assessment.violations) || throw(ValidityDomainError(assessment))
    return true
end

struct StudyDescriptor
    id::Symbol
    name::String
    domain::Symbol
    status::Symbol
    source_path::String

    function StudyDescriptor(id::Symbol, name::String, domain::Symbol, status::Symbol, source_path::String)
        status in STUDY_STATUSES || error("Unsupported study status $(status). Allowed statuses: $(join(STUDY_STATUSES, ", ")).")
        return new(id, name, domain, status, source_path)
    end
end

Base.@kwdef struct StudyRunRequest
    study::Symbol
    case_path::Union{Nothing,String} = nothing
    output_dir::Union{Nothing,String} = nothing
end

struct ResultQuantity
    key::Symbol
    value::Any
    unit::String
    base::Union{Nothing,String}
    description::String
end

function ResultQuantity(
    key::Symbol,
    value;
    unit::AbstractString,
    base::Union{Nothing,AbstractString} = nothing,
    description::AbstractString = "",
)
    return ResultQuantity(
        key,
        value,
        String(unit),
        base === nothing ? nothing : String(base),
        String(description),
    )
end

struct StudyAssumption
    key::Symbol
    value::Any
    description::String
end

function StudyAssumption(key::Symbol, value; description::AbstractString = "")
    return StudyAssumption(key, value, String(description))
end

struct StudyWarning
    code::Symbol
    severity::Symbol
    message::String
end

function StudyWarning(code::Symbol, message::AbstractString; severity::Symbol = :warning)
    severity in (:info, :warning, :error) || error("Unsupported warning severity $(severity).")
    return StudyWarning(code, severity, String(message))
end

Base.@kwdef mutable struct StudyResult
    study::Symbol
    status::Symbol
    quantities::Dict{Symbol,ResultQuantity} = Dict{Symbol,ResultQuantity}()
    assumptions::Vector{StudyAssumption} = StudyAssumption[]
    warnings::Vector{StudyWarning} = StudyWarning[]
    metadata::Dict{Symbol,Any} = Dict{Symbol,Any}()

    function StudyResult(
        study::Symbol,
        status::Symbol,
        quantities::Dict{Symbol,ResultQuantity},
        assumptions::Vector{StudyAssumption},
        warnings::Vector{StudyWarning},
        metadata::Dict{Symbol,Any},
    )
        status in RESULT_STATUSES || error("Unsupported result status $(status). Allowed statuses: $(join(RESULT_STATUSES, ", ")).")
        if status == :ok && any(w -> w.severity != :info, warnings)
            status = :warning
        end
        return new(study, status, quantities, assumptions, warnings, metadata)
    end
end

result_quantity(args...; kwargs...) = ResultQuantity(args...; kwargs...)
study_assumption(args...; kwargs...) = StudyAssumption(args...; kwargs...)
study_warning(args...; kwargs...) = StudyWarning(args...; kwargs...)

function study_result(
    study::Symbol;
    status::Symbol = :ok,
    quantities = ResultQuantity[],
    assumptions = StudyAssumption[],
    warnings = StudyWarning[],
    metadata = Dict{Symbol,Any}(),
)
    quantity_map = Dict{Symbol,ResultQuantity}()
    for quantity in quantities
        quantity_map[quantity.key] = quantity
    end
    return StudyResult(
        study = study,
        status = status,
        quantities = quantity_map,
        assumptions = collect(assumptions),
        warnings = collect(warnings),
        metadata = Dict{Symbol,Any}(metadata),
    )
end

function add_quantity!(result::StudyResult, quantity::ResultQuantity)
    result.quantities[quantity.key] = quantity
    return result
end

function add_assumption!(result::StudyResult, assumption::StudyAssumption)
    push!(result.assumptions, assumption)
    return result
end

function add_warning!(result::StudyResult, warning::StudyWarning)
    push!(result.warnings, warning)
    if result.status == :ok && warning.severity != :info
        result.status = :warning
    end
    return result
end

result_ok(result::StudyResult) = result.status in (:ok, :warning)
has_warnings(result::StudyResult) = any(w -> w.severity != :info, result.warnings)

function study_not_implemented(desc::StudyDescriptor)
    error("Study $(desc.id) ($(desc.name)) has status $(desc.status) and is not implemented as a runnable study yet.")
end

end

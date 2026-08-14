module TransformerApparatus

using ..StudyCore: ContractQuantity,
                   DynamicStateInventory,
                   FieldCoupledDetailed,
                   LegacyDetailed,
                   ModelValidityDomain,
                   NumericDomainBound,
                   ParameterProvenance,
                   PhysicalModelParameter,
                   ScientificModelContract,
                   SwitchingStateEquivalent
using ..TransformerParameters: MultiphaseTransformerParameterResult,
                               SaturableTransformerParameterResult
using ..NonlinearNetwork: AbstractNonlinearCurrentDevice,
                          NonlinearDeviceEventSurface,
                          PhysicalConstitutiveCurrent
import ..NonlinearNetwork: accept_nonlinear_device_state!,
                           apply_nonlinear_device_event!,
                           finish_nonlinear_device_step!,
                           nonlinear_current_jacobian!,
                           nonlinear_device_event_candidate_time,
                           nonlinear_device_event_surfaces,
                           nonlinear_device_event_value,
                           nonlinear_device_formulation,
                           nonlinear_device_provenance,
                           nonlinear_terminal_nodes,
                           prepare_nonlinear_device_step!
using LinearAlgebra
using SHA

include(joinpath(@__DIR__, "transformers", "connection_and_tiers.jl"))
include(joinpath(@__DIR__, "transformers", "magnetic_core.jl"))
include(joinpath(@__DIR__, "transformers", "apparatus_models.jl"))
include(joinpath(@__DIR__, "transformers", "event_contracts.jl"))
include(joinpath(@__DIR__, "transformers", "fixed_step_runtime.jl"))
include(joinpath(@__DIR__, "transformers", "events.jl"))
include(joinpath(@__DIR__, "transformers", "results.jl"))

end

module ModernMachines

using ..StudyCore: ContractQuantity,
                   DynamicStateInventory,
                   FieldCoupledDetailed,
                   ModelValidityDomain,
                   NumericDomainBound,
                   ParameterProvenance,
                   PhysicalModelParameter,
                   ScientificModelContract
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

include(joinpath(@__DIR__, "machines", "modern_machine_contracts.jl"))
include(joinpath(@__DIR__, "machines", "modern_machine_electromagnetics.jl"))
include(joinpath(@__DIR__, "machines", "modern_machine_controls_and_shaft.jl"))
include(joinpath(@__DIR__, "machines", "modern_machine_runtime.jl"))
include(joinpath(@__DIR__, "machines", "modern_machine_results.jl"))

end

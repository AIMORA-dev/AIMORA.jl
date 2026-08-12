module AIMORAAIMORAProjectExt

using AIMORA
using AIMORAProject

import AIMORA.NativeExtensions: ExtensionFailure,
                                ExtensionIdentity,
                                extension_identity,
                                resolve_extension

const PROJECT_SERVICE_TO_RUNTIME = Dict(
    AIMORAProject.ExtensionInitializationService => :initialize,
    AIMORAProject.ExtensionNonlinearCurrentService => :nonlinear_current,
    AIMORAProject.ExtensionJacobianService => :jacobian,
    AIMORAProject.ExtensionCompanionStampService => :companion_stamp,
    AIMORAProject.ExtensionStateAcceptanceService => :state_acceptance,
    AIMORAProject.ExtensionEventService => :event,
    AIMORAProject.ExtensionSampledTaskService => :sampled_task,
    AIMORAProject.ExtensionSourceService => :source,
    AIMORAProject.ExtensionOutputService => :output,
    AIMORAProject.ExtensionCheckpointService => :checkpoint,
    AIMORAProject.ExtensionReusableDefinitionService => :reusable_definition,
)

function extension_identity(declaration::AIMORAProject.ExtensionDeclaration)
    implementation = declaration.implementation
    return ExtensionIdentity(
        implementation.package_uuid,
        implementation.package.namespace.value,
        implementation.package.name.value,
        implementation.package.version,
        Symbol(implementation.symbol.value),
        declaration.api_version,
        implementation.content_hash.sha256,
    )
end

function _binding_failure(code::Symbol, declaration, message::String)
    throw(ExtensionFailure(code, :resolve_project_declaration, extension_identity(declaration), message))
end

function resolve_extension(
    registry::AIMORA.NativeExtensions.ExtensionRegistry,
    declaration::AIMORAProject.ExtensionDeclaration,
)
    registration = resolve_extension(registry, extension_identity(declaration))
    contract = registration.contract
    declaration.representation == AIMORAProject.InstantaneousEMT || _binding_failure(
        :unsupported_extension_representation,
        declaration,
        "project extension representation is not admitted by the runtime contract",
    )
    declaration.fidelity == AIMORAProject.SwitchingDetailed || _binding_failure(
        :unsupported_extension_fidelity,
        declaration,
        "project extension fidelity is not admitted by the runtime contract",
    )
    contract.representation === :instantaneous_emt || _binding_failure(
        :extension_representation_mismatch,
        declaration,
        "registered extension representation differs from the project declaration",
    )
    contract.fidelity === :switching_detailed || _binding_failure(
        :extension_fidelity_mismatch,
        declaration,
        "registered extension fidelity differs from the project declaration",
    )
    contract.terminal_count == length(declaration.terminals) || _binding_failure(
        :extension_terminal_count_mismatch,
        declaration,
        "registered extension terminal count differs from the project declaration",
    )
    project_services = Set(PROJECT_SERVICE_TO_RUNTIME[service] for service in declaration.services)
    project_services == Set(contract.services) || _binding_failure(
        :extension_service_mismatch,
        declaration,
        "registered extension services differ from the project declaration",
    )
    return registration
end

end

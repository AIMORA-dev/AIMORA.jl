module AIMORAAIMORAProjectExt

using AIMORA
using AIMORAProject

import AIMORA.NativeExtensions: ExtensionFailure,
                                ExtensionIdentity,
                                extension_identity,
                                resolve_extension
import AIMORA.EMTTaskPlatform: CarrierEMTTask,
                               ConverterControlEMTTask,
                               EMTLogicalTime,
                               EMTTaskEffect,
                               EMTTaskFamily,
                               EMTTaskSpec,
                               InterfaceEMTTask,
                               InvalidateEMTInterface,
                               InvalidateEMTOutput,
                               InvalidateEMTPowerHistory,
                               InvalidateEMTTopology,
                               MechanicalEMTTask,
                               ProtectionEMTTask,
                               SourceEMTTask,
                               ThermalEMTTask,
                               UserDefinedEMTTask,
                               emt_task_spec

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

const PROJECT_TASK_FAMILY_TO_ENGINE = Dict{AIMORAProject.ControlTaskFamily,EMTTaskFamily}(
    AIMORAProject.ProtectionControlTask => ProtectionEMTTask,
    AIMORAProject.CarrierControlTask => CarrierEMTTask,
    AIMORAProject.ConverterControlTask => ConverterControlEMTTask,
    AIMORAProject.MechanicalControlTask => MechanicalEMTTask,
    AIMORAProject.SourceControlTask => SourceEMTTask,
    AIMORAProject.ThermalControlTask => ThermalEMTTask,
    AIMORAProject.InterfaceControlTask => InterfaceEMTTask,
    AIMORAProject.UserDefinedControlTask => UserDefinedEMTTask,
)

const PROJECT_TASK_EFFECT_TO_ENGINE = Dict{AIMORAProject.ControlTaskInvalidation,EMTTaskEffect}(
    AIMORAProject.InvalidateControlPowerHistory => InvalidateEMTPowerHistory,
    AIMORAProject.InvalidateControlTopology => InvalidateEMTTopology,
    AIMORAProject.InvalidateControlInterface => InvalidateEMTInterface,
    AIMORAProject.InvalidateControlOutput => InvalidateEMTOutput,
)

function _project_logical_time(
    project::AIMORAProject.CanonicalProject,
    value::AIMORAProject.PhysicalValue{AIMORAProject.ScalarQuantity},
)
    quantity = AIMORAProject.convert_quantity(
        project.units,
        value.quantity,
        AIMORAProject.UnitId("s"),
    )
    rational = AIMORAProject.exact_rational(quantity.value)
    return EMTLogicalTime(rational.numerator, rational.denominator)
end

"""Resolve one inert Project task declaration into the dependency-light public EMT task contract."""
function emt_task_spec(
    project::AIMORAProject.CanonicalProject,
    declaration::AIMORAProject.ControlTaskDeclaration,
)
    return EMTTaskSpec(
        declaration.task.value,
        PROJECT_TASK_FAMILY_TO_ENGINE[declaration.family],
        _project_logical_time(project, declaration.epoch),
        _project_logical_time(project, declaration.period),
        _project_logical_time(project, declaration.phase),
        _project_logical_time(project, declaration.computational_delay);
        priority = declaration.priority,
        read_resources = getfield.(collect(declaration.read_resources), :value),
        write_resources = getfield.(collect(declaration.write_resources), :value),
        predecessors = getfield.(collect(declaration.predecessors), :value),
        effects = EMTTaskEffect[
            PROJECT_TASK_EFFECT_TO_ENGINE[effect] for effect in declaration.invalidations
        ],
    )
end

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

module AIMORA

include("solver_api/backend.jl")

# Public, dependency-light engineering core.
include("core/study.jl")
include("core/validation.jl")
include("core/inputs.jl")
include("core/tables.jl")
include("core/project.jl")
include("io/project_io.jl")
include("models/inverter.jl")
include("models/switch_detailed_vsc.jl")
include("models/nonlinear_network.jl")
include("extensions/native_components.jl")
include("models/inverter_assets.jl")
include("models/transformer_parameters.jl")
include("io/transformer_parameter_input.jl")
include("studies/transformer_parameters.jl")
include("io/transformer_parameter_report.jl")
include("RealtimeLoop.jl")
include("Figures.jl")
include("studies/catalog.jl")
include("studies/input_profiles.jl")
include("studies/power_flow.jl")
include("studies/short_circuit.jl")
include("studies/protection.jl")
include("studies/arc_flash.jl")

end

module AIMORA

export require_solver, solver_available, solver_status

const _SOLVER_SOURCE_DIR = joinpath(@__DIR__, "solvers")
const _REQUIRED_SOLVER_FILES = (
    "companion.jl",
    "nodal.jl",
    "timestep.jl",
    "over16_timestep_integration.jl",
)

if isdir(_SOLVER_SOURCE_DIR)
    Base.include_dependency(_SOLVER_SOURCE_DIR)
else
    Base.include_dependency(dirname(_SOLVER_SOURCE_DIR))
end

_solver_sources_present() = all(
    name -> isfile(joinpath(_SOLVER_SOURCE_DIR, name)),
    _REQUIRED_SOLVER_FILES,
)

const _SOLVER_LOADED = _solver_sources_present()

"""
Return `true` when this AIMORA module was loaded with the complete licensed
numerical component.
"""
solver_available() = _SOLVER_LOADED

"""
Report whether this installation is the public `:open_core` or a licensed
`:full_engine`, and identify how the numerical component was supplied.
"""
function solver_status()
    available = solver_available()
    return (
        available = available,
        mode = available ? :full_engine : :open_core,
        source = available ? :licensed_component : :not_installed,
    )
end

"""
Return `nothing` when the proprietary solver is available; otherwise throw an
error explaining how an authorized developer can initialize it.
"""
function require_solver()
    solver_available() && return nothing
    error(
        "The licensed AIMORA numerical component is not installed. " *
        "Use the installation instructions supplied with your distribution.",
    )
end

# Public, dependency-light engineering core.
include("core/study.jl")
include("core/validation.jl")
include("core/inputs.jl")
include("core/tables.jl")
include("core/project.jl")
include("io/project_io.jl")
include("models/inverter.jl")
include("models/switch_detailed_vsc.jl")
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

# The complete EMT runtime is loaded only in an authorized checkout. No solver
# implementation is committed to this public repository.
if solver_available()
    include("solvers/companion.jl")
    include("solvers/timestep.jl")
    include("models/branches.jl")
    include("models/sources.jl")
    include("models/switches.jl")
    include("models/nonlinear.jl")
    include("models/tacs.jl")
    include("models/lines.jl")
    include("models/machines.jl")
    include("solvers/over16_timestep_integration.jl")
    include("solvers/nodal.jl")
    include("io/deck_parser.jl")
    include("studies/cable_constants.jl")
    include("studies/line_constants.jl")
    include("studies/emt.jl")
    include("io/reports.jl")
end

end

module AIMORA

export require_solver, solver_available, solver_status

const _SOLVER_SOURCE_DIR = joinpath(@__DIR__, "julia", "solvers")
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
Return `true` when this AIMORA module was loaded with the complete proprietary
solver source component.
"""
solver_available() = _SOLVER_LOADED

function solver_status()
    available = solver_available()
    return (
        available = available,
        mode = available ? :full_engine : :open_core,
        source = available ? :private_git_submodule : :not_installed,
    )
end

function require_solver()
    solver_available() && return nothing
    error(
        "The proprietary AIMORA solver is not installed. Authorized developers " *
        "must initialize the private src/julia/solvers Git submodule.",
    )
end

# Public, dependency-light engineering core.
include("julia/core/study.jl")
include("julia/core/validation.jl")
include("julia/core/inputs.jl")
include("julia/core/tables.jl")
include("julia/core/project.jl")
include("julia/io/project_io.jl")
include("julia/models/inverter.jl")
include("julia/models/inverter_assets.jl")
include("julia/RealtimeLoop.jl")
include("julia/Figures.jl")
include("julia/studies/catalog.jl")
include("julia/studies/input_profiles.jl")
include("julia/studies/power_flow.jl")
include("julia/studies/short_circuit.jl")
include("julia/studies/protection.jl")
include("julia/studies/arc_flash.jl")

# The complete EMT runtime is loaded only in an authorized checkout. No solver
# implementation is committed to this public repository.
if solver_available()
    include("julia/solvers/companion.jl")
    include("julia/solvers/timestep.jl")
    include("julia/models/branches.jl")
    include("julia/models/sources.jl")
    include("julia/models/switches.jl")
    include("julia/models/nonlinear.jl")
    include("julia/models/tacs.jl")
    include("julia/models/lines.jl")
    include("julia/models/machines.jl")
    include("julia/solvers/over16_timestep_integration.jl")
    include("julia/solvers/nodal.jl")
    include("julia/io/deck_parser.jl")
    include("julia/studies/cable_constants.jl")
    include("julia/studies/line_constants.jl")
    include("julia/studies/emt.jl")
    include("julia/io/reports.jl")
end

end

module Machines
include(joinpath(@__DIR__, "machines", "synchronous_machine_core.jl"))
include(joinpath(@__DIR__, "machines", "universal_machine_runtime.jl"))
include(joinpath(@__DIR__, "machines", "machine_equations_and_terminal_state.jl"))

end

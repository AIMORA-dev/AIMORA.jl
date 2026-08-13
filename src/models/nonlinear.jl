module Nonlinear
include(joinpath(@__DIR__, "nonlinear", "transformer_models.jl"))
include(joinpath(@__DIR__, "nonlinear", "switching_resistor.jl"))
include(joinpath(@__DIR__, "nonlinear", "triggered_timed_resistance.jl"))
include(joinpath(@__DIR__, "nonlinear", "piecewise_nonlinear_inductor.jl"))
include(joinpath(@__DIR__, "nonlinear", "preprocessing.jl"))
include(joinpath(@__DIR__, "nonlinear", "network_devices.jl"))
include(joinpath(@__DIR__, "nonlinear", "inverse_columns_and_compensation.jl"))
include(joinpath(@__DIR__, "nonlinear", "power_semiconductor_fidelity.jl"))
include(joinpath(@__DIR__, "nonlinear", "power_semiconductor_switches.jl"))
include(joinpath(@__DIR__, "nonlinear", "power_semiconductor_bridges.jl"))
include(joinpath(@__DIR__, "nonlinear", "power_semiconductor_topologies.jl"))
include(joinpath(@__DIR__, "nonlinear", "simultaneous_zno_and_updates.jl"))

end

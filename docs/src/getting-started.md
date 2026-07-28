# Getting Started

## Public Open Core

```julia
using Pkg
Pkg.add(url = "https://github.com/AIMORA-dev/AIMORA.jl")

using AIMORA
AIMORA.solver_status()
AIMORA.StudyCatalog.available_studies()
```

Run the public package tests with:

```bash
julia --project=. test/runtests.jl
```

The inverter example is available without the private solver:

```bash
julia --project=. examples/inverter/run.jl
```

## Authorized Full Engine

```bash
git clone --recurse-submodules git@github.com:AIMORA-dev/AIMORA.jl.git
cd AIMORA.jl
julia --project=. test/runtests.jl
```

`AIMORA.solver_status().mode` must then be `:full_engine`. EMT, line-model, and
deck examples require this installation mode.

## Legacy Reference

Build and run historical reference decks through
`BPAEMTPReference.jl`. The Fortran executable is an external validation oracle,
not a production dependency.

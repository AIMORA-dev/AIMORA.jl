# AIMORA.jl

**Analytical Integration for Multiphysics Operations and Response Analysis**

AIMORA is a Julia-native power-system engineering platform. The public package
owns study contracts, project and asset schemas, validation, model metadata,
reporting utilities, and open engineering models. Its current production EMT
runtime is available to authorized installations through a separate
proprietary solver repository.

## Installation modes

### Public open core

```julia
using Pkg
Pkg.add(url = "https://github.com/AIMORA-dev/AIMORA.jl")

using AIMORA
AIMORA.solver_status()
```

The public package loads without private files. Project data, study input
profiles, asset tables, validation, inverter examples, and future study
interfaces remain available.

### Authorized full engine

```bash
git clone --recurse-submodules git@github.com:AIMORA-dev/AIMORA.jl.git
cd AIMORA.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```

The private Git submodule is mounted at `src/julia/solvers`. Its source remains
readable and testable on authorized development machines, but no solver blob
exists in this public repository's history.

## Current scientific scope

- Julia-only electromagnetic-transient execution for the accepted
  `aimora_bpa_emtp_replacement_v1` target set in authorized full-engine
  installations
- Typed project, scenario, study, result, validation, and asset-table APIs
- Open inverter model and fixed-step demonstration
- Overhead and cable line-constants implementations in the full engine
- Explicit architecture contracts for power flow, short circuit, protection,
  and arc flash; these studies are not yet claimed as implemented
- Optional CUDA fixed-admittance batching when the private solver and a
  functional CUDA device are present

CPU remains the default. Backend-neutral support for AMD, Intel, and other GPU
families is an architectural target, not a current compatibility claim.

## Package layout

```text
src/AIMORA.jl              package entrypoint
src/julia/core/            public study, model, table, and validation contracts
src/julia/models/          public equipment and component implementations
src/julia/studies/         public study workflows and declared roadmap
src/julia/io/              public input/output and reporting implementations
src/julia/solvers/         private Git submodule; not present in public history
examples/                  runnable public and authorized full-engine examples
test/                      public-core and private-integration package tests
```

## Related repositories

- [BPAEMTPReference.jl](https://github.com/AIMORA-dev/BPAEMTPReference.jl):
  compiled historical reference and Julia wrapper
- `AIMORASolvers.jl`: private solver source component
- `AIMORAValidation`: private cross-package qualification workspace
- [AIMORADocs](https://github.com/AIMORA-dev/AIMORADocs): unified documentation
- [AIMORACases.jl](https://github.com/AIMORA-dev/AIMORACases.jl): example and benchmark
  cases
- [AIMORACatalogs.jl](https://github.com/AIMORA-dev/AIMORACatalogs.jl): open
  catalog schemas and distributable model data

## Licence

Public AIMORA code is available under the MIT licence. The solver submodule is
proprietary and governed by its separate licence. Historical BPA source and
compiled-reference tooling live in their own repository with preserved
provenance.

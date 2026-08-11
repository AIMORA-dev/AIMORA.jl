# AIMORA.jl

**Analytical Integration for Multiphase Overvoltage and Response Analysis**

AIMORA is a Julia-native power-system engineering platform. The public package owns study contracts, project and asset schemas, validation, model metadata, reporting utilities, and open engineering models. Production numerical capabilities are supplied through separately licensed distributions.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/AIMORA-dev/AIMORA.jl")

using AIMORA
AIMORA.solver_status()
```

The public package loads independently. Project data, study input profiles, asset tables, validation, inverter examples, and future study interfaces remain available. A production-solve request without an explicitly registered backend returns `AIMORA.SolverUnavailableResult`. Installation and explicit activation instructions for licensed capabilities are provided with the corresponding distribution.

## Current scientific scope

- Julia-only electromagnetic-transient execution for the accepted `aimora_bpa_emtp_replacement_v1` target set in the production distribution
- Typed project, scenario, study, result, validation, and asset-table APIs
- Open inverter model and fixed-step demonstration
- Overhead and cable line-constants implementations in the full engine
- Explicit architecture contracts for power flow, short circuit, protection, and arc flash; these studies are not yet claimed as implemented
- Optional CUDA fixed-admittance batching when the numerical backend and a functional CUDA device are present

CPU remains the default. Backend-neutral support for AMD, Intel, and other GPU families is an architectural target, not a current compatibility claim.

## Package layout

```text
src/AIMORA.jl              package entrypoint
src/solver_api/             public backend contracts and typed availability results
src/core/            public study, model, table, and validation contracts
src/models/          public equipment and component implementations
src/studies/         public study workflows and declared roadmap
src/io/              public input/output and reporting implementations
test/                      public and production integration package tests
check.jl                   package structure and publication-boundary check
Makefile                   test and check commands
```

## Related repositories

- [AIMORAResources](https://github.com/AIMORA-dev/AIMORAResources): compiled historical reference, unified documentation, public cases, catalogues, independent reference models, report templates, teaching material, and provenance
- [AIMORAPlatform](https://github.com/AIMORA-dev/AIMORAPlatform): project, format, layout, service, visual, reporting, and symbol packages

## Licence

This repository's AIMORA-authored content is distributed under the PolyForm Noncommercial License 1.0.0. Research, education, personal study, public-interest noncommercial use, and other purposes permitted by that licence are free; commercial use requires a separate written agreement with Ahmed Elkholy <ahmed_elkholy@f-eng.tanta.edu.eg>. There is no licence key, telemetry, or technical feature restriction in this public package. Clearly identified third-party material retains its own terms, and copies received under an earlier licence retain those prior grants.

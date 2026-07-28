# AIMORA.jl Documentation

AIMORA is a Julia-native power-system engineering platform. The public package
owns typed project and study contracts, engineering validation, asset schemas,
open model implementations, and reporting boundaries.

Authorized installations add the proprietary numerical solver through the
private `src/julia/solvers` Git submodule. Compiled BPA EMTP remains a separate
external reference and is never called by the production Julia timestep loop.

## Current Boundary

The authorized full engine is accepted for the named
`aimora_bpa_emtp_replacement_v1` target set at 93/93. This is not a claim that
arbitrary future decks or planned non-EMT studies are implemented.

## Documentation Map

- [Getting Started](getting-started.md)
- [Architecture](architecture.md)
- [Studies](studies.md)

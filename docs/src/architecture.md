# Architecture

AIMORA separates public engineering ownership from proprietary numerical
implementation.

```text
src/AIMORA.jl
src/julia/core/       study, validation, project, and table contracts
src/julia/io/         project and result input/output
src/julia/models/     equipment and component models
src/julia/studies/    study workflows and declared maturity
src/julia/solvers/    private Git submodule
examples/             runnable public and authorized examples
test/                 open-core and private-integration package tests
```

The package entrypoint always loads the public core. When the complete private
solver checkout is present, it also loads the EMT network, parser, timestep,
line, machine, and report runtime.

Core engineering objects include:

- `Project`, `Case`, `Scenario`, `Revision`, and `StudySettings`
- `StudyDescriptor`, `StudyRunRequest`, and `StudyResult`
- `TableSchema`, `ColumnSpec`, `ValidationIssue`, and `ValidationResult`
- per-study required and optional input profiles
- open inverter and asset-table models

The GUI and catalogs must consume these same schemas. They must not create a
second, GUI-only definition of units, defaults, or required fields.

## Source Boundaries

- The public repository contains no solver blob in any commit.
- Authorized local clones restore the private solver through Git submodules.
- The Fortran reference lives in `BPAEMTPReference.jl`.
- Cross-package oracle comparison lives in the private validation workspace.
- CPU is the default. CUDA support is conditional on both the private solver
  and a functional device. Other GPU families remain future compatibility
  targets until validated.

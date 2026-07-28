# Studies

The platform grows from a shared physical model and study-specific parameter
facets. A study must validate its own required inputs rather than forcing every
asset to carry unrelated data.

## Catalog Status

| Study ID | Catalog status | Runtime availability |
|---|---|---|
| `:emt` | `:implemented` | Authorized full engine |
| `:line_constants` | `:implemented` | Authorized full engine |
| `:cable_constants` | `:implemented` | Authorized full engine |
| `:power_flow` | `:planned` | Not implemented |
| `:short_circuit` | `:planned` | Not implemented |
| `:protection` | `:planned` | Not implemented |
| `:arc_flash` | `:planned` | Not implemented |

The catalog status describes scientific implementation maturity. Installation
availability is reported separately by `AIMORA.solver_status()`.

## Shared Assets and Study Facets

A physical line, machine, transformer, or converter has one identity,
topology, rating, and provenance. Each study adds a typed parameter facet:

- power-flow impedance, limits, and controls;
- sequence data and grounding for fault studies;
- dynamic controls and states for RMS studies;
- geometry, frequency dependence, histories, and switching state for EMT.

Simpler data may be derived from detailed data only when a documented physical
conversion exists. AIMORA must never silently invent missing detailed
parameters.

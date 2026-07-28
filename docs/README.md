# AIMORA documentation

The package keeps a compact Markdown manual under `docs/src`.

- `index.md`: package boundary and navigation
- `getting-started.md`: public and authorized installation
- `architecture.md`: public/private source ownership
- `studies.md`: study maturity and model facets

The unified hosted manual is maintained in the separate public `AIMORADocs`
repository.

Run the dependency-light consistency check with:

```bash
julia --project=. docs/check.jl
```

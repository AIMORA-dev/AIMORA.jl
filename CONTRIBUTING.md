# Contributing

Public contributions should target study contracts, models, data schemas, documentation, examples, or tests that do not disclose proprietary solver source.

Before opening a pull request:

```bash
julia --project=. test/runtests.jl
git diff --check
```

Public pull requests and fork CI never receive access to the private solver. Solver changes are reviewed and tested in the private repository and the private integration workspace.

## Contribution licence

By submitting a contribution, you certify that you have the right to provide it and agree that you retain your copyright while granting Ahmed Elkholy a perpetual, worldwide, irrevocable, nonexclusive, royalty-free right to use, reproduce, modify, distribute, sublicense, and offer the contribution under this repository's PolyForm Noncommercial terms and under separately negotiated commercial agreements. Clearly identified third-party material is accepted only with compatible written redistribution and commercial-use rights; submitting it does not override its original licence.

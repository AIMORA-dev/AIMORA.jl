# Contributing

Public contributions should target study contracts, models, data schemas,
documentation, examples, or tests that do not disclose proprietary solver
source.

Before opening a pull request:

```bash
julia --project=. test/runtests.jl
git diff --check
```

Public pull requests and fork CI never receive access to the private solver.
Solver changes are reviewed and tested in the private repository and the
private integration workspace.

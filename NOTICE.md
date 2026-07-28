# AIMORA notices

- AIMORA's Julia production design is informed by and validated against the
  historical Bonneville Power Administration Electromagnetic Transients
  Program.
- The compiled Fortran program is an external reference oracle and is never
  part of the production Julia timestep loop.
- `src/solvers` is a separate proprietary Git repository. It is not
  licensed under the public repository's MIT licence.
- Product and study maturity statements apply only to named, tested target
  sets. Future studies listed in the catalog are not implied to be complete.

# Security

<!-- guidelines:start (shared development standard 0.1.0; edit it there, not here) -->
## Reporting a problem

Report a vulnerability privately through GitHub: **Security > Report a
vulnerability** on this repository. Please do not open a public issue for it.
You should hear back within a week.

## Supported versions

Only the latest release gets fixes.
<!-- guidelines:end -->

This is a geometry library: it reads coordinate tables and CSV files you give
it, takes windIO data as already-parsed mappings, and writes figures and mesh
files. It runs no network services.

## Things to know

- The kernel does not parse YAML itself; load windIO files with
  `yaml.safe_load` before passing them in.
- `GEOS_LIBRARY_PATH` loads a shared library from the path you set. Point it
  only at a GEOS build you trust.

# Float3109s.jl

Reference Julia implementation of an analytic floating-point geometry model for
parameterized P3109-style formats.

## Included

- `src/Float3109s.jl`
- `test/runtests.jl`
- `docs/make.jl`
- Documenter pages in `docs/src/`
- `.github/workflows/ci.yml`

## Usage

```julia
using Pkg
Pkg.develop(path=".")
Pkg.test()

using Float3109s
f = P3109Format(8, 3, 0, 1)
decode(f, 4)
format_summary(f)
```

## Notes

This package is intended as a reference and validation layer. It uses exact
rational arithmetic and explicit checks rather than machine-float shortcuts.

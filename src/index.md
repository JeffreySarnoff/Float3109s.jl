# Float3109s.jl

A reference Julia implementation of an analytic floating-point geometry model
for parameterized P3109-style formats.

## What this package does

Float3109s.jl provides exact, closed-form descriptions of every structural
property of a P3109 floating-point format.  Given a format's bitwidth `K` and
precision `P` — together with a signedness and domain choice — the package can:

- **Count** every class of value (normals, subnormals, infinities, ...)
- **Locate** the code point of any boundary value (smallest subnormal,
  largest normal, NaN, ±Inf, ...)
- **Decode** any code point to its exact rational value using `ClosedRational`
- **Enumerate** every representable value for small formats

All arithmetic is exact.  There is no rounding and no floating-point
approximation anywhere in the pipeline.

## Design intent

This package is structured as a **correctness reference** first.  It is meant
to serve as a test oracle and specification companion for the IEEE P3109
working group.  Optimized kernels can be layered on top while preserving these
routines as the ground truth.

## Quick start

```julia
using Float3109s

# An 8-bit, precision-4, signed, extended (has ±Inf) format
fmt = Format{is_signed, is_extended}(8, 4)

# Structural counts
nFiniteValuesOf(fmt)       # total finite values (including zero)
nPosNormalsOf(fmt)         # positive normal values
nBinadesOf(fmt)            # number of normal binades

# Boundary code points
cp_nan(fmt)                # code point of NaN
cp_pos_normal_max(fmt)     # code point of the largest positive normal

# Exact value decoding
FiniteValueOf(fmt, 1)            # exact rational value of code point 1
AllPositiveFiniteValuesOf(fmt)  # every positive finite value, in order
```

## Package structure

| Source file           | Purpose                                                             |
|:----------------------|:--------------------------------------------------------------------|
| `closedrationals.jl`  | `ClosedRational` exact rational type with NaN/±Inf support          |
| `types.jl`            | `Format` type, signedness/domain traits, accessors                  |
| `support.jl`          | `twopow`, `dyadic_twopow` -- exact power-of-two helpers            |
| `utils.jl`            | `sign_half_offset`, `sign_reduce`, `significand_scale`, `EpsilonOf` |
| `counts.jl`           | Closed-form counts of every value class                             |
| `codepoints.jl`       | Code-point locations for boundaries and specials                    |
| `floatvalues.jl`      | Boundary *values* (subnormal min/max, normal min/max)               |
| `values.jl`           | `FiniteValueOf` decoder, `AllFiniteValuesOf`, ordinal accessors     |

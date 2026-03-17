# Validation

Float3109s.jl includes test suites designed to verify the internal consistency
of the model across all four format variants and a wide range of `(K, P)`
combinations.

## Test files

| File | What it tests |
|:-----|:--------------|
| `test_counts.jl` | Partition identities, hand-computed values, cross-validation against exhaustive enumeration |
| `test_codepoints.jl` | Boundary locations, ordering, contiguity, sign operations, classification predicates |
| `test_floatvalues.jl` | Boundary values, cross-checks against `FiniteValueOf`, even spacing of subnormals |
| `test_values.jl` | `FiniteValueOf` decoding, monotonicity, negative symmetry, ULP doubling, enumeration lengths |

## Test strategy

### Format coverage

Every test section sweeps across multiple format variants:

- **All four trait combinations**: unsigned/signed crossed with finite/extended
- **Range of `(K, P)`**: typically `K = 2..10`, `P = 1..(K-1)`
- **Corner cases**: `P = 1` (no subnormals), `P = K` or `P = K-1` (minimal
  exponent), `K = 2` (minimal format)

### Partition invariants

The count tests verify that the partition identities documented in the
[Counts](@ref) section hold exactly for every valid `(K, P, S, D)` combination.
These algebraic checks catch any inconsistency between counting functions.

### Cross-validation against exhaustive enumeration

For small formats (`K <= 8`), the tests walk every code point, classify it
by category, and compare the tallies against the closed-form count functions.
This independently confirms that the formulas match the code-point layout.

### Hand-computed golden values

Several test sections verify specific small formats against manually derived
values.  These serve as regression anchors:

- `Format{UnsignedFormat, FiniteFormat}(3, 2)` — 8 code points, full value ladder
  `[0, 1/4, 1/2, 3/4, 1, 3/2, 2]`
- `Format{UnsignedFormat, FiniteFormat}(4, 3)` — 16 code points, 3 subnormals
- `Format{SignedFormat, ExtendedFormat}(4, 2)` — 16 code points, full positive and
  negative layout with ±Inf
- `Format{SignedFormat, FiniteFormat}(4, 2)` — 16 code points, symmetric signed
  layout without infinities

### Structural properties

Beyond partition checks, the tests verify:

- **Strict monotonicity** of positive-half values with increasing code point
- **Negative symmetry**: `FiniteValueOfOrdinalNeg(fmt, i) == -FiniteValueOfOrdinalPos(fmt, i)`
- **ULP doubling** between consecutive normal binades
- **Contiguity**: no gaps between subnormal and normal code-point ranges
- **Involution**: `cp_changesign(cp_changesign(cp)) == cp`
- **Dyadic rationality**: all decoded values have power-of-2 denominators

## Running the tests

```julia
using Pkg
Pkg.test("Float3109s")
```

Or run individual test files:

```julia
include("test/test_counts.jl")
include("test/test_codepoints.jl")
include("test/test_floatvalues.jl")
include("test/test_values.jl")
```

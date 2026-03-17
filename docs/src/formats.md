# Formats

## The `Format` type

Every format in Float3109s.jl is an instance of `Format{S,D}`, parameterized
by two trait types that encode signedness and domain at the type level.

```julia
struct Format{S<:Signedness, D<:Domain} <: FloatFormat
    K::Int   # bitwidth
    P::Int   # precision (including implicit leading bit)
end
```

### Signedness traits

| Trait | Meaning |
|:------|:--------|
| `is_unsigned` | No sign bit. All non-zero numerical values are positive. |
| `is_signed` | One sign bit. Positive and negative halves are symmetric. |

### Domain traits

| Trait | Meaning |
|:------|:--------|
| `is_finite` | No infinities. Every non-NaN code point is a finite number. |
| `is_extended` | Has +Inf (and -Inf for signed formats). |

### Constructing formats

Combine the traits as type parameters:

```julia
using Float3109s

# All four variants of an 8-bit, precision-4 format:
Format{UnsignedFormat, FiniteFormat}(8, 4)
Format{UnsignedFormat, ExtendedFormat}(8, 4)
Format{SignedFormat,   FiniteFormat}(8, 4)
Format{SignedFormat,   ExtendedFormat}(8, 4)
```

The constructor enforces `W >= 1` (the exponent field must have at least one
bit).  This means:

- Unsigned formats require `K >= P`
- Signed formats require `K >= P + 1`

### Querying traits

```julia
fmt = Format{SignedFormat, ExtendedFormat}(8, 4)

is_signed(fmt)     # true
is_unsigned(fmt)    # false
is_finite(fmt)      # false
is_extended(fmt)    # true
```

## Format accessors

These functions extract derived structural quantities from a format.

| Function | Returns | Formula |
|:---------|:--------|:--------|
| `BitwidthOf(fmt)` | Total bits `K` | `fmt.K` |
| `PrecisionOf(fmt)` | Precision `P` | `fmt.P` |
| `TrailingBitsOf(fmt)` | Trailing significand bits | `P - 1` |
| `SignBitsOf(fmt)` | Number of sign bits | `0` (unsigned) or `1` (signed) |
| `NonSignificantBitsOf(fmt)` | Non-significand bits | `K - P` |
| `ExponentBitsOf(fmt)` | Exponent field width `W` | `K - P + 1` (unsigned) or `K - P` (signed) |
| `ExponentBiasOf(fmt)` | Exponent bias `B` | `2^(K-P)` (unsigned) or `2^(K-P-1)` (signed) |
| `EpsilonOf(fmt)` | Machine epsilon | `2^(1 - P)` |

### Example

```julia
fmt = Format{SignedFormat, ExtendedFormat}(8, 4)

BitwidthOf(fmt)       # 8
PrecisionOf(fmt)      # 4
ExponentBitsOf(fmt)   # 4   (= 8 - 4)
ExponentBiasOf(fmt)   # 8   (= 2^3)
EpsilonOf(fmt)        # 0.125
```

## Internal helpers

These functions are not exported but are used throughout the package.

### `twopow(n)`

Computes `2^n` as `UInt128` for non-negative `n`, or as `BigFloat` for
negative `n`.  Used for code-point and structural calculations.

### `dyadic_twopow(n)`

Computes `2^n` as an exact `Rational{BigInt}`.  Works for both positive and
negative `n`.  Used wherever exact rational arithmetic is needed.

### `significand_scale(fmt)`

Returns `2^(P-1)`, the number of trailing-significand states per binade.
This is the boundary between subnormal and normal code points.

### `sign_half_offset(fmt)`

Returns `2^(K-1)` for signed formats and `0` for unsigned formats.
This is the code-point offset that separates the positive and negative halves.

### `sign_reduce(fmt, cp)`

Decomposes a code point into its sign and positive-half code point:

```julia
(s, cp_abs) = sign_reduce(fmt, cp)
# s = 0  →  positive (or zero)
# s = 1  →  negative
# cp_abs →  the corresponding positive-half code point
```

Throws `ArgumentError` if `cp` is the NaN midpoint in a signed format.

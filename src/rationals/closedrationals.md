# ExtendedRational

`ExtendedRational` is the exact rational number type used throughout Float3109s.jl.
It extends Julia's rational arithmetic to include the three IEEE-style special
values NaN, +Inf, and -Inf, forming a *closed* number system where every
operation returns a `ExtendedRational`.

## Canonical representation

Values are stored as `(num::BigInt, den::BigInt)` in reduced form:

| Value    | `num` | `den` | Display      |
|:---------|------:|------:|:-------------|
| Finite   | any   | > 0   | `num//den`   |
| Zero     | 0     | 1     | `0//1`       |
| +Inf     | 1     | 0     | `Inf`        |
| -Inf     | -1    | 0     | `-Inf`       |
| NaN      | 0     | 0     | `NaN`        |

The inner constructor enforces these invariants automatically:
- Finite values are GCD-reduced with `den > 0`.
- Any `(0, 0)` input becomes NaN.
- Any `(positive, 0)` input becomes +Inf.
- Any `(negative, 0)` input becomes -Inf.

## Construction

```julia
ExtendedRational(3, 4)          # 3//4
ExtendedRational(6, 8)          # 3//4 (auto-reduced)
ExtendedRational(-1, 0)         # -Inf
ExtendedRational(0, 0)          # NaN

ExtendedRational(1.5)           # 3//2 (from Float64)
ExtendedRational(Float32(0.1))  # exact Float32 rational
ExtendedRational(Inf)           # Inf
ExtendedRational(NaN)           # NaN

ExtendedRational(3//4)          # from Julia Rational

NaN(ExtendedRational)           # NaN
Inf(ExtendedRational)           # +Inf
NegInf(ExtendedRational)        # -Inf
```

Integer arguments of any `Integer` subtype are accepted and promoted to `BigInt`.

## Predicates

```julia
isnan(q)      # true for NaN (0//0)
isfinite(q)   # true for finite values (den != 0)
iszero(q)     # true for zero (0//1)
isone(q)      # true for one (1//1)
signbit(q)    # true for negative (num < 0)
```

## Arithmetic

All operations follow IEEE-style special-value semantics.

### Basic operations

```julia
+a, -a        # unary plus/minus
a + b         # addition
a - b         # subtraction
a * b         # multiplication
abs(a)        # absolute value
```

### Special-value rules

| Operation       | Result |
|:----------------|:-------|
| NaN op anything | NaN    |
| Inf + Inf       | Inf    |
| Inf + (-Inf)    | NaN    |
| Inf * 0         | NaN    |
| Inf * finite    | ±Inf   |

### Integer division

```julia
div(a, b)     # truncated toward zero
fld(a, b)     # floored toward -Inf
cld(a, b)     # ceiled toward +Inf
```

For finite `a` and infinite `b`:
- `div` returns 0
- `fld` returns 0 or -1 depending on signs
- `cld` returns 0 or 1 depending on signs

### Remainders

```julia
rem(a, b)     # a - div(a, b) * b  (sign follows a)
mod(a, b)     # a - fld(a, b) * b  (sign follows b)
```

### 1-based division

```julia
fld1(a, b)    # fld(a - 1, b) + 1
mod1(a, b)    # mod(a - 1, b) + 1   (result in [1, b])
```

### Tuple forms

```julia
divrem(a, b)  # (div(a, b), rem(a, b))
fldmod(a, b)  # (fld(a, b), mod(a, b))
fldmod1(a, b) # (fld1(a, b), mod1(a, b))
```

### GCD and LCM

For rational numbers `a/b` and `c/d`:

```julia
gcd(a, b)     # gcd(num_a, num_b) // lcm(den_a, den_b)
lcm(a, b)     # lcm(num_a, num_b) // gcd(den_a, den_b)
```

Special cases:
- `gcd(Inf, x) = abs(x)`, `gcd(Inf, Inf) = Inf`
- `lcm(0, x) = 0`, `lcm(Inf, x) = Inf`

## Display

```julia
julia> q = ExtendedRational(3, 4)
ExtendedRational(3//4)          # REPL (text/plain)

julia> print(q)
3//4                          # compact (IO)

julia> string(q)
"3//4"

julia> string(NaN(ExtendedRational))
"NaN"
```

## Accessors

```julia
numerator(q)    # BigInt numerator
denominator(q)  # BigInt denominator
Tuple(q)        # (num, den) as a tuple
```

## Why not `Rational{BigInt}`?

Julia's built-in `Rational{BigInt}` has no representation for NaN or ±Inf.
Since P3109 formats include these special values, `ExtendedRational` provides
a single unified type that can represent *every* code-point value exactly,
eliminating the need for separate marker types or `Union` dispatch.

# Concepts

This page describes the mathematical model behind Float3109s.jl.

## Format parameters

Every P3109 format is defined by four quantities:

| Symbol | Name | Meaning |
|:-------|:-----|:--------|
| `K`    | Bitwidth  | Total number of bits per code point |
| `P`    | Precision | Number of significand bits, including the implicit leading bit |
| `S`    | Signedness | `is_unsigned` or `is_signed` |
| `D`    | Domain     | `is_finite` (no infinities) or `is_extended` (has ±Inf) |

These are encoded as Julia type parameters on [`Format{S,D}`](@ref Format):

```julia
Format{SignedFormat, ExtendedFormat}(8, 4)   # K=8, P=4, signed, with infinities
Format{UnsignedFormat, FiniteFormat}(8, 4)   # K=8, P=4, unsigned, finite-only
```

## Derived quantities

From `K`, `P`, and the signedness, the exponent field width and bias follow:

| Quantity | Unsigned formula | Signed formula |
|:---------|:-----------------|:---------------|
| Sign bits `Σ` | `0` | `1` |
| Non-significant bits | `K − P` | `K − P` |
| Exponent bits `W`    | `K − P + 1` | `K − P` |
| Exponent bias `B`    | `2^(K − P)` | `2^(K − P − 1)` |

The structural constraint is `W ≥ 1`, which means:
- Unsigned: `K ≥ P`
- Signed: `K ≥ P + 1`

## Code-point space layout

The `2^K` code points are partitioned into distinct regions.  The layout
depends on signedness and domain.

### Unsigned format

```
0          1  ...  m−1        m  ...  cp_max−δ  ...  cp_max
│          │        │          │        │               │
zero    subnormals       normals      [+Inf]          NaN
```

- Code point `0` is zero.
- Code points `1` to `m − 1` (where `m = 2^(P−1)`) are positive subnormals.
- Code points `m` to `cp_max − δ` are positive normals, where `δ` accounts
  for +Inf (if extended) and NaN.
- The last code point `cp_max = 2^K − 1` is NaN.
- In extended formats, the code point just before NaN is +Inf.

### Signed format

```
Positive half:
0       1 ... m−1    m ... NaN−δ   [+Inf]   NaN

Negative half:
NaN+1  ...  NaN+m−1  ...  cp_max−δ  ...  [−Inf]  cp_max
│            │              │                │       │
neg subs   neg normals   neg normals      [−Inf]  (−Inf or last neg normal)
```

- The midpoint `2^(K−1)` is NaN.
- The lower half `[0, NaN)` holds zero and positive values.
- The upper half `(NaN, 2^K − 1]` holds negative values, with magnitude
  increasing with code point.
- In extended formats, `NaN − 1` is +Inf and `2^K − 1` is −Inf.

## Value classes

Every code point belongs to exactly one class:

| Class | Description |
|:------|:------------|
| **Zero** | The single code point `0`, representing the value `0`. |
| **Subnormal** | Values with magnitude below the normal range. Exist only when `P ≥ 2`. Evenly spaced with step size equal to the smallest representable positive value. |
| **Normal** | Values within the normal binade range. Organized into `2B − 1` binades, each containing `2^(P−1)` evenly spaced values. The ULP doubles between consecutive binades. |
| **Infinity** | `+Inf` and/or `−Inf`, present only in extended-domain formats. |
| **NaN** | Exactly one NaN per format. |

### Prenormals

The package uses the term **prenormal** to refer to zero and subnormals
collectively — the `2^(P−1)` code points in the range `[0, m)` of the
positive half.  For `P = 1` there are no subnormals; the only prenormal
is zero itself.

### Partition identities

The counting functions satisfy a set of exact partition identities.  These
are useful both as sanity checks and as a specification of the model:

```
nValuesOf        = nNumericalValuesOf + nNaNsOf
nNumericalValuesOf = nFiniteValuesOf + nInfsOf
nFiniteValuesOf  = nPrenormalsOf + nNormalsOf
nPrenormalsOf    = nZerosOf + nSubnormalsOf
nNonZeroFiniteValuesOf = nPosFiniteValuesOf + nNegFiniteValuesOf
nInfsOf          = nPosInfsOf + nNegInfsOf
```

## Binade structure

Normal values are organized into **binades** — groups of `2^(P−1)` consecutive
code points that share the same exponent.  Within each binade, values are
evenly spaced.  The step size (ULP) doubles from one binade to the next:

```
Binade k:   values are {2^k + t · 2^(k+1−P)}  for t = 0, 1, ..., 2^(P−1) − 1
ULP in binade k:  Δ_k = 2^(k+1−P)
```

The total number of binades is `2B − 1`, where `B` is the exponent bias.

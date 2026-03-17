# Code Points

The functions in `codepoints.jl` return the exact code-point locations of
boundary values and special values.  Code points are returned as the
narrowest unsigned integer type that fits `K` bits (e.g., `UInt8` for
`K <= 8`, `UInt16` for `K <= 16`, etc.).

## Range boundaries

```julia
cp_min(fmt)   # always 0 (the smallest code point)
cp_max(fmt)   # 2^K - 1 (the largest code point)
```

## Special code points

### NaN

```julia
cp_nan(fmt)
```

- **Unsigned**: NaN is the last code point, `2^K - 1`.
- **Signed**: NaN is the midpoint, `2^(K-1)`.

Every format has exactly one NaN.

### Infinities

```julia
cp_inf(fmt)      # +Inf code point, or `nothing` for finite formats
cp_neginf(fmt)   # -Inf code point, or `nothing` if absent
```

- **Finite-domain formats**: both return `nothing`.
- **Unsigned extended**: `cp_inf` returns `cp_nan - 1`; `cp_neginf` returns
  `nothing` (no negative infinity in unsigned formats).
- **Signed extended**: `cp_inf` returns `cp_nan - 1`; `cp_neginf` returns
  `cp_max` (the last code point).

### Zero

```julia
cp_zero(fmt)   # always 0
```

## Subnormal boundaries

Subnormals exist only when `P >= 2`.  When `P = 1`, these functions return
`nothing`.

### Positive subnormals

```julia
cp_pos_subnormal_min(fmt)   # 1   (or nothing if P=1)
cp_pos_subnormal_max(fmt)   # 2^(P-1) - 1   (or nothing if P=1)
```

The positive subnormals occupy a contiguous range from code point `1` to
`2^(P-1) - 1`, immediately below the first normal.

### Negative subnormals (signed only)

```julia
cp_neg_subnormal_min(fmt)   # cp_nan + 1   (or nothing)
cp_neg_subnormal_max(fmt)   # (or nothing for unsigned or P=1)
```

Unsigned formats always return `nothing`.  For signed formats with `P >= 2`,
the negative subnormals sit immediately above the NaN midpoint, mirroring
the positive subnormals.

## Normal boundaries

### Positive normals

```julia
cp_pos_normal_min(fmt)   # 2^(P-1)  (always present)
cp_pos_normal_max(fmt)   # depends on signedness and domain
```

The value of `cp_pos_normal_max` differs by format variant:

| Variant | `cp_pos_normal_max` |
|:--------|:--------------------|
| unsigned finite | `2^K - 2` |
| unsigned extended | `2^K - 3` |
| signed finite | `cp_nan - 1` |
| signed extended | `cp_nan - 2` |

### Negative normals (signed only)

```julia
cp_neg_normal_min(fmt)   # smallest-magnitude negative normal (or nothing)
cp_neg_normal_max(fmt)   # largest-magnitude negative normal (or nothing)
```

Unsigned formats return `nothing`.  For signed formats:

- `cp_neg_normal_min` is immediately after the last negative subnormal
  (or after NaN if `P = 1`).
- `cp_neg_normal_max` is `cp_max` for finite formats, or `cp_max - 1`
  for extended formats (since `cp_max` is -Inf).

## Contiguity

The code-point ranges are contiguous with no gaps:

```
cp_pos_subnormal_max + 1 == cp_pos_normal_min
cp_neg_subnormal_max + 1 == cp_neg_normal_min   (signed, P >= 2)
```

## Ordinal accessors

These functions return the code point of the `i`-th value (1-based) within
a class:

```julia
cp_ordinal_ith_pos_subnormal(fmt, i)   # code point of the i-th positive subnormal
cp_ordinal_ith_pos_normal(fmt, i)      # code point of the i-th positive normal
```

Both throw `ArgumentError` if `i` is out of range.  The unsafe variants
skip bounds checking:

```julia
unsafe_cp_ordinal_ith_pos_subnormal(fmt, i)
unsafe_cp_ordinal_ith_pos_normal(fmt, i)
```

## Classification predicates

```julia
cp_is_positive(fmt, cp)      # true if 0 < cp < cp_nan
cp_is_nonnegative(fmt, cp)   # true if 0 <= cp < cp_nan
cp_is_negative(fmt, cp)      # true if cp > cp_nan (signed only; always false for unsigned)
```

```julia
cp_positive_max(fmt)   # cp_nan - 1 (largest positive code point, may be +Inf)
```

## Sign operations (signed formats)

```julia
pos_cp_to_neg_cp(fmt, cp_pos)   # map a positive code point to its negative counterpart
neg_cp_to_pos_cp(fmt, cp_neg)   # map a negative code point to its positive counterpart
cp_changesign(fmt, cp)          # flip the sign; zero and NaN are fixed points
```

`cp_changesign` is an involution: applying it twice returns the original code
point.

## Example

```julia
fmt = Format{SignedFormat, ExtendedFormat}(4, 2)

cp_zero(fmt)              # 0x00
cp_pos_subnormal_min(fmt) # 0x01
cp_pos_subnormal_max(fmt) # 0x01   (only one subnormal when P=2)
cp_pos_normal_min(fmt)    # 0x02
cp_pos_normal_max(fmt)    # 0x06
cp_inf(fmt)               # 0x07
cp_nan(fmt)               # 0x08
cp_neg_subnormal_min(fmt) # 0x09
cp_neg_normal_max(fmt)    # 0x0e
cp_neginf(fmt)            # 0x0f
```

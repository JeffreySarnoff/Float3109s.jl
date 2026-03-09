# Values

The functions in `floatvalues.jl` and `values.jl` decode code points into
their exact rational values as `ClosedRational` values (see
[ClosedRational](@ref) for the type definition).

## Decoding a single code point

```julia
FiniteValueOf(fmt, cp) -> ClosedRational
```

Returns the exact value of code point `cp` in format `fmt`.

- `cp = 0` returns `0 // 1`.
- NaN, +Inf, and -Inf code points are not accepted (use `ValueOf` for those).
- Out-of-range code points throw `ArgumentError`.
- All finite values are returned as exact dyadic rationals (denominators are
  powers of 2).

### How decoding works

After sign reduction, the positive-half code point `cp_abs` is decoded as:

**Subnormal** (`cp_abs < 2^(P-1)`):

```
value = cp_abs * 2^(2 - P - B)
```

Subnormals are evenly spaced with a constant step size (the subnormal ULP).

**Normal** (`cp_abs >= 2^(P-1)`):

```
e = cp_abs ÷ 2^(P-1)       # stored exponent index
t = cp_abs mod 2^(P-1)     # trailing significand
k = e - B                  # unbiased exponent
value = 2^k + t * 2^(k+1-P)
```

For signed formats with `s = 1`, the final value is negated.

## Boundary values

These functions return the exact value at specific boundaries, or `nothing`
when the boundary does not exist (e.g., subnormal boundaries when `P = 1`).

### Subnormal boundaries

```julia
val_pos_subnormal_min(fmt)   # smallest positive subnormal value (nothing if P=1)
val_pos_subnormal_max(fmt)   # largest positive subnormal value (nothing if P=1)
```

### Normal boundaries

```julia
val_pos_normal_min(fmt)      # smallest positive normal = 2^(1 - B)
val_pos_normal_max(fmt)      # largest positive finite normal
```

`val_pos_normal_max` has different formulas depending on the format variant,
accounting for the code points consumed by +Inf and -Inf in extended formats.

### Ordinal subnormal values

```julia
val_ordinal_ith_pos_subnormal(fmt, i)      # value of the i-th positive subnormal (1-based)
val_cardinal_ith_pos_subnormal(fmt, cardinal)  # value of the cardinal-th subnormal (0-based)
```

Since subnormals are evenly spaced:

```
val_ordinal_ith_pos_subnormal(fmt, i) = i * 2^(2 - P - B)
```

## Enumerating all values

### All finite values

```julia
AllFiniteValuesOf(fmt) -> Vector{ClosedRational}
```

Returns every finite value in code-point order.  Zero is included; NaN and
infinities are excluded.  For signed formats, the order is:
`zero, positive values, negative values` (following code-point order).

### Positive finite values only

```julia
AllPositiveFiniteValuesOf(fmt) -> Vector{ClosedRational}
```

Returns every positive finite value in code-point order.  Zero is excluded.

The length of the returned vector matches `nPosFiniteValuesOf(fmt)`.

## Ordinal value accessors

These provide 1-based ordinal access to positive and negative values, ordered
by increasing magnitude:

```julia
FiniteValueOfOrdinalPos(fmt, i)   # i-th positive finite value (i=1 is smallest)
FiniteValueOfOrdinalNeg(fmt, i)   # i-th negative finite value (i=1 is closest to zero)
```

- `FiniteValueOfOrdinalPos(fmt, 1)` returns the smallest positive value.
- `FiniteValueOfOrdinalPos(fmt, nPosFiniteValuesOf(fmt))` returns the largest
  positive finite value.
- `FiniteValueOfOrdinalNeg` is only valid for signed formats; it throws
  `ArgumentError` for unsigned formats.

## Key properties

### Strict monotonicity

Positive-half values are strictly increasing with code point:

```
FiniteValueOf(fmt, cp) < FiniteValueOf(fmt, cp + 1)
```

for all consecutive positive finite code points.

### Negative symmetry (signed formats)

Every positive value has a negative counterpart of equal magnitude:

```
FiniteValueOfOrdinalNeg(fmt, i) == -FiniteValueOfOrdinalPos(fmt, i)
```

### ULP doubling

The step size doubles between consecutive normal binades:

```julia
fmt = Format{is_unsigned, is_finite}(8, 4)

v1 = FiniteValueOf(fmt, cp_pos_normal_min(fmt) + 1) - FiniteValueOf(fmt, cp_pos_normal_min(fmt))
# ULP in first binade

v2 = FiniteValueOf(fmt, cp_pos_normal_min(fmt) + 9) - FiniteValueOf(fmt, cp_pos_normal_min(fmt) + 8)
# ULP in second binade (8 values per binade when P=4)

v2 == 2 * v1   # true
```

## Example: complete value ladder

```julia
fmt = Format{is_unsigned, is_finite}(3, 2)

AllFiniteValuesOf(fmt)
# 7-element Vector{ClosedRational}:
#  0//1  1//4  1//2  3//4  1//1  3//2  2//1
```

Here `K=3`, `P=2`, so:
- 1 subnormal at `1/4`
- 3 binades of 2 normals each: `{1/2, 3/4}`, `{1, 3/2}`, `{2}`
  (the last binade has only one value because `cp = 6` is the last before NaN)

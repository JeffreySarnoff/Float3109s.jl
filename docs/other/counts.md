# Counts

The functions in `counts.jl` return exact, closed-form counts for every class
of value in a format.  All return non-negative integers.

## Total partition

These functions decompose the full `2^K` code-point space:

```julia
nValuesOf(fmt)                 # 2^K  (all code points)
nNaNsOf(fmt)                   # always 1
nNumericalValuesOf(fmt)        # 2^K - 1  (everything except NaN)
nZerosOf(fmt)                  # always 1
nNonZeroNumericalValuesOf(fmt) # numerical values minus zero
```

## Infinity counts

```julia
nInfsOf(fmt)      # 0 (finite), 1 (unsigned extended), 2 (signed extended)
nPosInfsOf(fmt)   # 0 or 1
nNegInfsOf(fmt)   # 0 or 1
```

| Format variant | `nInfsOf` | `nPosInfsOf` | `nNegInfsOf` |
|:---------------|:---------:|:------------:|:------------:|
| any finite     | 0         | 0            | 0            |
| unsigned extended | 1      | 1            | 0            |
| signed extended   | 2      | 1            | 1            |

## Finite value counts

```julia
nFiniteValuesOf(fmt)          # numerical values minus infinities
nNonZeroFiniteValuesOf(fmt)   # finite values minus zero
nNonNegFiniteValuesOf(fmt)    # positive finite + zero
```

## Positive / negative split

For signed formats, non-zero numerical values split evenly between positive
and negative:

```julia
nPosValuesOf(fmt)    # positive non-zero values (including +Inf if present)
nNegValuesOf(fmt)    # negative non-zero values (including -Inf if present)
nNonNegValuesOf(fmt) # positive values + zero
```

For unsigned formats, `nNegValuesOf` is always `0` and all non-zero numerical
values are positive.

The finite versions subtract out infinities:

```julia
nPosFiniteValuesOf(fmt)    # nPosValuesOf - nPosInfsOf
nNegFiniteValuesOf(fmt)    # nNegValuesOf - nNegInfsOf
nNonNegFiniteValuesOf(fmt) # nPosFiniteValuesOf + nZerosOf
```

## Subnormal and prenormal counts

**Prenormals** are the `2^(P-1)` code points in the positive-half range
`[0, m)` — that is, zero together with all subnormals:

```julia
nNonNegPrenormalsOf(fmt)   # 2^(P-1)   (or 1 when P=1)
nPrenormalsOf(fmt)         # unsigned: same as nNonNegPrenormals
                           # signed:   nNonNegPrenormals + (nNonNegPrenormals - 1)
```

Subnormals are prenormals minus zero:

```julia
nPosSubnormalsOf(fmt)   # 2^(P-1) - 1   (or 0 when P=1)
nNegSubnormalsOf(fmt)   # 0 (unsigned), equals nPosSubnormalsOf (signed)
nSubnormalsOf(fmt)      # nPosSubnormalsOf + nNegSubnormalsOf
```

When `P = 1`, there are no subnormals — the only prenormal is zero itself.

## Normal counts

```julia
nNormalsOf(fmt)      # nFiniteValuesOf - nPrenormalsOf
nPosNormalsOf(fmt)   # all normals (unsigned), or half (signed)
nNegNormalsOf(fmt)   # 0 (unsigned), equals nPosNormalsOf (signed)
```

## Binade count

```julia
nBinadesOf(fmt)   # 2B - 1, where B = ExponentBiasOf(fmt)
```

Each binade contains exactly `2^(P-1)` normal code points.

## Partition identities

The following equalities always hold and serve as self-consistency checks:

```
nValuesOf           == nNumericalValuesOf + nNaNsOf
nNumericalValuesOf  == nZerosOf + nNonZeroNumericalValuesOf
nNonZeroNumericalValuesOf == nNonZeroFiniteValuesOf + nInfsOf
nFiniteValuesOf     == nPrenormalsOf + nNormalsOf
nPrenormalsOf       == nZerosOf + nSubnormalsOf
nSubnormalsOf       == nPosSubnormalsOf + nNegSubnormalsOf
nNormalsOf          == nPosNormalsOf + nNegNormalsOf
nInfsOf             == nPosInfsOf + nNegInfsOf
nNonZeroFiniteValuesOf == nPosFiniteValuesOf + nNegFiniteValuesOf
nPosFiniteValuesOf  == nPosSubnormalsOf + nPosNormalsOf
nNegFiniteValuesOf  == nNegSubnormalsOf + nNegNormalsOf
```

## Example

```julia
fmt = Format{is_signed, is_extended}(8, 4)

nValuesOf(fmt)           # 256
nFiniteValuesOf(fmt)     # 253   (256 - 1 NaN - 2 Infs)
nPosSubnormalsOf(fmt)    # 7     (2^3 - 1)
nPosNormalsOf(fmt)       # 119
nBinadesOf(fmt)          # 15    (2*8 - 1)
```

# API Reference

Complete reference for all exported types and functions.

## Types and Traits

```@docs
Qx64
Format
Signedness
Domain
is_signed
is_unsigned
is_finite
is_extended
```

## Format Accessors

```@docs
BitwidthOf
PrecisionOf
TrailingBitsOf
SignBitsOf
NonSignificantBitsOf
ExponentBitsOf
ExponentBiasOf
EpsilonOf
```

## Counting Functions

### Total and Numerical

```@docs
nValuesOf
nNumericalValuesOf
nNonZeroNumericalValuesOf
nNaNsOf
nZerosOf
```

### Infinities

```@docs
nInfsOf
nPosInfsOf
nNegInfsOf
```

### Finite Values

```@docs
nFiniteValuesOf
nNonZeroFiniteValuesOf
nNonNegFiniteValuesOf
nPosFiniteValuesOf
nNegFiniteValuesOf
```

### Positive / Negative Split

```@docs
nPosValuesOf
nNegValuesOf
nNonNegValuesOf
```

### Prenormals and Subnormals

```@docs
nNonNegPrenormalsOf
nPrenormalsOf
nPosSubnormalsOf
nNegSubnormalsOf
nSubnormalsOf
```

### Normals and Binades

```@docs
nNormalsOf
nPosNormalsOf
nNegNormalsOf
nBinadesOf
```

## Code Point Functions

### Range and Specials

```@docs
cp_min
cp_max
cp_zero
cp_nan
cp_inf
cp_neginf
```

### Subnormal Boundaries

```@docs
cp_pos_subnormal_min
cp_pos_subnormal_max
cp_neg_subnormal_min
cp_neg_subnormal_max
```

### Normal Boundaries

```@docs
cp_pos_normal_min
cp_pos_normal_max
cp_neg_normal_min
cp_neg_normal_max
```

### Ordinal Accessors

```@docs
cp_ordinal_ith_pos_subnormal
cp_ordinal_ith_pos_normal
```

### Classification and Sign

```@docs
cp_positive_max
cp_is_positive
cp_is_nonnegative
cp_is_negative
cp_changesign
```

## Value Functions

### Boundary Values

```@docs
val_pos_subnormal_min
val_pos_subnormal_max
val_pos_normal_min
val_pos_normal_max
```

### Decoding

```@docs
ValueOf
FiniteValueOf
```

### Enumeration

```@docs
AllFiniteValuesOf
AllPositiveFiniteValuesOf
```

### Value Ordinal Accessors

```@docs
FiniteValueOfOrdinalPos
FiniteValueOfOrdinalNeg
```

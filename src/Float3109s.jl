"""
Float3109s.jl

Reference implementation of an analytic floating-point geometry model for
parameterized P3109-style formats.

The implementation emphasizes:

• exact arithmetic using `Rational{BigInt}`
• explicit structure for format parameters
• direct correspondence between formulas and code
• suitability as a correctness reference and test oracle
"""
module Float3109s

export Format, Signedness, Domain,
    is_signed, is_unsigned, is_finite, is_extended,
    ExtendedRational, NegInf,
    BitwidthOf, PrecisionOf, TrailingBitsOf, SignBitsOf,
    NonSignificantBitsOf, ExponentBitsOf, ExponentBiasOf,
    EpsilonOf,
    nBinadesOf,
    nValuesOf, nNumericalValuesOf, nNonZeroNumericalValuesOf,
    nNonNegValuesOf, nPosValuesOf, nNegValuesOf,
    nFiniteValuesOf, nNonZeroFiniteValuesOf,
    nNonNegFiniteValuesOf, nPosFiniteValuesOf, nNegFiniteValuesOf,
    nPrenormalsOf, nNonNegPrenormalsOf,
    nSubnormalsOf, nPosSubnormalsOf, nNegSubnormalsOf,
    nNormalsOf, nPosNormalsOf, nNegNormalsOf,
    nZerosOf, nInfsOf, nPosInfsOf, nNegInfsOf, nNaNsOf,
    cp_min, cp_max,
    cp_zero, cp_nan, cp_inf, cp_neginf,
    cp_pos_subnormal_min, cp_pos_subnormal_max, cp_neg_subnormal_min, cp_neg_subnormal_max,
    cp_pos_normal_min, cp_pos_normal_max, cp_neg_normal_min, cp_neg_normal_max,
    cp_ordinal_ith_pos_subnormal, cp_ordinal_ith_pos_normal,
    cp_positive_max, cp_is_positive, cp_is_nonnegative, cp_is_negative,
    cp_changesign,
    val_pos_subnormal_min, val_pos_subnormal_max,
    val_neg_subnormal_min, val_neg_subnormal_max,
    val_pos_normal_min, val_pos_normal_max,
    val_neg_normal_min, val_neg_normal_max,
    FiniteValueOf, FiniteValueOfOrdinalPos, FiniteValueOfOrdinalNeg,
    AllFiniteValuesOf, AllPositiveFiniteValuesOf

setprecision(BigFloat, 1024)

include("ExtendedRationals.jl")
include("support.jl")
include("types.jl")
include("utils.jl")
include("counts.jl")
include("codepoints.jl")
include("floatvalues.jl")
include("values.jl")

end # module

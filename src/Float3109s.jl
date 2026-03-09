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

include("closedrationals.jl")
include("support.jl")
include("types.jl")
include("utils.jl")
include("counts.jl")
include("codepoints.jl")
include("floatvalues.jl")
include("values.jl")

end # module

#=
export AbstractP3109Format,
    P3109Format,
    DecodedValue,
    FiniteValue,
    PosInfValue,
    NegInfValue,
    NaNValue,
    twopow,
    dyadic_twopow,
    exponent_width,
    significand_scale,
    bias,
    sign_half_offset,
    cp_zero,
    cp_sub_min,
    cp_sub_max,
    cp_nor_min,
    cp_nor_max,
    nan_cp,
    posinf_cp,
    neginf_cp,
    xmin_normal,
    xmin_subnormal,
    classify_special,
    sign_reduce,
    decode_subnormal_abs,
    normal_components,
    decode_normal_abs,
    decode,
    floorlog2_dyadic,
    encode_finite,
    n_finite,
    n_positive_finite,
    n_subnormals,
    n_normals,
    n_binades,
    ulp_in_binade,
    ulp_of_abs_value,
    exact_mean_step_normal_binades,
    exact_second_moment_step_normal_binades,
    exact_variance_step_normal_binades,
    quantile_step_normal_binades,
    value_step_count,
    roundtrip_decode_encode,
    enumerate_positive_finite_values,
    enumerate_positive_steps,
    verify_monotone_positive_half,
    exhaustive_roundtrip_check,
    exhaustive_monotonicity_check,
    format_summary

"""
Abstract supertype for floating-point formats.
"""
abstract type AbstractP3109Format <: AbstractFloat end

"""
    P3109Format(K, P, Σ, Δ)

Defines a floating-point format parameterized by

• `K` — total number of code-point bits
• `P` — precision including the hidden leading bit
• `Σ` — signedness flag (`0` unsigned, `1` signed)
• `Δ` — extended-range flag (`0` finite only, `1` infinities present)

The constructor validates the basic structural invariant

`W = K - P + 1 - Σ ≥ 1`

so that the exponent field width is well-defined.
"""
struct P3109Format <: AbstractP3109Format
    K::Int
    P::Int
    Σ::Int
    Δ::Int
    function P3109Format(K::Integer, P::Integer, Σ::Integer, Δ::Integer)
        K = Int(K)
        P = Int(P)
        Σ = Int(Σ)
        Δ = Int(Δ)

        Σ in (0, 1) || throw(ArgumentError("Σ must be 0 or 1"))
        Δ in (0, 1) || throw(ArgumentError("Δ must be 0 or 1"))
        K >= 1 || throw(ArgumentError("K must be ≥ 1"))
        P >= 1 || throw(ArgumentError("P must be ≥ 1"))

        W = K - P + 1 - Σ
        W >= 1 || throw(ArgumentError("W = K - P + 1 - Σ must be ≥ 1"))

        new(K, P, Σ, Δ)
    end
end

"""
Superclass for decoded values.
"""
abstract type DecodedValue end

"""
Represents a finite decoded value using exact rational arithmetic.
"""
struct FiniteValue <: DecodedValue
    value::Rational{BigInt}
end

"""
Positive infinity marker.
"""
struct PosInfValue <: DecodedValue end

"""
Negative infinity marker.
"""
struct NegInfValue <: DecodedValue end

"""
NaN marker.
"""
struct NaNValue <: DecodedValue end

"""
    twopow(n)

Compute `twopown` as a `UInt128`.

This is intended for structural quantities and code-point calculations where
the exponent is known to fit in `UInt128` shifting semantics.
"""
twopow(n::Integer) = UInt128(1) << Int(n)

"""
    dyadic_twopow(n)

Return the exact dyadic rational `twopown`.
"""
function dyadic_twopow(n::Integer)
    return n >= 0 ? (BigInt(1) << Int(n)) // 1 : 1 // (BigInt(1) << Int(-n))
end

"""
    exponent_width(f)

Exponent width in bits.
"""
exponent_width(f::P3109Format) = f.K - f.P + 1 - f.Σ

"""
    significand_scale(f)

Return `twopow(P-1)`, the number of trailing-significand states per binade.
"""
significand_scale(f::P3109Format) = twopow(f.P - 1)

"""
    bias(f)

Exponent bias `twopow(W-1)`.
"""
bias(f::P3109Format) = twopow(exponent_width(f) - 1)

"""
    sign_half_offset(f)

Offset between the positive and negative halves of the code-point space.
Returns `0` for unsigned formats.
"""
sign_half_offset(f::P3109Format) = f.Σ == 1 ? twopow(f.K - 1) : UInt128(0)

"""
    cp_zero(f)

Code point of zero.
"""
cp_zero(::P3109Format) = UInt128(0)

"""
    cp_sub_min(f)

Smallest positive subnormal code point.
"""
cp_sub_min(::P3109Format) = UInt128(1)

"""
    cp_sub_max(f)

Largest positive subnormal code point.
"""
cp_sub_max(f::P3109Format) = significand_scale(f) - 1

"""
    cp_nor_min(f)

Smallest positive normal code point.
"""
cp_nor_min(f::P3109Format) = significand_scale(f)

"""
    cp_nor_max(f)

Largest positive finite normal code point.
"""
cp_nor_max(f::P3109Format) = twopow(f.K - f.Σ) - 1 - UInt128(f.Δ * (2 - f.Σ))

"""
    nan_cp(f)

Code point used for NaN.
"""
nan_cp(f::P3109Format) = f.Σ == 0 ? twopow(f.K) - 1 : twopow(f.K - 1)

"""
    posinf_cp(f)

Code point used for positive infinity, or `nothing` if the format has no infinity.
"""
function posinf_cp(f::P3109Format)
    f.Δ == 1 || return nothing
    return f.Σ == 0 ? twopow(f.K) - 2 : twopow(f.K - 1) - 1
end

"""
    neginf_cp(f)

Code point used for negative infinity, or `nothing` if the format has no negative infinity.
"""
function neginf_cp(f::P3109Format)
    (f.Δ == 1 && f.Σ == 1) || return nothing
    return twopow(f.K) - 1
end

"""
    xmin_normal(f)

Smallest positive normal value.
"""
xmin_normal(f::P3109Format) = dyadic_twopow(1 - Int(bias(f)))

"""
    xmin_subnormal(f)

Smallest positive subnormal value.
"""
xmin_subnormal(f::P3109Format) = dyadic_twopow(2 - f.P - Int(bias(f)))

"""
    classify_special(f, cp)

Classify `cp` as one of

- `:zero`
- `:nan`
- `:posinf`
- `:neginf`
- `:finite`
"""
function classify_special(f::P3109Format, cp::Integer)
    cp == 0 && return :zero
    cp == Int(nan_cp(f)) && return :nan
    p = posinf_cp(f)
    n = neginf_cp(f)
    p !== nothing && cp == Int(p) && return :posinf
    n !== nothing && cp == Int(n) && return :neginf
    return :finite
end

"""
    sign_reduce(f, cp)

Reduce a code point to positive-half form.

Returns a named tuple `(s, cp_abs)` where

- `s` is the sign indicator (`0` or `1`)
- `cp_abs` is the positive-half code point

For signed formats, the midpoint is reserved for NaN and is rejected here.
"""
function sign_reduce(f::P3109Format, cp::Integer)
    if f.Σ == 0
        return (; s=0, cp_abs=cp)
    end
    H = Int(sign_half_offset(f))
    cp < H && return (; s=0, cp_abs=cp)
    cp > H && return (; s=1, cp_abs=cp - H)
    throw(ArgumentError("midpoint code point is reserved for NaN in signed formats"))
end

"""
    decode_subnormal_abs(f, cp_abs)

Decode a positive-half subnormal code point exactly.
"""
function decode_subnormal_abs(f::P3109Format, cp_abs::Integer)
    1 <= cp_abs <= Int(cp_sub_max(f)) || throw(ArgumentError("cp_abs not subnormal"))
    return BigInt(cp_abs) * xmin_subnormal(f)
end

"""
    normal_components(f, cp_abs)

Return the Euclidean decomposition of a positive-half normal code point.

The result is a named tuple `(e, t, k)` with

- `e` — stored exponent index
- `t` — trailing significand index
- `k` — unbiased binade index
"""
function normal_components(f::P3109Format, cp_abs::Integer)
    m = Int(significand_scale(f))
    Int(cp_nor_min(f)) <= cp_abs <= Int(cp_nor_max(f)) ||
        throw(ArgumentError("cp_abs not normal"))
    e = fld(cp_abs, m)
    t = mod(cp_abs, m)
    k = e - Int(bias(f))
    return (; e, t, k)
end

"""
    decode_normal_abs(f, cp_abs)

Decode a positive-half normal code point exactly.
"""
function decode_normal_abs(f::P3109Format, cp_abs::Integer)
    comps = normal_components(f, cp_abs)
    return dyadic_twopow(comps.k) + BigInt(comps.t) * dyadic_twopow(comps.k + 1 - f.P)
end

"""
    decode(f, cp)

Decode code point `cp` into a `DecodedValue`.

Finite values are returned as `FiniteValue` with exact rational payload.
Special values are returned as marker objects.
"""
function decode(f::P3109Format, cp::Integer)::DecodedValue
    0 <= cp <= Int(twopow(f.K) - 1) || throw(ArgumentError("cp out of range"))
    cls = classify_special(f, cp)
    cls === :zero && return FiniteValue(0 // 1)
    cls === :nan && return NaNValue()
    cls === :posinf && return PosInfValue()
    cls === :neginf && return NegInfValue()

    red = sign_reduce(f, cp)
    cp_abs = red.cp_abs
    val =
        if 1 <= cp_abs <= Int(cp_sub_max(f))
            decode_subnormal_abs(f, cp_abs)
        elseif Int(cp_nor_min(f)) <= cp_abs <= Int(cp_nor_max(f))
            decode_normal_abs(f, cp_abs)
        else
            throw(ArgumentError("finite cp_abs lies outside subnormal/normal ranges"))
        end

    return FiniteValue(red.s == 0 ? val : -val)
end

"""
    floorlog2_dyadic(x)

Return `floor(log2(x))` for a positive dyadic rational `x`.
"""
function floorlog2_dyadic(x::Rational{BigInt})
    x > 0 || throw(ArgumentError("x must be positive"))
    num = numerator(x)
    den = denominator(x)
    return ndigits(num, base=2) - ndigits(den, base=2)
end

"""
    encode_finite(f, x)

Encode a finite rational value `x` that is already known to be exactly representable
in `f`.

This function performs representability checks and throws `ArgumentError` if the
input is not exact for the format.
"""
function encode_finite(f::P3109Format, x::Rational{BigInt})::UInt128
    x == 0 // 1 && return 0

    s = x < 0 ? 1 : 0
    f.Σ == 0 && s == 1 && throw(ArgumentError("negative value not representable in unsigned format"))

    a = abs(x)
    H = sign_half_offset(f)

    if a < xmin_normal(f)
        scale = BigInt(1) << Int(f.P + Int(bias(f)) - 2)
        cp_abs = a * scale
        denominator(cp_abs) == 1 || throw(ArgumentError("subnormal value not exactly representable"))
        return UInt128(numerator(cp_abs)) + UInt128(s) * H
    else
        k = floorlog2_dyadic(a)
        scaled = a * dyadic_twopow(f.P - 1 - k)
        denominator(scaled) == 1 || throw(ArgumentError("normal value not exactly representable"))
        cp_abs = (BigInt(k + Int(bias(f)) - 1) << Int(f.P - 1)) + numerator(scaled)
        return UInt128(cp_abs) + UInt128(s) * H
    end
end

"""
    n_finite(f)

Total number of finite numeric values in the format.
"""
function n_finite(f::P3109Format)
    return (BigInt(1) << f.K) - (1 + f.Δ * (1 + f.Σ))
end

"""
    n_positive_finite(f)

Number of positive finite values.
"""
function n_positive_finite(f::P3109Format)
    return (BigInt(1) << (f.K - f.Σ)) - 1 - f.Δ * (2 - f.Σ)
end

"""
    n_subnormals(f)

Number of positive subnormal values.
"""
function n_subnormals(f::P3109Format)
    return (BigInt(1) << (f.P - 1)) - 1
end

"""
    n_normals(f)

Number of positive normal values.
"""
function n_normals(f::P3109Format)
    return n_positive_finite(f) - n_subnormals(f)
end

"""
    n_binades(f)

Number of positive normal binades.
"""
function n_binades(f::P3109Format)
    return 2 * BigInt(bias(f)) - 1
end

"""
    ulp_in_binade(f, k)

Exact ULP size for normal binade `k`.
"""
ulp_in_binade(f::P3109Format, k::Integer) = dyadic_twopow(k + 1 - f.P)

"""
    ulp_of_abs_value(f, a)

Return the exact ULP associated with positive magnitude `a`.
"""
function ulp_of_abs_value(f::P3109Format, a::Rational{BigInt})
    a > 0 || throw(ArgumentError("a must be positive"))
    if a < xmin_normal(f)
        return xmin_subnormal(f)
    else
        k = floorlog2_dyadic(a)
        return ulp_in_binade(f, k)
    end
end

"""
    exact_mean_step_normal_binades(f)

Exact mean step size across the normal binades of `f`, weighting each binade equally.
"""
function exact_mean_step_normal_binades(f::P3109Format)
    B = Int(bias(f))
    s = 0 // 1
    for k in (1-B):(B-1)
        s += dyadic_twopow(k + 1 - f.P)
    end
    return s // (2 * B - 1)
end

"""
    exact_second_moment_step_normal_binades(f)

Exact second moment of the step-size distribution across normal binades.
"""
function exact_second_moment_step_normal_binades(f::P3109Format)
    B = Int(bias(f))
    s = 0 // 1
    for k in (1-B):(B-1)
        Δ = dyadic_twopow(k + 1 - f.P)
        s += Δ^2
    end
    return s // (2 * B - 1)
end

"""
    exact_variance_step_normal_binades(f)

Exact variance of the step-size distribution across normal binades.
"""
function exact_variance_step_normal_binades(f::P3109Format)
    μ = exact_mean_step_normal_binades(f)
    m2 = exact_second_moment_step_normal_binades(f)
    return m2 - μ^2
end

"""
    quantile_step_normal_binades(f, p)

Return the step size associated with quantile `p ∈ [0,1]` of the normal-binade index set.
"""
function quantile_step_normal_binades(f::P3109Format, p::Real)
    0.0 <= p <= 1.0 || throw(ArgumentError("p must lie in [0,1]"))
    kmin = 1 - Int(bias(f))
    kmax = Int(bias(f)) - 1
    ks = collect(kmin:kmax)
    idx = clamp(ceil(Int, p * length(ks)), 1, length(ks))
    return dyadic_twopow(ks[idx] + 1 - f.P)
end

"""
    value_step_count(f, x1, x2)

Return the number of representable value steps between two exactly representable
finite values `x1 < x2`.
"""
function value_step_count(f::P3109Format, x1::Rational{BigInt}, x2::Rational{BigInt})
    x1 < x2 || throw(ArgumentError("require x1 < x2"))
    cp1 = encode_finite(f, x1)
    cp2 = encode_finite(f, x2)
    return Int(cp2 - cp1)
end

"""
    roundtrip_decode_encode(f, cp)

Check whether `encode_finite(decode(f, cp)) == cp` for a finite code point.
Returns `true` for special values by convention.
"""
function roundtrip_decode_encode(f::P3109Format, cp::Integer)
    dv = decode(f, cp)
    dv isa FiniteValue || return true
    return encode_finite(f, dv.value) == UInt128(cp)
end

"""
    enumerate_positive_finite_values(f)

Enumerate all positive-half finite decoded values, including zero, in code-point order.
"""
function enumerate_positive_finite_values(f::P3109Format)
    vals = Rational{BigInt}[]
    for cp in 0:Int(cp_nor_max(f))
        dv = decode(f, cp)
        dv isa FiniteValue || continue
        push!(vals, dv.value)
    end
    return vals
end

"""
    enumerate_positive_steps(f)

Enumerate successive positive-half step sizes in code-point order.
"""
function enumerate_positive_steps(f::P3109Format)
    vals = enumerate_positive_finite_values(f)
    steps = Rational{BigInt}[]
    for i in 1:(length(vals)-1)
        push!(steps, vals[i+1] - vals[i])
    end
    return steps
end

"""
    verify_monotone_positive_half(f)

Verify that decoded positive-half finite values are strictly increasing with code point.
"""
function verify_monotone_positive_half(f::P3109Format)
    prev = nothing
    for cp in 0:Int(cp_nor_max(f))
        dv = decode(f, cp)
        dv isa FiniteValue || continue
        if prev !== nothing && !(prev < dv.value)
            return false
        end
        prev = dv.value
    end
    return true
end

"""
    exhaustive_roundtrip_check(f)

Exhaustively verify finite decode/encode roundtrips across the full code-point space of `f`.

Returns `true` if all finite code points roundtrip successfully.
"""
function exhaustive_roundtrip_check(f::P3109Format)
    for cp in 0:Int(twopow(f.K) - 1)
        if !roundtrip_decode_encode(f, cp)
            return false
        end
    end
    return true
end

"""
    exhaustive_monotonicity_check(f)

Verify strict monotonicity on the positive half by exhaustive enumeration.
"""
exhaustive_monotonicity_check(f::P3109Format) = verify_monotone_positive_half(f)

"""
    format_summary(f)

Return a named tuple summarizing core structural properties of `f`.
"""
function format_summary(f::P3109Format)
    return (
        K=f.K,
        P=f.P,
        Σ=f.Σ,
        Δ=f.Δ,
        W=exponent_width(f),
        bias=bias(f),
        m=significand_scale(f),
        cp_nor_max=cp_nor_max(f),
        n_finite=n_finite(f),
        n_positive_finite=n_positive_finite(f),
        n_subnormals=n_subnormals(f),
        n_normals=n_normals(f),
        n_binades=n_binades(f),
    )
end

end # module
=#
# values.jl — decode the ith code point of any Format to its exact value

#=
"""
    AllHexStringValuesOf(fmt) -> Vector{Qx64}

Return the exact rational values as "%a" sprintf strings for every code point
in `fmt`, in code-point order.  Special code points map to the canonical forms:

Special code points map to the canonical string forms:

- zero  → `0//1`
- +Inf  → "Inf"`
- -Inf  → "-Inf"`
- NaN   → "NaN"`
"""
function AllHexStringValuesOf(@nospecialize(fmt::Format))
    vals = String[]
    sizehint!(vals, Int(nValuesOf(fmt)))%

    for cp in Int128(0):Int128(cp_max(fmt))
        hexstring = HexStringValueOf(fmt, cp)
        push!(vals, hexstring)
    end
    return vals
end
=#

"""
    HexStringValueOf(fmt, cp) -> Rational{BigInt}

Return the exact HexString ("%a" sprintf) value of code point `cp`
in format `fmt`.

Special code points map to the canonical string forms:

- zero  → `0//1`
- +Inf  → "Inf"`
- -Inf  → "-Inf"`
- NaN   → "NaN"`

For finite numerical code points the result is an exact dyadic rational.
"""
function HexStringValueOf(@nospecialize(fmt::Format), cp::Integer)
    0 <= cp <= cp_max(fmt) || throw(ArgumentError("code point $cp out of range [0, $(cp_max(fmt))]"))

    cp == cp_nan(fmt) && return "NaN"
    cp == cp_zero(fmt) && return hex_sprintf(Qx64(0, 1))

    inf_cp = cp_inf(fmt)
    inf_cp !== nothing && cp == inf_cp && return hex_sprintf(Qx64(1, 0))
    ninf_cp = cp_neginf(fmt)
    ninf_cp !== nothing && cp == ninf_cp && return "NegInf"

    red = sign_reduce(fmt, cp)
    val = _decode_positive_half(fmt, red.cp_abs)
    return hex_sprintf(red.s == 0 ? val : -val)
end

"""
    AllValuesOf(fmt) -> Vector{Qx64}

Return the exact rational values for every code point in `fmt`, in
code-point order.  Special code points map to the canonical `Qx64` forms:
"""
function AllValuesOf(@nospecialize(fmt::Format))
    vals = Qx64[]
    sizehint!(vals, Int(nValuesOf(fmt)))
    cp_of_nan = cp_nan(fmt)
    for cp in Int128(0):Int128(cp_max(fmt))
        if cp == cp_of_nan
            push!(vals, Qx64(1, 0))
            continue
        end
        push!(vals, ValueOf(fmt, cp))
    end
    return vals
end

"""
    ValueOf(fmt, cp) -> Qx64

Return the exact `Qx64` value of code point `cp` in format `fmt`.

Special code points map to the canonical `Qx64` forms:

- zero  → `0//1`
- +Inf  → `1//0`
- -Inf  → `-1//0`
- NaN   → `0//0`

For finite numerical code points the result is an exact dyadic rational.
"""
@inline function _decode_signed_finite_unchecked(@nospecialize(fmt::Format{SignedFormat,T}), cp::Integer) where T
    h = sign_half_offset(fmt)
    if cp < h
        return _decode_positive_half(fmt, cp)
    end
    return -_decode_positive_half(fmt, cp - h)
end

@inline function ValueOf(@nospecialize(fmt::Format{UnsignedFormat,FiniteFormat}), cp::Integer)
    cpmax = cp_max(fmt)
    0 <= cp <= cpmax || throw(ArgumentError("code point $cp out of range [0, $cpmax]"))
    cp == 0 && return zero(Qx64)
    cp == cpmax && return Qx64(NaN)
    _decode_positive_half(fmt, cp)
end

@inline function ValueOf(@nospecialize(fmt::Format{UnsignedFormat,ExtendedFormat}), cp::Integer)
    cpmax = cp_max(fmt)
    cpnan = cpmax
    cpinf = cpnan - one(typeofcp(fmt))
    0 <= cp <= cpmax || throw(ArgumentError("code point $cp out of range [0, $cpmax]"))
    cp == 0 && return zero(Qx64)
    cp == cpnan && return Qx64(NaN)
    cp == cpinf && return Qx64(Inf)
    _decode_positive_half(fmt, cp)
end

@inline function ValueOf(@nospecialize(fmt::Format{SignedFormat,FiniteFormat}), cp::Integer)
    cpmax = cp_max(fmt)
    cpnan = cp_nan(fmt)
    0 <= cp <= cpmax || throw(ArgumentError("code point $cp out of range [0, $cpmax]"))
    cp == 0 && return zero(Qx64)
    cp == cpnan && return Qx64(NaN)
    _decode_signed_finite_unchecked(fmt, cp)
end

@inline function ValueOf(@nospecialize(fmt::Format{SignedFormat,ExtendedFormat}), cp::Integer)
    cpmax = cp_max(fmt)
    cpnan = cp_nan(fmt)
    cpinf = cpnan - one(typeofcp(fmt))
    cpninf = cpmax
    0 <= cp <= cpmax || throw(ArgumentError("code point $cp out of range [0, $cpmax]"))
    cp == 0 && return zero(Qx64)
    cp == cpnan && return Qx64(NaN)
    cp == cpinf && return Qx64(Inf)
    cp == cpninf && return Qx64(-Inf)
    _decode_signed_finite_unchecked(fmt, cp)
end

"""
    FiniteValueOf(fmt, cp) -> Qx64

Return the exact rational value of a finite code point `cp` in format `fmt`.

The code point must be finite (not NaN, not ±Inf).  Zero is accepted.
"""
function FiniteValueOf(@nospecialize(fmt::Format), cp::Integer)
    0 <= cp <= cp_max(fmt) || throw(ArgumentError("code point $cp out of range [0, $(cp_max(fmt))]"))
    cp == cp_zero(fmt) && return zero(Qx64)

    red = sign_reduce(fmt, cp)
    val = _decode_positive_half(fmt, red.cp_abs)
    return red.s == 0 ? val : -val
end

"""
    _decode_positive_half(fmt, cp_abs) -> Qx64

Decode a positive-half code point `cp_abs` (after sign reduction) to its
exact rational value.  Caller must ensure `cp_abs > 0`.

Let `m = twopow(P-1)` (the significand scale) and `B = ExponentBias`.

Subnormal (`cp_abs < m`):  value = cp_abs · twopow(2 − P − B)
Normal    (`cp_abs ≥ m`):  e = cp_abs ÷ m,  t = cp_abs mod m,  k = e − B
                           value = twopowk + t · twopow(k + 1 − P)
"""
@inline function _dyadic_from_mantissa_shift(mantissa::Integer, shift::Integer)
    if shift >= 0
        return (BigInt(mantissa) << Int(shift)) // BigInt(1)
    end
    return BigInt(mantissa) // (BigInt(1) << Int(-shift))
end

function _decode_positive_half_legacy(@nospecialize(fmt::Format), cp_abs::Integer)
    P = PrecisionOf(fmt)
    B = BigInt(ExponentBiasOf(fmt))
    m = BigInt(significand_scale(fmt))
    ca = BigInt(cp_abs)

    if ca < m
        return Qx64(dyadic_twopow(2 - P - Int(B))) * Qx64(ca)
    end

    e = fld(ca, m)
    t = mod(ca, m)
    k = Int(e - B)
    return Qx64(dyadic_twopow(k)) + Qx64(t) * Qx64(dyadic_twopow(k + 1 - P))
end

function _decode_positive_half(@nospecialize(fmt::Format), cp_abs::Integer)
    P = PrecisionOf(fmt)
    Bu = ExponentBiasOf(fmt)
    mu = significand_scale(fmt)

    if Bu <= typemax(Int) && mu <= typemax(Int) && cp_abs <= typemax(Int)
        B = Int(Bu)
        m = Int(mu)
        ca = Int(cp_abs)

        if ca < m
            return _dyadic_from_mantissa_shift(ca, 2 - P - B)
        end

        e = ca >> (P - 1)
        t = ca & (m - 1)
        mantissa = m + t
        return _dyadic_from_mantissa_shift(mantissa, e - B + 1 - P)
    end

    _decode_positive_half_legacy(fmt, cp_abs)
end

# =========================================================================
# Convenience: decode every code point at once
# =========================================================================

"""
    AllFiniteValuesOf(fmt) -> Vector{Qx64}

Return the exact rational values for every finite numerical code point
in `fmt`, in code-point order.  Zero is included; NaN and ±Inf are excluded.
"""
function AllFiniteValuesOf(@nospecialize(fmt::Format))
    vals = Qx64[]
    sizehint!(vals, Int(nFiniteValuesOf(fmt)) + 1)

    nan_cp = cp_nan(fmt)
    inf_cp = cp_inf(fmt)
    ninf_cp = cp_neginf(fmt)

    for cp in BigInt(0):BigInt(cp_max(fmt))
        cp == nan_cp && continue
        inf_cp !== nothing && cp == inf_cp && continue
        ninf_cp !== nothing && cp == ninf_cp && continue
        push!(vals, FiniteValueOf(fmt, cp))
    end
    return vals
end

"""
    AllPositiveFiniteValuesOf(fmt) -> Vector{Qx64}

Return the exact rational values for every strictly positive
finite code point in `fmt`, in code-point order.  Zero is excluded.
"""
function AllPositiveFiniteValuesOf(@nospecialize(fmt::Format))
    vals = Qx64[]
    sizehint!(vals, Int(nPosFiniteValuesOf(fmt)))

    cpstart = BigInt(cp_zero(fmt)) + 1
    cpend = BigInt(cp_nan(fmt)) - 1

    inf_cp = cp_inf(fmt)
    if inf_cp !== nothing
        cpend = min(cpend, BigInt(inf_cp) - 1)
    end

    for cp in cpstart:cpend
        push!(vals, FiniteValueOf(fmt, cp))
    end
    return vals
end

"""
    AllNonnegativeFiniteValuesOf(fmt) -> Vector{Qx64}

Return the exact rational values for every nonnegative finite code point
in `fmt`, in code-point order.  Zero is included.
"""
function AllNonnegativeFiniteValuesOf(@nospecialize(fmt::Format))
    vals = Qx64[]
    sizehint!(vals, Int(nNonnegFiniteValuesOf(fmt)))

    cpstart = BigInt(0)
    cpend = BigInt(cp_nan(fmt)) - 1

    inf_cp = cp_inf(fmt)
    if inf_cp !== nothing
        cpend = min(cpend, BigInt(inf_cp) - 1)
    end

    for cp in cpstart:cpend
        push!(vals, FiniteValueOf(fmt, cp))
    end
    return vals
end

"""
    AllNegativeFiniteValuesOf(fmt) -> Vector{Qx64}

Return the exact rational values for every strictly negative
finite code point in `fmt`, in code-point order.
"""
function AllNegativeFiniteValuesOf(@nospecialize(fmt::Format))
    is_unsigned(fmt) && return nothing

    posvals = AllPositiveFiniteValuesOf(fmt)
    map(-, posvals)
end

# =========================================================================
# Ordinal value accessors (1-based)
# =========================================================================

"""
    FiniteValueOfOrdinalPos(fmt, i) -> Qx64

Return the exact value of the `i`-th positive finite value (1-based),
ordered by increasing magnitude.

`i = 1` is the smallest positive value, `i = nPosFiniteValuesOf(fmt)` is
the largest.
"""
function FiniteValueOfOrdinalPos(@nospecialize(fmt::Format), i::Integer)
    n = nPosFiniteValuesOf(fmt)
    1 <= i <= n || throw(ArgumentError("ordinal $i out of range [1, $n]"))
    cp = BigInt(cp_zero(fmt)) + BigInt(i)
    return FiniteValueOf(fmt, cp)
end

"""
    FiniteValueOfOrdinalNeg(fmt, i) -> Qx64

Return the exact value of the `i`-th negative finite value (1-based),
ordered by increasing magnitude (i.e. `i = 1` is the negative value
closest to zero).

Only valid for signed formats.
"""
function FiniteValueOfOrdinalNeg(@nospecialize(fmt::Format{SignedFormat,T}), i::Integer) where T
    n = nNegFiniteValuesOf(fmt)
    1 <= i <= n || throw(ArgumentError("ordinal $i out of range [1, $n]"))
    cp = BigInt(cp_nan(fmt)) + BigInt(i)
    return FiniteValueOf(fmt, cp)
end

FiniteValueOfOrdinalNeg(@nospecialize(fmt::Format{UnsignedFormat,T}), i::Integer) where T =
    throw(ArgumentError("unsigned formats have no negative values"))

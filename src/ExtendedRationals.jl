module ExtendedRationals

export ExtendedRational, NegInf

import Base: NaN, Inf,
    +, -, *, /, ÷, ^,
    ==, isequal, isless, <, <=, >, >=,
    div, rem, divrem,
    fld, mod, fldmod, fld1, mod1, fldmod1,
    cld, gcd, lcm,
    abs, inv,
    iszero, isone, isfinite, isnan, isinf, signbit,
    numerator, denominator, Tuple, string, show,
    hash, promote_rule, convert,
    zero, one, oneunit, float, Rational,
    real, imag, conj

"""
    ExtendedRational <: Real

Exact rational number type closed under IEEE-style special values.

Stored as a pair of `BigInt` fields `(num, den)` with canonical forms:

- finite: `num//den` with `gcd(num, den) == 1` and `den > 0`
- zero: `0//1`
- `+Inf`: `1//0`
- `-Inf`: `-1//0`
- `NaN`: `0//0`

This type is intended to behave as a well-integrated Julia `Real` subtype.
It provides canonical equality and hashing, total ordering via `isless`,
promotion/conversion with common numeric types, exact finite arithmetic, and
IEEE-style closure policies for `NaN` and infinities.
"""
struct ExtendedRational <: Real
    num::BigInt
    den::BigInt

    function ExtendedRational(num::BigInt, den::BigInt)
        if den == 0
            if num == 0
                return new(BigInt(0), BigInt(0))     # NaN
            elseif num > 0
                return new(BigInt(1), BigInt(0))     # +Inf
            else
                return new(BigInt(-1), BigInt(0))    # -Inf
            end
        end

        if num == 0
            return new(BigInt(0), BigInt(1))
        end

        g = gcd(num, den)
        n = div(num, g)
        d = div(den, g)
        if d < 0
            n = -n
            d = -d
        end
        return new(n, d)
    end
end

const NegInf = ExtendedRational(-1, 0)

# -----------------------------------------------------------------------------
# constructors and basic accessors
# -----------------------------------------------------------------------------

Base.numerator(q::ExtendedRational) = q.num
Base.denominator(q::ExtendedRational) = q.den
Base.Tuple(q::ExtendedRational) = (q.num, q.den)

ExtendedRational(nd::Tuple{<:Integer,<:Integer}) = ExtendedRational(BigInt(nd[1]), BigInt(nd[2]))
ExtendedRational(num::Integer, den::Integer) = ExtendedRational(BigInt(num), BigInt(den))
ExtendedRational(num::BigInt) = ExtendedRational(num, one(BigInt))
ExtendedRational(num::Integer) = ExtendedRational(BigInt(num))
ExtendedRational(x::Rational{<:Integer}) = ExtendedRational(numerator(x), denominator(x))

function ExtendedRational(x::AbstractFloat)
    if isnan(x)
        return NaN(ExtendedRational)
    elseif isinf(x)
        return signbit(x) ? NegInf(ExtendedRational) : Inf(ExtendedRational)
    else
        r = Rational{BigInt}(x)
        return ExtendedRational(numerator(r), denominator(r))
    end
end

# ExtendedRational(NaN) = ExtendedRational(0, 0)
# ExtendedRational(Inf) = ExtendedRational(1, 0)
# ExtendedRational(NegInf) = ExtendedRational(-1, 0)

# Base.NaN(::ExtendedRational)  = ExtendedRational(0,0)
# Base.Inf(::ExtendedRational)  = ExtendedRational(1,0)
# NegInf(::ExtendedRational)    = ExtendedRational(-1,0)

# -----------------------------------------------------------------------------
# predicates and simple structure
# -----------------------------------------------------------------------------

@inline Base.isnan(q::ExtendedRational) = q.den == 0 && q.num == 0
@inline Base.isfinite(q::ExtendedRational) = q.den != 0
@inline Base.isinf(q::ExtendedRational) = q.den == 0 && q.num != 0
@inline Base.iszero(q::ExtendedRational) = q.den != 0 && q.num == 0
@inline Base.isone(q::ExtendedRational) = q.den == 1 && q.num == 1
@inline Base.signbit(q::ExtendedRational) = q.num < 0

Base.zero(::Type{ExtendedRational}) = ExtendedRational(0)
Base.zero(::ExtendedRational) = zero(ExtendedRational)
Base.one(::Type{ExtendedRational}) = ExtendedRational(1)
Base.one(::ExtendedRational) = one(ExtendedRational)
Base.oneunit(::Type{ExtendedRational}) = one(ExtendedRational)
Base.oneunit(::ExtendedRational) = one(ExtendedRational)

Base.real(q::ExtendedRational) = q
Base.imag(::ExtendedRational) = zero(ExtendedRational)
Base.conj(q::ExtendedRational) = q

# -----------------------------------------------------------------------------
# promotion and conversion
# -----------------------------------------------------------------------------

Base.convert(::Type{ExtendedRational}, q::ExtendedRational) = q
Base.convert(::Type{ExtendedRational}, x::Integer) = ExtendedRational(x)
Base.convert(::Type{ExtendedRational}, x::Rational{<:Integer}) = ExtendedRational(x)
Base.convert(::Type{ExtendedRational}, x::AbstractFloat) = ExtendedRational(x)

Base.promote_rule(::Type{ExtendedRational}, ::Type{<:Integer}) = ExtendedRational
Base.promote_rule(::Type{ExtendedRational}, ::Type{<:Rational}) = ExtendedRational
Base.promote_rule(::Type{ExtendedRational}, ::Type{<:AbstractFloat}) = ExtendedRational

function Base.convert(::Type{Rational{BigInt}}, q::ExtendedRational)
    isfinite(q) || throw(DomainError(q, "cannot convert NaN or infinity to Rational{BigInt}"))
    return numerator(q) // denominator(q)
end

function Base.convert(::Type{T}, q::ExtendedRational) where {T<:AbstractFloat}
    if isnan(q)
        return T(NaN)
    elseif isinf(q)
        return signbit(q) ? T(-Inf) : T(Inf)
    else
        return T(numerator(q)) / T(denominator(q))
    end
end

Base.float(q::ExtendedRational) = Float64(q)

# -----------------------------------------------------------------------------
# equality, ordering, hashing
# -----------------------------------------------------------------------------

function Base.:(==)(a::ExtendedRational, b::ExtendedRational)
    isnan(a) || isnan(b) && return false
    return a.num == b.num && a.den == b.den
end

function Base.isequal(a::ExtendedRational, b::ExtendedRational)
    isnan(a) && isnan(b) && return true
    return a.num == b.num && a.den == b.den
end

function Base.hash(q::ExtendedRational, h::UInt)
    # Hash is consistent within this type; no cross-type isequal is defined,
    # so agreement with hash(::Rational) is not required.
    if isnan(q)
        return hash((:ExtendedRational, :NaN), h)
    elseif isinf(q)
        return hash((:ExtendedRational, signbit(q) ? -Inf : Inf), h)
    else
        return hash((:ExtendedRational, q.num, q.den), h)
    end
end

"""
A total order compatible with canonical representation.
Finite values are ordered by value, then `-Inf`, then `+Inf`, then `NaN`
(NaN is the maximum element of this total order).
"""
function Base.isless(a::ExtendedRational, b::ExtendedRational)
    if isnan(b)
        return !isnan(a)          # everything is less than NaN, except NaN itself
    elseif isnan(a)
        return false              # NaN is not less than anything
    elseif isinf(a) && isinf(b)
        return a.num < b.num      # -Inf < +Inf
    elseif isinf(a)
        return signbit(a)         # -Inf < finite, +Inf is not
    elseif isinf(b)
        return !signbit(b)        # finite < +Inf, finite is not < -Inf
    else
        return a.num * b.den < b.num * a.den
    end
end

Base.:<(a::ExtendedRational, b::ExtendedRational) = isless(a, b)
Base.:<=(a::ExtendedRational, b::ExtendedRational) = !isless(b, a)
Base.:>(a::ExtendedRational, b::ExtendedRational) = isless(b, a)
Base.:>=(a::ExtendedRational, b::ExtendedRational) = !isless(a, b)

# -----------------------------------------------------------------------------
# helpers for finite arithmetic
# -----------------------------------------------------------------------------

function _add_finite(a::ExtendedRational, b::ExtendedRational)
    g = gcd(a.den, b.den)
    ad = div(a.den, g)
    bd = div(b.den, g)
    n = a.num * bd + b.num * ad
    iszero(n) && return zero(ExtendedRational)
    gn = gcd(abs(n), g)
    return ExtendedRational(div(n, gn), ad * div(b.den, gn))
end

function _mul_finite(a::ExtendedRational, b::ExtendedRational)
    g1 = gcd(abs(a.num), b.den)
    g2 = gcd(abs(b.num), a.den)
    n1 = div(a.num, g1)
    d1 = div(a.den, g2)
    n2 = div(b.num, g2)
    d2 = div(b.den, g1)
    return ExtendedRational(n1 * n2, d1 * d2)
end

# -----------------------------------------------------------------------------
# unary ops and arithmetic
# -----------------------------------------------------------------------------

Base.abs(q::ExtendedRational) = isnan(q) ? NaN(q) : ExtendedRational(abs(q.num), q.den)
Base.:+(q::ExtendedRational) = q
Base.:-(q::ExtendedRational) = ExtendedRational(-q.num, q.den)

function Base.:+(a::ExtendedRational, b::ExtendedRational)
    isnan(a) || isnan(b) && return NaN(ExtendedRational)
    isinf(a) && isinf(b) && return signbit(a) == signbit(b) ? a : NaN(ExtendedRational)
    isinf(a) && return a
    isinf(b) && return b
    return _add_finite(a, b)
end

Base.:-(a::ExtendedRational, b::ExtendedRational) = a + (-b)

function Base.:*(a::ExtendedRational, b::ExtendedRational)
    isnan(a) || isnan(b) && return NaN(ExtendedRational)
    (isinf(a) && iszero(b)) || (isinf(b) && iszero(a)) && return NaN(ExtendedRational)
    (isinf(a) || isinf(b)) && return signbit(a) != signbit(b) ?
                                     NegInf(ExtendedRational) : Inf(ExtendedRational)
    return _mul_finite(a, b)
end

function Base.inv(q::ExtendedRational)
    isnan(q) && return NaN(ExtendedRational)
    # Canonical zero has num == 0 and den == 1, so signbit(q) == false.
    # There is no negative-zero representation; inv(0) is always +Inf.
    iszero(q) && return Inf(ExtendedRational)
    isinf(q) && return zero(ExtendedRational)
    return ExtendedRational(q.den, q.num)
end

Base.:/(a::ExtendedRational, b::ExtendedRational) = a * inv(b)
# Base.:(÷)(a::ExtendedRational, b::ExtendedRational) = div(a, b)

function Base.:^(q::ExtendedRational, n::Integer)
    n == 0 && return one(ExtendedRational)       # q^0 = 1 for all q (IEEE convention)
    n < 0 && return inv(q)^(-n)
    isnan(q) && return NaN(ExtendedRational)
    if isinf(q)
        # n > 0 is guaranteed here (n == 0 and n < 0 handled above).
        iseven(n) && return Inf(ExtendedRational)
        return signbit(q) ? NegInf(ExtendedRational) : Inf(ExtendedRational)
    end
    return ExtendedRational(q.num^n, q.den^n)
end

# -----------------------------------------------------------------------------
# quotient / remainder family
# -----------------------------------------------------------------------------

function Base.div(a::ExtendedRational, b::ExtendedRational)
    isnan(a) || isnan(b) || iszero(b) || isinf(a) && return NaN(ExtendedRational)
    isinf(b) && return zero(ExtendedRational)
    return ExtendedRational(div(a.num * b.den, a.den * b.num), 1)
end

function Base.fld(a::ExtendedRational, b::ExtendedRational)
    isnan(a) || isnan(b) || iszero(b) || isinf(a) && return NaN(ExtendedRational)
    # For finite a: a / ±Inf = 0 exactly (exact rational arithmetic, not a limit),
    # so floor(0) = 0 regardless of the signs of a and b.
    isinf(b) && return zero(ExtendedRational)
    return ExtendedRational(fld(a.num * b.den, a.den * b.num), 1)
end

function Base.cld(a::ExtendedRational, b::ExtendedRational)
    isnan(a) || isnan(b) || iszero(b) || isinf(a) && return NaN(ExtendedRational)
    # For finite a: a / ±Inf = 0 exactly, so ceil(0) = 0.
    isinf(b) && return zero(ExtendedRational)
    return ExtendedRational(cld(a.num * b.den, a.den * b.num), 1)
end

function Base.rem(a::ExtendedRational, b::ExtendedRational)
    isnan(a) || isnan(b) || iszero(b) || isinf(a) && return NaN(ExtendedRational)
    isinf(b) && return a
    return a - div(a, b) * b
end

function Base.mod(a::ExtendedRational, b::ExtendedRational)
    isnan(a) || isnan(b) || iszero(b) || isinf(a) && return NaN(ExtendedRational)
    isinf(b) && return a
    return a - fld(a, b) * b
end

Base.fld1(a::ExtendedRational, b::ExtendedRational) = fld(a - one(ExtendedRational), b) + one(ExtendedRational)
Base.mod1(a::ExtendedRational, b::ExtendedRational) = mod(a - one(ExtendedRational), b) + one(ExtendedRational)

Base.divrem(a::ExtendedRational, b::ExtendedRational) = (div(a, b), rem(a, b))
Base.fldmod(a::ExtendedRational, b::ExtendedRational) = (fld(a, b), mod(a, b))
Base.fldmod1(a::ExtendedRational, b::ExtendedRational) = (fld1(a, b), mod1(a, b))

# -----------------------------------------------------------------------------
# gcd / lcm
# -----------------------------------------------------------------------------

"""
    gcd(a::ExtendedRational, b::ExtendedRational)

For finite values uses the exact rational identity
`gcd(a/b, c/d) = gcd(a, c) / lcm(b, d)`.

For non-finite values a closure policy is used:
- `gcd(NaN, x) = NaN`
- `gcd(±Inf, ±Inf) = +Inf`
- `gcd(±Inf, x_finite) = abs(x_finite)`
"""
function Base.gcd(a::ExtendedRational, b::ExtendedRational)
    isnan(a) || isnan(b) && return NaN(ExtendedRational)
    isinf(a) && isinf(b) && return Inf(ExtendedRational)
    isinf(a) && return abs(b)
    isinf(b) && return abs(a)
    return ExtendedRational(gcd(a.num, b.num), lcm(a.den, b.den))
end

"""
    lcm(a::ExtendedRational, b::ExtendedRational)

For finite values uses the exact rational identity
`lcm(a/b, c/d) = lcm(a, c) / gcd(b, d)`.

For non-finite values a closure policy is used:
- `lcm(NaN, x) = NaN`
- `lcm(0, x) = 0` for finite `x`
- `lcm(±Inf, x)` is `+Inf` when `x ≠ 0`
"""
function Base.lcm(a::ExtendedRational, b::ExtendedRational)
    isnan(a) || isnan(b) && return NaN(ExtendedRational)
    iszero(a) || iszero(b) && return zero(ExtendedRational)
    isinf(a) || isinf(b) && return Inf(ExtendedRational)
    return ExtendedRational(lcm(a.num, b.num), gcd(a.den, b.den))
end

# -----------------------------------------------------------------------------
# display
# -----------------------------------------------------------------------------

function Base.show(io::IO, q::ExtendedRational)
    if isnan(q)
        print(io, "NaN")
    elseif isinf(q)
        print(io, signbit(q) ? "-Inf" : "Inf")
    else
        show(io, q.num)
        print(io, "//")
        show(io, q.den)
    end
end

# Delegate string() to show to avoid duplicating logic.
Base.string(q::ExtendedRational) = sprint(show, q)

function Base.show(io::IO, ::MIME"text/plain", q::ExtendedRational)
    if isnan(q)
        print(io, "ExtendedRational(NaN)")
    elseif isinf(q)
        print(io, signbit(q) ? "ExtendedRational(-Inf)" : "ExtendedRational(Inf)")
    else
        print(io, "ExtendedRational(")
        show(io, q.num)
        print(io, "//")
        show(io, q.den)
        print(io, ")")
    end
end

# -----------------------------------------------------------------------------
# BitLengths
# -----------------------------------------------------------------------------

"""
    bitlen_unsigned(x::Integer) -> Int

Number of bits needed to store nonnegative integer `x` as an unsigned value.
By convention, `bitlen_unsigned(0) == 1`.
"""
function bitlen_unsigned(x::Integer)::Int
    x < 0 && throw(DomainError(x, "unsigned bit length requires x ≥ 0"))
    x == 0 && return 1
    return ndigits(x; base=2)
end

"""
    bitlen_signed(x::Integer) -> Int

Minimum number of bits needed to store integer `x` in a signed two's-complement
range `[-2^(w-1), 2^(w-1)-1]`.

Uses the closed-form identity:
- `x ≥ 0`: `ndigits(x; base=2) + 1`  (one extra bit for the sign)
- `x < 0`: `ndigits(-x - 1; base=2) + 1`  (mirrors the positive branch via ~x)

Special case `x ∈ {0, -1}` both require 1 bit.
"""
function bitlen_signed(x::Integer)::Int
    # 1-bit signed range is [-1, 0]; handle both endpoints explicitly.
    (x == 0 || x == -1) && return 1
    x > 0 && return ndigits(x; base=2) + 1
    return ndigits(-x - 1; base=2) + 1   # x < -1
end

"""
    reduced_num_den(r::Rational) -> (p, q)

Return `(numerator(r), denominator(r))` for `r` in lowest terms.
Julia's `Rational` always stores a positive denominator, so no sign
normalisation is required here.
"""
function reduced_num_den(r::Rational)
    return (numerator(r), denominator(r))
end

"""
    required_bits_unsigned(r::Rational) -> Int

Minimum unsigned bit width sufficient to store the reduced numerator magnitude
and reduced denominator of `r` exactly, assuming separate unsigned magnitudes.
"""
function required_bits_unsigned(r::Rational)::Int
    p, q = reduced_num_den(r)
    return max(bitlen_unsigned(abs(p)), bitlen_unsigned(q))
end

"""
    required_bits_signed(r::Rational) -> Int

Minimum signed bit width sufficient to store the reduced numerator and reduced
denominator exactly in a signed integer type.
"""
function required_bits_signed(r::Rational)::Int
    p, q = reduced_num_den(r)
    return max(bitlen_signed(p), bitlen_signed(q))
end

"""
    required_bits(r::Rational) -> NamedTuple

Return both unsigned- and signed-storage requirements for exact representation
of `r` as a reduced rational.
"""
function required_bits(r::Rational)
    p, q = reduced_num_den(r)
    return (
        numerator=p,
        denominator=q,
        unsigned_bits=required_bits_unsigned(r),
        signed_bits=required_bits_signed(r),
    )
end

# -----------------------------------------------------------------------------
# Helpers for dyadic-rational (binary-expression) construction
# -----------------------------------------------------------------------------

"""
    pow2(n::Integer) -> BigInt

Exact `2^n` as a `BigInt`, for `n ≥ 0`.
"""
function pow2(n::Integer)::BigInt
    n < 0 && throw(DomainError(n, "pow2 requires n ≥ 0"))
    return big(1) << n
end

"""
    dyadic_rational(num_terms, den_terms) -> Rational{BigInt}

Construct an exact rational from sums of signed powers of two.

Each term is a pair `(s, k)` meaning `s * 2^k`, where `s` is typically `+1` or `-1`
and `k ≥ 0`.

Example:
    # (2^5 - 2^2 + 1) / (2^7 + 2^3)
    r = dyadic_rational([(1,5), (-1,2), (1,0)], [(1,7), (1,3)])
"""
function dyadic_rational(num_terms, den_terms)::Rational{BigInt}
    num = sum(big(s) * pow2(k) for (s, k) in num_terms)
    den = sum(big(s) * pow2(k) for (s, k) in den_terms)
    den == 0 && throw(DivideError())
    return num // den
end

"""
    required_bits_binary_form(num_terms, den_terms)

Return reduced rational and bit requirements for a binary-expression rational.
"""
function required_bits_binary_form(num_terms, den_terms)
    r = dyadic_rational(num_terms, den_terms)
    return (rational=r, required_bits(r)...)
end

#=
# -----------------------------------------------------------------------------
# Examples
# -----------------------------------------------------------------------------

# 2^(-n) = 1 / 2^n
n = 7
r1 = big(1) // pow2(n)
@show r1
@show required_bits(r1)
# unsigned_bits = n+1 = 8
# signed_bits   = n+2 = 9

# (2^n - 1) / (2^n + 1)
r2 = (pow2(n) - 1) // (pow2(n) + 1)
@show r2
@show required_bits(r2)

# A more general dyadic-looking example:
# (2^10 + 2^6) / 2^12 = 17/64 after reduction
r3info = required_bits_binary_form([(1,10), (1,6)], [(1,12)])
@show r3info
=#

end # module ExtendedRationals

# module ExtendedRationals

# export ExtendedRational, NegInf

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
    ExtendedRational{I<:Integer} <: Real

Exact rational number type closed under IEEE-style special values.

Stored as a pair of `I` fields `(num, den)` with canonical forms:

- finite: `num//den` with `gcd(num, den) == 1` and `den > 0`
- zero: `0//1`
- `+Inf`: `1//0`
- `-Inf`: `-1//0`
- `NaN`: `0//0`

The type parameter `I` specifies the integer type used for numerator and
denominator (e.g. `BigInt`, `Int64`, `Int128`).

This type is intended to behave as a well-integrated Julia `Real` subtype.
It provides canonical equality and hashing, total ordering via `isless`,
promotion/conversion with common numeric types, exact finite arithmetic, and
IEEE-style closure policies for `NaN` and infinities.
"""
struct ExtendedRational{I<:Integer} <: Real
    num::I
    den::I

    function ExtendedRational{I}(num::I, den::I) where {I<:Integer}
        if den == 0
            if num == 0
                return new{I}(zero(I), zero(I))     # NaN
            elseif num > 0
                return new{I}(one(I), zero(I))       # +Inf
            else
                return new{I}(-one(I), zero(I))      # -Inf
            end
        end

        if num == 0
            return new{I}(zero(I), one(I))
        end

        if signbit(den)
            num = -num
            den = -den
        end

        g = gcd(num, den)
        n = div(num, g)
        d = div(den, g)

        return new{I}(n, d)
    end
end

# -----------------------------------------------------------------------------
# constructors and basic accessors
# -----------------------------------------------------------------------------

Base.numerator(q::ExtendedRational) = q.num
Base.denominator(q::ExtendedRational) = q.den
Base.Tuple(q::ExtendedRational) = (q.num, q.den)

ExtendedRational{I}(nd::Tuple{<:Integer,<:Integer}) where {I<:Integer} = ExtendedRational{I}(I(nd[1]), I(nd[2]))
ExtendedRational{I}(num::Integer, den::Integer) where {I<:Integer} = ExtendedRational{I}(I(num), I(den))
ExtendedRational{I}(num::Integer) where {I<:Integer} = ExtendedRational{I}(I(num), one(I))
ExtendedRational{I}(x::Rational{<:Integer}) where {I<:Integer} = ExtendedRational{I}(I(numerator(x)), I(denominator(x)))

function ExtendedRational{I}(x::AbstractFloat) where {I<:Integer}
    if isnan(x)
        return NaN(ExtendedRational{I})
    elseif isinf(x)
        return signbit(x) ? NegInf(ExtendedRational{I}) : Inf(ExtendedRational{I})
    else
        r = Rational{I}(x)
        return ExtendedRational{I}(numerator(r), denominator(r))
    end
end

const NegInf = ExtendedRational{BigInt}(-1, 0)

# -----------------------------------------------------------------------------
# predicates and simple structure
# -----------------------------------------------------------------------------

@inline Base.isnan(q::ExtendedRational) = q.den == 0 && q.num == 0
@inline Base.isfinite(q::ExtendedRational) = q.den != 0
@inline Base.isinf(q::ExtendedRational) = q.den == 0 && q.num != 0
@inline Base.iszero(q::ExtendedRational) = q.den != 0 && q.num == 0
@inline Base.isone(q::ExtendedRational) = q.den == 1 && q.num == 1
@inline Base.signbit(q::ExtendedRational) = q.num < 0

Base.zero(::Type{ExtendedRational{I}}) where {I<:Integer} = ExtendedRational{I}(zero(I))
Base.zero(::ExtendedRational{I}) where {I<:Integer} = zero(ExtendedRational{I})
Base.one(::Type{ExtendedRational{I}}) where {I<:Integer} = ExtendedRational{I}(one(I))
Base.one(::ExtendedRational{I}) where {I<:Integer} = one(ExtendedRational{I})
Base.oneunit(::Type{ExtendedRational{I}}) where {I<:Integer} = one(ExtendedRational{I})
Base.oneunit(::ExtendedRational{I}) where {I<:Integer} = one(ExtendedRational{I})

Base.real(q::ExtendedRational) = q
Base.imag(::ExtendedRational{I}) where {I<:Integer} = zero(ExtendedRational{I})
Base.conj(q::ExtendedRational) = q

# -----------------------------------------------------------------------------
# promotion and conversion
# -----------------------------------------------------------------------------

Base.convert(::Type{ExtendedRational{I}}, q::ExtendedRational{I}) where {I<:Integer} = q
Base.convert(::Type{ExtendedRational{I}}, q::ExtendedRational) where {I<:Integer} = ExtendedRational{I}(I(q.num), I(q.den))
Base.convert(::Type{ExtendedRational{I}}, x::Integer) where {I<:Integer} = ExtendedRational{I}(x)
Base.convert(::Type{ExtendedRational{I}}, x::Rational{<:Integer}) where {I<:Integer} = ExtendedRational{I}(x)
Base.convert(::Type{ExtendedRational{I}}, x::AbstractFloat) where {I<:Integer} = ExtendedRational{I}(x)

Base.promote_rule(::Type{ExtendedRational{I}}, ::Type{<:Integer}) where {I<:Integer} = ExtendedRational{I}
Base.promote_rule(::Type{ExtendedRational{I}}, ::Type{<:Rational}) where {I<:Integer} = ExtendedRational{I}
Base.promote_rule(::Type{ExtendedRational{I}}, ::Type{<:AbstractFloat}) where {I<:Integer} = ExtendedRational{I}
Base.promote_rule(::Type{ExtendedRational{I1}}, ::Type{ExtendedRational{I2}}) where {I1<:Integer, I2<:Integer} = ExtendedRational{promote_type(I1, I2)}

function Base.convert(::Type{Rational{I}}, q::ExtendedRational) where {I<:Integer}
    isfinite(q) || throw(DomainError(q, "cannot convert NaN or infinity to Rational{$I}"))
    return I(numerator(q)) // I(denominator(q))
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

function _add_finite(a::ExtendedRational{I}, b::ExtendedRational{I}) where {I<:Integer}
    g = gcd(a.den, b.den)
    ad = div(a.den, g)
    bd = div(b.den, g)
    n = a.num * bd + b.num * ad
    iszero(n) && return zero(ExtendedRational{I})
    gn = gcd(abs(n), g)
    return ExtendedRational{I}(div(n, gn), ad * div(b.den, gn))
end

function _mul_finite(a::ExtendedRational{I}, b::ExtendedRational{I}) where {I<:Integer}
    g1 = gcd(abs(a.num), b.den)
    g2 = gcd(abs(b.num), a.den)
    n1 = div(a.num, g1)
    d1 = div(a.den, g2)
    n2 = div(b.num, g2)
    d2 = div(b.den, g1)
    return ExtendedRational{I}(n1 * n2, d1 * d2)
end

# -----------------------------------------------------------------------------
# unary ops and arithmetic
# -----------------------------------------------------------------------------

Base.abs(q::ExtendedRational{I}) where {I<:Integer} = isnan(q) ? NaN(q) : ExtendedRational{I}(abs(q.num), q.den)
Base.:+(q::ExtendedRational) = q
Base.:-(q::ExtendedRational{I}) where {I<:Integer} = ExtendedRational{I}(-q.num, q.den)

function Base.:+(a::ExtendedRational{I}, b::ExtendedRational{I}) where {I<:Integer}
    isnan(a) || isnan(b) && return NaN(ExtendedRational{I})
    isinf(a) && isinf(b) && return signbit(a) == signbit(b) ? a : NaN(ExtendedRational{I})
    isinf(a) && return a
    isinf(b) && return b
    return _add_finite(a, b)
end

Base.:-(a::ExtendedRational, b::ExtendedRational) = a + (-b)

function Base.:*(a::ExtendedRational{I}, b::ExtendedRational{I}) where {I<:Integer}
    isnan(a) || isnan(b) && return NaN(ExtendedRational{I})
    (isinf(a) && iszero(b)) || (isinf(b) && iszero(a)) && return NaN(ExtendedRational{I})
    (isinf(a) || isinf(b)) && return signbit(a) != signbit(b) ?
                                     NegInf(ExtendedRational{I}) : Inf(ExtendedRational{I})
    return _mul_finite(a, b)
end

function Base.inv(q::ExtendedRational{I}) where {I<:Integer}
    isnan(q) && return NaN(ExtendedRational{I})
    iszero(q) && return Inf(ExtendedRational{I})
    isinf(q) && return zero(ExtendedRational{I})
    return ExtendedRational{I}(q.den, q.num)
end

Base.:/(a::ExtendedRational, b::ExtendedRational) = a * inv(b)

function Base.:^(q::ExtendedRational{I}, n::Integer) where {I<:Integer}
    n == 0 && return one(ExtendedRational{I})       # q^0 = 1 for all q (IEEE convention)
    n < 0 && return inv(q)^(-n)
    isnan(q) && return NaN(ExtendedRational{I})
    if isinf(q)
        iseven(n) && return Inf(ExtendedRational{I})
        return signbit(q) ? NegInf(ExtendedRational{I}) : Inf(ExtendedRational{I})
    end
    return ExtendedRational{I}(q.num^n, q.den^n)
end

# -----------------------------------------------------------------------------
# quotient / remainder family
# -----------------------------------------------------------------------------

function Base.div(a::ExtendedRational{I}, b::ExtendedRational{I}) where {I<:Integer}
    isnan(a) || isnan(b) || iszero(b) || isinf(a) && return NaN(ExtendedRational{I})
    isinf(b) && return zero(ExtendedRational{I})
    return ExtendedRational{I}(div(a.num * b.den, a.den * b.num), one(I))
end

function Base.fld(a::ExtendedRational{I}, b::ExtendedRational{I}) where {I<:Integer}
    isnan(a) || isnan(b) || iszero(b) || isinf(a) && return NaN(ExtendedRational{I})
    isinf(b) && return zero(ExtendedRational{I})
    return ExtendedRational{I}(fld(a.num * b.den, a.den * b.num), one(I))
end

function Base.cld(a::ExtendedRational{I}, b::ExtendedRational{I}) where {I<:Integer}
    isnan(a) || isnan(b) || iszero(b) || isinf(a) && return NaN(ExtendedRational{I})
    isinf(b) && return zero(ExtendedRational{I})
    return ExtendedRational{I}(cld(a.num * b.den, a.den * b.num), one(I))
end

function Base.rem(a::ExtendedRational{I}, b::ExtendedRational{I}) where {I<:Integer}
    isnan(a) || isnan(b) || iszero(b) || isinf(a) && return NaN(ExtendedRational{I})
    isinf(b) && return a
    return a - div(a, b) * b
end

function Base.mod(a::ExtendedRational{I}, b::ExtendedRational{I}) where {I<:Integer}
    isnan(a) || isnan(b) || iszero(b) || isinf(a) && return NaN(ExtendedRational{I})
    isinf(b) && return a
    return a - fld(a, b) * b
end

Base.fld1(a::ExtendedRational, b::ExtendedRational) = fld(a - one(a), b) + one(a)
Base.mod1(a::ExtendedRational, b::ExtendedRational) = mod(a - one(a), b) + one(a)

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
function Base.gcd(a::ExtendedRational{I}, b::ExtendedRational{I}) where {I<:Integer}
    isnan(a) || isnan(b) && return NaN(ExtendedRational{I})
    isinf(a) && isinf(b) && return Inf(ExtendedRational{I})
    isinf(a) && return abs(b)
    isinf(b) && return abs(a)
    return ExtendedRational{I}(gcd(a.num, b.num), lcm(a.den, b.den))
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
function Base.lcm(a::ExtendedRational{I}, b::ExtendedRational{I}) where {I<:Integer}
    isnan(a) || isnan(b) && return NaN(ExtendedRational{I})
    iszero(a) || iszero(b) && return zero(ExtendedRational{I})
    isinf(a) || isinf(b) && return Inf(ExtendedRational{I})
    return ExtendedRational{I}(lcm(a.num, b.num), gcd(a.den, b.den))
end

# -----------------------------------------------------------------------------
# display
# -----------------------------------------------------------------------------

Base.string(q::ExtendedRational) = sprint(show, q)

function Base.show(io::IO, ::MIME"text/plain", q::ExtendedRational)
    if isnan(q)
        print(io, "NaN")
    elseif isinf(q)
        print(io, signbit(q) ? "-Inf" : "Inf")
    else
        print(io, "(")
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

# end # module ExtendedRationals

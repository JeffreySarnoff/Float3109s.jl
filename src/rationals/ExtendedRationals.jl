# module Qx64s

# export Qx64, NegInf

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
    Qx64{I<:Integer} <: Real

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
struct Qx64{I<:Integer} <: Real
    num::I
    den::I

    function Qx64{I}(num::I, den::I) where {I<:Integer}
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

Base.numerator(q::Qx64) = q.num
Base.denominator(q::Qx64) = q.den
Base.Tuple(q::Qx64) = (q.num, q.den)

Qx64{I}(nd::Tuple{<:Integer,<:Integer}) where {I<:Integer} = Qx64{I}(I(nd[1]), I(nd[2]))
Qx64{I}(num::Integer, den::Integer) where {I<:Integer} = Qx64{I}(I(num), I(den))
Qx64{I}(num::Integer) where {I<:Integer} = Qx64{I}(I(num), one(I))
Qx64{I}(x::Rational{<:Integer}) where {I<:Integer} = Qx64{I}(I(numerator(x)), I(denominator(x)))

function Qx64{I}(x::AbstractFloat) where {I<:Integer}
    if isnan(x)
        return NaN(Qx64{I})
    elseif isinf(x)
        return signbit(x) ? NegInf(Qx64{I}) : Inf(Qx64{I})
    else
        r = Rational{I}(x)
        return Qx64{I}(numerator(r), denominator(r))
    end
end

const NegInf = Qx64{BigInt}(-1, 0)

# -----------------------------------------------------------------------------
# predicates and simple structure
# -----------------------------------------------------------------------------

@inline Base.isnan(q::Qx64) = q.den == 0 && q.num == 0
@inline Base.isfinite(q::Qx64) = q.den != 0
@inline Base.isinf(q::Qx64) = q.den == 0 && q.num != 0
@inline Base.iszero(q::Qx64) = q.den != 0 && q.num == 0
@inline Base.isone(q::Qx64) = q.den == 1 && q.num == 1
@inline Base.signbit(q::Qx64) = q.num < 0

Base.zero(::Type{Qx64{I}}) where {I<:Integer} = Qx64{I}(zero(I))
Base.zero(::Qx64{I}) where {I<:Integer} = zero(Qx64{I})
Base.one(::Type{Qx64{I}}) where {I<:Integer} = Qx64{I}(one(I))
Base.one(::Qx64{I}) where {I<:Integer} = one(Qx64{I})
Base.oneunit(::Type{Qx64{I}}) where {I<:Integer} = one(Qx64{I})
Base.oneunit(::Qx64{I}) where {I<:Integer} = one(Qx64{I})

Base.real(q::Qx64) = q
Base.imag(::Qx64{I}) where {I<:Integer} = zero(Qx64{I})
Base.conj(q::Qx64) = q

# -----------------------------------------------------------------------------
# promotion and conversion
# -----------------------------------------------------------------------------

Base.convert(::Type{Qx64{I}}, q::Qx64{I}) where {I<:Integer} = q
Base.convert(::Type{Qx64{I}}, q::Qx64) where {I<:Integer} = Qx64{I}(I(q.num), I(q.den))
Base.convert(::Type{Qx64{I}}, x::Integer) where {I<:Integer} = Qx64{I}(x)
Base.convert(::Type{Qx64{I}}, x::Rational{<:Integer}) where {I<:Integer} = Qx64{I}(x)
Base.convert(::Type{Qx64{I}}, x::AbstractFloat) where {I<:Integer} = Qx64{I}(x)

Base.promote_rule(::Type{Qx64{I}}, ::Type{<:Integer}) where {I<:Integer} = Qx64{I}
Base.promote_rule(::Type{Qx64{I}}, ::Type{<:Rational}) where {I<:Integer} = Qx64{I}
Base.promote_rule(::Type{Qx64{I}}, ::Type{<:AbstractFloat}) where {I<:Integer} = Qx64{I}
Base.promote_rule(::Type{Qx64{I1}}, ::Type{Qx64{I2}}) where {I1<:Integer,I2<:Integer} = Qx64{promote_type(I1, I2)}

function Base.convert(::Type{Rational{I}}, q::Qx64) where {I<:Integer}
    isfinite(q) || throw(DomainError(q, "cannot convert NaN or infinity to Rational{$I}"))
    return I(numerator(q)) // I(denominator(q))
end

function Base.convert(::Type{T}, q::Qx64) where {T<:AbstractFloat}
    if isnan(q)
        return T(NaN)
    elseif isinf(q)
        return signbit(q) ? T(-Inf) : T(Inf)
    else
        return T(numerator(q)) / T(denominator(q))
    end
end

Base.float(q::Qx64) = Float64(q)

# -----------------------------------------------------------------------------
# equality, ordering, hashing
# -----------------------------------------------------------------------------

function Base.:(==)(a::Qx64, b::Qx64)
    isnan(a) || isnan(b) && return false
    return a.num == b.num && a.den == b.den
end

function Base.isequal(a::Qx64, b::Qx64)
    isnan(a) && isnan(b) && return true
    return a.num == b.num && a.den == b.den
end

function Base.hash(q::Qx64, h::UInt)
    if isnan(q)
        return hash((:Qx64, :NaN), h)
    elseif isinf(q)
        return hash((:Qx64, signbit(q) ? -Inf : Inf), h)
    else
        return hash((:Qx64, q.num, q.den), h)
    end
end

"""
A total order compatible with canonical representation.
Finite values are ordered by value, then `-Inf`, then `+Inf`, then `NaN`
(NaN is the maximum element of this total order).
"""
function Base.isless(a::Qx64, b::Qx64)
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

Base.:<(a::Qx64, b::Qx64) = isless(a, b)
Base.:<=(a::Qx64, b::Qx64) = !isless(b, a)
Base.:>(a::Qx64, b::Qx64) = isless(b, a)
Base.:>=(a::Qx64, b::Qx64) = !isless(a, b)

# -----------------------------------------------------------------------------
# helpers for finite arithmetic
# -----------------------------------------------------------------------------

function _add_finite(a::Qx64{I}, b::Qx64{I}) where {I<:Integer}
    g = gcd(a.den, b.den)
    ad = div(a.den, g)
    bd = div(b.den, g)
    n = a.num * bd + b.num * ad
    iszero(n) && return zero(Qx64{I})
    gn = gcd(abs(n), g)
    return Qx64{I}(div(n, gn), ad * div(b.den, gn))
end

function _mul_finite(a::Qx64{I}, b::Qx64{I}) where {I<:Integer}
    g1 = gcd(abs(a.num), b.den)
    g2 = gcd(abs(b.num), a.den)
    n1 = div(a.num, g1)
    d1 = div(a.den, g2)
    n2 = div(b.num, g2)
    d2 = div(b.den, g1)
    return Qx64{I}(n1 * n2, d1 * d2)
end

# -----------------------------------------------------------------------------
# unary ops and arithmetic
# -----------------------------------------------------------------------------

Base.abs(q::Qx64{I}) where {I<:Integer} = isnan(q) ? NaN(q) : Qx64{I}(abs(q.num), q.den)
Base.:+(q::Qx64) = q
Base.:-(q::Qx64{I}) where {I<:Integer} = Qx64{I}(-q.num, q.den)

function Base.:+(a::Qx64{I}, b::Qx64{I}) where {I<:Integer}
    isnan(a) || isnan(b) && return NaN(Qx64{I})
    isinf(a) && isinf(b) && return signbit(a) == signbit(b) ? a : NaN(Qx64{I})
    isinf(a) && return a
    isinf(b) && return b
    return _add_finite(a, b)
end

Base.:-(a::Qx64, b::Qx64) = a + (-b)

function Base.:*(a::Qx64{I}, b::Qx64{I}) where {I<:Integer}
    isnan(a) || isnan(b) && return NaN(Qx64{I})
    (isinf(a) && iszero(b)) || (isinf(b) && iszero(a)) && return NaN(Qx64{I})
    (isinf(a) || isinf(b)) && return signbit(a) != signbit(b) ?
                                     NegInf(Qx64{I}) : Inf(Qx64{I})
    return _mul_finite(a, b)
end

function Base.inv(q::Qx64{I}) where {I<:Integer}
    isnan(q) && return NaN(Qx64{I})
    iszero(q) && return Inf(Qx64{I})
    isinf(q) && return zero(Qx64{I})
    return Qx64{I}(q.den, q.num)
end

Base.:/(a::Qx64, b::Qx64) = a * inv(b)

function Base.:^(q::Qx64{I}, n::Integer) where {I<:Integer}
    n == 0 && return one(Qx64{I})       # q^0 = 1 for all q (IEEE convention)
    n < 0 && return inv(q)^(-n)
    isnan(q) && return NaN(Qx64{I})
    if isinf(q)
        iseven(n) && return Inf(Qx64{I})
        return signbit(q) ? NegInf(Qx64{I}) : Inf(Qx64{I})
    end
    return Qx64{I}(q.num^n, q.den^n)
end

# -----------------------------------------------------------------------------
# quotient / remainder family
# -----------------------------------------------------------------------------

function Base.div(a::Qx64{I}, b::Qx64{I}) where {I<:Integer}
    isnan(a) || isnan(b) || iszero(b) || isinf(a) && return NaN(Qx64{I})
    isinf(b) && return zero(Qx64{I})
    return Qx64{I}(div(a.num * b.den, a.den * b.num), one(I))
end

function Base.fld(a::Qx64{I}, b::Qx64{I}) where {I<:Integer}
    isnan(a) || isnan(b) || iszero(b) || isinf(a) && return NaN(Qx64{I})
    isinf(b) && return zero(Qx64{I})
    return Qx64{I}(fld(a.num * b.den, a.den * b.num), one(I))
end

function Base.cld(a::Qx64{I}, b::Qx64{I}) where {I<:Integer}
    isnan(a) || isnan(b) || iszero(b) || isinf(a) && return NaN(Qx64{I})
    isinf(b) && return zero(Qx64{I})
    return Qx64{I}(cld(a.num * b.den, a.den * b.num), one(I))
end

function Base.rem(a::Qx64{I}, b::Qx64{I}) where {I<:Integer}
    isnan(a) || isnan(b) || iszero(b) || isinf(a) && return NaN(Qx64{I})
    isinf(b) && return a
    return a - div(a, b) * b
end

function Base.mod(a::Qx64{I}, b::Qx64{I}) where {I<:Integer}
    isnan(a) || isnan(b) || iszero(b) || isinf(a) && return NaN(Qx64{I})
    isinf(b) && return a
    return a - fld(a, b) * b
end

Base.fld1(a::Qx64, b::Qx64) = fld(a - one(a), b) + one(a)
Base.mod1(a::Qx64, b::Qx64) = mod(a - one(a), b) + one(a)

Base.divrem(a::Qx64, b::Qx64) = (div(a, b), rem(a, b))
Base.fldmod(a::Qx64, b::Qx64) = (fld(a, b), mod(a, b))
Base.fldmod1(a::Qx64, b::Qx64) = (fld1(a, b), mod1(a, b))

# -----------------------------------------------------------------------------
# gcd / lcm
# -----------------------------------------------------------------------------

"""
    gcd(a::Qx64, b::Qx64)

For finite values uses the exact rational identity
`gcd(a/b, c/d) = gcd(a, c) / lcm(b, d)`.

For non-finite values a closure policy is used:
- `gcd(NaN, x) = NaN`
- `gcd(±Inf, ±Inf) = +Inf`
- `gcd(±Inf, x_finite) = abs(x_finite)`
"""
function Base.gcd(a::Qx64{I}, b::Qx64{I}) where {I<:Integer}
    isnan(a) || isnan(b) && return NaN(Qx64{I})
    isinf(a) && isinf(b) && return Inf(Qx64{I})
    isinf(a) && return abs(b)
    isinf(b) && return abs(a)
    return Qx64{I}(gcd(a.num, b.num), lcm(a.den, b.den))
end

"""
    lcm(a::Qx64, b::Qx64)

For finite values uses the exact rational identity
`lcm(a/b, c/d) = lcm(a, c) / gcd(b, d)`.

For non-finite values a closure policy is used:
- `lcm(NaN, x) = NaN`
- `lcm(0, x) = 0` for finite `x`
- `lcm(±Inf, x)` is `+Inf` when `x ≠ 0`
"""
function Base.lcm(a::Qx64{I}, b::Qx64{I}) where {I<:Integer}
    isnan(a) || isnan(b) && return NaN(Qx64{I})
    iszero(a) || iszero(b) && return zero(Qx64{I})
    isinf(a) || isinf(b) && return Inf(Qx64{I})
    return Qx64{I}(lcm(a.num, b.num), gcd(a.den, b.den))
end

# -----------------------------------------------------------------------------
# display
# -----------------------------------------------------------------------------

Base.string(q::Qx64) = sprint(show, q)

function Base.show(io::IO, ::MIME"text/plain", q::Qx64)
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

# end # module Qx64s

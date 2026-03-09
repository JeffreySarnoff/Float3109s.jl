import Base: NaN, Inf, +, -, *, ÷, /, ^,
    div, rem, divrem,
    fld, mod, fldmod, fld1, mod1, fldmod1,
    cld, gcd, lcm,
    iszero, isone, isfinite, isnan, signbit, convert,
    numerator, denominator, Tuple

"""
    ClosedRational

A type representing the set of large rational numbers enhanced to support NaN.

note: Julia abstract type Real supports NaN
"""

struct ClosedRational <: Real
    num::BigInt
    den::BigInt
    function ClosedRational(num::BigInt, den::BigInt)
        if den == 0
            if num == 0
                new(BigInt(0), BigInt(0))       # NaN
            elseif num > 0
                new(BigInt(1), BigInt(0))       # +Inf
            else
                new(BigInt(-1), BigInt(0))      # -Inf
            end
        else
            g = gcd(num, den)
            n = div(num, g)
            d = div(den, g)
            if d < 0
                n = -n
                d = -d
            end
            new(n, d)
        end
    end
end

numerator(q::ClosedRational) = q.num
denominator(q::ClosedRational) = q.den

Base.Tuple(q::ClosedRational) = (q.num, q.den)

ClosedRational(nd::Tuple{I,I}) where I<:BigInt = ClosedRational(nd...)
ClosedRational(nd::Tuple{I,I}) where I<:Integer = ClosedRational(BigInt(nd[1]), BigInt(nd[2]))

ClosedRational(num::I, den::I) where {I<:Integer} =
    ClosedRational(BigInt(num), BigInt(den))

ClosedRational(num::BigInt) = ClosedRational(num, one(BigInt))
ClosedRational(num::Integer) = ClosedRational(BigInt(num))

Base.NaN(::Type{ClosedRational}) = ClosedRational(0, 0)
Base.Inf(::Type{ClosedRational}) = ClosedRational(1, 0)
NegInf(::Type{ClosedRational}) = ClosedRational(-1, 0)

Base.NaN(q::ClosedRational) = Base.NaN(ClosedRational)
Base.Inf(q::ClosedRational) = Base.Inf(ClosedRational)
NegInf(q::ClosedRational) = NegInf(ClosedRational)

ClosedRational(x::AbstractFloat) =
    if isnan(x)
        ClosedRational(0, 0)  # NaN
    elseif isinf(x)
        if signbit(x)
            ClosedRational(-1, 0) # -Inf
        else
            ClosedRational(1, 0)  # +Inf
        end
    else
        r = Rational{BigInt}(x)
        ClosedRational(r.num, r.den)
    end

ClosedRational(x::Rational{I}) where I<:Integer = ClosedRational(x.num, x.den)

# predicates

Base.isnan(q::ClosedRational) = q.num == 0 && q.den == 0
Base.isfinite(q::ClosedRational) = q.den != 0
Base.iszero(q::ClosedRational) = q.num == 0 && q.den != 0
Base.isone(q::ClosedRational) = q.num == 1 && q.den == 1
Base.signbit(q::ClosedRational) = q.num < 0

# abs

Base.abs(q::ClosedRational) =
    isnan(q) ? NaN(ClosedRational) : ClosedRational(abs(q.num), q.den)

# unary +/-

Base.(+)(q::ClosedRational) = q
Base.(-)(q::ClosedRational) = ClosedRational(-q.num, q.den)

# addition

function +(a::ClosedRational, b::ClosedRational)
    if isnan(a) || isnan(b)
        return NaN(ClosedRational)
    elseif !isfinite(a) && !isfinite(b)
        return signbit(a) == signbit(b) ? a : NaN(ClosedRational)
    elseif !isfinite(a)
        return a
    elseif !isfinite(b)
        return b
    else
        return ClosedRational(a.num * b.den + b.num * a.den, a.den * b.den)
    end
end

# subtraction

function -(a::ClosedRational, b::ClosedRational)
    if isnan(a) || isnan(b)
        return NaN(ClosedRational)
    elseif !isfinite(a) && !isfinite(b)
        return signbit(a) != signbit(b) ? a : NaN(ClosedRational)
    elseif !isfinite(a)
        return a
    elseif !isfinite(b)
        return -b
    else
        return ClosedRational(a.num * b.den - b.num * a.den, a.den * b.den)
    end
end

# multiplication

function *(a::ClosedRational, b::ClosedRational)
    if isnan(a) || isnan(b)
        return NaN(ClosedRational)
    elseif !isfinite(a)
        iszero(b) ? NaN(ClosedRational) :
        ClosedRational(signbit(a) != signbit(b) ? -1 : 1, 0)
    elseif !isfinite(b)
        iszero(a) ? NaN(ClosedRational) :
        ClosedRational(signbit(a) != signbit(b) ? -1 : 1, 0)
    else
        return ClosedRational(a.num * b.num, a.den * b.den)
    end
end

# div (truncated toward zero), fld (toward -∞), cld (toward +∞)

function div(a::ClosedRational, b::ClosedRational)
    if isnan(a) || isnan(b) || !isfinite(a) || iszero(b)
        return NaN(ClosedRational)
    elseif !isfinite(b)
        return ClosedRational(0)
    else
        return ClosedRational(div(a.num * b.den, a.den * b.num), one(BigInt))
    end
end

function fld(a::ClosedRational, b::ClosedRational)
    if isnan(a) || isnan(b) || !isfinite(a) || iszero(b)
        return NaN(ClosedRational)
    elseif !isfinite(b)
        return ClosedRational(signbit(a) != signbit(b) ? -1 : 0)
    else
        return ClosedRational(fld(a.num * b.den, a.den * b.num), one(BigInt))
    end
end

function cld(a::ClosedRational, b::ClosedRational)
    if isnan(a) || isnan(b) || !isfinite(a) || iszero(b)
        return NaN(ClosedRational)
    elseif !isfinite(b)
        return ClosedRational(signbit(a) != signbit(b) ? 0 : 1)
    else
        return ClosedRational(cld(a.num * b.den, a.den * b.num), one(BigInt))
    end
end

# rem(a, b) = a - div(a, b) * b  (sign of a)
# mod(a, b) = a - fld(a, b) * b  (sign of b)

function rem(a::ClosedRational, b::ClosedRational)
    if isnan(a) || isnan(b) || !isfinite(a) || iszero(b)
        return NaN(ClosedRational)
    elseif !isfinite(b)
        return a
    else
        return a - div(a, b) * b
    end
end

function mod(a::ClosedRational, b::ClosedRational)
    if isnan(a) || isnan(b) || !isfinite(a) || iszero(b)
        return NaN(ClosedRational)
    elseif !isfinite(b)
        return a
    else
        return a - fld(a, b) * b
    end
end

# fld1(a, b) = fld(a-1, b) + 1   (1-based floored division)
# mod1(a, b) = mod(a-1, b) + 1   (result in [1, b] instead of [0, b))

function fld1(a::ClosedRational, b::ClosedRational)
    return fld(a - ClosedRational(1), b) + ClosedRational(1)
end

function mod1(a::ClosedRational, b::ClosedRational)
    return mod(a - ClosedRational(1), b) + ClosedRational(1)
end

divrem(a::ClosedRational, b::ClosedRational) = (div(a, b), rem(a, b))
fldmod(a::ClosedRational, b::ClosedRational) = (fld(a, b), mod(a, b))
fldmod1(a::ClosedRational, b::ClosedRational) = (fld1(a, b), mod1(a, b))

# gcd(a/b, c/d) = gcd(a,c) / lcm(b,d)
# lcm(a/b, c/d) = lcm(a,c) / gcd(b,d)

function gcd(a::ClosedRational, b::ClosedRational)
    if isnan(a) || isnan(b)
        return NaN(ClosedRational)
    elseif !isfinite(a) && !isfinite(b)
        return Inf(ClosedRational)
    elseif !isfinite(a)
        return abs(b)
    elseif !isfinite(b)
        return abs(a)
    else
        return ClosedRational(gcd(a.num, b.num), lcm(a.den, b.den))
    end
end

function lcm(a::ClosedRational, b::ClosedRational)
    if isnan(a) || isnan(b)
        return NaN(ClosedRational)
    elseif iszero(a) || iszero(b)
        return ClosedRational(0)
    elseif !isfinite(a) || !isfinite(b)
        return Inf(ClosedRational)
    else
        return ClosedRational(lcm(a.num, b.num), gcd(a.den, b.den))
    end
end

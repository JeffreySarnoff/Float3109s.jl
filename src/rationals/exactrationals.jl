module ExactRationals

export FastRational,
    DyadicRational,
    numerator,
    denominator,
    exponent,
    dyadic_denominator,
    isdyadic,
    to_fastrational,
    to_dyadic,
    dyadic_div,
    run_selftests,
    run_benchmarks

import Base: +, -, *, /, ==, isless, hash, show, zero, one,
    abs, sign, inv, convert, promote_rule

############################
# FastRational definition
############################

struct FastRational{T<:Integer} <: Real
    num::T
    den::T

    function FastRational{T}(n::T, d::T) where {T<:Integer}
        d == 0 && throw(DivideError())
        if n == 0
            return new{T}(zero(T), one(T))
        end
        if d < 0
            n = -n
            d = -d
        end
        g = gcd(n, d)
        new{T}(div(n, g), div(d, g))
    end
end

FastRational(n::T, d::T) where {T<:Integer} = FastRational{T}(n, d)
FastRational(n::Integer, d::Integer) = FastRational(promote(n, d)...)
FastRational(n::Integer) = FastRational(n, 1)

numerator(x::FastRational) = x.num
denominator(x::FastRational) = x.den

Base.show(io::IO, x::FastRational) = print(io, x.num, "//", x.den)

Base.zero(::Type{FastRational{T}}) where {T<:Integer} = FastRational{T}(0, 1)
Base.one(::Type{FastRational{T}}) where {T<:Integer} = FastRational{T}(1, 1)

Base.abs(x::FastRational{T}) where {T<:Integer} = FastRational{T}(abs(x.num), x.den)
Base.sign(x::FastRational) = sign(x.num)

Base.:(==)(x::FastRational, y::FastRational) = x.num == y.num && x.den == y.den
Base.hash(x::FastRational, h::UInt) = hash(x.den, hash(x.num, h))

Base.isless(x::FastRational, y::FastRational) = x.num * y.den < y.num * x.den

Base.:-(x::FastRational{T}) where {T<:Integer} = FastRational{T}(-x.num, x.den)

############################
# FastRational kernels
############################

@inline function Base.:+(x::FastRational{T}, y::FastRational{T}) where {T<:Integer}
    x.num == 0 && return y
    y.num == 0 && return x
    g = gcd(x.den, y.den)
    xd = div(x.den, g)
    yd = div(y.den, g)
    n = x.num * yd + y.num * xd
    n == 0 && return FastRational{T}(0, 1)
    d = x.den * yd
    g2 = gcd(n, g)
    FastRational{T}(div(n, g2), div(d, g2))
end

@inline function Base.:-(x::FastRational{T}, y::FastRational{T}) where {T<:Integer}
    x.num == 0 && return FastRational{T}(-y.num, y.den)
    y.num == 0 && return x
    g = gcd(x.den, y.den)
    xd = div(x.den, g)
    yd = div(y.den, g)
    n = x.num * yd - y.num * xd
    n == 0 && return FastRational{T}(0, 1)
    d = x.den * yd
    g2 = gcd(n, g)
    FastRational{T}(div(n, g2), div(d, g2))
end

@inline function Base.:*(x::FastRational{T}, y::FastRational{T}) where {T<:Integer}
    x.num == 0 && return zero(x)
    y.num == 0 && return zero(y)
    g1 = gcd(x.num, y.den)
    g2 = gcd(y.num, x.den)
    a = div(x.num, g1)
    d = div(y.den, g1)
    c = div(y.num, g2)
    b = div(x.den, g2)
    FastRational{T}(a * c, b * d)
end

@inline function Base.:/(x::FastRational{T}, y::FastRational{T}) where {T<:Integer}
    y.num == 0 && throw(DivideError())
    x.num == 0 && return zero(x)
    g1 = gcd(x.num, y.num)
    g2 = gcd(y.den, x.den)
    a = div(x.num, g1)
    c = div(y.num, g1)
    d = div(y.den, g2)
    b = div(x.den, g2)
    FastRational{T}(a * d, b * c)
end

############################
# DyadicRational definition
############################

struct DyadicRational{T<:Integer} <: Real
    num::T
    exp::Int

    function DyadicRational{T}(num::T, exp::Integer) where {T<:Integer}
        exp < 0 && throw(ArgumentError("exp must be ≥0"))
        if num == 0
            return new{T}(zero(T), 0)
        end
        tz = trailing_zeros(num)
        s = min(tz, exp)
        new{T}(num >> s, exp - s)
    end
end

DyadicRational(num::T, exp::Integer) where {T<:Integer} = DyadicRational{T}(num, exp)
DyadicRational(n::Integer) = DyadicRational(n, 0)

numerator(x::DyadicRational) = x.num
exponent(x::DyadicRational) = x.exp
dyadic_denominator(x::DyadicRational) = big(1) << x.exp
denominator(x::DyadicRational) = dyadic_denominator(x)

Base.show(io::IO, x::DyadicRational) =
    x.exp == 0 ? print(io, x.num) : print(io, x.num, "//2^", x.exp)

Base.zero(::Type{DyadicRational{T}}) where {T<:Integer} = DyadicRational{T}(0, 0)
Base.one(::Type{DyadicRational{T}}) where {T<:Integer} = DyadicRational{T}(1, 0)

Base.abs(x::DyadicRational{T}) where {T<:Integer} = DyadicRational{T}(abs(x.num), x.exp)
Base.sign(x::DyadicRational) = sign(x.num)

Base.:(==)(x::DyadicRational, y::DyadicRational) = x.num == y.num && x.exp == y.exp
Base.hash(x::DyadicRational, h::UInt) = hash(x.exp, hash(x.num, h))

function Base.isless(x::DyadicRational{T}, y::DyadicRational{T}) where {T<:Integer}
    if x.exp == y.exp
        x.num < y.num
    elseif x.exp < y.exp
        (x.num << (y.exp - x.exp)) < y.num
    else
        x.num < (y.num << (x.exp - y.exp))
    end
end

############################
# Dyadic kernels
############################

@inline function Base.:+(x::DyadicRational{T}, y::DyadicRational{T}) where {T<:Integer}
    x.num == 0 && return y
    y.num == 0 && return x
    if x.exp == y.exp
        DyadicRational{T}(x.num + y.num, x.exp)
    elseif x.exp < y.exp
        DyadicRational{T}((x.num << (y.exp - x.exp)) + y.num, y.exp)
    else
        DyadicRational{T}(x.num + (y.num << (x.exp - y.exp)), x.exp)
    end
end

@inline function Base.:-(x::DyadicRational{T}, y::DyadicRational{T}) where {T<:Integer}
    x.num == 0 && return DyadicRational{T}(-y.num, y.exp)
    y.num == 0 && return x
    if x.exp == y.exp
        DyadicRational{T}(x.num - y.num, x.exp)
    elseif x.exp < y.exp
        DyadicRational{T}((x.num << (y.exp - x.exp)) - y.num, y.exp)
    else
        DyadicRational{T}(x.num - (y.num << (x.exp - y.exp)), x.exp)
    end
end

@inline Base.:*(x::DyadicRational{T}, y::DyadicRational{T}) where {T<:Integer} =
    DyadicRational{T}(x.num * y.num, x.exp + y.exp)

############################
# Dyadic division
############################

function dyadic_div(x::DyadicRational{T}, y::DyadicRational{T}) where {T<:Integer}
    y.num == 0 && throw(DivideError())
    x.num == 0 && return zero(x)

    b = y.num
    tb = trailing_zeros(b)
    ob = b >> tb

    rem(x.num, ob) == 0 || throw(ArgumentError("exact quotient not dyadic"))

    qnum = div(x.num, ob)
    qexp = x.exp - y.exp - tb

    qexp >= 0 ? DyadicRational{T}(qnum, qexp) :
    DyadicRational{T}(qnum << (-qexp), 0)
end

Base.:/(x::DyadicRational, y::DyadicRational) = to_fastrational(x) / to_fastrational(y)

############################
# Conversions
############################

isdyadic(x::FastRational) = (x.den & (x.den - 1)) == 0
isdyadic(::DyadicRational) = true

function to_fastrational(x::DyadicRational)
    FastRational(big(x.num), big(1) << x.exp)
end

function to_dyadic(x::FastRational)
    isdyadic(x) || throw(ArgumentError("not dyadic"))
    DyadicRational(big(x.num), trailing_zeros(x.den))
end

############################
# Self tests
############################

function run_selftests()

    using Test

    @test FastRational(6, 8) == FastRational(3, 4)
    @test FastRational(3, 4) + FastRational(5, 6) == FastRational(19, 12)
    @test FastRational(3, 4) * FastRational(5, 6) == FastRational(5, 8)

    x = DyadicRational(3, 2)
    y = DyadicRational(5, 3)

    @test x + y == DyadicRational(11, 3)
    @test x * y == DyadicRational(15, 5)

    println("All tests passed")

end

############################
# Benchmarks
############################

function run_benchmarks(N=100_000)

    using Random

    xs = [FastRational(rand(-100:100), rand(1:100)) for i in 1:N]
    ys = [FastRational(rand(-100:100), rand(1:100)) for i in 1:N]

    t = @elapsed begin
        s = FastRational(0, 1)
        for i in 1:N
            s += xs[i] + ys[i]
        end
    end

    println("FastRational add loop: ", t, " seconds")

    dx = [DyadicRational(rand(-100:100), rand(0:8)) for i in 1:N]
    dy = [DyadicRational(rand(-100:100), rand(0:8)) for i in 1:N]

    t = @elapsed begin
        s = DyadicRational(0, 0)
        for i in 1:N
            s += dx[i] + dy[i]
        end
    end

    println("Dyadic add loop: ", t, " seconds")

end

end # module


























module ExactRationals

"""
    ExactRationals

A compact, performance-oriented exact arithmetic module for Julia v1.12.

It provides two exact real-number types:

- `FastRational{T}` for general exact rationals `num/den` with canonical normalization.
- `DyadicRational{T}` for exact dyadics `num/2^exp`, where normalization is done
  entirely with bit operations, eliminating general `gcd` from dyadic addition,
  subtraction, and multiplication.

The module also includes:

- conversion helpers between dyadics and general rationals,
- self-tests via `run_selftests()`, and
- simple microbenchmarks via `run_benchmarks()`.

This file is intentionally self-contained so it can be dropped into a project,
`include`d directly, or used as a starting point for a package.
"""
export FastRational,
    DyadicRational,
    numerator,
    denominator,
    exponent,
    dyadic_denominator,
    isdyadic,
    to_fastrational,
    to_dyadic,
    dyadic_div,
    run_selftests,
    run_benchmarks

import Base: +, -, *, /, ==, <, <=, >, >=, isless, hash, show, zero, one,
    abs, sign, inv, convert, promote_rule, promote, widen, float

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

@inline _checked_nonzero_den(d) = d == 0 && throw(DivideError())

@inline function _normalize_sign(n::T, d::T) where {T<:Integer}
    d < 0 ? (-n, -d) : (n, d)
end

@inline function _trailing_zeros_nonzero(x::T) where {T<:Integer}
    x == 0 && throw(ArgumentError("trailing zeros are undefined for zero in this helper"))
    return trailing_zeros(x)
end

@inline function _is_power_of_two_abs(x::T) where {T<:Integer}
    x == 0 && return false
    ax = abs(x)
    return (ax & (ax - one(T))) == 0
end

# -----------------------------------------------------------------------------
# FastRational
# -----------------------------------------------------------------------------

"""
    FastRational{T} <: Real

Exact rational number in canonical form.

Invariants:
- `den > 0`
- `gcd(num, den) == 1`
- zero is represented as `0//1`

Arithmetic is optimized with:
- denominator-gcd reduction for addition/subtraction,
- cross-cancellation for multiplication/division.
"""
struct FastRational{T<:Integer} <: Real
    num::T
    den::T

    function FastRational{T}(n::T, d::T) where {T<:Integer}
        _checked_nonzero_den(d)
        if n == 0
            return new{T}(zero(T), one(T))
        end
        n, d = _normalize_sign(n, d)
        g = gcd(n, d)
        return new{T}(div(n, g), div(d, g))
    end
end

FastRational(n::T, d::T) where {T<:Integer} = FastRational{T}(n, d)
FastRational(n::Integer, d::Integer) = begin
    T = promote_type(typeof(n), typeof(d))
    FastRational(convert(T, n), convert(T, d))
end
FastRational(n::Integer) = FastRational(n, 1)

@inline numerator(x::FastRational) = x.num
@inline denominator(x::FastRational) = x.den

Base.show(io::IO, x::FastRational) = print(io, x.num, "//", x.den)

Base.zero(::Type{FastRational{T}}) where {T<:Integer} = FastRational{T}(zero(T), one(T))
Base.one(::Type{FastRational{T}}) where {T<:Integer} = FastRational{T}(one(T), one(T))
Base.zero(x::FastRational{T}) where {T<:Integer} = zero(FastRational{T})
Base.one(x::FastRational{T}) where {T<:Integer} = one(FastRational{T})

Base.abs(x::FastRational{T}) where {T<:Integer} = FastRational{T}(abs(x.num), x.den)
Base.sign(x::FastRational) = sign(x.num)
Base.:-(x::FastRational{T}) where {T<:Integer} = FastRational{T}(-x.num, x.den)

Base.:(==)(x::FastRational, y::FastRational) = x.num == y.num && x.den == y.den
Base.hash(x::FastRational, h::UInt) = hash(x.den, hash(x.num, h))

@inline function Base.isless(x::FastRational, y::FastRational)
    return x.num * y.den < y.num * x.den
end
Base.:<(x::FastRational, y::FastRational) = isless(x, y)
Base.:<=(x::FastRational, y::FastRational) = !isless(y, x)
Base.:>(x::FastRational, y::FastRational) = isless(y, x)
Base.:>=(x::FastRational, y::FastRational) = !isless(x, y)

function Base.:+(x::FastRational{T}, y::FastRational{T}) where {T<:Integer}
    x.num == 0 && return y
    y.num == 0 && return x

    g = gcd(x.den, y.den)
    xd = div(x.den, g)
    yd = div(y.den, g)

    n = x.num * yd + y.num * xd
    n == 0 && return FastRational{T}(zero(T), one(T))

    # denominator = lcm(x.den, y.den) = x.den * (y.den / g)
    d = x.den * yd

    # Common factor can only remain between n and g here.
    g2 = gcd(n, g)
    return FastRational{T}(div(n, g2), div(d, g2))
end

function Base.:-(x::FastRational{T}, y::FastRational{T}) where {T<:Integer}
    x.num == 0 && return FastRational{T}(-y.num, y.den)
    y.num == 0 && return x

    g = gcd(x.den, y.den)
    xd = div(x.den, g)
    yd = div(y.den, g)

    n = x.num * yd - y.num * xd
    n == 0 && return FastRational{T}(zero(T), one(T))

    d = x.den * yd   # = lcm(x.den, y.den)

    g2 = gcd(n, g)
    n = div(n, g2)
    d = div(d, g2)

    return FastRational{T}(n, d)
end

function Base.:*(x::FastRational{T}, y::FastRational{T}) where {T<:Integer}
    x.num == 0 && return zero(x)
    y.num == 0 && return zero(y)

    g1 = gcd(x.num, y.den)
    g2 = gcd(y.num, x.den)

    a = div(x.num, g1)
    d = div(y.den, g1)
    c = div(y.num, g2)
    b = div(x.den, g2)

    return FastRational{T}(a * c, b * d)
end

function Base.inv(x::FastRational{T}) where {T<:Integer}
    x.num == 0 && throw(DivideError())
    return x.num > 0 ? FastRational{T}(x.den, x.num) : FastRational{T}(-x.den, -x.num)
end

function Base.:/(x::FastRational{T}, y::FastRational{T}) where {T<:Integer}
    y.num == 0 && throw(DivideError())

    g1 = gcd(x.num, y.num)
    g2 = gcd(y.den, x.den)

    a = div(x.num, g1)
    c = div(y.num, g1)
    d = div(y.den, g2)
    b = div(x.den, g2)

    return FastRational{T}(a * d, b * c)
end

Base.promote_rule(::Type{FastRational{T}}, ::Type{S}) where {T<:Integer,S<:Integer} =
    FastRational{promote_type(T, S)}

Base.convert(::Type{FastRational{T}}, n::Integer) where {T<:Integer} =
    FastRational{T}(convert(T, n), one(T))

Base.convert(::Type{FastRational{T}}, x::FastRational{S}) where {T<:Integer,S<:Integer} =
    FastRational{T}(convert(T, x.num), convert(T, x.den))

# Integer interactions
for op in (:+, :-, :*, :/)
    @eval Base.$op(x::FastRational, n::Integer) = Base.$op(x, convert(typeof(x), n))
    @eval Base.$op(n::Integer, x::FastRational) = Base.$op(convert(typeof(x), n), x)
end

# -----------------------------------------------------------------------------
# DyadicRational
# -----------------------------------------------------------------------------

"""
    DyadicRational{T} <: Real

Exact dyadic rational `num / 2^exp` in canonical form.

Invariants:
- `exp ≥ 0`
- if `num == 0`, representation is `0 / 2^0`
- otherwise `num` is odd

Normalization uses only bit operations:
- strip common factors of 2 from `num`
- decrease `exp` by the same amount

This eliminates general `gcd` from dyadic addition, subtraction, and multiplication.
"""
struct DyadicRational{T<:Integer} <: Real
    num::T
    exp::Int

    function DyadicRational{T}(num::T, exp::Integer) where {T<:Integer}
        exp < 0 && throw(ArgumentError("exp must be ≥ 0"))
        if num == 0
            return new{T}(zero(T), 0)
        end
        tz = _trailing_zeros_nonzero(num)
        s = min(Int(tz), Int(exp))
        return new{T}(num >> s, Int(exp) - s)
    end
end

DyadicRational(num::T, exp::Integer) where {T<:Integer} = DyadicRational{T}(num, exp)
DyadicRational(num::Integer) = DyadicRational(num, 0)

@inline numerator(x::DyadicRational) = x.num
@inline exponent(x::DyadicRational) = x.exp
@inline dyadic_denominator(x::DyadicRational) = big(1) << x.exp
@inline denominator(x::DyadicRational) = dyadic_denominator(x)

"""
    isdyadic(x::FastRational) -> Bool

Return `true` if the exact rational has denominator a power of two.
"""
@inline isdyadic(x::FastRational) = _is_power_of_two_abs(x.den)
@inline isdyadic(::DyadicRational) = true

Base.show(io::IO, x::DyadicRational) =
    x.exp == 0 ? print(io, x.num) : print(io, x.num, "//2^", x.exp)

Base.zero(::Type{DyadicRational{T}}) where {T<:Integer} = DyadicRational{T}(zero(T), 0)
Base.one(::Type{DyadicRational{T}}) where {T<:Integer} = DyadicRational{T}(one(T), 0)
Base.zero(x::DyadicRational{T}) where {T<:Integer} = zero(DyadicRational{T})
Base.one(x::DyadicRational{T}) where {T<:Integer} = one(DyadicRational{T})

Base.abs(x::DyadicRational{T}) where {T<:Integer} = DyadicRational{T}(abs(x.num), x.exp)
Base.sign(x::DyadicRational) = sign(x.num)
Base.:-(x::DyadicRational{T}) where {T<:Integer} = DyadicRational{T}(-x.num, x.exp)

Base.:(==)(x::DyadicRational, y::DyadicRational) = x.num == y.num && x.exp == y.exp
Base.hash(x::DyadicRational, h::UInt) = hash(x.exp, hash(x.num, h))

function Base.isless(x::DyadicRational{T}, y::DyadicRational{T}) where {T<:Integer}
    if x.exp == y.exp
        return x.num < y.num
    elseif x.exp < y.exp
        return (x.num << (y.exp - x.exp)) < y.num
    else
        return x.num < (y.num << (x.exp - y.exp))
    end
end
Base.:<(x::DyadicRational, y::DyadicRational) = isless(x, y)
Base.:<=(x::DyadicRational, y::DyadicRational) = !isless(y, x)
Base.:>(x::DyadicRational, y::DyadicRational) = isless(y, x)
Base.:>=(x::DyadicRational, y::DyadicRational) = !isless(x, y)

function Base.:+(x::DyadicRational{T}, y::DyadicRational{T}) where {T<:Integer}
    x.num == 0 && return y
    y.num == 0 && return x

    if x.exp == y.exp
        return DyadicRational{T}(x.num + y.num, x.exp)
    elseif x.exp < y.exp
        return DyadicRational{T}((x.num << (y.exp - x.exp)) + y.num, y.exp)
    else
        return DyadicRational{T}(x.num + (y.num << (x.exp - y.exp)), x.exp)
    end
end

Base.:-(x::DyadicRational{T}, y::DyadicRational{T}) where {T<:Integer} = x + (-y)

Base.:*(x::DyadicRational{T}, y::DyadicRational{T}) where {T<:Integer} =
    DyadicRational{T}(x.num * y.num, x.exp + y.exp)

"""
    dyadic_div(x, y) -> DyadicRational

Exact division within the dyadic domain. Throws if the exact quotient is not dyadic.

This stays `gcd`-free in the common binary-only case and is useful when exact closure
can be established structurally.
"""
function dyadic_div(x::DyadicRational{T}, y::DyadicRational{T}) where {T<:Integer}
    y.num == 0 && throw(DivideError())

    # x / y = (x.num / y.num) * 2^(y.exp - x.exp)
    # The quotient is dyadic iff the odd part of y.num divides x.num.
    a = x.num
    b = y.num

    # Remove all powers of two from b. They can be absorbed into the exponent shift.
    tb = _trailing_zeros_nonzero(b)
    oddb = b >> tb
    extra_shift = Int(tb) + y.exp - x.exp

    # Exact dyadic quotient requires oddb | a.
    rem(a, oddb) == 0 || throw(ArgumentError("exact quotient is not dyadic"))
    qnum = div(a, oddb)

    if extra_shift >= 0
        return DyadicRational{T}(qnum, extra_shift)
    else
        return DyadicRational{T}(qnum << (-extra_shift), 0)
    end
end

"""
    to_fastrational(x::DyadicRational) -> FastRational

Convert a dyadic into the corresponding general rational.
"""
function to_fastrational(x::DyadicRational{T}) where {T<:Integer}
    return FastRational(big(x.num), big(1) << x.exp)
end

"""
    to_dyadic(x::FastRational) -> DyadicRational

Convert a general rational to a dyadic when possible. Throws if the denominator is
not a power of two.
"""
function to_dyadic(x::FastRational)
    isdyadic(x) || throw(ArgumentError("rational is not dyadic"))
    exp = trailing_zeros(x.den)
    return DyadicRational(big(x.num), exp)
end

# Dyadic division by general `/` returns a general exact rational, since the dyadic
# domain is not closed under division.
Base.:/(x::DyadicRational, y::DyadicRational) = to_fastrational(x) / to_fastrational(y)

# Dyadic inverse returns a general exact rational for semantic correctness.
Base.inv(x::DyadicRational) = inv(to_fastrational(x))

Base.promote_rule(::Type{DyadicRational{T}}, ::Type{S}) where {T<:Integer,S<:Integer} =
    DyadicRational{promote_type(T, S)}

Base.convert(::Type{DyadicRational{T}}, n::Integer) where {T<:Integer} =
    DyadicRational{T}(convert(T, n), 0)

Base.convert(::Type{DyadicRational{T}}, x::DyadicRational{S}) where {T<:Integer,S<:Integer} =
    DyadicRational{T}(convert(T, x.num), x.exp)

for op in (:+, :-, :*)
    @eval Base.$op(x::DyadicRational, n::Integer) = Base.$op(x, convert(typeof(x), n))
    @eval Base.$op(n::Integer, x::DyadicRational) = Base.$op(convert(typeof(x), n), x)
end

# -----------------------------------------------------------------------------
# Float conversion helpers
# -----------------------------------------------------------------------------

"""
    Base.float(x::FastRational)
    Base.float(x::DyadicRational)

Convert to `Float64` for convenience.
"""
Base.float(x::FastRational) = float(x.num) / float(x.den)
Base.float(x::DyadicRational) = ldexp(float(x.num), -x.exp)

# -----------------------------------------------------------------------------
# Self-tests
# -----------------------------------------------------------------------------

"""
    run_selftests(; verbose=true)

Run a collection of unit tests for both exact types. Returns `true` on success.
"""
function run_selftests(; verbose::Bool=true)
    using Test

    verbose && println("Running ExactRationals self-tests ...")

    @testset "FastRational normalization" begin
        @test FastRational(6, 8) == FastRational(3, 4)
        @test FastRational(-6, -8) == FastRational(3, 4)
        @test FastRational(6, -8) == FastRational(-3, 4)
        @test FastRational(0, 17) == FastRational(0, 1)
    end

    @testset "FastRational arithmetic" begin
        a = FastRational(3, 4)
        b = FastRational(5, 6)
        @test a + b == FastRational(19, 12)
        @test a - b == FastRational(-1, 12)
        @test a * b == FastRational(5, 8)
        @test a / b == FastRational(9, 10)
        @test a + 1 == FastRational(7, 4)
        @test 2 * a == FastRational(3, 2)
    end

    @testset "Dyadic normalization" begin
        @test DyadicRational(12, 5) == DyadicRational(3, 3)
        @test DyadicRational(0, 99) == DyadicRational(0, 0)
        @test DyadicRational(-40, 4) == DyadicRational(-5, 1)
    end

    @testset "Dyadic arithmetic" begin
        x = DyadicRational(3, 2)   # 3/4
        y = DyadicRational(5, 3)   # 5/8
        @test x + y == DyadicRational(11, 3)
        @test x - y == DyadicRational(1, 3)
        @test x * y == DyadicRational(15, 5)
        @test x + 1 == DyadicRational(7, 2)
        @test 2 * x == DyadicRational(3, 1)
    end

    @testset "Dyadic conversions" begin
        x = DyadicRational(7, 4)
        r = to_fastrational(x)
        @test r == FastRational(big(7), big(16))
        @test to_dyadic(FastRational(7, 16)) == DyadicRational(big(7), 4)
        @test isdyadic(FastRational(3, 8))
        @test !isdyadic(FastRational(3, 10))
        @test_throws ArgumentError to_dyadic(FastRational(3, 10))
    end

    @testset "Dyadic exact division in the dyadic domain" begin
        x = DyadicRational(3, 4)   # 3/16
        y = DyadicRational(3, 6)   # 3/64
        @test dyadic_div(x, y) == DyadicRational(4, 0)
        @test_throws ArgumentError dyadic_div(DyadicRational(3, 3), DyadicRational(5, 2))
    end

    @testset "Cross-check against Base Rational" begin
        for a in -7:7, b in 1:7, c in -7:7, d in 1:7
            x = FastRational(a, b)
            y = FastRational(c, d)
            xr = a // b
            yr = c // d
            @test x + y == FastRational(numerator(xr + yr), denominator(xr + yr))
            @test x - y == FastRational(numerator(xr - yr), denominator(xr - yr))
            @test x * y == FastRational(numerator(xr * yr), denominator(xr * yr))
            c != 0 && @test x / y == FastRational(numerator(xr / yr), denominator(xr / yr))
        end
    end

    verbose && println("All self-tests passed.")
    return true
end

# -----------------------------------------------------------------------------
# Benchmarks
# -----------------------------------------------------------------------------

"""
    run_benchmarks(; N=100_000, T=Int64, verbose=true)

Run simple microbenchmarks comparing:
- Base `Rational{T}`
- `FastRational{T}`
- `DyadicRational{T}`

For richer results, install `BenchmarkTools.jl`.
Returns a named tuple of measured timings when possible.
"""
function run_benchmarks(; N::Int=100_000, T::Type{<:Integer}=Int64, verbose::Bool=true)
    verbose && println("Running ExactRationals benchmarks with N = $N ...")

    # Build reproducible, nontrivial workloads.
    rsA = [mod(37i + 11, 257) - 128 for i in 1:N]
    rsB = [mod(61i + 17, 251) + 1 for i in 1:N]
    rsC = [mod(29i + 7, 241) - 120 for i in 1:N]
    rsD = [mod(47i + 13, 239) + 1 for i in 1:N]

    base_x = Rational{T}[T(rsA[i]) // T(rsB[i]) for i in 1:N]
    base_y = Rational{T}[T(rsC[i]) // T(rsD[i]) for i in 1:N]

    fast_x = FastRational{T}[FastRational(T(rsA[i]), T(rsB[i])) for i in 1:N]
    fast_y = FastRational{T}[FastRational(T(rsC[i]), T(rsD[i])) for i in 1:N]

    dy_x = DyadicRational{T}[DyadicRational(T(2 * rsA[i] + 1), mod(i, 12)) for i in 1:N]
    dy_y = DyadicRational{T}[DyadicRational(T(2 * rsC[i] + 1), mod(i + 5, 12)) for i in 1:N]

    add_base() = begin
        s = zero(Rational{T})
        @inbounds for i in eachindex(base_x)
            s += base_x[i] + base_y[i]
        end
        s
    end

    add_fast() = begin
        s = zero(FastRational{T})
        @inbounds for i in eachindex(fast_x)
            s += fast_x[i] + fast_y[i]
        end
        s
    end

    add_dy() = begin
        s = zero(DyadicRational{T})
        @inbounds for i in eachindex(dy_x)
            s += dy_x[i] + dy_y[i]
        end
        s
    end

    mul_base() = begin
        s = one(Rational{T})
        @inbounds for i in 1:min(N, 20_000)
            s *= base_x[i] * base_y[i]
        end
        s
    end

    mul_fast() = begin
        s = one(FastRational{T})
        @inbounds for i in 1:min(N, 20_000)
            s *= fast_x[i] * fast_y[i]
        end
        s
    end

    mul_dy() = begin
        s = one(DyadicRational{T})
        @inbounds for i in 1:min(N, 20_000)
            s *= dy_x[i] * dy_y[i]
        end
        s
    end

    have_bt = false
    try
        @eval import BenchmarkTools
        have_bt = true
    catch
        verbose && println("BenchmarkTools not found; using Base.@elapsed fallback timings.")
    end

    if have_bt
        bt = BenchmarkTools
        t_add_base = minimum(bt.@benchmark add_base()).time
        t_add_fast = minimum(bt.@benchmark add_fast()).time
        t_add_dy = minimum(bt.@benchmark add_dy()).time
        t_mul_base = minimum(bt.@benchmark mul_base()).time
        t_mul_fast = minimum(bt.@benchmark mul_fast()).time
        t_mul_dy = minimum(bt.@benchmark mul_dy()).time
    else
        t_add_base = @elapsed add_base()
        t_add_fast = @elapsed add_fast()
        t_add_dy = @elapsed add_dy()
        t_mul_base = @elapsed mul_base()
        t_mul_fast = @elapsed mul_fast()
        t_mul_dy = @elapsed mul_dy()
    end

    results = (
        add=(base=t_add_base, fast=t_add_fast, dyadic=t_add_dy),
        mul=(base=t_mul_base, fast=t_mul_fast, dyadic=t_mul_dy),
    )

    if verbose
        unit = have_bt ? "ns (minimum)" : "s (@elapsed)"
        println("Benchmark summary [", unit, "]")
        println("  add: Base Rational = ", results.add.base,
            ", FastRational = ", results.add.fast,
            ", DyadicRational = ", results.add.dyadic)
        println("  mul: Base Rational = ", results.mul.base,
            ", FastRational = ", results.mul.fast,
            ", DyadicRational = ", results.mul.dyadic)
    end

    return results
end

# -----------------------------------------------------------------------------
# Script entry point
# -----------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    println("ExactRationals self-contained module")
    run_selftests()
    println()
    println("Example benchmark run:")
    run_benchmarks(N=20_000)
end

end # module
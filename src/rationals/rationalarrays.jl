module RationalArrays

export FastRational,
    RationalSoA,
    DyadicSoA,
    numerator,
    denominator,
    exponent,
    push!,
    normalize!,
    add!,
    sub!,
    mul!,
    div!,
    dyadic_add!,
    dyadic_sub!,
    dyadic_mul!,
    dyadic_normalize!,
    sum_rationals,
    run_selftests,
    run_benchmarks

import Base: getindex, setindex!, length, size, IndexStyle, push!, show

# =============================================================================
# Scalar exact rational
# =============================================================================

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

# =============================================================================
# Structure-of-arrays containers
# =============================================================================

struct RationalSoA{T<:Integer}
    nums::Vector{T}
    dens::Vector{T}
    function RationalSoA{T}(nums::Vector{T}, dens::Vector{T}) where {T<:Integer}
        length(nums) == length(dens) || throw(ArgumentError("length mismatch"))
        new{T}(nums, dens)
    end
end

RationalSoA{T}() where {T<:Integer} = RationalSoA{T}(T[], T[])
RationalSoA(nums::Vector{T}, dens::Vector{T}) where {T<:Integer} = RationalSoA{T}(nums, dens)

struct DyadicSoA{T<:Integer}
    nums::Vector{T}
    exps::Vector{Int}
    function DyadicSoA{T}(nums::Vector{T}, exps::Vector{Int}) where {T<:Integer}
        length(nums) == length(exps) || throw(ArgumentError("length mismatch"))
        new{T}(nums, exps)
    end
end

DyadicSoA{T}() where {T<:Integer} = DyadicSoA{T}(T[], Int[])
DyadicSoA(nums::Vector{T}, exps::Vector{Int}) where {T<:Integer} = DyadicSoA{T}(nums, exps)

@inline length(A::RationalSoA) = length(A.nums)
@inline size(A::RationalSoA) = (length(A),)
IndexStyle(::Type{<:RationalSoA}) = IndexLinear()

@inline length(A::DyadicSoA) = length(A.nums)
@inline size(A::DyadicSoA) = (length(A),)
IndexStyle(::Type{<:DyadicSoA}) = IndexLinear()

function getindex(A::RationalSoA{T}, i::Int) where {T}
    @boundscheck checkbounds(A.nums, i)
    return FastRational{T}(A.nums[i], A.dens[i])
end

function setindex!(A::RationalSoA{T}, x::FastRational{T}, i::Int) where {T}
    @boundscheck checkbounds(A.nums, i)
    A.nums[i] = x.num
    A.dens[i] = x.den
    return A
end

function push!(A::RationalSoA{T}, x::FastRational{T}) where {T}
    push!(A.nums, x.num)
    push!(A.dens, x.den)
    return A
end

@inline numerator(A::RationalSoA) = A.nums
@inline denominator(A::RationalSoA) = A.dens

@inline exponent(A::DyadicSoA) = A.exps
@inline numerator(A::DyadicSoA) = A.nums

function getindex(A::DyadicSoA{T}, i::Int) where {T}
    @boundscheck checkbounds(A.nums, i)
    return (A.nums[i], A.exps[i])
end

function setindex!(A::DyadicSoA{T}, x::Tuple{T,Int}, i::Int) where {T}
    @boundscheck checkbounds(A.nums, i)
    A.nums[i] = x[1]
    A.exps[i] = x[2]
    return A
end

function push!(A::DyadicSoA{T}, x::Tuple{T,Int}) where {T}
    push!(A.nums, x[1])
    push!(A.exps, x[2])
    return A
end

show(io::IO, A::RationalSoA) = print(io, "RationalSoA(", length(A), " elements)")
show(io::IO, A::DyadicSoA) = print(io, "DyadicSoA(", length(A), " elements)")

# =============================================================================
# Scalar kernels over raw numerators/denominators
# =============================================================================

@inline function add_numden(an::T, ad::T, bn::T, bd::T) where {T<:Integer}
    g = gcd(ad, bd)
    ad1 = div(ad, g)
    bd1 = div(bd, g)
    n = an * bd1 + bn * ad1
    n == 0 && return (zero(T), one(T))
    d = ad * bd1
    g2 = gcd(n, g)
    return (div(n, g2), div(d, g2))
end

@inline function sub_numden(an::T, ad::T, bn::T, bd::T) where {T<:Integer}
    g = gcd(ad, bd)
    ad1 = div(ad, g)
    bd1 = div(bd, g)
    n = an * bd1 - bn * ad1
    n == 0 && return (zero(T), one(T))
    d = ad * bd1
    g2 = gcd(n, g)
    return (div(n, g2), div(d, g2))
end

@inline function mul_numden(an::T, ad::T, bn::T, bd::T) where {T<:Integer}
    an == 0 && return (zero(T), one(T))
    bn == 0 && return (zero(T), one(T))
    g1 = gcd(an, bd)
    g2 = gcd(bn, ad)
    a = div(an, g1)
    d = div(bd, g1)
    b = div(bn, g2)
    c = div(ad, g2)
    return (a * b, c * d)
end

@inline function div_numden(an::T, ad::T, bn::T, bd::T) where {T<:Integer}
    bn == 0 && throw(DivideError())
    an == 0 && return (zero(T), one(T))
    g1 = gcd(an, bn)
    g2 = gcd(bd, ad)
    a = div(an, g1)
    b = div(bn, g1)
    d = div(bd, g2)
    c = div(ad, g2)
    nn = a * d
    dd = c * b
    if dd < 0
        nn = -nn
        dd = -dd
    end
    return (nn, dd)
end

# =============================================================================
# Bulk rational kernels
# =============================================================================

function add!(C::RationalSoA{T}, A::RationalSoA{T}, B::RationalSoA{T}) where {T<:Integer}
    n = length(A)
    length(B) == n || throw(ArgumentError("length mismatch"))
    length(C) == n || throw(ArgumentError("length mismatch"))
    @inbounds for i in 1:n
        cn, cd = add_numden(A.nums[i], A.dens[i], B.nums[i], B.dens[i])
        C.nums[i] = cn
        C.dens[i] = cd
    end
    return C
end

function sub!(C::RationalSoA{T}, A::RationalSoA{T}, B::RationalSoA{T}) where {T<:Integer}
    n = length(A)
    length(B) == n || throw(ArgumentError("length mismatch"))
    length(C) == n || throw(ArgumentError("length mismatch"))
    @inbounds for i in 1:n
        cn, cd = sub_numden(A.nums[i], A.dens[i], B.nums[i], B.dens[i])
        C.nums[i] = cn
        C.dens[i] = cd
    end
    return C
end

function mul!(C::RationalSoA{T}, A::RationalSoA{T}, B::RationalSoA{T}) where {T<:Integer}
    n = length(A)
    length(B) == n || throw(ArgumentError("length mismatch"))
    length(C) == n || throw(ArgumentError("length mismatch"))
    @inbounds for i in 1:n
        cn, cd = mul_numden(A.nums[i], A.dens[i], B.nums[i], B.dens[i])
        C.nums[i] = cn
        C.dens[i] = cd
    end
    return C
end

function div!(C::RationalSoA{T}, A::RationalSoA{T}, B::RationalSoA{T}) where {T<:Integer}
    n = length(A)
    length(B) == n || throw(ArgumentError("length mismatch"))
    length(C) == n || throw(ArgumentError("length mismatch"))
    @inbounds for i in 1:n
        cn, cd = div_numden(A.nums[i], A.dens[i], B.nums[i], B.dens[i])
        C.nums[i] = cn
        C.dens[i] = cd
    end
    return C
end

function normalize!(A::RationalSoA{T}) where {T<:Integer}
    @inbounds for i in eachindex(A.nums, A.dens)
        n = A.nums[i]
        d = A.dens[i]
        d == 0 && throw(DivideError())
        if n == 0
            A.nums[i] = zero(T)
            A.dens[i] = one(T)
            continue
        end
        if d < 0
            n = -n
            d = -d
        end
        g = gcd(n, d)
        A.nums[i] = div(n, g)
        A.dens[i] = div(d, g)
    end
    return A
end

function sum_rationals(A::RationalSoA{T}) where {T<:Integer}
    sn = zero(T)
    sd = one(T)
    @inbounds for i in eachindex(A.nums, A.dens)
        sn, sd = add_numden(sn, sd, A.nums[i], A.dens[i])
    end
    return FastRational{T}(sn, sd)
end

# =============================================================================
# Dyadic SoA kernels
# =============================================================================

@inline function _dyadic_normalize_pair(n::T, e::Int) where {T<:Integer}
    e < 0 && throw(ArgumentError("dyadic exponent must be nonnegative"))
    if n == 0
        return (zero(T), 0)
    end
    s = min(trailing_zeros(n), e)
    return (n >> s, e - s)
end

function dyadic_add!(C::DyadicSoA{T}, A::DyadicSoA{T}, B::DyadicSoA{T}) where {T<:Integer}
    n = length(A)
    length(B) == n || throw(ArgumentError("length mismatch"))
    length(C) == n || throw(ArgumentError("length mismatch"))
    @inbounds for i in 1:n
        an, ae = A.nums[i], A.exps[i]
        bn, be = B.nums[i], B.exps[i]
        if ae == be
            cn, ce = an + bn, ae
        elseif ae < be
            cn, ce = (an << (be - ae)) + bn, be
        else
            cn, ce = an + (bn << (ae - be)), ae
        end
        C.nums[i], C.exps[i] = _dyadic_normalize_pair(cn, ce)
    end
    return C
end

function dyadic_sub!(C::DyadicSoA{T}, A::DyadicSoA{T}, B::DyadicSoA{T}) where {T<:Integer}
    n = length(A)
    length(B) == n || throw(ArgumentError("length mismatch"))
    length(C) == n || throw(ArgumentError("length mismatch"))
    @inbounds for i in 1:n
        an, ae = A.nums[i], A.exps[i]
        bn, be = B.nums[i], B.exps[i]
        if ae == be
            cn, ce = an - bn, ae
        elseif ae < be
            cn, ce = (an << (be - ae)) - bn, be
        else
            cn, ce = an - (bn << (ae - be)), ae
        end
        C.nums[i], C.exps[i] = _dyadic_normalize_pair(cn, ce)
    end
    return C
end

function dyadic_mul!(C::DyadicSoA{T}, A::DyadicSoA{T}, B::DyadicSoA{T}) where {T<:Integer}
    n = length(A)
    length(B) == n || throw(ArgumentError("length mismatch"))
    length(C) == n || throw(ArgumentError("length mismatch"))
    @inbounds for i in 1:n
        cn, ce = _dyadic_normalize_pair(A.nums[i] * B.nums[i], A.exps[i] + B.exps[i])
        C.nums[i] = cn
        C.exps[i] = ce
    end
    return C
end

function dyadic_normalize!(A::DyadicSoA{T}) where {T<:Integer}
    @inbounds for i in eachindex(A.nums, A.exps)
        A.nums[i], A.exps[i] = _dyadic_normalize_pair(A.nums[i], A.exps[i])
    end
    return A
end

# =============================================================================
# Self-tests
# =============================================================================

function run_selftests(; verbose::Bool=true)
    using Test

    verbose && println("Running RationalArrays self-tests ...")

    @testset "FastRational basics" begin
        @test FastRational(6, 8) == FastRational(3, 4)
        @test numerator(FastRational(-6, -8)) == 3
        @test denominator(FastRational(-6, -8)) == 4
    end

    @testset "RationalSoA indexing and push!" begin
        A = RationalSoA{Int}()
        push!(A, FastRational(3, 4))
        push!(A, FastRational(5, 6))
        @test length(A) == 2
        @test A[1] == FastRational(3, 4)
        @test A[2] == FastRational(5, 6)
    end

    @testset "RationalSoA bulk kernels" begin
        A = RationalSoA(Int[1, 3, 5], Int[2, 4, 6])
        B = RationalSoA(Int[1, 1, 1], Int[3, 5, 7])
        C = RationalSoA(zeros(Int, 3), ones(Int, 3))

        add!(C, A, B)
        @test C[1] == FastRational(5, 6)
        @test C[2] == FastRational(19, 20)
        @test C[3] == FastRational(41, 42)

        sub!(C, A, B)
        @test C[1] == FastRational(1, 6)
        @test C[2] == FastRational(11, 20)
        @test C[3] == FastRational(29, 42)

        mul!(C, A, B)
        @test C[1] == FastRational(1, 6)
        @test C[2] == FastRational(3, 20)
        @test C[3] == FastRational(5, 42)

        div!(C, A, B)
        @test C[1] == FastRational(3, 2)
        @test C[2] == FastRational(15, 4)
        @test C[3] == FastRational(35, 6)
    end

    @testset "RationalSoA normalization and sum" begin
        A = RationalSoA(Int[2, -6, 0], Int[4, -8, 99])
        normalize!(A)
        @test A[1] == FastRational(1, 2)
        @test A[2] == FastRational(3, 4)
        @test A[3] == FastRational(0, 1)

        B = RationalSoA(Int[1, 1, 1], Int[2, 3, 6])
        @test sum_rationals(B) == FastRational(1, 1)
    end

    @testset "DyadicSoA kernels" begin
        A = DyadicSoA(Int[3, 5, 7], Int[2, 3, 1])
        B = DyadicSoA(Int[1, 3, 1], Int[2, 2, 1])
        C = DyadicSoA(zeros(Int, 3), zeros(Int, 3))

        dyadic_add!(C, A, B)
        @test C[1] == (1, 0)     # 3/4 + 1/4 = 1
        @test C[2] == (11, 3)    # 5/8 + 3/4 = 11/8
        @test C[3] == (4, 0)     # 7/2 + 1/2 = 4

        dyadic_sub!(C, A, B)
        @test C[1] == (1, 1)     # 1/2
        @test C[2] == (-1, 3)    # -1/8
        @test C[3] == (3, 0)     # 3

        dyadic_mul!(C, A, B)
        @test C[1] == (3, 4)
        @test C[2] == (15, 4)
        @test C[3] == (7, 2)
    end

    verbose && println("All self-tests passed.")
    return true
end

# =============================================================================
# Microbenchmarks
# =============================================================================

function run_benchmarks(; N::Int=100_000, T::Type{<:Integer}=Int64, verbose::Bool=true)
    verbose && println("Running RationalArrays microbenchmarks with N = $N ...")

    as = T[mod(37 * i + 11, 251) - 125 for i in 1:N]
    bs = T[mod(61 * i + 17, 241) + 1 for i in 1:N]
    cs = T[mod(29 * i + 7, 239) - 119 for i in 1:N]
    ds = T[mod(47 * i + 13, 233) + 1 for i in 1:N]

    A = RationalSoA(copy(as), copy(bs))
    normalize!(A)
    B = RationalSoA(copy(cs), copy(ds))
    normalize!(B)
    C = RationalSoA(zeros(T, N), ones(T, N))

    dyA = DyadicSoA(T[2 * mod(i, 127) + 1 for i in 1:N], Int[mod(i, 8) for i in 1:N])
    dyB = DyadicSoA(T[2 * mod(3i + 1, 127) + 1 for i in 1:N], Int[mod(i + 3, 8) for i in 1:N])
    dyC = DyadicSoA(zeros(T, N), zeros(Int, N))

    t_add = @elapsed add!(C, A, B)
    t_sub = @elapsed sub!(C, A, B)
    t_mul = @elapsed mul!(C, A, B)
    t_div = @elapsed div!(C, A, B)

    t_dadd = @elapsed dyadic_add!(dyC, dyA, dyB)
    t_dsub = @elapsed dyadic_sub!(dyC, dyA, dyB)
    t_dmul = @elapsed dyadic_mul!(dyC, dyA, dyB)

    results = (
        rational=(add=t_add, sub=t_sub, mul=t_mul, div=t_div),
        dyadic=(add=t_dadd, sub=t_dsub, mul=t_dmul),
    )

    if verbose
        println("RationalSoA timings (seconds): add=$(t_add), sub=$(t_sub), mul=$(t_mul), div=$(t_div)")
        println("DyadicSoA   timings (seconds): add=$(t_dadd), sub=$(t_dsub), mul=$(t_dmul)")
    end

    return results
end

end # module
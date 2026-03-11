module FastMachine
# This file owns the machine-native data model and fast kernels.

export Dyadic,
    DyadicBuffer,
    DyadicMatrix,
    dyadic,
    numerator,
    exponent,
    normalize,
    normalize!,
    isnormalized,
    lazy_add!,
    lazy_sub!,
    lazy_mul!

import Base: +, -, *, /, ==, isless, zero, one, show, length, getindex,
    setindex!, push!, size, axes, IndexStyle, eltype, iterate,
    hash, iszero, abs, sign, convert, promote_rule, float

struct Dyadic{T<:Integer} <: Real
    num::T
    exp::Int
    function Dyadic{T}(num::T, exp::Integer) where {T<:Integer}
        exp < 0 && throw(ArgumentError("exp must be nonnegative"))
        if num == 0
            return new{T}(zero(T), 0)
        end
        s = min(trailing_zeros(num), Int(exp))
        new{T}(num >> s, Int(exp) - s)
    end
end

dyadic(num::T, exp::Integer=0) where {T<:Integer} = Dyadic{T}(num, exp)
dyadic(num::Integer, exp::Integer=0) = begin
    T = promote_type(typeof(num), Int)
    Dyadic{T}(convert(T, num), exp)
end

@inline numerator(x::Dyadic) = x.num
@inline exponent(x::Dyadic) = x.exp

Base.zero(::Type{Dyadic{T}}) where {T<:Integer} = Dyadic{T}(zero(T), 0)
Base.one(::Type{Dyadic{T}}) where {T<:Integer} = Dyadic{T}(one(T), 0)
Base.iszero(x::Dyadic) = iszero(x.num)
Base.abs(x::Dyadic{T}) where {T<:Integer} = Dyadic{T}(abs(x.num), x.exp)
Base.sign(x::Dyadic) = sign(x.num)
Base.hash(x::Dyadic, h::UInt) = hash(x.exp, hash(x.num, h))
Base.:(==)(x::Dyadic, y::Dyadic) = x.num == y.num && x.exp == y.exp

function Base.isless(x::Dyadic{T}, y::Dyadic{T}) where {T<:Integer}
    if x.exp == y.exp
        x.num < y.num
    elseif x.exp < y.exp
        (x.num << (y.exp - x.exp)) < y.num
    else
        x.num < (y.num << (x.exp - y.exp))
    end
end

function Base.:+(x::Dyadic{T}, y::Dyadic{T}) where {T<:Integer}
    if x.exp == y.exp
        Dyadic{T}(x.num + y.num, x.exp)
    elseif x.exp < y.exp
        Dyadic{T}((x.num << (y.exp - x.exp)) + y.num, y.exp)
    else
        Dyadic{T}(x.num + (y.num << (x.exp - y.exp)), x.exp)
    end
end

function Base.:-(x::Dyadic{T}, y::Dyadic{T}) where {T<:Integer}
    if x.exp == y.exp
        Dyadic{T}(x.num - y.num, x.exp)
    elseif x.exp < y.exp
        Dyadic{T}((x.num << (y.exp - x.exp)) - y.num, y.exp)
    else
        Dyadic{T}(x.num - (y.num << (x.exp - y.exp)), x.exp)
    end
end

Base.:-(x::Dyadic{T}) where {T<:Integer} = Dyadic{T}(-x.num, x.exp)
Base.:*(x::Dyadic{T}, y::Dyadic{T}) where {T<:Integer} = Dyadic{T}(x.num * y.num, x.exp + y.exp)

# Machine-only division policy:
# do NOT provide exact generic division here unless you want Rational{BigInt}.
# Either omit / or define a checked/approximate variant.
# Best practice: leave exact division to exact_promoted.jl.

Base.convert(::Type{Dyadic{T}}, n::Integer) where {T<:Integer} =
    Dyadic{T}(convert(T, n), 0)

Base.convert(::Type{Dyadic{T}}, x::Dyadic{S}) where {T<:Integer,S<:Integer} =
    Dyadic{T}(convert(T, x.num), x.exp)

Base.promote_rule(::Type{Dyadic{T}}, ::Type{S}) where {T<:Integer,S<:Integer} =
    Dyadic{promote_type(T, S)}

Base.float(x::Dyadic) = ldexp(float(x.num), -x.exp)

@inline function normalize(num::T, exp::Int) where {T<:Integer}
    exp < 0 && throw(ArgumentError("exp must be nonnegative"))
    if num == 0
        return (zero(T), 0)
    end
    s = min(trailing_zeros(num), exp)
    return (num >> s, exp - s)
end

mutable struct DyadicBuffer{T<:Integer}
    nums::Vector{T}
    exps::Vector{Int}
    normalized::Bool
    function DyadicBuffer{T}(nums::Vector{T}, exps::Vector{Int}, normalized::Bool=false) where {T<:Integer}
        length(nums) == length(exps) || throw(ArgumentError("length mismatch"))
        new{T}(nums, exps, normalized)
    end
end

DyadicBuffer{T}() where {T<:Integer} = DyadicBuffer{T}(T[], Int[], true)

Base.length(B::DyadicBuffer) = length(B.nums)
@inline isnormalized(B::DyadicBuffer) = B.normalized
Base.getindex(B::DyadicBuffer{T}, i::Int) where {T<:Integer} = Dyadic{T}(B.nums[i], B.exps[i])

function normalize!(B::DyadicBuffer{T}) where {T<:Integer}
    @inbounds for i in eachindex(B.nums, B.exps)
        B.nums[i], B.exps[i] = normalize(B.nums[i], B.exps[i])
    end
    B.normalized = true
    return B
end

function lazy_add!(C::DyadicBuffer{T}, A::DyadicBuffer{T}, B::DyadicBuffer{T}) where {T<:Integer}
    n = length(A)
    length(B) == n || throw(ArgumentError("length mismatch"))
    length(C.nums) == n || throw(ArgumentError("output length mismatch"))
    @inbounds for i in 1:n
        an, ae = A.nums[i], A.exps[i]
        bn, be = B.nums[i], B.exps[i]
        if ae == be
            C.nums[i] = an + bn
            C.exps[i] = ae
        elseif ae < be
            C.nums[i] = (an << (be - ae)) + bn
            C.exps[i] = be
        else
            C.nums[i] = an + (bn << (ae - be))
            C.exps[i] = ae
        end
    end
    C.normalized = false
    return C
end

function lazy_sub!(C::DyadicBuffer{T}, A::DyadicBuffer{T}, B::DyadicBuffer{T}) where {T<:Integer}
    n = length(A)
    length(B) == n || throw(ArgumentError("length mismatch"))
    length(C.nums) == n || throw(ArgumentError("output length mismatch"))
    @inbounds for i in 1:n
        an, ae = A.nums[i], A.exps[i]
        bn, be = B.nums[i], B.exps[i]
        if ae == be
            C.nums[i] = an - bn
            C.exps[i] = ae
        elseif ae < be
            C.nums[i] = (an << (be - ae)) - bn
            C.exps[i] = be
        else
            C.nums[i] = an - (bn << (ae - be))
            C.exps[i] = ae
        end
    end
    C.normalized = false
    return C
end

function lazy_mul!(C::DyadicBuffer{T}, A::DyadicBuffer{T}, B::DyadicBuffer{T}) where {T<:Integer}
    n = length(A)
    length(B) == n || throw(ArgumentError("length mismatch"))
    length(C.nums) == n || throw(ArgumentError("output length mismatch"))
    @inbounds for i in 1:n
        C.nums[i] = A.nums[i] * B.nums[i]
        C.exps[i] = A.exps[i] + B.exps[i]
    end
    C.normalized = false
    return C
end

struct DyadicMatrix{T<:Integer} <: AbstractMatrix{Dyadic{T}}
    nums::Matrix{T}
    exps::Matrix{Int}
    function DyadicMatrix{T}(nums::Matrix{T}, exps::Matrix{Int}) where {T<:Integer}
        size(nums) == size(exps) || throw(ArgumentError("size mismatch"))
        new{T}(nums, exps)
    end
end

function DyadicMatrix(nums::AbstractMatrix{T}, exps::AbstractMatrix{<:Integer}) where {T<:Integer}
    DyadicMatrix{T}(Matrix{T}(nums), Int.(exps))
end

Base.size(A::DyadicMatrix) = size(A.nums)
Base.axes(A::DyadicMatrix) = axes(A.nums)
Base.IndexStyle(::Type{<:DyadicMatrix}) = IndexCartesian()
Base.eltype(::Type{DyadicMatrix{T}}) where {T<:Integer} = Dyadic{T}
Base.getindex(A::DyadicMatrix{T}, i::Int, j::Int) where {T<:Integer} = Dyadic{T}(A.nums[i, j], A.exps[i, j])

function Base.setindex!(A::DyadicMatrix{T}, x::Dyadic{T}, i::Int, j::Int) where {T<:Integer}
    A.nums[i, j] = x.num
    A.exps[i, j] = x.exp
    return A
end

end # module
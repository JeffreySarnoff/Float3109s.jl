module DyadicCASEngine

export Dyadic,
    DyadicBuffer,
    DyadicMatrix,
    dyadic,
    numerator,
    exponent,
    value,
    normalize,
    normalize!,
    isnormalized,
    lazy_add!,
    lazy_sub!,
    lazy_mul!,
    pushnormalized!,
    common_denominator_form,
    bareiss_det,
    bareiss_solve,
    crt_pair,
    crt_combine,
    reconstruct_rational,
    reconstruct_vector,
    reconstruct_matrix,
    solve_mod_prime,
    solve_via_modular_reconstruction

import Base: +, -, *, /, ==, isless, zero, one, show, length, getindex, setindex!,
    push!, size, axes, IndexStyle, eltype, iterate, hash, iszero, abs,
    sign, convert, promote_rule, float

# =============================================================================
# Dyadic scalar
# =============================================================================

"""
    Dyadic{T<:Integer} <: Real

Exact dyadic rational `num / 2^exp` with `exp ≥ 0`.

Canonical invariants:
- `exp ≥ 0`
- `num == 0  =>  exp == 0`
- `num != 0  =>  num` is odd

This is a natural exact representation for binary-scaled quantities, exact decoded
floating-point values, and power-of-two denominator arithmetic.
"""
struct Dyadic{T<:Integer} <: Real
    num::T
    exp::Int

    function Dyadic{T}(num::T, exp::Integer) where {T<:Integer}
        exp < 0 && throw(ArgumentError("exp must be nonnegative"))
        if num == 0
            return new{T}(zero(T), 0)
        end
        s = min(trailing_zeros(num), Int(exp))
        return new{T}(num >> s, Int(exp) - s)
    end
end

"""
    dyadic(num, [exp=0])

Construct a canonical dyadic value.
"""
dyadic(num::T, exp::Integer=0) where {T<:Integer} = Dyadic{T}(num, exp)
dyadic(num::Integer, exp::Integer=0) = begin
    T = promote_type(typeof(num), Int)
    Dyadic{T}(convert(T, num), exp)
end

"""
    numerator(x::Dyadic)

Return the stored integer numerator.
"""
@inline numerator(x::Dyadic) = x.num

"""
    exponent(x::Dyadic)

Return the stored denominator exponent.
"""
@inline exponent(x::Dyadic) = x.exp

"""
    value(x::Dyadic) -> Rational{BigInt}

Return the exact mathematical value.
"""
value(x::Dyadic) = BigInt(x.num) // (BigInt(1) << x.exp)

Base.show(io::IO, x::Dyadic) =
    x.exp == 0 ? print(io, x.num) : print(io, x.num, "//2^", x.exp)

Base.zero(::Type{Dyadic{T}}) where {T<:Integer} = Dyadic{T}(zero(T), 0)
Base.one(::Type{Dyadic{T}}) where {T<:Integer} = Dyadic{T}(one(T), 0)
Base.zero(x::Dyadic{T}) where {T<:Integer} = zero(Dyadic{T})
Base.one(x::Dyadic{T}) where {T<:Integer} = one(Dyadic{T})

Base.iszero(x::Dyadic) = iszero(x.num)
Base.abs(x::Dyadic{T}) where {T<:Integer} = Dyadic{T}(abs(x.num), x.exp)
Base.sign(x::Dyadic) = sign(x.num)
Base.hash(x::Dyadic, h::UInt) = hash(x.exp, hash(x.num, h))
Base.:(==)(x::Dyadic, y::Dyadic) = x.num == y.num && x.exp == y.exp

function Base.isless(x::Dyadic{T}, y::Dyadic{T}) where {T<:Integer}
    if x.exp == y.exp
        return x.num < y.num
    elseif x.exp < y.exp
        return (x.num << (y.exp - x.exp)) < y.num
    else
        return x.num < (y.num << (x.exp - y.exp))
    end
end

function Base.:+(x::Dyadic{T}, y::Dyadic{T}) where {T<:Integer}
    if x.exp == y.exp
        return Dyadic{T}(x.num + y.num, x.exp)
    elseif x.exp < y.exp
        return Dyadic{T}((x.num << (y.exp - x.exp)) + y.num, y.exp)
    else
        return Dyadic{T}(x.num + (y.num << (x.exp - y.exp)), x.exp)
    end
end

function Base.:-(x::Dyadic{T}, y::Dyadic{T}) where {T<:Integer}
    if x.exp == y.exp
        return Dyadic{T}(x.num - y.num, x.exp)
    elseif x.exp < y.exp
        return Dyadic{T}((x.num << (y.exp - x.exp)) - y.num, y.exp)
    else
        return Dyadic{T}(x.num - (y.num << (x.exp - y.exp)), x.exp)
    end
end

Base.:-(x::Dyadic{T}) where {T<:Integer} = Dyadic{T}(-x.num, x.exp)
Base.:*(x::Dyadic{T}, y::Dyadic{T}) where {T<:Integer} = Dyadic{T}(x.num * y.num, x.exp + y.exp)

"""
    /(x::Dyadic, y::Dyadic)

Exact division. Since dyadics are not closed under general division, this returns
a Julia `Rational{BigInt}`.
"""
Base.:/(x::Dyadic{T}, y::Dyadic{T}) where {T<:Integer} = value(x) / value(y)

Base.float(x::Dyadic) = ldexp(float(x.num), -x.exp)

Base.convert(::Type{Dyadic{T}}, n::Integer) where {T<:Integer} =
    Dyadic{T}(convert(T, n), 0)

Base.convert(::Type{Dyadic{T}}, x::Dyadic{S}) where {T<:Integer,S<:Integer} =
    Dyadic{T}(convert(T, x.num), x.exp)

Base.promote_rule(::Type{Dyadic{T}}, ::Type{S}) where {T<:Integer,S<:Integer} =
    Dyadic{promote_type(T, S)}

# Integer mixed-mode convenience.
Base.:+(x::Dyadic, n::Integer) = x + convert(typeof(x), n)
Base.:+(n::Integer, x::Dyadic) = convert(typeof(x), n) + x
Base.:-(x::Dyadic, n::Integer) = x - convert(typeof(x), n)
Base.:-(n::Integer, x::Dyadic) = convert(typeof(x), n) - x
Base.:*(x::Dyadic, n::Integer) = x * convert(typeof(x), n)
Base.:*(n::Integer, x::Dyadic) = convert(typeof(x), n) * x

# =============================================================================
# Lazy-normalization dyadic buffer
# =============================================================================

"""
    DyadicBuffer{T}

Structure-of-arrays storage for many dyadic values.

Fields:
- `nums`: numerators
- `exps`: denominator exponents
- `normalized`: whether all entries are known canonical

This supports lazy normalization: fast bulk kernels write mathematically correct raw
dyadic pairs and postpone canonicalization.
"""
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

"""
    isnormalized(B)

Return whether `B` is currently marked normalized.
"""
@inline isnormalized(B::DyadicBuffer) = B.normalized

"""
    normalize(num, exp)

Canonicalize a raw dyadic pair.
"""
@inline function normalize(num::T, exp::Int) where {T<:Integer}
    exp < 0 && throw(ArgumentError("exp must be nonnegative"))
    if num == 0
        return (zero(T), 0)
    end
    s = min(trailing_zeros(num), exp)
    return (num >> s, exp - s)
end

"""
    normalize!(B)

Canonicalize every element in `B` in place.
"""
function normalize!(B::DyadicBuffer{T}) where {T<:Integer}
    @inbounds for i in eachindex(B.nums, B.exps)
        B.nums[i], B.exps[i] = normalize(B.nums[i], B.exps[i])
    end
    B.normalized = true
    return B
end

"""
    pushnormalized!(B, x)

Append a canonical dyadic entry.
"""
function pushnormalized!(B::DyadicBuffer{T}, x::Dyadic{T}) where {T<:Integer}
    push!(B.nums, x.num)
    push!(B.exps, x.exp)
    B.normalized &= true
    return B
end

Base.push!(B::DyadicBuffer{T}, x::Dyadic{T}) where {T<:Integer} = pushnormalized!(B, x)
Base.getindex(B::DyadicBuffer{T}, i::Int) where {T<:Integer} = Dyadic{T}(B.nums[i], B.exps[i])

"""
    lazy_add!(C, A, B)

Elementwise dyadic addition without per-element normalization.
"""
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

"""
    lazy_sub!(C, A, B)

Elementwise dyadic subtraction without per-element normalization.
"""
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

"""
    lazy_mul!(C, A, B)

Elementwise dyadic multiplication without per-element normalization.
"""
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

# =============================================================================
# DyadicMatrix with AbstractMatrix-style API
# =============================================================================

"""
    DyadicMatrix{T} <: AbstractMatrix{Dyadic{T}}

Matrix of exact dyadic values stored in structure-of-arrays form:
- `nums[i,j]`
- `exps[i,j]`

This is designed as a staging type for exact binary-scaled linear algebra.
"""
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

Base.getindex(A::DyadicMatrix{T}, i::Int, j::Int) where {T<:Integer} =
    Dyadic{T}(A.nums[i, j], A.exps[i, j])

function Base.setindex!(A::DyadicMatrix{T}, x::Dyadic{T}, i::Int, j::Int) where {T<:Integer}
    A.nums[i, j] = x.num
    A.exps[i, j] = x.exp
    return A
end

function Base.iterate(A::DyadicMatrix, state...)
    iter = iterate(CartesianIndices(A), state...)
    iter === nothing && return nothing
    I, st = iter
    return (A[I], st)
end

Base.getindex(A::DyadicMatrix, I::CartesianIndex) = A[I.I...]

function Base.show(io::IO, A::DyadicMatrix)
    m, n = size(A)
    println(io, "DyadicMatrix($m×$n)")
    for i in 1:m
        for j in 1:n
            print(io, A[i, j], j == n ? "" : "  ")
        end
        i < m && println(io)
    end
end

"""
    common_denominator_form(A)

Return `(M, E)` such that every entry of `A` is represented as `M[i,j] / 2^E`.
This converts a dyadic matrix into an integer matrix plus one shared scale factor.
"""
function common_denominator_form(A::DyadicMatrix{T}) where {T<:Integer}
    E = maximum(A.exps)
    M = Matrix{BigInt}(undef, size(A)...)
    @inbounds for i in axes(A.nums, 1), j in axes(A.nums, 2)
        M[i, j] = BigInt(A.nums[i, j]) << (E - A.exps[i, j])
    end
    return M, E
end

# =============================================================================
# Fraction-free Bareiss elimination
# =============================================================================

"""
    bareiss_det(M::AbstractMatrix{<:Integer}) -> BigInt

Compute the determinant using fraction-free Bareiss elimination.
"""
function bareiss_det(M::AbstractMatrix{<:Integer})
    n, m = size(M)
    n == m || throw(ArgumentError("matrix must be square"))
    A = Matrix{BigInt}(M)
    n == 0 && return BigInt(1)

    prev = BigInt(1)
    signflip = false

    for k in 1:n-1
        if A[k, k] == 0
            pivot_row = findfirst(r -> A[r, k] != 0, k+1:n)
            pivot_row === nothing && return BigInt(0)
            r = pivot_row
            A[k, :], A[r, :] = copy(A[r, :]), copy(A[k, :])
            signflip = !signflip
        end

        pivot = A[k, k]
        for i in k+1:n
            for j in k+1:n
                A[i, j] = div(pivot * A[i, j] - A[i, k] * A[k, j], prev)
            end
        end
        for i in k+1:n
            A[i, k] = 0
        end
        prev = pivot
    end

    return signflip ? -A[n, n] : A[n, n]
end

"""
    bareiss_det(A::DyadicMatrix) -> Rational{BigInt}

Determinant of a dyadic matrix by clearing a common power-of-two denominator and
applying Bareiss to the resulting integer matrix.
"""
function bareiss_det(A::DyadicMatrix)
    M, E = common_denominator_form(A)
    n = size(M, 1)
    d = bareiss_det(M)
    return d // (BigInt(1) << (E * n))
end

"""
    bareiss_solve(M, b) -> Vector{Rational{BigInt}}

Exact solution of `M*x = b` using fraction-free elimination followed by exact back
substitution.
"""
function bareiss_solve(M::AbstractMatrix{<:Integer}, b::AbstractVector{<:Integer})
    n, m = size(M)
    n == m || throw(ArgumentError("matrix must be square"))
    length(b) == n || throw(ArgumentError("dimension mismatch"))

    A = Matrix{BigInt}(undef, n, n + 1)
    for i in 1:n
        for j in 1:n
            A[i, j] = BigInt(M[i, j])
        end
        A[i, n+1] = BigInt(b[i])
    end

    prev = BigInt(1)
    signflip = false

    for k in 1:n-1
        if A[k, k] == 0
            pivot_row = findfirst(r -> A[r, k] != 0, k+1:n)
            pivot_row === nothing && throw(ArgumentError("singular matrix"))
            r = pivot_row
            A[k, :], A[r, :] = copy(A[r, :]), copy(A[k, :])
            signflip = !signflip
        end

        pivot = A[k, k]
        for i in k+1:n
            for j in k+1:n+1
                A[i, j] = div(pivot * A[i, j] - A[i, k] * A[k, j], prev)
            end
        end
        for i in k+1:n
            A[i, k] = 0
        end
        prev = pivot
    end

    detA = signflip ? -A[n, n] : A[n, n]
    detA == 0 && throw(ArgumentError("singular matrix"))

    x = Vector{Rational{BigInt}}(undef, n)
    for i in n:-1:1
        s = A[i, n+1] // BigInt(1)
        for j in i+1:n
            s -= A[i, j] * x[j]
        end
        x[i] = s / A[i, i]
    end
    return x
end

"""
    bareiss_solve(A::DyadicMatrix, bnums, bexps)

Solve a dyadic linear system exactly.
"""
function bareiss_solve(A::DyadicMatrix{T}, bnums::Vector{T}, bexps::Vector{Int}) where {T<:Integer}
    n, m = size(A)
    n == m || throw(ArgumentError("matrix must be square"))
    (length(bnums) == n && length(bexps) == n) || throw(ArgumentError("dimension mismatch"))

    M, EA = common_denominator_form(A)
    EB = maximum(bexps)
    b = Vector{BigInt}(undef, n)
    @inbounds for i in 1:n
        b[i] = BigInt(bnums[i]) << (EB - bexps[i])
    end

    if EA >= EB
        rhs = [bi << (EA - EB) for bi in b]
        return bareiss_solve(M, rhs)
    else
        x = bareiss_solve(M, b)
        return [xi / (BigInt(1) << (EB - EA)) for xi in x]
    end
end

# =============================================================================
# CRT and rational reconstruction
# =============================================================================

"""
    crt_pair(a, m, b, n)

Combine two congruences with coprime moduli.
"""
function crt_pair(a::Integer, m::Integer, b::Integer, n::Integer)
    gcd(m, n) == 1 || throw(ArgumentError("moduli must be coprime"))
    g, u, _ = gcdx(m, n)
    g == 1 || throw(ArgumentError("moduli must be coprime"))
    t = mod((b - a) * u, n)
    mn = BigInt(m) * BigInt(n)
    x = mod(BigInt(a) + BigInt(m) * t, mn)
    return x, mn
end

"""
    crt_combine(residues, moduli)

Combine many congruences by repeated CRT.
"""
function crt_combine(residues::AbstractVector{<:Integer}, moduli::AbstractVector{<:Integer})
    length(residues) == length(moduli) || throw(ArgumentError("length mismatch"))
    x = BigInt(mod(residues[1], moduli[1]))
    m = BigInt(moduli[1])
    for i in 2:length(residues)
        x, m = crt_pair(x, m, BigInt(mod(residues[i], moduli[i])), BigInt(moduli[i]))
    end
    return x, m
end

"""
    reconstruct_rational(x, m)

Reconstruct a small exact rational from a residue class.
"""
function reconstruct_rational(x::Integer, m::Integer)
    x = BigInt(mod(x, m))
    m = BigInt(m)

    r0, r1 = m, x
    t0, t1 = BigInt(0), BigInt(1)
    bound = isqrt(m ÷ 2)

    while abs(r1) > bound
        q = div(r0, r1)
        r0, r1 = r1, r0 - q * r1
        t0, t1 = t1, t0 - q * t1
    end

    t1 == 0 && throw(ArgumentError("reconstruction failed"))
    p, q = r1, t1
    if q < 0
        p = -p
        q = -q
    end
    g = gcd(p, q)
    p = div(p, g)
    q = div(q, g)
    (abs(p) <= bound && abs(q) <= bound) || throw(ArgumentError("reconstruction failed"))
    return p // q
end

reconstruct_vector(xs::AbstractVector{<:Integer}, m::Integer) =
    [reconstruct_rational(x, m) for x in xs]

function reconstruct_matrix(X::AbstractMatrix{<:Integer}, m::Integer)
    R = Matrix{Rational{BigInt}}(undef, size(X)...)
    for i in axes(X, 1), j in axes(X, 2)
        R[i, j] = reconstruct_rational(X[i, j], m)
    end
    return R
end

# =============================================================================
# Modular solve
# =============================================================================

"""
    solve_mod_prime(M, b, p)

Solve `M*x = b (mod p)` over the prime field `ℤ/pℤ`.
"""
function solve_mod_prime(M::AbstractMatrix{<:Integer}, b::AbstractVector{<:Integer}, p::Integer)
    n, m = size(M)
    n == m || throw(ArgumentError("matrix must be square"))
    length(b) == n || throw(ArgumentError("dimension mismatch"))

    A = Matrix{Int}(undef, n, n + 1)
    for i in 1:n
        for j in 1:n
            A[i, j] = mod(M[i, j], p)
        end
        A[i, n+1] = mod(b[i], p)
    end

    row = 1
    for col in 1:n
        pivot = 0
        for r in row:n
            if A[r, col] % p != 0
                pivot = r
                break
            end
        end
        pivot == 0 && throw(ArgumentError("singular modulo $p"))

        if pivot != row
            A[row, :], A[pivot, :] = copy(A[pivot, :]), copy(A[row, :])
        end

        g, invp, _ = gcdx(A[row, col], p)
        g == 1 || throw(ArgumentError("pivot not invertible modulo $p"))
        invp = mod(invp, p)

        for j in col:n+1
            A[row, j] = mod(A[row, j] * invp, p)
        end

        for r in 1:n
            r == row && continue
            f = A[r, col]
            if f != 0
                for j in col:n+1
                    A[r, j] = mod(A[r, j] - f * A[row, j], p)
                end
            end
        end

        row += 1
    end

    x = Vector{Int}(undef, n)
    for i in 1:n
        x[i] = mod(A[i, n+1], p)
    end
    return x
end

"""
    solve_via_modular_reconstruction(M, b; primes=...)

Solve an exact rational linear system by modular solves, CRT, and rational
reconstruction.
"""
function solve_via_modular_reconstruction(
    M::AbstractMatrix{<:Integer},
    b::AbstractVector{<:Integer};
    primes=[101, 103, 107, 109, 113, 127, 131],
)
    n, m = size(M)
    n == m || throw(ArgumentError("matrix must be square"))
    length(b) == n || throw(ArgumentError("dimension mismatch"))

    residues = [Int[] for _ in 1:n]
    mods = Int[]

    for p in primes
        xmod = solve_mod_prime(M, b, p)
        push!(mods, p)
        for i in 1:n
            push!(residues[i], xmod[i])
        end
    end

    xcrt = Vector{BigInt}(undef, n)
    modulus = BigInt(1)
    for i in 1:n
        xcrt[i], modulus = crt_combine(residues[i], mods)
    end

    return reconstruct_vector(xcrt, modulus)
end

end # module

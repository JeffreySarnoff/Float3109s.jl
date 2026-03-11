module ExactPromoted

using ..FastMachine: Dyadic, DyadicMatrix, numerator, exponent

export exact_value,
    exact_div,
    common_denominator_form_big,
    bareiss_det_exact,
    bareiss_solve_exact,
    crt_pair_big,
    crt_combine_big,
    reconstruct_rational,
    reconstruct_vector,
    reconstruct_matrix,
    solve_mod_prime,
    solve_via_modular_reconstruction

"""
    exact_value(x::Dyadic) -> Rational{BigInt}

Return the exact value of a dyadic as a `Rational{BigInt}`.

Promotion rule:
- always promote numerator and denominator to `BigInt`
- because `2^exp` may not fit in the input integer type
"""
exact_value(x::Dyadic) = BigInt(numerator(x)) // (BigInt(1) << exponent(x))

"""
    exact_div(x::Dyadic, y::Dyadic) -> Rational{BigInt}

Exact division of dyadics. This promotes because dyadics are not closed under
general division.
"""
exact_div(x::Dyadic, y::Dyadic) = exact_value(x) / exact_value(y)

"""
    common_denominator_form_big(A::DyadicMatrix) -> (Matrix{BigInt}, Int)

Convert `A` to `(M, E)` such that `A[i,j] = M[i,j] / 2^E`.

Promotion rule:
- always promote lifted matrix entries to `BigInt`
- because shifts may exceed machine precision
"""
function common_denominator_form_big(A::DyadicMatrix{T}) where {T<:Integer}
    E = maximum(A.exps)
    M = Matrix{BigInt}(undef, size(A)...)
    @inbounds for i in axes(A.nums, 1), j in axes(A.nums, 2)
        M[i, j] = BigInt(A.nums[i, j]) << (E - A.exps[i, j])
    end
    return M, E
end

"""
    bareiss_det_exact(M::AbstractMatrix{<:Integer}) -> BigInt

Exact determinant via fraction-free Bareiss. Always promotes to `BigInt`.
"""
function bareiss_det_exact(M::AbstractMatrix{<:Integer})
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
    bareiss_det_exact(A::DyadicMatrix) -> Rational{BigInt}

Exact determinant of a dyadic matrix. Promotes through `BigInt`.
"""
function bareiss_det_exact(A::DyadicMatrix)
    M, E = common_denominator_form_big(A)
    n = size(M, 1)
    d = bareiss_det_exact(M)
    return d // (BigInt(1) << (E * n))
end

"""
    bareiss_solve_exact(M, b) -> Vector{Rational{BigInt}}

Exact linear solve by fraction-free elimination. Promotes to `BigInt`.
"""
function bareiss_solve_exact(M::AbstractMatrix{<:Integer}, b::AbstractVector{<:Integer})
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

function crt_pair_big(a::Integer, m::Integer, b::Integer, n::Integer)
    gcd(m, n) == 1 || throw(ArgumentError("moduli must be coprime"))
    g, u, _ = gcdx(m, n)
    g == 1 || throw(ArgumentError("moduli must be coprime"))
    t = mod((b - a) * u, n)
    mn = BigInt(m) * BigInt(n)
    x = mod(BigInt(a) + BigInt(m) * t, mn)
    return x, mn
end

function crt_combine_big(residues::AbstractVector{<:Integer}, moduli::AbstractVector{<:Integer})
    length(residues) == length(moduli) || throw(ArgumentError("length mismatch"))
    x = BigInt(mod(residues[1], moduli[1]))
    m = BigInt(moduli[1])
    for i in 2:length(residues)
        x, m = crt_pair_big(x, m, BigInt(mod(residues[i], moduli[i])), BigInt(moduli[i]))
    end
    return x, m
end

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

function solve_via_modular_reconstruction(M::AbstractMatrix{<:Integer}, b::AbstractVector{<:Integer};
    primes=[101, 103, 107, 109, 113, 127, 131])
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
        xcrt[i], modulus = crt_combine_big(residues[i], mods)
    end

    return reconstruct_vector(xcrt, modulus)
end

end
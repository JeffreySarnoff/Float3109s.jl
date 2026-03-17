if !isdefined(@__MODULE__, :FormatTrait)
    include("traits.jl")
end

setprecision(BigFloat, 1024)
const DyadicSignificandInteger = Int32
const DyadicExponentInteger = Int32

const CPOINT = UInt64
const FPVALUE = BigFloat

# =========================================================================
# Formats
# =========================================================================

abstract type FloatFormat <: AbstractFloat end

"""
    Format{S<:Signedness, D<:Domain}(K, P)

A P3109 floating-point format with bitwidth `K` and precision `P`.

The type parameters encode signedness (`is_signed` / `is_unsigned`) and
domain (`is_finite` / `is_extended`) at the type level.
"""
struct Format{S<:Signedness,D<:Domain} <: FloatFormat
    K::Int
    P::Int
end

function Format(K::Integer, P::Integer, S::Type{Sgn}, D::Type{Dom}) where {Sgn<:Signedness,Dom<:Domain}
    Format{S,D}(K, P)
end

"""Return `true` if `fmt` is unsigned."""
is_unsigned(fnt::Format{UnsignedFormat,T}) where T = true
is_unsigned(fnt::Format{SignedFormat,T}) where T = false

"""Return `true` if `fmt` is signed."""
is_signed(fnt::Format{UnsignedFormat,T}) where T = false
is_signed(fnt::Format{SignedFormat,T}) where T = true

"""Return `true` if `fmt` is finite-only (no infinities)."""
is_finite(fnt::Format{T,FiniteFormat}) where T = true
is_finite(fnt::Format{T,ExtendedFormat}) where T = false

"""Return `true` if `fmt` is extended (has ±Inf)."""
is_extended(fnt::Format{T,FiniteFormat}) where T = false
is_extended(fnt::Format{T,ExtendedFormat}) where T = true

is_unsigned_finite(fmt::Format{UnsignedFormat,FiniteFormat}) = true
is_unsigned_extended(fmt::Format{UnsignedFormat,ExtendedFormat}) = true
is_signed_finite(fmt::Format{SignedFormat,FiniteFormat}) = true
is_signed_extended(fmt::Format{SignedFormat,ExtendedFormat}) = true

is_unsigned_finite(@nospecialize fmt::Format{<:Signedness,<:Domain}) = false
is_unsigned_extended(@nospecialize fmt::Format{<:Signedness,<:Domain}) = false
is_signed_finite(@nospecialize fmt::Format{<:Signedness,<:Domain}) = false
is_signed_extended(@nospecialize fmt::Format{<:Signedness,<:Domain}) = false

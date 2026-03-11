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


# =========================================================================
# Accessors
# =========================================================================

"""Total number of bits per code point."""
BitwidthOf(@nospecialize fmt::Format) = fmt.K

"""Precision: number of significand bits including the implicit leading bit."""
PrecisionOf(@nospecialize fmt::Format) = fmt.P

"""Number of trailing significand bits (`P - 1`)."""
TrailingBitsOf(@nospecialize fmt::Format) = fmt.P - 1

"""Number of sign bits: `0` for unsigned, `1` for signed."""
SignBitsOf(@nospecialize fmt::Format) = 0 + is_signed(fmt)

"""Number of non-significand bits (`K - P`)."""
NonSignificantBitsOf(@nospecialize fmt::Format) = BitwidthOf(fmt) - PrecisionOf(fmt)

"""Exponent field width in bits. Unsigned: `K - P + 1`. Signed: `K - P`."""
ExponentBitsOf(@nospecialize fmt::Format{UnsignedFormat,T}) where T = NonSignificantBitsOf(fmt) + 1
ExponentBitsOf(@nospecialize fmt::Format{SignedFormat,T}) where T = NonSignificantBitsOf(fmt)

"""Exponent bias `B`. Unsigned: `2^(K-P)`. Signed: `2^(K-P-1)`."""
ExponentBiasOf(@nospecialize fmt::Format{UnsignedFormat,T}) where T = UInt128(1) << NonSignificantBitsOf(fmt)
ExponentBiasOf(@nospecialize fmt::Format{SignedFormat,T}) where T = UInt128(1) << (NonSignificantBitsOf(fmt) - 1)

"""Exponent unbiased minimum value"""
ExponentMinOf(@nospecialize fmt::Format) = 1 - ExponentBiasOf(fmt)

"""Exponent unbiased maximum value"""
function ExponentMaxOf(@nospecialize fmt::Format{UnsignedFormat,T}) where T
    if PrecisionOf(fmt) > 2
        ExponentBiasOf(fmt) - 1
    else
        adjust = (2 - PrecisionOf(fmt)) * is_extended(fmt) + (1 - PrecisionOf(fmt) - is_extended(fmt))
        twopow(BitwidthOf(fmt) - PrecisionOf(fmt)) - 1 - adjust
    end
end

function ExponentMaxOf(@nospecialize fmt::Format{SignedFormat,T}) where T
    if PrecisionOf(fmt) > 1
        ExponentBiasOf(fmt) - 1
    else
        adjust = (2 - PrecisionOf(fmt)) * is_extended(fmt)
        twopow(BitwidthOf(fmt) - PrecisionOf(fmt)) - 1 - adjust
    end
end

"""Exponent spread `E` (number of distinct exponent values)."""
ExponentSpreadOf(@nospecialize fmt::Format) = ExponentMaxOf(fmt) - ExponentMinOf(fmt) + 1


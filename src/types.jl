setprecision(BigFloat, 1024)
const CPOINT = UInt64
const FPVALUE = BigFloat

# =========================================================================
# Traits for Signedness and Domain
# =========================================================================

"""Abstract supertype for signedness traits (`is_signed`, `is_unsigned`)."""
abstract type Signedness end

"""Trait type indicating a signed format (one sign bit)."""
struct is_signed <: Signedness end

"""Trait type indicating an unsigned format (no sign bit)."""
struct is_unsigned <: Signedness end

"""Test whether `x` is the `is_signed` trait."""
is_signed(x::Signedness) = x isa is_signed

"""Test whether `x` is the `is_unsigned` trait."""
is_unsigned(x::Signedness) = x isa is_unsigned

"""Abstract supertype for domain traits (`is_finite`, `is_extended`)."""
abstract type Domain end

"""Trait type indicating a finite-only format (no infinities)."""
struct is_finite <: Domain end

"""Trait type indicating an extended format (has ±Inf code points)."""
struct is_extended <: Domain end

"""Test whether `x` is the `is_finite` trait."""
is_finite(x::Domain) = x isa is_finite

"""Test whether `x` is the `is_extended` trait."""
is_extended(x::Domain) = x isa is_extended

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
is_unsigned(fnt::Format{is_unsigned,T}) where T = true
is_unsigned(fnt::Format{is_signed,T}) where T = false

"""Return `true` if `fmt` is signed."""
is_signed(fnt::Format{is_unsigned,T}) where T = false
is_signed(fnt::Format{is_signed,T}) where T = true

"""Return `true` if `fmt` is finite-only (no infinities)."""
is_finite(fnt::Format{T,is_finite}) where T = true
is_finite(fnt::Format{T,is_extended}) where T = false

"""Return `true` if `fmt` is extended (has ±Inf)."""
is_extended(fnt::Format{T,is_finite}) where T = false
is_extended(fnt::Format{T,is_extended}) where T = true


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
ExponentBitsOf(@nospecialize fmt::Format{is_unsigned,T}) where T = NonSignificantBitsOf(fmt) + 1
ExponentBitsOf(@nospecialize fmt::Format{is_signed,T}) where T = NonSignificantBitsOf(fmt)

"""Exponent bias `B`. Unsigned: `2^(K-P)`. Signed: `2^(K-P-1)`."""
ExponentBiasOf(@nospecialize fmt::Format{is_unsigned,T}) where T = UInt128(1) << NonSignificantBitsOf(fmt)
ExponentBiasOf(@nospecialize fmt::Format{is_signed,T}) where T = UInt128(1) << (NonSignificantBitsOf(fmt) - 1)

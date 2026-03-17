#=
Prompt 1: Formal Mathematical Specification Generator
You are a numerical analysis researcher and floating-point systems architect
specializing in custom floating-point algebraic systems. Your task is to help me
develop a complete formal algebra for P3109 floating-point formats. 
I would like you to construct a rigorous mathematical specification 
that defines the algebraic structure of the format family. 
Include: (1) formal definitions of the format parameters (K, P, signedness, domain), 
(2) the code-point space and mapping to numerical values, 
(3) definitions of subnormal, normal, zero, infinities, and NaN, (
(4) algebraic operations on code points and on decoded values,
(5) closure properties, 
(6) ordering relations, 
(7) binade structure and exponent lattice, 
(8) exact rational representation of values, and 
(9) proofs or derivations for key identities. 
Present the result as a well-structured technical document (20–40 pages equivalent) 
with equations, definitions, and theorems. Provide the output in LaTeX source format
and also produce a downloadable compiled PDF of the algebra specification and also
produce a downloadable markdown versionof the technical cocument using Unicode
and not using latex macros in the md file.
Write in a formal academic tone suitable for researchers 
in numerical computing and computer arithmetic.

Prompt 2: Implementation-Oriented Algebra and Reference Library
Imagine you are a computer arithmetic engineer designing reference implementations 
for novel floating-point systems. Your task is to help me derive a 
complete algebraic framework for P3109 floating-point formats and translate it
into an implementable reference library. I would like the output to include: 
(1) a formal algebra describing the structure of the format family, 
(2) definitions for code-point arithmetic and value decoding, 
(3) rules governing subnormal, normal, infinity, and NaN representations, 
(4) algebraic operations such as sign transformations, ordinal indexing, exponent binades, and rational decoding, 
(5) algorithms for encoding and decoding values, and 
(6) enumerators for all values in a format. 
Provide well-documented source code (Julia v1.12) implementing the algebra. 
Include unit tests, examples, and documentation. 
Package the final deliverable as a downloadable archive (ZIP) 
containing the library, documentation, and example scripts. \
Use a technical yet practical tone aimed at engineers implementing floating-point systems.

Prompt 3: Computational Algebra Generator with Dataset Export
Assume the role of an expert in computer arithmetic, algebraic number systems, 
and floating-point standardization. Your task is to help me construct a 
complete algebraic model of P3109 floating-point formats and generate machine-readable
outputs describing the system. 
The model should define: 
(1) the parameterized format Format(K,P,S,D), 
(2) the code-point algebra and transformations, 
(3) exact rational value decoding for every code point, 
(4) classification of values (zero, subnormal, normal, infinity, NaN), 
(5) algebraic relations between binades and significands, 
(6) symmetry relations between positive and negative halves, and 
(7) counting formulas for all value classes. 
Use precise mathematical notation and step-by-step derivations. 
Then generate exportable datasets containing the complete value sets and 
structural properties for example formats. 
Provide the outputs as downloadable files (CSV, JSON, and Markdown documentation). 
Present the explanation in a clear, research-grade technical tone 
for mathematicians and systems designers.
=#


setprecision(BigFloat, 1024)

"""
    twopow(n)

Compute `twopow` (2^n) as a `UInt128`.

This is intended for structural quantities and code-point calculations where
the exponent is known to fit in `UInt128` shifting semantics.
"""
twopow(n::Integer) = n >= 0 ? UInt128(1) << Int(n) : BigFloat(1) / (UInt128(1) << Int(-n))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

abstract type FormatTrait end
abstract type Signedness <: FormatTrait end
abstract type Domain <: FormatTrait end

struct Signedness{T} end

const UnsignedFormat = Signedness{0}
const SignedFormat = Signedness{1}
const IsUnsigned = UnsignedFormat()
const IsSigned = SignedFormat()

is_unsigned(x::UnsignedFormat) = true
is_unsigned(x::SignedFormat) = false
is_signed(x::UnsignedFormat) = false
is_signed(x::SignedFormat) = true

struct Domain{T} end

const FiniteFormat = Domain{0}
const ExtendedFormat = Domain{1}
const IsFinite = FiniteFormat()
const IsExtended = ExtendedFormat()

is_finite(x::FiniteFormat) = true
is_finite(x::ExtendedFormat) = false
is_extended(x::FiniteFormat) = false
is_extended(x::ExtendedFormat) = true

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

"""
    sign_half_offset(f)

Offset between the positive and negative halves of the code-point space.
Returns `0` for unsigned formats.
"""
sign_half_offset(@nospecialize fmt::Format{UnsignedFormat,T}) where T = 0
sign_half_offset(@nospecialize fmt::Format{SignedFormat,T}) where T = twopow(BitwidthOf(fmt) - 1)

"""Signed offset for code-point arithmetic: `0` (unsigned) or `1 - 2^(K-1)` (signed)."""
sign_offset(@nospecialize fmt::Format{UnsignedFormat,T}) where T = 0
sign_offset(@nospecialize fmt::Format{SignedFormat,T}) where T = 1 - twopow(BitwidthOf(fmt) - 1)

"""
    sign_reduce(fmt, cp)

Reduce a code point to positive-half form.

Returns a named tuple `(s, cp_abs)` where

- `s` is the sign indicator (`0` or `1`)
- `cp_abs` is the positive-half code point

For signed formats, the midpoint is reserved for NaN and is rejected here.
"""
sign_reduce(@nospecialize(fmt::Format{UnsignedFormat,T}), cp::Integer) where T = (s=0, cp_abs=cp)

function sign_reduce(@nospecialize(fmt::Format{SignedFormat,T}), cp::Integer) where T
    H = Int(sign_half_offset(fmt))
    cp < H && return (; s=0, cp_abs=cp)
    cp > H && return (; s=1, cp_abs=cp - H)
    throw(ArgumentError("midpoint code point is reserved for NaN in signed formats"))
end

"""
    significand_scale(fmt)

Return `twopow(P-1)`, the number of trailing-significand states per binade.
"""
significand_scale(fmt::Format) = twopow(TrailingBitsOf(fmt))


"""Machine epsilon `2^(1-P)` for the format."""
EpsilonOf(@nospecialize fmt::Format) = twopow(1 - PrecisionOf(fmt))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

"""Number of NaN code points (always 1)."""
nNaNsOf(@nospecialize fmt::Format) = 1

"""Number of zero code points (always 1)."""
nZerosOf(@nospecialize fmt::Format) = 1

"""Total number of infinity code points: 0 (finite), 1 (unsigned extended), 2 (signed extended)."""
nInfsOf(@nospecialize fmt::Format{T,FiniteFormat}) where T = 0
nInfsOf(@nospecialize fmt::Format{UnsignedFormat,ExtendedFormat}) = 1
nInfsOf(@nospecialize fmt::Format{SignedFormat,ExtendedFormat}) = 2

"""Number of positive infinity code points: 0 or 1."""
nPosInfsOf(@nospecialize fmt::Format{T,FiniteFormat}) where T = 0
nPosInfsOf(@nospecialize fmt::Format{UnsignedFormat,ExtendedFormat}) = 1
nPosInfsOf(@nospecialize fmt::Format{SignedFormat,ExtendedFormat}) = 1

"""Number of negative infinity code points: 0 or 1."""
nNegInfsOf(@nospecialize fmt::Format{T,FiniteFormat}) where T = 0
nNegInfsOf(@nospecialize fmt::Format{UnsignedFormat,ExtendedFormat}) = 0
nNegInfsOf(@nospecialize fmt::Format{SignedFormat,ExtendedFormat}) = 1

"""Total number of code points: `2^K`."""
nValuesOf(@nospecialize fmt::Format) = twopow(BitwidthOf(fmt))

"""Number of numerical (non-NaN) code points: `2^K - 1`."""
nNumericalValuesOf(@nospecialize fmt::Format) = nValuesOf(fmt) - 1

"""Number of non-zero numerical code points."""
nNonZeroNumericalValuesOf(@nospecialize fmt::Format) = nNumericalValuesOf(fmt) - nZerosOf(fmt)

"""Number of finite values (numerical minus infinities), including zero."""
nFiniteValuesOf(@nospecialize fmt::Format) = nNumericalValuesOf(fmt) - nInfsOf(fmt)

"""Number of non-zero finite values."""
nNonZeroFiniteValuesOf(@nospecialize fmt::Format) = nFiniteValuesOf(fmt) - nZerosOf(fmt)

"""Number of negative non-zero values (including -Inf if present). Always 0 for unsigned."""
nNegValuesOf(@nospecialize fmt::Format{UnsignedFormat,T}) where T = 0
nNegValuesOf(@nospecialize fmt::Format{SignedFormat,T}) where T = nNonZeroNumericalValuesOf(fmt) ÷ 2

"""Number of positive non-zero values (including +Inf if present)."""
nPosValuesOf(@nospecialize fmt::Format{UnsignedFormat,T}) where T = nNonZeroNumericalValuesOf(fmt)
nPosValuesOf(@nospecialize fmt::Format{SignedFormat,T}) where T = nNonZeroNumericalValuesOf(fmt) ÷ 2

"""Number of non-negative values (positive values plus zero)."""
nNonNegValuesOf(@nospecialize fmt::Format) = nPosValuesOf(fmt) + nZerosOf(fmt)

"""Number of negative finite values."""
nNegFiniteValuesOf(@nospecialize fmt::Format) = nNegValuesOf(fmt) - nNegInfsOf(fmt)

"""Number of positive finite values."""
nPosFiniteValuesOf(@nospecialize fmt::Format) = nPosValuesOf(fmt) - nPosInfsOf(fmt)

"""Number of non-negative finite values (positive finite plus zero)."""
nNonNegFiniteValuesOf(@nospecialize fmt::Format) = nPosFiniteValuesOf(fmt) + nZerosOf(fmt)

#

"""Number of non-negative prenormal code points (zero + positive subnormals): `2^(P-1)`, or 1 when `P=1`."""
nNonNegPrenormalsOf(@nospecialize fmt::Format) = PrecisionOf(fmt) == 1 ? 1 : twopow(PrecisionOf(fmt) - 1)

"""Total number of prenormal code points (zero + all subnormals)."""
nPrenormalsOf(@nospecialize fmt::Format{UnsignedFormat,T}) where T = PrecisionOf(fmt) == 1 ? 1 : nNonNegPrenormalsOf(fmt)
nPrenormalsOf(@nospecialize fmt::Format{SignedFormat,T}) where T = PrecisionOf(fmt) == 1 ? 1 : nNonNegPrenormalsOf(fmt) + (nNonNegPrenormalsOf(fmt) - 1)

"""Number of positive subnormal values: `2^(P-1) - 1`, or 0 when `P=1`."""
nPosSubnormalsOf(@nospecialize fmt::Format) = PrecisionOf(fmt) == 1 ? 0 : nNonNegPrenormalsOf(fmt) - 1

"""Number of negative subnormal values. Always 0 for unsigned."""
nNegSubnormalsOf(@nospecialize fmt::Format{UnsignedFormat,T}) where T = 0
nNegSubnormalsOf(@nospecialize fmt::Format{SignedFormat,T}) where T = PrecisionOf(fmt) == 1 ? 0 : nPosSubnormalsOf(fmt)

"""Total number of subnormal values (positive + negative)."""
nSubnormalsOf(@nospecialize fmt::Format) = PrecisionOf(fmt) == 1 ? 0 : nPosSubnormalsOf(fmt) + nNegSubnormalsOf(fmt)

#

"""Total number of normal values (finite minus prenormals)."""
nNormalsOf(@nospecialize fmt::Format) = nFiniteValuesOf(fmt) - nPrenormalsOf(fmt)

"""Number of positive normal values."""
nPosNormalsOf(@nospecialize fmt::Format{UnsignedFormat,T}) where T = nNormalsOf(fmt)
nPosNormalsOf(@nospecialize fmt::Format{SignedFormat,T}) where T = nNormalsOf(fmt) ÷ 2

"""Number of negative normal values. Always 0 for unsigned."""
nNegNormalsOf(@nospecialize fmt::Format{UnsignedFormat,T}) where T = 0
nNegNormalsOf(@nospecialize fmt::Format{SignedFormat,T}) where T = nPosNormalsOf(fmt)

#

"""Number of normal binades: `2B - 1` where `B = ExponentBiasOf(fmt)`."""
nBinadesOf(@nospecialize fmt::Format) = 2 * ExponentBiasOf(fmt) - 1

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

"""Return the narrowest unsigned integer type that fits `K` bits."""
@inline function typeofcp(fmt::Format)
    K = BitwidthOf(fmt)
    if K <= 8
        UInt8
    elseif K <= 16
        UInt16
    elseif K <= 32
        UInt32
    elseif K <= 64
        UInt64
    elseif K <= 128
        UInt128
    else
        error("Unsupported bitwidth")
    end
end

"""Smallest code point (always 0)."""
cp_min(fmt::Format) = zero(typeofcp(fmt))

"""Largest code point: `2^K - 1`."""
cp_max(fmt::Format) = typeofcp(fmt)(twopow(BitwidthOf(fmt)) - 1)

"""Code point of zero (always 0)."""
cp_zero(fmt::Format) = zero(typeofcp(fmt))

"""Code point of NaN. Unsigned: `cp_max`. Signed: midpoint `2^(K-1)`."""
cp_nan(@nospecialize fmt::Format{UnsignedFormat,T}) where T = cp_max(fmt)
cp_nan(@nospecialize fmt::Format{SignedFormat,T}) where T = typeofcp(fmt)(sign_half_offset(fmt))

"""Code point of +Inf, or `nothing` for finite-domain formats."""
cp_inf(@nospecialize fmt::Format{UnsignedFormat,ExtendedFormat}) = cp_nan(fmt) - one(typeofcp(fmt))
cp_inf(@nospecialize fmt::Format{SignedFormat,ExtendedFormat}) = cp_nan(fmt) - one(typeofcp(fmt))
cp_inf(@nospecialize fmt::Format{UnsignedFormat,FiniteFormat}) = nothing
cp_inf(@nospecialize fmt::Format{SignedFormat,FiniteFormat}) = nothing

"""Alias for [`cp_inf`](@ref)."""
cp_posinf = cp_inf

"""Code point of -Inf, or `nothing` if absent (unsigned or finite-domain)."""
cp_neginf(@nospecialize fmt::Format{UnsignedFormat,ExtendedFormat}) = nothing
cp_neginf(@nospecialize fmt::Format{SignedFormat,ExtendedFormat}) = cp_max(fmt)
cp_neginf(@nospecialize fmt::Format{UnsignedFormat,FiniteFormat}) = nothing
cp_neginf(@nospecialize fmt::Format{SignedFormat,FiniteFormat}) = nothing

"""
    cp_pos_subnormal_min(fmt)

the code point for the positive subnormal of least magnitude, or `nothing` if P == 1
"""
cp_pos_subnormal_min(fmt::Format) = PrecisionOf(fmt) > 1 ? one(typeofcp(fmt)) : nothing

"""
    cp_pos_subnormal_max(fmt)

the code point for the positive subnormal of greatest magnitude, or `nothing` if P == 1
"""
cp_pos_subnormal_max(fmt::Format) = PrecisionOf(fmt) > 1 ? typeofcp(fmt)(twopow(TrailingBitsOf(fmt)) - 1) : nothing

"""
    cp_neg_subnormal_min(fmt)

the code point for the negative subnormal of least magnitude, or `nothing`
"""
cp_neg_subnormal_min(@nospecialize fmt::Format{UnsignedFormat,T}) where T = nothing
cp_neg_subnormal_min(@nospecialize fmt::Format{SignedFormat,T}) where T =
    PrecisionOf(fmt) > 1 ? cp_nan(fmt) + one(typeofcp(fmt)) : nothing

"""
    cp_neg_subnormal_max(fmt)

the code point for the negative subnormal of greatest magnitude, or `nothing`
"""
cp_neg_subnormal_max(@nospecialize fmt::Format{UnsignedFormat,T}) where T = nothing
cp_neg_subnormal_max(@nospecialize fmt::Format{SignedFormat,T}) where T =
    PrecisionOf(fmt) > 1 ? typeofcp(fmt)(twopow(BitwidthOf(fmt) - 1) + twopow(PrecisionOf(fmt) - 1) - 1) : nothing

"""
    unsafe_cp_ordinal_ith_pos_subnormal(fmt, i)

the code point for the (1-based) ith positive subnormal value — no bounds check
"""
unsafe_cp_ordinal_ith_pos_subnormal(fmt::Format, i::Integer) = typeofcp(fmt)(i)

"""
    cp_ordinal_ith_pos_subnormal(fmt, i)

the code point for the (1-based) ith positive subnormal value
"""
function cp_ordinal_ith_pos_subnormal(fmt::Format, i::Integer)
    imax = nPosSubnormalsOf(fmt)
    1 <= i <= imax || throw(ArgumentError("i ($i) must be in 1..$imax"))
    unsafe_cp_ordinal_ith_pos_subnormal(fmt, i)
end

"""
    cp_pos_normal_min(fmt)

the code point for the positive normal of least magnitude
"""
cp_pos_normal_min(fmt::Format) = typeofcp(fmt)(twopow(PrecisionOf(fmt) - 1))

"""
    cp_pos_normal_max(fmt)

the code point for the positive normal of greatest magnitude
"""
cp_pos_normal_max(@nospecialize fmt::Format{UnsignedFormat,FiniteFormat}) = typeofcp(fmt)(twopow(BitwidthOf(fmt)) - 2)
cp_pos_normal_max(@nospecialize fmt::Format{UnsignedFormat,ExtendedFormat}) = typeofcp(fmt)(twopow(BitwidthOf(fmt)) - 3)
cp_pos_normal_max(@nospecialize fmt::Format{SignedFormat,FiniteFormat}) = cp_nan(fmt) - one(typeofcp(fmt))
cp_pos_normal_max(@nospecialize fmt::Format{SignedFormat,ExtendedFormat}) = cp_nan(fmt) - 2 * one(typeofcp(fmt))

"""
    cp_neg_normal_min(fmt)

the code point for the negative normal of least magnitude, or `nothing`
"""
cp_neg_normal_min(@nospecialize fmt::Format{UnsignedFormat,T}) where T = nothing
function cp_neg_normal_min(@nospecialize fmt::Format{SignedFormat,T}) where T
    if PrecisionOf(fmt) > 1
        cp_neg_subnormal_max(fmt) + one(typeofcp(fmt))
    else
        cp_nan(fmt) + one(typeofcp(fmt))
    end
end

"""
    cp_neg_normal_max(fmt)

the code point for the negative normal of greatest magnitude, or `nothing`
"""
cp_neg_normal_max(@nospecialize fmt::Format{UnsignedFormat,T}) where T = nothing
cp_neg_normal_max(@nospecialize fmt::Format{SignedFormat,ExtendedFormat}) = cp_max(fmt) - one(typeofcp(fmt))
cp_neg_normal_max(@nospecialize fmt::Format{SignedFormat,FiniteFormat}) = cp_max(fmt)

"""
    unsafe_cp_ordinal_ith_pos_normal(fmt, i)

the code point for the (1-based) ith positive normal value — no bounds check
"""
unsafe_cp_ordinal_ith_pos_normal(fmt::Format, i::Integer) = typeofcp(fmt)(twopow(PrecisionOf(fmt) - 1) + (i - 1))

"""
    cp_ordinal_ith_pos_normal(fmt, i)

the code point for the (1-based) ith positive normal value
"""
function cp_ordinal_ith_pos_normal(fmt::Format, i::Integer)
    imax = nPosNormalsOf(fmt)
    1 <= i <= imax || throw(ArgumentError("i ($i) must be in 1..$imax"))
    unsafe_cp_ordinal_ith_pos_normal(fmt, i)
end

"""
    cp_positive_max(fmt)

code point of the largest positive value (may be infinity)
"""
cp_positive_max(fmt::Format) = cp_nan(fmt) - one(typeofcp(fmt))

"""Return `true` if `cp` is a positive code point (between zero and NaN, exclusive)."""
function cp_is_positive(@nospecialize(fmt::Format), cp::Integer)
    cp_zero(fmt) < cp < cp_nan(fmt)
end

"""Return `true` if `cp` is a non-negative code point (zero or positive)."""
function cp_is_nonnegative(@nospecialize(fmt::Format), cp::Integer)
    cp_zero(fmt) <= cp < cp_nan(fmt)
end

"""Return `true` if `cp` is a negative code point. Always `false` for unsigned formats."""
cp_is_negative(@nospecialize(fmt::Format{UnsignedFormat,T}), cp::Integer) where T = false
function cp_is_negative(@nospecialize(fmt::Format{SignedFormat,T}), cp::Integer) where T
    cp > cp_nan(fmt)
end

"""
    pos_cp_to_neg_cp(fmt, cp_pos)

given fmt, cp(+value), find cp(-value)
"""
function pos_cp_to_neg_cp(@nospecialize(fmt::Format{SignedFormat,T}), cp_pos::Integer) where T
    cp_is_positive(fmt, cp_pos) || throw(ArgumentError("cp_pos ($cp_pos) must be a positive code point"))
    unsafe_pos_cp_to_neg_cp(fmt, cp_pos)
end

"""Map a positive code point to its negative counterpart — no bounds check."""
unsafe_pos_cp_to_neg_cp(@nospecialize(fmt::Format{SignedFormat,T}), cp_pos::Integer) where T =
    twopow(BitwidthOf(fmt) - 1) + cp_pos

unsafe_pos_cp_to_neg_cp(@nospecialize(fmt::Format{UnsignedFormat,T}), cp_pos::Integer) where T = nothing

"""
    neg_cp_to_pos_cp(fmt, cp_neg)

given fmt, cp(-value), find cp(+value)
"""
function neg_cp_to_pos_cp(@nospecialize(fmt::Format{SignedFormat,T}), cp_neg::Integer) where T
    cp_is_negative(fmt, cp_neg) || throw(ArgumentError("cp_neg ($cp_neg) must be a negative code point"))
    unsafe_neg_cp_to_pos_cp(fmt, cp_neg)
end

"""Map a negative code point to its positive counterpart — no bounds check."""
unsafe_neg_cp_to_pos_cp(@nospecialize(fmt::Format{SignedFormat,T}), cp_neg::Integer) where T =
    cp_neg - twopow(BitwidthOf(fmt) - 1)

unsafe_neg_cp_to_pos_cp(@nospecialize(fmt::Format{UnsignedFormat,T}), cp_neg::Integer) where T = nothing

"""
    cp_changesign(fmt, cp)

given fmt, cp, find corresponding sign-changed code point
"""
function cp_changesign(@nospecialize(fmt::Format{SignedFormat,T}), cp::Integer) where T
    if cp == cp_zero(fmt) || cp == cp_nan(fmt)
        return cp
    end
    if cp_is_positive(fmt, cp)
        unsafe_pos_cp_to_neg_cp(fmt, cp)
    else
        unsafe_neg_cp_to_pos_cp(fmt, cp)
    end
end

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

"""
    val_pos_subnormal_min(fmt)

Exact value of the smallest positive subnormal, or `nothing` if `P=1`.
"""
val_pos_subnormal_min(fmt::Format) = PrecisionOf(fmt) > 1 ? twopow(2 - PrecisionOf(fmt) - ExponentBiasOf(fmt)) : nothing

"""
    val_pos_subnormal_max(fmt)

Exact value of the largest positive subnormal, or `nothing` if `P=1`.
"""
function val_pos_subnormal_max(fmt::Format)
    P = PrecisionOf(fmt)
    B = ExponentBiasOf(fmt)
    P == 1 && return nothing
    twopow(1 - B) * twopow(1 - twopow(1 - P))
end

"""Exact value of the `i`-th positive subnormal (1-based), or `nothing` if `P=1`."""
function val_ordinal_ith_pos_subnormal(fmt::Format, i::Integer)
    P = PrecisionOf(fmt)
    P == 1 && return nothing
    B = ExponentBiasOf(fmt)
    i * twopow(2 - P - B)
end

"""Exact value of the `cardinal`-th positive subnormal (0-based), or `nothing` if `P=1`."""
function val_cardinal_ith_pos_subnormal(fmt::Format, cardinal::Integer)
    P = PrecisionOf(fmt)
    P == 1 && return nothing
    B = ExponentBiasOf(fmt)
    (cardinal + 1) * twopow(2 - P - B)
end

"""
    val_pos_normal_min(fmt)

Exact value of the smallest positive normal: `2^(1 - B)`.
"""
val_pos_normal_min(fmt::Format) = twopow(1 - ExponentBiasOf(fmt))

"""
    val_pos_normal_max(fmt)

Exact value of the largest positive finite normal.
"""
val_pos_normal_max(fmt::Format{UnsignedFormat,FiniteFormat}) =
    (twopow(twopow(ExponentBiasOf(fmt)) - 1)) * (2 - twopow(1 - PrecisionOf(fmt)))

val_pos_normal_max(fmt::Format{UnsignedFormat,ExtendedFormat}) =
    (twopow(twopow(ExponentBiasOf(fmt)) - 1)) * (2 - 3 * twopow(1 - PrecisionOf(fmt)))

val_pos_normal_max(fmt::Format{SignedFormat,FiniteFormat}) =
    (twopow(twopow(ExponentBiasOf(fmt)) - 1)) * (2 - twopow(1 - PrecisionOf(fmt)))

val_pos_normal_max(fmt::Format{SignedFormat,ExtendedFormat}) =
    (twopow(twopow(ExponentBiasOf(fmt)) - 1)) * (2 - twopow(2 - PrecisionOf(fmt)))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



"""
    ValueOf(fmt, cp) -> Qx64

Return the exact `Qx64` value of code point `cp` in format `fmt`.

Special code points map to the canonical `Qx64` forms:

- zero  → `0//1`
- +Inf  → `1//0`
- -Inf  → `-1//0`
- NaN   → `0//0`

For finite numerical code points the result is an exact dyadic rational.
"""
function ValueOf(@nospecialize(fmt::Format), cp::Integer)
    0 <= cp <= cp_max(fmt) || throw(ArgumentError("code point $cp out of range [0, $(cp_max(fmt))]"))

    cp == cp_zero(fmt) && return zero(Qx64)
    cp == cp_nan(fmt) && return NaN(Qx64)

    inf_cp = cp_inf(fmt)
    inf_cp !== nothing && cp == inf_cp && return Inf(Qx64)
    ninf_cp = cp_neginf(fmt)
    ninf_cp !== nothing && cp == ninf_cp && return NegInf(Qx64)

    red = sign_reduce(fmt, cp)
    val = _decode_positive_half(fmt, red.cp_abs)
    return red.s == 0 ? val : -val
end

"""
    FiniteValueOf(fmt, cp) -> Qx64

Return the exact rational value of a finite code point `cp` in format `fmt`.

The code point must be finite (not NaN, not ±Inf).  Zero is accepted.
"""
function FiniteValueOf(@nospecialize(fmt::Format), cp::Integer)
    0 <= cp <= cp_max(fmt) || throw(ArgumentError("code point $cp out of range [0, $(cp_max(fmt))]"))
    cp == cp_zero(fmt) && return zero(Qx64)

    red = sign_reduce(fmt, cp)
    val = _decode_positive_half(fmt, red.cp_abs)
    return red.s == 0 ? val : -val
end

"""
    _decode_positive_half(fmt, cp_abs) -> Qx64

Decode a positive-half code point `cp_abs` (after sign reduction) to its
exact rational value.  Caller must ensure `cp_abs > 0`.

Let `m = twopow(P-1)` (the significand scale) and `B = ExponentBias`.

Subnormal (`cp_abs < m`):  value = cp_abs · twopow(2 − P − B)
Normal    (`cp_abs ≥ m`):  e = cp_abs ÷ m,  t = cp_abs mod m,  k = e − B
                           value = twopowk + t · twopow(k + 1 − P)
"""
function _decode_positive_half(@nospecialize(fmt::Format), cp_abs::Integer)
    P = PrecisionOf(fmt)
    B = BigInt(ExponentBiasOf(fmt))
    m = BigInt(significand_scale(fmt))      # twopow(P-1)
    ca = BigInt(cp_abs)

    if ca < m
        # subnormal
        return Qx64(dyadic_twopow(2 - P - Int(B))) * Qx64(ca)
    else
        # normal: binade decomposition
        e = fld(ca, m)
        t = mod(ca, m)
        k = Int(e - B)
        return Qx64(dyadic_twopow(k)) + Qx64(t) * Qx64(dyadic_twopow(k + 1 - P))
    end
end

# =========================================================================
# Convenience: decode every code point at once
# =========================================================================

"""
    AllFiniteValuesOf(fmt) -> Vector{Qx64}

Return the exact rational values for every finite numerical code point
in `fmt`, in code-point order.  Zero is included; NaN and ±Inf are excluded.
"""
function AllFiniteValuesOf(@nospecialize(fmt::Format))
    vals = Qx64[]
    sizehint!(vals, Int(nFiniteValuesOf(fmt)) + 1)

    nan_cp = cp_nan(fmt)
    inf_cp = cp_inf(fmt)
    ninf_cp = cp_neginf(fmt)

    for cp in BigInt(0):BigInt(cp_max(fmt))
        cp == nan_cp && continue
        inf_cp !== nothing && cp == inf_cp && continue
        ninf_cp !== nothing && cp == ninf_cp && continue
        push!(vals, FiniteValueOf(fmt, cp))
    end
    return vals
end

"""
    AllPositiveFiniteValuesOf(fmt) -> Vector{Qx64}

Return the exact rational values for every positive finite code point
in `fmt`, in code-point order.  Zero is excluded.
"""
function AllPositiveFiniteValuesOf(@nospecialize(fmt::Format))
    vals = Qx64[]
    sizehint!(vals, Int(nPosFiniteValuesOf(fmt)))

    cpstart = BigInt(cp_zero(fmt)) + 1
    cpend = BigInt(cp_nan(fmt)) - 1

    inf_cp = cp_inf(fmt)
    if inf_cp !== nothing
        cpend = min(cpend, BigInt(inf_cp) - 1)
    end

    for cp in cpstart:cpend
        push!(vals, FiniteValueOf(fmt, cp))
    end
    return vals
end

# =========================================================================
# Ordinal value accessors (1-based)
# =========================================================================

"""
    FiniteValueOfOrdinalPos(fmt, i) -> Qx64

Return the exact value of the `i`-th positive finite value (1-based),
ordered by increasing magnitude.

`i = 1` is the smallest positive value, `i = nPosFiniteValuesOf(fmt)` is
the largest.
"""
function FiniteValueOfOrdinalPos(@nospecialize(fmt::Format), i::Integer)
    n = nPosFiniteValuesOf(fmt)
    1 <= i <= n || throw(ArgumentError("ordinal $i out of range [1, $n]"))
    cp = BigInt(cp_zero(fmt)) + BigInt(i)
    return FiniteValueOf(fmt, cp)
end

"""
    FiniteValueOfOrdinalNeg(fmt, i) -> Qx64

Return the exact value of the `i`-th negative finite value (1-based),
ordered by increasing magnitude (i.e. `i = 1` is the negative value
closest to zero).

Only valid for signed formats.
"""
function FiniteValueOfOrdinalNeg(@nospecialize(fmt::Format{SignedFormat,T}), i::Integer) where T
    n = nNegFiniteValuesOf(fmt)
    1 <= i <= n || throw(ArgumentError("ordinal $i out of range [1, $n]"))
    cp = BigInt(cp_nan(fmt)) + BigInt(i)
    return FiniteValueOf(fmt, cp)
end

FiniteValueOfOrdinalNeg(@nospecialize(fmt::Format{UnsignedFormat,T}), i::Integer) where T =
    throw(ArgumentError("unsigned formats have no negative values"))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>




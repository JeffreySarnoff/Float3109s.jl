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


#=
        prior realization

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

=#
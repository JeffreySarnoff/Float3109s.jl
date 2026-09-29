export FINITE, EXTENDED,
    is_finite, is_extended, DomainOf

public FiniteFormat, ExtendedFormat, Domain

struct FiniteFormat end
struct ExtendedFormat end

is_finite(::Type{FiniteFormat}) = true
is_finite(::Type{ExtendedFormat}) = false
is_extended(::Type{FiniteFormat}) = false
is_extended(::Type{ExtendedFormat}) = true

const FINITE = FiniteFormat()
const EXTENDED = ExtendedFormat()

is_finite(::FiniteFormat) = true
is_finite(::ExtendedFormat) = false
is_extended(::FiniteFormat) = false
is_extended(::ExtendedFormat) = true

const Domain = Union{FiniteFormat, ExtendedFormat}

DomainOf(x::Domain) = x

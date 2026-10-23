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
const DomainType = Union{Type{FiniteFormat}, Type{ExtendedFormat}}
const Domains = Union{DomainType, Domain}

DomainOf(x::Domain) = x

format_char(::Type{FiniteFormat}) = FiniteFormatChar
format_char(::Type{ExtendedFormat}) = ExtendedFormatChar
format_char(::FiniteFormat) = FiniteInstanceChar
format_char(::ExtendedFormat) = ExtendedInstanceChar

Base.print(io::IO, x::Domains) =
    print(io, format_char(x))

Base.show(io::IO, x::Domains) =
    print(io, x)

Base.show(io::IO, ::MIME"text/plain", x::Domains) =
    print(io, x)

Base.string(::Type{FiniteFormat}) = FiniteStr
Base.string(::Type{ExtendedFormat}) = ExtendedStr
Base.string(::FiniteFormat) = FiniteFormatStr
Base.string(::ExtendedFormat) = ExtendedFormatStr

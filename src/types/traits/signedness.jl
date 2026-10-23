export UNSIGNED, SIGNED,
    is_unsigned, is_signed, SignednessOf

public UnsignedFormat, SignedFormat, Signedness

struct UnsignedFormat end
struct SignedFormat end

is_unsigned(::Type{UnsignedFormat}) = true
is_unsigned(::Type{SignedFormat}) = false
is_signed(::Type{UnsignedFormat}) = false
is_signed(::Type{SignedFormat}) = true

const UNSIGNED = UnsignedFormat()
const SIGNED = SignedFormat()

is_unsigned(::UnsignedFormat) = true
is_unsigned(::SignedFormat) = false
is_signed(::UnsignedFormat) = false
is_signed(::SignedFormat) = true

const Signedness = Union{UnsignedFormat, SignedFormat}
const SignednessType = Union{Type{UnsignedFormat}, Type{SignedFormat}}
const Signednesses = Union{SignednessType, Signedness}

SignednessOf(x::Signedness) = x

format_char(::Type{UnsignedFormat}) = UnsignedFormatChar
format_char(::Type{SignedFormat}) = SignedFormatChar
format_char(::UnsignedFormat) = UnsignedInstanceChar
format_char(::SignedFormat) = SignedInstanceChar

Base.print(io::IO, x::Signednesses) =
    print(io, format_char(x))

Base.show(io::IO, x::Signednesses) =
    print(io, x)

Base.show(io::IO, ::MIME"text/plain", x::Signednesses) =
    print(io, x)

Base.string(::Type{UnsignedFormat}) = UnsignedStr
Base.string(::Type{SignedFormat}) = SignedStr
Base.string(::UnsignedFormat) = UnsignedFormatStr
Base.string(::SignedFormat) = SignedFormatStr

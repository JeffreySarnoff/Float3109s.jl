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

SignednessOf(x::Signedness) = x


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


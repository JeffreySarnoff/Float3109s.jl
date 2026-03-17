
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

Base.eps(@nospecialize fmt::Format) = EpsilonOf(fmt)
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

"""Code point of one, negone."""
cp_one(fmt::Format{UnsignedFormat,T}) where {T} = twopow(BitwidthOf(fmt)) - 1
cp_one(fmt::Format{SignedFormat,T}) where {T} = twopow(BitwidthOf(fmt) - 2)

cp_negone(fmt::Format{UnsignedFormat,T}) where {T} = nothing
cp_negone(fmt::Format{SignedFormat,T}) where {T} = 3 * cp_one(fmt)

"""Code point of one-half, negone-half."""
cp_onehalf(fmt::Format{UnsignedFormat,T}) where {T} = twopow(GBitwidthOf(fmt) - 1) - twopow(PrecisionOf(fmt) - 1)
cp_onehalf(fmt::Format{SignedFormat,T}) where {T} = twopow(BitwidthOf(fmt) - 2) - twopow(PrecisionOf(fmt) - 1)

cp_negonehalf(fmt::Format{UnsignedFormat,T}) where {T} = nothing
cp_negonehalf(fmt::Format{SignedFormat,T}) where {T} = 3 * twopow(BitwidthOf(fmt) - 2) - twopow(PrecisionOf(fmt) - 1)

"""Code point of two, negtwo."""
cp_two(fmt::Format{UnsignedFormat,T}) where {T} = twopow(BitwidthOf(fmt) - 1) + twopow(PrecisionOf(fmt) - 1)
cp_two(fmt::Format{SignedFormat,T}) where {T} = twopow(BitwidthOf(fmt) - 2) + twopow(PrecisionOf(fmt) - 1)

cp_negtwo(fmt::Format{UnsignedFormat,T}) where {T} = nothing
cp_negtwo(fmt::Format{SignedFormat,T}) where {T} = 3 * twopow(BitwidthOf(fmt) - 2) + twopow(PrecisionOf(fmt) - 1)


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


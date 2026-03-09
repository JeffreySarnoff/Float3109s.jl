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
val_pos_normal_max(fmt::Format{is_unsigned,is_finite}) =
    (twopow(twopow(ExponentBiasOf(fmt)) - 1)) * (2 - twopow(1 - PrecisionOf(fmt)))

val_pos_normal_max(fmt::Format{is_unsigned,is_extended}) =
    (twopow(twopow(ExponentBiasOf(fmt)) - 1)) * (2 - 3 * twopow(1 - PrecisionOf(fmt)))

val_pos_normal_max(fmt::Format{is_signed,is_finite}) =
    (twopow(twopow(ExponentBiasOf(fmt)) - 1)) * (2 - twopow(1 - PrecisionOf(fmt)))

val_pos_normal_max(fmt::Format{is_signed,is_extended}) =
    (twopow(twopow(ExponentBiasOf(fmt)) - 1)) * (2 - twopow(2 - PrecisionOf(fmt)))

function val_pos_normal_max(fmt::Format{is_unsigned,is_extended})
    P = PrecisionOf(fmt)
    if P >= 3
        (twopow(P) - 3) * twopow(ExponentBiasOf(fmt) - P)
    elseif P == 2
        3 * twopow(ExponentBiasOf(fmt) - 3)
    else # P == 1
        twopow(ExponentBiasOf(fmt) - 3)
    end
end




import Base: eps

"""
    twopow(n)

Compute `twopown` as a `UInt128`.

This is intended for structural quantities and code-point calculations where
the exponent is known to fit in `UInt128` shifting semantics.
"""
twopow(n::Integer) = n >= 0 ? UInt128(1) << Int(n) : BigFloat(1) / (UInt128(1) << Int(-n))

"""
    dyadic_twopow(n)

Return the exact dyadic rational `twopown`.
"""
function dyadic_twopow(n::Integer)
    return n >= 0 ? (BigInt(1) << Int(n)) // 1 : 1 // (BigInt(1) << Int(-n))
end

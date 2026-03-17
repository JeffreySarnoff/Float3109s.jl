import Base: eps

"""
    twopow(n)

Compute `twopown` as a `Int128`.

This is intended for structural quantities and code-point calculations where
the exponent is known to fit in `Int128` shifting semantics.
"""
twopow(n::Integer) = n >= 0 ? Int128(1) << Int(n) : Int128(1) // (Int128(1) << Int(-n))

"""
    dyadic_twopow(n)

Return the exact dyadic rational `twopown`.
"""
function dyadic_twopow(n::Integer)
    return n >= 0 ? (BigInt(1) << Int(n)) // 1 : 1 // (BigInt(1) << Int(-n))
end

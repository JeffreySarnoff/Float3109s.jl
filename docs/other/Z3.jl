import z3

# 1. Define the Bit-Vector Code Point
K = 8
E_bits, M_bits = 4, 3
bias = 7
c = z3.BitVec('c', K)

# 2. Extract Fields
sign = z3.Extract(K-1, K-1, c)
exp  = z3.Extract(K-2, M_bits, c)
mant = z3.Extract(M_bits-1, 0, c)

# 3. Identify the Unique P3109 NaN (0x80 for K=8)
is_nan = (c == z3.BitVecVal(2**(K-1), K))

# Convert to integer sorts for math
sign_int = z3.BV2Int(sign)
exp_int  = z3.BV2Int(exp)
mant_int = z3.BV2Int(mant)

# Unroll power of 2 for the exponent (since E is small)
def pow2(e_val):
    min_e = -bias + 1
    max_e = (2**E_bits - 1) - bias
    res = z3.RealVal(1.0)
    for e in range(min_e, max_e + 1):
        pow_val = z3.RealVal(2**e) if e >= 0 else z3.RealVal(1.0 / (2**-e))
        res = z3.If(e_val == e, pow_val, res)
    return res

# Calculate numeric values
sign_factor = z3.If(sign_int == 1, z3.RealVal(-1), z3.RealVal(1))
mant_frac   = z3.ToReal(mant_int) / z3.RealVal(2**M_bits)

# Normal and Subnormal equations
normal_val    = sign_factor * pow2(exp_int - bias) * (z3.RealVal(1) + mant_frac)
subnormal_val = sign_factor * pow2(1 - bias) * mant_frac

# Final decoded value
is_zero_or_subnormal = (exp_int == 0)
decode_val = z3.If(is_zero_or_subnormal, subnormal_val, normal_val)

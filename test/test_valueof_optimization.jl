using Test
using Float3109s

function all_formats_types(K, P)
    fmts = Format[]
    for S in (UnsignedFormat, SignedFormat), D in (FiniteFormat, ExtendedFormat)
        Σ = S === SignedFormat ? 1 : 0
        W = K - P + 1 - Σ
        W >= 1 || continue
        push!(fmts, Format{S,D}(K, P))
    end
    fmts
end

function legacy_decode_positive_half(fmt::Format, cp_abs::Integer)
    P = PrecisionOf(fmt)
    B = BigInt(ExponentBiasOf(fmt))
    m = BigInt(significand_scale(fmt))
    ca = BigInt(cp_abs)

    if ca < m
        return Qx64(dyadic_twopow(2 - P - Int(B))) * Qx64(ca)
    end

    e = fld(ca, m)
    t = mod(ca, m)
    k = Int(e - B)
    Qx64(dyadic_twopow(k)) + Qx64(t) * Qx64(dyadic_twopow(k + 1 - P))
end

function legacy_valueof(fmt::Format, cp::Integer)
    0 <= cp <= cp_max(fmt) || throw(ArgumentError("code point $cp out of range [0, $(cp_max(fmt))]"))

    cp == cp_zero(fmt) && return zero(Qx64)
    cp == cp_nan(fmt) && return Qx64(NaN)

    inf_cp = cp_inf(fmt)
    inf_cp !== nothing && cp == inf_cp && return Qx64(Inf)
    ninf_cp = cp_neginf(fmt)
    ninf_cp !== nothing && cp == ninf_cp && return Qx64(-Inf)

    red = sign_reduce(fmt, cp)
    val = legacy_decode_positive_half(fmt, red.cp_abs)
    red.s == 0 ? val : -val
end

@inline function capture_result(f)
    try
        return (ok=true, value=f(), errtype=nothing)
    catch err
        return (ok=false, value=nothing, errtype=typeof(err))
    end
end

@testset "ValueOf optimized path matches legacy behavior (K <= 10)" begin
    for K in 2:10, P in 1:(K-1), fmt in all_formats_types(K, P)
        cpmax = Int(cp_max(fmt))
        for cp in 0:cpmax
            got = capture_result(() -> ValueOf(fmt, cp))
            ref = capture_result(() -> legacy_valueof(fmt, cp))

            @test got.ok == ref.ok
            if got.ok
                @test numerator(got.value) == numerator(ref.value)
                @test denominator(got.value) == denominator(ref.value)
            else
                @test got.errtype == ref.errtype
            end
        end
    end
end

@testset "_decode_positive_half optimized path matches legacy formula (K <= 10)" begin
    for K in 2:10, P in 1:(K-1), fmt in all_formats_types(K, P)
        cpnan = Int(cp_nan(fmt))
        cpinf = cp_inf(fmt)

        for cp_abs in 1:(cpnan - 1)
            cpinf !== nothing && cp_abs == Int(cpinf) && continue
            got = _decode_positive_half(fmt, cp_abs)
            ref = legacy_decode_positive_half(fmt, cp_abs)
            @test numerator(got) == numerator(ref)
            @test denominator(got) == denominator(ref)
        end
    end
end

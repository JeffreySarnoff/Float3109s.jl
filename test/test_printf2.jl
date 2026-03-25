using Test
using Float3109s

const P2 = Float3109s.Printf2
const FMTA = P2.Format("%a")

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

@testset "Printf2 normalize leading-one matches legacy" begin
    manual = [
        "0x1p+0",
        "0x1.8p+1",
        "-0x1.8p+1",
        "0x0.2p+0",
        "0x2.fp+3",
        "Inf",
        "NaN",
    ]
    for s in manual
        @test P2._normalize_a_string_to_leading_one(s) ==
              P2._normalize_a_string_to_leading_one_legacy(s)
    end

    for x in (BigFloat(0.0), BigFloat(-0.0), BigFloat(1.0), BigFloat(-1.0), BigFloat(0.125), BigFloat(3.5))
        s = P2.format(FMTA, x)
        @test P2._normalize_a_string_to_leading_one(s) ==
              P2._normalize_a_string_to_leading_one_legacy(s)
    end
end

@testset "Printf2 normalize leading-zero matches legacy" begin
    manual = [
        "0x1p+0",
        "0x1.8p+1",
        "-0x1.8p+1",
        "0x0.2p+0",
        "0x2.fp+3",
        "Inf",
        "NaN",
    ]
    for s in manual
        @test P2._normalize_a_string_to_leading_zero(s) ==
              P2._normalize_a_string_to_leading_zero_legacy(s)
    end

    for x in (BigFloat(0.0), BigFloat(-0.0), BigFloat(1.0), BigFloat(-1.0), BigFloat(0.125), BigFloat(3.5))
        s = P2.format(FMTA, x)
        @test P2._normalize_a_string_to_leading_zero(s) ==
              P2._normalize_a_string_to_leading_zero_legacy(s)
    end
end

@testset "sprintf2a matches legacy normalization pipeline" begin
    for K in 2:8, P in 1:(K-1), fmt in all_formats_types(K, P)
        for v in AllFiniteValuesOf(fmt)
            x = BigFloat(v)
            raw = P2.format(FMTA, x)
            @test P2.sprintf2a(x; subnormal=false) == P2._normalize_a_string_to_leading_one_legacy(raw)
            @test P2.sprintf2a(x; subnormal=true) == P2._normalize_a_string_to_leading_zero_legacy(raw)
        end
    end
end


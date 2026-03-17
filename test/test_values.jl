using Test
using Float3109s

# =========================================================================
# Helper: build all four format variants for a given (K, P)
# =========================================================================

function all_formats(K, P)
    fmts = []
    for (S, D) in ((is_unsigned, is_finite), (is_unsigned, is_extended),
        (is_signed, is_finite), (is_signed, is_extended))
        Σ = S === is_signed ? 1 : 0
        W = K - P + 1 - Σ
        W >= 1 || continue
        push!(fmts, Format{S,D}(K, P))
    end
    fmts
end

# =========================================================================
# §1  ValueOf — zero always decodes to 0//1
# =========================================================================

@testset "ValueOf(cp=0) == 0 for all formats" begin
    for K in 2:10, P in 1:(K-1)
        for fmt in all_formats(K, P)
            @test ValueOf(fmt, 0) == Rational{BigInt}(0)
        end
    end
end

# =========================================================================
# §2  ValueOf — NaN and infinity code points throw
# =========================================================================

@testset "ValueOf throws on NaN code point" begin
    for K in 2:8, P in 1:(K-1)
        for fmt in all_formats(K, P)
            @test_throws ArgumentError ValueOf(fmt, cp_nan(fmt))
        end
    end
end

@testset "ValueOf throws on +Inf code point" begin
    for K in 2:8, P in 1:K
        K - P + 1 >= 1 || continue
        ue = Format{UnsignedFormat,ExtendedFormat}(K, P)
        @test_throws ArgumentError ValueOf(ue, cp_inf(ue))
    end

    for K in 3:8, P in 1:(K-1)
        se = Format{SignedFormat,ExtendedFormat}(K, P)
        @test_throws ArgumentError ValueOf(se, cp_inf(se))
    end
end

@testset "ValueOf throws on -Inf code point" begin
    for K in 3:8, P in 1:(K-1)
        se = Format{SignedFormat,ExtendedFormat}(K, P)
        @test_throws ArgumentError ValueOf(se, cp_neginf(se))
    end
end

@testset "ValueOf throws on out-of-range code points" begin
    for K in 2:6, P in 1:(K-1)
        for fmt in all_formats(K, P)
            @test_throws ArgumentError ValueOf(fmt, -1)
            @test_throws ArgumentError ValueOf(fmt, Int(cp_max(fmt)) + 1)
        end
    end
end

# =========================================================================
# §3  ValueOf — first positive value is val_pos_subnormal_min (P≥2)
#              or val_pos_normal_min (P=1)
# =========================================================================

@testset "ValueOf(cp=1) is the smallest positive value" begin
    for K in 2:8, P in 1:(K-1)
        for fmt in all_formats(K, P)
            v = ValueOf(fmt, 1)
            @test v > 0
            if PrecisionOf(fmt) > 1
                @test v == val_pos_subnormal_min(fmt)
            else
                @test v == val_pos_normal_min(fmt)
            end
        end
    end
end

# =========================================================================
# §4  _decode_positive_half — subnormals are evenly spaced
# =========================================================================

@testset "Subnormal values are evenly spaced" begin
    for K in 3:8, P in 2:(K-1)
        for fmt in all_formats(K, P)
            nsub = nPosSubnormalsOf(fmt)
            nsub >= 2 || continue
            v1 = ValueOf(fmt, 1)
            v2 = ValueOf(fmt, 2)
            ulp = v2 - v1
            for i in 1:nsub
                @test ValueOf(fmt, i) == i * ulp
            end
        end
    end
end

# =========================================================================
# §5  _decode_positive_half — normal values: binade structure
# =========================================================================

@testset "Normal values within a binade are evenly spaced" begin
    for K in 3:8, P in 2:(K-1)
        for fmt in all_formats(K, P)
            m = Int(significand_scale(fmt))  # values per binade
            cp_nmin = Int(cp_pos_normal_min(fmt))
            cp_nmax = Int(cp_pos_normal_max(fmt))

            # check first binade
            if cp_nmin + m - 1 <= cp_nmax
                v_start = ValueOf(fmt, cp_nmin)
                v_next = ValueOf(fmt, cp_nmin + 1)
                ulp = v_next - v_start
                for j in 0:(m-1)
                    cp = cp_nmin + j
                    cp <= cp_nmax || break
                    @test ValueOf(fmt, cp) == v_start + j * ulp
                end
            end
        end
    end
end

@testset "ULP doubles between consecutive binades" begin
    for K in 4:8, P in 2:(K-2)
        for fmt in all_formats(K, P)
            m = Int(significand_scale(fmt))
            cp_nmin = Int(cp_pos_normal_min(fmt))
            cp_nmax = Int(cp_pos_normal_max(fmt))

            # need at least 2 binades
            cp_nmin + 2 * m - 1 <= cp_nmax || continue

            # ULP in binade 0
            ulp0 = ValueOf(fmt, cp_nmin + 1) - ValueOf(fmt, cp_nmin)
            # ULP in binade 1
            ulp1 = ValueOf(fmt, cp_nmin + m + 1) - ValueOf(fmt, cp_nmin + m)

            @test ulp1 == 2 * ulp0
        end
    end
end

# =========================================================================
# §6  Strict monotonicity of positive-half values
# =========================================================================

@testset "Positive-half values strictly increase K ≤ 8" begin
    for K in 2:8, P in 1:(K-1)
        for fmt in all_formats(K, P)
            vals = AllPositiveFiniteValuesOf(fmt)
            for i in 2:length(vals)
                @test vals[i] > vals[i-1]
            end
        end
    end
end

# =========================================================================
# §7  Negative values are negations of positive values (signed)
# =========================================================================

@testset "Negative values mirror positive values" begin
    for K in 3:7, P in 1:(K-1)
        for D in (is_finite, is_extended)
            fmt = Format{SignedFormat,D}(K, P)
            npos = nPosFiniteValuesOf(fmt)
            nneg = nNegFiniteValuesOf(fmt)
            @test npos == nneg

            for i in 1:npos
                vpos = ValueOfOrdinalPos(fmt, i)
                vneg = ValueOfOrdinalNeg(fmt, i)
                @test vneg == -vpos
            end
        end
    end
end

# =========================================================================
# §8  ValueOfOrdinalPos — first and last match boundary values
# =========================================================================

@testset "ValueOfOrdinalPos(1) is smallest positive value" begin
    for K in 2:8, P in 1:(K-1)
        for fmt in all_formats(K, P)
            v = ValueOfOrdinalPos(fmt, 1)
            @test v == ValueOf(fmt, 1)
            @test v > 0
        end
    end
end

@testset "ValueOfOrdinalPos(n) is largest positive finite value" begin
    for K in 2:8, P in 1:(K-1)
        for fmt in all_formats(K, P)
            n = nPosFiniteValuesOf(fmt)
            v = ValueOfOrdinalPos(fmt, n)
            @test v == val_pos_normal_max(fmt)
        end
    end
end

@testset "ValueOfOrdinalPos out-of-range throws" begin
    for K in 3:6, P in 2:(K-1)
        for fmt in all_formats(K, P)
            n = nPosFiniteValuesOf(fmt)
            @test_throws ArgumentError ValueOfOrdinalPos(fmt, 0)
            @test_throws ArgumentError ValueOfOrdinalPos(fmt, n + 1)
        end
    end
end

# =========================================================================
# §9  ValueOfOrdinalNeg — unsigned throws, signed works
# =========================================================================

@testset "ValueOfOrdinalNeg throws for unsigned" begin
    for K in 2:6, P in 1:K
        K - P + 1 >= 1 || continue
        for D in (is_finite, is_extended)
            fmt = Format{UnsignedFormat,D}(K, P)
            @test_throws ArgumentError ValueOfOrdinalNeg(fmt, 1)
        end
    end
end

@testset "ValueOfOrdinalNeg out-of-range throws for signed" begin
    for K in 3:6, P in 1:(K-1)
        for D in (is_finite, is_extended)
            fmt = Format{SignedFormat,D}(K, P)
            n = nNegFiniteValuesOf(fmt)
            @test_throws ArgumentError ValueOfOrdinalNeg(fmt, 0)
            @test_throws ArgumentError ValueOfOrdinalNeg(fmt, n + 1)
        end
    end
end

@testset "ValueOfOrdinalNeg increasing magnitude" begin
    for K in 3:7, P in 1:(K-1)
        for D in (is_finite, is_extended)
            fmt = Format{SignedFormat,D}(K, P)
            n = nNegFiniteValuesOf(fmt)
            n >= 2 || continue
            for i in 2:n
                # magnitude increases
                @test abs(ValueOfOrdinalNeg(fmt, i)) > abs(ValueOfOrdinalNeg(fmt, i - 1))
            end
        end
    end
end

# =========================================================================
# §10  AllFiniteValuesOf — length matches count
# =========================================================================

@testset "AllFiniteValuesOf length == nFiniteValuesOf" begin
    for K in 2:7, P in 1:(K-1)
        for fmt in all_formats(K, P)
            @test length(AllFiniteValuesOf(fmt)) == nFiniteValuesOf(fmt)
        end
    end
end

@testset "AllPositiveFiniteValuesOf length == nPosFiniteValuesOf" begin
    for K in 2:7, P in 1:(K-1)
        for fmt in all_formats(K, P)
            @test length(AllPositiveFiniteValuesOf(fmt)) == nPosFiniteValuesOf(fmt)
        end
    end
end

# =========================================================================
# §11  AllFiniteValuesOf — all values are exact dyadic rationals
# =========================================================================

@testset "All values are dyadic rationals (denominator is power of 2)" begin
    for K in 2:7, P in 1:(K-1)
        for fmt in all_formats(K, P)
            for v in AllFiniteValuesOf(fmt)
                d = denominator(v)
                # d must be a power of 2
                @test d > 0
                @test d & (d - 1) == 0  # power-of-2 check
            end
        end
    end
end

# =========================================================================
# §12  AllPositiveFiniteValuesOf — all positive, no zero
# =========================================================================

@testset "AllPositiveFiniteValuesOf contains only positive values" begin
    for K in 2:7, P in 1:(K-1)
        for fmt in all_formats(K, P)
            for v in AllPositiveFiniteValuesOf(fmt)
                @test v > 0
            end
        end
    end
end

# =========================================================================
# §13  Hand-computed full value ladder: unsigned finite K=3 P=2
# =========================================================================

@testset "Full ladder: unsigned finite K=3 P=2" begin
    fmt = Format{UnsignedFormat,FiniteFormat}(3, 2)
    @test AllFiniteValuesOf(fmt) == Rational{BigInt}[0, 1//4, 1//2, 3//4, 1, 3//2, 2]
    @test AllPositiveFiniteValuesOf(fmt) == Rational{BigInt}[1//4, 1//2, 3//4, 1, 3//2, 2]
end

@testset "Full ladder: unsigned extended K=3 P=2" begin
    fmt = Format{UnsignedFormat,ExtendedFormat}(3, 2)
    # cp 6 = +Inf, so finite values are cp 0..5
    @test AllFiniteValuesOf(fmt) == Rational{BigInt}[0, 1//4, 1//2, 3//4, 1, 3//2]
    @test AllPositiveFiniteValuesOf(fmt) == Rational{BigInt}[1//4, 1//2, 3//4, 1, 3//2]
end

# =========================================================================
# §14  Hand-computed full value ladder: signed extended K=4 P=2
# =========================================================================

@testset "Full ladder: signed extended K=4 P=2" begin
    # K=4, P=2, W=2, B=2
    # Positive: cp 1=sub, cp 2..6=normals
    # Negative: cp 9=sub, cp 10..14=normals
    fmt = Format{SignedFormat,ExtendedFormat}(4, 2)

    pos_vals = AllPositiveFiniteValuesOf(fmt)
    @test length(pos_vals) == nPosFiniteValuesOf(fmt)

    # all positive values are strictly increasing
    for i in 2:length(pos_vals)
        @test pos_vals[i] > pos_vals[i-1]
    end

    # first positive = val_pos_subnormal_min
    @test pos_vals[1] == val_pos_subnormal_min(fmt)

    # last positive = val_pos_normal_max
    @test pos_vals[end] == val_pos_normal_max(fmt)

    # All finite values include zero and both halves
    all_vals = AllFiniteValuesOf(fmt)
    @test all_vals[1] == 0  # zero first
    @test length(all_vals) == nFiniteValuesOf(fmt)
end

# =========================================================================
# §15  ValueOf decode/encode roundtrip: exhaustive for small formats
# =========================================================================

@testset "Decode roundtrip: every finite cp decodes to a unique value K ≤ 7" begin
    for K in 2:7, P in 1:(K-1)
        for fmt in all_formats(K, P)
            vals = AllFiniteValuesOf(fmt)
            # no duplicates (except zero appears once)
            nonzero = filter(v -> v != 0, vals)
            @test length(nonzero) == length(unique(nonzero))
        end
    end
end

# =========================================================================
# §16  Signed format: AllFiniteValues has zero in the middle
# =========================================================================

@testset "Signed AllFiniteValuesOf: negative < 0 < positive" begin
    for K in 3:7, P in 2:(K-1)
        for D in (is_finite, is_extended)
            fmt = Format{SignedFormat,D}(K, P)
            vals = AllFiniteValuesOf(fmt)

            # zero should be present exactly once
            @test count(v -> v == 0, vals) == 1

            # values before zero are non-negative (code-point order: 0 is first)
            @test vals[1] == 0

            # positive values come next, then negative
            pos_count = count(v -> v > 0, vals)
            neg_count = count(v -> v < 0, vals)
            @test pos_count == nPosFiniteValuesOf(fmt)
            @test neg_count == nNegFiniteValuesOf(fmt)
        end
    end
end

# =========================================================================
# §17  ValueOf at binade boundaries: value equals exact power of 2
# =========================================================================

@testset "First normal in each binade is a power of 2" begin
    for K in 3:8, P in 2:(K-1)
        for fmt in all_formats(K, P)
            m = Int(significand_scale(fmt))
            cp_nmin = Int(cp_pos_normal_min(fmt))
            cp_nmax = Int(cp_pos_normal_max(fmt))
            B = Int(ExponentBiasOf(fmt))

            binade = 0
            cp = cp_nmin
            while cp <= cp_nmax
                v = ValueOf(fmt, cp)
                k = binade + 1 - B
                expected = Rational{BigInt}(BigInt(1) << (k + Int(B) - 1)) // (BigInt(1) << (Int(B) - 1))
                # simpler: v should be a power of 2
                num = numerator(v)
                den = denominator(v)
                @test num & (num - 1) == 0  # numerator is power of 2
                @test den & (den - 1) == 0  # denominator is power of 2
                binade += 1
                cp += m
            end
        end
    end
end

# =========================================================================
# §18  Large format: K=16 structural checks
# =========================================================================

@testset "K=16 P=8: basic structural sanity" begin
    for (S, D) in ((is_unsigned, is_finite), (is_signed, is_extended))
        Σ = S === is_signed ? 1 : 0
        fmt = Format{S,D}(16, 8)

        # ValueOf at boundaries
        @test ValueOf(fmt, 0) == 0
        @test ValueOf(fmt, 1) > 0
        @test ValueOf(fmt, 1) == val_pos_subnormal_min(fmt)

        cp_nm = cp_pos_normal_min(fmt)
        @test ValueOf(fmt, cp_nm) == val_pos_normal_min(fmt)

        # Smallest normal > largest subnormal
        cp_sm = cp_pos_subnormal_max(fmt)
        @test ValueOf(fmt, cp_nm) > ValueOf(fmt, cp_sm)

        # Normal min is a power of 2
        v = ValueOf(fmt, cp_nm)
        num = numerator(v)
        den = denominator(v)
        @test num & (num - 1) == 0
        @test den & (den - 1) == 0
    end
end

# =========================================================================
# §19  P=1 formats: no subnormals, all normals
# =========================================================================

@testset "P=1: ValueOf(1) == val_pos_normal_min" begin
    for K in 2:8
        for fmt in all_formats(K, 1)
            @test ValueOf(fmt, 1) == val_pos_normal_min(fmt)
        end
    end
end

@testset "P=1: every positive value is a power of 2" begin
    for K in 2:7
        for fmt in all_formats(K, 1)
            pos = AllPositiveFiniteValuesOf(fmt)
            for v in pos
                num = numerator(v)
                den = denominator(v)
                @test num == 1 || (num & (num - 1) == 0)
                @test den & (den - 1) == 0
            end
        end
    end
end

# =========================================================================
# §20  Minimal formats
# =========================================================================

@testset "Minimal: unsigned finite K=2 P=1" begin
    fmt = Format{UnsignedFormat,FiniteFormat}(2, 1)
    # 4 cps: 0=zero, 1,2=normals, 3=NaN
    # B=2, normals: cp1 → 2^(1-2)=1/2, cp2 → 2^(2-2)=1
    @test AllFiniteValuesOf(fmt) == Rational{BigInt}[0, 1//2, 1]
end

@testset "Minimal: unsigned finite K=2 P=2" begin
    fmt = Format{UnsignedFormat,FiniteFormat}(2, 2)
    # 4 cps: 0=zero, 1=subnormal, 2=normal, 3=NaN
    # B=1, sub: 1*2^(2-2-1)=1/2, normal: cp2 e=1,k=0: 2^0=1
    @test AllFiniteValuesOf(fmt) == Rational{BigInt}[0, 1//2, 1]
end

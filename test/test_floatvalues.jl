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
# §1  val_pos_subnormal_min — P=1 returns nothing, P≥2 returns exact value
# =========================================================================

@testset "val_pos_subnormal_min: P=1 → nothing" begin
    for K in 2:8
        for fmt in all_formats(K, 1)
            @test val_pos_subnormal_min(fmt) === nothing
        end
    end
end

@testset "val_pos_subnormal_min: smallest positive value" begin
    # For P≥2, val_pos_subnormal_min = 2^(2 - P - B)
    # This should equal ValueOf(fmt, 1) — the first positive code point's value
    for K in 3:8, P in 2:(K-1)
        for fmt in all_formats(K, P)
            vmin = val_pos_subnormal_min(fmt)
            @test vmin !== nothing
            @test vmin > 0
            # cross-check with ValueOf at code point 1
            @test vmin == ValueOf(fmt, 1)
        end
    end
end

@testset "val_pos_subnormal_min hand-computed: unsigned finite K=3 P=2" begin
    # B = ExponentBias = 2^(K-P) = 2^1 = 2
    # val = 2^(2-2-2) = 2^(-2) = 1/4
    fmt = Format{UnsignedFormat,FiniteFormat}(3, 2)
    @test val_pos_subnormal_min(fmt) == 1 // 4
end

@testset "val_pos_subnormal_min hand-computed: unsigned finite K=4 P=3" begin
    # B = 2^(4-3) = 2, val = 2^(2-3-2) = 2^(-3) = 1/8
    fmt = Format{UnsignedFormat,FiniteFormat}(4, 3)
    @test val_pos_subnormal_min(fmt) == 1 // 8
end

# =========================================================================
# §2  val_pos_subnormal_max — boundary between subnormal and normal
# =========================================================================

@testset "val_pos_subnormal_max: P=1 → nothing" begin
    for K in 2:8
        for fmt in all_formats(K, 1)
            @test val_pos_subnormal_max(fmt) === nothing
        end
    end
end

@testset "val_pos_subnormal_max: cross-check with ValueOf" begin
    # Should equal ValueOf at the last subnormal code point
    for K in 3:8, P in 2:(K-1)
        for fmt in all_formats(K, P)
            vmax = val_pos_subnormal_max(fmt)
            @test vmax !== nothing
            cp_sub_last = cp_pos_subnormal_max(fmt)
            @test vmax == ValueOf(fmt, cp_sub_last)
        end
    end
end

@testset "val_pos_subnormal_max < val_pos_normal_min" begin
    for K in 3:8, P in 2:(K-1)
        for fmt in all_formats(K, P)
            vsmax = val_pos_subnormal_max(fmt)
            vnmin = val_pos_normal_min(fmt)
            @test vsmax < vnmin
        end
    end
end

@testset "val_pos_subnormal_max hand-computed: unsigned finite K=3 P=2" begin
    # Only 1 subnormal at cp=1, val = 1/4
    # So subnormal_max = subnormal_min = 1/4
    fmt = Format{UnsignedFormat,FiniteFormat}(3, 2)
    @test val_pos_subnormal_max(fmt) == 1 // 4
end

@testset "val_pos_subnormal_max hand-computed: unsigned finite K=4 P=3" begin
    # 3 subnormals at cp 1,2,3 with values 1/8, 2/8, 3/8
    # max = 3/8
    fmt = Format{UnsignedFormat,FiniteFormat}(4, 3)
    @test val_pos_subnormal_max(fmt) == 3 // 8
end

# =========================================================================
# §3  val_pos_normal_min — smallest normal value
# =========================================================================

@testset "val_pos_normal_min: equals 2^(1-B)" begin
    for K in 2:10, P in 1:(K-1)
        for fmt in all_formats(K, P)
            B = ExponentBiasOf(fmt)
            expected = twopow(1 - B)
            @test val_pos_normal_min(fmt) == expected
        end
    end
end

@testset "val_pos_normal_min: cross-check with ValueOf" begin
    for K in 3:8, P in 1:(K-1)
        for fmt in all_formats(K, P)
            cp_nm = cp_pos_normal_min(fmt)
            @test val_pos_normal_min(fmt) == ValueOf(fmt, cp_nm)
        end
    end
end

@testset "val_pos_normal_min hand-computed: unsigned finite K=3 P=2" begin
    # B=2, val = 2^(1-2) = 1/2
    fmt = Format{UnsignedFormat,FiniteFormat}(3, 2)
    @test val_pos_normal_min(fmt) == 1 // 2
end

# =========================================================================
# §4  val_pos_normal_max — largest positive finite value
# =========================================================================

@testset "val_pos_normal_max: cross-check with ValueOf" begin
    for K in 3:8, P in 1:(K-1)
        for fmt in all_formats(K, P)
            cp_nm = cp_pos_normal_max(fmt)
            @test val_pos_normal_max(fmt) == ValueOf(fmt, cp_nm)
        end
    end
end

@testset "val_pos_normal_max hand-computed: unsigned finite K=3 P=2" begin
    # Normals: cp 2=1/2, 3=3/4, 4=1, 5=3/2, 6=2  → max = 2
    # Formula: (2^P - 1) * 2^(B - P) = (4-1)*2^(2-2) = 3*1 = 3
    # Wait — let me re-derive.  B=2, P=2.
    # cp 6: e = 6÷2 = 3, t = 0, k = 3-2 = 1, val = 2^1 + 0 = 2
    # Formula gives (2^2 - 1)*2^(2-2) = 3*1 = 3  — that's wrong!
    # Actually the last normal cp is 6 and its value is 2.  Let's just check ValueOf.
    fmt = Format{UnsignedFormat,FiniteFormat}(3, 2)
    @test val_pos_normal_max(fmt) == ValueOf(fmt, cp_pos_normal_max(fmt))
end

@testset "val_pos_normal_max: unsigned finite vs signed finite same (K,P)" begin
    # For finite formats: formula is the same for unsigned and signed
    for K in 3:8, P in 2:(K-1)
        uf = Format{UnsignedFormat,FiniteFormat}(K, P)
        sf = Format{SignedFormat,FiniteFormat}(K, P)
        # They may have different exponent biases, so values differ.
        # But each should match its own ValueOf.
        @test val_pos_normal_max(uf) == ValueOf(uf, cp_pos_normal_max(uf))
        @test val_pos_normal_max(sf) == ValueOf(sf, cp_pos_normal_max(sf))
    end
end

@testset "val_pos_normal_max: extended < finite (same K,P)" begin
    # Extended formats sacrifice one code point for +Inf, so max finite is smaller
    for K in 3:8, P in 2:(K-1)
        for S in (is_unsigned, is_signed)
            Σ = S === is_signed ? 1 : 0
            K - P + 1 - Σ >= 1 || continue
            ff = Format{S,FiniteFormat}(K, P)
            fe = Format{S,ExtendedFormat}(K, P)
            @test val_pos_normal_max(fe) < val_pos_normal_max(ff)
        end
    end
end

# =========================================================================
# §5  val_pos_normal_max — all four variants hand-computed K=4 P=2
# =========================================================================

@testset "val_pos_normal_max all variants K=4 P=2" begin
    # unsigned finite K=4 P=2: W=3, B=2^2=4
    # cp_pos_normal_max = 2^4 - 2 = 14
    # e=14÷2=7, t=0, k=7-4=3, val = 2^3 = 8
    uf = Format{UnsignedFormat,FiniteFormat}(4, 2)
    @test val_pos_normal_max(uf) == ValueOf(uf, cp_pos_normal_max(uf))

    # unsigned extended K=4 P=2: cp_pos_normal_max = 2^4 - 3 = 13
    # e=13÷2=6, t=1, k=6-4=2, val = 2^2 + 1*2^(2+1-2) = 4+2 = 6
    ue = Format{UnsignedFormat,ExtendedFormat}(4, 2)
    @test val_pos_normal_max(ue) == ValueOf(ue, cp_pos_normal_max(ue))

    # signed finite K=4 P=2: W=1, B=2^0=1
    # cp_nan = 8, cp_pos_normal_max = 8-1 = 7
    # e=7÷2=3, t=1, k=3-1=2, val = 2^2 + 1*2^(2+1-2) = 4+2 = 6
    sf = Format{SignedFormat,FiniteFormat}(4, 2)
    @test val_pos_normal_max(sf) == ValueOf(sf, cp_pos_normal_max(sf))

    # signed extended K=4 P=2: cp_pos_normal_max = 8-2 = 6
    # e=6÷2=3, t=0, k=3-1=2, wait B=2^(K-P-1)=2^1=2, k=3-2=1
    # val = 2^1 + 0 = 2
    se = Format{SignedFormat,ExtendedFormat}(4, 2)
    @test val_pos_normal_max(se) == ValueOf(se, cp_pos_normal_max(se))
end

# =========================================================================
# §6  val_pos_normal_max — P=1 edge case (all normals, no trailing bits)
# =========================================================================

@testset "val_pos_normal_max P=1: each binade has exactly one value" begin
    for K in 2:8
        for (S, D) in ((is_unsigned, is_finite), (is_unsigned, is_extended),
            (is_signed, is_finite), (is_signed, is_extended))
            Σ = S === is_signed ? 1 : 0
            W = K - 1 + 1 - Σ
            W >= 1 || continue
            fmt = Format{S,D}(K, 1)
            @test val_pos_normal_max(fmt) == ValueOf(fmt, cp_pos_normal_max(fmt))
        end
    end
end

# =========================================================================
# §7  val_ordinal_ith_pos_subnormal — ordinal accessor
# =========================================================================

@testset "val_ordinal_ith_pos_subnormal: P=1 → nothing" begin
    for K in 2:8
        for fmt in all_formats(K, 1)
            @test val_ordinal_ith_pos_subnormal(fmt, 1) === nothing
        end
    end
end

@testset "val_ordinal_ith_pos_subnormal: matches ValueOf" begin
    for K in 3:8, P in 2:(K-1)
        for fmt in all_formats(K, P)
            nsub = nPosSubnormalsOf(fmt)
            for i in 1:nsub
                @test val_ordinal_ith_pos_subnormal(fmt, i) == ValueOf(fmt, i)
            end
        end
    end
end

@testset "val_ordinal_ith_pos_subnormal: evenly spaced" begin
    # Subnormals are evenly spaced: val(i) = i * ULP_sub
    for K in 3:8, P in 2:(K-1)
        for fmt in all_formats(K, P)
            nsub = nPosSubnormalsOf(fmt)
            nsub >= 2 || continue
            v1 = val_ordinal_ith_pos_subnormal(fmt, 1)
            v2 = val_ordinal_ith_pos_subnormal(fmt, 2)
            ulp_sub = v2 - v1
            for i in 2:nsub
                vi = val_ordinal_ith_pos_subnormal(fmt, i)
                @test vi == i * ulp_sub
            end
        end
    end
end

# =========================================================================
# §8  val_cardinal_ith_pos_subnormal — 0-based ordinal
# =========================================================================

@testset "val_cardinal_ith_pos_subnormal: cardinal = ordinal - 1" begin
    for K in 3:8, P in 2:(K-1)
        for fmt in all_formats(K, P)
            nsub = nPosSubnormalsOf(fmt)
            for i in 1:nsub
                @test val_cardinal_ith_pos_subnormal(fmt, i - 1) == val_ordinal_ith_pos_subnormal(fmt, i)
            end
        end
    end
end

# =========================================================================
# §9  Negative value boundaries via symmetry
# =========================================================================

@testset "val_neg = -val_pos for signed formats (via ValueOf)" begin
    for K in 3:7, P in 2:(K-1)
        for D in (is_finite, is_extended)
            fmt = Format{SignedFormat,D}(K, P)

            # smallest negative subnormal = -val_pos_subnormal_min
            vmin = val_pos_subnormal_min(fmt)
            @test ValueOf(fmt, cp_neg_subnormal_min(fmt)) == -vmin

            # largest negative subnormal = -val_pos_subnormal_max
            vmax = val_pos_subnormal_max(fmt)
            @test ValueOf(fmt, cp_neg_subnormal_max(fmt)) == -vmax

            # smallest negative normal = -val_pos_normal_min
            @test ValueOf(fmt, cp_neg_normal_min(fmt)) == -val_pos_normal_min(fmt)

            # largest negative normal = -val_pos_normal_max
            @test ValueOf(fmt, cp_neg_normal_max(fmt)) == -val_pos_normal_max(fmt)
        end
    end
end

# =========================================================================
# §10  Full value ladder: hand-computed for Format{unsigned,finite}(3,2)
# =========================================================================

@testset "Complete value ladder: unsigned finite K=3 P=2" begin
    fmt = Format{UnsignedFormat,FiniteFormat}(3, 2)
    # cp 0: 0, cp 1: 1/4, cp 2: 1/2, cp 3: 3/4, cp 4: 1, cp 5: 3/2, cp 6: 2
    expected = [0 // 1, 1 // 4, 1 // 2, 3 // 4, 1 // 1, 3 // 2, 2 // 1]
    actual = AllFiniteValuesOf(fmt)
    @test actual == expected
end

@testset "Complete value ladder: unsigned finite K=4 P=3" begin
    fmt = Format{UnsignedFormat,FiniteFormat}(4, 3)
    # B=2, subnormal ULP = 2^(2-3-2) = 2^(-3) = 1/8
    # cp 0: 0
    # cp 1: 1/8,  cp 2: 2/8=1/4,  cp 3: 3/8  (subnormals)
    # cp 4: e=1,k=-1: 1/2           cp 5: 1/2+1/8=5/8
    #   cp 6: 1/2+2/8=3/4           cp 7: 1/2+3/8=7/8
    # cp 8: e=2,k=0: 1              cp 9: 1+1/4=5/4
    #   cp 10: 1+2/4=3/2            cp 11: 1+3/4=7/4
    # cp 12: e=3,k=1: 2             cp 13: 2+1/2=5/2
    #   cp 14: 2+2/2=3
    # cp 15: NaN
    expected = Rational{BigInt}[
        0, 1//8, 1//4, 3//8,
        1//2, 5//8, 3//4, 7//8,
        1, 5//4, 3//2, 7//4,
        2, 5//2, 3
    ]
    actual = AllFiniteValuesOf(fmt)
    @test actual == expected
end

# =========================================================================
# §11  Consistency: val functions agree with each other
# =========================================================================

@testset "val_pos_subnormal_min == val_ordinal_ith_pos_subnormal(1)" begin
    for K in 3:8, P in 2:(K-1)
        for fmt in all_formats(K, P)
            @test val_pos_subnormal_min(fmt) == val_ordinal_ith_pos_subnormal(fmt, 1)
        end
    end
end

@testset "val_pos_subnormal_max == val_ordinal_ith_pos_subnormal(nPosSubnormals)" begin
    for K in 3:8, P in 2:(K-1)
        for fmt in all_formats(K, P)
            n = nPosSubnormalsOf(fmt)
            @test val_pos_subnormal_max(fmt) == val_ordinal_ith_pos_subnormal(fmt, n)
        end
    end
end

# =========================================================================
# §12  val_pos_normal_max monotonically increases with K (fixed P, S, D)
# =========================================================================

@testset "val_pos_normal_max increases with K" begin
    for P in 2:4
        for (S, D) in ((is_unsigned, is_finite), (is_signed, is_extended))
            prev = nothing
            for K in (P+1):10
                Σ = S === is_signed ? 1 : 0
                K - P + 1 - Σ >= 1 || continue
                fmt = Format{S,D}(K, P)
                v = val_pos_normal_max(fmt)
                if prev !== nothing
                    @test v > prev
                end
                prev = v
            end
        end
    end
end

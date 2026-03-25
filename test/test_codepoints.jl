using Test
using Float3109s

# =========================================================================
# Helper: build all four format variants for a given (K, P)
# =========================================================================

function all_formats(K, P)
    fmts = []
    for (S, D) in ((UnsignedFormat, FiniteFormat), (UnsignedFormat, ExtendedFormat),
        (SignedFormat, FiniteFormat), (SignedFormat, ExtendedFormat))
        Σ = S === SignedFormat ? 1 : 0
        W = K - P + 1 - Σ
        W >= 1 || continue
        push!(fmts, Format{S,D}(K, P))
    end
    fmts
end

# =========================================================================
# §1  typeofcp — correct UInt type for each bitwidth
# =========================================================================

@testset "typeofcp returns correct unsigned integer type" begin
    for K in 1:8
        fmt = Format{UnsignedFormat,FiniteFormat}(K, 1)
        @test cp_zero(fmt) isa UInt8
        @test cp_max(fmt) isa UInt8
    end
    for K in 9:16
        W = K - 1 + 1
        W >= 1 || continue
        fmt = Format{UnsignedFormat,FiniteFormat}(K, 1)
        @test cp_zero(fmt) isa UInt16
        @test cp_max(fmt) isa UInt16
    end
    # K=32, P=1 → UInt32
    fmt32 = Format{UnsignedFormat,FiniteFormat}(32, 1)
    @test cp_zero(fmt32) isa UInt32
    @test cp_max(fmt32) isa UInt32

    # K=64, P=1 → UInt64
    fmt64 = Format{UnsignedFormat,FiniteFormat}(64, 1)
    @test cp_zero(fmt64) isa UInt64
end

# =========================================================================
# §2  cp_min, cp_max, cp_zero — basic range
# =========================================================================

@testset "cp_min / cp_max / cp_zero basic properties" begin
    for K in 2:10
        for fmt in all_formats(K, 1)
            @test cp_min(fmt) == 0
            @test cp_zero(fmt) == 0
            @test cp_max(fmt) == twopow(K) - 1
            @test cp_min(fmt) <= cp_max(fmt)
        end
    end
end

# =========================================================================
# §3  cp_nan — location depends on signedness
# =========================================================================

@testset "cp_nan location" begin
    # unsigned: NaN is at cp_max (last code point)
    for K in 2:10, P in 1:K
        K - P + 1 >= 1 || continue
        for D in (FiniteFormat, ExtendedFormat)
            fmt = Format{UnsignedFormat,D}(K, P)
            @test cp_nan(fmt) == cp_max(fmt)
        end
    end

    # signed: NaN is at the midpoint (sign_half_offset)
    for K in 3:10, P in 1:(K-1)
        for D in (FiniteFormat, ExtendedFormat)
            fmt = Format{SignedFormat,D}(K, P)
            @test cp_nan(fmt) == twopow(K - 1)
        end
    end
end

# =========================================================================
# §4  cp_inf / cp_neginf — presence and location
# =========================================================================

@testset "cp_inf / cp_neginf for finite-domain → nothing" begin
    for K in 3:10, P in 1:(K-1)
        for S in (UnsignedFormat, SignedFormat)
            Σ = S === SignedFormat ? 1 : 0
            K - P + 1 - Σ >= 1 || continue
            fmt = Format{S,FiniteFormat}(K, P)
            @test cp_inf(fmt) === nothing
            @test cp_neginf(fmt) === nothing
        end
    end
end

@testset "cp_inf for unsigned extended" begin
    for K in 2:10, P in 1:K
        K - P + 1 >= 1 || continue
        fmt = Format{UnsignedFormat,ExtendedFormat}(K, P)
        @test cp_inf(fmt) == cp_nan(fmt) - 1
        @test cp_neginf(fmt) === nothing  # unsigned has no -Inf
    end
end

@testset "cp_inf / cp_neginf for signed extended" begin
    for K in 3:10, P in 1:(K-1)
        fmt = Format{SignedFormat,ExtendedFormat}(K, P)
        @test cp_inf(fmt) == cp_nan(fmt) - 1
        @test cp_neginf(fmt) == cp_max(fmt)
    end
end

@testset "cp_inf and cp_neginf are distinct from each other and from NaN" begin
    for K in 3:10, P in 1:(K-1)
        fmt = Format{SignedFormat,ExtendedFormat}(K, P)
        @test cp_inf(fmt) != cp_nan(fmt)
        @test cp_neginf(fmt) != cp_nan(fmt)
        @test cp_inf(fmt) != cp_neginf(fmt)
    end
end

# =========================================================================
# §5  Hand-computed special code points for Format{unsigned, finite}(3, 2)
# =========================================================================

@testset "Hand-computed: unsigned finite K=3 P=2 specials" begin
    # 8 code points: 0=zero, 1=sub, 2..6=normals, 7=NaN
    fmt = Format{UnsignedFormat,FiniteFormat}(3, 2)
    @test cp_zero(fmt) == 0x00
    @test cp_nan(fmt) == 0x07
    @test cp_inf(fmt) === nothing
    @test cp_neginf(fmt) === nothing
    @test cp_min(fmt) == 0x00
    @test cp_max(fmt) == 0x07
end

@testset "Hand-computed: unsigned extended K=3 P=2 specials" begin
    # 8 code points: 0=zero, 1=sub, 2..5=normals, 6=+Inf, 7=NaN
    fmt = Format{UnsignedFormat,ExtendedFormat}(3, 2)
    @test cp_nan(fmt) == 0x07
    @test cp_inf(fmt) == 0x06
    @test cp_neginf(fmt) === nothing
end

@testset "Hand-computed: signed extended K=4 P=2 specials" begin
    # 16 code points: 0=zero, 1=sub, 2..6=normals, 7=+Inf
    # 8=NaN, 9=neg sub, 10..14=neg normals, 15=-Inf
    fmt = Format{SignedFormat,ExtendedFormat}(4, 2)
    @test cp_zero(fmt) == 0x00
    @test cp_nan(fmt) == 0x08
    @test cp_inf(fmt) == 0x07
    @test cp_neginf(fmt) == 0x0f
end

# =========================================================================
# §6  Subnormal code point boundaries
# =========================================================================

@testset "cp_pos_subnormal_min / max — P=1 gives nothing" begin
    for K in 2:10
        for fmt in all_formats(K, 1)
            @test cp_pos_subnormal_min(fmt) === nothing
            @test cp_pos_subnormal_max(fmt) === nothing
        end
    end
end

@testset "cp_pos_subnormal_min / max — P≥2" begin
    for K in 3:10, P in 2:(K-1)
        for fmt in all_formats(K, P)
            @test cp_pos_subnormal_min(fmt) == 1
            @test cp_pos_subnormal_max(fmt) == twopow(P - 1) - 1
            # subnormal range is contiguous
            @test cp_pos_subnormal_max(fmt) == cp_pos_normal_min(fmt) - 1
        end
    end
end

@testset "cp_neg_subnormal — unsigned gives nothing" begin
    for K in 3:10, P in 2:K
        K - P + 1 >= 1 || continue
        for D in (FiniteFormat, ExtendedFormat)
            fmt = Format{UnsignedFormat,D}(K, P)
            @test cp_neg_subnormal_min(fmt) === nothing
            @test cp_neg_subnormal_max(fmt) === nothing
        end
    end
end

@testset "cp_neg_subnormal — signed, P=1 gives nothing" begin
    for K in 3:10
        for D in (FiniteFormat, ExtendedFormat)
            fmt = Format{SignedFormat,D}(K, 1)
            @test cp_neg_subnormal_min(fmt) === nothing
            @test cp_neg_subnormal_max(fmt) === nothing
        end
    end
end

@testset "cp_neg_subnormal — signed, P≥2" begin
    for K in 4:10, P in 2:(K-1)
        for D in (FiniteFormat, ExtendedFormat)
            fmt = Format{SignedFormat,D}(K, P)
            nsmin = cp_neg_subnormal_min(fmt)
            nsmax = cp_neg_subnormal_max(fmt)
            @test nsmin !== nothing
            @test nsmax !== nothing
            # least-magnitude negative subnormal is right after NaN
            @test nsmin == cp_nan(fmt) + 1
            # range is contiguous: nsmax == cp_neg_normal_min - 1
            @test nsmax == cp_neg_normal_min(fmt) - 1
        end
    end
end

# =========================================================================
# §7  Normal code point boundaries
# =========================================================================

@testset "cp_pos_normal_min equals twopow(P-1)" begin
    for K in 2:10, P in 1:(K-1)
        for fmt in all_formats(K, P)
            @test cp_pos_normal_min(fmt) == twopow(P - 1)
        end
    end
end

@testset "cp_pos_normal_max hand-computed K=3 P=2" begin
    # unsigned finite: normals at cp 2..6 → max = 6 = 2^3 - 2
    @test cp_pos_normal_max(Format{UnsignedFormat,FiniteFormat}(3, 2)) == 6

    # unsigned extended: normals at cp 2..5 → max = 5 = 2^3 - 3
    @test cp_pos_normal_max(Format{UnsignedFormat,ExtendedFormat}(3, 2)) == 5
end

@testset "cp_pos_normal_max hand-computed signed K=4 P=2" begin
    # signed finite: normals at cp 2..7 → max = 7 = NaN-1
    @test cp_pos_normal_max(Format{SignedFormat,FiniteFormat}(4, 2)) == 7

    # signed extended: normals at cp 2..6 → max = 6 = NaN-2
    @test cp_pos_normal_max(Format{SignedFormat,ExtendedFormat}(4, 2)) == 6
end

@testset "cp_neg_normal_min / max — unsigned gives nothing" begin
    for K in 2:10, P in 1:K
        K - P + 1 >= 1 || continue
        for D in (FiniteFormat, ExtendedFormat)
            fmt = Format{UnsignedFormat,D}(K, P)
            @test cp_neg_normal_min(fmt) === nothing
            @test cp_neg_normal_max(fmt) === nothing
        end
    end
end

@testset "cp_neg_normal_max for signed formats" begin
    for K in 3:10, P in 1:(K-1)
        sf = Format{SignedFormat,FiniteFormat}(K, P)
        @test cp_neg_normal_max(sf) == cp_max(sf)

        se = Format{SignedFormat,ExtendedFormat}(K, P)
        @test cp_neg_normal_max(se) == cp_max(se) - 1  # last cp is -Inf
    end
end

# =========================================================================
# §8  Ordinal subnormal/normal code point accessors
# =========================================================================

@testset "cp_ordinal_ith_pos_subnormal sequential" begin
    for K in 3:8, P in 2:(K-1)
        for fmt in all_formats(K, P)
            nsub = nPosSubnormalsOf(fmt)
            nsub > 0 || continue
            for i in 1:nsub
                cp = cp_ordinal_ith_pos_subnormal(fmt, i)
                @test cp == i  # subnormals are at code points 1, 2, ..., nsub
            end
            # out-of-range throws
            @test_throws ArgumentError cp_ordinal_ith_pos_subnormal(fmt, 0)
            @test_throws ArgumentError cp_ordinal_ith_pos_subnormal(fmt, nsub + 1)
        end
    end
end

@testset "cp_ordinal_ith_pos_normal sequential" begin
    for K in 3:8, P in 1:(K-1)
        for fmt in all_formats(K, P)
            nnorm = nPosNormalsOf(fmt)
            nnorm > 0 || continue
            cp_first = cp_pos_normal_min(fmt)
            for i in 1:min(nnorm, 10)  # test first 10
                cp = cp_ordinal_ith_pos_normal(fmt, i)
                @test cp == cp_first + (i - 1)
            end
            # last one
            @test cp_ordinal_ith_pos_normal(fmt, nnorm) == cp_pos_normal_max(fmt)
            # out-of-range throws
            @test_throws ArgumentError cp_ordinal_ith_pos_normal(fmt, 0)
            @test_throws ArgumentError cp_ordinal_ith_pos_normal(fmt, nnorm + 1)
        end
    end
end

# =========================================================================
# §9  cp_positive_max
# =========================================================================

@testset "cp_positive_max == cp_nan - 1" begin
    for K in 2:10, P in 1:(K-1)
        for fmt in all_formats(K, P)
            @test cp_positive_max(fmt) == cp_nan(fmt) - 1
        end
    end
end

# =========================================================================
# §10  cp_is_positive / cp_is_nonnegative / cp_is_negative
# =========================================================================

@testset "cp_is_positive / nonnegative / negative classification" begin
    for K in 3:7, P in 1:(K-1)
        for fmt in all_formats(K, P)
            @testset "K=$K P=$P $(typeof(fmt))" begin
                # zero is nonneg but not positive
                @test cp_is_nonnegative(fmt, cp_zero(fmt))
                @test !cp_is_positive(fmt, cp_zero(fmt))

                # NaN is not positive and not nonnegative
                @test !cp_is_positive(fmt, cp_nan(fmt))
                @test !cp_is_nonnegative(fmt, cp_nan(fmt))

                # first positive subnormal (or normal if P=1)
                cp1 = PrecisionOf(fmt) > 1 ? cp_pos_subnormal_min(fmt) : cp_pos_normal_min(fmt)
                @test cp_is_positive(fmt, cp1)
                @test cp_is_nonnegative(fmt, cp1)
                @test !cp_is_negative(fmt, cp1)
            end
        end
    end
end

@testset "cp_is_negative — unsigned always false" begin
    for K in 2:8, P in 1:K
        K - P + 1 >= 1 || continue
        for D in (FiniteFormat, ExtendedFormat)
            fmt = Format{UnsignedFormat,D}(K, P)
            for cp in 0:Int(cp_max(fmt))
                @test cp_is_negative(fmt, cp) == false
            end
        end
    end
end

@testset "cp_is_negative — signed, only code points > NaN" begin
    for K in 3:7, P in 1:(K-1)
        for D in (FiniteFormat, ExtendedFormat)
            fmt = Format{SignedFormat,D}(K, P)
            nan = cp_nan(fmt)
            for cp in 0:Int(cp_max(fmt))
                if cp > nan
                    @test cp_is_negative(fmt, cp)
                else
                    @test !cp_is_negative(fmt, cp)
                end
            end
        end
    end
end

# =========================================================================
# §11  pos_cp_to_neg_cp / neg_cp_to_pos_cp roundtrip (signed only)
# =========================================================================

@testset "pos ↔ neg code point roundtrip (signed)" begin
    for K in 3:7, P in 1:(K-1)
        for D in (FiniteFormat, ExtendedFormat)
            fmt = Format{SignedFormat,D}(K, P)
            # walk all positive code points
            for cp in 1:(Int(cp_nan(fmt))-1)
                cp_is_positive(fmt, cp) || continue
                # skip +Inf code point
                inf_val = cp_inf(fmt)
                inf_val !== nothing && cp == inf_val && continue

                neg = pos_cp_to_neg_cp(fmt, cp)
                @test cp_is_negative(fmt, neg)
                back = neg_cp_to_pos_cp(fmt, neg)
                @test back == cp
            end
        end
    end
end

@testset "pos_cp_to_neg_cp rejects non-positive" begin
    fmt = Format{SignedFormat,FiniteFormat}(4, 2)
    @test_throws ArgumentError pos_cp_to_neg_cp(fmt, 0)            # zero
    @test_throws ArgumentError pos_cp_to_neg_cp(fmt, cp_nan(fmt))  # NaN
end

# =========================================================================
# §12  cp_changesign — involution and fixed points (signed only)
# =========================================================================

@testset "cp_changesign is an involution" begin
    for K in 3:7, P in 1:(K-1)
        for D in (FiniteFormat, ExtendedFormat)
            fmt = Format{SignedFormat,D}(K, P)
            for cp in 0:Int(cp_max(fmt))
                cp == cp_nan(fmt) && (
                    @test cp_changesign(fmt, cp) == cp;
                    continue
                )
                cp == cp_zero(fmt) && (
                    @test cp_changesign(fmt, cp) == cp;
                    continue
                )
                flipped = cp_changesign(fmt, cp)
                @test cp_changesign(fmt, flipped) == cp
            end
        end
    end
end

@testset "cp_changesign fixed points: zero and NaN" begin
    for K in 3:8, P in 1:(K-1)
        for D in (FiniteFormat, ExtendedFormat)
            fmt = Format{SignedFormat,D}(K, P)
            @test cp_changesign(fmt, cp_zero(fmt)) == cp_zero(fmt)
            @test cp_changesign(fmt, cp_nan(fmt)) == cp_nan(fmt)
        end
    end
end

@testset "cp_changesign swaps positive ↔ negative" begin
    for K in 3:7, P in 1:(K-1)
        for D in (FiniteFormat, ExtendedFormat)
            fmt = Format{SignedFormat,D}(K, P)
            for cp in 1:(Int(cp_nan(fmt))-1)
                @test cp_is_positive(fmt, cp)
                flipped = cp_changesign(fmt, cp)
                @test cp_is_negative(fmt, flipped)
            end
        end
    end
end

# =========================================================================
# §13  Code point ordering: subnormals < normals < specials
# =========================================================================

@testset "Code point ordering within positive half" begin
    for K in 3:8, P in 2:(K-1)
        for fmt in all_formats(K, P)
            @testset "K=$K P=$P $(typeof(fmt))" begin
                # zero < first subnormal < last subnormal < first normal < last normal
                @test cp_zero(fmt) < cp_pos_subnormal_min(fmt)
                @test cp_pos_subnormal_min(fmt) <= cp_pos_subnormal_max(fmt)
                @test cp_pos_subnormal_max(fmt) < cp_pos_normal_min(fmt)
                @test cp_pos_normal_min(fmt) <= cp_pos_normal_max(fmt)
                @test cp_pos_normal_max(fmt) < cp_nan(fmt)

                # for extended formats, +Inf sits between last normal and NaN
                if is_extended(fmt)
                    @test cp_pos_normal_max(fmt) < cp_inf(fmt)
                    @test cp_inf(fmt) < cp_nan(fmt)
                end
            end
        end
    end
end

@testset "Code point ordering within negative half (signed, P≥2)" begin
    for K in 4:8, P in 2:(K-1)
        for D in (FiniteFormat, ExtendedFormat)
            fmt = Format{SignedFormat,D}(K, P)
            @testset "K=$K P=$P $(D)" begin
                # NaN < neg sub min < neg sub max < neg normal min < neg normal max
                @test cp_nan(fmt) < cp_neg_subnormal_min(fmt)
                @test cp_neg_subnormal_min(fmt) <= cp_neg_subnormal_max(fmt)
                @test cp_neg_subnormal_max(fmt) < cp_neg_normal_min(fmt)
                @test cp_neg_normal_min(fmt) <= cp_neg_normal_max(fmt)

                # for extended, -Inf is last
                if is_extended(fmt)
                    @test cp_neg_normal_max(fmt) < cp_neginf(fmt)
                    @test cp_neginf(fmt) == cp_max(fmt)
                end
            end
        end
    end
end

# =========================================================================
# §14  Contiguity: no gaps between subnormal max and normal min
# =========================================================================

@testset "Positive subnormal/normal boundary is contiguous" begin
    for K in 3:10, P in 2:(K-1)
        for fmt in all_formats(K, P)
            @test cp_pos_subnormal_max(fmt) + 1 == cp_pos_normal_min(fmt)
        end
    end
end

@testset "Negative subnormal/normal boundary is contiguous (signed, P≥2)" begin
    for K in 4:10, P in 2:(K-1)
        for D in (FiniteFormat, ExtendedFormat)
            fmt = Format{SignedFormat,D}(K, P)
            @test cp_neg_subnormal_max(fmt) + 1 == cp_neg_normal_min(fmt)
        end
    end
end

# =========================================================================
# §15  Count consistency: number of code points in each range matches counts.jl
# =========================================================================

@testset "Code point range sizes match count functions K ≤ 8" begin
    for K in 3:8, P in 2:(K-1)
        for fmt in all_formats(K, P)
            @testset "K=$K P=$P $(typeof(fmt))" begin
                # positive subnormals
                nsub = Int(cp_pos_subnormal_max(fmt)) - Int(cp_pos_subnormal_min(fmt)) + 1
                @test nsub == nPosSubnormalsOf(fmt)

                # positive normals
                nnor = Int(cp_pos_normal_max(fmt)) - Int(cp_pos_normal_min(fmt)) + 1
                @test nnor == nPosNormalsOf(fmt)

                # negative (signed only)
                if is_signed(fmt) && P > 1
                    nsub_neg = Int(cp_neg_subnormal_max(fmt)) - Int(cp_neg_subnormal_min(fmt)) + 1
                    @test nsub_neg == nNegSubnormalsOf(fmt)

                    nnor_neg = Int(cp_neg_normal_max(fmt)) - Int(cp_neg_normal_min(fmt)) + 1
                    @test nnor_neg == nNegNormalsOf(fmt)
                end
            end
        end
    end
end

# =========================================================================
# §16  Hand-computed exhaustive walkthrough: signed extended K=4 P=2
# =========================================================================

@testset "Exhaustive K=4 P=2 signed extended code point map" begin
    fmt = Format{SignedFormat,ExtendedFormat}(4, 2)

    # positive half
    @test cp_zero(fmt) == 0
    @test cp_pos_subnormal_min(fmt) == 1
    @test cp_pos_subnormal_max(fmt) == 1   # twopow(1)-1 = 1
    @test cp_pos_normal_min(fmt) == 2       # twopow(1) = 2
    @test cp_pos_normal_max(fmt) == 6       # NaN(8) - 2
    @test cp_inf(fmt) == 7
    @test cp_nan(fmt) == 8

    # negative half
    @test cp_neg_subnormal_min(fmt) == 9    # NaN + 1
    @test cp_neg_subnormal_max(fmt) == 9    # 8 + 1 = 9
    @test cp_neg_normal_min(fmt) == 10
    @test cp_neg_normal_max(fmt) == 14      # cp_max - 1 = 14
    @test cp_neginf(fmt) == 15
    @test cp_max(fmt) == 15
end

# =========================================================================
# §17  Hand-computed exhaustive walkthrough: unsigned finite K=4 P=3
# =========================================================================

@testset "Exhaustive K=4 P=3 unsigned finite code point map" begin
    # K=4, P=3, W=4-3+1=2, ExponentBias=2^(4-3)=2
    # significand_scale = twopow(2) = 4
    # 16 code points: 0=zero, 1..3=subnormals, 4..14=normals, 15=NaN
    fmt = Format{UnsignedFormat,FiniteFormat}(4, 3)

    @test cp_zero(fmt) == 0
    @test cp_pos_subnormal_min(fmt) == 1
    @test cp_pos_subnormal_max(fmt) == 3   # twopow(2)-1
    @test cp_pos_normal_min(fmt) == 4       # twopow(2)
    @test cp_pos_normal_max(fmt) == 14      # 2^4 - 2
    @test cp_nan(fmt) == 15

    @test nPosSubnormalsOf(fmt) == 3
    @test nPosNormalsOf(fmt) == 11
end

# =========================================================================
# §18  Minimal format edge cases
# =========================================================================

@testset "Minimal: unsigned finite K=2 P=1 (no subnormals)" begin
    # 4 code points: 0=zero, 1,2=normals, 3=NaN
    fmt = Format{UnsignedFormat,FiniteFormat}(2, 1)
    @test cp_zero(fmt) == 0
    @test cp_pos_subnormal_min(fmt) === nothing
    @test cp_pos_subnormal_max(fmt) === nothing
    @test cp_pos_normal_min(fmt) == 1   # twopow(0) = 1
    @test cp_pos_normal_max(fmt) == 2   # 2^2 - 2
    @test cp_nan(fmt) == 3
end

@testset "Minimal: unsigned finite K=2 P=2 (1 binade)" begin
    # 4 code points: 0=zero, 1=subnormal, 2=normal, 3=NaN
    fmt = Format{UnsignedFormat,FiniteFormat}(2, 2)
    @test cp_pos_subnormal_min(fmt) == 1
    @test cp_pos_subnormal_max(fmt) == 1   # twopow(1)-1
    @test cp_pos_normal_min(fmt) == 2
    @test cp_pos_normal_max(fmt) == 2       # 2^2-2
    @test cp_nan(fmt) == 3
end

@testset "Minimal: signed finite K=3 P=1" begin
    # 8 code points: 0=zero, 1,2,3=pos normals, 4=NaN, 5,6,7=neg normals
    fmt = Format{SignedFormat,FiniteFormat}(3, 1)
    @test cp_zero(fmt) == 0
    @test cp_pos_subnormal_min(fmt) === nothing
    @test cp_pos_normal_min(fmt) == 1
    @test cp_pos_normal_max(fmt) == 3   # NaN - 1
    @test cp_nan(fmt) == 4
    @test cp_neg_subnormal_min(fmt) === nothing
    @test cp_neg_normal_min(fmt) == 5   # NaN + 1
    @test cp_neg_normal_max(fmt) == 7   # cp_max
end


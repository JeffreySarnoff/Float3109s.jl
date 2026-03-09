using Test
using Float3109s

# =========================================================================
# Helper: build all four format variants for a given (K, P)
# =========================================================================

function all_formats(K, P)
    (
        uf = Format{is_unsigned, is_finite}(K, P),
        ue = Format{is_unsigned, is_extended}(K, P),
        sf = Format{is_signed,   is_finite}(K, P),
        se = Format{is_signed,   is_extended}(K, P),
    )
end

# =========================================================================
# §1  Universal invariants — must hold for every valid (K, P) combination
# =========================================================================

@testset "Universal invariants across format variants" begin
    # Sweep over a range of (K, P) including edge cases
    for K in 3:12
        for P in 1:(K - 1)   # need W ≥ 1 for signed: K - P + 1 - 1 ≥ 1 → P ≤ K - 1
            K - P + 1 >= 1 || continue  # unsigned constraint

            fmts = all_formats(K, P)

            for fmt in fmts
                @testset "K=$K P=$P $(typeof(fmt))" begin
                    # --- total partition ---
                    @test nValuesOf(fmt) == twopow(BitwidthOf(fmt))
                    @test nValuesOf(fmt) == nNumericalValuesOf(fmt) + nNaNsOf(fmt)
                    @test nNumericalValuesOf(fmt) == nZerosOf(fmt) + nNonZeroNumericalValuesOf(fmt)
                    @test nNonZeroNumericalValuesOf(fmt) == nNonZeroFiniteValuesOf(fmt) + nInfsOf(fmt)

                    # --- finite partition ---
                    @test nFiniteValuesOf(fmt) == nNonZeroFiniteValuesOf(fmt) + nZerosOf(fmt)
                    @test nNonZeroFiniteValuesOf(fmt) == nPosFiniteValuesOf(fmt) + nNegFiniteValuesOf(fmt)
                    @test nNonNegFiniteValuesOf(fmt) == nPosFiniteValuesOf(fmt) + nZerosOf(fmt)

                    # --- sign split of non-zero numerical ---
                    @test nNonZeroNumericalValuesOf(fmt) == nPosValuesOf(fmt) + nNegValuesOf(fmt)

                    # --- infinity partition ---
                    @test nInfsOf(fmt) == nPosInfsOf(fmt) + nNegInfsOf(fmt)

                    # --- finite = prenormals + normals ---
                    @test nFiniteValuesOf(fmt) == nPrenormalsOf(fmt) + nNormalsOf(fmt)

                    # --- prenormals = zero + subnormals ---
                    @test nPrenormalsOf(fmt) == nZerosOf(fmt) + nSubnormalsOf(fmt)

                    # --- subnormal split ---
                    @test nSubnormalsOf(fmt) == nPosSubnormalsOf(fmt) + nNegSubnormalsOf(fmt)

                    # --- normal split ---
                    @test nNormalsOf(fmt) == nPosNormalsOf(fmt) + nNegNormalsOf(fmt)

                    # --- pos finite = pos subnormals + pos normals ---
                    @test nPosFiniteValuesOf(fmt) == nPosSubnormalsOf(fmt) + nPosNormalsOf(fmt)

                    # --- neg finite = neg subnormals + neg normals ---
                    @test nNegFiniteValuesOf(fmt) == nNegSubnormalsOf(fmt) + nNegNormalsOf(fmt)

                    # --- non-negativity of all counts ---
                    @test nNaNsOf(fmt) >= 0
                    @test nZerosOf(fmt) >= 0
                    @test nInfsOf(fmt) >= 0
                    @test nPosInfsOf(fmt) >= 0
                    @test nNegInfsOf(fmt) >= 0
                    @test nValuesOf(fmt) >= 0
                    @test nNumericalValuesOf(fmt) >= 0
                    @test nFiniteValuesOf(fmt) >= 0
                    @test nPosValuesOf(fmt) >= 0
                    @test nNegValuesOf(fmt) >= 0
                    @test nPrenormalsOf(fmt) >= 0
                    @test nSubnormalsOf(fmt) >= 0
                    @test nPosSubnormalsOf(fmt) >= 0
                    @test nNegSubnormalsOf(fmt) >= 0
                    @test nNormalsOf(fmt) >= 0
                    @test nPosNormalsOf(fmt) >= 0
                    @test nNegNormalsOf(fmt) >= 0
                    @test nBinadesOf(fmt) >= 1
                end
            end
        end
    end
end

# =========================================================================
# §2  Unsigned-specific invariants
# =========================================================================

@testset "Unsigned format invariants" begin
    for K in 3:10, P in 1:K
        K - P + 1 >= 1 || continue

        for fmt in (Format{is_unsigned, is_finite}(K, P),
                    Format{is_unsigned, is_extended}(K, P))
            @testset "K=$K P=$P $(typeof(fmt))" begin
                @test nNegValuesOf(fmt) == 0
                @test nNegFiniteValuesOf(fmt) == 0
                @test nNegSubnormalsOf(fmt) == 0
                @test nNegNormalsOf(fmt) == 0
                @test nNegInfsOf(fmt) == 0
                @test nPosValuesOf(fmt) == nNonZeroNumericalValuesOf(fmt)
            end
        end
    end
end

# =========================================================================
# §3  Signed-specific invariants
# =========================================================================

@testset "Signed format invariants" begin
    for K in 3:10, P in 1:(K - 1)
        K - P >= 1 || continue

        for fmt in (Format{is_signed, is_finite}(K, P),
                    Format{is_signed, is_extended}(K, P))
            @testset "K=$K P=$P $(typeof(fmt))" begin
                # signed formats have symmetric positive/negative counts
                @test nPosValuesOf(fmt) == nNegValuesOf(fmt)
                @test nPosSubnormalsOf(fmt) == nNegSubnormalsOf(fmt)
                @test nPosNormalsOf(fmt) == nNegNormalsOf(fmt)
                @test nPosFiniteValuesOf(fmt) == nNegFiniteValuesOf(fmt)
            end
        end
    end
end

# =========================================================================
# §4  Finite vs Extended domain invariants
# =========================================================================

@testset "Finite-domain formats have no infinities" begin
    for K in 3:10, P in 1:(K - 1)
        for fmt in (Format{is_unsigned, is_finite}(K, P),
                    Format{is_signed, is_finite}(K, P))
            @test nInfsOf(fmt) == 0
            @test nPosInfsOf(fmt) == 0
            @test nNegInfsOf(fmt) == 0
        end
    end
end

@testset "Extended-domain infinity counts" begin
    for K in 3:10, P in 1:K
        K - P + 1 >= 1 || continue
        ue = Format{is_unsigned, is_extended}(K, P)
        @test nInfsOf(ue) == 1
        @test nPosInfsOf(ue) == 1
        @test nNegInfsOf(ue) == 0
    end

    for K in 3:10, P in 1:(K - 1)
        se = Format{is_signed, is_extended}(K, P)
        @test nInfsOf(se) == 2
        @test nPosInfsOf(se) == 1   # only +Inf is positive
        @test nNegInfsOf(se) == 1   # only -Inf is negative
    end
end

# =========================================================================
# §5  P = 1 corner case: no subnormals
# =========================================================================

@testset "P=1 → no subnormals" begin
    for K in 2:10
        K >= 2 || continue
        for (S, D) in ((is_unsigned, is_finite), (is_unsigned, is_extended),
                        (is_signed, is_finite),  (is_signed, is_extended))
            # check W constraint
            Σ = S === is_signed ? 1 : 0
            W = K - 1 + 1 - Σ
            W >= 1 || continue

            fmt = Format{S, D}(K, 1)
            @test nPosSubnormalsOf(fmt) == 0
            @test nNegSubnormalsOf(fmt) == 0
            @test nSubnormalsOf(fmt) == 0
            @test nNonNegPrenormalsOf(fmt) == 1   # just zero
            @test nPrenormalsOf(fmt) == 1
        end
    end
end

# =========================================================================
# §6  Hand-computed values for specific small formats
# =========================================================================

@testset "Hand-computed: Format{unsigned, finite}(3, 2)" begin
    # K=3, P=2: 8 code points, W=2, bias=2
    # cp 0=zero, 1=subnormal, 2..6=normals, 7=NaN
    fmt = Format{is_unsigned, is_finite}(3, 2)

    @test nValuesOf(fmt) == 8
    @test nNaNsOf(fmt) == 1
    @test nNumericalValuesOf(fmt) == 7
    @test nZerosOf(fmt) == 1
    @test nNonZeroNumericalValuesOf(fmt) == 6
    @test nInfsOf(fmt) == 0
    @test nFiniteValuesOf(fmt) == 7
    @test nNonZeroFiniteValuesOf(fmt) == 6
    @test nPosValuesOf(fmt) == 6
    @test nNegValuesOf(fmt) == 0
    @test nPosFiniteValuesOf(fmt) == 6
    @test nNegFiniteValuesOf(fmt) == 0

    @test nNonNegPrenormalsOf(fmt) == 2   # twopow(1) = 2
    @test nPrenormalsOf(fmt) == 2          # zero + 1 subnormal
    @test nPosSubnormalsOf(fmt) == 1
    @test nNegSubnormalsOf(fmt) == 0
    @test nSubnormalsOf(fmt) == 1

    @test nNormalsOf(fmt) == 5
    @test nPosNormalsOf(fmt) == 5
    @test nNegNormalsOf(fmt) == 0

    @test nBinadesOf(fmt) == 3  # 2*2-1 = 3
end

@testset "Hand-computed: Format{unsigned, extended}(3, 2)" begin
    # Same as above but cp 6=+Inf, so only 5 positive finite values
    fmt = Format{is_unsigned, is_extended}(3, 2)

    @test nValuesOf(fmt) == 8
    @test nNumericalValuesOf(fmt) == 7
    @test nInfsOf(fmt) == 1
    @test nPosInfsOf(fmt) == 1
    @test nNegInfsOf(fmt) == 0
    @test nFiniteValuesOf(fmt) == 6
    @test nNonZeroFiniteValuesOf(fmt) == 5
    @test nPosFiniteValuesOf(fmt) == 5

    @test nPosSubnormalsOf(fmt) == 1
    @test nPosNormalsOf(fmt) == 4
end

@testset "Hand-computed: Format{signed, finite}(4, 2)" begin
    # K=4, P=2, W=4-2+1-1=2, bias=2
    # 16 code points: positive half 0..7, negative half 9..15
    # cp 0=zero, 1=pos subnormal, 2..7=pos normals (6 normals)
    # cp 8=NaN
    # cp 9=neg subnormal, 10..15=neg normals (6 normals)
    fmt = Format{is_signed, is_finite}(4, 2)

    @test nValuesOf(fmt) == 16
    @test nNaNsOf(fmt) == 1
    @test nNumericalValuesOf(fmt) == 15
    @test nNonZeroNumericalValuesOf(fmt) == 14
    @test nInfsOf(fmt) == 0
    @test nFiniteValuesOf(fmt) == 15
    @test nNonZeroFiniteValuesOf(fmt) == 14

    @test nPosValuesOf(fmt) == 7
    @test nNegValuesOf(fmt) == 7
    @test nPosFiniteValuesOf(fmt) == 7
    @test nNegFiniteValuesOf(fmt) == 7

    @test nPosSubnormalsOf(fmt) == 1
    @test nNegSubnormalsOf(fmt) == 1
    @test nSubnormalsOf(fmt) == 2
    @test nPosNormalsOf(fmt) == 6
    @test nNegNormalsOf(fmt) == 6
    @test nNormalsOf(fmt) == 12
end

@testset "Hand-computed: Format{signed, extended}(4, 2)" begin
    # K=4, P=2, W=2, bias=2
    # cp 0=zero, 1=pos sub, 2..6=pos normals, 7=+Inf
    # cp 8=NaN
    # cp 9=neg sub, 10..14=neg normals, 15=-Inf
    fmt = Format{is_signed, is_extended}(4, 2)

    @test nValuesOf(fmt) == 16
    @test nNaNsOf(fmt) == 1
    @test nNumericalValuesOf(fmt) == 15
    @test nNonZeroNumericalValuesOf(fmt) == 14
    @test nInfsOf(fmt) == 2
    @test nPosInfsOf(fmt) == 1
    @test nNegInfsOf(fmt) == 1
    @test nFiniteValuesOf(fmt) == 13
    @test nNonZeroFiniteValuesOf(fmt) == 12

    @test nPosValuesOf(fmt) == 7
    @test nNegValuesOf(fmt) == 7
    @test nPosFiniteValuesOf(fmt) == 6
    @test nNegFiniteValuesOf(fmt) == 6

    @test nPosSubnormalsOf(fmt) == 1
    @test nNegSubnormalsOf(fmt) == 1
    @test nSubnormalsOf(fmt) == 2
    @test nPosNormalsOf(fmt) == 5
    @test nNegNormalsOf(fmt) == 5
    @test nNormalsOf(fmt) == 10
end

# =========================================================================
# §7  Minimal formats — smallest valid parameter combinations
# =========================================================================

@testset "Minimal format: unsigned finite K=2 P=1" begin
    # K=2, P=1: W=2, bias=2, 4 code points
    # cp 0=zero, 1,2=normals, 3=NaN
    fmt = Format{is_unsigned, is_finite}(2, 1)

    @test nValuesOf(fmt) == 4
    @test nNumericalValuesOf(fmt) == 3
    @test nFiniteValuesOf(fmt) == 3
    @test nPosSubnormalsOf(fmt) == 0     # P=1 → no subnormals
    @test nPrenormalsOf(fmt) == 1         # just zero
    @test nNormalsOf(fmt) == 2
    @test nBinadesOf(fmt) == 3            # 2*2-1 = 3
end

@testset "Minimal format: unsigned finite K=2 P=2" begin
    # K=2, P=2: W=1, ExponentBias=1 (unsigned: 2^NonSig = 2^0 = 1)
    # 4 code points, 1 trailing bit
    # cp 0=zero, 1=subnormal, 2=normal, 3=NaN
    fmt = Format{is_unsigned, is_finite}(2, 2)

    @test nValuesOf(fmt) == 4
    @test nNumericalValuesOf(fmt) == 3
    @test nPosSubnormalsOf(fmt) == 1
    @test nPosNormalsOf(fmt) == 1
    @test nBinadesOf(fmt) == 1   # 2*1-1 = 1
end

@testset "Minimal format: signed finite K=3 P=1" begin
    # K=3, P=1, W=3-1+1-1=2, bias=2
    # 8 code points: pos half 0..3, neg half 5..7
    # cp 0=zero, 1,2,3=pos normals, 4=NaN, 5,6,7=neg normals
    fmt = Format{is_signed, is_finite}(3, 1)

    @test nValuesOf(fmt) == 8
    @test nNumericalValuesOf(fmt) == 7
    @test nPosSubnormalsOf(fmt) == 0
    @test nNegSubnormalsOf(fmt) == 0
    @test nPosNormalsOf(fmt) == 3
    @test nNegNormalsOf(fmt) == 3
end

# =========================================================================
# §8  Larger formats — structural sanity checks
# =========================================================================

@testset "8-bit formats (like IEEE-like tiny floats)" begin
    # K=8, P=4: similar to E4M3 style
    for (S, D) in ((is_unsigned, is_finite), (is_unsigned, is_extended),
                    (is_signed, is_finite),   (is_signed, is_extended))
        Σ = S === is_signed ? 1 : 0
        W = 8 - 4 + 1 - Σ
        W >= 1 || continue
        fmt = Format{S, D}(8, 4)

        @test nValuesOf(fmt) == 256
        @test nNaNsOf(fmt) == 1
        @test nZerosOf(fmt) == 1
    end
end

@testset "16-bit format counts" begin
    # K=16, P=8: moderate-size format
    fmt_uf = Format{is_unsigned, is_finite}(16, 8)
    fmt_se = Format{is_signed, is_extended}(16, 8)

    @test nValuesOf(fmt_uf) == 65536
    @test nValuesOf(fmt_se) == 65536

    # unsigned has more positive values than signed
    @test nPosValuesOf(fmt_uf) > nPosValuesOf(fmt_se)

    # signed extended has exactly 2 infinities
    @test nInfsOf(fmt_se) == 2
end

# =========================================================================
# §9  Cross-validation with exhaustive enumeration (small formats only)
# =========================================================================

@testset "Cross-validate counts against exhaustive enumeration K ≤ 8" begin
    for K in 2:8
        for P in 1:K
            for (S, D) in ((is_unsigned, is_finite), (is_unsigned, is_extended),
                            (is_signed, is_finite),  (is_signed, is_extended))
                Σ = S === is_signed ? 1 : 0
                W = K - P + 1 - Σ
                W >= 1 || continue

                fmt = Format{S, D}(K, P)

                @testset "Enum K=$K P=$P $(S) $(D)" begin
                    # count by walking all code points
                    n_finite_enum = 0
                    n_pos_finite_enum = 0
                    n_neg_finite_enum = 0
                    n_zero_enum = 0
                    n_pos_sub_enum = 0
                    n_neg_sub_enum = 0

                    nan_cp_val = cp_nan(fmt)
                    inf_cp_val = cp_inf(fmt)
                    ninf_cp_val = cp_neginf(fmt)

                    for cp in BigInt(0):BigInt(cp_max(fmt))
                        # skip specials
                        cp == nan_cp_val && continue
                        inf_cp_val !== nothing && cp == inf_cp_val && continue
                        ninf_cp_val !== nothing && cp == ninf_cp_val && continue

                        val = ValueOf(fmt, cp)
                        n_finite_enum += 1
                        if val == 0
                            n_zero_enum += 1
                        elseif val > 0
                            n_pos_finite_enum += 1
                            # subnormal if cp_abs < significand_scale
                            cp_abs = cp < BigInt(sign_half_offset(fmt)) ? cp : cp - BigInt(sign_half_offset(fmt))
                            if cp_abs < significand_scale(fmt)
                                n_pos_sub_enum += 1
                            end
                        else
                            n_neg_finite_enum += 1
                            # negative subnormal
                            cp_abs = cp - BigInt(sign_half_offset(fmt))
                            if cp_abs < significand_scale(fmt)
                                n_neg_sub_enum += 1
                            end
                        end
                    end

                    @test nFiniteValuesOf(fmt) == n_finite_enum
                    @test nZerosOf(fmt) == n_zero_enum
                    @test nPosFiniteValuesOf(fmt) == n_pos_finite_enum
                    @test nNegFiniteValuesOf(fmt) == n_neg_finite_enum
                    @test nPosSubnormalsOf(fmt) == n_pos_sub_enum
                    @test nNegSubnormalsOf(fmt) == n_neg_sub_enum
                end
            end
        end
    end
end

# =========================================================================
# §10  nBinadesOf consistency
# =========================================================================

@testset "nBinadesOf equals number of normal binades" begin
    for K in 3:10, P in 1:(K - 1)
        for (S, D) in ((is_unsigned, is_finite), (is_signed, is_extended))
            Σ = S === is_signed ? 1 : 0
            W = K - P + 1 - Σ
            W >= 1 || continue

            fmt = Format{S, D}(K, P)
            B = ExponentBiasOf(fmt)

            @test nBinadesOf(fmt) == 2 * B - 1

            # each binade has twopow(P-1) normal code points
            if P >= 2
                @test nPosNormalsOf(fmt) == nBinadesOf(fmt) * twopow(P - 1)
            else
                # P=1: each binade has 1 normal
                @test nPosNormalsOf(fmt) == nBinadesOf(fmt)
            end
        end
    end
end

# =========================================================================
# §11  High-precision formats — large P relative to K
# =========================================================================

@testset "High precision: P close to K" begin
    # K=8, P=7 unsigned → W = 8-7+1 = 2
    fmt = Format{is_unsigned, is_finite}(8, 7)
    @test nPosSubnormalsOf(fmt) == 63   # twopow(6) - 1
    @test nBinadesOf(fmt) == 3           # 2*2-1

    # K=8, P=8 unsigned → W = 8-8+1 = 1
    fmt = Format{is_unsigned, is_finite}(8, 8)
    @test nPosSubnormalsOf(fmt) == 127  # twopow(7) - 1
    @test nBinadesOf(fmt) == 1           # 2*1-1
    @test nPosNormalsOf(fmt) == 127      # 1 binade × 2^7 code points... wait
    # 1 binade with twopow(7) = 128 entries but one is taken by the binade start
    # Actually nPosNormals = nFinite - nPrenormals for unsigned
end

# =========================================================================
# §12  Wide exponent formats — large W, small P
# =========================================================================

@testset "Wide exponent: P=1 various K" begin
    for K in 2:10
        fmt = Format{is_unsigned, is_finite}(K, 1)
        @test nPosSubnormalsOf(fmt) == 0
        @test nSubnormalsOf(fmt) == 0
        @test nPrenormalsOf(fmt) == 1   # just zero
        @test nNormalsOf(fmt) == nFiniteValuesOf(fmt) - 1
        # All finite non-zero values are normal
        @test nPosNormalsOf(fmt) == nNonZeroFiniteValuesOf(fmt)
    end
end

@testset "Wide exponent: P=2 various K" begin
    for K in 2:10
        fmt = Format{is_unsigned, is_finite}(K, 2)
        @test nPosSubnormalsOf(fmt) == 1   # twopow(1) - 1
        @test nNonNegPrenormalsOf(fmt) == 2
    end
end

# =========================================================================
# §13  Comparison: finite vs extended for same (K, P)
# =========================================================================

@testset "Finite has more finite values than extended" begin
    for K in 3:10, P in 1:(K - 1)
        uf = Format{is_unsigned, is_finite}(K, P)
        ue = Format{is_unsigned, is_extended}(K, P)
        sf = Format{is_signed, is_finite}(K, P)
        se = Format{is_signed, is_extended}(K, P)

        @test nFiniteValuesOf(uf) > nFiniteValuesOf(ue)
        @test nFiniteValuesOf(sf) > nFiniteValuesOf(se)

        # difference equals number of infinities
        @test nFiniteValuesOf(uf) - nFiniteValuesOf(ue) == nInfsOf(ue)
        @test nFiniteValuesOf(sf) - nFiniteValuesOf(se) == nInfsOf(se)
    end
end

# =========================================================================
# §14  Comparison: unsigned vs signed for same (K, P, D)
# =========================================================================

@testset "Unsigned has more positive values than signed" begin
    for K in 3:10, P in 1:(K - 1)
        for D in (is_finite, is_extended)
            uf = Format{is_unsigned, D}(K, P)
            sf = Format{is_signed, D}(K, P)

            @test nPosFiniteValuesOf(uf) > nPosFiniteValuesOf(sf)
        end
    end
end

# =========================================================================
# §15  AllFiniteValuesOf / AllPositiveFiniteValuesOf cross-check
# =========================================================================

@testset "Count functions match length of enumerated value lists K ≤ 7" begin
    for K in 2:7
        for P in 1:K
            for (S, D) in ((is_unsigned, is_finite), (is_unsigned, is_extended),
                            (is_signed, is_finite),  (is_signed, is_extended))
                Σ = S === is_signed ? 1 : 0
                W = K - P + 1 - Σ
                W >= 1 || continue

                fmt = Format{S, D}(K, P)
                @testset "Enum-len K=$K P=$P $(S) $(D)" begin
                    all_fin = AllFiniteValuesOf(fmt)
                    @test length(all_fin) == nFiniteValuesOf(fmt)

                    all_pos = AllPositiveFiniteValuesOf(fmt)
                    @test length(all_pos) == nPosFiniteValuesOf(fmt)
                end
            end
        end
    end
end

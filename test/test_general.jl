
@testset "P3109Format constructor" begin
    f = P3109Format(8, 3, 0, 1)
    @test f.K == 8
    @test f.P == 3
    @test f.Σ == 0
    @test f.Δ == 1

    @test_throws ArgumentError P3109Format(1, 4, 0, 0)
    @test_throws ArgumentError P3109Format(8, 3, 2, 0)
    @test_throws ArgumentError P3109Format(8, 3, 0, 2)
end

@testset "Basic structural quantities" begin
    f = P3109Format(8, 3, 0, 1)
    @test exponent_width(f) == 6
    @test significand_scale(f) == UInt128(4)
    @test bias(f) == UInt128(32)
    @test sign_half_offset(f) == UInt128(0)
end

@testset "Special code points" begin
    fu = P3109Format(8, 3, 0, 1)
    fs = P3109Format(8, 3, 1, 1)

    @test nan_cp(fu) == UInt128(255)
    @test posinf_cp(fu) == UInt128(254)
    @test neginf_cp(fu) === nothing

    @test nan_cp(fs) == UInt128(128)
    @test posinf_cp(fs) == UInt128(127)
    @test neginf_cp(fs) == UInt128(255)
end

@testset "Decode monotonicity and roundtrip" begin
    f = P3109Format(6, 3, 0, 0)
    @test verify_monotone_positive_half(f)
    @test exhaustive_roundtrip_check(f)
    @test exhaustive_monotonicity_check(f)
end

@testset "Known unsigned Binary3p2 finite-only ladder" begin
    f = P3109Format(3, 2, 0, 0)
    vals = Rational{BigInt}[]
    for cp in 0:Int(cp_nor_max(f))
        dv = decode(f, cp)
        dv isa FiniteValue || continue
        push!(vals, dv.value)
    end
    @test vals == Rational{BigInt}[0//1, 1//4, 1//2, 3//4, 1//1, 3//2, 2//1, 3//1]
end

@testset "Statistics are nonnegative" begin
    f = P3109Format(8, 3, 0, 0)
    μ = exact_mean_step_normal_binades(f)
    m2 = exact_second_moment_step_normal_binades(f)
    σ2 = exact_variance_step_normal_binades(f)

    @test μ > 0
    @test m2 > 0
    @test σ2 >= 0
end

@testset "Exhaustive small-format sweep K ≤ 8" begin
    for K in 2:8
        for P in 1:K
            for Σ in 0:1, Δ in 0:1
                W = K - P + 1 - Σ
                W >= 1 || continue
                f = P3109Format(K, P, Σ, Δ)
                @test exhaustive_roundtrip_check(f)
                @test exhaustive_monotonicity_check(f)
            end
        end
    end
end

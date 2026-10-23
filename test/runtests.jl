using Float3109s
using Test

@testset "Float3109s.jl" begin
    @testset "Binary parameter bounds" begin
        F = Float3109s
        for D in (F.FiniteFormat, F.ExtendedFormat)
            for S in (F.UnsignedFormat, F.SignedFormat)
                x = F.Binary(3, 1, S, D)
                @test F.Binary(3, 1, S(), D()) === x
                @test F.BitwidthOf(x) === F.IntFormat(3)
                @test F.PrecisionOf(x) === F.IntFormat(1)
                @test F.SignednessOf(x) === S()
                @test F.DomainOf(x) === D()
                @test F.Binary(3, 2, S, D) isa F.Binary
                for (k, p) in ((2, 1), (0, 1), (-1, 1), (3, 0), (3, -1), (3, 4))
                    @test_throws ArgumentError F.Binary(k, p, S, D)
                    @test_throws ArgumentError F.Binary{F.IntFormat(k), F.IntFormat(p), S, D}()
                end
                @test_throws ArgumentError F.Binary{3.0, F.IntFormat(1), S, D}()
                @test_throws ArgumentError F.Binary{F.IntFormat(3), 1.0, S, D}()
            end
            @test F.Binary(3, 3, F.UnsignedFormat, D) isa F.Binary
            @test_throws ArgumentError F.Binary(3, 3, F.SignedFormat, D)
        end
    end
end

@testset "Datum constructors and trait dispatch" begin
    F = Float3109s
    for S in (F.UnsignedFormat, F.SignedFormat), D in (F.FiniteFormat, F.ExtendedFormat)
        format_type = F.BinaryFormat(5, 2, S, D)
        format = F.Binary(format_type)
        B = typeof(format)
        for code in (UInt8(0x11), UInt16(0x11))
            T = F.Datum{B, typeof(code)}
            x = T(code)
            @test F.Datum(format_type, code) === x
            @test F.Datum(B, code) === x
            @test F.Datum(format, code) === x
            @test (@inferred F.codepoint(x)) === code
            for arg in (F.AbstractDatum{B}, T, x)
                @test (@inferred F.BitwidthOf(arg)) === F.IntFormat(5)
                @test (@inferred F.PrecisionOf(arg)) === F.IntFormat(2)
                @test (@inferred F.is_unsigned(arg)) === (S === F.UnsignedFormat)
                @test (@inferred F.is_signed(arg)) === (S === F.SignedFormat)
                @test (@inferred F.is_finite(arg)) === (D === F.FiniteFormat)
                @test (@inferred F.is_extended(arg)) === (D === F.ExtendedFormat)
                @test (@inferred F.FormatOf(arg)) === F.FormatOf(format)
            end
            for type in (F.AbstractDatum{B}, T)
                @test (@inferred F.SignednessOf(type)) === S
                @test (@inferred F.DomainOf(type)) === D
            end
            @test (@inferred F.SignednessOf(x)) === S()
            @test (@inferred F.DomainOf(x)) === D()
            @test x == F.Datum(format, code)
            @test isequal(x, F.Datum(format, code))
            @test hash(x) == hash(F.Datum(format, code))
            @test !(x == F.Datum(format, code + one(code)))
        end
    end
    @test_throws ArgumentError F.Datum(F.BinaryFormat(2, 1, F.UnsignedFormat, F.FiniteFormat), 0x01)
end

@testset "Binary inference and trait validation" begin
    F = Float3109s
    for S in (F.UnsignedFormat, F.SignedFormat), D in (F.FiniteFormat, F.ExtendedFormat)
        T = F.Binary{Int16(5), Int16(2), S, D}
        x = @inferred T()
        @test (@inferred F.SignednessOf(T)) === S
        @test (@inferred F.DomainOf(T)) === D
        @test (@inferred F.SignednessOf(x)) === S()
        @test (@inferred F.DomainOf(x)) === D()
        for arg in (T, x)
            @test (@inferred F.BitwidthOf(arg)) === Int16(5)
            @test (@inferred F.PrecisionOf(arg)) === Int16(2)
            @test (@inferred F.is_unsigned(arg)) === (S === F.UnsignedFormat)
            @test (@inferred F.is_signed(arg)) === (S === F.SignedFormat)
            @test (@inferred F.is_finite(arg)) === (D === F.FiniteFormat)
            @test (@inferred F.is_extended(arg)) === (D === F.ExtendedFormat)
        end
        @test F.Binary(Int32(5), Int8(2), S(), D()) === x
        @test_throws InexactError F.Binary(32768, 2, S, D)
    end
    @test_throws ArgumentError F.Binary(5, 2, F.Signedness, F.FiniteFormat)
    @test_throws ArgumentError F.Binary(5, 2, F.SignedFormat, F.Domain)
end

@testset "Signedness accessor display" begin
    F = Float3109s
    for (S, type_str, instance_str) in (
        (F.UnsignedFormat, F.UnsignedStr, F.UnsignedFormatStr),
        (F.SignedFormat, F.SignedStr, F.SignedFormatStr),
    ), D in (F.FiniteFormat, F.ExtendedFormat)
        T = F.Binary{Int16(5), Int16(2), S, D}
        @test F.SignednessOf(T) === S
        @test F.SignednessOf(T()) === S()
        for (tag, expected) in ((F.SignednessOf(T), type_str), (F.SignednessOf(T()), instance_str))
            @test string(tag) === expected
            @test sprint(print, tag) == expected
            @test repr(MIME"text/plain"(), tag) == expected
        end
    end
end

@testset "Domain accessor display" begin
    F = Float3109s
    for (D, type_str, instance_str) in (
        (F.FiniteFormat, F.FiniteStr, F.FiniteFormatStr),
        (F.ExtendedFormat, F.ExtendedStr, F.ExtendedFormatStr),
    ), S in (F.UnsignedFormat, F.SignedFormat)
        T = F.Binary{Int16(5), Int16(2), S, D}
        @test F.DomainOf(T) === D
        @test F.DomainOf(T()) === D()
        for (tag, expected) in ((F.DomainOf(T), type_str), (F.DomainOf(T()), instance_str))
            @test string(tag) === expected
            @test sprint(print, tag) == expected
            @test repr(tag) == expected
            @test repr(MIME"text/plain"(), tag) == expected
        end
    end
end

abstract type BinaryFormat{K, P, Σ<:Signedness, Δ<:Domain} end

BitwidthOf(::Type{<:BinaryFormat{K}}) where {K} = K
PrecisionOf(::Type{<:BinaryFormat{K,P}}) where {K,P} = P
SignednessOf(::Type{<:BinaryFormat{K,P,S}}) where {K,P,S} = S
DomainOf(::Type{<:BinaryFormat{K,P,S,D}}) where {K,P,S,D} = D

is_unsigned(::Type{<:BinaryFormat{K,P,S,D}}) where {K,P,S,D} =
    is_unsigned(S)
is_signed(::Type{<:BinaryFormat{K,P,S,D}}) where {K,P,S,D} =
    is_signed(S)
is_finite(::Type{<:BinaryFormat{K,P,S,D}}) where {K,P,S,D} =
    is_finite(D)
is_extended(::Type{<:BinaryFormat{K,P,S,D}}) where {K,P,S,D} =
    is_extended(D)

BitwidthOf(::BinaryFormat{K}) where {K} = K
PrecisionOf(::BinaryFormat{K,P}) where {K,P} = P
SignednessOf(::BinaryFormat{K,P,S}) where {K,P,S} = S()
DomainOf(::BinaryFormat{K,P,S,D}) where {K,P,S,D} = D()

is_unsigned(::BinaryFormat{K,P,S,D}) where {K,P,S,D} =
    is_unsigned(S)
is_signed(::BinaryFormat{K,P,S,D}) where {K,P,S,D} =
    is_signed(S)
is_finite(::BinaryFormat{K,P,S,D}) where {K,P,S,D} =
    is_finite(D)
is_extended(::BinaryFormat{K,P,S,D}) where {K,P,S,D} =
    is_extended(D)

struct Binary{K, P, Σ, Δ} <: BinaryFormat{K, P, Σ, Δ}
    function Binary{K, P, Σ, Δ}() where {
        K, P, Σ<:Signedness, Δ<:Domain
    }
        K isa IntFormat && P isa IntFormat ||
            throw(ArgumentError("K and P must be values of type $IntFormat"))
        K > 2 ||
            throw(ArgumentError("K ($K) must be greater than 2"))
        P > 0 ||
            throw(ArgumentError("P ($P) must be greater than 0"))

        isconcretetype(Σ) ||
            throw(ArgumentError("Σ must be a concrete signedness tag type"))
        isconcretetype(Δ) ||
            throw(ArgumentError("Δ must be a concrete domain tag type"))

        if is_unsigned(Σ)
            P <= K ||
                throw(ArgumentError("P ($P) must not exceed K ($K) for unsigned formats"))
        else
            P < K ||
                throw(ArgumentError("P ($P) must be less than K ($K) for signed formats"))
        end

        new{K, P, Σ, Δ}()
    end
end

Binary(k::Integer, p::Integer, σ::Signedness, δ::Domain) =
    Binary(k, p, typeof(σ), typeof(δ))

function Binary(k::Integer, p::Integer, ::Type{S}, ::Type{D}) where {
    S<:Signedness, D<:Domain
}
    K = IntFormat(k)
    P = IntFormat(p)
    Binary{K, P, S, D}()
end

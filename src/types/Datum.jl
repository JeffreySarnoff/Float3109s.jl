export AbstractDatum, Datum,
    codepoint,
    BitwidthOf, PrecisionOf, SignednessOf, DomainOf,
    is_unsigned, is_signed, is_finite, is_extended

public FormatOf


abstract type AbstractDatum{F<:Binary} <: AbstractFloat end

BitwidthOf(::Type{<:AbstractDatum{F}}) where {F} = BitwidthOf(F)
PrecisionOf(::Type{<:AbstractDatum{F}}) where {F} = PrecisionOf(F)
SignednessOf(::Type{<:AbstractDatum{F}}) where {F} = SignednessOf(F)
DomainOf(::Type{<:AbstractDatum{F}}) where {F} = DomainOf(F)

is_unsigned(::Type{<:AbstractDatum{F}}) where {F} = is_unsigned(F)
is_signed(::Type{<:AbstractDatum{F}}) where {F} = is_signed(F)
is_finite(::Type{<:AbstractDatum{F}}) where {F} = is_finite(F)
is_extended(::Type{<:AbstractDatum{F}}) where {F} = is_extended(F)

BitwidthOf(::AbstractDatum{F}) where {F} = BitwidthOf(F)
PrecisionOf(::AbstractDatum{F}) where {F} = PrecisionOf(F)
SignednessOf(::AbstractDatum{F}) where {F} = SignednessOf(F)()
DomainOf(::AbstractDatum{F}) where {F} = DomainOf(F)()

is_unsigned(::AbstractDatum{F}) where {F} = is_unsigned(F)
is_signed(::AbstractDatum{F}) where {F} = is_signed(F)
is_finite(::AbstractDatum{F}) where {F} = is_finite(F)
is_extended(::AbstractDatum{F}) where {F} = is_extended(F)

struct Datum{ F<:Binary, U<:UIntDatum } <: AbstractDatum{F}
    code::U
end

Datum(f::Binary, code::UIntDatum) = Datum{typeof(f), typeof(code)}(code)

function Datum(T::Type{<:BinaryFormat}, code::UIntDatum)
    Datum(Binary(T), code)
end

codepoint(@nospecialize x::Datum) = x.code

@inline function Base.:(==)(x::Datum{F}, y::Datum{G}) where {F,G}
    return F === G && x.code == y.code
end

@inline function Base.isequal(x::Datum{F}, y::Datum{G}) where {F,G}
    return F === G && isequal(x.code, y.code)
end

@inline function Base.hash(x::Datum{F}, h::UInt) where {F}
    return hash(x.code, hash(F, h))
end

FormatOf(::Type{<:AbstractDatum{F}}) where {F} = supertype(F)
FormatOf(::AbstractDatum{F}) where {F} = supertype(F)

#=
julia> sf52x11 = Datum{sf52type, UInt8}( 0x11 )
Datum{Binary{5, 2, 𝙎, 𝙁}, UInt8}(0x11)
=#

export AbstractDatum, Datum
    BitwidthOf, PrecisionOf, SignednessOf, DomainOf,
    is_unsigned, is_signed, is_finite, is_extended

public FormatOf


abstract type AbstractDatum{F<:Binary} <: AbstractFloat end

struct Datum{ F<:Binary, U<:UIntDatum } <: AbstractDatum{F}
    code::U
end


BitwidthOf(::Datum{F, U}) where {F, U} = BitwidthOf(F)
PrecisionOf(::Datum{F, U}) where {F, U} = PrecisionOf(F)
SignednessOf(::Datum{F, U}) where {F, U} = SignednessOf(F)
DomainOf(::Datum{F, U}) where {F, U} = DomainOf(F)

is_unsigned(::Datum{F, U}) where {F, U} = is_unsigned(F)
is_signed(::Datum{F, U}) where {F, U} = is_signed(F)
is_finite(::Datum{F, U}) where {F, U} = is_finite(F)
is_extended(::Datum{F, U}) where {F, U} = is_extended(F)

FormatOf(x::Datum{F,U}) where {F<:Binary, U<:UIntDatum} =
    supertype(F)

#=
julia> sf52x11 = Datum{sf52type, UInt8}( 0x11 )
Datum{Binary{5, 2, 𝙎, 𝙁}, UInt8}(0x11)

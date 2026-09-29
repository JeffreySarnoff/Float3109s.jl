abstract type AbstractBinary{F<:Binary} <: AbstractFloat end

struct Datum{ F<:Binary, U<:UIntDatum } <: AbstractBinary{F}
    code::U
end

#=
julia> sf52x11 = Datum{sf52type, UInt8}( 0x11 )
Datum{Binary{5, 2, 𝙎, 𝙁}, UInt8}(0x11)
=#
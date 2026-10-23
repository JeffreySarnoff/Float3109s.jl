module Float3109s

include("constants.jl")

include("types/traits/signedness.jl")
include("types/traits/domain.jl")

include("types/Binary.jl")
include("types/Datum.jl")

# tiny trial

sf52type = Binary{5%IntFormat, 2%IntFormat, SignedFormat, FiniteFormat}
sf52 = sf52type()
sf52b = Binary(5%IntFormat, 2%IntFormat, SIGNED, FINITE)
sf52 == sf52b

ue43type = Binary{4%IntFormat, 3%IntFormat, UnsignedFormat, ExtendedFormat}
ue43 = ue43type()
ue43b = Binary(4%IntFormat, 3%IntFormat, UNSIGNED, EXTENDED)
ue43 == ue43b

println(string("\t type", "\t\t => \t", "struct"))
println(string("\t BinaryFormat", "\t => \t", "Bunary","\n"))
println(string("    ", sf52type, "\t =>     ", sf52))
println(string("    ", ue43type, "\t =>     ", ue43))

end  # Float3109s

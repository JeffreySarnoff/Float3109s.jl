module Float3109s

include("constants.jl")

include("types/traits/signedness.jl")
include("types/traits/domain.jl")

include("types/Binary.jl")
include("types/Datum.jl")




sf52type = Binary{5%Int16, 2%Int16, SignedFormat, FiniteFormat}
sf52 = sf52type()

ue43type = Binary{4%Int16, 3%Int16, UnsignedFormat, ExtendedFormat}
ue43 = Binary(4%Int16, 3%Int16, UNSIGNED, EXTENDED)


(sf52type, sf52, ue43type, ue43)

end  # Float3109s

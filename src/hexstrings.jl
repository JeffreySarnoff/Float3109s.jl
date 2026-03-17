function AllHexStringValuesOf(fmt::Format{signedness,domain}) where {signedness,domain}
    values = AllValuesOf(fmt)
    map(hex_sprintf, values)
end

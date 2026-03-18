const FmtKinds = (;
    uf=(UnsignedFormat, FiniteFormat),
    ue=(UnsignedFormat, ExtendedFormat),
    sf=(SignedFormat, FiniteFormat),
    se=(SignedFormat, ExtendedFormat))

codetype(K) = K <= 8 ? UInt8 : K <= 16 ? UInt16 : K <= 32 ? UInt32 : UInt64
hex_codepoint(K, x) = K <= 8 ? @sprintf("0x%02x", x) : @sprintf("0x%04x", x)

for K in MinK:MaxK
    codepoints = map(x -> hex_codepoint(K, x), collect(0:twopow(K)-1))

    subdir = string("K", K)
    subpath = joinpath(HexStringsBase, subdir)
    for P in 1:K
        subsubdir = string("P", P)
        subsubpath = joinpath(subpath, subsubdir)
        basename = string("Binary", K, "p", P)
        fournames = (;
            uf=string(basename, "uf", Ext), ue=string(basename, "ue", Ext),
            sf=string(basename, "sf", Ext), se=string(basename, "se", Ext))
        fourpaths = (;
            uf=joinpath(subsubpath, fournames.uf), ue=joinpath(subsubpath, fournames.ue),
            sf=joinpath(subsubpath, fournames.sf), se=joinpath(subsubpath, fournames.se))

        tables = []

        for fmtkind in keys(fourpaths)
            signedness, domain = FmtKinds[Symbol(fmtkind)]
            if P < K
                fmt = Format{signedness,domain}(K, P)
            else
                fmt = Format{UnsignedFormat,domain}(K, P)
            end
            hexstring_values = AllHexStringValuesOf(fmt)
            localtable = columntable((; codepoint=codepoints, hexstring=hexstring_values))
            # writetable(fourpaths[Symbol(fmtkind)], localtable)
            CSV.write(fourpaths[Symbol(fmtkind)], localtable)
        end
    end
end

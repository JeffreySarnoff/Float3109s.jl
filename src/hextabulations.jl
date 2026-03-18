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
    subpath = aspath(HexStringsBase, subdir)
    for P in 1:K-1
        subsubdir = string("P", P)
        subsubpath = aspath(subpath, subsubdir)
        basename = string("Binary", K, "p", P)
        fournames = (;
            uf=string(basename, "uf", Ext), ue=string(basename, "ue", Ext),
            sf=string(basename, "sf", Ext), se=string(basename, "se", Ext))
        fourpaths = (;
            uf=aspath(subsubpath, fournames.uf), ue=aspath(subsubpath, fournames.ue),
            sf=aspath(subsubpath, fournames.sf), se=aspath(subsubpath, fournames.se))

        for fmtkind in keys(fourpaths)
            signedness, domain = FmtKinds[Symbol(fmtkind)]
            fmt = Format{signedness,domain}(K, P)
            hexstring_values = AllHexStringValuesOf(fmt)
            subnormal_icons = AllSubnormalIconsOf(fmt)
            localtable = columntable((; codepoint=codepoints, sprintf_a=hexstring_values, subnormal=subnormal_icons))
            CSV.write(fourpaths[Symbol(fmtkind)], localtable)
        end
    end

    for P in K:K
        subsubdir = string("P", P)
        subsubpath = aspath(subpath, subsubdir)
        basename = string("Binary", K, "p", P)

        twonames = (;
            uf=string(basename, "uf", Ext), ue=string(basename, "ue", Ext))
        twopaths = (;
            uf=aspath(subsubpath, twonames.uf), ue=aspath(subsubpath, twonames.ue))

        for fmtkind in keys(twopaths)
            signedness, domain = FmtKinds[Symbol(fmtkind)]
            fmt = Format{signedness,domain}(K, P)
            subnormal_icons = AllSubnormalIconsOf(fmt)
            hexstring_values = AllHexStringValuesOf(fmt)
            localtable = columntable((; codepoint=codepoints, sprintf_a=hexstring_values, subnormal=subnormal_icons))
            CSV.write(twopaths[Symbol(fmtkind)], localtable)
        end
    end
end

const SubnormalIcon = '↯'

function AllSubnormalIconsOf(fmt::Format)
    subnormal_icons = fill(' ', nValuesOf(fmt))
    PrecsionOf(fmt) == 1 && return subnormal_icons

    npos_subnormals = nPosSubnormalsOf(fmt)
    nneg_subnormals = is_signed(fmt) ? nNegSubnormalsOf(fmt) : 0

    cpnan = cp_nan(fmt)

    cp_firstpos_subnormal = 0x01
    cp_lastpos_subnormal = npos_subnormals

    cp_firstneg_subnormal = is_signed(fmt) ? cpnan + 0x01 : -1
    cp_lastneg_subnormal = is_signed(fmt) ? cpnan + nneg_subnormals : -1

    subnormal_icons[cp_firstpos_subnormal:cp_lastpos_subnormal] .= SubnormalIcon
    if is_signed(fmt)
        subnormal_icons[cp_firstneg_subnormal:cp_lastneg_subnormal] .= SubnormalIcon
    end

    subnormal_icons
end
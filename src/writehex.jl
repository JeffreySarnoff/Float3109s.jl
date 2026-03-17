const HexBase = joinpath(P3109Base, "HexStrings/")

for k in 3:16
    mkpath(joinpath(HexBase, "K$k"))
    for p in 1:k
        mkpath(joinpath(HexBase, "K$k", "P$p"))
    end
end


function save_uf_table(ufpath, uftable)
    CSV.write(ufpath, uftable)
end

function get_uf_table(codepoints, K, P, cp_nan)
    values =
        table = Vector{Tuple{UInt64,String}}()
    for cp in UInt64(0):UInt64(cp_nan - 1)
        finite_vals = AllFiniteValuesOf(Format{UnsignedFormat,FiniteFormat}(K, P))
        finite_strs = map(hex_sprintf, finite_vals)
        strs = push!(finite_strs, "NaN")
        vals = finite_valspush!(table, (cp, hex_sprintf(val)))
    end
    return columntable((; codepoint=vals, hexstring=strs))
end

for k in 3:16
    finalcodepoint = 0x0002^k - 0x0001
    codepoints = collect(k <= 8 ? collect(0x00:UInt8(finalcodepoint)) : collect(0x0000:UInt16(finalcodepoint)))
    cp_unsigned_nan = finalcodepoint
    cp_signed_nan = finalcodepoint >> 0x01
    K = k
    for p in 1:(k-1)
        P = p
        currpath = joinpath(HexBase, "K$k", "P$p")
        filebasename = "Binary$(k)p($p)"
        ufpath = joinpath(currpath, string(filebasename, "uf", ".csv"))
        uepath = joinpath(currpath, string(filebasename, "ue", ".csv"))
        sfpath = joinpath(currpath, string(filebasename, "sf", ".csv"))
        sepath = joinpath(currpath, string(filebasename, "se", ".csv"))

        uftable = invoke_latest(get_uf_table(codepoints, K, P, cp_unsigned_nan))
        save_uf_table(ufpath, uftable)
        uetable = get_ue_table(codepoints, K, P, cp_unsigned_nan)
        sftable = get_sf_table(codepoints, K, P, cp_signed_nan)
        setable = get_se_table(codepoints, K, P, cp_signed_nan)
        save_ue_table(uepath, uetable)
        save_sf_table(sfpath, sftable)
        save_se_table(sepath, setable)
    end
    K = k
    P = k
    currpath = joinpath(HexBase, "K$k", "P$k")
    filebasename = "Binary$(k)p($k)"
    ufpath = joinpath(currpath, string(filebasename, "uf", ".csv"))
    uepath = joinpath(currpath, string(filebasename, "ue", ".csv"))
    uftable = get_uf_table(codepoints, K, P, cp_unsigned_nan)
    uetable = get_ue_table(K, P, cp_unsigned_nan)
    save_uf_table(ufpath, uftable)
    save_ue_table(uepath, uetable)
end


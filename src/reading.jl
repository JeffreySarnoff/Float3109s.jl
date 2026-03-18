function read_hexstring_csv(K, P, signedness::Sugnedness, domain::Domain)
    basedir = HexStringsBase
    subdir = string("K", K)
    subsubdir = string("P", P)
    csvdir = aspath(basedir, subdir, subsubdir)

    filebase = string("Binary", K, "p", P)
    suffix = string(signedness == UnsignedFormat ? "u" : "s", domain == FiniteFormat ? "f" : "e")
    filename = string(filebase, suffix, Ext)

    csvfile = aspath(csvdir, filename)

    fmt = Format{signedness,domain}(K, P)
    path = aspath(HexStringsBase, string("K", K), string("P", P), string("Binary", K, "p", P, signedness == UnsignedFormat ? "u" : "s", domain == FiniteFormat ? "f" : "e", Ext))
    df = CSV.File(path) |> DataFrame
    return df
end

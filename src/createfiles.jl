# P3109Base = s"C:/git/P3109/"
const HexStringsBase = aspath(P3109Base, "HexStrings")

for K in MinK:MaxK
    subdir = string("K", K)
    subpath = aspath(HexStringsBase, subdir)
    if !isdir(subpath)
        @info "Creating directory $subpath"
        mkpath(subpath)
    end
end

for K in MinK:MaxK
    subdir = string("K", K)
    subpath = aspath(HexStringsBase, subdir)
    for P in 1:K-1
        subsubdir = string("P", P)
        subsubpath = aspath(subpath, subsubdir)
        if !isdir(subsubpath)
            @info "Creating directory $subsubpath"
            mkpath(subsubpath)
        end
        basename = string("Binary", K, "p", P)
        fournames = (; uf=string(basename, "uf", Ext), ue=string(basename, "ue", Ext),
            sf=string(basename, "sf", Ext), se=string(basename, "se", Ext))
        fourpaths = (; uf=aspath(subsubpath, fournames.uf), ue=aspath(subsubpath, fournames.ue),
            sf=aspath(subsubpath, fournames.sf), se=aspath(subsubpath, fournames.se))
        for path in values(fourpaths)
            if !isfile(path)
                @info "Creating file $path"
                open(path, "w") do io
                    write(io, "")
                end
            end
        end
    end
    for P in K:K
        subsubdir = string("P", P)
        subsubpath = aspath(subpath, subsubdir)
        if !isdir(subsubpath)
            @info "Creating directory $subsubpath"
            mkpath(subsubpath)
        end
        basename = string("Binary", K, "p", P)
        twonames = (; uf=string(basename, "uf", Ext), ue=string(basename, "ue", Ext))
        twopaths = (; uf=aspath(subsubpath, twonames.uf), ue=aspath(subsubpath, twonames.ue))
        for path in values(twopaths)
            if !isfile(path)
                @info "Creating file $path"
                open(path, "w") do io
                    write(io, "")
                end
            end
        end
    end
end


const DecFmtKinds = (;
    uf=(UnsignedFormat, FiniteFormat),
    ue=(UnsignedFormat, ExtendedFormat),
    sf=(SignedFormat, FiniteFormat),
    se=(SignedFormat, ExtendedFormat))

const DecSubnormalMarker = "!"
const Log2DenomColumn = :log2_denom

const StrBig = Union{String,BigInt,Int64}
const VecStrBig = Vector{StrBig}
const UINT = Union{Int,UInt8,UInt16}
const VecUI = Vector{UINT}
const VecStr = Vector{String}

function DataFrameOf(codepoints::Vector{String}, numerators::Vector{Union{String,BigInt}}, log2denoms::Vector{Union{BigInt,String}}, subnormal_markers::Vector{String})
    DataFrame(
        codepoint=codepoints,
        numerator=numerators,
        log2_denominator=log2denoms,
        subnormal=subnormal_markers
    )
end

dec_codepoint(K, x) = K <= 8 ? @sprintf("0x%02x", x) : @sprintf("0x%04x", x)

function AllSubnormalMarkersOf(fmt::Format)
    markers = fill(" ", Int(nValuesOf(fmt)))

    cp_pos_min = cp_pos_subnormal_min(fmt)
    cp_pos_max = cp_pos_subnormal_max(fmt)
    if cp_pos_min !== nothing && cp_pos_max !== nothing
        markers[Int(cp_pos_min)+1:Int(cp_pos_max)+1] .= DecSubnormalMarker
    end

    if is_signed(fmt)
        cp_neg_min = cp_neg_subnormal_min(fmt)
        cp_neg_max = cp_neg_subnormal_max(fmt)
        if cp_neg_min !== nothing && cp_neg_max !== nothing
            markers[Int(cp_neg_min)+1:Int(cp_neg_max)+1] .= DecSubnormalMarker
        end
    end

    markers
end

function NumeratorAndLog2DenomColumnsOf(fmt::Format)
    n = Int(nValuesOf(fmt))
    numerators = Vector{Union{String,BigInt}}(undef, n)
    log2denoms = Vector{Union{String,BigInt}}(undef, n)

    cpnan = cp_nan(fmt)
    cpinf = cp_inf(fmt)
    cpneginf = cp_neginf(fmt)

    for cp in Int128(0):Int128(cp_max(fmt))
        idx = Int(cp) + 1

        if cp == cpnan
            numerators[idx] = "NaN"
            log2denoms[idx] = zero(BigInt)
            continue
        elseif cp == cpinf
            numerators[idx] = "Inf"
            log2denoms[idx] = zero(BigInt)
            continue
        elseif cp == cpneginf
            numerators[idx] = "NegInf"
            log2denoms[idx] = zero(BigInt)
            continue
        end

        q = FiniteValueOf(fmt, cp)
        den = denominator(q)

        numerators[idx] = numerator(q)
        log2denoms[idx] = ispow2(den) ? BigInt(trailing_zeros(den)) : error("Denominator (fmt=$(fmt), codepoint=$(idx-1)) is not a power of 2")
    end

    return (; numerators, log2denoms)
end

for K in MinK:MaxK
    codepoints = map(x -> dec_codepoint(K, x), collect(0:twopow(K)-1))

    subdir = string("K", K)
    subpath = aspath(RationalsBase, subdir)
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
            path = fourpaths[fmtkind]
            signedness, domain = DecFmtKinds[Symbol(fmtkind)]
            fmt = Format{signedness,domain}(K, P)
            comps = NumeratorAndLog2DenomColumnsOf(fmt)
            subnormal_markers = AllSubnormalMarkersOf(fmt)
            #localtable = NamedTuple{
            #    (:codepoint, :numerator, :log2_denominator, :subnormal)
            #}((codepoints, comps.numerators, comps.log2denoms, subnormal_markers))
            df = DataFrameOf(codepoints, comps.numerators, comps.log2denoms, subnormal_markers)
            CSV.write(path, df)
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
            path = twopaths[fmtkind]
            signedness, domain = DecFmtKinds[Symbol(fmtkind)]
            fmt = Format{signedness,domain}(K, P)
            comps = NumeratorAndLog2DenomColumnsOf(fmt)
            subnormal_markers = AllSubnormalMarkersOf(fmt)
            #localtable = NamedTuple{
            #    (:codepoint, :numerator, Log2DenomColumn, :subnormal)
            #}((codepoints, comps.numerators, comps.log2denoms, subnormal_markers))
            # CSV.write(twopaths[Symbol(fmtkind)], localtable)
            df = DataFrameOf(codepoints, comps.numerators, comps.log2denoms, subnormal_markers)
            CSV.write(path, df)
        end
    end
end

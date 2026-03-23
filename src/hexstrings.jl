hex_sprintf(x) = @sprintf2("%a", BigFloat(x))

function AllHexStringValuesOf(fmt::Format{signedness,domain}) where {signedness,domain}
    values = AllValuesOf(fmt)
    hex_values = Vector{String}(undef, length(values))

    cp_pos_min = cp_pos_subnormal_min(fmt)
    cp_pos_max = cp_pos_subnormal_max(fmt)
    cp_neg_min = cp_neg_subnormal_min(fmt)
    cp_neg_max = cp_neg_subnormal_max(fmt)

    cp_pos_min_i = cp_pos_min === nothing ? nothing : Int128(cp_pos_min)
    cp_pos_max_i = cp_pos_max === nothing ? nothing : Int128(cp_pos_max)
    cp_neg_min_i = cp_neg_min === nothing ? nothing : Int128(cp_neg_min)
    cp_neg_max_i = cp_neg_max === nothing ? nothing : Int128(cp_neg_max)

    for (i, value) in enumerate(values)
        cp = Int128(i - 1)
        is_subnormal = false

        if cp_pos_min_i !== nothing && cp_pos_max_i !== nothing
            is_subnormal = cp_pos_min_i <= cp <= cp_pos_max_i
        end
        if !is_subnormal && cp_neg_min_i !== nothing && cp_neg_max_i !== nothing
            is_subnormal = cp_neg_min_i <= cp <= cp_neg_max_i
        end

        hex_values[i] = sprintf2a(BigFloat(value); subnormal=is_subnormal)
    end

    hex_values
end

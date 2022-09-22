is_nonnegative_integers(x::Integer) = x ≥ 0
is_nonnegative_integers(x) = false

function fill_missing_bins(df::DataFrame, index::Symbol, fill_pairs...)
    empty_df = DataFrame()
    empty_df[!, index] = 1:maximum(df[!, index])
    filled_df = leftjoin(empty_df, df, on=index)
    for (k, v) in fill_pairs
        filled_df[:, k] .= coalesce.(filled_df[!, k], v)
    end
    sort!(filled_df, index)
    return filled_df
end

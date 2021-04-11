@inline function combine_fuzzysets(fs::AbstractMatrix{T}, op_ratio) where {T}
    return op_ratio .* fuzzyset_union(fs) .+
           (one(T) - op_ratio) .* fuzzyset_intersection(fs)
end

@inline fuzzyset_union(fs::AbstractMatrix) = fs .+ fs' .- (fs .* fs')
@inline fuzzyset_intersection(fs::AbstractMatrix) = fs .* fs'

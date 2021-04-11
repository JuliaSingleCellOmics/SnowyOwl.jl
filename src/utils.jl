using StatsBase

function argsort(x::AbstractVector)
    rank = StatsBase.ordinalrank(x)
    invpermute!(collect(1:length(x)), rank)
end

function argsort(X::AbstractMatrix; dims=1)
    if dims == 1
        idxs = [argsort(view(X, :, j)) for j = 1:size(X,2)]
        return hcat(idxs...)
    elseif dims == 2
        idxs = [argsort(view(X, i, :))' for i = 1:size(X,1)]
        return vcat(idxs...)
    else
        throw(ArgumentError("Invalid dims, dims only accepts 1 or 2."))
    end
end

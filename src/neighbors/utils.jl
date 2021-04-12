function argsort(x::AbstractVector)
    rank = StatsBase.ordinalrank(x)
    invpermute!(collect(1:length(x)), rank)
end

function argsort(X::AbstractMatrix; dims=1)
    res = similar(X, Int)
    if dims == 1
        for j = 1:size(X,2)
            view(res, :, j) .= argsort(view(X, :, j))
        end
    elseif dims == 2
        for i = 1:size(X,1)
            view(res, i, :) .= argsort(view(X, i, :))
        end
    else
        throw(ArgumentError("Invalid dims, dims only accepts 1 or 2."))
    end
    return res
end

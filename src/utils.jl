sortperm(V::AbstractMatrix; kwargs...) = sortperm!(similar(V, Int), V; kwargs...)

function sortperm!(X::AbstractMatrix{<:Integer}, V::AbstractMatrix; dims::Int=1, kwargs...)
    if dims == 1
        for j in 1:size(X,2)
            sortperm!(view(X, :, j), view(V, :, j), kwargs...)
        end
    elseif dims == 2
        for i in 1:size(X,1)
            sortperm!(view(X, i, :), view(V, i, :), kwargs...)
        end
    else
        throw(ArgumentError("Invalid dims, dims only accepts 1 or 2."))
    end

    return X
end

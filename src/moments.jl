normalize_indicator(i::AbstractArray) = i ./ sum(i)
normalize_indicator(i::AbstractArray, dims) = i ./ sum(i, dims=dims)

to_indicator_matrix(i::AbstractVector) = diagm(i)

function union_diagonal!(A::AbstractMatrix{T}) where {T}
    for i in 1:min(size(A)...)
        A[i, i] = one(T)
    end
    return A
end

_generate_slices(n, k) = Base.ntuple(i->(i==n) ? k : Colon(), 2)

"""
    graph_filter(ind, order)

Generate a graph filter of `order`-th order moment depending on the indicator `ind`.

## Arguments
- `ind`: the indicator vector or matrix with values of 0s or 1s.
- `order`: the order of moment for a graph filter. If `order=1`, generates a first moment graph filter,
which is the mean operator. If `order=2`, generates a second moment graph filter, which is the covariance operator.
"""
function graph_filter(ind::AbstractVector, order::Integer)
    if order == 1
        return sparse(normalize_indicator(ind))
    elseif order == 2
        adj = to_indicator_matrix(ind)
        return sparse(normalize_indicator(adj))
    else
        throw(ArgumentError("order other than 1 and 2 is not supported while get $order."))
    end
end

function graph_filter(adj::AbstractMatrix, order::Integer; dims::Integer=1, self::Bool=true)
    if order == 1
        self && union_diagonal!(adj)
        return sparse(normalize_indicator(adj, dims))
    elseif order == 2
        n = dims%2 + 1
        blocks = [graph_filter(vec(adj[_generate_slices(n, k)...]), 2) for k in 1:size(adj, n)]
        return cat(blocks...,dims=(1,2))
    else
        throw(ArgumentError("order other than 1 and 2 is not supported while get $order."))
    end
end

function moment(A::AbstractMatrix, X::AbstractMatrix, order::Integer; dims::Integer=1)
    if order == 1
        return first_moment(A, X; dims=dims)
    elseif order == 2
        return second_moment(A, X; dims=dims)
    else
        throw(ArgumentError("order other than 1 and 2 is not supported while get $order."))
    end
end

function first_moment(A::AbstractMatrix, X::AbstractMatrix; dims::Integer=1)
    adj = A .> 0
    if dims == 1
        G = graph_filter(adj, 1; dims=1)
        return G*X
    elseif dims == 2
        G = graph_filter(adj, 1; dims=2)
        return X*G
    else
        throw(ArgumentError("dims other than 1 and 2 is not supported while get $dims."))
    end
end

# second_moment(A::AbstractMatrix, X::AbstractMatrix; dims::Integer=1) = second_moment(A, X, X; dims=dims)

# function second_moment(A::AbstractMatrix, X::AbstractMatrix, Y::AbstractMatrix; dims::Integer=1)
#     # x' * G_2 * y - (G_1' * x)*(G_1' * y)
#     adj = A .> 0
#     if dims == 1
#         G = graph_filter(adj, 2; dims=1)
#         return X'*G*Y - first_moment(A, X; dims=1)*first_moment(A, Y; dims=1)'
#     elseif dims == 2
#         G = graph_filter(adj, 2; dims=2)
#         return vec(X')'*G*vec(Y') - first_moment(A, X; dims=2)*first_moment(A, Y; dims=2)'
#     else
#         throw(ArgumentError("dims other than 1 and 2 is not supported while get $dims."))
#     end
# end

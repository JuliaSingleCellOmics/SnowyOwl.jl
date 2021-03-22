using LinearAlgebra: diagm
using SparseArrays: sparse

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

function normalize_neighbors(C::AbstractMatrix)
    neighbor_graph = C .> 0
    union_diagonal!(neighbor_graph)
    neighbor_graph .= normalize_indicator(neighbor_graph, 2)
    return neighbor_graph
end

function moment(A::AbstractArray, order::Integer)
    if order == 1
        return first_moment(A)
    elseif order == 2
        return second_moment(A)
    else
        throw(ArgumentError("order other than 1 and 2 is not supported while get $order."))
    end
end

function first_moment(A::AbstractArray)
    
end

function second_moment(A::AbstractArray)
    
end

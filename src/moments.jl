using LinearAlgebra: diagm

normalize_indicator(i::AbstractArray) = i ./ sum(i)
normalize_indicator(i::AbstractArray, dims) = i ./ sum(i, dims=dims)

to_indicator_matrix(i::AbstractVector) = diagm(i)

function union_diagonal!(A::AbstractMatrix{T}) where {T}
    for i in 1:min(size(A)...)
        A[i, i] = one(T)
    end
    return A
end

function normalize_neighbors(C::AbstractMatrix)
    neighbor_graph = C .> 0
    union_diagonal!(neighbor_graph)
    neighbor_graph ./= sum(neighbor_graph, dims=2)
    return neighbor_graph
end

function moment(A::AbstractArray, order::Integer)
    
end

function first_moment(A::AbstractArray)
    
end

function second_moment(A::AbstractArray)
    
end

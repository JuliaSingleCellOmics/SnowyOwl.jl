using LinearAlgebra: diagm

normalize_indicator(i::AbstractArray) = i ./ sum(i)
to_indicator_matrix(i::AbstractVector) = diagm(i)

function normalize_neighbors(C::AbstractMatrix)
    neighbor_graph = C .> 0
    for i in 1:min(size(neighbor_graph)...)
        neighbor_graph[i, i] = true
    end
    neighbor_graph ./= sum(neighbor_graph, dims=2)
    return neighbor_graph
end

function moment(A::AbstractArray, order::Integer)
    
end

function first_moment(A::AbstractArray)
    
end

function second_moment(A::AbstractArray)
    
end

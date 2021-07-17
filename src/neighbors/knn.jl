knn_search(dist_mat, k, metric::Symbol) = knn_search(dist_mat, k, Val(metric))

"""
    knn_search(dist_mat, k, :precomputed) -> knns, dists
Find the `k` nearest neighbors of each point in a precomputed distance
matrix.
"""
knn_search(dist_mat, k, ::Val{:precomputed}) = _knn_from_dists(dist_mat, k)

"""
    knn_search(X, k, metric) -> knns, dists
Find the `k` nearest neighbors of each point in `X` by `metric`.
"""
function knn_search(X,
                    k,
                    metric::SemiMetric)
    if size(X, 2) < 4096
        return knn_search(X, k, metric, Val(:pairwise))
    else
        return knn_search(X, k, metric, Val(:approximate))
    end
end

# compute all pairwise distances
# return the nearest k to each point v, other than v itself
function knn_search(X::AbstractMatrix{S},
                    k,
                    metric,
                    ::Val{:pairwise}) where {S <: Real}
    num_points = size(X, 2)
    dist_mat = Array{S}(undef, num_points, num_points)
    pairwise!(dist_mat, metric, X, dims=2)
    # all_dists is symmetric distance matrix
    return _knn_from_dists(dist_mat, k)
end

# find the approximate k nearest neighbors using NNDescent
function knn_search(X::AbstractMatrix{S},
                    k,
                    metric,
                    ::Val{:approximate}) where {S <: Real}
    knngraph = nndescent(X, k, metric)
    return knn_matrices(knngraph)
end

function _knn_from_dists(dist_mat::AbstractMatrix{S}, k) where {S <: Real}
    knns_ = [partialsortperm(view(dist_mat, :, i), 2:k+1) for i in 1:size(dist_mat, 1)]
    dists_ = [dist_mat[:, i][knns_[i]] for i in eachindex(knns_)]
    knns = hcat(knns_...)::Matrix{Int}
    dists = hcat(dists_...)::Matrix{S}
    return knns, dists
end

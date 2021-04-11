function spectral_embedding(adjm, dim)
    # Deal with multiple connected components
    # need metric=:euclidean
    L = normalized_laplacian(adjm, Float64)
    U, Σ, V = svd(L, full=true)
    # eigs = eigen(L, 2:dim+1)
    # return eigs.vectors
end

function optimize_embedding()

end

function smooth_knn_dist(distm, k, n_iter=64, local_connectivity=1.0, bandwidth=1.0)

end

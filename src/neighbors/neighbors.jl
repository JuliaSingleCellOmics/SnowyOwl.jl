# neighbors = Neighbors(adata)
# neighbors.compute_neighbors(
#     n_pcs = None if use_rep == "X" else n_pcs,
#     use_rep = "X" if adata.n_vars < 50 or n_pcs == 0 else "X_pca",
#     random_state=0,
# )

# """
# "rapids" methods for cuML - RAPIDS Machine Learning Library is not supported.
# """
# function neighborhood_graph!(prof::Profile, n_neighbors::Integer=30, method::Symbol=:umap;
#                              use_knn=true, metric::Metric=Euclidean(), use_highly_variable=true)
#     X = _choose_representation(prof, use_rep=use_rep, n_pcs=n_pcs)
#     neighborhood_graph(X, n_neighbors, method; use_knn=use_knn, metric=metric)
#     # adata.obsp["distances"] = neighbors.distances
#     # adata.obsp["connectivities"] = neighbors.connectivities
#     # adata.uns["neighbors"]["connectivities_key"] = "connectivities"
#     # adata.uns["neighbors"]["distances_key"] = "distances"
# end

"""
https://github.com/theislab/scanpy/blob/c488909a54e9ab1462186cca35b537426e4630db/scanpy/neighbors/__init__.py#L720
"""
function neighborhood_graph(X::AbstractMatrix, n_neighbors::Integer=30, method=:umap; use_knn=True, metric::Metric=Euclidean())
    @assert method == :umap && use_knn "UMAP method supports only use_knn=true."
    @assert method in (:umap, :gaussian) "Unsupported method of $method, only support (:umap, :gaussian)."

    use_dense_distance = (metric isa Euclidean && size(X,2) < 8192) || !use_knn
    neighbor_search(X, n_neighbors, metric, Val(use_dense_distance))

    # if !use_dense_distance
    #     # we need self._distances also for method == 'gauss' if we didn't use dense distances
    #     self._distances, self._connectivities = compute_connectivities(X, Val(:umap), graph)
    #     # knn_indices,
    #     # knn_distances,
    #     # self._adata.shape[0],
    #     # self.n_neighbors
    # else
    #     self._distances, self._connectivities = compute_connectivities(X, Val(method), graph)
    # end

    # self._number_connected_components = connected_components(connectivities)
end

"""
exact knn graph
https://github.com/theislab/scanpy/blob/c488909a54e9ab1462186cca35b537426e4630db/scanpy/neighbors/__init__.py#L773
"""
function neighbor_search(X::AbstractMatrix, n_neighbors, metric::Metric, use_dense_distance::Val{true})
    graph = nndescent(X, n_neighbors, metric)
    indices, distances = knn_matrices(graph)

    dist = pairwise(metric, X, dims=2)
    use_knn && (dist = sparse(dist, n_neighbors))
    return dist
end

"""
approximate knn graph
https://github.com/theislab/scanpy/blob/c488909a54e9ab1462186cca35b537426e4630db/scanpy/neighbors/__init__.py#L786
# approximate nearest neighbors for non-euclidean case
"""
function neighbor_search(X::AbstractMatrix, n_neighbors, metric::Metric, use_dense_distance::Val{false})
    if size(X, 2) < 4096
        X = pairwise(metric, X, dims=2)
        knn_indices, knn_distances = knn_graph(X, Val(:precomputed), n_neighbors)
    else
        knn_indices, knn_distances = knn_graph(X, metric, n_neighbors)
    end
    return dist
end

"""
https://github.com/lmcinnes/umap/blob/ae5255be571c4a90bd93611a07612b4565c7984d/umap/umap_.py#L310
"""
function knn_graph(X::AbstractMatrix{T}, metric::Val{:precomputed}, n_neighbors::Integer=30) where T
    # Compute indices of n nearest neighbors
    knn_indices = SnowyOwl.argsort(X, dims=1)[2:(n_neighbors+1), :]
    # Compute the nearest neighbor distances
    knn_dists = sort(X, dims=1)[2:(n_neighbors+1), :]
    # Prune any nearest neighbours that are infinite distance apart.
    knn_dists[isinf.(knn_dists)] .= -one(T)

    return knn_indices, knn_dists
end

# """
# https://github.com/lmcinnes/umap/blob/ae5255be571c4a90bd93611a07612b4565c7984d/umap/umap_.py#L323
# """
# function knn_graph(X::AbstractMatrix, metric::Metric=Euclidean(), n_neighbors::Integer=30)
#     n_obs = size(X, 2)

#     n_trees = min(64, 5 + int(round(n_obs^0.5 / 20.0)))
#     n_iters = max(5, round(Int, log2(n_obs)))

#     knn_search_index = NNDescent(
#         X,
#         n_neighbors=n_neighbors,
#         metric=metric,
#         n_trees=n_trees,
#         n_iters=n_iters,
#         max_candidates=60,
#         low_memory=low_memory,
#     )
#     knn_indices, knn_dists = knn_search_index.neighbor_graph

#     return knn_indices, knn_dists
# end

# function compute_connectivities(X::AbstractMatrix, method::Val{:umap}, graph)
#     # knn_indices,
#     # knn_dists,
#     # n_obs,
#     # n_neighbors,
#     # set_op_mix_ratio=1.0,
#     # local_connectivity=1.0,

#     """\
#     This is from umap.fuzzy_simplicial_set [McInnes18]_.
#     Given a set of data X, a neighborhood size, and a measure of distance
#     compute the fuzzy simplicial set (here represented as a fuzzy graph in
#     the form of a sparse matrix) associated to the data. This is done by
#     locally approximating geodesic distance at each point, creating a fuzzy
#     simplicial set for each such point, and then combining all the local
#     fuzzy simplicial sets into a global one via a fuzzy union.
#     """

#     #     from umap.umap_ import fuzzy_simplicial_set

#     # X = coo_matrix(([], ([], [])), shape=(n_obs, 1))
#     # connectivities = fuzzy_simplicial_set(
#     #     X,
#     #     n_neighbors,
#     #     None,
#     #     None,
#     #     knn_indices=knn_indices,
#     #     knn_dists=knn_dists,
#     #     set_op_mix_ratio=set_op_mix_ratio,
#     #     local_connectivity=local_connectivity,
#     # )

#     # distances = _get_sparse_matrix_from_indices_distances_umap(
#     #     knn_indices, knn_dists, n_obs, n_neighbors
#     # )

#     # return distances, connectivities.tocsr()

# end

# function compute_connectivities(X::AbstractMatrix, method::Val{:gaussian}, graph)
#     # overwrite the umap connectivities if method is 'gauss'

#     # self is class Neighbors
#     # density_normalize=True

#     # # init distances
#     # if self.knn:
#     #     Dsq = self._distances.power(2)
#     #     indices, distances_sq = _get_indices_distances_from_sparse_matrix(
#     #         Dsq, self.n_neighbors
#     #     )
#     # else:
#     #     Dsq = np.power(self._distances, 2)
#     #     indices, distances_sq = _get_indices_distances_from_dense_matrix(
#     #         Dsq, self.n_neighbors
#     #     )

#     # # exclude the first point, the 0th neighbor
#     # indices = indices[:, 1:]
#     # distances_sq = distances_sq[:, 1:]

#     # # choose sigma, the heuristic here doesn't seem to make much of a difference,
#     # # but is used to reproduce the figures of Haghverdi et al. (2016)
#     # if self.knn:
#     #     # as the distances are not sorted
#     #     # we have decay within the n_neighbors first neighbors
#     #     sigmas_sq = np.median(distances_sq, axis=1)
#     # else:
#     #     # the last item is already in its sorted position through argpartition
#     #     # we have decay beyond the n_neighbors neighbors
#     #     sigmas_sq = distances_sq[:, -1] / 4
#     # sigmas = np.sqrt(sigmas_sq)

#     # # compute the symmetric weight matrix
#     # if not issparse(self._distances):
#     #     Num = 2 * np.multiply.outer(sigmas, sigmas)
#     #     Den = np.add.outer(sigmas_sq, sigmas_sq)
#     #     W = np.sqrt(Num / Den) * np.exp(-Dsq / Den)
#     #     # make the weight matrix sparse
#     #     if not self.knn:
#     #         mask = W > 1e-14
#     #         W[~mask] = 0
#     #     else:
#     #         # restrict number of neighbors to ~k
#     #         # build a symmetric mask
#     #         mask = np.zeros(Dsq.shape, dtype=bool)
#     #         for i, row in enumerate(indices):
#     #             mask[i, row] = True
#     #             for j in row:
#     #                 if i not in set(indices[j]):
#     #                     W[j, i] = W[i, j]
#     #                     mask[j, i] = True
#     #         # set all entries that are not nearest neighbors to zero
#     #         W[~mask] = 0
#     # else:
#     #     W = (
#     #         Dsq.copy()
#     #     )  # need to copy the distance matrix here; what follows is inplace
#     #     for i in range(len(Dsq.indptr[:-1])):
#     #         row = Dsq.indices[Dsq.indptr[i] : Dsq.indptr[i + 1]]
#     #         num = 2 * sigmas[i] * sigmas[row]
#     #         den = sigmas_sq[i] + sigmas_sq[row]
#     #         W.data[Dsq.indptr[i] : Dsq.indptr[i + 1]] = np.sqrt(num / den) * np.exp(
#     #             -Dsq.data[Dsq.indptr[i] : Dsq.indptr[i + 1]] / den
#     #         )
#     #     W = W.tolil()
#     #     for i, row in enumerate(indices):
#     #         for j in row:
#     #             if i not in set(indices[j]):
#     #                 W[j, i] = W[i, j]
#     #     W = W.tocsr()

#     # self._connectivities = W

# end

# "analyze the connected components of a sparse graph"
# connected_components(::DenseArray) = 1
# function connected_components(conn::SparseMatrixCSC)
#     graph = SimpleGraph(conn .> 0)
#     cc = LightGraphs.connected_components(graph)
#     n_cc = length(cc)
#     return n_cc, cc
# end

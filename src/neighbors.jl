using Distances

# use_highly_variable=True,
# metric_kwds=None,
# num_threads=-1,
# copy=False,


# neighbors = Neighbors(adata)
# neighbors.compute_neighbors(
#     n_neighbors=30,
#     knn=True,
#     n_pcs = None if use_rep == "X" else n_pcs,
#     method="umap",
#     use_rep = "X" if adata.n_vars < 50 or n_pcs == 0 else "X_pca",
#     random_state=0,
#     metric="euclidean",
#     metric_kwds=metric_kwds,
#     write_knn_indices=True,
# )

# adata.obsp["distances"] = neighbors.distances
# adata.obsp["connectivities"] = neighbors.connectivities
# adata.uns["neighbors"]["connectivities_key"] = "connectivities"
# adata.uns["neighbors"]["distances_key"] = "distances"

function neighborhood_graph!(X, method; use_knn=True, metric::Metric=Euclidean)
    # method in ['umap', 'gauss', 'rapids']
    # cuML - RAPIDS Machine Learning Library
    # R = pairwise(dist, X, dims=2)

    # neighbor search
    use_dense_dist = (metric == Euclidean && X.shape[0] < 8192) || !use_knn
    if use_dense_dist
        _distances = pairwise_distances(X, metric=metric)
        knn_indices, knn_distances = _get_indices_distances_from_dense_matrix(_distances, n_neighbors)
        self._distances = use_knn ? _get_sparse_matrix_from_indices_distances_numpy(knn_indices, knn_distances, X.shape[0], n_neighbors) : _distances
    elseif method == "rapids"
        knn_indices, knn_distances = compute_neighbors_rapids(X, n_neighbors)
    else # non-euclidean case and approximate nearest neighbors
        if X.shape[0] < 4096
            X = pairwise_distances(X, metric=metric)
            metric = "precomputed"
        end
        knn_indices, knn_distances, forest = compute_neighbors_umap(X, n_neighbors, random_state, metric=metric, metric_kwds=metric_kwds)
    end

    # compute connectivities
    if !use_dense_dist || method in ["umap", "rapids"]
        # we need self._distances also for method == 'gauss' if we didn't use dense distances
        self._distances, self._connectivities = _compute_connectivities_umap(
            knn_indices,
            knn_distances,
            self._adata.shape[0],
            self.n_neighbors,
        )
    end

    if method == "gauss"  # overwrite the umap connectivities if method is 'gauss'
        self._compute_connectivities_diffmap()
    end

    if issparse(self._connectivities)
        self._connected_components = connected_components(self._connectivities)  # Analyze the connected components of a sparse graph
        self._number_connected_components = self._connected_components[0]
    else
        self._number_connected_components = 1
    end
end


def compute_neighbors(
    self,
    n_neighbors: int = 30,
    knn: bool = True,
    n_pcs: Optional[int] = None,
    use_rep: Optional[str] = None,
    method: _Method = 'umap', # in {'umap', 'gauss', 'rapids'}
    random_state: AnyRandom = 0,
    write_knn_indices: bool = False,
    metric: _Metric = 'euclidean',
    metric_kwds: Mapping[str, Any] = MappingProxyType({}),
):
    self._rp_forest = None
    X = _choose_representation(self._adata, use_rep=use_rep, n_pcs=n_pcs)
    # neighbor search
    # computed connectivities



def compute_neighbors_umap(
    X: Union[np.ndarray, csr_matrix],
    n_neighbors: int,
    random_state: AnyRandom = None,
    metric: Union[_Metric, _MetricFn] = 'euclidean',
    metric_kwds: Mapping[str, Any] = MappingProxyType({}),
    angular: bool = False,
    verbose: bool = False,
):
    from umap.umap_ import nearest_neighbors # umap 0.5.0

    knn_indices, knn_dists, forest = nearest_neighbors(
        X,
        n_neighbors,
        random_state=random_state,
        metric=metric,
        metric_kwds=metric_kwds,
        angular=angular,
        verbose=verbose,
    )

    return knn_indices, knn_dists, forest


def compute_neighbors_rapids(X: np.ndarray, n_neighbors: int):
    from cuml.neighbors import NearestNeighbors

    nn = NearestNeighbors(n_neighbors=n_neighbors)
    X_contiguous = np.ascontiguousarray(X, dtype=np.float32)
    nn.fit(X_contiguous)
    knn_distsq, knn_indices = nn.kneighbors(X_contiguous)
    return knn_indices, np.sqrt(knn_distsq)  # cuml uses sqeuclidean metric so take sqrt


def _get_indices_distances_from_dense_matrix(D, n_neighbors: int):
    sample_range = np.arange(D.shape[0])[:, None]
    indices = np.argpartition(D, n_neighbors - 1, axis=1)[:, :n_neighbors]
    indices = indices[sample_range, np.argsort(D[sample_range, indices])]
    distances = D[sample_range, indices]
    return indices, distances

def _get_sparse_matrix_from_indices_distances_numpy(indices, distances, n_obs, n_neighbors):
    n_nonzero = n_obs * n_neighbors
    indptr = np.arange(0, n_nonzero + 1, n_neighbors)
    D = csr_matrix(
        (distances.ravel(), indices.ravel(), indptr),
        shape=(n_obs, n_obs),
    )
    D.eliminate_zeros()
    return D


def _compute_connectivities_diffmap(self, density_normalize=True):
    # init distances
    if self.knn:
        Dsq = self._distances.power(2)
        indices, distances_sq = _get_indices_distances_from_sparse_matrix(
            Dsq, self.n_neighbors
        )
    else:
        Dsq = np.power(self._distances, 2)
        indices, distances_sq = _get_indices_distances_from_dense_matrix(
            Dsq, self.n_neighbors
        )

    # exclude the first point, the 0th neighbor
    indices = indices[:, 1:]
    distances_sq = distances_sq[:, 1:]

    # choose sigma, the heuristic here doesn't seem to make much of a difference,
    # but is used to reproduce the figures of Haghverdi et al. (2016)
    if self.knn:
        # as the distances are not sorted
        # we have decay within the n_neighbors first neighbors
        sigmas_sq = np.median(distances_sq, axis=1)
    else:
        # the last item is already in its sorted position through argpartition
        # we have decay beyond the n_neighbors neighbors
        sigmas_sq = distances_sq[:, -1] / 4
    sigmas = np.sqrt(sigmas_sq)

    # compute the symmetric weight matrix
    if not issparse(self._distances):
        Num = 2 * np.multiply.outer(sigmas, sigmas)
        Den = np.add.outer(sigmas_sq, sigmas_sq)
        W = np.sqrt(Num / Den) * np.exp(-Dsq / Den)
        # make the weight matrix sparse
        if not self.knn:
            mask = W > 1e-14
            W[~mask] = 0
        else:
            # restrict number of neighbors to ~k
            # build a symmetric mask
            mask = np.zeros(Dsq.shape, dtype=bool)
            for i, row in enumerate(indices):
                mask[i, row] = True
                for j in row:
                    if i not in set(indices[j]):
                        W[j, i] = W[i, j]
                        mask[j, i] = True
            # set all entries that are not nearest neighbors to zero
            W[~mask] = 0
    else:
        W = (
            Dsq.copy()
        )  # need to copy the distance matrix here; what follows is inplace
        for i in range(len(Dsq.indptr[:-1])):
            row = Dsq.indices[Dsq.indptr[i] : Dsq.indptr[i + 1]]
            num = 2 * sigmas[i] * sigmas[row]
            den = sigmas_sq[i] + sigmas_sq[row]
            W.data[Dsq.indptr[i] : Dsq.indptr[i + 1]] = np.sqrt(num / den) * np.exp(
                -Dsq.data[Dsq.indptr[i] : Dsq.indptr[i + 1]] / den
            )
        W = W.tolil()
        for i, row in enumerate(indices):
            for j in row:
                if i not in set(indices[j]):
                    W[j, i] = W[i, j]
        W = W.tocsr()

    self._connectivities = W

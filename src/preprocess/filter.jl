"""
Filter out cells which are not satisfy specific criteria.
"""
function filter_cells!(p::AnnotatedProfile; min_counts::Real=0, max_counts::Real=maximum(sum(p.data,dims=1)),
                       min_genes::Real=0, max_genes::Real=nrow(p))
      check_cell_args(min_counts, max_counts, min_genes, max_genes)
      return _filter_cells!(p, min_counts, max_counts, min_genes, max_genes)
end

function filter_cells(p::AnnotatedProfile; min_counts::Real=0, max_counts::Real=maximum(sum(p.data,dims=1)),
                      min_genes::Real=0, max_genes::Real=nrow(p))
      check_cell_args(min_counts, max_counts, min_genes, max_genes)
      return _filter_cells!(copy(p), min_counts, max_counts, min_genes, max_genes)
end

function _filter_cells!(p::AnnotatedProfile, min_counts::Real, max_counts::Real, min_genes::Real, max_genes::Real)
      expressed_genes = sum(p.data .!= 0, dims=1)
      bool_idx = min_genes .<= expressed_genes .<= max_genes
      library_size = sum(p.data, dims=1)
      bool_idx .&= min_counts .<= library_size .<= max_counts
      filter_layers!(p, obs_idx=vec(bool_idx))
      idx = findall(vec(bool_idx))
      p.data = p.data[:, idx]
      p.obs = p.obs[idx, :]
      return p
end

function check_cell_args(min_counts::Real, max_counts::Real, min_genes::Real, max_genes::Real)
      (min_counts < 0) && throw(ArgumentError("`min_counts` must not lower than zero but got $min_counts."))
      (min_genes < 0) && throw(ArgumentError("`min_genes` must not lower than zero but got $min_genes."))
      (max_genes < min_genes) &&
            throw(ArgumentError("`max_genes` must not lower than minimum number of genes ($min_genes) but got $max_genes."))
      (max_counts < min_counts) &&
            throw(ArgumentError("`max_counts` must not lower than minimum count ($min_counts) but got $max_counts."))
end

"""
Filter out genes which are not satisfy specific criteria.
"""
function filter_genes!(p::AnnotatedProfile; min_counts::Real=0, max_counts::Real=maximum(sum(p.data,dims=2)),
                       min_cells::Real=0, max_cells::Real=ncol(p))
      check_genes_args(min_counts, max_counts, min_cells, max_cells)
      return _filter_genes!(p, min_counts, max_counts, min_cells, max_cells)
end

function filter_genes(p::AnnotatedProfile; min_counts::Real=0, max_counts::Real=maximum(sum(p.data,dims=2)),
                      min_cells::Real=0, max_cells::Real=ncol(p))
      check_genes_args(min_counts, max_counts, min_cells, max_cells)
      return _filter_genes!(copy(p), min_counts, max_counts, min_cells, max_cells)
end

function _filter_genes!(p::AnnotatedProfile, min_counts::Real, max_counts::Real, min_cells::Real, max_cells::Real)
      expressed_cells = sum(p.data .!= 0, dims=2)
      bool_idx = min_cells .<= expressed_cells .<= max_cells
      expressed_counts = sum(p.data, dims=2)
      bool_idx .&= min_counts .<= expressed_counts .<= max_counts
      filter_layers!(p, var_idx=vec(bool_idx))
      idx = findall(vec(bool_idx))
      p.data = p.data[idx, :]
      p.var = p.var[idx, :]
      return p
end

function check_genes_args(min_counts::Real, max_counts::Real, min_cells::Real, max_cells::Real)
      (min_counts < 0) && throw(ArgumentError("`min_counts` must not lower than zero but got $min_counts."))
      (min_cells < 0) && throw(ArgumentError("`min_cells` must not lower than zero but got $min_cells."))
      (max_cells < min_cells) &&
            throw(ArgumentError("`max_cells` must not lower than minimum number of cells ($min_cells) but got $max_cells."))
      (max_counts < min_counts) &&
            throw(ArgumentError("`max_counts` must not lower than minimum count ($min_counts) but got $max_counts."))
end

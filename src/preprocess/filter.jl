function filter_cells!(p::Profile; min_counts::Real=0, min_genes::Real=0,
                       max_counts::Real=maximum(sum(p.data,dims=2)),
                       max_genes::Real=nrow(p))
      (min_counts < 0) && throw(ArgumentError("`min_counts` must not lower than zero but got $min_counts."))
      (min_genes < 0) && throw(ArgumentError("`min_genes` must not lower than zero but got $min_genes."))
      num_genes = nrow(p)
      (max_genes > num_genes) &&
            throw(ArgumentError("`max_genes` must not larger than max number of genes ($num_genes) but got $max_genes."))
      mtx_max_counts = maximum(sum(p.data, dims=2))
      (max_counts > mtx_max_counts) &&
            throw(ArgumentError("`max_counts` must not larger than max count of matrix ($mtx_max_counts) but got $max_counts."))
      
      expressed_genes = sum(p.data .!= 0, dims=2)
      idx = min_genes .<= expressed_genes .<= max_genes
      library_size = sum(p.data, dims=2)
      idx .&= min_counts .<= library_size .<= max_counts
      idx = findall(vec(idx))
      p.data = p.data[:, idx]
      p.obs = p.obs[idx, :]
      filter_layers!(p, obs_idx=idx)
      p
end

function filter_genes!(p::Profile; min_counts::Real=0, min_cells::Real=0,
                       max_counts::Real=maximum(sum(p.data,dims=1)),
                       max_cells::Real=ncol(p))
      (min_counts < 0) && throw(ArgumentError("`min_counts` must not lower than zero but got $min_counts."))
      (min_cells < 0) && throw(ArgumentError("`min_cells` must not lower than zero but got $min_cells."))
      num_cells = ncol(p)
      (max_cells > num_cells) &&
            throw(ArgumentError("`max_cells` must not larger than max number of cells ($num_cells) but got $max_cells."))
      mtx_max_counts = maximum(sum(p.data, dims=1))
      (max_counts > mtx_max_counts) &&
            throw(ArgumentError("max_counts must not larger than max count of matrix ($mtx_max_counts) but got $max_counts."))
      
      expressed_cells = sum(p.data .!= 0, dims=1)
      idx = min_cells .<= expressed_cells .<= max_cells
      expressed_counts = sum(p.data, dims=1)
      idx .&= min_counts .<= expressed_counts .<= max_counts
      idx = findall(vec(idx))
      p.data = p.data[idx, :]
      p.var = p.var[idx, :]
      filter_layers!(p, var_idx=idx)
      p
end

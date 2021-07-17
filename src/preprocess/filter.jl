function filter_cells!(p::Profile; min_counts::Real=0, min_genes::Real=0,
                       max_counts::Real=maximum(sum(p.data,dims=1)),
                       max_genes::Real=nrow(p))
      num_genes = nrow(p)
      mtx_max_counts = maximum(sum(p.data, dims=1))
      # TODO: can be refactored to more informative error message
      if min_counts < 0 || min_genes < 0 || max_genes > num_genes ||
            max_counts > mtx_max_counts
            throw(ArgumentError("""min_counts and min_genes must not lower than zero,
                                 max_counts must not larger than max count of matrix and
                                 max_genes must not larger than max number of genes.
                                 """))
      end
      expressed_genes = sum(p.data .!= 0, dims=1)
      idx = min_genes .<= expressed_genes .<= max_genes
      library_size = sum(p.data, dims=1)
      idx .&= min_counts .<= library_size .<= max_counts
      idx = findall(vec(idx))
      p.data = p.data[:, idx]
      p.obs = p.obs[idx, :]
      p
end

function filter_genes!(p::Profile; min_counts::Real=0, min_cells::Real=0,
                       max_counts::Real=maximum(sum(p.data,dims=2)),
                       max_cells::Real=ncol(p))
      num_cells = ncol(p)
      mtx_max_counts = maximum(sum(p.data, dims=2))
      if min_counts < 0 || min_cells < 0 || max_cells > num_cells ||
            max_counts > mtx_max_counts
            throw(ArgumentError("""min_counts and min_cells must not lower than zero,
                                 max_counts must not larger than max count of matrix and
                                 max_cells must not larger than max number of cells.
                                 """))
      end
      expressed_cells = sum(p.data .!= 0, dims=2)
      idx = min_cells .<= expressed_cells .<= max_cells
      expressed_counts = sum(p.data, dims=2)
      idx .&= min_counts .<= expressed_counts .<= max_counts
      idx = findall(vec(idx))
      p.data = p.data[idx, :]
      p.var = p.var[idx, :]
      p
end

"""
    quality_control_metrics(p; kwargs...)
    quality_control_metrics(X, obs, var; obsindex, varindex, kwargs...)

Calculate quality control metrics.

# Arguments

- `p::Profile`: The profile object to calculate on.
- `X::AbstractMatrix`: The count matrix to calculate on.
- `obs::DataFrame`: The feature information matrix with cell information.
- `var::DataFrame`: The feature information matrix with gene information.

# Common Keyword arguments

- `obsindex::Symbol=:barcode`: The index of dataframe `obs`.
- `varindex::Symbol=:gene_symbols`: The index of dataframe `var`.
- `use_log1p::Bool`: Computing log1p-transformed quality control metric.
- `qc_vars::AbstractVector{String}=["mt"]`: The boolean mask for identifying variables you could calculate on.
- `percent_top=nothing`: The proportions of top genes to calculate on.

# Example

See also [`quality_control_metrics!`](@ref) for non-inplace operation.
"""
function quality_control_metrics end

"""
    quality_control_metrics!(p; kwargs...)
    quality_control_metrics!(X, obs, var; kwargs...)

Calculate quality control metrics.

# Arguments

- `p::AnnotatedProfile`: The profile object to calculate on.
- `X::AbstractMatrix`: The count matrix to calculate on.
- `obs::DataFrame`: The feature information matrix with cell information.
- `var::DataFrame`: The feature information matrix with gene information.

# Common Keyword arguments

- `use_log1p::Bool`: Computing log1p-transformed quality control metric.
- `qc_vars::AbstractVector{String}=["mt"]`: The boolean mask for identifying variables you could calculate on.
- `percent_top=nothing`: The proportions of top genes to calculate on.

# Example

See also [`quality_control_metrics`](@ref) for non-inplace operation.
"""
function quality_control_metrics! end


quality_control_metrics(p::AnnotatedProfile; kwargs...) = quality_control_metrics!(copy(p); kwargs...)

function quality_control_metrics!(p::AnnotatedProfile; kwargs...)
    omic = p.omics[:RNA]
    obs_cols = obsnames(p)
    var_cols = names(omic.var)

    @info "Calculating quality control metrics:"
    obs_result = describe_obs!(p.RNA.count, p.obs, p.RNA.var; kwargs...)
    var_result = describe_var!(p.RNA.count, p.RNA.var; kwargs...)
    @info "  => generated qc metrics for $(setdiff(names(obs_result), obs_cols)) on observations"
    @info "  => generated qc metrics for $(setdiff(names(var_result), var_cols)) on features or variables"

    # omic.pipeline[:qc] = Dict(:qc_vars => qc_vars)
    @info "  => added :qc to pipeline in RNA"

    return p
end

function quality_control_metrics!(X::AbstractMatrix, obs::DataFrame, var::DataFrame; kwargs...)
    describe_obs!(X, obs, var; kwargs...)
    describe_var!(X, var; kwargs...)
    return obs, var
end

function quality_control_metrics(X::AbstractMatrix, obs::DataFrame, var::DataFrame;
    obsindex::Symbol=:barcode, varindex::Symbol=:gene_symbols, kwargs...)
    obs_metrics = DataFrame()
    var_metrics = DataFrame()
    obs_metrics[!, obsindex] = obs[!, obsindex]
    var_metrics[!, varindex] = var[!, varindex]
    describe_obs!(X, obs_metrics, var; kwargs...)
    describe_var!(X, var_metrics; kwargs...)
    return obs_metrics, var_metrics
end

function describe_obs!(X::AbstractMatrix, obs::DataFrame, var::DataFrame; use_log1p::Bool=false,
    qc_vars::AbstractVector{String}=["mt"], percent_top=nothing)

    obs[!, :n_genes_by_counts] = vec(count(X .!= 0, dims=1))
    obs[!, :total_counts] = vec(sum(X, dims=1))
    
    if use_log1p
        obs[!, :log1p_n_genes_by_counts] = log1p.(obs[!, :n_genes_by_counts])
        obs[!, :log1p_total_counts] = log1p.(obs[!, :total_counts])
    end

    # TODO: implement percent_top

    for qc_var in qc_vars
        colname = Symbol("total_counts_$qc_var")
        obs[!, colname] = vec(sum(X[var[!, qc_var], :], dims=1))
        if use_log1p
            obs[!, Symbol("log1p_total_counts_$qc_var")] = log1p.(obs[!, colname])
        end
        obs[!, Symbol("pct_total_counts_$qc_var")] = ( obs[!, colname] ./ obs[!, :total_counts] ) .* 100
    end

    # TODO: implement percent_top in qc_var

    return obs
end

function describe_var!(X::AbstractMatrix, var::DataFrame; use_log1p::Bool=false)
    var[!, :n_cells_by_counts] = vec(count(X .!= 0, dims=2))
    var[!, :mean_counts] = vec(mean(X, dims=2))
    var[!, :pct_dropout_by_counts] = (1 .- var[!, :n_cells_by_counts] ./ size(X, 1)) .* 100
    var[!, :total_counts] = vec(sum(X, dims=2))
    
    if use_log1p
        var[!, :log1p_mean_counts] = log1p.(var[!, :mean_counts])
        var[!, :log1p_total_counts] = log1p.(var[!, :total_counts])
    end

    return var
end

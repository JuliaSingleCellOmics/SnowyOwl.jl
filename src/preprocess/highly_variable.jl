"""
    highly_variable_genes(prof, method; omicsname=:RNA, layer=:count, kwargs...)
    highly_variable_genes(X, var, method; varname=:gene_symbols, kwargs...)

Calculate highly variable genes and return a new `DataFrame` with column :highly_variable,
which selects highly variable genes from given gene set. Additional information including
`:means`, `:dispersions` and `:dispersions_norm` will be added to returned `DataFrame`.

# Arguments

- `prof::AnnotatedProfile`: The profile object to calculate on.
- `X::AbstractMatrix`: The count matrix to calculate on.
- `var::DataFrame`: The feature information matrix with gene information.
- `method::HighlyVariableMethod`: Method to calculate highly variable genes, available for
    `CellRangerHVG`, `SeuratHVG` and `Seuratv3HVG`.

# Specific keyword arguments

- `omicsname::Symbol`: The `OmicsProfile` specified to calculate on.
- `layer::Symbol`: The layer specified to calculate on.
- `varname::Symbol`: The variable name to be specified as identifier for genes.

# Common keyword arguments

- `ntop_genes::Int=-1`: Number of top variable genes to be selected. Specify `-1` to switch
    to selection by mean and dispersion. Available for `CellRangerHVG` and `SeuratHVG` methods.
- `min_disp::Real=0.5`: Minimum dispersion for selecting highly variable genes. Available
    for `CellRangerHVG` and `SeuratHVG` methods.
- `max_disp::Real=Inf`: Maximum dispersion for selecting highly variable genes. Available
    for `CellRangerHVG` and `SeuratHVG` methods.
- `min_mean::Real=0.0125`: Minimum mean for selecting highly variable genes. Available for
    `CellRangerHVG` and `SeuratHVG` methods.
- `max_mean::Real=3.`: Maximum mean for selecting highly variable genes. Available for
    `CellRangerHVG` and `SeuratHVG` methods.

# Examples

See also [`highly_variable_genes!`](@ref) for inplace operation.
"""
function highly_variable_genes end

"""
    highly_variable_genes!(prof, method; omicsname=:RNA, layer=:count, kwargs...)
    highly_variable_genes!(X, var, method; kwargs...)

Calculate highly variable genes and modify `var` directly by adding column :highly_variable,
which selects highly variable genes from given gene set. Additional information including
`:means`, `:dispersions` and `:dispersions_norm` will be added to `var`.

# Arguments

- `prof::AnnotatedProfile`: The profile object to calculate on.
- `X::AbstractMatrix`: The count matrix to calculate on.
- `var::DataFrame`: The feature information matrix with gene information.
- `method`: Method to calculate highly variable genes, available for `CellRangerHVG`,
    `SeuratHVG` and `Seuratv3HVG`.

# Specific keyword arguments

- `omicsname::Symbol`: The `OmicsProfile` specified to calculate on.
- `layer::Symbol`: The layer specified to calculate on.

# Common keyword arguments

- `ntop_genes::Int=-1`: Number of top variable genes to be selected. Specify `-1` to switch
    to selection by mean and dispersion. Available for `CellRangerHVG` and `SeuratHVG` methods.
- `min_disp::Real=0.5`: Minimum dispersion for selecting highly variable genes. Available
    for `CellRangerHVG` and `SeuratHVG` methods.
- `max_disp::Real=Inf`: Maximum dispersion for selecting highly variable genes. Available
    for `CellRangerHVG` and `SeuratHVG` methods.
- `min_mean::Real=0.0125`: Minimum mean for selecting highly variable genes. Available for
    `CellRangerHVG` and `SeuratHVG` methods.
- `max_mean::Real=3.`: Maximum mean for selecting highly variable genes. Available for
    `CellRangerHVG` and `SeuratHVG` methods.

# Examples

See also [`highly_variable_genes`](@ref) for non-inplace operation.
"""
function highly_variable_genes! end


highly_variable_genes(p::AnnotatedProfile, method::HighlyVariableMethod=Seuratv3HVG(); kwargs...) =
    highly_variable_genes!(deepcopy(p), method; kwargs...)

function highly_variable_genes!(p::AnnotatedProfile, method::HighlyVariableMethod=Seuratv3HVG();
                                omicsname::Symbol=:RNA, layer::Symbol=:count, kwargs...)
    omic = p.omics[omicsname]
    @assert haskey(omic.layers, layer) "$layer not found in layers."
    X = getlayer(omic, layer)
    if method isa SeuratHVG && haskey(omic.pipeline, :log1p)
        X .*= log(omic.pipeline[:log1p][:base])
    end

    @info "Calculating highly variable genes using $method method over $omicsname.$layer:"
    cols = names(omic.var)
    hvg_result = highly_variable_genes!(X, omic.var, method; kwargs...)
    @info "  => generated statistics for $(setdiff(names(hvg_result), cols))"

    omic.pipeline[:hvg] = Dict(:method => Symbol(method), :layer => layer)
    @info "  => added :hvg to pipeline in $omicsname"

    return p
end

function highly_variable_genes(X::AbstractMatrix, var::DataFrame, method::HighlyVariableMethod;
                               varname::Symbol=:gene_symbols, kwargs...)
    df = DataFrame()
    df[!, varname] = var[!, varname]
    return highly_variable_genes!(copy(X), df, method; kwargs...)
end


function highly_variable_genes!(df::DataFrame, X::AbstractMatrix, ::Seuratv3HVG;
        ntop_genes::Integer=2000, span::Real=0.3, check_values::Bool=true)
    error("This method is not implemented.")
end

function highly_variable_genes!(X::AbstractMatrix, var::DataFrame, ::SeuratHVG;
    ntop_genes::Int=-1, min_disp=0.5, max_disp=Inf, min_mean=0.0125, max_mean=3.,
    nbins::Int=20)
    X = expm1.(X)
    μ, σ² = mean_and_var(X, 2, corrected=true)
    μ, σ² = vec(μ), vec(σ²)

    # compute dispersion
    replace!(μ, 0. => 1e-12)
    dispersion = replace!(σ² ./ μ, 0. => NaN)
    dispersion .= log.(dispersion)
    μ .= log1p.(μ)

    var[!, :means] = μ
    var[!, :dispersions] = dispersion
    h = fit(Histogram, μ; nbins=nbins)
    var[!, :mean_bin] = map(x -> StatsBase.binindex(h, x), μ)
    gdf = groupby(var, :mean_bin)
    disp_bin = combine(gdf, [:dispersions => mean, :dispersions => std])
    disp_bin = fill_missing_bins(disp_bin, :mean_bin, :dispersions_mean => 0,
                                 :dispersions_std => NaN)
    # retrieve those genes that have nan std, these are the ones where only a single gene
    # fell in the bin and implicitly set them to have a normalized disperion of 1
    one_gene_per_bin = isnan.(disp_bin.dispersions_std)
    gen_indices = findall(one_gene_per_bin[var.mean_bin])
    if length(gen_indices) > 0
        msg = "Gene indices $gen_indices fell into a single bin: their normalized dispersion" *
            "was set to 1.\nDecreasing `n_bins` will likely avoid this effect."
        @debug msg
    end
    disp_bin[one_gene_per_bin, :dispersions_std] .= disp_bin[one_gene_per_bin, :dispersions_mean]
    disp_bin[one_gene_per_bin, :dispersions_mean] .= 0
    disp_norm = normalize_by_bins(var.dispersions, disp_bin.dispersions_mean,
                          disp_bin.dispersions_std, var.mean_bin)
    var[!, :dispersions_norm] = disp_norm

    var[!, :highly_variable] = select_ntop_genes(μ, disp_norm, ntop_genes, nrow(var),
        min_mean, max_mean, min_disp, max_disp)
    select!(var, Not(:mean_bin))
    return var
end

function highly_variable_genes!(X::AbstractMatrix, var::DataFrame, ::CellRangerHVG;
    ntop_genes::Int=-1, min_disp=0.5, max_disp=Inf, min_mean=0.0125, max_mean=3.)
    μ, σ² = mean_and_var(X, 2, corrected=true)
    μ, σ² = vec(μ), vec(σ²)

    # compute dispersion
    replace!(μ, 0. => 1e-12)
    dispersion = σ² ./ μ

    var[!, :means] = μ
    var[!, :dispersions] = dispersion
    bins = [-Inf, percentile(μ, 10:5:100)..., Inf]
    h = fit(Histogram, μ, bins)
    var[!, :mean_bin] = map(x -> StatsBase.binindex(h, x), μ)
    gdf = groupby(var, :mean_bin)
    disp_bin = combine(gdf, [:dispersions => median, :dispersions => mad])
    disp_bin = fill_missing_bins(disp_bin, :mean_bin, :dispersions_median => 0,
                                 :dispersions_mad => NaN)

    disp_norm = normalize_by_bins(var.dispersions, disp_bin.dispersions_median,
                                        disp_bin.dispersions_mad, var.mean_bin)
    var[!, :dispersions_norm] = disp_norm

    var[!, :highly_variable] = select_ntop_genes(μ, disp_norm, ntop_genes, nrow(var),
        min_mean, max_mean, min_disp, max_disp)
    select!(var, Not(:mean_bin))
    return var
end

normalize_by_bins(x::AbstractVector, center, disp, bins) = @. (x - center[bins]) / disp[bins]

function select_ntop_genes(μ::AbstractVector, disp_norm::AbstractVector, ntop_genes::Int,
        total_genes::Int, min_mean::Real, max_mean::Real, min_disp::Real, max_disp::Real)
    if ntop_genes < 0
        return select_ntop_genes(μ, disp_norm, min_mean, max_mean, min_disp, max_disp)
    else
        if ntop_genes > total_genes
            @info "n_top_genes ($ntop_genes) is greater than number of gene, returning all genes."
            ntop_genes = total_genes
        end
        return select_ntop_genes(disp_norm, ntop_genes)
    end
end

function select_ntop_genes(μ::AbstractVector, disp_norm::AbstractVector, min_mean::Real,
        max_mean::Real, min_disp::Real, max_disp::Real)
    disp_norm = replace(disp_norm, NaN => 0)
    return (min_mean .< μ .< max_mean) .& (min_disp .< disp_norm .< max_disp)
end

function select_ntop_genes(disp_norm::AbstractVector, ntop_genes::Int)
    filtered_disp_norm = disp_norm[.!isnan.(disp_norm)]
    sort!(filtered_disp_norm, rev=true)
    if ntop_genes > length(filtered_disp_norm)
        @warn "`n_top_genes` > number of normalized dispersions, returning all genes with normalized dispersions."
        ntop_genes = length(filtered_disp_norm)
    end
    cutoff = filtered_disp_norm[ntop_genes]
    gene_subset = replace(disp_norm, NaN => 0) .>= cutoff
    @debug "$ntop_genes top genes correspond to a normalized dispersion cutoff of $cutoff."
    return gene_subset
end

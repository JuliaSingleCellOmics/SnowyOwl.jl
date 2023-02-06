"""
    project(prof, method; omicsname=:RNA, layer=:count, kwargs...)
    project(X, method; kwargs...)

Project the count matrix from original space to low-dimensional space.

# Arguments

- `prof::AnnotatedProfile`: The profile object to calculate on.
- `X::AbstractMatrix`: The count matrix to calculate on. The dimension of `X` should be
    (`nfeature`, `nsample`) where `nfeature` is number of features and `nsample` is number
    of samples.
- `method::AnalysisMethod`: Method to project count matrix, available for `PCAMethod()` and
    `UMAPMethod()`.

# Specific keyword arguments

- `omicsname::Symbol`: The `OmicsProfile` specified to calculate on.
- `layer::Symbol`: The layer specified to calculate on.

# Common keyword arguments

- `dims::Int=50`: The dimension of target space. Should lower or equal to the dimension of
    feature from count matrix which is `dims` <= `nfeature`.
- `n_neighbors::Int=20`: The number of local neighborhood points used in approximations of
    manifold structure. It is only effective when `method = UMAPMethod()`.
- `min_dist::Real=0.5`: This controls how tightly the embedding is allowed compress points
    together. Larger values prevent points from packing together and will prone to the
    preservation of global structure instead. Smaller values result in a more clustered/clumped
    embedding and prone to the preservation of local structure. It is only effective when
    `method = UMAPMethod()`.

# Examples

See also [`project!`](@ref) for inplace operation.
"""
function project end

"""
    project!(prof, method; omicsname=:RNA, layer=:count, kwargs...)
    project!(X, method; kwargs...)

Project the data matrix from original space to (usually) low-dimensional space.

# Arguments

- `prof::AnnotatedProfile`: The profile object to calculate on.
- `X::AbstractMatrix`: The count matrix to calculate on. The dimension of `X` should be
    (`nfeature`, `nsample`) where `nfeature` is number of features and `nsample` is number
    of samples.
- `method::AnalysisMethod`: Method to project count matrix, available for `PCAMethod()` and
    `UMAPMethod()`.

# Specific keyword arguments

- `omicsname::Symbol`: The `OmicsProfile` specified to calculate on.
- `layer::Symbol`: The layer specified to calculate on.

# Common keyword arguments

- `dims::Int=50`: The dimension of target space. Should lower or equal to the dimension of
    feature from count matrix which is `dims` <= `nfeature`.
- `n_neighbors::Int=20`: The number of local neighborhood points used in approximations of
    manifold structure. It is only effective when `method = UMAPMethod()`.
- `min_dist::Real=0.5`: This controls how tightly the embedding is allowed compress points
    together. Larger values prevent points from packing together and will prone to the
    preservation of global structure instead. Smaller values result in a more clustered/clumped
    embedding and prone to the preservation of local structure. It is only effective when
    `method = UMAPMethod()`.

# Examples

See also [`project`](@ref) for non-inplace operation.
"""
function project! end

project(p::AnnotatedProfile, method::AnalysisMethod; kwargs...) =
    project!(deepcopy(p), method; kwargs...)

function project(X::AbstractMatrix, method::AnalysisMethod; dims::Int=50, kwargs...)
    nfeature, nsample = size(X)
    @assert dims ≤ nfeature "the feature dimension of count matrix must be lower or equal to dims, but got $nfeature"
    return project!(similar(X, dims, nsample), X, method; dims=dims, kwargs...)
end

function project!(p::AnnotatedProfile, method::AnalysisMethod; dims::Int=50,
                  omicsname::Symbol=:RNA, layer::Symbol=:count, kwargs...)
    omic = p.omics[omicsname]
    @assert haskey(omic.layers, layer) "$layer not found in layers."
    X = getlayer(omic, layer)

    @info "Projecting $omicsname.$layer using $method method:"
    ncell = size(X, 2)
    Y = similar(X, dims, ncell)
    embed, model = project!(Y, X, method; dims=dims, return_model=true, kwargs...)
    @info "  => projected $layer"

    setpipeline!(omic, Symbol(method), Dict(pairs(kwargs)))
    saveembeddings!(p, method, omicsname, layer, embed)
    savemodel!(p, method, omicsname, layer, model)
    @info "  => added $method to pipeline in $omicsname"

    return p
end

function project!(Y::AbstractMatrix, X::AbstractMatrix, ::PCAMethod; dims::Int=50,
                  return_model::Bool=false)
    @assert size(Y) == (dims, size(X, 2))
    model = fit(PCA, X; maxoutdim=dims, pratio=1.)
    Y .= MultivariateStats.predict(model, X)
    return return_model ? (Y, model) : Y
end

function project!(Y::AbstractMatrix, X::AbstractMatrix, ::UMAPMethod; dims::Int=2,
                  n_neighbors::Int=20, min_dist::Real=0.5, return_model::Bool=false)
    @assert size(Y) == (dims, size(X, 2))
    Y .= UMAP.umap(X, dims, n_neighbors=n_neighbors, min_dist=min_dist)
    model = nothing  # TODO: fetch model from UMAP package
    return return_model ? (Y, model) : Y
end

function saveembeddings!(p::AnnotatedProfile, method::AnalysisMethod, omicsname::Symbol,
                         layer::Symbol, embed)
    project_name = Symbol(omicsname, "_", layer, "_", lowercase(repr(method)))
    p.obsm[project_name] = embed
    return p
end

function savemodel!(p::AnnotatedProfile, ::PCAMethod, omicsname::Symbol, layer::Symbol, model)
    omic = p.omics[omicsname]
    omic.varm[:PCs] = model.proj
    return p
end

function savemodel!(p::AnnotatedProfile, ::UMAPMethod, omicsname::Symbol, layer::Symbol, model)
    # TODO: decouple information from UMAP struct and save them
    return p
end

"""
    project(prof, method; omicsname=:RNA, layer=:count, kwargs...)
    project(X, method; kwargs...)

Project the count matrix from original space to low-dimensional space.

# Arguments

- `prof::AnnotatedProfile`: The profile object to calculate on.
- `X::AbstractMatrix`: The count matrix to calculate on. The dimension of `X` should be
    (`nfeature`, `nsample`) where `nfeature` is number of features and `nsample` is number
    of samples.
- `method::Symbol`: Method to project count matrix, available for `:pca` and `:umap`.

# Specific keyword arguments

- `omicsname::Symbol`: The `OmicsProfile` specified to calculate on.
- `layer::Symbol`: The layer specified to calculate on.

# Common keyword arguments

- `dims::Int=50`: The dimension of target space. Should lower or equal to the dimension of
    feature from count matrix which is `dims` <= `nfeature`.

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
- `method::Symbol`: Method to project count matrix, available for `:pca` and `:umap`.

# Specific keyword arguments

- `omicsname::Symbol`: The `OmicsProfile` specified to calculate on.
- `layer::Symbol`: The layer specified to calculate on.

# Common keyword arguments

- `dims::Int=50`: The dimension of target space. Should lower or equal to the dimension of
    feature from count matrix which is `dims` <= `nfeature`.

# Examples

See also [`project`](@ref) for non-inplace operation.
"""
function project! end

project(p::AnnotatedProfile, method::Symbol; kwargs...) =
    project!(deepcopy(p), method; kwargs...)

function project(X::AbstractMatrix, method::Symbol; dims::Int=50, kwargs...)
    nfeature, nsample = size(X)
    @assert dims ≤ nfeature "the feature dimension of count matrix must be lower or equal to dims, but got $nfeature"
    return project!(similar(X, dims, nsample), X, method; dims=dims, kwargs...)
end

function project!(p::AnnotatedProfile, method::Symbol; omicsname::Symbol=:RNA, layer::Symbol=:count, kwargs...)
    omic = p.omics[omicsname]
    @assert haskey(omic.layers, layer) "$layer not found in layers."
    X = getlayer(omic, layer)

    @info "Projecting $omicsname.$layer using $method method:"
    project!(X, method; kwargs...)
    @info "  => projected $layer"

    setpipeline!(omic, method, Dict(pairs(kwargs)))
    @info "  => added $method to pipeline in $omicsname"

    return p
end

project!(Y::AbstractMatrix, X::AbstractMatrix, method::Symbol; kwargs...) =
    project!(Y, X, Val(method); kwargs...)

function project!(Y::AbstractMatrix, X::AbstractMatrix, ::Val{:pca}; dims::Int=50)
    @assert size(Y) == (dims, size(X, 2))
    model = fit(PCA, X; maxoutdim=dims, pratio=1.)
    Y .= MultivariateStats.transform(model, X)
    return Y
end

function project!(Y::AbstractMatrix, X::AbstractMatrix, ::Val{:umap}; dims::Int=2, n_neighbors=20, min_dist=0.5)
    @assert size(Y) == (dims, size(X, 2))
    Y .= UMAP.umap(X, dims, n_neighbors=n_neighbors, min_dist=min_dist)
    return Y
end

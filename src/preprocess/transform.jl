"""
    normalize(prof, method; omicsname=:RNA, layer=:count, kwargs...)
    normalize(X, method; kwargs...)

Normalize counts per cell.

# Arguments

- `prof::AnnotatedProfile`: The profile object to calculate on.
- `X::AbstractMatrix`: The count matrix to calculate on.
- `method`: Method to calculate highly variable genes, available for `LogNormalization`,
    `RelativeNormalization`, `CenteredLogRatioNormalization` and `CustomNormalization`.

# Specific keyword arguments

- `omicsname::Symbol`: The `OmicsProfile` specified to calculate on.
- `layer::Symbol`: The layer specified to calculate on.

# Common keyword arguments

- `scaled_size::Real=1e4`: Scaled total counts for each cell over all genes.

# Examples

See also [`normalize!`](@ref) for inplace operation.
"""
function normalize end

"""
    normalize!(prof, method; omicsname=:RNA, layer=:count, kwargs...)
    normalize!(X, method; kwargs...)

Normalize counts per cell.

# Arguments

- `prof::AnnotatedProfile`: The profile object to calculate on.
- `X::AbstractMatrix`: The count matrix to calculate on.
- `method`: Method to calculate highly variable genes, available for `LogNormalization`,
    `RelativeNormalization`, `CenteredLogRatioNormalization` and `CustomNormalization`.

# Specific keyword arguments

- `omicsname::Symbol`: The `OmicsProfile` specified to calculate on.
- `layer::Symbol`: The layer specified to calculate on.

# Common keyword arguments

- `scaled_size::Real=1e4`: Scaled total counts for each cell over all genes.

# Examples

See also [`normalize`](@ref) for non-inplace operation.
"""
function normalize! end

normalize(p::AnnotatedProfile, method::NormalizationMethod=LogNormalization(); kwargs...) =
    normalize!(deepcopy(p), method; kwargs...)

normalize(X::AbstractMatrix, method::NormalizationMethod; kwargs...) =
    normalize!(copy(X), method; kwargs...)

function normalize!(p::AnnotatedProfile, method::NormalizationMethod=LogNormalization();
                    omicsname::Symbol=:RNA, layer::Symbol=:count, kwargs...)
    omic = p.omics[omicsname]
    @assert haskey(omic.layers, layer) "$layer not found in layers."
    X = getlayer(omic, layer)

    @info "Normalizing $omicsname.$layer using $method method:"
    normalize!(X, method; kwargs...)
    @info "  => normalized $omicsname.$layer"

    setpipeline!(omic, :normalize, Dict(:method => Symbol(method), :layer => layer))
    @info "  => added :normalize to pipeline in $omicsname"

    return p
end

normalize!(X::AbstractMatrix, method::NormalizationMethod; kwargs...) =
    normalize!(X, method; kwargs...)

function normalize!(X::AbstractMatrix, ::LogNormalization; scaled_size::Real=1e4)
    factors = scaled_size ./ sum(X, dims=1)
    X .= @. log1p(X * factors)
    return X
end

function normalize!(X::AbstractMatrix, ::RelativeNormalization; scaled_size::Real=1e4)
    factors = scaled_size ./ sum(X, dims=1)
    X .*= factors
    return X
end

function normalize!(X::AbstractMatrix, ::CenteredLogRatioNormalization)
    factor = sum(log1p.(X[X .> 0])) / length(X)
    X .= @. log1p(X / exp(factor))
    return X
end

function normalize!(X::AbstractMatrix, ::CustomNormalization; scaled_size::Real=1e4,
                    dims::Int=1)
    error("This method is not implemented.")
end

"""
    logarithmize(prof; omicsname=:RNA, layer=:count, kwargs...)
    logarithmize(X; kwargs...)

Logarithmize the data matrix.

# Arguments

- `prof::AnnotatedProfile`: The profile object to calculate on.
- `X::AbstractMatrix`: The count matrix to calculate on.

# Specific keyword arguments

- `omicsname::Symbol`: The `OmicsProfile` specified to calculate on.
- `layer::Symbol`: The layer specified to calculate on.

# Common keyword arguments

- `base::Real=ℯ`: Base of the logarithm.

# Examples

See also [`logarithmize!`](@ref) for inplace operation.
"""
function logarithmize end

"""
    logarithmize!(prof; omicsname=:RNA, layer=:count, kwargs...)
    logarithmize!(X; kwargs...)

Logarithmize the data matrix.

# Arguments

- `prof::AnnotatedProfile`: The profile object to calculate on.
- `X::AbstractMatrix`: The count matrix to calculate on.

# Specific keyword arguments

- `omicsname::Symbol`: The `OmicsProfile` specified to calculate on.
- `layer::Symbol`: The layer specified to calculate on.

# Common keyword arguments

- `base::Real=ℯ`: Base of the logarithm.

# Examples

See also [`logarithmize`](@ref) for non-inplace operation.
"""
function logarithmize! end

logarithmize(p::AnnotatedProfile; kwargs...) = logarithmize!(deepcopy(p); kwargs...)

logarithmize(X::AbstractMatrix; kwargs...) = logarithmize!(copy(X); kwargs...)

function logarithmize!(p::AnnotatedProfile; omicsname::Symbol=:RNA, layer::Symbol=:count,
                       kwargs...)
    omic = p.omics[omicsname]
    @assert haskey(omic.layers, layer) "$layer not found in layers."
    X = getlayer(omic, layer)

    @info "Logarithmizing $omicsname.$layer:"
    logarithmize!(X; kwargs...)
    @info "  => calculated log1p for $layer"

    base = haskey(kwargs, :base) ? kwargs[:base] : ℯ
    setpipeline!(omic, :log1p, Dict{Symbol,Any}(:base => base))
    @info "  => added :log1p to pipeline in $omicsname"

    return p
end

function logarithmize!(X::AbstractMatrix; base::Real=ℯ)
    X .= @. log1p(X) / log(base)
    return X
end

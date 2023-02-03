"""
    normalize(prof, method; omicsname=:RNA, layer=:count, kwargs...)
    normalize(X, method; kwargs...)

Normalize counts per cell.

# Arguments

- `prof::AnnotatedProfile`: The profile object to calculate on.
- `X::AbstractMatrix`: The count matrix to calculate on.
- `method`: Method to calculate highly variable genes, available for `:lognormalize`,
    `:relative`, `:clr` and `:custom`.

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
- `method`: Method to calculate highly variable genes, available for `:lognormalize`,
    `:relative`, `:clr` and `:custom`.

# Specific keyword arguments

- `omicsname::Symbol`: The `OmicsProfile` specified to calculate on.
- `layer::Symbol`: The layer specified to calculate on.

# Common keyword arguments

- `scaled_size::Real=1e4`: Scaled total counts for each cell over all genes.

# Examples

See also [`normalize`](@ref) for non-inplace operation.
"""
function normalize! end

normalize(p::AnnotatedProfile, method::Symbol=:lognormalize; kwargs...) =
    normalize!(deepcopy(p), method; kwargs...)

normalize(X::AbstractMatrix, method::Symbol; kwargs...) =
    normalize(copy(X), Val(method); kwargs...)

function normalize!(p::AnnotatedProfile, method::Symbol=:lognormalize;
                    omicsname::Symbol=:RNA, layer::Symbol=:count, kwargs...)
    omic = p.omics[omicsname]
    @assert haskey(omic.layers, layer) "$layer not found in layers."
    X = getlayer(omic, layer)

    @info "Normalizing $omicsname.$layer using $method method:"
    normalize!(X, Val(method); kwargs...)
    @info "  => normalized $omicsname.$layer"

    omic.pipeline[:normalize] = Dict(:method => method, :layer => layer)
    @info "  => added :normalize to pipeline in $omicsname"

    return p
end

normalize!(X::AbstractMatrix, method::Symbol; kwargs...) =
    normalize!(X, Val(method); kwargs...)

function normalize!(X::AbstractMatrix, method::Val{:lognormalize}; scaled_size::Real=1e4)
    colsums = sum(X, dims=1)
    X .= @. log1p(X) / colsums * scaled_size
    return X
end

"""
Normalize count data to relative counts per cell by dividing by the total counts per cell.
"""
function normalize!(X::AbstractMatrix, method::Val{:relative}; scaled_size::Real=1e4)
    X .= @. X / sum(X, dims=1) * scaled_size
    return X
end

# function normalize!(X::AbstractMatrix, method::Val{:clr}; scaled_size::Real=1e4, dims::Int=1)

# end

# function normalize!(X::AbstractMatrix, method::Val{:custom}; scaled_size::Real=1e4, dims::Int=1)

# end

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

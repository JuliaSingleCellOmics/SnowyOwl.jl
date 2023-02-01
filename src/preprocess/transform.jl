# normalize_total(p::AnnotatedProfile, library_size::Real=1e6)

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

"""
    scale(p; kwargs...)
    scale(X; kwargs...)

Scale each gene to zero mean and unit variance.

# Arguments

- `p::AnnotatedProfile`: The profile object to calculate on.
- `X::AbstractMatrix`: The count matrix to calculate on.

# Common Keyword arguments

- `max_value::Float64`: Clip to this value after scaling.

# Example

See also [`scale!`](@ref) for inplace operation.
"""
function scale end

"""
    scale!(p; kwargs...)
    scale!(X; kwargs...)

Scale each gene to zero mean and unit variance.

# Arguments

- `p::AnnotatedProfile`: The profile object to calculate on.
- `X::AbstractMatrix`: The count matrix to calculate on.

# Common Keyword arguments

- `max_value::Float64`: Clip to this value after scaling.

# Example

See also [`scale`](@ref) for non-inplace operation.
"""
function scale! end


scale(p::AnnotatedProfile; kwargs...) = scale!(copy(p); kwargs...)

function scale!(p::AnnotatedProfile; max_value=nothing)
    @info "Scaling each gene:"
    X = p.RNA.count
    X = mapcols!(zscore, X)

    if max_value !== nothing
        @info "Clipping at max_value $(max_value)"
        X[X.>=max_value] .= max_value

    return p
end

function scale!(X::AbstractMatrix; max_value=nothing)
    @info "Scaling each gene:"
    X = mapcols!(zscore, X)

    if max_value !== nothing
        @info "Clipping at max_value $(max_value)"
        X[X.>=max_value] .= max_value
    
    return X
end

function scale(X::AbstractMatrix; max_value=nothing)
    @info "Scaling each gene:"
    X = mapcols(zscore, X)

    if max_value !== nothing
        @info "Clipping at max_value $(max_value)"
        X[X.>=max_value] .= max_value

    return X
end

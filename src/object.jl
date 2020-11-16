mutable struct Profile{T<:AbstractMatrix,S<:AbstractMatrix}
    data::T
    obs::DataFrame
    var::DataFrame
    layers::Dict{Symbol,S}
    pipeline::OrderedDict

    function Profile(data::T, obs::DataFrame, var::DataFrame) where {T<:AbstractMatrix}
        r, c = size(data)
        @assert nrow(obs) == c
        @assert nrow(var) == r
        new{T,S}(data, obs, var, Dict{Symbol,Matrix}(), OrderedDict{Symbol,Dict}())
    end
end

obsnames(p::Profile) = names(p.obs)
varnames(p::Profile) = names(p.var)
layernames(p::Profile) = keys(p.layers)
nrow(p::Profile) = size(p.data, 1)
ncol(p::Profile) = size(p.data, 2)

function Base.show(io::IO, p::Profile)
    println(io, "Profile(n_var × n_obs = ", nrow(p), " × ", ncol(p), ")")
    isempty(p.obs) || println(io, "    obs: ", join(obsnames(p), ", "))
    isempty(p.var) || println(io, "    var: ", join(varnames(p), ", "))
    isempty(p.layers) || println(io, "    layers: ", join(layernames(p), ", "))
    isempty(p.pipeline) || println(io, "    pipeline: ", join(keys(p.pipeline), ", "))
end

Base.maximum(p::Profile) = maximum(p.data)
Base.minimum(p::Profile) = minimum(p.data)

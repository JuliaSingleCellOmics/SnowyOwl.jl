mutable struct Profile{T<:AbstractMatrix,S<:AbstractMatrix}
    data::T
    obs::DataFrame
    var::DataFrame
    layers::Dict{Symbol,S}
    pipeline::OrderedDict
end

function Profile(data::T, obs::DataFrame, var::DataFrame) where {T<:AbstractMatrix}
    r, c = size(data)
    @assert nrow(obs) == c
    @assert nrow(var) == r
    Profile{T,Matrix}(data, obs, var, Dict{Symbol,Matrix}(), OrderedDict{Symbol,Dict}())
end

obsnames(p::Profile) = names(p.obs)
varnames(p::Profile) = names(p.var)
layernames(p::Profile) = keys(p.layers)
nrow(p::Profile) = size(p.data, 1)
ncol(p::Profile) = size(p.data, 2)
nvar(p::Profile) = size(p.data, 1)
nobs(p::Profile) = size(p.data, 2)

function Base.show(io::IO, p::Profile)
    println(io, "Profile(n_var × n_obs = ", nrow(p), " × ", ncol(p), ")")
    isempty(p.obs) || println(io, "    obs: ", join(obsnames(p), ", "))
    isempty(p.var) || println(io, "    var: ", join(varnames(p), ", "))
    isempty(p.layers) || println(io, "    layers: ", join(layernames(p), ", "))
    isempty(p.pipeline) || println(io, "    pipeline: ", join(keys(p.pipeline), ", "))
end

Base.maximum(p::Profile) = maximum(p.data)
Base.minimum(p::Profile) = minimum(p.data)

Base.size(p::Profile) = size(p.data)
Base.axes(p::Profile) = axes(p.data)

function Base.getindex(p::Profile, inds...)
    p_ = Profile(getindex(p.data, inds...),
                 getindex(p.obs, inds[2], :),
                 getindex(p.var, inds[1], :))
    for (k, v) in p.layers
        p_.layers[k] = getindex(v, inds...)
    end
    p_.pipeline = copy(p.pipeline)
    p_
end

Base.setproperty!(prof::Profile, name::Symbol, x) = setproperty!(prof, Val(name), x)
Base.setproperty!(prof::Profile, ::Val{S}, x) where S = setfield!(prof, S, x)

function Base.setproperty!(prof::Profile, ::Val{:obs}, x)
    @assert nrow(x) == size(prof.data, 2)
    setfield!(prof, :obs, x)
end

function Base.setproperty!(prof::Profile, ::Val{:var}, x)
    @assert nrow(x) == size(prof.data, 1)
    setfield!(prof, :var, x)
end

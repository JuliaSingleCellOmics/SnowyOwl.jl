mutable struct Profile{T<:AbstractMatrix,S}
    data::T
    var::DataFrame
    obs::DataFrame
    layers::Dict{Symbol,S}
    pipeline::OrderedDict
end

function Profile(data::M, var::DataFrame, obs::DataFrame; T=float(eltype(data))) where {M<:AbstractMatrix}
    @assert (nrow(var), nrow(obs)) == size(data)
    S = Union{Matrix{T},SparseMatrixCSC{T,UInt32}}
    data = T.(data)
    layers = Dict{Symbol,S}()
    pipeline = OrderedDict{Symbol,Dict}()
    Profile{typeof(data),S}(data, copy(var), copy(obs), layers, pipeline)
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

Base.copy(p::Profile) = Profile(copy(p.data), copy(p.var), copy(p.obs), copy(p.layers), copy(p.pipeline))

Base.maximum(p::Profile) = maximum(p.data)
Base.minimum(p::Profile) = minimum(p.data)

Base.size(p::Profile) = size(p.data)
Base.axes(p::Profile) = axes(p.data)

function Base.getindex(p::Profile, inds...)
    p_ = Profile(getindex(p.data, inds[1], inds[2]),
                 getindex(p.var, inds[1], :),
                 getindex(p.obs, inds[2], :))
    for (k, v) in p.layers
        if size(v) == size(p_.data)
            p_.layers[k] = getindex(v, inds[1], inds[2])
        end
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

Base.filter(x::Pair{Symbol,T}, prof::Profile) where {T} = filter!(x, copy(prof))

function Base.filter!(x::Pair{Symbol,T}, prof::Profile) where {T}
    col, f = x
    sel = f.(prof.var[:,col])
    filter!(x, prof.var)
    prof.data = prof.data[sel, :]
    filter_layers!(prof, var_idx=sel)
    prof
end

function filter_layers!(prof::Profile; var_idx=(:), obs_idx=(:))
    if isempty(prof.layers)
        for k in keys(prof.layers)
            if size(prof.layers[k]) == size(prof.data)
                prof.layers[k] = prof.layers[k][var_idx, obs_idx]
            end
        end
    end
    prof
end
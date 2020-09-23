struct Profile{T<:AbstractMatrix}
    data::T
    obs::DataFrame
    var::DataFrame

    function Profile(data, obs, var)
        r, c = size(data)
        @assert nrow(obs) == c
        @assert nrow(var) == r
        new(data, obs, var)
    end
end

# function Base.show(::Profile)
    
# end
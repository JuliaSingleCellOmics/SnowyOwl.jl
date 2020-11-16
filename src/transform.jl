normalize_total(p::Profile, library_size::Real=1e6)

"""
Logarithmize the data matrix.
"""
function log1p(p::Profile; base::Real=ℯ)
    log1p.(p.data) ./ log(base)
end

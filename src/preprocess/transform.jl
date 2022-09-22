# normalize_total(p::AnnotatedProfile, library_size::Real=1e6)

"""
Logarithmize the data matrix.
"""
function log1p!(p::AnnotatedProfile; base::Real=ℯ)
    p.data .= @. log1p(p.data) / log(base)
    p.pipeline[:log1p] = Dict{Symbol,Any}(:base => base)
    return p
end

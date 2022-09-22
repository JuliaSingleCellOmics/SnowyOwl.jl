module Neighbors

using SparseArrays
using LinearAlgebra

export
    #moments
    normalize_indicator,
    to_indicator_matrix,
    union_diagonal!,
    graph_filter,
    moment,
    first_moment

# neighborhood_graph!()

include("moments.jl")

end

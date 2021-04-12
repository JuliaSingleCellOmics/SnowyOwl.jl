module Neighbors

using SparseArrays
using LinearAlgebra

using Distances
using Graphs
using NearestNeighborDescent
using StatsBase

export
    #moments
    normalize_indicator,
    to_indicator_matrix,
    union_diagonal!,
    graph_filter,
    moment,
    first_moment,

    # neighbors
    neighborhood_graph


include("utils.jl")
include("moments.jl")
include("neighbors.jl")

end

module SnowyOwl

using Reexport

include("datasets.jl")
include("preprocess/Preprocess.jl")
include("analysis/Analysis.jl")
include("neighbors/Neighbors.jl")
include("plots/Plots.jl")

function __init__()
    Dataset.__init__pbmc3k()
end

end

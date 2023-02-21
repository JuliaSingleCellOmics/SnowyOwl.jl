module Preprocess

using DataFrames
using OmicsProfiles
using StatsBase

export
    # filter
    filter_cells!,
    filter_cells,
    filter_genes!,
    filter_genes,

    # qc
    quality_control_metrics!,
    quality_control_metrics,

    # transform
    NormalizationMethod,
    LogNormalization,
    RelativeNormalization,
    CenteredLogRatioNormalization,
    CustomNormalization,
    normalize,
    normalize!,
    logarithmize,
    logarithmize!,

    # highly_variable
    highly_variable_genes,
    highly_variable_genes!

include("methodtype.jl")
include("utils.jl")
include("filter.jl")
include("qc.jl")
include("transform.jl")
include("highly_variable.jl")

end

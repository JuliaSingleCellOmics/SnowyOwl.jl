module SnowyOwl

using SparseArrays
using DataFrames

export
    # io
    read_mtx,
    read_features,
    read_barcodes,
    read_genes,
    read_cells,

    # object
    Profile,
    obsnames,
    varnames,
    layernames,
    nrow,
    ncol,
    nvar,
    nobs,

    # datasets
    load_pbmc68k,

    # filter
    filter_cells!,
    filter_genes!,

    # model
    unspliced,
    spliced,
    mRNA,

    #moments
    normalize_indicator,
    to_indicator_matrix,
    union_diagonal!,
    graph_filter,
    moment,
    first_moment

include("io.jl")
include("object.jl")
include("datasets.jl")
include("preprocess/filter.jl")
include("model.jl")
include("moments.jl")

end

module SnowyOwl

using DataStructures: OrderedDict
using SparseArrays
using PyCall
using CSV, DataFrames, CodecZlib, Mmap
using JLD2

import DataFrames: nrow, ncol

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
    filter_genes!

include("io.jl")
include("object.jl")
include("datasets.jl")
include("filter.jl")

end

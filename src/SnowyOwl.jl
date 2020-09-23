module SnowyOwl

using SparseArrays
using PyCall
using CSV, DataFrames, CodecZlib, Mmap

export
    read_mtx,
    read_features,
    read_barcodes,
    read_genes,
    read_cells

include("io.jl")

end

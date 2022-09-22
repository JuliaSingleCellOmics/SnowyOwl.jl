module SnowyOwl

using LinearAlgebra
using SparseArrays

using CSV, DataFrames
using CodecZlib, Mmap
using JLD2
using DataStructures: OrderedDict
using OmicsProfiles
using StatsBase

import DataFrames: nrow, ncol

const DEFAULT_FEATURE_COLS = [:ensembleid, :genesymbol, :type]
const DEFAULT_BARCODE_COLS = [:barcode]
const FEATURE_COLS = [:featurekey, :featurename, :featuretype, :chromosome, :featurestart, :featureend, :isgene, :genus_species]

export
    # io
    read_mtx,
    read_features,
    read_barcodes,
    read_genes,
    read_cells,

    # datasets
    load_pbmc68k,

    # highly_variable
    highly_variable_genes,
    highly_variable_genes!,

    #moments
    normalize_indicator,
    to_indicator_matrix,
    union_diagonal!,
    graph_filter,
    moment,
    first_moment

include("io.jl")
include("datasets.jl")
# include("preprocess/filter.jl")
include("highly_variable.jl")
include("moments.jl")
include("utils.jl")

end

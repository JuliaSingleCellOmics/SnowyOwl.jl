module SnowyOwl

using Reexport

@reexport using OmicsProfiles

import DataFrames: nrow, ncol

const DEFAULT_FEATURE_COLS = [:ensembleid, :genesymbol, :type]
const DEFAULT_BARCODE_COLS = [:barcode]
const FEATURE_COLS = [:featurekey, :featurename, :featuretype, :chromosome, :featurestart, :featureend, :isgene, :genus_species]


include("rw.jl")
include("datasets.jl")
include("preprocess/Preprocess.jl")
include("neighbors/Neighbors.jl")
include("plots.jl")

end

using SnowyOwl
using OmicsProfiles
using LinearAlgebra
using SparseArrays
using DataFrames
using Distributions
using StatsBase
using Test
using OmicsProfiles

tests = [
    # "filter",
    "qc",
    "highly_variable",
    "moments",
]

@testset "SnowyOwl.jl" begin
    for t in tests
        include("$(t).jl")
    end
end

using SnowyOwl
using OmicsProfiles
using LinearAlgebra
using SparseArrays
using DataFrames
using Distributions
using StatsBase
using Test
using OmicsProfiles
using Plots

const TEST_PATH = @__DIR__

tests = [
    # "filter",
    "qc",
    "transform",
    "highly_variable",
    "moments",
    "plots",
]

@testset "SnowyOwl.jl" begin
    for t in tests
        include("$(t).jl")
    end
end

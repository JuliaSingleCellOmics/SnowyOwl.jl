using SnowyOwl
using LinearAlgebra
using SparseArrays
using DataFrames
using Distributions
using Test

tests = [
    # "filter",
    "highly_variable",
    "moments",
]

@testset "SnowyOwl.jl" begin
    for t in tests
        include("$(t).jl")
    end
end

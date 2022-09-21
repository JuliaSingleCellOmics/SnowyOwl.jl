using SnowyOwl
using DataFrames
using LinearAlgebra
using SparseArrays
using Test

tests = [
    # "filter",
    "moments",
]

@testset "SnowyOwl.jl" begin
    for t in tests
        include("$(t).jl")
    end
end

using SnowyOwl
using DataFrames
using LinearAlgebra
using SparseArrays
using Test

tests = [
    "object",
    # "io"
    "filter",
    "moments",
    "qc",
]

@testset "SnowyOwl.jl" begin
    for t in tests
        include("$(t).jl")
    end
end

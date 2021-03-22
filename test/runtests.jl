using SnowyOwl
using DataFrames
using Test

tests = [
    "object",
    # "io"
    "moments",
]

@testset "SnowyOwl.jl" begin
    for t in tests
        include("$(t).jl")
    end
end

using SnowyOwl
using Test

tests = [
    "io"
]

@testset "SnowyOwl.jl" begin
    for t in tests
        include("$(t).jl")
    end
end

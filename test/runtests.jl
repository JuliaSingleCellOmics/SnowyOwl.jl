using SnowyOwl
using OmicsProfiles
using LinearAlgebra
using SparseArrays
using DataFrames
using Distributions
using StatsBase
using CUDA
using Plots
using Test

const TEST_PATH = @__DIR__

ENV["DATADEPS_ALWAYS_ACCEPT"] = true

cuda_tests = [
    "cuda",
]

tests = [
    "utils",
    "datasets",
    # "filter",
    "qc",
    "transform",
    "highly_variable",
    "moments",
    "project",
    "plots",
]

if CUDA.functional()
    CUDA.allowscalar(false)
    append!(tests, cuda_tests)
else
    @warn "CUDA unavailable, not testing GPU support"
end

@testset "SnowyOwl.jl" begin
    for t in tests
        include("$(t).jl")
    end
end

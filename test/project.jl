@testset "project" begin
    nfeature = 10
    nsample = 1000

    X = rand(nfeature, nsample)
    Y = SnowyOwl.Analysis.project(X, :pca; dims=5)
    @test size(Y) == (5, nsample)

    Y = SnowyOwl.Analysis.project(X, :umap; dims=2)
    @test size(Y) == (2, nsample)
end

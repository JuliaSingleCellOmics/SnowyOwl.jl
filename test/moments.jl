@testset "moments" begin
    ind = rand([0,1], 10)
    @test sum(SnowyOwl.Neighbors.normalize_indicator(ind)) ≈ 1.

    Ind = diagm(ind)
    @test sum(SnowyOwl.Neighbors.normalize_indicator(Ind)) ≈ 1.

    @test SnowyOwl.Neighbors.to_indicator_matrix(ind) == Ind

    Ind = rand([0,1], 10, 10)
    @test sum(SnowyOwl.Neighbors.normalize_indicator(Ind, 1), dims=1) ≈ ones(Int64, 1, 10)
    @test sum(SnowyOwl.Neighbors.normalize_indicator(Ind, 2), dims=2) ≈ ones(Int64, 10)

    @test diag(SnowyOwl.Neighbors.union_diagonal!(Ind)) == ones(Int64, 10)

    ind = rand([0,1], 10)

    G = SnowyOwl.Neighbors.graph_filter(ind, 1)
    @test G isa SparseVector
    @test sum(G .!= 0) == sum(ind .!= 0)
    @test sum(G) ≈ 1.

    G = SnowyOwl.Neighbors.graph_filter(ind, 2)
    @test G isa SparseMatrixCSC
    @test sum(G .!= 0) == sum(ind .!= 0)
    @test sum(G) ≈ 1.

    adj = rand([0,1], 10, 10)

    G = SnowyOwl.Neighbors.graph_filter(adj, 1, dims=1)
    @test G isa SparseMatrixCSC
    @test sum(G .!= 0, dims=1) == sum(adj .!= 0, dims=1)
    @test sum(G, dims=1) ≈ ones(Int64, 1, 10)

    G = SnowyOwl.Neighbors.graph_filter(adj, 1, dims=2)
    @test G isa SparseMatrixCSC
    @test sum(G .!= 0, dims=2) == sum(adj .!= 0, dims=2)
    @test sum(G, dims=2) ≈ ones(Int64, 10)

    adj = [1 0 0; 1 0 1; 1 1 1]
    blocks = [diagm([1/3, 1/3, 1/3]), diagm([0, 0, 1]), diagm([0, 1/2, 1/2])]
    @test SnowyOwl.Neighbors.graph_filter(adj, 2, dims=1) == cat(blocks...,dims=(1,2))

    blocks = [diagm([1, 0, 0]), diagm([1/2, 0, 1/2]), diagm([1/3, 1/3, 1/3])]
    @test SnowyOwl.Neighbors.graph_filter(adj, 2, dims=2) == cat(blocks...,dims=(1,2))

    A = rand([0,1], 5, 5)
    X = rand(5, 10)
    @test size(SnowyOwl.Neighbors.moment(A, X, 1, dims=1)) == size(X)
    @test size(SnowyOwl.Neighbors.moment(A, Matrix(X'), 1, dims=2)) == size(X')
    # @test size(SnowyOwl.Neighbors.moment(A, X, 2, dims=1)) == (size(X, 1), 5*10)
    # @test size(SnowyOwl.Neighbors.moment(A, Matrix(X'), 2, dims=2)) == (5*10, size(X, 1))
end

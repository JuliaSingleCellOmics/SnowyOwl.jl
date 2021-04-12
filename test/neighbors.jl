@testset "neighbors" begin
    @testset "knn_graph" begin
        X = [5 1 0 3 9 25;
             2 9 8 4 20 1]
        n_neighbors = 3

        D = pairwise(Euclidean(), X, dims=2)
        knn_indices, knn_dists = knn_graph(D, Val(:precomputed), n_neighbors)
        @test knn_indices == [4 3 2 1 2 1;
                              3 4 4 3 3 4;
                              2 1 1 2 4 5]
        @test knn_dists == sort(D, dims=1)[2:(n_neighbors+1), :]

    end

    # @testset "connected_components" begin

    # end
end
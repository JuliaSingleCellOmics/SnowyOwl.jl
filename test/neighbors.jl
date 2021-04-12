@testset "neighbors" begin
    @testset "neighborhood_graph" begin
        X = [5 2 7 7 1 8;
             1 9 3 4 3 6;
             0 8 6 Inf 11 2]
        n_neighbors = 3
        
        knn_indices, knn_dists = neighborhood_graph(X, Val(:umap), Val(:precomputed), n_neighbors)
        @test knn_indices == [5 2 1; 1 3 5; 1 6 3]
        @test knn_dists == [1 2 5;
                            1 3 3;
                            0 2 6]
        
    end
    
    # @testset "connected_components" begin

    # end
end
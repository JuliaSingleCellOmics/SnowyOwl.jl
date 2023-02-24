@testset "utils" begin
    X = [5 2 7 7;
         1 9 3 4;
         0 8 6 10]

    @test SnowyOwl.sortperm(X, dims=1) == [3 1 2 2;
                                          2 3 3 1;
                                          1 2 1 3]
    @test SnowyOwl.sortperm(X, dims=2) == [2 1 3 4;
                                          1 3 4 2;
                                          1 3 2 4]

    @test_throws ArgumentError SnowyOwl.sortperm(X, dims=3)
end

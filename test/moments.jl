@testset "moments" begin
    ind = rand([0,1], 10)
    @test sum(normalize_indicator(ind)) ≈ 1.

    Ind = diagm(ind)
    @test sum(normalize_indicator(Ind)) ≈ 1.

    @test to_indicator_matrix(ind) == Ind

    Ind = rand([0,1], 10, 10)
    @test sum(normalize_indicator(Ind, 2), dims=2) ≈ ones(Int64, 10)
end

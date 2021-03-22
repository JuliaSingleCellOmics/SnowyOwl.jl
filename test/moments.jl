@testset "moments" begin
    ind = rand([0,1], 10)    
    @test sum(normalize_indicator(ind)) ≈ 1.
end

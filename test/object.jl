@testset "object" begin
    r, c = (100, 500)
    data = rand([0, 1], r, c)
    obs = DataFrame(A=rand(c), B=rand(c))
    var = DataFrame(C=rand(r), D=rand(r))
    @test_throws AssertionError Profile(data, var, obs)

    prof = Profile(data, obs, var)
    @test obsnames(prof) == ["A", "B"]
    @test varnames(prof) == ["C", "D"]
    @test nrow(prof) == r
    @test ncol(prof) == c
    @test nvar(prof) == r
    @test nobs(prof) == c
    @test maximum(prof) == 1
    @test minimum(prof) == 0
    @test size(prof) == (r, c)
    @test axes(prof) == (Base.OneTo(r), Base.OneTo(c))

    prof2 = prof[1:50, :]
    @test prof2.data == prof.data[1:50, :]

    idx = rand([false,true], r)
    prof3 = prof[idx, :]
    @test prof3.data == prof.data[idx, :]
    @test prof3.obs == prof.obs
    @test prof3.var == prof.var[idx, :]
end

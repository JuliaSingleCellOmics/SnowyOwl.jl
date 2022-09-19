@testset "qc" begin
    r, c = (100, 500)
    data = rand(0:10, r, c)
    obs = DataFrame(index=1:c, A=rand(c), B=rand(c))
    var = DataFrame(index=1:r, C=rand(r), D=rand(r), mt=rand(Bool,r))
    prof = Profile(data, var, obs)

    quality_control_metrics!(prof)
    @test prof.obs[!, :n_genes_by_counts] == count(data .!= 0, dims=1)
    @test prof.obs[!, :total_counts] == sum(data, dims=1)
    @test prof.obs[!, :log1p_n_genes_by_counts] = log1p.(count(data .!= 0, dims=1))

    prof = Profile(data, var, obs)
    prof2 = copy(prof)
    @test prof2 !== prof
    @test prof2.data == prof.data
    @test prof2.var == prof.var
    @test prof2.obs == prof.obs
    
    prof2 = quality_control_metrics(prof)

end
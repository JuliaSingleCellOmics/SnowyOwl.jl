@testset "qc" begin
    r, c = (100, 500)
    data = rand(0:10, r, c)
    obs = DataFrame(index=1:c, A=rand(c), B=rand(c))
    var = DataFrame(index=1:r, C=rand(r), D=rand(r), mt=rand(Bool,r))
    prof = Profile(data, var, obs)

    prof1 = copy(prof)
    quality_control_metrics!(prof, log1p=true)
    @test prof1.obs[!, :n_genes_by_counts] == count(data .!= 0, dims=1)
    @test prof1.obs[!, :total_counts] == sum(data, dims=1)
    @test prof1.obs[!, :log1p_n_genes_by_counts] == log1p.(count(data .!= 0, dims=1))
    @test prof1.obs[!, :log1p_total_counts] == log1p.(sum(data, dims=1))
    @test prof1.obs[!, :total_counts_mt] == sum(data[:, var[!, :mt]])
    @test prof1.obs[!, :log1p_total_counts_mt] == log1p.(sum(data[:, var[!, :mt]]))
    @test prof1.obs[!, :pct_total_counts_mt] == ( sum(data[:, var[!, :mt]]) ./ sum(data, dims=1) ) .* 100
    @test prof1.var[!, :n_cells_by_counts] == count(data .!= 0, dims=2)
    @test prof1.var[!, :mean_counts] == mean(data, dims=2)
    @test prof1.var[!, :pct_dropout_by_counts] == (1 .- count(data .!= 0, dims=2) ./ size(data, 1)) .* 100
    @test prof1.var[!, :total_counts] == sum(data, dims=2)
    @test prof1.var[!, :log1p_mean_counts] == log1p.(mean(data, dims=2))
    @test prof1.var[!, :log1p_total_counts] == log1p.(sum(data, dims=2))

    prof2 = copy(prof)
    quality_control_metrics(prof2, log1p=true)
    @test prof2 == prof1

    quality_control_metrics!(data, obs, var, log1p=true)
    @test obs[!, :n_genes_by_counts] == count(data .!= 0, dims=1)
    @test obs[!, :total_counts] == sum(data, dims=1)
    @test obs[!, :log1p_n_genes_by_counts] == log1p.(count(data .!= 0, dims=1))
    @test obs[!, :log1p_total_counts] == log1p.(sum(data, dims=1))
    @test obs[!, :total_counts_mt] == sum(data[:, var[!, :mt]])
    @test obs[!, :log1p_total_counts_mt] == log1p.(sum(data[:, var[!, :mt]]))
    @test obs[!, :pct_total_counts_mt] == ( sum(data[:, var[!, :mt]]) ./ sum(data, dims=1) ) .* 100
    @test var[!, :n_cells_by_counts] == count(data .!= 0, dims=2)
    @test var[!, :mean_counts] == mean(data, dims=2)
    @test var[!, :pct_dropout_by_counts] == (1 .- count(data .!= 0, dims=2) ./ size(data, 1)) .* 100
    @test var[!, :total_counts] == sum(data, dims=2)
    @test var[!, :log1p_mean_counts] == log1p.(mean(data, dims=2))
    @test var[!, :log1p_total_counts] == log1p.(sum(data, dims=2))

    obs_metrics, var_metrics = quality_control_metrics(data, obs, var, log1p=true)
    @test obs_metrics == obs
    @test var_metrics == var
end
@testset "qc" begin
    r, c = (100, 500)
    data = rand(0:10, r, c)
    obs = DataFrame(index=1:c, A=rand(c), B=rand(c))
    var = DataFrame(index=1:r, C=rand(r), D=rand(r), mt=rand(Bool,r))
    prof = Profile(data, :RNA, var, obs, varindex=:index, obsindex=:index)

    obs_metrics, var_metrics = quality_control_metrics(data, obs, var, obsname=:index, varname=:index, use_log1p=true)
    @test obs_metrics[!, :n_genes_by_counts] == vec(count(data .!= 0, dims=1))
    @test obs_metrics[!, :total_counts] == vec(sum(data, dims=1))
    @test obs_metrics[!, :log1p_n_genes_by_counts] == vec(log1p.(count(data .!= 0, dims=1)))
    @test obs_metrics[!, :log1p_total_counts] == vec(log1p.(sum(data, dims=1)))
    @test obs_metrics[!, :total_counts_mt] == vec(sum(data[var[!, :mt], :], dims=1))
    @test obs_metrics[!, :log1p_total_counts_mt] == vec(log1p.(sum(data[var[!, :mt], :], dims=1)))
    @test obs_metrics[!, :pct_total_counts_mt] == vec((sum(data[var[!, :mt], :], dims=1) ./ sum(data, dims=1)) .* 100)
    @test var_metrics[!, :n_cells_by_counts] == vec(count(data .!= 0, dims=2))
    @test var_metrics[!, :mean_counts] == vec(mean(data, dims=2))
    @test var_metrics[!, :pct_dropout_by_counts] == vec((1 .- count(data .!= 0, dims=2) ./ size(data, 1)) .* 100)
    @test var_metrics[!, :total_counts] == vec(sum(data, dims=2))
    @test var_metrics[!, :log1p_mean_counts] == vec(log1p.(mean(data, dims=2)))
    @test var_metrics[!, :log1p_total_counts] == vec(log1p.(sum(data, dims=2)))

    quality_control_metrics!(data, obs, var, use_log1p=true)
    @test obs[!, names(obs_metrics)] == obs_metrics
    @test var[!, names(var_metrics)] == var_metrics

    prof1 = quality_control_metrics(prof, use_log1p=true)
    @test prof1.obs == prof.obs
    @test prof1.RNA.var == prof.RNA.var

    quality_control_metrics!(prof1, use_log1p=true)
    @test prof1.obs == prof.obs
    @test prof1.RNA.var == prof.RNA.var
end

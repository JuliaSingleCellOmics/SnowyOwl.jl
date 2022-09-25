@testset "qc" begin
    r, c = (100, 500)
    data = rand(0:10, r, c)
    obs = DataFrame(index=1:c, A=rand(c), B=rand(c))
    var = DataFrame(index=1:r, C=rand(r), D=rand(r), mt=rand(Bool,r))
    prof = Profile(data, :RNA, var, obs, varindex=:index, obsindex=:index)

    prof1 = copy(prof)
    quality_control_metrics!(prof1, use_log1p=true)
    @test prof1.obs[!, :n_genes_by_counts] == vec(count(data .!= 0, dims=1))
    @test prof1.obs[!, :total_counts] == vec(sum(data, dims=1))
    @test prof1.obs[!, :log1p_n_genes_by_counts] == vec(log1p.(count(data .!= 0, dims=1)))
    @test prof1.obs[!, :log1p_total_counts] == vec(log1p.(sum(data, dims=1)))
    @test prof1.obs[!, :total_counts_mt] == vec(sum(data[var[!, :mt], :], dims=1))
    @test prof1.obs[!, :log1p_total_counts_mt] == vec(log1p.(sum(data[var[!, :mt], :], dims=1)))
    @test prof1.obs[!, :pct_total_counts_mt] == vec((sum(data[var[!, :mt], :], dims=1) ./ sum(data, dims=1)) .* 100)
    @test prof1.RNA.var[!, :n_cells_by_counts] == vec(count(data .!= 0, dims=2))
    @test prof1.RNA.var[!, :mean_counts] == vec(mean(data, dims=2))
    @test prof1.RNA.var[!, :pct_dropout_by_counts] == vec((1 .- count(data .!= 0, dims=2) ./ size(data, 1)) .* 100)
    @test prof1.RNA.var[!, :total_counts] == vec(sum(data, dims=2))
    @test prof1.RNA.var[!, :log1p_mean_counts] == vec(log1p.(mean(data, dims=2)))
    @test prof1.RNA.var[!, :log1p_total_counts] == vec(log1p.(sum(data, dims=2)))

    prof2 = quality_control_metrics(prof, use_log1p=true)
    @test prof2 == prof1

    quality_control_metrics!(data, obs, var, use_log1p=true)
    @test obs[!, :n_genes_by_counts] == vec(count(data .!= 0, dims=1))
    @test obs[!, :total_counts] == vec(sum(data, dims=1))
    @test obs[!, :log1p_n_genes_by_counts] == vec(log1p.(count(data .!= 0, dims=1)))
    @test obs[!, :log1p_total_counts] == vec(log1p.(sum(data, dims=1)))
    @test obs[!, :total_counts_mt] == vec(sum(data[var[!, :mt], :], dims=1))
    @test obs[!, :log1p_total_counts_mt] == vec(log1p.(sum(data[var[!, :mt], :], dims=1)))
    @test obs[!, :pct_total_counts_mt] == vec((sum(data[var[!, :mt], :], dims=1) ./ sum(data, dims=1)) .* 100)
    @test var[!, :n_cells_by_counts] == vec(count(data .!= 0, dims=2))
    @test var[!, :mean_counts] == vec(mean(data, dims=2))
    @test var[!, :pct_dropout_by_counts] == vec((1 .- count(data .!= 0, dims=2) ./ size(data, 1)) .* 100)
    @test var[!, :total_counts] == vec(sum(data, dims=2))
    @test var[!, :log1p_mean_counts] == vec(log1p.(mean(data, dims=2)))
    @test var[!, :log1p_total_counts] == vec(log1p.(sum(data, dims=2)))

    obs_metrics, var_metrics = quality_control_metrics(data, obs, var, obsname=:index, varname=:index, use_log1p=true)
    @test obs_metrics == obs[!, [:index, :n_genes_by_counts, :total_counts, :log1p_n_genes_by_counts, :log1p_total_counts, :total_counts_mt, :log1p_total_counts_mt, :pct_total_counts_mt]]
    @test var_metrics == var[!, [:index, :n_cells_by_counts, :mean_counts, :pct_dropout_by_counts, :total_counts, :log1p_mean_counts, :log1p_total_counts]]
end

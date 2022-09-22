@testset "filter" begin
    r, c = (100, 500)
    data = rand(0:10, r, c)
    obs = DataFrame(A=rand(c), B=rand(c))
    var = DataFrame(C=rand(r), D=rand(r))
    layer_a = rand(r, c)
    layer_b = rand(c, r)

    prof = Profile(data, var, obs)
    prof.layers[:a] = copy(layer_a)
    prof.layers[:b] = copy(layer_b)
    
    min_genes = 90
    expressed_genes = vec(sum(data .!= 0, dims=1))
    SnowyOwl.Preprocess.filter_cells!(prof; min_genes=min_genes)
    @test prof.data == data[:, min_genes .<= expressed_genes]
    @test prof.obs == obs[min_genes .<= expressed_genes, :]
    @test prof.var == var
    @test prof.layers[:a] == layer_a[:, min_genes .<= expressed_genes]
    @test prof.layers[:b] == layer_b
    @test_throws ArgumentError SnowyOwl.Preprocess.filter_cells!(prof; min_genes=-1)
    @test_throws ArgumentError SnowyOwl.Preprocess.filter_cells!(prof; max_genes=-1)
    @test_throws ArgumentError SnowyOwl.Preprocess.filter_cells!(prof; min_counts=-1)
    @test_throws ArgumentError SnowyOwl.Preprocess.filter_cells!(prof; max_counts=-1)

    prof = Profile(data, var, obs)
    prof.layers[:a] = copy(layer_a)
    prof.layers[:b] = copy(layer_b)
    
    min_cells = 300
    expressed_cells = vec(sum(data .!= 0, dims=2))
    SnowyOwl.Preprocess.filter_genes!(prof, min_cells=min_cells)
    @test prof.data == data[min_cells .<= expressed_cells, :]
    @test prof.obs == obs
    @test prof.var == var[min_cells .<= expressed_cells, :]
    @test prof.layers[:a] == layer_a[min_cells .<= expressed_cells, :]
    @test prof.layers[:b] == layer_b
    @test_throws ArgumentError SnowyOwl.Preprocess.filter_genes!(prof; min_cells=-1)
    @test_throws ArgumentError SnowyOwl.Preprocess.filter_genes!(prof; max_cells=-1)
    @test_throws ArgumentError SnowyOwl.Preprocess.filter_genes!(prof; min_counts=-1)
    @test_throws ArgumentError SnowyOwl.Preprocess.filter_genes!(prof; max_counts=-1)
end

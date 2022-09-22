@testset "highly_variable" begin
    ngenes, ncells = (100, 500)
    nhvgs = 10
    X = rand(NegativeBinomial(1, 0.8), ngenes, ncells)
    var = DataFrame(gene_symbols=1:ngenes, A=rand(ngenes))

    min_disp, max_disp = 0.5, Inf
    min_mean, max_mean = 0.0125, 3.
    hvg = highly_variable_genes(X, var, :seurat; varname=:gene_symbols,
                                min_disp=min_disp, max_disp=max_disp,
                                min_mean=min_mean, max_mean=max_mean)
    top_genes = (min_mean .< hvg.means .< max_mean) .&
                (min_disp .< hvg.dispersions_norm .< max_disp)
    @test names(hvg) == ["gene_symbols", "means", "dispersions", "dispersions_norm",
                         "highly_variable"]
    @test all(top_genes .== hvg.highly_variable)

    highly_variable_genes!(X, var, :cellranger; ntop_genes=nhvgs)
    @test names(var) == ["gene_symbols", "A", "means", "dispersions", "dispersions_norm",
                         "highly_variable"]
    @test count(var.highly_variable) == nhvgs
end

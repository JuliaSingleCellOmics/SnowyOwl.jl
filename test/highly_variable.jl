@testset "highly_variable" begin
    ngenes, ncells = (100, 500)
    nhvgs = 10
    X = rand(NegativeBinomial(1, 0.8), ngenes, ncells)
    var = DataFrame(gene_symbols=1:ngenes, A=rand(ngenes))

    min_disp, max_disp = 0.5, Inf
    min_mean, max_mean = 0.0125, 3.
    hvg = SnowyOwl.Preprocess.highly_variable_genes(X, var, :seurat; varname=:gene_symbols,
                                                    min_disp=min_disp, max_disp=max_disp,
                                                    min_mean=min_mean, max_mean=max_mean)
    top_genes = (min_mean .< hvg.means .< max_mean) .&
                (min_disp .< hvg.dispersions_norm .< max_disp)
    @test names(hvg) == ["gene_symbols", "means", "dispersions", "dispersions_norm",
                         "highly_variable"]
    @test all(top_genes .== hvg.highly_variable)

    SnowyOwl.Preprocess.highly_variable_genes!(X, var, :cellranger; ntop_genes=nhvgs)
    @test names(var) == ["gene_symbols", "A", "means", "dispersions", "dispersions_norm",
                         "highly_variable"]
    @test count(var.highly_variable) == nhvgs

    obs = DataFrame(barcode=1:ncells)
    var = DataFrame(gene_symbols=1:ngenes, A=rand(ngenes))

    prof = Profile(X, :RNA, var, obs; varindex=:gene_symbols, obsindex=:barcode)
    prof2 = SnowyOwl.Preprocess.highly_variable_genes(prof, :cellranger; layer=:count,
                                                      min_disp=min_disp, max_disp=max_disp,
                                                      min_mean=min_mean, max_mean=max_mean)
    top_genes = (min_mean .< prof2.RNA.var.means .< max_mean) .&
                (min_disp .< prof2.RNA.var.dispersions_norm .< max_disp)
    @test names(prof2.RNA.var) == ["gene_symbols", "A", "means", "dispersions",
                                   "dispersions_norm", "highly_variable"]
    @test all(top_genes .== prof2.RNA.var.highly_variable)

    SnowyOwl.Preprocess.highly_variable_genes!(prof, :seurat; omicsname=:RNA, layer=:count,
                                               ntop_genes=nhvgs)
    @test names(prof.RNA.var) == ["gene_symbols", "A", "means", "dispersions",
                                  "dispersions_norm", "highly_variable"]
    @test count(prof.RNA.var.highly_variable) == nhvgs
    @test haskey(prof.RNA.pipeline, :hvg)
end

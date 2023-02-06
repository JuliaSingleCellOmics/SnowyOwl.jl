@testset "project" begin
    ngene = 10
    ncell = 1000
    X = rand(ngene, ncell)

    Y = SnowyOwl.Analysis.project(X, SnowyOwl.Analysis.PCAMethod(); dims=5)
    @test size(Y) == (5, ncell)

    newY = similar(X, 5, ncell)
    SnowyOwl.Analysis.project!(newY, X, SnowyOwl.Analysis.PCAMethod(); dims=5)
    @test newY == Y

    Y = SnowyOwl.Analysis.project(X, SnowyOwl.Analysis.UMAPMethod(); dims=2)
    @test size(Y) == (2, ncell)

    obs = DataFrame(barcode=1:ncell)
    var = DataFrame(gene_symbols=1:ngene)
    prof = Profile(X, :RNA, var, obs; varindex=:gene_symbols, obsindex=:barcode)

    prof2 = SnowyOwl.Analysis.project(prof, SnowyOwl.Analysis.PCAMethod();
                                      dims=5, omicsname=:RNA, layer=:count)
    @test size(prof2.obsm[:RNA_count_pca]) == (5, ncell)
    @test haskey(prof2.RNA.pipeline, :pca)

    SnowyOwl.Analysis.project!(prof, SnowyOwl.Analysis.PCAMethod();
                               dims=5, omicsname=:RNA, layer=:count)
    @test haskey(prof.RNA.pipeline, :pca)
end

@testset "transform" begin
    @testset "logarithmize" begin
        ngenes, ncells = (100, 500)
        X = fill(2., ngenes, ncells)

        logX = SnowyOwl.Preprocess.logarithmize(X)
        @test logX == fill(log1p(2.), ngenes, ncells)

        SnowyOwl.Preprocess.logarithmize!(X; base=2)
        @test X == fill(log1p(2.)/log(2.), ngenes, ncells)

        obs = DataFrame(barcode=1:ncells)
        var = DataFrame(gene_symbols=1:ngenes)

        X = fill(2., ngenes, ncells)
        prof = Profile(X, :RNA, var, obs; varindex=:gene_symbols, obsindex=:barcode)
        prof2 = SnowyOwl.Preprocess.logarithmize(prof; layer=:count)
        @test prof2.RNA.count == fill(log1p(2.), ngenes, ncells)

        SnowyOwl.Preprocess.logarithmize!(prof; omicsname=:RNA, layer=:count,
                                                   base=2)
        @test prof.RNA.count == fill(log1p(2.)/log(2.), ngenes, ncells)
        @test haskey(prof.RNA.pipeline, :log1p)
    end
end

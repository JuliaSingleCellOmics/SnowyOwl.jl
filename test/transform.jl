@testset "transform" begin
    @testset "normalize" begin
        ngenes, ncells = (100, 500)
        X = rand(ngenes, ncells) * 20

        @testset "relative normalization" begin
            normX = SnowyOwl.Preprocess.normalize(X, SnowyOwl.Preprocess.RelativeNormalization())
            library_size = sum(normX, dims=1)
            @test all(library_size .≈ 1e4)
        end

        @testset "log normalization" begin
            normX = SnowyOwl.Preprocess.normalize(X, SnowyOwl.Preprocess.LogNormalization())
            library_size = sum(expm1.(normX), dims=1)
            @test all(library_size .≈ 1e4)
        end

        @testset "clr normalization" begin
            normX = SnowyOwl.Preprocess.normalize(X, SnowyOwl.Preprocess.CenteredLogRatioNormalization())
            library_size = sum(expm1.(normX), dims=1)
            @test_skip all(library_size .≈ 1e4) # TODO: not sure how to test
        end

        SnowyOwl.Preprocess.normalize!(X, SnowyOwl.Preprocess.RelativeNormalization())
        library_size = sum(X, dims=1)
        @test all(library_size .≈ 1e4)

        obs = DataFrame(barcode=1:ncells)
        var = DataFrame(gene_symbols=1:ngenes)

        X = rand(ngenes, ncells) * 20
        prof = Profile(X, :RNA, var, obs; varindex=:gene_symbols, obsindex=:barcode)
        prof2 = SnowyOwl.Preprocess.normalize(prof, SnowyOwl.Preprocess.RelativeNormalization())
        library_size = sum(prof2.RNA.count, dims=1)
        @test all(library_size .≈ 1e4)

        SnowyOwl.Preprocess.normalize!(prof, SnowyOwl.Preprocess.RelativeNormalization())
        library_size = sum(prof.RNA.count, dims=1)
        @test all(library_size .≈ 1e4)
        @test haskey(prof.RNA.pipeline, :normalize)
    end

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

@testset "cuda" begin
    ngenes, ncells = (100, 500)

    @testset "normalize" begin
        X = rand(ngenes, ncells) * 20 |> cu

        @testset "relative normalization" begin
            normX = SnowyOwl.Preprocess.normalize(X, SnowyOwl.Preprocess.RelativeNormalization())
            library_size = collect(sum(normX, dims=1))
            @test all(library_size .≈ 1e4)
        end

        @testset "log normalization" begin
            normX = SnowyOwl.Preprocess.normalize(X, SnowyOwl.Preprocess.LogNormalization())
            library_size = collect(sum(expm1.(normX), dims=1))
            @test all(library_size .≈ 1e4)
        end

        @testset "clr normalization" begin
            normX = SnowyOwl.Preprocess.normalize(X, SnowyOwl.Preprocess.CenteredLogRatioNormalization())
            library_size = collect(sum(expm1.(normX), dims=1))
            @test_skip all(library_size .≈ 1e4) # TODO: not sure how to test
        end
    end

    @testset "logarithmize" begin
        X = fill(2., ngenes, ncells) |> cu

        logX = SnowyOwl.Preprocess.logarithmize(X)
        @test collect(logX) == fill(log1p(2.f0), ngenes, ncells)
    end

    @testset "highly_variable" begin
        nhvgs = 10
        X = Float32.(rand(NegativeBinomial(1, 0.8), ngenes, ncells)) |> cu
        var = DataFrame(gene_symbols=1:ngenes, A=rand(ngenes))

        @testset "Seurat method" begin
            min_disp, max_disp = 0.5, Inf
            min_mean, max_mean = 0.0125, 3.
            hvg = SnowyOwl.Preprocess.highly_variable_genes(X, var, SnowyOwl.Preprocess.SeuratHVG();
                                                            varname=:gene_symbols,
                                                            min_disp=min_disp, max_disp=max_disp,
                                                            min_mean=min_mean, max_mean=max_mean)
            top_genes = (min_mean .< hvg.means .< max_mean) .&
                        (min_disp .< hvg.dispersions_norm .< max_disp)
            @test names(hvg) == ["gene_symbols", "means", "dispersions", "dispersions_norm",
                                 "highly_variable"]
            @test all(top_genes .== hvg.highly_variable)
        end

        @testset "cell ranger method" begin
            hvg = SnowyOwl.Preprocess.highly_variable_genes(X, var, SnowyOwl.Preprocess.CellRangerHVG();
                                                            ntop_genes=nhvgs)
            @test names(hvg) == ["gene_symbols", "means", "dispersions", "dispersions_norm",
                                 "highly_variable"]
            @test count(hvg.highly_variable) == nhvgs
        end
    end
end

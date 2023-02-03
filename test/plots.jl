@testset "plots" begin
    @testset "scatter" begin
        x = collect(1:10)
        y = collect(1:10)
        z = collect(1:10)

        p = SnowyOwl.Plots.scatterplot(x, y)
        @test p isa Plots.Plot
        p = SnowyOwl.Plots.scatterplot(x, y, z)
        @test p isa Plots.Plot

        p = SnowyOwl.Plots.pca(x, y)
        @test p isa Plots.Plot
        p = SnowyOwl.Plots.pca(x, y, z)
        @test p isa Plots.Plot

        p = SnowyOwl.Plots.umap(x, y)
        @test p isa Plots.Plot
        p = SnowyOwl.Plots.umap(x, y, z)
        @test p isa Plots.Plot
    end

    @testset "save" begin
        x = collect(1:10)
        y = collect(1:10)

        p = SnowyOwl.Plots.scatterplot(x, y)

        filepath = joinpath(TEST_PATH, "test_plot")
        SnowyOwl.Plots.save(p, filepath, nothing)
        @test !isfile(filepath * ".png")

        SnowyOwl.Plots.save(p, filepath, :png)
        @test isfile(filepath * ".png")
        rm(filepath * ".png")

        SnowyOwl.Plots.save(p, filepath, [:svg, :png])
        @test isfile(filepath * ".png")
        @test isfile(filepath * ".svg")
        rm(filepath * ".png")
        rm(filepath * ".svg")
    end
end

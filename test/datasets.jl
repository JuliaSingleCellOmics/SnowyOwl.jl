@testset "datasets" begin
    DATADEPS_DIR = joinpath(homedir(), ".julia", "datadeps")
    
    @testset "PBMC3k" begin
        @test SnowyOwl.Dataset.pbmc3k_folder() == joinpath(DATADEPS_DIR, "PBMC3k", "filtered_gene_bc_matrices", "hg19")
    end
end

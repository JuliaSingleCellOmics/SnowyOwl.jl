@testset "io" begin
    folder = joinpath(pkgdir(SnowyOwl), "HumanColonicMesenchymeIBD",
                      "5e5a9f34-734a-4501-920c-3ecf5e2e2f36.mtx")
    mtx = read_mtx(joinpath(folder, "matrix.mtx.gz"))
    features = read_genes(joinpath(folder, "genes.tsv.gz"))
    cells = read_cells(joinpath(folder, "cells.tsv.gz"))
end
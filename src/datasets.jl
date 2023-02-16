module Dataset

using DataDeps
using DataFrames
using JLD2

export load_pbmc68k

function __init__pbmc3k()
    DEPNAME = "PBMC3k"

    register(DataDep(
            DEPNAME,
            """
            Dataset: The 3k Peripheral Blood Mononuclear Cells (PBMC)
            Authors: 10X Genomics
            Website: https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k/
            """,
            "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
            "847d6ebd9a1ec9a768f2be7e40ca42cbfe75ebeb6d76a4c24167041699dc28b5",
            post_fetch_method = DataDeps.unpack
        ))
end

pbmc3k_folder() = datadep"PBMC3k/filtered_gene_bc_matrices/hg19"

function load_pbmc68k(local_path)
    @load local_path data
    expr_mat = data[:all_data][Symbol("17820")][:hg19]

    obs = DataFrame(barcode=expr_mat[:barcodes])
    var = DataFrame(gene_id=expr_mat[:genes], gene_symbol=expr_mat[:gene_symbols])
    prof = Profile(SparseMatrixCSC(expr_mat[:mat]'), :RNA, var, obs; varindex=gene_id, obsindex=barcode)
    prof
end

end

function load_pbmc68k(local_path)
    @load local_path data
    expr_mat = data[:all_data][Symbol("17820")][:hg19]

    obs = DataFrame(barcode=expr_mat[:barcodes])
    var = DataFrame(gene_id=expr_mat[:genes], gene_symbol=expr_mat[:gene_symbols])
    prof = Profile(SparseMatrixCSC(expr_mat[:mat]'), :RNA, var, obs; varindex=gene_id, obsindex=barcode)
    prof
end

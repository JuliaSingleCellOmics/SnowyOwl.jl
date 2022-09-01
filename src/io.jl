function read_mtx(filename::String)
    py"""
    from scipy.io import mmread
    mtx = mmread($filename)
    """
    sparse(py"mtx.row" .+ 1, py"mtx.col" .+ 1, py"mtx.data")
end

function read_tsv_gz(filename::String, header::Vector{Symbol}=Symbol[])
    gz_unzipper = transcode(GzipDecompressor, Mmap.mmap(filename))
    header = isempty(header) ? 1 : header
    return CSV.File(gz_unzipper; delim='\t', header=header) |> DataFrame
end

read_features(filename::String, names=DEFAULT_FEATURE_COLS) = read_tsv_gz(filename, names)
read_barcodes(filename::String, names=DEFAULT_BARCODE_COLS) = read_tsv_gz(filename, names)

## Human Cell Atlas (HCA)
read_genes(filename::String) = read_tsv_gz(filename)
read_cells(filename::String) = read_tsv_gz(filename)

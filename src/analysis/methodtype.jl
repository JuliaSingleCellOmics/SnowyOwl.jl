abstract type AnalysisMethod end

struct PCAMethod <: AnalysisMethod end

Base.show(io::IO, ::PCAMethod) = print(io, "PCA")

Base.Symbol(::PCAMethod) = :pca

struct UMAPMethod <: AnalysisMethod end

Base.show(io::IO, ::UMAPMethod) = print(io, "UMAP")

Base.Symbol(::UMAPMethod) = :umap

abstract type NormalizationMethod end

struct LogNormalization <: NormalizationMethod end

Base.show(io::IO, ::LogNormalization) = print(io, "log normalization")

Base.Symbol(::LogNormalization) = :lognormalize

"""
    RelativeNormalization()

Normalize count data to relative counts per cell by dividing by the total counts per cell.
"""
struct RelativeNormalization <: NormalizationMethod end

Base.show(io::IO, ::RelativeNormalization) = print(io, "relative normalization")

Base.Symbol(::RelativeNormalization) = :relative

struct CenteredLogRatioNormalization <: NormalizationMethod end

Base.show(io::IO, ::CenteredLogRatioNormalization) = print(io, "centered log ratio normalization")

Base.Symbol(::CenteredLogRatioNormalization) = :clr

struct CustomNormalization <: NormalizationMethod end

Base.show(io::IO, ::CustomNormalization) = print(io, "custom normalization")

Base.Symbol(::CustomNormalization) = :custom

abstract type HighlyVariableMethod end

"""
    CellRangerHVG()

Calculate highly variable features using cell ranger method. The top `ntop_genes` most
variable genes were identified based on their mean and dispersion (variance/mean).
Genes were grouped based on their mean expression. Within each bin, z score is calculated
from median absolute deviation of all genes within the bin. Genes are identified whose
z-scores were higher than `min_disp`.

Zheng et al. (2017), Massively parallel digital transcriptional profiling of single cells, Nature Communications.
"""
struct CellRangerHVG <: HighlyVariableMethod end

Base.show(io::IO, ::CellRangerHVG) = print(io, "cell ranger")

Base.Symbol(::CellRangerHVG) = :cellranger

"""
    SeuratHVG()

Calculate highly variable features using Seurat method. The mean and dispersion (variance/mean)
are calculated for each gene across all cells. Genes were grouped based on their mean
expression. Within each bin, z-score is calculated from the dispersion of all genes within
the bin. Genes are identified whose z-scores were higher than `min_disp`.

Satija et al. (2015), Spatial reconstruction of single-cell gene expression data, Nature Biotechnology.
"""
struct SeuratHVG <: HighlyVariableMethod end

Base.show(io::IO, ::SeuratHVG) = print(io, "seurat")

Base.Symbol(::SeuratHVG) = :seurat

"""
    Seuratv3HVG()

Calculate highly variable features using Seurat v3 method.

Stuart et al. (2019), Comprehensive integration of single-cell data Cell.
"""
struct Seuratv3HVG <: HighlyVariableMethod end

Base.show(io::IO, ::Seuratv3HVG) = print(io, "seurat v3")

Base.Symbol(::Seuratv3HVG) = :seuratv3

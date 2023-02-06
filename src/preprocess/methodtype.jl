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

struct CellRangerHVG <: HighlyVariableMethod end

Base.show(io::IO, ::CellRangerHVG) = print(io, "cell ranger")

Base.Symbol(::CellRangerHVG) = :cellranger

struct SeuratHVG <: HighlyVariableMethod end

Base.show(io::IO, ::SeuratHVG) = print(io, "seurat")

Base.Symbol(::SeuratHVG) = :seurat

struct Seuratv3HVG <: HighlyVariableMethod end

Base.show(io::IO, ::Seuratv3HVG) = print(io, "seurat v3")

Base.Symbol(::Seuratv3HVG) = :seuratv3

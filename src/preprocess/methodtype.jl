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

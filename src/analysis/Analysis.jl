module Analysis

using MultivariateStats
using OmicsProfiles
using UMAP

export
    # methodtype
    AnalysisMethod,
    PCAMethod,
    UMAPMethod,

    # project
    project!,
    project

include("methodtype.jl")
include("project.jl")

end

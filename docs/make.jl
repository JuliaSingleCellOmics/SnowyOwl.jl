using SnowyOwl
using Documenter

makedocs(;
    modules=[SnowyOwl],
    authors="Yueh-Hua Tu",
    repo="https://github.com/JuliaSingleCellOmics/SnowyOwl.jl/blob/{commit}{path}#L{line}",
    sitename="SnowyOwl.jl",
    clean = false,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliasinglecellomics.github.io/SnowyOwl.jl/stable/",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaSingleCellOmics/SnowyOwl.jl",
    target = "build",
)

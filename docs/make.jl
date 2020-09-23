using SnowyOwl
using Documenter

makedocs(;
    modules=[SnowyOwl],
    authors="Yueh-Hua Tu",
    repo="https://github.com/yuehhua/SnowyOwl.jl/blob/{commit}{path}#L{line}",
    sitename="SnowyOwl.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://yuehhua.github.io/SnowyOwl.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/yuehhua/SnowyOwl.jl",
)

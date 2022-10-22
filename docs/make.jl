using UnfoldSim
using Documenter

DocMeta.setdocmeta!(UnfoldSim, :DocTestSetup, :(using UnfoldSim); recursive=true)

makedocs(;
    modules=[UnfoldSim],
    authors="Luis Lips, Benedikt Ehinger, Judith Schepers",
    repo="https://github.com/behinger/UnfoldSim.jl/blob/{commit}{path}#{line}",
    sitename="UnfoldSim.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://behinger.github.io/UnfoldSim.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/behinger/UnfoldSim.jl",
    devbranch="main",
)

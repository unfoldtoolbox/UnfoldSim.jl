using UnfoldSim
using Documenter
using Glob
using Literate


GENERATED = joinpath(@__DIR__, "src", "literate")
for subfolder ∈ ["explanations","HowTo","tutorials","reference"]
    local SOURCE_FILES = Glob.glob(subfolder*"/*.jl", GENERATED)
    #config=Dict(:repo_root_path=>"https://github.com/unfoldtoolbox/UnfoldSim")
    foreach(fn -> Literate.markdown(fn, GENERATED*"/"*subfolder), SOURCE_FILES)

end


DocMeta.setdocmeta!(UnfoldSim, :DocTestSetup, :(using UnfoldSim); recursive=true)

makedocs(;
    modules=[UnfoldSim],
authors="Luis Lips, Benedikt Ehinger, Judith Schepers",
    repo="https://github.com/behinger/UnfoldSim.jl/blob/{commit}{path}#{line}",
    sitename="UnfoldSim.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://unfoldtoolbox.github.io/UnfoldSim.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials"=>[
                "Quickstart" => "literate/tutorials/quickstart.md",
                "Simulate ERPs" => "literate/tutorials/simulateERP.md",
                "Poweranalysis" => "literate/tutorials/poweranalysis.md",
        ],
        "Reference"=>[
                "NoiseTypes" => "./literate/reference/noisetypes.md",
                "ComponentBasisTypes" => "./literate/reference/basistypes.md",
        ],
        "HowTo" => [
                "New Experimental Design" => "./literate/HowTo/newDesign.md",
        ],
        "DocStrings" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/unfoldtoolbox/UnfoldSim.jl",
    devbranch="main",
)

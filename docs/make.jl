using UnfoldSim
using Documenter
using Glob
using Literate


GENERATED = joinpath(@__DIR__, "src", "generated")
SOURCE = joinpath(@__DIR__, "literate")

for subfolder ∈ ["explanations", "HowTo", "tutorials", "reference"]
    local SOURCE_FILES = Glob.glob(subfolder * "/*.jl", SOURCE)
    #config=Dict(:repo_root_path=>"https://github.com/unfoldtoolbox/UnfoldSim")
    foreach(fn -> Literate.markdown(fn, GENERATED * "/" * subfolder), SOURCE_FILES)

end


DocMeta.setdocmeta!(UnfoldSim, :DocTestSetup, :(using UnfoldSim); recursive = true)

makedocs(;
    modules = [UnfoldSim],
    authors = "Luis Lips, Benedikt Ehinger, Judith Schepers",
    #repo="https://github.com/unfoldtoolbox/UnfoldSim.jl/blob/{commit}{path}#{line}",
    repo = Documenter.Remotes.GitHub("unfoldtoolbox", "UnfoldSim.jl"),
    sitename = "UnfoldSim.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://unfoldtoolbox.github.io/UnfoldSim.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Quickstart" => "generated/tutorials/quickstart.md",
            "Simulate ERPs" => "generated/tutorials/simulateERP.md",
            "Poweranalysis" => "generated/tutorials/poweranalysis.md",
            "Multi-subject simulation" => "generated/tutorials/multisubject.md",
        ],
        "Reference" => [
            "Overview: Toolbox Functions" => "./generated/reference/overview.md",
            "Overview: NoiseTypes" => "./generated/reference/noisetypes.md",
            "Overview: OnsetTypes" => "./generated/reference/onsettypes.md",
            "Overview: Components (EEG, fMRI, Pupil)" => "./generated/reference/basistypes.md",
        ],
        "HowTo" => [
            "Define a new, (imbalanced) design" => "./generated/HowTo/newDesign.md",
            "Repeating a design" => "./generated/HowTo/repeatTrials.md",
            "Define a new duration & jitter component" => "./generated/HowTo/newComponent.md",
            "Generate multi channel data" => "./generated/HowTo/multichannel.md",
            "Use predefined design / onsets data" => "./generated/HowTo/predefinedData.md",
        ],
        "API / DocStrings" => "api.md",
    ],
)

deploydocs(;
    repo = "github.com/unfoldtoolbox/UnfoldSim.jl",
    devbranch = "main",
    versions = "v#.#",
    push_preview = true,
)

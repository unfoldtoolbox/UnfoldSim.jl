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
    authors = "Judith Schepers, Luis Lips, Maanik Marathe, Benedikt Ehinger",
    #repo="https://github.com/unfoldtoolbox/UnfoldSim.jl/blob/{commit}{path}#{line}",
    repo = Documenter.Remotes.GitHub("unfoldtoolbox", "UnfoldSim.jl"),
    sitename = "UnfoldSim.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://unfoldtoolbox.github.io/UnfoldSim.jl",
        edit_link = "main",
        sidebar_sitename = false,
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Quickstart" => "generated/tutorials/quickstart.md",
            "Simulate event-related potentials (ERPs)" => "generated/tutorials/simulateERP.md",
            "Power analysis" => "generated/tutorials/poweranalysis.md",
            "Multi-subject simulation" => "generated/tutorials/multisubject.md",
        ],
        "Reference" => [
            "Overview of functionality" => "./generated/reference/overview.md",
            "Overview: Experimental design types" => "./generated/reference/designtypes.md",
            "Overview: Basis function (component) types" => "./generated/reference/basistypes.md",
            "Overview: Onset types" => "./generated/reference/onsettypes.md",
            "Overview: Noise types" => "./generated/reference/noisetypes.md",
        ],
        "HowTo" => [
            "Define a new (imbalanced) design" => "./generated/HowTo/newDesign.md",
            "Get multiple trials with identical subject/item combinations" => "./generated/HowTo/repeatTrials.md",
            "Define a new component (with variable duration and shift)" => "./generated/HowTo/newComponent.md",
            "Generate multi channel data" => "./generated/HowTo/multichannel.md",
            "Use existing experimental designs & onsets in the simulation" => "./generated/HowTo/predefinedData.md",
            "Get ground truth via EffectsDesign" => "./generated/HowTo/getGroundTruth.md",
        ],
        "API / Docstrings" => "api.md",
    ],
)

# deploydocs(;
#     repo = "github.com/unfoldtoolbox/UnfoldSim.jl",
#     #devbranch = "main",
#     #versions = "v#.#",
#     push_preview = true,
# )

deploydocs(;
    repo = "github.com/unfoldtoolbox/UnfoldSim.jl",
    devbranch = "main",
    #versions = ["stable" => "v^", "v#.#.#"],
    push_preview = true,
)

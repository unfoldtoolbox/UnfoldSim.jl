
using Test
using Random
using StableRNGs
using Statistics
using LinearAlgebra
using MixedModelsSim
using DataFrames
using Distributions # For LogNormal function

function gen_debug_design(; n_subjects = 20, n_item = 100)
    # define design parameters
    item_btwn = Dict(:stimType => ["A", "B"])

    # instantiate the design
    return MultiSubjectDesign(;
        n_subjects = n_subjects,
        n_items = n_item,
        items_between = item_btwn,
    )
end

function gen_debug_component()
    basisfunction = [zeros(10) ones(10) zeros(10)]
    formula = @formula(0 ~ 1 + stimType + (1 + stimType | subj))
    contrasts = Dict(:stimType => DummyCoding())
    β = [2.0, 0.5]
    σ_ranef = Dict(:subj => [1, 0.0])


    # instantiate the component(s)
    return MixedModelComponent(;
        basis = basisfunction,
        formula = formula,
        β = β,
        σs = σ_ranef,
        contrasts = contrasts,
    )
    #return MixedModelComponent(basisfunction, formula, contrasts, β, σ_ranef)
end

function gen_debug_simulation(;
    design = gen_debug_design(),
    component = gen_debug_component(),
    noisetype = PinkNoise(),
    onset = UniformOnset(; width = 0, offset = 100),
)
    return Simulation(design, component, onset, noisetype)
end

# # Get ground truth via EffectsDesign

# Usually, to test a method, you want to compare your results to a known ground truth. In UnfoldSim you can obtain your ground truth via the `EffectsDesign`.
# Doing it this way let's you marginalize any effects/ variables of your original design. You can find more on what marginalized effects are here in the [Unfold.jl documentation](https://unfoldtoolbox.github.io/Unfold.jl/dev/generated/HowTo/effects/)

# ## Setup
using UnfoldSim
using Unfold
using CairoMakie
using Random

# ## Simulation
# First let's make up a SingleSubject simulation

# !!! note
#     Getting a ground truth for a MultiSubjectDesign is not implemented yet

design =
    SingleSubjectDesign(;
        conditions = Dict(
            :condition => ["bike", "face"],
            :continuous => range(0, 5, length = 10),
        ),
    ) |> x -> RepeatDesign(x, 100);

# **n170** has a condition effect, faces are more negative than bikes
n1 = LinearModelComponent(;
    basis = n170(),
    formula = @formula(0 ~ 1 + condition),
    β = [5, 3],
);
# **p300** has a continuous effect, higher continuous values will result in larger P300's.
# We include both a linear and a quadratic effect of the continuous variable.
p3 = LinearModelComponent(;
    basis = p300(),
    formula = @formula(0 ~ 1 + continuous + continuous^2),
    β = [5, 1, 0.2],
);

components = [n1, p3]
data, evts = simulate(
    MersenneTwister(1),
    design,
    components,
    UniformOnset(; width = 0, offset = 1000),
    PinkNoise(),
);

# ## GroundTruthDesign
# To marginalize effects we first have to specify an effects dictionary and subsequently hand this dict plus the original design to `EffectsDesign()`

effects_dict = Dict(:condition => ["bike", "face"])

effects_design = EffectsDesign(design, effects_dict)

# !!! note
#     We only specified the condition levels here, by default every unspecified variable will be set to a "typical" (i.e. the mean) value.

# And finally we can simulate our ground truth ERP with marginalized effects

gt_data, gt_events = simulate(
    MersenneTwister(1),
    effects_design,
    components,
    UniformOnset(; width = 0, offset = 1000),
    NoNoise(),
    return_epoched = true,
);

# ## Compare with Unfold.jl results

m = fit(
    UnfoldModel,
    Dict(
        Any => (
            @formula(0 ~ 1 + condition + spl(continuous, 4)),
            firbasis(τ = [-0.1, 1], sfreq = 100, name = "basis"),
        ),
    ),
    evts,
    data,
);

eff = effects(effects_dict, m)

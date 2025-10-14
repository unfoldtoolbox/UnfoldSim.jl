# # Simulate ground truth marginalized Effects

# Often when testing some algorithm, we want to compare our results to a known ground truth. In the case of marginalized effects via the `Unfold.effects`/ `Effects.jl` interface, we can do this using an `EffectsDesign`.
# You can find more on what marginalized effects are here in the [Unfold.jl documentation](https://unfoldtoolbox.github.io/Unfold.jl/dev/generated/HowTo/effects/)

# ### Setup
# ```@raw html
# <details>
# <summary>Click to expand</summary>
# ```
## Load required packages
using UnfoldSim
using Unfold
using CairoMakie
using UnfoldMakie
using Random
# ```@raw html
# </details >
# ```
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
    ) |> x -> RepeatDesign(x, 5);

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

# ## Simulate marginalized effects directly
# To marginalize effects we first have to specify an effects dictionary and subsequently hand this dict plus the original design to `EffectsDesign()`

effects_dict = Dict(:condition => ["bike", "face"])

effects_design = EffectsDesign(design, effects_dict)

# !!! note
#     We only specified the condition levels here, by default every unspecified variable will be set to a "typical" (i.e. the mean) value.

# And finally we can simulate our ground truth marginal effects

gt_data, gt_events = simulate(
    MersenneTwister(1),
    effects_design,
    components,
    NoOnset(),
    NoNoise(),
    return_epoched = true,
);
@show gt_events

# Additionally, we can get the simulated effects into a tidy dataframe using Unfold's `result_to_table`.
# Note that the data has to be reshaped into a channel X times X predictor form. (In our one channel example `size(gt_data) = (45,2)`, missing the channel dimension)

g = reshape(gt_data, 1, size(gt_data)...)
times = range(1, size(gt_data, 1));
gt_effects = Unfold.result_to_table([g], [gt_events], [times], ["effects"])
first(gt_effects, 5)


# ## Compare with Unfold.jl results

m = fit(
    UnfoldModel,
    [
        Any => (
            @formula(0 ~ 1 + condition + spl(continuous, 4)),
            firbasis(τ = [-0.1, 1], sfreq = 100, name = "basis"),
        ),
    ],
    evts,
    data,
);

ef = effects(effects_dict, m);

# !!! note
#      The ground truth is shorter because the ground truth typically returns values between `[0 maxlength(components)]`, whereas in our unfold-model we included a baseline period of 0.1s.
#      If you want to actually compare results with the ground truth, you could either us `UnfoldSim.pad_array()` or set the Unfold modelling window to `τ=[0,1]`

gt_effects.type .= "UnfoldSim effects"
ef.type .= "Unfold effects"

gt_effects.time = gt_effects.time ./ 100 .- 1 / 100
ef.continuous .= 2.5 # needed to be able to easily merge the two dataframes
comb = vcat(gt_effects, ef)
plot_erp(comb; mapping = (; color = :type, col = :condition))

# The simulated ground truth marginal effects, and the fitted marginal effects look similar as expected, but the fitted has some additional noise because of finite data (also as expected).
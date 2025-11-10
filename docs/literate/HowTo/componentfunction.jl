# # [Define design-dependent component basis functions](@id componentfunction)
# Here you will learn how to specify a component basis that uses a function that depends on the `design` and returns per-event basis vectors, instead of the same basis vector for all events.


# ### Setup
# ```@raw html
# <details>
# <summary>Click to expand</summary>
# ```
## Load required packages
using UnfoldSim
using Unfold
using Random
using DSP
using CairoMakie, UnfoldMakie
# ```@raw html
# </details >
# ```


sfreq = 100;

# ## Design
# Let's generate a design with a categorical effect and a continuous duration effect
design = UnfoldSim.SingleSubjectDesign(;
    conditions = Dict(
        :category => ["dog", "cat"],
        :duration => Int.(round.(20 .+ rand(100) .* sfreq)),
    ),
);


# Instead of defining a "boring" vector basis function e.g. `[0,0,1,2,3,3,2,1,0,0,0]`, let's use a function - in our case, a Hanning window with the size depending on the experimental design's duration.
# !!! important
#     Two things have to be taken care of:
#     1. in case a rng is required to e.g. generate the design, or your basis function depends on it, you have to specify a two-argument basis function: `(rng,design)->...`
#     2. a `maxlength` in samples has to be specified via a tuple `(function,maxlength)`

mybasisfun = design -> hanning.(generate_events(design).duration)
signal = LinearModelComponent(;
    basis = (mybasisfun, 100),
    formula = @formula(0 ~ 1 + category),
    β = [1, 0.5],
);

erp = UnfoldSim.simulate_component(MersenneTwister(1), signal, design);

# After simulation, we are ready to plot it. We expect that the length of the simulated responses is scaled by the design's duration. To show it more effectively, we sort by duration.
##---
f = Figure()
df = UnfoldMakie.eeg_array_to_dataframe(erp')
df.duration = repeat(generate_events(design).duration, inner = size(erp, 1))
df.category = repeat(generate_events(design).category, inner = size(erp, 1))
plot_erp!(
    f[1, 1],
    df,
    mapping = (; group = :group => nonnumeric, col = :category), #  color = :duration, fails right nowUnfoldMakie#353
    layout = (; legend_position = :left),
    colorbar = (; label = "Duration"),
)
plot_erpimage!(
    f[2, 1],
    erp,
    sortvalues = generate_events(design).duration,
    layout = (; legend_position = :bottom),
)
f

# The scaling by the two `condition` effect levels and the modified event duration by the `duration` are clearly visible.

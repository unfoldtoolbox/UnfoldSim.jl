# # Component Functions
# HowTo put arbitrary functions into components

using UnfoldSim
using Unfold
using Random
using DSP
using CairoMakie, UnfoldMakie

sfreq = 100;

# ## Design
# Let's generate a design with a categorical effect and a continuous duration effect
design = UnfoldSim.SingleSubjectDesign(;
    conditions = Dict(
        :category => ["dog", "cat"],
        :duration => Int.(round.(20 .+ rand(100) .* sfreq)),
    ),
);


# Instead of defining a boring vector basis function e.g. `[0,0,1,2,3,3,2,1,0,0,0]`, let's use function, generating random values for now.
# !!! important
#     because any function depending on `design` can be used, two things have to be taken care of:
#     
#     1. in case a random component exist, specify a `RNG`, the basis might be evaluated multiple times inside `simulate`
#     2. a `maxlength` has to be specified via a tuple `(function.maxlength)``
mybasisfun = design -> hanning.(generate_events(design).duration)
signal = LinearModelComponent(;
    basis = (mybasisfun, 100),
    formula = @formula(0 ~ 1 + category),
    β = [1, 0.5],
);

erp = UnfoldSim.simulate_component(MersenneTwister(1), signal, design);


# Finally, let's plot it, sorted by duration
plot_erpimage(erp, sortvalues = generate_events(design).duration)

# # Quickstart

# To get started with data simulation, the user needs to provide four ingredients: an experimental design, an event basis function (component), an inter-onset distribution and a noise specification.

# !!! tip
#       Use `subtypes(AbstractNoise)` (or `subtypes(AbstractComponent)` etc.) to find already implemented building blocks.

# ## Specify the simulation ingredients

# ### Setup
# ```@raw html
# <details>
# <summary>Click to expand</summary>
# ```
## Load required packages
using UnfoldSim
using Random # to get an RNG
using CairoMakie # for plotting

# ```@raw html
# </details >
# ```

# ### Experimental Design
# Define a 1 x 2 design with 20 trials. That is, one condition (`cond_A`) with two levels.
design =
    SingleSubjectDesign(; conditions = Dict(:cond_A => ["level_A", "level_B"])) |>
    x -> RepeatDesign(x, 10);

# ### Event basis function (Component)
# Define a simple component and ground truth simulation formula. Akin to ERP components, we call one simulation signal a component.
# 
# !!! note 
#        You could easily specify multiple components by providing a vector of components, which are automatically added at the same onsets. This procedure simplifies to generate some response that is independent of simulated condition, whereas other depends on it.
signal = LinearModelComponent(;
    basis = [0, 0, 0, 0.5, 1, 1, 0.5, 0, 0],
    formula = @formula(0 ~ 1 + cond_A),
    β = [1, 0.5],
);

# ### Onsets and Noise
# We will start with a uniform (but overlapping, `offset` < `length(signal.basis)`) inter-onset distribution.
onset = UniformOnset(; width = 20, offset = 4);

# And we will use some noise
noise = PinkNoise(; noiselevel = 0.2);

# ## Combine & Generate
# Finally, we will combine all ingredients and simulate some data.
data, events = simulate(MersenneTwister(1), design, signal, onset, noise);
# `data` is a `n-sample` Vector (but could be a Matrix for e.g. `MultiSubjectDesign` or epoched data).

# `events` is a DataFrame that contains a column `latency` with the onsets of events (in samples).

# ## Plot them!
lines(data; color = "black")
vlines!(events.latency; color = ["orange", "teal"][1 .+ (events.cond_A.=="level_B")])

current_axis().title = "Simulated data"
current_axis().xlabel = "Time [samples]"
current_axis().ylabel = "Amplitude [μV]"

current_figure()

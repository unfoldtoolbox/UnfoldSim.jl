using UnfoldSim
using Random
using CairoMakie

# !!! Tipp
#        Use `subtypes(AbstractNoise)` (or `subtypes(AbstractComponent)` etc.) to find already implemented building blocks

# #### "Experimental" Design
# Define a 1 x 2 design with 20 trials. That is, one condition (`condaA`) with two levels.
design = SingleSubjectDesign(;
        n_trials=10,
        conditions=Dict(:condA=>["levelA","levelB"])
        );

# #### Component / Signal
# Define a simple component and ground truth simulation formula. Akin to ERP components, we call one simulation signal a component.
# 
# !!! highlight
#        You could easily specify multiple components by providing a vector of components, which are automatically added at the same onsets. This procedure simplifies to generate some response that is independent of simulated condition, whereas other depends on it.
signal = LinearModelComponent(;
        basis=[0,0,0,0.5,1,1,0.5,0,0],
        formula = @formula(0~1+condA),
        β = [1,0.5]
        );

# #### Onsets and Noise
# We will start with a uniform (but overlapping, `offset` < `length(signal.basis)`) onset-distribution
onset = UniformOnset(;width=20,offset=4);

# And we will use some noise
noise = PinkNoise(;noiselevel=0.2);

# #### Combine & Generate
# We will put it all together in one `Simulation` type
simulation = Simulation(design, signal,  onset, noise);

# finally, we will simulate some data
data,events = simulate(MersenneTwister(1),simulation);
# Data is a `n-sample` Vector (but could be a Matrix for e.g. `MultiSubjectDesign`).

# events is a DataFrame that contains a column `latency` with the onsets of events.

# ## Plot them!
lines(data;color="black")
vlines!(events.latency;color=["orange","teal"][1 .+ (events.condA.=="levelB")])
current_figure()
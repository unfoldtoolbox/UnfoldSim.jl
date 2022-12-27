using UnfoldSim
using Random
using CairoMakie

# !!! Tipp
#   Use `subtypes(AbstractNoise)`` (or AbstractComponent etc.) to find already implemented buildingblocks

# Define a 1 x 2 design with 20 trials
design = SingleSubjectDesign(;
        n_trials=10,
        conditions=Dict(:condA=>["levelA","levelB"])
        );

# Define a simple component + ground truth formula
signal = LinearModelComponent(;
        basis=[0,0,0,0.5,1,1,0.5,0,0],
        formula = @formula(0~1+condA),
        β = [1,0.5]
        );

# start with uniform (but overlapping, `offset` < `length(signal.basis)` onsets
onset = UniformOnset(;width=20,offset=4);

# and some noise
noise = PinkNoise(;noiselevel=0.2);

# put it all together
simulation = Simulation(design, signal,  onset, noise);

# ## Generate data
data,events = simulate(MersenneTwister(1),simulation);

# ## Plot them!
lines(data;color="black")
vlines!(events.latency;color=["orange","teal"][1 .+ (events.condA.=="levelB")])
current_figure()
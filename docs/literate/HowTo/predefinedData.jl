# # Use existing experimental designs & onsets in the simulation

# Let's say you want to use the events data frame (containing the levels of the experimental variables and the event onsets (latencies)) from a previous study in your simulation.

# ### Setup
# ```@raw html
# <details>
# <summary>Click to expand</summary>
# ```
## Load required packages
using UnfoldSim
using DataFrames
using Random
using CairoMakie # for plotting
# ```@raw html
# </details >
# <br />
# ```

# From a previous study, we (somehow, e.g. by using [pyMNE.jl](https://unfoldtoolbox.github.io/Unfold.jl/dev/HowTo/pymne/)) imported an event data frame like this:
my_events = DataFrame(:condition => [:A, :B, :B, :A, :A], :latency => [7, 13, 22, 35, 41])

# To use exactly these values, we can generate a new `AbstractDesign`, which will always return this event dataframe
struct MyManualDesign <: AbstractDesign
    my_events::Any
end
UnfoldSim.generate_events(d::MyManualDesign) = deepcopy(d.my_events) ## generate function which is called internally in UnfoldSim
UnfoldSim.size(d::MyManualDesign) = size(d.my_events, 1); ## necessary function to tell what the dimensionality of the experimental design is
# !!! note
#     Note the `UnfoldSim.generate_events` which tells Julia to "overload" the `generate_events` function as defined in UnfoldSim.


# Next we generate a `MyManualDesign`
mydesign = MyManualDesign(my_events);

# We could already use this "solo" and simulate some data, for example:
signal = LinearModelComponent(;
    basis = [1, 1, 0.5, 0, 0],
    formula = @formula(0 ~ 1 + condition),
    β = [1, 0.5],
);

data, events =
    simulate(MersenneTwister(1), mydesign, signal, UniformOnset(; width = 10, offset = 5))
lines(data) # plotting 
vlines!(my_events.latency; linestyle = :dash)
current_figure()
# Looks good, but the events don't match our custom onsets yet.

# ## Custom Timings
# Finally, we want to use our custom timings as well. For this we define a new `AbstractOnset`. Again, it simply returns our manually provided latencies
struct MyManualOnset <: AbstractOnset end
UnfoldSim.simulate_onsets(rng, onset::MyManualOnset, simulation::Simulation) =
    generate_events(simulation.design).latency
# !!! hint
#     This is a bit of a trick, it relies that `MyManualOnset` is always used in combination with `MyManualDesign`. You could of course repeat the structure from `MyManualDesign` also for `MyManualOnset` and have an explicit field in the structure containing the onsets.

# And that's it
data, events = simulate(MersenneTwister(1), mydesign, signal, MyManualOnset())
lines(data) # plotting 
vlines!(my_events.latency, linestyle = :dash)
current_figure()
# now everything matches, lovely!
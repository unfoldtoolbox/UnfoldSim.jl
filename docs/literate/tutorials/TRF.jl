using UnfoldSim
using CairoMakie
using Random

# # Temporal response functions
# So far, we simulated event-related potentials. In this tutorial, we will switch gears a bit, and simulate a continuous response to a continuous feature vector. 
#
# One good example is the visual response to a grayscale circle, which continuously changes its luminance (this goes back to VESPA, described in the ground breaking [Lalor 2006 paper](10.1016/j.neuroimage.2006.05.054). The brain will react to this luminance change continuously. TRFs typically describe the process to recover the brain response (in terms of a filter response). Here we set out to simulate such a dataset first and foremost.
#
# ## Simulation
# we start with the simplest possible design, one condition
design = SingleSubjectDesign(conditions=Dict(:dummy=>["A"]));

# next we define the convolution kernel that the feature signal should be convolved with (= the brain response == the TRF)
brain_response = [LinearModelComponent(basis=p100(),formula=@formula(0~1),β=[1]),
                  LinearModelComponent(basis=n170(),formula=@formula(0~1),β=[1]),
                  LinearModelComponent(basis=p300(),formula=@formula(0~1),β=[1])];

# !!! hint
#     For a fancy way to write the above code, you could use `LinearModelComponent.([p100,n170,p300],@formula(0~1),Ref([1]))`, notice the `.(...)` indicating broadcasting

# Now we can simulate our feature signal, 10_000 samples of random gray scale values
feature = rand(10_000)

# Next we have to nest the response in a `TRFComponent` and add ou
trf_component = TRFComponent(brain_response,feature)

# Finally, when simulating, we have only a single "event" (that is, TRF-) onset, the first sample. Therefore, we use  `TRFOnset` to indicate this.
dat,evts = simulate(design,trf_component,UnfoldSim.TRFOnset())
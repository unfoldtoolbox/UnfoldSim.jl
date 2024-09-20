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
feature = rand(1_000)

# Next we have to nest the response in a `TRFComponent` and add ou
trf_component = TRFComponent(brain_response,feature);

# Finally, when simulating, we have only a single "event" (that is, TRF-) onset, the first sample. Therefore, we use  `TRFOnset` to indicate this.
dat,evts = simulate(design,trf_component,UnfoldSim.TRFOnset());

# Let's plot the feature signal and the TRF response
f,ax,h =lines(dat)
lines!(feature)
f

# ## Multivariate TRFs
# Now TRFs can depend on multiple variables e.g. the luminance and the size of the circle.
feature_luminance = rand(1_000)
feature_size = rand(1_000)

# We could call the `simulate` function twice and simply add the results:
dat_l,_ = simulate(design,TRFComponent(brain_response,feature_luminance),UnfoldSim.TRFOnset());
dat_r,_ = simulate(design,TRFComponent(brain_response,feature_size),UnfoldSim.TRFOnset());

# let's plot and compare to the previous plot
f,ax,h = lines(dat_l .+ dat_r)
lines!(dat)
f
# as visible, the blue line (luminance+size) has ~double the amplitude. This is because we simulated  two brain responses and simply added them.

# A bit more convenient way is to do the following
dat_combined,_ = simulate(design,[TRFComponent(brain_response,feature_size),TRFComponent(brain_response,feature_luminance)],UnfoldSim.TRFOnset());
f,ax,h = lines(dat_combined)
lines!(dat_l .+ dat_r)
f
# where you can see that the results are equivalent.

# ## Another cool feature is to modulate the feature vector based on the design
design_mod = SingleSubjectDesign(conditions=Dict(:condition=>["A","B"]));

# Let's take only a single component for simplicity. Note how the `formula` has been changed. The β allows to control the amplitude. In this linear model component, the default contrast-function is `Dummy` (also known as `Reference` coding), which means, the second beta indicated a "difference"
brain_response_mod = LinearModelComponent(basis=p100(),formula=@formula(0~1+condition),β=[1,1]);

# let's simulate another feature signal, but this time, we simulate a Matrix.
feature_mod = rand(1000,2)
feature_mod[:,2] .= 0
feature_mod[500:600,2] .= 1;

# to better understand how our (experimental) feature now looks like, let's visualize it
series(feature_mod')
# As visible, the first column has again a random signal, indicating e.g. luminance changes. The second temporal feature indicates some offset (a colorchange?) between 500 and 600 samples.


dat_mod,_ = simulate(design_mod,TRFComponent([brain_response_mod],feature_mod),UnfoldSim.TRFOnset());
lines(dat_mod)

# !!! hint
#     Excourse: now you rightfully might ask why the jump is to >10 a.u.? The reason is that you are effectively convolving a feature-signal above 1 (feature_mod * (β_0 + β_1) = 1 * (1+1)), and with a kernel with maximum = 1, these 2's add up. The kernel has a rough width of ~5 which results in additional 5*2 => 10 per affected sample

# ## Combination with a event-related design
# Right now there is no easy interface to do this. You have to simulate a TRF signal, and an rERP signal and then add the two signals.

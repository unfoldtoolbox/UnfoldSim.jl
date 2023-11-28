
using UnfoldSim
using UnfoldMakie
using CairoMakie
using DataFrames
using Random


# ## Specifying a design

# We are using a one-level design for testing here.
design = SingleSubjectDesign(conditions=Dict(:condA=>["levelA"]))

# Next we generate two simple components at two different times without any formula attached (we have a single condition anyway)
c = LinearModelComponent(;basis=p100(),formula = @formula(0~1),β = [1]);
c2 = LinearModelComponent(;basis=p300(),formula = @formula(0~1),β = [1]);


# ## The multichannel component
# next similar to the nested design above, we can nest the component in a `MultichannelComponent`. We could either provide the projection marix manually, e.g.:
mc = UnfoldSim.MultichannelComponent(c, [1,2,-1,3,5,2.3,1])

# or maybe more convenient: use the pair-syntax: Headmodel=>Label which makes use of a headmodel (HaRTmuT is currently easily available in UnfoldSim)
hart = headmodel(type="hartmut")
mc = UnfoldSim.MultichannelComponent(c, hart=>"Left Postcentral Gyrus")
mc2 = UnfoldSim.MultichannelComponent(c2, hart=>"Right Occipital Pole")

# !!! hint
#     You could also specify a noise-specific component which is applied prior to projection & summing with other components
# 
# finally we need to define the onsets of the signal
onset = UniformOnset(;width=20,offset=4);

# ## Simulation

# Now as usual we simulate data. Inspecting data shows our result is now indeed ~230 Electrodes large! Nice!
data,events = simulate(MersenneTwister(1),design, [mc,mc2],  onset, PinkNoise(noiselevel=0.05)) 
size(data)


# !!! hint
#     The noise declared in the `simulate` function is added after mixing to channels, each channel receives independent noise. It is also possible to add noise to each individual component+source prior to projection. This would introduce correlated noise.
#
# ## Plotting
# Let's plot using Butterfly & Topoplot
# first we convert the electrodes to positions usable in TopoPlots.jl
pos3d = hart.electrodes["pos"];
pos2d = to_positions(pos3d')
pos2d = [Point2f(p[1]+0.5,p[2]+0.5) for p in pos2d];

# now plot!
f = Figure()
df = DataFrame(:estimate => data[:],:channel => repeat(1:size(data,1),outer=size(data,2)),:time => repeat(1:size(data,2),inner=size(data,1)))
plot_butterfly!(f[1,1:2],df;positions=pos2d)
plot_topoplot!(f[2,1],df[df.time .== 28,:];positions=pos2d,visual=(;enlarge=0.5,label_scatter=false),axis=(;limits=((0,1),(0,0.9))))
plot_topoplot!(f[2,2],df[df.time .== 48,:];positions=pos2d,visual=(;enlarge=0.5,label_scatter=false),axis=(;limits=((0,1),(0,0.9))))
f


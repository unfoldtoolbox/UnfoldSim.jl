# UnfoldSim

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://unfoldtoolbox.github.io/UnfoldSim.jl/dev/)
[![Build Status](https://github.com/unfoldtoolbox/UnfoldSim.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/unfoldtoolbox/UnfoldSim.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/v/UnfoldSim.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/unfoldtoolbox/UnfoldSim.jl)
[![DOI](https://zenodo.org/badge/413455526.svg)](https://zenodo.org/badge/latestdoi/413455526)




A package to simulate single timeseries model-based ERPs, fMRI activity, pupil dilation etc.
If you have one channel, it is a timeseries of (overlapping) event-related activity and some noise - you might have fun here!

Many Tutorials, Guides, How-Tos and References available in the documentation!

![grafik](https://user-images.githubusercontent.com/10183650/213565922-90feec23-3b51-4602-b50c-31561dbfc261.png)



# Quickstart
```julia
data,evts = UnfoldSim.predef_eeg(;n_trials=20,noiselevel=0.8)
```
Produces continuous "EEG" with PinkNoise and some Overlap between 20 events

# Slightly longer
```julia
# start by defining the design / event-table
design = SingleSubjectDesign(;conditions=Dict(:condA=>["levelA","levelB"])) |> d->RepeatDesign(d,10);
# next define a ground-truth signal + relation to events/design with Wilkinson Formulas
signal = LinearModelComponent(;
        basis=[0,0,0,0.5,1,1,0.5,0,0],
        formula = @formula(0~1+condA),
        β = [1,0.5]
        );
# finally, define some Onset Distribution and Noise, and simulate!
data,events = simulate(Random.MersenneTwister(1),design, signal,  UniformOnset(;offset=5,width=4), PinkNoise());        
```
All components (design, components, onsets, noise) can be easily modified and you can simply plugin your own!

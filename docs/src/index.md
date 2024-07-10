```@meta
CurrentModule = UnfoldSim
```

# UnfoldSim.jl

Documentation for [UnfoldSim.jl](https://github.com/unfoldtoolbox/UnfoldSim.jl). UnfoldSim.jl is a Julia package for simulating multivariate timeseries data with a special focus on EEG data.

## Start simulating timeseries
We offer some predefined (EEG) signals, check them out!

For instance a P1/N170/P300 complex (containing three typical ERP components).
```@example main
using UnfoldSim
using CairoMakie # plotting

data, evts = UnfoldSim.predef_eeg(; n_repeats = 1, noiselevel = 0.8)

lines(data; color = "black")
vlines!(evts.latency; color = ["orange", "teal"][1 .+ (evts.condition.=="car")])

current_figure()
```

## Or simulate epoched data directly
```@example main

data, evts = UnfoldSim.predef_eeg(; n_repeats = 20, noiselevel = 0.8, return_epoched = true)
heatmap(data[:, sortperm(evts, [:condition, :continuous])])

```

```@meta
CurrentModule = UnfoldSim
```

# UnfoldSim

Documentation for [UnfoldSim](https://github.com/behinger/UnfoldSim.jl).

## Start simulating timeseries
We offer some predefined signals, check them out!

For instance an P1/N170/P300 complex.
```@example main
using UnfoldSim
using CairoMakie
data,evts = UnfoldSim.predef_eeg(;n_repeats=1,noiselevel=0.8)

lines(data;color="black")
vlines!(evts.latency;color=["orange","teal"][1 .+ (evts.condition .=="car")])

current_figure()
```

## Or simulate epoched data directly
```@example main

data,evts = UnfoldSim.predef_eeg(;n_repeats=20,noiselevel=0.8,return_epoched=true)
heatmap(data[:,sortperm(evts,[:condition,:continuous])])

```

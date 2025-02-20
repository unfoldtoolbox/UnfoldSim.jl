```@meta
CurrentModule = UnfoldSim
```

# UnfoldSim.jl Documentation

Welcome to [UnfoldSim.jl](https://github.com/unfoldtoolbox/UnfoldSim.jl): a Julia package for simulating multivariate timeseries data with a focus on EEG, especially event-related potentials (ERPs). 
The user provides four ingredients: 1) an experimental design, with both categorical and continuous variables, 2) event basis functions specified via linear or hierarchical models, 3) an inter-event onset distribution, and 4) a noise specification.

```@raw html
<div style="width:60%; margin: auto;">

<img src="https://github.com/unfoldtoolbox/UnfoldSim.jl/blob/assets/docs/src/assets/UnfoldSim_features_animation.gif?raw=true"/>
</div>
```

## Key features
- **Modularity:** Choose or implement different designs, components, onset distributions or noise types
- **Multi-subject** & complex experimental designs
- **Multi-channel** via EEG-forward models
- **Continuous or epoched** data with potentially **overlapping** signals
- Potential support for **other modalities**, e.g. single-voxel fMRI or pupil dilation

## Installation
```julia-repl
julia> using Pkg; Pkg.add("UnfoldSim")
```
For more detailed instructions please refer to [Installing Julia & UnfoldSim.jl](@ref).

## Usage example
#### Start simulating time series data
We offer some predefined (EEG) signals, check them out!
For instance, a P1/N170/P300 complex (containing three typical ERP components).

```@example quickstart
using UnfoldSim
using CairoMakie # plotting #hide
data, events = UnfoldSim.predef_eeg(; n_repeats = 1, noiselevel = 0.8)

lines(data;color="black")
vlines!(evts.latency;color=["orange","teal"][1 .+ (evts.condition .=="car")])

current_figure() #hide
```

#### Or simulate epoched data directly
```@example quickstart

data,evts = UnfoldSim.predef_eeg(;n_repeats=20,noiselevel=0.8,return_epoched=true)
heatmap(data[:,sortperm(evts,[:condition,:continuous])])

```

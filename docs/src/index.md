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

f = Figure(size = (900, 400)) #hide
ax = Axis( #hide
    f[1, 1], #hide
    title="Simulated EEG data", #hide
    titlesize=18, #hide
    xlabel="Time [samples]", #hide
    ylabel="Amplitude [µV]", #hide
    xlabelsize=16, #hide
    ylabelsize=16, #hide
    xgridvisible=false, #hide
    ygridvisible=false, #hide
) #hide

colors = ifelse.(events.condition .== "face", :orange , :teal) #hide

lines!(data; color="black") #hide
vlines!(events.latency; color=colors, label=events.condition) #hide

xlims!(current_axis(), [0, 400]) #hide

current_figure() #hide
```

#### Or simulate epoched data directly
```@example quickstart

data, events = UnfoldSim.predef_eeg(; n_repeats = 20, noiselevel = 0.8, return_epoched = true)

f2 = Figure(size = (900, 400)) #hide
ax = Axis( #hide
    f2[1, 1], #hide
    title="ERP image", #hide
    titlesize=18, #hide
    xlabel="Time [samples]", #hide
    ylabel="Trials", #hide
    xlabelsize=16, #hide
    ylabelsize=16, #hide
    xgridvisible=false, #hide
    ygridvisible=false, #hide
) #hide
hm = heatmap!(data[:, sortperm(events, [:condition, :continuous])]) #hide
Colorbar(f2[:, end+1], hm, label = "Amplitude [µV]") #hide

current_figure() #hide
```

## Where to start: Learning roadmap
#### 1. First steps
📌 Goal: Learn about the simulation workflow and run your first simulation\
🔗 [Quickstart](@ref) | [Simulate event-related potentials (ERPs)](@ref)

#### 2. Intermediate topics
📌 Goal: Learn about multi-subject and multi-channel simulations \
🔗 [Multi-subject simulation](@ref) | [Generate multi channel data](@ref)

#### 3. Advanced topics
📌 Goal: Learn how to implement your own experimental designs, components, etc.\
🔗 [Define a new (imbalanced) design](@ref) | [Define a new component (with variable duration and shift)](@ref) | [Use existing experimental designs & onsets in the simulation](@ref)


## Statement of need
EEG researchers often analyze data containing (temporally) overlapping events (e.g. stimulus onset and button press, or consecutive eye-fixations), non-linear effects, and complex experimental designs. For a multitude of reasons, we often need to simulate such kinds of data: Simulated EEG data is useful to test preprocessing and analysis tools, validate statistical methods, illustrate conceptual issues, test toolbox functionalities, and find limitations of traditional analysis workflows. For instance, such simulation tools allow for testing the assumptions of new analysis algorithms and testing their robustness against any violation of these assumptions.

```@raw html
<!---
Note: The statement of need is also used in the `README.md`. Make sure that they are synchronized.
-->
```

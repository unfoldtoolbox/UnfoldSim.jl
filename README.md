# UnfoldSim

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://s-ccs.github.io/UnfoldSim.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://s-ccs.github.io/UnfoldSim.jl/dev/)
[![Build Status](https://github.com/s-ccs/UnfoldSim.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/s-ccs/UnfoldSim.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/s-ccs/UnfoldSim.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/s-ccs/UnfoldSim.jl)


A package to simulate single timeseries model-based ERPs, fMRI activity, pupil dilation etc.
If you have one channel, it is a timeseries of (overlapping) event-related activity and some noise - you might have fun here!

**Note: Not released yet**

# Usage

1) Define the experiment design
```julia
# define design parameters
n_subj = 2
n_item = 2
btwn_subj = nothing
btwn_item = Dict("stimType" => ["A", "B"])
both_win = nothing

# instantiate the design
design = ExperimentDesign(n_subj, n_item, btwn_subj, btwn_item, both_win)
```

2) Define the ERP component(s)
```julia
# define component parameters
basisfunction = padarray(hanning(25), (-50, 50), 0)
formula = @formula(dv ~ 1 + stimType + (1 + stimType | subj) + (1 | item))
contrasts = Dict(:stimType => DummyCoding())
β = [2.0, 0.5]
σ_ranef = Dict(:subj => create_re(0.2, 0.0), :item=>create_re(0.0))
σ_res = 0.0001

# instantiate the component(s)
p100 = Component(basisfunction, formula, contrasts, β, σ_ranef, σ_res)
```

3) Set a seed
```julia
rng = MersenneTwister(2)
```

4) Put all together into a simulation
```julia
# some more parameters

noisetype = PinkNoise(;noiselevel=2)
onset = UniformOnset(;offset=200,width=50) # 200-250 samples distance between events

simulation = Simulation(design, [p100],  onset,noisetype)
```

5) Simulate EEG data
```julia
eeg, onsets = simulate(rng, simulation);
```

5.1) Convert to eeg, onsets to data, evts (similar to Unfold)
```julia
data, evts = UnfoldSim.convert(eeg, onsets, design)
```

6) Simulate only ERPs (without noise)
```julia
erps = simulate_erps(rng, design, [p100]);
```

using UnfoldSim
using CairoMakie
using DSP
using StableRNGs

# ## Basistypes
# There are several basis types directly implemented. They can be easily used for the `components`.
#
# !!! note
#       You can use any arbitrary shape defined by yourself! We often make use of `hanning(50)` from the DSP.jl package.

# ## EEG
# By default, the EEG bases assume a sampling rate of 100, which can easily be changed by e.g. p100(;sfreq=300)
f = Figure()
ax = f[1,1] = Axis(f)
for b in [p100,n170,p300,n400]
    lines!(ax,b(),label=string(b))
    scatter!(ax,b(),label=string(b))
end
axislegend(ax,merge=true)
f

# ## fMRI
# default hrf TR is 1. Get to know all your favourite shapes!
##--
f = Figure()
plotConfig = (:peak=>1:3:10,
             :psUnder=>10:5:30,
             :amplitude=>2:5,
             :shift=>0:3:10,
             :peak_width => 0.1:0.5:1.5,
             :psUnder_width => 0.1:0.5:1.5,
             )

for (ix,pl) = enumerate(plotConfig)
    col = (ix-1)%3 +1
    row = Int(ceil(ix/3))
    
    ax = f[row,col] = Axis(f)
    cfg = collect(pl)
    for k = cfg[2]    
        lines!(ax,UnfoldSim.hrf(;TR=0.1,(cfg[1]=>k,)...),label=string(k))
    end
    
    axislegend(string(cfg[1]);merge=true,)
end
f

# ## Pupil
# We use the simplified PuRF from Hoeks & Levelt, 1993. Note that https://www.science.org/doi/10.1126/sciadv.abi9979 show some evidence in their supplementary material, that the convolution model is not fully applicable.
f = Figure()
plotConfig = (:n=>5:3:15,
             :tmax=>0.5:0.2:1.1,
             )

for (ix,pl) = enumerate(plotConfig)
    ax = f[1,ix] = Axis(f)
    cfg = collect(pl)
    for k = cfg[2]    
        lines!(ax,UnfoldSim.PuRF(;(cfg[1]=>k,)...),label=string(k))
    end
    
    axislegend(string(cfg[1]);merge=true,)
end
f
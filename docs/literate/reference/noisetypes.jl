using UnfoldSim
using CairoMakie
using DSP
using StableRNGs
import StatsBase.autocor
# ## What's the noise?
# There are several noise-types directly implemented. Here is a comparison

f = Figure()
ax_sig = f[1,1:2] = Axis(f;title="1.000 samples of noise")
ax_spec = f[2,1] = Axis(f;title="Welch Periodigram")
ax_auto = f[2,2] = Axis(f;title="Autocorrelogram (every 10th lag)")
for n = [PinkNoise RedNoise WhiteNoise NoNoise ExponentialNoise]

    ## generate
    noisevec = gen_noise(StableRNG(1),n(),10000)

    ## plot 1000 samples
    lines!(ax_sig,noisevec[1:1000];label=string(n))

    ## calc spectrum
    perio = welch_pgram(noisevec)

    ## plot spectrum
    lines!(ax_spec,freq(perio),log10.(power(perio)))
    
    lags = 0:10:500
    autocor_vec = autocor(noisevec,lags)
    lines!(ax_auto,lags,autocor_vec)
    
end
f[1:2,3] = Legend(f,ax_sig,"NoiseType")
f


# !!! hint
#    We recommed for smaller signals the `ExponentialNoise`, maybe with a removed DC offset or a HighPass filter. For long signals, this Noise requires lot's of memory though. maybe Pinknoise is a better choice then. 
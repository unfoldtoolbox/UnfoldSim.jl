using UnfoldSim
using CairoMakie
using DSP
using StableRNGs
# ## What's the noise?
# There are several noise-types directly implemented. Here is a comparison

f = Figure()
ax1 = f[1,1] = Axis(f)
ax2 = f[2,1] = Axis(f)
for n = [PinkNoise RedNoise WhiteNoise NoNoise]
    noisevec = gen_noise(StableRNG(1),n(),10000)
    lines!(ax1,noisevec[1:1000];label=string(n))
    perio = welch_pgram(noisevec)
    lines!(ax2,freq(perio),log10.(power(perio)))
    
end
f[1:2,2] = Legend(f,ax1,"NoiseType")

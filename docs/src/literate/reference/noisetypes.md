```@meta
EditURL = "<unknown>/docs/src/literate/reference/noisetypes.jl"
```

````@example noisetypes
using UnfoldSim
using CairoMakie
using DSP
using StableRNGs
import StatsBase.autocor
````

## What's the noise?
There are several noise-types directly implemented. Here is a comparison

````@example noisetypes
##---
````

---

````@example noisetypes
f = Figure()
ax_sig = f[1,1:2] = Axis(f;title="1.000 samples of noise")
ax_spec = f[2,1] = Axis(f;title="Welch Periodigram")
ax_auto = f[2,2] = Axis(f;title="Autocorrelogram (every 10th lag)")
for n = [PinkNoise RedNoise WhiteNoise NoNoise ExponentialNoise]
````

generate

````@example noisetypes
    noisevec = gen_noise(StableRNG(1),n(),10000)
````

plot 1000 samples

````@example noisetypes
    lines!(ax_sig,noisevec[1:1000];label=string(n))
````

calc spectrum

````@example noisetypes
    perio = welch_pgram(noisevec)
````

plot spectrum

````@example noisetypes
    lines!(ax_spec,freq(perio),log10.(power(perio)))

    lags = 1:10:500
    autocor_vec = autocor(noisevec,lags)
    lines!(ax_auto,lags,autocor_vec)

end
f[1:2,3] = Legend(f,ax1,"NoiseType")
f
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


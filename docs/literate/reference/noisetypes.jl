# # Overview: Noise types

# There are different types of noise signals which differ in their power spectra.
# If you are not familiar with different types/colors of noise yet, have a look at the [colors of noise Wikipedia page](https://en.wikipedia.org/wiki/Colors_of_noise).

# There are several noise types directly implemented in UnfoldSim.jl. Here is a comparison:


using UnfoldSim
using CairoMakie
using DSP
using StableRNGs
import StatsBase.autocor

f = Figure()
ax_sig =
    f[1, 1:3] =
        Axis(f; title = "1.000 samples of noise", xlabel = "Time", ylabel = "Amplitude")
ax_spec =
    f[2, 1:2] = Axis(
        f;
        title = "Welch Periodogram",
        xlabel = "Normalized frequency",
        ylabel = "log(Power)",
    )
ax_auto =
    f[2, 3:4] = Axis(
        f;
        title = "Autocorrelogram (every 10th lag)",
        xlabel = "Lag",
        ylabel = "Autocorrelation",
    )
for n in [PinkNoise RedNoise WhiteNoise NoNoise ExponentialNoise]

    ## generate
    noisevec = simulate_noise(StableRNG(1), n(), 10000)

    ## plot 1000 samples
    lines!(ax_sig, noisevec[1:1000]; label = string(n))

    ## calc spectrum
    perio = welch_pgram(noisevec)

    ## plot spectrum
    lines!(ax_spec, freq(perio), log10.(power(perio)))

    lags = 0:10:500
    autocor_vec = autocor(noisevec, lags)
    lines!(ax_auto, lags, autocor_vec)

end
f[1, 4] = Legend(f, ax_sig, "Noise type", tellheight = true)
f


# !!! hint
#     From a theoretical point, `ExponentialNoise` seems to be the best fit for the AR spectrum of EEG signals. PinkNoise seems to be the most common choice in research papers.

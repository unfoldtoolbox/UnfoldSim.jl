#----------------
# Types
#-------------


"""
    PinkNoise <: AbstractNoise
Generate Pink Noise using the SignalAnalysis.jl implementation.

`noiselevel` is used to scale the noise
"""
@with_kw struct PinkNoise <: AbstractNoise
    noiselevel = 1
    func = SignalAnalysis.PinkGaussian
end


"""
    RedNoise <: AbstractNoise
Generate Red Noise using the SignalAnalysis.jl implementation.

`noiselevel` is used to scale the noise
"""
@with_kw struct RedNoise <: AbstractNoise
    noiselevel = 1
    func = SignalAnalysis.RedGaussian
end


"""
    WhiteNoise <: WhiteNoise
    noiselevel = 1
    imfilter = 0
Generate White Noise using `randn` - thus Gaussian noise.
`noiselevel` is used to scale the noise

Using `imfilter` > 0 it is possible to smooth the noise using Image.imfilter.
"""
@with_kw struct WhiteNoise <: AbstractNoise
    noiselevel = 1
    imfilter = 0
end

"""
    RealisticNoise <: AbstractNoise
Not implemented - planned to use Artifacts.jl to provide real EEG data to add.
"""
@with_kw struct RealisticNoise <: AbstractNoise
    noiselevel = 1
end

"""
    NoNoise <: AbstractNoise

Return zeros instead of noise.
"""
struct NoNoise <: AbstractNoise end

"""
    AutoRegressiveNoise <: AbstractNoise
Not implemented
"""
struct AutoRegressiveNoise <: AbstractNoise end



""" 
    ExponentialNoise <: AbstractNoise

Noise with exponential decay in AR spectrum. Implements the algorithm (3) from Markus Deserno 2002: https://www.cmu.edu/biolphys/deserno/pdf/corr_gaussian_random.pdf
"How to generate exponentially correlated Gaussian random numbers"

The noise has std of 1 and mean of 0 (over many samples)

The factor "τ" defines the decay over samples. 

`noiselevel` is used to scale the noise

"""
@with_kw struct ExponentialNoise <: AbstractNoise
    noiselevel = 1
    τ = 1000
end


"""
    simulate_noise(rng, t::Union{PinkNoise,RedNoise}, n::Int)

Generate Pink or Red Noise using the `SignalAnalysis.jl` implementation.
"""
function simulate_noise(rng, t::Union{PinkNoise,RedNoise}, n::Int)
    return t.noiselevel .* rand(rng, t.func(n, 1.0))
end

"""
    simulate_noise(rng, t::NoNoise, n::Int)
Return zeros instead of noise.
"""
function simulate_noise(rng, t::NoNoise, n::Int)
    return zeros(n)
end


function simulate_noise(rng, t::WhiteNoise, n::Int)
    noisevector = randn(rng, n)
    if !isnothing(t.imfilter)
        noisevector = imfilter(noisevector, Kernel.gaussian((t.imfilter,)))
    end
    return t.noiselevel .* noisevector
end



function simulate_noise(rng, t::RealisticNoise, n::Int)
    error("not implemented")
    return 0
end


function simulate_noise(rng, t::ExponentialNoise, n::Int)

    τ = t.τ


    @assert τ > 0
    f = exp(-1 / τ)
    #s = 0:n-1
    g = randn(rng, n)
    r = similar(g)
    r[1] = g[1]
    for n = 1:length(g)-1
        r[n+1] = f * r[n] + sqrt(1 - f^2) * g[n+1]
    end




    return t.noiselevel .* r
end

"""
    add_noise!(rng, noisetype::AbstractNoise, signal)

Generate and add noise to a data matrix.
Assumes that the signal can be linearized, that is, that the noise is stationary
"""
function add_noise!(rng, noisetype::AbstractNoise, signal)

    # generate noise
    noise = simulate_noise(deepcopy(rng), noisetype, length(signal))

    noise = reshape(noise, size(signal))

    # add noise to data
    signal .+= noise

end



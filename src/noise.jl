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

Noise with exponential decay in AR spectrum.

`noiselevel` is used to scale the noise

!!! warning
    With the current implementation we try to get exponential decay over the whole autoregressive (AR) spectrum, which is N-Samples (the total number of samples) long. This involves the inversion of a cholesky matrix of size NxN matrix, which will need lots of RAM for non-trivial problems.
"""
@with_kw struct ExponentialNoise <: AbstractNoise
    noiselevel = 1
    ν = 1.5 # exponential factor of AR decay "nu"
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

    function exponential_correlation(x; nu = 1, length_ratio = 1)
        # Author: Jaromil Frossard
        # generate exponential function
        R = length(x) * length_ratio
        return exp.(-3 * (x / R) .^ nu)
    end

    Σ = Symmetric(Circulant(exponential_correlation([0:1:(n-1);], nu = t.ν)), :L)

    # cholesky(Σ) is n x n diagonal, lots of RAM :S
    return t.noiselevel .* 10 .* (randn(rng, n)'*cholesky(Σ).U)[1, :]
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

#----------------
# Types
#-------------

@with_kw struct PinkNoise <: AbstractNoise
    noiselevel = 1
    func = SignalAnalysis.PinkGaussian
end


@with_kw struct RedNoise <: AbstractNoise
    noiselevel = 1
    func = SignalAnalysis.RedGaussian
end


@with_kw struct WhiteNoise <: AbstractNoise
    noiselevel = 1
    imfilter = 0
end


@with_kw struct RealisticNoise <: AbstractNoise
    noiselevel = 1
end

struct NoNoise <: AbstractNoise end

struct AutoRegressiveNoise <: AbstractNoise end



""" 
    Noise with exponential decay in AR spectrum
    !!! warning
        Current implementation: cholesky of NxN matrix needs to be calculated, might need lot's of RAM

"""

@with_kw struct ExponentialNoise <: AbstractNoise
    noiselevel = 1
    ν = 1.5 # exponential factor of AR decay "nu"
end


"""
    simulate_noise(t::Union{PinkNoise, RedNoise}, n::Int)

Generate noise of a given type t and length n
"""
function simulate_noise(rng, t::Union{PinkNoise,RedNoise}, n::Int)
    return t.noiselevel .* rand(rng, t.func(n, 1.0))
end

function simulate_noise(rng, t::NoNoise, n::Int)
    return zeros(n)
end

"""
    simulate_noise(t::WhiteNoise, n::Int)

Generate noise of a given type t and length n
"""
function simulate_noise(rng, t::WhiteNoise, n::Int)
    noisevector = randn(rng, n)
    if !isnothing(t.imfilter)
        noisevector = imfilter(noisevector, Kernel.gaussian((t.imfilter,)))
    end
    return t.noiselevel .* noisevector
end


"""
    simulate_noise(t::RealisticNoise, n::Int)

Generate noise of a given type t and length n
"""
function simulate_noise(rng, t::RealisticNoise, n::Int)
    error("not implemented")
    return 0
end


function simulate_noise(rng, t::ExponentialNoise, n::Int)

    function exponentialCorrelation(x; nu = 1, length_ratio = 1)
        # Author: Jaromil Frossard
        # generate exponential function
        R = length(x) * length_ratio
        return exp.(-3 * (x / R) .^ nu)
    end

    Σ = Symmetric(Circulant(exponentialCorrelation([0:1:(n-1);], nu = t.ν)), :L)

    # cholesky(Σ) is n x n diagonal, lots of RAM :S
    return t.noiselevel .* 10 .* (randn(rng, n)'*cholesky(Σ).U)[1, :]
end

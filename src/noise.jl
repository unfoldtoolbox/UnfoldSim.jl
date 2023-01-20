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
    gen_noise(t::Union{PinkNoise, RedNoise}, n::Int)

Generate noise of a given type t and length n
"""
function gen_noise(rng, t::Union{PinkNoise, RedNoise}, n::Int)
    return t.noiselevel .*rand(rng, t.func(n, 1.0))
end

function gen_noise(rng,t::NoNoise,n::Int)
    return zeros(n)
end

"""
    gen_noise(t::WhiteNoise, n::Int)

Generate noise of a given type t and length n
"""
function gen_noise(rng, t::WhiteNoise, n::Int)
    noisevector = randn(rng, n)
    if !isnothing(t.imfilter)
        noisevector = imfilter(noisevector, Kernel.gaussian((t.imfilter,)))
    end
    return t.noiselevel .*noisevector
end


"""
    gen_noise(t::RealisticNoise, n::Int)

Generate noise of a given type t and length n
"""
function gen_noise(rng, t::RealisticNoise, n::Int)
    error("not implemented")
    return 0
end


function gen_noise(rng,t::ExponentialNoise, n::Int)
    # helper functions from Jaromil Frossard
    function circulant(x)
        # Author: Jaromil Frossard
        # returns a symmetric matrix where X was circ-shifted.
        lx = length(x)
        ids = [1:1:(lx-1);]
        a = Array{Float64,2}(undef, lx, lx)
        for i = 1:length(x)
            if i == 1
                a[i, :] = x
            else
                a[i, :] = vcat(x[i], a[i-1, ids])
            end
        end
        return Symmetric(a)
    end
    
    
    function exponentialCorrelation(x; nu = 1, length_ratio = 1)
        # Author: Jaromil Frossard
        # generate exponential function
        R = length(x) * length_ratio
        return exp.(-3 * (x / R) .^ nu)
    end
    
    Σ = circulant(exponentialCorrelation([0:1:(n-1);], nu = t.ν))
    Ut = LinearAlgebra.cholesky(Σ).U'
    return t.noiselevel .* 10 .* (randn(rng, n)'*Ut')[1, :]    
end
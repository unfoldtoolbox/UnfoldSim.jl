abstract type Noise end


@with_kw struct PinkNoise <: Noise
    func = SignalAnalysis.PinkGaussian
    noiselevel = 1
end


@with_kw struct RedNoise <: Noise
    func = SignalAnalysis.RedGaussian
    noiselevel = 1
end


@with_kw struct WhiteNoise <: Noise
    func = randn
    noiselevel = 1
    imfilter = 5
end


@with_kw struct RealisticNoise <: Noise 
    func = error
    noiselevel = 1
end


struct NoNoise <: Noise end


"""
    gen_noise(t::Union{PinkNoise, RedNoise}, n::Int)

Generate noise of a given type t and length n
"""
function gen_noise(rng, t::Union{PinkNoise, RedNoise}, n::Int)
    return t.noiselevel .* rand(rng, t.func(n, 1.0))
end

function gen_noise(rng,t::NoNoise,n::Int)
    return zeros(n)
end

"""
    gen_noise(t::WhiteNoise, n::Int)

Generate noise of a given type t and length n
"""
function gen_noise(rng, t::WhiteNoise, n::Int)
    noisevector = t.noiseLevel .* randn(rng, n)
    if !isnothing(t.imfilter)
        noisevector = imfilter(noisevector, Kernel.gaussian((t.imfilter,)))
    end
    return noisevector
end


"""
    gen_noise(t::RealisticNoise, n::Int)

Generate noise of a given type t and length n
"""
function gen_noise(rng, t::RealisticNoise, n::Int)
    return 0
end

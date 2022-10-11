abstract type Noise end


@with_kw struct PinkNoise <: Noise
    func = SignalAnalysis.PinkGaussian
end


@with_kw struct RedNoise <: Noise
    func = SignalAnalysis.RedGaussian
end


struct WhiteNoise <: Noise end


struct RealisticNoise <: Noise end


"""
    gen_noise(t::Union{PinkNoise, RedNoise}, n::Int)

Generate noise of a given type t and length n
"""
function gen_noise(rng::MersenneTwister, t::Union{PinkNoise, RedNoise}, n::Int)
    return rand(rng, noise.func(n, 1.0))
end


"""
    gen_noise(t::WhiteNoise, n::Int)

Generate noise of a given type t and length n
"""
function gen_noise(rng::MersenneTwister, t::WhiteNoise, n::Int)
    return randn(rng, n)
end


"""
    gen_noise(t::RealisticNoise, n::Int)

Generate noise of a given type t and length n
"""
function gen_noise(rng::MersenneTwister, t::RealisticNoise, n::Int)
    return 0
end

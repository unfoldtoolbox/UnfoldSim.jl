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
    noisevector = t.noiselevel .* randn(rng, n)
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
    error("not implemented")
    return 0
end

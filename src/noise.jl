#----------------
# Noise types
#----------------

"""
    PinkNoise <: AbstractNoise

A type for generating Pink Noise using the `SignalAnalysis.jl` implementation.

Tip: To manually create noise samples use the [`simulate_noise`](@ref) function.

# Fields
- `noiselevel = 1` (optional): Factor that is used to scale the noise.
- `func = SignalAnalysis.PinkGaussian` (optional): Function that is used to create the noise samples.
    This field is for internal use and should not typically be modified directly by the user.
    Changes to this field may result in unexpected behavior.

# Examples
```julia-repl
julia> noise = PinkNoise(noiselevel = 3)
PinkNoise
  noiselevel: Int64 3
  func: UnionAll

julia> using StableRNGs

julia> simulate_noise(StableRNG(1), noise, 5)
5-element Vector{Float64}:
 2.578878369756878
 3.4972108606501786
 2.878568584946028
 2.2725654770788384
 3.5291669151888683
```

See also [`RedNoise`](@ref), [`WhiteNoise`](@ref), [`ExponentialNoise`](@ref), [`NoNoise`](@ref).
"""
@with_kw struct PinkNoise <: AbstractNoise
    noiselevel = 1
    func = SignalAnalysis.PinkGaussian
end


"""
    RedNoise <: AbstractNoise

A type for generating Red Noise using the `SignalAnalysis.jl` implementation.

Tip: To manually create noise samples use the [`simulate_noise`](@ref) function.

# Fields
- `noiselevel = 1` (optional): Factor that is used to scale the noise.
- `func = SignalAnalysis.RedGaussian` (optional): Function that is used to create the noise samples.
    This field is for internal use and should not typically be modified directly by the user.
    Changes to this field may result in unexpected behavior.

# Examples
```julia-repl
julia> noise = RedNoise(noiselevel = 2)
RedNoise
  noiselevel: Int64 2
  func: UnionAll

julia> using StableRNGs

julia> simulate_noise(StableRNG(2), noise, 3)
3-element Vector{Float64}:
 -0.34153942884005967
 -0.4651387715669636
 -0.4951538876376382
```

See also [`PinkNoise`](@ref), [`WhiteNoise`](@ref), [`ExponentialNoise`](@ref), [`NoNoise`](@ref).
"""
@with_kw struct RedNoise <: AbstractNoise
    noiselevel = 1
    func = SignalAnalysis.RedGaussian
end


"""
    WhiteNoise <: AbstractNoise

A type for generating White Noise using `randn` - thus Gaussian noise.

Using `imfilter` > 0 it is possible to smooth the noise using Image.imfilter.
Tip: To manually create noise samples use the [`simulate_noise`](@ref) function.

# Fields
- `noiselevel = 1` (optional): Factor that is used to scale the noise.
- `imfilter = 0` (optional):

# Examples
```julia-repl
julia> noise = WhiteNoise()
WhiteNoise
  noiselevel: Int64 1
  imfilter: Int64 0

julia> using StableRNGs

julia> simulate_noise(StableRNG(1), noise, 3)
3-element Vector{Float64}:
 -0.5325200748641231
  0.098465514284785
  0.7528865221245234
```

See also [`PinkNoise`](@ref), [`RedNoise`](@ref), [`ExponentialNoise`](@ref), [`NoNoise`](@ref).
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

A type for simulations without noise; return zeros instead of noise.

Tip: To manually create noise samples use the [`simulate_noise`](@ref) function.

# Examples
```julia-repl
julia> noise = NoNoise()
NoNoise()

julia> using StableRNGs

julia> simulate_noise(StableRNG(1), noise, 3)
3-element Vector{Float64}:
 0.0
 0.0
 0.0
```
See also [`PinkNoise`](@ref), [`RedNoise`](@ref), [`ExponentialNoise`](@ref), [`WhiteNoise`](@ref).
"""
struct NoNoise <: AbstractNoise end

"""
    AutoRegressiveNoise <: AbstractNoise

Not implemented.
"""
struct AutoRegressiveNoise <: AbstractNoise end

""" 
    ExponentialNoise <: AbstractNoise

Type for generating noise with exponential decay in AR spectrum.

Tip: To manually create noise samples use the [`simulate_noise`](@ref) function.

!!! warning
    With the current implementation we try to get exponential decay over the whole autoregressive (AR) spectrum, which is N samples (the total number of samples in the signal) long. This involves the inversion of a Cholesky matrix of size NxN matrix, which will need lots of RAM for non-trivial problems.

# Fields
- `noiselevel = 1` (optional): Factor that is used to scale the noise.
- `ν = 1.5 ` (optional): Exponential factor of AR decay "nu".

# Examples
```julia-repl
julia> noise = ExponentialNoise()
ExponentialNoise
  noiselevel: Int64 1
  ν: Float64 1.5

julia> using StableRNGs

julia> simulate_noise(StableRNG(1), noise, 5)
5-element Vector{Float64}:
  -5.325200748641231
  -3.437402125380177
   2.7852625669058884
  -1.5381022393382109
 -14.818799857226612
```
See also [`PinkNoise`](@ref), [`RedNoise`](@ref), [`NoNoise`](@ref), [`WhiteNoise`](@ref).
"""
@with_kw struct ExponentialNoise <: AbstractNoise
    noiselevel = 1
    ν = 1.5 # exponential factor of AR decay "nu"
end

#-----------------------------
# Noise simulation functions
#-----------------------------

"""
    simulate_noise(rng, t::AbstractNoise, n::Int)

Generate noise samples of the given type `t`.

# Arguments
- `rng::AbstractRNG`: Random number generator (RNG) to make the process reproducible.
- `t::AbstractNoise`: Instance of a noise type e.g. `PinkNoise()`.
- `n::Int`: The number of noise samples that should be generated.

# Returns
- `Vector{Float64}`: Vector of length `n` containing the noise samples.

# Examples
```julia-repl
# Here we use White Noise as an example but it works in the same way for the other noise types.
julia> noise = WhiteNoise()
WhiteNoise
  noiselevel: Int64 1
  imfilter: Int64 0

julia> using StableRNGs

julia> simulate_noise(StableRNG(1), noise, 3)
3-element Vector{Float64}:
 -0.5325200748641231
  0.098465514284785
  0.7528865221245234
```
"""
function simulate_noise end


function simulate_noise(rng, t::Union{PinkNoise,RedNoise}, n::Int)
    return t.noiselevel .* rand(rng, t.func(n, 1.0))
end


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

Generate and add noise to a signal.

Assumes that the signal can be linearized, that is, that the noise is stationary.
"""
function add_noise!(rng, noisetype::AbstractNoise, signal)

    # generate noise
    noise = simulate_noise(deepcopy(rng), noisetype, length(signal))

    noise = reshape(noise, size(signal))

    # add noise to data
    signal .+= noise

end

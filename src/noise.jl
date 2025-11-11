#----------------
# Noise types
#----------------

"""
    PinkNoise <: AbstractNoise

A type for generating Pink Noise using the `SignalAnalysis.jl` implementation.

The noise values are sampled from a standard normal distribution 𝒩(μ=0, σ=1).
That means that ~95% of the values are between ~-2 and ~2 (with `noiselevel = 1`).\\
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

The noise values are sampled from a standard normal distribution 𝒩(μ=0, σ=1).
That means that ~95% of the values are between ~-2 and ~2 (with `noiselevel = 1`).\\
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

The noise values are sampled from a standard normal distribution 𝒩(μ=0, σ=1).
That means that ~95% of the values are between ~-2 and ~2 (with `noiselevel = 1`).\\
Tip: To manually create noise samples use the [`simulate_noise`](@ref) function.

# Fields
- `noiselevel = 1` (optional): Factor that is used to scale the noise.
- `imfilter = nothing` (optional): Use `imfilter > 0` to smooth the noise using `Image.imfilter` with a Gaussian kernel with `σ = imfilter`.

# Examples
```julia-repl
julia> noise = WhiteNoise()
WhiteNoise
  noiselevel: Int64 1
  imfilter: Nothing nothing

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
    imfilter = nothing
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

Noise with exponential decay in AR spectrum. Implements the algorithm (3) from Markus Deserno 2002: https://www.cmu.edu/biolphys/deserno/pdf/corr_gaussian_random.pdf
"How to generate exponentially correlated Gaussian random numbers"

The noise has std of 1 and mean of 0 (over many samples)

The factor "τ" defines the decay over samples. 

# Fields
- `noiselevel = 1` (optional): Factor that is used to scale the noise.
- `τ`: Exponential factor of AR decay "tau" in samples. Must be > 0. This factor should depend on your sampling rate, therefore no default is provided.

# Examples
```julia-repl
julia> noise = ExponentialNoise(τ = 10)
julia> using StableRNGs
julia> simulate_noise(StableRNG(1), noise, 5)
5-element Vector{Float64}:
 -0.5325200748641231
 -0.43992168173929114
 -0.07751069370021019
 -0.42927675219497446
 -1.2458281425035913
```
See also [`PinkNoise`](@ref), [`RedNoise`](@ref), [`NoNoise`](@ref), [`WhiteNoise`](@ref).
"""
@with_kw struct ExponentialNoise <: AbstractNoise
    noiselevel = 1
    τ::Any
end

#-----------------------------
# Noise simulation functions
#-----------------------------

"""
    simulate_noise(rng, t::AbstractNoise, n::Int, [Simulation::Simulation])
    simulate_noise(rng, t::AbstractNoise, n::Tuple, [Simulation::Simulation])

Generate noise samples of the given type `t`.

For details, see the documentation of the individual noise types.
Use `subtypes(AbstractNoise)` for a list of the implemented noise types.

# Arguments
- `rng::AbstractRNG`: Random number generator (RNG) to make the process reproducible.
- `t::AbstractNoise`: Instance of a noise type e.g. `PinkNoise()`.
- `n::Int`: The number of noise samples that should be generated. If a tuple is provided, the `prod(n)` samples are generated. Future usage might generate multi-dimensional noise more directly.
- `simulation::Simulation` (optional): Currently not used, but provided for future extensions where the noise generation might depend on the simulation design.
# Returns
- `Vector`: Vector of length `n` containing the noise samples.

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

    τ = t.τ


    @assert τ > 0
    f = exp(-1 / τ)
    #s = 0:n-1
    g = randn(rng, n)
    r = similar(g)
    r[1] = g[1]
    for n = 1:(length(g)-1)
        r[n+1] = f * r[n] + sqrt(1 - f^2) * g[n+1]
    end



simulate_noise(rng::AbstractRNG, t::AbstractNoise, signal, Simulation::Simulation) =
    simulate_noise(rng, t, prod(signal))


    return t.noiselevel .* r
end

"""
    add_noise!(rng::AbstractRNG, noisetype::AbstractNoise, signal[, simulation::Simulation])

Add simulated noise to `signal` in-place.

Arguments
- `rng`: random number generator.
- `noisetype`: noise descriptor (e.g. `WhiteNoise`, `PinkNoise`, `ExponentialNoise`, `NoNoise`).
- `signal`: array-like container to mutate; generated noise is reshaped to `size(signal)`.
- `simulation` (optional): reserved for future use.

Returns
- `nothing` (mutates `signal`).

Notes
- Internally calls `simulate_noise(deepcopy(rng), ...)`.
"""

#    add_noise!(rng, noisetype, signal)
function add_noise!(
    rng::AbstractRNG,
    noisetype::AbstractNoise,
    signal,
    simulation::Simulation,
)

    # generate noise
    noise = simulate_noise(deepcopy(rng), noisetype, size(signal), simulation)

    noise = reshape(noise, size(signal))

    # add noise to data
    signal .+= noise

end

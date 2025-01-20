#----------------
# Types
#----------------

"""
    UniformOnset <: AbstractOnset

Provide a Uniform Distribution for the inter-event distances (in samples).

Tip: To manually generate inter-event distance samples use the [`simulate_interonset_distances`](@ref) function.

# Fields
- `width = 50` (optional): Width of the uniform distribution (=> the "jitter"). Since the lower bound is 0, `width` is also the upper bound.
- `offset = 0` (optional): The minimal distance between events. The maximal distance is `offset + width`.

# Examples
```julia-repl
julia> onset_distribution = UniformOnset(width = 25, offset = 5)
UniformOnset
  width: Int64 25
  offset: Int64 5
```

See also [`LogNormalOnset`](@ref), [`NoOnset`](@ref).
"""
@with_kw struct UniformOnset <: AbstractOnset
    width = 50 # how many samples jitter?
    offset = 0 # minimal offset?
end
"""
    LogNormalOnset <: AbstractOnset

Log-normal inter-event distances (in samples) using the `Distributions.jl` truncated LogNormal distribution ([code and mathematical reference](https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.LogNormal)).

Be careful with large `μ` and `σ` values, as they are on logscale. σ>8 can quickly give you out-of-memory sized signals! \\
Tip: To manually generate inter-event distance samples use the [`simulate_interonset_distances`](@ref) function.

# Fields
- `μ`: The mean of the log-transformed variable (the log-normal random variable's logarithm follows a normal distribution).
- `σ`: The standard deviation of the log-transformed variable.
- `offset = 0` (optional): The minimal distance between events.
- `truncate_upper = nothing` (optional): Upper limit (in samples) at which the distribution is truncated.

# Examples
```julia-repl
julia> onset_distribution = LogNormalOnset(3, 0.25, 10, 25)
LogNormalOnset
  μ: Int64 3
  σ: Float64 0.25
  offset: Int64 10
  truncate_upper: Int64 25
```

See also [`UniformOnset`](@ref), [`NoOnset`](@ref).
"""
@with_kw struct LogNormalOnset <: AbstractOnset
    μ::Any  # mean
    σ::Any  # variance
    offset = 0 # additional offset
    truncate_upper = nothing # truncate at some sample?
end

"""
    NoOnset <: AbstractOnset

For cases where the user wants to simulate epoched data without any overlap between consecutive events.

# Examples
```julia-repl
julia> onset_distribution = NoOnset()
NoOnset()
```

See also [`UniformOnset`](@ref), [`LogNormalOnset`](@ref).
"""
struct NoOnset <: AbstractOnset end


#-----------------------------
# Onset simulation functions
#-----------------------------

"""
    simulate_interonset_distances(rng, onset::AbstractOnset, design::AbstractDesign)

Generate the inter-onset distance vector by sampling from the respective distribution (in samples).

# Arguments
- `rng`: Random number generator (RNG) to make the process reproducible.
- `onset::AbstractOnset`: Inter-onset distance distribution to sample from.
- `design::AbstractDesign`: Experimental design with conditions and covariates.

# Returns
- `Vector{Int64}`: Inter-onset distances in samples. Note that these are distances between onsets and no latencies.

# Examples
```julia-repl
# Create an experimental design
julia> design_single = SingleSubjectDesign(;
           conditions = Dict(
               :stimulus_type => ["natural", "artificial"],
               :contrast_level => range(0, 1, length = 3),
           ),
       );

# Create an inter-onset distance distribution
julia> onset_distribution = LogNormalOnset(3, 0.5, 5, nothing);

julia> using StableRNGs

julia> simulate_interonset_distances(StableRNG(1), onset_distribution, design_single)
6-element Vector{Int64}:
 20
 26
 34
 18
 12
 23
```

See also [`simulate_onsets`](@ref).
"""
function simulate_interonset_distances end


function simulate_interonset_distances(rng, onset::UniformOnset, design::AbstractDesign)
    return Int.(
        round.(rand(deepcopy(rng), onset.offset:(onset.offset+onset.width), size(design)))
    )
end

function simulate_interonset_distances(rng, onset::LogNormalOnset, design::AbstractDesign)
    s = size(design)
    fun = LogNormal(onset.μ, onset.σ)
    if !isnothing(onset.truncate_upper)
        fun = truncated(fun; upper = onset.truncate_upper)
    end
    return Int.(round.(onset.offset .+ rand(deepcopy(rng), fun, s)))
end



"""
    simulate_onsets(rng, onset::AbstractOnset, simulation::Simulation)

Call `simulate_interonset_distances` to generate distances between events and then add them up to generate the actual latencies in samples.

Please note that this function is mainly for internal use in the context of `simulate` function calls. \n
Also note that the accumulation of onsets starts at 1 to avoid indexing problems in the case that the first sampled onset is 0.

# Arguments
- `rng`: Random number generator (RNG) to make the process reproducible.
- `onset::AbstractOnset`: Inter-onset distance distribution which is passed to `simulate_interonset_distances`.
- `simulation::Simulation`: Simulation object which contains design, component(s), inter-onset distance distribution and noise.

# Returns

# Examples
```julia-repl
# Create Simulation object
julia> design_single = SingleSubjectDesign(;
           conditions = Dict(
               :stimulus_type => ["natural", "artificial"],
               :contrast_level => range(0, 1, length = 3),
           ),
       );

julia> p1_component = LinearModelComponent(; basis = p100(), formula = @formula(0 ~ 1), β = [5]);

julia> simulation = Simulation(design_single, p1_component, UniformOnset(), NoNoise());

julia> using StableRNGs

# Simulate onsets for this simulation
julia> simulate_onsets(StableRNG(1), simulation.onset, simulation)
6-element Vector{Int64}:
  20
  70
  97
 110
 150
 182
```

See also [`simulate_interonset_distances`](@ref).
"""
function simulate_onsets(rng, onset::AbstractOnset, simulation::Simulation)

    # sample different onsets
    onsets = simulate_interonset_distances(rng, onset, simulation.design)

    if maximum(onsets) > 10000
        @warn "Maximum of inter-event distances was $(maximum(onsets)) - are you sure this is what you want?"
    end
    # accumulate them
    onsets_accum = accumulate(+, onsets, dims = 1, init = 1)

    return onsets_accum
end

#----------------
# Types
#---------------

"""
    struct UniformOnset <: AbstractOnset
Provide a Uniform Distribution of the inter-event-distances.
`width`  is the width of the uniform distribution (=> the jitter). Since the lower bound is 0, `width` is also the upper bound.
`offset` is the minimal distance. The maximal distance is `offset + width`.
"""
@with_kw struct UniformOnset <: AbstractOnset
    width = 50 # how many samples jitter?
    offset = 0 # minimal offset?
end
"""
    @with_kw struct LogNormalOnset <: AbstractOnset
Log-normal inter-event distances using the `Distributions.jl` truncated LogNormal distribution.

Be careful with large `μ` and `σ` values, as they are on logscale. σ>8 can quickly give you out-of-memory sized signals!
"""
@with_kw struct LogNormalOnset <: AbstractOnset
    μ::Any  # mean
    σ::Any  # variance
    offset = 0 # additional offset
    truncate_upper = nothing # truncate at some sample?
end

"""
    struct NoOnset <: AbstractOnset end
In the case that the user directly wants no overlap to be simulated (=> epoched data).
"""
struct NoOnset <: AbstractOnset end

"""
    simulate_interonset_distances(rng, onset::UniformOnset, design::AbstractDesign)
    simulate_interonset_distances(rng, onset::LogNormalOnset, design::AbstractDesign)
Generate the inter-event-onset vector in samples (returns Int).
"""

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
# main call from `simulation`
"""
function simulate_onsets(rng, onset::AbstractOnset, simulation::Simulation)

    # sample different onsets
    onsets = simulate_interonset_distances(rng, onset, simulation.design)

    # accumulate them
    onsets_accum = accumulate(+, onsets, dims = 1, init = 1)

    return onsets_accum
end

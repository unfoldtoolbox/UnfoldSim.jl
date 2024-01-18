#----------------
# Types
#---------------

@with_kw struct UniformOnset <: AbstractOnset
    width = 50 # how many samples jitter?
    offset = 0 # minimal offset?
end
@with_kw struct LogNormalOnset <: AbstractOnset
    μ::Any  # mean
    σ::Any  # variance
    offset = 0 # additional offset
    truncate_upper = nothing # truncate at some sample?
end

# In the case that the user directly wants the erps/epoched data (no overlap) 
struct NoOnset <: AbstractOnset end

#-------------
function rand_onsets(rng, onset::UniformOnset, design::AbstractDesign)
    return Int.(
        round.(rand(deepcopy(rng), onset.offset:(onset.offset+onset.width), size(design)))
    )
end

function rand_onsets(rng, onset::LogNormalOnset, design::AbstractDesign)
    s = size(design)
    fun = LogNormal(onset.μ, onset.σ)
    if !isnothing(onset.truncate_upper)
        fun = truncated(fun; upper = onset.truncate_upper)
    end
    return Int.(round.(onset.offset .+ rand(deepcopy(rng), fun, s)))
end


# main call from `simulation`
function generate(rng, onset::AbstractOnset, simulation::Simulation)

	# sample different onsets
	onsets = rand_onsets(rng, onset, simulation.design)

    # accumulate them
    onsets_accum = accumulate(+, onsets, dims = 1)

    return onsets_accum
end

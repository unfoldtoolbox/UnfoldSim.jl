#----------------
# Types
#---------------

@with_kw struct UniformOnset<:AbstractOnset
    width=50 # how many samples jitter?
    offset=0 # minimal offset?
end
@with_kw struct LogNormalOnset<:AbstractOnset
    μ  # mean
    σ  # variance
    offset = 0 # additional offset
    truncate_upper = nothing # truncate at some sample?
end

#-------------
function rand_onsets(rng,onset::UniformOnset,design::AbstractDesign)
    return rand(deepcopy(rng), onset.offset:(onset.offset + onset.width), size(design))
end

function rand_onsets(rng,onset::LogNormalOnset,design::AbstractDesign)
    size = size(design)
    fun = LogNormal(onset.μ,onset.σ)
    if !isnothing(onset.truncate_upper)
        fun = truncated(fun;upper=onset.truncate_upper)
    end
    return onset.offset .+ rand(deepcopy(rng), fun, size)
end


# main call from `simulation`
function generate(rng,onset::AbstractOnset,simulation::Simulation)
    
	# sample different onsets
	onsets = rand_onsets(rng,onset,simulation.design)
	
    # accumulate them
	onsets_accum = accumulate(+, onsets, dims=1)
	return onsets_accum
end


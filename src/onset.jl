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
function rand_onsets(rng,onset::UniformOnset,n_item,n_subj,simulation)
    return rand(deepcopy(rng), onset.offset:(onset.offset + onset.width), (n_item, n_subj))
end

function rand_onsets(rng,onset::LogNormalOnset,n_item,n_subj,simulation)
    fun = LogNormal(onset.μ,onset.σ)
    if !isnothing(onset.truncate_upper)
        fun = truncated(fun;upper=onset.truncate_upper)
    end
    return onset.offset .+ rand(deepcopy(rng), fun, (n_item, n_subj))
end


# main call from `simulation`
function gen_onsets(rng,simulation::Simulation)
    (;n_item, n_subj) = simulation.design

	# sample different onsets
	onsets = rand_onsets(rng,simulation.onset,n_item,n_subj,simulation)
	
    # accumulate them
	onsets_accum = accumulate(+, onsets, dims=1)
	return onsets_accum
end


abstract type Onset end

@with_kw struct UniformOnset<:Onset
    width=50 # how many samples jitter?
    offset=0 # minimal offset?

end


function gen_onsets(rng,simulation::Simulation, onset::UniformOnset)
    (;n_item, n_subj) = simulation.design
    
	
	# sample different onsets
	onsets = rand(deepcopy(rng), onset.offset:(onset.offset + onset.width), (n_item, n_subj))
	
    # accumulate them
	onsets_accum = accumulate(+, onsets, dims=1)
	return onsets_accum
end


#@with_kw struct LogNormalOnset<:Onset

#end
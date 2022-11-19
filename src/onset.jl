abstract type Onset end

@with_kw struct UniformOnset<:Onset
    width=50 # how many samples jitter?
    offset=0 # minimal offset?

end


function gen_onsets(rng,simulation::Simulation, onset::UniformOnset)
    (;n_item, n_subj) = simulation.design
    
	
	# min + max fixationlen
	min_fixationlen = onset.offset + maxlength(simulation.components) - onset.width÷2
	max_fixationlen = onset.offset + maxlength(simulation.components) + onset.width÷2

	# sample different fixationlens & combine them with epoch lens
	onsets = rand(deepcopy(rng), min_fixationlen:max_fixationlen, (n_item, n_subj))
	
	onsets_accum = accumulate(+, onsets, dims=1)
	return onsets_accum
end


#@with_kw struct LogNormalOnset<:Onset

#end
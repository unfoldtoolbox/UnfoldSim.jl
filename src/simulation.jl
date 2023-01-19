
# helper to move input ::Component to ::Vector{Component}
Simulation(design::AbstractDesign,component::AbstractComponent,onset::AbstractOnset,noisetype::AbstractNoise) = Simulation(design,[component],onset,noisetype)

"""
Simulate eeg data given a simulation design, effect sizes and variances
"""
simulate(rng,design, signal,  onset, noise;kwargs...) = simulate(rng,Simulation(design, signal,  onset, noise);kwargs...)
function simulate(rng, simulation::Simulation;return_epoched::Bool=false)
	
	# unpacking fields
	(; design, components, onset, noisetype) = simulation 

	# create epoch data / erps
	erps = simulate(deepcopy(rng), components,simulation)

	if !return_epoched
		# we only need to simulate onsets & pull everything together, if we 
		# want a continuous EEG 	
		
		onsets = generate(deepcopy(rng),onset,simulation)
		
		# XXX todo: Separate Subjects in Time by adding offset to onsets!!

		# combine erps with onsets
		max_length = Int(ceil(maximum(onsets))) .+ maxlength(components)
		
		n_subj = length(size(design))==1 ? 1 : size(design)[2]
		n_trial = size(design)[1]
		eeg = zeros(max_length,n_subj)
		# not all designs have multiple subjects
		for s in 1:n_subj
			for i in 1:n_trial
				one_onset = onsets[CartesianIndex(i, s)]
				eeg[one_onset:one_onset+maxlength(components)-1,s] .+= @view erps[:, (s-1)*n_trial+i]
			end
		end	
	else
		eeg = erps
		onsets = [] # this is still a bit ugly
	end
	add_noise!(rng,noisetype,eeg)

	return convert(eeg,onsets,design;reshape=!return_epoched)

end


"""
Simulates erp data given the specified parameters 
"""
function simulate(rng, components::Vector{<:AbstractComponent},simulation::Simulation)

	epoch_data = zeros(maxlength(components), length(simulation.design))

	# Simulate each component
	for c in components
		# add them up


		epoch_data[1:length(c),:] += simulate(rng,c,simulation)
	end
	return epoch_data
end


function add_noise!(rng,noisetype::AbstractNoise,eeg)

	# generate noise
	noise = gen_noise(deepcopy(rng), noisetype, length(eeg))
	
	noise = reshape(noise, size(eeg))
	
	# add noise to data
	eeg .+= noisetype.noiselevel .* noise
	
end


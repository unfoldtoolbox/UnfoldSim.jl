
# helper to move input ::Component to ::Vector{Component}
Simulation(design::AbstractDesign,component::AbstractComponent,onset::AbstractOnset,noisetype::AbstractNoise) = Simulation(design,[component],onset,noisetype)

# by default no noise
Simulation(design::AbstractDesign,component,onset::AbstractOnset) = Simulation(design,component,onset,NoNoise())

"""
Simulate eeg data given a simulation design, effect sizes and variances

make use of `return_epoched=true` to skip the Onset-calculation + conversion to continuous data and get the epoched data directly
"""

simulate(rng,design::AbstractDesign, signal,  onset::AbstractOnset, noise::AbstractNoise;kwargs...) = simulate(rng,Simulation(design, signal,  onset, noise);kwargs...)
function simulate(rng, simulation::Simulation;return_epoched::Bool=false)
	
	# unpacking fields
	(; design, components, onset, noisetype) = simulation 

	# create epoch data / erps
	erps = simulate(deepcopy(rng), components,simulation)
#	@debug size(erps)

	n_subj = length(size(design))==1 ? 1 : size(design)[2]
	n_trial = size(design)[1]
	n_ch = n_channels(components)


	if !return_epoched
		# we only need to simulate onsets & pull everything together, if we 
		# want a continuous EEG 	
		
		onsets = generate(deepcopy(rng),onset,simulation)
		
		# XXX todo: Separate Subjects in Time by adding offset to onsets!!

		# combine erps with onsets
		maxlen = maxlength(components)
		max_length = Int(ceil(maximum(onsets))) .+ maxlen
		
	
		eeg = zeros(n_ch,max_length,n_subj)

		# not all designs have multiple subjects
		if n_subj == 1; eeg = dropdims(eeg,dims=3); end
		
		# not all designs have multiple channels
		if n_ch == 1; eeg = dropdims(eeg,dims=1); end

		
		for e in 1:n_ch
			for s in 1:n_subj
				for i in 1:n_trial
					one_onset = onsets[CartesianIndex(i, s)]
					adderp!(eeg,erps,e,s,one_onset:one_onset+maxlen-1,(s-1)*n_trial+i)
				end
			end	
		end
	else
		eeg = erps
		onsets = [] # this is still a bit ugly

	end


	add_noise!(deepcopy(rng),noisetype,eeg)

	# create events data frame
    events = UnfoldSim.generate(design)

	# save the onsets in the events df
	events.latency = onsets[:,]

	return eeg, events

end


"""
Helper function to add inplace the erps to the EEG, but for both 2D (1 channel) and 3D (X channel case)
"""
function adderp!(eeg,erps::Vector,e,s,tvec,erpvec)
	@views eeg[tvec] .+=  erps[:, erpvec]
end
function adderp!(eeg,erps::Matrix,e,s,tvec,erpvec)
	@views eeg[tvec,s] .+=  erps[:, erpvec]
end
function adderp!(eeg,erps::AbstractArray,e,s,tvec,erpvec)
	@views eeg[e,tvec,s] .+=  erps[e,:, erpvec]
end



"""
Simulates erp data given the specified parameters 
"""
function simulate(rng, components::Vector{<:AbstractComponent},simulation::Simulation)
	if n_channels(components) > 1
		epoch_data = zeros(n_channels(components),maxlength(components), length(simulation.design))
	else
		epoch_data = zeros(maxlength(components), length(simulation.design))
	end

	for c in components
		simulateandadd!(epoch_data,c,simulation,rng)
	end
	return epoch_data
end
function simulateandadd!(epoch_data::AbstractMatrix,c,simulation,rng)
	@debug "matrix"
	@views epoch_data[1:length(c),:] .+= simulate(rng,c,simulation)
end
function simulateandadd!(epoch_data::AbstractArray,c,simulation,rng)
	@debug "3D Array"
		@views epoch_data[:,1:length(c),:] .+= simulate(rng,c,simulation)
end




function add_noise!(rng,noisetype::AbstractNoise,eeg)

	# generate noise
	noise = gen_noise(deepcopy(rng), noisetype, length(eeg))
	
	noise = reshape(noise, size(eeg))
	
	# add noise to data
	eeg .+= noise
	
end


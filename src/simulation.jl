
# helper to move input ::Component to ::Vector{Component}
Simulation(design::AbstractDesign,component::AbstractComponent,onset::AbstractOnset,noisetype::AbstractNoise) = Simulation(design,[component],onset,noisetype)

"""
Simulate eeg data given a simulation design, effect sizes and variances
"""
simulate(rng,design, signal,  onset, noise) = simulate(rng,Simulation(design, signal,  onset, noise))
function simulate(rng, simulation::Simulation)
	
	# unpacking fields
	(; design, components, onset, noisetype) = simulation 

	# create epoch data / erps
	erps = simulate(deepcopy(rng), components,simulation)

	onsets = generate(deepcopy(rng),onset,simulation)
	
	# XXX todo: Separate Subjects in Time by adding offset to onsets!!

	# combine erps with onsets
	max_length = Int(ceil(maximum(onsets))) .+ maxlength(components)
	#eeg_continuous = Array{Float64,2}(0,max_length,n_subj)
	
	n_subj = length(size(design))==1 ? 1 : size(design)[2]
	n_trial = size(design)[1]
	eeg_continuous = zeros(max_length,n_subj)
	# not all designs have multiple subjects
	for s in 1:n_subj
		for i in 1:n_trial
			one_onset = onsets[CartesianIndex(i, s)]
			eeg_continuous[one_onset:one_onset+maxlength(components)-1,s] .+= @view erps[:, (s-1)*n_trial+i]
		end
	end	

	add_noise!(rng,noisetype,eeg_continuous)
	
	return convert(eeg_continuous,onsets,design)

end


"""
Simulates erp data given the specified parameters 
"""
function simulate(rng, components::Vector{<:AbstractComponent},simulation::Simulation)

	epoch_data = zeros(maxlength(components), length(simulation.design))

	# Simulate each component
	for c in components
		# add them up


		epoch_data += simulate(rng,c,simulation)
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


"""
Function to convert output similar to unfold (data, evts)
"""
function convert(eeg, onsets, design)
	data = eeg[:,]
	evt = UnfoldSim.generate(design)
	
	evt.latency = (onsets' .+ range(0,size(eeg,2)-1).*size(eeg,1) )'[:,]

	if :d ∈	names(evt)
    	select!(evt, Not([:dv]))
	end

	return data,evt
	
end


"""
Pads array with specified value, length
"""
function padarray(arr, len, val)
	for l in len
		pad = fill(val, abs(l))
		arr = l > 0 ? vcat(arr, pad) : vcat(pad, arr)
	end
	return arr
end

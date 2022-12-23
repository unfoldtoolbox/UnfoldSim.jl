
# helper to move input ::Component to ::Vector{Component}
Simulation(design::AbstractDesign,component::AbstractComponent,onset::AbstractOnset,noisetype::AbstractNoise) = Simulation(design,[component],onset,noisetype)

"""
Simulate eeg data given a simulation design, effect sizes and variances
"""
function simulate(rng, simulation::Simulation)
	
	# unpacking fields
	(; design, components, onset, noisetype) = simulation 
	(;n_item, n_subj) = design

	# create epoch data / erps
	erps = simulate(deepcopy(rng), design, components)

	onsets = gen_onsets(deepcopy(rng),simulation)
	
	# combine erps with onsets
	max_length = maximum(onsets) .+ maxlength(components)
	#eeg_continuous = Array{Float64,2}(0,max_length,n_subj)
	eeg_continuous = zeros(max_length,n_subj)
	for s in 1:n_subj
		for i in 1:n_item
			one_onset = onsets[CartesianIndex(i, s)]
			eeg_continuous[one_onset:one_onset+maxlength(components)-1,s] .+= @view erps[:, (s-1)*n_item+i]
		end
	end	

	add_noise!(rng,eeg_continuous,noisetype)
	
	return eeg_continuous, onsets
end


function add_noise!(rng,eeg,noisetype)

	# generate noise
	noise = gen_noise(deepcopy(rng), noisetype, length(eeg))
	
	noise = reshape(noise, size(eeg))
	
	# add noise to data
	eeg .+= noisetype.noiselevel .* noise
	
end

"""
Simulates erp data given the specified parameters 
"""
function simulate(rng, design::AbstractDesign, components::Vector{<:AbstractComponent})

	epoch_data = Array{Float64}(Int(length(components)), dims(design))

	# Simulate each component
	for c in components
		# add them up
		epoch_data += simulate(rng,c,design)
	end
	return epoch_data
end


"""
Function to convert output similar to unfold (data, evts)
"""
function convert(eeg, onsets, design)
	data = eeg[:,]
	evt = UnfoldSim.generate(design)
	
	evt.latency = (onsets' .+ range(0,size(eeg,2)-1).*size(eeg,1) )'[:,]
	
	rename!(evt,:subj => :subject)

    select!(evt, Not([:dv]))

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

# helper to move input ::Component to ::Vector{Component}
Simulation(design::ExperimentDesign,component::Component,onset::AbstractOnset,noisetype::AbstractNoise) = Simulation(design,[component],onset,noisetype)

"""
Simulate eeg data given a simulation design, effect sizes and variances
"""
function simulate(rng, simulation)
	
	# unpacking fields
	(; design, components, onset, noisetype) = simulation 
	(;n_item, n_subj) = design

	# create epoch data / erps
	erps = simulate_erps(deepcopy(rng), design, components)

	onsets = gen_onsets(rng,simulation)
	
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
	eeg += noisetype.noiselevel .* noise
end

"""
Simulates erp data given the specified parameters 
"""
function simulate_erps(rng, design, components)

	# unpacking fields
	(; n_subj, n_item) = design

	epoch_data = []

	# for each components
	for (; basis, formula, contrasts, β, σ_ranef, σ_res) in components

		# create model
		m = MixedModels.MixedModel(formula, generate(design), contrasts=contrasts)

		# limit runtime (in seconds)
		m.optsum.maxtime = 1
	
		# fit mixed model to experiment design and dummy data
		refit!(m, progress=false)

		# empty epoch data
		epoch_data_component = zeros(Int(length(basis)), n_subj*n_item)

		# residual variance for lmm
		σ_lmm = σ_res # 0.0001
			
		# iterate over each timepoint
		for t in eachindex(basis)

			# select weight from basis
			b = basis[t]
			
			# update random effects parametes of model
			if σ_ranef !== nothing
				k = (collect(keys(σ_ranef))...,)
				v = b .* (collect(values(σ_ranef))...,) ./ σ_lmm

				namedre = NamedTuple{k}(v)
				
				MixedModelsSim.update!(m; namedre...)
			end

			# simulate with new parameters
			simulate!(deepcopy(rng), m, β = b .* [β...], σ = σ_lmm)

			# save data to array
			epoch_data_component[t, :] = m.y
		end

		push!(epoch_data, epoch_data_component)
	end

	epoch_data = +(epoch_data...)

	return epoch_data
end


"""
Function to convert output similar to unfold (data, evts)
"""
function convert(eeg, onsets, design)
	# data 
	data = eeg[:,]

	# unpack
	(;n_subj) = design

	# generate design data frame
	ed = generate(design)

	# create & fill event data frame
	evts = DataFrame()
	a = collect(0:n_subj-1) * size(eeg, 1)
	evts.latency = (onsets' .+ a)'[:,]
	#evts.type .= "sim"
	insertcols!(evts, :type=>"sim")
	evts.trialnum = 1:size(evts, 1)
	for i in Set(ed.stimType)
		evts[!, Symbol("cond"*i)] = [(i == d ? 1 : 0) for d in ed.stimType]
	end
	evts.cond =ed.stimType
	evts.stimulus = ed.item
	evts.subject = ed.subj
	evts.urlatency = onsets[:,]

	return data, evts
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
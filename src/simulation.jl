"""
Simulation
"""
struct Simulation
	design
	components
	epochlen
	fixationlen
	noisetype
	noiselevel
end


"""
Simulate eeg data given a simulation design, effect sizes and variances
"""
function simulate(rng, simulation)
	
	# unpacking fields
	(; design, components, epochlen, fixationlen, noisetype, noiselevel) = simulation 
	(;n_item, n_subj) = design

	# create epoch data / erps
	erps = simulate_erps(deepcopy(rng), design, components)

	# deviation from average fixation samples
	std_fixationlen = 25 
	
	# min + max fixationlen
	min_fixationlen = fixationlen - std_fixationlen
	max_fixationlen = fixationlen + std_fixationlen

	# max length eeg 
	upperbound = (epochlen + max_fixationlen) * n_item + (epochlen + max_fixationlen)

	# sample different fixationlens & combine them with epoch lens
	fixationlens = rand(deepcopy(rng), min_fixationlen:max_fixationlen, (n_item, n_subj))
	epochlens = [0, repeat([epochlen], n_item-1)...]
	onsets = accumulate(+, (fixationlens .+ epochlens), dims=1)

	# one hot encoding of onsets
	onsets_1hot = hcat([[Int.(i in col) for i in 1:upperbound] for col in eachcol(onsets)]...)

	# combine erps with onsets
	e = []
	for s in 1:n_subj
		e_s = zeros(upperbound)
		for i in 1:n_item
			onset = onsets[CartesianIndex(i, s)]
			e_s[onset:onset+epochlen-1] = erps[:, (s-1)*n_item+i]
		end
		push!(e, e_s)
	end
	eeg = hcat(e...)
	
	# generate noise
	noise = gen_noise(deepcopy(rng), noisetype, length(eeg))
	noise = imfilter(noise, Kernel.gaussian((5,)))
	noise = reshape(noise, size(eeg))
	
	# add noise to data
	eeg = eeg + noiselevel .* noise

	return eeg, onsets
end


"""
Simulates erp data given the specified parameters 
"""
function simulate_erps(rng, design, components)

	# unpacking fields
	(; n_subj, n_item) = design

	epoch_data = []

	# for each components
	for (;basis, formula, contrasts, β, σ_ranef, σ_res) in components
	
		# fit mixed model to experiment design and dummy data
		m = MixedModels.fit(
			MixedModels.MixedModel, 
			formula, 
			generate(design), 
			contrasts=contrasts
		)

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
	evts.type .= "sim"
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
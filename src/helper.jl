"""
Pads array with specified value, length
padarray(arr, len, val)
"""
padarray(arr::Vector,len::Tuple,val) = padarray(padarray(arr,len[1],val),len[2],val)
function padarray(arr::Vector, len::Int, val)
		pad = fill(val, abs(len))
		arr = len > 0 ? vcat(arr, pad) : vcat(pad, arr)
	return arr
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



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
function convert(eeg, onsets, design;reshape=true)
	evt = UnfoldSim.generate(design)
	if reshape
		data = eeg[:,]
		evt.latency = (onsets' .+ range(0,size(eeg,2)-1).*size(eeg,1) )'[:,]
	else
		data = eeg
	end

	if :d ∈	names(evt)
    	select!(evt, Not([:dv]))
	end

	return data,evt
	
end

"""
Takes an array of 'm' target coordinate arrays (1-by-3) and a matrix (n-by-3) of all available positions, and returns an array of size 'm' containing the indices of the respective items in 'pos' that are nearest to each of the target coordinates.
"""
function closest_srcs(coords_list, pos) 
	out = [];
	s = size(pos);
	dist = zeros(s[1]);
	diff = zeros(s[2]);
	for coords in coords_list
		for i=1:s[1] 
			for j=1:s[2]
				diff[j] = pos[i,j] - coords[j];
			end
			dist[i] = norm(diff);
		end
		push!(out,findmin(dist)[2]) #retain only the index of the minimum difference.
	end
	return out
end

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


# TODO: Transfer function to Unfold.jl
"""
Function to convert output similar to unfold (data, evts)
"""
function convert(eeg, onsets, design,n_ch,;reshape=true)
	evt = UnfoldSim.generate(design)
	@debug size(eeg)
	if reshape
		n_subj = length(size(design))==1 ? 1 : size(design)[2]

		if n_ch == 1
		data = eeg[:,]
		
		evt.latency = (onsets' .+ range(0,size(eeg,2)-1).*size(eeg,1) )'[:,]
		elseif n_subj == 1
			data = eeg
			@debug size(onsets)
			evt.latency = onsets
		else # multi subject + multi channel
			data = eeg[:,:,]
			evt.latency = (onsets' .+ range(0,size(eeg,3)-1).*size(eeg,2) )'[:,]
		end
	else
		data = eeg
	end

	return data,evt
	
end


"""
	closest_src(coords_list::AbstractVector{<:AbstractVector}, pos)
	closest_src(coords::Vector{<:Real}, pos) 
	
Takes an array of 'm' target coordinate vector (size 3) (or vector of vectors) and a matrix (n-by-3) of all available positions, and returns an array of size 'm' containing the indices of the respective items in 'pos' that are nearest to each of the target coordinates.
"""
closest_src(coords_list::AbstractVector{<:AbstractVector}, pos)  = closest_src.(coords_list, Ref(pos))

function closest_src(coords::Vector{<:Real}, pos) 
	s = size(pos);
	dist = zeros(s[1]);
	diff = zeros(s[2]);
		for i=1:s[1] 
			for j=1:s[2]
				diff[j] = pos[i,j] - coords[j];
			end
			dist[i] = norm(diff);
		end
	return findmin(dist)[2]
	
	
end

"""
	closest_src(head::Hartmut,label::String)
Returns src-`ix` of the Headmodel `Hartmut` which is closest to the average of the `label`.

!!! important
	We use the average in eucledean space, but the cortex is a curved surface. In most cases they will not overlap. Ideally we would calculate the average on the surface, but this is a bit more complex to do (you'd need to calculate the vertices etc.)

```julia
hartmut = headmodel()
pos = closest_src(hartmut=>"Left Middle Temporal Gyrus, posterior division")
```
"""
function closest_src(head::Hartmut,label::String)

	pos = head.cortical["pos"]
	ix = findall(head.cortical["label"] .== label)
	@assert sum(ix)>0 """could not find label $label in hartmut.cortical["label"] - try unique(hartmut.cortical["label"]) for a list"""

	ix = UnfoldSim.closest_src(mean(pos[ix,:],dims=1)[1,:],pos)
	return ix
end
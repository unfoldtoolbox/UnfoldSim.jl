"""
<<<<<<< HEAD
Pads array with specified value, length
pad_array(arr, len, val)
"""
pad_array(arr::Vector, len::Tuple, val) =
    pad_array(pad_array(arr, len[1], val), len[2], val)
=======
    pad_array(arr::Vector, len::Int, val)    
    pad_array(arr::Vector, len::Tuple, val)

Pads the input array `arr` with a specified value `val` either before or after the existing elements, based on the sign of `len`.

# Arguments
- `arr::Vector`: The input array to be padded.
- `len::Union{Int, Tuple}`: The number of times that `val` should be added.
    If `len` is negative, the values are added before the existing array. Otherwise, they are added after the existing array.
    If `len` is a tuple, `pad_array` is called twice which enables padding the array before and after in one function call.
- `val`: The value to be used for padding.

# Returns
- `Vector`: Padded vector with the new length `length(arr) + sum(abs.(len))`.

# Examples
```julia-repl
# Create an array that will be padded
julia> my_array = rand(5)
5-element Vector{Float64}:
 0.0420017254437951
 0.19179144973603235
 0.5388760239550549
 0.6973699906283798
 0.9966598131018376

 # Pad the array with zeros before the original array
julia> pad_array(my_array, -2, 0)
7-element Vector{Float64}:
 0.0
 0.0
 0.0420017254437951
 0.19179144973603235
 0.5388760239550549
 0.6973699906283798
 0.9966598131018376

# Pad the array with the value 5 before and after the original array
julia> pad_array(my_array, (-2, 1), 5)
8-element Vector{Float64}:
 5.0
 5.0
 0.0420017254437951
 0.19179144973603235
 0.5388760239550549
 0.6973699906283798
 0.9966598131018376
 5.0
```
"""
pad_array(arr::Vector, len::Tuple, val) =
    pad_array(pad_array(arr, len[1], val), len[2], val)

>>>>>>> main
function pad_array(arr::Vector, len::Int, val)
    pad = fill(val, abs(len))
    arr = len > 0 ? vcat(arr, pad) : vcat(pad, arr)
    return arr
end


<<<<<<< HEAD

"""
Obsolete - # TODO: Transfer function to Unfold.jl

Function to convert output similar to unfold (data, events)
"""
function convert(eeg, onsets, design, n_chan; reshape = true)
    events = UnfoldSim.generate_events(design)
    @debug size(eeg)
    if reshape
        n_subjects = length(size(design)) == 1 ? 1 : size(design)[2]

        if n_chan == 1
            data = eeg[:,]

            events.latency = (onsets' .+ range(0, size(eeg, 2) - 1) .* size(eeg, 1))'[:,]
        elseif n_subjects == 1
            data = eeg
            @debug size(onsets)
            events.latency = onsets
        else # multi subject + multi channel
            data = eeg[:, :]
            events.latency = (onsets' .+ range(0, size(eeg, 3) - 1) .* size(eeg, 2))'[:,]
        end
    else
        data = eeg
    end

    return data, events

end


"""
	closest_src(coords_list::AbstractVector{<:AbstractVector}, pos)
	closest_src(coords::Vector{<:Real}, pos) 
	
Takes an array of 'm' target coordinate vector (size 3) (or vector of vectors) and a matrix (n-by-3) of all available positions, and returns an array of size 'm' containing the indices of the respective items in 'pos' that are nearest to each of the target coordinates.
"""
closest_src(coords_list::AbstractVector{<:AbstractVector}, pos) =
    closest_src.(coords_list, Ref(pos))

=======
"""
	closest_src(coords::Vector{<:Real}, pos)

Return the index of the position that is closest to the target coordinates (using the Euclidean distance).

# Arguments
- `coords::Vector{<:Real}`: Target coordinate vector (typically with length 3).
- `pos`: Matrix that contains all possible positions.
    Its dimensions are the number of positions `n` times number of entries in the position coordinate vectors (usually 3) i.e. `n x 3`.

# Returns
- `Int`: Row index of the position in `pos` that is closest to `coords`.

# Examples
```julia-repl
julia> positions = [0 0 0; 2 2 2; 3 3 3; 4 4 4]
4×3 Matrix{Int64}:
 0  0  0
 2  2  2
 3  3  3
 4  4  4

julia> target_coordinates = [0.5, 0, 0.5]
3-element Vector{Float64}:
 0.5
 0.0
 0.5

 julia> UnfoldSim.closest_src(target_coordinates, positions)
1
```
"""
>>>>>>> main
function closest_src(coords::Vector{<:Real}, pos)
    s = size(pos)
    dist = zeros(s[1])
    diff = zeros(s[2])
    for i ∈ 1:s[1]
        for j ∈ 1:s[2]
            diff[j] = pos[i, j] - coords[j]
        end
        dist[i] = norm(diff)
    end
    return findmin(dist)[2]


end

"""
	closest_src(coords_list::AbstractVector{<:AbstractVector}, pos)

When applied to a vector of target coordinate vectors, return the index of the closest position in `pos` for each of the targets. 

# Returns
- `Vector`: Index of the closest position in `pos` for each vector in `coords_list`.

# Examples
```julia-repl
julia> positions = [0 0 0; 2 2 2; 3 3 3; 4 4 4];
julia> target_coordinate_vectors = [[0.5, 0, 0.5], [10, 1, 8]];

julia> UnfoldSim.closest_src(target_coordinate_vectors, positions)
2-element Vector{Int64}:
 1
 4
```
"""
closest_src(coords_list::AbstractVector{<:AbstractVector}, pos) =
    closest_src.(coords_list, Ref(pos))

"""
    closest_src(head::Hartmut, label::String)

Return the index of the Harmut model source point that is closest (using Euclidean distance) to the mean position of the source points matching the given `label`.

!!! important
	We use the average in Euclidean space, but the cortex is a curved surface. In most cases they will not overlap. Ideally we would calculate the average on the surface, but this is a bit more complex to do (you'd need to calculate the vertices etc.).

# Arguments
- `head::Hartmut`: Headmodel of type `Hartmut`.
- `label::String`: Label of the source point(s) of interest.

# Returns
- `Int`: Index of the closest Hartmut source point.

# Examples
```julia-repl
julia> hartmut = Hartmut();
Please cite: HArtMuT: Harmening Nils, Klug Marius, Gramann Klaus and Miklody Daniel - 10.1088/1741-2552/aca8ce

julia> label = "Right Cingulate Gyrus, posterior division";

julia> UnfoldSim.closest_src(hartmut, label)
1875
```
"""
function closest_src(head::Hartmut, label::String)

    pos = head.cortical["pos"]
    ix = findall(head.cortical["label"] .== label)
    @assert sum(ix) > 0 """Could not find label $label in hartmut.cortical["label"] - try unique(hartmut.cortical["label"]) for a list."""

    ix = UnfoldSim.closest_src(mean(pos[ix, :], dims = 1)[1, :], pos)
    return ix
end


<<<<<<< HEAD


# One channel case

"""
    epoch(data::AbstractVector, args...; kwargs...)
    epoch(
        data::AbstractArray{T,2},
        events,
        τ::Tuple{Number,Number},
        sfreq;
        eventtime::Symbol = :latency,
    ) where {T<:Union{Missing,Number}}
Helper function to epoch data.

Adapted from Unfold.jl: https://github.com/unfoldtoolbox/Unfold.jl/blob/b3a21c2bb7e93d2f45ec64b0197f4663a6d7939a/src/utilities.jl#L40

"""
function epoch(data::AbstractVector, args...; kwargs...)
    data_r = reshape(data, (1, :))
    ep = epoch(data_r, args...; kwargs...)
    return dropdims(ep; dims = 1)
end
=======
"""
    epoch(
        data::AbstractArray{T,2},
        events,
        τ::Tuple{Number,Number},
        sfreq;
        eventtime::Symbol = :latency,
    ) where {T<:Union{Missing,Number}}

Helper function to segment continuous `data` into epochs based on the given `events` and the time window `τ`.
>>>>>>> main

Adapted from Unfold.jl: https://github.com/unfoldtoolbox/Unfold.jl/blob/b3  a21c2bb7e93d2f45ec64b0197f4663a6d7939a/src/utilities.jl#L40

# Arguments
- `data::AbstractArray{T,2}`: Continuous data with the dimensions `channels x continuous_time`.
- `events`: Events data frame (based on the design) with latencies.
- `τ::Tuple{Number,Number}`: Time window for epoching in s.
- `sfreq`: Sampling frequency in Hz.

# Keyword arguments
- `eventtime::Symbol = :latency`: The name of the column in `events` that contains the latencies.

# Returns
- `Array`: Epoched data with the dimensions `channels x times x event` (or `times x event` for the single channel case).

# Examples
```julia-repl
# One channel example
julia> data, events = UnfoldSim.predef_eeg();

julia> size(data), size(events)
((120199,), (2000, 3))

julia> UnfoldSim.epoch(data, events, (-0.2, 1), 100)
121×2000 Matrix{Float64}:
  0.114127   -0.105347     …  -0.125485   1.6383
  0.128198    0.0474874        0.0112935  1.28122
  0.0547917  -0.0832977       -0.126181   0.850062
  0.0992842  -0.230224        -0.0449072  0.496583
  0.024461   -0.175023        -0.0223837  0.170389
  0.165133   -0.000793527  …  -0.0278197  0.104454
  ⋮                        ⋱              
 -0.362249    2.91297          0.546699   0.0
 -0.199265    2.58394          0.171159   0.0
 -0.184075    2.34611         -0.0269841  0.0
 -0.13901     2.11971         -0.0552873  0.0
 -0.0674085   1.74561      …  -0.187959   0.0
```
"""
function epoch(
    data::AbstractArray{T,2},
    events,
    τ::Tuple{Number,Number},
    sfreq;
    eventtime::Symbol = :latency,
) where {T<:Union{Missing,Number}}
    # data: channels x times

    # partially taken from EEG.jl

    n_epochs = size(events, 1)

    times = range(τ[1], stop = τ[2], step = 1 ./ sfreq)
    length_epochs = length(times)
    n_chan = size(data, 1)
    epochs = Array{T}(undef, Int(n_chan), Int(length_epochs), Int(n_epochs))


    # User feedback
    @debug "Creating epochs: $n_chan x $length_epochs x $n_epochs"

    for si ∈ 1:size(events, 1)
        # d_start and d_end are the start and end of the epoch (in samples) in the data
        d_start = Int(round(events[si, eventtime]) + times[1] .* sfreq)
        d_end = Int(round(events[si, eventtime]) + times[end] .* sfreq)

        # e_start and e_end are the start and end within the epoch (in samples)
        e_start = 1
        e_end = length_epochs
        #println("d: $(size(data)),e: $(size(epochs)) | $d_start,$d_end,$e_start,$e_end | $(events[si,eventtime])")

        # Case that the start of the epoch is before the start of the data/recording (e.g. if the start is before i.e. negative relative to the event)
        if d_start < 1
            #@warn "d_start $d_start"
            e_start = e_start + (-d_start + 1)
            d_start = 1
        end
        # Case that the end of the epoch is after the end of the data/recording
        if d_end > size(data, 2)
            e_end = e_end - (d_end - size(data, 2))
            d_end = size(data, 2)
        end
        #println("d: $(size(data)),e: $(size(epochs)) | $d_start,$d_end,$e_start,$e_end | $(events[si,eventtime])")
        epochs[:, e_start:e_end, si] = data[:, d_start:d_end]
    end
    return epochs
end

"""
   epoch(data::AbstractVector, args...; kwargs...)

   One channel case. The data is reshaped and then passed to the multi-channel epoch method.
"""
function epoch(data::AbstractVector, args...; kwargs...)
    data_r = reshape(data, (1, :))
    ep = epoch(data_r, args...; kwargs...)
    return dropdims(ep; dims = 1)
end

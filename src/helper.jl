"""
Pads array with specified value, length
pad_array(arr, len, val)
"""
pad_array(arr::Vector, len::Tuple, val) =
    pad_array(pad_array(arr, len[1], val), len[2], val)
function pad_array(arr::Vector, len::Int, val)
    pad = fill(val, abs(len))
    arr = len > 0 ? vcat(arr, pad) : vcat(pad, arr)
    return arr
end



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
	closest_src(head::Hartmut,label::String)
Returns src-`ix` of the Headmodel `Hartmut` which is closest to the average of the `label`.

!!! important
	We use the average in eucledean space, but the cortex is a curved surface. In most cases they will not overlap. Ideally we would calculate the average on the surface, but this is a bit more complex to do (you'd need to calculate the vertices etc.)

```julia
hartmut = headmodel()
pos = closest_src(hartmut=>"Left Middle Temporal Gyrus, posterior division")
```
"""
function closest_src(head::Hartmut, label::String)

    pos = head.cortical["pos"]
    ix = findall(head.cortical["label"] .== label)
    @assert sum(ix) > 0 """could not find label $label in hartmut.cortical["label"] - try unique(hartmut.cortical["label"]) for a list"""

    ix = UnfoldSim.closest_src(mean(pos[ix, :], dims = 1)[1, :], pos)
    return ix
end




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
Helper function to epoch data

Adapted from Unfold.jl: https://github.com/unfoldtoolbox/Unfold.jl/blob/b3a21c2bb7e93d2f45ec64b0197f4663a6d7939a/src/utilities.jl#L40

"""
function epoch(data::AbstractVector, args...; kwargs...)
    data_r = reshape(data, (1, :))
    ep = epoch(data_r, args...; kwargs...)
    return dropdims(ep; dims = 1)
end

function epoch(
    data::AbstractArray{T,2},
    events,
    τ::Tuple{Number,Number},
    sfreq;
    eventtime::Symbol = :latency,
) where {T<:Union{Missing,Number}}
    # data: channels x times

    # partial taken from EEG.jl

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

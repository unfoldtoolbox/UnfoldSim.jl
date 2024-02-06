# helper to move input ::Component to ::Vector{Component}
Simulation(
    design::AbstractDesign,
    component::AbstractComponent,
    onset::AbstractOnset,
    noisetype::AbstractNoise,
) = Simulation(design, [component], onset, noisetype)

"""
    simulate(
    rng,
    design::AbstractDesign,
    signal,
    onset::AbstractOnset,
    noise::AbstractNoise = NoNoise();
    return_epoched=false,
    )

Main simulation function, given `Design`, [Array of] `Component`, `Onset` and optional [`Noise`], returns continuous or epoched EEG data.

## optional
- `return_epoched` (Bool, default: `false`):  Skip the Onset-calculation and conversion to continuous data and return the epoched data directly (see also remarks below).

## Return 
Depending on the design, the components and on `return_epoched`, the output can be a 1-D, 2-D, 3-D or 4-D Array. For example, a 4-D Array would have the dimensions `channels x time x trials x subjects`

## Notes
Some remarks to how the noise is added:
  - If `return_epoched = true` and `onset =NoOnset()` the noise is added to the epoched data matrix
  - If `onset` is not `NoOnset`, a continuous eeg signal is created and the noise is added to this i.e. this means that the noise won't be the same as in the `onset = NoOnset()` case even if `return_epoched = true`.
  - The case `return_epoched = false` and `onset = NoOnset()` is not possible and therefore covered by an assert statement

"""

simulate(
    rng,
    design::AbstractDesign,
    signal,
    onset::AbstractOnset,
    noise::AbstractNoise = NoNoise();
    kwargs...,
) = simulate(rng, Simulation(design, signal, onset, noise); kwargs...)

function simulate(rng, simulation::Simulation; return_epoched::Bool = false)
    (; design, components, onset, noisetype) = simulation

    # equivalent to !(isa(onset,NoOnset) && return_epoched == false)
    @assert !isa(onset, NoOnset) || return_epoched == true "It is not possible to get continuous data without specifying a specific onset distribution. Please either specify an onset distribution (other than `NoOnset`) or set `return_epoched = true` to get epoched data without overlap."

    # create epoch data / erps
    erps = simulate_responses(deepcopy(rng), components, simulation)

    # create events data frame
    events = UnfoldSim.generate_events(design)

    if isa(onset, NoOnset)
        # reshape the erps such that the last dimension is split in two dimensions (trials per subject and subject)
        # such that the resulting dimensions are dimensions: channels x times x trials x subjects
        # TODO: This assumes a balanced design, but create_continuous_eeg also assumes this, so we should be fine ;)
        size_erps = size(erps)
        eeg = reshape(erps, size_erps[1:end-1]..., size(design)...)
    else # if there is an onset distribution given the next step is to create a continuous EEG
        eeg, latencies = create_continuous_eeg(rng, erps, simulation)
        events.latency = latencies
    end

    add_noise!(deepcopy(rng), noisetype, eeg)

    # In case the data should be epoched & onset distribution is given i.e. the signals might be overlapping
    if return_epoched && !isa(onset, NoOnset)

        # use epoch function to epoch the continuous (possibly overlapping) signal
        if length(size(design)) == 1 # if there is only one subject
            eeg = epoch(eeg, events, (0, maxlength(components) - 1), 1)
        else # multi-subject case
            evt_epoch = groupby(events, :subject) |> collect
            # Epoch data per subject
            # Note: Ref() is needed to prevent broadcasting of τ and sfreq (due to applying epoch elementwise)
            eeg =
                epoch.(
                    eachslice(eeg, dims = length(size(eeg))),
                    evt_epoch,
                    Ref((0, maxlength(components) - 1)),
                    Ref(1),
                )

            # TODO: This assumes a balanced design, but create_continuous_eeg also assumes this, so we should be fine ;)
            # Concatenate the epoched data of all subjects again.
            eeg = reshape(reduce(hcat, vec.(eeg)), size(eeg[1])..., length(eeg))
            #cat(eeg..., dims = length(size(erps)) + 1) #TODO: find a way to use always the right dims
        end


    end
    return eeg, events

end

function create_continuous_eeg(rng, erps, simulation)

    (; design, components, onset, noisetype) = simulation

    n_subj = length(size(design)) == 1 ? 1 : size(design)[2]
    n_trial = size(design)[1]
    n_ch = n_channels(components)

    # we only need to simulate onsets & pull everything together, if we 
    # want a continuous EEG 	
    onsets = simulate_onsets(deepcopy(rng), onset, simulation)

    # flatten onsets (since subjects are concatenated in the events df)
    latencies = onsets[:,]

    # combine erps with onsets
    max_length_component = maxlength(components)
    max_length_continuoustime = Int(ceil(maximum(onsets))) .+ max_length_component


    eeg = zeros(n_ch, max_length_continuoustime, n_subj)

    for e = 1:n_ch
        for s = 1:n_subj
            for i = 1:n_trial
                one_onset = onsets[CartesianIndex(i, s)]
                adderp!(
                    eeg,
                    erps,
                    e,
                    s,
                    one_onset:one_onset+max_length_component-1,
                    (s - 1) * n_trial + i,
                )
            end
        end
    end

    # not all designs have multiple subjects
    if n_subj == 1
        eeg = dropdims(eeg, dims = 3)
    end

    # not all designs have multiple channels
    if n_ch == 1
        eeg = dropdims(eeg, dims = 1)
    end

    return eeg, latencies
end


"""
Helper function to add inplace the erps to the EEG, but for both 2D (1 channel) and 3D (X channel case)
"""
function adderp!(eeg, erps::Vector, e, s, tvec, erpvec)
    @views eeg[e, tvec, s] .+= erps[:, erpvec]
end
function adderp!(eeg, erps::Matrix, e, s, tvec, erpvec)#
    @views eeg[e, tvec, s] .+= erps[:, erpvec]
end
function adderp!(eeg, erps::AbstractArray, e, s, tvec, erpvec)
    @views eeg[e, tvec, s] .+= erps[e, :, erpvec]
end

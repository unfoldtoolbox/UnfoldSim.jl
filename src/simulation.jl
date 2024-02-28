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

Main simulation function, given `Design`, [Array of] `Component`, `Onset` and optional [`Noise`], returns continuous or epoched signal.

## optional
- `return_epoched` (Bool, default: `false`):  Skip the Onset-calculation and conversion to continuous data and return the epoched data directly (see also remarks below).

## Return 
Depending on the design, the components and on `return_epoched`, the output can be a 1-D, 2-D, 3-D or 4-D Array. For example, a 4-D Array would have the dimensions `channels x time x trials x subjects`

## Notes
Some remarks to how the noise is added:
  - If `return_epoched = true` and `onset =NoOnset()` the noise is added to the epoched data matrix
  - If `onset` is not `NoOnset`, a continuous signal is created and the noise is added to this i.e. this means that the noise won't be the same as in the `onset = NoOnset()` case even if `return_epoched = true`.
  - The case `return_epoched = false` and `onset = NoOnset()` is not possible and therefore covered by an assert statement
  - `simulate(rng,design::SequenceDesign,...)`
    If no `design.rng` was defined for `SequenceDesign`, we replace it with the `simulation`-function call `rng` object

"""


function simulate(args...; kwargs...)
    @warn "No random generator defined, used the default (`Random.MersenneTwister(1)`) with a fixed seed. This will always return the same results and the user is strongly encouraged to provide their own random generator!"
    simulate(MersenneTwister(1), args...; kwargs...)
end

function simulate(
    rng::AbstractRNG,
    design::AbstractDesign,
    signal,
    onset::AbstractOnset,
    noise::AbstractNoise = NoNoise();
    kwargs...,
)

    if is_SequenceDesign(design)
        design = sequencedesign_add_rng(rng, design)
    end
    simulate(rng, Simulation(design, signal, onset, noise); kwargs...)
end


sequencedesign_add_rng(rng, design::AbstractDesign) = design
sequencedesign_add_rng(rng, design::RepeatDesign) =
    RepeatDesign(sequencedesign_add_rng(rng, design.design), design.repeat)
sequencedesign_add_rng(rng, design::SequenceDesign) =
    isnothing(design.rng) ?
    SequenceDesign(design.design, design.sequence, design.sequencelength, rng) : design


is_SequenceDesign(d::AbstractDesign) = false
is_SequenceDesign(d::RepeatDesign) = is_SequenceDesign(d.design)
is_SequenceDesign(d::SequenceDesign) = true


function simulate(rng::AbstractRNG, simulation::Simulation; return_epoched::Bool = false)
    (; design, components, onset, noisetype) = simulation

    # equivalent to !(isa(onset,NoOnset) && return_epoched == false)
    @assert !isa(onset, NoOnset) || return_epoched == true "It is not possible to get continuous data without specifying a specific onset distribution. Please either specify an onset distribution (other than `NoOnset`) or set `return_epoched = true` to get epoched data without overlap."

    # create epoch data / responses
    responses = simulate_responses(deepcopy(rng), components, simulation)

    # create events data frame
    events = UnfoldSim.generate_events(design)

    if isa(onset, NoOnset)
        # reshape the responses such that the last dimension is split in two dimensions (trials per subject and subject)
        # such that the resulting dimensions are dimensions: channels x times x trials x subjects
        # TODO: This assumes a balanced design, but create_continuous_signal also assumes this, so we should be fine ;)
        size_responses = size(responses)
        signal = reshape(responses, size_responses[1:end-1]..., size(design)...)
    else # if there is an onset distribution given the next step is to create a continuous signal
        signal, latencies = create_continuous_signal(rng, responses, simulation)
        events.latency = latencies
    end

    add_noise!(deepcopy(rng), noisetype, signal)

    # In case the data should be epoched & onset distribution is given i.e. the signals might be overlapping
    if return_epoched && !isa(onset, NoOnset)

        # use epoch function to epoch the continuous (possibly overlapping) signal
        if length(size(design)) == 1 # if there is only one subject
            signal = epoch(signal, events, (0, maxlength(components) - 1), 1)
        else # multi-subject case
            events_epoch = groupby(events, :subject) |> collect
            # Epoch data per subject
            # Note: Ref() is needed to prevent broadcasting of τ and sfreq (due to applying epoch elementwise)
            signal =
                epoch.(
                    eachslice(signal, dims = length(size(signal))),
                    events_epoch,
                    Ref((0, maxlength(components) - 1)),
                    Ref(1),
                )

            # TODO: This assumes a balanced design, but create_continuous_signal also assumes this, so we should be fine ;)
            # Concatenate the epoched data of all subjects again.
            signal = reshape(reduce(hcat, vec.(signal)), size(signal[1])..., length(signal))
            #cat(signal..., dims = length(size(responses)) + 1) #TODO: find a way to use always the right dims
        end


    end
    return signal, events

end


"""
    create_continuous_signal(rng, responses, simulation)
Based on the responses and simulation parameters, simulate onset latencies and add together a continuous signal.


"""
function create_continuous_signal(rng, responses, simulation)

    (; design, components, onset, noisetype) = simulation

    n_subjects = length(size(design)) == 1 ? 1 : size(design)[2]
    n_trials = size(design)[1]
    n_chan = n_channels(components)

    # we only need to simulate onsets & pull everything together, if we 
    # want a continuous signal 	
    onsets = simulate_onsets(deepcopy(rng), onset, simulation)

    # flatten onsets (since subjects are concatenated in the events df)
    latencies = onsets[:,]

    # combine responses with onsets
    max_length_component = maxlength(components)
    offset_range = maxoffset(simulation.components) - minoffset(simulation.components)
    max_length_continuoustime =
        Int(ceil(maximum(onsets))) .+ max_length_component .+ offset_range


    signal = zeros(n_chan, max_length_continuoustime, n_subjects)
    @debug size(signal), offset_range
    for e = 1:n_chan
        for s = 1:n_subjects
            for i = 1:n_trials
                one_onset = onsets[CartesianIndex(i, s)]
                add_responses!(
                    signal,
                    responses,
                    e,
                    s,
                    one_onset+minoffset(simulation.components):one_onset+max_length_component-1+maxoffset(
                        simulation.components,
                    ),
                    (s - 1) * n_trials + i,
                )
            end
        end
    end

    # not all designs have multiple subjects
    if n_subjects == 1
        signal = dropdims(signal, dims = 3)
    end

    # not all designs have multiple channels
    if n_chan == 1
        signal = dropdims(signal, dims = 1)
    end

    return signal, latencies
end


"""
    add_responses!(signal, responses::Vector, e, s, tvec, erpvec)
    add_responses!(signal, responses::Matrix, e, s, tvec, erpvec)
    add_responses!(signal, responses::AbstractArray, e, s, tvec, erpvec)
Helper function to add inplace the responses to the signal, but for both 2D (1 channel) and 3D (X channel case).
"""
function add_responses!(signal, responses::Vector, e, s, tvec, erpvec)
    @views signal[e, tvec, s] .+= responses[:, erpvec]
end
function add_responses!(signal, responses::Matrix, e, s, tvec, erpvec)#
    #    @debug size(signal), size(responses), e, s, size(tvec), size(erpvec)
    #@debug tvec, erpvec
    @views signal[e, tvec, s] .+= responses[:, erpvec]
end
function add_responses!(signal, responses::AbstractArray, e, s, tvec, erpvec)
    @views signal[e, tvec, s] .+= responses[e, :, erpvec]
end

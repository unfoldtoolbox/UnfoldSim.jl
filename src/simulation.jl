# helper to move input ::Component to ::Vector{Component}
Simulation(
    design::AbstractDesign,
    component::AbstractComponent,
    onset::AbstractOnset,
    noisetype::AbstractNoise,
) = Simulation(design, [component], onset, noisetype)


function simulate(design::AbstractDesign, signal, onset::AbstractOnset, args...; kwargs...)
    @warn "No random generator defined, used the default (`Random.MersenneTwister(1)`) with a fixed seed. This will always return the same results and the user is strongly encouraged to provide their own random generator!"
    simulate(MersenneTwister(1), design, signal, onset, args...; kwargs...)
end

"""
    simulate(
    rng::AbstractRNG,
    design::AbstractDesign,
    signal,
    onset::AbstractOnset,
    noise::AbstractNoise = NoNoise();
    return_epoched = false,
    )

	simulate(
    design::AbstractDesign,
    signal,
    onset::AbstractOnset,
    noise::AbstractNoise = NoNoise();
    return_epoched = false,
    )

Return continuous or epoched signal, given `Design`, [Array of] `Component`, `Onset` and 
optional [`Noise`] and [`RNG`]. Main simulation function.

# Arguments
- `design::AbstractDesign`: Desired experimental design.
- `signal`: `Component` for the desired signal.
- `onset::AbstractOnset`: Desired onset.
- `noise::AbstractNoise = NoNoise()` (optional): Desired noise.
- `rng::AbstractRNG` (optional): Random number generator, important to ensure reproducibility.

# Keyword arguments
- `return_epoched = false`: Skip the Onset-calculation and conversion to continuous data 
    and return the epoched data directly (see also Notes below).

# Returns
- `signal` : Generated signal. Depending on the design, on the components and on 
    `return_epoched`, the output can be a 1-D, 2-D, 3-D or 4-D Array. 
    For example, a 4-D Array would have the dimensions `channels x time x trials x subjects`.
- `events`: Generated events.

# Examples
Adapted from the quickstart tutorial in the UnfoldSim docs.
```julia-repl
julia> using UnfoldSim

julia> using Random # to get an RNG

julia> design =
    SingleSubjectDesign(; conditions = Dict(:cond_A => ["level_A", "level_B"])) |>
    x -> RepeatDesign(x, 10);

julia> signal = LinearModelComponent(; 
    basis = [0, 0, 0, 0.5, 1, 1, 0.5, 0, 0],
    formula = @formula(0 ~ 1 + cond_A),
    β = [1, 0.5],
);

julia> onset = UniformOnset(; width = 20, offset = 4);

julia> noise = PinkNoise(; noiselevel = 0.2);

julia> data, events = simulate(MersenneTwister(1), design, signal, onset, noise)
([-0.045646938524459196, 0.15784406738265955, 0.012640319497460443, 0.026669512219327673, 0.15329144053662508, 0.06412654786607011, -0.16766777448918685, -0.08012027590515228, 0.0020515088981202137, -0.24874482217391175  …  0.24621397283439814, 0.1710771262918883, -0.01527736524528042, 0.7639978745937471, 1.5600315092771557, 1.624219837479329, 1.2889713838347956, 0.26819223928179, 0.16535758767503866, 0.21291936972924855], 20×2 DataFrame
 Row │ cond_A   latency 
     │ String   Int64   
─────┼──────────────────
   1 │ level_A       18
   2 │ level_B       39
   3 │ level_A       45
   4 │ level_B       49
  ⋮  │    ⋮        ⋮
  18 │ level_B      258
  19 │ level_A      281
  20 │ level_B      303

julia> data1, events1 = simulate(design, signal, onset, noise)
┌ Warning: No random generator defined, used the default (`Random.MersenneTwister(1)`) with a fixed seed. This will always return the same results and the user is strongly encouraged to provide their own random generator!
```

# Notes
Some remarks to how the noise is added:
- If `return_epoched = true` and `onset = NoOnset()` the noise is added to the epoched data matrix.
- If `onset` is not `NoOnset`, a continuous signal is created and the noise is added to this 
    i.e. this means that the noise won't be the same as in the `onset = NoOnset()` case even if `return_epoched = true`.
- The case `return_epoched = false` and `onset = NoOnset()` is not possible and therefore 
    covered by an assert statement.
"""
simulate(
    rng::AbstractRNG,
    design::AbstractDesign,
    signal,
    onset::AbstractOnset,
    noise::AbstractNoise = NoNoise();
    kwargs...,
) = simulate(rng, Simulation(design, signal, onset, noise); kwargs...)


function simulate(rng::AbstractRNG, simulation::Simulation; return_epoched::Bool = false)
    (; design, components, onset, noisetype) = simulation

    # equivalent to !(isa(onset,NoOnset) && return_epoched == false)
    @assert !isa(onset, NoOnset) || return_epoched == true "It is not possible to get continuous data without specifying a specific onset distribution. Please either specify an onset distribution (other than `NoOnset`) or set `return_epoched = true` to get epoched data without overlap."

    # create epoch data / responses
    responses = simulate_responses(deepcopy(rng), components, simulation)

    # create events data frame
    events = UnfoldSim.generate_events(deepcopy(rng), design)

    if isa(onset, NoOnset)
        # reshape the responses such that the last dimension is split in two dimensions (trials per subject and subject)
        # such that the resulting dimensions are dimensions: channels x times x trials x subjects
        # TODO: This assumes a balanced design, but create_continuous_signal also assumes this, so we should be fine ;)
        size_responses = size(responses)
        signal = reshape(responses, size_responses[1:end-1]..., size(design)...)
    else # if there is an onset distribution given the next step is to create a continuous signal
        signal, latencies = create_continuous_signal(deepcopy(rng), responses, simulation)
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
    max_length_continuoustime = Int(ceil(maximum(onsets))) .+ max_length_component


    signal = zeros(n_chan, max_length_continuoustime, n_subjects)

    for e = 1:n_chan
        for s = 1:n_subjects
            for i = 1:n_trials
                one_onset = onsets[CartesianIndex(i, s)]
                add_responses!(
                    signal,
                    responses,
                    e,
                    s,
                    one_onset:one_onset+max_length_component-1,
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
    @views signal[e, tvec, s] .+= responses[:, erpvec]
end
function add_responses!(signal, responses::AbstractArray, e, s, tvec, erpvec)
    @views signal[e, tvec, s] .+= responses[e, :, erpvec]
end

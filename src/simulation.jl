# helper to move input ::Component to ::Vector{Component}
Simulation(
    design::AbstractDesign,
    component::AbstractComponent,
    onset::AbstractOnset,
    noisetype::AbstractNoise,
) = Simulation(design, [component], onset, noisetype)


function simulate(
    design::AbstractDesign,
    components,
    onset::AbstractOnset,
    args...;
    kwargs...,
)
    @warn "No random generator defined, used the default (`Random.MersenneTwister(1)`) with a fixed seed. This will always return the same results and the user is strongly encouraged to provide their own random generator!"
    simulate(MersenneTwister(1), design, components, onset, args...; kwargs...)
end

"""
    simulate(
    rng::AbstractRNG,
    design::AbstractDesign,
    components,
    onset::AbstractOnset,
    noise::AbstractNoise = NoNoise();
    return_epoched = false,
    )

	simulate(
    design::AbstractDesign,
    components,
    onset::AbstractOnset,
    noise::AbstractNoise = NoNoise();
    return_epoched = false,
    )

Simulate continuous or epoched signal, given `design`, [Array of] `component`, `onset` and 
optional `noise` and `rng`. Main simulation function.

# Arguments
- `rng::AbstractRNG` (optional): Random number generator, important to ensure reproducibility.
- `design::AbstractDesign`: Desired experimental design.
- `components`: `Component`(s) for the desired signal.
- `onset::AbstractOnset`: Desired inter-onset distance distribution.
- `noise::AbstractNoise = NoNoise()` (optional): Desired noise.

# Keyword arguments
- `return_epoched::Bool = false`: If set to `true` epoched data is returned, otherwise a continuous signal is returned (see also Notes below).

# Returns
- `(signal, events)::Tuple{Array, DataFrame}`:
    - `signal` : Generated signal. Depending on the design, on the components and on 
        `return_epoched`, the output can be a 1-D, 2-D, 3-D or 4-D Array. 
        For example, a 4-D Array would have the dimensions `channels x time x trials x subjects`.
    - `events`: Generated events data frame with latencies.

# Examples
Adapted from the [quickstart tutorial](https://unfoldtoolbox.github.io/UnfoldSim.jl/stable/generated/tutorials/quickstart/) in the UnfoldSim docs.
```julia-repl
julia> using Random # to get an RNG

julia> design =
    SingleSubjectDesign(; conditions = Dict(:cond_A => ["level_A", "level_B"])) |>
    x -> RepeatDesign(x, 10);

julia> component = LinearModelComponent(; 
    basis = [0, 0, 0, 0.5, 1, 1, 0.5, 0, 0],
    formula = @formula(0 ~ 1 + cond_A),
    β = [1, 0.5],
);

julia> onset = UniformOnset(; width = 20, offset = 4);

julia> noise = PinkNoise(; noiselevel = 0.2);

# Variant 1: Use a custom RNG.
julia> data, events = simulate(MersenneTwister(2), design, component, onset, noise);

julia> data
293-element Vector{Float64}:
 -0.013583193323430123
  0.09159433856866195
  ⋮
 -0.25190584567097907
 -0.20179992275876316

julia> events
20×2 DataFrame
 Row │ cond_A   latency 
     │ String   Int64   
─────┼──────────────────
   1 │ level_A        9
   2 │ level_B       20
   3 │ level_A       27
   4 │ level_B       37
  ⋮  │    ⋮        ⋮
  18 │ level_B      257
  19 │ level_A      271
  20 │ level_B      284
         13 rows omitted

# Variant 2: Without specifying an RNG, MersenneTwister(1) will be used for the simulation.
julia> data1, events1 = simulate(design, component, onset, noise);
┌ Warning: No random generator defined, used the default (`Random.MersenneTwister(1)`) with a fixed seed. This will always return the same results and the user is strongly encouraged to provide their own random generator!
```

# Notes
Some remarks on how the noise is added:
- If `return_epoched = true` and `onset = NoOnset()` the noise is added to the epoched data matrix.
- If `onset` is not `NoOnset`, a continuous signal is created and the noise is added to this 
    i.e. this means that the noise won't be the same as in the `onset = NoOnset()` case even if `return_epoched = true`.
- The case `return_epoched = false` and `onset = NoOnset()` is not possible and therefore 
    covered by an assert statement.

Additional remarks on the overlap of adjacent signals when `return_epoched = true`:
- If `onset = NoOnset()` there will not be any overlapping signals in the data because the onset calculation and conversion to a continuous signal is skipped.
- If an inter-onset distance distribution is given, a continuous signal(potentially with overlap) is constructed and partitioned into epochs afterwards.
"""
simulate(
    rng::AbstractRNG,
    design::AbstractDesign,
    components,
    onset::AbstractOnset,
    noise::AbstractNoise = NoNoise();
    kwargs...,
) = simulate(rng, Simulation(design, components, onset, noise); kwargs...)


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

"""TODO docstring
s: Vector of AbstractContinuousSignal and/or AbstractNoise 
NOTE:
1. It is assumed that `AbstractContinuousSignal.controlsignal` 
    is defined from the start of the simulation. If you require the artifact to start at
    a different time, please take care to pre-pad the controlsignal with zeros.
2. EEG simulation components must be multichannel, and it is assumed that
    the number and position of EEG channels are the same as in the artifact head model. 
"""
function simulate(rng::AbstractRNG,d::AbstractDesign,c::AbstractComponent,o::AbstractOnset,s::AbstractVector)
    @assert all(x -> x isa AbstractContinuousSignal || x isa AbstractNoise, s) "Artifact-related inputs should all be of type AbstractContinuousSignal or AbstractNoise"
    sim = Simulation(d, c, o, NoNoise()) # since generated controlsignal might depend on some aspect of the simulation

    # Current Assumptions:
    # events from artifacts - not considered/handled right now
    # EEG component is assumed to be multichannel and having the same number and positions of channels as the artifact.
    # noise is considered to be independent of EEG/artifacts.
    # controlsignals are generated only for AbstractContinuousSignal and not AbstractNoise.
    # currently any number of noise components and/or artifacts are allowed. e.g. we can simulate RedNoise as well as PinkNoise and multiple EyeMovements each starting from t=0.
    # eeg and artifact signals are simulated separately and then added together, padding to the larger of the two. Then noise is added to this signal.
    # multiple EyeMovement artifacts are currently allowed. (TODO: check if it makes sense to restrict this. This could be useful if different Eyemovements have different offsets?)

    println("Simulating EEG with no noise...")
    eeg_signal,evts = simulate(rng,d,c,o,NoNoise());
    @show(size(eeg_signal))

    println("Generating controlsignals...")
    controlsignal = generate_controlsignal.(deepcopy(rng),s,Ref(sim)) 
    # Vector (n_feat) of Matrix (__ x time) -> each kind of artifact could have a different shape of controlsignal, so keep them as elements of the vector rather than combining to a matrix and losing information of which row(s) corresp. to which artifact  
    @show size.(controlsignal) 
    
    println("Simulating continuous signals...")
    artifact_signal = simulate_continuoussignal.(deepcopy(rng),s,controlsignal,Ref(sim)); #TODO handle events: right now for simplicity assume no events are being returned
    
    println("Removing empty artifact signals...")
    filter!(signal -> size(signal) != (0,0), artifact_signal) # remove the empty arrays (e.g. from AbstractNoise, PowerLineNoise)

    combined_signals = [[eeg_signal] ; artifact_signal]

    row_counts = [size(mat, 1) for mat in combined_signals];
    if length(unique(row_counts)) > 1
        @warn "Simulated EEG and artifacts do not have the same number of channels: $(row_counts)"
    end

    max_cols = maximum([size(mat, 2) for mat in combined_signals])
    for i in 1:length(combined_signals)
        mat = combined_signals[i]
        if size(mat, 2) < max_cols
            combined_signals[i] = hcat(mat, zeros(size(mat, 1), max_cols - size(mat, 2)))
        end
    end

    sum_signals = reduce(+, combined_signals) # sum of eeg and artifacts

    println("Adding noise...")    
    add_noise!.(rng,[x for x in s if x isa AbstractNoise],Ref(sum_signals))

    return combined_signals, sum_signals, evts
end

"""
    create_continuous_signal(rng, responses, simulation)

Simulate onset latencies and add together a continuous signal, based on the given responses and simulation parameters. Helper function.

# Arguments
- `rng`: Random number generator, important to ensure reproducibility.
- `responses`: Responses to be combined with the given onsets.
- `simulation`: Simulation parameters, including design, components, onsets, and noisetype.

# Returns
- `(signal, latencies)::Tuple{Array, Array}`:
    - `signal` contains the generated signal. Has the dimensions `channels x continuous_time x subjects`.
    - `latencies` contains the onset latencies.

# Examples
```julia-repl
julia> using StableRNGs # to get an RNG

julia> design = SingleSubjectDesign(; conditions = Dict(:cond => ["natural", "artificial"]));

julia> c1 = LinearModelComponent(; basis = p100(), formula = @formula(0 ~ 1 + cond), β = [1, 0.5]);

julia> c2 = LinearModelComponent(; basis = p300(), formula = @formula(0 ~ 1), β = [2]);

julia> simulation = Simulation(design, [c1, c2], UniformOnset(; width = 0, offset = 30), PinkNoise());

julia> responses = simulate_responses(StableRNG(1), [c1, c2], simulation)
45×2 Matrix{Float64}:
 0.0        0.0
 0.0        0.0
 ⋮          
 0.0233794  0.0233794
 0.0        0.0

julia> signal, latencies = UnfoldSim.create_continuous_signal(StableRNG(1), responses, simulation);

julia> signal
106-element Vector{Float64}:
 0.0
 0.0
 0.0
 ⋮
 0.023379444289913343
 0.0
 0.0

 julia> latencies
2-element Vector{Int64}:
 31
 61
```
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
    add_responses!(signal, responses::Vector, e, s, tvec, trial_idx)
    add_responses!(signal, responses::Matrix, e, s, tvec, trial_idx)
    add_responses!(signal, responses::AbstractArray, e, s, tvec, trial_idx)

Add (in-place) the given `responses` to the `signal`, for both 2D (1 channel) and 3D (X channel case). Helper function.

# Arguments
- `signal`: Continuous EEG signal to be modified in place. Has the dimensions `channels x continuous_time x subjects`.
- `responses::Union{Vector, Matrix, AbstractArray}`: Responses to be added. In the multi-channel case, the dimensions are `channels x maxlength(components) x length(simulation.design)`, else `maxlength(components) x length(simulation.design)`.  The data for all the subjects and their respective trials is concatenated.
- `e`: Index of the channel (in `signal`) for which to add the response.
- `s`: Index of the subject (in `signal`) for which to add the response.
- `tvec`: Time points (indices in `signal`) at which to add the response.
- `trial_idx`: Index of the particular trial (in `responses`) from where the response is to be added.

# Returns
- `Nothing`: `signal` is modified in-place.

# Examples
```julia-repl
julia> signal, responses, tvec = zeros(5,15,2), ones(5,6), 1:5;

julia> UnfoldSim.add_responses!(signal, responses, 1, 2, tvec, 5);

julia> signal
5×15×2 Array{Float64, 3}:
[:, :, 1] =
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

[:, :, 2] =
 1.0  1.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
function add_responses!(signal, responses::Vector, e, s, tvec, trial_idx)
    @views signal[e, tvec, s] .+= responses[:, trial_idx]
end
function add_responses!(signal, responses::Matrix, e, s, tvec, trial_idx)
    @views signal[e, tvec, s] .+= responses[:, trial_idx]
end
function add_responses!(signal, responses::AbstractArray, e, s, tvec, trial_idx)
    @views signal[e, tvec, s] .+= responses[e, :, trial_idx]
end

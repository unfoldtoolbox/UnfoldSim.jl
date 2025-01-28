# Here we define some predefined simulations. Convenient if you just want to have a quick simulation :)

# If no RNG is given, create one.
predef_2x2(; kwargs...) = predef_2x2(MersenneTwister(1); kwargs...) # without rng always call same one
predef_eeg(; kwargs...) = predef_eeg(MersenneTwister(1); kwargs...) # without rng always call same one
predef_eeg(nsubjects::Int; kwargs...) = predef_eeg(MersenneTwister(1), nsubjects; kwargs...) # without rng always call same one


"""
    predef_eeg(; <keyword arguments>)
    predef_eeg(rng; <keyword arguments>)
    predef_eeg(rng, n_subjects; <keyword arguments>)

Simulate data for a P1/N1/P3 component complex.

Note that this function is mainly used for demonstration and internal testing purposes or whenever a quick simulation is needed. \n
In case `n_subjects` is defined - `MixedModelComponents` are generated (multi-subject simulation), else `LinearModelComponents` (single-subject simulation). \n
The most used keyword argument is: `return_epoched = true` which returns already epoched data. If you want epoched data without overlap, specify `onset = NoOnset()` and `return_epoched = true`.

# Keyword arguments

## Design
- `n_repeats = 100`: Number of times the experimental design is repeated. Only used in the single-subject case.
- `event_order_function = shuffle`: Random trial order. Use `event_order_function = (rng, x) -> x` to deactivate.
- `conditions = Dict(:condition => ["car", "face"], :continuous => range(-5, 5, length = 10))`: Conditions and covariates used in this predefined design.

## Component / Signal
- `sfreq = 100`: Sampling frequency.
- `p1 = (p100(; sfreq = sfreq), @formula(0 ~ 1), [5], Dict())`: P1 with amplitude 5; no effects.
- `n1 = (n170(; sfreq = sfreq), @formula(0 ~ 1 + condition), [5,-3], Dict())`: N1 with amplitude 5, dummy-coded condition effect (levels "car", "face") of -3.
- `p3 = (p300(; sfreq = sfreq), @formula(0 ~ 1 + continuous), [5,1], Dict())`: P3 with amplitude 5, continuous effect range [-5,5] with slope 1.

## Onset
- `overlap = (0.5,0.2)`, # convenient parameterization for the default `onset::UniformOnset`. (`offset`, `width`) in seconds. If you do not want any overlap, either use `onset=NoOnset()`, or put the offset to a value larger than the maximum used component length, e.g. `overlap=(1,0.2)`. Put the `width` to `0` to have no jitter between events.
- `onset = UniformOnset(; offset = sfreq * 0.5 * overlap[1], width = sfreq * 0.5 * overlap[2])`, 

## Noise
- `noiselevel = 0.2`.
- `noise = PinkNoise(; noiselevel = noiselevel)`.

## Other parameters
- `return_epoched = false`: If true, already epoched data is returned. Otherwise, continuous data is returned.

# Returns
- `(data, events)::Tuple{Array, DataFrame}`: 
    - `data` contains the simulated EEG data.
      Its dimensionality depends on the status of `return_epoched` and whether it's a single- or multisubject simulation (1D, 2D or 3D array):
      `continuous_time`, `continuous_time x n_subjects`, `times x size(design)[1]` or `times x size(design)[1] x n_subjects`. 

    - `events` contains the event combinations based on the experimental design.
      Each row corresponds to one combination of condition/covariate levels which is often equivalent to one stimulus or trial.


# Examples
```julia-repl
julia> data, events = UnfoldSim.predef_eeg();

julia> events
2000×3 DataFrame
  Row │ continuous  condition  latency 
      │ Float64     String     Int64   
──────┼────────────────────────────────
    1 │   2.77778   car             62
    2 │  -5.0       face           132
  ⋮   │     ⋮           ⋮         ⋮
 1999 │  -0.555556  face        120096
 2000 │  -2.77778   car         120154

julia> data
120199-element Vector{Float64}:
  0.31631798033146774
  0.40338935529989906
  0.46409775558165056
  0.5862082040156747
  ⋮
 -0.1879589005111152
 -0.3163314509311509
 -0.22230944464885682
 -0.01320095208877194
```

See also [`predef_2x2`](@ref).
"""
function predef_eeg(
    rng;
    # design
    n_repeats = 100,
    event_order_function = shuffle,

    # component / signal
    sfreq = 100,
    p1 = (p100(; sfreq = sfreq), @formula(0 ~ 1), [5], Dict()),
    n1 = (n170(; sfreq = sfreq), @formula(0 ~ 1 + condition), [5, -3], Dict()),
    p3 = (p300(; sfreq = sfreq), @formula(0 ~ 1 + continuous), [5, 1], Dict()),
    kwargs...,
)

    design =
        SingleSubjectDesign(;
            conditions = Dict(
                :condition => ["car", "face"],
                :continuous => range(-5, 5, length = 10),
            ),
            event_order_function = event_order_function,
        ) |> x -> RepeatDesign(x, n_repeats)

    return predef_eeg(rng, design, LinearModelComponent, [p1, n1, p3]; sfreq, kwargs...)
end

function predef_eeg(
    rng::AbstractRNG,
    design::AbstractDesign,
    T::Type{<:AbstractComponent},
    comps;
    sfreq = 100,
    # noise
    noiselevel = 0.2,
    noise = PinkNoise(; noiselevel = noiselevel),

    # onset
    overlap = (0.5, 0.2),
    onset = UniformOnset(; offset = sfreq * overlap[1], width = sfreq * overlap[2]), #put offset to 1 for no overlap. put width to 0 for no jitter
    kwargs...,
)

    components = []
    for c in comps
        append!(components, [T(c...)])
    end
    return simulate(rng, design, components, onset, noise; kwargs...)
end

function predef_eeg(
    rng::AbstractRNG,
    n_subjects;
    # design
    n_items = 100,
    event_order_function = shuffle,
    conditions = Dict(
        :condition => ["car", "face"],
        :continuous => range(-5, 5, length = 10),
    ),
    # component / signal
    sfreq = 100,
    p1 = (;
        s = p100(; sfreq = sfreq),
        f = @formula(0 ~ 1 + (1 | subject) + (1 | item)),
        β = [5],
        σs = Dict(:subject => [1], :item => [1]),
        c = Dict(),
    ),
    n1 = (;
        s = n170(; sfreq = sfreq),
        f = @formula(
            0 ~ 1 + condition + (1 + condition | subject) + (1 + condition | item)
        ),
        β = [5, -3],
        σs = Dict(:subject => [1, 1], :item => [0.5, 0.5]),
        c = Dict(),
    ),
    p3 = (;
        s = p300(; sfreq = sfreq),
        f = @formula(
            0 ~ 1 + continuous + (1 + continuous | subject) + (1 + continuous | item)
        ),
        β = [5, 1],
        σs = Dict(:subject => [1, 1], :item => [0.5, 0.5]),
        c = Dict(),
    ),
    kwargs...,
)


    design = MultiSubjectDesign(;
        n_subjects = n_subjects,
        n_items = n_items,
        items_between = conditions,
        event_order_function = event_order_function,
    )

    return predef_eeg(rng, design, MixedModelComponent, [p1, n1, p3]; sfreq, kwargs...)

end
"""

    predef_2x2(rng::AbstractRNG; <keyword arguments>)

Simulate data for a 2x2 design i.e. a design with two conditions with two levels each.

Note that this function is mainly used for demonstration and internal testing purposes or whenever a quick simulation is needed. \n

The most used keyword argument is: `return_epoched = true` which returns already epoched data. If you want epoched data without overlap, specify `onset = NoOnset()` and `return_epoched = true`. \n
Be careful if you modify `n_items` with `n_subjects = 1`, `n_items` has to be a multiple of 4 (or your equivalent conditions factorial, e.g. all combinations length).

In difference to `predef_EEG`, `predef_2x2` is sample based (no sampling rate to be specified), and also has a 2x2 design, instead of a 2-categorical, 1-continuous design.

# Keyword arguments

## Design
- `n_items = 100`: Number of items.
- `n_subjects = 1`: Number of subjects.
- `conditions = Dict(:A => ["a_small","a_big"], :B => ["b_tiny","b_large"])`: Experimental conditions with their levels.
- `event_order_function = shuffle`: Random trial order.

## Component / Signal
- `signalsize = 100`: Length of simulated hanning window.
- `basis = hanning(signalsize)`: The actual "function". `signalsize` is only used here.
- `β = [1, -0.5, .5, +1]`: The parameters for the fixed effects.
- `σs = Dict(:subject => [1, 0.5, 0.5, 0.5],:item => [1])`: Only needed in n_subjects >= 2 cases; specifies the random effects.
- `contrasts = Dict(:A => EffectsCoding(), :B => EffectsCoding())`: Effect coding by default.
- `formula = n_subjects == 1 ? @formula(0 ~ 1 + A*B) : @formula(dv ~ 1 + A*B + (A*B|subject) + (1|item))`: Model formula with interaction.

## Onset
- `overlap = (0.5,0.2)`,
- `onset = UniformOnset(; offset = signalsize * overlap[1], width = signalsize * overlap[2])`, # Put offset to 1 for no overlap. put width to 0 for no jitter

## Noise
- `noiselevel = 0.2`.
- `noise = PinkNoise(; noiselevel = noiselevel)`.

## Other parameters
- `return_epoched = false`: If true, already epoched data is returned. Otherwise, continuous data is returned.

# Returns
- `(data, events)::Tuple{Array, DataFrame}`: 
    - `data` contains the simulated EEG data.
      Its dimensionality depends on the status of `return_epoched` and whether n_subjects >=2:
      `continuous_time`, `continuous_time x n_subjects`, `times x size(design)[1]` or `times x size(design)[1] x n_subjects`. 

    - `events` contains the event combinations based on the experimental design.
      Each row corresponds to one combination of condition levels which is often equivalent to one stimulus or trial.

# Examples
```julia-repl
julia> data, events = UnfoldSim.predef_2x2();

julia> events
100×3 DataFrame
 Row │ A        B        latency 
     │ String   String   Int64   
─────┼───────────────────────────
   1 │ a_small  b_large       62
   2 │ a_big    b_tiny       132
  ⋮  │    ⋮        ⋮        ⋮
  99 │ a_big    b_large     5883
 100 │ a_big    b_tiny      5935
                  96 rows omitted

julia> data
6035-element Vector{Float64}:
  0.32384561187165956
  0.4108799249488322
  0.4715514814540277
  0.5936253009785152
  ⋮
  0.047664408295942984
  0.051193074728432035
 -0.08951617593287822
 -0.14456000097460356
```

See also [`predef_eeg`](@ref).
"""
function predef_2x2(
    rng::AbstractRNG;
    # design
    n_items = 100,
    n_subjects = 1,
    conditions = Dict(:A => ["a_small", "a_big"], :B => ["b_tiny", "b_large"]),
    event_order_function = shuffle,

    # component / signal
    signalsize = 100,
    basis = hanning(signalsize), # the component "function"
    β = [1, -0.5, 0.5, +1],
    σs = Dict(:subject => [1, 0.5, 0.5, 0.5], :item => [1]), # only in n_subjects>2 case
    contrasts = Dict(:A => EffectsCoding(), :B => EffectsCoding()),
    formula = n_subjects == 1 ? @formula(0 ~ 1 + A * B) :
              @formula(dv ~ 1 + A * B + (A * B | subject) + (1 | item)),

    # noise
    noiselevel = 0.2,
    noise = PinkNoise(; noiselevel = noiselevel),

    # onset
    overlap = (0.5, 0.2),
    onset = UniformOnset(;
        offset = signalsize * overlap[1],
        width = signalsize * overlap[2],
    ), #put offset to 1 for no overlap. put width to 0 for no jitter
    kwargs...,
)

    if n_subjects == 1
        design =
            SingleSubjectDesign(;
                conditions = conditions,
                event_order_function = event_order_function,
            ) |> x -> RepeatDesign(x, n_items ./ length(x))

        signal = LinearModelComponent(;
            basis = basis,
            formula = formula,
            β = β,
            contrasts = contrasts,
        )
    else
        design = MultiSubjectDesign(;
            n_subjects = n_subjects,
            n_items = n_items,
            items_between = conditions,
            event_order_function = event_order_function,
        )
        signal = MixedModelComponent(;
            basis = basis,
            formula = formula,
            β = β,
            σs = σs,
            contrasts = contrasts,
        )

    end

    data, events = simulate(rng, design, signal, onset, noise; kwargs...)
    return data, events
end

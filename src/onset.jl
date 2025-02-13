#----------------
# Types
#---------------

"""
    struct UniformOnset <: AbstractOnset
Provide a Uniform Distribution of the inter-event-distances.
`width`  is the width of the uniform distribution (=> the jitter). Since the lower bound is 0, `width` is also the upper bound.
`offset` is the minimal distance. The maximal distance is `offset + width`.

For a more advanced parameter specification, see `UniformOnsetFormula``, which allows to specify the onset-parameters depending on the `Design` employed via a linear regression model
"""
@with_kw struct UniformOnset <: AbstractOnset
    width = 50 # how many samples jitter?
    offset = 0 # minimal offset?
end
"""
    @with_kw struct LogNormalOnset <: AbstractOnset
Log-normal inter-event distances using the `Distributions.jl` truncated LogNormal distribution.

Be careful with large `μ` and `σ` values, as they are on logscale. σ>8 can quickly give you out-of-memory sized signals!

For a more advanced parameter specification, see `LogNormalOnsetFormula, which allows to specify the onset-parameters depending on the `Design` employed via linear regression model
"""
@with_kw struct LogNormalOnset <: AbstractOnset
    μ::Any  # mean
    σ::Any  # variance
    offset = 0 # additional offset
    truncate_upper = nothing # truncate at some sample?
end

"""
    struct NoOnset <: AbstractOnset end
In the case that the user directly wants no overlap to be simulated (=> epoched data).
"""
struct NoOnset <: AbstractOnset end

"""
    simulate_interonset_distances(rng, onset::UniformOnset, design::AbstractDesign)
    simulate_interonset_distances(rng, onset::LogNormalOnset, design::AbstractDesign)
    simulate_interonset_distances(rng, onset::UniformOnsetFormula, design::AbstractDesign)
    simulate_interonset_distances(rng, onset::LogNormalOnsetFormula, design::AbstractDesign)
Generate the inter-event-onset vector in samples (returns Int).
"""

function simulate_interonset_distances(rng, onset::UniformOnset, design::AbstractDesign)
    return Int.(
        round.(rand(deepcopy(rng), onset.offset:(onset.offset+onset.width), size(design)))
    )
end

function simulate_interonset_distances(rng, onset::LogNormalOnset, design::AbstractDesign)
    s = size(design)
    fun = LogNormal(onset.μ, onset.σ)
    if !isnothing(onset.truncate_upper)
        fun = truncated(fun; upper = onset.truncate_upper)
    end
    return Int.(round.(onset.offset .+ rand(deepcopy(rng), fun, s)))
end


#function simulate_interonset_distances(rng, onset::AbstractOnset,design::)


contains_design(d::AbstractDesign, target::Type) = false
contains_design(d::Union{RepeatDesign,SequenceDesign,SubselectDesign}, target::Type) =
    d.design isa target ? true : contains_design(d.design, target)


"""
    simulate_onsets(rng, onset::AbstractOnset, simulation::Simulation)
Call `simulate_interonset_distances` to generate distances between events and then add them up to generate the actual latencies in samples.
# main call from `simulation`
"""
function simulate_onsets(rng, onset::AbstractOnset, simulation::Simulation)

    # sample different onsets
    onsets = simulate_interonset_distances(rng, onset, simulation.design)


    if contains_design(simulation.design, SequenceDesign)
        currentsequence = sequencestring(deepcopy(rng), simulation.design)
        if !isnothing(findfirst("_", currentsequence))

            @assert currentsequence[end] == '_' "the blank-indicator '_' has to be the last sequence element"
            df = generate_events(simulation.design)
            nrows_df = size(df, 1)
            stepsize = length(currentsequence) - 1
            # add to every stepsize onset the maxlength of the response
            #@debug onsets[stepsize:stepsize:end]
            @debug stepsize
            onsets[stepsize+1:stepsize:end] .+= 2 .* maxlength(simulation.components)
            #@debug onsets[stepsize:stepsize:end]
        end
    end

    if maximum(onsets) > 10000
        @warn "Maximum of inter-event-distances was $(maximum(onsets)) - are you sure this is what you want?"
    end
    # accumulate them
    onsets_accum = accumulate(+, onsets, dims = 1, init = 1)

    return onsets_accum
end

"""
    UniformOnsetFormula <: AbstractOnset
provide a Uniform Distribution of the inter-event-distances, but with regression formulas.
This is helpful if your overlap/event-distribution should be dependend on some condition, e.g. more overlap in cond='A' than cond='B'.

**width**

        -`width_formula`: choose a formula depending on your `Design`, default `@formula(0~1)`
        -`width_β`: Choose a `Vector` of betas, number needs to fit the formula chosen, no default.
        -`width_contrasts` (optional): Choose a contrasts-`Dict`ionary according to the StatsModels specifications, default `Dict()``
    
**offset** is the minimal distance. The maximal distance is `offset + width`.

        -`offset_formula`: choose a formula depending on your `design`, default `@formula(0~1)``
        -`offset_β`: Choose a `Vector` of betas, number needs to fit the formula chosen, default `[0]`
        -`offset_contrasts` (optional): Choose a contrasts-`Dict`ionary according to the StatsModels specifications, default `Dict()`

See `UniformOnset` for a simplified version without linear regression specifications
"""
@with_kw struct UniformOnsetFormula <: AbstractOnset
    width_formula = @formula(0 ~ 1)
    width_β::Vector
    width_contrasts::Dict = Dict()
    offset_formula = @formula(0 ~ 1)
    offset_β::Vector = [0]
    offset_contrasts::Dict = Dict()
end


function simulate_interonset_distances(rng, o::UniformOnsetFormula, design::AbstractDesign)
    events = generate_events(design)
    widths =
        UnfoldSim.generate_designmatrix(o.width_formula, events, o.width_contrasts) *
        o.width_β
    offsets =
        UnfoldSim.generate_designmatrix(o.offset_formula, events, o.offset_contrasts) *
        o.offset_β

    return Int.(
        round.(reduce(vcat, rand.(deepcopy(rng), range.(offsets, offsets .+ widths), 1)))
    )
end


"""
    LogNormalOnsetFormula <: AbstractOnset
provide a LogNormal Distribution of the inter-event-distances, but with regression formulas.
This is helpful if your overlap/event-distribution should be dependend on some condition, e.g. more overlap in cond='A' than cond='B'.

**μ**

        -`μ_formula`: choose a formula depending on your `Design`, default `@formula(0~1)`
        -`μ_β`: Choose a `Vector` of betas, number needs to fit the formula chosen, default `[0]`
        -`μ_contrasts` (optional): Choose a contrasts-`Dict`ionary according to the StatsModels specifications, default `Dict()``
   
        -`σ_formula`: choose a formula depending on your `Design`, default `@formula(0~1)`
        -`σ_β`: Choose a `Vector` of betas, number needs to fit the formula chosen, default `[0]`
        -`σ_contrasts` (optional): Choose a contrasts-`Dict`ionary according to the StatsModels specifications, default `Dict()``
    
**offset** is the minimal distance. The maximal distance is `offset + width`.

        -`offset_formula`: choose a formula depending on your `design`, default `@formula(0~1)``
        -`offset_β`: Choose a `Vector` of betas, number needs to fit the formula chosen, default `[0]`
        -`offset_contrasts` (optional): Choose a contrasts-`Dict`ionary according to the StatsModels specifications, default `Dict()`

`truncate_upper` - truncate at some sample, default nothing

See `LogNormalOnset` for a simplified version without linear regression specifications
"""
@with_kw struct LogNormalOnsetFormula <: AbstractOnset
    μ_formula = @formula(0 ~ 1)
    μ_β::Vector
    μ_contrasts::Dict = Dict()
    σ_formula = @formula(0 ~ 1)
    σ_β::Vector
    σ_contrasts::Dict = Dict()
    offset_formula = @formula(0 ~ 1)
    offset_β::Vector = [0]
    offset_contrasts::Dict = Dict()
    truncate_upper = nothing # truncate at some sample?
end

function simulate_interonset_distances(
    rng,
    o::LogNormalOnsetFormula,
    design::AbstractDesign,
)
    events = generate_events(design)


    μs = UnfoldSim.generate_designmatrix(o.μ_formula, events, o.μ_contrasts) * o.μ_β
    σs = UnfoldSim.generate_designmatrix(o.σ_formula, events, o.σ_contrasts) * o.σ_β
    offsets =
        UnfoldSim.generate_designmatrix(o.offset_formula, events, o.offset_contrasts) *
        o.offset_β


    funs = LogNormal.(μs, σs)
    if !isnothing(o.truncate_upper)
        funs = truncated.(funs; upper = o.truncate_upper)
    end
    #@debug reduce(hcat, rand.(deepcopy(rng), funs, 1))
    return Int.(round.(offsets .+ reduce(vcat, rand.(deepcopy(rng), funs, 1))))
end
"""
    SequenceOnset <: AbstractOnset
A struct that defines as one argument the used Onsets for a [`SequenceDesign`](@ref).

All fields are mandatory. Works best with [`SequenceDesign`](@ref).

# Fields
- `onset::Dict`: for each Sequence event a Onset is defined.

# Examples
```julia-repl
sequence_onset = SequenceOnset(
    Dict('S'=>UniformOnset(width=0,offset=85*fs/100),
         'C'=>(DriftOnset(), UniformOnset(width=0, offset=150)),
         'R'=>UniformOnset(width=0,offset=120*fs/100)))
```
"""
struct SequenceOnset <: AbstractOnset
    onset::Dict
end

"""
    SequenceOnset <: AbstractOnset
A struct that defines the used Onsets for a [`Drift_Component`](@ref).

All fields are mandatory. Works best with [`Drift_Component`](@ref).

# Fields
- `onset::Dict`: onset for the Drift_Component.

# Examples
```julia-repl
sequence_onset = SequenceOnset(
    Dict('S'=>UniformOnset(width=0,offset=85*fs/100),
         'C'=>(DriftOnset(), UniformOnset(width=0, offset=150)),
         'R'=>UniformOnset(width=0,offset=120*fs/100)))
```
"""
struct DriftOnset{T} <: AbstractOnset
    onset::T
end
DriftOnset() = DriftOnset(UniformOnset(width=0,offset=0))

"""
    UnfoldSim.simulate_interonset_distances(rng, onset::AbstractOnset, design::AbstractDesign, components::AbstractComponent)

Base function to simulate Abstract Onset between components in an [`SequenceDesign`](@ref).

# Arguments
- `rng::StableRNG`: Random seed to ensure reproducibility.
- `onset::AbstractOnset`: Onset of type AbstractOnset which defines how the onset is created.
- `design::AbstractDesign`: Design for which the onsets are simulated.
- `components::AbstractComponent`: The Component for which the onset is simulated.

# Returns
- `simulate_interonset_distances`: function call.
"""
UnfoldSim.simulate_interonset_distances(rng, onset::AbstractOnset, design::AbstractDesign, components::AbstractComponent) = UnfoldSim.simulate_interonset_distances(rng,onset,design)

"""
    UnfoldSim.simulate_interonset_distances(rng, onset::DriftOnset, design::AbstractDesign, components::AbstractComponent)

Generates list of onsets for multiple [`Drift_Component`](@ref) in an [`SequenceDesign`](@ref).

# Arguments
- `rng::StableRNG`: Random seed to ensure reproducibility.
- `onset::DriftOnset`: DriftOnset defines to create onsets for a [`Drift_Component`](@ref).
- `design::AbstractDesign`: Design for which the onsets are simulated.
- `components::AbstractComponent`: The Component for which the onset is simulated.

# Returns
- `Vector{Float64}`: the generated onsets for the drift components in the SequenceDesign.
"""
function UnfoldSim.simulate_interonset_distances(rng, onset::DriftOnset, design::AbstractDesign, components::AbstractComponent)
    rts = calculate_response_times_for_ssm(deepcopy(rng), components, design)
    return Int.(round.(rts))
end

"""
    UnfoldSim.simulate_interonset_distances(rng, onset::Tuple{DriftOnset, UniformOnset}, design::AbstractDesign, components::AbstractComponent)

Generates list of onsets for multiple [`Drift_Component`](@ref) in an [`SequenceDesign`](@ref) and possibility to ad an [`UniformOnset`](@ref).

# Arguments
- `rng::StableRNG`: Random seed to ensure reproducibility.
- `onset::Tuple{DriftOnset, UniformOnset}`: DriftOnset defines to create onsets for a [`Drift_Component`](@ref) on top with an [`UniformOnset`](@ref).
- `design::AbstractDesign`: Design for which the onsets are simulated.
- `components::AbstractComponent`: The Component for which the onset is simulated.

# Returns
- `Vector{Float64}`: the generated onsets for the drift components in the SequenceDesign.
"""
function UnfoldSim.simulate_interonset_distances(rng, onset::Tuple{DriftOnset, UniformOnset}, design::AbstractDesign, components::AbstractComponent)
    rts = calculate_response_times_for_ssm(deepcopy(rng), components, design)
    jitter = Int.(
        round.(rand(deepcopy(rng), onset[2].offset:(onset[2].offset+onset[2].width), size(design)))
    )
    rts = Int.(round.(rts))
    rts = rts .+ jitter
    return rts
end

"""
    UnfoldSim.simulate_interonset_distances(rng, onset::Char, design::AbstractDesign, components::AbstractComponent)

Generates list of onsets for the end of a sequence in an [`SequenceDesign`](@ref).

# Arguments
- `rng::StableRNG`: Random seed to ensure reproducibility.
- `onset::Char`: Defines to simulates onsets at the end of a sequence.
- `design::AbstractDesign`: Design for which the onsets are simulated.
- `components::AbstractComponent`: The Component for which the onset is simulated.

# Returns
- `Vector{Float64}`: the generated onsets for the end of a sequence in the SequenceDesign.
"""
function UnfoldSim.simulate_interonset_distances(rng, onset::Char, design::AbstractDesign, components::AbstractComponent)
    @assert onset == '_'
    df = generate_events(rng, design)
    nrows_df = Int(size(df, 1))
    onsets = repeat([UnfoldSim.maxlength([components])],nrows_df)
    return onsets
end

"""
    UnfoldSim.simulate_onsets(rng, onset::SequenceOnset, simulation::Simulation)

Generates list of onsets for all events of an [`SequenceDesign`](@ref), how to simulate the onsets is defined in the [`SequenceOnset`](@ref).

# Arguments
- `rng::StableRNG`: Random seed to ensure reproducibility.
- `onset::SequenceOnset`: onset definition for each event in the sequence design.
- `simulation::Simulation`: Simulation which contains the design and other elements for the experiment.

# Returns
- `Vector{Float64}`: the generated onsets for all events in the SequenceDesign.
"""
function UnfoldSim.simulate_onsets(rng, onset::SequenceOnset, simulation::Simulation)
    @assert isa(simulation.design.design, SequenceDesign)
    events = generate_events(deepcopy(rng), simulation.design)
    onset_map = Dict()
    onset_counter = Dict()
    for k in keys(onset.onset)
        sub_design = UnfoldSim.SubselectDesign(simulation.design, k)
        onsets_for_k = simulate_interonset_distances(deepcopy(rng), onset.onset[k], sub_design, simulation.components[k][1])
        onset_map[k] = onsets_for_k
        onset_counter[k] = 1
    end
    final_onsets = []
    for (i, evt_k) in enumerate(events.event)
        push!(final_onsets, onset_map[evt_k][onset_counter[evt_k]])
        onset_counter[evt_k] += 1
    end
    final_onsets = vcat(final_onsets[end], final_onsets[1:end-1])
    if maximum(final_onsets) > 10000
        @warn "Maximum of inter-event-distances was $(maximum(final_onsets)) - are you sure this is what you want?"
    end
    onsets_accum = accumulate(+, final_onsets, dims = 1, init = 1)
    return onsets_accum

end
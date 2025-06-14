#----------------
# Types
#----------------

"""
    UniformOnset <: AbstractOnset

Provide a Uniform Distribution for the inter-event distances (in samples).

Tip: To manually generate inter-event distance samples use the [`simulate_interonset_distances`](@ref) function.

# Fields
- `width = 50` (optional): Width of the uniform distribution (=> the "jitter"). Since the lower bound is 0, `width` is also the upper bound.
- `offset = 0` (optional): The minimal distance between events. The maximal distance is `offset + width`.

# Examples
```julia-repl
julia> onset_distribution = UniformOnset(width = 25, offset = 5)
UniformOnset
  width: Int64 25
  offset: Int64 5
```

See also [`LogNormalOnset`](@ref), [`NoOnset`](@ref).
"""
@with_kw struct UniformOnset <: AbstractOnset
    width = 50 # how many samples jitter?
    offset = 0 # minimal offset?
end
"""
    LogNormalOnset <: AbstractOnset

Log-normal inter-event distances (in samples) using the `Distributions.jl` truncated LogNormal distribution ([code and mathematical reference](https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.LogNormal)).

Be careful with large `μ` and `σ` values, as they are on logscale. σ>8 can quickly give you out-of-memory sized signals! \\
Tip: To manually generate inter-event distance samples use the [`simulate_interonset_distances`](@ref) function.

# Fields
- `μ`: The mean of the log-transformed variable (the log-normal random variable's logarithm follows a normal distribution).
- `σ`: The standard deviation of the log-transformed variable.
- `offset = 0` (optional): The minimal distance between events.
- `truncate_upper = nothing` (optional): Upper limit (in samples) at which the distribution is truncated.

# Examples
```julia-repl
julia> onset_distribution = LogNormalOnset(3, 0.25, 10, 25)
LogNormalOnset
  μ: Int64 3
  σ: Float64 0.25
  offset: Int64 10
  truncate_upper: Int64 25
```

See also [`UniformOnset`](@ref), [`NoOnset`](@ref).
"""
@with_kw struct LogNormalOnset <: AbstractOnset
    μ::Any  # mean
    σ::Any  # variance
    offset = 0 # additional offset
    truncate_upper = nothing # truncate at some sample?
end

"""
    NoOnset <: AbstractOnset

For cases where the user wants to simulate epoched data without any overlap between consecutive events.

# Examples
```julia-repl
julia> onset_distribution = NoOnset()
NoOnset()
```

See also [`UniformOnset`](@ref), [`LogNormalOnset`](@ref).
"""
struct NoOnset <: AbstractOnset end


#-----------------------------
# Onset simulation functions
#-----------------------------

"""
    simulate_interonset_distances(rng, onset::AbstractOnset, design::AbstractDesign)

Generate the inter-onset distance vector by sampling from the respective distribution (in samples).

# Arguments
- `rng`: Random number generator (RNG) to make the process reproducible.
- `onset::AbstractOnset`: Inter-onset distance distribution to sample from.
- `design::AbstractDesign`: Experimental design with conditions and covariates.

# Returns
- `Vector{Integer}`: Inter-onset distances in samples. Note that these are distances between onsets and no latencies.

# Examples
```julia-repl
# Create an experimental design
julia> design_single = SingleSubjectDesign(;
           conditions = Dict(
               :stimulus_type => ["natural", "artificial"],
               :contrast_level => range(0, 1, length = 3),
           ),
       );

# Create an inter-onset distance distribution
julia> onset_distribution = LogNormalOnset(3, 0.5, 5, nothing);

julia> using StableRNGs

julia> simulate_interonset_distances(StableRNG(1), onset_distribution, design_single)
6-element Vector{Int64}:
 20
 26
 34
 18
 12
 23
```

See also [`simulate_onsets`](@ref).
"""
function simulate_interonset_distances end


function simulate_interonset_distances(rng, onset::UniformOnset, design::AbstractDesign)
    return Int.(
        round.(
            rand(
                deepcopy(rng),
                onset.offset:(onset.offset+onset.width),
                size(deepcopy(rng), design),
            )
        )
    )
end

function simulate_interonset_distances(rng, onset::LogNormalOnset, design::AbstractDesign)
    s = size(deepcopy(rng), design)
    fun = LogNormal(onset.μ, onset.σ)
    if !isnothing(onset.truncate_upper)
        fun = truncated(fun; upper = onset.truncate_upper)
    end
    return Int.(round.(onset.offset .+ rand(deepcopy(rng), fun, s)))
end


contains_design(d::AbstractDesign, target::Type) = false
contains_design(d::Union{RepeatDesign,SequenceDesign,SubselectDesign}, target::Type) =
    d.design isa target ? true : contains_design(d.design, target)


"""
    simulate_onsets(rng, onset::AbstractOnset, simulation::Simulation)

Call `simulate_interonset_distances` to generate distances between events and then add them up to generate the actual latencies in samples.

Please note that this function is mainly for internal use in the context of `simulate` function calls. \n
Also note that the accumulation of onsets starts at 1 to avoid indexing problems in the case that the first sampled onset is 0.

In case of a SequenceDesign with a '_' no-overlap indicator, we use twice the `maxlength(components)` as the distance following that sequence character.

# Arguments
- `rng`: Random number generator (RNG) to make the process reproducible.
- `onset::AbstractOnset`: Inter-onset distance distribution which is passed to `simulate_interonset_distances`.
- `simulation::Simulation`: Simulation object which contains design, component(s), inter-onset distance distribution and noise.

# Returns

# Examples
```julia-repl
# Create Simulation object
julia> design_single = SingleSubjectDesign(;
           conditions = Dict(
               :stimulus_type => ["natural", "artificial"],
               :contrast_level => range(0, 1, length = 3),
           ),
       );

julia> p1_component = LinearModelComponent(; basis = p100(), formula = @formula(0 ~ 1), β = [5]);

julia> simulation = Simulation(design_single, p1_component, UniformOnset(), NoNoise());

julia> using StableRNGs

# Simulate onsets for this simulation
julia> simulate_onsets(StableRNG(1), simulation.onset, simulation)
6-element Vector{Int64}:
  20
  70
  97
 110
 150
 182
```

See also [`simulate_interonset_distances`](@ref).
"""
function simulate_onsets(rng, onset::AbstractOnset, simulation::Simulation)

    # sample different onsets
    onsets = simulate_interonset_distances(rng, onset, simulation.design)


    if contains_design(simulation.design, SequenceDesign)
        currentsequence = sequencestring(deepcopy(rng), simulation.design)
        if !isnothing(findfirst("_", currentsequence))

            @assert currentsequence[end] == '_' "the blank-indicator '_' has to be the last sequence element"
            df = generate_events(deepcopy(rng), simulation.design)
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
        @warn "Maximum of inter-event distances was $(maximum(onsets)) - are you sure this is what you want?"
    end
    # accumulate them
    onsets_accum = accumulate(+, onsets, dims = 1, init = 1)
    onsets_accum = onsets_accum .- minoffset(simulation.components)

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
    events = generate_events(deepcopy(rng), design)
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
    events = generate_events(deepcopy(rng), design)


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
Allows to specify inter-onset-distance functions per event.

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
SequenceOnset(args::Pair...) = SequenceOnset(Dict(args...))
"""
    DriftOnset<: AbstractOnset
A type that defines the onsets for a [`DriftComponent`](@ref).

Works best with [`DriftComponent`](@ref).

# Fields
- `onset::Dict`: onset for the DriftComponent.

# Examples
```julia-repl
drift_onset = DriftOnset()
```
"""
struct DriftOnset{T} <: AbstractOnset
    onset::T
end
DriftOnset() = DriftOnset(UniformOnset(width = 0, offset = 0))

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

Generates list of onsets for multiple [`DriftComponent`](@ref) in an [`SequenceDesign`](@ref), one for each sequence. Onsets are rounded, returned as `Int`

# Arguments
- `rng::StableRNG`: Random seed to ensure reproducibility.
- `onset::DriftOnset`: DriftOnset defines to create onsets for a [`DriftComponent`](@ref).
- `design::AbstractDesign`: Design for which the onsets are simulated.
- `components::AbstractComponent`: The Component for which the onset is simulated.

# Returns
- `Vector{Int}`: the generated onsets for the drift components in the SequenceDesign.
"""
function UnfoldSim.simulate_interonset_distances(rng, onset::DriftOnset, design::AbstractDesign, components::AbstractComponent)
    rts = calculate_response_times_for_ssm(deepcopy(rng), components, design)
    return Int.(round.(rts))
end

"""
    UnfoldSim.simulate_interonset_distances(rng, onset::Tuple{DriftOnset, UniformOnset}, design::AbstractDesign, components::AbstractComponent)

Generates list of onsets for multiple [`DriftComponent`](@ref) in an [`SequenceDesign`](@ref) and possibility to ad an [`UniformOnset`](@ref).

# Arguments
- `rng::StableRNG`: Random seed to ensure reproducibility.
- `onset::Tuple{DriftOnset, UniformOnset}`: DriftOnset defines to create onsets for a [`DriftComponent`](@ref) on top with an [`UniformOnset`](@ref).
- `design::AbstractDesign`: Design for which the onsets are simulated.
- `components::AbstractComponent`: The Component for which the onset is simulated.

# Returns
- `Vector{Float64}`: the generated onsets for the drift components in the SequenceDesign.
"""
function UnfoldSim.simulate_interonset_distances(rng, onset::Tuple, design::AbstractDesign, components::AbstractComponent)
    return Int.(round.(reduce(.+, simulate_interonset_distances.(deepcopy(rng), onset, Ref(design), Ref(components)))))
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
    onsets = repeat([UnfoldSim.maxlength([components])], nrows_df)
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
    @assert isa(simulation.design.design, SequenceDesign) "`SequenceOnset`is currently only compatible with a `SequenceDesign`"
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
        @warn "Number of simulated inter-event-distances was $(maximum(final_onsets)) - are you sure this is what you want?"
    end
    onsets_accum = accumulate(+, final_onsets, dims = 1, init = 1)
    return onsets_accum

end

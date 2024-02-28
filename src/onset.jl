#----------------
# Types
#---------------

"""
    struct UniformOnset <: AbstractOnset
Provide a Uniform Distribution of the inter-event-distances.
`width`  is the width of the uniform distribution (=> the jitter). Since the lower bound is 0, `width` is also the upper bound.
`offset` is the minimal distance. The maximal distance is `offset + width`.

For a more advanced parameter specification, see `FormulaUniformOnset``, which allows to specify the onset-parameters depending on the `Design` employed via a linear regression model
"""
@with_kw struct UniformOnset <: AbstractOnset
    width = 50 # how many samples jitter?
    offset = 0 # minimal offset?
end
"""
    @with_kw struct LogNormalOnset <: AbstractOnset
Log-normal inter-event distances using the `Distributions.jl` truncated LogNormal distribution.

Be careful with large `μ` and `σ` values, as they are on logscale. σ>8 can quickly give you out-of-memory sized signals!

For a more advanced parameter specification, see `FormulaLogNormalOnset, which allows to specify the onset-parameters depending on the `Design` employed via linear regression model
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
    simulate_interonset_distances(rng, onset::FormulaUniformOnset, design::AbstractDesign)
    simulate_interonset_distances(rng, onset::FormulaLogNormalOnset, design::AbstractDesign)
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
            onsets[stepsize+1:stepsize:end] .= 2 .* maxlength(simulation.components)
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
    FormulaUniformOnset <: AbstractOnset
provide a Uniform Distribution of the inter-event-distances, but with regression formulas.
This is helpful if your overlap/event-distribution should be dependend on some condition, e.g. more overlap in cond='A' than cond='B'.

    - `width`:  is the width of the uniform distribution (=> the jitter). Since the lower bound is 0, `width` is also the upper bound.
        -`width_formula`: choose a formula depending on your `Design`
        -`width_β`: Choose a vector of betas, number needs to fit the formula chosen
        -`width_contrasts` (optional): Choose a contrasts-dictionary according to the StatsModels specifications
    `offset` is the minimal distance. The maximal distance is `offset + width`.
        -`offset_formula`: choose a formula depending on your `design`
        -`offset_β`: Choose a vector of betas, number needs to fit the formula chosen
        -`offset_contrasts` (optional): Choose a contrasts-dictionary according to the StatsModels specifications

See `UniformOnset` for a simplified version without linear regression specifications
"""
@with_kw struct FormulaUniformOnset <: AbstractOnset
    width_formula = @formula(0 ~ 1)
    width_β::Vector = [50]
    width_contrasts::Dict = Dict()
    offset_formula = @formula(0 ~ 1)
    offset_β::Vector = [0]
    offset_contrasts::Dict = Dict()
end


function simulate_interonset_distances(rng, o::FormulaUniformOnset, design::AbstractDesign)
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


@with_kw struct FormulaLogNormalOnset <: AbstractOnset
    μ_formula = @formula(0 ~ 1)
    μ_β::Vector = [0]
    μ_contrasts::Dict = Dict()
    σ_formula = @formula(0 ~ 1)
    σ_β::Vector = [0]
    σ_contrasts::Dict = Dict()
    offset_formula = @formula(0 ~ 1)
    offset_β::Vector = [0]
    offset_contrasts::Dict = Dict()
    truncate_upper = nothing # truncate at some sample?
end

function simulate_interonset_distances(
    rng,
    o::FormulaLogNormalOnset,
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
        fun = truncated.(fun; upper = o.truncate_upper)
    end
    return Int.(round.(offsets .+ rand.(deepcopy(rng), funs, 1)))
end



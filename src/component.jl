"""
A component that adds a hierarchical relation between parameters according to a LMM defined via MixedModels.jl

- `basis`: an object, if accessed, provides a 'basis-function', e.g. `hanning(40)`, this defines the response at a single event. It will be weighted by the model-prediction
- `formula`: Formula-Object in the style of MixedModels.jl e.g. `@formula 0~1+cond + (1|subject)` - left-handside is ignored
- `β` Vector of betas, must fit the formula
- `σs` Dict of random effect variances, e.g. `Dict(:subject=>[0.5,0.4])` or to specify correlationmatrix `Dict(:subject=>[0.5,0.4,I(2,2)],...)`. Technically, this will be passed to MixedModels.jl `create_re` function, which creates the θ matrices.
- `contrasts`: Dict in the style of MixedModels.jl. Default is empty.

All arguments can be named, in that case `contrasts` is optional

Works best with `MultiSubjectDesign`
```julia
MixedModelComponent(;
    basis=hanning(40),
    formula=@formula(0~1+cond+(1+cond|subject)),
    β = [1.,2.],
    σs= Dict(:subject=>[0.5,0.4]),
    contrasts=Dict(:cond=>EffectsCoding())
)

```
"""
@with_kw struct MixedModelComponent <: AbstractComponent
    basis::Any
    formula::Any # e.g. 0~1+cond 
    β::Vector
    σs::Dict # Dict(:subject=>[0.5,0.4]) or to specify correlationmatrix Dict(:subject=>[0.5,0.4,I(2,2)],...)
    contrasts::Dict = Dict()
end

"""
A multiple regression component for one subject

- `basis`: an object, if accessed, provides a 'basis-function', e.g. `hanning(40)`, this defines the response at a single event. It will be weighted by the model-prediction
- `formula`: StatsModels Formula-Object  `@formula 0~1+cond` (left side must be 0)
- `β` Vector of betas, must fit the formula
- `contrasts`: Dict. Default is empty, e.g. `Dict(:condA=>EffectsCoding())`

All arguments can be named, in that case `contrasts` is optional

Works best with `SingleSubjectDesign`
```julia
LinearModelComponent(;
    basis=hanning(40),
    formula=@formula(0~1+cond),
    β = [1.,2.],
    contrasts=Dict(:cond=>EffectsCoding())
)

```
"""
@with_kw struct LinearModelComponent <: AbstractComponent
    basis::Any
    formula::Any # e.g. 0~1+cond - left side must be "0"
    β::Vector
    contrasts::Dict = Dict()
end


"""
Wrapper for an `AbstractComponent` to project it to multiple target-channels via `projection`. optional adds `noise` to the source prior to projection.
"""
@with_kw struct MultichannelComponent <: AbstractComponent
    component::AbstractComponent
    projection::AbstractVector
    noise::AbstractNoise # optional
end

MultichannelComponent(c::AbstractComponent, p) =
    MultichannelComponent(c::AbstractComponent, p, NoNoise())

function MultichannelComponent(
    component::AbstractComponent,
    projection::Pair{<:AbstractHeadmodel,String},
    noise::AbstractNoise,
)
    ix = closest_src(projection[1], projection[2])
    mg = magnitude(projection[1])
    return MultichannelComponent(component, mg[:, ix], noise)
end
Base.length(c::MultichannelComponent) = length(c.component)

"""
    n_channels(c::AbstractComponent)
Return the number of channels. By default = 1.
"""
n_channels(c::AbstractComponent) = 1

"""
    n_channels(c::MultichannelComponent)
For `MultichannelComponent` return the length of the projection vector.

"""
n_channels(c::MultichannelComponent) = length(c.projection)

"""
For a vector of `MultichannelComponent`s, return the first but asserts all are of equal length.
"""
function n_channels(c::Vector{<:AbstractComponent})
    all_channels = n_channels.(c)
    @assert length(unique(all_channels)) == 1 "Error - projections of different channels cannot be different from eachother"
    return all_channels[1]
end

"""
    simulate_component(rng,c::MultichannelComponent,design::AbstractDesign)
Return the projection of a component from source to "sensor" space.
"""
function simulate_component(rng, c::MultichannelComponent, design::AbstractDesign)
    y = simulate_component(rng, c.component, design)

    for trial = 1:size(y, 2)
        y[:, trial] .= y[:, trial] .+ simulate_noise(rng, c.noise, size(y, 1))
    end

    y_proj = kron(y, c.projection)
    return reshape(y_proj, length(c.projection), size(y)...)
end



Base.length(c::AbstractComponent) = length(c.basis)

"""
    maxlength(c::Vector{<:AbstractComponent}) = maximum(length.(c))
maximum of individual component lengths
"""
maxlength(c::Vector{<:AbstractComponent}) = maximum(length.(c))

"""
    simulate_component(rng, c::AbstractComponent, simulation::Simulation)
By default call `simulate_component` with `(::Abstractcomponent,::AbstractDesign)` instead of the whole simulation. This allows users to provide a hook to do something completely different :)
"""
simulate_component(rng, c::AbstractComponent, simulation::Simulation) =
    simulate_component(rng, c, simulation.design)

"""
    simulate_component(rng, c::AbstractComponent, simulation::Simulation)
Generate a linear model design matrix, weight it by c.β and multiply the result with the given basis vector.

julia> c = UnfoldSim.LinearModelComponent([0,1,1,0],@formula(0~1+cond),[1,2],Dict())
julia> design = MultiSubjectDesign(;n_subjects=2,n_items=50,items_between=(;:cond=>["A","B"]))
julia> simulate_component(StableRNG(1),c,design)
"""
function simulate_component(rng, c::LinearModelComponent, design::AbstractDesign)
    events = generate_events(design)

    # special case, intercept only 
    # https://github.com/JuliaStats/StatsModels.jl/issues/269
    if c.formula.rhs == ConstantTerm(1)
        X = ones(nrow(events), 1)
    else
        if isempty(c.contrasts)
            m = StatsModels.ModelFrame(c.formula, events)
        else
            m = StatsModels.ModelFrame(c.formula, events; contrasts = c.contrasts)
        end
        X = StatsModels.modelmatrix(m)
    end
    y = X * c.β
    return y' .* c.basis
end
"""
    simulate_component(rng, c::MixedModelComponent, design::AbstractDesign)
Generates a MixedModel and simulates data according to c.β and c.σs.

A trick is used to remove the Normal-Noise from the MixedModel which might lead to rare numerical instabilities. Practically, we upscale the σs by factor 10000, and provide a σ=0.0001. Internally this results in a normalization where the response scale is 10000 times larger than the noise.

Currently, it is not possible to use a different basis for fixed and random effects, but a code-stub exists (it is slow though).

julia> design = MultiSubjectDesign(;n_subjects=2,n_items=50,items_between=(;:cond=>["A","B"]))
julia> c = UnfoldSim.MixedModelComponent([0.,1,1,0],@formula(0~1+cond+(1|subject)),[1,2],Dict(:subject=>[2],),Dict())
julia> simulate(StableRNG(1),c,design)

"""
function simulate_component(rng, c::MixedModelComponent, design::AbstractDesign)
    events = generate_events(design)

    # add the mixed models lefthandside
    lhs_column = :tmp_dv
    @assert string(lhs_column) ∉ names(events) "Error: Wow you are unlucky, we have to introduce a temporary lhs-symbol which we name ``:tmp_dv` - you seem to have a condition called `:tmp_dv` in your dataset as well. Please rename it!"
    f = FormulaTerm(Term(:tmp_dv), c.formula.rhs)
    events[!, lhs_column] .= 0

    # create dummy
    if isempty(c.contrasts)
        m = MixedModels.MixedModel(f, events)
    else
        m = MixedModels.MixedModel(f, events; contrasts = c.contrasts)
    end


    # empty epoch data
    epoch_data_component = zeros(Int(length(c.basis)), length(design))

    # residual variance for lmm
    σ_lmm = 0.0001
    if 1 == 1
        named_random_effects = weight_σs(c.σs, 1.0, σ_lmm)
        θ = createθ(m; named_random_effects...)
        @debug named_random_effects, θ, m.θ
        try
            simulate!(deepcopy(rng), m.y, m; β = c.β, σ = σ_lmm, θ = θ)
        catch e
            if isa(e, DimensionMismatch)
                @warn "Most likely your σs's do not match the formula!"
            elseif isa(e, ArgumentError)
                @warn "Most likely your β's do not match the formula!"
            end
            rethrow(e)
        end

        # save data to array
        #@show size(m.y)
        #@show size(c.basis)


        epoch_data_component = kron(c.basis, m.y')


    else
        # iterate over each timepoint
        for t in eachindex(c.basis)

            # select weight from basis
            # right now, it is the same, but maybe changein thefuture?
            basis_β = c.basis[t]
            basis_σs = c.basis[t]


            # weight random effects by the basis function
            named_random_effects = weight_σs(c.σs, basis_σs, σ_lmm)

            θ = createθ(m; named_random_effects...)


            # simulate with new parameters; will update m.y
            simulate!(deepcopy(rng), m.y, m; β = basis_β .* c.β, σ = σ_lmm, θ = θ)

            # save data to array
            epoch_data_component[t, :] = m.y
        end
    end
    return epoch_data_component

end


"""
Weights a σs Dict for MixedModels.jl by a Float64

Finally sales it by σ_lmm, as a trick to simulate noise-free LMMs

I anticipate a function
    `function weight_σs(σs::Dict,b_σs::Dict,σ_lmm::Float64)`
where each σs entry can be weighted individually
"""
function weight_σs(σs::Dict, b_σs::Float64, σ_lmm::Float64)
    #k = (collect(keys(σs))...,)
    #val = values(σs)

    keys = Symbol[]
    vals = LowerTriangular[]

    for (k, v) in σs

        scale = (x) -> b_σs ./ σ_lmm .* x

        if v[end] isa Matrix
            v = create_re.(scale(v[1:end-1])...; corrmat = v[end])
        else
            v = create_re.(scale(v)...;)
        end

        push!(keys, k)
        push!(vals, v)
    end

    named_random_effects = NamedTuple(keys .=> vals)

    return named_random_effects
end

#----

"""
    simulate_responses(
        rng,
        components::Vector{<:AbstractComponent},
        simulation::Simulation)
Simulate multiple component responses and accumulates them on a per-event basis.
"""
function simulate_responses(
    rng,
    components::Vector{<:AbstractComponent},
    simulation::Simulation,
)
    if n_channels(components) > 1
        epoch_data =
            zeros(n_channels(components), maxlength(components), length(simulation.design))
    else
        epoch_data = zeros(maxlength(components), length(simulation.design))
    end

    for c in components
        simulate_and_add!(epoch_data, c, simulation, rng)
    end
    return epoch_data
end


"""
    simulate_and_add!(epoch_data::AbstractMatrix, c, simulation, rng)
    simulate_and_add!(epoch_data::AbstractArray, c, simulation, rng)
Helper function to call `simulate_component` and add it to a provided Array.
"""
function simulate_and_add!(epoch_data::AbstractMatrix, c, simulation, rng)
    @debug "matrix"
    @views epoch_data[1:length(c), :] .+= simulate_component(rng, c, simulation)
end
function simulate_and_add!(epoch_data::AbstractArray, c, simulation, rng)
    @debug "3D Array"
    @views epoch_data[:, 1:length(c), :] .+= simulate_component(rng, c, simulation)
end

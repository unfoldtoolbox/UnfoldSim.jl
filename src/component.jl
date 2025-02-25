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
# backwards compatability after introducing the `offset` field`
@with_kw struct MixedModelComponent <: AbstractComponent
    basis::Any
    formula::Any # e.g. 0~1+cond 
    β::Vector
    σs::Dict # Dict(:subject=>[0.5,0.4]) or to specify correlationmatrix Dict(:subject=>[0.5,0.4,I(2,2)],...)
    contrasts::Dict = Dict()
    offset::Int = 0
end
MixedModelComponent(basis, formula, β, σs, contrasts) =
    MixedModelComponent(basis, formula, β, σs, contrasts, 0)
"""
A multiple regression component for one subject

- `basis`: an object, if accessed, provides a 'basis-function', e.g. `hanning(40)`, this defines the response at a single event. It will be weighted by the model-prediction. Can also be a tuple `(fun::Function,maxlength::Int)` with a function `fun` that either generates a matrix `size = (maxlength,size(design,1))` or a vector of vectors. If a larger matrix is generated, it is automatically cutoff at `maxlength`
- `formula`: StatsModels Formula-Object  `@formula 0~1+cond` (left side must be 0)
- `β` Vector of betas, must fit the formula
- `contrasts`: Dict. Default is empty, e.g. `Dict(:condA=>EffectsCoding())`
- `offset`: Int. Default is 0. Can be used to shift the basis function in time
All arguments can be named, in that case `contrasts` and  `offset` are optional.

Works best with `SingleSubjectDesign`
```julia
# use a hanning window of size 40 as the component basis
LinearModelComponent(;
    basis=hanning(40),
    formula=@formula(0~1+cond),
    β = [1.,2.],
    contrasts=Dict(:cond=>EffectsCoding())
)

# define a function returning random numbers as the component basis
maxlength = 15
my_signal_function = d->rand(StableRNG(1),maxlength,length(d))
LinearModelComponent(;
    basis=(my_signal_function,maxlength),
    formula=@formula(0~1),
    β = [1.],
)

```
"""
@with_kw struct LinearModelComponent <: AbstractComponent
    basis::Union{Tuple{Function,Int},Array}
    formula::FormulaTerm # e.g. 0~1+cond - left side must be "0"
    β::Vector
    contrasts::Dict = Dict()
    offset::Int = 0
    function LinearModelComponent(basis, formula, β, contrasts, offset)
        @assert isa(basis, Tuple{Function,Int}) ".basis needs to be an `::Array` or a `Tuple(function::Function,maxlength::Int)`"
        @assert basis[2] > 0 "`maxlength` needs to be longer than 0"
        new(basis, formula, β, contrasts, offset)
    end
    LinearModelComponent(basis::Array, formula, β, contrasts, offset) =
        new(basis, formula, β, contrasts, offset)
end

# backwards compatability after introducing the `offset` field
LinearModelComponent(basis, formula, β, contrasts) =
    LinearModelComponent(basis, formula, β, contrasts, 0)
"""
    offset(AbstractComponent)

Should the `basis` be shifted? Returns c.offset for most components, if not implemented for a type, returns 0. Can be positive or negative, but has to be Integer
"""
offset(c::AbstractComponent)::Int = 0
offset(c::LinearModelComponent)::Int = c.offset
offset(c::MixedModelComponent)::Int = c.offset

maxoffset(c::Vector{<:AbstractComponent}) = maximum(offset.(c))
maxoffset(d::Dict{<:Char,<:Vector{<:AbstractComponent}}) = maximum(maxoffset.(values(d)))
minoffset(c::Vector{<:AbstractComponent}) = minimum(offset.(c))
minoffset(d::Dict{<:Char,<:Vector{<:AbstractComponent}}) = minimum(minoffset.(values(d)))



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
    @assert length(unique(all_channels)) == 1 "Error - projections of different components have to be of the same output (=> channel) dimension"
    return all_channels[1]
end

function n_channels(components::Dict)
    all_channels = [n_channels(c) for c in values(components)]
    @assert length(unique(all_channels)) == 1 "Error - projections of different components have to be of the same output (=> channel) dimension"
    return all_channels[1]

end
"""
    simulate_component(rng,c::MultichannelComponent,design::AbstractDesign)
Return the projection of a component from source to "sensor" space.
"""
function simulate_component(rng, c::MultichannelComponent, design::AbstractDesign)
    y = simulate_component(rng, c.component, design)

    for trial = 1:size(y, 2)
        y[:, trial] .= y[:, trial] .+ simulate_noise(deepcopy(rng), c.noise, size(y, 1))
    end

    y_proj = kron(y, c.projection)
    return reshape(y_proj, length(c.projection), size(y)...)
end

"""
    basis(c::AbstractComponent)

returns the basis of the component (typically `c.basis`)
"""
basis(c::AbstractComponent) = c.basis

"""
    basis(c::AbstractComponent,design)
evaluates the basis, if basis is a vector, directly returns it. if basis is a tuple `(f::Function,maxlength::Int)`, evaluates the function with input `design`. Cuts the resulting vector or Matrix at `maxlength`
"""
basis(c::AbstractComponent, design) = basis(basis(c), design)


basis(b::AbstractVector, design) = b
function basis(basis::Tuple{Function,Int}, design)
    f = basis[1]
    maxlength = basis[2]
    basis_out = f(design)
    if isa(basis_out, AbstractVector{<:AbstractVector}) || isa(basis_out, AbstractMatrix)
        if isa(basis_out, AbstractMatrix)
            l = size(basis_out, 2)
        else
            l = length(basis_out) # vector of vector case
        end
        @assert l == size(generate_events(design))[1] "Component basis function needs to either return a Vector of vectors or a Matrix with dim(2) == size(design,1) [l / $(size(design,1))], or a Vector of Vectors with length(b) == size(design,1) [$l / $(size(design,1))]. "
    end
    limit_basis(basis_out, maxlength)
end


function limit_basis(b::AbstractVector{<:AbstractVector}, maxlength)

    # first cut off maxlength
    b = limit_basis.(b, maxlength)
    # now fill up with 0's
    Δlengths = maxlength .- length.(b)

    b = pad_array.(b, Δlengths, 0)
    basis_out = reduce(hcat, b)


    return basis_out
end
limit_basis(b::AbstractVector{<:Number}, maxlength) = b[1:min(length(b), maxlength)]
limit_basis(b::AbstractMatrix, maxlength) = b[1:min(length(b), maxlength), :]

Base.length(c::AbstractComponent) = isa(basis(c), Tuple) ? basis(c)[2] : length(basis(c))



"""
    maxlength(c::Vector{<:AbstractComponent}) = maximum(length.(c))
    maxlength(components::Dict) 
maximum of individual component lengths
"""
maxlength(c::Vector{<:AbstractComponent}) = maximum(length.(c))


maxlength(components::Dict) = maximum([maximum(length.(c)) for c in values(components)])
"""
    simulate_component(rng, c::AbstractComponent, simulation::Simulation)
By default call `simulate_component` with `(::Abstractcomponent,::AbstractDesign)` instead of the whole simulation. This allows users to provide a hook to do something completely different :)
"""
simulate_component(rng, c::AbstractComponent, simulation::Simulation) =
    simulate_component(rng, c, simulation.design)

"""
    simulate_component(rng, c::LinearModelComponent, design::AbstractDesign)
Generate a linear model design matrix, weight it by c.β and multiply the result with the given basis vector.

julia> c = UnfoldSim.LinearModelComponent([0,1,1,0],@formula(0~1+cond),[1,2],Dict())
julia> design = MultiSubjectDesign(;n_subjects=2,n_items=50,items_between=(;:cond=>["A","B"]))
julia> simulate_component(StableRNG(1),c,design)
"""
function simulate_component(rng, c::LinearModelComponent, design::AbstractDesign)
    events = generate_events(deepcopy(rng), design)
    X = generate_designmatrix(c.formula, events, c.contrasts)
    y = X * c.β

    return y' .* basis(c, design)
end


"""
Helper function to generate a designmatrix from formula, events and contrasts.
"""
function generate_designmatrix(formula, events, contrasts)
    # special case, intercept only 
    # https://github.com/JuliaStats/StatsModels.jl/issues/269
    if formula.rhs == ConstantTerm(1)
        X = ones(nrow(events), 1)
    else
        if isempty(contrasts)
            m = StatsModels.ModelFrame(formula, events)
        else
            m = StatsModels.ModelFrame(formula, events; contrasts = contrasts)
        end
        X = StatsModels.modelmatrix(m)
    end
    return X
end

"""
    simulate_component(rng, c::MixedModelComponent, design::AbstractDesign)
Generates a MixedModel and simulates data according to c.β and c.σs.

A trick is used to remove the Normal-Noise from the MixedModel which might lead to rare numerical instabilities. Practically, we upscale the σs by factor 10000, and provide a σ=0.0001. Internally this results in a normalization where the response scale is 10000 times larger than the noise.

Currently, it is not possible to use a different basis for fixed and random effects, but a code-stub exists (it is slow though).

- `return_parameters` (Bool,false) - can be used to return the per-event parameters used to weight the basis function. Sometimes useful to see what is simulated

julia> design = MultiSubjectDesign(;n_subjects=2,n_items=50,items_between=(;:cond=>["A","B"]))
julia> c = UnfoldSim.MixedModelComponent([0.,1,1,0],@formula(0~1+cond+(1|subject)),[1,2],Dict(:subject=>[2],),Dict())
julia> simulate_component(StableRNG(1),c,design)

"""
function simulate_component(
    rng,
    c::MixedModelComponent,
    design::AbstractDesign;
    return_parameters = false,
)
    events = generate_events(deepcopy(rng), design)

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
    #epoch_data_component = zeros(Int(length(c.basis)), length(design))

    # residual variance for lmm
    σ_lmm = 0.0001

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

    @debug size(basis(c, design))
    # in case the parameters are of interest, we will return those, not them weighted by basis
    b = return_parameters ? [1.0] : basis(c, design)
    @debug :b, typeof(b), size(b), :m, size(m.y')
    if isa(b, AbstractMatrix)
        epoch_data_component = ((m.y' .* b))
    else
        epoch_data_component = kron(b, m.y')
    end
    return epoch_data_component
    #=
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
    =#
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
    epoch_data = init_epoch_data(components, simulation.design)
    simulate_responses!(rng, epoch_data, components, simulation)
    return epoch_data
end

function simulate_responses!(
    rng,
    epoch_data::AbstractArray,
    components::Vector,
    simulation::Simulation,
)
    for c in components
        simulate_and_add!(epoch_data, c, simulation, deepcopy(rng)) # TODO: `deepcopy` necessary?
    end
    return epoch_data
end
function init_epoch_data(components, design)
    max_offset = maxoffset(components)
    min_offset = minoffset(components)
    range_offset = (max_offset - min_offset)
    if n_channels(components) > 1
        epoch_data = zeros(
            n_channels(components),
            maxlength(components) + range_offset,
            length(design),
        )
    else
        epoch_data = zeros(maxlength(components) + range_offset, length(design))
    end
    return epoch_data
end

function simulate_responses(rng, event_component_dict::Dict, s::Simulation)
    #@debug rng.state
    epoch_data = init_epoch_data(event_component_dict, s.design)
    #@debug rng.state
    evts = generate_events(s.design)
    #@debug rng.state
    @debug size(epoch_data), size(evts)
    multichannel = n_channels(event_component_dict) > 1
    for key in keys(event_component_dict)
        if key == '_'
            continue
        end
        s_key = Simulation(
            s.design |> x -> SubselectDesign(x, key),
            event_component_dict,
            s.onset,
            s.noisetype,
        )
        ix = evts.event .== key
        if multichannel
            simulate_responses!(
                rng,
                @view(epoch_data[:, :, ix]),
                event_component_dict[key],
                s_key,
            )
        else
            #@debug sum(ix), size(simulate_responses(rng, event_component_dict[key], s_key)), key
            simulate_responses!(
                rng,
                @view(epoch_data[:, ix]),
                event_component_dict[key],
                s_key,
            )
        end
    end
    return epoch_data
end

"""
    simulate_and_add!(epoch_data::AbstractMatrix, c, simulation, rng)
    simulate_and_add!(epoch_data::AbstractArray, c, simulation, rng)
Helper function to call `simulate_component` and add it to a provided Array.
"""
function simulate_and_add!(
    epoch_data::AbstractMatrix,
    component::AbstractComponent,
    simulation,
    rng,
)
    @debug "matrix"

    off = offset(component) - minoffset(simulation.components)


    @views epoch_data[1+off:length(component)+off, :] .+=
        simulate_component(rng, component, simulation)
end
function simulate_and_add!(
    epoch_data::AbstractArray,
    component::AbstractComponent,
    simulation,
    rng,
)
    @debug "3D Array"
    off = offset(component) - minoffset(simulation.components)
    @views epoch_data[:, 1+off:length(component)+off, :] .+=
        simulate_component(rng, component, simulation)
end



"""
    Drift_Component <: AbstractComponent
A component that adds an evidence accumulation signal according to a sequential sampling model selected in the field model_type.

All fields are mandatory. Works best with [`SequenceDesign`](@ref).

# Fields
- `basis`: an object, if accessed, provides a 'basis-function', e.g. `hanning(40)`, this defines the response at a single event. It will be weighted by the model-prediction. Can also be a tuple `(fun::Function,maxlength::Int)` with a function `fun` that either generates a matrix `size = (maxlength,size(design,1))` or a vector of vectors. If a larger matrix is generated, it is automatically cutoff at `maxlength`
- `time_vec`: Vector which defines the length of the signal
- `Δt`: Float that indicates the timesteps used to generate the time_vec
- `model_type`: Model struct which defines the model to use to generate the traces, e.g. `KellyModel`
- `model_parameters`: Dict. Containing the parameters for the simulation model specified in model_type.

# Examples
```julia-repl
# use the KellyModel and its default parameters to simulate traces from 0:1/500:1.0
fs=500
Δt = 1/fs;
tEnd = 1.0
time_vec = 0:Δt:tEnd;
model_parameter = create_kelly_parameters_dict(KellyModel())
Drift_Component(
    simulate_component,
    time_vec,
    Δt,
    KellyModel,
    model_parameter)
```

See also [`create_kelly_parameters_dict`](@ref), [`KellyModel`](@ref).
"""
struct Drift_Component <: AbstractComponent
    basis::Any
    time_vec::Any
    Δt::Any
    model_type::Any
    model_parameters::Dict
end
"""
    simulate_component(rng, c::Drift_Component, design::AbstractDesign)

Generate evidence accumulation traces by using the model and its parameters specified in the component c.

# Returns
- `Matrix{Float64}`: Simulated component for each event in the events data frame. The output dimensions are `length(c.time_vec) x size(events, 1)`.

# Examples
```julia-repl
# use the KellyModel and its default parameters to simulate traces from 0:1/500:1.0
julia> model_parameter = create_kelly_parameters_dict(KellyModel());

julia> c = Drift_Component(simulate_component, 0:1/500:1.0, 1/500, KellyModel, model_parameter);

julia> design_single = SingleSubjectDesign(conditions = Dict(:drift_rate => [0.5, 0.8], :condition => [1]));

julia> design_seq = SequenceDesign(design_single,"SCR_");

julia> simulate_component(StableRNG(1),c,design_seq)
501x6 Matrix{Float64}:
0.0  0.0  0.0  0.0  0.0  0.0
0.0  0.0  0.0  0.0  0.0  0.0
⋮                        ⋮
```
"""
function UnfoldSim.simulate_component(rng, c::Drift_Component, design::AbstractDesign)
    _, traces = trace_sequential_sampling_model(deepcopy(rng), c, design)
    return traces
end
"""
    calculate_response_times_for_ssm(rng, component::Drift_Component, design::AbstractDesign)

Generate response times of the evidence accumulation by using the model and its parameters specified in the component.

Using same rng as in simulate_component(rng, component::Drift_Component, design::AbstractDesign) to ensure that the response times match the generated traces. 

# Arguments
- `rng::StableRNG`: Random seed to ensure the same traces are created as in the use of the [`simulate_component`](@ref) function.
- `component::Drift_Component`: Component to specify the model and its parameters to simulate the evidence accumulation.
- `design::UnfoldSim.SubselectDesign`: Subselection of the Sequence Design to ensure we only generate rt for the drift component events.

# Returns
- `Vector{Float64}`: Simulated response time for each event in the events data frame. The output dimension is `size(events, 1)`.

```julia-repl
# use the KellyModel and its default parameters to simulate traces from 0:1/500:1.0
julia> model_parameter = create_kelly_parameters_dict(KellyModel());

julia> c = Drift_Component(simulate_component, 0:1/500:1.0, 1/500, KellyModel, model_parameter);

julia> design_single = SingleSubjectDesign(conditions = Dict(:drift_rate => [0.5, 0.8], :condition => [1]));

julia> design_seq = SequenceDesign(design_single,"SCR_");

julia> calculate_response_times_for_ssm(StableRNG(1),c,design_seq)
2-element Vector{Float64}:
 260.70134768436486
 360.1329203034039
```
"""
function calculate_response_times_for_ssm(rng, component::Drift_Component, design::UnfoldSim.SubselectDesign)
    rts, _ = trace_sequential_sampling_model(deepcopy(rng), component, design)
    return rts
end
Base.length(c::Drift_Component) = length(c.time_vec)
"""
    get_model_parameter(rng, evt, d::Dict)

Construct a parameter dictionary of parameter names and values.
"""
function get_model_parameter(rng, evt, d::Dict)
    result_parameter = Dict()
    for key in keys(d)
        result_parameter[key] = get_model_parameter(rng, evt, d[key])
    end
    return result_parameter
end
"""
    get_model_parameter(rng, evt, val)

Recursive default call return the model parameter as it is.
"""
get_model_parameter(rng, evt, val) = val
"""
    get_model_parameter(rng, evt, val::String)
Recursive call if the model parameter specified as string from the conditions in the events data frame.
"""
get_model_parameter(rng, evt, val::String) = evt[val]

"""
    MixedModelComponent <: AbstractComponent

A component that adds a hierarchical relation between parameters according to a Linear Mixed Model (LMM) defined via `MixedModels.jl`.

All fields can be named. Works best with [`MultiSubjectDesign`](@ref).

# Fields
- `basis::Any`: an object, if accessed, provides a 'basis function', e.g. `hanning(40)::Vector`, this defines the response at a single event. It will be weighted by the model prediction. Future versions will allow for functions, as of v0.3 this is restricted to array-like objects
- `formula::Any`: Formula-object in the style of MixedModels.jl e.g. `@formula 0 ~ 1 + cond + (1|subject)`. The left-hand side is ignored.
- `β::Vector` Vector of betas (fixed effects), must fit the formula.
- `σs::Dict` Dict of random effect variances, e.g. `Dict(:subject => [0.5, 0.4])` or to specify correlation matrix `Dict(:subject=>[0.5,0.4,I(2,2)],...)`. Technically, this will be passed to the MixedModels.jl `create_re` function, which creates the θ matrices.
- `contrasts::Dict` (optional): Dict in the style of MixedModels.jl. Determines which coding scheme to use for which categorical variables. Default is empty which corresponds to dummy coding. For more information see <https://juliastats.org/StatsModels.jl/stable/contrasts>.

# Examples
```julia-repl
julia> MixedModelComponent(;
           basis = hanning(40),
           formula = @formula(0 ~ 1 + cond + (1 + cond|subject)),
           β = [1., 2.],
           σs= Dict(:subject => [0.5, 0.4]),
           contrasts=Dict(:cond => EffectsCoding())
       )
MixedModelComponent
  basis: Array{Float64}((40,)) [0.0, 0.006474868681043577, 0.02573177902642726, 0.0572719871733951, 0.10027861829824952, 0.1536378232452003, 0.21596762663442215, 0.28565371929847283, 0.3608912680417737, 0.43973165987233853  …  0.43973165987233853, 0.3608912680417737, 0.28565371929847283, 0.21596762663442215, 0.1536378232452003, 0.10027861829824952, 0.0572719871733951, 0.02573177902642726, 0.006474868681043577, 0.0]
  formula: StatsModels.FormulaTerm{StatsModels.ConstantTerm{Int64}, Tuple{StatsModels.ConstantTerm{Int64}, StatsModels.Term, StatsModels.FunctionTerm{typeof(|), Vector{StatsModels.AbstractTerm}}}}
  β: Array{Float64}((2,)) [1.0, 2.0]
  σs: Dict{Symbol, Vector{Float64}}
  contrasts: Dict{Symbol, EffectsCoding}
```

See also [`LinearModelComponent`](@ref), [`MultichannelComponent`](@ref).
"""
@with_kw struct MixedModelComponent <: AbstractComponent
    basis::Any
    formula::Any # e.g. 0~1+cond 
    β::Vector
    σs::Dict # Dict(:subject=>[0.5,0.4]) or to specify correlationmatrix Dict(:subject=>[0.5,0.4,I(2,2)],...)
    contrasts::Dict = Dict()
end

"""
    LinearModelComponent <: AbstractComponent

A multiple regression component for one subject.

All fields can be named. Works best with [`SingleSubjectDesign`](@ref).

# Fields
- `basis::Any`: an object, if accessed, provides a 'basis function', e.g. `hanning(40)::Vector`, this defines the response at a single event. It will be weighted by the model prediction. Future versions will allow for functions, as of v0.3 this is restricted to array-like objects
- `formula::Any`: StatsModels `formula` object, e.g.  `@formula 0 ~ 1 + cond` (left-hand side must be 0).
- `β::Vector` Vector of betas/coefficients, must fit the formula.
- `contrasts::Dict` (optional): Determines which coding scheme to use for which categorical variables. Default is empty which corresponds to dummy coding.
     For more information see <https://juliastats.org/StatsModels.jl/stable/contrasts>.

# Examples
```julia-repl
julia> LinearModelComponent(;
           basis = hanning(40),
           formula = @formula(0 ~ 1 + cond),
           β = [1., 2.],
           contrasts = Dict(:cond => EffectsCoding())
       )
LinearModelComponent
  basis: Array{Float64}((40,)) [0.0, 0.006474868681043577, 0.02573177902642726, 0.0572719871733951, 0.10027861829824952, 0.1536378232452003, 0.21596762663442215, 0.28565371929847283, 0.3608912680417737, 0.43973165987233853  …  0.43973165987233853, 0.3608912680417737, 0.28565371929847283, 0.21596762663442215, 0.1536378232452003, 0.10027861829824952, 0.0572719871733951, 0.02573177902642726, 0.006474868681043577, 0.0]
  formula: StatsModels.FormulaTerm{StatsModels.ConstantTerm{Int64}, Tuple{StatsModels.ConstantTerm{Int64}, StatsModels.Term}}
  β: Array{Float64}((2,)) [1.0, 2.0]
  contrasts: Dict{Symbol, EffectsCoding}
```

See also [`MixedModelComponent`](@ref), [`MultichannelComponent`](@ref).
"""
@with_kw struct LinearModelComponent <: AbstractComponent
    basis::Any
    formula::Any # e.g. 0 ~ 1 + cond - left side must be "0"
    β::Vector
    contrasts::Dict = Dict()
end


"""
    MultichannelComponent <: AbstractComponent

Projects a `AbstractComponent` to multiple "channels" via the `projection` vector.
    
Optionally, `noise` can be added to the source prior to projection.
By default a `MultichannelComponent` can be constructed using one of the following options for `projection`:
- `projection::AbstractVector`: Directly pass a custom projection vector.
- `projection::Pair{<:AbstractHeadmodel,String}`: Generate a projection vector by specifying which headmodel to use and which sources should be active.

# Fields
- `component::AbstractComponent`: The component that should be projected to the sensors.
- `projection::AbstractVector` or `projection::Pair{<:AbstractHeadmodel,String}`: Vector `p` that projects the (source) component `c[t]` (where `t` is time) to the sensors `s`.
    The length of `p` equals the number of sensors `s`. Typically, it is a slice of the leadfield matrix. `out[s,t] = p[s]*c[t]`.
- `noise::AbstractNoise` (optional): Noise added in the source space. Default is `NoNoise`.

# Examples
```julia-repl
# Variant 1: Specify the projection vector manually
julia> c1 = LinearModelComponent(; basis = p100(), formula = @formula(0 ~ 1), β = [1]);

julia> mc1 = UnfoldSim.MultichannelComponent(c, [1, 2, -1, 3, 5, 2.3, 1])
MultichannelComponent
  component: LinearModelComponent
  projection: Array{Float64}((7,)) [1.0, 2.0, -1.0, 3.0, 5.0, 2.3, 1.0]
  noise: NoNoise NoNoise()

# Variant 2: Use a headmodel and specify a source
julia> c2 = LinearModelComponent(; basis = p300(), formula = @formula(0 ~ 1), β = [1]);

julia> hart = headmodel(type = "hartmut");
Please cite: HArtMuT: Harmening Nils, Klug Marius, Gramann Klaus and Miklody Daniel - 10.1088/1741-2552/aca8ce

julia> mc2 = UnfoldSim.MultichannelComponent(c2, hart => "Right Occipital Pole")
MultichannelComponent
  component: LinearModelComponent
  projection: Array{Float64}((227,)) [-0.03461859471337842, -0.04321094803502425, 0.0037088347968313525, -0.014722528968861278, -0.0234889834534478, 0.02731807504242923, 0.038863688452528036, 0.1190531258070562, -0.09956890221613562, -0.0867729334438599  …  0.37435404409695094, -0.020863789022627935, 0.25627478723535513, -0.05777985212119245, 0.37104376432271147, -0.19446620423767172, 0.2590764703721097, -0.12923837607416555, 0.1732886690359311, 0.4703016561960567]
  noise: NoNoise NoNoise()
```

See also [`LinearModelComponent`](@ref), [`MixedModelComponent`](@ref).
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

Base.length(c::AbstractComponent) = length(c.basis)
Base.length(c::MultichannelComponent) = length(c.component)

"""
    n_channels(c::AbstractComponent)

Return the number of channels for the given component `c`. By default = 1.
"""
n_channels(c::AbstractComponent) = 1

"""
    n_channels(c::MultichannelComponent)

For `MultichannelComponent` return the length of the projection vector.

"""
n_channels(c::MultichannelComponent) = length(c.projection)

"""
    n_channels(c::Vector{<:AbstractComponent})

For a vector of `MultichannelComponent`s, return the number of channels for the first component but assert all are of equal length.
"""
function n_channels(c::Vector{<:AbstractComponent})
    all_channels = n_channels.(c)
    @assert length(unique(all_channels)) == 1 "Error - projections of different channels cannot be different from each other."
    return all_channels[1]
end

"""
    maxlength(c::Vector{<:AbstractComponent}) = maximum(length.(c))

Return the maximum of the individual components' lengths.
"""
maxlength(c::Vector{<:AbstractComponent}) = maximum(length.(c))

"""
    simulate_component(rng, c::AbstractComponent, simulation::Simulation)

By default call `simulate_component` with `(rng, c::Abstractcomponent, design::AbstractDesign)` instead of the whole simulation. This function exist solely to provide a "hook" if for a custom component something else than the design is necessary, e.g. a dependency on the onsets, noise or similar.
"""
simulate_component(rng, c::AbstractComponent, simulation::Simulation) =
    simulate_component(rng, c, simulation.design)

"""
    simulate_component(rng, c::LinearModelComponent, design::AbstractDesign)

Generate a linear model design matrix, weight it by the coefficients `c.β` and multiply the result with the given basis vector.

# Returns
- `Matrix{Float64}`: Simulated component for each event in the events data frame. The output dimensions are `length(c.basis) x length(design)`.

# Examples
```julia-repl
julia> design = SingleSubjectDesign(; conditions = Dict(:cond => ["natural", "artificial"]));

julia> c = UnfoldSim.LinearModelComponent([0, 1, 1, 0], @formula(0 ~ 1 + cond), [1, 2], Dict());

julia> using StableRNGs

julia> simulate_component(StableRNG(1), c, design)
4×2 Matrix{Float64}:
 0.0  0.0
 3.0  1.0
 3.0  1.0
 0.0  0.0
```
"""
function simulate_component(rng, c::LinearModelComponent, design::AbstractDesign)
    events = generate_events(deepcopy(rng), design)

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
    simulate_component(rng, c::MixedModelComponent, design::AbstractDesign, return_parameters = false)

Generate a MixedModel and simulate data according to the given parameters `c.β` and `c.σs`.


# Keyword arguments
- `return_parameters::Bool = false`: Can be used to return the per-event parameters used to weight the basis function. Sometimes useful to inspect what is simulated.

# Returns
- `Matrix{Float64}`: Simulated component for each event in the events data frame. The output dimensions are `length(c.basis) x length(design)`.

# Notes
1) MixedModels/Sim does not allow simulation of data without white noise of the residuals. Because we want our own noise, we use the following trick to remove the MixedModels-Noise:
Practically, we upscale the specified `σs` by factor 10_000, and request a white-noise-level of `σ = 0.0001`.
Internally in MixedModels/Sim, `σs` are relative to `σ`, and thus are normalized correctly, while keeping the noise 10_000 times smaller than the random effects

We cannot exclude that this trick runs into strange numerical issues if the random effect `σs` are very large compared to the fixed effects.

2) Currently, it is not possible to use a different basis for fixed and random effects. If this is needed, some code-scaffold is available but commented out at the moment and requires a bit of implementation work.



# Examples
```julia-repl
julia> design = MultiSubjectDesign(; n_subjects = 2, n_items = 6, items_between = Dict(:cond => ["A", "B"]));

julia> c = UnfoldSim.MixedModelComponent([0, 1, 1, 0], @formula(0 ~ 1 + cond + (1|subject)), [1, 2], Dict(:subject => [2],), Dict());

julia> using StableRNGs

julia> simulate_component(StableRNG(1), c, design)
4×12 Matrix{Float64}:
 -0.0      -0.0       -0.0      -0.0       -0.0     -0.0       0.0       0.0      0.0       0.0      0.0       0.0
 -2.70645  -0.706388  -2.70632  -0.706482  -2.7066  -0.706424  0.325722  2.32569  0.325627  2.32564  0.325468  2.32565
 -2.70645  -0.706388  -2.70632  -0.706482  -2.7066  -0.706424  0.325722  2.32569  0.325627  2.32564  0.325468  2.32565
 -0.0      -0.0       -0.0      -0.0       -0.0     -0.0       0.0       0.0      0.0       0.0      0.0       0.0

julia> simulate_component(StableRNG(1), c, design, return_parameters = true)
1×12 Matrix{Float64}:
 -2.70645  -0.706388  -2.70632  -0.706482  -2.7066  -0.706424  0.325722  2.32569  0.325627  2.32564  0.325468  2.32565
```
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

    # in case the parameters are of interest, we will return those, not them weighted by basis
    epoch_data_component = kron(return_parameters ? [1.0] : c.basis, m.y')
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
    simulate_component(rng, c::MultichannelComponent, design::AbstractDesign)

Return the projection of a `MultichannelComponent c` from "source" to "sensor" space.

# Returns
- `Array{Float64,3}`: Projected simulated component for each event in the events data frame. The output dimensions are `length(c.projection) x length(c.basis) x length(design)`.

# Examples
```julia-repl
julia> design = SingleSubjectDesign(; conditions = Dict(:cond => ["natural", "artificial"]));

julia> c = LinearModelComponent(; basis = p100(), formula = @formula(0 ~ 1 + cond), β = [1, 0.5]);
julia> mc = UnfoldSim.MultichannelComponent(c, [1, 2, -1, 3, 5, 2.3, 1], PinkNoise());

julia> using StableRNGs

julia> simulate_component(StableRNG(1), mc, design)
7×15×2 Array{Float64, 3}:
[:, :, 1] =
  0.859626   1.16574   0.959523   0.757522   1.17639   1.65156   1.3742    1.76706   2.76971   2.0306    1.17429   1.00922   1.09519   0.754659   2.25662
  1.71925    2.33147   1.91905    1.51504    2.35278   3.30312   2.7484    3.53412   5.53942   4.0612    2.34858   2.01845   2.19039   1.50932    4.51324
 -0.859626  -1.16574  -0.959523  -0.757522  -1.17639  -1.65156  -1.3742   -1.76706  -2.76971  -2.0306   -1.17429  -1.00922  -1.09519  -0.754659  -2.25662
  2.57888    3.49721   2.87857    2.27257    3.52917   4.95469   4.1226    5.30118   8.30913   6.09179   3.52287   3.02767   3.28558   2.26398    6.76985
  4.29813    5.82868   4.79761    3.78761    5.88194   8.25781   6.871     8.8353   13.8485   10.153     5.87145   5.04612   5.47597   3.7733    11.2831
  1.97714    2.68119   2.2069     1.7423     2.70569   3.79859   3.16066   4.06424   6.37033   4.67037   2.70087   2.32121   2.51894   1.73572    5.19022
  0.859626   1.16574   0.959523   0.757522   1.17639   1.65156   1.3742    1.76706   2.76971   2.0306    1.17429   1.00922   1.09519   0.754659   2.25662

[:, :, 2] =
  0.859626   1.16574   0.959523   0.757522   1.17639   1.65156   1.31571   1.56047   2.39471   1.54567   0.689367   0.634223   0.888605   0.69617   2.25662
  1.71925    2.33147   1.91905    1.51504    2.35278   3.30312   2.63142   3.12094   4.78942   3.09135   1.37873    1.26845    1.77721    1.39234   4.51324
 -0.859626  -1.16574  -0.959523  -0.757522  -1.17639  -1.65156  -1.31571  -1.56047  -2.39471  -1.54567  -0.689367  -0.634223  -0.888605  -0.69617  -2.25662
  2.57888    3.49721   2.87857    2.27257    3.52917   4.95469   3.94713   4.68142   7.18413   4.63702   2.0681     1.90267    2.66582    2.08851   6.76985
  4.29813    5.82868   4.79761    3.78761    5.88194   8.25781   6.57855   7.80236  11.9735    7.72837   3.44684    3.17112    4.44303    3.48085  11.2831
  1.97714    2.68119   2.2069     1.7423     2.70569   3.79859   3.02613   3.58909   5.50783   3.55505   1.58554    1.45871    2.04379    1.60119   5.19022
  0.859626   1.16574   0.959523   0.757522   1.17639   1.65156   1.31571   1.56047   2.39471   1.54567   0.689367   0.634223   0.888605   0.69617   2.25662
```
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
    weight_σs(σs::Dict, b_σs::Float64, σ_lmm::Float64)

Weight a `σs` Dict for MixedModels.jl by `b_σs`, a scaling factor typically from a `basis`.

Finally scales it again by `σ_lmm`, as a trick to simulate noise-free LMMs (see `MixedModelsComponent`)

In the future, we anticipate a function
    `function weight_σs(σs::Dict,b_σs::Dict,σ_lmm::Float64)`
where each `σs` entry can be weighted individually by a matching `b_σs`, but it is not implemented.

# Arguments
- `σs::Dict` = a Dict of random effects as output of MixedModels.create_re
- `b_σs::Float64` = a scaling factor, typically one entry of a basis function from a component
- `σ_lmm::Float64` = a scaling factor to simulate near-zero noise LMMs
# Returns
    `NamedTuple` of the weighted random effects

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

Simulate multiple component responses and accumulate them on a per-event basis.

# Returns
- `epoch_data`: `Matrix` (or `Array` in the multi-channel case) of combined simulated components.
    The output dimensions are `maxlength(components) x length(simulation.design)` for single-channel components and
    `n_channels(components) x maxlength(components) x length(simulation.design)` for multi-channel components.

# Examples
```julia-repl
julia> design = SingleSubjectDesign(; conditions = Dict(:cond => ["natural", "artificial"]));

julia> c1 = LinearModelComponent(; basis = p100(), formula = @formula(0 ~ 1 + cond), β = [1, 0.5]);
julia> c2 = LinearModelComponent(; basis = p300(), formula = @formula(0 ~ 1), β = [2]);

julia> simulation = Simulation(design, [c1, c2], UniformOnset(; width = 0, offset = 30), PinkNoise());

julia> using StableRNGs

julia> simulate_responses(StableRNG(1), [c1, c2], simulation)
45×2 Matrix{Float64}:
 0.0        0.0
 0.0        0.0
 0.0        0.0
 0.0        0.0
 0.0        0.0
 ⋮          
 0.352614   0.352614
 0.203907   0.203907
 0.0924246  0.0924246
 0.0233794  0.0233794
 0.0        0.0
```
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
        simulate_and_add!(epoch_data, c, simulation, deepcopy(rng)) # TODO: `deepcopy` necessary?
    end
    return epoch_data
end


"""
    simulate_and_add!(epoch_data::AbstractMatrix, c, simulation, rng)
    simulate_and_add!(epoch_data::AbstractArray, c, simulation, rng)

Helper function to call `simulate_component` and add it to a provided Array `epoch_data`.
"""
function simulate_and_add!(epoch_data::AbstractMatrix, c, simulation, rng)
    @debug "matrix"
    @views epoch_data[1:length(c), :] .+= simulate_component(rng, c, simulation)
end
function simulate_and_add!(epoch_data::AbstractArray, c, simulation, rng)
    @debug "3D Array"
    @views epoch_data[:, 1:length(c), :] .+= simulate_component(rng, c, simulation)
end

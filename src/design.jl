""" 
    MultiSubjectDesign <: AbstractDesign

A type for specifying the experimental design for multiple subjects (based on the given random-effects structure).

### Fields
- `n_subjects`::Int -> number of subjects
- `n_items`::Int -> number of items (sometimes ≈trials)
- `subjects_between` = Dict{Symbol,Vector} -> effects between subjects, e.g. young vs old 
- `items_between` = Dict{Symbol,Vector} -> effects between items, e.g. natural vs artificial images, (but shown to all subjects if not specified also in `subjects_between`)
- `both_within` = Dict{Symbol,Vector}	-> effects completly crossed
- `event_order_function` = `(rng,x)->x`; # can be used to sort, or e.g. `(rng,x)->shuffle(rng,x)` (or shorter just `event_order_function=shuffle`)

Tip: Check the resulting dataframe using `generate_events(design)`


### Example
```julia
# declaring same condition both sub-between and item-between results in a full between subject/item design
design = MultiSubjectDesign(;
		n_items = 10,
		n_subjects = 30,
		subjects_between = Dict(:cond => ["levelA", "levelB"]),
		items_between = Dict(:cond => ["levelA", "levelB"]),
		);
```
See also [`SingleSubjectDesign`](@ref), [`RepeatDesign`](@ref)
"""
@with_kw struct MultiSubjectDesign <: AbstractDesign
    n_subjects::Int
    n_items::Int
    subjects_between::Dict{Symbol,Vector} = Dict()
    items_between::Dict{Symbol,Vector} = Dict()
    both_within::Dict{Symbol,Vector} = Dict()
    event_order_function = (rng, x) -> x
end


"""
    SingleSubjectDesign <: AbstractDesign

A type for specifying the experimental for a single subject (based on the given conditions).

### Fields
- `conditions = Dict{Symbol,Vector}`` of conditions, e.g. `Dict(:A=>["a_small","a_big"],:B=>["b_tiny","b_large"])`
- `event_order_function` = (rng::AbstractRNG,x::DataFrame)->x; # can be used to sort by specifying `sort`, or shuffling by providing `shuffle`, or custom functions following the interface `(rng,x)->my_shuffle(rng,x)`

Number of trials / rows in `generate_events(design)` depend on the full factorial of your `conditions`.

To increase the number of repetitions simply use `RepeatDesign(SingleSubjectDesign(...),5)`

If conditions are omitted (or set to `nothing`), a single trial is simulated with a column `:dummy` and content `:dummy` - this is for convenience.

Tip: Check the resulting dataframe using `generate_events(design)`

### Example
```julia
design = SingleSubjectDesign(;
    conditions = Dict(
        :stimulus_type => ["natural", "artificial"],
        :contrast_level => range(0, 1, length = 5),
);
```
See also [`MultiSubjectDesign`](@ref), [`RepeatDesign`](@ref)
"""
@with_kw struct SingleSubjectDesign <: AbstractDesign
    conditions::Dict{Symbol,Vector} = Dict()
    event_order_function = (rng, x) -> x
end


""" Returns dimension of experiment design"""
size(design::MultiSubjectDesign) = (design.n_items, design.n_subjects)
size(design::SingleSubjectDesign) = (*(length.(values(design.conditions))...),)

"""
Generates full-factorial DataFrame of design.conditions


Afterwards applies design.event_order_function.

If conditions is `nothing`, a single trial is simulated with a column `:dummy` and content `:dummy` - this is for convenience.


julia> d = SingleSubjectDesign(;conditions= Dict(:A=>nlevels(5),:B=>nlevels(2)))
julia> generate_events(d)
"""
function generate_events(rng::AbstractRNG, design::SingleSubjectDesign)
    if isempty(design.conditions)
        events = DataFrame(:dummy => [:dummy])
    else
        # we get a Dict(:A=>["1","2"],:B=>["3","4"]), but needed a list
        # of named tuples for MixedModelsSim.factorproduct function.
        events =
            factorproduct(((; k => v) for (k, v) in pairs(design.conditions))...) |>
            DataFrame
    end
    # by default does nothing
    return apply_event_order_function(design.event_order_function, rng, events)

end

"""
    apply_event_order_function(fun::Function, rng::AbstractRNG, events::DataFrame)

apply fun(rng,events), raise an error if function is wrongly defined. Convenience function to not repeat the error handling at multiple places.
"""
function apply_event_order_function(fun, rng, events)
    try
        return fun(rng, events)

    catch e
        error(
            "Problem in `event_order_function` - Make sure the function allows for two inputs: `(rng::AbstractRNG,x::DataFrame)`",
        )
    end
end
"""
    generate_events([rng::AbstractRNG],design::MultiSubjectDesign)
Generate full factorial Dataframe according to MixedModelsSim.jl 's `simdat_crossed` function.
Note: n_items = you can think of it as `trials` or better, as `stimuli`.

Note: No condition can be named `dv` which is used internally in MixedModelsSim / MixedModels as a dummy left-side

Afterwards applies `design.event_order_function`. Could be used to duplicate trials, sort, subselect etc.

Finally it sorts by `:subject`

julia> d = MultiSubjectDesign(;n_subjects = 10,n_items=20,both_within= Dict(:A=>nlevels(5),:B=>nlevels(2)))
julia> generate_events(d)
"""

function generate_events(rng::AbstractRNG, design::MultiSubjectDesign)

    # check that :dv is not in any condition
    allconditions = [design.subjects_between, design.items_between, design.both_within]

    @assert all(isempty.(allconditions)) ||
            :dv ∉ keys(merge(allconditions[.!isempty.(allconditions)]...)) "due to technical limitations in MixedModelsSim.jl, `:dv` cannot be used as a factorname"

    data = DataFrame(
        MixedModelsSim.simdat_crossed(
            design.n_subjects,
            design.n_items,
            subj_btwn = isempty(design.subjects_between) ? nothing :
                        design.subjects_between,
            item_btwn = isempty(design.items_between) ? nothing : design.items_between,
            both_win = isempty(design.both_within) ? nothing : design.both_within,
        ),
    )
    rename!(data, :subj => :subject)
    select!(data, Not(:dv)) # remove the default column from MixedModelsSim.jl - we don't need it in UnfoldSim.jl
    # by default does nothing
    data = apply_event_order_function(design.event_order_function, rng, data)

    # sort by subject
    data = sort!(data, (order(:subject)))

    return data

end

generate_events(design::AbstractDesign) = generate_events(MersenneTwister(1), design)

# length is the same of all dimensions
length(design::AbstractDesign) = *(size(design)...)


# ----

"""
    RepeatDesign{T} <: AbstractDesign
Repeat a design DataFrame multiple times to mimick repeatedly recorded trials.

```julia
designOnce = MultiSubjectDesign(;
		n_items=2,
		n_subjects = 2,
		subjects_between =Dict(:cond=>["levelA","levelB"]),
		items_between =Dict(:cond=>["levelA","levelB"]),
		);

design = RepeatDesign(designOnce,4);
```
See also [`SingleSubjectDesign`](@ref), [`MultiSubjectDesign`](@ref)
"""
@with_kw struct RepeatDesign{T} <: AbstractDesign
    design::T
    repeat::Int = 1
end

"""
    UnfoldSim.generate_events(design::RepeatDesign{T})

In a repeated design, iteratively calls the underlying {T} Design and concatenates. In case of MultiSubjectDesign, sorts by subject.
"""
function UnfoldSim.generate_events(rng::AbstractRNG, design::RepeatDesign)
    df = map(x -> generate_events(rng, design.design), 1:design.repeat) |> x -> vcat(x...)
    if isa(design.design, MultiSubjectDesign)
        sort!(df, [:subject])
    end
    return df

end
Base.size(design::RepeatDesign{MultiSubjectDesign}) =
    size(design.design) .* (design.repeat, 1)
Base.size(design::RepeatDesign{SingleSubjectDesign}) = size(design.design) .* design.repeat

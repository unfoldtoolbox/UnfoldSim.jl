""" 
    MultiSubjectDesign

- `n_subjects`::Int -> number of subjects
- `n_items`::Int -> number of items (sometimes ≈trials)
- `subjects_between` = Dict{Symbol,Vector} -> effects between subjects, e.g. young vs old 
- `items_between` = Dict{Symbol,Vector} -> effects between items, e.g. natural vs artificial images, (but shown to all subjects if not specified also in `subjects_between`)
- `both_within` = Dict{Symbol,Vector}	-> effects completly crossed
- `event_order_function` = `x->x`; # can be used to sort, or e.g. `x->shuffle(MersenneTwister(42),x)` - be sure to fix/update the rng accordingly!!

tipp: check the resulting dataframe using `generate_events(design)`



```julia
# declaring same condition both sub-between and item-between results in a full between subject/item design
design = MultiSubjectDesign(;
		n_items = 10,
		n_subjects = 30,
		subjects_between = Dict(:cond => ["levelA", "levelB"]),
		items_between = Dict(:cond => ["levelA", "levelB"]),
		);
```
"""
@with_kw struct MultiSubjectDesign <: AbstractDesign
    n_subjects::Int
    n_items::Int
    subjects_between::Dict{Symbol,Vector} = Dict()
    items_between::Dict{Symbol,Vector} = Dict()
    both_within::Dict{Symbol,Vector} = Dict()
    event_order_function = x -> x # can be used to sort, or x->shuffle(rng,x)
end


"""
- conditions = Dict{Symbol,Vector} of conditions, e.g. `Dict(:A=>["a_small","a_big"],:B=>["b_tiny","b_large"])`
- `event_order_function` = x->x; # can be used to sort, or x->shuffle(MersenneTwister(42),x) - be sure to fix/update the rng accordingly!!

Number of trials / rows in `generate_events(design)` depend on the full factorial of your `conditions`.

To increase the number of repetitions simply use `RepeatDesign(SingleSubjectDesign(...),5)`

If conditions are omitted (or set to `nothing`), a single trial is simulated with a column `:dummy` and content `:dummy` - this is for convenience.

tipp: check the resulting dataframe using `generate_events(design)`
"""
@with_kw struct SingleSubjectDesign <: AbstractDesign
    conditions::Dict{Symbol,Vector} = Dict()
    event_order_function = x -> x
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
function generate_events(design::SingleSubjectDesign)
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
    return design.event_order_function(events)
end

"""
    generate_events(design::MultiSubjectDesign)
Generate full factorial Dataframe according to MixedModelsSim.jl 's `simdat_crossed` function.
Note: n_items = you can think of it as `trials` or better, as `stimuli`.

Note: No condition can be named `dv` which is used internally in MixedModelsSim / MixedModels as a dummy left-side

Afterwards applies `design.event_order_function``.  Could be used to duplicate trials, sort, subselect etc.

Finally it sorts by `:subject`

julia> d = MultiSubjectDesign(;n_subjects = 10,n_items=20,both_within= Dict(:A=>nlevels(5),:B=>nlevels(2)))
julia> generate_events(d)
"""
function generate_events(design::MultiSubjectDesign)

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
    data = design.event_order_function(data)

    # sort by subject
    data = sort!(data, (order(:subject)))

    return data

end


# length is the same of all dimensions
length(design::AbstractDesign) = *(size(design)...)



# ----

"""
    RepeatDesign{T}
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
"""
@with_kw struct RepeatDesign{T} <: AbstractDesign
    design::T
    repeat::Int = 1
end


function check_sequence(s::String)
    blankfind = findall('_', s)
    @assert length(blankfind) <= 1 && (length(blankfind) == 0 || length(s) == blankfind[1]) "the blank-indicator '_' has to be the last sequence element"
    return s
end
@with_kw struct SequenceDesign{T} <: AbstractDesign
    design::T
    sequence::String = ""
    sequencelength::Int = 0
    rng = nothing

    SequenceDesign{T}(d, s, sl, r) where {T<:AbstractDesign} =
        new(d, check_sequence(s), sl, r)
end

SequenceDesign(design, sequence) = SequenceDesign(design = design, sequence = sequence)

generate_events(design::SequenceDesign{MultiSubjectDesign}) = error("not yet implemented")


generate_events(rng, design::AbstractDesign) = generate_events(design)
generate_events(design::SequenceDesign) = generate_events(deepcopy(design.rng), design)

function generate_events(rng, design::SequenceDesign)
    df = generate_events(design.design)
    nrows_df = size(df, 1)

    rng = if isnothing(rng)
        @warn "Could not (yet) find an rng for `SequenceDesign` - ignore this message if you called `generate_events` yourself, be worried if you called `simulate` and still see this message. Surpress this message by defining the `rng` when creating the `SequenceDesign`"
        MersenneTwister(1)
    else
        rng
    end
    #   @debug design.sequence
    currentsequence = sequencestring(rng, design.sequence)
    #    @debug currentsequence
    currentsequence = replace(currentsequence, "_" => "")
    df = repeat(df, inner = length(currentsequence))

    df.event .= repeat(collect(currentsequence), nrows_df)

    return df

end



get_rng(design::AbstractDesign) = nothing
get_rng(design::SequenceDesign) = design.rng

"""
    generate_events(rng,design::RepeatDesign{T})

In a repeated design, iteratively calls the underlying {T} Design and concatenates. In case of MultiSubjectDesign, sorts by subject.
"""
function generate_events(design::RepeatDesign)
    design = deepcopy(design)
    df =
        map(x -> generate_events(get_rng(design.design), design.design), 1:design.repeat) |>
        x -> vcat(x...)
    if isa(design.design, MultiSubjectDesign)
        sort!(df, [:subject])
    end
    return df

end


"""
Internal helper design to subset a sequence design in its individual components
"""
struct SubselectDesign{T} <: AbstractDesign
    design::T
    key::Char
end

function generate_events(design::SubselectDesign)
    return subset(generate_events(design.design), :event => x -> x .== design.key)
end


Base.size(design::RepeatDesign{MultiSubjectDesign}) =
    size(design.design) .* (design.repeat, 1)
Base.size(design::RepeatDesign{SingleSubjectDesign}) = size(design.design) .* design.repeat
#Base.size(design::SequenceDesign) =
#size(design.design) .* length(replace(design.sequence, "_" => "",r"\{.*\}"=>""))

#Base.size(design::) = size(design.design) .* design.repeat

# No way to find out what size it is without actually generating first...
Base.size(
    design::Union{<:SequenceDesign,<:SubselectDesign,<:RepeatDesign{<:SequenceDesign}},
) = size(generate_events(design), 1)
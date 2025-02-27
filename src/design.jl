# Design types

"""
    SingleSubjectDesign <: AbstractDesign

A type for specifying the experimental design for a single subject (based on the given conditions).

Tip: Check the resulting dataframe using the [`generate_events`](@ref) function. \\
The number of trials/rows in the output of `generate_events([rng, ]design)` depends on the full factorial of your `conditions`. \\
To increase the number of repetitions, e.g. by 5, simply use `RepeatDesign(SingleSubjectDesign(...), 5)`. \\
If conditions are omitted (or set to `nothing`), a single trial is simulated with a column `:dummy` and content `:dummy` - this is for convenience.

# Fields
- `conditions::Dict{Symbol,Vector}`: Experimental conditions, e.g. `Dict(:A => ["a_small","a_big"], :B => ["b_tiny","b_large"])`.
- `event_order_function = (rng, x) -> x`: Can be used to sort by specifying `sort`, or shuffling by providing `shuffle`, or custom functions following the interface `(rng, x) -> my_shuffle(rng,x)`.
    The default is the identify function, i.e. not changing the order of the events.

# Examples
```julia-repl
julia> using StableRNGs # For using the `generate_events` function in a reproducible way

julia> design =
           SingleSubjectDesign(;
               conditions = Dict(
                   :stimulus_type => ["natural", "artificial"],
                   :contrast_level => range(0, 1, length = 5),
               ),
           )
SingleSubjectDesign
  conditions: Dict{Symbol, Vector}
  event_order_function: #10 (function of type UnfoldSim.var"#10#14")

julia> generate_events(StableRNG(1), design)
10×2 DataFrame
 Row │ contrast_level  stimulus_type 
     │ Float64         String        
─────┼───────────────────────────────
   1 │           0.0   natural
   2 │           0.25  natural
   3 │           0.5   natural
  ⋮  │       ⋮               ⋮
   8 │           0.5   artificial
   9 │           0.75  artificial
  10 │           1.0   artificial
                       4 rows omitted
```

See also [`MultiSubjectDesign`](@ref), [`RepeatDesign`](@ref)
"""
@with_kw struct SingleSubjectDesign <: AbstractDesign
    conditions::Dict{Symbol,Vector} = Dict()
    event_order_function = (rng, x) -> x
end

""" 
    MultiSubjectDesign <: AbstractDesign

A type for specifying the experimental design for multiple subjects (based on the given random-effects structure).

Tip: Check the resulting dataframe using the [`generate_events`](@ref) function. \\
Please note that the number of items `n_items` has to be a multiple of the number of between-item levels. The sample applies for `n_subjects` and the number of between-subject levels.

# Fields
- `n_subjects::Int`: Number of subjects.
- `n_items::Int`: Number of items/stimuli (sometimes ≈ trials).
- `subjects_between::Dict{Symbol,Vector}`: Effects between subjects, e.g. young vs old.
- `items_between::Dict{Symbol,Vector}`: Effects between items, e.g. natural vs artificial images, (but shown to all subjects if not specified also in `subjects_between`).
- `both_within::Dict{Symbol,Vector}`: Effects completely crossed i.e. conditions/covariates that are both within-subject and within-item.
- `event_order_function = (rng, x) -> x`: Can be used to sort, or shuffle the events e.g. `(rng, x) -> shuffle(rng, x)` (or shorter just `event_order_function = shuffle`).
    The default is the identify function, i.e. not changing the order of the events.

# Examples
```julia-repl
# Declaring the same condition both between-subject and between-item results in a full between-subject/item design.
julia> design = MultiSubjectDesign(;
                       n_items = 10,
                       n_subjects = 30,
                       subjects_between = Dict(:cond => ["levelA", "levelB"]),
                       items_between = Dict(:cond => ["levelA", "levelB"]),
                       )
MultiSubjectDesign
  n_subjects: Int64 30
  n_items: Int64 10
  subjects_between: Dict{Symbol, Vector}
  items_between: Dict{Symbol, Vector}
  both_within: Dict{Symbol, Vector}
  event_order_function: #3 (function of type UnfoldSim.var"#3#7")

julia> generate_events(StableRNG(1), design)
150×3 DataFrame
 Row │ subject  cond    item   
     │ String   String  String 
─────┼─────────────────────────
   1 │ S01      levelA  I01
   2 │ S01      levelA  I03
   3 │ S01      levelA  I05
  ⋮  │    ⋮       ⋮       ⋮
 148 │ S30      levelB  I06
 149 │ S30      levelB  I08
 150 │ S30      levelB  I10
               144 rows omitted
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
    RepeatDesign{T} <: AbstractDesign

Repeat a design (and the corresponding events DataFrame) multiple times to mimick repeatedly recorded trials.

Tip: Check the resulting dataframe using the [`generate_events`](@ref) function.
Please note that when using an `event_order_function`(e.g. `shuffle`) in a `RepeatDesign`, the corresponding RNG is shared across repetitions and not deep-copied for each repetition.
As a result, the order of events will differ for each repetition.

# Fields
- `design::T`: The experimental design that should be repeated.
- `repeat::Int = 1`: The number of repetitions.

# Examples
```julia-repl
julia> using StableRNGs # For using the `generate_events` function in a reproducible way

julia> design_once =
           SingleSubjectDesign(;
               conditions = Dict(
                   :stimulus_type => ["natural", "artificial"],
                   :contrast_level => range(0, 1, length = 2),
               ),
               event_order_function = shuffle,
           );

julia> generate_events(StableRNG(1), design_once)
4×2 DataFrame
 Row │ contrast_level  stimulus_type 
     │ Float64         String        
─────┼───────────────────────────────
   1 │            1.0  natural
   2 │            1.0  artificial
   3 │            0.0  natural
   4 │            0.0  artificial

julia> design = RepeatDesign(design_once, 2)
RepeatDesign{SingleSubjectDesign}
  design: SingleSubjectDesign
  repeat: Int64 2

julia> generate_events(StableRNG(1), design)
8×2 DataFrame
 Row │ contrast_level  stimulus_type 
     │ Float64         String        
─────┼───────────────────────────────
   1 │            1.0  natural
   2 │            1.0  artificial
   3 │            0.0  natural
   4 │            0.0  artificial
   5 │            1.0  artificial
   6 │            0.0  natural
   7 │            1.0  natural
   8 │            0.0  artificial
```

See also [`SingleSubjectDesign`](@ref), [`MultiSubjectDesign`](@ref)
"""
@with_kw struct RepeatDesign{T} <: AbstractDesign
    design::T
    repeat::Int = 1
end

#----
# Design helper functions

"Return the dimensions of the experiment design."
size(design::MultiSubjectDesign) = (design.n_items, design.n_subjects)
size(design::SingleSubjectDesign) = (*(length.(values(design.conditions))...),)

Base.size(design::RepeatDesign{MultiSubjectDesign}) =
    size(design.design) .* (design.repeat, 1)
Base.size(design::RepeatDesign{SingleSubjectDesign}) = size(design.design) .* design.repeat

"Length is the product of all dimensions and equals the number of events in the corresponding events dataframe."
length(design::AbstractDesign) = *(size(design)...)

"""
    apply_event_order_function(fun, rng, events)

Apply `fun(rng, events)`, raise an error if function `fun` is wrongly defined. Convenience function to not repeat the error handling at multiple places.
"""
function apply_event_order_function(fun, rng, events)
    try
        return fun(rng, events)

    catch e
        error(
            "Problem in `event_order_function` - Make sure the function allows for two inputs: `(rng::AbstractRNG, x::DataFrame)`",
        )
    end
end

#----
# `generate_events` functions for all design types

"""
    generate_events(design::AbstractDesign)
    generate_events(rng::AbstractRNG, design::AbstractDesign)

Generate a full-factorial events DataFrame based on the experimental conditions and covariates defined in the design.

# Arguments
- `design::AbstractDesign`: Experimental design for which the events DataFrame should be created.
- `rng::AbstractRNG` (optional): Random number generator (RNG) to make the process reproducible. If none is given, `MersenneTwister(1)` will be used.

# Returns
- `DataFrame`: Each row corresponds to one combination of condition/covariate levels which is often equivalent to one stimulus or trial.
"""
function generate_events end

# If no rng is given, create one.
generate_events(design::AbstractDesign) = generate_events(MersenneTwister(1), design)

"""
    generate_events(design::SingleSubjectDesign)
    generate_events(rng::AbstractRNG, design::SingleSubjectDesign)

Generate the events DataFrame based on `design.conditions` and afterwards apply `design.event_order_function`.

If `design.conditions` is `nothing`, a single trial is simulated with a column `:dummy` and content `:dummy` - this is for convenience.

# Examples
```julia-repl
julia> using Random # for shuffling
julia> using StableRNGs

julia> design = SingleSubjectDesign(;
           conditions = Dict(:A => ["small", "large"], :B => range(1, 5, length = 3)),
           event_order_function = shuffle,
       );

# Variant 1: Without specifying an RNG, MersenneTwister(1) will be used for the shuffling specified as `event_order_function`.
julia> generate_events(design)
6×2 DataFrame
 Row │ A       B       
     │ String  Float64 
─────┼─────────────────
   1 │ large       5.0
   2 │ large       1.0
   3 │ large       3.0
   4 │ small       1.0
   5 │ small       3.0
   6 │ small       5.0

# Variant 2: Use a custom RNG.
julia> generate_events(StableRNG(1), design)
6×2 DataFrame
 Row │ A       B       
     │ String  Float64 
─────┼─────────────────
   1 │ large       5.0
   2 │ large       3.0
   3 │ small       1.0
   4 │ small       3.0
   5 │ large       1.0
   6 │ small       5.0
```
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
    generate_events(design::MultiSubjectDesign)
    generate_events(rng::AbstractRNG, design::MultiSubjectDesign)

Generate the events Dataframe according to `MixedModelsSim.jl`'s `simdat_crossed` function.

Afterwards apply `design.event_order_function` and finally sort by `:subject`. \\
Note: No condition can be named `dv` which is used internally in MixedModelsSim / MixedModels as a dummy left-side

# Examples
```julia-repl
julia> using Random # for shuffling
julia> using StableRNGs

julia> design = MultiSubjectDesign(;
                       n_items = 4,
                       n_subjects = 5,
                       both_within = Dict(:condition => ["red", "green"]),
                       event_order_function = shuffle,
                       );

julia> generate_events(StableRNG(1), design)
40×3 DataFrame
 Row │ subject  item    condition 
     │ String   String  String    
─────┼────────────────────────────
   1 │ S1       I2      red
   2 │ S1       I2      green
   3 │ S1       I3      green
  ⋮  │    ⋮       ⋮         ⋮
  38 │ S5       I3      red
  39 │ S5       I2      red
  40 │ S5       I4      green
                   34 rows omitted
```
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

# ----


function check_sequence(s::String)
    blankfind = findall('_', s)
    @assert length(blankfind) <= 1 && (length(blankfind) == 0 || length(s) == blankfind[1]) "the blank-indicator '_' has to be the last sequence element"
    return s
end


"""
    SequenceDesign{T} <: AbstractDesign
Enforce a sequence of events for each entry of a provided `AbstractDesign`.
The sequence string can contain any number of `char`, but the `_` character is used to indicate a break between events without any overlap.

It is also possible to define variable length sequences using `{}`. For example, `A{10,20}` would result in a sequence of 10 to 20 `A`'s.

Another variable sequence is defined using `[]`. For example, `S[ABC]` would result in any one sequence `SA`, `SB`, `SC`.

Important: The exact same variable sequence is used for current rows of a design. Only, if you later nest in a `RepeatDesign` then each `RepeatDesign` repetition will gain a new variable sequence. If you need imbalanced designs, please refer to the `ImbalancedDesign` tutorial

```julia
design = SingleSubjectDesign(conditions = Dict(:condition => ["one", "two"]))
design = SequenceDesign(design, "SCR_")
```
Would result in a `generate_events(design)`
```repl
6×2 DataFrame
 Row │ condition  event 
     │ String     Char  
─────┼──────────────────
   1 │ one        S
   2 │ one        C
   3 │ one        R
   4 │ two        S
   5 │ two        C
   6 │ two        R
```

## Example for Sequence -> Repeat vs. Repeat -> Sequence

### Sequence -> Repeat 
```julia
design = SingleSubjectDesign(conditions = Dict(:condition => ["one", "two"]))
design = SequenceDesign(design, "[AB]")
design = RepeatDesign(design,2)
generate_events(design)
```


```repl
4×2 DataFrame
 Row │ condition  event 
     │ String     Char  
─────┼──────────────────
   1 │ one        A
   2 │ two        A
   3 │ one        B
   4 │ two        B
```
Sequence -> Repeat: a sequence design is repeated, then for each repetition a sequence is generated and applied. Events have different values

### Repeat -> Sequence
```julia
design = SingleSubjectDesign(conditions = Dict(:condition => ["one", "two"]))
design = RepeatDesign(design,2)
design = SequenceDesign(design, "[AB]")
generate_events(design)
```

```repl
4×2 DataFrame
 Row │ condition  event 
     │ String     Char  
─────┼──────────────────
   1 │ one        A
   2 │ two        A
   3 │ one        A
   4 │ two        A
```
Repeat -> Sequence: the design is first repeated, then for that design one sequence generated and applied. All events are the same


See also [`SingleSubjectDesign`](@ref), [`MultiSubjectDesign`](@ref), [`RepeatDesign`](@ref)
"""
@with_kw struct SequenceDesign{T} <: AbstractDesign
    design::T
    sequence::String = ""
    sequencelength::Int = 0
    SequenceDesign{T}(d, s, sl) where {T<:AbstractDesign} = new(d, check_sequence(s), sl)
end

SequenceDesign(design, sequence) = SequenceDesign(design = design, sequence = sequence)

generate_events(rng, design::SequenceDesign{MultiSubjectDesign}) =
    error("not yet implemented")


generate_events(rng, design::AbstractDesign) = generate_events(design)

function generate_events(rng, design::SequenceDesign)
    df = generate_events(deepcopy(rng), design.design)
    nrows_df = size(df, 1)

    #   @debug design.sequence
    currentsequence = sequencestring(rng, design.sequence)
    #    @debug currentsequence
    currentsequence = replace(currentsequence, "_" => "")
    df = repeat(df, inner = length(currentsequence))

    df.event .= repeat(collect(currentsequence), nrows_df)

    return df

end


"""
    
    UnfoldSim.generate_events([rng::AbstractRNG, ]design::RepeatDesign{T})

For a `RepeatDesign`, iteratively call `generate_events` for the underlying {T} design and concatenate the results.

In case of `MultiSubjectDesign`, sort by subject. \\
Please note that when using an `event_order_function`(e.g. `shuffle`) in a `RepeatDesign`, the corresponding RNG is shared across repetitions and not deep-copied for each repetition.
As a result, the order of events will differ for each repetition.
"""

function UnfoldSim.generate_events(rng::AbstractRNG, design::RepeatDesign)
    df = map(x -> generate_events(rng, design.design), 1:design.repeat) |> x -> vcat(x...)

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

function generate_events(rng, design::SubselectDesign)
    return subset(generate_events(rng, design.design), :event => x -> x .== design.key)
end



# --- 
# Effects

"""
    EffectsDesign <: AbstractDesign
Design to obtain ground truth simulation.

## Fields
- `design::AbstractDesign`
   The design of your (main) simulation.
- `effects_dict::Dict`
   Effects.jl style dictionary specifying variable effects. See also [Unfold.jl marginalized effects](https://unfoldtoolbox.github.io/Unfold.jl/stable/generated/HowTo/effects/)
"""
struct EffectsDesign <: AbstractDesign
    design::AbstractDesign
    effects_dict::Dict
end
EffectsDesign(design::MultiSubjectDesign, effects_dict::Dict) = error("not yet implemented")
UnfoldSim.size(t::EffectsDesign) = size(generate_events(t), 1)

"""
    expand_grid(design)

calculate all possible combinations of the key/value pairs of the design-dict. Copied from Effects.jl
"""
function expand_grid(design)
    colnames = tuple(Symbol.(keys(design))...)
    rowtab = NamedTuple{colnames}.(Base.Iterators.product(values(design)...))

    return DataFrame(vec(rowtab))
end

typical_value(v::Vector{<:Number}) = [mean(v)]
typical_value(v) = unique(v)

"""
    UnfoldSim.generate_events(rng,design::EffectsDesign)

Generates events to simulate marginalized effects using an Effects.jl reference-grid dictionary. Every covariate that is in the `EffectsDesign` but not in the `effects_dict` will be set to a `typical_value` (i.e. the mean)

```julia
# Example

effects_dict = Dict{Symbol,Union{<:Number,<:String}}(:conditionA=>[0,1])
SingleSubjectDesign(...) |> x-> EffectsDesign(x,effects_dict)
```
"""
function UnfoldSim.generate_events(rng, t::EffectsDesign)
    effects_dict = Dict{Any,Any}(t.effects_dict)
    #effects_dict = t.effects_dict
    current_design = generate_events(deepcopy(rng), t.design)
    to_be_added = setdiff(names(current_design), string.(keys(effects_dict)))
    for tba in to_be_added
        effects_dict[tba] = typical_value(current_design[:, tba])
    end
    return expand_grid(effects_dict)
end


#Base.size(design::SequenceDesign) =
#size(design.design) .* length(replace(design.sequence, "_" => "",r"\{.*\}"=>""))

#Base.size(design::) = size(design.design) .* design.repeat

# ---  
# Size for Sequence design
# No way to find out what size it is without actually generating first...
Base.size(
    design::Union{<:SequenceDesign,<:SubselectDesign,<:RepeatDesign{<:SequenceDesign}},
) = size(generate_events(design), 1)

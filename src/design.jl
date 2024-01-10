""" 
- n_subjects::Int -> number of subjects
- n_items::Int -> number of items (sometimes ≈trials)
- subjects_between = nothing -> effects between subjects, e.g. young vs old 
- items_between = nothing -> effects between items, e.g. natural vs artificial images, but shown to all subjects
- both_within = nothing	-> effects completly crossed
- tableModifyFun = x->x; # can be used to sort, or x->shuffle(MersenneTwister(42),x) - be sure to fix/update the rng accordingly!!

tipp: check the resulting dataframe using `generate(design)`

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
    subjects_between = nothing
    items_between = nothing
    both_within = nothing
    tableModifyFun = x -> x # can be used to sort, or x->shuffle(rng,x)
end


"""
- conditions = Dict of conditions, e.g. `Dict(:A=>["a_small","a_big"],:B=>["b_tiny","b_large"])`
- tableModifyFun = x->x; # can be used to sort, or x->shuffle(MersenneTwister(42),x) - be sure to fix/update the rng accordingly!!

Number of trials / rows in `generate(design)` depend on the full factorial of your `conditions`.

To increase the number of repetitions simply use `RepeatDesign(SingleSubjectDesign(...),5)`

tipp: check the resulting dataframe using `generate(design)`
"""
@with_kw struct SingleSubjectDesign <: AbstractDesign
    conditions = nothing
    tableModifyFun = x -> x
end


""" Returns dimension of experiment design"""
size(expdesign::MultiSubjectDesign) = (expdesign.n_items, expdesign.n_subjects)
size(expdesign::SingleSubjectDesign) = (*(length.(values(expdesign.conditions))...),)

"""
Generates full-factorial DataFrame of expdesign.conditions

Afterwards applies expdesign.tableModifyFun.

julia> d = SingleSubjectDesign(;conditions= Dict(:A=>nlevels(5),:B=>nlevels(2)))
julia> generate(d)
"""
function generate(expdesign::SingleSubjectDesign)
    # we get a Dict(:A=>["1","2"],:B=>["3","4"]), but needed a list
    # of named tuples for MixedModelsSim.factorproduct function.
    evts =
        factorproduct(((; k => v) for (k, v) in pairs(expdesign.conditions))...) |>
        DataFrame

    # by default does nothing
    return expdesign.tableModifyFun(evts)
end

"""
Generates full factorial Dataframe according to MixedModelsSim.jl 's simdat_crossed function
Note: n_items = you can think of it as `trials` or better, as stimuli

Note: No condition can be named `dv` which is used internally in MixedModelsSim / MixedModels as a dummy left-side

Afterwards applies expdesign.tableModifyFun.  Could be used to duplicate trials, sort, subselect etc.

Finally it sorts by `:subject`

julia> d = MultiSubjectDesign(;n_subjects = 10,n_items=20,both_within= Dict(:A=>nlevels(5),:B=>nlevels(2)))
julia> generate(d)
"""
function generate(expdesign::MultiSubjectDesign)
    #generate(expdesign::AbstractDesign) = generate(MersenneTwister(1),expdesign)

    # check that :dv is not in any condition
    allconditions =
        [expdesign.subjects_between, expdesign.items_between, expdesign.both_within]
    @assert :dv ∉ keys(merge(allconditions[.!isnothing.(allconditions)]...)) "due to technical limitations in MixedModelsSim.jl, `:dv` cannot be used as a factorname"


	data = DataFrame(
		MixedModelsSim.simdat_crossed(
			expdesign.n_subjects,
			expdesign.n_items,
			subj_btwn = expdesign.subjects_between,
			item_btwn = expdesign.items_between,
			both_win = expdesign.both_within,
		),
	)
	rename!(data, :subj => :subject)
	select!(data, Not(:dv)) # remove the default column from MixedModelsSim.jl - we don't need it in UnfoldSim.jl
	# by default does nothing
	data = expdesign.tableModifyFun(data)

	# sort by subject
	data = sort!(data, (order(:subject)))

	return data

end


# length is the same of all dimensions
length(expdesign::AbstractDesign) = *(size(expdesign)...)



# ----

"""
repeat a design DataFrame multiple times to mimick repeatedly recorded trials

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

function UnfoldSim.generate(design::RepeatDesign)
	df = map(x -> generate(design.design), 1:design.repeat) |> x -> vcat(x...)
	if isa(design.design, MultiSubjectDesign)
		sort!(df, [:subject])
	end
	return df

end
Base.size(design::RepeatDesign{MultiSubjectDesign}) =
	size(design.design) .* (design.repeat, 1)
Base.size(design::RepeatDesign{SingleSubjectDesign}) = size(design.design) .* design.repeat

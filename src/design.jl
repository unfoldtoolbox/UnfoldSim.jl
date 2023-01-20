""" 
n_subjects::Int -> number of subjects
n_items::Int -> number of items (=trials)
subjects_between = nothing -> effects between subjects, e.g. young vs old 
items_between = nothing -> effects between items, e.g. natural vs artificial images, but shown to all subjects
both_within = nothing	-> effects completly crossed
tableModifyFun = x->x; # can be used to sort, or x->shuffle(MersenneTwister(42),x) - be sure to fix/update the rng accordingly!!

tipp: check the resulting dataframe using `generate(design)`
"""
@with_kw struct MultiSubjectDesign <: AbstractDesign
    n_subjects::Int
    n_items::Int
    subjects_between = nothing
    items_between = nothing
    both_within = nothing
    tableModifyFun = x->x; # can be used to sort, or x->shuffle(rng,x)
end


"""

n_trials::Int -> number of trials
conditions = Dict of conditions, e.g. `Dict(:A=>["a_small","a_big"],:B=>["b_tiny","b_large"])`
tableModifyFun = x->x; # can be used to sort, or x->shuffle(MersenneTwister(42),x) - be sure to fix/update the rng accordingly!!

tipp: check the resulting dataframe using `generate(design)`
"""
@with_kw struct SingleSubjectDesign <: AbstractDesign
	n_trials::Int
	conditions = nothing
	tableModifyFun = x->x;
end


""" Returns dimension of experiment design"""
size(expdesign::MultiSubjectDesign) = (expdesign.n_items, expdesign.n_subjects)
size(expdesign::SingleSubjectDesign) = (expdesign.n_trials,)

"""
Generates full factorial DataFrame of expdesign.conditions x expdesign.n_trials.
Afterwards applies expdesign.tableModifyFun.

julia> d = SingleSubjectDesign(;n_trials = 10,conditions= Dict(:A=>nlevels(5),:B=>nlevels(2)))
julia> generate(d)
"""
function generate(expdesign::SingleSubjectDesign)
	# we get a Dict(:A=>["1","2"],:B=>["3","4"]), but needed a list
	# of named tuples for MixedModelsSim.factorproduct function.
	evts = factorproduct(
		(;trial=1:expdesign.n_trials),
		((;k=>v) for (k,v) in pairs(expdesign.conditions))...) |> DataFrame
	select!(evts,Not(:trial))
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
	data = DataFrame(
		MixedModelsSim.simdat_crossed(
			expdesign.n_subjects, 
			expdesign.n_items, 
			subj_btwn=expdesign.subjects_between, 
			item_btwn=expdesign.items_between, 
			both_win=expdesign.both_within
		)
	)
	rename!(data,:subj => :subject)
	
	# by default does nothing
	data = expdesign.tableModifyFun(data)
	
	# sort by subject
	data = sort!(data,(order(:subject)))

	return data
	
end


# length is the same of all dimensions
length(expdesign::AbstractDesign) = *(size(expdesign)...)
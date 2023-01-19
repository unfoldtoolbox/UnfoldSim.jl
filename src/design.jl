""" 
n_subj::Int -> number of subjects
n_item::Int -> number of items (=trials)
subj_btwn = nothing -> effects between subjects, e.g. young vs old 
item_btwn = nothing -> effects between items, e.g. natural vs artificial images, but shown to all subjects
both_win = nothing	-> effects completly crossed
tableModifyFun = x->x; # can be used to sort, or x->shuffle(MersenneTwister(42),x) - be sure to fix/update the rng accordingly!!

tipp: check the resulting dataframe using `generate(design)`
"""
@with_kw struct MultiSubjectDesign <: AbstractDesign
    n_subj::Int
    n_item::Int
    subj_btwn = nothing
    item_btwn = nothing
    both_win = nothing
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
size(expdesign::MultiSubjectDesign) = (expdesign.n_item, expdesign.n_subj)
size(expdesign::SingleSubjectDesign) = (expdesign.n_trials,)
"""
    generate(expdesign::ExperimentDesign)

Generates experiment design data frame
"""
function generate(expdesign::SingleSubjectDesign)
	tmpdesign = MultiSubjectDesign(;
	n_subj=2,
	n_item=expdesign.n_trials,
	item_btwn=expdesign.conditions,
	tableModifyFun=expdesign.tableModifyFun)
	evts = generate(tmpdesign)
	return evts[evts.subject .== "S1",Not(["subject","dv","item"])]
end
function generate(expdesign::MultiSubjectDesign)
	#generate(expdesign::AbstractDesign) = generate(MersenneTwister(1),expdesign)
	data = DataFrame(
		MixedModelsSim.simdat_crossed(
			expdesign.n_subj, 
			expdesign.n_item, 
			subj_btwn=expdesign.subj_btwn, 
			item_btwn=expdesign.item_btwn, 
			both_win=expdesign.both_win
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
"""
Experiment Design
"""
@with_kw struct MultiSubjectDesign <: AbstractDesign
    n_subj::Int
    n_item::Int
    subj_btwn = nothing
    item_btwn = nothing
    both_win = nothing
    tableModifyFun = x->x; # can be used to sort, or x->permute(rng,x)
end

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
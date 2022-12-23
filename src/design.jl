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



""" Returns dimension of experiment design"""
dims(expdesign::MultiSubjectDesign) = expdesign.n_subj * expdesign.n_item


"""
    generate(expdesign::ExperimentDesign)

Generates experiment design data frame
"""

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
	# by default does nothing
	data = expdesign.tableModifyFun(data)
	
	# sort by subject
	data = sort!(data,(order(:subj)))

	return data
	
end

"""
Experiment Design
"""
mutable struct ExperimentDesign
    n_subj::Int
    n_item::Int
    subj_btwn::Any
    item_btwn::Any
    both_win::Any
end


""" Returns dimension of experiment design"""
dims(expdesign::ExperimentDesign) = expdesign.n_subj * expdesign.n_item


"""
    generate(expdesign::ExperimentDesign)

Generates experiment design data frame
"""
function generate(expdesign::ExperimentDesign)
	return sort!(DataFrame(
		MixedModelsSim.simdat_crossed(
			expdesign.n_subj, 
			expdesign.n_item, 
			subj_btwn=expdesign.subj_btwn, 
			item_btwn=expdesign.item_btwn, 
			both_win=expdesign.both_win
		)
	))
end

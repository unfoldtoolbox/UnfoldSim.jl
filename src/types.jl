abstract type AbstractOnset end
abstract type AbstractNoise end

# find other types in onset.jl and noise.jl
"""
Experiment Design
"""
struct ExperimentDesign
    n_subj::Int
    n_item::Int
    subj_btwn
    item_btwn
    both_win
end




"""
ERP Component
"""
struct Component
	basis
	formula
	contrasts
	β
	σ_ranef
	σ_res
end

struct Simulation
	design::ExperimentDesign
	components::Vector{Component}
	onset::AbstractOnset
	noisetype::AbstractNoise
end


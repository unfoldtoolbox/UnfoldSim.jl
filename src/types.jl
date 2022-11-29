abstract type AbstractOnset end
abstract type AbstractNoise end

abstract type AbstractDesign end

abstract type AbstractComponent end

# find other types in onset.jl and noise.jl
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




"""
ERP Component
"""
struct Component <: AbstractComponent
	basis
	formula
	contrasts
	β
	σ_ranef
	σ_res
end

struct Simulation
	design::AbstractDesign
	components::Vector{AbstractComponent}
	onset::AbstractOnset
	noisetype::AbstractNoise
end


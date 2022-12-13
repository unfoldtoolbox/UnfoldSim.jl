
using Test
using Random
using StableRNGs
using Statistics

function gen_debug_design()
	# define design parameters
	n_subj = 20
	n_item = 100

	btwn_item = Dict("stimType" => ["A", "B"])
		
	# instantiate the design
	return MultiSubjectDesign(;n_subj=n_subj, n_item=n_item, btwn_item = btwn_item)
end

function gen_debug_component()
    basisfunction = [zeros(10) ones(10) zeros(10)]
    formula = @formula(dv ~ 1 + stimType + (1 + stimType | subj))
    contrasts = Dict(:stimType => DummyCoding())
    β = [2.0, 0.5]
    σ_ranef = Dict(:subj => create_re(1, 0.0))
    σ_res = 0.0001
    
    # instantiate the component(s)
    return Component(basisfunction, formula, contrasts, β, σ_ranef, σ_res)
end

function gen_debug_simulation(;design = gen_debug_design(),component = gen_debug_component(),noisetype = PinkNoise(),onset = UniformOnset(;width=0,offset=100))
	return Simulation(design, component, onset,noisetype)
end

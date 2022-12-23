"""
MixedModelComponent
works best with MultiSubjectDesign
"""
@with_kw struct MixedModelComponent <: AbstractComponent
    basis
    formula # e.g. dv~1+cond - left side must be "dv"
    β::Vector
    σs::Dict # Dict(:subj=>[0.5,0.4]) or to specify correlationmatrix Dict(:subj=>[0.5,0.4,I(2,2)],...)
    contrasts::Dict = Dict()
end

@with_kw struct LinearModelComponent <: AbstractComponent
    basis
    formula # e.g. 0~1+cond - left side must be "0"
    β::Vector
    contrasts::Dict = Dict()
end


Base.length(c::AbstractComponent) = length(c.basis)
maxlength(c::Vector{AbstractComponent}) = maximum(length.(c))

"""
simulate a linearModel

julia> c = UnfoldSim.LinearModelComponent([0,1,1,0],@formula(0~1+cond),[1,2],Dict())
julia> design = MultiSubjectDesign(;n_subj=2,n_item=50,item_btwn=(;:cond=>["A","B"]))
julia> simulate(StableRNG(1),c,design)
"""
function simulate(rng,c::LinearModelComponent,design::AbstractDesign)
    evts = generate(design)
    if isempty(c.contrasts)
        m = StatsModels.ModelFrame(c.formula, evts)
    else
        m = StatsModels.ModelFrame(c.formula, evts;contrasts=c.contrasts)
    end
    X = StatsModels.modelmatrix(m)    
    y = X * c.β
    return y' .* c.basis
end
"""
simulate MixedModelComponent

julia> design = MultiSubjectDesign(;n_subj=2,n_item=50,item_btwn=(;:cond=>["A","B"]))
julia> c = UnfoldSim.MixedModelComponent([0.,1,1,0],@formula(dv~1+cond+(1|subj)),[1,2],Dict(:subj=>[2],),Dict())
julia> simulate(StableRNG(1),c,design)

"""
function simulate(rng,c::MixedModelComponent,design::AbstractDesign)
	evts = generate(design)

	# create dummy
    if isempty(c.contrasts)
        m = MixedModels.MixedModel(c.formula, evts)
    else
	    m = MixedModels.MixedModel(c.formula, evts; contrasts=c.contrasts)
    end

	# limit runtime (in seconds)
	m.optsum.maxtime = 1
	m.optsum.feval = 1 # also do only a single step; 

	# fit mixed model to experiment design and dummy data
	refit!(m, progress=false) # to initialize it
	

	# empty epoch data
	epoch_data_component = zeros(Int(length(c.basis)), dims(design))

	# residual variance for lmm
	σ_lmm = 	0.0001
		
	# iterate over each timepoint
	for t in eachindex(c.basis)

			# select weight from basis
            # right now, it is the same, but maybe changein thefuture?
			basis_β  = c.basis[t]
			basis_σs = c.basis[t]
			
			
            # weight random effects by the basis function
            namedre = weight_σs(c.σs,basis_σs,σ_lmm)
            MixedModelsSim.update!(m; namedre...)
        
			# simulate with new parameters; will update m.y
			simulate!(
				deepcopy(rng), # same parameter for each timepont
				m, 
				β = basis_β .* c.β, # weight the beta by the basisfunction
				σ = σ_lmm # no noise
			)

			# save data to array
			epoch_data_component[t, :] = m.y
		end

		return epoch_data_component
	
end


"""
Weights a σs Dict for MixedModels.jl by a Float64

Finally sales it by σ_lmm, as a trick to simulate noise-free LMMs

I anticipate a function
    `function weight_σs(σs::Dict,b_σs::Dict,σ_lmm::Float64)`
where each σs entry can be weighted individually
"""
function weight_σs(σs::Dict,b_σs::Float64,σ_lmm::Float64)
    #k = (collect(keys(σs))...,)
    #val = values(σs)

    keys = Symbol[]
    vals = LowerTriangular[]

    for (k,v) in σs

        scale = (x)-> b_σs./σ_lmm .* x
        
        if v[end] isa Matrix
            v = create_re.(scale(v[1:end-1])...;corrmat=v[end])    
        else
            v = create_re.(scale(v)...;)
        end
        
        push!(keys,k)
        push!(vals,v)
    end
    
    namedre = NamedTuple(keys.=>vals)
    
   return namedre
end

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

Base.length(c::Component) = length(c.basis)
maxlength(c::Vector{Component}) = maximum(length.(c))


"""
simulate component?
"""

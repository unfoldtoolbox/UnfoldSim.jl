
Base.length(c::Component) = length(c.basis)
maxlength(c::Vector{Component}) = maximum(length.(c))


"""
simulate component?
"""

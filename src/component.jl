
Base.length(c::Component) = length(c.basis)
maxlength(c::Vector{AbstractComponent}) = maximum(length.(c))


"""
simulate component?
"""

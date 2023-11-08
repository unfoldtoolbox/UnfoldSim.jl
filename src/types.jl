abstract type AbstractOnset end
abstract type AbstractNoise end

abstract type AbstractDesign end

abstract type AbstractComponent end

# find other types in onset.jl and noise.jl
# and in design.jl and component.jl

abstract type AbstractHeadmodel end


struct Simulation
	design::AbstractDesign
	components::Vector{AbstractComponent}
	onset::AbstractOnset
	noisetype::AbstractNoise
end


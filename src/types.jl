abstract type AbstractOnset end
abstract type AbstractNoise end

abstract type AbstractDesign end

abstract type AbstractComponent end

# find other types in onset.jl and noise.jl
# and in design.jl and component.jl

abstract type AbstractHeadmodel end


struct Simulation
    design::AbstractDesign
    components::Union{
        <:Dict{<:Char,<:Vector{<:AbstractComponent}},
        <:Vector{<:AbstractComponent},
    }
    onset::AbstractOnset
    noisetype::AbstractNoise
end


Simulation(
    design::AbstractDesign,
    components::Dict{<:Char,<:Vector},
    onset::AbstractOnset,
    noisetype::AbstractNoise,
) = Simulation(design, Dict{Char,Vector{<:AbstractComponent}}(components), onset, noisetype)

"""
"""
abstract type AbstractOnset end
"""
"""
abstract type AbstractNoise end
"""
"""
abstract type AbstractDesign end
"""
"""
abstract type AbstractComponent end

# find other types in onset.jl and noise.jl
# and in design.jl and component.jl
"""
"""
abstract type AbstractHeadmodel end

"""
    Simulation

A type to store all "ingredients" for a simulation including their parameters.

Can either be created by the user or will be created automatically when calling the [`simulate`](@ref) function with the required "ingredients".
Tip: Use the `subtypes` function to get an overview of the implemented "ingredients", e.g. `subtypes(AbstractDesign)`.

# Fields
- `design::AbstractDesign`: Experimental design.
- `components::Vector{AbstractComponent}`: Response function(s) for the events (e.g. the ERP shape in EEG).
- `onset::AbstractOnset`: Inter-onset distance distribution.
- `noisetype::AbstractNoise`: Noise type.

# Examples
```julia-repl
julia> design = SingleSubjectDesign(; conditions = Dict(:stimulus_type => ["natural", "artificial"]));

julia> component = LinearModelComponent(;
           basis = p100(),
           formula = @formula(0 ~ 1 + stimulus_type),
           β = [2, 3],
       );

julia> onset = UniformOnset();

julia> noise = PinkNoise();

julia> simulation = Simulation(design, component, onset, noise);

julia> using StableRNGs

julia> data, events = simulate(StableRNG(1), simulation);

julia> events
2×2 DataFrame
 Row │ stimulus_type  latency 
     │ String         Int64   
─────┼────────────────────────
   1 │ natural             20
   2 │ artificial          70

julia> data
85-element Vector{Float64}:
  0.8596261232522926
  1.1657369535500595
  0.9595228616486761
  ⋮
  0.9925202143746904
  0.2390652543395527
 -0.11672523788068771
```
"""
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

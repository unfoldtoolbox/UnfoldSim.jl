abstract type AbstractOnset end
abstract type AbstractNoise end

abstract type AbstractDesign end

abstract type AbstractComponent end

# find other types in onset.jl and noise.jl
# and in design.jl and component.jl

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
    components::Vector{AbstractComponent}
    onset::AbstractOnset
    noisetype::AbstractNoise
end


# Eye Movement
abstract type AbstractControlSignal{T} end

struct HREFCoordinates{T} <: AbstractControlSignal{T}
    val::Matrix{T}
end

struct GazeDirectionVectors{T} <: AbstractControlSignal{T}
    val::Matrix{T} # size 3 x n_timepoints
end



# Artifacts
"""
TODO docstring
`controlsignal`: Defines the control signal, always starting from the first time point. 
"""
abstract type AbstractContinuousSignal end
# in future we may want to allow generating the noise on only a particular channel.

@with_kw struct EyeMovement{T} <: AbstractContinuousSignal
    controlsignal::T
    headmodel
    eye_model::String = "crd"
    # events # <-- from realdata (or from controlsignal?) or passed in by user. will be added into the events dataframe returned by simulation function
end

struct TRF <: AbstractContinuousSignal
    controlsignal
    # TBD
end

# harmonics are always weighted the same relative to each other 
@with_kw struct PowerLineNoise <: AbstractContinuousSignal
    controlsignal = nothing
    base_freq::Float64 = 50
    harmonics::Array{Int64} = [1 3 5]
    weights_harmonics::Array{Float64} = ones(length(harmonics))
    sampling_rate::Float64 = 1000
end

@with_kw struct UserDefinedContinuousSignal <: AbstractContinuousSignal
    controlsignal = nothing
    signal::Array
end

@with_kw struct DCDriftNoise <: AbstractContinuousSignal
    controlsignal = nothing
    scaling_factor = 1
end
@with_kw struct ARDriftNoise <: AbstractContinuousSignal
    controlsignal = nothing
    σ
end
@with_kw struct LinearDriftNoise <: AbstractContinuousSignal
    controlsignal = nothing
    scaling_factor = 1
end

@with_kw struct DriftNoise <: AbstractContinuousSignal
    ar::ARDriftNoise = ARDriftNoise(σ=1)
    linear::LinearDriftNoise = LinearDriftNoise()
    dc::DCDriftNoise = DCDriftNoise()
end

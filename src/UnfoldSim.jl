module UnfoldSim

   using DSP
   using Random
   using DataFrames
   using Parameters
   using StatsModels
   using MixedModels
   using ImageFiltering
   using MixedModelsSim
   using SignalAnalysis

   import Base.length

   include("design.jl")
   include("component.jl")
   include("noise.jl")
   include("simulation.jl")
    include("onset.jl")

   export 
    @formula,
    EffectsCoding, 
    ExperimentDesign,
    generate,
    dims,
    Component,
    PinkNoise,
    RedNoise,
    WhiteNoise,
    NoNoise,
    gen_noise,
    Simulation,
    simulate,
    simulate_erps,
    padarray,
    convert,
    UniformOnset


    export create_re
end

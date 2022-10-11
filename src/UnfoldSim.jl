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

   include("design.jl")
   include("component.jl")
   include("noise.jl")
   include("simulation.jl")


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
    gen_noise,
    Simulation,
    simulate,
    simulate_erps,
    padarray,
    convert
end
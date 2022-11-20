module UnfoldSim

   using DSP
   using Random
   using DataFrames
   using Distributions # for LogNormal Onset
   using Parameters
   using StatsModels
   using MixedModels
   using ImageFiltering # for Noise-filter (can be replaced maybe?)
   using MixedModelsSim
   using SignalAnalysis 

   import Base.length
   include("types.jl")
   include("design.jl")
   include("component.jl")
   include("noise.jl")
   include("simulation.jl")
   include("onset.jl")

   
   # statsmodels re-export
   export @formula,DummyCoding,EffectsCoding
   # mixedModels re-export
   export create_re
   # main types
   export ExperimentDesign,Simulation, Component
   # noise functions

   export PinkNoise,RedNoise,WhiteNoise,NoNoise #,RealNoise (not implemented yet)
   # UnfoldSim functions 
   export simulate, simulate_erps,gen_noise,padarray,convert

   # export Offsets
   export UniformOnset,LogNormalOnset
   
end

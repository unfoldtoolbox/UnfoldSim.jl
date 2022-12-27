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
   using LinearAlgebra

   import Base.length
   import Base.size
   include("types.jl")
   include("design.jl")
   include("component.jl")
   include("noise.jl")
   include("simulation.jl")
   include("onset.jl")

   export AbstractComponent,AbstractNoise,AbstactOnset,AbstractDesign
   # statsmodels re-export
   export @formula,DummyCoding,EffectsCoding
   # mixedModels re-export
   export create_re
   
   # main types
   export Simulation
   
   # component types
   export MixedModelComponent,LinearModelComponent

   # export designs
   export MultiSubjectDesign,SingleSubjectDesign

   # noise functions
   export PinkNoise,RedNoise,WhiteNoise,NoNoise,ExponentialNoise #,RealNoise (not implemented yet)
   # UnfoldSim functions 
   export simulate, gen_noise
   
   # utilities
   export padarray,convert

   # export Offsets
   export UniformOnset,LogNormalOnset
   
end

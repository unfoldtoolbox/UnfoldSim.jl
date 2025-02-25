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
using ToeplitzMatrices # for AR Expo. Noise "Circulant"
using StatsModels
using HDF5, Artifacts, FileIO
using Automa # for sequence

using LinearAlgebra # headmodel

using SequentialSamplingModels # for SequentialSamplingModels

import DSP.hanning
import Base.length
import Base.size
import Base.show
include("types.jl")
include("design.jl")
include("component.jl")
include("noise.jl")
include("simulation.jl")
include("onset.jl")
include("predefinedSimulations.jl")
include("headmodel.jl")
include("helper.jl")
include("sequence.jl")
include("bases.jl")
include("sequentialSamplingModelSimulation.jl")

export size, length
export AbstractComponent, AbstractNoise, AbstractOnset, AbstractDesign
# statsmodels re-export
export @formula, DummyCoding, EffectsCoding
# mixedModels re-export
export create_re

# main types
export Simulation

# component types
export MixedModelComponent, LinearModelComponent

# export designs
export MultiSubjectDesign, SingleSubjectDesign, RepeatDesign, SequenceDesign

# noise functions
export PinkNoise, RedNoise, WhiteNoise, NoNoise, ExponentialNoise #,RealNoise (not implemented yet)

# UnfoldSim functions 
export simulate,
    simulate_responses,
    simulate_interonset_distances,
    simulate_component,
    simulate_onsets,
    simulate_noise,
    generate_events

# utilities
export pad_array, convert

# export Offsets
export UniformOnset, LogNormalOnset, NoOnset, UniformOnsetFormula, LogNormalOnsetFormula

# re-export StatsModels
export DummyCoding, EffectsCoding

# export bases
export p100, n170, p300, n400, hrf, PuRF

# headmodel
export AbstractHeadmodel, Hartmut, headmodel, leadfield, orientation, magnitude

# multichannel
export MultichannelComponent

# evidence accumulation
export DriftComponent, DriftOnset, SequenceOnset, KellyModel
end

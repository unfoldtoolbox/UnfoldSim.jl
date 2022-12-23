using UnfoldSim
using CairoMakie
using DSP
using StableRNGs
using DataFrames
import StatsBase.autocor
using Parameters
# ## Define a new Design
# A design specifies how much data is generated, and how the event-table(s)
# should be generated. Implemented examples are `MultiSubjectDesign`, (t.b.d.) `SingleSubjectdesign`
#
# We need 4 things, a `struct<:AbstractDesign`, a `dims` function and a `generate` function
# 1) We need a `ImbalanceSubjectDesign` struct
#
@with_kw struct ImbalanceSubjectDesign <: UnfoldSim.AbstractDesign
    nTrials::Int
    balance::Float64 = 0.5 # default balanced    
end
##
# 2) we need a `dims(design::ImbalanceSubjectDesign)` function to tell how many events we will have
dims(design::ImbalanceSubjectDesign) = design.nTrials

# We need a type `generate(design::ImbalanceSubjectDesign)` 
function generate(design::ImbalanceSubjectDesign)
    nA = Int(round.(design.nTrials .* balance))
    nB = Int(round.(design.nTrials .* (1-balance)))
    @assert nA + nB  ≈ design.nTrials
    levels = vcat(repeat(["levelA"],nA),repeat(["levelB"],nB))
    return DataFrame(Dict(:condition=>levels))
end

# now we can test the function
design = ImbalanceSubjectDesign(;nTrials=6,balance=0.2)
generate(design)

# !!! important
#   it is the users task to ensure that each run is reproducible. So if you have a random process (e.g. shuffling), be sure to 
#   safe a RNG object in your struct and use it in your generate function.
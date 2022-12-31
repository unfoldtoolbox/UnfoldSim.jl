using UnfoldSim
using StableRNGs
using DataFrames
using Parameters
# ## Define a new Design
# A design specifies how much data is generated, and how the event-table(s)
# should be generated. Already implemented examples are `MultiSubjectDesign` and `SingleSubjectdesign`
#
# We need 3 things for a new design: a `struct<:AbstractDesign`, a `size` and a `generate` function
#
# #### 1) `type`
# We need a `ImbalanceSubjectDesign` struct. You are free to implement it as you wish, as long as the other two functions are implemented
#
@with_kw struct ImbalanceSubjectDesign <: UnfoldSim.AbstractDesign
    nTrials::Int
    balance::Float64 = 0.5 # default balanced    
end;

# #### 2) `size`
# we need a `size(design::ImbalanceSubjectDesign)` function to tell how many events we will have.
# This is used at different places, e.g. in the Default onset implementation

## note the trailling , to make it a Tuple
size(design::ImbalanceSubjectDesign) = (design.nTrials,);

# #### 3) `generate`
# We need a type `generate(design::ImbalanceSubjectDesign)` function. This function should return the actual table as a `DataFrame`
function generate(design::ImbalanceSubjectDesign)
    nA = Int(round.(design.nTrials .* design.balance))
    nB = Int(round.(design.nTrials .* (1-design.balance)))
    @assert nA + nB  ≈ design.nTrials
    levels = vcat(repeat(["levelA"],nA),repeat(["levelB"],nB))
    return DataFrame(Dict(:condition=>levels))
end;

# Finally, we can test the function and see whether it returns a Design-DataFrame as we requested
design = ImbalanceSubjectDesign(;nTrials=6,balance=0.2)
generate(design)

# !!! important
#       it is the users task to ensure that each run is reproducible. So if you have a random process (e.g. shuffling), be sure to 
#       safe a RNG object in your struct and use it in your generate function.
# # Define a new (imbalanced) design

# A design specifies how much data is generated, and how the event-table(s)
# should be generated. Already implemented examples are `MultiSubjectDesign` and `SingleSubjectDesign`.


# We need 3 things for a new design: a `struct<:AbstractDesign`, a `size` and a `generate_events` function.


# ### Setup
# ```@raw html
# <details>
# <summary>Click to expand</summary>
# ```
using UnfoldSim
using StableRNGs
using DataFrames
using Parameters
using Random

# ```@raw html
# </details>
# <br />
# ```

# #### 1) `type`
# We need a `ImbalanceSubjectDesign` struct. You are free to implement it as you wish, as long as the other two functions are implemented
#
@with_kw struct ImbalanceSubjectDesign <: UnfoldSim.AbstractDesign
    n_trials::Int
    balance::Float64 = 0.5 # default balanced    
end;

# #### 2) `size`
# we need a `size(design::ImbalanceSubjectDesign)` function to tell how many events we will have.
# This is used at different places, e.g. in the Default onset implementation

## note the trailing , to make it a Tuple
UnfoldSim.size(design::ImbalanceSubjectDesign) = (design.n_trials,);

# #### 3) `generate_events`
# We need a type `generate_events(rng::AbstractRNG, design::ImbalanceSubjectDesign)` function. This function should return the actual table as a `DataFrame`
function UnfoldSim.generate_events(rng::AbstractRNG, design::ImbalanceSubjectDesign)
    nA = Int(round.(design.n_trials .* design.balance))
    nB = Int(round.(design.n_trials .* (1 - design.balance)))
    @assert nA + nB ≈ design.n_trials
    levels = vcat(repeat(["levelA"], nA), repeat(["levelB"], nB))
    return DataFrame(Dict(:condition => levels))
end;

# Finally, we can test the function and see whether it returns a Design-DataFrame as we requested
design = ImbalanceSubjectDesign(; n_trials = 6, balance = 0.2)
generate_events(design)

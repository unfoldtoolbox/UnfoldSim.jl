# # Define a new component (with variable duration and shift)

# We want a new component that changes its duration and shift depending on a column in the event design. This is somewhat already implemented in the HRF + Pupil bases.
# !!! hint
#     if you are just interested to use duration-dependency in your simulation, check out the component-function tutorial


# ### Setup
# ```@raw html
# <details>
# <summary>Click to expand</summary>
# ```
using UnfoldSim
import UnfoldSim.simulate_component
using Unfold
using Random
using DSP
using CairoMakie, UnfoldMakie

sfreq = 100;
# ```@raw html
# </details >
# ```

# ## Design
# Let's generate a design with two columns, shift + duration
design = UnfoldSim.SingleSubjectDesign(;
    conditions = Dict(
        :shift => rand(100) .* sfreq / 5,
        :duration => 20 .+ rand(100) .* sfreq / 5,
    ),
)


# ## Implement a new AbstractComponent
# We also need a new `AbstractComponent`
struct TimeVaryingComponent <: AbstractComponent
    basisfunction::Any
    maxlength::Any
end

# We have to define the length of a component 
Base.length(c::TimeVaryingComponent) = length(c.maxlength)

# While we could have put the TimeVaryingComponent.basisfunction directly into the simulate function, I thought this is a bit more modular
function UnfoldSim.simulate_component(rng, c::TimeVaryingComponent, design::AbstractDesign)
    evts = generate_events(design)
    return c.basisfunction(evts, c.maxlength)
end

# finally, the actual function that does the shifting + duration
function basis_shiftduration(evts, maxlength)
    basis = hanning.(Int.(round.(evts.duration))) ## hanning as long as duration
    if "shift" ∈ names(evts)
        basis = pad_array.(basis, Int.(round.(.-evts.shift)), 0) ## shift by adding 0 in front
    end
    ## we should make sure that all bases have maxlength by appending / truncating
    difftomax = maxlength .- length.(basis)
    if any(difftomax .< 0)
        @warn "basis longer than max length in at least one case. either increase maxlength or redefine function. Trying to truncate the basis"
        basis[difftomax.>0] = pad_array.(basis[difftomax.>0], difftomax[difftomax.>0], 0)
        return [b[1:maxlength] for b in basis]
    else
        return pad_array.(basis, difftomax, 0)
    end
end

# ## Simulate data with the new component type
erp = UnfoldSim.simulate_component(
    MersenneTwister(1),
    TimeVaryingComponent(basis_shiftduration, 50),
    design,
)
plot_erpimage(hcat(erp...), sortvalues = generate_events(design).shift)

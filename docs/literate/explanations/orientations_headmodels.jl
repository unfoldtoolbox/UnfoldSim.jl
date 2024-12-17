# # On orientations in Leadfields for simulation

# This is a somewhat loose organisation of some thoughts we had regarding leadfields, simulations with a focus on orientation of sources

# ### Setup
# ```@raw html
# <details>
# <summary>Click to expand</summary>
# ```
## Load required packages

using CairoMakie
using UnfoldMakie
using UnfoldSim
using LinearAlgebra: norm
# ```@raw html
# </details >
# ```

# ### The leadfield
hart = headmodel()
pos = to_positions(hart.electrodes["pos"]')
L = leadfield(hart)
size(L)

# The leadfield describes the connection between source-points and electrodes. It contains the precalculated forward simulation of the electrical potentials.
# As seen from the size, we actually have three leadfields, one for `x`, one for `y `and one for the `z` direction.
#
# The question we will explore is: Which direction should we simulate. 
#
# Let's first plot the three orientations so we know what we are talking about
f = Figure()
ix = 50
[
    plot_topoplot!(
        f[1, k],
        L[:, ix, k];
        positions = pos,
        visual = (; label_scatter = false),
    ) for k = 1:3
]
f

# These are the simulated source activatios of the source-point `$(ix)`, if the underlying dipolar source would be oriented in x,y, or z direction.

# But often, we do not want that, but want to simulate the "best" direction. But how to get that?

# Luckily, often headmodellers provide orientations orthogonal to the cortex, in direction of the pyramidal cells - a physiological plausible prior.
# This is the default we use 

m_perp = magnitude(hart)
plot_topoplot!(
    f[2, 2],
    m_perp[:, ix];
    positions = pos,
    visual = (; label_scatter = false),
    axis = (; title = "perpendicular"),
)
f

# We could also calculate the vector norm:

m_norm = mapslices(norm, L; dims = (3))
ix = 50
plot_topoplot!(
    f[2, 1],
    m_norm[:, ix, 1];
    positions = pos,
    visual = (; label_scatter = false),
    axis = (; title = "norm"),
)
f

# other options are sum, or maximum, or ...
plot_topoplot!(
    f[2, 3],
    sum(L[:, ix, :], dims = 2)[:, 1];
    positions = pos,
    visual = (; label_scatter = false),
    axis = (; title = "sum"),
)
plot_topoplot!(
    f[3, 3],
    maximum(L[:, ix, :], dims = 2)[:, 1];
    positions = pos,
    visual = (; label_scatter = false),
    axis = (; title = "max"),
)
f

# Which is the best? We don't know! If you have good ideas & your own headmodels please tell us. For now we recommend just to use "perpendicular" :)
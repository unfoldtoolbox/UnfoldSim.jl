# # Power analysis
# For a power analysis, we will repeatedly simulate data, and check whether we can find a significant effect.
# We perform the power analysis on epoched data.

# ### Setup
# ```@raw html
# <details>
# <summary>Click to expand</summary>
# ```
## Load required packages
using UnfoldSim
using Statistics
using HypothesisTests
using DataFrames
using Random
# ```@raw html
# </details >
# <br />
# ```

# ## Simulation loop
pvals = fill(NaN, 100)
@time for seed in eachindex(pvals)
    ## Simulate data of 30 subjects
    data, evts = UnfoldSim.predef_2x2(
        MersenneTwister(seed);
        n_subjects = 20, ## 30 subjects
        overlap = (1, 0), ## deactivate overlap
        noiselevel = 10,  ## add more noise to make it more challenging
        return_epoched = true, ## saves us the epoching step
    )


    ## take the mean over a pre-specified timewindow
    evts.y = dropdims(mean(data[40:60, :, :], dims = 1), dims = (1))[:]

    ## extract the two levels of condition A
    evts_reduced = combine(groupby(evts, [:subject, :A]), :y => mean)
    y_big = evts_reduced[evts_reduced.A.=="a_big", :y_mean]
    y_small = evts_reduced[evts_reduced.A.=="a_small", :y_mean]

    ## calculate a one-sided t-test
    pvals[seed] = pvalue(OneSampleTTest(y_big, y_small))
end

# Let's calculate the power
power = mean(pvals .< 0.05) * 100

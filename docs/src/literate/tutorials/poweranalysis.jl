using UnfoldSim
using Unfold
using Statistics
using HypothesisTests
using DataFrames
using Random

# ## Simple Poweranalysis Script
# For a power analysis, we will repeatedly simulate data, and check whether we can find a significant effect.
# 

pvals = fill(NaN,100)
for seed = eachindex(pvals)
    # Simulate data of 30 subjects
    data,evts = UnfoldSim.predef_2x2(MersenneTwister(seed);
                n_subjects=30, # 30 subjects
                overlap=(1,0), # deactivate overlap
                noiselevel=4,  # add more noise to make it more challenging
                )

    # cut the continuous simulated signal to epochs
    data_e,times = Unfold.epoch(;data=data,tbl=evts,τ=(0.4,0.6),sfreq=100)

    # take the mean over a pre-specified timewindow
    evts.y = dropdims(mean(data_e,dims=2),dims=(1,2))
    
    # extract the two levels of condition A
    evts_reduced = combine(groupby(evts,[:subject,:A]),:y=>mean)
    y_big = evts_reduced[evts_reduced.A .=="a_big",:y_mean]
    y_small = evts_reduced[evts_reduced.A .=="a_small",:y_mean]

    # calculate a one-sided t-test
    pvals[seed] = pvalue(OneSampleTTest(y_big,y_small))
end
power = mean(pvals .<0.05)*100;

# the calculated power is $power%%
using UnfoldSim
using Unfold
using CairoMakie
using UnfoldMakie
using DataFrames

# # Multi-Subject simulation
# Similar to the single subject case, multi-subject simulation depends on:
# - `Design` (typically a `MultiSubjectDesign`)
# - `Components` (typically a `MixedModelComponent`)
# - `Onset` (any)
# - `Noise` (any)

# ## Design

# Our first design should be 20 subjects, with 3 items each. Any individual image is shown  only either as large or small, thus we choose `items_between`. 
design = MultiSubjectDesign(
    n_subjects = 20,
    n_items = 4,
    items_between = Dict(:condition => ["large", "small"]),
)


# ## Between, within?
# In the beginning, the distinction between `between-items`, `between-subjects` and `within-subjects`, `within-items` and `both-between`, `both-within` feels daunting.
# 
# We base our terminology on `MixedModelsSim` which uses the following definitions:
# - `subjects_between` -> effects between subjects, e.g. young vs old
# - `items_between` -> effects between items, e.g. natural vs artificial images, (but shown to all subjects if not specified in subjects_between as well)
# - `both_within` -> effects completly crossed, e.g. word vs. scramble, where the "original" word is the item, and shown to all subjects 


# ## Components
# For multi-subject, similar to the `LinearModelComponent` specified before, we have to define the fixed effect `β`, the model parameters that are applied to all subjects.
β = [1, 2] # 1 = intercept, 2 = difference between large and small

# In addition, we have to provide random effects `σs`, which define the spread (and  correlation) of the subjects around the fixed effects, foreach parameter
σs = Dict(
    :subject => [0.5, 1], # we have more spread in the condition-effect
    :item => [1], # the item-variability is higher than the subject-variability 
)

# now we are ready to build it together
signal = MixedModelComponent(;
    basis = UnfoldSim.hanning(50),
    formula = @formula(0 ~ 1 + condition + (1 + condition | subject) + (1 | item)),
    β = β,
    σs = σs,
    contrasts = Dict(:condition => EffectsCoding()), # we highly recommend specifying your contrasts, by Default its Dummy/ReferenceCoding with alphabetically sorted levels (relying 100% on StatsModels.jl)
)

# and simulate!
data, evts = simulate(design, signal, NoOnset(), NoNoise(), return_epoched = true);

# We get data with 50 samples (our `basis` from above), with `2` trials/items and 20 subjects. We get items and subjects separately because we chose no-overlap (via `NoOnset`) and `return_epoched=true``.
size(data)

first(evts, 5)

# Finally, let's plot the data
f = Figure()

for k = 1:4
    series(
        f[1, k],
        data[:, k, :]',
        solid_color = :black,
        axis = (; limits = ((0, 50), (-4, 5))),
    )
    Label(f[1, k, Top()], text = "Item:" * evts[k, :item] * ", c:" * evts[k, :condition])
end
f

# Some remarks on interpreting the plot:
# - The β main-effect of small / Large (#1 and #3 vs. #2 and #4) is clearly visible. 
# - The variability between subjects, is the variability between the individual curves.
# - The item effect shows up e.g. that #2 vs. #4 column show different values.


# # Continuous Signals / Overlap
# Let's continue our tutorial and simulate overlapping signals instead.
#
# We replace the `NoOnset` with an `UniformOnset` between 20 and 70 sampples after each event.` We further remove the `return_epoched`, because we want to have continuous data for now.
data, evts = simulate(design, signal, UniformOnset(offset = 20, width = 50), NoNoise());
size(data)

# The data is now $size(data,1) x $size(data,2), with the first dimension being continuous data, and the latter still the subjects.
series(data', solid_color = :black)

# Each line is one subject, and it looks a bit unstructured, because the event-onsets are of course random or each subject.
# !!! note
#     All subjects have the same sequence of trials, if you need to change this, specify a `event_order_function` in the `MultiSubjectDesign`


# # Analyzing these data with Unfold.jl
# We will analyze these data using the `Unfold.jl` toolbox. While preliminary support for deconvolution (overlap correction) for mixed models is available, here we will simply apply it separately to each timepoint, following the MassUnivariate approach.
data, evts = simulate(
    design,
    signal,
    UniformOnset(offset = 20, width = 50),
    NoNoise();
    return_epoched = true,
);
size(data)

# For Unfold.jl, we have to reshape the data, so that all subjects are concatenated.
data = reshape(data, size(data, 1), :)
times = range(0, 1, length = size(data, 1))
m = fit(
    UnfoldModel,
    @formula(0 ~ 1 + condition + (1 | item) + (1 + condition | subject)),
    evts,
    data,
    times,
)
plot_erp(coeftable(m))#, mapping = (; col = :group)) # FIXME facetting by col currently broken, waiting for new UnfoldMakie release!

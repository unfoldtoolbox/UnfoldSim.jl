# # Overview: Experimental design types

# The experimental design specifies the experimental conditions and other variables that are supposed to have an influence on the simulated data.
# Currently, there are three types of designs implemented: `SingleSubjectDesign`, `MultiSubjectDesign` and `RepeatDesign`.

# ### Setup
# ```@raw html
# <details>
# <summary>Click to expand</summary>
# ```
## Load required packages
using UnfoldSim
using Random
# ```@raw html
# </details >
# ```

# ## Single-subject designs
# As the name suggests, the `SingleSubjectDesign` type can be used to specify the experimental design for a single subject. Using the `conditions` arguments,
# the user can specify all relevant conditions or predictors and their levels or value range. 

# The current implementation assumes a full factorial design (also called fully crossed design)
# in which each level of a factor occurs with each level of the other factors. Moreover, in the current implementation, there is exactly one instance of each of these factor combinations. 

# Example:
design_single = SingleSubjectDesign(;
    conditions = Dict(
        :stimulus_type => ["natural", "artificial"],
        :contrast_level => range(0, 1, length = 3),
    ),
);

# In order to inspect the design, we can use the `generate_events` function to create an event table based on the design we specified.
generate_events(design_single)

# To change the order of the trials e.g. to sort or shuffle them, one can use the `event_order_function` argument.
# Example: Randomize the order of trials
design_single_shuffled = SingleSubjectDesign(;
    conditions = Dict(
        :stimulus_type => ["natural", "artificial"],
        :contrast_level => range(0, 1, length = 3),
    ),
    event_order_function = shuffle,
);
# ```@raw html
# <details>
# <summary>Click to expand event table </summary>
# ```
generate_events(design_single_shuffled)
# ```@raw html
# </details >
# ```

# ## Multi-subject designs
# The `MultiSubjectDesign` type can be used to simulate data for an experiment with multiple subjects. Internally, it uses the [MixedModelsSim.jl package](https://github.com/RePsychLing/MixedModelsSim.jl).
# One needs to specify the number of subjects `n_subjects` and the number of items `n_items` i.e. stimuli.
# In addition, one needs to decide for every experimental factor whether it should be between- or within-subject (and item).

# !!! note 
#     For factors that are not listed in `items_between` it is assumed that they vary within-item (accordingly for `subjects_between`).

design_multi = MultiSubjectDesign(
    n_subjects = 6,
    n_items = 4,
    items_between = Dict(:colour => ["red", "blue"]),
    subjects_between = Dict(:age_group => ["young", "old"]),
    both_within = Dict(:luminance => range(0, 1, length = 3)),
);

# ```@raw html
# <details>
# <summary>Click to expand event table </summary>
# ```
generate_events(design_multi)
# ```@raw html
# </details >
# <br />
# ```

# As with the `SingleSubjectDesign` one can use the `event_order_function` argument to determine the order of events/trials.

# !!! important
#     The number of subjects/items has to be a divisor of the number of factor level combinations, i.e. it is assumed that the design is balanced
#     which means that there is an equal number of observations for all possible factor level combinations.

# ## Repeat designs
# The `RepeatDesign` type is a functionality to encapsulate single- or multi-subject designs. It allows to repeat a generated event table multiple times.
# In other words, the `RepeatDesign` type allows to have multiple instances of the same item/subject/factor level combination.

# Example:
# Assume, we have the following single-subject design from above:
# ```@raw html
# <details>
# <summary>Click to expand event table </summary>
# ```
generate_events(design_single)
# ```@raw html
# </details >
# <br />
# ```


# But instead of having only one instance of the factor combinations e.g. `stimulus_type`: `natural` and `contrast_level`: `0`, we will repeat the design three times such that there are three occurrences of each combination.
design_repeated = RepeatDesign(design_single, 3);
generate_events(design_repeated)

# [Here](@ref howto_repeat_design) one can find another example of how to repeat design entries for multi-subject designs.

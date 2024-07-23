using UnfoldSim


# # [Repeating Design entries](@id howto_repeat_design)
# Sometimes we want to repeat a design, that is, have multiple trials with identical values, but it is not always straight forward to implement. 
# For instance, there is no way to easily modify `MultiSubjectDesign` to have multiple identical subject/item combinations,
# without doing awkward repetitions of condition-levels or something.

# If you struggle with this problem `RepeatDesign` is an easy tool for you:

designOnce = MultiSubjectDesign(;
    n_items = 2,
    n_subjects = 2,
    subjects_between = Dict(:cond => ["levelA", "levelB"]),
    items_between = Dict(:cond => ["levelA", "levelB"]),
);

design = RepeatDesign(designOnce, 4);
generate_events(design)

# As you can see, the design was simply repeated.

# !!! note
#       If you implemented your own `AbstractDesign`, you need to define the size function accordingly. E.g.:
#       `Base.size(design::RepeatDesign{SingleSubjectDesign}) = size(design.design).*design.repeat`

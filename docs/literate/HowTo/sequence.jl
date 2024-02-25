using UnfoldSim
using CairoMakie
using StableRNGs

# ## Stimulus - Response design

# let's say we want to simulate a stimulus response, followed by a button press response. 
# First we generate the minimal design of the experiment by specifying our conditins (a one-condition-two-levels design in our case)
design = SingleSubjectDesign(conditions = Dict(:condition => ["one", "two"]))
generate_events(design)
# next we use the `SequenceDesign` and nest our initial design in it. "SR_" is code for an "S" event and an "R" event - only single letter events are supported! The `_` is a signal for the Onset generator to generate a bigger pause - no overlap between adjacend `SR` pairs
design = SequenceDesign(design, "SRR_")
generate_events(design)
# The main thing that happened is that the design was repeated for every event (each 'letter') of the sequence, and an `eventtype` column was added.
# !!! hint
#     more advaned sequences exist, like "SR{1,3}", or "A[BC]" etc.

# Finally, let's repeat the design 2 times
design = RepeatDesign(design, 2)
generate_events(design)

# This results in 20 trials that nicely follow our sequence

# Next we have to specify for both events `S` and `R` what the responses should look like.

p1 = LinearModelComponent(;
    basis = p100(),
    formula = @formula(0 ~ 1 + condition),
    β = [1, 0.5],
);

n1 = LinearModelComponent(;
    basis = n170(),
    formula = @formula(0 ~ 1 + condition),
    β = [1, 0.5],
);
p3 = LinearModelComponent(;
    basis = UnfoldSim.hanning(Int(0.5 * 100)), # sfreq = 100 for the other bases
    formula = @formula(0 ~ 1 + condition),
    β = [1, 2],
);

components = Dict('S' => [p1, n1], 'R' => [p3])

data, evts = simulate(
    StableRNG(1),
    design,
    components,
    UniformOnset(offset = 10, width = 100),
    NoNoise(),
)
lines(data)
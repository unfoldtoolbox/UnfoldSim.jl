using Base: add_sum
using UnfoldSim
using CairoMakie
using StableRNGs

# ## Stimulus - Response design

# Let's say we want to simulate a stimulus response, followed by a button press response. 
# 
# First we generate the minimal design of the experiment by specifying our conditins (a one-condition-two-levels design in our case)
design = SingleSubjectDesign(conditions = Dict(:condition => ["one", "two"]))
generate_events(design)
# Next we use the `SequenceDesign` and nest our initial design in it. "`SR_`" is code for an "`S`" event and an "`R`" event - only single letter events are supported! The "`_`" is a signal for the Onset generator to generate a bigger pause - no overlap between adjacend "`SR`" pairs
design = SequenceDesign(design, "SR_")
generate_events(StableRNG(1), design)
# The main thing that happened is that the design was repeated for every event (each 'letter') of the sequence, and an `eventtype` column was added.
# !!! hint
#     more advaned sequences are possible as well, like "SR{1,3}", or "A[BC]". Infinite sequences are **not** possible like "AB*". 

# Finally, let's repeat the current design 4 times
design = RepeatDesign(design, 4)
generate_events(StableRNG(1), design)

# This results in 16 trials that nicely follow our sequence

# !!! hint
#     There is a difference between `SequenceDesign(RepeatDesign)` and `RepeatDesign(SequenceDesign)` for variable sequences e.g. "A[BC]", where in the former case,  one sequence is drawn e.g. "AC" and applied to all repeated rows, in the latter, one sequence for each repeat is drawn.


# Next we have to specify for both events `S` and `R` what the responses should look like.
p1 = LinearModelComponent(;
    basis = p100(),
    formula = @formula(0 ~ 1 + condition),
    β = [1, 0.5],
)

n1 = LinearModelComponent(;
    basis = n170(),
    formula = @formula(0 ~ 1 + condition),
    β = [1, 0.5],
)

p3 = LinearModelComponent(;
    basis = UnfoldSim.hanning(Int(0.5 * 100)), # sfreq = 100 for the other bases
    formula = @formula(0 ~ 1 + condition),
    β = [1, 0],
)

resp = LinearModelComponent(;
    basis = UnfoldSim.hanning(Int(0.5 * 100)), # sfreq = 100 for the other bases
    formula = @formula(0 ~ 1 + condition),
    β = [1, 2],
    offset = -10,
)
nothing ## hide


# We combine them into a dictionary with a sequence-`Char` as key and simulate
components = Dict('S' => [p1, n1, p3], 'R' => [resp])

data, evts = simulate(
    StableRNG(1),
    design,
    components,
    UniformOnset(offset = 40, width = 10),
    NoNoise(),
)
nothing ## hide

# Finally we can plot the results
lines(data)
vlines!(evts.latency[evts.event.=='S'], color = (:darkblue, 0.5))
vlines!(evts.latency[evts.event.=='R'], color = (:darkred, 0.5))
xlims!(0, 500)
current_figure()

# As visible, the `R` response always follows the `S` response. Due to the "`_`" we have large breaks between the individual sequences.
# # Simulate event-related potentials (ERPs)

# One subfield of EEG research focuses on so-called event-related potentials (ERPs) which are defined as brain responses time-locked to a certain event e.g. stimulus onset.
# The waveform of an ERP usually consists of multiple ERP components which denote the peaks and troughs of the waveform.

# ERP components are characterized (and named) by their timing relative to the event, their polarity (positive or negative) and their scalp topography. 
# For example, the N170 describes a negative deflection which occurrs roughly 170 ms after the onset of (certain) visual stimuli.
# Often, researchers are interested how a component (e.g. its amplitude or timing) changes depending on certain experimental factors. 
# For example, N170 has been shown to be related to face processing and its amplitude is modulated by whether the stimulus is a face or an object e.g. a car.
# ([Source](https://neuraldatascience.io/7-eeg/components.html))

# Here we will learn how to simulate a typical ERP complex with P100, N170, P300.

# ### Setup
# ```@raw html
# <details>
# <summary>Click to expand</summary>
# ```
## Load required packages
using UnfoldSim
using CairoMakie
using Random
using Unfold
using UnfoldMakie
# ```@raw html
# </details >
# ```

# ## Simulation
# Let's grab a `SingleSubjectDesign` and add a continuous predictor
design =
    SingleSubjectDesign(;
        conditions = Dict(
            :condition => ["car", "face"],
            :continuous => range(0, 5, length = 10),
        ),
    ) |> x -> RepeatDesign(x, 100);

# Let's make use of the prespecified basis functions, but use different formulas + parameters for each!

# **p100** is unaffected by our design and has amplitude of 5
p1 = LinearModelComponent(; basis = p100(), formula = @formula(0 ~ 1), β = [5]);

# **n170** has a condition effect, faces are more negative than cars
n1 = LinearModelComponent(;
    basis = n170(),
    formula = @formula(0 ~ 1 + condition),
    β = [5, 3],
);
# **p300** has a continuous effect, higher continuous values will result in larger P300's.
# We include both a linear and a quadratic effect of the continuous variable.
p3 = LinearModelComponent(;
    basis = p300(),
    formula = @formula(0 ~ 1 + continuous + continuous^2),
    β = [5, 1, 0.2],
);

# Now we can simply combine the components and simulate 
components = [p1, n1, p3]
data, evts = simulate(
    MersenneTwister(1),
    design,
    components,
    UniformOnset(; width = 0, offset = 1000),
    PinkNoise(),
);


# ## Analysis
# Let's check that everything worked out well, by using Unfold
m = fit(
    UnfoldModel,
    Dict(
        Any => (
            @formula(0 ~ 1 + condition + spl(continuous, 4)),
            firbasis(τ = [-0.1, 1], sfreq = 100, name = "basis"),
        ),
    ),
    evts,
    data,
);

# first the "pure" beta/linear regression parameters
plot_erp(
    coeftable(m);
    axis = (
        title = "Estimated regression parameters",
        xlabel = "Time [s]",
        ylabel = "Amplitude [μV]",
    ),
)

# and now beautifully visualized as marginal betas / predicted ERPs
f = plot_erp(
    effects(Dict(:condition => ["car", "face"], :continuous => 0:0.5:5), m);
    axis = (
        title = "Predicted event-related potential (ERP)",
        xlabel = "Time [s]",
        ylabel = "Amplitude [μV]",
    ),
    mapping = (:color => :continuous, linestyle = :condition, group = :continuous),
    legend = (; valign = :top, halign = :right, tellwidth = false),
    categorical_color = false,
);

## Workaround to separate legend and colorbar (will be fixed in a future UnfoldMakie version)
legend = f.content[2]
f[:, 1] = legend
current_figure()

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
using UnfoldSim # For simulation
using Random # For randomization
using StableRNGs # To get an RNG
using Unfold # For analysis
using CairoMakie # For plotting
using UnfoldMakie # For plotting
# ```@raw html
# </details >
# ```

# ## Specify the simulation "ingredients"
# 1\. We specify an **experimental design** with one subject in two experimental conditions including a continuous variable with 10 values.
# To mimic randomization in an experiment, we shuffle the trials using the `event_order_function` argument. To generate more trials we repeat the design 100 times which results in 2000 trials in total.
design =
    SingleSubjectDesign(;
        conditions = Dict(
            :condition => ["car", "face"],
            :continuous => range(0, 5, length = 10),
        ),
        event_order_function = x -> shuffle(StableRNG(1), x),
    ) |> x -> RepeatDesign(x, 100);

# The `generate_events` function can be used to create an events data frame from the specified experimental design.
events_df = generate_events(design);
first(events_df, 5)
# Above you can see the first five rows extracted from the events data frame representing the experimental design. Each row corresponds to one event.
# The columns *continuous* and *condition* display the levels of the predictor variables for the specific event.

# 2\. Next, we create a signal consisting of two different **components**.
# For the first component, we use the prespecified N170 base with an intercept of 5µV and a condition effect of 3&nbsp;µV for the “face/car” condition i.e. faces will have a more negative signal than cars.
# For the second component, we use the prespecified P300 base and include a linear and a quadratic effect of the continuous variable: the larger the value of the continuous variable, the larger the simulated potential.

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

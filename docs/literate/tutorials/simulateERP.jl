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

# ### 1. Experimental design
# We specify an experimental design with one subject in two experimental conditions including a continuous variable with 10 values.
# To mimic randomization in an experiment, we shuffle the trials using the `event_order_function` argument. To generate more trials we repeat the design 100 times which results in 2000 trials in total.
design =
    SingleSubjectDesign(;
        conditions = Dict(
            :condition => ["car", "face"],
            :continuous => range(0, 5, length = 10),
        ),
        event_order_function = shuffle,
    ) |> x -> RepeatDesign(x, 100);

# The `generate_events` function can be used to create an events data frame from the specified experimental design.
events_df = generate_events(StableRNG(1), design);
first(events_df, 5)
# Above you can see the first five rows extracted from the events data frame representing the experimental design. Each row corresponds to one event.
# The columns *continuous* and *condition* display the levels of the predictor variables for the specific event.

# ### 2. Event basis functions (Components)
# Next, we create a signal consisting of three different components.

# For the first component, we use the prespecified **P100** base which will be unaffected by our design and has an amplitude of 5 µV.
p1 = LinearModelComponent(; basis = p100(), formula = @formula(0 ~ 1), β = [5]);

# For the second component, we use the prespecified **N170** base with an intercept of 5 µV and a condition effect of 3 µV for the “face/car” condition i.e. faces will have a more negative signal than cars.
n1 = LinearModelComponent(;
    basis = n170(),
    formula = @formula(0 ~ 1 + condition),
    β = [5, 3],
);
# For the third component, we use the prespecified **P300** base and include a linear and a quadratic effect of the continuous variable: the larger the value of the continuous variable, the larger the simulated potential.
p3 = LinearModelComponent(;
    basis = p300(),
    formula = @formula(0 ~ 1 + continuous + continuous^2),
    β = [5, 1, 0.2],
);

# ### 3. Inter-onset distribution
# In the next step, we specify an inter-onset distribution, in this case, a uniform distribution with an inter-event distance of exactly 200 samples.
onset = UniformOnset(; width = 0, offset = 200);

# ### 4. Noise specification
# As the last ingredient, we specify the **noise**, in this case, Pink noise.
noise = PinkNoise(; noiselevel = 2);

# ## Simulate data

# Finally, we combine all the ingredients and simulate data. To make the simulation reproducible, one can specify a random generator.

## Combine the components in a vector
components = [p1, n1, p3]

## Simulate data
eeg_data, events_df = simulate(StableRNG(1), design, components, onset, noise);

# To inspect the simulated data we will visualize the first 1400 samples.

# ```@raw html
# <details>
# <summary>Click to show the code for the figure below</summary>
# ```
f = Figure(size = (1000, 400))
ax = Axis(
    f[1, 1],
    title = "Simulated EEG data",
    titlesize = 18,
    xlabel = "Time [samples]",
    ylabel = "Amplitude [µV]",
    xlabelsize = 16,
    ylabelsize = 16,
    xgridvisible = false,
    ygridvisible = false,
)

n_samples = 1400
lines!(eeg_data[1:n_samples]; color = "black")
v_lines = [
    vlines!(
        [r["latency"]];
        color = ["orange", "teal"][1+(r["condition"]=="car")],
        label = r["condition"],
    ) for r in
    filter(:latency => x -> x < n_samples, events_df)[:, ["latency", "condition"]] |>
    eachrow
]
xlims!(ax, 0, n_samples)
axislegend("Event onset"; unique = true);
# ```@raw html
# </details >
# ```
current_figure()

# The vertical lines denote the event onsets and their colour represents the respective condition i.e. car or face.

# ## Validate the simulation results
# To validate the simulation results, we use the [`Unfold.jl` package](https://github.com/unfoldtoolbox/Unfold.jl/) to fit an Unfold regression model to the simulated data and examine the estimated regression parameters and marginal effects. For the formula, we include a categorical predictor for *condition* and a non-linear predictor (based on splines) for *continuous*.
m = fit(
    UnfoldModel,
    [
        Any => (
            @formula(0 ~ 1 + condition + spl(continuous, 4)),
            firbasis(τ = [-0.1, 1], sfreq = 100, name = "basis"),
        ),
    ],
    events_df,
    eeg_data,
);

# To inspect the modelling results we will visualize the model coefficient estimates together with the estimated marginal effects i.e. the predicted ERPs.

# ```@raw html
# <details>
# <summary>Click to show the code for the figure below</summary>
# ```
## Create a data frame with the model coefficients and extract the coefficient names
coefs = coeftable(m)
coefnames = unique(coefs.coefname)

f2 = Figure(size = (1000, 400))
ga = f2[1, 1] = GridLayout()
gb = f2[1, 2] = GridLayout()

## Plot A: Estimated regression parameters
ax_A = Axis(
    ga[1, 1],
    title = "Estimated regression parameters",
    titlegap = 12,
    xlabel = "Time [s]",
    ylabel = "Amplitude [μV]",
    xlabelsize = 16,
    ylabelsize = 16,
    xgridvisible = false,
    ygridvisible = false,
)

for coef in coefnames
    estimate = filter(:coefname => ==(coef), coefs)

    lines!(ax_A, estimate.time, estimate.estimate, label = coef)
end
axislegend("Coefficient", framevisible = false)
hidespines!(ax_A, :t, :r)

## Plot B: Marginal effects
plot_B = plot_erp!(
    gb,
    effects(Dict(:condition => ["car", "face"], :continuous => 0:0.5:5), m);
    mapping = (; color = :continuous, linestyle = :condition, group = :continuous),
    legend = (; valign = :top, halign = :right),
    axis = (
        title = "Marginal effects",
        titlegap = 12,
        xlabel = "Time [s]",
        ylabel = "Amplitude [μV]",
        xlabelsize = 16,
        ylabelsize = 16,
        xgridvisible = false,
        ygridvisible = false,
    ),
)

## Add letter labels to the plots
## Adapted from: https://docs.makie.org/stable/tutorials/layout-tutorial/
for (label, layout) in zip(["A", "B"], [ga, gb])
    Label(
        layout[1, 1, TopLeft()],
        label,
        fontsize = 26,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right,
    )
end
# ```@raw html
# </details >
# ```

current_figure()

# In subplot A, one can see the model coefficient estimates and as intended the first component is unaffected by the experimental design,
# there is a condition effect in the second component and an effect of the continuous variable on the third component.
# The relation between the levels of the continuous variable and the scaling of the third component is even clearer visible in subplot B which depicts the estimated marginal effects of the predictors
# which are obtained by evaluating the estimated function at specific values of the continuous variable. 

# As shown in this example, `UnfoldSim.jl` and `Unfold.jl` can be easily combined to investigate the effects of certain features,
# e.g. the type of noise or its intensity on the analysis result and thereby assess the robustness of the analysis.

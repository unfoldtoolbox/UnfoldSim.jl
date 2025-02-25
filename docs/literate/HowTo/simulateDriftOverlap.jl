# # Simulate an Evidence Accumulation Overlap and Deconvolution

# We want to use the SequenceDesign with three components to simulate an full evidence accumulation process with a Overlap between the evidence accumulation and the response signal. Afterwards we use the Deconvolution Method from Unfold to unfold the signals again.



# ### Setup
# ```@raw html
# <details>
# <summary>Click to expand</summary>
# ```
using UnfoldSim
using Unfold
using StableRNGs
using CairoMakie, UnfoldMakie
using SequentialSamplingModels

fs = 500
Δt = 1 / fs; # time step
tEnd = 1.0 # trial Duration
time_vec = 0:Δt:tEnd; # time base - let's make it a typical stimulus duration

p3 = LinearModelComponent(;
    basis = UnfoldSim.hanning(Int(0.7 * fs)),
    formula = @formula(0 ~ 1 + condition),
    β = [0.8, 0],
)
resp = LinearModelComponent(;
    basis = UnfoldSim.hanning(Int(0.7 * fs)),
    formula = @formula(0 ~ 1 + condition),
    β = [0.5, 0],
)

# ```@raw html
# </details >
# ```
# ## Design
# Let's generate a single design with two different drift_rates as condition, so we can use them later for the Simulation
design_single =
    SingleSubjectDesign(conditions = Dict(:drift_rate => [3, 5], :condition => [1]))
design_seq = SequenceDesign(design_single, "SCR_")
design_rep = RepeatDesign(design_seq, 500)

# ## Create the Drift Component using the special KellyModel for the evidence accumulation
# As parameter for the models drift_rate we use the drift_rate specified in the design by directly naming the condition as string.
v = "drift_rate" # Drift Rate from the design
# Now we retrieve the models default parameters, but we change the drift_rate and
kelly_model = KellyModel(drift_rate = v, motor_onset = 0.4, event_onset = 0.2)
kelly_model_parameter = UnfoldSim.create_kelly_parameters_dict(kelly_model)
drift = Drift_Component(simulate_component, time_vec, Δt, KellyModel, kelly_model_parameter)
components = Dict('S' => [p3], 'C' => [drift], 'R' => [resp])
# ## Create the Onset Component to simulate the overlap
# To simulate an overlap between the drift_component in 'C' and the response component in 'R'. We have to specify the UniformOnset of the 'C' Component therefore with an negative offset to produce the overlap.
seq_onset = SequenceOnset(
    Dict(
        'S' => UniformOnset(width = 0, offset = 1.2 * fs),
        'C' => (DriftOnset(), UniformOnset(width = 0, offset = -100)),
        'R' => UniformOnset(width = 0, offset = 2 * fs),
    ),
)

# ## Simulate data with the created design, components and onsets
data, evts = UnfoldSim.simulate(
    StableRNG(12),
    design_rep,
    components,
    seq_onset,
    NoNoise(), # PinkNoise(noiselevel=1)
)
# ## Plot EEG to see the overlap as it appears
lines(data)
vlines!(evts.latency[evts.event.=='S'], color = (:green, 0.5))
vlines!(evts.latency[evts.event.=='C'], color = (:blue, 0.5))
vlines!(evts.latency[evts.event.=='R'], color = (:orange, 0.5))
CairoMakie.xlims!(0, 2000)
current_figure()

# ## Plot ERP with the overlap in the 'R' Component
evts.event = string.(evts.event)
data_epochs, times_epoch =
    Unfold.epoch(data = data, tbl = evts, τ = (-0.1, 0.9), sfreq = fs);
f = @formula(0 ~ 1)
m = fit(
    UnfoldModel,
    ["S" => (f, times_epoch), "R" => (f, times_epoch), "C" => (f, times_epoch)],
    evts,
    data_epochs,
);
results = coeftable(m)
plot_erp(results; mapping = (; color = :eventname))

# ## Deconvolution of the overlap
fir = firbasis(τ = (-0.1, 0.9), sfreq = fs)
evts2 = deepcopy(evts)
evts2.event = string.(evts2.event)
m = fit(
    UnfoldModel,
    ["S" => (f, fir), "R" => (f, deepcopy(fir)), "C" => (f, deepcopy(fir))],
    evts2,
    data,
);
results = coeftable(m)
plot_erp(results; mapping = (; color = :eventname))
# # Simulate an Evidence Accumulation EEG

# We want to use the SequenceDesign with a stimulus component, a evidence accumulation component (DriftComponent) and a response component to simulate an full evidence accumulation process.



# ### Setup
# ```@raw html
# <details>
# <summary>Click to expand</summary>
# ```
using UnfoldSim
using Unfold
using Random
using CairoMakie, UnfoldMakie

fs = 500
Δt = 1/fs; # time step
tEnd = 1.0 # trial Duration
time_vec = 0:Δt:tEnd; # time base - let's make it a typical stimulus duration
# ```@raw html
# </details >
# ```

# ## Design
# Let's generate a single design with two different drift_rates as condition
design_single = SingleSubjectDesign(
    conditions = Dict(
        :drift_rate => [0.5, 0.8],
        :condition => [1])
    )
# Now we convert the SingleSubjectDesign to an SequenceDesign with a Sequence of S: stimulus, C: component, R: response, _: gap
design_seq = SequenceDesign(
    design_single,
    "SCR_"
    )
# On top we can now repeat the Sequence to have multiple trials. For this experiment we use 100 replicas.
design_rep = RepeatDesign(
    design_seq,
    100
    )

# ## Create the 3 components we want to use
# First the stimulus component as simple p3
p3 = LinearModelComponent(;
    basis = UnfoldSim.hanning(Int(0.5 * fs)),
    formula = @formula(0 ~ 1 + condition),
    β = [1.0, 0],
)
# Second a response component
resp = LinearModelComponent(;
    basis = UnfoldSim.hanning(Int(0.5 * fs)),
    formula = @formula(0 ~ 1 + condition),
    β = [0.5, 0],
)
# Third we create our Drift_Component which implies the evidence accumulation. For the Drift_Component we can choose between different Models which are used for the Simulation. As parameter for the models drift_rate we use the drift_rate specified in the design by directly naming the condition.
v = "drift_rate" # Drift Rate
a = 2.0 # The amount of information that is considered for a decision.
t = 0.10 # The duration for a non-decisional processes (encoding and response execution).
z = 0.50 # An indicator of an an initial bias towards a decision.
ddm_parameter = Dict(:ν=>v, :α=>a, :z=>z, :τ=>t)
drift = Drift_Component(
    simulate_component,
    time_vec,
    Δt,
    DDM,
    ddm_parameter)
# As last step we have to specify the components as a Dict connection the components with the events of the design.
components = Dict('S' => [p3], 'C' => [drift], 'R' => [resp])

# ## Create the Onset Component for the design
# For the stimulus and response we use a simple Uniform onset with a distance of 250 and 300. The Drift_Component uses the special DriftOnset in combination with an UniformOnset. The DriftOnset returns the response time after the model reaches the decision threshold.
# Important is to know that, the onset for each component defines the onset for the next. So in this case: S->C, C->R, R->S.
seq_onset = SequenceOnset(
    Dict('S'=>UniformOnset(width=0,offset=250),
         'C'=>(DriftOnset(), UniformOnset(width=0, offset=150)),
         'R'=>UniformOnset(width=0,offset=300)))

# ## Simulate data with the created design, components and onsets
data, evts = UnfoldSim.simulate(
    StableRNG(12),
    design_rep,
    components,
    seq_onset,
    NoNoise() # PinkNoise(noiselevel=1)
)
lines(data)
vlines!(evts.latency[evts.event.=='S'], color = (:green, 0.5))
vlines!(evts.latency[evts.event.=='C'], color = (:blue, 0.5))
vlines!(evts.latency[evts.event.=='R'], color = (:orange, 0.5))
CairoMakie.xlims!(0, 2000)
current_figure()

# ## Plot ERP of simulated EEG
evts.event = string.(evts.event)
data_epochs, times_epoch = Unfold.epoch(data = data, tbl = evts, τ = (0, 1.0), sfreq = fs);
f = @formula(0 ~ 1)
m = fit(UnfoldModel, ["S"=>(f,times_epoch),"R"=>(f,times_epoch),"C"=>(f,times_epoch)], evts, data_epochs);
results = coeftable(m)
plot_erp(results;mapping=(;color=:eventname))
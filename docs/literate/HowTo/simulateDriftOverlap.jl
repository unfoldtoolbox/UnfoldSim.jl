# # Simulate an Evidence Accumulation Overlap and Deconvolution

# We want to use the SequenceDesign with three components to simulate an full evidence accumulation process with a Overlap between the evidence accumulation and the response signal. Afterwards we use the Deconvolution Method from Unfold to unfold the signals again.



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

design_single = SingleSubjectDesign()
design_seq = SequenceDesign(design_single, "SCR_")
design_rep = RepeatDesign(design_seq, 500)
p3 = LinearModelComponent(;
    basis = UnfoldSim.hanning(Int(0.5 * fs)),
    formula = @formula(0 ~ 1 + condition),
    β = [1.0, 0],
)
resp = LinearModelComponent(;
    basis = UnfoldSim.hanning(Int(0.5 * fs)),
    formula = @formula(0 ~ 1 + condition),
    β = [0.5, 0],
)
v = 0.6 # Drift Rate
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
components = Dict('S' => [p3], 'C' => [drift], 'R' => [resp])

# ```@raw html
# </details >
# ```

# ## Create the Onset Component to simulate the overlap
# To simulate an overlap between the drift_component in 'C' and the response component in 'R'. We have to specify the UniformOnset of the 'C' Component therefore with an negative offset to produce the overlap.
seq_onset = SequenceOnset(
    Dict('S'=>UniformOnset(width=0,offset=250),
         'C'=>(DriftOnset(), UniformOnset(width=0, offset=-150)),
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

# ## Plot ERP with the overlap in the 'R' Component
evts.event = string.(evts.event)
data_epochs, times_epoch = Unfold.epoch(data = data, tbl = evts, τ = (0, 1.0), sfreq = fs);
f = @formula(0 ~ 1)
m = fit(UnfoldModel, ["S"=>(f,times_epoch),"R"=>(f,times_epoch),"C"=>(f,times_epoch)], evts, data_epochs);
results = coeftable(m)
plot_erp(results;mapping=(;color=:eventname))

# ## Deconvolution of the overlap
fir = firbasis(τ=(0,1.0),sfreq=fs)
evts2 = deepcopy(evts)
evts2.event = string.(evts2.event)
m = fit(UnfoldModel, ["S"=>(f,fir),"R"=>(f,deepcopy(fir)),"C"=>(f,deepcopy(fir))], evts2, data);
results = coeftable(m)
plot_erp(results;mapping=(;color=:eventname))
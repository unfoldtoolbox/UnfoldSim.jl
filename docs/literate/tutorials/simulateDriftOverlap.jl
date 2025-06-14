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
max_length = tEnd / Δt

c_stimulus = LinearModelComponent(;
    basis = ones(Int(0.1 * fs)),
    formula = @formula(0 ~ 1 + condition),
    β = [1, 0],
)
c_response = LinearModelComponent(;
    basis = ones(Int(0.2 * fs)),
    formula = @formula(0 ~ 1 + condition),
    β = [-0.5, 0],
)

# ```@raw html
# </details >
# ```
# ## Design
# Let's generate a single design with two drift_rates
design_single =
    SingleSubjectDesign(conditions = Dict(:drift_rate => [1, 2], :condition => [1]))
design_seq = SequenceDesign(design_single, "SCR_")
design_rep = RepeatDesign(design_seq, 100) # 100 SCR-repeats

# # Simplified S-C-R drift-model with LBA
# ## LBA
# The following parameters are described in the [SequentialSamplingModels.LBA model documentation](https://itsdfish.github.io/SequentialSamplingModels.jl/dev/lba/)
v = [0.8] ## Drift Rate, needs to be in [ ] for compatability with `SequentialSamplingModels.jl``
A = 0.2 ## The starting point is sampled uniformly between `[0, A]`
k = 0.4 ## Difference of A with the bound/threshold ($A + k = b$), where b is the decision threshold.

lba_parameter = Dict(:ν => v, :A => A, :k => k, :τ => 0.0) # we didnt specify σ thus it is 1 by default,

# !!! important
# The non-decision time has to be specified here as "τ=>0." because it is not defined whether non-decision time applies to stimulus encoding or response execution. If you have a non-0 value here, the events following the DriftComponent will be shifted by τ.
drift = DriftComponent(max_length, fs, LBA, lba_parameter)

# Let's have a quick look at the simulation
series(
    hcat(simulate_component.(StableRNG.(1:10), Ref(drift), Ref(design_single))...)',
    solid_color = (:black, 0.5),
)
# Super! Our start values are sampled between [0, 0.2] as indicated, and the boundary is at A+k = 0.2 + 0.4 = 0.6

# ### Simulate EEG data 
components = Dict('S' => [c_stimulus], 'C' => [drift], 'R' => [c_response])

# Next we define a SequenceOnset ingredient, this allows us fine-grained control on when events should happen.

seq_onset = SequenceOnset(
    Dict(
        'S' => UniformOnset(width = 0, offset = 0.2 * fs), # distance S<->C - we could choose 0.1*fs to immediately start after the stimulus component
        'C' => DriftOnset(), ## automatically determines which is the respective DriftComponent
        'R' => UniformOnset(width = 0, offset = 2 * fs), # distance R<->S
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
f_data, _ax, _h = lines(data[1:5000])
vlines!(evts.latency, color = Int.(evts.event), alpha = 0.3)
xlims!(700, 3500)
f_data

# # Complex  S-C-R drift-model ala `Kelly et al.`

# In addition to changing the model to a more realistic version, we will also make use of the `drift_rate` specification in our `design``

v = "drift_rate" ## take the drift rate directly from the design

# The KellyModel has _many_ default parameters, let's grab them and immediately change a few

kelly_model =
    UnfoldSim.KellyModel(drift_rate = v, motor_delay = 0.4, sensor_encoding_delay = 0.2)



% parameters:
% Urgency:
Z = 0.3; % starting level at evidence onset 
U = 1; % slope (buildup rate) of dynamic urgency at the motor level (prop. to bound /sec)
Sz = 0.4; % start point variability at the motor level applied independently to both sides - full width of uniform distribution
Su = 0.17 % between-trial var in urgency rate

d = [1 2]; % drift rate 
s = 0.5; % gaussian noise /sec at accumulator level
evonT = 0.2; % time at which sensory evidence starts impacting on the sensory accumulator
%accT = 0.2; %  Accumulation onset time. Let's say same as evidence encoding onset.
mT = 0.1; %  motor time
st = 0.1; %  motor time variability
Terz = 0.1; % variability applied to sensory encoding delays

postAccDurM = 0.1; % mean duration of post decision accumulation
postAccDurS = 0.2; % range of uniform dist of post-decision accumulation duration

CPPrampdownDur = 0.1; % in Sec; CPP will linearly ramp back down after the accumulation stops


kelly_model = UnfoldSim.KellyModel(
    drift_rate = v,
    motor_delay = 0.1,
    motor_delay_varibility = 0.1,
    sensor_encoding_delay_variability = 0.1,
    sensor_encoding_delay = 0.1,
    accumulative_level_noise=0.5,
    urgency = 1,
    urgency_variability = 0.17,
    post_accumulation_duration = 0.1,
    post_accumulation_duration_variability = 0.2,
    boundary = 0.7,
    ramp_down_duration = 0.1
)
kelly_model_parameter = UnfoldSim.create_kelly_parameters_dict(kelly_model)
drift = DriftComponent(max_length, fs, UnfoldSim.KellyModel, kelly_model_parameter)
components = Dict('S' => [c_stimulus], 'C' => [drift], 'R' => [c_response])

# We are ready to simulate!
data, evts = UnfoldSim.simulate(
    StableRNG(12),
    design_rep,
    components,
    seq_onset,
    NoNoise(), # PinkNoise(noiselevel=1)
)
# ## Plot EEG to see the overlap as it appears
f_data, _ax, _h = lines(data[1:5000])
vlines!(evts.latency, color = Int.(evts.event), alpha = 0.3)
xlims!(700, 3500)
f_data

# !!! note
# There are varying gaps between the stimulus offset <-> driftstart and driftend <-> response start. This is because the `KellyModel` by default has a `sensor_encoding_delay` and a  `motor_delay` (including a `motor_delay_variability` etc.)

# # ERP & Deconvolution 
# Let's calculate traditional ERPs and deconvolved ones!
evts.event = string.(evts.event) # due to a Unfold-bug we have to convert :S
data_epochs, times_epoch =
    Unfold.epoch(data = data, tbl = evts, τ = (-1.3, 1.3), sfreq = fs);
f = @formula(0 ~ 1 + drift_rate)
m = fit(UnfoldModel, ["S" => (f, times_epoch), "R" => (f, times_epoch)], evts, data_epochs);
plot_erp(
    effects(Dict(:drift_rate => [2, 3]), m);
    mapping = (;
        color = :drift_rate => nonnumeric,
        col = :eventname => UnfoldMakie.sorter(["S", "R"]),
    ),
)

# ## Deconvolution of the overlap
fir_stim = firbasis(τ = (-0.1, 1.3), sfreq = fs)
fir_resp = firbasis(τ = (-1, 0.3), sfreq = fs)

m = fit(UnfoldModel, ["S" => (f, fir_stim), "R" => (f, fir_resp)], evts, data);
plot_erp(
    effects(Dict(:drift_rate => [2, 3]), m);
    mapping = (;
        color = :drift_rate => nonnumeric,
        col = :eventname => UnfoldMakie.sorter(["S", "R"]),
    ),
)
# Deconvolution sucessfully removed the stimulus and response overlaps from our estimates, but it kept the drift-diffusion aspect nearly untouched. It is important to stress the "nearly", because depending on the `delay` parameters, this can potentially lead to problematic situation where it is unclear whether the drift activity should be assigned to S or R. In most realistic situations we found to be safe.
# Further note that the drift-diffusion amplitude can be reduced because it is split up between `S` and `R` events.`
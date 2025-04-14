"""
    KellyModel

A advanced drift diffusion Model which can be used to simulate evidence accumulation.

All fields can be named. Is used with [`DriftComponent`](@ref).

# Fields T::Union{Real,String}
- `drift_rate::T`: defines the amount of evidence accumulated per time step. (roughly the steepness of the trace)
- `sensor_encoding_delay::T`: constant event onset delay in seconds. (mimics sensory evidence)
- `sensor_encoding_delay_variability::T`: Normal variability in the delay of the event onset in seconds, added ontop of `event_onset`. (mimics sensory encoding delay)
- `event_onset_distribution::Normal`: Normal distribution for sensor encoding delay
- `accumulative_level_noise::T`: sigma of the normal-distributed noise added to the accumulation process.
- `boundary::T`: the threshold of evidence needed to make a decision.
- `motor_delay::T`: fixed delay between boundary reached and response time in seconds. (mimics motor time)
- `motor_delay_variability::T`: variability in delay between boundary reached and response time in seconds. (mimics different reaction times of participants)
- `motor_onset_distribution::Normal`: Normal distribution for motor delay
- `post_accumulation_duration::T`: fixed time the accumulation process resumes after boundary reached in seconds. (mimics evidence overshoot)
- `post_accumulation_duration_variability::T`: variability in time the accumulation process resumes after boundary reached in seconds.
- `post_accumulation_distribution::Normal`: Normal distribution for post accumulation duration
- `ramp_down_duration::T`: duration (in s) of the post accumulation ramp down process.

# Examples
```julia-repl
julia> KellyModel();
KellyModel{Float64}(6.0, 0.2, 0.4, 0.5, 1.0, 0.1, 0.4, 0.1, 0.2, 0.1)
```

See also [`LinearModelComponent`](@ref), [`MultichannelComponent`](@ref).
"""
Base.@kwdef mutable struct KellyModel <: SequentialSamplingModels.SSM2D
    drift_rate::Union{Real,String} = 6.0                    # drift rate
    sensor_encoding_delay::Union{Real,String} = 0.2                   # onset (sensory evidence)
    sensor_encoding_delay_variability::Union{Real,String} = 0.1         # sensory encoding delay
    event_onset_distribution::Normal = Normal(sensor_encoding_delay, sensor_encoding_delay_variability) # Normal distribution for sensor encoding delay
    accumulative_level_noise::Union{Real,String} = 0.5      # accumulation level noise
    boundary::Union{Real,String} = 1.0                      # boundary height
    motor_delay::Union{Real,String} = 0.4                   # motor onset
    motor_delay_variability::Union{Real,String} = 0.1                   # motor delay
    motor_onset_distribution::Normal = Normal(motor_delay, motor_delay_variability) # Normal distribution for motor delay
    post_accumulation_duration::Union{Real,String} = 0.1    # mean post-decision duration
    post_accumulation_duration_variability::Union{Real,String} = 0.001  # variability post-decision
    post_accumulation_distribution::Normal = Normal(post_accumulation_duration, post_accumulation_duration_variability) # Normal distribution for post accumulation duration
    ramp_down_duration::Union{Real,String} = 0.1            # CPPrampdown duration
end

"""
    create_kelly_parameters_dict(model::KellyModel)

Convert a `KellyModel` instance into a dictionary containing its parameters.

# Arguments
- `model::KellyModel`: The Kelly model instance whose parameters will be extracted.

# Returns
- `Dict{Symbol, Any}`: A dictionary where the keys are the parameter names as symbols, and the values are the corresponding parameter values.

# Examples
```julia-repl
julia> model = KellyModel(drift_rate=5.5, boundary=1.2)
KellyModel(5.5, 0.2, 0.4, 0.5, 1.2, 0.1, 0.4, 0.1, 0.2, 0.1)

julia> create_kelly_parameters_dict(model)
Dict{Symbol, Any}(:drift_rate => 5.5, :event_onset => 0.2, :sensor_encoding_delay => 0.4, 
                  :accumulative_level_noise => 0.5, :boundary => 1.2, :motor_onset => 0.1, 
                  :motor_delay => 0.4, :post_accumulation_duration => 0.1, 
                  :post_accumulation_duration_variability => 0.2, :ramp_down_duration => 0.1)
"""
function create_kelly_parameters_dict(model::KellyModel)
    return Dict(name => getfield(model, name) for name in fieldnames(typeof(model)))
end


"""
    KellyModel_simulate_cpp(rng, model::KellyModel, time_vec, Δt)

Generate a single response time, and evidence trace of an evidence accumulation process using the Kelly model.

# Arguments
- `rng::StableRNG`: Random seed to ensure the same traces are created for reconstruction.
- `model::KellyModel`:  specifies the model and its parameters to simulate the evidence accumulation.
- `time_vec::StepRangeLen`: range of time steps for which the evidence is accumulated.
- `Δt::Float64`: size of the time steps

# Returns
- `Float64`: Simulated response time for this trial.
- `Vector{Float64}`: evidence values over time. The output dimension is `length(time_vec)`.

# Examples
```julia-repl
# use the KellyModel and its default parameters to simulate traces from 0:1/500:1.0
julia> KellyModel_simulate_cpp(StableRNG(1), KellyModel(), 0:1/500:1.0, 1/500)
(260.70134768436486, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
```
"""
function KellyModel_simulate_cpp(rng, model::KellyModel, time_vec, Δt)
    evidence = zeros(length(time_vec));
    evidence[time_vec .>= rand(rng,model.event_onset_distribution)] .= 1; 
    startAccT = time_vec[findfirst(evidence .== 1)];

    noise = vcat(zeros(
        sum(time_vec .< startAccT)),
        randn(rng,sum(time_vec .>=  startAccT)) .*  model.accumulative_level_noise .*sqrt(Δt));

    cum_evidence = cumsum(evidence .* model.drift_rate .* Δt .+ noise); # This is the cumulative differential evidence, just as in a 1d DDM. 
    
    
    # terminate the decision process on boundary crossing, record threshold-crossing samplepoint:
    cum_evidence = abs.(cum_evidence);
    dti = findfirst(cum_evidence .> model.boundary); # finding the sample point of threshold crossing of each, then will pick the earlier as the winner
    if isnothing(dti)  # Check if no crossing was found
        dti = length(time_vec)  # Set to the last time step
    end
    # now record RT in sec after adding motor time, with variability
    rt = time_vec[dti] + rand(rng,model.motor_onset_distribution)

    # now make the CPP peak and go down linearly after a certain amount of post-dec accum time for this trial:
    post_acc_duration = rand(rng,model.post_accumulation_distribution)
    # so post_acc_duration is the post accumulation duration time, where the accumulation spikes over the threshold

    # acc_stop_index is the accumulation Stop index which is the index from the time Vector where the accumulation really stops
    acc_stop_index = dti + (post_acc_duration÷Δt)|>Int;
    # Take the absolute value of the accumulations
    cum_evidence = abs.(cum_evidence)
    if acc_stop_index < length(time_vec)
        nT = length(time_vec)
        tmp = cum_evidence[acc_stop_index] .- (1:(nT-acc_stop_index)) .*cum_evidence[acc_stop_index] .* (Δt ./ model.ramp_down_duration)
        cum_evidence[(acc_stop_index+1):end] .= max.(Ref(0), tmp);
    end
    return rt / Δt, cum_evidence[1:end]
end


"""
    simulate_drift_component(rng, component::DriftComponent, design::AbstractDesign)

Generate response times and evidence Vectors of an given [`AbstractDesign`](@ref) with a [`DriftComponent`](@ref) which contains the model used for the simulation.

# Arguments
- `rng::StableRNG`: Random seed to ensure the same traces are created for reconstruction.
- `component::DriftComponent`: Component to specify the model and its parameters to simulate the evidence accumulation.
- `design::AbstractDesign`: design of the experiment preferable SequenceDesign.

# Returns
- `Vector{Float64}`: Simulated response times for the trials.
- `Matrix{Float64}`: evidence values over time for each trial. The output dimensions are `c.max_length x size(events, 1)`.

# Examples
```julia-repl
julia> model_parameter = Dict(:motor_onset => 0.4, :event_onset => 0.2);

julia> c = DriftComponent(500, 500, KellyModel, model_parameter);

julia> design_single = SingleSubjectDesign(conditions = Dict(:drift_rate => [0.5, 0.8], :condition => [1]));

julia> design_seq = SequenceDesign(design_single,"SCR_");

julia> simulate_component(StableRNG(1),c,design_seq)
Vector{Float64}, 501x6 Matrix{Float64}:
([96.65745162948949, 273.7368235451535, 271.86040880709123, 128.41057786118193, 342.35208862144276, 237.14773586760617], [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0])
```
"""
function simulate_drift_component(rng, component::DriftComponent, design::AbstractDesign)
    events = generate_events(deepcopy(rng), design)
    traces = Matrix{Float64}(undef, component.max_length, size(events, 1))
    rts = Vector{Float64}(undef, size(events, 1))
    for (i, evt) in enumerate(eachrow(events))
        parameters = get_model_parameter(rng, evt, component.model_parameters)
        model = component.model_type(; (key => parameters[key] for key in keys(parameters))...)
        rt, evidence = SSM_Simulate(rng, model, component.sfreq, component.max_length)

        rts[i] = rt
        traces[:, i] = evidence[1:component.max_length]
    end
    return rts, traces
end

"""
    SSM_Simulate(rng, model::SequentialSamplingModels.SSM2D, sfreq, max_length)

Generate response time and evidence Vector of component.max_length by using a SequentialSamplingModels.SSM2D as DDM or LBA for the simulation.

# Arguments
- `rng::StableRNG`: Random seed to ensure the same traces are created for reconstruction.
- `model::SequentialSamplingModels.SSM2D`: SequentialSamplingModel to simulate the evidence and response time.
- `sfreq::Real`: sample frequency used to simulate the signal.
- `max_length::Int`: maximum length of the simulated signal.

# Returns
- `Float64`: Simulated response time for the trial.
- `Vector{Float64}`: evidence values over time. The output dimension is `component.max_length`.

# Examples
```julia-repl
julia> model = DDM()
julia> SSM_Simulate(StableRNG(1), model, 500, 500)
Float64, Vector{Float64}:
(96.65745162948949, [0.0 0.0 … 0.0 0.0])
```
"""
function SSM_Simulate(rng, model::SequentialSamplingModels.SSM2D, sfreq, max_length)
    Δt = 1 / sfreq
    time_steps, evidence = SequentialSamplingModels.simulate(rng, model; Δt)
    if !(evidence isa Vector{Float64})
        evidence = hcat(evidence...)
        evidence = collect(vec(evidence))
    end
    # Store results for this trial
    rt = (time_steps[end] + model.τ) / Δt
    if length(evidence) < max_length
        final_value = 0
        append!(evidence, fill(final_value, max_length - length(evidence)))
    else
        rt = max_length + (model.τ / Δt)
        evidence = evidence[1:max_length]
    end
    return rt, evidence
end

"""
    SSM_Simulate(rng, model::KellyModel, sfreq, max_length)

Generate response time and evidence Vector of max_length by using the Kelly Model for the simulation.

# Arguments
- `rng::StableRNG`: Random seed to ensure the same traces are created for reconstruction.
- `model::KellyModel`: SequentialSamplingModel to simulate the evidence and response time.
- `sfreq::Real`: sample frequency used to simulate the trace.
- `max_length::Int`: maximum length of the simulated trace.

# Returns
- `Float64`: Simulated response time for the trial.
- `Vector{Float64}`: evidence values over time. The output dimension is `c.max_length`.

# Examples
```julia-repl
julia> model = KellyModel()

julia> SSM_Simulate(StableRNG(1), model, 500, 500)
Float64, Vector{Float64}:
(96.65745162948949, [0.0 0.0 … 0.0 0.0])
```
"""
function SSM_Simulate(rng, model::KellyModel, sfreq, max_length)
    Δt = 1/sfreq
    time_vec = 0:Δt:max_length*Δt
    rt, evidence = KellyModel_simulate_cpp(rng, model, time_vec, Δt)
    if length(evidence) < max_length
        final_value = 0
        append!(evidence, fill(final_value, max_length - length(evidence)))
    end
    evidence = evidence[1:max_length]
    return rt, evidence
end
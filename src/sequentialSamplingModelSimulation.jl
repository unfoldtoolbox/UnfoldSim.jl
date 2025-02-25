"""
    KellyModel

A advanced drift diffusion Model which can be used to simulate evidence accumulation.

All fields can be named. Is used with [`DriftComponent`](@ref).

# Fields T::Real
- `drift_rate::T`: defines the amount of evidence accumulated per time step. (steepness of the trace)
- `event_onset::T`: event onset delay in seconds. (mimics sensory evidence)
- `sensor_encoding_delay::T`: variability in the delay of the event onset in seconds. (mimics sensory encoding delay)
- `accumulative_level_noise::T`: some noise added to the accumulation process.
- `boundary::T`: the threshold of evidence needed to make a decision.
- `motor_onset::T`: fixed delay between boundary reached and response time in seconds. (mimics motor time)
- `motor_delay::T`: variability in delay between boundary reached and response time in seconds. (mimics different reaction times of participants)
- `post_accumulation_duration_mean::T`: fixed time the accumulation process resumes after boundary reached in seconds. (mimics evidence overshoot)
- `post_accumulation_duration_variability::T`: variability in time the accumulation process resumes after boundary reached in seconds. (mimics diff of participants)
- `CPPrampdownDur::T`: duration post accumulation process ramp down.

# Examples
```julia-repl
julia> KellyModel();
KellyModel{Float64}(6.0, 0.2, 0.4, 0.5, 1.0, 0.1, 0.4, 0.1, 0.2, 0.1)
```

See also [`LinearModelComponent`](@ref), [`MultichannelComponent`](@ref).
"""
mutable struct KellyModel
    drift_rate::Union{Real,String} # drift rate
    event_onset::Union{Real,String} # onset(sensory evidence)
    sensor_encoding_delay::Union{Real,String} # var(sensory encoding delay)
    accumulative_level_noise::Union{Real,String} # accum level noise
    boundary::Union{Real,String} # boundaryary height
    motor_onset::Union{Real,String} # onset(motor)
    motor_delay::Union{Real,String} # var(motor)
    post_accumulation_duration_mean::Union{Real,String} # mean(post decision)
    post_accumulation_duration_variability::Union{Real,String} # var(post decision)
    CPPrampdownDur::Union{Real,String} # CPPrampdown duration

    # Constructor with default values
    function KellyModel(;
        drift_rate = 6.0,
        event_onset = 0.2,
        sensor_encoding_delay = 0.1,
        accumulative_level_noise = 0.5,
        boundary = 1.0,
        motor_onset = 0.4,
        motor_delay = 0.1,
        post_accumulation_duration_mean = 0.1,
        post_accumulation_duration_variability = 0.2,
        CPPrampdownDur = 0.1,
    )
        return new(drift_rate, event_onset, sensor_encoding_delay, accumulative_level_noise, boundary, motor_onset,
            motor_delay, post_accumulation_duration_mean, post_accumulation_duration_variability,
            CPPrampdownDur)
    end
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
                  :motor_delay => 0.4, :post_accumulation_duration_mean => 0.1, 
                  :post_accumulation_duration_variability => 0.2, :CPPrampdownDur => 0.1)
"""
function create_kelly_parameters_dict(model::KellyModel)
    return Dict(name => getfield(model, name) for name in fieldnames(typeof(model)))
end


"""
    KellyModel_simulate_cpp(rng, model::KellyModel, time_vec, Δt)

Generate response time and evidence of an evidence accumulation process using the Kelly model.

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
    evidence[time_vec .>= (model.event_onset+(rand(rng) -.5)*model.sensor_encoding_delay)] .= 1; 
    startAccT = time_vec[findfirst(evidence .== 1)];

    noise = vcat(zeros(
        sum(time_vec .< startAccT)),
        randn(rng,sum(time_vec .>=  startAccT)) .*  model.accumulative_level_noise .*sqrt(Δt));
    ev=evidence;
    ev[time_vec .< startAccT] .= 0; # set to zero all of the evidence before accT
    cum_evidence = cumsum(ev .* model.drift_rate .* Δt .+ noise,dims=1); # This is the cumulative differential evidence, just as in a 1d DDM. 
    
    
    # terminate the decision process on boundary crossing, record threshold-crossing samplepoint:
    cum_evidence = abs.(cum_evidence);
    dti = findfirst(cum_evidence .> model.boundary); # finding the sample point of threshold crossing of each, then will pick the earlier as the winner
    if isnothing(dti)  # Check if no crossing was found
        dti = length(time_vec)  # Set to the last time step
    end
    # now record RT in sec after adding motor time, with variability
    rt = time_vec[dti] + model.motor_onset + (rand(rng) - 0.5) * model.motor_delay

    # now make the CPP peak and go down linearly after a certain amount of post-dec accum time for this trial:
    post_acc_duration = model.post_accumulation_duration_mean .+ model.post_accumulation_duration_variability .* rand(rng);
    # so post_acc_duration is the post accumulation duration time, where the accumulation spikes over the threshold

    # acc_stop_index is the accumulation Stop index which is the index from the time Vector where the accumulation really stops
    acc_stop_index = dti + (post_acc_duration÷Δt)|>Int;
    # Take the absolute value of the accumulations
    cum_evidence = abs.(cum_evidence)
    if acc_stop_index < length(time_vec)
        nT = length(time_vec)
        tmp = cum_evidence[acc_stop_index] .- (1:(nT-acc_stop_index)) .*cum_evidence[acc_stop_index] .* (Δt ./ model.CPPrampdownDur)
        cum_evidence[(acc_stop_index+1):end] .= max.(Ref(0), tmp);
    end
    return rt / Δt, cum_evidence[1:end]
end


"""
    trace_sequential_sampling_model(rng, component::DriftComponent, design::AbstractDesign)

Generate response times and evidence Vectors of an given [`AbstractDesign`](@ref) with a [`DriftComponent`](@ref) which contains the model used for the simulation.

# Arguments
- `rng::StableRNG`: Random seed to ensure the same traces are created for reconstruction.
- `component::DriftComponent`: Component to specify the model and its parameters to simulate the evidence accumulation.
- `design::AbstractDesign`: design of the experiment preferable SequenceDesign.

# Returns
- `Vector{Float64}`: Simulated response times for the trials.
- `Matrix{Float64}`: evidence values over time for each trial. The output dimensions are `length(c.time_vec) x size(events, 1)`.

# Examples
```julia-repl
julia> model_parameter = Dict(:motor_onset => 0.4, :event_onset => 0.2);

julia> c = DriftComponent(0:1/500:1.0, 1/500, KellyModel, model_parameter);

julia> design_single = SingleSubjectDesign(conditions = Dict(:drift_rate => [0.5, 0.8], :condition => [1]));

julia> design_seq = SequenceDesign(design_single,"SCR_");

julia> simulate_component(StableRNG(1),c,design_seq)
Vector{Float64}, 501x6 Matrix{Float64}:
([96.65745162948949, 273.7368235451535, 271.86040880709123, 128.41057786118193, 342.35208862144276, 237.14773586760617], [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0])
```
"""
function trace_sequential_sampling_model(rng, component::DriftComponent, design::AbstractDesign)
    events = generate_events(deepcopy(rng), design)
    traces = Matrix{Float64}(undef, length(component.time_vec), size(events, 1))
    rts = Vector{Float64}(undef, size(events, 1))
    for (i, evt) in enumerate(eachrow(events))
        parameters = get_model_parameter(rng, evt, component.model_parameters)
        model = component.model_type(; (key => parameters[key] for key in keys(parameters))...)
        rt, evidence = SSM_Simulate(rng, model, component.Δt, component.time_vec)

        rts[i] = rt
        traces[:, i] = evidence[1:length(component.time_vec)]
    end
    return rts, traces
end

"""
    SSM_Simulate(rng, model::KellyModel, Δt, time_vec)

Generate response time and evidence Vector of length(time_vec) by using the Kelly Model for the simulation.

# Arguments
- `rng::StableRNG`: Random seed to ensure the same traces are created for reconstruction.
- `model::KellyModel`: SequentialSamplingModel to simulate the evidence and response time.
- `time_vec::StepRangeLen`: range of time steps for which the evidence is accumulated.
- `Δt::Float64`: size of the time steps.

# Returns
- `Float64`: Simulated response time for the trial.
- `Vector{Float64}`: evidence values over time. The output dimension is `length(c.time_vec)`.

# Examples
```julia-repl
julia> model = KellyModel()

julia> SSM_Simulate(StableRNG(1), model, 1/500, 0:1/500:1.0)
Float64, Vector{Float64}:
(96.65745162948949, [0.0 0.0 … 0.0 0.0])
```
"""
function SSM_Simulate(rng, model::KellyModel, Δt, time_vec)
    max_steps = length(time_vec)
    rt, evidence = KellyModel_simulate_cpp(rng, model, time_vec, Δt)
    if length(evidence) < max_steps
        final_value = 0
        append!(evidence, fill(final_value, max_steps - length(evidence)))
    end
    evidence = evidence[1:max_steps]
    return rt, evidence
end

"""
    SSM_Simulate(rng, model::LBA, Δt, time_vec)

Generate response time and evidence Vector of length(time_vec) by using the LBA for the simulation.

# Arguments
- `rng::StableRNG`: Random seed to ensure the same traces are created for reconstruction.
- `model::LBA`: SequentialSamplingModel to simulate the evidence and response time.
- `time_vec::StepRangeLen`: range of time steps for which the evidence is accumulated.
- `Δt::Float64`: size of the time steps.

# Returns
- `Float64`: Simulated response time for the trial.
- `Vector{Float64}`: evidence values over time. The output dimension is `length(c.time_vec)`.

# Examples
```julia-repl
julia> model = LBA()

julia> SSM_Simulate(StableRNG(1), model, 1/500, 0:1/500:1.0)
Float64, Vector{Float64}:
(96.65745162948949, [0.0 0.0 … 0.0 0.0])
```
"""
function SSM_Simulate(rng, model::LBA, Δt, time_vec)
    max_steps = length(time_vec)
    time_steps, evidence = SequentialSamplingModels.simulate(rng, model; Δt)
    evidence = hcat(evidence...)
    evidence = collect(vec(evidence))
    # Store results for this trial
    rt = (time_steps[end] + model.τ) / Δt
    if length(evidence) < max_steps
        final_value = 0
        append!(evidence, fill(final_value, max_steps - length(evidence)))
    else
        rt = (time_vec[end] + model.τ) / Δt
        evidence = evidence[1:max_steps]
    end
    return rt, evidence
end

"""
    SSM_Simulate(rng, model::DDM, Δt, time_vec)

Generate response time and evidence Vector of length(time_vec) by using the DDM for the simulation.

# Arguments
- `rng::StableRNG`: Random seed to ensure the same traces are created for reconstruction.
- `model::DDM`: SequentialSamplingModel to simulate the evidence and response time.
- `time_vec::StepRangeLen`: range of time steps for which the evidence is accumulated.
- `Δt::Float64`: size of the time steps.

# Returns
- `Float64`: Simulated response time for the trial.
- `Vector{Float64}`: evidence values over time. The output dimension is `length(c.time_vec)`.

# Examples
```julia-repl
julia> model = DDM()

julia> SSM_Simulate(StableRNG(1), model, 1/500, 0:1/500:1.0)
Float64, Vector{Float64}:
(96.65745162948949, [0.0 0.0 … 0.0 0.0])
```
"""
function SSM_Simulate(rng, model::DDM, Δt, time_vec)
    max_steps = length(time_vec)
    time_steps, evidence = SequentialSamplingModels.simulate(rng, model; Δt)
    evidence = evidence .- model.α * model.z

    # Store results for this trial
    rt = (time_steps[end] + model.τ) / Δt
    if length(evidence) < max_steps
        final_value = 0
        append!(evidence, fill(final_value, max_steps - length(evidence)))
    end
    evidence = evidence[1:max_steps]
    return rt, evidence
end
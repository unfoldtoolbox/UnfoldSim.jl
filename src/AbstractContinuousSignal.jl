"""
controlsignal: 
2 x n_timepoints Matrix containing x,y angles specifying n gaze direction vectors (measured from center gaze reference)
TODO docstring
"""
# -> also takes the simulation object since it might influence the generated controlsignal 
# always take controlsignal only of Gaze Direction Vectors: Matrix of size 3 x n_timepoints 
function simulate_continuoussignal(rng::AbstractRNG, s::EyeMovement, controlsignal::AbstractMatrix, sim::Simulation)
    headmodel = s.headmodel
    eye_model = s.eye_model
    return simulate_eyemovement(headmodel,controlsignal, eye_model)
end

"""
controlsignal: weights to each channel x timepoint? -> user can decide this themselves.
for now as an initial step it can directly be the simulated eeg itself.
#TODO docstring
"""
function simulate_continuoussignal(rng::AbstractRNG, s::PowerLineNoise, controlsignal::AbstractArray, sim::Simulation;)
    base_freq = s.base_freq
    harmonics = s.harmonics
    sampling_rate = s.sampling_rate
    weights_harmonics = s.weights_harmonics

    n_samples = size(controlsignal,2)
    k = 0:1:n_samples-1
    # TODO add check for nyquist criterion? sampling rate & freq -> warn or error?

    harmonics_signals = [sin.(2 * pi * (base_freq.*h)/sampling_rate .* k) for h in harmonics].*weights_harmonics
    # TODO: use controlsignal as weights for each channel x timepoint?
    
    return reduce(+,harmonics_signals)

end

# AbstractNoise: just returns empty Array - the noise simulation is handled separately in artifact-aware simulate()
function simulate_continuoussignal(rng::AbstractRNG, s::AbstractNoise, controlsignal::AbstractArray, sim::Simulation)
    return zeros(Float64, 0, 0)
end

# function simulate_continuoussignal(rng, s::TRF, controlsignal::)
#     # TBD
# end


# generate_controlsignal - returns controlsignal for the specified type of AbstractContinuousSignal
# -> also takes the simulation object since it might influence the generated controlsignal (e.g. generate blinks right after event A -> controlsignal depends on the design) 

# for EyeMovement, always return GazeDirectionVectors

function generate_controlsignal(rng::AbstractRNG, cs::GazeDirectionVectors, sim::Simulation)
    @assert size(cs.val)[1] == 3 "Please make sure gaze data has the shape 3 x n_timepoints."
    return cs.val
end

function generate_controlsignal(rng::AbstractRNG, cs::HREFCoordinates, sim::Simulation)
    return reduce(hcat,gazevec_from_angle_3d.(cs.val[1,:],cs.val[2,:])) # always return a 3 x time_points matrix
end

function generate_controlsignal(rng::AbstractRNG, cs::EyeMovement, sim::Simulation)
    return generate_controlsignal(rng, cs.controlsignal, sim)
end

# PowerLineNoise: returns empty Array for now.
# later, maybe weights for each channel x timepoint
function generate_controlsignal(rng::AbstractRNG, cs::PowerLineNoise, sim::Simulation)
    return zeros(Float64, 0, 0)
end

# AbstractNoise: returns empty Array for now, however in future the controlsignal for AbstractNoise may depend on the Simulation.
function generate_controlsignal(rng::AbstractRNG, cs::AbstractNoise, sim::Simulation)
    return zeros(Float64, 0, 0)
end

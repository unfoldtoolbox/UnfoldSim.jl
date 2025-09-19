"""
controlsignal: 
2 x n_timepoints Matrix containing x,y angles specifying n gaze direction vectors (measured from center gaze reference)
TODO docstring
"""
# -> also takes the simulation object since it might influence the generated controlsignal 
# always take controlsignal only of GazeDirectionVectors
function simulate_continuoussignal(rng::AbstractRNG, s::EyeMovement, controlsignal::GazeDirectionVectors, sim::Simulation)
    headmodel = s.headmodel
    return simulate_eyemovement(headmodel,controlsignal.coords; eye_model=s.eye_model)
end

function simulate_continuoussignal(rng::AbstractRNG, s::PowerLineNoise, controlsignal::AbstractArray, sim::Simulation;)
    base_freq = s.base_freq
    harmonics = s.harmonics
    sampling_rate = s.sampling_rate

    k = 0:1:size(controlsignal)-1 # assumes controlsignal is just 1D 
    # TODO add check for nyquist criterion? sampling rate & freq -> warn or error?

    n_samples = sampling_rate*length(controlsignal) + 1 # TODO check if we need the +1 or not depending on UnfoldSim conventions
    # TODO generate sinusoids at each harmonic

    # not sure yet if controlsignal should be multidimensional (e.g. different magnitudes for different channels => n_chan x time) 
    # or if we just take the simplest case: the same vector will be used for all channels.
    # weight the values at each point, if we need to have different relative strengths of harmonics or if we want to switch on/off the PLN
end

# AbstractNoise: doesn't simulate anything, returns empty Array. Since noise simulation will be handled separately in artifact-aware simulate()
function simulate_continuoussignal(rng::AbstractRNG, s::AbstractNoise, controlsignal::AbstractArray, sim::Simulation)
    
end

# function simulate_continuoussignal(rng, s::TRF, controlsignal::)
#     # TBD
# end



# sketching possibilities for generate_controlsignal

# generate_controlsignal - returns controlsignal for the specified type of AbstractContinuousSignal
# -> also takes the simulation object since it might influence the generated controlsignal (e.g. generate blinks right after event A -> controlsignal depends on the design) 

# for EyeMovement, always return GazeDirectionVectors

function generate_controlsignal(rng::AbstractRNG, cs::GazeDirectionVectors, sim::Simulation)
    return cs.val
end

function generate_controlsignal(rng::AbstractRNG, cs::HREFCoordinates, sim::Simulation)
    return gazevec_from_angle_3d.(cs.val[1,:],cs.val[2,:])
end

function generate_controlsignal(rng::AbstractRNG, s::EyeMovement, sim::Simulation)
    return generate_controlsignal(rng, s.controlsignal, sim)
end

# AbstractNoise: returns empty Array. Since noise simulation will be handled separately in artifact-aware simulate()
function generate_controlsignal(rng::AbstractRNG, s::AbstractNoise, sim::Simulation)
    return []
end

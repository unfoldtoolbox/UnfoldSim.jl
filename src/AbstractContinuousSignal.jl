"""
controlsignal: 
2 x n_timepoints Matrix containing x,y angles specifying n gaze direction vectors (measured from center gaze reference)
TODO docstring
"""
# -> also takes the simulation object since it might influence the generated controlsignal 
# always take controlsignal only of Gaze Direction Vectors: Matrix of size 3 x n_timepoints 
function simulate_continuoussignal(rng::AbstractRNG, s::EyeMovement, controlsignal::AbstractMatrix, sim::Simulation)
    @show typeof(s), objectid(rng)
    println(rand(rng))
    headmodel = s.headmodel
    eye_model = s.eye_model
    return simulate_eyemovement(headmodel,controlsignal, eye_model)
end

"""
controlsignal: actual weights (channel x time) to be applied to the simulated single-channel PLN

controlsignal ---> n_channels x n_timepoints
#TODO docstring
"""
function simulate_continuoussignal(rng::AbstractRNG, s::PowerLineNoise, controlsignal::AbstractMatrix, sim::Simulation)
    base_freq = s.base_freq
    harmonics = s.harmonics
    sampling_rate = s.sampling_rate
    weights_harmonics = s.weights_harmonics

    n_samples = size(controlsignal)[end]
    k = 0:1:n_samples-1
    # TODO add check for nyquist criterion? sampling rate & freq -> warn or error?

    harmonics_signals = [sin.(2 * pi * (base_freq.*h)/sampling_rate .* k) for h in harmonics].*weights_harmonics
    
    # @show typeof(harmonics_signals[1]), size(controlsignal), size(weights_harmonics)
    return reduce(+,harmonics_signals)' .*controlsignal
end


function simulate_continuoussignal(rng::AbstractRNG, s::AbstractNoise, controlsignal::AbstractArray, sim::Simulation)
    return reshape(simulate_noise(rng,s,length(controlsignal)),size(controlsignal)) .* controlsignal
end


# generate_controlsignal - returns controlsignal for the specified type of AbstractContinuousSignal
# -> also takes the simulation object since it might influence the generated controlsignal (e.g. generate blinks right after event A -> controlsignal depends on the design) 

# for EyeMovement, always return a matrix containing gaze direction vectors

function generate_controlsignal(rng::AbstractRNG, cs::GazeDirectionVectors, sim::Simulation)
    @assert size(cs.val)[1] == 3 "Please make sure gaze data has the shape 3 x n_timepoints."
    return cs.val
end

function generate_controlsignal(rng::AbstractRNG, cs::HREFCoordinates, sim::Simulation)
    return reduce(hcat,gazevec_from_angle_3d.(cs.val[1,:],cs.val[2,:])) # always return a 3 x time_points matrix
end

function generate_controlsignal(rng::AbstractRNG, cs::AbstractContinuousSignal, sim::Simulation)
    return generate_controlsignal(deepcopy(rng), cs.controlsignal, sim)
end

function generate_controlsignal(rng::AbstractRNG, cs::PowerLineNoise, sim::Simulation)
    return nothing
end

# AbstractNoise: returns nothing for now, however in future the controlsignal for AbstractNoise may depend on the Simulation.
function generate_controlsignal(rng::AbstractRNG, cs::AbstractNoise, sim::Simulation)
    return nothing
end

# identity function - in case the user already specified the final controlsignal while creating the artifact
function generate_controlsignal(rng::AbstractRNG, cs::AbstractMatrix, sim::Simulation)
    return cs
end


# # for drifts / allowing user to pass in a vector of values for individual channels
# # Dict: channel => weight array
# function generate_controlsignal(rng::AbstractRNG, cs::Dict{Int,Array{Int}}, sim::Simulation)
#     # place the given weight arrays together into a common zero matrix at the specified channels.
#     # return controlsignal matrix
# end

# for drifts / allowing user to pass in a function for individual channels
# function generate_controlsignal(rng::AbstractRNG, cs::Dict{Int,function}, sim::Simulation)
#     # generate the weight arrays from the function
#     # place them into a common zero matrix at the specified channels.
#     # return controlsignal matrix
# end


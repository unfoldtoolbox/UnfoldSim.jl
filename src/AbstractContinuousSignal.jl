abstract type AbstractContinuousSignal end

struct EyeMovement <: AbstractContinuousSignal
    controlsignal
    headmodel
    # events # <-- from realdata (or from controlsignal?) or passed in by user. will be added into the events dataframe returned by simulation function
end

struct TRF <: AbstractContinuousSignal
    controlsignal
    # TBD
end

struct PowerLineNoise <: AbstractContinuousSignal
    controlsignal
    base_freq::Float64 = 50
    harmonics::Array{Int} = [1 3 5]
    sampling_rate::Float64 = 500
end

"""
controlsignal: 
2 x n_timepoints Matrix containing x,y angles specifying n gaze direction vectors (measured from center gaze reference)
TODO docstring
"""
function simulate_continuoussignal(rng::AbstractRNG,s::EyeMovement,controlsignal::HREFCoordinates, headmodel::AbstractHeadmodel, eye_model="crd")
    return simulate_eyemovement(headmodel,gazevec_from_angle_3d.(controlsignal[1,:],controlsignal[2,:]); eye_model=eye_model)
end

function simulate_continuoussignal(rng::AbstractRNG,s::PowerLineNoise,controlsignal::AbstractArray, base_freq::Float64, harmonics::Array{Int} = [1 3 5], sampling_rate::Float64 = 500.0)
    k = 0:1:size(controlsignal)-1 # assumes controlsignal is just 1D 
    # TODO add check for nyquist criterion? sampling rate & freq -> warn or error?

    n_samples = sampling_rate*length(controlsignal) + 1 # TODO check if we need the +1 or not depending on UnfoldSim conventions
    # TODO generate sinusoids at each harmonic

    # not sure yet if controlsignal should be multidimensional (e.g. different magnitudes for different channels => n_chan x time) 
    # or if we just take the simplest case: the same vector will be used for all channels.
    # weight the values at each point, if we need to have different relative strengths of harmonics or if we want to switch on/off the PLN
end

# function simulate_continuoussignal(rng, s::TRF, controlsignal::)
#     # TBD
# end



# sketching possibilities for generate_controlsignal

function generate_controlsignal(rng::AbstractRNG, c::GazeDirectionVectors)
    return c
end

function generate_controlsignal(rng::AbstractRNG, c::HREFCoordinates)
    return gazevec_from_angle_3d.(c[1,:],c[2,:])
end

abstract type AbstractControlSignal end

struct HREFCoordinates <: AbstractControlSignal end

struct GazeDirectionVectors <: AbstractControlSignal end
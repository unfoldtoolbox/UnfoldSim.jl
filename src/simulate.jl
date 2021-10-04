mutable struct UnfoldSim

    coef # event x coefficients
    noise # either vector of event, or function
    overlap # event x onset
    signal # event x function


end
function simulate(ufS::UnfoldSim)
    s = genEmptySignals(ufS) # an object that automatically extends itself by 0s - not sure we need that
    s = distributeEvents(s,ufS.overlap)
    s = scaleSticks(s,ufS.coef)
    s = convolve(s,ufS.signal)
    s = addNoise(s,ufS.noise)
    s = sumSignals(s)
end

# returns a function to generate an event
function signal()
    
    
end
    
function kernel(n,peakTime,peakWidth)
    #n = 100
    #peakTime = 20
    #peakWidth = 6
    sig = zeros(n)
    kern = DSP.Windows.hanning(peakWidth)
    sig[(peakTime-Int(ceil(peakWidth/ 2))):(peakTime+Int(floor(peakWidth/2))-1)] = kern
    return sig
end
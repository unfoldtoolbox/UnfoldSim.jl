# here we define some commonly used basis for simulation

## EEG
# P100, N170, P300, N400
p100(;fs=100) = hanning(0.1,0.1,fs)
p300(;fs=100) = hanning(0.3,0.3,fs)
n170(;fs=100) = -hanning(0.15,0.17,fs)
n400(;fs=100) = -hanning(0.4,0.4,fs)

"""
generate a hanning window

duration: in s
offset: in s, defines hanning peak
fs: sampling rate in Hz
"""
function DSP.hanning(duration,offset,fs)
    signal = hanning(Int(round(duration*fs)))
    return padarray(signal,-Int(round(offset*fs/2)),0)
end

## pupil

function PuRF(;n = 10.1,tmax = 0.93,fs=100)
    t = (0:1/fs:3*tmax)
    @show t
    return PuRF(t,n,tmax)./PuRF(tmax,n,tmax)
end

function PuRF(t,n,tmax)
    return t.^n .* exp.(-n .* t ./ (tmax))
end
## fMRI
"""
Generate a HRF kernel. 

TR = 1/fs
default parameters taken from SPM

Code adapted from Unfold.jl
"""
function hrf(;TR=1,peak=6.,psUnder=16,length = 32.,peak_width=1.,psUnder_width=1,amplitude=6,shift=0)

    # code adapted from SPM12b
    mt = Int(max(length,32)) # simulation maxtime
    dt = TR / mt
    
    #duration = 1.0 / mt
    
    
    box_mt = ones(mt)
    u = range(0, stop = ceil(length / dt)) .- shift / dt

    
    # Note the inverted scale parameter compared to SPM12.
    g1 = Gamma(peak ./ peak_width, peak_width ./ dt)
    #spm_Gpdf(u,p(2)/p(4),dt/p(4))
    g2 = Gamma(psUnder ./ psUnder_width, psUnder_width./ dt)
    # g1 - g2/p(5);
    hrf = pdf.(g1, u) .- pdf.(g2, u) ./ amplitude

    # hrf = hrf([0:floor(p(7)/RT)]*fMRI_T + 1);
    hrf = hrf ./ sum(hrf)
    hrf = conv(box_mt, hrf)'
    
    hrf = hrf[((range(0, stop = Int(floor(length ./ TR)))*mt)).+1]
    return hrf
end
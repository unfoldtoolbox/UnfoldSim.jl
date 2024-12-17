# here we define some commonly used basis for simulation

## EEG
"""
    p100(;sfreq=100)
Generator for Hanning window, peak at 100ms, width 100ms, at kwargs `sfreq` (default 100). Returns a vector.
"""
p100(; sfreq = 100) = hanning(0.1, 0.1, sfreq)

"""
    p300(;sfreq=100)
Generator for Hanning window, peak at 300ms, width 300ms, at kwargs `sfreq` (default 100). Returns a vector.
"""
p300(; sfreq = 100) = hanning(0.3, 0.3, sfreq)

"""
    n170(;sfreq=100)
Generator for Hanning window, negative (!) peak at 170ms, width 150ms, at kwargs `sfreq` (default 100). Returns a vector.
"""
n170(; sfreq = 100) = -hanning(0.15, 0.17, sfreq)

"""
    n400(;sfreq=100)
Generator for Hanning window, negative (!) peak at 400ms, width 400ms, at kwargs `sfreq` (default 100). Returns a vector.
"""
n400(; sfreq = 100) = -hanning(0.4, 0.4, sfreq)

"""
generate a hanning window

width: in s
offset: in s, defines hanning peak, must be > round(width/2)
sfreq: sampling rate in Hz
"""
function DSP.hanning(width, offset, sfreq)
    width = width * sfreq
    offset = offset * sfreq
    signal = hanning(Int(round(width)))
    pad_by = Int(round(offset - length(signal) / 2))

    pad_by < 0 ? error("offset has to be > round(width/2)") : ""
    return pad_array(signal, -pad_by, 0)
end

## pupil
"""
    PuRF()
Default generator for PuRF Pupil Response Function.
"""
function PuRF(; n = 10.1, tmax = 0.93, sfreq = 100)
    t = (0:1/sfreq:3*tmax)
    return PuRF(t, n, tmax) ./ PuRF(tmax, n, tmax)
end

function PuRF(t, n, tmax)
    return t .^ n .* exp.(-n .* t ./ (tmax))
end
## fMRI
"""
Generate a HRF kernel. 

TR = 1/sfreq
default parameters taken from SPM

Code adapted from Unfold.jl
"""
function hrf(;
    TR = 1,
    peak = 6.0,
    post_undershoot = 16,
    length = 32.0,
    peak_width = 1.0,
    post_undershoot_width = 1,
    amplitude = 6,
    shift = 0,
)

    # code adapted from SPM12b
    mt = Int(max(length, 32)) # simulation maxtime
    dt = TR / mt

    #duration = 1.0 / mt


    box_mt = ones(mt)
    u = range(0, stop = ceil(length / dt)) .- shift / dt


    # Note the inverted scale parameter compared to SPM12.
    g1 = Gamma(peak ./ peak_width, peak_width ./ dt)
    #spm_Gpdf(u,p(2)/p(4),dt/p(4))
    g2 = Gamma(post_undershoot ./ post_undershoot_width, post_undershoot_width ./ dt)
    # g1 - g2/p(5);
    hrf = pdf.(g1, u) .- pdf.(g2, u) ./ amplitude

    # hrf = hrf([0:floor(p(7)/RT)]*fMRI_T + 1);
    hrf = hrf ./ sum(hrf)
    hrf = conv(box_mt, hrf)'

    hrf = hrf[((range(0, stop = Int(floor(length ./ TR)))*mt)).+1]
    return hrf
end

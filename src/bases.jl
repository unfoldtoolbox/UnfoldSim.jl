# Here we define some commonly used basis functions for simulation

## EEG
"""
<<<<<<< HEAD
    p100(;sfreq=100)
Generator for Hanning window, peak at 100ms, width 100ms, at kwargs `sfreq` (default 100). Returns a vector.
=======
    p100(; sfreq = 100)

Generate a Hanning window mimicking a P100 EEG component with a peak at 100ms and a width of 100ms.

# Keyword arguments
- `sfreq = 100`: Sampling frequency in Hz.

# Returns
- `Vector`: Contains a shifted (i.e. zero-padded) hanning window.

# Examples
```julia-repl
julia> p100(; sfreq = 200)
30-element Vector{Float64}:
 0.0
 0.0
 0.0
 0.0
 0.0
 ⋮
 0.37725725642960045
 0.22652592093878665
 0.10542974530180327
 0.02709137914968268
 0.0
```

See also [`p300`](@ref), [`n170`](@ref), [`n400`](@ref), [`hanning`](@ref).
>>>>>>> main
"""
p100(; sfreq = 100) = hanning(0.1, 0.1, sfreq)

"""
<<<<<<< HEAD
    p300(;sfreq=100)
Generator for Hanning window, peak at 300ms, width 300ms, at kwargs `sfreq` (default 100). Returns a vector.
=======
    p300(; sfreq = 100)

Generate a Hanning window mimicking a P300 EEG component with a peak at 300ms and a width of 300ms.

# Keyword arguments
- `sfreq = 100`: Sampling frequency in Hz.

# Returns
- `Vector`: Contains a shifted (i.e. zero-padded) hanning window.

# Examples
```julia-repl
julia> p300(; sfreq = 150)
67-element Vector{Float64}:
 0.0
 0.0
 0.0
 0.0
 0.0
 ⋮
 0.07937323358440934
 0.04518400232274078
 0.02025351319275137
 0.005089279059533658
 0.0
```

See also [`p100`](@ref), [`n170`](@ref), [`n400`](@ref), [`hanning`](@ref).
>>>>>>> main
"""
p300(; sfreq = 100) = hanning(0.3, 0.3, sfreq)

"""
<<<<<<< HEAD
    n170(;sfreq=100)
Generator for Hanning window, negative (!) peak at 170ms, width 150ms, at kwargs `sfreq` (default 100). Returns a vector.
=======
    n170(; sfreq = 100)

Generate a Hanning window mimicking an N170 EEG component with a negative (!) peak at 170ms and a width of 150ms.

# Keyword arguments
- `sfreq = 100`: Sampling frequency in Hz.

# Returns
- `Vector`: Contains a shifted (i.e. zero-padded) hanning window.

# Examples
```julia-repl
julia> n170(; sfreq = 120)
28-element Vector{Float64}:
 -0.0
 -0.0
 -0.0
 -0.0
 -0.0
  ⋮
 -0.45386582026834904
 -0.2771308221117309
 -0.1304955413896705
 -0.03376388529782215
 -0.0
```
See also [`p100`](@ref), [`p300`](@ref), [`n400`](@ref), [`hanning`](@ref).
>>>>>>> main
"""
n170(; sfreq = 100) = -hanning(0.15, 0.17, sfreq)

"""
<<<<<<< HEAD
    n400(;sfreq=100)
Generator for Hanning window, negative (!) peak at 400ms, width 400ms, at kwargs `sfreq` (default 100). Returns a vector.
=======
    n400(; sfreq = 100)

Generate a Hanning window mimicking an N400 EEG component with a negative (!) peak at 400ms and a width of 400ms.

# Keyword arguments
- `sfreq = 100`: Sampling frequency in Hz.

# Returns
- `Vector`: Contains a shifted (i.e. zero-padded) hanning window.

# Examples
```julia-repl
julia> n400(; sfreq = 250)
150-element Vector{Float64}:
 -0.0
 -0.0
 -0.0
 -0.0
 -0.0
  ⋮
 -0.016025649301821876
 -0.009035651368646647
 -0.00402259358460233
 -0.0010066617640578368
 -0.0
```
See also [`p100`](@ref), [`p300`](@ref), [`n170`](@ref), [`hanning`](@ref).
>>>>>>> main
"""
n400(; sfreq = 100) = -hanning(0.4, 0.4, sfreq)

"""
    hanning(duration, offset, sfreq)

Generate a (potentially shifted) hanning window with a certain duration. 

Note: This function extends the `DSP.hanning` function using multiple dispatch.

# Arguments
- `duration`: in s.
- `offset`: in s, defines hanning peak i.e. shift of the hanning window.
- `sfreq`: Sampling rate in Hz.

# Returns
- `Vector`: Contains a shifted (i.e. zero-padded) hanning window.

# Examples
```julia-repl
julia> UnfoldSim.hanning(0.1, 0.3, 100)
25-element Vector{Float64}:
 0.0
 0.0
 0.0
 0.0
 0.0
 ⋮
 0.9698463103929542
 0.75
 0.41317591116653485
 0.11697777844051105
 0.0
```
"""
function DSP.hanning(duration, offset, sfreq)
    signal = hanning(Int(round(duration * sfreq)))
    return pad_array(signal, -Int(round(offset * sfreq / 2)), 0)
end

<<<<<<< HEAD
## pupil
"""
    PuRF()
Default generator for PuRF Pupil Response Function.
=======

## Pupil
"""
    PuRF(; n = 10.1, tmax = 0.93, sfreq = 100)

Default generator for PuRF Pupil Response Function. The canonical PRF is a gamma function and implemented according to Denison 2020 equation (2) going back to Hoeks & Levelt, 1993.

The pupil response is evaluated at t = 0:1/sfreq:3*tmax. The response is normalized by the peak-maximum at tmax, thus typically a pupil-response of 1 is returned (disregarding numerical issues).


# Keyword arguments:
- `n = 10.1`: shape parameter
- `tmax = 0.93`: peak maximum
- `sfreq = 100`: sampling frequency

# Returns
- `Vector`: canonical pupil response with length(0:1/sfreq:3*tmax) entries.

# Examples
```julia-repl
julia> PuRF(; n = 5)
280-element Vector{Float64}:
 0.0
 2.0216617815131253e-8
 6.130689396024061e-7
 4.4118063684811444e-6
 1.761817666835793e-5
 ⋮
 0.012726253506722554
 0.012280989091455786
 0.011850525657416842
 0.011434405338911133
 0.011032182932283816
```
>>>>>>> main
"""
function PuRF(; n = 10.1, tmax = 0.93, sfreq = 100)
    t = (0:1/sfreq:3*tmax)
    return _PuRF(t, n, tmax) ./ _PuRF(tmax, n, tmax)
end

function _PuRF(t, n, tmax)
    return t .^ n .* exp.(-n .* t ./ (tmax))
end


## fMRI
"""
    hrf(;
    TR = 1,
    peak = 6.0,
    post_undershoot = 16,
    length = 32.0,
    peak_width = 1.0,
    post_undershoot_width = 1,
    amplitude = 6,
    shift = 0)

Generate a parameterized BOLD haemodynamic response function (HRF) kernel based on gamma-functions.

Implementation and default parameters were taken from the SPM-toolbox.

Note: TR = 1/sfreq

# Keyword arguments
- `TR = 1`: repetition time, 1/sfreq.
- `length = 32.0`: total length of the kernel in seconds.
- `amplitude = 6`: maximal amplitude.
- `peak = 6.0`: peak timing.
- `peak_width = 1.0`: width of the peak.
- `post_undershoot = 16`: post-undershoot timing.
- `post_undershoot_width = 1`: post-undershoot width.
- `shift = 0`: shift the whole HRF.

# Returns
- `Vector`: HRF BOLD response.

# Examples
```julia-repl
julia> hrf()
33-element Vector{Float64}:
  0.0
  0.0007715994433635659
  0.019784004131204957
  0.08202939459091822
  0.158157713522699
  ⋮
 -0.0006784790038792572
 -0.00042675060451877775
 -0.000263494738348278
 -0.00015990722628360688
 -9.548780093799345e-5
```
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

    # Code adapted from SPM12b
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

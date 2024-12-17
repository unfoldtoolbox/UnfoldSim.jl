"""
    Hartmut <: AbstractHeadmodel

Type for accessing the HArtMuT model (Harmening et al., 2022) which is a head model including not only brain sources but also muscular and ocular sources.

# Fields
- `artefacual::Any`: Dict with artefactual sources (e.g. muscular and ocular) containing label, leadfield, orientation and position.
- `cortical::Any`:  Dict with cortical sources containing label, leadfield, orientation and position.
- `electrodes::Any`: Dict with electrode labels and their positions in 3D space.

# Examples
```julia-repl
julia> h = Hartmut()
Please cite: HArtMuT: Harmening Nils, Klug Marius, Gramann Klaus and Miklody Daniel - 10.1088/1741-2552/aca8ce
HArtMuT-Headmodel
227 electrodes:  (AF3,AF3h...) - hartmut.electrodes
2004 source points: (Left Middle Temporal Gyrus, posterior division,...) - hartmut.cortical
4260 source points: (Muscle_DepressorLabii_left,...) - hartmut.artefactual

In addition to UnfoldSim.jl please cite:
HArtMuT: Harmening Nils, Klug Marius, Gramann Klaus and Miklody Daniel - 10.1088/1741-2552/aca8ce

# Access artefactual sources
julia> h.artefactual
Dict{String, Any} with 4 entries:
  "label"       => ["Muscle_DepressorLabii_left", "Muscle_DepressorLabii_left", "Muscle_DepressorLabii_left", "Muscle_DepressorLabii_left", "Muscle_DepressorLabii_left", "Muscle_DepressorLabii_left", "Muscle_DepressorLabii_left", "Mu…
  "leadfield"   => [-0.00148682 -0.00229078 … -0.307231 -0.321305; 0.0141537 0.0134523 … -0.141177 -0.140583; … ; 0.108747 0.1091 … 0.106614 0.117126; 0.0257134 0.0248186 … 0.00105064 -0.00433068;;; 0.0898798 0.0891799 … 2.02466 1.81…
  "orientation" => [0.54244 0.482395 -0.687789; 0.546949 0.507161 -0.666059; … ; 0.193693 -0.979684 0.0519768; 0.327641 -0.944196 -0.0338583]
  "pos"         => [-29.3728 58.6061 -120.829; -28.4183 59.0185 -121.448; … ; -21.8088 70.5968 -28.7679; -21.1545 70.1282 -31.4184]
```
"""
struct Hartmut <: AbstractHeadmodel
    artefactual::Any
    cortical::Any
    electrodes::Any
end

function Hartmut() # Outer constructor
    println("""Please cite: $(hartmut_citation())""")
    path = joinpath(artifact"hartmut", "hartmut.h5")
    h = h5open(path)


    weirdchan = ["Nk1", "Nk2", "Nk3", "Nk4"]
    ## getting index of these channels from imported hartmut model data, exclude them in the topoplot
    remove_indices = findall(l -> l ∈ weirdchan, h["electrodes"] |> read |> x -> x["label"])

    function select_channel(x)

        if "leadfield" ∈ keys(x)
            x["leadfield"] = x["leadfield"][Not(remove_indices), :, :] .* 10e3 # this scaling factor seems to generate potentials with +-1 as max
        else
            x["label"] = x["label"][Not(remove_indices)]
            pos3d = x["pos"][Not(remove_indices), :]
            pos3d = pos3d ./ (4 * maximum(pos3d, dims = 1))
            x["pos"] = pos3d
        end
        return x
    end
    headmodel = Hartmut(
        h["artefacts"] |> read |> select_channel,
        h["cortical"] |> read |> select_channel,
        h["electrodes"] |> read |> select_channel,
    )
end

function Base.show(io::IO, h::Hartmut)
    src_label = h.cortical["label"]
    ele_label = h.electrodes["label"]
    art_label = h.artefactual["label"]

    println(
        io,
        """HArtMuT-Headmodel
        $(length(ele_label)) electrodes:  ($(ele_label[1]),$(ele_label[2])...) - hartmut.electrodes
        $(length(src_label)) source points: ($(src_label[1]),...) - hartmut.cortical
        $(length(art_label)) source points: ($(art_label[1]),...) - hartmut.artefactual

        In addition to UnfoldSim.jl please cite:
        $(hartmut_citation())
        """,
    )
end

"Return the citation string for HArtMuT."
hartmut_citation() =
    "HArtMuT: Harmening Nils, Klug Marius, Gramann Klaus and Miklody Daniel - 10.1088/1741-2552/aca8ce"

"""
    leadfield(hart::Hartmut; type = "cortical")

Return the leadfield for the (cortical or artefacual) sources of the HArtMuT model.

# Keyword arguments
- `type = "cortical"`: Defines whether the "cortical" or "artefactual" leadfield should be returned.

# Returns
- `Array{Float64, 3}`: Leadfield values i.e. how much each source contributes to the potential measured at each electrode.
    The output dimensions are `electrodes x sources x spatial dimension`.

# Examples
```julia-repl
julia> h = Hartmut();
Please cite: HArtMuT: Harmening Nils, Klug Marius, Gramann Klaus and Miklody Daniel - 10.1088/1741-2552/aca8ce

julia> lf = leadfield(h);

julia> size(lf)
(227, 2004, 3)

# Access the leadfield for one spatial dimension
julia> lf[:,:,1]
227×2004 Matrix{Float64}:
  0.151429    0.0284341  …  -0.119927     -0.158114
  0.1732      0.0441432     -0.110515     -0.165316
  0.249592    0.10857       -0.000593027  -0.0206122
  0.245206    0.104854      -0.0200251    -0.055392
  0.126496    0.0118467     -0.128248     -0.146107
  0.253092    0.111688   …   0.02277       0.0184387
  ⋮                      ⋱                
  0.20306     0.138837       0.133486      0.18334
 -0.689154   -1.00904       -0.108276     -0.105398
  0.192729    0.148448       0.0924756     0.142322
 -1.26181    -1.59936       -0.140598     -0.133054
  0.213982    0.145953   …   0.115515      0.170698
  0.0731569   0.0794415      0.0485819     0.107929
```
"""
leadfield(hart::Hartmut; type = "cortical") =
    type == "cortical" ? hart.cortical["leadfield"] : hart.artefactual["leadfield"]

"""
    orientation(hart::Hartmut; type = "cortical")

Return the orientations of the (cortical or artefacual) sources of the HArtMuT model. 

The norm of the orientation vectors is 1 and the values are between -1 and 1.

# Keyword arguments
- `type = "cortical"`: Defines whether the "cortical" or "artefactual" orientations should be returned.

# Returns
- `Matrix{Float64}`: Orientations in 3D space. The output dimensions are `sources x spatial dimensions`.

# Examples
```julia-repl
julia> h = Hartmut();
Please cite: HArtMuT: Harmening Nils, Klug Marius, Gramann Klaus and Miklody Daniel - 10.1088/1741-2552/aca8ce

julia> orientation(h)
2004×3 Matrix{Float64}:
 -0.921919   0.364292   0.131744
 -0.900757  -0.415843  -0.125345
 -0.954087   0.117479  -0.27553
 -0.814613  -0.55344    0.17352
 -0.790526   0.276849  -0.546281
  ⋮                    
 -0.962905   0.20498   -0.17549
 -0.828358   0.557468  -0.0552498
 -0.963785  -0.265607  -0.0239074
 -0.953909   0.203615   0.22045
 -0.75762    0.128027   0.640016
```
"""
orientation(hart::Hartmut; type = "cortical") =
    type == "cortical" ? hart.cortical["orientation"] : hart.artefactual["orientation"]


@deprecate headmodel() Hartmut()

"""
    magnitude(headmodel::AbstractHeadmodel)

Extract magnitude of the orientation-including leadfield.

By default use the orientation specified in the headmodel
Fallback: along the third dimension using `norm` - the maximal projection

# Returns
- `Matrix{Float64}`: ... The output dimensions are `electrodes x sources`.
"""
magnitude(headmodel::AbstractHeadmodel) = magnitude(leadfield(headmodel))

"""
    magnitude(headmodel::Hartmut; type = "perpendicular")

Extract magnitude of 3-orientation-leadfield, `type` (default: "perpendicular") => uses the provided source-point orientations - otherwise falls back to `norm`.

# Keyword arguments
- `type = "perpendicular"`: 

# Examples
```julia-repl
```
"""
magnitude(headmodel::Hartmut; type = "perpendicular") =
    type == "perpendicular" ? magnitude(leadfield(headmodel), orientation(headmodel)) :
    magnitude(leadfield(headmodel))

"""
    magnitude(lf::AbstractArray{T,3}, orientation::AbstractArray{T,2}) where {T<:Real}

Return the magnitude along an orientation of the leadfield.

# Arguments
- `lf::AbstractArray{T,3}`:
- `orientation::AbstractArray{T,2}`:

# Examples
```julia-repl
```
"""
function magnitude(lf::AbstractArray{T,3}, orientation::AbstractArray{T,2}) where {T<:Real}
    si = size(lf)
    magnitude = fill(NaN, si[1:2])
    for e = 1:si[1]
        for s = 1:si[2]
            magnitude[e, s] = lf[e, s, :]' * orientation[s, :]
        end
    end
    return magnitude
end

"""
    magnitude(lf::AbstractArray{T,3}) where {T<:Real}

If orientation is not specified, returns the maximal magnitude (norm of leadfield).

# Arguments:
- `lf::AbstractArray{T,3}`:

# Examples
```julia-repl
```
"""
function magnitude(lf::AbstractArray{T,3}) where {T<:Real}
    si = size(lf)
    magnitude = fill(NaN, si[1:2])
    for e = 1:si[1]
        for s = 1:si[2]
            magnitude[e, s] = norm(lf[e, s, :])
        end
    end
    return magnitude
end

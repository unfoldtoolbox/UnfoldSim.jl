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

julia> leadfield(h)
227×2004×3 Array{Float64, 3}:
[:, :, 1] =
  0.151429    0.0284341   0.0381777  -0.031533    0.173189    0.00485958   0.0490966   0.0707025   0.00818434  …  -0.029391   -0.480143  -0.140076   -0.405818  -0.0802413  -0.480738    -0.205534   -0.119927     -0.158114
  0.1732      0.0441432   0.0518934  -0.0101328   0.181351    0.0204942    0.0610806   0.0848026   0.0261072      -0.0259066  -0.48539   -0.0838597  -0.423522  -0.0628599  -0.48615     -0.163832   -0.110515     -0.165316
  0.249592    0.10857     0.115235    0.10436     0.211371    0.10154      0.125745    0.145194    0.106947        0.0252987  -0.142679   0.373272   -0.20579    0.0663492  -0.243599     0.199307   -0.000593027  -0.0206122
  0.245206    0.104854    0.110464    0.0929579   0.209626    0.0941604    0.119648    0.140814    0.102867        0.0164751  -0.218093   0.317887   -0.261457   0.0519611  -0.294938     0.155672   -0.0200251    -0.055392
  0.126496    0.0118467   0.0241024  -0.0535207   0.164124   -0.011041     0.037146    0.0555699  -0.00916642     -0.0319598  -0.465564  -0.190625   -0.380701  -0.0963106  -0.467464    -0.241251   -0.128248     -0.146107
  ⋮                                                           ⋮                                                ⋱   ⋮                                                         ⋮                                     
 -0.689154   -1.00904    -0.640231   -0.803683   -0.669248   -0.989008    -0.924508   -0.780335   -0.659238       -0.134393   -0.13584   -0.253791   -0.010676  -0.259598    0.00974083  -0.609214   -0.108276     -0.105398
  0.192729    0.148448    0.131563    0.19185     0.132415    0.164303     0.14867     0.155177    0.157946        0.120672    0.261752   0.189576    0.130332   0.300065    0.144634     0.553918    0.0924756     0.142322
 -1.26181    -1.59936    -1.57901    -1.32513    -1.15289    -1.44583     -1.5914     -1.34047    -0.99371        -0.124783   -0.225718  -0.321099   -0.058764  -0.275479   -0.0440558   -0.679299   -0.140598     -0.133054
  0.213982    0.145953    0.135357    0.185493    0.155274    0.159222     0.152356    0.160551    0.154056    …   0.123235    0.303627   0.301454    0.13918    0.309234    0.149528     0.621073    0.115515      0.170698
  0.0731569   0.0794415   0.0814297   0.190411    0.0333774   0.140429     0.105861    0.086699    0.0747807       0.0283169   0.180094  -0.0449803   0.137933  -0.0166636   0.118638     0.0344215   0.0485819     0.107929

[:, :, 2] =Return the leadfield for the (cortical or artefacual) sources of the HArtMuT model.
  0.591764     0.471349    0.573102    0.478932   0.607152     0.502843    0.472437    0.608999    0.526104     0.540295   …   0.364791    0.712709   0.532255    1.05677    0.335067    1.09307     0.600539   0.418957    0.782385
  0.543546     0.436331    0.534807    0.414284   0.559653     0.459095    0.432672    0.542914    0.50076      0.480347       0.363619    0.768387   0.564811    1.11043    0.336465    1.14023     0.624134   0.423507    0.826305
  0.332102     0.262255    0.310747    0.14202    0.347224     0.244054    0.230912    0.255358    0.285465     0.227095       0.385679    0.720729   0.58951     0.970872   0.320419    0.961721    0.640017   0.377979    0.777559
  0.349808     0.276847    0.333381    0.162083   0.364274     0.262838    0.249247    0.277439    0.309421     0.245992       0.379013    0.773062   0.606652    1.03547    0.324358    1.02911     0.653921   0.391519    0.820855
  0.64474      0.506963    0.609941    0.549661   0.657857     0.54828     0.513286    0.680338    0.544963     0.605602       0.366573    0.650619   0.495855    0.9933     0.332749    1.03381     0.573857   0.412086    0.733588
  ⋮                                                            ⋮                                                           ⋱   ⋮                                                         ⋮                                 
 -0.996154    -1.84228    -1.38123    -0.458889  -1.13876     -1.25072    -1.09428    -1.2008     -0.926254    -1.02257       -0.175818   -0.195926  -0.168892   -0.210932  -0.078534   -0.195056   -0.301431  -0.0385821  -0.249137
 -0.00458421  -0.046744   -0.042647   -0.246742  -0.00330122  -0.108935   -0.0846609  -0.169111   -0.0787297   -0.129365      -0.0836064  -0.196058  -0.1195     -0.270287  -0.138095   -0.274416   -0.126704  -0.113781   -0.259266
 -1.0759       0.284423   -0.749177    0.720501  -1.28191      0.537338    0.0481887   0.0145133   0.212858    -0.366328       0.0101172  -0.10103   -0.118584   -0.107945   0.0493172  -0.0989716  -0.133399   0.0586076  -0.0935456
  0.0553323    0.0195319   0.0230949  -0.172482   0.0608517   -0.0368073  -0.0236916  -0.0873639  -0.00752649  -0.0633672  …   0.10244    -0.098749  -0.0492562  -0.192389  -0.0118357  -0.212741    0.046137  -0.034545   -0.107372
 -0.211627    -0.225928   -0.270457   -0.436503  -0.182126    -0.349264   -0.288703   -0.396514   -0.194699    -0.302807      -0.246299   -0.326628  -0.188054   -0.289232  -0.204373   -0.304826   -0.347969  -0.281474   -0.296293

[:, :, 3] =
  0.271254   0.0957242   0.235673    0.0724304    0.229273    0.12759     0.186508     0.194947    0.0488128   …   0.0745411   0.500774   -0.124186     0.310166    0.239008    0.443472     0.0841403   0.480781   0.126728
  0.296812   0.117088    0.257683    0.0778473    0.254908    0.136939    0.199821     0.220435    0.0896733       0.100678    0.585576   -0.120079     0.372598    0.261439    0.534567     0.0981147   0.528447   0.159212
  0.237229   0.0763321   0.164717   -0.0145197    0.169283    0.0370844   0.118791     0.151308    0.108431        0.0869617   0.506642   -0.136966     0.218055    0.2576      0.414411     0.116566    0.447633   0.0846618
  0.260593   0.0959073   0.193354    0.00676149   0.197061    0.062048    0.141588     0.177994    0.128168        0.111609    0.591143   -0.128143     0.30101     0.276576    0.51247      0.123914    0.504287   0.126872
  0.232089   0.0660315   0.201022    0.0601517    0.188268    0.110154    0.164513     0.157126   -0.00060752      0.0441638   0.406901   -0.128795     0.24055     0.213069    0.344352     0.0689211   0.422925   0.0893109
  ⋮                                                           ⋮                                                ⋱   ⋮                                                            ⋮                                  
  0.104551   0.306119   -0.0314079   1.13434      0.182324    1.31835     0.410955     0.660495   -1.02091         0.0726301  -0.061279    0.0754638    0.0190903  -0.021282   -0.00128803  -0.0624448   0.180966   0.0313662
  0.121412   0.0182608   0.0259512  -0.0729851    0.0625917  -0.0685532   0.00716254   0.0612821   0.0831421       0.105609   -0.0598804   0.0601159   -0.0583101   0.0379451  -0.0286753    0.0313028   0.112618  -0.0552693
  0.148321   0.356874    0.458193    1.4662       0.310998    1.50553     0.81209      0.960325   -1.27739         0.0422736  -0.0393025   0.0572244    0.0378037   0.032761    0.0159193   -0.048176    0.199276   0.0319228
  0.130557   0.0182334   0.0383301  -0.0806093    0.066031   -0.0653128   0.0160211    0.0613658   0.0816704   …   0.0812635  -0.0336645   0.0369646   -0.0670275   0.097711   -0.0188773    0.0565507   0.119389  -0.0709868
 -0.151815  -0.261587   -0.311147   -0.239542    -0.201341   -0.353552   -0.280478    -0.273946   -0.243196       -0.324282   -0.356306    0.00222309  -0.240861   -0.328941   -0.257704    -0.114429   -0.411469  -0.198202
```
"""
leadfield(hart::Hartmut; type = "cortical") =
    type == "cortical" ? hart.cortical["leadfield"] : hart.artefactual["leadfield"]

"""
    orientation(hart::Hartmut; type = "cortical")

Return the orientations of the (cortical or artefacual) sources of the HArtMuT model.

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

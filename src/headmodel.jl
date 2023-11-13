

struct Hartmut  <: AbstractHeadmodel
    artefactual
    cortical
    electrodes
end

function Base.show(io::IO,h::Hartmut)
    src_label = h.cortical["label"]
    ele_label = h.electrodes["label"]
    art_label = h.artefactual["label"]
    
    println(io,"""HArtMuT-Headmodel
    $(length(ele_label)) electrodes:  ($(ele_label[1]),$(ele_label[2])...) - hartmut.electrodes
    $(length(src_label)) source points: ($(src_label[1]),...) - hartmut.cortical
    $(length(art_label)) source points: ($(art_label[1]),...) - hartmut.artefactual
    
    In addition to UnfoldSim.jl please cite:
    $(hartmut_citation())
    """)
end

"Returns citation-string for HArtMuT"
hartmut_citation() = "HArtMuT: Harmening Nils, Klug Marius, Gramann Klaus and Miklody Daniel - 10.1088/1741-2552/aca8ce"

"Returns the leadfield"
leadfield(hart::Hartmut;type="cortical") = type == "cortical" ? hart.cortical["leadfield"] : hart.artefactual["leadfield"]
orientation(hart::Hartmut;type="cortical") = type == "cortical" ? hart.cortical["orientation"] : hart.artefactual["orientation"]


"""
Load a headmodel, using Artifacts.jl automatically downloads the required files

Currently only `type="hartmut"` is implemented
"""
function headmodel(;type="hartmut")
    if type == "hartmut"
        println("""Please cite: $(hartmut_citation())""")
        path = joinpath(artifact"hartmut", "hartmut.h5")
        h = h5open(path)

               
        weirdchan = ["Nk1","Nk2","Nk3","Nk4"]
        ## getting index of these channels from imported hartmut model data, exclude them in the topoplot
        remove_indices = findall(l -> l ∈ weirdchan, h["electrodes"]|>read|>x->x["label"]);
        
        function sel_chan(x)
            
            if "leadfield" ∈     keys(x)
                x["leadfield"] = x["leadfield"][Not(remove_indices),:,:].*10e3 # this scaling factor seems to generate potentials with +-1 as max
            else
                x["label"] = x["label"][Not(remove_indices)]
                pos3d= x["pos"][Not(remove_indices),:]
                pos3d = pos3d ./ (4*maximum(pos3d,dims=1))
                x["pos"]  = pos3d
            end
            return x
        end
        headmodel = Hartmut(
                    h["artefacts"]|>read|>sel_chan,
                    h["cortical"]|>read|>sel_chan,
                    h["electrodes"]|>read|>sel_chan,
                    )
    else
        error("unknown headmodel. currently only 'hartmut' allowed")
    end
    
    return headmodel
end

"""
Extracts magnitude of the orientation-including leadfield.

By default uses the orientation specified in the headmodel

Fallback: along the third dimension using `norm` - the maximal projection
"""
magnitude(headmodel::AbstractHeadmodel) = magnitude(leadfield(headmodel))

"""
Extract magnitude of 3-orientation-leadfield, 
`type` (default: "perpendicular") => uses the provided source-point orientations - otherwise falls back to `norm`
"""
magnitude(headmodel::Hartmut;type="perpendicular") = type == "perpendicular" ? magnitude(leadfield(headmodel),orientation(headmodel)) : magnitude(leadfield(headmodel))


function magnitude(lf::AbstractArray{T,3},orientation::AbstractArray{T,2}) where {T<:Real}
    si = size(lf);
    magnitude = fill(NaN,si[1:2]);
    for e=1:si[1] 
    	for s=1:si[2]
	    	magnitude[e,s] = lf[e,s,:]' * orientation[s,:]
	    end
    end
    return magnitude
end


function magnitude(lf::AbstractArray{T,3}) where {T<:Real}
    si = size(lf);
    magnitude = fill(NaN,si[1:2]);
    for e=1:si[1] 
    	for s=1:si[2]
	    	magnitude[e,s] = norm(lf[e,s,:]);
	    end
    end
    return magnitude
end

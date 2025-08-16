#=
TODOs

- check all the packages used and whether they are actually needed
- do we need ___?
    - read_new_hartmut()
    - pos2dfrom3d() (it's technically just 2 lines. see unfoldsim docs multichannel example)

- docstrings for all functions
=#



# For helper function to read in the large model
# using HDF5
# using DataFrames

# for angle conversion - need to check if we still need this or if the simple angle conversion util code will suffice
# using CartesianFromSpherical


"""
TODO docstring
Construct the gaze vector in 3-D world coordinates, given the gaze angle (in degrees) as measured from the center gaze (looking straight ahead) in the world x-y plane. Z-component is always 0.  
"""
function gazevec_from_angle(angle_deg)
	# just x,y plane for now. gaze angle measured from front neutral gaze, not from x-axis
	return [sind(angle_deg) cosd(angle_deg) 0]
end

"""
TODO docstring
Calculate the gaze vector in 3-D world coordinates, given the horizontal and vertical angles (in degrees) from center gaze direction i.e. looking straight ahead.
"""
function gazevec_from_angle_3d(angle_H, angle_V)
	# angles measured from center gaze position => use complementary angle for θ 
	return Array{Float32}(CartesianFromSpherical()(Spherical(1, deg2rad(90-angle_H), deg2rad(angle_V))))
end


"""
    read_eyemodel(;p::String = "HArtMuT_NYhead_extra_eyemodel.mat")

Read and return the HArtMuT eye model from the given file path.

# Fields
- `p::String` (optional): Path of the file from which to read the eye model. 

# Returns
- `model`: Eye model read from the specified file path.
"""
function read_eyemodel(;p::String = "HArtMuT_NYhead_extra_eyemodel.mat")
    file = matopen(p);
    hartmut = read(file,"eyemodel")
    close(file)
    hartmut["label"] = hartmut["labels"]
    delete!(hartmut,"labels")
    hartmut
end

"""
TODO docstring

For each of a given set of labels, find the indices of all the sources in the given model which have that label.

Arguments:
- `headmodel`: Head model containing the sources in which to search for labels
- `labels::Vector{String}`: Labels to search for. 

Returns:
- `labelsourceindices::Dict`: Dict with the keys being the label and the values being the indices of sources with that label in the given model. 
"""
function hart_indices_from_labels(headmodel,labels::Vector{String}=["dummy"])
    labelsourceindices = Dict()
    for l in labels
        labelsourceindices[l] = findall(k->occursin(l,k),headmodel["label"][:])
    end
    return labelsourceindices
end


"""
TODO docstring
Calculate orientations from the given positions with respect to the reference. `direction` is by default "away" from the reference point; set `towards` to calculate orientations away from the reference point.
"""
function calc_orientations(reference, positions; towards=false, away=true)
    if towards
        orientation_vecs = reference .- positions
    else
        orientation_vecs = positions .- reference
    end
    return orientation_vecs ./ norm.(eachrow(orientation_vecs))
end


"""
TODO docstring
Calculate the angle between two 3-D vectors.
"""
function angle_between(a,b) 
    return acosd.(dot(a, b)/(norm(a)*norm(b)))
end

# """
# TODO docstring
# For the given gaze direction, return a weight value for each source point in the model: 1 for a cornea source, 0 for retina sources and all other points in the model. 

# # Arguments
# - `model`: Head model.
# - `sim_idx`: Indices of the (eye) source points for which weights should be calculated.
# - `gazedir`: The gaze direction vector.
# - `max_cornea_angle_deg`: Angle (in degrees) made by the radial vector containing the cornea point farthest from the cornea center with the radial vector containing the cornea center.

# # Returns
# - `eyeweights`: Weights of the source points: 1 for a cornea source, 0 for retina sources and all other points in the model.
# """
# function weights_from_gazedir(model, sim_idx, gazedir, max_cornea_angle_deg)
# 	eyeweights = zeros(size(model["pos"])[1]) # all sources other than those defined by sim_idx will be set to zero magnitude 
# 	eyeweights[sim_idx] .= mapslices(x -> is_corneapoint(x,gazedir,max_cornea_angle_deg), model["orientation"][sim_idx,:],dims=2)
# 	return eyeweights
# end


"""
    is_corneapoint(orientation::Vector{Float64}, gazedir::Vector{Float64}, max_cornea_angle_deg::Float64)


"""
function is_corneapoint(orientation::Vector{Float64}, gazedir::Vector{Float64}, max_cornea_angle_deg::Float64)
	if(angle_between(orientation,gazedir)<=max_cornea_angle_deg)
		return 1
	else 
		return -1
	end
end


"""
TODO docstring
"""
function ensemble_leadfield(eyemodel, sim_idx::Vector{Int}, gazedir::Vector{Float64}, max_cornea_angle_deg::Float64)
	mag_model = magnitude(eyemodel["leadfield"],eyemodel["orientation"])
	source_weights = zeros(size(eyemodel["pos"])[1]) # all sources other than those defined by sim_idx will be set to zero magnitude 
	source_weights[sim_idx] .= mapslices(x -> is_corneapoint(x,gazedir,max_cornea_angle_deg), eyemodel["orientation"][sim_idx,:],dims=2)
	
	weighted_sum = sum(mag_model[:,idx].* source_weights[idx] for idx in sim_idx,dims=2)
	return weighted_sum
end


"""
TODO docstring
Read the eye model from a mat file, remove channels Nk1-Nk4, and calculate left/right eye indices, eye centers, orientations and add them to the eyemodel as properties. 
Return imported model, labelsourceindices of eyemodel
"""
function import_eyemodel(; labels=[
    r"EyeRetina_Choroid_Sclera_left$" # match the end of the search term, or else it also matches "leftright" sources.
    ,r"EyeRetina_Choroid_Sclera_right$"
    ,"EyeCornea_left" # eyemodel - only one set of sources per eye, and the labels are different, compared to the full HArtMuT model
    ,"EyeCornea_right"
    ,"EyeCenter_left"
    ,"EyeCenter_right"
], modelpath="HArtMuT_NYhead_extra_eyemodel_hull_mesh8_2025-03-01.mat"
)   
    eyemodel = read_eyemodel(; p=modelpath)
    remove_indices = [164, 165, 166, 167] # since eyemodel structure doesn't exactly correspond to the main hartmut mat structure expected by the read_new_hartmut function, just get the indices of the electrodes that it drops & drop the same indices from eyemodel directly 
    eyemodel["leadfield"] = eyemodel["leadfield"][Not(remove_indices), :, :]

    lsi_eyemodel = hart_indices_from_labels(eyemodel,labels) 

    # add electrode positions to eyemodel - import small headmodel temporarily since the eyemodel does not contain them
	hart_small = Hartmut()
	pos3d = hart_small.electrodes["pos"]
    eyemodel["electrodes"] = deepcopy(hart_small.electrodes)
    eyemodel["electrodes"]["pos2d"] = pos2dfrom3d(pos3d)

    # add indices and positions of certain sets of sources, to make them easier to access
    eyemodel["eyeleft_idx"] = [ lsi_eyemodel["EyeCornea_left"] ; lsi_eyemodel[r"EyeRetina_Choroid_Sclera_left$"] ] # eyemodel_left_idx
    eyemodel["eyeright_idx"] = [ lsi_eyemodel["EyeCornea_right"] ; lsi_eyemodel[r"EyeRetina_Choroid_Sclera_right$"] ] # eyemodel_right_idx
    eyemodel["eyecenter_left_pos"] = eyemodel["pos"][lsi_eyemodel["EyeCenter_left"],:] # eye_center_L
    eyemodel["eyecenter_right_pos"] = eyemodel["pos"][lsi_eyemodel["EyeCenter_right"],:] # eye_center_R
    eyemodel["eyecenter_left_idx"] = lsi_eyemodel["EyeCenter_left"]
    eyemodel["eyecenter_right_idx"] = lsi_eyemodel["EyeCenter_right"]

    # add orientations (by default all set to zero)
	eyemodel["orientation"] = zeros(size(eyemodel["pos"]))

    return eyemodel, lsi_eyemodel
end


"""
TODO docstring
Calculate the sum of leadfields of just the source points at the given indices. 
"""
function selected_sources_eeg(eyemodel,src_idx,equiv_orientations)
	# equiv_ori_model = model["orientation"]
	# equiv_ori_model[idx,:] .= equiv_orientations[1:length(idx),:]
	mag_eyemodel_equiv = magnitude(eyemodel["leadfield"],eyemodel["orientation"])
	mag = sum(mag_eyemodel_equiv[:,ii] for ii in src_idx)
end


function find_avg_gazedir()
    # finding cornea centers and approximate gaze direction (as mean of cornea orientations)
    cornea_center_R = Statistics.mean(eyemodel["pos"][lsi_eyemodel["EyeCornea_right"],:],dims=1)
    cornea_center_L = Statistics.mean(eyemodel["pos"][lsi_eyemodel["EyeCornea_left"],:],dims=1)
    gazedir_R = mean(eyemodel["orientation"][lsi_eyemodel["EyeCornea_right"],:], dims=1).*10
    gazedir_L = mean(eyemodel["orientation"][lsi_eyemodel["EyeCornea_left"],:], dims=1).*10
end

#TODO auto-find cornea max. angle from center and cornea point positions


"""
TODO docstring
"""
function simulate_eyemovement(eyemodel, gazevectors::Vector{Vector{Float64}}; time::Vector{Number}, samplingrate::Number, crd=false, ensemble=false)
    if (!crd && !ensemble)
        @warn "Neither CRD nor Ensemble model selected. No simulation performed."
        return []
    end
    
    eyemodel, lsi_eyemodel = import_eyemodel()

    # leadfields: matrix with dimensions [electrodes x n_gazepoints] 
    leadfields = zeros(size(eyemodel["pos"])[1],length(gazevectors))

    if crd
        # select eye centers as source points; set orientation = gazevector; simulate
        for ix in eachindex(gazevectors)
            eyemodel["orientation"][eyemodel["eyecenter_left_idx"]] = gazevectors[ix]
            eyemodel["orientation"][eyemodel["eyecenter_right_idx"]] = gazevectors[ix]
            leadfields_crd[:,ix] = selected_sources_eeg(eyemodel,[eyemodel["eyecenter_left_idx"] eyemodel["eyecenter_right_idx"]],gazevectors[ix])
        end

    elseif ensemble
        # select cornea and retina sources to simulate; set orientation; simulate EEG for each of the gaze vectors

        sim_srcs_idx = [eyemodel["eyeleft_idx"]; eyemodel["eyeright_idx"]] # indices in eyemodel, of the points which we want to use to simulate data, i.e. retina & cornea. used for ensemble simulation. 

        # calculate eyeball source orientations away from center, later use negative weightage for the points where the source dipole points towards the center (i.e. the retina points).
        eyemodel["orientation"][eyemodel["eyeleft_idx"],:] = calc_orientations(eyemodel["eyecenter_left_pos"], eyemodel["pos"][eyemodel["eyeleft_idx"],:]; towards=true)
        eyemodel["orientation"][eyemodel["eyeright_idx"],:] = calc_orientations(eyemodel["eyecenter_right_pos"], eyemodel["pos"][eyemodel["eyeright_idx"],:]; towards=true) 

        for ix in eachindex((gazevectors))
            # for each gazepoint, calculate the leadfield and store it in the corresponding column
            leadfields[:,ix] = ensemble_leadfield(eyemodel, sim_srcs_idx, gazevectors[ix], 54.0384)
        end

    end

    # placeholder: matching the simulated leadfields with time points?

    leadfields
end
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
Calculate orientations from the given positions with respect to the reference. `direction` is by default "towards" the reference point; if a different value is given, orientations are calculated away from the reference point.
"""
function calc_orientations(reference, positions; direction="towards")
    if direction == "towards"
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
Simulate eye movement start to end using ensemble method 
"""
function sim_gazevecs_ensemble(gazevectors::Vector{Vector{Float64}}) # TODO rename this and turn it into a proper function. For now, moved here from pluto notebook
# simulate for gazevectors - ensemble method

    sim_srcs_idx = [eyemodel["eyeleft_idx"]; eyemodel["eyeright_idx"]] # indices in eyemodel, of the points which we want to use to simulate data, i.e. retina & cornea. used for ensemble simulation. 

    # leadfields: matrix of dimensions [electrodes x n_gazepoints] 
	leadfields_ensemble = zeros(227,length(gazevectors))
	for ix in eachindex((gazevectors)) # 1:length(gazevectors) #TODO delete this comment once it's proven to work 
		# for each gazepoint, calculate the leadfield and store it in the corresponding column
		leadfields_ensemble[:,ix] = ensemble_leadfield(eyemodel, sim_srcs_idx,
			gazevectors[ix], 54.0384) #leadfield_from_gazedir
	end
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

    # labels for which to find source indices
    # labels = [
    #             r"EyeRetina_Choroid_Sclera_left$" # match the end of the search term, or else it also matches "leftright" sources.
	# 			,r"EyeRetina_Choroid_Sclera_right$"
	# 			,"EyeCornea_left" # new spherical eyemodel - only one set of sources per eye and the labels are different
	# 			,"EyeCornea_right" # same as above
	# 			,"EyeCenter_right"
	# 			,"EyeCenter_left"
	# 		]
   
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
"""
function simulate_crd(eyemodel, gazevectors::Vector{Vector{Float64}}) # TODO add type & dimensions
    leadfields_crd = zeros(227,length(gazevectors)) # TODO replace hardcoded 227 with number of channels in model?
	eyecenter_idx = [lsi_eyemodel["EyeCenter_left"][1],lsi_eyemodel["EyeCenter_right"][1]] #TODO remove duplication of calculating center indices
	

    #original from notebook
	# for ix in eachindex(gazevectors)
	# 	# calculate leadfield due to 2 individual CRD sources and add
	# 	leadfields_crd[:,ix] = (
	# 		leadfield_specific_sources_orientations(deepcopy(eyemodel),eyecenter_idx,gazevectors[ix])
	# 	)
	# end

    # not yet tested
    leadfields_crd[:,ix] = (
			leadfield_specific_sources_orientations(deepcopy(eyemodel),eyecenter_idx,gazevectors[ix])
		for ix in eachindex(gazevectors) )

    return leadfields_crd
end


"""
TODO docstring
"""
function leadfield_specific_sources_orientations(model,idx,equiv_orientations)
	# take just a selected subset of points in the model, along with a new orientation for those points, and calculate the sum of leadfields of just these points with the given orientation. 
    equiv_ori_model = model["orientation"]
    
    for ii in idx
        equiv_ori_model[ii,:] = equiv_orientations[:]
    end
    
    mag_eyemodel_equiv = magnitude(eyemodel["leadfield"],equiv_ori_model)
    mag = sum(mag_eyemodel_equiv[:,ii] for ii in idx)
    return mag
end


"""
TODO docstring
Take just a selected subset of points in the model, along with new orientations for those points, 
and calculate the sum of leadfields of just these points with the new orientations. 
"""
function equiv_dipole_mag(model,idx,equiv_orientations)
	equiv_ori_model = model["orientation"]
	equiv_ori_model[idx,:] .= equiv_orientations[1:length(idx),:]
	mag_eyemodel_equiv = magnitude(eyemodel["leadfield"],equiv_ori_model)
	mag = sum(mag_eyemodel_equiv[:,ii] for ii in idx)
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
    # leadfields: matrix of dimensions [electrodes x n_gazepoints] 
    leadfields = zeros(size(eyemodel["pos"])[1],length(gazevectors))

    if crd
        # select eye centers as points
        # set orientation = gazedir
        leadfields = simulate_crd(eyemodel,gazevectors)
    elseif ensemble
        # select cornea and retina sources to simulate; set orientation; simulate EEG for each of the gaze vectors

        sim_srcs_idx = [eyemodel["eyeleft_idx"]; eyemodel["eyeright_idx"]] # indices in eyemodel, of the points which we want to use to simulate data, i.e. retina & cornea. used for ensemble simulation. 

        # calculate eyeball source orientations away from center, later use negative weightage for the points where the source dipole points towards the center (i.e. the retina points).
        eyemodel["orientation"][eyemodel["eyeleft_idx"],:] = calc_orientations(eyemodel["eyecenter_left_pos"], eyemodel["pos"][eyemodel["eyeleft_idx"],:]; direction="away")
        eyemodel["orientation"][eyemodel["eyeright_idx"],:] = calc_orientations(eyemodel["eyecenter_right_pos"], eyemodel["pos"][eyemodel["eyeright_idx"],:]; direction="away") 

        # leadfields: matrix with dimensions [electrodes x n_gazepoints]
        for ix in eachindex((gazevectors))
            # for each gazepoint, calculate the leadfield and store it in the corresponding column
            leadfields[:,ix] = ensemble_leadfield(eyemodel, sim_srcs_idx, gazevectors[ix], 54.0384)
        end

    end

    # placeholder: matching the simulated leadfields with time points?

    leadfields
end
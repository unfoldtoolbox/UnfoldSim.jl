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
	# just x,y plane. gaze angle measured from front neutral gaze, not from x-axis
	return [sind(angle_deg) cosd(angle_deg) 0]
end

"""
TODO docstring
Calculate the gaze vector in 3-D world coordinates, given the horizontal and vertical angles (in degrees) from center gaze direction i.e. looking straight ahead.
"""
function gazevec_from_angle_3d(angle_H, angle_V)
	# angles measured from center gaze position => use complementary angle for θ; 
    # ϕ is already measured upwards from the x-y plane, according to the conventions of CoordinateTransformations.jl
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
Calculate orientation vectors from the given positions with respect to the reference. The vectors are normalized to have length 1. 
`direction` is by default "away" from the reference point; set `direction="towards"` to calculate orientations towards the reference point.
"""
function calc_orientations(reference, positions; direction::String="away")
    if direction=="towards"
        orientation_vecs = reference .- positions
    else 
        orientation_vecs = positions .- reference
    end
    return orientation_vecs ./ norm.(eachrow(orientation_vecs))
end


"""
TODO docstring
Return the angle (in degrees) between two 3-D vectors in Cartesian coordinates.
"""
function angle_between(a,b) 
    return acosd.(dot(a, b)/(norm(a)*norm(b)))
end


"""
    is_corneapoint(orientation::Vector{Float64}, gazedir::Vector{Float64}, max_cornea_angle_deg::Float64)

Return 1 if the angle between the given orientation vector and the gaze direction vector is less than or equal to the given maximum cornea angle; else return -1.

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
	
	weighted_sum = sum(mag_model[:,idx].* source_weights[idx] for idx in sim_idx)
	return weighted_sum
end


"""
TODO docstring 
Single gaze direction. Assumes orientations are already set in eyemodel
"""
function generate_eyegaze_eeg(eyemodel, src_idx::Vector{Int}, weights::Vector{Int})

    mag = magnitude(eyemodel["leadfield"],eyemodel["orientation"])
    signal = sum(mag[:,idx].* weights[idx] for idx in src_idx)

    return signal
end


"""
Simulate the eye movement and return the corresponding scalp potentials at each of the channels, 
given a head model, an array of gaze direction vectors defining the eye movement, and optionally an eye-model type. 

# Arguments
- `headmodel`: Head model to be used.
- `gazevectors::Vector{Vector{Float64}}`: Vector of gaze direction vectors (in Cartesian coordinates) to define the eye movement trajectory. The vector should point from the eyes towards the point at which the gaze is directed. 
- `eye_model::String` (optional): Choice of model to use for the eye. Options available are "crd" (default) and "ensemble".

# Returns
- 
    
TODO docstring; type for headmodel; auto-conversion from x,y angles to gaze direction vectors
"""
function simulate_eyemovement(headmodel, gazevectors::GazeDirectionVectors; eye_model::String="crd")
    # when gaze direction vector changes, 
    # CRD: orientation changes, weight stays the same (=1 for all points)
    # Ensemble: orientation stays the same, weight changes (=1 for cornea points, -1 for retina points, calculated based on gaze direction)

    # leadfields: matrix with dimensions [electrodes x n_gazepoints] 
    leadfields = zeros(size(headmodel["pos"])[1],length(gazevectors))

    # weights: matrix of dimensions [n_channel x n_gazevec]

    if eye_model=="crd"
        weights = ones(size(headmodel["pos"])[1])
        src_idx = [headmodel["eyecenter_left_idx"] headmodel["eyecenter_right_idx"]]
        # select eye centers as source points; set orientation = gazevector
        for ix in eachindex(gazevectors)
            headmodel["orientation"][headmodel["eyecenter_left_idx"]] = gazevectors[ix]
            headmodel["orientation"][headmodel["eyecenter_right_idx"]] = gazevectors[ix]
            # or just use src_idx - will help later if we want to simulate just for one eye or so

            leadfields[:,ix] = generate_eyegaze_eeg(headmodel, src_idx, weights)
        end

    elseif eye_model=="ensemble"
        weights = zeros(size(headmodel["pos"])[1])

        # select cornea and retina sources to simulate; set orientation; simulate EEG for each of the gaze vectors

        # indices in headmodel, of the points which we want to use to simulate data, i.e. retina & cornea. used for ensemble simulation. 
        src_idx = [headmodel["eyeleft_idx"]; headmodel["eyeright_idx"]]
        
        
        for ix in eachindex((gazevectors))
            # weight changes depending on retina/cornea type based on current gazevector
            weights[src_idx] .= mapslices(x -> is_corneapoint(x,gazedir,max_cornea_angle_deg), headmodel["orientation"][src_idx,:])
            # change the above line to separately calculate weights for left and right eye if giving gazepoints and calculating gaze vectors separately. or just run simulate_eyemovement twice, once with gazevectors based on left eye and once with right eye? 

            leadfields[:,ix] = generate_eyegaze_eeg(headmodel, src_idx, weights)
        end

    else
        @error "Neither CRD nor Ensemble model selected. No eye movement simulation performed."
    end

    return leadfields
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

    # calculate eye surface sources' orientations away from the respective eye centers
    # later use negative weightage for the points where the source dipole points towards the center (i.e. the retina points).  
    eyemodel["orientation"][eyemodel["eyeleft_idx"],:] = calc_orientations(eyemodel["eyecenter_left_pos"], eyemodel["pos"][eyemodel["eyeleft_idx"],:])
    eyemodel["orientation"][eyemodel["eyeright_idx"],:] = calc_orientations(eyemodel["eyecenter_right_pos"], eyemodel["pos"][eyemodel["eyeright_idx"],:]) 

    return eyemodel
end

"""
TODO docstring

"""
function find_avg_gazedir()
    # finding cornea centers and approximate gaze direction (as mean of cornea orientations)
    cornea_center_R = Statistics.mean(eyemodel["pos"][lsi_eyemodel["EyeCornea_right"],:],dims=1)
    cornea_center_L = Statistics.mean(eyemodel["pos"][lsi_eyemodel["EyeCornea_left"],:],dims=1)
    gazedir_R = mean(eyemodel["orientation"][lsi_eyemodel["EyeCornea_right"],:], dims=1).*10
    gazedir_L = mean(eyemodel["orientation"][lsi_eyemodel["EyeCornea_left"],:], dims=1).*10
end

#TODO auto-find cornea max. angle from center and cornea point positions


function a_z_simulation()
    # import hartmut model - modified with new eye points
    eyemodel::AbstractHeadmodel = import_eyemodel()

    # import href gaze coordinates
    sample_data = example_data_eyemovements()
    href_trajectory::HREFCoordinates = sample_data[1:2,:]

    # setup basic ingredients for simulate
    #       design, .... , noise

        # Ques: should it be possible to simulate just for eyemovement, without specifying any design etc? Just passing in the eye trajectory?

    # call simulate with href coords
    simulate(d, c, o, [EyeMovement(href_trajectory, eyemodel, nothing) NoNoise()]) # TODO clarify!

    # plot simulated data
end
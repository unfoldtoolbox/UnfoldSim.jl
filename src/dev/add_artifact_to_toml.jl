using Pkg.Artifacts, ArtifactsUtils

# Get path to the Artifact.toml
artifacts_toml = joinpath(@__DIR__, "Artifacts.toml")

# Check if noise artefact exists in Artefact.toml and get hash
noise_hash = artifact_hash("RealNoise", artifacts_toml)

# Check if has does not exist (i.e. Artifact not added) or if artifact doesn't exist on disk;
# create/ add to .toml if fullfilled
if isnothing(noise_hash) || !artifact_exists(noise_hash)
    add_artifact!(
                     "Artifacts.toml",
                     "RealNoise",
                     "https://github.com/cormullion/juliamono/releases/download/v0.030/JuliaMono.tar.gz", # Change this line to actual dataset
                     force=true,
                     lazy=true, # If lazy is set to true, even if download information is available, this artifact will not be downloaded until it is accessed via the artifact"name" syntax, or ensure_artifact_installed() is called upon it.
                    )
end

#=
For the specific use case of using artifacts that were previously bound, 
we have the shorthand notation artifact"name" which will automatically 
search for the Artifacts.toml file contained within the current package, look up the given artifact by name,
install it if it is not yet installed, then return the path to that given artifact.
=#


# The following can be used to install the artifact/ ensure it is installed
#import Pkg; Pkg.ensure_artifact_installed("RealNoise", "Artifacts.toml")

# The below can be used to get the path to the artifact file
# According to docs this should also install the Artifact if it is not already installed
# artifact"RealNoise"
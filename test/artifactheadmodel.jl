using Base: make_atomicreplace
using MAT
using HDF5


dlpath = download(
    "https://github.com/harmening/HArtMuT/raw/main/HArtMuTmodels/HArtMuT_NYhead_small.mat",
);
file = matopen(dlpath);
hartmut = read(file, "HArtMuT");
close(file)

fname = tempname(); # temporary file

fid = h5open(fname, "w")


label = string.(hartmut["electrodes"]["label"][:, 1])
chanpos = Float64.(hartmut["electrodes"]["chanpos"])

cort_label = string.(hartmut["cortexmodel"]["labels"][:, 1])
cort_orient = hartmut["cortexmodel"]["orientation"]
cort_leadfield = hartmut["cortexmodel"]["leadfield"]
cort_pos = hartmut["cortexmodel"]["pos"]


art_label = string.(hartmut["artefactmodel"]["labels"][:, 1])
art_orient = hartmut["artefactmodel"]["orientation"]
art_leadfield = hartmut["artefactmodel"]["leadfield"]
art_pos = hartmut["artefactmodel"]["pos"]

e = create_group(fid, "electrodes")
e["label"] = label
e["pos"] = chanpos
c = create_group(fid, "cortical")
c["label"] = cort_label
c["orientation"] = cort_orient
c["leadfield"] = cort_leadfield
c["pos"] = cort_pos
a = create_group(fid, "artefacts")
a["label"] = art_label
a["orientation"] = art_orient
a["leadfield"] = art_leadfield
a["pos"] = art_pos

close(fid)
mkdir(joinpath(tempdir(), "artifact"))
mv(fname, joinpath(tempdir(), "artifact", "hartmut.h5"))
using ArtifactUtils
artifact_id = artifact_from_directory("/tmp/artifact")
gist = upload_to_gist(artifact_id)
add_artifact!("Artifacts.toml", "hartmut", gist)

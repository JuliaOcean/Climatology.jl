using Pkg.Artifacts
using OceanStateEstimation

# This is the path to the Artifacts.toml we will manipulate
p=dirname(pathof(OceanStateEstimation))
artifact_toml = joinpath(p, "../Artifacts.toml")
#artifact_toml = "Artifacts.toml"

ECCOclim_hash = artifact_hash("ECCOclim", artifact_toml)
# If the name was not bound, or the hash it was bound to does not exist, create it!
if ECCOclim_hash == nothing || !artifact_exists(ECCOclim_hash)
    # create_artifact() returns the content-hash of the artifact directory once we're finished creating it
    ECCOclim_hash = create_artifact() do artifact_dir
        cp(joinpath(p,"../examples/nctiles_climatology.csv"), joinpath(artifact_dir, "nctiles_climatology.csv"))
    end

    # Now bind that hash within our `Artifacts.toml`.
    bind_artifact!(artifact_toml, "ECCOclim", ECCOclim_hash)
end

OCCAclim_hash = artifact_hash("OCCAclim", artifact_toml)
# If the name was not bound, or the hash it was bound to does not exist, create it!
if OCCAclim_hash == nothing || !artifact_exists(OCCAclim_hash)
    OCCAclim_hash = create_artifact() do artifact_dir
        cp(joinpath(p,"../examples/OCCA_climatology.csv"), joinpath(artifact_dir, "OCCA_climatology.csv"))
    end
    # Now bind that hash within our `Artifacts.toml`.
    bind_artifact!(artifact_toml, "OCCAclim", OCCAclim_hash)
end

# Get the path of the iris dataset, either newly created or previously generated.
# this should be something like `~/.julia/artifacts/dbd04e28be047a54fbe9bf67e934be5b5e0d357a`
ECCOclim_path = artifact_path(ECCOclim_hash)
OCCAclim_path = artifact_path(OCCAclim_hash)

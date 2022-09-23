module ScratchSpace

using Downloads, Scratch

# This will be filled in inside `__init__()`
path = ""

# Downloads a resource, stores it within path
function download_dataset(url)
    fname = joinpath(path, basename(url))
    if !isfile(fname)
        Downloads.download(url, fname)
    end
    return fname
end

function __init__()
    global path = @get_scratch!("downloaded_files")
end

end

module downloads

using OceanStateEstimation: pkg_pth
using OceanStateEstimation: ScratchSpace
using Artifacts, LazyArtifacts, Downloads, Tar, CodecZlib
using Statistics, FortranFiles, MeshArrays, MITgcmTools

##

artifact_toml = joinpath(pkg_pth, "../Artifacts.toml")

MITPROFclim_hash = artifact_hash("MITPROFclim", artifact_toml)
MITPROFclim_path = artifact_path(MITPROFclim_hash)*"/"

"""
    MITPROFclim_download()

Download lazy artifact to `MITPROFclim_path`.
"""   
MITPROFclim_download() = artifact"MITPROFclim"

CBIOMESclim_hash = artifact_hash("CBIOMESclim", artifact_toml)
CBIOMESclim_path = artifact_path(CBIOMESclim_hash)*"/"

"""
    CBIOMESclim_download()

Download lazy artifact to `CBIOMESclim_path`.
"""   
CBIOMESclim_download() = artifact"CBIOMESclim"

ECCOdiags_hash = artifact_hash("ECCOdiags", artifact_toml)
ECCOdiags_path = artifact_path(ECCOdiags_hash)*"/"

"""
    ECCOdiags_download()

Download lazy artifact to `ECCOdiags_path`.
"""    
ECCOdiags_download() = artifact"ECCOdiags"

## Dataverse Donwloads

"""
    get_from_dataverse(lst::String,nam::String,pth::String)

```
using OceanStateEstimation, CSV, DataFrames
pth=dirname(pathof(OceanStateEstimation))
lst=joinpath(pth,"../examples/OCCA_climatology.csv")
nams=CSV.read(lst,DataFrame)
[get_from_dataverse(lst,nam,ScratchSpace.path) for nam in nams.name[:]]
```
"""
function get_from_dataverse(lst::String,nam::String,pth::String)
    lists=dataverse_lists(lst)
    ii = findall([occursin("$nam", lists.name[i]) for i=1:length(lists.ID)])
    for i in ii 
        nam1=ScratchSpace.download_dataset(lists.URL[i])
        if length(ii)>1
            !isdir(joinpath(pth,nam)) ? mkdir(joinpath(pth,nam)) : nothing
            nam2=joinpath(pth,nam,lists.name[i])
            mv(nam1,nam2)
        else
            nam2=joinpath(pth,lists.name[i])
            mv(nam1,nam2)
        end
    end
end

get_from_dataverse(nam::String,pth::String) = get_from_dataverse("../examples/nctiles_climatology.csv",nam,pth)

"""
    dataverse_lists(lst::String)

Read and derive lists (ID,name,URL) from csv file (ID,name) and return as tuple
```
lists=dataverse_lists(lst)
```
"""
function dataverse_lists(lst::String)
    tmp=readlines(lst)
    ID=[parse(Int,tmp[j][1:findfirst(isequal(','),tmp[j])-1]) for j=2:length(tmp)]
    name=[tmp[j][findfirst(isequal(','),tmp[j])+1:end] for j=2:length(tmp)]
    tmp="https://dataverse.harvard.edu/api/access/datafile/"
    URL=[tmp*"$(ID[j])" for j=1:length(ID)]
    return (ID=ID,name=name,URL=URL)
end

"""
    get_ecco_files(γ::gcmgrid,v::String,t=1)

```
using MeshArrays, OceanStateEstimation
γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
tmp=get_ecco_files(γ,"oceQnet")
```
"""
function get_ecco_files(γ::gcmgrid,v::String,t=1)
    get_ecco_variable_if_needed(v)
    return read_nctiles(joinpath(ScratchSpace.path,"$v/$v"),"$v",γ,I=(:,:,t))
end

"""
    get_ecco_variable_if_needed(v::String)

Download ECCO output for variable `v` to scratch space if needed
"""
function get_ecco_variable_if_needed(v::String)
    lst=joinpath(pkg_pth,"../examples/nctiles_climatology.csv")
    !isdir(joinpath(ScratchSpace.path,v)) ? get_from_dataverse(lst,v,ScratchSpace.path) : nothing
end

"""
    get_ecco_velocity_if_needed()

Download ECCO output for `u,v,w` to scratch space if needed
"""
function get_ecco_velocity_if_needed()
    get_ecco_variable_if_needed("UVELMASS")
    get_ecco_variable_if_needed("VVELMASS")
    get_ecco_variable_if_needed("WVELMASS")
end

"""
    get_occa_variable_if_needed(v::String)

Download OCCA output for variable `v` to scratch space if needed
"""
function get_occa_variable_if_needed(v::String)
    lst=joinpath(pkg_pth,"../examples/OCCA_climatology.csv")
    fil=joinpath(ScratchSpace.path,v*".0406clim.nc")
    !isfile(fil) ? get_from_dataverse(lst,v,ScratchSpace.path) : nothing
end

"""
    get_occa_velocity_if_needed()

Download OCCA output for `u,v,w` to scratch space if needed
"""
function get_occa_velocity_if_needed()
    nams = ("DDuvel","DDvvel","DDwvel","DDtheta","DDsalt")
    [get_occa_variable_if_needed(nam) for nam in nams]
    "done"
end

## zenodo.org Downloads

"""
    ECCOdiags_add(nam::String)

Add data to the ECCOdiags_path folder. Known options for `nam` include 
"release1", "release3", "release4", and "interp_coeffs". Note that 
"release2" is the estimate that's readily donwloaded by ECCOdiags_download().
"""
function ECCOdiags_add(nam::String)
    if nam=="release1"
        url="https://zenodo.org/record/6123262/files/ECCOv4r1_analysis.tar.gz?download=1"
        fil="ECCOv4r1_analysis.tar.gz"
    elseif nam=="release3"
        url="https://zenodo.org/record/6123288/files/ECCOv4r3_analysis.tar.gz?download=1"
        fil="ECCOv4r3_analysis.tar.gz"
    elseif nam=="release4"
        url="https://zenodo.org/record/6123127/files/ECCOv4r4_analysis.tar.gz?download=1"
        fil="ECCOv4r4_analysis.tar.gz"
    elseif nam=="interp_coeffs"
        url="https://zenodo.org/record/5784905/files/interp_coeffs_halfdeg.jld2?download=1"
        fil="interp_coeffs_halfdeg.jld2"
    else
        println("unknown release name")
        url=missing
        fil=missing
    end
    if (!ismissing(url))&&(nam=="interp_coeffs")
        Downloads.download(url,joinpath(ECCOdiags_path,fil);timeout=60000.0)
    elseif (!ismissing(url))&&(!isdir(joinpath(ECCOdiags_path,fil)[1:end-7]))
        println("downloading "*nam*" ... started")
        Downloads.download(url,joinpath(ECCOdiags_path,fil);timeout=60000.0)
        tmp_path=open(joinpath(ECCOdiags_path,fil)) do io
            Tar.extract(CodecZlib.GzipDecompressorStream(io))
        end
        mv(joinpath(tmp_path,fil[1:end-7]),joinpath(ECCOdiags_path,fil[1:end-7]))
        rm(joinpath(ECCOdiags_path,fil))
        println("downloading "*nam*" ... done")
    end
end

end

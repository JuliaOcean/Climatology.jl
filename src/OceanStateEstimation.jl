module OceanStateEstimation

using Statistics, Pkg.Artifacts, Downloads
using FortranFiles, MeshArrays, MITgcmTools
export dataverse_lists, get_from_dataverse, get_ecco_files
export get_occa_velocity_if_needed, get_ecco_velocity_if_needed
export get_ecco_variable_if_needed
export ECCOclim_path, OCCAclim_path, MITPROFclim_path

##

p=dirname(pathof(OceanStateEstimation))
artifact_toml = joinpath(p, "../Artifacts.toml")
ECCOclim_hash = artifact_hash("ECCOclim", artifact_toml)
ECCOclim_path = artifact_path(ECCOclim_hash)*"/"
OCCAclim_hash = artifact_hash("OCCAclim", artifact_toml)
OCCAclim_path = artifact_path(OCCAclim_hash)*"/"
MITPROFclim_hash = artifact_hash("MITPROFclim", artifact_toml)
MITPROFclim_path = artifact_path(MITPROFclim_hash)*"/"

##

"""
    get_from_dataverse(lst::String,nam::String,pth::String)

```
using OceanStateEstimation, CSV, DataFrames
lst=joinpath(dirname(pathof(OceanStateEstimation)),"../examples/OCCA_climatology.csv")
nams = CSV.File(lst) |> DataFrame!
nams = nams.name[:]
[get_from_dataverse(lst,nam,OCCAclim_path) for nam in nams]
```
"""
function get_from_dataverse(lst::String,nam::String,pth::String)
    lists=dataverse_lists(lst)
    ii = findall([occursin("$nam", lists.name[i]) for i=1:length(lists.ID)])
    !isdir("$pth"*"$nam") ? mkdir("$pth"*"$nam") : nothing
    for i in ii
        nam1=Downloads.download(lists.URL[i])
        nam2=joinpath("$pth"*"$nam/",lists.name[i])
        mv(nam1,nam2)
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
γ=GridSpec("LatLonCap","./")
tmp=get_ecco_files(γ,"oceQnet")
```
"""
function get_ecco_files(γ::gcmgrid,v::String,t=1)
    pth=artifact_path(ECCOclim_hash)*"/"    
    lst=joinpath(dirname(pathof(OceanStateEstimation)),"../examples/nctiles_climatology.csv")
    !isdir("$pth"*v) ? get_from_dataverse(lst,v,ECCOclim_path) : nothing
    #return read_nctiles(ECCOclim_path*"$v/$v","$v",γ,I=(:,:,:,t))
    return read_nctiles(ECCOclim_path*"$v/$v","$v",γ,I=(:,:,t))
end

"""
    get_ecco_variable_if_needed(v::String)

Download ECCO output for variable `v` to `ECCOclim_path` if needed
"""
function get_ecco_variable_if_needed(v::String)
    p=dirname(pathof(OceanStateEstimation))
    lst=joinpath(p,"../examples/nctiles_climatology.csv")
    pth=ECCOclim_path
    !isdir(pth*v) ? get_from_dataverse(lst,v,pth) : nothing
end

"""
    get_ecco_velocity_if_needed()

Download ECCO output for `u,v,w` to `ECCOclim_path` if needed
"""
function get_ecco_velocity_if_needed()
    get_ecco_variable_if_needed("UVELMASS")
    get_ecco_variable_if_needed("VVELMASS")
    get_ecco_variable_if_needed("WVELMASS")
end

"""
    get_occa_velocity_if_needed()

Download `MITgcm` transport output to `OCCAclim_path` if needed
"""
function get_occa_velocity_if_needed()
    p=dirname(pathof(OceanStateEstimation))
    lst=joinpath(p,"../examples/OCCA_climatology.csv")
    pth=OCCAclim_path
    nams = ("DDuvel.0406clim.nc","DDvvel.0406clim.nc","DDwvel.0406clim.nc","DDtheta.0406clim.nc","DDsalt.0406clim.nc")
    if !isfile("$pth"*"DDuvel.0406clim.nc") 
        tmp=joinpath(pth,"tmp/")
        !isdir(tmp) ? mkdir(tmp) : nothing
        [get_from_dataverse(lst,nam,tmp) for nam in nams]
        [mv(joinpath(tmp,nam,nam),joinpath(pth,nam)) for nam in nams]
    end
end

end # module

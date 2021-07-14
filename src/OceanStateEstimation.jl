module OceanStateEstimation

using Statistics, Pkg.Artifacts, Downloads
using FortranFiles, MeshArrays, MITgcmTools
export dataverse_lists, get_from_dataverse, get_ecco_files
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

end # module

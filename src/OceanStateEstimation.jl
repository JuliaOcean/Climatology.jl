module OceanStateEstimation

using Statistics, Pkg.Artifacts   
using FortranFiles, MeshArrays, MITgcmTools
export dataverse_lists, get_from_dataverse, get_ecco_files

##

p=dirname(pathof(OceanStateEstimation))
artifact_toml = joinpath(p, "../Artifacts.toml")
ECCOclim_hash = artifact_hash("ECCOclim", artifact_toml)
OCCAclim_hash = artifact_hash("OCCAclim", artifact_toml)

fil1=joinpath(p,"../examples/OCCA_climatology.csv")
fil2=joinpath(artifact_path(OCCAclim_hash)*"/", "OCCA_climatology.csv")
!isfile(fil2) ? cp(fil1,fil2) : nothing

fil1=joinpath(p,"../examples/nctiles_climatology.csv")
fil2=joinpath(artifact_path(ECCOclim_hash)*"/", "nctiles_climatology.csv")
!isfile(fil2) ? cp(fil1,fil2) : nothing

##

"""
    get_from_dataverse(lst::String,nam::String,pth::String)

```
using OceanStateEstimation, CSV, DataFrames, Pkg.Artifacts
pth=artifact_path(OceanStateEstimation.OCCAclim_hash)*"/"
lst=joinpath(pth,"OCCA_climatology.csv")
nams = CSV.File(lst) |> DataFrame!
nams = nams.name[:]
[get_from_dataverse(lst,nam,pth) for nam in nams]
```
"""
function get_from_dataverse(lst::String,nam::String,pth::String)
    lists=dataverse_lists(lst)
    ii = findall([occursin("$nam", lists.name[i]) for i=1:length(lists.ID)])
    !isdir("$pth"*"$nam") ? mkdir("$pth"*"$nam") : nothing
    for i in ii
        nam1=download(lists.URL[i]);
        nam2=joinpath("$pth"*"$nam/",lists.name[i])
        run(`mv $nam1 $nam2`);
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
    lst=joinpath(pth,"nctiles_climatology.csv")
    !isdir("$pth"*v) ? get_from_dataverse(lst,v,pth) : nothing
    #return read_nctiles("$pth"*"$v/$v","$v",γ,I=(:,:,:,t))
    return read_nctiles("$pth"*"$v/$v","$v",γ,I=(:,:,t))
end

end # module

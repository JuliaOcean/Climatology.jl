module OceanStateEstimation

using Statistics
using FortranFiles, MeshArrays, MITgcmTools
export get_from_dataverse, get_ecco_files

##

"""
    get_from_dataverse(lst::String,nam::String,pth::String)

```
using OceanStateEstimation, CSV, DataFrames
lst=joinpath(dirname(pathof(OceanStateEstimation)),"../examples/OCCA_climatology.csv")
pth=joinpath(dirname(pathof(OceanStateEstimation)),"../examples/OCCA_climatology/")
nams = CSV.File(lst) |> DataFrame!
nams = nams.name[:]
[get_from_dataverse(lst,nam,pth) for nam in nams]
```
"""
function get_from_dataverse(lst::String,nam::String,pth::String)
    tmp=readlines(lst)
    ID=[parse(Int,tmp[j][1:findfirst(isequal(','),tmp[j])-1]) for j=2:length(tmp)]
    name=[tmp[j][findfirst(isequal(','),tmp[j])+1:end] for j=2:length(tmp)]
    ii = findall([occursin("$nam", name[i]) for i=1:length(ID)])
    !isdir("$pth"*"$nam") ? mkdir("$pth"*"$nam") : nothing
    for i in ii
        id1=ID[i]
        nam1=name[i]
        nam2=joinpath("$pth"*"$nam/",nam1)
        run(`wget --content-disposition https://dataverse.harvard.edu/api/access/datafile/$id1`);
        run(`mv $nam1 $nam2`);
    end
end

get_from_dataverse(nam::String,pth::String) = get_from_dataverse("../examples/nctiles_climatology.csv",nam,pth)

"""
    get_ecco_files(γ::gcmgrid,v::String,t=1)

```
using MeshArrays, OceanStateEstimation
γ=GridSpec("LatLonCap","./")
tmp=get_ecco_files(γ,"oceQnet")
```
"""
function get_ecco_files(γ::gcmgrid,v::String,t=1)
    lst=joinpath(dirname(pathof(OceanStateEstimation)),"../examples/nctiles_climatology.csv")
    pth=joinpath(dirname(pathof(OceanStateEstimation)),"../examples/nctiles_climatology/")
    !isdir("$pth") ? mkdir("$pth") : nothing
    !isdir("$pth"*v) ? get_from_dataverse(lst,v,pth) : nothing
    #return read_nctiles("$pth"*"$v/$v","$v",γ,I=(:,:,:,t))
    return read_nctiles("$pth"*"$v/$v","$v",γ,I=(:,:,t))
end

end # module

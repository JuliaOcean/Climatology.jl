module OceanStateEstimation

using CSV, DataFrames, Statistics
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
    tmp = CSV.File(lst) |> DataFrame!
    ii = findall([occursin("$nam", tmp[i,:name]) for i=1:size(tmp,1)])
    !isdir("$pth"*"$nam") ? mkdir("$pth"*"$nam") : nothing
    for i in ii
        ID=tmp[i,:ID]
        nam1=tmp[i,:name]
        nam2=joinpath("$pth"*"$nam/",nam1)
        run(`wget --content-disposition https://dataverse.harvard.edu/api/access/datafile/$ID`);
        run(`mv $nam1 $nam2`);
    end
end

get_from_dataverse(nam::String,pth::String) = get_from_dataverse("../examples/nctiles_climatology.csv",nam,pth)

"""
    get_ecco_files(v="oceQnet",t)

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

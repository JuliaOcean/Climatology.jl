module OceanStateEstimation

using CSV, DataFrames, Statistics
using FortranFiles, MeshArrays, MITgcmTools
export get_grid_if_needed, get_ecco_files, get_from_dataverse


##

"""
    get_ecco_files(v="oceQnet",t)

```
using MeshArrays, OceanStateEstimation
get_grid_if_needed()
pth="../examples/GRID_LLC90/"
γ=GridSpec("LatLonCap",pth)
Γ=GridLoad(γ)
tmp=get_ecco_files(γ)
```
"""
function get_ecco_files(γ::gcmgrid,v="oceQnet",t=1)
    pth="../examples/nctiles_climatology/"
    !isdir("$pth") ? mkdir("$pth") : nothing
    !isdir("$pth"*v) ? get_from_dataverse(v,pth) : nothing
    #return read_nctiles("$pth"*"$v/$v","$v",γ,I=(:,:,:,t))
    return read_nctiles("$pth"*"$v/$v","$v",γ,I=(:,:,t))
end

function get_from_dataverse(nam::String,pth::String)
    tmp = CSV.File("nctiles_climatology.csv") |> DataFrame!
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

function get_grid_if_needed()
    if !isdir("../examples/GRID_LLC90")
        run(`git clone https://github.com/gaelforget/GRID_LLC90 ../examples/GRID_LLC90`)
    end
end

end # module

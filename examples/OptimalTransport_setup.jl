
"""
object: setup to compute optimal transport between model and/or satellite climatologies
date: 2021/10/28
author: Gaël Forget

- examples/CBIOMES_climatology_compare.jl
"""

import OceanStateEstimation, NCTiles
using Statistics, LinearAlgebra, JLD2

## load files

fil_out=joinpath(OceanStateEstimation.CBIOMESclim_path,"CBIOMES-global-alpha-climatology.nc")
nc=NCTiles.NCDataset(fil_out,"r")
lon=nc["lon"][:]
lat=nc["lat"][:]
uni=nc["Chl"].attrib["units"]

## region and base distance (Cost) definition

i1=findall( (lon.>-180.0).*(lon.<-120.0) )
j1=findall( (lat.>-20.0).*(lat.<50.0) )

## main arrays
Chl_from_Mod=nc["Chl"][i1,j1,:]
fil_sat="examples_climatology_prep/gridded_geospatial_montly_clim_360_720_ver_0_2.nc"
Chl_from_Sat=NCTiles.NCDataset(fil_sat,"r")["Chl"][i1,j1,:]

## cost matrix

if !isfile("examples_EMD_paper_exploration/example_Cost.jld2")
    #this only needs to be done one
    #C = [[i,j] for i in i1, j in j1]
    C = [[lon[i],lat[j]] for i in i1, j in j1]
    C=C[:]

    gcdist(lo1,lo2,la1,la2) = acos(sind(la1)*sind(la2)+cosd(la1)*cosd(la2)*cosd(lo1-lo2))

    #C=[gcdist(C[i][1],C[j][1],C[i][2],C[j][2]) for i in 1:length(C), j in 1:length(C)]

    nx=length(C)
    Cost=zeros(nx,nx)
    for i in 1:length(C), j in 1:length(C)
        i!==j ? Cost[i,j]=gcdist(C[i][1],C[j][1],C[i][2],C[j][2]) : nothing
    end

    @save "examples_EMD_paper_exploration/example_Cost.jld2" Cost
end
Cost=load("examples_EMD_paper_exploration/example_Cost.jld2")["Cost"]
println("reusing Cost matrix computed previously\n")

## helper functions

function preprocess_Chl(a,b)
    k=findall(ismissing.(a).|ismissing.(b));
    a[k].=0.0; b[k].=0.0;
    k=findall((a.<0).|(b.<0));
    a[k].=0.0; b[k].=0.0;
    k=findall(isnan.(a).|isnan.(b));
    a[k].=0.0; b[k].=0.0;

    M=0.1
    k=findall((a.>M).|(b.>M));
    a[findall(a.>M)].=M;
    b[findall(b.>M)].=M;

    a=Float64.(a); a=a/sum(a)
    b=Float64.(b); b=b/sum(b)

    a,b
end

## 

function export_zm()
    M=NaN*zeros(140,12)
    S=NaN*zeros(140,12)
    for t in 1:12
        a=Chl_from_Mod[:,:,t][:]
        b=Chl_from_Sat[:,:,t][:]
        a,b=preprocess_Chl(a,b)
        M[:,t]=sum(reshape(a,(120,140)),dims=1)[:]
        S[:,t]=sum(reshape(b,(120,140)),dims=1)[:]
    end
    (M,S)
end


"""
object: compute optimal transport between model and/or satellite climatologies
date: 2021/10/28
author: Gaël Forget

- examples/CBIOMES_climatology_compare.jl
"""

using OptimalTransport, Statistics, LinearAlgebra, JLD2

## load files

fil_out=joinpath(CBIOMESclim_path,"CBIOMES-global-alpha-climatology.nc")
nc=NCTiles.NCDataset(fil_out,"r")
lon=nc["lon"][:]
lat=nc["lat"][:]
uni=nc["Chl"].attrib["units"]
Chl_from_Mod=nc["Chl"][:]

fil_sat="examples/gridded_geospatial_montly_clim_360_720_ver_0_2.nc"
Chl_from_Sat=NCTiles.NCDataset(fil_sat,"r")["Chl"][:]

## region and base distance (Cost) definition

i1=findall( (lon.>-180.0).*(lon.<-120.0) )
j1=findall( (lat.>-20.0).*(lat.<50.0) )

if !isfile("examples/example_Cost.jld2")
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

    @save "examples/example_Cost.jld2" Cost
end
Cost=load("examples/example_Cost.jld2")["Cost"]
println("reusing Cost matrix computed previously\n")

## helper functions

function preprocess_Chl(a,b)
    k=findall(xor.(ismissing.(a),ismissing.(b)));
    a[k].=0.0; b[k].=0.0;
    k=findall(xor.(a.<0,b.<0));
    a[k].=0.0; b[k].=0.0;
    k=findall(xor.(isnan.(a),isnan.(b)));
    a[k].=0.0; b[k].=0.0;

    M=0.1
    k=findall(xor.(a.>M,b.>M));
    a[findall(a.>M)].=M;
    b[findall(b.>M)].=M;

    a=Float64.(a); a=a/sum(a)
    b=Float64.(b); b=b/sum(b)

    a,b
end

## 

a=mean(Chl_from_Sat[i1,j1,:],dims=3)[:]
b=mean(Chl_from_Mod[i1,j1,:],dims=3)[:]
a,b=preprocess_Chl(a,b)

ε=0.1
G=sinkhorn(a,b, Cost, ε)
#
bb=G'*a #predict b, step1
bb=bb[:]*sum(b)/sum(bb) #predict b, step2
#
OptCost=sinkhorn2(a,b, Cost, ε, plan=G) #compute optimal cost
#
cc=dot(G, Cost) #compute optimal cost, directly

#@time sinkhorn2(a,b, Cost, ε);
#@time sinkhorn(a,b, Cost, ε);
#@time sinkhorn_divergence(a,b, Cost, ε);

## loop over months

if false
d=zeros(12,12)
t=[time()]
for i in 1:12, j in 1:12
    if j!==i
        a=Chl_from_Mod[i1,j1,i][:];
        b=Chl_from_Sat[i1,j1,j][:];
        a,b=preprocess_Chl(a,b)

        ε = 0.1
        d[i,j] = sinkhorn2(a,b, Cost, ε)

        #timing
        δt=time()-t[1]
        t[1]=time()
        println("$i,$j,$(δt)")
    end
end
end
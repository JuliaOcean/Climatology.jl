import Climatology: read_IAP, file_IAP, write_H_to_T

import NCDatasets, MeshArrays
using DataStructures: OrderedDict
using MeshArrays: gridmask, Integration

"""
    file_IAP(path,y,m)
"""
file_IAP(path,y,m)=begin
    mm=(m<10 ? "0$m" : "$m")        
    joinpath(path,"IAPv4_Temp_monthly_1_6000m_year_$(y)_month_$(mm).nc")
end

"""
    read_IAP(F,var,tim,tmp=[])

```
using Climatology, NCDatasets, MeshArrays

p0="IAPv4_IAP_Temperature_gridded_1month_netcdf/monthly/"
fil=Climatology.file_IAP(p0,"2023","12")
depth=Dataset(fil)["depth_std"][:]
temp=Climatology.read_IAP(fil,"temp",1,[])
mask=1.0*(!ismissing).(temp)

G=Gris_simple.GridLoad_lonlatdep(depth,mask)
tmp=zeros(G.XC.grid)*ones(length(depth))
Climatology.read_IAP(fil,"temp",1,tmp)
```
"""
function read_IAP(F,var,tim,tmp=[])
    fil=F

    ds=Dataset(fil)
    temp=permutedims(ds[var][:,:,:],(2,3,1))
    close(ds)
    temp[findall(ismissing.(temp))].=0
    temp[findall(isnan.(temp))].=0
    if !isempty(tmp)
        tmp.=read(Float32.(temp),tmp)
        tmp
    else
        temp
    end
end

"""
    write_H_to_T(file::String,M::gridmask,G::NamedTuple,H::Array)

Write `H / Integration.volumes(M,G)` to file.

```
using Climatology, NCDatasets, MeshArrays

G=MeshArrays.Grids_simple.GridLoad_lonlatdep(depth,mask)
M=Integration.define_sums(grid=G,regions=(10,5))
H=ones(length(M.names),length(M.depths),3)
V=Integration.volumes(M,G)
Climatology.write_H_to_T(tempname()*".nc",M,G,H,V)
```
"""
function write_H_to_T(file::String,M::gridmask,G::NamedTuple,H::Array,V::Array)
    nb,nz,nt=size(H)
    inv_vol=1.0./V
    #inv_vol[V.==0].=0
    pos=gridpos(M,(10,5))
    arr2d=zeros(36,32)
    arr3d=zeros(36,32,nz)
    arr4d=zeros(36,32,nz,nt)
    println(nz)

    ds = Dataset(file,"c")
    defDim(ds,"lon",36); defDim(ds,"lat",32); 
    defDim(ds,"dep",size(H,2)); defDim(ds,"tim",size(H,3));
    ds.attrib["title"] = "this is a test file"

    dlo=10; dla=5;
    lons=collect(-180:dlo:180); lons=0.5*(lons[1:end-1]+lons[2:end])
    lats=[-90 ; -75:dla:75 ; 90]; lats=0.5*(lats[1:end-1]+lats[2:end])
    vlo=defVar(ds,"lon",lons,("lon",), attrib = OrderedDict("units" => "degree", "long_name" => "Longitude"))
    vla=defVar(ds,"lat",lats,("lat",), attrib = OrderedDict("units" => "degree", "long_name" => "Latitude"))

    v1 = defVar(ds,"volume",Float32,("lon","lat","dep"), attrib = OrderedDict("units" => "m^3",))

    arr3d.=0 
    for ii in 1:nb
    i,j=pos[ii]
    [arr3d[i,j,k]=V[ii,k] for k in 1:nz]
    end

    v1[:,:,:] = arr3d

    v = defVar(ds,"temperature",Float32,("lon","lat","dep","tim"), attrib = OrderedDict("units" => "degree Celsius",))
    v.attrib["comments"] = "this is a string attribute with Unicode Ω ∈ ∑ ∫ f(x) dx"

    arr4d.=0
    for t in 1:nt
    for ii in 1:nb
    i,j=pos[ii]
    [arr4d[i,j,k,t]=H[ii,k,t]*inv_vol[ii,k] for k in 1:nz]
    end
    end

    v[:,:,:,:] = arr4d

    close(ds)
    file
end

"""
    gridpos(M::gridmask,res::Tuple)

```
gridpos(M,(10,5))
```
"""
gridpos(M::gridmask,res::Tuple)=begin
    n=length(M.names)
    allpos=fill((0,0),n)
    for i in 1:n
        t1=split(M.names[i],"Nto")
        t2=split(t1[2],"N_")
        t3=split(t2[2],"Eto")
        t4=split(t3[2],"E")
        tt=[t1[1] t2[1] t3[1] t4[1]]
        tt=parse.(Ref(Int),tt)

        dlo=res[1]; dla=res[2]
        lons=collect(-180:dlo:180)
        lats=[-90 ; -75:dla:75 ; 90]
        thispos=(findall(lons.==tt[3])[1],findall(lats.==tt[1])[1])

        allpos[i]=thispos
    end
    allpos
end

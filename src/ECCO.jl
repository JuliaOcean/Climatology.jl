
module ECCO

using Pkg, DataFrames
import Climatology: pkg_pth

"""
    ECCO.standard_analysis_setup(pth0::String)

Create temporary run folder `pth` where data folder `pth0` will be linked. 

Data folder `pth0` should be the path to ECCO data.

For example:

```
using Climatology, Pkg
pth=ECCO.standard_analysis_setup(ScratchSpaces.ECCO)
```

The `Project.toml` file found in `pth` provides an environment ready for `ECCO` analyses. 

This environment can be activated and instantiated:

```
Pkg.activate(pth)
Pkg.instantiate()
```
"""
function standard_analysis_setup(pth0="",sol0="")
	
	#1. setup run folder and create link to ECCO data folder
	pth=joinpath(tempdir(),"ECCO_diags_dev"); 
	!isdir(pth) ? mkdir(pth) : nothing
    if in(sol0,["r2","r3","r4","r5"])
        pth1=joinpath(pth,"ECCOv4"*sol0)
    else
        pth1=joinpath(pth,sol0)
    end

    !isdir(pth1) ? mkdir(pth1) : nothing
	link0=joinpath(pth1,"diags")
	!isfile(link0)&& !islink(link0)&& !isempty(pth0) ? symlink(pth0,link0) : nothing
	
	#2. copy Project.toml to run folder
	tmp0=pkg_pth
	tmp1=joinpath(tmp0,"..","examples","ECCO","ECCO_standard_Project.toml")
	tmp2=joinpath(pth,"Project.toml")
	!isfile(tmp2) ? cp(tmp1,tmp2) : nothing
		
	return pth1
end


add_diag!(list,file=tempname(),name="variable",units="unknown",dims=("time",)) = begin
    append!(list,DataFrame("file"=>file,"name"=>name,"units"=>units,"dims"=>dims))
end

#time series
function diagnostics_set1(path_in=".")
    list=DataFrame("file"=>String[],"name"=>String[],"units"=>String[],"dims"=>Tuple[])
    add_diag!(list,joinpath(path_in,"THETA_glo3d","glo3d.jld2"),"temperature_global","degreeC",("time",))
    add_diag!(list,joinpath(path_in,"THETA_glo2d","glo2d.jld2"),"temperature_global_level","degreeC",("depth","time"))
    add_diag!(list,joinpath(path_in,"SALT_glo3d","glo3d.jld2"),"salinity_global","PSS",("time",))
    add_diag!(list,joinpath(path_in,"SALT_glo2d","glo2d.jld2"),"salinity_global_level","PSS",("depth","time"))
    add_diag!(list,joinpath(path_in,"trsp","trsp.jld2"),"volume_transport","m3/s",("section","depth","time"))
    add_diag!(list,joinpath(path_in,"MHT","MHT.jld2"),"meridional_heat_transport","PW",("latMT","time"))
    add_diag!(list,joinpath(path_in,"THETA_zonmean","zonmean.jld2"),"temperature_zonal_level","degreeC",("latZM","depth","time"))
    add_diag!(list,joinpath(path_in,"SALT_zonmean","zonmean.jld2"),"salinity_zonal_level","PSS",("latZM","depth","time"))
    add_diag!(list,joinpath(path_in,"MXLDEPTH_zonmean2d","zonmean2d.jld2"),"MLD_zonal","m",("latZM","time"))
    add_diag!(list,joinpath(path_in,"SSH_zonmean2d","zonmean2d.jld2"),"SSH_zonal","m",("latZM","time"))
    add_diag!(list,joinpath(path_in,"SIarea_zonmean2d","zonmean2d.jld2"),"SIarea_zonal","nondimensional",("latZM","time"))
    add_diag!(list,joinpath(path_in,"overturn","overturn.jld2"),"overturn","m3/s",("latMT","depth","time"))
    list
end

#2d climatologies on ECCO's LLC90 grid
function diagnostics_set2(path_in=".")
    list=DataFrame("file"=>String[],"name"=>String[],"units"=>String[],"dims"=>Tuple[])
    add_diag!(list,joinpath(path_in,"BSF_clim","BSF.jld2"),"BSF_clim","m3/s",("time",))
    add_diag!(list,joinpath(path_in,"MXLDEPTH_clim","MXLDEPTH.jld2"),"MXLDEPTH_clim","m",("time",))
    add_diag!(list,joinpath(path_in,"SIarea_clim","SIarea.jld2"),"SIarea_clim","nondimensional",("time",))
    add_diag!(list,joinpath(path_in,"SSH_clim","SSH.jld2"),"SSH_clim","m",("time",))
    list
end

#3d climatologies on ECCO's LLC90 grid
function diagnostics_set3(path_in=".")
    list=DataFrame("file"=>String[],"name"=>String[],"units"=>String[],"dims"=>Tuple[])
    add_diag!(list,joinpath(path_in,"THETA_clim","THETA_k01.jld2"),"THETA_clim","degreeC",("time",))
    add_diag!(list,joinpath(path_in,"SALT_clim","SALT_k01.jld2"),"SALT_clim","PSS",("time",))
    list
end

end

##

module ECCO_helpers

using MeshArrays, TOML, JLD2, Glob
import Climatology: read_Dataset

"""
    parameters(P0,params)

Prepare parameter NamedTuple for use in `ECCO_diagnostics.driver`.

`P1=parameters(P0,p)` 

is faster than e.g. `parameters(pth,"r2",p)` as grid, etc get copied from `P0` to `P1`.
"""
function parameters(P,params)
    calc=params.calc
    nam=params.nam
    kk=params.lev

    pth_out=dirname(P.pth_out)
    if sum(calc.==("overturn","MHT","trsp"))==0
        pth_out=joinpath(pth_out,nam*"_"*calc)
    else
        pth_out=joinpath(pth_out,calc)
    end    

    return (pth_in=P.pth_in,pth_out=pth_out,list_steps=P.list_steps,nt=P.nt,
    calc=calc,nam=nam,kk=kk,sol=P.sol,γ=P.γ,Γ=P.Γ,LC=P.LC)
end

"""
    parameters(pth0::String,sol0::String,params)

Prepare parameter NamedTuple for use in `ECCO_diagnostics.driver`.

For example, to compute zonal mean temperatures at level 5:

```
p=(calc = "zonmean", nam = "THETA", lev = 5)
pth=ECCO.standard_analysis_setup(ScratchSpaces.ECCO)
P0=ECCO_helpers.parameters(pth,"r2",p)
```

or, from a predefined list:

```
list0=ECCO_helpers.standard_list_toml("")
pth=ECCO.standard_analysis_setup(ScratchSpaces.ECCO)
P1=ECCO_helpers.parameters(pth,"r2",list0[1])
```
"""
function parameters(pth0::String,sol0::String,params)

    calc=params.calc
    nam=params.nam
    kk=params.lev
    if in(sol0,["r2","r3","r4","r5"])
      sol="ECCOv4"*sol0*"_analysis"
      pth_in=joinpath(pth0,"ECCOv4"*sol0,"diags")
    else
      sol=sol0*"_analysis"
      pth_in=joinpath(pth0,sol0,"diags")
    end

    !ispath(pth_in) ? pth_in=joinpath(pth0,"diags") : nothing
    list_steps=list_time_steps(pth_in)

    if sol0=="r1"||sol0=="r2"
        fil=joinpath(pth_in,"THETA","THETA.0001.nc")
        if isfile(fil)
            nt=read_Dataset(fil) do ds
                data = length(ds["tim"][:])
            end
        else
            nt=12
        end
    elseif sol0=="r3"
        nt=288
    elseif sol0=="r4"
        nt=312
    elseif sol0=="r5"
        nt=336
    else
        nt=length(list_steps)
    end

    pth_out=joinpath(pth0,sol)

    if sum(calc.==("overturn","MHT","trsp"))==0
        pth_out=joinpath(pth_out,nam*"_"*calc)
    else
        pth_out=joinpath(pth_out,calc)
    end    

    γ,Γ,LC=GridLoad_Plus()

    P=(pth_in=pth_in,pth_out=pth_out,list_steps=list_steps,nt=nt,
    calc=calc,nam=nam,kk=kk,sol=sol,γ=γ,Γ=Γ,LC=LC)
end

#STATE/state_3d_set1.0000241020.meta
#    'THETA   ' 'SALT    ' 'DRHODR  '
#TRSP/trsp_3d_set1.0000241020.meta
#    'UVELMASS' 'VVELMASS' 'WVELMASS' 'GM_PsiX ' 'GM_PsiY '
#TRSP/trsp_3d_set3.0000241020.meta
#    'DFxE_TH ' 'DFyE_TH ' 'ADVx_TH ' 'ADVy_TH ' 'DFxE_SLT' 'DFyE_SLT' 'ADVx_SLT' 'ADVy_SLT'

function list_time_steps(pth_in)
    println(pth_in)
    if !isempty(glob("STATE/state_3d_set1*.data",pth_in))
        list=basename.(glob("STATE/state_3d_set1*.data",pth_in)) 
    elseif !isempty(glob("state_3d_set1*.data",pth_in))
        list=basename.(glob("state_3d_set1*.data",pth_in))
    else
        list=[]
    end
    return list
end

nansum(x) = sum(filter(!isnan,x))
nansum(x,y) = mapslices(nansum,x,dims=y)

function GridLoad_Plus()
    G=GridLoad(ID=:LLC90,option=:light)
    γ=G.XC.grid
    nr=length(G.RC)

    hFacC=GridLoadVar("hFacC",γ)
    hFacW=GridLoadVar("hFacW",γ)
    hFacS=GridLoadVar("hFacS",γ)
    mskC=hFacC./hFacC

    tmp=[nansum(mskC[i,j].*G.RAC[i]) for j in 1:nr, i in eachindex(G.RAC)]
    tot_RAC=nansum(tmp,2)
    tmp=[nansum(hFacC[i,j].*G.RAC[i].*G.DRF[j]) for j in 1:nr, i in eachindex(G.RAC)]
    tot_VOL=nansum(tmp,2)

    G=merge(G,(hFacC=hFacC,hFacW=hFacW,hFacS=hFacS,mskC=mskC,tot_RAC=tot_RAC,tot_VOL=tot_VOL))
	LC=LatitudeCircles(-89.0:89.0,G)
    return γ,G,LC
end

import Base:push!

function push!(allcalc::Vector{String},allnam::Vector{String},allkk::Vector{Int};
    calc="unknown",nam="unknown",kk=1)
    push!(allcalc,calc)
    push!(allnam,nam)
    push!(allkk,kk)
end

function standard_list_toml(fil)
    
    allcalc=String[]
    allnam=String[]
    allkk=Int[]

    push!(allcalc,allnam,allkk;calc="trsp")
    push!(allcalc,allnam,allkk;calc="MHT")
    push!(allcalc,allnam,allkk;calc="zonmean2d",nam="SIarea")
    push!(allcalc,allnam,allkk;calc="zonmean2d",nam="MXLDEPTH")
    push!(allcalc,allnam,allkk;calc="zonmean2d",nam="SSH")
    push!(allcalc,allnam,allkk;calc="zonmean",nam="THETA")
    push!(allcalc,allnam,allkk;calc="glo2d",nam="THETA")
    push!(allcalc,allnam,allkk;calc="glo3d",nam="THETA")
    push!(allcalc,allnam,allkk;calc="zonmean",nam="SALT")
    push!(allcalc,allnam,allkk;calc="glo2d",nam="SALT")
    push!(allcalc,allnam,allkk;calc="glo3d",nam="SALT")
    push!(allcalc,allnam,allkk;calc="overturn")
    [push!(allcalc,allnam,allkk;calc="clim",nam="THETA",kk=kk) for kk in [1 10 20 29 38 44]]
    [push!(allcalc,allnam,allkk;calc="clim",nam="SALT",kk=kk) for kk in [1 10 20 29 38 44]]
    push!(allcalc,allnam,allkk;calc="clim",nam="SSH")
    push!(allcalc,allnam,allkk;calc="clim",nam="MXLDEPTH")
    push!(allcalc,allnam,allkk;calc="clim",nam="SIarea")
    push!(allcalc,allnam,allkk;calc="clim",nam="BSF")

    tmp1=Dict("calc"=>allcalc,"nam"=>allnam,"kk"=>allkk)
    if !isempty(fil)
        open(fil, "w") do io
            TOML.print(io, tmp1)
        end
    end

    out=[(calc=allcalc[i],nam=allnam[i],lev=allkk[i]) for i in 1:length(allcalc)]

    return out
end

##

function transport_lines()
    lonPairs=[]    
    latPairs=[]    
    namPairs=[]    

    push!(lonPairs,[-173 -164]); push!(latPairs,[65.5 65.5]); push!(namPairs,"Bering Strait");
    push!(lonPairs,[-5 -5]); push!(latPairs,[34 40]); push!(namPairs,"Gibraltar");
    push!(lonPairs,[-81 -77]); push!(latPairs,[28 26]); push!(namPairs,"Florida Strait");
    push!(lonPairs,[-81 -79]); push!(latPairs,[28 22]); push!(namPairs,"Florida Strait W1");
    push!(lonPairs,[-76 -76]); push!(latPairs,[21 8]); push!(namPairs,"Florida Strait S1");
    push!(lonPairs,[-77 -77]); push!(latPairs,[26 24]); push!(namPairs,"Florida Strait E1");
    push!(lonPairs,[-77 -77]); push!(latPairs,[24 22]); push!(namPairs,"Florida Strait E2");
    push!(lonPairs,[-65 -50]); push!(latPairs,[66 66]); push!(namPairs,"Davis Strait");
    push!(lonPairs,[-35 -20]); push!(latPairs,[67 65]); push!(namPairs,"Denmark Strait");
    push!(lonPairs,[-16 -7]); push!(latPairs,[65 62.5]); push!(namPairs,"Iceland Faroe");
    push!(lonPairs,[-6.5 -4]); push!(latPairs,[62.5 57]); push!(namPairs,"Faroe Scotland");
    push!(lonPairs,[-4 8]); push!(latPairs,[57 62]); push!(namPairs,"Scotland Norway");
    push!(lonPairs,[-68 -63]); push!(latPairs,[-54 -66]); push!(namPairs,"Drake Passage");
    push!(lonPairs,[103 103]); push!(latPairs,[4 -1]); push!(namPairs,"Indonesia W1");
    push!(lonPairs,[104 109]); push!(latPairs,[-3 -8]); push!(namPairs,"Indonesia W2");
    push!(lonPairs,[113 118]); push!(latPairs,[-8.5 -8.5]); push!(namPairs,"Indonesia W3");
    push!(lonPairs,[118 127 ]); push!(latPairs,[-8.5 -15]); push!(namPairs,"Indonesia W4");
    push!(lonPairs,[127 127]); push!(latPairs,[-25 -68]); push!(namPairs,"Australia Antarctica");
    push!(lonPairs,[38 46]); push!(latPairs,[-10 -22]); push!(namPairs,"Madagascar Channel");
    push!(lonPairs,[46 46]); push!(latPairs,[-22 -69]); push!(namPairs,"Madagascar Antarctica");
    push!(lonPairs,[20 20]); push!(latPairs,[-30 -69.5]); push!(namPairs,"South Africa Antarctica");
    push!(lonPairs,[-76 -72]); push!(latPairs,[21 18.5]); push!(namPairs,"Florida Strait E3");
    push!(lonPairs,[-72 -72]); push!(latPairs,[18.5 10]); push!(namPairs,"Florida Strait E4");

    lonPairs,latPairs,namPairs
end

function transport_lines(Γ,pth_trsp)
    mkdir(pth_trsp)
    lonPairs,latPairs,namPairs=transport_lines()
    for ii in 1:length(lonPairs)
        lons=Float64.(lonPairs[ii])
        lats=Float64.(latPairs[ii])
        name=namPairs[ii]
        Trsct=Transect(name,lons,lats,Γ,format=:NamedTuple)
        jldsave(joinpath(pth_trsp,"$(Trsct.name).jld2"),
            tabC=Trsct.tabC,tabW=Trsct.tabW,tabS=Trsct.tabS); 
    end
    return true
end

function reload_transport_lines(pth_trsp)
    list_trsp=readdir(pth_trsp)
    ntr=length(list_trsp)
    TR=[load(joinpath(pth_trsp,list_trsp[itr])) for itr in 1:ntr]
    return list_trsp,MeshArrays.Dict_to_NamedTuple.(TR),ntr
end

end #module ECCO_helpers

## generic read function

module ECCO_io

using MeshArrays
import Climatology: read_nctiles_alias, read_Dataset, read_mdsio_alias

"""
    read_monthly(P,nam,t)

Read record `t` for variable `nam` from file locations specified via parameters `P`.

The method used to read `nam` is selected based on `nam`'s value. Methods include:

- `read_monthly_default`
- `read_monthly_SSH`
- `read_monthly_MHT`
- `read_monthly_BSF`
"""
function read_monthly(P,nam,t) 
    if nam=="SSH"
        read_monthly_SSH(P,t)
    elseif nam=="MHT"
        read_monthly_MHT(P,t)
    elseif nam=="BSF"
        read_monthly_BSF(P,t)
    else
        read_monthly_default(P,nam,t)
    end
end

function read_monthly_SSH(P,t)
    (; Γ) = P
    ETAN=read_monthly_default(P,"ETAN",t)
    sIceLoad=read_monthly_default(P,"sIceLoad",t)
    (ETAN+sIceLoad/1029.0)*Γ.mskC[:,1]
end

function read_monthly_MHT(P,t)
    (; Γ) = P

    U=read_monthly_default(P,"ADVx_TH",t)
    V=read_monthly_default(P,"ADVy_TH",t)
    U=U+read_monthly_default(P,"DFxE_TH",t)
    V=V+read_monthly_default(P,"DFyE_TH",t)

    [U[i][findall(isnan.(U[i]))].=0.0 for i in eachindex(U)]
    [V[i][findall(isnan.(V[i]))].=0.0 for i in eachindex(V)]

    Tx=0.0*U[:,1]
    Ty=0.0*V[:,1]
    [Tx=Tx+U[:,z] for z=1:nr]
    [Ty=Ty+V[:,z] for z=1:nr]

    return Tx,Ty
end

function read_monthly_BSF(P,t)
    (; Γ) = P

    U=read_monthly_default(P,"UVELMASS",t)
    V=read_monthly_default(P,"VVELMASS",t)
    MeshArrays.UVtoTransport!(U,V,Γ)
    
    nz=size(Γ.hFacC,2)
    μ=Γ.mskC[:,1]
    Tx=0.0*U[:,1]
	Ty=0.0*V[:,1]
	for z=1:nz
		Tx=Tx+U[:,z]
		Ty=Ty+V[:,z]
	end

    #convergence & land mask
    TrspCon=μ.*convergence(Tx,Ty)

    #scalar potential
    TrspPot=ScalarPotential(TrspCon)

    #Divergent transport component
    (TxD,TyD)=gradient(TrspPot,Γ)
    TxD=TxD.*Γ.DXC
    TyD=TyD.*Γ.DYC

    #Rotational transport component
    TxR = Tx-TxD
    TyR = Ty-TyD

    #vector Potential
    TrspPsi=VectorPotential(TxR,TyR,Γ)

    GC.gc()

    return TrspPsi
end

function read_monthly_default(P,nam,t)
    (; pth_in, sol, list_steps, γ) = P

    var_list3d=("THETA","SALT","UVELMASS","VVELMASS",
                "ADVx_TH","ADVy_TH","DFxE_TH","DFyE_TH")
    if ispath(joinpath(pth_in,"STATE"))
      mdsio_list3d=("STATE/state_3d_set1","STATE/state_3d_set1",
        "TRSP/trsp_3d_set1","TRSP/trsp_3d_set1","TRSP/trsp_3d_set2",
        "TRSP/trsp_3d_set2","TRSP/trsp_3d_set2","TRSP/trsp_3d_set2")
    else
      mdsio_list3d=("state_3d_set1","state_3d_set1",
        "trsp_3d_set1","trsp_3d_set1","trsp_3d_set2",
        "trsp_3d_set2","trsp_3d_set2","trsp_3d_set2")
    end

    var_list2d=("MXLDEPTH","SIarea","sIceLoad","ETAN")
    if ispath(joinpath(pth_in,"STATE"))
      mdsio_list2d=("STATE/state_2d_set1","STATE/state_2d_set1",
        "STATE/state_2d_set1","STATE/state_2d_set1")
    else
      mdsio_list2d=("state_2d_set1","state_2d_set1","state_2d_set1","state_2d_set1")
    end

    if (sol=="ECCOv4r1_analysis")||(sol=="ECCOv4r2_analysis")||(sol=="ECCOv4r3_analysis")
        nct_path=joinpath(pth_in,nam)
        try
            if sum(var_list3d.==nam)==1
                tmp=read_nctiles_alias(nct_path,nam,γ,I=(:,:,:,t))
            else
                tmp=read_nctiles_alias(nct_path,nam,γ,I=(:,:,t))
            end
        catch
            error("failed: call to `read_nctiles`
            This method is provided by `MITgcm.jl`
            and now activated by `using MITgcm` ")
        end
    elseif (sol=="ECCOv4r4_analysis")
        y0=Int(floor((t-1)/12))+1992
        m0=mod1(t,12)
        nct_path=joinpath(pth_in,nam,string(y0))
        m0<10 ? fil=nam*"_$(y0)_0$(m0).nc" : fil=nam*"_$(y0)_$(m0).nc"
        tmp0=read_Dataset(joinpath(nct_path,fil))[nam]
        til0=Tiles(γ,90,90)
        if sum(var_list3d.==nam)==1
            tmp=MeshArray(γ,γ.ioPrec,nr)
            for i in 1:13, k in 1:50
              ff=til0[i].face
              ii=collect(til0[i].i)
              jj=collect(til0[i].j)
              tmp[ff,k][ii,jj]=tmp0[:,:,i,k,1]
            end
            tmp
        else
            tmp=MeshArray(γ,γ.ioPrec)
            for i in 1:13
              ff=til0[i].face
              ii=collect(til0[i].i)
              jj=collect(til0[i].j)
              tmp[ff][ii,jj]=tmp0[:,:,i,1]
            end
            tmp
        end
    else
      if !isempty(findall(var_list3d.==nam))
        fil=mdsio_list3d[ findall(var_list3d.==nam)[1] ]
        fil1=joinpath(pth_in,fil*list_steps[t][14:end])
	tmp=read_mdsio_alias(fil1,Symbol(nam))
        tmp=P.Γ.mskC*read(tmp,γ)     
      else
        fil=mdsio_list2d[ findall(var_list2d.==nam)[1] ]
        tmp=read_mdsio_alias(joinpath(pth_in,fil*list_steps[t][14:end]),Symbol(nam))
        tmp=P.Γ.mskC[:,1]*read(tmp,P.Γ.XC)
      end
    end
end

end #module ECCO_io

##

module ECCO_diagnostics

using SharedArrays, Distributed, Printf, JLD2, MeshArrays
import Climatology: ECCO_io, ECCO_helpers

"""
List of variables derived in this module:

- climatologies
- global means
- zonal means
- geographic maps
- transect transports
- MOC, MHT

Sample workflow:

```
## Setup Computation Parameters
@everywhere sol0="r2"
@everywhere nam="THETA"
@everywhere calc="clim"
@everywhere kk=1

## Preliminary Steps
@everywhere include("ECCO_pkg_grid_etc.jl")
@everywhere pth_in,pth_out,pth_tmp,sol,nt,list_steps=ECCO_path_etc(sol0,calc,nam)
!isdir(pth_out) ? mkdir(pth_out) : nothing
!isdir(pth_tmp) ? mkdir(pth_tmp) : nothing

## Main Computation
include("ECCO_standard_analysis.jl")
```
"""

## climatological mean

function comp_clim(P,tmp_m,tmp_s1,tmp_s2,m)
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol, γ, Γ) = P

    nm=length(m:12:nt)
    tmp_m[:,:,m].=0.0
    tmp_s1[:,:,m].=0.0
    tmp_s2[:,:,m].=0.0
    for t in m:12:nt
        tmp=ECCO_io.read_monthly(P,nam,t)
        ndims(tmp)>1 ? tmp=tmp[:,kk] : nothing
        tmp_m[:,:,m]=tmp_m[:,:,m]+1.0/nm*γ.write(tmp)
        tmp_s1[:,:,m]=tmp_s1[:,:,m]+γ.write(tmp)
        tmp_s2[:,:,m]=tmp_s2[:,:,m]+γ.write(tmp).^2
    end
end

function main_clim(P)
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol, γ, Γ) = P

    tmp_s1 = SharedArray{Float64}(γ.ioSize...,12)
    tmp_s2 = SharedArray{Float64}(γ.ioSize...,12)
    tmp_m = SharedArray{Float64}(γ.ioSize...,12)

    tmp=ECCO_io.read_monthly(P,nam,1)
    ndims(tmp)>1 ? nz=size(tmp,2) : nz=1
    nz==1 ? kk=1 : nothing
    nz>1 ? suff=Printf.@sprintf("_k%02d",kk) : suff=""

    @sync @distributed for m in 1:12
        comp_clim(P,tmp_m,tmp_s1,tmp_s2,m)
        GC.gc()
    end

    tmp0=read(tmp_m[:],γ)

    tmp=1.0/nt*sum(tmp_s1,dims=3)
    tmp1=read(tmp[:],Γ.XC)

    tmp=1/nt*sum(tmp_s2,dims=3)-tmp.^2
    tmp[findall(tmp.<0.0)].=0.0
    tmp=sqrt.(nt/(nt-1)*tmp)
    tmp2=read(tmp[:],Γ.XC)

    fil_out=joinpath(pth_out,nam*suff*".jld2")
    save(fil_out,"mean",tmp1,"std",tmp2,"mon",tmp0)

    return true
end

##

nansum(x) = sum(filter(!isnan,x))
nansum(x,y) = mapslices(nansum,x,dims=y)

## global mean

function comp_glo(P,glo,t)
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol, Γ) = P
    nr=length(Γ.DRF)

    tmp=ECCO_io.read_monthly(P,nam,t)
    if calc=="glo2d"
        tmp=[nansum(tmp[i,j].*Γ.RAC[i]) for j in 1:nr, i in eachindex(Γ.RAC)]
    else
        tmp=[nansum(tmp[i,j].*Γ.hFacC[i,j].*Γ.RAC[i]*Γ.DRF[j]) for j in 1:nr, i in eachindex(Γ.RAC)]
    end
    glo[:,t]=nansum(tmp,2)
end
    
function main_glo(P)
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol, Γ) = P
    nr=length(Γ.DRF)

    glo = SharedArray{Float64}(nr,nt)
    @sync @distributed for t in 1:nt
        comp_glo(P,glo,t)
        GC.gc()
    end

    if calc=="glo2d"
        tmp=[glo[r,t]/Γ.tot_RAC[r] for t in 1:nt, r in 1:nr]
    else
        tmp=[nansum(glo[:,t])/nansum(Γ.tot_VOL) for t in 1:nt]
    end
    save_object(joinpath(pth_out,calc*".jld2"),collect(tmp))
end

##

function comp_msk0(P,msk0,zm0,l)
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol, γ, Γ) = P
    nr=length(Γ.DRF)

    lats=load(joinpath(pth_out,calc*"_lats.jld2"),"single_stored_object")
    dlat=lats[2]-lats[1]
    la0=lats[l]-dlat/2
    la1=lats[l]+dlat/2
    if la1<0.0
        msk=1.0*(Γ.YC.>=la0)*(Γ.YC.<la1)
    elseif la0>0.0
        msk=1.0*(Γ.YC.>la0)*(Γ.YC.<=la1)
    else
        msk=1.0*(Γ.YC.>=la0)*(Γ.YC.<=la1)
    end
    msk[findall(msk.==0.0)].=NaN;
    msk0[:,:,l]=write(msk*Γ.RAC)

    tmp2=[nansum(Γ.mskC[i,j].*msk[i].*Γ.RAC[i]) for j in 1:nr, i in eachindex(Γ.RAC)]
    zm0[l,:]=1.0 ./nansum(tmp2,2)
end

function zmsum!(tmp1,tmp,msk,idx)
   tmp1.=0.0
   for j in 1:length(tmp1)
     for i in 1:length(idx)
       tmp1[j]+=tmp[idx[i],j]*msk[idx[i]]
     end
   end
end

function comp_zonmean(P,zm,t,msk0,zm0)
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol, γ, Γ) = P
    nl=size(msk0,3)
    idx0=[findall(msk0[:,:,l].>0) for l in 1:nl]
    comp_zonmean(P,zm,t,msk0,zm0,idx0)
end

function comp_zonmean(P,zm,t,msk0,zm0,idx0)
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol, γ, Γ) = P
    nr=length(Γ.DRF)

    lats=load(joinpath(pth_out,calc*"_lats.jld2"),"single_stored_object")
    nl=length(lats)

    tmp=write(ECCO_io.read_monthly(P,nam,t))
    tmp[findall(isnan.(tmp))].=0.0
    tmp1=zeros(nr)
    for l in 1:nl
        zmsum!(tmp1,tmp,msk0[:,:,l],idx0[l])
        zm[l,:,t]=tmp1.*zm0[l,:]
    end
end

function comp_zonmean2d(P,zm,t,msk0,zm0)
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol, γ, Γ) = P

    lats=load(joinpath(pth_out,calc*"_lats.jld2"),"single_stored_object")
    nl=length(lats)

    tmp=ECCO_io.read_monthly(P,nam,t)
    for l in 1:nl
        mskrac=read(msk0[:,:,l],γ)
        tmp1=[nansum(tmp[i].*mskrac[i]) for i in eachindex(Γ.RAC)]
        zm[l,t]=nansum(tmp1)*zm0[l,1]
    end
end

function main_zonmean(P)
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol, γ, Γ) = P
    nr=length(Γ.DRF)

    dlat=2.0
    lats=(-90+dlat/2:dlat:90-dlat/2)
    save_object(joinpath(pth_out,calc*"_lats.jld2"),collect(lats))
    nl=length(lats)

    msk0 = SharedArray{Float64}(γ.ioSize...,nl)
    zm0 = SharedArray{Float64}(nl,nr)
    @sync @distributed for l in 1:nl
        comp_msk0(P,msk0,zm0,l)
    end
    save_object(joinpath(pth_out,calc*"_zm0.jld2"),collect(zm0))
    save_object(joinpath(pth_out,calc*"_msk0.jld2"),collect(msk0))

    #to speed up main loop, reuse:
    #- precomputed msk*RAC once and for all
    #- precomputed 1.0./nansum(tmp2,2)

    msk0=load(joinpath(pth_out,calc*"_msk0.jld2"),"single_stored_object")
    zm0=load(joinpath(pth_out,calc*"_zm0.jld2"),"single_stored_object")
    idx0=[findall(msk0[:,:,l].>0) for l in 1:nl]

    if (calc=="zonmean")
        zm = SharedArray{Float64}(nl,nr,nt)
        @sync @distributed for t in 1:nt
            comp_zonmean(P,zm,t,msk0,zm0,idx0)
            GC.gc()
        end
    else
        zm = SharedArray{Float64}(nl,nt)
        @sync @distributed for t in 1:nt
            comp_zonmean2d(P,zm,t,msk0,zm0)
            GC.gc()
	end
    end
    save_object(joinpath(pth_out,calc*".jld2"),collect(zm))

    return true
end

##

function comp_overturn(P,ov,t)
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol, LC, Γ) = P

    nr=length(Γ.DRF)
    nl=length(LC)

    U=ECCO_io.read_monthly(P,"UVELMASS",t)
    V=ECCO_io.read_monthly(P,"VVELMASS",t)
    MeshArrays.UVtoTransport!(U,V,Γ)

    UV=Dict("U"=>0*U[:,1],"V"=>0*V[:,1],"dimensions"=>["x","y"])

    #integrate across latitude circles
    for z=1:nr
        UV["U"].=U[:,z]
        UV["V"].=V[:,z]
        [ov[l,z,t]=ThroughFlow(UV,LC[l],Γ) for l=1:nl]
    end
    #integrate from bottom
    ov[:,:,t]=reverse(cumsum(reverse(ov[:,:,t],dims=2),dims=2),dims=2)
    #
    true
end

function main_overturn(P)  
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol, LC, Γ) = P

    nr=length(Γ.DRF)
    nl=length(LC)

    ov = SharedArray{Float64}(nl,nr,nt)
    @sync @distributed for t in 1:nt
        comp_overturn(P,ov,t)
        GC.gc()
    end
    
    save_object(joinpath(pth_out,calc*".jld2"),collect(ov))
	"Done with overturning"
end

##

function comp_MHT(P,MHT,t)
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol, LC, Γ) = P

    nr=length(Γ.DRF)
    nl=length(LC)

    U=ECCO_io.read_monthly(P,"ADVx_TH",t)+ECCO_io.read_monthly(P,"DFxE_TH",t)
    V=ECCO_io.read_monthly(P,"ADVy_TH",t)+ECCO_io.read_monthly(P,"DFyE_TH",t)

    [U[i][findall(isnan.(U[i]))].=0.0 for i in eachindex(U)]
    [V[i][findall(isnan.(V[i]))].=0.0 for i in eachindex(V)]
    Tx=0.0*U[:,1]
    Ty=0.0*V[:,1]
    [Tx=Tx+U[:,z] for z=1:nr]
    [Ty=Ty+V[:,z] for z=1:nr]

    UV=Dict("U"=>Tx,"V"=>Ty,"dimensions"=>["x","y"])
    [MHT[l,t]=1e-15*4e6*ThroughFlow(UV,LC[l],Γ) for l=1:nl]
end

function main_MHT(P)  
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol, LC) = P

    nl=length(LC)
    MHT = SharedArray{Float64}(nl,nt)
    @sync @distributed for t in 1:nt
        comp_MHT(P,MHT,t)
        GC.gc()
    end
    save_object(joinpath(pth_out,calc*".jld2"),collect(MHT))
    "Done with MHT"
end

##

function comp_trsp(P,trsp,t)
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol, Γ) = P

    U=ECCO_io.read_monthly(P,"UVELMASS",t)
    V=ECCO_io.read_monthly(P,"VVELMASS",t)
    MeshArrays.UVtoTransport!(U,V,Γ)

    UV=Dict("U"=>0*U[:,1],"V"=>0*V[:,1],"dimensions"=>["x","y"])

    pth_trsp=joinpath(pth_out,"..","ECCO_transport_lines")
    list_trsp,msk_trsp,ntr=ECCO_helpers.reload_transport_lines(pth_trsp)

    #integrate across transport lines
    for z=1:length(Γ.DRF)
        UV["U"].=U[:,z]
        UV["V"].=V[:,z]
        [trsp[itr,z,t]=ThroughFlow(UV,msk_trsp[itr],Γ) for itr=1:ntr]
    end
end

function main_trsp(P) 
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol) = P

    list_trsp=readdir(joinpath(pth_out,"..","ECCO_transport_lines"))
    ntr=length(list_trsp)
 
    nr=length(P.Γ.DRF)
    trsp = SharedArray{Float64}(ntr,nr,nt)
    @sync @distributed for t in 1:nt
        comp_trsp(P,trsp,t)
        GC.gc()
    end
    
    trsp=[(nam=list_trsp[itr],val=trsp[itr,:,:]) for itr=1:ntr]
    save_object(joinpath(pth_out,calc*".jld2"),collect(trsp))
	"Done with transports"
end

"""
    driver(P)

Call main computation loop as specified by parameters `P`.

The main computation loop choice depends on the `P` parameter values. Methods include:

- `main_clim`
- `main_glo`
- `main_zonmean`
- `main_overturn`
- `main_MHT`
- `main_trsp`
"""
function driver(P)
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol) = P

    if calc=="clim"
        main_clim(P)
    elseif (calc=="glo2d")||(calc=="glo3d")
        main_glo(P)
    elseif (calc=="zonmean")||(calc=="zonmean2d")
        main_zonmean(P)
    elseif (calc=="overturn")
        main_overturn(P)
    elseif (calc=="MHT")
        main_MHT(P)
    elseif (calc=="trsp")
        main_trsp(P)
    else
        println("unknown calc")
    end
end

end #module ECCO_diagnostics

##

module ECCO_procs

using JLD2, MeshArrays, DataDeps, Statistics, Climatology, TOML
import Climatology: ECCOdiag

function longname(n)
	if occursin("_k",n)
		ln=split(n,"_k")[1]*" at level "*split(n,"_k")[2]
	else
		ln=n
	end
	occursin("BSF",ln) ? ln=replace(ln, "BSF" => "Horizontal Streamfunction (m3/s)") : nothing
	occursin("MXLDEPTH",ln) ? ln=replace(ln, "MXLDEPTH" => "Mixed Layer Depth (m)") : nothing
	occursin("SIarea",ln) ? ln=replace(ln, "SIarea" => "Ice Concentration (0 to 1)") : nothing
	occursin("SSH",ln) ? ln=replace(ln, "SSH" => "Free Surface Height (m)") : nothing
	occursin("THETA",ln) ? ln=replace(ln, "THETA" => "Potential Temperature (degree C)") : nothing
	occursin("SALT",ln) ? ln=replace(ln, "SALT" => "Salinity (psu)") : nothing
	return ln
end

function climatology_files(pth_out)
	list_clim=readdir(pth_out)
	kk=findall(occursin.(Ref("clim"),list_clim))
	list_clim=list_clim[kk]
	clim_files=[]
	for ii in 1:length(list_clim)
		tmp=joinpath.(Ref(list_clim[ii]),readdir(joinpath(pth_out,list_clim[ii])))
		[push!(clim_files,i) for i in tmp]
	end
	clim_files
end

##

function years_min_max(sol)
	year0=1992
	year1=2011
	if occursin("ECCOv4r3",sol)
		year1=2015
	elseif occursin("ECCOv4r4",sol)
		year1=2017
	elseif occursin("ECCOv4r5",sol)
		year1=2019
	elseif occursin("OCCA2HR1",sol)
		year0=1980
		year1=2024
	elseif occursin("OCCA2HR2",sol)
		year0=1960
		year1=2024
	end
	return year0,year1
end

##

function parameters()

	γ=GridSpec(ID=:LLC90)
	Γ=GridLoad(γ;option="full")
	#LC=LatitudeCircles(-89.0:89.0,Γ)

	μ = land_mask(Γ)
    λ = interpolation_setup()
	path0=ECCOdiags_add("OCCA2HR1")

    tmp=load(ECCOdiag(path=path0,name="trsp"))
	ntr=length(tmp)
	list_trsp=[vec(tmp)[i].nam for i in 1:ntr] 
	list_trsp=[i[1:end-5] for i in list_trsp]

	pth_colors=joinpath(dirname(pathof(Climatology)),"..","examples","ECCO")	
	clim_colors1=TOML.parsefile(joinpath(pth_colors,"clim_colors1.toml"))
	clim_colors2=TOML.parsefile(joinpath(pth_colors,"clim_colors2.toml"))


	clim_files=climatology_files(path0)
	clim_name=[split(basename(f),'.')[1] for f in clim_files]
	clim_longname=longname.(clim_name) 

	#"Done with listing solutions, file names, color codes"
	(γ=γ,Γ=Γ,λ=λ,μ=μ,list_trsp=list_trsp,
	clim_colors1=clim_colors1,clim_colors2=clim_colors2,
	clim_files=clim_files,clim_name=clim_name,clim_longname=clim_longname)
end

##

function glo(pth_out,nam,k,year0,year1)
    nam_full=nam*(k>0 ? "_glo2d" : "_glo3d")
    tmp=load(ECCOdiag(path=pth_out,name=nam_full))

	occursin("THETA",nam) ? ln=longname("THETA") : ln=longname("SALT")
	if k>0
		nt=Int(length(tmp[:])./50.0)
		tmp=reshape(tmp,(nt,50))
		tmp=tmp[:,k]
		occursin("THETA",fil) ? rng=[18.0,19.0] : rng=[34.65,34.80]
		txt=ln*" -- level $(k)" 
		k>1 ? rng=[extrema(tmp)...] : nothing
	else
		nt=length(tmp[:])
		occursin("THETA",nam) ? rng=[3.5,3.65] : rng=[34.724,34.728]
		txt=ln
	end

	x=vec(0.5:nt)
	x=year0 .+ x./12.0

	(y=tmp,txt=txt,rng=rng,x=x)
end

function map(nammap,P,statmap,timemap,pth_out)
	ii=findall(P.clim_longname.==nammap)[1]
	nam=P.clim_name[ii]; file=nam*".jld2"
    nam_full=split(nam,"_")[1]*"_clim"
    tmp=load(ECCOdiag(path=pth_out,name=nam_full),file=file,variable=statmap)
	tmp=(statmap!=="mon" ? tmp : tmp[:,timemap])

	DD=Interpolate(P.μ*tmp,P.λ.f,P.λ.i,P.λ.j,P.λ.w)
	DD=reshape(DD,size(P.λ.lon))
	#DD[findall(DD.==0.0)].=NaN
	statmap=="std" ? rng=P.clim_colors2[nam] : rng=P.clim_colors1[nam]
	levs=rng[1] .+collect(0.0:0.05:1.0)*(rng[2]-rng[1])

	ttl=P.clim_longname[ii]
	(λ=P.λ,field=DD,levels=levs,title=ttl)
end

function TimeLat(namzm,pth_out,year0,year1,cmap_fac,k_zm,P)
	fn(x)=transpose(x);
	if namzm=="MXLDEPTH"
		levs=(0.0:50.0:400.0); cm=:turbo
		dlat=2.0; y=vec(-90+dlat/2:dlat:90-dlat/2)
		nam=namzm*"_zonmean2d"
	elseif namzm=="SIarea"
		levs=(0.0:0.1:1.0); cm=:turbo
		dlat=2.0; y=vec(-90+dlat/2:dlat:90-dlat/2)
		nam=namzm*"_zonmean2d"
	elseif namzm=="THETA"
		levs=(-2.0:2.0:34.0); cm=:turbo
		dlat=2.0; y=vec(-90+dlat/2:dlat:90-dlat/2)
		nam=namzm*"_zonmean"
	elseif namzm=="SALT"
		levs=(32.6:0.2:36.2); cm=:turbo
		dlat=2.0; y=vec(-90+dlat/2:dlat:90-dlat/2)
		nam=namzm*"_zonmean"
	elseif (namzm=="ETAN")||(namzm=="SSH")
		levs=10*(-0.15:0.02:0.15); cm=:turbo
		dlat=2.0; y=vec(-90+dlat/2:dlat:90-dlat/2)
		nam=namzm*"_zonmean2d"
	else
		levs=missing
        nam="missing"
	end

	tmp=load(ECCOdiag(path=pth_out,name=nam))
	if length(size(tmp))==3
		z=fn(tmp[:,k_zm,:])
		x=vec(0.5:size(tmp,3))
		addon1=" at $(Int(round(P.Γ.RC[k_zm])))m "
	else
		z=fn(tmp[:,:])
		x=vec(0.5:size(tmp,2))
		addon1=""
	end

	x=year0 .+ x./12.0
	ttl="$(longname(namzm)) : Zonal Mean $(addon1)"
	(x=x,y=y,z=z,levels=cmap_fac*levs,title=ttl,ylims=(-90.0,90.0),year0=year0,year1=year1)
end

function TimeLatAnom(namzmanom2d,pth_out,year0,year1,cmap_fac,k_zm2d,l0,l1,P)
	namzm=namzmanom2d
	if namzm=="MXLDEPTH"
		levs=(-100.0:25.0:100.0)/2.0; fn=transpose; cm=:turbo
		nam=namzm*"_zonmean2d"
	elseif namzm=="SIarea"
		levs=(-0.5:0.1:0.5)/5.0; fn=transpose; cm=:turbo
		nam=namzm*"_zonmean2d"
	elseif namzm=="THETA"
		levs=(-2.0:0.25:2.0)/5.0; fn=transpose; cm=:turbo
		nam=namzm*"_zonmean"
	elseif namzm=="SALT"
		levs=(-0.5:0.1:0.5)/5.0; fn=transpose; cm=:turbo
		nam=namzm*"_zonmean"
	elseif (namzm=="ETAN")||(namzm=="SSH")
		levs=(-0.5:0.1:0.5)/2.0; fn=transpose; cm=:turbo
        nam=namzm*"_zonmean2d"
    else
		fn=transpose
		levs=missing
        nam="missing"
    end

	tmp=load(ECCOdiag(path=pth_out,name=nam))
	if length(size(tmp))==3
		z=fn(tmp[:,k_zm2d,:])
		x=vec(0.5:size(tmp,3)); 
		addon1=" -- at $(Int(round(P.Γ.RC[k_zm2d])))m "
	else
		z=fn(tmp[:,:])
		x=vec(0.5:size(tmp,2)); 
		addon1=""
	end

	dlat=2.0; y=vec(-90+dlat/2:dlat:90-dlat/2)
	nt=size(z,1)

	m0=(1992-year0)*12
	
	if true
		#a. subtract monthly mean
		ref1="1992-2011 monthy mean"
		for m in 1:12
			zmean=vec(mean(z[m0+m:12:m0+240,:],dims=1))
			[z[t,:]=z[t,:]-zmean for t in m:12:nt]
		end
	else
		#b. subtract time mean
		ref1="1992-2011 annual mean"
		zmean=vec(mean(z[m0+1:m0+240,:],dims=1))
		[z[t,:]=z[t,:]-zmean for t in 1:nt]
	end

	x=1992.0-m0/12.0 .+ x./12.0
	ttl="$(longname(namzm)) -- minus $(ref1) $(addon1)"

	(x=x,y=y,z=z,levels=cmap_fac*levs,title=ttl,ylims=(y[l0],y[l1]),year0=year0,year1=year1)
end

fn_DepthTime(x)=transpose(x)	

function DepthTime(namzmanom,pth_out,facA,l_Tzm,year0,year1,k0,k1,P)
if namzmanom=="THETA"
	levs=(-3.0:0.4:3.0)/8.0; cm=:turbo
elseif namzmanom=="SALT"
	levs=(-0.5:0.1:0.5)/10.0;cm=:turbo
else
	levs=missing;
end
nam_full=namzmanom*"_zonmean"
tmp=load(ECCOdiag(path=pth_out,name=nam_full))

dlat=2.0
lats=(-90+dlat/2:dlat:90-dlat/2)

z=fn_DepthTime(tmp[l_Tzm,:,:])
addon1=" -- at $(lats[l_Tzm])N "
x=vec(0.5:size(tmp,3)); 
y=vec(P.Γ.RC)
nt=size(tmp,3)

#a. subtract monthly mean
ref1="1992-2011 monthy mean"
m0=(1992-year0)*12
for m in 1:12
	zmean=vec(mean(z[m0+m:12:m0+240,:],dims=1))
	[z[t,:]=z[t,:]-zmean for t in m:12:nt]
end
#b. subtract time mean
#ref1="1992-2011 annual mean"
#zmean=vec(mean(z[1:240,:],dims=1))
#[z[t,:]=z[t,:]-zmean for t in 1:nt]

x=year0 .+ x./12.0
ttl="$(longname(namzmanom)) -- minus $(ref1) $(addon1)"

(x=x,y=y,z=z,levels=facA*levs,title=ttl,ylims=(P.Γ.RC[k1],P.Γ.RC[k0]),year0=year0,year1=year1)

end

end #module ECCO_procs

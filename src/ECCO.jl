
module ECCO

using Pkg
import OceanStateEstimation: pkg_pth

"""
    ECCO.standard_analysis_setup(pth0)

Create a run folder where folder `pth0` will be linked. 

Folder `pth0` (of type `String`) should be the path to the user's ECCO data folder.

```
using Revise
using OceanStateEstimation, Pkg

#pth0="nctiles_monthly" #edit path as needed
pth=ECCO.standard_analysis_setup(pth0)

Pkg.activate(pth)
Pkg.instantiate()

include("ECCO_standard_loop.jl")
```
"""
function standard_analysis_setup(pth0)
	println(pth0)
	
	#1. setup run folder and create link to ECCO data folder
	pth=joinpath(tempdir(),"ECCO_diags_dev"); 
	!isdir(pth) ? mkdir(pth) : nothing
	pth1=joinpath(pth,"ECCOv4r2")
	!isdir(pth1) ? mkdir(pth1) : nothing
	link0=joinpath(pth1,"nctiles_monthly")
	!isfile(link0)&& !islink(link0) ? symlink(pth0,link0) : nothing
	
	#2. copy Project.toml to run folder
	tmp0=pkg_pth
	tmp1=joinpath(tmp0,"..","examples","ECCO","ECCO_standard_Project.toml")
	tmp2=joinpath(pth,"Project.toml")
	!isfile(tmp2) ? cp(tmp1,tmp2) : nothing
		
	return pth
end


end

##

module ECCO_helper_functions

using MeshArrays, TOML, JLD2

function parameters(pth0::String,sol0::String,list0,gg)

    calc=list0["calc"][gg]
    nam=list0["nam"][gg]
    kk=list0["kk"][gg]
    sol="ECCOv4"*sol0*"_analysis"

    if sol0=="r1"||sol0=="r2"
        fil=isfile(joinpath(pth0,"THETA","THETA.0001.nc"))
        if isfile(fil)
            nt=NCDataset(fil) do ds
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
    end

    if sol0!=="r5"
        pth_in=joinpath(pth0,"ECCOv4"*sol0,"nctiles_monthly")
    else
        pth_in=joinpath(pth0,"ECCOv4"*sol0,"diags")
    end
    list_steps=list_time_steps(pth_in)

    pth_out=joinpath(pth0,sol)

    if sum(calc.==("overturn","MHT","trsp"))==0
        pth_out=joinpath(pth_out,nam*"_"*calc)
    else
        pth_out=joinpath(pth_out,calc)
    end    

    P=(pth_in=pth_in,pth_out=pth_out,list_steps=list_steps,nt=nt,
    calc=calc,nam=nam,kk=kk,sol=sol)
end

#STATE/state_3d_set1.0000241020.meta
#    'THETA   ' 'SALT    ' 'DRHODR  '
#TRSP/trsp_3d_set1.0000241020.meta
#    'UVELMASS' 'VVELMASS' 'WVELMASS' 'GM_PsiX ' 'GM_PsiY '
#TRSP/trsp_3d_set3.0000241020.meta
#    'DFxE_TH ' 'DFyE_TH ' 'ADVx_TH ' 'ADVy_TH ' 'DFxE_SLT' 'DFyE_SLT' 'ADVx_SLT' 'ADVy_SLT'

function list_time_steps(pth_in)
    if isdir(joinpath(pth_in,"STATE"))
        list=readdir(joinpath(pth_in,"STATE"))
        list=list[findall(occursin.(Ref("state_3d_set1"),list))]
        list=list[findall(occursin.(Ref("data"),list))]
        [list[i][15:end] for i in 1:length(list)]
    else
        list=[]
    end
    return list
end

function get_nr_nt(pth_in,nam)
    nct_path=joinpath(pth_in,nam)
    lst=readdir(nct_path)
    lst=lst[findall(occursin.(Ref(".nc"),lst))]
    fil1=joinpath(nct_path,lst[1])
    ds=NCTiles.NCDatasets.Dataset(fil1)
    siz=size(ds[nam])
    siz[end-1],siz[end]
end

nansum(x) = sum(filter(!isnan,x))
nansum(x,y) = mapslices(nansum,x,dims=y)

function GridLoad_Main()

    γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
    nr=50
    XC=GridLoadVar("XC",γ)
    YC=GridLoadVar("YC",γ)
    RAC=GridLoadVar("RAC",γ)
    hFacC=GridLoadVar("hFacC",γ)
    hFacW=GridLoadVar("hFacW",γ)
    hFacS=GridLoadVar("hFacS",γ)
    DRF=GridLoadVar("DRF",γ)
    DXC=GridLoadVar("DXC",γ)
    DYC=GridLoadVar("DYC",γ)
    DXG=GridLoadVar("DXG",γ)
    DYG=GridLoadVar("DYG",γ)

    mskC=hFacC./hFacC

    tmp=hFacC./hFacC
    tmp=[nansum(tmp[i,j].*RAC[i]) for j in 1:nr, i in eachindex(RAC)]
    tot_RAC=nansum(tmp,2)

    tmp=[nansum(hFacC[i,j].*RAC[i].*DRF[j]) for j in 1:nr, i in eachindex(RAC)]
    tot_VOL=nansum(tmp,2)

    Γ=(XC=XC,YC=YC,RAC=RAC,DXG=DXG,DYG=DYG,DXC=DXC,DYC=DYC,DRF=DRF,
        hFacC=hFacC,hFacW=hFacW,hFacS=hFacS,
        mskC=mskC,tot_RAC=tot_RAC,tot_VOL=tot_VOL)

	LC=LatitudeCircles(-89.0:89.0,Γ)

    return γ,Γ,LC
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
    push!(allcalc,allnam,allkk;calc="clim",nam="BSF")
    push!(allcalc,allnam,allkk;calc="clim",nam="SSH")
    push!(allcalc,allnam,allkk;calc="clim",nam="MXLDEPTH")
    push!(allcalc,allnam,allkk;calc="clim",nam="SIarea")

    tmp1=Dict("calc"=>allcalc,"nam"=>allnam,"kk"=>allkk)
    open(fil, "w") do io
        TOML.print(io, tmp1)
    end

end

##

#temporary fix:
γ,Γ,LC=GridLoad_Main()

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

function transport_lines(pth_trsp)
    mkdir(pth_trsp)
    lonPairs,latPairs,namPairs=transport_lines()
    for ii in 1:length(lonPairs)
        lons=Float64.(lonPairs[ii])
        lats=Float64.(latPairs[ii])
        name=namPairs[ii]
        Trsct=Transect(name,lons,lats,Γ)
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

end #module ECCO_helper_functions

## generic read function

module ECCO_read_monthly

using MeshArrays, MITgcmTools

#temporary fix:
import OceanStateEstimation: ECCO_helper_functions
γ,Γ,LC=ECCO_helper_functions.GridLoad_Main()

function main(sol,nam,t,list_steps; pth_in="", pth_out="") 
    if nam=="SSH"
        read_monthly_SSH(sol,nam,t,list_steps; pth_in=pth_in, pth_out=pth_out)
    elseif nam=="MHT"
        read_monthly_MHT(sol,nam,t,list_steps; pth_in=pth_in, pth_out=pth_out)
    elseif nam=="BSF"
        read_monthly_BSF(sol,nam,t,list_steps; pth_in=pth_in, pth_out=pth_out)
    else
        read_monthly_default(sol,nam,t,list_steps; pth_in=pth_in, pth_out=pth_out)
    end
end

function read_monthly_SSH(sol,nam,t,list_steps; pth_in="", pth_out="")
    ETAN=read_monthly_default(sol,"ETAN",t,list_steps; pth_in=pth_in, pth_out=pth_out)
    sIceLoad=read_monthly_default(sol,"sIceLoad",t,list_steps; pth_in=pth_in, pth_out=pth_out)
    (ETAN+sIceLoad/1029.0)*Γ.mskC[:,1]
end

function read_monthly_MHT(sol,nam,t,list_steps; pth_in="", pth_out="")
    U=read_monthly_default(sol,"ADVx_TH",t,list_steps; pth_in=pth_in, pth_out=pth_out)
    V=read_monthly_default(sol,"ADVy_TH",t,list_steps; pth_in=pth_in, pth_out=pth_out)
    U=U+read_monthly_default(sol,"DFxE_TH",t,list_steps; pth_in=pth_in, pth_out=pth_out)
    V=V+read_monthly_default(sol,"DFyE_TH",t,list_steps; pth_in=pth_in, pth_out=pth_out)

    [U[i][findall(isnan.(U[i]))].=0.0 for i in eachindex(U)]
    [V[i][findall(isnan.(V[i]))].=0.0 for i in eachindex(V)]

    Tx=0.0*U[:,1]
    Ty=0.0*V[:,1]
    [Tx=Tx+U[:,z] for z=1:nr]
    [Ty=Ty+V[:,z] for z=1:nr]

    return Tx,Ty
end

function read_monthly_BSF(sol,nam,t,list_steps; pth_in="", pth_out="")
    U=read_monthly_default(sol,"UVELMASS",t,list_steps; pth_in=pth_in, pth_out=pth_out)
    V=read_monthly_default(sol,"VVELMASS",t,list_steps; pth_in=pth_in, pth_out=pth_out)
    (Utr,Vtr)=UVtoTransport(U,V,Γ)
    
    nz=size(Γ.hFacC,2)
    μ=Γ.mskC[:,1]
    Tx=0.0*Utr[:,1]
	Ty=0.0*Vtr[:,1]
	for z=1:nz
		Tx=Tx+Utr[:,z]
		Ty=Ty+Vtr[:,z]
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

    return TrspPsi
end

function read_monthly_default(sol,nam,t,list; pth_in="", pth_out="")
    var_list3d=("THETA","SALT","UVELMASS","VVELMASS",
                "ADVx_TH","ADVy_TH","DFxE_TH","DFyE_TH")
    mdsio_list3d=("STATE/state_3d_set1","STATE/state_3d_set1",
        "TRSP/trsp_3d_set1","TRSP/trsp_3d_set1","TRSP/trsp_3d_set2",
        "TRSP/trsp_3d_set2","TRSP/trsp_3d_set2","TRSP/trsp_3d_set2")

    var_list2d=("MXLDEPTH","SIarea","sIceLoad","ETAN")
    mdsio_list2d=("STATE/state_2d_set1","STATE/state_2d_set1",
                  "STATE/state_2d_set1","STATE/state_2d_set1")

    if (sol=="ECCOv4r1_analysis")||(sol=="ECCOv4r2_analysis")||(sol=="ECCOv4r3_analysis")
        nct_path=joinpath(pth_in,nam)
        if sum(var_list3d.==nam)==1
            tmp=read_nctiles(nct_path,nam,γ,I=(:,:,:,t))
        else
            tmp=read_nctiles(nct_path,nam,γ,I=(:,:,t))
        end
    elseif (sol=="ECCOv4r4_analysis")
        y0=Int(floor((t-1)/12))+1992
        m0=mod1(t,12)
        nct_path=joinpath(pth_in,nam,string(y0))
        m0<10 ? fil=nam*"_$(y0)_0$(m0).nc" : fil=nam*"_$(y0)_$(m0).nc"
        tmp0=NCTiles.NCDatasets.Dataset(joinpath(nct_path,fil))[nam]
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
        ii=findall(var_list3d.==nam)[1];
        fil=mdsio_list3d[ii]
        meta=read_meta(joinpath(pth_in,fil*".0000241020.meta"))
        kk=findall(vec(meta.fldList).==nam)[1]
        tmp=read_mdsio(joinpath(pth_in,fil*list[t][14:end]))[:,:,:,kk]
        tmp=Γ.mskC*read(tmp,γ)     
      else
        ii=findall(var_list2d.==nam)[1];
        fil=mdsio_list2d[ii]
        meta=read_meta(joinpath(pth_in,fil*".0000241020.meta"))
        kk=findall(vec(meta.fldList).==nam)[1]
        tmp=read_mdsio(joinpath(pth_in,fil*list[t][14:end]))[:,:,kk]
        tmp=Γ.mskC[:,1]*read(tmp,Γ.XC)
      end
    end
end

end #module ECCO_read_monthly

##

module ECCO_diagnostics

using SharedArrays, Distributed, Printf, JLD2, MeshArrays
import OceanStateEstimation: ECCO_read_monthly

#temporary fix:
read_monthly=ECCO_read_monthly.main

#temporary fix:
import OceanStateEstimation: ECCO_helper_functions
γ,Γ,LC=ECCO_helper_functions.GridLoad_Main()
nr=length(Γ.DRF)
nl=length(LC)
#nt=12

list_time_steps=ECCO_helper_functions.list_time_steps

#list_trsp,msk_trsp,ntr=ECCO_helper_functions.reload_transport_lines(pth0)

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

function comp_clim(sol,nam,tmp_m,tmp_s1,tmp_s2,kk,m; pth_in="", pth_out="")
    list_steps=list_time_steps(pth_in)
    nm=length(m:12:nt)
    tmp_m[:,:,m].=0.0
    tmp_s1[:,:,m].=0.0
    tmp_s2[:,:,m].=0.0
    for t in m:12:nt
        tmp=read_monthly(sol,nam,t,list_steps; pth_in=pth_in, pth_out=pth_out)
        ndims(tmp)>1 ? tmp=tmp[:,kk] : nothing
        tmp_m[:,:,m]=tmp_m[:,:,m]+1.0/nm*γ.write(tmp)
        tmp_s1[:,:,m]=tmp_s1[:,:,m]+γ.write(tmp)
        tmp_s2[:,:,m]=tmp_s2[:,:,m]+γ.write(tmp).^2
    end
end

function main_clim(P)
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol) = P

    list_steps=list_time_steps(pth_in)

    tmp_s1 = SharedArray{Float64}(γ.ioSize...,12)
    tmp_s2 = SharedArray{Float64}(γ.ioSize...,12)
    tmp_m = SharedArray{Float64}(γ.ioSize...,12)

    tmp=read_monthly(sol,nam,1,list_steps; pth_in=pth_in, pth_out=pth_out)
    ndims(tmp)>1 ? nz=size(tmp,2) : nz=1
    nz==1 ? kk=1 : nothing
    nz>1 ? suff=Printf.@sprintf("_k%02d",kk) : suff=""

    @sync @distributed for m in 1:12
        comp_clim(sol,nam,tmp_m,tmp_s1,tmp_s2,kk,m; pth_in=pth_in, pth_out=pth_out)
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

function comp_glo(sol,nam,calc,glo,t; pth_in="", pth_out="")
    list_steps=list_time_steps(pth_in)
    tmp=read_monthly(sol,nam,t,list_steps; pth_in=pth_in, pth_out=pth_out)
    if calc=="glo2d"
        tmp=[nansum(tmp[i,j].*Γ.RAC[i]) for j in 1:nr, i in eachindex(Γ.RAC)]
    else
        tmp=[nansum(tmp[i,j].*Γ.hFacC[i,j].*Γ.RAC[i]*Γ.DRF[j]) for j in 1:nr, i in eachindex(Γ.RAC)]
    end
    glo[:,t]=nansum(tmp,2)
end
    
function main_glo(P)
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol) = P

    glo = SharedArray{Float64}(nr,nt)
    @sync @distributed for t in 1:nt
        comp_glo(sol,nam,calc,glo,t; pth_in=pth_in, pth_out=pth_out)
    end

    if calc=="glo2d"
        tmp=[glo[r,t]/Γ.tot_RAC[r] for t in 1:nt, r in 1:nr]
    else
        tmp=[nansum(glo[:,t])/nansum(Γ.tot_VOL) for t in 1:nt]
    end
    save_object(joinpath(pth_out,calc*".jld2"),collect(tmp))
end

##

function comp_msk0(msk0,zm0,l,calc; pth_in="", pth_out="")
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

function comp_zonmean(sol,nam,calc,zm,t,msk0,zm0; pth_in="", pth_out="")
    lats=load(joinpath(pth_out,calc*"_lats.jld2"),"single_stored_object")
    nl=length(lats)
    list_steps=list_time_steps(pth_in)
    tmp=read_monthly(sol,nam,t,list_steps; pth_in=pth_in, pth_out=pth_out)
    for l in 1:nl
        mskrac=read(msk0[:,:,l],γ)
        tmp1=[nansum(tmp[i,j].*mskrac[i]) for j in 1:nr, i in eachindex(Γ.RAC)]
        zm[l,:,t]=nansum(tmp1,2).*zm0[l,:]
    end
end

function comp_zonmean2d(sol,nam,calc,zm,t,msk0,zm0; pth_in="", pth_out="")
    lats=load(joinpath(pth_out,calc*"_lats.jld2"),"single_stored_object")
    nl=length(lats)
    list_steps=list_time_steps(pth_in)
    tmp=read_monthly(sol,nam,t,list_steps; pth_in=pth_in, pth_out=pth_out)
    for l in 1:nl
        mskrac=read(msk0[:,:,l],γ)
        tmp1=[nansum(tmp[i].*mskrac[i]) for i in eachindex(Γ.RAC)]
        zm[l,t]=nansum(tmp1)*zm0[l,1]
    end
end

function main_zonmean(P)
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol) = P

    dlat=2.0
    lats=(-90+dlat/2:dlat:90-dlat/2)
    save_object(joinpath(pth_out,calc*"_lats.jld2"),collect(lats))
    nl=length(lats)

    msk0 = SharedArray{Float64}(γ.ioSize...,nl)
    zm0 = SharedArray{Float64}(nl,nr)
    @sync @distributed for l in 1:nl
        comp_msk0(msk0,zm0,l,calc; pth_in=pth_in, pth_out=pth_out)
    end
    save_object(joinpath(pth_out,calc*"_zm0.jld2"),collect(zm0))
    save_object(joinpath(pth_out,calc*"_msk0.jld2"),collect(msk0))

    #to speed up main loop, reuse:
    #- precomputed msk*RAC once and for all
    #- precomputed 1.0./nansum(tmp2,2)

    msk0=load(joinpath(pth_out,calc*"_msk0.jld2"),"single_stored_object")
    zm0=load(joinpath(pth_out,calc*"_zm0.jld2"),"single_stored_object")

    if (calc=="zonmean")
        zm = SharedArray{Float64}(nl,nr,nt)
        @sync @distributed for t in 1:nt
            comp_zonmean(sol,nam,calc,zm,t,msk0,zm0; pth_in=pth_in, pth_out=pth_out)
        end
    else
        zm = SharedArray{Float64}(nl,nt)
        @sync @distributed for t in 1:nt
            comp_zonmean2d(sol,nam,calc,zm,t,msk0,zm0; pth_in=pth_in, pth_out=pth_out)
        end
    end
    save_object(joinpath(pth_out,calc*".jld2"),collect(zm))

    return true
end

##

function comp_overturn(sol,ov,t; pth_in="", pth_out="")
    list_steps=list_time_steps(pth_in)
    U=read_monthly(sol,"UVELMASS",t,list_steps; pth_in=pth_in, pth_out=pth_out)
    V=read_monthly(sol,"VVELMASS",t,list_steps; pth_in=pth_in, pth_out=pth_out)
    (Utr,Vtr)=UVtoTransport(U,V,Γ)
    #integrate across latitude circles
    for z=1:nr
        UV=Dict("U"=>Utr[:,z],"V"=>Vtr[:,z],"dimensions"=>["x","y"])
        [ov[l,z,t]=ThroughFlow(UV,LC[l],Γ) for l=1:nl]
    end
    #integrate from bottom
    ov[:,:,t]=reverse(cumsum(reverse(ov[:,:,t],dims=2),dims=2),dims=2)
    #
    true
end

function main_overturn(P)  
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol) = P

    ov = SharedArray{Float64}(nl,nr,nt)
    @sync @distributed for t in 1:nt
        comp_overturn(sol,ov,t; pth_in=pth_in, pth_out=pth_out)
    end
    
    save_object(joinpath(pth_out,calc*".jld2"),collect(ov))
	"Done with overturning"
end

##

function comp_MHT(sol,nam,MHT,t; pth_in="", pth_out="")
    list_steps=list_time_steps(pth_in)
    U=read_monthly(sol,"ADVx_TH",t,list_steps; pth_in=pth_in, pth_out=pth_out)
    V=read_monthly(sol,"ADVy_TH",t,list_steps; pth_in=pth_in, pth_out=pth_out)
    U=U+read_monthly(sol,"DFxE_TH",t,list_steps; pth_in=pth_in, pth_out=pth_out)
    V=V+read_monthly(sol,"DFyE_TH",t,list_steps; pth_in=pth_in, pth_out=pth_out)

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
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol) = P

    MHT = SharedArray{Float64}(nl,nt)
    @sync @distributed for t in 1:nt
        comp_MHT(sol,nam,MHT,t; pth_in=pth_in, pth_out=pth_out)
    end
    save_object(joinpath(pth_out,calc*".jld2"),collect(MHT))
	"Done with MHT"
end

##

function comp_trsp(sol,nam,trsp,t; pth_in="", pth_out="")
    list_steps=list_time_steps(pth_in)
    U=read_monthly(sol,"UVELMASS",t,list_steps; pth_in=pth_in, pth_out=pth_out)
    V=read_monthly(sol,"VVELMASS",t,list_steps; pth_in=pth_in, pth_out=pth_out)
    (Utr,Vtr)=UVtoTransport(U,V,Γ)

    pth_trsp=joinpath(pth_out,"..","ECCO_transport_lines")
    #list_trsp=readdir(pth_trsp)
    #ntr=length(list_trsp)
    list_trsp,msk_trsp,ntr=ECCO_helper_functions.reload_transport_lines(pth_trsp)

    #integrate across transport lines
    for z=1:nr
        UV=Dict("U"=>Utr[:,z],"V"=>Vtr[:,z],"dimensions"=>["x","y"])
        [trsp[itr,z,t]=ThroughFlow(UV,msk_trsp[itr],Γ) for itr=1:ntr]
    end
end

function main_trsp(P) 
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol) = P

    list_trsp=readdir(joinpath(pth_out,"..","ECCO_transport_lines"))
    ntr=length(list_trsp)
 
    trsp = SharedArray{Float64}(ntr,nr,nt)
    @sync @distributed for t in 1:nt
        comp_trsp(sol,nam,trsp,t; pth_in=pth_in, pth_out=pth_out)
    end
    
    trsp=[(nam=list_trsp[itr],val=trsp[itr,:,:]) for itr=1:ntr]
    save_object(joinpath(pth_out,calc*".jld2"),collect(trsp))
	"Done with transports"
end

function main_function(P)
    (; pth_in, pth_out, list_steps, nt, calc, nam, kk, sol) = P

    println("starting calc,sol,nam=$(calc),$(sol),$(nam) ...")
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
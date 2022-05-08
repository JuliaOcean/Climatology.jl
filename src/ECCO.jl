
module ECCO

using Pkg
import OceanStateEstimation: pkg_pth

"""
    ECCO.standard_analysis_setup(pth0)

Create a run folder where folder `pth0` will be linked. 

Folder `pth0` (of type `String`) should be the path to the user's ECCO data folder.

```
using OceanStateEstimation, Pkg

pth0="nctiles_monthly" #edit path as needed
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

#temporary fix:
import OceanStateEstimation: ECCO_helper_functions
γ,Γ,LC=ECCO_helper_functions.GridLoad_Main()

function path_etc(pth0::String,sol0::String,calc::String,nam::String)
    sol="ECCOv4"*sol0*"_analysis"

    if sol0=="r1"||sol0=="r2"
        nt=240
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
        pth_tmp=joinpath(pth_out,nam*"_"*calc)
    else
        pth_tmp=joinpath(pth_out,calc)
    end    

    return pth_in,pth_out,pth_tmp,nt,list_steps
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

    return Tx,Ty
end

function read_monthly_BSF(sol,nam,t,list_steps; pth_in="", pth_out="")
    U=read_monthly(sol,"UVELMASS",t,list_steps; pth_in=pth_in, pth_out=pth_out)
    V=read_monthly(sol,"VVELMASS",t,list_steps; pth_in=pth_in, pth_out=pth_out)
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
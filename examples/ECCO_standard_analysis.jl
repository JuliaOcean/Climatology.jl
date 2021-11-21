
"""
List of variables derived in this notebook:

- climatologies
- global means
- zonal means
- geographic maps
- transect transports
- MOC, MHT

Sample workflow:

```
#@everywhere nam="MXLDEPTH"
#@everywhere nam="THETA"
@everywhere calc="trsp"
@everywhere sol0="r2"
include("ECCO_standard_analysis.jl")
```
"""

@everywhere using MeshArrays, MITgcmTools, NCTiles
@everywhere using JLD2, UUIDs, Unitful
@everywhere using Distributed, SharedArrays

## Preliminary Steps

#STATE/state_3d_set1.0000241020.meta
#    'THETA   ' 'SALT    ' 'DRHODR  '
#TRSP/trsp_3d_set1.0000241020.meta
#    'UVELMASS' 'VVELMASS' 'WVELMASS' 'GM_PsiX ' 'GM_PsiY '
#TRSP/trsp_3d_set3.0000241020.meta
#    'DFxE_TH ' 'DFyE_TH ' 'ADVx_TH ' 'ADVy_TH ' 'DFxE_SLT' 'DFyE_SLT' 'ADVx_SLT' 'ADVy_SLT'

@everywhere sol="ECCOv4"*sol0*"_analysis"
@everywhere γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
if sol0!=="r5"
  @everywhere pth_in="ECCOv4"*sol0*"/nctiles_monthly/"
else
  @everywhere pth_in="ECCOv4"*sol0*"/diags/"
end
@everywhere pth_out=joinpath("ECCOv4"*sol0*"/",sol)
#@everywhere pth_out=joinpath(tempdir(),"ECCO_diags",sol)

!isdir(pth_out) ? mkdir(pth_out) : nothing

@everywhere function list_time_steps(pth_in)
    list=readdir(joinpath(pth_in,"STATE"))
    list=list[findall(occursin.(Ref("state_3d_set1"),list))]
    list=list[findall(occursin.(Ref("data"),list))]
    [list[i][15:end] for i in 1:length(list)]
end

if sol0=="r5"
 @everywhere list=list_time_steps(pth_in)
else
 @everywhere list=[]
end

@everywhere function get_nr_nt(pth_in,nam)
    nct_path=joinpath(pth_in,nam)
    lst=readdir(nct_path)
    lst=lst[findall(occursin.(Ref(".nc"),lst))]
    fil1=joinpath(nct_path,lst[1])
    ds=NCTiles.NCDatasets.Dataset(fil1)
    siz=size(ds[nam])
    siz[end-1],siz[end]
end

#@everywhere (nr,nt)=get_nr_nt(pth_in,"THETA")
@everywhere if sol=="ECCOv4r2_analysis"
    nt=240
elseif sol=="ECCOv4r3_analysis"
    nt=288
elseif sol=="ECCOv4r4_analysis"
    nt=312
elseif sol=="ECCOv4r5_analysis"
    nt=336
end
@everywhere nr=50

@everywhere function GridLoadVar(nam,γ)
    pc=fill(0.5,2); pg=fill(0.0,2); pu=[0.,0.5]; pv=[0.5,0.];
    list_n=("XC","XG","YC","YG","RAC","RAW","RAS","RAZ","DXC","DXG","DYC","DYG","Depth","AngleCS","AngleSN")
    list_u=(u"°",u"°",u"°",u"°",u"m^2",u"m^2",u"m^2",u"m^2",u"m",u"m",u"m",u"m",u"m",1.0,1.0)
    list_p=(pc,pg,pc,pg,pc,pu,pv,pg,pu,pv,pv,pu,pc,pc,pc)
    #
    list3d_n=("hFacC","hFacS","hFacW");
    list3d_u=(1.0,1.0,1.0)
    list3d_p=(fill(0.5,3),[0.,0.5,0.5],[0.5,0.,0.5])
    #
    list1d_n=("DRC","DRF","RC","RF")

    if sum(nam.==list_n)==1
        ii=findall(nam.==list_n)[1]
        m=varmeta(list_u[ii],list_p[ii],missing,list_n[ii],list_n[ii])
        tmp1=γ.read(joinpath(γ.path,list_n[ii]*".data"),MeshArray(γ,γ.ioPrec;meta=m))
    elseif sum(nam.==list1d_n)==1
        fil=joinpath(γ.path,nam*".data")
        γ.ioPrec==Float64 ? reclen=8 : reclen=4
        n3=Int64(stat(fil).size/reclen)

        fid = open(fil)
        tmp1 = Array{γ.ioPrec,1}(undef,n3)
        read!(fid,tmp1)
        tmp1 = hton.(tmp1)
    elseif sum(nam.==list3d_n)==1
        ii=findall(nam.==list3d_n)[1]
        m=varmeta(list3d_u[ii],list3d_p[ii],missing,list3d_n[ii],list3d_n[ii]);
        tmp1=γ.read(joinpath(γ.path,list3d_n[ii]*".data"),MeshArray(γ,γ.ioPrec,nr;meta=m))
    else
        tmp1=missing
    end
end

@everywhere XC=GridLoadVar("XC",γ)
@everywhere YC=GridLoadVar("YC",γ)
@everywhere RAC=GridLoadVar("RAC",γ)
@everywhere hFacC=GridLoadVar("hFacC",γ)
@everywhere DRF=GridLoadVar("DRF",γ)

@everywhere nansum(x) = sum(filter(!isnan,x))
@everywhere nansum(x,y) = mapslices(nansum,x,dims=y)

@everywhere mskC=hFacC./hFacC

tmp=hFacC./hFacC
tmp=[nansum(tmp[i,j].*RAC[i]) for j in 1:nr, i in eachindex(RAC)]
tot_RAC=nansum(tmp,2)

tmp=[nansum(hFacC[i,j].*RAC[i].*DRF[j]) for j in 1:nr, i in eachindex(RAC)]
tot_VOL=nansum(tmp,2)

## select variable and calculation

@everywhere pth_tmp=joinpath(pth_out,nam*"_"*calc)
!isdir(pth_tmp) ? mkdir(pth_tmp) : nothing

## generic read function

@everywhere function read_monthly(sol,nam,t,list)
    var_list3d=("THETA","SALT","UVELMASS","VVELMASS",
                "ADVx_TH","ADVy_TH","DFxE_TH","DFyE_TH")
    mdsio_list3d=("STATE/state_3d_set1","STATE/state_3d_set1",
        "TRSP/trsp_3d_set1","TRSP/trsp_3d_set1","TRSP/trsp_3d_set2",
        "TRSP/trsp_3d_set2","TRSP/trsp_3d_set2","TRSP/trsp_3d_set2")

    if (sol=="ECCOv4r2_analysis")||(sol=="ECCOv4r3_analysis")
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
        meta=read_meta(joinpath(pth_in,"STATE/state_3d_set1.0000241020.meta"))
        kk=findall(vec(meta.fldList).==nam)[1]
        tmp=read_mdsio(joinpath(pth_in,"STATE/state_3d_set1."*list[t]))[:,:,:,kk]
        tmp=mskC*read(tmp,γ)     
    end
end

## climatological mean

if calc=="clim"
    @distributed for m in 1:12
        nm=length(m:12:nt)
        tmp=1/nm*read_monthly(sol,nam,m,list)
        for t in m+12:12:nt
            tmp=tmp+1/nm*read_monthly(sol,nam,t,list)
        end
        save_object(joinpath(pth_tmp,"$m.jld2"),tmp)
    end
end

## global mean

if (calc=="glo2d")||(calc=="glo3d")
    glo = SharedArray{Float64}(nr,nt)

    @sync @distributed for t in 1:nt

        tmp=read_monthly(sol,nam,t,list)
        if calc=="glo2d"
            tmp=[nansum(tmp[i,j].*RAC[i]) for j in 1:nr, i in eachindex(RAC)]
        else
            tmp=[nansum(tmp[i,j].*hFacC[i,j].*RAC[i]*DRF[j]) for j in 1:nr, i in eachindex(RAC)]
        end
        glo[:,t]=nansum(tmp,2)
    end

    if calc=="glo2d"
        tmp=[glo[r,t]/tot_RAC[r] for t in 1:nt, r in 1:nr]
    else
        tmp=[nansum(glo[:,t])/nansum(tot_VOL) for t in 1:nt]
    end
    save_object(joinpath(pth_tmp,calc*".jld2"),tmp)
end

##

#same as LatitudeCircles? split out of LatitudeCircles?
@everywhere function edge_mask(mskCint::MeshArray)
    mskCint=1.0*mskCint

    #treat the case of blank tiles:
    #mskCint[findall(RAC.==0)].=NaN
    
    mskCplus=exchange(mskCint)

    #edge tracer mask
    mskCedge=similar(mskCint)
    for i in eachindex(mskCedge)
        tmp1=mskCplus[i]
        tmp2=tmp1[2:end-1,1:end-2]+tmp1[2:end-1,3:end]+
            tmp1[1:end-2,2:end-1]+tmp1[3:end,2:end-1]
        mskCedge[i]=1.0*(tmp2.>0).*(tmp1[2:end-1,2:end-1].==0.0)
    end

    #edge velocity mask:
    mskWedge=similar(mskCint)
    mskSedge=similar(mskCint)
    for i in eachindex(mskCedge)
        mskWedge[i]=mskCplus[i][2:end-1,2:end-1] - mskCplus[i][1:end-2,2:end-1]
        mskSedge[i]=mskCplus[i][2:end-1,2:end-1] - mskCplus[i][2:end-1,1:end-2]
    end

    #treat the case of blank tiles:
    #mskCedge[findall(isnan.(mskCedge))].=0.0
    #mskWedge[findall(isnan.(mskWedge))].=0.0
    #mskSedge[findall(isnan.(mskSedge))].=0.0

    return mskCedge,mskWedge,mskSedge
end

##

if (calc=="zonmean")||(calc=="zonmean2d")
    dlat=2.0
    lats=(-90+dlat/2:dlat:90-dlat/2)
    save_object(joinpath(pth_tmp,calc*"_lats.jld2"),lats)
    nl=length(lats)

    msk0 = SharedArray{Float64}(γ.ioSize...,nl)
    zm0 = SharedArray{Float64}(nl,nr)
    tmp0=hFacC./hFacC
    @sync @distributed for l in 1:nl
        la0=lats[l]-dlat/2
        la1=lats[l]+dlat/2
        if la1<0.0
            msk=1.0*(YC.>=la0)*(YC.<la1)
        elseif la0>0.0
            msk=1.0*(YC.>la0)*(YC.<=la1)
        else
            msk=1.0*(YC.>=la0)*(YC.<=la1)
        end
        msk[findall(msk.==0.0)].=NaN;
        msk0[:,:,l]=write(msk*RAC)

        tmp2=[nansum(tmp0[i,j].*msk[i].*RAC[i]) for j in 1:nr, i in eachindex(RAC)]
        zm0[l,:]=1.0 ./nansum(tmp2,2)
    end
    save_object(joinpath(pth_tmp,calc*"_zm0.jld2"),zm0)
    save_object(joinpath(pth_tmp,calc*"_msk0.jld2"),msk0)

    #to speed up main loop, reuse:
    #- precomputed msk*RAC once and for all
    #- precomputed 1.0./nansum(tmp2,2)

    @everywhere msk0=load(joinpath(pth_tmp,calc*"_msk0.jld2"),"single_stored_object")
    @everywhere zm0=load(joinpath(pth_tmp,calc*"_zm0.jld2"),"single_stored_object")

    if (calc=="zonmean")
        zm = SharedArray{Float64}(nl,nr,nt)
        @sync @distributed for t in 1:nt
            tmp=read_monthly(sol,nam,t,list)
            for l in 1:nl
                mskrac=read(msk0[:,:,l],γ)
                tmp1=[nansum(tmp[i,j].*mskrac[i]) for j in 1:nr, i in eachindex(RAC)]
                zm[l,:,t]=nansum(tmp1,2).*zm0[l,:]
            end
        end
    else
        zm = SharedArray{Float64}(nl,nt)
        @sync @distributed for t in 1:nt
            if nam=="SSH"
                #tmp=read_monthly_ssh(sol,t,list)
                ETAN=read_monthly(sol,"ETAN",t,list)
                sIceLoad=read_monthly(sol,"sIceLoad",t,list)
                tmp=(ETAN+sIceLoad/1029.0)*mskC[:,1]
            else
                tmp=read_monthly(sol,nam,t,list)
            end
            for l in 1:nl
                mskrac=read(msk0[:,:,l],γ)
                tmp1=[nansum(tmp[i].*mskrac[i]) for i in eachindex(RAC)]
                zm[l,t]=nansum(tmp1)*zm0[l,1]
            end
        end
    end
    save_object(joinpath(pth_tmp,calc*".jld2"),zm)

    #levs=collect(-1.0:-1.0:-50.0)
    #heatmap(lats,kk,yearN-year1,colorrange=(-1.0,1.0),colormap="jet")
end

##

if (calc=="overturn")
    @everywhere Γ=GridLoad(γ,option="full")
	@everywhere LC=LatitudeCircles(-89.0:89.0,Γ)
    @everywhere nl=length(LC)
	
    ov = SharedArray{Float64}(nl,nr,nt)
    @sync @distributed for t in 1:nt
        U=read_monthly(sol,"UVELMASS",t,list)
        V=read_monthly(sol,"VVELMASS",t,list)
        (Utr,Vtr)=UVtoTransport(U,V,Γ)
        #integrate across latitude circles
        for z=1:nr
            UV=Dict("U"=>Utr[:,z],"V"=>Vtr[:,z],"dimensions"=>["x","y"])
            [ov[l,z,t]=ThroughFlow(UV,LC[l],Γ) for l=1:nl]
        end
        #integrate from bottom
        ov[:,:,t]=reverse(cumsum(reverse(ov[:,:,t],dims=2),dims=2),dims=2);
    end
    
    save_object(joinpath(pth_tmp,calc*".jld2"),ov)
	"Done with overturning"
end

if (calc=="MHT")
    @everywhere Γ=GridLoad(γ,option="full")
    @everywhere LC=LatitudeCircles(-89.0:89.0,Γ)
    @everywhere nl=length(LC)

    MHT = SharedArray{Float64}(nl,nt)
    @sync @distributed for t in 1:nt
        U=read_monthly(sol,"ADVx_TH",t,list)
        V=read_monthly(sol,"ADVy_TH",t,list)
        U=U+read_monthly(sol,"DFxE_TH",t,list)
        V=V+read_monthly(sol,"DFyE_TH",t,list)

        [U[i][findall(isnan.(U[i]))].=0.0 for i in eachindex(U)]
        [V[i][findall(isnan.(V[i]))].=0.0 for i in eachindex(V)]

        Tx=0.0*U[:,1]
        Ty=0.0*V[:,1]
        [Tx=Tx+U[:,z] for z=1:nr]
        [Ty=Ty+V[:,z] for z=1:nr]

        UV=Dict("U"=>Tx,"V"=>Ty,"dimensions"=>["x","y"])
        [MHT[l,t]=1e-15*4e6*ThroughFlow(UV,LC[l],Γ) for l=1:nl]
    end
    
    save_object(joinpath(pth_tmp,calc*".jld2"),MHT)
	"Done with MHT"
end

##

if (calc=="trsp")
    @everywhere Γ=GridLoad(γ,option="full")
    @everywhere pth_trsp=joinpath(tempdir(),"ECCO_transport_lines")
    @everywhere list_trsp=readdir(pth_trsp)
    @everywhere ntr=length(list_trsp)
    @everywhere TR=[load(joinpath(pth_trsp,list_trsp[itr])) for itr in 1:ntr]
    @everywhere TR=MeshArrays.Dict_to_NamedTuple.(TR)
    
    trsp = SharedArray{Float64}(ntr,nr,nt)
    @sync @distributed for t in 1:nt
        U=read_monthly(sol,"UVELMASS",t,list)
        V=read_monthly(sol,"VVELMASS",t,list)
        (Utr,Vtr)=UVtoTransport(U,V,Γ)
        #integrate across transport lines
        for z=1:nr
            UV=Dict("U"=>Utr[:,z],"V"=>Vtr[:,z],"dimensions"=>["x","y"])
            [trsp[itr,z,t]=ThroughFlow(UV,TR[itr],Γ) for itr=1:ntr]
        end
    end
    
    trsp=[(nam=list_trsp[itr],val=trsp[itr,:,:]) for itr=1:ntr]
    save_object(joinpath(pth_tmp,calc*".jld2"),trsp)
	"Done with transports"
end

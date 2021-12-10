
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
@everywhere using JLD2, UUIDs, Unitful, Printf
@everywhere using Distributed, SharedArrays

@everywhere include("ECCO_helper_functions.jl")
@everywhere γ,Γ=GridLoadSome()
@everywhere nr=length(Γ.DRF)
@everywhere LC=LatitudeCircles(-89.0:89.0,Γ)
@everywhere nl=length(LC)
@everywhere list_trsp,msk_trsp=GridLoad_TR()
@everywhere ntr=length(msk_trsp)

## Preliminary Steps

@everywhere sol="ECCOv4"*sol0*"_analysis"

@everywhere if sol=="ECCOv4r2_analysis"
    nt=240
elseif sol=="ECCOv4r3_analysis"
    nt=288
elseif sol=="ECCOv4r4_analysis"
    nt=312
elseif sol=="ECCOv4r5_analysis"
    nt=336
end

if sol0!=="r5"
  @everywhere pth_in="ECCOv4"*sol0*"/nctiles_monthly/"
else
  @everywhere pth_in="ECCOv4"*sol0*"/diags/"
end
@everywhere list_steps=list_time_steps(pth_in)

@everywhere pth_out=joinpath("ECCO_diags",sol)
!isdir(pth_out) ? mkdir(pth_out) : nothing

## select variable and calculation

@everywhere if sum(calc.==("overturn","MHT","trsp"))==0
    pth_tmp=joinpath(pth_out,nam*"_"*calc)
else
    pth_tmp=joinpath(pth_out,calc)
end    
!isdir(pth_tmp) ? mkdir(pth_tmp) : nothing

## climatological mean

if calc=="clim"
    tmp_s1 = SharedArray{Float64}(γ.ioSize...,12)
    tmp_s2 = SharedArray{Float64}(γ.ioSize...,12)
    tmp_m = SharedArray{Float64}(γ.ioSize...,12)

    tmp=read_monthly(sol,nam,1,list_steps)
    ndims(tmp)>1 ? nz=size(tmp,2) : nz=1
    nz==1 ? kk=1 : nothing
    nz>1 ? suff=Printf.@sprintf("_k%02d",kk) : suff=""

    @sync @distributed for m in 1:12
        nm=length(m:12:nt)
        tmp_m[:,:,m].=0.0
        tmp_s1[:,:,m].=0.0
        tmp_s2[:,:,m].=0.0
        for t in m:12:nt
            tmp=read_monthly(sol,nam,t,list_steps)
            ndims(tmp)>1 ? tmp=tmp[:,kk] : nothing
            tmp_m[:,:,m]=tmp_m[:,:,m]+1.0/nm*γ.write(tmp)
            tmp_s1[:,:,m]=tmp_s1[:,:,m]+γ.write(tmp)
            tmp_s2[:,:,m]=tmp_s2[:,:,m]+γ.write(tmp).^2
        end
        #save_object(joinpath(pth_tmp,nam*suff*Printf.@sprintf("_m%02d.jld2",m)),tmp_m)
    end

    tmp0=read(tmp_m[:],γ)

    tmp=1.0/nt*sum(tmp_s1,dims=3)
    tmp1=read(tmp[:],Γ.XC)

    tmp=1/nt*sum(tmp_s2,dims=3)-tmp.^2
    tmp[findall(tmp.<0.0)].=0.0
    tmp=sqrt.(nt/(nt-1)*tmp)
    tmp2=read(tmp[:],Γ.XC)

    fil_out=joinpath(pth_tmp,nam*suff*".jld2")
    save(fil_out,"mean",tmp1,"std",tmp2,"mon",tmp0)
end

## global mean

if (calc=="glo2d")||(calc=="glo3d")
    
    @everywhere function comp_glo(glo,t)
        tmp=read_monthly(sol,nam,t,list_steps)
        if calc=="glo2d"
            tmp=[nansum(tmp[i,j].*Γ.RAC[i]) for j in 1:nr, i in eachindex(Γ.RAC)]
        else
            tmp=[nansum(tmp[i,j].*Γ.hFacC[i,j].*Γ.RAC[i]*Γ.DRF[j]) for j in 1:nr, i in eachindex(Γ.RAC)]
        end
        glo[:,t]=nansum(tmp,2)
    end
    
    glo = SharedArray{Float64}(nr,nt)
    @sync @distributed for t in 1:nt
        comp_glo(glo,t)
    end

    if calc=="glo2d"
        tmp=[glo[r,t]/Γ.tot_RAC[r] for t in 1:nt, r in 1:nr]
    else
        tmp=[nansum(glo[:,t])/nansum(Γ.tot_VOL) for t in 1:nt]
    end
    save_object(joinpath(pth_tmp,calc*".jld2"),tmp)
end

##

if (calc=="zonmean")||(calc=="zonmean2d")
    dlat=2.0
    lats=(-90+dlat/2:dlat:90-dlat/2)
    save_object(joinpath(pth_tmp,calc*"_lats.jld2"),lats)
    nl=length(lats)

    @everywhere function comp_msk0(msk0,zm0,l)
        lats=load(joinpath(pth_tmp,calc*"_lats.jld2"),"single_stored_object")
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

    msk0 = SharedArray{Float64}(γ.ioSize...,nl)
    zm0 = SharedArray{Float64}(nl,nr)
    @sync @distributed for l in 1:nl
        comp_msk0(msk0,zm0,l)
    end
    save_object(joinpath(pth_tmp,calc*"_zm0.jld2"),zm0)
    save_object(joinpath(pth_tmp,calc*"_msk0.jld2"),msk0)

    #to speed up main loop, reuse:
    #- precomputed msk*RAC once and for all
    #- precomputed 1.0./nansum(tmp2,2)

    @everywhere msk0=load(joinpath(pth_tmp,calc*"_msk0.jld2"),"single_stored_object")
    @everywhere zm0=load(joinpath(pth_tmp,calc*"_zm0.jld2"),"single_stored_object")

    @everywhere function comp_zonmean(zm,t)
        lats=load(joinpath(pth_tmp,calc*"_lats.jld2"),"single_stored_object")
        nl=length(lats)
        tmp=read_monthly(sol,nam,t,list_steps)
        for l in 1:nl
            mskrac=read(msk0[:,:,l],γ)
            tmp1=[nansum(tmp[i,j].*mskrac[i]) for j in 1:nr, i in eachindex(Γ.RAC)]
            zm[l,:,t]=nansum(tmp1,2).*zm0[l,:]
        end
    end

    @everywhere function comp_zonmean2d(zm,t)
        lats=load(joinpath(pth_tmp,calc*"_lats.jld2"),"single_stored_object")
        nl=length(lats)
        tmp=read_monthly(sol,nam,t,list_steps)
        for l in 1:nl
            mskrac=read(msk0[:,:,l],γ)
            tmp1=[nansum(tmp[i].*mskrac[i]) for i in eachindex(Γ.RAC)]
            zm[l,t]=nansum(tmp1)*zm0[l,1]
        end
    end

    if (calc=="zonmean")
        zm = SharedArray{Float64}(nl,nr,nt)
        @sync @distributed for t in 1:nt
            comp_zonmean(zm,t)
        end
    else
        zm = SharedArray{Float64}(nl,nt)
        @sync @distributed for t in 1:nt
            comp_zonmean2d(zm,t)
        end
    end
    save_object(joinpath(pth_tmp,calc*".jld2"),zm)

    #levs=collect(-1.0:-1.0:-50.0)
    #heatmap(lats,kk,yearN-year1,colorrange=(-1.0,1.0),colormap="jet")
end

##

if (calc=="overturn")
    @everywhere function comp_overturn(ov,t)
        U=read_monthly(sol,"UVELMASS",t,list_steps)
        V=read_monthly(sol,"VVELMASS",t,list_steps)
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

    ov = SharedArray{Float64}(nl,nr,nt)
    @sync @distributed for t in 1:nt
        comp_overturn(ov,t)
    end
    
    save_object(joinpath(pth_tmp,calc*".jld2"),ov)
	"Done with overturning"
end

if (calc=="MHT")
    @everywhere function comp_MHT(MHT,t)
        U=read_monthly(sol,"ADVx_TH",t,list_steps)
        V=read_monthly(sol,"ADVy_TH",t,list_steps)
        U=U+read_monthly(sol,"DFxE_TH",t,list_steps)
        V=V+read_monthly(sol,"DFyE_TH",t,list_steps)

        [U[i][findall(isnan.(U[i]))].=0.0 for i in eachindex(U)]
        [V[i][findall(isnan.(V[i]))].=0.0 for i in eachindex(V)]

        Tx=0.0*U[:,1]
        Ty=0.0*V[:,1]
        [Tx=Tx+U[:,z] for z=1:nr]
        [Ty=Ty+V[:,z] for z=1:nr]

        UV=Dict("U"=>Tx,"V"=>Ty,"dimensions"=>["x","y"])
        [MHT[l,t]=1e-15*4e6*ThroughFlow(UV,LC[l],Γ) for l=1:nl]
    end

    MHT = SharedArray{Float64}(nl,nt)
    @sync @distributed for t in 1:nt
        comp_MHT(MHT,t)
    end
    
    save_object(joinpath(pth_tmp,calc*".jld2"),MHT)
	"Done with MHT"
end

##

if (calc=="trsp")
    @everywhere function comp_trsp(trsp,t)
        U=read_monthly(sol,"UVELMASS",t,list_steps)
        V=read_monthly(sol,"VVELMASS",t,list_steps)
        (Utr,Vtr)=UVtoTransport(U,V,Γ)
        #integrate across transport lines
        for z=1:nr
            UV=Dict("U"=>Utr[:,z],"V"=>Vtr[:,z],"dimensions"=>["x","y"])
            [trsp[itr,z,t]=ThroughFlow(UV,msk_trsp[itr],Γ) for itr=1:ntr]
        end
    end

    trsp = SharedArray{Float64}(ntr,nr,nt)
    @sync @distributed for t in 1:nt
        comp_trsp(trsp,t)
    end
    
    trsp=[(nam=list_trsp[itr],val=trsp[itr,:,:]) for itr=1:ntr]
    save_object(joinpath(pth_tmp,calc*".jld2"),trsp)
	"Done with transports"
end

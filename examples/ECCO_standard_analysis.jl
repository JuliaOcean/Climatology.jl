
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
@everywhere include("ECCO_helper_functions.jl")

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
    @distributed for m in 1:12
        nm=length(m:12:nt)
        tmp=1/nm*read_monthly(sol,nam,m,list_steps)
        for t in m+12:12:nt
            tmp=tmp+1/nm*read_monthly(sol,nam,t,list_steps)
        end
        save_object(joinpath(pth_tmp,"$m.jld2"),tmp)
    end
end

## global mean

if (calc=="glo2d")||(calc=="glo3d")
    glo = SharedArray{Float64}(nr,nt)

    @sync @distributed for t in 1:nt

        tmp=read_monthly(sol,nam,t,list_steps)
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
            tmp=read_monthly(sol,nam,t,list_steps)
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
                ETAN=read_monthly(sol,"ETAN",t,list_steps)
                sIceLoad=read_monthly(sol,"sIceLoad",t,list_steps)
                tmp=(ETAN+sIceLoad/1029.0)*mskC[:,1]
            else
                tmp=read_monthly(sol,nam,t,list_steps)
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
        U=read_monthly(sol,"UVELMASS",t,list_steps)
        V=read_monthly(sol,"VVELMASS",t,list_steps)
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
        U=read_monthly(sol,"UVELMASS",t,list_steps)
        V=read_monthly(sol,"VVELMASS",t,list_steps)
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
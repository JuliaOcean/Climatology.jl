
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

@everywhere function comp_clim(tmp_m,tmp_s1,tmp_s2,kk,m)
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
end

function main_clim(sol,nam,kk=1)

    tmp_s1 = SharedArray{Float64}(γ.ioSize...,12)
    tmp_s2 = SharedArray{Float64}(γ.ioSize...,12)
    tmp_m = SharedArray{Float64}(γ.ioSize...,12)

    tmp=read_monthly(sol,nam,1,list_steps)
    ndims(tmp)>1 ? nz=size(tmp,2) : nz=1
    nz==1 ? kk=1 : nothing
    nz>1 ? suff=Printf.@sprintf("_k%02d",kk) : suff=""

    @sync @distributed for m in 1:12
        comp_clim(tmp_m,tmp_s1,tmp_s2,kk,m)
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

    return true
end

if calc=="clim"
    main_clim(sol,nam,kk)
end

## global mean

@everywhere function comp_glo(glo,t)
    tmp=read_monthly(sol,nam,t,list_steps)
    if calc=="glo2d"
        tmp=[nansum(tmp[i,j].*Γ.RAC[i]) for j in 1:nr, i in eachindex(Γ.RAC)]
    else
        tmp=[nansum(tmp[i,j].*Γ.hFacC[i,j].*Γ.RAC[i]*Γ.DRF[j]) for j in 1:nr, i in eachindex(Γ.RAC)]
    end
    glo[:,t]=nansum(tmp,2)
end
    
function main_glo(sol,nam)

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

if (calc=="glo2d")||(calc=="glo3d")
    main_glo(sol,nam)
end

##

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

function main_zonmean(sol,nam)
    dlat=2.0
    lats=(-90+dlat/2:dlat:90-dlat/2)
    save_object(joinpath(pth_tmp,calc*"_lats.jld2"),lats)
    nl=length(lats)

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

    return true
end

if (calc=="zonmean")||(calc=="zonmean2d")
    main_zonmean(sol,nam)
end

##

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

function main_overturn(sol,nam)  
    ov = SharedArray{Float64}(nl,nr,nt)
    @sync @distributed for t in 1:nt
        comp_overturn(ov,t)
    end
    
    save_object(joinpath(pth_tmp,calc*".jld2"),ov)
	"Done with overturning"
end

if (calc=="overturn")
    main_overturn(sol,nam)
end

##

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

function main_MHT(sol,nam)  
    MHT = SharedArray{Float64}(nl,nt)
    @sync @distributed for t in 1:nt
        comp_MHT(MHT,t)
    end
    save_object(joinpath(pth_tmp,calc*".jld2"),MHT)
	"Done with MHT"
end

if (calc=="MHT")
    main_MHT(sol,nam)
end

##

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

function main_trsp(sol,nam)  
    trsp = SharedArray{Float64}(ntr,nr,nt)
    @sync @distributed for t in 1:nt
        comp_trsp(trsp,t)
    end
    
    trsp=[(nam=list_trsp[itr],val=trsp[itr,:,:]) for itr=1:ntr]
    save_object(joinpath(pth_tmp,calc*".jld2"),trsp)
	"Done with transports"
end

if (calc=="trsp")
    main_trsp(sol,nam)
end

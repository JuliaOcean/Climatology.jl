
"""
object: load model and satellite climatologies and cross-compare
date: 2021/10/28
author: Gaël Forget

- CBIOMES-global-alpha-climatology.nc 
- examples/gridded_geospatial_montly_clim_360_720_ver_0_2.nc
- examples/gridded_darwin_montly_clim_360_720_ver_0_2_6.nc 
"""

using OceanStateEstimation, JLD2, NCTiles, Statistics
import CairoMakie as Mkie

##

selectPlot=1
println("\n selectPlot=$(selectPlot)\n")

##

fil_out=joinpath(CBIOMESclim_path,"CBIOMES-global-alpha-climatology.nc")
nc=NCTiles.NCDataset(fil_out,"r")
lon=nc["lon"][:]
lat=nc["lat"][:]
uni=nc["Chl"].attrib["units"]
Chl_from_Mod=nc["Chl"][:]

##

Chl_from_Rrs=load(joinpath(tempdir(),"Chl_from_Rrs_clim.jld2"))["clim_Chl"]
fil_sat="examples/gridded_geospatial_montly_clim_360_720_ver_0_2.nc"
Chl_from_Sat=NCTiles.NCDataset(fil_sat,"r")["Chl"][:]
#fil_="examples/gridded_darwin_montly_clim_360_720_ver_0_2_6.nc"
#Chl_from_Bror=NCTiles.NCDataset(fil_Bror,"r")["Chl"][:]

## plotting functions

function PacificMapFig(val,lev,ttl)
	fig = Mkie.Figure(size = (600,400), backgroundcolor = :grey95, fontsize=12)
	ax = Mkie.Axis(fig[1,1], title=ttl*" units: $(uni)",xlabel="longitude",ylabel="latitude")
	_,hm1=PacificMap!(ax,val,lev,ttl)
	Mkie.Colorbar(fig[1,2], hm1, height = Mkie.Relative(0.65))
	fig
end

function PacificMap!(ax,val,lev,ttl)
	hm1=Mkie.contourf!(ax,[lon.-360.0;lon],lat,[val;val], colormap = Mkie.cgrad(:Spectral, length(lev), categorical = true), 
    tickfont = (4, :black), levels = lev)
	#
	lon1=[-180.0 -120.0 -120.0 -180.0 -180.0]
	lat1=[-20.0 -20.0 50.0 50.0 -20.0]
	Mkie.lines!(ax,lon1[:],lat1[:],color=:gray25,linewidth=3.0)
	#
	Mkie.xlims!(ax, (-200.0, -100.0)); Mkie.ylims!(ax, (-30.0, 60.0))
	#
	ax,hm1
end

## helper functions

function regional_mean(val)
    val[findall(ismissing.(val))].=NaN
    lon1=[-180.0 -120.0 -120.0 -180.0 -180.0]
	lat1=[-20.0 -20.0 50.0 50.0 -20.0]

	i1=findall( (lon.>-180.0).*(lon.<=-120.0) )
	j1=findall( (lat.>-20.0).*(lat.<=50.0) )

	val_region=val[i1,j1]
	msk_region=1.0 .+ 0.0*val_region
	val_region[findall(isnan.(val_region))].=0.0
	msk_region[findall(isnan.(msk_region))].=0.0

	sum(val_region,dims=(1,2))[1]./sum(msk_region,dims=(1,2))[1]
end

function regional_median(val)
    val[findall(ismissing.(val))].=NaN
    lon1=[-180.0 -120.0 -120.0 -180.0 -180.0]
	lat1=[-20.0 -20.0 50.0 50.0 -20.0]

	i1=findall( (lon.>-180.0).*(lon.<=-120.0) )
	j1=findall( (lat.>-20.0).*(lat.<=50.0) )

	val_region=val[i1,j1]
	median(val_region[findall((!isnan).(val_region))])
end

##

fig=[Mkie.Figure()]
if selectPlot==1 #annual mean difference
	x=dropdims(mean(Chl_from_Mod,dims=3)-mean(Chl_from_Sat,dims=3),dims=3); ttl="Chl_from_Mod-Chl_from_Sat, annual mean"
	##x=dropdims(mean(Chl_from_Rrs,dims=3)-mean(Chl_from_Sat,dims=3),dims=3); ttl="Chl_from_Rrs-Chl_from_Sat, annual mean"
	##x=dropdims(mean(Chl_from_Mod,dims=3)-mean(Chl_from_Rrs,dims=3),dims=3); ttl="Chl_from_Mod-Chl_from_Rrs, annual mean"
	lev=2.0*(-0.1:0.02:0.1) .- regional_median(x)
	fig[1]=PacificMapFig(x,lev,ttl)
elseif selectPlot==2 #monthly mean difference
	m=1
	x=Chl_from_Mod[:,:,m]-Chl_from_Sat[:,:,m]; 
	ttl="Chl_from_Mod-Chl_from_Sat, month $m"
	lev=2.0*(-0.1:0.02:0.1) .- regional_mean(x)
	fig[1]=PacificMapFig(x,lev,ttl)
end

## see below for selectPlot = 3,4 

if selectPlot<3
	fig[1]
elseif selectPlot==3 #12 monthly misfit maps
	FIG = Mkie.Figure(size = (900,1200), backgroundcolor = :grey95, fontsize=12)
	f(x)=Int(ceil(x/3.0))
	g(x)=Int(x-3*(f(x)-1.0))
	for m in 1:12
		ax = Mkie.Axis(FIG[f(m),g(m)], title="month $m ; units: $(uni)",xlabel="longitude",ylabel="latitude")
		x=Chl_from_Mod[:,:,m]-Chl_from_Sat[:,:,m]
		#x=Chl_from_Rrs[:,:,m]-Chl_from_Sat[:,:,m]
		#x=Chl_from_Mod[:,:,m]-Chl_from_Rrs[:,:,m]
		lev=2.0*(-0.1:0.02:0.1) .- regional_median(x)
		PacificMap!(ax,x,lev,ttl)
	end
	FIG
elseif selectPlot==4 #12 monthly anomaly map
	FIG = Mkie.Figure(size = (900,1200), backgroundcolor = :grey95, fontsize=12)
	f(x)=Int(ceil(x/3.0))
	g(x)=Int(x-3*(f(x)-1.0))
	for m in 1:12
		ax = Mkie.Axis(FIG[f(m),g(m)], title="month $m ; units: $(uni)",xlabel="longitude",ylabel="latitude")
		x=dropdims(Chl_from_Mod[:,:,m]-mean(Chl_from_Mod,dims=3),dims=3); ttl="model anomalies"
		#x=dropdims(Chl_from_Sat[:,:,m]-mean(Chl_from_Sat,dims=3),dims=3); ttl="sat. anomalies"
		#x=dropdims(Chl_from_Rrs[:,:,m]-mean(Chl_from_Rrs,dims=3),dims=3); ttl="model(rrs) anomalies"
		lev=2.0*(-0.1:0.02:0.1) .- regional_median(x)
		PacificMap!(ax,x,lev,ttl)
	end
	FIG
end

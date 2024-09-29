
module ClimatologyMakieExt

	using Makie, Climatology
	import Climatology: Statistics, RollingFunctions, plot_examples, load
	import Climatology: ECCOdiag, SSTdiag, SeaLevelAnomaly
	import Statistics: mean
	import Makie: plot
	import RollingFunctions: runmean

	## 1. ECCO

	function plot(x::ECCOdiag)
		if !isempty(x.options)
			o=x.options
			if string(o.plot_type)=="ECCO_map"
				map(ECCO_procs.map(o.nammap,o.P,o.statmap,o.timemap,x.path))
			elseif string(o.plot_type)=="ECCO_TimeLat"
				nam=split(x.name,"_")[1]
				TimeLat(ECCO_procs.TimeLat(nam,x.path,o.year0,o.year1,o.cmap_fac,o.k,o.P); years_to_display=o.years_to_display)
			elseif string(o.plot_type)=="ECCO_TimeLatAnom"
				nam=split(x.name,"_")[1]
				TimeLat(ECCO_procs.TimeLatAnom(nam,x.path,o.year0,o.year1,o.cmap_fac,o.k,o.l0,o.l1,o.P); years_to_display=o.years_to_display)
			elseif string(o.plot_type)=="ECCO_DepthTime"
				nam=split(x.name,"_")[1]
				DepthTime(ECCO_procs.DepthTime(nam,x.path,o.facA,o.l,o.year0,o.year1,o.k0,o.k1,o.P); years_to_display=o.years_to_display)
			elseif string(o.plot_type)=="ECCO_GlobalMean"
				gl1=ECCO_procs.glo(x.path,x.name,o.k,o.year0,o.year1)
				glo(gl1,o.year0,o.year1; years_to_display=o.years_to_display)
			elseif x.name=="OHT"&&string(o.plot_type)=="ECCO_OHT1"
				OHT(x.path)
			elseif x.name=="overturn"&&string(o.plot_type)=="ECCO_Overturn1"
				figov1(x.path,o.kk,o.low1,o.year0,o.year1; years_to_display=o.years_to_display)
			elseif x.name=="overturn"&&string(o.plot_type)=="ECCO_Overturn2"
				figov2(x.path,o.grid)
			elseif x.name=="trsp"&&string(o.plot_type)=="ECCO_Transports"
				transport(o.namtrs,o.ncols,x.path,o.list_trsp,o.year0,o.year1,years_to_display=o.years_to_display)
			else
				println("unknown option (b)")	
			end
		else
			println("unknown option (a)")
		end
	end

	##

	to_range!(DD,levs::Tuple) = to_range!(DD,range(levs[1],levs[2],length=10))

	function to_range!(DD,levs)
		DD[findall(DD.<=levs[1])].=levs[1]+(levs[2]-levs[1])/100
		DD[findall(DD.>=levs[end])].=levs[end]-(levs[end]-levs[end-1])/100
	end

#	years_to_display=(1960,2023)
	years_to_display=(1980,2024)

	function axtr1(ax,namtr,pth_out,list_trsp,year0,year1;years_to_display=years_to_display)
		itr=findall(list_trsp.==namtr)[1]
		tmp=vec(load(ECCOdiag(path=pth_out,name="trsp")))[itr]
		
		nt=size(tmp.val,2)
		x=vec(0.5:nt)
	
		txt=tmp.nam[1:end-5]
		val=1e-6*vec(sum(tmp.val,dims=1)[:])
		valsmo = runmean(val, 12)
	
		x=vec(0.5:nt)
		x=year0 .+ x./12.0

		hm1=lines!(ax,x,val,label="ECCO estimate")
		valsmo[1:5].=NaN
		valsmo[end-4:end].=NaN
		lines!(ax,x,valsmo,linewidth=4.0,color=:red)
		xlims!(ax,years_to_display)
	end

	function transport(namtrs,ncols,pth_out,list_trsp,year0,year1;years_to_display=years_to_display)
		if ncols > 1
			fig1 = Figure(size = (2000,1000),markersize=0.1)
		else
			fig1 = Figure(size = (900,400),markersize=0.1)
		end
		for na in 1:length(namtrs)
			txt=namtrs[na]
			jj=div.(na,ncols,RoundUp)
			kk=na-(jj.-1)*ncols
			ax1 = Axis(fig1[jj,kk], title=" $txt (in Sv)",
				xticks=(year0:4:year1),ylabel="transport, in Sv")
			axtr1(ax1,namtrs[na],pth_out,list_trsp,year0,year1,years_to_display=years_to_display)
			#ylims!(ax1,rng)
		end
		fig1
	end

	function figov1(pth_out,kk,low1,year0,year1;years_to_display=years_to_display)
		tmp=-1e-6*load(ECCOdiag(path=pth_out,name="overturn"))
	
		nt=size(tmp,3)
		x=vec(0.5:nt)
		x=year0 .+ x./12.0
		lats=vec(-89.0:89.0)

		fig1 = Figure(size = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1],ylabel="Sv",
			title="Global Overturning, in Sv, at kk=$(kk)",
			xticks=(year0:4:year1))
		for ll in 115:10:145
			ov=tmp[ll,kk,:]
			ov=runmean(ov, 12)
			ov[1:5].=NaN
			ov[end-4:end].=NaN
			hm1=lines!(x,ov,label="$(lats[ll])N")
		end
		xlims!(ax1,years_to_display)
		ylims!(ax1,(5,20))
		low1!="auto" ? ylims!(ax1,(low1,20.0)) : nothing
		fig1[1, 2] = Legend(fig1, ax1, "estimate", framevisible = false)
	
		fig1
	end

	function figov2(pth_out,Γ; ClipToRange=true)
		tmp=-1e-6*load(ECCOdiag(path=pth_out,name="overturn"))
		ovmean=dropdims(mean(tmp[:,:,1:240],dims=3),dims=3)

		x=vec(-89.0:89.0); y=reverse(vec(Γ.RF[1:end-1])); #coordinate variables
		z=reverse(ovmean,dims=2); z[z.==0.0].=NaN

		levs=(-40.0:5.0:40.0)
		ClipToRange ? to_range!(z,levs) : nothing
	
		fig1 = Figure(size = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1], title="Meridional Overturning Streamfunction (in Sv, time mean)",
				xlabel="latitude",ylabel="depth (in m)")
		hm1=contourf!(ax1,x,y,z,levels=levs)
		Colorbar(fig1[1,2], hm1, height = Relative(0.65))
		fig1
	end

	function OHT(pth_out)
		tmp=load(ECCOdiag(path=pth_out,name="MHT"))
		MT=vec(mean(tmp[:,1:240],dims=2))

		x=vec(-89.0:89.0)
		fig1 = Figure(size = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1], title="Northward Heat Transport (in PW, time mean)",
			xticks=(-90.0:10.0:90.0),yticks=(-2.0:0.25:2.0),
			xlabel="latitude",ylabel="Transport (in PW)")
		hm1=lines!(x,MT)
		ylims!(ax1,(-2.0,2.0))
		fig1
	end

	function glo(gl1,year0,year1;years_to_display=years_to_display)
		ttl="Global Mean $(gl1.txt)"
		zlb=gl1.txt
		rng=gl1.rng

		if false
			fac=4e6*1.335*10^9*10^9/1e21
			ttl="Ocean Heat Uptake (Zetta-Joules)"
			zlb="Zetta-Joules"
			rng=(-100.0,300.0)
			y=fac*(gl1.y.-gl1.y[1])
		else
			y=gl1.y
		end

		fig1 = Figure(size = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1], title=ttl,
			xticks=collect(year0:4:year1),ylabel=zlb)
		hm1=lines!(ax1,gl1.x,y)
		xlims!(ax1,years_to_display)
		ylims!(ax1,rng)
		fig1
	end

	function DepthTime(XYZ; ClipToRange=true, years_to_display=years_to_display)
		ClipToRange ? to_range!(XYZ.z,XYZ.levels) : nothing
		fig1 = Figure(size = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1], title=XYZ.title,
			xticks=collect(XYZ.year0:4:XYZ.year1))
		hm1=contourf!(ax1,XYZ.x,XYZ.y,XYZ.z,levels=XYZ.levels,colormap=:turbo)
		Colorbar(fig1[1,2], hm1, height = Relative(0.65))
		haskey(XYZ,:years_to_display) ? xlims!(ax1,XYZ.years_to_display) : xlims!(ax1,years_to_display)
		ylims!(ax1,XYZ.ylims)
		fig1
	end

	function TimeLat(XYZ; ClipToRange=true, years_to_display=years_to_display)
		ClipToRange ? to_range!(XYZ.z,XYZ.levels) : nothing
		fig1 = Figure(size = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1], title=XYZ.title,
			xticks=collect(XYZ.year0:4:XYZ.year1),yticks=collect(-90.0:20.0:90.0),ylabel="latitude")
		hm1=contourf!(ax1,XYZ.x,XYZ.y,XYZ.z,levels=XYZ.levels,colormap=:turbo)
		Colorbar(fig1[1,2], hm1, height = Relative(0.65))
		xlims!(ax1,years_to_display)
		ylims!(ax1,XYZ.ylims...)
		fig1
	end

	function map(X; ClipToRange=true)
		ClipToRange ? to_range!(X.field,X.levels) : nothing
		fig = Figure(size = (900,600), backgroundcolor = :grey95)
		ax = Axis(fig[1,1], title=X.title,xlabel="longitude",ylabel="latitude")
		hm1=contourf!(ax,X.λ.lon[:,1],X.λ.lat[1,:],X.field,levels=X.levels,colormap=:turbo)
		Colorbar(fig[1,2], hm1, height = Relative(0.65))
		fig	
	end

	## 2. OISST

	function plot(x::SSTdiag)
		if !isempty(x.options)
			o=x.options
			if string(o.plot_type)=="map_base"
				fig,ax,im=SST_plots.map_base()
				fig
			elseif string(o.plot_type)=="local_and_global"
				SST_plots.local_and_global(o.ts,o.ts_global,o.kdf0)
			elseif string(o.plot_type)=="by_year"
				SST_plots.by_year(o.ts)
			elseif string(o.plot_type)=="by_time"
				SST_plots.by_time(o.ts,show_anom=o.show_anom,show_clim=o.show_clim)
			elseif string(o.plot_type)=="TimeLat"
				SST_plots.TimeLat(o.ts,o.zm,o.title)
			elseif string(o.plot_type)=="MHW"
				SST_plots.MHW(o.ts)
			elseif string(o.plot_type)=="map"
				SST_plots.plot_sst_map(o.to_map)
			else
				error("unknown plot_type")
			end
		else
			error("unknown options")
		end
	end
	
	module SST_plots

	using Makie
	import Climatology: load, Statistics, SSTdiag
	import Climatology: MeshArrays, DataDeps
	import Statistics: median
	
	#
	
	function by_time(ts; show_anom = true, show_clim=true)
		tim=collect(1:length(ts.sst))/365.25 .+ 1982
		f,a=lines(tim,ts.sst,label="SST",linewidth=4)
		show_clim ? lines!(a,tim,ts.clim,color=:orange,label="seasonal climatology",linewidth=1) : nothing
		show_anom ? lines!(a,tim,ts.anom,color=:red,label="SST - seasonal cycle") : nothing
		a.title=ts.title
		xlims!(1982,2024)
		axislegend(a,position=:rb)    
		f
	end
	
	function by_year(ts)
		f,a,l=lines(ts.sst[1:365],color=:gray)
		[lines!(ts.sst[ (1:365) .+ 365*(y-1)] ,color=:gray) for y in 2:length(1982:2022)]
		lines!(ts.sst[ 365*(2023-1982):365*(2024-1982)],color=:orange)
		lines!(ts.sst[ 365*(2024-1982):end],color=:red,linewidth=2)
		for y in 2021:2022
			tt1=vec(1:365) .+(y-1982)*365; lines!(ts.sst[tt1],color=:blue)
		end
		a.title="SST year by year (red=2024, orange=2023, blue=2021:2022)"
		f
	end
	
	#
		
	function to_range!(DD,levs)
		DD[findall(DD.<=levs[1])].=levs[1]+(levs[2]-levs[1])/100
		DD[findall(DD.>=levs[end])].=levs[end]-(levs[end]-levs[end-1])/100
	end
	
	function TimeLat(list,zm,ttl; 
		ClipToRange=true, year0=1982, year1=2024, lat0=-90, lat1=90)
		x=collect(1:length(list.year))/365.25 .+ 1982
		dy=Int(180/size(zm,1))
		y=collect(-90+dy/2:dy:90-dy/2)
		z=permutedims(zm)
		levs=(-2.0:0.25:2.0)/5.0
	
		ClipToRange ? to_range!(z,levs) : nothing
		fig1 = Figure(resolution = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1], title=ttl,
		xticks=collect(year0:4:year1),yticks=collect(-90.0:20.0:90.0),ylabel="latitude")
		hm1=contourf!(ax1,x[1:7:end],y,z[1:7:end,:],levels=levs,colormap=:curl)
		Colorbar(fig1[1,2], hm1, height = Relative(0.65))
		xlims!(ax1,year0,year1)
		ylims!(ax1,lat0,lat1)
		fig1
	end
	
	#
	
	function lowres_scatter(kdf,fig=[],ax=[]; input=[])
		(i,j)=([x.i for x in kdf],[x.j for x in kdf])
		(ii,jj)=(10*i.-5,10*j.-95)
		if isa(fig,Array)
			f,a=scatter(ii,jj,color=input,markersize=10)
			c=(:blue,:red)
		else
			(f,a)=(ax,fig)
			c=(:skyblue,:pink)
		end
		text!(a,ii.+1,jj,text=string.(i),fontsize=11,color=c[1])
		text!(a,ii.+1,jj.-3,text=string.(j),fontsize=11,color=c[2])
		f		
	end
	
	#
	
	function local_and_global(ts,ts_global,kdf0)
		tim=collect(1:length(ts.anom))/365.25 .+ 1982
		fig,ax,li=lines(tim,ts.anom .-median(ts.anom),label="local")
		lines!(tim,ts_global.anom .-median(ts_global.anom),label="global")
		ax.title="local and global SST anomalies"
		xlims!(1982,2024)
		ylims!(-2.5,2.5)
		axislegend(ax,position = :rb)
		fig
	end
	
	function map_base()
		earth_jpg=joinpath(MeshArrays.mydatadep("basemap_jpg1"),
		"Blue_Marble_Next_Generation_%2B_topography_%2B_bathymetry.jpg") 
		
		earth_img=load(earth_jpg)
		earth_img=reverse(permutedims(earth_img),dims=2)
		earth_img=circshift(earth_img,(1800,0))
	
		#fig = Figure(resolution = (1200, 800)) #, backgroundcolor = :grey80)
		fig=with_theme(Figure,theme_light())
		ax = Axis(fig[1, 1])
#		im=image!(ax, -0.05 .. 359.95, -89.95 .. 89.95, 0.5 .+0.5*Gray.(earth_img))
		im=image!(ax, -0.05 .. 359.95, -89.95 .. 89.95, earth_img)
		hidedecorations!(ax)
	
		fig,ax,im
	end
	
	##
	
	function MHW(ts,ttl="SST anomaly with extreme warm periods in red")
		x=ts.sst-ts.clim
		y=fill(:blue,size(x))
		y[findall(x.>=ts.high)].=:red
		tim=collect(1:length(ts.sst))/365.25 .+ 1982
		
		fig,ax,li=lines(tim,x,color=y)
		xlims!(1982,2024)
		ax.title=ttl
		fig
	end
	
	function plot_sst_map(to_map)
		fig=plot(SSTdiag(options=(plot_type=:map_base,)))
		ax=current_axis()
		hm=heatmap!(ax,to_map.lon,to_map.lat,to_map.field,colormap=to_map.colormap,colorrange=to_map.colorrange)
		to_map.showgrid ? lowres_scatter(ax) : nothing
		scatter!(ax,to_map.lon1,to_map.lat1,marker=:circle,color=:blue,markersize=30)
		scatter!(ax,to_map.lon1,to_map.lat1,marker=:x,color=:yellow,markersize=15)
		Colorbar(fig[1, 2],hm)
		ax.title=to_map.title
		fig
	end
	
end
	
## 3. SeaLevelAnomaly

function plot(x::SeaLevelAnomaly)
	SLA_PLOTS.default_plot(x)
end

module SLA_PLOTS

using Makie
import Climatology: SeaLevelAnomaly, SLA_MAIN, Statistics
import Statistics: mean

## Satellite

"""
    default_plot(b::SeaLevelAnomaly; dates=[], kwargs...)
	
```
using Climatology
sla=make_plot(SeaLevelAnomaly(),:sla_podaac)
plot(sla)
```
"""
default_plot(b::SeaLevelAnomaly) = begin
	fig,_,_=prep_movie(b.data[1]; b.options...)
	fig
end

function prep_movie(ds; topo=[], colormap=:PRGn, color=:black, 
	time=1, dates=[], resolution = (600, 400))
	lon=ds["lon"][:]
	lat=ds["lat"][:]
	store=ds["SLA"][:,:,:]

	nt=size(store,3)
	kk=findall((!isnan).(store[:,:,end]))

	n=Observable(time)
	SLA=@lift(store[:,:,$n])
	SLA2=@lift($(SLA).-mean($(SLA)[kk]))

	fig=Figure(size=resolution,fontsize=11)
	ax=Axis(fig[1,1])
    hm=heatmap!(lon,lat,SLA2,colorrange=0.25.*(-1.0,1.0),colormap=colormap)

	if !isempty(topo)
		lon[1]>0.0 ? lon_off=360.0 : lon_off=0.0
		contour!(lon_off.+topo.lon,topo.lat,topo.z,levels=-300:100:300,color=color,linewidth=1)
		contour!(lon_off.+topo.lon,topo.lat,topo.z,levels=-2500:500:-500,color=color,linewidth=0.25)
		contour!(lon_off.+topo.lon,topo.lat,topo.z,levels=-6000:1000:-3000,color=color,linewidth=0.1)
	end

	lon0=minimum(lon)+(maximum(lon)-minimum(lon))/20.0
	lat0=maximum(lat)-(maximum(lat)-minimum(lat))/10.0
	
	if isempty(dates)
		println("no date")
	else
	    dtxt=@lift(string(dates[$n]))
		text!(lon0,lat0,text=dtxt,color=:blue2,fontsize=14,font = :bold)	
	end
	
	Colorbar(fig[1,2],hm)

	fig,n,nt
end

function make_movie(ds,tt; framerate = 90, dates=[])
	fig,n,nt=prep_movie(ds,dates=dates)
    record(fig,tempname()*".mp4", tt; framerate = framerate) do t
        n[] = t
    end
end

end

##

end

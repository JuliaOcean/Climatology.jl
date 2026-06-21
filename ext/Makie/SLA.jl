
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

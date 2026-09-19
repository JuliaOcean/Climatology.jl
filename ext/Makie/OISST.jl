
	function plot(x::SSTdiag)
		isempty(x.options) && error("unknown options")
		o = x.options
		pt = string(o.plot_type)
		if pt=="map_base"
			fig,_,_ = SST_plots.map_base()
			fig
		elseif pt=="local_and_global"
			SST_plots.local_and_global(x)
		elseif pt=="by_year"
			SST_plots.by_year(x)
		elseif pt=="by_time"
			SST_plots.by_time(x)
		elseif pt=="TimeLat"
			SST_plots.TimeLat(x)
		elseif pt=="MHW"
			SST_plots.MHW(x)
		elseif pt=="map"
			SST_plots.plot_sst_map(x)
		else
			error("unknown plot_type")
		end
	end

module SST_plots

using Makie
import Climatology: load, Statistics, SSTdiag
import Climatology: MeshArrays, DataDeps
import Statistics: median

import Climatology: getopt

# unified title resolution: explicit SSTdiag.name wins, else caller-supplied default
title_or(X::SSTdiag, default::AbstractString="") =
    (X.name in ("", "unknown") ? default : X.name)

#

function by_time(X::SSTdiag)
    o = X.options
    ts = o.timeseries
    show_anom = getopt(o,:show_anom,true)
    show_clim = getopt(o,:show_clim,true)
    year0,year1 = getopt(o,:period,(1982,2024))

    tim = collect(1:length(ts.sst))/365.25 .+ 1982
    f,a = lines(tim,ts.sst,label="SST",linewidth=4)
    show_clim ? lines!(a,tim,ts.clim,color=:orange,label="seasonal climatology",linewidth=1) : nothing
    show_anom ? lines!(a,tim,ts.anom,color=:red,label="SST - seasonal cycle") : nothing
    a.title = title_or(X, ts.title)
    xlims!(year0,year1)
    axislegend(a,position=:rb)
    f
end

function by_year(X::SSTdiag)
    ts = X.options.timeseries
    f,a,l = lines(ts.sst[1:365],color=:gray)
    [lines!(ts.sst[ (1:365) .+ 365*(y-1)] ,color=:gray) for y in 2:length(1982:2022)]
    lines!(ts.sst[ 365*(2023-1982):365*(2024-1982)],color=:orange)
    lines!(ts.sst[ 365*(2024-1982):end],color=:red,linewidth=2)
    for y in 2021:2022
        tt1=vec(1:365) .+(y-1982)*365; lines!(ts.sst[tt1],color=:blue)
    end
    a.title = title_or(X,"SST year by year (red=2024, orange=2023, blue=2021:2022)")
    f
end

#

function to_range!(DD,levs)
    DD[findall(DD.<=levs[1])].=levs[1]+(levs[2]-levs[1])/100
    DD[findall(DD.>=levs[end])].=levs[end]-(levs[end]-levs[end-1])/100
end

function TimeLat(X::SSTdiag)
    o = X.options
    list = o.timeseries
    zm   = o.zonal_mean
    year0,year1 = getopt(o,:period,(1982,2024))
    lat0,lat1   = getopt(o,:ylims,(-90,90))
    clip_to_range = getopt(o,:clip_to_range,true)

    x = collect(1:length(list.year))/365.25 .+ 1982
    dy = Int(180/size(zm,1))
    y = collect(-90+dy/2:dy:90-dy/2)
    z = permutedims(zm)
    levs = (-2.0:0.25:2.0)/5.0

    clip_to_range ? to_range!(z,levs) : nothing
    fig1 = Figure(resolution = (900,400),markersize=0.1)
    ax1 = Axis(fig1[1,1], title=title_or(X,"OISST anomaly"),
        xticks=collect(year0:4:year1),yticks=collect(-90.0:20.0:90.0),ylabel="latitude")
    hm1 = contourf!(ax1,x[1:7:end],y,z[1:7:end,:],levels=levs,colormap=:curl)
    Colorbar(fig1[1,2], hm1, height = Relative(0.65))
    xlims!(ax1,year0,year1)
    ylims!(ax1,lat0,lat1)
    fig1
end

#

function lowres_scatter(kdf,fig=[],ax=[]; input=[])
    (i,j) = ([x.i for x in kdf],[x.j for x in kdf])
    (ii,jj) = (10*i.-5,10*j.-95)
    if isa(fig,Array)
        f,a = scatter(ii,jj,color=input,markersize=10)
        c=(:blue,:red)
    else
        (f,a) = (ax,fig)
        c=(:skyblue,:pink)
    end
    text!(a,ii.+1,jj,text=string.(i),fontsize=11,color=c[1])
    text!(a,ii.+1,jj.-3,text=string.(j),fontsize=11,color=c[2])
    f
end

function local_and_global(X::SSTdiag)
    o = X.options
    ts        = o.timeseries
    ts_global = o.timeseries_global
    year0,year1 = getopt(o,:period,(1982,2024))
    ylim0,ylim1 = getopt(o,:ylims,(-2.5,2.5))

    tim = collect(1:length(ts.anom))/365.25 .+ 1982
    fig,ax,li = lines(tim,ts.anom .-median(ts.anom),label="local")
    lines!(tim,ts_global.anom .-median(ts_global.anom),label="global")
    ax.title = title_or(X,"local and global SST anomalies")
    xlims!(year0,year1)
    ylims!(ylim0,ylim1)
    axislegend(ax,position = :rb)
    fig
end

function map_base()
    earth_jpg = joinpath(MeshArrays.mydatadep("basemap_jpg1"),
        "Blue_Marble_Next_Generation_+_topography_+_bathymetry.jpg")
    earth_img = load(earth_jpg)
    earth_img = reverse(permutedims(earth_img),dims=2)
    earth_img = circshift(earth_img,(1800,0))

    fig = with_theme(Figure,theme_light())
    ax = Axis(fig[1, 1])
    im = image!(ax, -0.05 .. 359.95, -89.95 .. 89.95, earth_img)
    hidedecorations!(ax)
    fig,ax,im
end

##

function MHW(X::SSTdiag)
    o = X.options
    ts = o.timeseries
    year0,year1 = getopt(o,:period,(1982,2024))

    x = ts.sst-ts.clim
    y = fill(:blue,size(x))
    y[findall(x.>=ts.high)].=:red
    tim = collect(1:length(ts.sst))/365.25 .+ 1982

    fig,ax,li = lines(tim,x,color=y)
    xlims!(year0,year1)
    ax.title = title_or(X,"SST anomaly with extreme warm periods in red")
    fig
end

function plot_sst_map(X::SSTdiag)
    md = X.options.map_data
    fig,ax,_ = map_base()
    hm = heatmap!(ax,md.lon,md.lat,md.field,colormap=md.colormap,colorrange=md.colorrange)
    md.showgrid ? lowres_scatter(ax) : nothing
    scatter!(ax,md.lon1,md.lat1,marker=:circle,color=:blue,markersize=30)
    scatter!(ax,md.lon1,md.lat1,marker=:x,color=:yellow,markersize=15)
    Colorbar(fig[1, 2],hm)
    ax.title = title_or(X,"SST map")
    fig
end

end

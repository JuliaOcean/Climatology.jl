module plots

using CairoMakie, Statistics, FileIO, Colors, Downloads

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

function save_fig(fig,trigger=true; file="")
    isempty(file) ? fil=tempname()*".png" : fil=joinpath(tempdir(),file)
    save(fil,fig)
    println(fil)
    fig
end

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
	earth_jpg=joinpath(tempdir(),"Blue_Marble.jpg")
	url="https://upload.wikimedia.org/wikipedia/commons/5/56/Blue_Marble_Next_Generation_%2B_topography_%2B_bathymetry.jpg"
	!isfile(earth_jpg) ? Downloads.download(url,earth_jpg) : nothing

	earth_img=load(earth_jpg)
    earth_img=reverse(permutedims(earth_img),dims=2)
	earth_img=circshift(earth_img,(1800,0))

    #fig = Figure(resolution = (1200, 800)) #, backgroundcolor = :grey80)
	fig=with_theme(Figure,theme_light())
    ax = Axis(fig[1, 1])
    im=image!(ax, -0.05 .. 359.95, -89.95 .. 89.95, 0.5 .+0.5*Gray.(earth_img))
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

end

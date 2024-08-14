
using GLMakie, PlutoUI, Printf, JLD2, MeshArrays, Colors, DataDeps

#for source in glob("files*","/Users/gforget/mywork/data/OISST/")
#symlink(source,basename(source))
#end
 
module inc
	using NCDatasets, Statistics, DataFrames, CSV, Glob, Dates
	include("sst_module.jl")
end

(df,gdf,kdf)=inc.coarse_grain.lowres_read(;path="files_csv/");
list=inc.sst_files.read_files_list(path="files_csv/")[1:length(unique(df.t)),:];
G=inc.coarse_grain.grid(list.fil[1]);

tt1=Observable(length(list.fil))
nt=length(list.fil)

function plot!(tt1,fig) 
	fil1=@lift(list.fil[$tt1])
	fil2=@lift(isfile($fil1) ? $fil1 : $fil1[1:end-3]*"_preliminary.nc")

	sst1 = @lift(inc.NCDataset($fil2,"r") do ds
		ds["sst"][:,:,1,1]
	end)

	anom1 = @lift(inc.NCDataset($fil2,"r") do ds
		ds["anom"][:,:,1,1]
	end)

	println(fil2)
	ttl=@lift(@sprintf("Date : %04i / %02i",list.year[$tt1],list.month[$tt1]))
	
	fig
	ax=current_axis(fig)
	hm=heatmap!(ax,G.lon,G.lat,anom1,colormap=:curl,colorrange=6 .*(-1.0,1.0))
	ct=contour!(ax,G.lon,G.lat,sst1,levels=-2:2:40,color=:black,linewidth=1)
	#ct=contour!(ax,G.lon,G.lat,sst1,levels=-2:2:40,color=:white,linestyle = :dash,LineWidth=2)
	#ct=contour!(ax,G.lon,G.lat,sst1,levels=-2:2:40,color=:white,LineWidth=0.5)
	text!(80,-85,text="Data : NOAA OISST")
	text!(10,-85,text=ttl)
	Colorbar(fig[1, 2],hm)
end

lon,lat,earth_img=demo.get_basemap()
fig1=plot_examples(:basemap,lon,lat,0.5 .+0.5*Gray.(earth_img))
plot!(tt1,fig1)

#save(fig1,file="sst_map_dev1.png")

# ╔═╡ 15cff610-fcf7-4365-8c72-358d04df4604
if false
	#1:2:365*4
	fil=joinpath(tempdir(),"sst_v5.mp4")
	record(fig1,fil, 1:nt; framerate = 15) do t 
		tt1[] = t
	end
end

fig1
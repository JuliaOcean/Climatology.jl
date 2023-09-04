### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 91f04e7e-4645-11ec-2d30-ddd4d9932541
begin	
	using Pkg; Pkg.activate()
	using OceanStateEstimation, MeshArrays, ClimateModels, MITgcmTools
	using JLD2, PlutoUI, Glob, TOML, Statistics
	"Done with packages"
end

# ╔═╡ a468baa1-2e5b-40ce-b33c-2e275d720c8e
module plots

	using CairoMakie, JLD2, RollingFunctions, Statistics

	to_range!(DD,levs::Tuple) = to_range!(DD,range(levs[1],levs[2],length=10))

	function to_range!(DD,levs)
		DD[findall(DD.<=levs[1])].=levs[1]+(levs[2]-levs[1])/100
		DD[findall(DD.>=levs[end])].=levs[end]-(levs[end]-levs[end-1])/100
	end

	function axtr1(ax,namtr,pth_out,list_trsp)
		fil_trsp=joinpath(pth_out,"trsp/trsp.jld2")

		itr=findall(list_trsp.==namtr)[1]
		tmp=vec(load(fil_trsp,"single_stored_object"))[itr]
		
		nt=size(tmp.val,2)
		x=vec(0.5:nt)
	
		txt=tmp.nam[1:end-5]
		val=1e-6*vec(sum(tmp.val,dims=1)[:])
		valsmo = runmean(val, 12)
	
		x=vec(0.5:nt)
		x=1992.0 .+ x./12.0

		hm1=lines!(ax,x,val,label="ECCO estimate")
		valsmo[1:5].=NaN
		valsmo[end-4:end].=NaN
		lines!(ax,x,valsmo,linewidth=4.0,color=:red)
		xlims!(ax,(1992.0,2020.0))
	end

	function transport(namtrs,ncols,pth_out,list_trsp)
		if ncols > 1
			fig1 = Figure(resolution = (2000,1000),markersize=0.1)
		else
			fig1 = Figure(resolution = (900,400),markersize=0.1)
		end
		for na in 1:length(namtrs)
			txt=namtrs[na]
			jj=div.(na,ncols,RoundUp)
			kk=na-(jj.-1)*ncols
			ax1 = Axis(fig1[jj,kk], title=" $txt (in Sv)",
				xticks=(1992.0:4:2020.0),ylabel="transport, in Sv")
			axtr1(ax1,namtrs[na],pth_out,list_trsp)
		end
		#ylims!(ax1,rng)
		fig1
	end

	function figov1(pth_out,kk,low1)
		fil=joinpath(pth_out,"overturn/overturn.jld2")
		tmp=-1e-6*load(fil,"single_stored_object")
	
		nt=size(tmp,3)
		x=vec(0.5:nt)
		x=1992.0 .+ x./12.0
		lats=vec(-89.0:89.0)

		fig1 = Figure(resolution = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1],ylabel="Sv",
			title="Global Overturning, in Sv, at kk=$(kk)",
			xticks=(1992.0:4:2020.0))
		for ll in 115:10:145
			ov=tmp[ll,kk,:]
			ov=runmean(ov, 12)
			ov[1:5].=NaN
			ov[end-4:end].=NaN
			hm1=lines!(x,ov,label="$(lats[ll])N")
		end
		xlims!(ax1,(1992.0,2020.0))
		low1!="auto" ? ylims!(ax1,(low1,20.0)) : nothing
		fig1[1, 2] = Legend(fig1, ax1, "estimate", framevisible = false)

	
		fig1
	end

	function figov2(pth_out,Γ; ClipToRange=true)
		fil=joinpath(pth_out,"overturn/overturn.jld2")
		tmp=-1e-6*load(fil,"single_stored_object")
		
		ovmean=dropdims(mean(tmp[:,:,1:240],dims=3),dims=3)
			
		x=vec(-89.0:89.0); y=reverse(vec(Γ.RF[1:end-1])); #coordinate variables
		z=reverse(ovmean,dims=2); z[z.==0.0].=NaN

		levs=(-40.0:5.0:40.0)
		ClipToRange ? to_range!(z,levs) : nothing
	
		fig1 = Figure(resolution = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1], title="Meridional Overturning Streamfunction (in Sv, 92-11)",
				xlabel="latitude",ylabel="depth (in m)")
		hm1=contourf!(ax1,x,y,z,levels=levs,clims=extrema(levs))
		Colorbar(fig1[1,2], hm1, height = Relative(0.65))
		fig1
	end

	function OHT(pth_out)
		fil=joinpath(pth_out,"MHT/MHT.jld2")
		tmp=load(fil,"single_stored_object")
		MT=vec(mean(tmp[:,1:240],dims=2))
	
		x=vec(-89.0:89.0)
		fig1 = Figure(resolution = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1], title="Northward Heat Transport (in PW, 92-11)",
			xticks=(-90.0:10.0:90.0),yticks=(-2.0:0.25:2.0),
			xlabel="latitude",ylabel="Transport (in PW)")
		hm1=lines!(x,MT)
		ylims!(ax1,(-2.0,2.0))
		fig1
	end

	function glo(gl1)
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

		fig1 = Figure(resolution = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1], title=ttl,
			xticks=collect(1992.0:4:2020.0),ylabel=zlb)
		hm1=lines!(ax1,gl1.x,y)
		xlims!(ax1,(1992.0,2020.0))
		ylims!(ax1,rng)
		fig1
	end

	function DepthTime(x,y,z,levs,ttl,RC1,RC0; ClipToRange=true)
		ClipToRange ? to_range!(z,levs) : nothing
		fig1 = Figure(resolution = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1], title=ttl,
			xticks=collect(1992.0:4:2020.0))
		hm1=contourf!(ax1,x,y,z,levels=levs,colormap=:turbo)
		Colorbar(fig1[1,2], hm1, height = Relative(0.65))
		xlims!(ax1,1992.0,2020.0)
		ylims!(ax1,RC1,RC0)
		
		fig1
	end

	function TimeLat(x,y,z,levs,ttl,y0,y1; ClipToRange=true)
		ClipToRange ? to_range!(z,levs) : nothing
		fig1 = Figure(resolution = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1], title=ttl,
			xticks=collect(1992.0:4:2020.0),yticks=collect(-90.0:20.0:90.0),ylabel="latitude")
		hm1=contourf!(ax1,x,y,z,levels=levs,colormap=:turbo)
		Colorbar(fig1[1,2], hm1, height = Relative(0.65))
		xlims!(ax1,1992.0,2020.0)
		ylims!(ax1,y0,y1)
		fig1
	end

	function map(λ,DD,levs,ttl; ClipToRange=true)
		ClipToRange ? to_range!(DD,levs) : nothing
		fig = Figure(resolution = (900,600), backgroundcolor = :grey95)
		ax = Axis(fig[1,1], title=ttl,xlabel="longitude",ylabel="latitude")
		hm1=contourf!(ax,λ.lon[:,1],λ.lat[1,:],DD,levels=levs,colormap=:turbo)
		Colorbar(fig[1,2], hm1, height = Relative(0.65))
		fig	
	end

	function save_fig(fig,trigger)
		fil=tempname()*".png"
		plots.save(fil,fig)
		fil
	end
end

# ╔═╡ 6f721618-d955-4c51-ba44-2873f8609831
PlutoUI.TableOfContents()

# ╔═╡ 63b0b781-c6b0-46a1-af06-a228af8211dc
md"""#  Standard Views of The Ocean State


!!! introduction
	This [Julia](https://julialang.org) [notebook](https://github.com/fonsp/Pluto.jl) let's you explore [ECCO](https://ecco-group.org) ocean state estimates interactively -- [ECCO version 4](https://doi.org/10.5194/gmd-8-3071-2015) [releases 1 to 4](https://ecco-group.org/products.htm) initially. 

!!! note
    - In you are viewing a live version of the notebook, plots will update according to drop down menus as seen in this [video demo](https://youtu.be/UEmBnzspSRg). Directions to run the notebook via [Pluto.jl](https://github.com/fonsp/Pluto.jl), are provided at the bottom of the page. 
    - If instead you are viewing the static html version hosted online, then this interactivity is disabled.
"""

# ╔═╡ 8c4093d7-30aa-4ebe-a429-5d2c2f72fdc3
md"""## Climatology Map"""

# ╔═╡ 1df3bd3c-1396-4cd0-bfd2-3a05dec68261
md"""## Zonal Mean vs Time"""

# ╔═╡ bb3b3089-ab83-4683-9cf0-860a55a9af97
begin
	k_zm_select = @bind k_zm PlutoUI.Slider(1:50, default=1, show_value=true)
	namzm_select = @bind namzm PlutoUI.Select(["MXLDEPTH","THETA","SALT","SSH","SIarea"])
	
	md"""Select a quantity and plot it as a function of time and latitude.
	
	- variable for zonal mean vs time : $(namzm_select)
	- level for zonal mean vs time : $(k_zm_select)
	
	!!! note
	    Choosing a `level` only has an effect if the selected variable, $(namzm), is three-dimensional.

	"""
end

# ╔═╡ 31e97c10-69a6-4074-8b51-89d845620548
md"""## Zonal Mean vs Time (anomalies)"""

# ╔═╡ 5b21c86e-1d75-4510-b474-97ac33fcb271
begin
	namzmanom2d_select = @bind namzmanom2d Select(["MXLDEPTH","SIarea","SSH","THETA","SALT"],default="SALT")
	k_zm2d_select = @bind k_zm2d PlutoUI.Slider(1:50,show_value=true)
	cmap_fac_select = @bind cmap_fac Select(vec([0.05 0.1 0.25 0.5 0.75 1.0 1.5 2.0 5.0]), default=1.0)
	l0_select = @bind l0 PlutoUI.Slider(1:90;default=1, show_value=true)
	l1_select = @bind l1 PlutoUI.Slider(1:90;default=90, show_value=true)

	#cmap_fac_select = @bind cmap_fac Select(string.([0.05 0.1 0.25 0.5 0.75 1.0 1.5 2.0 5.0])[:])
	#cmap_fac_select = @bind cmap_fac Select([1 2])
md"""Select a quantity and plot it's anomaly as a function of time and latitude.

- variable for zonal mean anomaly vs time : $(namzmanom2d_select)
- depth level for zonal mean vs time : $(k_zm2d_select)
- latitude index, min : $(l0_select)
- latitude index, max : $(l1_select)
- scaling factor for color range : $(cmap_fac_select)

!!! note
	Choosing a `level` only has an effect if the selected variable, $(namzmanom2d), is three-dimensional.
"""
end

# ╔═╡ 7dbbb44c-22db-4c30-b71c-58fbab3f78b6
md"""## Depth vs Time (Anomalies)"""

# ╔═╡ 302c84ce-c39d-456b-b748-e3f5ddec0eda
begin
	namzmanom_select = @bind namzmanom Select(["THETA","SALT"])
	l_zm_select = @bind l_Tzm PlutoUI.Slider(8:5:90;default=28,show_value=true)
	k0_select = @bind k0 PlutoUI.Slider(1:50;default=1, show_value=true)
	k1_select = @bind k1 PlutoUI.Slider(1:50;default=30, show_value=true)
	facA_select = @bind facA Select(vec([0.05 0.1 0.25 0.5 0.75 1.0 1.5 2.0 5.0]), default=1.0)

	md"""Settings:
	
	- variable for depth vs time anomaly : $(namzmanom_select)
	- latitude index for depth vs time : $(l_zm_select)
	- top depth level : $(k0_select)
	- bottom depth level : $(k1_select)
	- scaling factor for color range : $(facA_select)
	"""
end

# ╔═╡ 6d9cdea2-272b-4953-8af7-1425817f07a2
md"""## Global Means"""

# ╔═╡ 92d1fc2f-9bdc-41dc-af49-9412f931d882
begin
	ngl1_select = @bind ngl1 Select(["THETA","SALT"];default="THETA")
	kgl1_select = @bind kgl1 PlutoUI.Slider(0:50;default=0, show_value=true)

	md"""Settings
	
	- variable for global mean vs time : $(ngl1_select)
	- depth index, k, for time series : $(kgl1_select)
	  - _(k=0 for volume average; k>0 for area average at level k)_
	"""
end

# ╔═╡ d9c2d8a0-4e5b-4fb5-84cd-c7c989608af5
md"""## Ocean Transports

Here we look at a few aspects of the ocean circulation as estimated in ECCO.

### Meridional Heat Transport
"""

# ╔═╡ c2cd21d9-3fe7-42ec-b6a8-ce34d0770d63
md"""### Overturning Streamfunction"""

# ╔═╡ 215cf4aa-e540-4882-8a33-b7976a6e1b04
md"""### Overturning Time Series"""

# ╔═╡ 7a9269b9-b7aa-4dec-bc86-636a0be6ad01
begin
	ktr1_select = @bind ktr1 PlutoUI.Slider(1:50;default=29, show_value=true)
	low1_select = @bind low1 Select(["auto",-10.0,0.0,5.0,10.0];default="auto")
	
	md"""Seetings:
	
	- level for overturning vs time : $(ktr1_select)
	- lower limit for overturning plot : $(low1_select)
	"""
end

# ╔═╡ ac1cb355-5d59-4d98-9b0a-181a89625b21
md"""### Transport Across One Section"""

# ╔═╡ 0b8ce7b9-8f41-451f-9ec5-5bff418bcafb
md"""### Transport Across Multiple Sections"""

# ╔═╡ 0f308191-13ca-4056-a85f-3a0061958e28
begin
	space = html"<br><br>"
	md"""## Appendices
	$(space)
	"""
end

# ╔═╡ 84e9d2b7-b596-4663-81e0-9bcb172dbdfc
Pkg.status()

# ╔═╡ 6100db56-1681-4f9f-bc53-3792da75e85e
begin
	ECCOdiags_add("release5")
	interpolation_setup()
	"Done with downloading files"
end

# ╔═╡ 64cd25be-2875-4249-b59c-19dcda28a127
begin
	pth=MeshArrays.GRID_LLC90
	γ=GridSpec("LatLonCap",pth)
	Γ=GridLoad(γ;option="full")
	#LC=LatitudeCircles(-89.0:89.0,Γ)

	μ = land_mask(Γ)

	λ_file = joinpath(tempdir(),"interp_coeffs_halfdeg.jld2")
	if !isfile(λ_file)
		lon=[i for i=-179.75:0.5:179.75, j=-89.75:0.5:89.75]
		lat=[j for i=-179.75:0.5:179.75, j=-89.75:0.5:89.75]		
		(f,i,j,w)=InterpolationFactors(Γ,vec(lon),vec(lat))
		jldsave(λ_file; lon=lon, lat=lat, f=f, i=i, j=j, w=w)
	end
	λ = interpolation_setup(λ_file)
	
	"Done with ECCO grid and interpolation"
end

# ╔═╡ a522d3ef-1c94-4eb4-87bc-355965d2ac4a
begin
	sol_list=glob("ECCOv4r*_analysis",ScratchSpaces.ECCO)
	sol_list=[basename(i) for i in sol_list]

	fil_trsp=joinpath(ScratchSpaces.ECCO,"ECCOv4r2_analysis/trsp/trsp.jld2")
	ntr=length(load(fil_trsp,"single_stored_object"))
	list_trsp=[vec(load(fil_trsp,"single_stored_object"))[i].nam for i in 1:ntr] 
	list_trsp=[i[1:end-5] for i in list_trsp]

	pth_colors=joinpath(dirname(pathof(OceanStateEstimation)),"..","examples","ECCO")	
	clim_colors1=TOML.parsefile(joinpath(pth_colors,"clim_colors1.toml"))
	clim_colors2=TOML.parsefile(joinpath(pth_colors,"clim_colors2.toml"))

	function longname(n)
		if occursin("_k",n)
			ln=split(n,"_k")[1]*" at level "*split(n,"_k")[2]
		else
			ln=n
		end
		occursin("BSF",ln) ? ln=replace(ln, "BSF" => "Horizontal Streamfunction (m3/s)") : nothing
		occursin("MXLDEPTH",ln) ? ln=replace(ln, "MXLDEPTH" => "Mixed Layer Depth (m)") : nothing
		occursin("SIarea",ln) ? ln=replace(ln, "SIarea" => "Ice Concentration (0 to 1)") : nothing
		occursin("SSH",ln) ? ln=replace(ln, "SSH" => "Free Surface Height (m)") : nothing
		occursin("THETA",ln) ? ln=replace(ln, "THETA" => "Potential Temperature (degree C)") : nothing
		occursin("SALT",ln) ? ln=replace(ln, "SALT" => "Salinity (psu)") : nothing
		return ln
	end

	function climatology_files(pth_out)
		list_clim=readdir(pth_out)
		kk=findall(occursin.(Ref("clim"),list_clim))
		list_clim=list_clim[kk]
		clim_files=[]
		for ii in 1:length(list_clim)
			tmp=joinpath.(Ref(list_clim[ii]),readdir(joinpath(pth_out,list_clim[ii])))
			[push!(clim_files,i) for i in tmp]
		end
		clim_files
	end
	
	pth_tmp01=joinpath(ScratchSpaces.ECCO,"ECCOv4r2_analysis")
	clim_files=climatology_files(pth_tmp01)
	clim_name=[split(basename(f),'.')[1] for f in clim_files]
	clim_longname=longname.(clim_name) 

	"Done with listing solutions, file names, color codes"
end

# ╔═╡ 17fc2e78-628e-4082-8191-adf07abcc3ff
begin
	nammap_select = @bind nammap Select(clim_longname)
	statmap_select = @bind statmap Select(["mean","std","mon"])	
	timemap_select = @bind timemap Select(1:12)
	md"""
	- file for time mean map : $(nammap_select)
	- choice of statistic for time mean map : $(statmap_select)
	- (optional) if `mon` was selected then show month # : $(timemap_select)
	"""
end

# ╔═╡ 5d320375-0a3c-4197-b35d-f6610173329d
begin
	function glo(pth_out,nam,k=0)
		if k>0
			fil=fil=joinpath(pth_out,nam*"_glo2d/glo2d.jld2")
		else
			fil=joinpath(pth_out,nam*"_glo3d/glo3d.jld2")
		end
		tmp=vec(load(fil,"single_stored_object"))
		occursin("THETA",fil) ? ln=longname("THETA") : ln=longname("SALT")
		if k>0
			nt=Int(length(tmp[:])./50.0)
		    tmp=reshape(tmp,(nt,50))
			tmp=tmp[:,k]
			occursin("THETA",fil) ? rng=(18.0,19.0) : rng=(34.65,34.80)
			txt=ln*" -- level $(k)" 
			k>1 ? rng=extrema(tmp) : nothing
		else
			nt=length(tmp[:])
			occursin("THETA",fil) ? rng=(3.55,3.65) : rng=(34.724,34.728)
			txt=ln
		end

		x=vec(0.5:nt)
		x=1992.0 .+ x./12.0

		(y=tmp,txt=txt,rng=rng,x=x)
	end
end

# ╔═╡ aa340276-cfed-4f0d-a2f1-e6cc18c0bba8
begin
	ntr1_select = @bind ntr1 Select(list_trsp)
	
	md"""Settings:
	
	- transect for transport vs time : $(ntr1_select)	
	"""
end

# ╔═╡ 8b286e86-692f-419c-83c1-f9120e4e35de
begin
	ntr2_select = @bind namtrs MultiCheckBox(list_trsp; orientation=:row, select_all=true, default=[list_trsp[1],list_trsp[2]])
	
	md"""Select Sections:
	
$(ntr2_select)	

!!! note
    The layout of this multiple-panel display should update as you select and unselect sections.
	
	"""
end

# ╔═╡ c46f0656-3627-448b-a779-dad2d980e3cf
md"""### Select Solution
$(space)
Changing solution will update all plots.
"""

# ╔═╡ 8fced956-e527-4ed0-94d4-321368f09773
begin
	sol_select = @bind sol Select(sol_list,default="ECCOv4r2_analysis")
	md"""select a solution : $(sol_select)"""
end

# ╔═╡ 79a9794e-85c6-400e-8b44-3742b56544a2
begin
	pth_out=joinpath(ScratchSpaces.ECCO,sol)
	md"""### Input Data Files
$(space)

Here we read and display results from previous computation that derived transports and other quantities like zonal means from gridded model output. Plots include interpolation from model grid to regular grid.
	
Folder name : $(pth_out)
	"""
end

# ╔═╡ 8563e63d-0096-49f0-8368-e32c4457f5a3
with_terminal() do
	fil_list=readdir(pth_out)
	println.(fil_list)
	"Subfolders And files list:"
end

# ╔═╡ 935cb17d-07b3-4c0c-b863-448ab327d57b
begin
	bind_SaveAllPlots = @bind SaveAllPlots PlutoUI.Button("Save Plots")
	
	md"""### Save All Plots
	$(space)
	All plots will be saved in the folder listed below.
	$(space)
	$(bind_SaveAllPlots)
	"""
end

# ╔═╡ 0a956b36-9306-42e2-a296-3a1840a4cf5b
begin
	MC=ModelConfig(model="ECCO_plots")
	pathof(MC)
end

# ╔═╡ 4d8aa01d-09ef-4f0b-bc7e-16b9ca71a884
begin
	function data_map(nammap)
		ii=findall(clim_longname.==nammap)[1]
		nam=clim_name[ii]
		
		fil=joinpath(pth_out,clim_files[ii])
		if statmap!=="mon"
			tmp=load(fil,statmap)
		else
			tmp=load(fil,statmap)[:,timemap]
		end
	
		DD=Interpolate(μ*tmp,λ.f,λ.i,λ.j,λ.w)
		DD=reshape(DD,size(λ.lon))
		#DD[findall(DD.==0.0)].=NaN
		statmap=="std" ? rng=clim_colors2[nam] : rng=clim_colors1[nam]
		levs=rng[1] .+collect(0.0:0.05:1.0)*(rng[2]-rng[1])
	
		ttl=clim_longname[ii]
		return λ,DD,levs,ttl
	end
	
	MC.outputs[:map]=plots.map(data_map(nammap)...)
end

# ╔═╡ 39ca358a-6e4b-45ed-9ccb-7785884a9868
begin
	pth_out

	function data_TimeLat(namzm)
		fn(x)=transpose(x);
		if namzm=="MXLDEPTH"
			levs=(0.0:50.0:400.0); cm=:turbo
			dlat=2.0; y=vec(-90+dlat/2:dlat:90-dlat/2)
			fil=joinpath(pth_out,namzm*"_zonmean2d/zonmean2d.jld2")
		elseif namzm=="SIarea"
			levs=(0.0:0.1:1.0); cm=:turbo
			dlat=2.0; y=vec(-90+dlat/2:dlat:90-dlat/2)
			fil=joinpath(pth_out,namzm*"_zonmean2d/zonmean2d.jld2")
		elseif namzm=="THETA"
			levs=(-2.0:2.0:34.0); cm=:turbo
			dlat=2.0; y=vec(-90+dlat/2:dlat:90-dlat/2)
			fil=joinpath(pth_out,namzm*"_zonmean/zonmean.jld2")
		elseif namzm=="SALT"
			levs=(32.6:0.2:36.2); cm=:turbo
			dlat=2.0; y=vec(-90+dlat/2:dlat:90-dlat/2)
			fil=joinpath(pth_out,namzm*"_zonmean/zonmean.jld2")
		elseif (namzm=="ETAN")||(namzm=="SSH")
			levs=10*(-0.15:0.02:0.15); cm=:turbo
			dlat=2.0; y=vec(-90+dlat/2:dlat:90-dlat/2)
			fil=joinpath(pth_out,namzm*"_zonmean2d/zonmean2d.jld2")
		else
			levs=missing
		end
	
		tmp=load(fil,"single_stored_object")
		if length(size(tmp))==3
			z=fn(tmp[:,k_zm,:])
			x=vec(0.5:size(tmp,3))
			addon1=" at $(Int(round(Γ.RC[k_zm])))m "
		else
			z=fn(tmp[:,:])
			x=vec(0.5:size(tmp,2))
			addon1=""
		end
	
		x=1992.0 .+ x./12.0
		ttl="$(longname(namzm)) : Zonal Mean $(addon1)"
		return x,y,z,cmap_fac*levs,ttl,-90.0,90.0
	end
	
	MC.outputs[:TimeLat]=plots.TimeLat(data_TimeLat(namzm)...)
end

# ╔═╡ 2d819d3e-f62e-4a73-b51c-0e1204da2369
begin
	pth_out

	function data_TimeLatAnom(namzmanom2d)
		namzm=namzmanom2d
		if namzm=="MXLDEPTH"
			levs=(-100.0:25.0:100.0)/2.0; fn=transpose; cm=:turbo
			fil=joinpath(pth_out,namzm*"_zonmean2d/zonmean2d.jld2")
		elseif namzm=="SIarea"
			levs=(-0.5:0.1:0.5)/5.0; fn=transpose; cm=:turbo
			fil=joinpath(pth_out,namzm*"_zonmean2d/zonmean2d.jld2")
		elseif namzm=="THETA"
			levs=(-2.0:0.25:2.0)/5.0; fn=transpose; cm=:turbo
			fil=joinpath(pth_out,namzm*"_zonmean/zonmean.jld2")
		elseif namzm=="SALT"
			levs=(-0.5:0.1:0.5)/5.0; fn=transpose; cm=:turbo
			fil=joinpath(pth_out,namzm*"_zonmean/zonmean.jld2")
		elseif (namzm=="ETAN")||(namzm=="SSH")
			levs=(-0.5:0.1:0.5)/2.0; fn=transpose; cm=:turbo
			fil=joinpath(pth_out,namzm*"_zonmean2d/zonmean2d.jld2")
		else
			fn=transpose
			levs=missing
		end
	
		tmp=load(fil,"single_stored_object")
		if length(size(tmp))==3
			z=fn(tmp[:,k_zm2d,:])
			x=vec(0.5:size(tmp,3)); 
			addon1=" -- at $(Int(round(Γ.RC[k_zm2d])))m "
		else
			z=fn(tmp[:,:])
			x=vec(0.5:size(tmp,2)); 
			addon1=""
		end
	
		dlat=2.0; y=vec(-90+dlat/2:dlat:90-dlat/2)
		nt=size(z,1)
	
		if true
			#a. subtract monthly mean
			ref1="1992-2011 monthy mean"
			for m in 1:12
				zmean=vec(mean(z[m:12:240,:],dims=1))
				[z[t,:]=z[t,:]-zmean for t in m:12:nt]
			end
		else
			#b. subtract time mean
			ref1="1992-2011 annual mean"
			zmean=vec(mean(z[1:240,:],dims=1))
			[z[t,:]=z[t,:]-zmean for t in 1:nt]
		end
	
		x=1992.0 .+ x./12.0
		ttl="$(longname(namzm)) -- minus $(ref1) $(addon1)"
	
		return x,y,z,cmap_fac*levs,ttl,y[l0],y[l1]
	end

	MC.outputs[:TimeLatAnom]=plots.TimeLat(data_TimeLatAnom(namzmanom2d)...)
end

# ╔═╡ 3f73757b-bab9-4d72-9fff-8884e96e76cd
begin
	pth_out
	
	fn_DepthTime(x)=transpose(x)	

	function data_DepthTime(namzmanom)
	if namzmanom=="THETA"
		levs=(-3.0:0.4:3.0)/8.0; cm=:turbo
		fil=joinpath(pth_out,namzmanom*"_zonmean/zonmean.jld2")
	elseif namzmanom=="SALT"
		levs=(-0.5:0.1:0.5)/10.0;cm=:turbo
		fil=joinpath(pth_out,namzmanom*"_zonmean/zonmean.jld2")
	else
		levs=missing;
	end

	dlat=2.0
	lats=(-90+dlat/2:dlat:90-dlat/2)

	tmp=load(fil,"single_stored_object")
	z=fn_DepthTime(tmp[l_Tzm,:,:])
	addon1=" -- at $(lats[l_Tzm])N "
	x=vec(0.5:size(tmp,3)); 
	y=vec(Γ.RC)
	nt=size(tmp,3)
	
	#a. subtract monthly mean
	ref1="1992-2011 monthy mean"
	for m in 1:12
		zmean=vec(mean(z[m:12:240,:],dims=1))
		[z[t,:]=z[t,:]-zmean for t in m:12:nt]
	end
	#b. subtract time mean
	#ref1="1992-2011 annual mean"
	#zmean=vec(mean(z[1:240,:],dims=1))
	#[z[t,:]=z[t,:]-zmean for t in 1:nt]

	x=1992.0 .+ x./12.0
	ttl="$(longname(namzmanom)) -- minus $(ref1) $(addon1)"

	return x,y,z,facA*levs,ttl,Γ.RC[k1],Γ.RC[k0]
	end

	MC.outputs[:DepthTime]=plots.DepthTime(data_DepthTime(namzmanom)...)
end

# ╔═╡ 16fd6241-8ec1-449d-93ac-ef84c8325867
begin
	save_global=true
	gl1=glo(pth_out,ngl1,kgl1)
	MC.outputs[:global]=plots.glo(gl1)
end

# ╔═╡ a19561bb-f9d6-4f05-9696-9b69bba024fc
MC.outputs[:OHT]=plots.OHT(pth_out)

# ╔═╡ 594c8843-f03f-4230-bdba-a943d535524d
MC.outputs[:overturning]=plots.figov2(pth_out,Γ)

# ╔═╡ 88e85850-b09d-4f46-b104-3489ffe63fa0
MC.outputs[:overturnings]=plots.figov1(pth_out,ktr1,low1)

# ╔═╡ f5e41a76-e56c-4889-821a-68abcb5a72c8
begin
	save_transport=true
	MC.outputs[:transport]=plots.transport([ntr1],1,pth_out,list_trsp)
end

# ╔═╡ 8702a6cf-69de-4e9c-8e77-81f39b55efc7
begin
		#namtrs=[ntr1,ntr1,ntr1,ntr1]
		ncols=Int(floor(sqrt(length(namtrs))))
		MC.outputs[:transports]=plots.transport(namtrs,ncols,pth_out,list_trsp)
end

# ╔═╡ 1fb8f44b-d6f7-4539-8459-fdae07bb6a58
begin
	isempty(MC.outputs) ? listSave=["N/A"] : listSave=collect(keys(MC.outputs))
	bind_SaveOnePlot = @bind SaveOnePlot PlutoUI.Button("Save Plot")
	bind_SelectePlot = @bind SelectPlot PlutoUI.Select(listSave)
	
	md"""### Save Plot
	$(space)
	The selected plot is saved in file listed below.
	$(space)
	$(bind_SaveOnePlot)
	$(bind_SelectePlot)
	"""
end

# ╔═╡ 60b73ddd-ec82-4a52-87b6-058728f150a4
md""" $(bind_SaveOnePlot) for $(bind_SelectePlot) in $(sol_select)"""

# ╔═╡ 0477e49b-d8b2-4308-b692-cadcdfe28892
md""" $(bind_SaveOnePlot) for $(bind_SelectePlot) in $(sol_select)"""

# ╔═╡ 22faa18e-cdf9-411f-8ddb-5b779e44db01
md""" $(bind_SaveOnePlot) for $(bind_SelectePlot) in $(sol_select)"""

# ╔═╡ 3088bca4-0db3-4e4d-a7e5-8afb0f356271
md""" $(bind_SaveOnePlot) for $(bind_SelectePlot) in $(sol_select)"""

# ╔═╡ e88a17f0-5e42-4d0b-8253-e83cabfec4d2
md""" $(bind_SaveOnePlot) for $(bind_SelectePlot) in $(sol_select)"""

# ╔═╡ b55432ac-4960-4983-8330-4ea957a05eee
md""" $(bind_SaveOnePlot) for $(bind_SelectePlot) in $(sol_select)"""

# ╔═╡ 347dc728-2224-4e91-9d7b-45badef8f9a0
md""" $(bind_SaveOnePlot) for $(bind_SelectePlot) in $(sol_select)"""

# ╔═╡ 53069bcc-9b28-40bf-9053-4ec0c6099611
md""" $(bind_SaveOnePlot) for $(bind_SelectePlot) in $(sol_select)"""

# ╔═╡ edf6e079-9aad-4969-b6e3-06dd45b99d68
md""" $(bind_SaveOnePlot) for $(bind_SelectePlot) in $(sol_select)"""

# ╔═╡ 339c792e-7ef1-4554-9f12-d616bc9a7e5b
md""" $(bind_SaveOnePlot) for $(bind_SelectePlot) in $(sol_select)"""

# ╔═╡ c1c4e3c1-d581-4103-8f46-e555c64df86a
let
	SelectPlot!=="N/A" ? fil=plots.save_fig(MC.outputs[SelectPlot],SaveOnePlot) : fil="N/A"
	md"""File = $(fil)"""
end

# ╔═╡ c6ca87f7-fa0d-4cb5-9050-5204f43e0d69
begin
	SaveAllPlots
	MC.outputs
end

# ╔═╡ ff40a006-915a-4d35-847f-5f10085f60a2
begin
	SaveAllPlots
	
	!isdir(pathof(MC)) ? setup(MC) : nothing
	p=joinpath(pathof(MC),"plots")
	!isdir(p) ? mkdir(p) : nothing

	listplots=("overturning","overturnings","transport","transports","OHT",
		"global","DepthTime","TimeLat","TimeLatAnom","map")
	if !isempty(MC.outputs)
		[plots.save(joinpath(p,f*".png"),MC.outputs[Symbol(f)]) for f in listplots]
	end

	display(readdir(p))
end

# ╔═╡ 77339a25-c26c-4bfe-84ee-15274389619f
md""" ### User Directions
$(space)

!!! summary
    Running this notebook on a local computer requires [downloading julia](https://julialang.org/downloads/) (version 1.7 and above), if not already done, and then one can proceed as show below. `Code for steps 2 to 4` is given first in the `grey box`. These commands should be executed in the Julia terminal window (the `REPL`) after installing `julia` in step 1.

```
using Pluto

import OceanStateEstimation
OceanStateEstimation.ECCOdiags_add("release5")
OceanStateEstimation.interpolation_setup()

#optional : ~250M each
OceanStateEstimation.ECCOdiags_add("release4")
OceanStateEstimation.ECCOdiags_add("release3")
OceanStateEstimation.ECCOdiags_add("release2")
OceanStateEstimation.ECCOdiags_add("release1")

Pluto.run()
```

1. [start julia](https://docs.julialang.org/en/v1/manual/getting-started/)
1. download input files (_incl. in code shown above_)
1. [add Pluto](https://github.com/fonsp/Pluto.jl) using [Pkg.jl](https://pkgdocs.julialang.org/v1/getting-started/) (_incl. in code shown above_)
1. [start Pluto](https://github.com/fonsp/Pluto.jl/wiki/🔎-Basic-Commands-in-Pluto) (_incl. in code shown above_)
1. Once Pluto opens, in your web browser, paste the [notebook url](https://raw.githubusercontent.com/gaelforget/OceanStateEstimation.jl/master/examples/ECCO/ECCO_standard_plots.jl) (not the julia code) and click open.

**At first, it may take a couple minutes** for the whole suite of plots to get display, as code compilation takes place (on the fly). Afterwards, since compiled code is available, things should update much faster when using the drop down menus. Another point to note is that notebooks like this one can be executed in various ways outside of Pluto as well.

For more on the underlying software and additional notebooks like this, take a look at the list below.

- [Julia](https://www.julialang.org) ; [docs here](https://docs.julialang.org/en/v1/) and many other resources @ <https://www.julialang.org>
- [Pluto.jl](https://github.com/fonsp/Pluto.jl) ; video from [JuliaCon 2020](https://www.youtube.com/watch?v=IAF8DjrQSSk) that introduced `Pluto.jl`
- [OceanStateEstimation.jl](https://gaelforget.github.io/OceanStateEstimation.jl/dev/) ; where this notebook is located along with the pre-processing loop.
- [MeshArrays.jl](https://juliaclimate.github.io/MeshArrays.jl/dev/) ; array framework able to represents climate model grids like the ones used in ECCO.
- [MITgcmTools.jl](https://gaelforget.github.io/MITgcmTools.jl/dev/) ; functions read/write MITgcm files (inputs & outputs) in various formats.
- [JuliaClimate Notebooks](https://juliaclimate.github.io/GlobalOceanNotebooks/) ; a longer series of notebooks that often use ECCO as an example.
- [OceanObs 2020 workshop](https://github.com/JuliaOcean/JuliaOceanSciencesMeeting2020) ; Julia (language) users and tools for oceanography.
- [JuliaCon 2021 workshop](https://github.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl) ; Modeling Marine Ecosystems At Multiple Scales Using Julia.
"""

# ╔═╡ Cell order:
# ╟─6f721618-d955-4c51-ba44-2873f8609831
# ╟─63b0b781-c6b0-46a1-af06-a228af8211dc
# ╟─8c4093d7-30aa-4ebe-a429-5d2c2f72fdc3
# ╟─4d8aa01d-09ef-4f0b-bc7e-16b9ca71a884
# ╟─17fc2e78-628e-4082-8191-adf07abcc3ff
# ╟─60b73ddd-ec82-4a52-87b6-058728f150a4
# ╟─1df3bd3c-1396-4cd0-bfd2-3a05dec68261
# ╟─39ca358a-6e4b-45ed-9ccb-7785884a9868
# ╟─bb3b3089-ab83-4683-9cf0-860a55a9af97
# ╟─0477e49b-d8b2-4308-b692-cadcdfe28892
# ╟─31e97c10-69a6-4074-8b51-89d845620548
# ╟─2d819d3e-f62e-4a73-b51c-0e1204da2369
# ╟─5b21c86e-1d75-4510-b474-97ac33fcb271
# ╟─22faa18e-cdf9-411f-8ddb-5b779e44db01
# ╟─7dbbb44c-22db-4c30-b71c-58fbab3f78b6
# ╟─3f73757b-bab9-4d72-9fff-8884e96e76cd
# ╟─302c84ce-c39d-456b-b748-e3f5ddec0eda
# ╟─3088bca4-0db3-4e4d-a7e5-8afb0f356271
# ╟─6d9cdea2-272b-4953-8af7-1425817f07a2
# ╟─5d320375-0a3c-4197-b35d-f6610173329d
# ╟─16fd6241-8ec1-449d-93ac-ef84c8325867
# ╟─92d1fc2f-9bdc-41dc-af49-9412f931d882
# ╟─e88a17f0-5e42-4d0b-8253-e83cabfec4d2
# ╟─d9c2d8a0-4e5b-4fb5-84cd-c7c989608af5
# ╟─a19561bb-f9d6-4f05-9696-9b69bba024fc
# ╟─b55432ac-4960-4983-8330-4ea957a05eee
# ╟─c2cd21d9-3fe7-42ec-b6a8-ce34d0770d63
# ╟─594c8843-f03f-4230-bdba-a943d535524d
# ╟─347dc728-2224-4e91-9d7b-45badef8f9a0
# ╟─215cf4aa-e540-4882-8a33-b7976a6e1b04
# ╟─88e85850-b09d-4f46-b104-3489ffe63fa0
# ╟─7a9269b9-b7aa-4dec-bc86-636a0be6ad01
# ╟─53069bcc-9b28-40bf-9053-4ec0c6099611
# ╟─ac1cb355-5d59-4d98-9b0a-181a89625b21
# ╟─f5e41a76-e56c-4889-821a-68abcb5a72c8
# ╟─aa340276-cfed-4f0d-a2f1-e6cc18c0bba8
# ╟─edf6e079-9aad-4969-b6e3-06dd45b99d68
# ╟─0b8ce7b9-8f41-451f-9ec5-5bff418bcafb
# ╟─8702a6cf-69de-4e9c-8e77-81f39b55efc7
# ╟─8b286e86-692f-419c-83c1-f9120e4e35de
# ╟─339c792e-7ef1-4554-9f12-d616bc9a7e5b
# ╟─0f308191-13ca-4056-a85f-3a0061958e28
# ╟─84e9d2b7-b596-4663-81e0-9bcb172dbdfc
# ╟─91f04e7e-4645-11ec-2d30-ddd4d9932541
# ╟─6100db56-1681-4f9f-bc53-3792da75e85e
# ╟─64cd25be-2875-4249-b59c-19dcda28a127
# ╟─a522d3ef-1c94-4eb4-87bc-355965d2ac4a
# ╟─a468baa1-2e5b-40ce-b33c-2e275d720c8e
# ╟─c46f0656-3627-448b-a779-dad2d980e3cf
# ╟─8fced956-e527-4ed0-94d4-321368f09773
# ╟─79a9794e-85c6-400e-8b44-3742b56544a2
# ╟─8563e63d-0096-49f0-8368-e32c4457f5a3
# ╟─1fb8f44b-d6f7-4539-8459-fdae07bb6a58
# ╟─c1c4e3c1-d581-4103-8f46-e555c64df86a
# ╟─935cb17d-07b3-4c0c-b863-448ab327d57b
# ╟─0a956b36-9306-42e2-a296-3a1840a4cf5b
# ╟─c6ca87f7-fa0d-4cb5-9050-5204f43e0d69
# ╟─ff40a006-915a-4d35-847f-5f10085f60a2
# ╟─77339a25-c26c-4bfe-84ee-15274389619f

### A Pluto.jl notebook ###
# v0.19.25

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
	using MeshArrays, OceanStateEstimation, ClimateModels
	using JLD2, PlutoUI, Glob, TOML, Statistics
	ECCOdiags_add("release2")
	ECCOdiags_add("interp_coeffs")
	space = html"<br><br>"
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
md"""#  Ocean State Estimate : Standard Plots


!!! introduction
	This [Julia](https://julialang.org) [notebook](https://github.com/fonsp/Pluto.jl) let's you explore [ECCO](https://ecco-group.org) ocean state estimates interactively -- [ECCO version 4](https://doi.org/10.5194/gmd-8-3071-2015) [releases 1 to 4](https://ecco-group.org/products.htm) initially. 

!!! note
    - In you are viewing a live version of the notebook, plots will update according to drop down menus as seen in this [video demo](https://youtu.be/UEmBnzspSRg). Directions to run the notebook via [Pluto.jl](https://github.com/fonsp/Pluto.jl), are provided at the bottom of the page. 
    - If instead you are viewing the static html version hosted online, then this interactivity is disabled.
"""

# ╔═╡ c46f0656-3627-448b-a779-dad2d980e3cf
md"""## Select Solution

Changing solution will update all plots.
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
md"""## Appendices"""

# ╔═╡ 64cd25be-2875-4249-b59c-19dcda28a127
begin
	pth=MeshArrays.GRID_LLC90
	γ=GridSpec("LatLonCap",pth)
	Γ=GridLoad(γ;option="full")
	#LC=LatitudeCircles(-89.0:89.0,Γ)

	function setup_interp(Γ)
		μ =Γ.hFacC[:,1]
		μ[findall(μ.>0.0)].=1.0
		μ[findall(μ.==0.0)].=NaN
	
		if !isfile(joinpath(ScratchSpaces.ECCO,"interp_coeffs_halfdeg.jld2"))
			lon=[i for i=-179.75:0.5:179.75, j=-89.75:0.5:89.75]
			lat=[j for i=-179.75:0.5:179.75, j=-89.75:0.5:89.75]
			
			(f,i,j,w)=InterpolationFactors(Γ,vec(lon),vec(lat))
			jldsave(joinpath(ScratchSpaces.ECCO,"interp_coeffs_halfdeg.jld2"); 
				lon=lon, lat=lat, f=f, i=i, j=j, w=w, μ=μ)
		end
	
		λ = load(joinpath(ScratchSpaces.ECCO,"interp_coeffs_halfdeg.jld2"))
		λ = MeshArrays.Dict_to_NamedTuple(λ)
	end
	
	λ=setup_interp(Γ)
	
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

# ╔═╡ 8fced956-e527-4ed0-94d4-321368f09773
begin
	sol_select = @bind sol Select(sol_list,default="ECCOv4r2_analysis")
	md"""select a solution : $(sol_select)"""
end

# ╔═╡ 0477e49b-d8b2-4308-b692-cadcdfe28892
md"""select a solution : $(sol_select)"""

# ╔═╡ 22faa18e-cdf9-411f-8ddb-5b779e44db01
md"""Select a solution : $(sol_select)"""

# ╔═╡ 3088bca4-0db3-4e4d-a7e5-8afb0f356271
md"""select a solution : $(sol_select)"""

# ╔═╡ e88a17f0-5e42-4d0b-8253-e83cabfec4d2
md"""select a solution : $(sol_select)"""

# ╔═╡ b55432ac-4960-4983-8330-4ea957a05eee
md"""select a solution : $(sol_select)"""

# ╔═╡ 347dc728-2224-4e91-9d7b-45badef8f9a0
md"""select a solution : $(sol_select)"""

# ╔═╡ 53069bcc-9b28-40bf-9053-4ec0c6099611
md"""select a solution : $(sol_select)"""

# ╔═╡ edf6e079-9aad-4969-b6e3-06dd45b99d68
md"""select a solution : $(sol_select)"""

# ╔═╡ 339c792e-7ef1-4554-9f12-d616bc9a7e5b
md"""select a solution : $(sol_select)"""

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

# ╔═╡ 79a9794e-85c6-400e-8b44-3742b56544a2
begin
	pth_out=joinpath(ScratchSpaces.ECCO,sol)
	md"""### Input Data Files

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
	
	md"""### Save All Plots at Once
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

# ╔═╡ 1fb8f44b-d6f7-4539-8459-fdae07bb6a58
begin
	bind_SaveOnePlot = @bind SaveOnePlot PlutoUI.Button("Save Plot")
	bind_SelectePlot = @bind SelectPlot PlutoUI.Select(collect(keys(MC.outputs)))
	
	md"""### Save Plot
	The selected plot is saved in file listed below.
	$(space)
	$(bind_SaveOnePlot)
	$(bind_SelectePlot)
	"""
end

# ╔═╡ c1c4e3c1-d581-4103-8f46-e555c64df86a
let
	fil=plots.save_fig(MC.outputs[SelectPlot],SaveOnePlot)
	md"""File = $(fil)"""
end

# ╔═╡ 4d8aa01d-09ef-4f0b-bc7e-16b9ca71a884
let
	ii=findall(clim_longname.==nammap)[1]
	nam=clim_name[ii]
	
	fil=joinpath(pth_out,clim_files[ii])
	if statmap!=="mon"
		tmp=load(fil,statmap)
	else
		tmp=load(fil,statmap)[:,timemap]
	end

	DD=Interpolate(λ.μ*tmp,λ.f,λ.i,λ.j,λ.w)
	DD=reshape(DD,size(λ.lon))
	#DD[findall(DD.==0.0)].=NaN
	statmap=="std" ? rng=clim_colors2[nam] : rng=clim_colors1[nam]
	levs=rng[1] .+collect(0.0:0.05:1.0)*(rng[2]-rng[1])

	ttl=clim_longname[ii]

	MC.outputs[:map]=plots.map(λ,DD,levs,ttl)
end

# ╔═╡ 39ca358a-6e4b-45ed-9ccb-7785884a9868
begin
	pth_out

	if namzm=="MXLDEPTH"
		levs=(0.0:50.0:400.0); fn(x)=transpose(x); cm=:turbo
		dlat=2.0; y=vec(-90+dlat/2:dlat:90-dlat/2)
		fil=joinpath(pth_out,namzm*"_zonmean2d/zonmean2d.jld2")
	elseif namzm=="SIarea"
		levs=(0.0:0.1:1.0); fn(x)=transpose(x); cm=:turbo
		dlat=2.0; y=vec(-90+dlat/2:dlat:90-dlat/2)
		fil=joinpath(pth_out,namzm*"_zonmean2d/zonmean2d.jld2")
	elseif namzm=="THETA"
		levs=(-2.0:2.0:34.0); fn(x)=transpose(x); cm=:turbo
		dlat=2.0; y=vec(-90+dlat/2:dlat:90-dlat/2)
		fil=joinpath(pth_out,namzm*"_zonmean/zonmean.jld2")
	elseif namzm=="SALT"
		levs=(32.6:0.2:36.2); fn(x)=transpose(x); cm=:turbo
		dlat=2.0; y=vec(-90+dlat/2:dlat:90-dlat/2)
		fil=joinpath(pth_out,namzm*"_zonmean/zonmean.jld2")
	elseif (namzm=="ETAN")||(namzm=="SSH")
		levs=10*(-0.15:0.02:0.15); fn(x)=transpose(x); cm=:turbo
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
	
	MC.outputs[:TimeLat]=plots.TimeLat(x,y,z,cmap_fac*levs,ttl,-90.0,90.0)
end

# ╔═╡ 2d819d3e-f62e-4a73-b51c-0e1204da2369
let
	pth_out

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

	MC.outputs[:TimeLatAnom]=plots.TimeLat(x,y,z,cmap_fac*levs,ttl,y[l0],y[l1])
end

# ╔═╡ 3f73757b-bab9-4d72-9fff-8884e96e76cd
let
	pth_out
	
	fn(x)=transpose(x)	
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
	z=fn(tmp[l_Tzm,:,:])
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
	
	MC.outputs[:DepthTime]=plots.DepthTime(x,y,z,facA*levs,ttl,Γ.RC[k1],Γ.RC[k0])
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

!!! summary
    Running this notebook on a local computer requires [downloading julia](https://julialang.org/downloads/) (version 1.7 and above), if not already done, and then one can proceed as show below. `Code for steps 2 to 4` is given first in the `grey box`. These commands should be executed in the Julia terminal window (the `REPL`) after installing `julia` in step 1.

```
using Pluto

import OceanStateEstimation
OceanStateEstimation.ECCOdiags_download()
OceanStateEstimation.ECCOdiags_add("interp_coeffs")

#optional : ~250M each
OceanStateEstimation.ECCOdiags_add("release4")
OceanStateEstimation.ECCOdiags_add("release3")
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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
ClimateModels = "f6adb021-9183-4f40-84dc-8cea6f651bb0"
Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
MeshArrays = "cb8c808f-1acf-59a3-9d2b-6e38d009f683"
OceanStateEstimation = "891f6deb-a4f5-4bc5-a2e3-1e8f649cdd2c"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
RollingFunctions = "b0e4dd01-7b14-53d8-9b45-175a3e362653"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
TOML = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[compat]
CairoMakie = "~0.10.4"
ClimateModels = "~0.2.15"
Glob = "~1.3.0"
JLD2 = "~0.4.29"
MeshArrays = "~0.2.32"
OceanStateEstimation = "~0.3.1"
PlutoUI = "~0.7.49"
RollingFunctions = "~0.7.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0-rc2"
manifest_format = "2.0"
project_hash = "7cd5bff773ece4a9f139df426b13573105a7588b"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "16b6dbc4cf7caee4e1e75c49485ec67b667098a0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.AccurateArithmetic]]
deps = ["LinearAlgebra", "Random", "VectorizationBase"]
git-tree-sha1 = "07af26e8d08c211ef85918f3e25d4c0990d20d70"
uuid = "22286c92-06ac-501d-9306-4abd417d9753"
version = "0.3.8"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cc37d689f599e8df4f464b2fa3870ff7db7492ef"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.1"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "38911c7737e123b28182d89027f4216cfc8a9da7"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.3"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e5f08b5689b1aad068e01751889f2f615c7db36d"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.29"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["Printf", "ScanByte", "TranscodingStreams"]
git-tree-sha1 = "d50976f217489ce799e366d9561d56a98a30d7fe"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.2"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "1dd4d9f5beebac0c03446918741b1a03dc5e5788"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.6"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CFTime]]
deps = ["Dates", "Printf"]
git-tree-sha1 = "ed2e76c1c3c43fd9d0cb9248674620b29d71f2d1"
uuid = "179af706-886a-5703-950a-314cd64e0468"
version = "0.1.2"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "Static"]
git-tree-sha1 = "2c144ddb46b552f72d7eafe7cc2f50746e41ea21"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.2"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "SnoopPrecompile", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "c700cce799b51c9045473de751e9319bdd1c6e94"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.9"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["Base64", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "SHA", "SnoopPrecompile"]
git-tree-sha1 = "2aba202861fd2b7603beb80496b6566491229855"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.10.4"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CatViews]]
deps = ["Random", "Test"]
git-tree-sha1 = "23d1f1e10d4e24374112fcf800ac981d14a54b24"
uuid = "81a5f4ea-a946-549a-aa7e-2a7f63a27d31"
version = "1.0.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ClimateModels]]
deps = ["CSV", "DataFrames", "Dates", "Downloads", "Git", "NetCDF", "OrderedCollections", "Pkg", "Statistics", "Suppressor", "TOML", "Test", "UUIDs"]
git-tree-sha1 = "81a64dad5c7e72992fcd4af4d14b1496bab5247b"
uuid = "f6adb021-9183-4f40-84dc-8cea6f651bb0"
version = "0.2.15"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "70232f82ffaab9dc52585e0dd043b5e0c6b714f1"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.12"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "be6ab11021cd29f0344d5c4357b163af05a48cba"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.21.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "89a9db8d28102b094992472d333674bd1a83ce2a"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.1"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    IntervalSetsExt = "IntervalSets"
    StaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "aa51303df86f8626a962fccb878430cdb0a97eee"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.5.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DiskArrays]]
deps = ["OffsetArrays"]
git-tree-sha1 = "6d5b17f84cffa44a1933b50997927fe109039abc"
uuid = "3c3547ce-8d99-4f5e-a174-61eb10b00ae3"
version = "0.3.10"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "49eba9ad9f7ead780bfb7ee319f962c811c6d3b2"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.8"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "13027f188d26206b9e7b863036f87d2f2e7d013a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.87"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "f9818144ce7c8c41edf5c4c179c684d92aa4d9fe"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.6.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "fc86b4fd3eff76c3ce4f5e96e2fdfa6282722885"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.0.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FortranFiles]]
git-tree-sha1 = "f8cec967f151a65f03afd826650c6e91d8b1da16"
uuid = "c58ffaec-ab22-586d-bfc5-781a99fd0b10"
version = "0.6.0"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "cabd77ab6a6fdff49bfd24af2ebe76e6e018a2b4"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.0.0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "38a92e40157100e796690421e34a11c107205c86"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "1cd7f0af1aa58abc02ea1d872953a97359cb87fa"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.4"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "0eb6de0b312688f852f347171aba888658e29f20"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "303202358e38d2b01ba46844b92e48a3c238fd9e"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.6"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Git]]
deps = ["Git_jll"]
git-tree-sha1 = "51764e6c2e84c37055e846c516e9015b4a291c7d"
uuid = "d7ba0133-e1db-5d97-8f8c-041e4b3a1eb2"
version = "1.3.0"

[[deps.Git_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "LibCURL_jll", "Libdl", "Libiconv_jll", "OpenSSL_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "d8be4aab0f4e043cc40984e9097417307cce4c03"
uuid = "f8c6e375-362e-5223-8a59-34ff63f689eb"
version = "2.36.1+2"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Glob]]
git-tree-sha1 = "97285bbd5230dd766e9ef6749b80fc617126d496"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.1"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "678d136003ed5bceaab05cf64519e3f956ffa4ba"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.9.1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "4cc2bb72df6ff40b055295fdef6d92955f9dede8"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.2+2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "734fd90dd2f920a2f1921d5388dcebe805b262dc"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.14"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "432b5b03176f8182bd6841fbfc42c718506a2d5f"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.15"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "c54b581a83008dc7f292e205f4c409ab5caa0f04"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.10"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "342f789fd041a55166764c351da1710db97ce0e0"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.6"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "36cbaebed194b292590cba2593da27b34763804a"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.8"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3d09a9f60edf77f8a4d99f9e015e8fbf9989605d"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.7+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0cb9352ef2e01574eeebdb102948a58740dcaf83"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2023.1.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "16c0cc91853084cb5f58a78bd209513900206ce6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.4"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "Requires", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "42c17b18ced77ff0be65957a591d34f4ed57c631"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.31"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "106b6aa272f294ba47e96bd3acbabdc0407b5c60"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.2"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.KahanSummation]]
git-tree-sha1 = "6292e7878fe190651e74148edb11356dbbc2e194"
uuid = "8e2b3108-d4c1-50be-a7a2-16352aec75c3"
version = "0.3.1"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "4a9513ad756e712177bd342ba6c022b515ed8d76"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.6"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "88b8f66b604da079a627b6fb2860d3704a6729a1"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.14"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDTypes", "SLEEFPirates", "SnoopPrecompile", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "defbfba8ddbccdc8ca3edb4a96a6d6fd3cd33ebd"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.157"

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.LoopVectorization.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MITgcmTools]]
deps = ["ClimateModels", "CodecZlib", "DataFrames", "Dates", "Downloads", "FortranFiles", "MeshArrays", "OrderedCollections", "Printf", "Scratch", "SparseArrays", "Suppressor", "Tar", "UUIDs"]
git-tree-sha1 = "57a9bb10a502f06e2727118682addf075c1405d0"
uuid = "62725fbc-3a66-4df3-9000-e33e85b3a198"
version = "0.2.5"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Makie]]
deps = ["Animations", "Base64", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "InteractiveUtils", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MakieCore", "Markdown", "Match", "MathTeXEngine", "MiniQhull", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "Printf", "Random", "RelocatableFolders", "Setfield", "Showoff", "SignedDistanceFields", "SnoopPrecompile", "SparseArrays", "StableHashTraits", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun"]
git-tree-sha1 = "74657542dc85c3b72b8a5a9392d57713d8b7a999"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.19.4"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "9926529455a331ed73c19ff06d16906737a876ed"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.6.3"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Match]]
git-tree-sha1 = "1d9bc5c1a6e7ee24effb93f175c9342f9154d97f"
uuid = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
version = "1.2.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test", "UnicodeFun"]
git-tree-sha1 = "8f52dbaa1351ce4cb847d95568cb29e62a307d93"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.5.6"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.MeshArrays]]
deps = ["CatViews", "Dates", "LazyArtifacts", "NearestNeighbors", "Pkg", "Printf", "SparseArrays", "Statistics", "Unitful"]
git-tree-sha1 = "93c6743e68d0f3d29cab2d368c1178fe87327599"
uuid = "cb8c808f-1acf-59a3-9d2b-6e38d009f683"
version = "0.2.36"

[[deps.MiniQhull]]
deps = ["QhullMiniWrapper_jll"]
git-tree-sha1 = "9dc837d180ee49eeb7c8b77bb1c860452634b0d1"
uuid = "978d7f02-9e05-4691-894f-ae31a51d76ca"
version = "0.4.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NCDatasets]]
deps = ["CFTime", "DataStructures", "Dates", "NetCDF_jll", "NetworkOptions", "Printf"]
git-tree-sha1 = "fe130b7201b7fd908d950076dbfc0671270894c5"
uuid = "85f8d34a-cbdd-5861-8df4-14fed0d494ab"
version = "0.12.13"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "2c3726ceb3388917602169bed973dbc97f1b51a8"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.13"

[[deps.NetCDF]]
deps = ["DiskArrays", "Formatting", "NetCDF_jll"]
git-tree-sha1 = "328178762645783b20495d408ab317b4c2d25b1a"
uuid = "30363a11-5582-574a-97bb-aa9a979735b9"
version = "0.11.7"

[[deps.NetCDF_jll]]
deps = ["Artifacts", "HDF5_jll", "JLLWrappers", "LibCURL_jll", "Libdl", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "072f8371f74c3b9e1b26679de7fbf059d45ea221"
uuid = "7243133f-43d8-5620-bbf4-c2c921802cf3"
version = "400.902.5+1"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "5ae7ca23e13855b3aba94550f26146c01d259267"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OceanStateEstimation]]
deps = ["CodecZlib", "Distributed", "Downloads", "FortranFiles", "JLD2", "MITgcmTools", "MeshArrays", "NCDatasets", "Pkg", "Printf", "Scratch", "SharedArrays", "Statistics", "TOML", "Tar"]
git-tree-sha1 = "359da75707d37eed1a6dee10de90d3e7e3356125"
uuid = "891f6deb-a4f5-4bc5-a2e3-1e8f649cdd2c"
version = "0.3.2"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "82d7c9e310fe55aa54996e6f7f94674e2a38fcb4"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.9"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "a4ca623df1ae99d09bc9868b008262d0c0ac1e4f"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9ff31d101d987eb9d66bd8b176ac7c277beccd09"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.20+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "67eae2738d63117a196f497d7db789821bce61d1"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.17"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "f809158b27eba0c18c269cf2a2be6ed751d3e81d"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.17"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "ec3edfe723df33528e085e632414499f26650501"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.0"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "84a314e3926ba9ec66ac097e3635e270986b0f10"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.9+0"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "240d7170f5ffdb285f9427b92333c3463bf65bf6"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.1"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "bc2bda41d798c2e66e7c44a11007bb329b15941b"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.0.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "548793c7859e28ef026dba514752275ee871169f"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QhullMiniWrapper_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Qhull_jll"]
git-tree-sha1 = "607cf73c03f8a9f83b36db0b86a3a9c14179621f"
uuid = "460c41e3-6112-5d7f-b78c-b6823adb3f2d"
version = "1.0.0+1"

[[deps.Qhull_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "238dd7e2cc577281976b9681702174850f8d4cbc"
uuid = "784f63db-0788-585a-bace-daefebcd302b"
version = "8.0.1001+0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "6ec7ac8412e83d57e313393220879ede1740f9ee"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "6d7bb727e76147ba18eed998700998e17b8e4911"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.4"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.RollingFunctions]]
deps = ["AccurateArithmetic", "KahanSummation", "LinearAlgebra", "LoopVectorization", "Statistics", "StatsBase"]
git-tree-sha1 = "32e8f64eab12eccd7ef87916c1a050291f122ac4"
uuid = "b0e4dd01-7b14-53d8-9b45-175a3e362653"
version = "0.7.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "8b20084a97b004588125caebf418d8cab9e393d1"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.4"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "cda0aece8080e992f6370491b08ef3909d1c04e7"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.38"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "2436b15f376005e8790e318329560dcc67188e84"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.3"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "77d3c4726515dca71f6d80fbb5e251088defe305"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.18"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableHashTraits]]
deps = ["CRC32c", "Compat", "Dates", "SHA", "Tables", "TupleTools", "UUIDs"]
git-tree-sha1 = "0b8b801b8f03a329a4e86b44c5e8a7d7f4fe10a3"
uuid = "c5dd0088-6c3f-4803-b00e-f31a60c170fa"
version = "0.3.1"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "08be5ee09a7632c32695d954a602df96a877bf0d"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.6"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "Requires", "SnoopPrecompile", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "33040351d2403b84afce74dae2e22d3f5b18edcb"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.4.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "b8d897fe7fa688e93aef573711cb207c08c9e11e"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.19"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "521a0e828e98bb69042fec1809c1b5a680eb7389"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.15"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.Suppressor]]
git-tree-sha1 = "c6ed566db2fe3931292865b966d6d140b7ef32a9"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "c97f60dd4f2331e1a495527f80d242501d2f9865"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.1"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "8621f5c499a8aa4aa970b1ae381aae0ef1576966"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.4"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.TupleTools]]
git-tree-sha1 = "3c712976c47707ff893cf6ba4354aa14db1d8938"
uuid = "9d95972d-f1c8-5527-a6e0-b4b365fa01f6"
version = "1.3.0"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "bb37ed24f338bc59b83e3fc9f32dd388e5396c53"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.12.4"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "b182207d4af54ac64cbc71797765068fdeff475d"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.64"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.4.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

# ╔═╡ Cell order:
# ╟─6f721618-d955-4c51-ba44-2873f8609831
# ╟─63b0b781-c6b0-46a1-af06-a228af8211dc
# ╟─c46f0656-3627-448b-a779-dad2d980e3cf
# ╟─8fced956-e527-4ed0-94d4-321368f09773
# ╟─1fb8f44b-d6f7-4539-8459-fdae07bb6a58
# ╟─c1c4e3c1-d581-4103-8f46-e555c64df86a
# ╟─8c4093d7-30aa-4ebe-a429-5d2c2f72fdc3
# ╟─4d8aa01d-09ef-4f0b-bc7e-16b9ca71a884
# ╟─17fc2e78-628e-4082-8191-adf07abcc3ff
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
# ╟─91f04e7e-4645-11ec-2d30-ddd4d9932541
# ╟─64cd25be-2875-4249-b59c-19dcda28a127
# ╟─a522d3ef-1c94-4eb4-87bc-355965d2ac4a
# ╟─a468baa1-2e5b-40ce-b33c-2e275d720c8e
# ╟─79a9794e-85c6-400e-8b44-3742b56544a2
# ╟─8563e63d-0096-49f0-8368-e32c4457f5a3
# ╟─935cb17d-07b3-4c0c-b863-448ab327d57b
# ╟─0a956b36-9306-42e2-a296-3a1840a4cf5b
# ╟─c6ca87f7-fa0d-4cb5-9050-5204f43e0d69
# ╟─ff40a006-915a-4d35-847f-5f10085f60a2
# ╟─77339a25-c26c-4bfe-84ee-15274389619f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

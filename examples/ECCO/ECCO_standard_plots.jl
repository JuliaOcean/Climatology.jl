### A Pluto.jl notebook ###
# v0.18.0

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
	using MeshArrays, OceanStateEstimation
	using JLD2, PlutoUI, Glob
	using TOML, Statistics, RollingFunctions
	import CairoMakie as Mkie
	"Done with packages"
end

# ╔═╡ 63b0b781-c6b0-46a1-af06-a228af8211dc
md"""#  ECCO Solutions : Standard Plots


!!! introduction
	This notebook let's you explore and compare ocean state estimates from the [ECCO version 4](https://doi.org/10.5194/gmd-8-3071-2015) series ([releases 1 to 4](https://ecco-group.org/products.htm), for now, with _release 5_ coming soon) using [Julia](https://julialang.org). If you are running a live version of the notebook via Pluto.jl (`instructions provided at bottom of this page`) then the various plots will update when you use the drow down menus, as seen in [this video](https://youtu.be/UEmBnzspSRg). If instead you are viewing the html page hosted online then this interactivity is disabled.
"""

# ╔═╡ 6f721618-d955-4c51-ba44-2873f8609831
PlutoUI.TableOfContents()

# ╔═╡ c46f0656-3627-448b-a779-dad2d980e3cf
md"""## Select a Solution"""

# ╔═╡ 8c4093d7-30aa-4ebe-a429-5d2c2f72fdc3
md"""## Climatology Maps"""

# ╔═╡ bb3b3089-ab83-4683-9cf0-860a55a9af97
begin
	k_zm_select = @bind k_zm Slider(1:50, default=1, show_value=true)
	namzm_select = @bind namzm Select(["MXLDEPTH","THETA","SALT","SSH","SIarea"])
	
	md"""## Zonal Means 
	
	Here we select a quantity and plot it as a function of time and latitude.
	
	- variable for zonal mean vs time : $(namzm_select)
	- level for zonal mean vs time : $(k_zm_select)
	
	!!! note
	    Choosing a `level` only has an effect if the selected variable, $(namzm), is three-dimensional.

	"""
end

# ╔═╡ 5b21c86e-1d75-4510-b474-97ac33fcb271
begin
	namzmanom2d_select = @bind namzmanom2d Select(["MXLDEPTH","SIarea","SSH","THETA","SALT"])
	k_zm2d_select = @bind k_zm2d Slider(1:50,show_value=true)
	cmap_fac_select = @bind cmap_fac Select(vec([0.05 0.1 0.25 0.5 0.75 1.0 1.5 2.0 5.0]), default=1.0)
	l0_select = @bind l0 Slider(1:90;default=1, show_value=true)
	l1_select = @bind l1 Slider(1:90;default=90, show_value=true)

	#cmap_fac_select = @bind cmap_fac Select(string.([0.05 0.1 0.25 0.5 0.75 1.0 1.5 2.0 5.0])[:])
	#cmap_fac_select = @bind cmap_fac Select([1 2])
md"""## Zonal Mean Anomalies

Here we select a quantity and plot it's anomaly as a function of time and latitude.

- variable for zonal mean anomaly vs time : $(namzmanom2d_select)
- depth level for zonal mean vs time : $(k_zm2d_select)
- latitude index, min : $(l0_select)
- latitude index, max : $(l1_select)
- scaling factor for color range : $(cmap_fac_select)

!!! note
	Choosing a `level` only has an effect if the selected variable, $(namzmanom2d), is three-dimensional.
"""
end

# ╔═╡ 302c84ce-c39d-456b-b748-e3f5ddec0eda
begin
	namzmanom_select = @bind namzmanom Select(["THETA","SALT"])
	l_zm_select = @bind l_Tzm Slider(8:5:90;default=28,show_value=true)
	k0_select = @bind k0 Slider(1:50;default=1, show_value=true)
	k1_select = @bind k1 Slider(1:50;default=30, show_value=true)
	facA_select = @bind facA Select(vec([0.05 0.1 0.25 0.5 0.75 1.0 1.5 2.0 5.0]), default=1.0)

	md"""### Depth vs Time Anomalies
	
	- variable for depth vs time anomaly : $(namzmanom_select)
	- latitude index for depth vs time : $(l_zm_select)
	- top depth level : $(k0_select)
	- bottom depth level : $(k1_select)
	- scaling factor for color range : $(facA_select)
	"""
end

# ╔═╡ 92d1fc2f-9bdc-41dc-af49-9412f931d882
begin
	ngl1_select = @bind ngl1 Select(["THETA","SALT"];default="THETA")
	kgl1_select = @bind kgl1 Slider(0:50;default=0, show_value=true)

	md"""## Global Means
	
	- variable for global mean vs time : $(ngl1_select)
	- depth index, k, for time series : $(kgl1_select)
	
	!!! note
	    k=0 for volume integral; k>0 for level temperature
	"""
end

# ╔═╡ d9c2d8a0-4e5b-4fb5-84cd-c7c989608af5
md"""## Ocean Transports

Here we look at a few aspects of the ocean circulation as estimated in ECCO.

### Meridional Heat Transport
"""

# ╔═╡ c2cd21d9-3fe7-42ec-b6a8-ce34d0770d63
md"""### Overturning Streamfunction"""

# ╔═╡ 7a9269b9-b7aa-4dec-bc86-636a0be6ad01
begin
	ktr1_select = @bind ktr1 Slider(1:50;default=29, show_value=true)
	low1_select = @bind low1 Select(["auto",-10.0,0.0,5.0,10.0];default="auto")
	
	md"""### Overturning Time Series
	- level for overturning vs time : $(ktr1_select)
	- lower limit for overturning plot : $(low1_select)
	"""
end

# ╔═╡ 0f308191-13ca-4056-a85f-3a0061958e28
md"""## Appendices"""

# ╔═╡ 64cd25be-2875-4249-b59c-19dcda28a127
begin
	pth=MeshArrays.GRID_LLC90
	γ=GridSpec("LatLonCap",pth)
	Γ=GridLoad(γ;option="full")
	#LC=LatitudeCircles(-89.0:89.0,Γ)
	"Done with ECCO grid"
end

# ╔═╡ 963c0bcf-5804-47a5-940e-68f348db95ea
begin
	function setup_interp(Γ)
		μ =Γ.hFacC[:,1]
		μ[findall(μ.>0.0)].=1.0
		μ[findall(μ.==0.0)].=NaN
	
		if !isfile(joinpath(ECCOdiags_path,"interp_coeffs_halfdeg.jld2"))
			lon=[i for i=-179.75:0.5:179.75, j=-89.75:0.5:89.75]
			lat=[j for i=-179.75:0.5:179.75, j=-89.75:0.5:89.75]
			
			(f,i,j,w)=InterpolationFactors(Γ,vec(lon),vec(lat))
			jldsave(joinpath(ECCOdiags_path,"interp_coeffs_halfdeg.jld2"); 
				lon=lon, lat=lat, f=f, i=i, j=j, w=w, μ=μ)
		end
	
		λ = load(joinpath(ECCOdiags_path,"interp_coeffs_halfdeg.jld2"))
		λ = MeshArrays.Dict_to_NamedTuple(λ)
	end
	
	λ=setup_interp(Γ)
	"Done with interpolation coefficients"	
end

# ╔═╡ a522d3ef-1c94-4eb4-87bc-355965d2ac4a
begin
	sol_list=glob("ECCOv4r?_analysis",ECCOdiags_path)
	sol_list=[basename(i) for i in sol_list]

	fil_trsp=joinpath(OceanStateEstimation.ECCOdiags_path,"ECCOv4r2_analysis/trsp/trsp.jld2")
	ntr=length(load(fil_trsp,"single_stored_object"))
	list_trsp=[vec(load(fil_trsp,"single_stored_object"))[i].nam for i in 1:ntr] 
	list_trsp=[i[1:end-5] for i in list_trsp]

	pth_colors=joinpath(dirname(pathof(OceanStateEstimation)),"..","examples","ECCO")	
	clim_colors1=TOML.parsefile(joinpath(pth_colors,"clim_colors1.toml"))
	clim_colors2=TOML.parsefile(joinpath(pth_colors,"clim_colors2.toml"))

	"Done with listing solutions, file names, color codes"
end

# ╔═╡ 8fced956-e527-4ed0-94d4-321368f09773
begin
	sol_select = @bind sol Select(sol_list,default="ECCOv4r2_analysis")
	md"""select a solution : $(sol_select)"""
end

# ╔═╡ b6978bbe-7de6-4827-8b0e-43ceb70e0b52
md"""select a solution : $(sol_select)"""

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

# ╔═╡ aa340276-cfed-4f0d-a2f1-e6cc18c0bba8
begin
	ntr1_select = @bind ntr1 Select(list_trsp)
	
	md"""### Transport Across One Section
	
	- transect for transport vs time : $(ntr1_select)	
	"""
end

# ╔═╡ 8b286e86-692f-419c-83c1-f9120e4e35de
begin
	ntr2_select = @bind namtrs MultiCheckBox(list_trsp; orientation=:row, select_all=true, default=[list_trsp[1],list_trsp[2]])
	
	md"""### Transport Across Multiple Sections

!!! note
    The layout of this multiple-panel display should update as you select and unselect sections.
	
$(ntr2_select)	
	"""
end

# ╔═╡ 758e43bf-e60b-42a3-a3e0-5b928c330892
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

# ╔═╡ 17fc2e78-628e-4082-8191-adf07abcc3ff
begin
	pth_tmp1=joinpath(OceanStateEstimation.ECCOdiags_path,"ECCOv4r2_analysis")
	clim_files=climatology_files(pth_tmp1)	
	nammap_select = @bind nammap Select(clim_files)
	statmap_select = @bind statmap Select(["mean","std","mon"])	
	timemap_select = @bind timemap Select(1:12)
	md"""
	- file for time mean map : $(nammap_select)
	- choice of statistic for time mean map : $(statmap_select)
	- (optional) if `mon` was selected then show month # : $(timemap_select)
	"""
end

# ╔═╡ 79a9794e-85c6-400e-8b44-3742b56544a2
begin
	pth_out=joinpath(OceanStateEstimation.ECCOdiags_path,sol)
	md"""## ECCO Files

Here we read and display results from a previous computation that derived transports and other quantities like zonal means from the gridded model output. Plotting includes interpolation from model grid to regular grid.
	
Folder name : $(pth_out)
	"""
end

# ╔═╡ 4d8aa01d-09ef-4f0b-bc7e-16b9ca71a884
let
	fil=joinpath(pth_out,nammap)
	if statmap!=="mon"
		tmp=load(fil,statmap)
	else
		tmp=load(fil,statmap)[:,timemap]
	end

	DD=Interpolate(λ.μ*tmp,λ.f,λ.i,λ.j,λ.w)
	DD=reshape(DD,size(λ.lon))
	#DD[findall(DD.==0.0)].=NaN
	kk3=basename(nammap)[1:end-5]
	statmap=="std" ? rng=clim_colors2[kk3] : rng=clim_colors1[kk3]
	levs=rng[1] .+collect(0.0:0.05:1.0)*(rng[2]-rng[1])

	fig = Mkie.Figure(resolution = (900,600), backgroundcolor = :grey95)
	ax = Mkie.Axis(fig[1,1], title="file is "*basename(nammap),xlabel="longitude",ylabel="latitude")
	hm1=Mkie.contourf!(ax,λ.lon[:,1],λ.lat[1,:],DD,levels=levs,colormap=:turbo)
	Mkie.Colorbar(fig[1,2], hm1, height = Mkie.Relative(0.65))
	fig	
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

	fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
	if !ismissing(levs)
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
		ax1 = Mkie.Axis(fig1[1,1], title="Zonal Mean $(namzm)$(addon1)",ylabel="latitude",
			xticks=collect(1992.0:4:2021.0),yticks=collect(-90.0:20.0:90.0))
		hm1 = Mkie.contourf!(ax1,x,y,z,levels=levs,clims=extrema(levs),colormap=cm)
		Mkie.xlims!(ax1,1992.0,2021.0)
		Mkie.Colorbar(fig1[1,2], hm1, height = Mkie.Relative(0.65))
	end
	
	fig1
end

# ╔═╡ 2d819d3e-f62e-4a73-b51c-0e1204da2369
let
	pth_out
	fn(x)=transpose(x)

	namzm=namzmanom2d
	if namzm=="MXLDEPTH"
		levs=(-100.0:25.0:100.0)/2.0; fn(x)=transpose(x); cm=:turbo
		fil=joinpath(pth_out,namzm*"_zonmean2d/zonmean2d.jld2")
	elseif namzm=="SIarea"
		levs=(-0.5:0.1:0.5)/5.0; fn(x)=transpose(x); cm=:turbo
		fil=joinpath(pth_out,namzm*"_zonmean2d/zonmean2d.jld2")
	elseif namzm=="THETA"
		levs=(-2.0:0.25:2.0)/5.0; fn(x)=transpose(x); cm=:turbo
		fil=joinpath(pth_out,namzm*"_zonmean/zonmean.jld2")
	elseif namzm=="SALT"
		levs=(-0.5:0.1:0.5)/5.0; fn(x)=transpose(x); cm=:turbo
		fil=joinpath(pth_out,namzm*"_zonmean/zonmean.jld2")
	elseif (namzm=="ETAN")||(namzm=="SSH")
		levs=(-0.5:0.1:0.5)/2.0; fn(x)=transpose(x); cm=:turbo
		fil=joinpath(pth_out,namzm*"_zonmean2d/zonmean2d.jld2")
	else
		levs=missing
	end

	tmp=load(fil,"single_stored_object")
	if length(size(tmp))==3
		z=fn(tmp[:,k_zm2d,:])
		x=vec(0.5:size(tmp,3)); 
		addon1=" at $(Int(round(Γ.RC[k_zm2d])))m "
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
	fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
	ax1 = Mkie.Axis(fig1[1,1], title="Anomaly of $(namzm)$(addon1) ; deviation from $(ref1)",
		xticks=collect(1992.0:4:2021.0),yticks=collect(-90.0:20.0:90.0),ylabel="latitude")
	hm1=Mkie.contourf!(ax1,x,y,z,levels=cmap_fac*levs,colormap=:turbo)
	Mkie.Colorbar(fig1[1,2], hm1, height = Mkie.Relative(0.65))
	Mkie.xlims!(ax1,1992.0,2021.0)
	Mkie.ylims!(ax1,y[l0],y[l1])
	fig1
end

# ╔═╡ 3f73757b-bab9-4d72-9fff-8884e96e76cd
let
	pth_out
	fn(x)=transpose(x)
	if namzmanom=="THETA"
		levs=(-3.0:0.4:3.0)/8.0; fn(x)=transpose(x); cm=:turbo
		fil=joinpath(pth_out,namzmanom*"_zonmean/zonmean.jld2")
	elseif namzmanom=="SALT"
		levs=(-0.5:0.1:0.5)/10.0; fn(x)=transpose(x); cm=:turbo
		fil=joinpath(pth_out,namzmanom*"_zonmean/zonmean.jld2")
	else
		levs=missing
	end

	dlat=2.0
	lats=(-90+dlat/2:dlat:90-dlat/2)

	tmp=load(fil,"single_stored_object")
	z=fn(tmp[l_Tzm,:,:])
	addon1=" at $(lats[l_Tzm])N "
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

	fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
	ax1 = Mkie.Axis(fig1[1,1], title="Anomaly of $(namzmanom)$(addon1) ; deviation from $(ref1)"
		,xticks=collect(1992.0:4:2021.0))
	hm1=Mkie.contourf!(ax1,x,y,z,levels=facA*levs,colormap=:turbo)
	Mkie.Colorbar(fig1[1,2], hm1, height = Mkie.Relative(0.65))
	Mkie.xlims!(ax1,1992.0,2021.0)
	Mkie.ylims!(ax1,Γ.RC[k1],Γ.RC[k0])
	
	fig1
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
		if k>0
			nt=Int(length(tmp[:])./50.0)
		    tmp=reshape(tmp,(nt,50))
			tmp=tmp[:,k]
			occursin("THETA",fil) ? rng=(18.0,19.0) : rng=(34.5,35.0)
			occursin("THETA",fil) ? txt="Temperature (degree C, level $(k))" : txt="Salinity (psu), level $(k)"
			k>1 ? rng=extrema(tmp) : nothing
		else
			nt=length(tmp[:])
			occursin("THETA",fil) ? rng=(3.55,3.65) : rng=(34.72,34.73)
			occursin("THETA",fil) ? txt="Temperature  (degree C)" : txt="Salinity (psu)"
		end

		x=vec(0.5:nt)
		x=1992.0 .+ x./12.0

		(y=tmp,txt=txt,rng=rng,x=x)
	end

	function onegloplot(gl1)
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

		fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
		ax1 = Mkie.Axis(fig1[1,1], title=ttl,
			xticks=collect(1992.0:4:2021.0),ylabel=zlb)
		hm1=Mkie.lines!(ax1,gl1.x,y)
		Mkie.xlims!(ax1,(1992.0,2021.0))
		Mkie.ylims!(ax1,rng)
		fig1
	end

	gl1=glo(pth_out,ngl1,kgl1)
	glfig1=onegloplot(gl1)
end

# ╔═╡ a19561bb-f9d6-4f05-9696-9b69bba024fc
let
	fil=joinpath(pth_out,"MHT/MHT.jld2")
	tmp=load(fil,"single_stored_object")
	MT=vec(mean(tmp[:,1:240],dims=2))

	x=vec(-89.0:89.0)
	fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
	ax1 = Mkie.Axis(fig1[1,1], title="Northward Heat Transport (in PW, 92-11)",
		xticks=(-90.0:10.0:90.0),yticks=(-2.0:0.25:2.0),
		xlabel="latitude",ylabel="Transport (in PW)")
	hm1=Mkie.lines!(x,MT)
	Mkie.ylims!(ax1,(-2.0,2.0))
	fig1
end

# ╔═╡ 12790dfb-5806-498b-8a08-3bfea0dac6a6
let
	fil=joinpath(pth_out,"overturn/overturn.jld2")
	tmp=-1e-6*load(fil,"single_stored_object")
	
	ovmean=dropdims(mean(tmp[:,:,1:240],dims=3),dims=3)
		
	x=vec(-89.0:89.0); y=reverse(vec(Γ.RF[1:end-1])); #coordinate variables
	z=reverse(ovmean,dims=2); z[z.==0.0].=NaN

	fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
	ax1 = Mkie.Axis(fig1[1,1], title="Meridional Overturning Streamfunction (in Sv, 92-11)",
			xlabel="latitude",ylabel="depth (in m)")
	hm1=Mkie.contourf!(ax1,x,y,z,levels=(-40.0:5.0:40.0),clims=(-40,40))
	Mkie.Colorbar(fig1[1,2], hm1, height = Mkie.Relative(0.65))
	fig1

end

# ╔═╡ 88e85850-b09d-4f46-b104-3489ffe63fa0
begin	
	function figov1(pth_out,kk=29)
		fil=joinpath(pth_out,"overturn/overturn.jld2")
		tmp=-1e-6*load(fil,"single_stored_object")
	
		nt=size(tmp,3)
		x=vec(0.5:nt)
		x=1992.0 .+ x./12.0
		lats=vec(-89.0:89.0)

		fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
		ax1 = Mkie.Axis(fig1[1,1],ylabel="Sv",
			title="Global Overturning, in Sv, at kk=$(kk)",
			xticks=(1992.0:4:2021.0))
		for ll in 115:10:145
			ov=tmp[ll,kk,:]
			ov=runmean(ov, 12)
			hm1=Mkie.lines!(x,ov,label="$(lats[ll])N")
		end
		Mkie.xlims!(ax1,(1992.0,2021.0))
		low1!="auto" ? Mkie.ylims!(ax1,(low1,20.0)) : nothing
		fig1[1, 2] = Mkie.Legend(fig1, ax1, "estimate", framevisible = false)

	
		fig1
	end

	figov1(pth_out,ktr1)
end

# ╔═╡ a468baa1-2e5b-40ce-b33c-2e275d720c8e
begin
	function axtr1(ax,namtr)
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

		hm1=Mkie.lines!(ax,x,val,label="ECCO estimate")
		Mkie.lines!(ax,x,valsmo,linewidth=4.0,color=:red)
		Mkie.xlims!(ax,(1992.0,2021.0))
	end

	function transport_plot(namtrs,ncols)
		if ncols > 1
			fig1 = Mkie.Figure(resolution = (2000,1000),markersize=0.1)
		else
			fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
		end
		for na in 1:length(namtrs)
			txt=namtrs[na]
			jj=div.(na,ncols,RoundUp)
			kk=na-(jj.-1)*ncols
			ax1 = Mkie.Axis(fig1[jj,kk], title=" $txt (in Sv)",
				xticks=(1992.0:4:2021.0),ylabel="transport, in Sv")
			axtr1(ax1,namtrs[na])
		end
		#Mkie.ylims!(ax1,rng)
		fig1
	end
end

# ╔═╡ 030dab23-18ed-4e1e-9074-4da8bb9e3ee8
transport_plot([ntr1],1)

# ╔═╡ 8702a6cf-69de-4e9c-8e77-81f39b55efc7
begin
		#namtrs=[ntr1,ntr1,ntr1,ntr1]
		ncols=Int(floor(sqrt(length(namtrs))))
		transport_plot(namtrs,ncols)
end

# ╔═╡ 8563e63d-0096-49f0-8368-e32c4457f5a3
with_terminal() do
	fil_list=readdir(pth_out)
	println.(fil_list)
	"Subfolders And files list:"
end

# ╔═╡ 77339a25-c26c-4bfe-84ee-15274389619f
md""" ## Directions

!!! summary
    Running this notebook on a local computer requires [downloading julia](https://julialang.org/downloads/) (version 1.7 and above), if not already done, and then one can proceed as show below. `Code for steps 2 to 4` is given first in the grey text box. These command should be executed in the Julia terminal window (the `REPL`) after installing `julia` in step 1.

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
Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
MeshArrays = "cb8c808f-1acf-59a3-9d2b-6e38d009f683"
OceanStateEstimation = "891f6deb-a4f5-4bc5-a2e3-1e8f649cdd2c"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
RollingFunctions = "b0e4dd01-7b14-53d8-9b45-175a3e362653"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
TOML = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[compat]
CairoMakie = "~0.7.3"
Glob = "~1.3.0"
JLD2 = "~0.4.21"
MeshArrays = "~0.2.31"
OceanStateEstimation = "~0.2.2"
PlutoUI = "~0.7.34"
RollingFunctions = "~0.6.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.AWS]]
deps = ["Base64", "Compat", "Dates", "Downloads", "GitHub", "HTTP", "IniFile", "JSON", "MbedTLS", "Mocking", "OrderedCollections", "Retry", "Sockets", "URIs", "UUIDs", "XMLDict"]
git-tree-sha1 = "07d944e4d9946c2061f97c1564d1b7ae8ea8f189"
uuid = "fbe9abb3-538b-5e4e-ba9e-bc94f4f92ebc"
version = "1.61.0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "745233d77146ad221629590b6d82fe7f1ddb478f"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "4.0.3"

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

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Blosc]]
deps = ["Blosc_jll"]
git-tree-sha1 = "310b77648d38c223d947ff3f50f511d08690b8d5"
uuid = "a74b3585-a348-5f62-a45c-50e91977d574"
version = "0.7.3"

[[deps.Blosc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Lz4_jll", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "91d6baa911283650df649d0aea7c28639273ae7b"
uuid = "0b7ba130-8d10-5ba8-a3d6-c5182647fed9"
version = "1.21.1+0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[deps.CFTime]]
deps = ["Dates", "Printf"]
git-tree-sha1 = "bca6cb6ee746e6485ca4535f6cc29cf3579a0f20"
uuid = "179af706-886a-5703-950a-314cd64e0468"
version = "0.1.1"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "49f14b6c56a2da47608fe30aed711b5882264d7a"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.9.11"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["Base64", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "SHA", "StaticArrays"]
git-tree-sha1 = "b1d884ee7dae11985314192270eb5762b9ed09ae"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.7.3"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.CatViews]]
deps = ["Random", "Test"]
git-tree-sha1 = "23d1f1e10d4e24374112fcf800ac981d14a54b24"
uuid = "81a5f4ea-a946-549a-aa7e-2a7f63a27d31"
version = "1.0.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f9982ef575e19b0e5c7a98c6e75ee496c0f73a93"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.12.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.ClimateModels]]
deps = ["AWS", "CFTime", "CSV", "DataFrames", "Dates", "Downloads", "Git", "NetCDF", "OrderedCollections", "Pkg", "Statistics", "Suppressor", "TOML", "Test", "UUIDs", "Zarr"]
git-tree-sha1 = "7e9ebc6da08f73d6fc16e7805245ec03b7efc3fd"
uuid = "f6adb021-9183-4f40-84dc-8cea6f651bb0"
version = "0.2.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "3f1f500312161f1ae067abe07d13b40f78f32e07"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.8"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "ae02104e835f219b8930c7664b8012c93475c340"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.2"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiskArrays]]
git-tree-sha1 = "8bd8656d6cec3e13885a822761e8ddfc3b6b2130"
uuid = "3c3547ce-8d99-4f5e-a174-61eb10b00ae3"
version = "0.3.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "38012bf3553d01255e83928eec9c998e19adfddf"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.48"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "d7ab55febfd0907b285fbf8dc0c73c0825d9d6aa"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.3.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ae13fcbc7ab8f16b0856729b050ef0c446aa3492"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.4+0"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.EzXML]]
deps = ["Printf", "XML2_jll"]
git-tree-sha1 = "0fa3b52a04a4e210aeb1626def9c90df3ae65268"
uuid = "8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615"
version = "1.1.0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "463cb335fa22c4ebacfd1faba5fde14edb80d96c"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.5"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "80ced645013a5dbdc52cf70329399c35ce007fae"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.13.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "04d13bfa8ef11720c24e4d840c0033d145537df7"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.17"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "deed294cde3de20ae0b2e0355a6c4e1c6a5ceffc"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.8"

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
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "StaticArrays"]
git-tree-sha1 = "770050893e7bc8a34915b4b9298604a3236de834"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.9.5"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Git]]
deps = ["Git_jll"]
git-tree-sha1 = "d7bffc3fe097e9589145493c08c41297b457e5d0"
uuid = "d7ba0133-e1db-5d97-8f8c-041e4b3a1eb2"
version = "1.2.1"

[[deps.GitHub]]
deps = ["Base64", "Dates", "HTTP", "JSON", "MbedTLS", "Sockets", "SodiumSeal", "URIs"]
git-tree-sha1 = "056781ae7b953289778408b136f8708a46837979"
uuid = "bc5e4493-9b4d-5f90-b8aa-2b2bcaad7a26"
version = "5.7.2"

[[deps.Git_jll]]
deps = ["Artifacts", "Expat_jll", "Gettext_jll", "JLLWrappers", "LibCURL_jll", "Libdl", "Libiconv_jll", "OpenSSL_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "6e93d42b97978709e9c941fa43d0f01701f0d290"
uuid = "f8c6e375-362e-5223-8a59-34ff63f689eb"
version = "2.34.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Glob]]
git-tree-sha1 = "4df9f7e06108728ebf00a0a11edee4b29a482bb2"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "1c5a84319923bea76fa145d49e93aa4394c73fc2"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "70938436e2720e6cb8a7f2ca9f1bbdbf40d7f5d0"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.6.4"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "fd83fa0bde42e01952757f01149dd968c06c4dba"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.0+1"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "9a5c62f231e5bba35695a20988fc7cd6de7eeb5a"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.3"

[[deps.ImageIO]]
deps = ["FileIO", "JpegTurbo", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "464bdef044df52e6436f8c018bea2d48c40bb27b"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.1"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "61feba885fac3a407465726d0c330b3055df897f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b15fc0a95c564ca2e0a7ae12c1f095ca848ceb31"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.5"

[[deps.IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

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
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "28b114b3279cdbac9a61c57b3e6548a572142b34"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.21"

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
git-tree-sha1 = "a77b273f1ddec645d1b7c4fd5fb98c8f90ad10a5"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "591e8dc09ad18386189610acafb970032c519707"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.3"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LRUCache]]
git-tree-sha1 = "d64a0aff6691612ab9fb0117b0995270871c5dfc"
uuid = "8ac3fa9e-de4c-5943-b1dc-09c6b5f20637"
version = "1.3.0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

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
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

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
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Lz4_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "5d494bc6e85c4c9b626ee0cab05daa4085486ab1"
uuid = "5ced341a-0733-55b8-9ab6-a4889d929147"
version = "1.9.3+0"

[[deps.MITgcmTools]]
deps = ["Artifacts", "ClimateModels", "DataFrames", "Dates", "LazyArtifacts", "MeshArrays", "NetCDF", "OrderedCollections", "Printf", "SparseArrays", "Suppressor", "UUIDs"]
git-tree-sha1 = "bc351e1452ebffa346997f809e164815c73c1a67"
uuid = "62725fbc-3a66-4df3-9000-e33e85b3a198"
version = "0.2.1"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "5455aef09b40e5020e1520f551fa3135040d4ed0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+2"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Makie]]
deps = ["Animations", "Base64", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Distributions", "DocStringExtensions", "FFMPEG", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MakieCore", "Markdown", "Match", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "Printf", "Random", "RelocatableFolders", "Serialization", "Showoff", "SignedDistanceFields", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "UnicodeFun"]
git-tree-sha1 = "475b854bff7867c37687d65f7b9498401ac6536d"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.16.4"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "c5fb1bfac781db766f9e4aef96adc19a729bc9b2"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.2.1"

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
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test"]
git-tree-sha1 = "70e733037bbf02d691e78f95171a1fa08cdc6332"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.2.1"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.MeshArrays]]
deps = ["CatViews", "Dates", "Downloads", "LazyArtifacts", "NearestNeighbors", "Pkg", "Printf", "SparseArrays", "Statistics", "Unitful"]
git-tree-sha1 = "c5b9b98540a900934d9b531f8a153ce49016cec7"
uuid = "cb8c808f-1acf-59a3-9d2b-6e38d009f683"
version = "0.2.31"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "29714d0a7a8083bba8427a4fbfb00a540c681ce7"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.3"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "16baacfdc8758bc374882566c9187e785e85c2f0"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.9"

[[deps.NetCDF]]
deps = ["DiskArrays", "Formatting", "NetCDF_jll"]
git-tree-sha1 = "8aca90c1a0f6c3bb65be11d307d8aad6b200087b"
uuid = "30363a11-5582-574a-97bb-aa9a979735b9"
version = "0.11.4"

[[deps.NetCDF_jll]]
deps = ["Artifacts", "HDF5_jll", "JLLWrappers", "LibCURL_jll", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Pkg", "Zlib_jll", "nghttp2_jll"]
git-tree-sha1 = "0cf4d1bf2ef45156aed85c9ac5f8c7e697d9288c"
uuid = "7243133f-43d8-5620-bbf4-c2c921802cf3"
version = "400.702.400+0"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore"]
git-tree-sha1 = "18efc06f6ec36a8b801b23f076e3c6ac7c3bf153"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Observables]]
git-tree-sha1 = "fe29afdef3d0c4a8286128d4e45cc50621b1e43d"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.4.0"

[[deps.OceanStateEstimation]]
deps = ["Artifacts", "CodecZlib", "Downloads", "FortranFiles", "LazyArtifacts", "MITgcmTools", "MeshArrays", "Statistics", "Tar"]
git-tree-sha1 = "cae706d0b95a1e762bdd01b67ff4965559bd918f"
uuid = "891f6deb-a4f5-4bc5-a2e3-1e8f649cdd2c"
version = "0.2.2"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

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
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "ee26b350276c51697c9c2d88a072b339f9f03d73"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.5"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "eb4dbb8139f6125471aa3da98fb70f02dc58e49c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.14"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "1155f6f937fa2b94104162f01fa400e192e4272f"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.4.2"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a121dfbba67c94a5bec9dde613c3d0cbcf3a12b"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.3+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "13468f237353112a01b2d6b32f3d0f80219944aa"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.2"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "6f1b25e8ea06279b5689263cc538f51331d7ca17"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.3"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8979e9802b4ac3d58c503a20f2824ad67f9074dd"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.34"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "db3a23166af8aebf4db5ef87ac5b00d36eb771e2"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "afadeba63d90ff223a6a48d2009434ecee2ec9e8"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.1"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "01d341f502250e81f6fec0afe662aa861392a3aa"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.2"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Retry]]
git-tree-sha1 = "41ac127cd281bb33e42aba46a9d3b25cd35fc6d5"
uuid = "20febd7b-183b-5ae2-ac4a-720e7ce64774"
version = "0.4.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.RollingFunctions]]
deps = ["LinearAlgebra", "Statistics", "StatsBase", "Test"]
git-tree-sha1 = "cdf9158377f81470b1b73c630d0853a3ec0c7445"
uuid = "b0e4dd01-7b14-53d8-9b45-175a3e362653"
version = "0.6.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SIMD]]
git-tree-sha1 = "39e3df417a0dd0c4e1f89891a281f82f5373ea3b"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.0"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "9cc2955f2a254b18be655a4ee70bc4031b2b189e"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "6a2f7d70512d205ca8c7ee31bfa9f142fe74310c"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.12"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

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

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SodiumSeal]]
deps = ["Base64", "Libdl", "libsodium_jll"]
git-tree-sha1 = "80cef67d2953e33935b41c6ab0a178b9987b1c99"
uuid = "2133526b-2bfb-4018-ac12-889fb3908a75"
version = "0.1.1"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "8d0c8e3d0ff211d9ff4a0c2307d876c99d10bdf1"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.2"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "00b725fffc9a7e9aac8850e4ed75b4c1acbe8cd2"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.5.5"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "95c6a5d0e8c69555842fc4a927fc485040ccc31c"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.5"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f35e1879a71cca95f4826a14cdbf0b9e253ed918"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.15"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "d21f2c564b21a202f4677c0fba5b5ee431058544"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.4"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "991d34bbff0d9125d93ba15887d6594e8e84b305"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.5.3"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "b649200e887a487468b71821e2644382699f1b0f"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.11.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "c69f9da3ff2f4f02e811c3323c22e5dfcb584cfa"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.1"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XMLDict]]
deps = ["EzXML", "IterTools", "OrderedCollections"]
git-tree-sha1 = "d9a3faf078210e477b291c79117676fca54da9dd"
uuid = "228000da-037f-5747-90a9-8195ccbf91a5"
version = "0.4.1"

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

[[deps.Zarr]]
deps = ["AWS", "Blosc", "CodecZlib", "DataStructures", "Dates", "DiskArrays", "HTTP", "JSON", "LRUCache", "OffsetArrays", "Pkg", "URIs"]
git-tree-sha1 = "47a53313f4493879345f217e20146c4c40774520"
uuid = "0a941bbe-ad1d-11e8-39d9-ab76183a1d99"
version = "0.7.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.isoband_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "a1ac99674715995a536bbce674b068ec1b7d893d"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.2+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78736dab31ae7a53540a6b752efc61f77b304c5b"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.8.6+1"

[[deps.libsodium_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "848ab3d00fe39d6fbc2a8641048f8f272af1c51e"
uuid = "a9144af2-ca23-56d9-984f-0d03f7b5ccf8"
version = "1.0.20+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

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
# ╟─63b0b781-c6b0-46a1-af06-a228af8211dc
# ╟─6f721618-d955-4c51-ba44-2873f8609831
# ╟─c46f0656-3627-448b-a779-dad2d980e3cf
# ╟─8fced956-e527-4ed0-94d4-321368f09773
# ╟─8c4093d7-30aa-4ebe-a429-5d2c2f72fdc3
# ╟─17fc2e78-628e-4082-8191-adf07abcc3ff
# ╟─4d8aa01d-09ef-4f0b-bc7e-16b9ca71a884
# ╟─b6978bbe-7de6-4827-8b0e-43ceb70e0b52
# ╟─bb3b3089-ab83-4683-9cf0-860a55a9af97
# ╟─39ca358a-6e4b-45ed-9ccb-7785884a9868
# ╟─0477e49b-d8b2-4308-b692-cadcdfe28892
# ╟─5b21c86e-1d75-4510-b474-97ac33fcb271
# ╟─2d819d3e-f62e-4a73-b51c-0e1204da2369
# ╟─22faa18e-cdf9-411f-8ddb-5b779e44db01
# ╟─302c84ce-c39d-456b-b748-e3f5ddec0eda
# ╟─3f73757b-bab9-4d72-9fff-8884e96e76cd
# ╟─3088bca4-0db3-4e4d-a7e5-8afb0f356271
# ╟─92d1fc2f-9bdc-41dc-af49-9412f931d882
# ╟─5d320375-0a3c-4197-b35d-f6610173329d
# ╟─e88a17f0-5e42-4d0b-8253-e83cabfec4d2
# ╟─d9c2d8a0-4e5b-4fb5-84cd-c7c989608af5
# ╟─a19561bb-f9d6-4f05-9696-9b69bba024fc
# ╟─b55432ac-4960-4983-8330-4ea957a05eee
# ╟─c2cd21d9-3fe7-42ec-b6a8-ce34d0770d63
# ╟─12790dfb-5806-498b-8a08-3bfea0dac6a6
# ╟─347dc728-2224-4e91-9d7b-45badef8f9a0
# ╟─7a9269b9-b7aa-4dec-bc86-636a0be6ad01
# ╟─88e85850-b09d-4f46-b104-3489ffe63fa0
# ╟─53069bcc-9b28-40bf-9053-4ec0c6099611
# ╟─aa340276-cfed-4f0d-a2f1-e6cc18c0bba8
# ╟─030dab23-18ed-4e1e-9074-4da8bb9e3ee8
# ╟─edf6e079-9aad-4969-b6e3-06dd45b99d68
# ╟─8b286e86-692f-419c-83c1-f9120e4e35de
# ╟─8702a6cf-69de-4e9c-8e77-81f39b55efc7
# ╟─339c792e-7ef1-4554-9f12-d616bc9a7e5b
# ╟─0f308191-13ca-4056-a85f-3a0061958e28
# ╟─91f04e7e-4645-11ec-2d30-ddd4d9932541
# ╟─64cd25be-2875-4249-b59c-19dcda28a127
# ╟─963c0bcf-5804-47a5-940e-68f348db95ea
# ╟─a522d3ef-1c94-4eb4-87bc-355965d2ac4a
# ╟─758e43bf-e60b-42a3-a3e0-5b928c330892
# ╟─a468baa1-2e5b-40ce-b33c-2e275d720c8e
# ╟─79a9794e-85c6-400e-8b44-3742b56544a2
# ╟─8563e63d-0096-49f0-8368-e32c4457f5a3
# ╟─77339a25-c26c-4bfe-84ee-15274389619f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

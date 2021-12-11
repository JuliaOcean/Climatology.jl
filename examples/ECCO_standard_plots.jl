### A Pluto.jl notebook ###
# v0.17.2

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
	using Pkg
	Pkg.activate()
	
	using JLD2, MeshArrays, OceanStateEstimation, PlutoUI
	using TOML, Statistics, RollingFunctions
	import CairoMakie as Mkie
	"Done with packages"
end

# ╔═╡ 63b0b781-c6b0-46a1-af06-a228af8211dc
md"""# ECCO.v4 Standard Analysis

Explore and compare ocean state estimates from the [ECCO version 4](https://doi.org/10.5194/gmd-8-3071-2015) framework ([release 1 to 5](https://ecco-group.org/products.htm), currently) using [Julia](https://julialang.org). 

- [MeshArrays.jl](https://juliaclimate.github.io/MeshArrays.jl/dev/)
- [MITgcmTools.jl](https://github.com/gaelforget/MITgcmTools.jl)
- [JuliaClimate Notebooks](https://juliaclimate.github.io/GlobalOceanNotebooks/)
- <https://youtu.be/UEmBnzspSRg>

"""

# ╔═╡ 6f721618-d955-4c51-ba44-2873f8609831
PlutoUI.TableOfContents()

# ╔═╡ bb3b3089-ab83-4683-9cf0-860a55a9af97
begin
	k_zm_select = @bind k_zm Slider(1:50, default=1, show_value=true)
	namzm_select = @bind namzm Select(["MXLDEPTH","THETA","SALT","SSH","SIarea"])
	
	md"""## Zonal Means 
	
	Here we select a quantity and plot it vs time and latitude.
	
	- variable for zonal mean vs time : $(namzm_select)
	- level for zonal mean vs time : $(k_zm_select)
	
	_note : choosing level only has an effect if $(namzm) is a three-dimensional variable._

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

Select a variable for zonal mean vs time : $(namzmanom2d_select)

- depth level for zonal mean vs time : $(k_zm2d_select)
- latitude index, min : $(l0_select)
- latitude index, max : $(l1_select)
- scaling factor for color range : $(cmap_fac_select)

!!! note
	Choosing level may only take effect if a three-dimensional variable was selected.
"""
end

# ╔═╡ 302c84ce-c39d-456b-b748-e3f5ddec0eda
begin
	namzmanom_select = @bind namzmanom Select(["THETA","SALT"])
	l_zm_select = @bind l_Tzm Slider(8:5:90;default=28,show_value=true)
	k0_select = @bind k0 Slider(1:50;default=1, show_value=true)
	k1_select = @bind k1 Slider(1:50;default=30, show_value=true)
	facA_select = @bind facA Select(vec([0.05 0.1 0.25 0.5 0.75 1.0 1.5 2.0 5.0]), default=1.0)

	md"""### Depth vs Time Plot
	
	- variable for zonal mean vs time : $(namzmanom_select)
	- latitude index for depth vs time : $(l_zm_select)
	- top depth level : $(k0_select)
	- bottom depth level : $(k1_select)
	- scaling factor for color range : $(facA_select)
	"""
end

# ╔═╡ 92d1fc2f-9bdc-41dc-af49-9412f931d882
begin
	ngl1_select = @bind ngl1 Select(["THETA","SALT"];default="THETA")
	kgl1_select = @bind kgl1 Slider(0:1;default=0, show_value=true)

	md"""## Global Means
	
	- variable for global mean vs time : $(ngl1_select)
	- latitude index for depth vs time : $(kgl1_select)
	"""
end

# ╔═╡ d9c2d8a0-4e5b-4fb5-84cd-c7c989608af5
md"""## Transports"""

# ╔═╡ 7a9269b9-b7aa-4dec-bc86-636a0be6ad01
begin
	ktr1_select = @bind ktr1 Slider(1:50;default=29, show_value=true)
	
	md"""
	- level for overturning vs time : $(ktr1_select)
	"""
end

# ╔═╡ 8fced956-e527-4ed0-94d4-321368f09773
begin
	sol_select = @bind sol Select(["ECCOv4r2_analysis","ECCOv4r3_analysis",
									"ECCOv4r4_analysis","ECCOv4r5_analysis"],default="ECCOv4r2_analysis")
	md"""select a solution : $(sol_select)"""
end

# ╔═╡ c46f0656-3627-448b-a779-dad2d980e3cf
md""" select a solution : $(sol_select)"""

# ╔═╡ 0477e49b-d8b2-4308-b692-cadcdfe28892
md"""select a solution : $(sol_select)"""

# ╔═╡ 22faa18e-cdf9-411f-8ddb-5b779e44db01
md"""Select a solution : $(sol_select)"""

# ╔═╡ e88a17f0-5e42-4d0b-8253-e83cabfec4d2
md"""select a solution : $(sol_select)"""

# ╔═╡ 53069bcc-9b28-40bf-9053-4ec0c6099611
md"""select a solution : $(sol_select)"""

# ╔═╡ 79a9794e-85c6-400e-8b44-3742b56544a2
pth_out=joinpath(OceanStateEstimation.ECCOdiags_path,sol)

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
			occursin("THETA",fil) ? txt="SST  (degree C)" : txt="SSS (psu)"
		else
			nt=length(tmp[:])
			occursin("THETA",fil) ? rng=(3.55,3.65) : rng=(34.72,34.73)
			occursin("THETA",fil) ? txt="Temperature  (degree C)" : txt="Salinity (psu)"
		end
		(y=tmp,txt=txt,rng=rng,x=vec(0.5:nt))
	end

	function onegloplot(gl1)
		fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
		ax1 = Mkie.Axis(fig1[1,1], title="Global Mean $(gl1.txt)",
			xticks=(12:24:336),xlabel="latitude",ylabel="$(gl1.txt)")
		hm1=Mkie.lines!(ax1,gl1.x,gl1.y)
		Mkie.xlims!(ax1,(0.0,336.0))
		Mkie.ylims!(ax1,gl1.rng)
		fig1
	end

	gl1=glo(pth_out,ngl1,kgl1)
	onegloplot(gl1)
end

# ╔═╡ a19561bb-f9d6-4f05-9696-9b69bba024fc
let
	fil=joinpath(pth_out,"MHT/MHT.jld2")
	tmp=load(fil,"single_stored_object")
	MT=vec(mean(tmp[:,1:240],dims=2))

	x=vec(-89.0:89.0)
	fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
	ax1 = Mkie.Axis(fig1[1,1], title="Northward Heat Transport (in PW, 92-11)",
					xticks=(-90.0:10.0:90.0),yticks=(-2.0:0.25:2.0))
	hm1=Mkie.lines!(x,MT,xlabel="latitude",ylabel="Transport (in PW)",label="ECCO estimate")
	Mkie.ylims!(ax1,(-2.0,2.0))
	fig1
end

# ╔═╡ 88e85850-b09d-4f46-b104-3489ffe63fa0
begin	
	function figov1(pth_out,kk=29)
		fil=joinpath(pth_out,"overturn/overturn.jld2")
		tmp=-1e-6*load(fil,"single_stored_object")
	
		nt=size(tmp,3)
		x=vec(0.5:nt)
		lats=vec(-89.0:89.0)

		fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
		ax1 = Mkie.Axis(fig1[1,1],xlabel="month",ylabel="Sv",
		title="Global Overturning, in Sv, at kk=$(kk)",xticks=(12:24:336))
		for ll in 115:10:145
			ov=tmp[ll,kk,:]
			ov=runmean(ov, 12)
			hm1=Mkie.lines!(x,ov,label="$(lats[ll])N")
		end
		Mkie.xlims!(ax1,(0.0,336.0))
		#Mkie.ylims!(ax1,rng)
		fig1[1, 2] = Mkie.Legend(fig1, ax1, "estimate", framevisible = false)

	
		fig1
	end

	figov1(pth_out,ktr1)
end

# ╔═╡ 8563e63d-0096-49f0-8368-e32c4457f5a3
readdir(pth_out)

# ╔═╡ 0f308191-13ca-4056-a85f-3a0061958e28
md"""## Appendices"""

# ╔═╡ 64cd25be-2875-4249-b59c-19dcda28a127
begin
	pth=MeshArrays.GRID_LLC90
	γ=GridSpec("LatLonCap",pth)
	Γ=GridLoad(γ;option="full")
	#LC=LatitudeCircles(-89.0:89.0,Γ)
	"Done with grid"
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
		
		ax1 = Mkie.Axis(fig1[1,1], title="Zonal Mean $(namzm)$(addon1)")
		hm1 = Mkie.contourf!(ax1,x,y,z,levels=levs,clims=extrema(levs),colormap=cm)
		Mkie.xlims!(ax1,0.0,336.0)
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
	
	fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
	ax1 = Mkie.Axis(fig1[1,1], title="Anomaly of $(namzm)$(addon1) ; deviation from $(ref1)")
	hm1=Mkie.contourf!(ax1,x,y,z,levels=cmap_fac*levs,colormap=:turbo)
	Mkie.Colorbar(fig1[1,2], hm1, height = Mkie.Relative(0.65))
	Mkie.xlims!(ax1,0.0,336.0)
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
	
	fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
	ax1 = Mkie.Axis(fig1[1,1], title="Anomaly of $(namzmanom)$(addon1) ; deviation from $(ref1)")
	hm1=Mkie.contourf!(ax1,x,y,z,levels=facA*levs,colormap=:turbo)
	Mkie.Colorbar(fig1[1,2], hm1, height = Mkie.Relative(0.65))
	Mkie.xlims!(ax1,0.0,336.0)
	Mkie.ylims!(ax1,Γ.RC[k1],Γ.RC[k0])
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
	ax1 = Mkie.Axis(fig1[1,1], title="Meridional Overturning Streamfunction (in Sv, 92-11)")
	hm1=Mkie.contourf!(ax1,x,y,z,levels=(-40.0:5.0:40.0),clims=(-40,40))
	Mkie.Colorbar(fig1[1,2], hm1, height = Mkie.Relative(0.65))
	fig1

end

# ╔═╡ 963c0bcf-5804-47a5-940e-68f348db95ea
begin
	function setup_interp(Γ)
		μ =Γ.hFacC[:,1]
		μ[findall(μ.>0.0)].=1.0
		μ[findall(μ.==0.0)].=NaN
	
		if !isfile(joinpath(tempdir(),"interp_coeffs_halfdeg.jld2"))
			lon=[i for i=-179.75:0.5:179.75, j=-89.75:0.5:89.75]
			lat=[j for i=-179.75:0.5:179.75, j=-89.75:0.5:89.75]
			
			(f,i,j,w)=InterpolationFactors(Γ,vec(lon),vec(lat))
			jldsave(joinpath(tempdir(),"interp_coeffs_halfdeg.jld2"); 
				lon=lon, lat=lat, f=f, i=i, j=j, w=w, μ=μ)
		end
	
		λ = load(joinpath(tempdir(),"interp_coeffs_halfdeg.jld2"))
		λ = MeshArrays.Dict_to_NamedTuple(λ)
	end
	
	λ=setup_interp(Γ)
	"Done with interpolation coefficients"	
end

# ╔═╡ a522d3ef-1c94-4eb4-87bc-355965d2ac4a
begin
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
	
	clim_colors1=TOML.parsefile("clim_colors1.toml")
	clim_colors2=TOML.parsefile("clim_colors2.toml")

	fil_trsp=joinpath(OceanStateEstimation.ECCOdiags_path,"ECCOv4r2_analysis/trsp/trsp.jld2")
	ntr=length(load(fil_trsp,"single_stored_object"))
	list_trsp=[vec(load(fil_trsp,"single_stored_object"))[i].nam for i in 1:ntr] 

	"Done with listing files"
end

# ╔═╡ 17fc2e78-628e-4082-8191-adf07abcc3ff
begin
	pth_tmp1=joinpath("ECCO_diags","ECCOv4r2_analysis")
	clim_files=climatology_files(pth_tmp1)	
	nammap_select = @bind nammap Select(clim_files)
	statmap_select = @bind statmap Select(["mean","std","mon"])
	timemap_select = @bind timemap Select(1:12)
	md"""## Climatology Maps

	- file for time mean map : $(nammap_select)
	- statistics for time mean map : $(statmap_select)
	- (optional) month for time mean map : $(timemap_select)
	
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
	ax = Mkie.Axis(fig[1,1], title=nammap,xlabel="longitude",ylabel="latitude")
	hm1=Mkie.contourf!(ax,λ.lon[:,1],λ.lat[1,:],DD,levels=levs,colormap=:turbo)
	Mkie.Colorbar(fig[1,2], hm1, height = Mkie.Relative(0.65))
	fig	
end

# ╔═╡ aa340276-cfed-4f0d-a2f1-e6cc18c0bba8
begin
	ntr1_select = @bind ntr1 Select(list_trsp)
	
	md"""### Transport Across One Section
	
	- transect for transport vs time : $(ntr1_select)	
	"""
end

# ╔═╡ 57d01a67-01c7-4d61-93c7-737ef2cbb6a9
begin
	function figtr1(namtr)
		itr=findall(list_trsp.==namtr)[1]
		fil_trsp=joinpath(pth_out,"trsp/trsp.jld2")
		tmp=vec(load(fil_trsp,"single_stored_object"))[itr]
		
		nt=size(tmp.val,2)
		x=vec(0.5:nt)
	
		txt=tmp.nam[1:end-5]
		val=1e-6*vec(sum(tmp.val,dims=1)[:])
		valsmo = runmean(val, 12)
	
		x=vec(0.5:nt)
		fig1 = Mkie.Figure(resolution = (900,400),markersize=0.1)
		ax1 = Mkie.Axis(fig1[1,1], title=" $txt (in Sv)",
			xticks=(12:24:336),xlabel="latitude",ylabel="transport, in Sv")
		hm1=Mkie.lines!(x,val,label="ECCO estimate")
		Mkie.lines!(x,valsmo,linewidth=4.0,color=:red)
		Mkie.xlims!(ax1,(0.0,336.0))
		#Mkie.ylims!(ax1,rng)
		fig1
	end
	figtr1(ntr1)
end

# ╔═╡ 8b286e86-692f-419c-83c1-f9120e4e35de
begin
	ntr2_select = @bind namtrs MultiCheckBox(list_trsp; orientation=:row, select_all=true, default=[list_trsp[1],list_trsp[2]])
	
	md"""### Transport Across Several Sections
	
	$(ntr2_select)	
	"""
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

		hm1=Mkie.lines!(ax,x,val,label="ECCO estimate")
		Mkie.lines!(ax,x,valsmo,linewidth=4.0,color=:red)
		Mkie.xlims!(ax,(0.0,336.0))
	end

	function figtr2(namtrs,ncols)
		fig1 = Mkie.Figure(resolution = (2000,1000),markersize=0.1)
		for na in 1:length(namtrs)
			txt=namtrs[na][1:end-5]
			jj=div.(na,ncols,RoundUp)
			kk=na-(jj.-1)*ncols
			ax1 = Mkie.Axis(fig1[jj,kk], title=" $txt (in Sv)",
				xticks=(12:24:336),xlabel="latitude",ylabel="transport, in Sv")
			axtr1(ax1,namtrs[na])
		end
		#Mkie.ylims!(ax1,rng)
		fig1
	end
end

# ╔═╡ 8702a6cf-69de-4e9c-8e77-81f39b55efc7
begin
		#namtrs=[ntr1,ntr1,ntr1,ntr1]
		ncols=Int(floor(sqrt(length(namtrs))))
		ff=figtr2(namtrs,ncols)
end

# ╔═╡ Cell order:
# ╟─63b0b781-c6b0-46a1-af06-a228af8211dc
# ╟─6f721618-d955-4c51-ba44-2873f8609831
# ╟─c46f0656-3627-448b-a779-dad2d980e3cf
# ╟─17fc2e78-628e-4082-8191-adf07abcc3ff
# ╟─4d8aa01d-09ef-4f0b-bc7e-16b9ca71a884
# ╟─bb3b3089-ab83-4683-9cf0-860a55a9af97
# ╟─39ca358a-6e4b-45ed-9ccb-7785884a9868
# ╟─0477e49b-d8b2-4308-b692-cadcdfe28892
# ╟─5b21c86e-1d75-4510-b474-97ac33fcb271
# ╟─2d819d3e-f62e-4a73-b51c-0e1204da2369
# ╟─22faa18e-cdf9-411f-8ddb-5b779e44db01
# ╟─302c84ce-c39d-456b-b748-e3f5ddec0eda
# ╟─3f73757b-bab9-4d72-9fff-8884e96e76cd
# ╟─92d1fc2f-9bdc-41dc-af49-9412f931d882
# ╟─5d320375-0a3c-4197-b35d-f6610173329d
# ╟─e88a17f0-5e42-4d0b-8253-e83cabfec4d2
# ╟─d9c2d8a0-4e5b-4fb5-84cd-c7c989608af5
# ╟─12790dfb-5806-498b-8a08-3bfea0dac6a6
# ╟─a19561bb-f9d6-4f05-9696-9b69bba024fc
# ╟─7a9269b9-b7aa-4dec-bc86-636a0be6ad01
# ╟─88e85850-b09d-4f46-b104-3489ffe63fa0
# ╟─53069bcc-9b28-40bf-9053-4ec0c6099611
# ╟─aa340276-cfed-4f0d-a2f1-e6cc18c0bba8
# ╟─57d01a67-01c7-4d61-93c7-737ef2cbb6a9
# ╟─8b286e86-692f-419c-83c1-f9120e4e35de
# ╟─8702a6cf-69de-4e9c-8e77-81f39b55efc7
# ╟─8fced956-e527-4ed0-94d4-321368f09773
# ╟─79a9794e-85c6-400e-8b44-3742b56544a2
# ╟─8563e63d-0096-49f0-8368-e32c4457f5a3
# ╟─0f308191-13ca-4056-a85f-3a0061958e28
# ╟─64cd25be-2875-4249-b59c-19dcda28a127
# ╟─91f04e7e-4645-11ec-2d30-ddd4d9932541
# ╟─963c0bcf-5804-47a5-940e-68f348db95ea
# ╟─a522d3ef-1c94-4eb4-87bc-355965d2ac4a
# ╟─a468baa1-2e5b-40ce-b33c-2e275d720c8e

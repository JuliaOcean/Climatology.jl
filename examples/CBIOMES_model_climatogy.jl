### A Pluto.jl notebook ###
# v0.16.4

using Markdown
using InteractiveUtils

# ╔═╡ 4e55fa2e-3193-11ec-3a50-2defc49cc54d
begin
	using Pkg
	Pkg.activate()
	using NCTiles, MITgcmTools, MeshArrays, OceanColorData, OceanStateEstimation
	using PlutoUI, JLD2, CairoMakie
	import CairoMakie as Mkie
	"Done with packages"
end

# ╔═╡ 979a1b66-0aaa-4b8f-868c-85db76d3fb99
md"""# Interpolated Climatology Workflow

In this notebook we illustrate how `MeshArrays.jl` and `NCTiles.jl` can be used to:

- start from native-grid output from the [CBIOMES-global](https://cbiomes.readthedocs.io/en/latest/)
- read from file, derive, and interapolate new variables
- write results to CF-compliant NetCDF files
- plot result maps interactively (TBD)

`Author: Gaël Forget`

`Date: 2021/10/25`
"""

# ╔═╡ aa41a861-7eca-4dda-92e8-bd8e96f2a21e
PlutoUI.TableOfContents()

# ╔═╡ 5fcb77a8-b234-45d9-a4bd-02bd4eb8fbd9
md"""## Typical Workflow

1. encode functions like `▶a(x,t)`, `▶b(x,t)` to derive target variables
   - input: `x` is path to input variable folder (`nctiles`); `t"` is time record
   - output: one two-dimensional field (the derivation end-result)
1. loop over target variables
   - output result to temporary storage (`jld2` files)
   - create and instantiate new netcdf files (`nc` files)
1. combine all `nc` files into one

!!! note
	- the first time, this notebook generates `half_degree_grid_coeffs.jld2` which may add $O(1')$ runtime
    - derivation methods ▶c, to ▶g illustrate different input/output setups (e.g. N inputs to M outputs)
    - to create and instantiate `nc` files, we define `NCTiles.NCT` instances and simply call `NCTiles.write()`
    - could use infinite time dimension and increment file instead to avoid reading 12 fields at once
"""

# ╔═╡ ad9d39a2-65b3-4fd9-aecd-f360ac34df0f
md"""## Configuration"""

# ╔═╡ d7d0fc22-0bd7-4c4b-9249-8084d75bb33e
md"""## Derivation Methods

The functions provided below derive the target variables listed hereafter.

| Variables [1/3] | Variables [2/3]  | Variables [3/3]  |
|:----------------|:----------------:|-----------------:|
| lon             | lat              |                  |
| SST             | SSS              | Chl              |
| MLD             | EuphoticDepth    | OceanDepth       |
| WindSpeed       | PAR              | TKE              |
| Rrs412          | Rrs443           | Rrs490           |
| Rrs510          | Rrs555           | Rrs670           |

_(note: could generate via `PrettyTables.jl` and `with__terminal` instead)_

### Interpolation, Subsampling
"""

# ╔═╡ 7a1058e4-3813-4b90-a9ee-ad5ff4eb6d6f
begin
	pth0="201805-CBIOMES-climatology/nctiles/"
	γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
	Γ=GridLoad(γ; option="full");
	"Done with reading grid files"
end

# ╔═╡ 3bc4c68a-360d-4dae-bf03-a630d5bc7f24
begin
	lon=[i for i=-179.75:0.5:179.75, j=-89.75:0.5:89.75]
	lat=[j for i=-179.75:0.5:179.75, j=-89.75:0.5:89.75]	
	fil_coeffs=joinpath(tempdir(),"half_degree_grid_coeffs.jld2")
	if !isfile(fil_coeffs)
		(f,i,j,w,j_f,j_x,j_y)=InterpolationFactors(Γ,vec(lon),vec(lat))
		jldsave(fil_coeffs;f,i,j,w,j_f,j_x,j_y)
	end
	C=load(fil_coeffs)
	"Done with interpolation coefficients"
end

# ╔═╡ 6eb8c05c-2624-4535-906b-3e38d502c997
begin	
	function ▶a(x::String,t::Int) #Γ,lon,C from calling scope
		y=basename(x)
		X=read_nctiles(x,y,γ;I=(:,:,1,t))[:,1]
		X[findall(Γ.hFacC[:,1] .== 0.)] .= NaN

		X=Interpolate(X,C["f"],C["i"],C["j"],C["w"])
		reshape(X,size(lon))
	end

	function ▶b(x::String,t::Int) #Γ,lon,C from calling scope
		y=basename(x)
		X=read_nctiles(x,y,γ;I=(:,:,t))[:,1]
		X[findall(Γ.hFacC[:,1] .== 0.)] .= NaN

		X=Interpolate(X,C["f"],C["i"],C["j"],C["w"])
		reshape(X,size(lon))
	end
	
	(▶a,▶b)
end

# ╔═╡ 2965efbb-a5b7-47eb-8a3b-0e8024d4b15f
md"""### Ocean Bottom Depth"""

# ╔═╡ d0c9aace-6fd1-4f04-92dc-4eede6cd47eb
function ▶c(x::String,t::Int) #Γ,lon,C from calling scope
	X=Interpolate(Γ.Depth,C["f"],C["i"],C["j"],C["w"])
	reshape(X,size(lon))
end

# ╔═╡ b8c2d4e6-7cbe-467c-bc67-4dc14b7eaf72
md"""### Remotely Sensed Reflectances"""

# ╔═╡ 9fbe6a83-9ed4-4e13-ae9b-d2742d60cc09
begin
	wvbd_out=Float64.([412, 443, 490, 510, 555, 670])
	wvbd_in=Float64.([400,425,450,475,500,525,550,575,600,625,650,675,700])
	Rirr=[["Rirr00$i" for i in 1:9];["Rirr0$i" for i in 10:13]]

	function ▶d(x::String,t::Int) #Γ,lon,C from calling scope
		tmp=fill(0.0,size(lon)...,13)
		Rrs=fill(0.0,size(lon)...,6)
		Chla=fill(0.0,size(lon)...)

		for ii in 1:length(Rirr)
			x = joinpath(pth0,"IrradianceReflectance",Rirr[ii],Rirr[ii])
			X=read_nctiles(x,Rirr[ii],γ;I=(:,:,t))[:,1]
			X[findall(Γ.hFacC[:,1] .== 0.)] .= NaN
			X=Interpolate(X,C["f"],C["i"],C["j"],C["w"])
			tmp[:,:,ii]=reshape(X,size(lon))
		end

		for i in 1:size(tmp,1), j in 1:size(tmp,2)
			Rrs[i,j,:].=RemotelySensedReflectance(tmp[i,j,:],wvbd_in,wvbd_out)
			Chla[i,j]=RrsToChla(Rrs[i,j,:])
		end
		
		Rrs,Chla
	end
end

# ╔═╡ 9828a5be-83f9-49f8-8f26-6f65451d021b
md"""### Chl (all phyto-plankton, 0-50m)""" 

# ╔═╡ ad30d090-7e49-4ad6-adbd-bcb77cbcb2f1
begin
	w50m=similar(Γ.hFacC); w50m.=0.0
	tot_w50m=similar(Γ.hFacC[:,1]); tot_w50m.=0.0
	
	for f in eachindex(w50m)
		ni=Γ.hFacC.fSize[f[1]][1]
		nj=Γ.hFacC.fSize[f[1]][2]
		for i in 1:ni, j in 1:nj
			Γ.RC[f[2]].>-50.0 ? w50m[f][i,j]=Γ.hFacC[f][i,j]*Γ.DRF[f[2]] : nothing
			Γ.RC[f[2]].>-50.0 ? tot_w50m[f[1]][i,j]+=Γ.hFacC[f][i,j]*Γ.DRF[f[2]] : nothing
		end
	end

	for f in eachindex(w50m)
		w50m[f]=w50m[f]./tot_w50m[f[1]]
	end
	
	w50m[findall(isnan.(w50m))]=0.0
	
	"Done with ▶e configuration"
end

# ╔═╡ a009b061-839c-4820-bfe8-d73dcf909dfa
begin
	function ▶e(x::String,t::Int) #Γ,lon,C from calling scope
		Chla=fill(0.0,size(lon)...)

		for ii in 1:35
			ii<10 ? nam="Chl0$(ii)" : nam="Chl$(ii)"
			x = joinpath(pth0,"Chlorophyll",nam,nam)
			X=read_nctiles(x,nam,γ;I=(:,:,:,t))[:,1]
			
			Y=read( sum(write(X).*write(w50m),dims=3) , Γ.XC)
			Y[findall(Γ.hFacC[:,1] .== 0.)] .= NaN
			Y=Interpolate(Y,C["f"],C["i"],C["j"],C["w"])
			Y=reshape(Y,size(lon))
			
			Chla=Chla+Y
		end
		
		Chla
	end
end

# ╔═╡ a3e1631a-56cc-428f-82ad-96e20d8e70d5
md"""### Surface PAR, Euphotic Depth"""

# ╔═╡ abb34421-c340-4f22-a75f-4a23d6af03dc
function ▶f(x::String,t::Int) #Γ,lon,C from calling scope
	y=basename(x)
	X=read_nctiles(x,y,γ;I=(:,:,1,t))[:,1]
	X[findall(Γ.hFacC[:,1] .== 0.)] .= NaN
	
	#convert from uEin/m^2/s to Ein/m^2/day
	X=86400.0/1000000.0*X

	X=Interpolate(X,C["f"],C["i"],C["j"],C["w"])
	reshape(X,size(lon))
end

# ╔═╡ fc858f0d-fed8-4f8d-a816-46b007157207
function ▶g(x::String,t::Int) #Γ,lon,C from calling scope
	y=basename(x)
	X=read_nctiles(x,y,γ;I=(:,:,:,t))
	X=write(X)
	
	ni=Γ.XC.grid.ioSize[1]
	nj=Γ.XC.grid.ioSize[2]
	euphotic_depth=fill(0.0,ni,nj)
	for i in 1:ni, j in 1:nj
		tmp1=X[i,j,:]; tmp2=0.01*tmp1[1]
		if tmp1[1]>0.0
			k1=minimum(findall(tmp1.<tmp2))
			tmp3=(tmp1[k1-1]-tmp2)/((tmp1[k1-1]-tmp1[k1]))
			euphotic_depth[i,j]=(1-tmp3)*Γ.RF[k1-1]+tmp3*Γ.RF[k1]
		end
	end	
	X=-read(euphotic_depth,Γ.XC)

	X[findall(Γ.hFacC[:,1] .== 0.)] .= NaN
	X=Interpolate(X,C["f"],C["i"],C["j"],C["w"])
	reshape(X,size(lon))
end

# ╔═╡ ec8fb91c-1b11-4a81-8fac-08a5f3b4e9b3
md"""## Main Computation Loops

### Individual Variables
"""

# ╔═╡ 3d541a75-7d88-4343-aed1-c4639f67208d
begin
	🏁1=true
	PlutoUI.with_terminal() do
		input=["THETA","SALT","MXLDEPTH","EXFwspee","GGL90TKE"]
		output=["SST","SSS","MLD","WindSpeed","TKE"]
		▶▶ = [▶a,▶a,▶b,▶b,▶a]
		pth=fill("PhysicalOceanography",5) 
		
		for ii in 1:length(input)
			X = joinpath(pth0,pth[ii],input[ii],input[ii])
			▶▶▶ = ▶▶[ii]
			clim=fill(0.0,size(lon)...,12)
			[clim[:,:,t]=▶▶▶(X,t) for t in 1:12]

			fil2=joinpath(tempdir(),output[ii]*"_clim.jld2")
			jldsave(fil2;clim)

			println("Done with pre-processing $(input[ii]) the 12 months")
		end
	end
end

# ╔═╡ 7e912200-f852-457d-8fea-c05997f8ae1d
begin
	🏁2=true
	PlutoUI.with_terminal() do
		input=["OceanDepth","Chl","PARF","PARF"]
		output=["OceanDepth","Chl","PAR","EuphoticDepth"]
		▶▶ = [▶c,▶e,▶f,▶g]
		pth=[fill("PhysicalOceanography",1);"Chlorophyll";fill("IrradianceReflectance",2)] 
		
		for ii in 1:length(input)
			X = joinpath(pth0,pth[ii],input[ii],input[ii])
			▶▶▶ = ▶▶[ii]
			clim=fill(0.0,size(lon)...,12)
			[clim[:,:,t]=▶▶▶(X,t) for t in 1:12]

			fil2=joinpath(tempdir(),output[ii]*"_clim.jld2")
			jldsave(fil2;clim)

			println("Done with pre-processing $(input[ii]) the 12 months")
		end
	end
end

# ╔═╡ 9006655e-3b58-47d3-9baf-6b287e95ece1
md""" ### Remotely Sensed Reflectances"""

# ╔═╡ cffa9cb4-8ddf-4fb4-9874-cd1817fa2aad
🏁3 = let
	clim_Chl=fill(0.0,size(lon)...,12)
	clim_Rrs=[fill(0.0,size(lon)...,12) for j in 1:13]

	for t in 1:12
		(Rrs,Chl)=▶d("Chl",t)
		clim_Chl[:,:,t].=Chl
		[clim_Rrs[j][:,:,t]=Rrs[:,:,j] for j in 1:6]
	end

	fil_Chl=joinpath(tempdir(),"Chl_from_Rrs_clim.jld2")
	jldsave(fil_Chl;clim_Chl)

	fil_Rrs=joinpath(tempdir(),"Rrs_clim.jld2")
	jldsave(fil_Rrs;clim_Rrs)

end

# ╔═╡ 7015ddff-66b3-4c2c-bf9a-ccfc18c3693e
md""" ### Chl Concentration"""

# ╔═╡ 6e389492-a103-419e-8513-9f22f1be4989
🏁4 = let
	clim=fill(0.0,size(lon)...,12)
	[clim[:,:,t]=▶e("Chl",t) for t in 1:12]

	fil_Chl=joinpath(tempdir(),"Chl_clim.jld2")
	jldsave(fil_Chl;clim)
end

# ╔═╡ 4d8ded5f-cfd1-45ac-815c-b0358806e846
md"""## File Creation Loops

### Individual Variables
"""

# ╔═╡ 99fc0197-02d5-411b-a703-375fe71ca983
📁1 = let
	🏁1, 🏁2, 🏁4
	meta0=[(name="SST",units="degC",long_name="Sea Surface Temperature (0-10m)")
		(name="SSS",units="psu",long_name="Sea Surface Salinity (0-10m)")
		(name="MLD",units="m",long_name="Mixed Layer Depth")
		(name="WindSpeed",units="m/s",long_name="Surface Wind Speed")
		(name="TKE",units="m^2/s^2",long_name="Turbulent Kinetic Energy (0-10m, from GGL90)")
		(name="OceanDepth",units="m",long_name="Bathymetry")
		(name="Chl",units="mg Chl",long_name="Chorophyll a (0-50m)")
		(name="PAR",units="m",long_name="PAR at Sea Surface")
		(name="EuphoticDepth",units="m",long_name="Euphotic Depth")]

	for ii in 1:length(meta0)
		nam=meta0[ii].name
		data=load(joinpath(tempdir(),nam*"_clim.jld2"))["clim"]
		outputfile=joinpath(tempdir(),nam*"_clim.nc")

		README = ["Source: Gael Forget","Product: CBIOMES-global climatology","Version: alpha"]
		meta=(outputfile=outputfile , README=README , yearrange=(1992,2011))
		meta=merge(meta0[ii],meta)
		
		nct=NCT(lon,lat,data,meta)
		write(nct)
	end	
end

# ╔═╡ b59d16e0-477c-44f3-8e0a-6fcedadd0a23
md"""### Remotely Sensed Reflectances"""

# ╔═╡ d30cf2f5-5918-458c-be0c-900700642974
📁2 = let
	🏁3

	meta1=[(name="Chl",units="mg/m^3",long_name="Chlorophyll a Concentration")
		(name="Rrs412",units="sr-1",long_name="Remote Sensing Reflectance at 412 nm")
		(name="Rrs443",units="sr-1",long_name="Remote Sensing Reflectance at 443 nm")
		(name="Rrs490",units="sr-1",long_name="Remote Sensing Reflectance at 490 nm")
		(name="Rrs510",units="sr-1",long_name="Remote Sensing Reflectance at 510 nm")
		(name="Rrs555",units="sr-1",long_name="Remote Sensing Reflectance at 555 nm")
		(name="Rrs670",units="sr-1",long_name="Remote Sensing Reflectance at 670 nm")]

	for ii in 1:length(meta1)
		nam=meta1[ii].name
		if ii==1
			data=load(joinpath(tempdir(),"Chl_from_Rrs_clim.jld2"))["clim_Chl"]
			outputfile=joinpath(tempdir(),"Chl_from_Rrs_clim.nc")
		else
			data=load(joinpath(tempdir(),"Rrs_clim.jld2"))["clim_Rrs"][ii-1]
			outputfile=joinpath(tempdir(),nam*"_clim.nc")
		end

		README = ["Source: Gael Forget","Product: CBIOMES-global climatology","Version: alpha"]
		meta=(outputfile=outputfile , README=README , yearrange=(1992,2011))
		meta=merge(meta1[ii],meta)
		
		nct=NCT(lon,lat,data,meta)
		write(nct)
	end	
end

# ╔═╡ 4d29a7f6-6071-4dd3-b422-3d9c513b0d86
md"""### Combine NetCDF files"""

# ╔═╡ 8b202092-62fd-46a1-a184-5a7a3cf40c46
"""
    CBIOMES_combine_files(fil_out::String)

Take all files generated earlier and combine them into one.
"""
function CBIOMES_combine_files(fil_out::String)

list_in=["Chl","EuphoticDepth","MLD","OceanDepth","PAR",
    "Rrs412","Rrs443","Rrs490","Rrs510","Rrs555","Rrs670",
    "SSS","SST","TKE","WindSpeed"]

!isfile(fil_out) ? cp(joinpath(tempdir(),"SST_clim.nc"),fil_out) : nothing

for ii in list_in
    ds = NCTiles.Dataset(fil_out,"a")
    if !haskey(ds,ii)
        fil_in=joinpath(tempdir(),ii*"_clim.nc")
        ds_in = NCTiles.Dataset(fil_in,"r")
        tmp11=ds_in[ii][:]
        u=ds_in[ii].attrib["units"]
        ln=ds_in[ii].attrib["long_name"]
        close(ds_in)

        ##

        v = NCTiles.defVar(ds,ii,Float64,("lon","lat","t"), 
        attrib = Dict("units" => u, "long_name" => ln))
        v[:] = tmp11
    end
    close(ds)
end

ds = NCTiles.Dataset(fil_out,"a")
ds.attrib["description,1"]="Source: Gael Forget"
ds.attrib["description,2"]="Product: CBIOMES-global climatology"
ds.attrib["description,3"]="Version: alpha"
if haskey(ds.attrib,"A")
    NCTiles.delete!(ds.attrib,"description")
    NCTiles.delete!(ds.attrib,"A")
    NCTiles.delete!(ds.attrib,"B")
end
close(ds)

"all set"
end

# ╔═╡ 0d28aa82-fa6f-466d-ba08-25016bce90a6
📁3 = let
	📁1, 📁2
	
	fil_out=joinpath(tempdir(),"CBIOMES-global-alpha-climatology.nc")
	rm(fil_out,force=true)
	CBIOMES_combine_files(fil_out)
end

# ╔═╡ 12a4b433-3722-4ebc-b75f-027539d979b0
md"""## Appendices


1. single variable test and comparison with earlier-code result
1. inspect a sample input file, and retrieve meta-data from it

"""

# ╔═╡ 8edb7bc3-c140-4e32-b626-d6333167ea06
let
	fil0="PhysicalOceanography/THETA/THETA"
	fil1=joinpath(pth0,fil0)

	SST=▶a(fil1,1)

	fig = Mkie.Figure(resolution = (900,600), backgroundcolor = :grey95)
	ax = Mkie.Axis(fig[1,1], title="SST in degC",xlabel="longitude",ylabel="latitude")
	hm1=Mkie.contourf!(ax,lon[:,1],lat[1,:],SST,levels = -2.0:2.0:34.0, tickfont = (4, :black))
	xlims!(ax, (-180.0, 180.0)); ylims!(ax, (-90.0, 90.0))
	Mkie.Colorbar(fig[1,2], hm1, height = Mkie.Relative(0.65))

md"""### Single Month Test

!!! note
	Plot below shows SST for January after interpolation to the half-degree resolution grid.

$fig
"""
end

# ╔═╡ d5ed2bcc-8034-4994-aaa7-2c5fbdfda37e
let
	🏁1
	
	t=11
	clim=load(joinpath(tempdir(),"SST_clim.jld2"))["clim"]
	pml=NCTiles.ncread("gridded_darwin_montly_clim_360_720_ver_0_2_6.nc","SST")

	RMSE=dropdims(sqrt.(sum((clim-pml).^2,dims=3)/12),dims=3)

	fig = Mkie.Figure(resolution = (900,600), backgroundcolor = :grey95)
	ax = Mkie.Axis(fig[1,1], title="log10(RMSE)",xlabel="longitude",ylabel="latitude")
	hm1=Mkie.contourf!(ax,lon[:,1],lat[1,:],log10.(RMSE), tickfont = (4, :black))
	xlims!(ax, (-180.0, 180.0)); ylims!(ax, (-90.0, 90.0))
	Mkie.Colorbar(fig[1,2], hm1, height = Mkie.Relative(0.65))

	md"""
	!!! note
	    Plot below shows room mean square difference from an earlier version generated with different software.
	
	$(fig)
	"""
end

# ╔═╡ 620ceef5-99f9-48d2-9079-3f60109d2005
begin
	📁3	
	fil_out=joinpath(tempdir(),"CBIOMES-global-alpha-climatology.nc")
	NCTiles.NCDataset(fil_out,"r")
end

# ╔═╡ 1318a437-81b2-4982-8f58-9273a5cc6df9
begin
	fil0="PhysicalOceanography/THETA/THETA"
	NCTiles.NCDataset(joinpath(pth0,fil0*".0058.nc"),"r")
end

# ╔═╡ 81f1ace0-cf04-4690-82c4-0b54219a0c1d
NCTiles.NCDataset(fil_out,"r")

# ╔═╡ 9cd5900d-32f3-4a07-a4a2-071a0615c527
let
	ncvars,ncdims,fileatts = readncfile(joinpath(pth0,fil0*".0058.nc"));
	#PlutoUI.with_terminal() do
	println("\n============ A look at "*fil0*".0058.nc using readncfile ========= \n\n")
	println("- Variables defined in ncvars: \n ")
	println.(keys(ncvars))
	println("\n")
	println("- Summary for some variables: \n ")
	show(ncvars["tim"])
	println("\n")
	show(ncvars["THETA"])
	println("\n")
end

# ╔═╡ Cell order:
# ╟─979a1b66-0aaa-4b8f-868c-85db76d3fb99
# ╟─aa41a861-7eca-4dda-92e8-bd8e96f2a21e
# ╟─5fcb77a8-b234-45d9-a4bd-02bd4eb8fbd9
# ╟─ad9d39a2-65b3-4fd9-aecd-f360ac34df0f
# ╟─4e55fa2e-3193-11ec-3a50-2defc49cc54d
# ╟─d7d0fc22-0bd7-4c4b-9249-8084d75bb33e
# ╟─7a1058e4-3813-4b90-a9ee-ad5ff4eb6d6f
# ╟─3bc4c68a-360d-4dae-bf03-a630d5bc7f24
# ╟─6eb8c05c-2624-4535-906b-3e38d502c997
# ╟─2965efbb-a5b7-47eb-8a3b-0e8024d4b15f
# ╟─d0c9aace-6fd1-4f04-92dc-4eede6cd47eb
# ╟─b8c2d4e6-7cbe-467c-bc67-4dc14b7eaf72
# ╟─9fbe6a83-9ed4-4e13-ae9b-d2742d60cc09
# ╟─9828a5be-83f9-49f8-8f26-6f65451d021b
# ╟─ad30d090-7e49-4ad6-adbd-bcb77cbcb2f1
# ╟─a009b061-839c-4820-bfe8-d73dcf909dfa
# ╟─a3e1631a-56cc-428f-82ad-96e20d8e70d5
# ╟─abb34421-c340-4f22-a75f-4a23d6af03dc
# ╟─fc858f0d-fed8-4f8d-a816-46b007157207
# ╟─ec8fb91c-1b11-4a81-8fac-08a5f3b4e9b3
# ╟─3d541a75-7d88-4343-aed1-c4639f67208d
# ╟─7e912200-f852-457d-8fea-c05997f8ae1d
# ╟─9006655e-3b58-47d3-9baf-6b287e95ece1
# ╟─cffa9cb4-8ddf-4fb4-9874-cd1817fa2aad
# ╟─7015ddff-66b3-4c2c-bf9a-ccfc18c3693e
# ╟─6e389492-a103-419e-8513-9f22f1be4989
# ╟─4d8ded5f-cfd1-45ac-815c-b0358806e846
# ╟─99fc0197-02d5-411b-a703-375fe71ca983
# ╟─b59d16e0-477c-44f3-8e0a-6fcedadd0a23
# ╟─d30cf2f5-5918-458c-be0c-900700642974
# ╟─4d29a7f6-6071-4dd3-b422-3d9c513b0d86
# ╟─8b202092-62fd-46a1-a184-5a7a3cf40c46
# ╟─0d28aa82-fa6f-466d-ba08-25016bce90a6
# ╟─12a4b433-3722-4ebc-b75f-027539d979b0
# ╟─8edb7bc3-c140-4e32-b626-d6333167ea06
# ╟─d5ed2bcc-8034-4994-aaa7-2c5fbdfda37e
# ╟─620ceef5-99f9-48d2-9079-3f60109d2005
# ╟─1318a437-81b2-4982-8f58-9273a5cc6df9
# ╟─81f1ace0-cf04-4690-82c4-0b54219a0c1d
# ╟─9cd5900d-32f3-4a07-a4a2-071a0615c527

### A Pluto.jl notebook ###
# v0.17.1

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

# ╔═╡ d123161e-49f1-11ec-1c1b-51871624545d
begin
	using Pkg
	Pkg.activate()
	
	using JLD2, MeshArrays, OceanStateEstimation, PlutoUI, Statistics
	import CairoMakie as Mkie
	"Done with packages"
end

# ╔═╡ 0ca1e0f5-9656-44d5-99b7-6d8f2a58f43a
readdir(joinpath(tempdir(),"ECCO_transport_lines"))

# ╔═╡ cc3e9d0c-8b71-432b-b68d-a8b832ca5f26
md"""# Tansport Lines For ECCO Analysis"""

# ╔═╡ 39924391-38ce-46a1-877f-80a7975340a0
begin
	pth=MeshArrays.GRID_LLC90
	γ=GridSpec("LatLonCap",pth)
	Γ=GridLoad(γ;option="full")
	#LC=LatitudeCircles(-89.0:89.0,Γ)
	"Done with grid"
end

# ╔═╡ 26299781-c3f7-493e-8ae8-e11aa2eb726e
md"""## Interactive Visualization"""

# ╔═╡ 2507c335-2886-42b5-b6b8-63279a2d60fe
begin
	lonPairs=[]    
	latPairs=[]    
	namPairs=[]    
	
	push!(lonPairs,[-173 -164]); push!(latPairs,[65.5 65.5]); push!(namPairs,"Bering Strait");
	push!(lonPairs,[-5 -5]); push!(latPairs,[34 40]); push!(namPairs,"Gibraltar");
	push!(lonPairs,[-81 -77]); push!(latPairs,[28 26]); push!(namPairs,"Florida Strait");
	push!(lonPairs,[-81 -79]); push!(latPairs,[28 22]); push!(namPairs,"Florida Strait W1");
	push!(lonPairs,[-76 -76]); push!(latPairs,[21 8]); push!(namPairs,"Florida Strait S1");
	push!(lonPairs,[-77 -77]); push!(latPairs,[26 24]); push!(namPairs,"Florida Strait E1");
	push!(lonPairs,[-77 -77]); push!(latPairs,[24 22]); push!(namPairs,"Florida Strait E2");
	push!(lonPairs,[-65 -50]); push!(latPairs,[66 66]); push!(namPairs,"Davis Strait");
	push!(lonPairs,[-35 -20]); push!(latPairs,[67 65]); push!(namPairs,"Denmark Strait");
	push!(lonPairs,[-16 -7]); push!(latPairs,[65 62.5]); push!(namPairs,"Iceland Faroe");
	push!(lonPairs,[-6.5 -4]); push!(latPairs,[62.5 57]); push!(namPairs,"Faroe Scotland");
	push!(lonPairs,[-4 8]); push!(latPairs,[57 62]); push!(namPairs,"Scotland Norway");
	push!(lonPairs,[-68 -63]); push!(latPairs,[-54 -66]); push!(namPairs,"Drake Passage");
	push!(lonPairs,[103 103]); push!(latPairs,[4 -1]); push!(namPairs,"Indonesia W1");
	push!(lonPairs,[104 109]); push!(latPairs,[-3 -8]); push!(namPairs,"Indonesia W2");
	push!(lonPairs,[113 118]); push!(latPairs,[-8.5 -8.5]); push!(namPairs,"Indonesia W3");
	push!(lonPairs,[118 127 ]); push!(latPairs,[-8.5 -15]); push!(namPairs,"Indonesia W4");
	push!(lonPairs,[127 127]); push!(latPairs,[-25 -68]); push!(namPairs,"Australia Antarctica");
	push!(lonPairs,[38 46]); push!(latPairs,[-10 -22]); push!(namPairs,"Madagascar Channel");
	push!(lonPairs,[46 46]); push!(latPairs,[-22 -69]); push!(namPairs,"Madagascar Antarctica");
	push!(lonPairs,[20 20]); push!(latPairs,[-30 -69.5]); push!(namPairs,"South Africa Antarctica");
	push!(lonPairs,[-76 -72]); push!(latPairs,[21 18.5]); push!(namPairs,"Florida Strait E3");
	push!(lonPairs,[-72 -72]); push!(latPairs,[18.5 10]); push!(namPairs,"Florida Strait E4");
end

# ╔═╡ 21d8bccb-a065-4965-806d-5757f8177dbb
@bind ii Select(1:length(namPairs))

# ╔═╡ 94b8ba05-dfb8-4075-a260-7968e8fdd78f
md"""## Main Computation Loop"""

# ╔═╡ 6cc62cf0-cb30-4d93-aad6-2ab16f60f95f
let
	pth_trsp=joinpath(tempdir(),"ECCO_transport_lines")
	!isdir(pth_trsp) ? mkdir(pth_trsp) : nothing
	
	for ii in 1:length(lonPairs)
		lons=Float64.(lonPairs[ii])
		lats=Float64.(latPairs[ii])
		name=namPairs[ii]
		Trsct=Transect(name,lons,lats,Γ)
		jldsave(joinpath(tempdir(),"ECCO_transport_lines","$(Trsct.name).jld2"),
			tabC=Trsct.tabC,tabW=Trsct.tabW,tabS=Trsct.tabS); 
	end
end

# ╔═╡ 25144e1b-21fc-4cc9-b63d-7b26eab1a673
md"""## Underlying Functions"""

# ╔═╡ f140908c-a8cd-4dd9-8244-ba5d08689559
#bring end points to the equator -> define 3D rotation matrix
function rotate_points(lons,lats)
	#... and of end points
	x0=cosd.(lats).*cosd.(lons)
	y0=cosd.(lats).*sind.(lons)
	z0=sind.(lats)

	#get the rotation matrix:
	#1) rotate around x axis to put first point at z=0
	theta=atan(-z0[1],y0[1])
	R1=[[1;0;0] [0;cos(theta);sin(theta)] [0;-sin(theta);cos(theta)]]
	tmp0=[x0;y0;z0]; tmp1=R1*tmp0; x1=tmp1[1,:]; y1=tmp1[2,:]; z1=tmp1[3,:]
	x0=x1; y0=y1; z0=z1
	#2) rotate around z axis to put first point at y=0
	theta=atan(x0[1],y0[1])
	R2=[[cos(theta);sin(theta);0] [-sin(theta);cos(theta);0] [0;0;1]]
	tmp0=transpose([x0 y0 z0]); tmp1=R2*tmp0; x1=tmp1[1,:]; y1=tmp1[2,:]; z1=tmp1[3,:]
	x0=x1; y0=y1; z0=z1
	#3) rotate around y axis to put second point at z=0
	theta=atan(-z0[2],-x0[2])
	R3=[[cos(theta);0;-sin(theta)] [0;1;0] [sin(theta);0;cos(theta)]]
	tmp0=transpose([x0 y0 z0]); tmp1=R3*tmp0; x1=tmp1[1,:]; y1=tmp1[2,:]; z1=tmp1[3,:]
	x0=x1; y0=y1; z0=z1

	x0,y0,z0,R3*R2*R1
end

# ╔═╡ 2b7d7df1-7309-4645-ad52-c46e29c7d610
function rotate_XCYC(Γ,R)
	#3D carthesian coordinates:
	lon=Γ.XC; lat=Γ.YC
	x=cosd.(lat)*cosd.(lon)
	y=cosd.(lat)*sind.(lon)
	z=sind.(lat)

	#rotation R:
	tmpx=γ.write(x); tmpy=γ.write(y); tmpz=γ.write(z)
	tmp1=findall((!isnan).(tmpx))
	tmpx2=tmpx[tmp1]; tmpy2=tmpy[tmp1]; tmpz2=tmpz[tmp1]
	tmp2=[tmpx2';tmpy2';tmpz2']
	
	tmp3=R*tmp2

	tmpx2=tmp3[1,:]; tmpy2=tmp3[2,:]; tmpz2=tmp3[3,:]
	tmpx[tmp1]=tmpx2; tmpy[tmp1]=tmpy2; tmpz[tmp1]=tmpz2
	x=γ.read(tmpx,Γ.XC); y=γ.read(tmpy,Γ.XC); z=γ.read(tmpz,Γ.XC)
	
	x,y,z
end

# ╔═╡ 963e421c-43fb-43d3-b667-1b9912f940b8
begin
	lons=Float64.(lonPairs[ii])
	lats=Float64.(latPairs[ii])
	name=namPairs[ii]

	x0,y0,z0,R=rotate_points(lons,lats)
	x,y,z=rotate_XCYC(Γ,R)
	mskCint=1.0*(z.>0)
	mskCedge,mskWedge,mskSedge=MeshArrays.edge_mask(mskCint)
	mskCedge,mskWedge,mskSedge=MeshArrays.shorter_paths!((x,y,z),(x0,y0,z0),(mskCedge,mskWedge,mskSedge))
	"Done with transport line masks"
end

# ╔═╡ 5e3fb2fa-4e6a-4cb3-b76f-d02011948530
function edge_mask(mskCint::MeshArray)
    mskCint=1.0*mskCint

    #treat the case of blank tiles:
    #mskCint[findall(RAC.==0)].=NaN
    
    mskCplus=exchange(mskCint)

    #edge tracer mask
    mskCedge=similar(mskCint)
    for i in eachindex(mskCedge)
        tmp1=mskCplus[i]
        tmp2=tmp1[2:end-1,1:end-2]+tmp1[2:end-1,3:end]+
            tmp1[1:end-2,2:end-1]+tmp1[3:end,2:end-1]
        mskCedge[i]=1.0*(tmp2.>0).*(tmp1[2:end-1,2:end-1].==0.0)
    end

    #edge velocity mask:
    mskWedge=similar(mskCint)
    mskSedge=similar(mskCint)
    for i in eachindex(mskCedge)
        mskWedge[i]=mskCplus[i][2:end-1,2:end-1] - mskCplus[i][1:end-2,2:end-1]
        mskSedge[i]=mskCplus[i][2:end-1,2:end-1] - mskCplus[i][2:end-1,1:end-2]
    end

    #treat the case of blank tiles:
    #mskCedge[findall(isnan.(mskCedge))].=0.0
    #mskWedge[findall(isnan.(mskWedge))].=0.0
    #mskSedge[findall(isnan.(mskSedge))].=0.0

    return mskCedge,mskWedge,mskSedge
end

# ╔═╡ d5e26bc2-5b15-458d-98eb-0fc5d4977b9d
function shorter_paths!(xyz,xyz0,msk_in)
	(x0,y0,z0)=xyz0[:]
	(x,y,z)=xyz[:]

	#split in two segments:
	theta=zeros(2)
	theta[1]=atan(y0[1],x0[1])
	theta[2]=atan(y0[2],x0[2])

	γ=msk_in[1],grid
	tmpx=γ.write(x); tmpy=γ.write(y); tmpz=γ.write(z);
	tmptheta=atan.(tmpy,tmpx)
	if theta[2]<0;
		tmp00=findall(tmptheta.<=theta[2])
		tmptheta[tmp00].=tmptheta[tmp00]+2*pi
		theta[2]=theta[2]+2*pi
	end

	msk_out=[]
	for kk in 1:3
		#select field to treat:
		mm=msk_in[kk]
		#select the shorther segment:
		tmpm=γ.write(mm)
       if theta[2]-theta[1]<=pi
            tmpm[findall( (tmptheta.>theta[2]).|(tmptheta.<theta[1]) )].=0.0
        else
            tmpm[findall( (tmptheta.<=theta[2]).&(tmptheta.>=theta[1]) )].=0.0
        end
        mm=γ.read(tmpm,mm);
        #store result:
		push!(msk_out,mm)
    end
	msk_out[:]
end

# ╔═╡ 75c98143-dd36-4446-adfa-c440fdf8aab1
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

# ╔═╡ ed36a2a5-44ea-43a7-a3bd-13f234b6580d
begin	
	λ=setup_interp(Γ)
	
	DD=Interpolate(λ.μ*mskCedge,λ.f,λ.i,λ.j,λ.w)
	DD=reshape(DD,size(λ.lon))
	#DD[findall(DD.==0.0)].=NaN
	
	fig = Mkie.Figure(resolution = (900,600), backgroundcolor = :grey95)
	ax = Mkie.Axis(fig[1,1], title="mskCdge",xlabel="longitude",ylabel="latitude")
	hm1=Mkie.heatmap!(ax,λ.lon[:,1],λ.lat[1,:],DD) #,colorrange=(-1.0,1.0).*0.2)
	Mkie.scatter!(ax,lons[:],lats[:],color=:red)
	#Mkie.Colorbar(fig[1,2], hm1, height = Mkie.Relative(0.65))
	fig
end

# ╔═╡ c70e3a7a-5512-4076-9aae-2d8bd86367e7
#joinpath(tempdir(),"ECCO_transport_lines","$(name).jld2")
        function MskToTab(msk::MeshArray)
          n=Int(sum(msk .!= 0)); k=0
          tab=Array{Int,2}(undef,n,4)
          for i=1:msk.grid.nFaces
            a=msk.f[i]
            b=findall( a .!= 0)
            for ii in eachindex(b)
              k += 1
              tab[k,:]=[i,b[ii][1],b[ii][2],a[b[ii]]]
            end
          end
          return tab
        end


# ╔═╡ Cell order:
# ╠═0ca1e0f5-9656-44d5-99b7-6d8f2a58f43a
# ╟─cc3e9d0c-8b71-432b-b68d-a8b832ca5f26
# ╠═d123161e-49f1-11ec-1c1b-51871624545d
# ╟─39924391-38ce-46a1-877f-80a7975340a0
# ╟─26299781-c3f7-493e-8ae8-e11aa2eb726e
# ╟─2507c335-2886-42b5-b6b8-63279a2d60fe
# ╟─21d8bccb-a065-4965-806d-5757f8177dbb
# ╟─ed36a2a5-44ea-43a7-a3bd-13f234b6580d
# ╟─963e421c-43fb-43d3-b667-1b9912f940b8
# ╟─94b8ba05-dfb8-4075-a260-7968e8fdd78f
# ╠═6cc62cf0-cb30-4d93-aad6-2ab16f60f95f
# ╟─25144e1b-21fc-4cc9-b63d-7b26eab1a673
# ╟─f140908c-a8cd-4dd9-8244-ba5d08689559
# ╟─2b7d7df1-7309-4645-ad52-c46e29c7d610
# ╟─5e3fb2fa-4e6a-4cb3-b76f-d02011948530
# ╟─d5e26bc2-5b15-458d-98eb-0fc5d4977b9d
# ╟─75c98143-dd36-4446-adfa-c440fdf8aab1
# ╟─c70e3a7a-5512-4076-9aae-2d8bd86367e7
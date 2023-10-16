
## 

module procs

using JLD2, MeshArrays, Statistics, OceanStateEstimation, Glob, TOML

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

##

function years_min_max(sol)
	year0=1992
	year1=2011
	if occursin("r3",sol)
		year1=2015
	elseif occursin("r4",sol)
		year1=2017
	elseif occursin("r5",sol)
		year1=2019
	elseif occursin("43y",sol)
		year0=1982
		year1=2023
	end
	return year0,year1
end

##

function parameters()

	pth=MeshArrays.GRID_LLC90
	γ=GridSpec("LatLonCap",pth)
	Γ=GridLoad(γ;option="full")
	#LC=LatitudeCircles(-89.0:89.0,Γ)

	##

	#ECCOdiags_add("release5")
	interpolation_setup()
	μ = land_mask(Γ)

	λ_file = joinpath(tempdir(),"interp_coeffs_halfdeg.jld2")
	if !isfile(λ_file)
		lon=[i for i=-179.75:0.5:179.75, j=-89.75:0.5:89.75]
		lat=[j for i=-179.75:0.5:179.75, j=-89.75:0.5:89.75]		
		(f,i,j,w)=InterpolationFactors(Γ,vec(lon),vec(lat))
		jldsave(λ_file; lon=lon, lat=lat, f=f, i=i, j=j, w=w)
	end
	λ = interpolation_setup(λ_file)

	##
	
	sol_list=glob("ECCOv4*_analysis",ScratchSpaces.ECCO)
	sol_list=[basename(i) for i in sol_list]

	fil_trsp=joinpath(ScratchSpaces.ECCO,"ECCOv4r2_analysis/trsp/trsp.jld2")
	ntr=length(load(fil_trsp,"single_stored_object"))
	list_trsp=[vec(load(fil_trsp,"single_stored_object"))[i].nam for i in 1:ntr] 
	list_trsp=[i[1:end-5] for i in list_trsp]

	pth_colors=joinpath(dirname(pathof(OceanStateEstimation)),"..","examples","ECCO")	
	clim_colors1=TOML.parsefile(joinpath(pth_colors,"clim_colors1.toml"))
	clim_colors2=TOML.parsefile(joinpath(pth_colors,"clim_colors2.toml"))


	pth_tmp01=joinpath(ScratchSpaces.ECCO,"ECCOv4r2_analysis")
	clim_files=climatology_files(pth_tmp01)
	clim_name=[split(basename(f),'.')[1] for f in clim_files]
	clim_longname=longname.(clim_name) 

	#"Done with listing solutions, file names, color codes"
	(γ=γ,Γ=Γ,λ=λ,μ=μ,sol_list=sol_list,list_trsp=list_trsp,
	clim_colors1=clim_colors1,clim_colors2=clim_colors2,
	clim_files=clim_files,clim_name=clim_name,clim_longname=clim_longname)
end

##

function glo(pth_out,nam,k,year0,year1)
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
	x=year0 .+ x./12.0

	(y=tmp,txt=txt,rng=rng,x=x)
end

function map(nammap,P,statmap,timemap,pth_out)
	ii=findall(P.clim_longname.==nammap)[1]
	nam=P.clim_name[ii]
	
	fil=joinpath(pth_out,P.clim_files[ii])
	if statmap!=="mon"
		tmp=load(fil,statmap)
	else
		tmp=load(fil,statmap)[:,timemap]
	end

	DD=Interpolate(P.λ.μ*tmp,P.λ.f,P.λ.i,P.λ.j,P.λ.w)
	DD=reshape(DD,size(P.λ.lon))
	#DD[findall(DD.==0.0)].=NaN
	statmap=="std" ? rng=P.clim_colors2[nam] : rng=P.clim_colors1[nam]
	levs=rng[1] .+collect(0.0:0.05:1.0)*(rng[2]-rng[1])

	ttl=P.clim_longname[ii]
	return P.λ,DD,levs,ttl
end

function TimeLat(namzm,pth_out,year0,year1,cmap_fac)
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

	x=year0 .+ x./12.0
	ttl="$(longname(namzm)) : Zonal Mean $(addon1)"
	return x,y,z,cmap_fac*levs,ttl,-90.0,90.0,year0,year1
end

function TimeLatAnom(namzmanom2d,pth_out,year0,year1,cmap_fac,k_zm2d,l0,l1,P)
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
		addon1=" -- at $(Int(round(P.Γ.RC[k_zm2d])))m "
	else
		z=fn(tmp[:,:])
		x=vec(0.5:size(tmp,2)); 
		addon1=""
	end

	dlat=2.0; y=vec(-90+dlat/2:dlat:90-dlat/2)
	nt=size(z,1)

	m0=(1992-year0)*12
	
	if true
		#a. subtract monthly mean
		ref1="1992-2011 monthy mean"
		for m in 1:12
			zmean=vec(mean(z[m0+m:12:m0+240,:],dims=1))
			[z[t,:]=z[t,:]-zmean for t in m:12:nt]
		end
	else
		#b. subtract time mean
		ref1="1992-2011 annual mean"
		zmean=vec(mean(z[m0+1:m0+240,:],dims=1))
		[z[t,:]=z[t,:]-zmean for t in 1:nt]
	end

	x=1992.0-m0/12.0 .+ x./12.0
	ttl="$(procs.longname(namzm)) -- minus $(ref1) $(addon1)"

	return x,y,z,cmap_fac*levs,ttl,y[l0],y[l1],year0,year1
end

fn_DepthTime(x)=transpose(x)	

function DepthTime(namzmanom,pth_out,facA,l_Tzm,year0,year1,k0,k1,P)
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
y=vec(P.Γ.RC)
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

x=year0 .+ x./12.0
ttl="$(longname(namzmanom)) -- minus $(ref1) $(addon1)"

return x,y,z,facA*levs,ttl,P.Γ.RC[k1],P.Γ.RC[k0]
end

end

##

module plots

	using CairoMakie, JLD2, RollingFunctions, Statistics

	to_range!(DD,levs::Tuple) = to_range!(DD,range(levs[1],levs[2],length=10))

	function to_range!(DD,levs)
		DD[findall(DD.<=levs[1])].=levs[1]+(levs[2]-levs[1])/100
		DD[findall(DD.>=levs[end])].=levs[end]-(levs[end]-levs[end-1])/100
	end

	function axtr1(ax,namtr,pth_out,list_trsp,year0,year1)
		fil_trsp=joinpath(pth_out,"trsp/trsp.jld2")

		itr=findall(list_trsp.==namtr)[1]
		tmp=vec(load(fil_trsp,"single_stored_object"))[itr]
		
		nt=size(tmp.val,2)
		x=vec(0.5:nt)
	
		txt=tmp.nam[1:end-5]
		val=1e-6*vec(sum(tmp.val,dims=1)[:])
		valsmo = runmean(val, 12)
	
		x=vec(0.5:nt)
		x=year0 .+ x./12.0

		hm1=lines!(ax,x,val,label="ECCO estimate")
		valsmo[1:5].=NaN
		valsmo[end-4:end].=NaN
		lines!(ax,x,valsmo,linewidth=4.0,color=:red)
		xlims!(ax,(year0,year1))
	end

	function transport(namtrs,ncols,pth_out,list_trsp,year0,year1)
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
				xticks=(year0:4:year1),ylabel="transport, in Sv")
			axtr1(ax1,namtrs[na],pth_out,list_trsp,year0,year1)
		end
		#ylims!(ax1,rng)
		fig1
	end

	function figov1(pth_out,kk,low1,year0,year1)
		fil=joinpath(pth_out,"overturn/overturn.jld2")
		tmp=-1e-6*load(fil,"single_stored_object")
	
		nt=size(tmp,3)
		x=vec(0.5:nt)
		x=year0 .+ x./12.0
		lats=vec(-89.0:89.0)

		fig1 = Figure(resolution = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1],ylabel="Sv",
			title="Global Overturning, in Sv, at kk=$(kk)",
			xticks=(year0:4:year1))
		for ll in 115:10:145
			ov=tmp[ll,kk,:]
			ov=runmean(ov, 12)
			ov[1:5].=NaN
			ov[end-4:end].=NaN
			hm1=lines!(x,ov,label="$(lats[ll])N")
		end
		xlims!(ax1,(year0,year1))
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

	function glo(gl1,year0,year1)
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
			xticks=collect(year0:4:year1),ylabel=zlb)
		hm1=lines!(ax1,gl1.x,y)
		xlims!(ax1,(year0,year1))
		ylims!(ax1,rng)
		fig1
	end

	function DepthTime(x,y,z,levs,ttl,RC1,RC0,year0,year1; ClipToRange=true)
		ClipToRange ? to_range!(z,levs) : nothing
		fig1 = Figure(resolution = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1], title=ttl,
			xticks=collect(year0:4:year1))
		hm1=contourf!(ax1,x,y,z,levels=levs,colormap=:turbo)
		Colorbar(fig1[1,2], hm1, height = Relative(0.65))
		xlims!(ax1,year0,year1)
		ylims!(ax1,RC1,RC0)
		
		fig1
	end

	function TimeLat(x,y,z,levs,ttl,y0,y1,year0,year1; ClipToRange=true)
		ClipToRange ? to_range!(z,levs) : nothing
		fig1 = Figure(resolution = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1], title=ttl,
			xticks=collect(year0:4:year1),yticks=collect(-90.0:20.0:90.0),ylabel="latitude")
		hm1=contourf!(ax1,x,y,z,levels=levs,colormap=:turbo)
		Colorbar(fig1[1,2], hm1, height = Relative(0.65))
		xlims!(ax1,year0,year1)
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

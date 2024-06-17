
module ClimatologyMakieExt

	using Makie, Climatology
	import Climatology: plot_examples, load, mean, runmean

	function plot_examples(ID=Symbol,stuff...)
        if ID==:ECCO_map
			map(stuff...)
		elseif ID==:ECCO_TimeLat
			TimeLat(stuff...)
		elseif ID==:ECCO_TimeLatAnom
			TimeLat(stuff...)
		elseif ID==:ECCO_DepthTime
			DepthTime(stuff...)
		elseif ID==:ECCO_GlobalMean
			glo(stuff...)
		elseif ID==:ECCO_OHT
			OHT(stuff...)
		elseif ID==:ECCO_Overturn1
			figov1(stuff...)
		elseif ID==:ECCO_Overturn2
			figov2(stuff...)
		elseif ID==:ECCO_Transports
			transport(stuff...)
		else
			println("unknown plot ID")
		end
	end

	##

	to_range!(DD,levs::Tuple) = to_range!(DD,range(levs[1],levs[2],length=10))

	function to_range!(DD,levs)
		DD[findall(DD.<=levs[1])].=levs[1]+(levs[2]-levs[1])/100
		DD[findall(DD.>=levs[end])].=levs[end]-(levs[end]-levs[end-1])/100
	end

	years_to_display=(1979,2025)

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
		xlims!(ax,years_to_display)
	end

	function transport(namtrs,ncols,pth_out,list_trsp,year0,year1)
		if ncols > 1
			fig1 = Figure(size = (2000,1000),markersize=0.1)
		else
			fig1 = Figure(size = (900,400),markersize=0.1)
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

		fig1 = Figure(size = (900,400),markersize=0.1)
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
		xlims!(ax1,years_to_display)
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
	
		fig1 = Figure(size = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1], title="Meridional Overturning Streamfunction (in Sv, time mean)",
				xlabel="latitude",ylabel="depth (in m)")
		hm1=contourf!(ax1,x,y,z,levels=levs)
		Colorbar(fig1[1,2], hm1, height = Relative(0.65))
		fig1
	end

	function OHT(pth_out)
		fil=joinpath(pth_out,"MHT/MHT.jld2")
		tmp=load(fil,"single_stored_object")
		MT=vec(mean(tmp[:,1:240],dims=2))

		x=vec(-89.0:89.0)
		fig1 = Figure(size = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1], title="Northward Heat Transport (in PW, time mean)",
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

		fig1 = Figure(size = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1], title=ttl,
			xticks=collect(year0:4:year1),ylabel=zlb)
		hm1=lines!(ax1,gl1.x,y)
		xlims!(ax1,years_to_display)
		ylims!(ax1,rng)
		fig1
	end

	function DepthTime(x,y,z,levs,ttl,RC1,RC0,year0,year1; ClipToRange=true)
		ClipToRange ? to_range!(z,levs) : nothing
		fig1 = Figure(size = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1], title=ttl,
			xticks=collect(year0:4:year1))
		hm1=contourf!(ax1,x,y,z,levels=levs,colormap=:turbo)
		Colorbar(fig1[1,2], hm1, height = Relative(0.65))
		xlims!(ax1,years_to_display)
		ylims!(ax1,RC1,RC0)
		
		fig1
	end

	function TimeLat(x,y,z,levs,ttl,y0,y1,year0,year1; ClipToRange=true)
		ClipToRange ? to_range!(z,levs) : nothing
		fig1 = Figure(size = (900,400),markersize=0.1)
		ax1 = Axis(fig1[1,1], title=ttl,
			xticks=collect(year0:4:year1),yticks=collect(-90.0:20.0:90.0),ylabel="latitude")
		hm1=contourf!(ax1,x,y,z,levels=levs,colormap=:turbo)
		Colorbar(fig1[1,2], hm1, height = Relative(0.65))
		xlims!(ax1,years_to_display)
		ylims!(ax1,y0,y1)
		fig1
	end

	function map(λ,DD,levs,ttl; ClipToRange=true)
		ClipToRange ? to_range!(DD,levs) : nothing
		fig = Figure(size = (900,600), backgroundcolor = :grey95)
		ax = Axis(fig[1,1], title=ttl,xlabel="longitude",ylabel="latitude")
		hm1=contourf!(ax,λ.lon[:,1],λ.lat[1,:],DD,levels=levs,colormap=:turbo)
		Colorbar(fig[1,2], hm1, height = Relative(0.65))
		fig	
	end

end

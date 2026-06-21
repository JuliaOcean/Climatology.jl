module ERA5_plot

using Makie
using RollingFunctions, Statistics

import Climatology: read_bulk_formulae, SurfaceFluxDiag

#"spfh","tmp2m_degC","wspeed"
#"pres","rain","d2m",
#"u10m","v10m",
#"ustr","vstr",
function plot_bulk_formulae(df)
    fig=Figure(size=(600,900))
    lst=["dlw","dsw","hl","hs","qnet"]
    for v in 1:length(lst)
        vv=lst[v]
        ax=Axis(fig[v,1],title=vv)
        lines!(df[!,vv])
    end
    fig
end

function plot_bulk_formulae(fil::String)
	df=read_bulk_formulae(fil)
	plot_bulk_formulae(df)
end

#"spfh","tmp2m_degC","wspeed"
#"pres","rain","d2m",
#"u10m","v10m",
#"ustr","vstr",
function plot_surface_balance(df,tim,sst)
    fig=Figure(size=(1500,900),fontsize=24)
	ax=Axis(fig[1,1],title="temperature",xlabel="day since Jan. 1",ylabel="degree C")
	lines!(tim,df[!,"tmp2m_degC"],label="tmp2m_degC",linewidth=2)
	lines!(tim,sst,label="sst",linewidth=4,color=:red)
    axislegend(ax,position = :lt)
	ax=Axis(fig[1,2],title="radiative components",xlabel="day since Jan. 1",ylabel="W/m2")
    lst=["lw","sw","dlw","dsw","ulw","usw"]
	[lines!(tim,rnmn(df[!,vv],24),label=vv) for vv in lst]
    axislegend(ax,position = :lb, orientation = :horizontal); ylims!(ax,(-300,500))
	ax=Axis(fig[2,2],title="Qnet & radiative components",xlabel="day since Jan. 1",ylabel="W/m2")
    lst=["lw","sw","qnet"]
	[lines!(tim,rnmn(df[!,vv],24),label=vv) for vv in lst]
    axislegend(ax,position = :lb, orientation = :horizontal); ylims!(ax,(-500,300))
	ax=Axis(fig[2,1],title="Qnet & turbulent components",xlabel="day since Jan. 1",ylabel="W/m2")
    lst=["hl","hs","qnet"]
	[lines!(tim,rnmn(df[!,vv],24),label=vv) for vv in lst]
    axislegend(ax,position = :lt, orientation = :horizontal); ylims!(ax,(-400,400))
	fig
end

function plot_Qnet_cumsum(df,tim,sst)
    fig=Figure(size=(500,300),fontsize=11)
	ax=Axis(fig[1,1],title="cumulated(Qnet')",xlabel="day since Jan. 1",ylabel="non-dimensional")
	z=rnmn(df[!,"qnet"],24)
	z=cumsum(z .-mean(z))
	z=z./sqrt(mean(z.^2))
    lines!(tim,z,linewidth=4)
	fig
end

function rnmn(hourly,n=24)
	result = rollmean(hourly,n)
	result = [fill(result[1],Int(n/2))
	0.5*(result[1:end-1]+result[2:end])
	fill(result[end],Int(n/2))]
end

end

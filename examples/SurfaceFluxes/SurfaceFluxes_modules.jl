
module ERA5_read

using MeshArrays, NCDatasets, AirSeaFluxes, DataFrames, Interpolations
import CSV, Statistics

function ij(lon0,lat0,lon,lat)
    tmp1=(lon.-lon0).^2 .+ (lat.-lat0)'.^2
    (ii,jj)=findall(tmp1.==minimum(tmp1))[1].I
end


##

Rdry=287.0597 ; Rvap=461.5250 ; a1=611.21 ; a3=17.502 ; a4=32.19 ; T0=273.16
#Calculation of E saturation water vapour from Teten's formula
E(dtas)=a1*exp(a3*(dtas-T0)/(dtas-a4))
#Calculation of saturation specific humidity at 2m qsat  (equal to huss)
qsat(ps,E)=(Rdry/Rvap)*E/(ps-((1-Rdry/Rvap)*E))

wspeed(u10,v10)=sqrt(u10^2+v10^2)

##

function read_from_nc(fil::String,ii,jj)

list_in=["dlw","dsw","pres","rain","d2m","tmp2m_degC","u10m","ustr","v10m","vstr"]#,"wspeed"];
list_ds=["msdwlwrf","msdwswrf","sp","tp","d2m","t2m","u10","metss","v10","mntss"]#,"..."]

offset=zeros(12)
offset[6]=-273.15
factor=ones(12)
factor[1]=-1.0
factor[2]=-1.0
factor[4]=1/3600
factor[7]=1.0
factor[8]=-1.0
factor[9]=1.0
factor[10]=-1.0

df=DataFrame()
for vv in 1:length(list_ds)
    v_ds=list_ds[vv]
    tmp=Dataset(fil)[v_ds][ii,jj,:]
    v_in=list_in[vv]
    df[!,v_in]=offset[vv].+factor[vv]*tmp
end


df.spfh=[qsat(df.pres[i],E(df.d2m[i])) for i in eachindex(df.pres)]
df.wspeed=wspeed.(df.u10m,df.v10m)

df

end

##

"""
    read_one_year(year,ii,jj)

```
year0=2023
lon0=205; lat0=45;
(ii,jj)=ij(lon0,lat0,lon,lat)

df=read_one_year(year0,ii,jj)

fil="ERA5_lon"*string(lon0)*"_lat"*string(lat0)*"_year"*string(year0)*".csv"
CSV.write(joinpath(tempdir(),fil),df)
```
"""
function read_one_year(year,ii,jj)
  path_to_data="/Volumes/My Passport/data/ERA5"
  df=DataFrame()
  for m in 1:12
    mm=(m>9 ? "" : "0")*string(m)
    fil=joinpath(path_to_data,string(year),"ERA5_$(year)_$(mm).nc")
    append!(df,read_all(fil,ii,jj))
  end
  df
end

##

stefanBoltzmann = 5.670e-8 #[J*K^-4*m^-2*s^-1]
albedo=0.06

upsw(dsw)=albedo*abs(dsw)
uplw(sst)=stefanBoltzmann*(sst+273.15)^4

function interpolate_sst(sst,tim)
	xs = 0.5:364.5
	interp_linear = linear_interpolation(xs, sst, extrapolation_bc=Line())
	interp_linear(tim)
end

function surface_balance(df,sst)
	fluxes=bulkformulae.(df.tmp2m_degC.+273.16,df.spfh,df.wspeed,sst)
	df.hl=[x.hl for x in fluxes]
	df.hs=[x.hs for x in fluxes]
	df.evap=[x.evap for x in fluxes]

	stefanBoltzmann = 5.670e-8 #[J*K^-4*m^-2*s^-1]
	albedo=0.06

	upsw(dsw)=albedo*abs(dsw)
	uplw(sst)=stefanBoltzmann*(sst+273.15)^4

	df.dlw=abs.(df.dlw);
	df.dsw=abs.(df.dsw);
	df.ulw=uplw.(sst);
	df.usw=upsw.(df.dsw);

	df.lw=df.dlw-df.ulw;
	df.sw=df.dsw-df.usw;

	df.qnet=df.hl+df.hs+df.lw+df.sw

	df
end

end

##

module ERA5_plot

import ..ERA5_read: read_from_nc, ij, NCDatasets, bulkformulae
import .NCDatasets: Dataset

using CairoMakie
using RollingFunctions, Statistics

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

function plot_bulk_formulae(fil::String,sst0=15)
	lon=Dataset(fil)["longitude"][:]
	lat=Dataset(fil)["latitude"][:]

	lon0=205; lat0=45;
	(ii,jj)=ij(lon0,lat0,lon,lat)

	#t2m_mean=mean(Dataset(fil)["t2m"],dims=3)
	t2m=Dataset(fil)["t2m"][ii,jj,:]

	df=read_from_nc(fil,ii,jj)
	sst=fill(15.0,length(df.dlw))

	fluxes=bulkformulae.(df.tmp2m_degC.+273.16,df.spfh,df.wspeed,sst)
	df.hl=[x.hl for x in fluxes]
	df.hs=[x.hs for x in fluxes]
	df.evap=[x.evap for x in fluxes]
	df.qnet=df.hl+df.hs-df.dlw-df.dsw

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

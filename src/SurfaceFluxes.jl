
module ERA5

using MeshArrays, DataFrames
import CSV, Statistics
import Climatology: read_Dataset
import Climatology.downloads: ERA5_OISST_download

##

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

"""
	read_lonlat(; path_to_data="ERA5_data")

```
using NCDatasets, Climatology
ERA5.read_lonlat()
```
"""
function read_lonlat(; path_to_data="ERA5_data")
	fil=joinpath(path_to_data,"2023/ERA5_2023_01.nc")
	lon=read_Dataset(fil)["longitude"][:]
	lat=read_Dataset(fil)["latitude"][:]
	lon,lat
end
##

"""
    read_from_nc(fil::String, ii, jj)

Read ERA5 reanalysis variables at grid point `(ii,jj)` from NetCDF file
`fil`, apply unit conversions, and derive specific humidity and wind
speed.

Reads 10 ERA5 dataset variables (`list_ds`, ECMWF short names) and
renames/rescales them to `list_in` (this package's working names), via
parallel `offset`/`factor` lookup tables indexed by position:

| `list_in`      | `list_ds`  | offset  | factor    | notes                          |
|:----------------|:-------------|:----------|:------------|:---------------------------------|
| `dlw`            | `msdwlwrf`    | 0         | `-1.0`        | sign flip                        |
| `dsw`            | `msdwswrf`    | 0         | `-1.0`        | sign flip                        |
| `pres`           | `sp`          | 0         | `1.0`         |                                   |
| `rain`           | `tp`          | 0         | `1/3600`      | per-hour rate from per-hour total|
| `d2m`            | `d2m`         | 0         | `1.0`         | dewpoint, Kelvin                 |
| `tmp2m_degC`     | `t2m`         | `-273.15` | `1.0`         | Kelvin → Celsius                 |
| `u10m`           | `u10`         | 0         | `1.0`         |                                   |
| `ustr`           | `metss`       | 0         | `-1.0`        | sign flip                        |
| `v10m`           | `v10`         | 0         | `1.0`         |                                   |
| `vstr`           | `mntss`       | 0         | `-1.0`        | sign flip                        |

then derives:
- `spfh`: saturation specific humidity at 2 m, via `qsat(pres, E(d2m))`
  (Teten's formula, see `E`/`qsat`);
- `wspeed`: 10 m wind speed, via `wspeed.(u10m, v10m)`.

Returns a `DataFrame` with one row per time record in `fil`, columns
`list_in` plus `spfh`, `wspeed`.

!!! note
    `offset`/`factor` are sized `12`, but only the first `10` entries
    (matching `list_in`/`list_ds`'s length) are ever indexed. Worth
    double-checking this isn't hiding an intended 11th/12th variable
    dropped from `list_in`/`list_ds` without trimming the tables.
"""
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
    tmp=read_Dataset(fil)[v_ds][ii,jj,:]
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
import Climatology.ERA5: read_lonlat, read_one_year, ij

lon,lat=read_lonlat()

year0=2023
lon0=205; lat0=45;
(ii,jj)=ij(lon0,lat0,lon,lat)

df=read_one_year(year0,ii,jj)

fil="ERA5_lon"*string(lon0)*"_lat"*string(lat0)*"_year"*string(year0)*".csv"
CSV.write(joinpath(tempdir(),fil),df)

using CairoMakie
da=SurfaceFluxDiag((default=true,),df)
plot(da)
```
"""
function read_one_year(year,ii,jj; path_to_data="ERA5_data")
  df=DataFrame()
  for m in 1:12
    mm=(m>9 ? "" : "0")*string(m)
    fil=joinpath(path_to_data,string(year),"ERA5_$(year)_$(mm).nc")
    append!(df,read_from_nc(fil,ii,jj))
  end
  df
end

function read_bulk_formulae(fil::String)
	lon=read_Dataset(fil)["longitude"][:]
	lat=read_Dataset(fil)["latitude"][:]

	lon0=205; lat0=45;
	(ii,jj)=ij(lon0,lat0,lon,lat)

	#t2m_mean=mean(Dataset(fil)["t2m"],dims=3)
	t2m=read_Dataset(fil)["t2m"][ii,jj,:]

	df=read_from_nc(fil,ii,jj)
	sst=fill(15.0,length(df.dlw))

	fluxes=bulkformulae.(df.tmp2m_degC.+273.16,df.spfh,df.wspeed,sst)
	df.hl=[x.hl for x in fluxes]
	df.hs=[x.hs for x in fluxes]
	df.evap=[x.evap for x in fluxes]
	df.qnet=df.hl+df.hs-df.dlw-df.dsw

	df
end

"""
	read_sample(path=ERA5_OISST_download())

"""
function read_sample(path=ERA5_OISST_download())
	fil_ERA5="ERA5_lon205_lat45_year2023.csv"
	fil_OISST="OISST_lon205_lat45_year2023.csv"
	ERA5=CSV.read(joinpath(path,fil_ERA5),DataFrames.DataFrame)
	OISST=CSV.read(joinpath(path,fil_OISST),DataFrames.DataFrame)
	tim=(1:365*24)./24 .-0.5/24
	(tim=tim,OISST=OISST,ERA5=ERA5)
end

##

using AirSeaFluxes, Interpolations

stefanBoltzmann = 5.670e-8 #[J*K^-4*m^-2*s^-1]
albedo=0.06

upsw(dsw)=albedo*abs(dsw)
uplw(sst)=stefanBoltzmann*(sst+273.15)^4

function interpolate_sst(sst,tim)
	xs = 0.5:364.5
	interp_linear = linear_interpolation(xs, sst, extrapolation_bc=Line())
	interp_linear(tim)
end

"""
    surface_balance(df, sst)

Compute the full air-sea surface heat budget for `df` given sea surface
temperature(s) `sst`, mutating `df` in place with the added flux columns
and returning it.

Requires `df` to already have `tmp2m_degC`, `spfh`, `wspeed` (e.g. from
[`read_from_nc`](@ref)) and `dlw`, `dsw` (downward long-/shortwave, sign
convention: made positive-down via `abs.(...)` here regardless of input
sign).

Computes, via `AirSeaFluxes.bulkformulae` (turbulent fluxes) plus
Stefan-Boltzmann/albedo (radiative fluxes):
- `hl`, `hs`, `evap`: latent heat flux, sensible heat flux, evaporation,
  from `bulkformulae(tmp2m_degC + 273.16, spfh, wspeed, sst)` (temperature
  converted to Kelvin);
- `ulw = σ(sst+273.15)^4`: upward longwave (Stefan-Boltzmann,
  `σ = 5.670e-8` W·m⁻²·K⁻⁴);
- `usw = albedo * |dsw|`: upward (reflected) shortwave, `albedo = 0.06`;
- `lw = dlw - ulw`, `sw = dsw - usw`: net longwave/shortwave;
- `qnet = hl + hs + lw + sw`: total net surface heat flux into the
  ocean.

!!! note
    Locally redefines `upsw`/`uplw` closures identical to the
    module-level `upsw`/`uplw` functions defined earlier in this file —
    harmless shadowing, but redundant.
"""
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

SST_demo_path=joinpath(tempdir(),"demo_OISST")

##

module SST_FILES

using Printf, DataFrames, CSV, Dates, Glob
import Climatology: read_Dataset, SST_demo_path

read_files_list(;path=SST_demo_path,file="oisst_whole_file_list.csv",add_ymd=true) = begin
    if add_ymd
        add_to_table(CSV.read(joinpath(path,file),DataFrame))
    else
        CSV.read(joinpath(path,file),DataFrame)
    end
end

function add_to_table(list)
    ymd!(list)
    list.t=collect(1:length(list.day))
    list
end


"""
    file_lists(path="")

Create file lists and output to csv.

- `whole_file_list.csv` : all files through today's date
- `to_get_file_list.csv` : files that remain to download

Sample file names :

```
url="https://www.ncei.noaa.gov/thredds/dodsC/OisstBase/NetCDF/V2.1/AVHRR/198201/oisst-avhrr-v02r01.19820101.nc"
url="https://www.ncei.noaa.gov/thredds/fileServer/OisstBase/NetCDF/V2.1/AVHRR/198201/oisst-avhrr-v02r01.19820101.nc"
```
"""
function file_lists(;path=tempname())
    #url0="https://www.ncei.noaa.gov/thredds/fileServer/OisstBase/NetCDF/V2.1/AVHRR/"
    url0="https://noaa-cdr-sea-surface-temp-optimum-interpolation-pds.s3.amazonaws.com/data/v2.1/avhrr/"

    !ispath(path) ? mkdir(path) : nothing
    
    ndays=( today()-Date(1982,1,1) ).value
    file_list=DataFrame(fil=String[],url=String[],todo=Bool[])
    for t in 1:ndays
        dd=Date(1982,1,1)+Dates.Day(t-1)
        y=year(dd)
        m=month(dd)
        d=day(dd)
        url=@sprintf "%s%04i%02i%s%04i%02i%02i.nc" url0 y m "/oisst-avhrr-v02r01." y m d
        fil=@sprintf "%s/%04i%02i%s%04i%02i%02i.nc" path y m "/oisst-avhrr-v02r01." y m d
        push!(file_list,(fil=fil,url=url,todo=!isfile(fil)))
    end

    fil1=joinpath(path,"oisst_whole_file_list.csv")
    CSV.write(fil1,file_list)
    fil2=joinpath(path,"oisst_to_get_file_list.csv")
    CSV.write(fil2,file_list[file_list.todo,:])
    
    return fil1,fil2    
end

"""
    SST_FILES.ersst_file_lists(; path=SST_demo_path, dir=joinpath(path,"files_ersst"), last_month=nothing)

Build ERSST (Extended Reconstructed SST, monthly, 1854–present) file
lists, mirroring [`file_lists`](@ref) (which targets daily OISST files
instead): writes `"ersst_whole_file_list.csv"` (all months through
`last_month`) and `"ersst_to_get_file_list.csv"` (months not yet
downloaded) to `path`. The listed local file names are `dir/ersst.v5.YYYYMM.nc`.

ERSST is monthly and a month is published only after it ends (there are no
`_preliminary` files, unlike OISST). By default `last_month` is the previous
calendar month; pass any `Date` in the desired last month to override.
If that month is not yet online, [`SST_processing.download_ersst_files`](@ref)
skips it and it stays in the to-get list.

Source URL pattern:
`https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/netcdf/ersst.v5.YYYYMM.nc`.

Returns `(fil1, fil2)`, the paths to the whole-list and to-get-list CSVs.
"""
function ersst_file_lists(;path=SST_demo_path,dir=joinpath(path,"files_ersst"),last_month=nothing)
    url0="https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/netcdf/"

    # ERSST is monthly and published after month end (no preliminary files):
    # by default list through the previous calendar month.
    last_month = isnothing(last_month) ? firstdayofmonth(today())-Dates.Month(1) : firstdayofmonth(last_month)

    file_list=DataFrame(fil=String[],url=String[],todo=Bool[])
    for dd in Date(1854,1,1):Dates.Month(1):last_month
        y=year(dd)
        m=month(dd)
        url=@sprintf "%s%s%04i%02i.nc" url0 "ersst.v5." y m
        fil=joinpath(dir,@sprintf("ersst.v5.%04i%02i.nc",y,m))
        push!(file_list,(fil=fil,url=url,todo=!isfile(fil)))
    end

    fil1=joinpath(path,"ersst_whole_file_list.csv")
    CSV.write(fil1,file_list)
    fil2=joinpath(path,"ersst_to_get_file_list.csv")
    CSV.write(fil2,file_list[file_list.todo,:])
    
    return fil1,fil2    
end

"""
    test_files(list,ii=[])

Test whether all downloaded files are valid.

```
list=CSV.read("oisst_whole_file_list.csv",DataFrame)
list_pb=sst_files.test_files(list)
[Downloads.download(r.url,r.fil) for r in eachrow(list[list_pb,:])]
```
"""
function test_files(list,ii=[]; print_fails=false)
    test=zeros(1,length(list.fil))
    isempty(ii) ? jj=collect(1:length(list.fil)) : jj=ii
    for f in jj
       try
        ds=read_Dataset(list.fil[f])
        close(ds)
       catch e
        print_fails ? println(basename(list.fil[f])) : nothing
        test[f]=1
       end
    end
    return [i[2] for i in findall(test.==1)]
end

function ymd(f)
	tmp=split(f,".")[end-1]
	parse.(Int,[tmp[1:4] tmp[5:6] tmp[7:8]])
end

function ymd!(d::DataFrame)
	tmp=ymd.(d.fil)
	d[!, :year]=[a[1] for a in tmp]
	d[!, :month]=[a[2] for a in tmp]
	d[!, :day]=[a[3] for a in tmp]
	d
end

"""
    SST_FILES.monthlymean(gdf, m; path0=pwd(), varname="sst")

Compute the mean of `varname` across all files in group `m` of grouped
file list `gdf` (e.g. grouped by calendar month), reading each file
relative to `path0`.

Used by `SST_processing.monthly_climatology` to average all years' daily
files for a given calendar month into one climatological monthly field.
"""
function monthlymean(gdf,m;path0=pwd(),varname="sst")
    list=joinpath.(path0,gdf[m].fil)
    ds=read_Dataset(list[1])
    tmp=0*ds[varname][:,:,1,1]
    [tmp.+=read_Dataset(f)[varname][:,:,1,1] for f in list]
    tmp./length(list)
end

###

read_lon_lat(fil) = begin
    lon=read_Dataset(fil)["lon"][:]
    lat=read_Dataset(fil)["lat"][:]
    lon,lat
end

###

"""
    SST_FILES.read_map(; variable="anom", file="", file_climatology="")

Read a single day's OISST field from `file` (falling back to
`file[1:end-3]*"_preliminary.nc"` if `file` itself doesn't exist — OISST
publishes near-real-time files under a `_preliminary` suffix before the
finalized file is available).

`variable` selects what to return:
- `"sst"`: the raw SST field;
- `"anom"` (default): the precomputed anomaly field as stored in `file`;
- `"anom_recompute"`: SST minus the corresponding calendar month's field
  from `file_climatology` (read at `file`'s own month, `mon_sst`, derived
  from `file`'s date-stamped name) — i.e. an anomaly recomputed against a
  caller-supplied climatology rather than using the file's own stored
  anomaly.

Returns the selected 2D field.
"""
function read_map(;variable="anom",file="",file_climatology="")
	(year_sst,mon_sst,day_sst)=ymd(file)	
    isfile(file) ? fil_sst1=file : fil_sst1=file[1:end-3]*"_preliminary.nc"

    ds= read_Dataset(fil_sst1)
    sst=ds["sst"][:,:,1,1]
    anom = ds["anom"][:,:,1,1]
    close(ds)

    x = if variable=="anom_recompute"
    	sst_clim = read_Dataset(file_climatology)["sst"][:,:,mon_sst]
        sst-sst_clim
    elseif variable=="anom"
        anom
    else
        sst
    end

	x
end

"""
    SST_FILES.map_data(; variable="anom", file="", file_climatology="")

Figure-free accessor for an OISST map: reads one daily field (via
[`read_map`](@ref)) together with its coordinates and default display
settings, without constructing a Figure.

`variable` is one of `"sst"`, `"anom"` (NCEI anomaly stored in `file`) or
`"anom_recompute"` (anomaly against `file_climatology`, which is then
required). Any other value throws an error (unlike `read_map`, which treats
unknown values as `"sst"`).

Returns a `NamedTuple` with:

- `field`: 2D array (`lon` × `lat`, `Union{Missing,Float32}`, `missing` over land/ice)
- `lon`, `lat`: coordinate vectors (regular 0.25° grid, `lon` in `[0,360)`)
- `colorrange`, `colormap`: default display settings for `variable`
- `title`: default title (same wording as the SST notebook)
- `metadata`: `(variable, date, file, units, source, climatology)`

`plot(SSTdiag(path, title, :map; map_data=...))` additionally expects
`lon1`, `lat1` (marked point) and `showgrid`; add them with `merge`:

```julia
md = SST_FILES.map_data(variable="sst", file=fil)
plot(SSTdiag("", md.title, :map; map_data=merge(md, (lon1=205.0, lat1=45.0, showgrid=false))))
```
"""
function map_data(;variable="anom",file="",file_climatology="")
    variable in ("sst","anom","anom_recompute") ||
        error("unknown variable \"$(variable)\" (use \"sst\", \"anom\" or \"anom_recompute\")")
    variable=="anom_recompute" && isempty(file_climatology) &&
        error("variable=\"anom_recompute\" requires file_climatology")

    (y,m,d)=ymd(file)
    ymd_str=@sprintf "%04i%02i%02i" y m d
    fil=(isfile(file) ? file : file[1:end-3]*"_preliminary.nc")

    field=read_map(variable=variable,file=file,file_climatology=file_climatology)
    lon,lat=read_lon_lat(fil)

    (colorrange,colormap,title)=if variable=="anom"
        (4 .*(-1.0,1.0),:curl,"SST anomaly (NCEI) for time "*ymd_str)
    elseif variable=="anom_recompute"
        (4 .*(-1.0,1.0),:curl,"SST anomaly (recomputed) for $(y) / $(m) / $(d)")
    else
        ((-1.0,31.0),:thermal,"SST for time "*ymd_str)
    end

    metadata=(variable=variable,date=Date(y,m,d),file=fil,units="degC",
        source="NOAA OISST v2.1 AVHRR daily, 0.25 degree",
        climatology=(variable=="anom_recompute" ? file_climatology : nothing))

    (field=field,lon=lon,lat=lat,colorrange=colorrange,colormap=colormap,title=title,metadata=metadata)
end

"""
    SST_FILES.local_days(path)

Read-only scan of the daily OISST files that exist under `path`, laid out as
`path/YYYYMM/oisst-avhrr-v02r01.YYYYMMDD[_preliminary].nc`. Nothing is
downloaded or written, and names that do not match (other folders, stray
files) are ignored. Errors if `path` does not exist.

Returns a `NamedTuple` with one entry per day on disk, sorted by date:

- `dates`: `Vector{Date}`
- `files`: the *nominal final* file name of each day (`path/YYYYMM/oisst-avhrr-v02r01.YYYYMMDD.nc`),
  which may not exist; [`map_data`](@ref) falls back to `_preliminary` by itself
- `actual`: a file that does exist (the final one if both exist)
- `preliminary`: `true` where only a `_preliminary` file exists

```julia
d = SST_FILES.local_days("/path/to/OISST/files")
SST_FILES.map_data(variable="sst", file=d.files[end])
```
"""
function local_days(path)
    final=Dict{Date,String}(); prelim=Dict{Date,String}()
    for d in readdir(path)
        occursin(r"^\d{6}$",d) || continue
        dd=joinpath(path,d)
        isdir(dd) || continue
        for f in readdir(dd)
            m=match(r"^oisst-avhrr-v02r01\.(\d{8})(_preliminary)?\.nc$",f)
            m===nothing && continue
            t=Date(m.captures[1],dateformat"yyyymmdd")
            (m.captures[2]===nothing ? final : prelim)[t]=joinpath(dd,f)
        end
    end
    dates=sort!(collect(union(keys(final),keys(prelim))))
    nominal=[joinpath(path,Dates.format(t,"yyyymm"),
             "oisst-avhrr-v02r01.$(Dates.format(t,"yyyymmdd")).nc") for t in dates]
    actual=[haskey(final,t) ? final[t] : prelim[t] for t in dates]
    (dates=dates,files=nominal,actual=actual,preliminary=[!haskey(final,t) for t in dates])
end

"""
    SST_FILES.ersst_map_data(t::Date; path)

Figure-free data for the ERSST v5 monthly `ssta` map of the month of `t`.
Reads `path/ersst.v5.YYYYMM.nc` (see [`SST_FILES.ersst_file_lists`](@ref)).

If that month is not on disk, the newest local month at or before it is used
(ERSST months appear only after they end, so a month newer than the newest file
maps to the newest file). `metadata.fallback` is then `true` and `title` says
so. A month older than the first file on disk maps to the first file.
Nothing is downloaded or written. Errors if `path` has no `ersst.v5.YYYYMM.nc`
files.

Returns a `NamedTuple` with `field` (`lon` × `lat`, `Float64`, `NaN` where
missing), `lon`, `lat`, `colorrange`, `colormap`, `title` and `metadata`
(`requested`, `file`, `fallback`, `units`, `source`, `variable`).

```julia
md = SST_FILES.ersst_map_data(Date(2020,6,15); path="/path/to/OISST/files_ersst")
plot(SSTdiag(path="", name=md.title, options=(plot_type=:ersst_map, map_data=md)))
```
"""
function ersst_map_data(t::Date;path)
    want="ersst.v5.$(Dates.format(t,"yyyymm")).nc"
    avail=sort(filter(f->occursin(r"^ersst\.v5\.\d{6}\.nc$",f),readdir(path)))
    isempty(avail) && error("No ersst.v5.YYYYMM.nc files in $(path)")
    fil=avail[max(searchsortedlast(avail,want),1)]    # newest file at or before `want`, else the first one

    ds=read_Dataset(joinpath(path,fil))
    lon=ds["lon"][:]; lat=ds["lat"][:]
    a=ds["ssta"][:,:,1,1]
    close(ds)
    field=[ismissing(v) ? NaN : Float64(v) for v in a]

    (field=field,lon=lon,lat=lat,colorrange=4 .*(-1.0,1.0),colormap=:curl,
     title=(fil==want ? fil : "$(fil) (requested $(want) not on disk)"),
     metadata=(requested=t,file=joinpath(path,fil),fallback=(fil!=want),units="degC",
               source="NOAA ERSST v5 monthly",variable="ssta"))
end

CG = parentmodule(@__MODULE__).SST_coarse_grain   # sibling module, looked up at run time

"""
    SST_FILES.resolve_request(request; path, fil="lowres_oisst_sst_10.0.csv")

Resolve `request = (product=:OISST,)` into the data bundle the OISST explorer
needs, reading only: the daily file list under `path` (see
[`SST_FILES.local_days`](@ref)), the coarse-grained table `fil` (see
[`SST_coarse_grain.lowres_read`](@ref)), the file list CSV
(see [`SST_FILES.read_files_list`](@ref)) and the coarse grid of the newest
day (see [`SST_coarse_grain.grid`](@ref)). Nothing is downloaded or written.
Errors if `path` has no OISST daily files, or if `fil` is missing.

The product is `:OISST` only. The ERSST month lookup is `SST_FILES.ersst_map_data`
and `SST_FILES.resolve_request` is not needed for it.

Returns a `NamedTuple`:

- `path`: the folder given
- `days`: the `local_days` result (`dates`, `files`, `actual`, `preliminary`)
- `df`, `gdf`, `kdf`: the coarse table, its groups by `(i,j)`, and their keys
- `list`: the file list (`read_files_list`, read only)
- `G`: the coarse grid from the newest day on disk
- `info`: the `table_info(df; list, newest)` summary
- `metadata`: `product`, `source`, `path`, `table` (file name), `newest` (`Date`)

```julia
b = SST_FILES.resolve_request((product=:OISST,); path=oisst_files)
b.info.ncells
```
"""
function resolve_request(request; path, fil="lowres_oisst_sst_10.0.csv")
    request.product==:OISST || error("Only product :OISST is handled here; got $(repr(request.product)).")
    isdir(path) || error("path is not a directory: $(repr(path))")

    days=local_days(path)
    isempty(days.dates) && error("No oisst-avhrr-v02r01.YYYYMMDD.nc files under $(path)/YYYYMM/")
    isfile(joinpath(path,fil)) || error("Coarse table not found: $(joinpath(path,fil))")

    (df,gdf,kdf)=CG.lowres_read(path=path,fil=fil)
    list=read_files_list(path=path)
    newest=last(days.dates)
    info=CG.table_info(df; list=list, newest=newest)
    G=CG.grid(last(days.actual))

    (path=path,days=days,df=df,gdf=gdf,kdf=kdf,list=list,G=G,info=info,
     metadata=(product=:OISST,source="NOAA OISST daily",path=path,
               table=fil,newest=newest))
end

end


## 

module SST_coarse_grain

using Statistics, DataFrames, CSV, Glob, Dates
import Climatology: read_Dataset, SST_demo_path

"""
    SST_coarse_grain.areamean(arr, ii, jj, dnl)

Mean of `arr` over the `dnl × dnl` block of native-resolution cells
corresponding to coarse-grid cell `(ii,jj)`, skipping `missing` values.
"""
@inline areamean(arr,ii,jj,dnl) = 
    mean(skipmissing(
        arr[(ii-1)*dnl.+collect(1:dnl),(jj-1)*dnl.+collect(1:dnl)]
        ))


"""
    SST_coarse_grain.indices(list, dlon=10.0)

Determine which coarse-grid cells (at `dlon`-degree resolution) contain
valid (non-`NaN`) ocean data, using the first file in `list` as a
representative sample.

Returns `(i=ii[kk], j=jj[kk], k=kk)`: the coarse-grid `i`/`j` indices of
valid cells, and their linear index `k` into the full coarse grid — used
by [`calc_zm`](@ref) and `SST_processing.coarse_grain` to avoid computing
or storing land/all-`NaN` cells.
"""
function indices(list,dlon=10.0)
    dnl=Int(dlon/0.25)
    nnl=Int(720/dnl)

    fil=(isfile(list.fil[1]) ? list.fil[1] : list.fil[1][1:end-3]*"_preliminary.nc")
    println(fil)
    arr=read_Dataset(fil)["sst"][:,:]

    ii=[ii for ii in 1:nnl*2, jj in 1:nnl]
    jj=[jj for ii in 1:nnl*2, jj in 1:nnl]    
    tmp=[areamean(arr,ii,jj,dnl) for ii in 1:nnl*2, jj in 1:nnl]
    kk=findall((!isnan).(tmp))
    (i=ii[kk],j=jj[kk],k=kk)
end

"""
    grid(fil)

Return `(lon=lon,lat=lat,msk=msk,area=area)` based on `fil`.
"""
function grid(fil)
    fil=(isfile(fil) ? fil : fil[1:end-3]*"_preliminary.nc")

    ds=read_Dataset(fil)
    lon=ds["lon"][:]
    lat=ds["lat"][:]
    msk=ds["sst"][:,:]
    msk[ismissing.(msk)].=NaN
    msk=1 .+ 0*msk[:,:]
    area=[cellarea(lon0,lon0+0.25,lat0,lat0+0.25) for lon0 in 0:0.25:360-0.25, lat0 in -90:0.25:90-0.25]
    close(ds)
    (lon=lon,lat=lat,msk=msk,area=area)
end

"""
    cellarea(lon0,lon1,lat0,lat1)

[source](https://gis.stackexchange.com/questions/29734/how-to-calculate-area-of-1-x-1-degree-cells-in-a-raster)

As a consequence of a theorem of Archimedes, the area of a cell spanning longitudes l0 to l1 (l1 > l0) and latitudes f0 to f1 (f1 > f0) is

```(sin(f1) - sin(f0)) * (l1 - l0) * R^2```

where

- l0 and l1 are expressed in radians (not degrees or whatever).
- l1 - l0 is calculated modulo 2*pi (e.g., -179 - 181 = 2 degrees, not -362 degrees).
- R is the authalic Earth radius, almost exactly 6371 km.

!!! note
    As a quick check, the entire globe area can be computed by letting `l1 - l0 = 2pi`, `f1 = pi/2`, `f0 = -pi/2`. The result is `4 * Pi * R^2`.
"""
function cellarea(lon0,lon1,lat0,lat1)
    EarthRadius = 6371.0
    #f0=20; f1=21; l0=349; l1=350;
    f0=-90; f1=90; l0=0; l1=360;
    1e6 * (sind(lat1) - sind(lat0)) * mod1(deg2rad(lon1 - lon0),2pi) * EarthRadius^2
end

@inline nansum(x) = sum(filter(!isnan,x))
@inline nansum(x,y) = mapslices(nansum,x,dims=y)

"""
    SST_coarse_grain.areaintegral(arr, i::Int, j::Int, G::NamedTuple, dnl)

Area-weighted sum of `arr` over the `dnl × dnl` block of native-resolution
cells corresponding to coarse-grid cell `(i,j)`, weighted by the native
grid's mask `G.msk` and cell area `G.area`. Used by [`calc_zm`](@ref) to
build per-latitude-band area weights.
"""
@inline areaintegral(arr,i::Int,j::Int,G::NamedTuple,dnl) = begin
    ii=(i-1)*dnl.+collect(1:dnl)
    jj=(j-1)*dnl.+collect(1:dnl)
    nansum(arr[ii,jj].*G.msk[ii,jj].*G.area[ii,jj])
end

"""
    SST_coarse_grain.calc_zm(G::NamedTuple, df, dnl=missing)

Compute a zonal-mean (latitude-band) time series from a coarse-grained,
long-format SST table `df` (as produced by `SST_processing.coarse_grain`
+ `SST_coarse_grain.lowres_read`, with columns `i`, `j`, `t`, `sst` for
coarse-grid longitude/latitude indices, time index, and SST value).

`G` is the coarse grid `NamedTuple` from `grid` (`lon`, `lat`, `msk`,
`area`). `dnl` is the coarse-graining factor in grid cells (e.g.
`dlon/0.25` for a `dlon`-degree coarse cell); if not given, defaults to
the equivalent of `dlon = 10.0`.

For each coarse latitude-band index `k` (from `minimum(df.j)` to
`maximum(df.j)`), computes the area-weighted mean SST across all
longitude cells in that band, for every time step (grouped via
`groupby(df, :t)`). Returns an `(nlat, ntime)` array `arr`, with rows for
latitude bands outside `[minimum(df.j), maximum(df.j)]` left as `NaN`.

!!! note
    The local variable computed from `dnl` when it isn't given is
    currently unused — `dnl` itself (`missing`, in that case) is passed
    directly to `areaintegral` regardless. Worth checking the
    `dnl=missing` default path actually behaves as intended.
"""
function calc_zm(G::NamedTuple,df,dnl=missing)
    gdf_tim=groupby(df, :t)
    arr=NaN*zeros(maximum(df.j),length(gdf_tim))

    dn=if isempty(dnl)
        dlon=10.0
        Int(dlon/0.25)
    else
        dnl
    end

    for k in minimum(df.j):maximum(df.j)
        area_tmp=[areaintegral(G.msk,x.i,x.j,G,dnl) for x in eachrow(gdf_tim[1])]
        area_tmp[gdf_tim[1].j.!==k].=0
        tmp1=[sum(tmp1.sst[:].*area_tmp)/sum(area_tmp) for tmp1 in gdf_tim]
        arr[k,:].=tmp1
    end
    return arr
end

"""
    lowres_merge(;path=SST_demo_path,variable="sst")

Merge all files found in chosen path.
"""
function merge_files(;path=SST_demo_path,variable="sst",dlon=10.0)
    path0=dirname(file_root(path=path,variable=variable))
    file_list=glob("$(variable)_lowres*csv",path0)

    df=DataFrame(i=Int[],j=Int[],t=Int[],sst=Float32[])
    [lowres_append!(df,f) for f in file_list]
    CSV.write(joinpath(path,"lowres_oisst_$(variable)_$(dlon).csv"),df)
end

function lowres_append!(df,f)
    tmp=CSV.read(f,DataFrame)
    tmp.t.=parse(Int,split(basename(f),"_")[end][1:8])
    append!(df,tmp)
    return tmp
end

file_root(;path=SST_demo_path,variable="sst") = joinpath(path,"$(variable)_lowres_files","$(variable)_lowres_")

"""
    lowres_read(;path=SST_demo_path,fil="lowres_oisst_sst_10.0.csv")

Read `sst_lowres.csv`
"""
function lowres_read(;path=SST_demo_path,fil="lowres_oisst_sst_10.0.csv")
    fil=joinpath(path,fil)
    df=CSV.read(fil,DataFrame)
    gdf=groupby(df, [:i, :j])
    kdf=keys(gdf)
    return (df,gdf,kdf)
end

"""
    SST_coarse_grain.lowres_index(lon0, lat0, kdf)

Find the index into grouped-keys `kdf` (coarse-grid `(i,j)` pairs) whose
cell center is nearest to `(lon0,lat0)`.

See also [`lowres_position`](@ref) (the inverse: index → coordinates).
"""
function lowres_index(lon0,lat0,kdf)
    (i,j)=([x.i for x in kdf],[x.j for x in kdf])
	dx=Int(360/maximum(i))
    (ii,jj)=(dx*i.-dx/2,dx*j.-dx/2 .-90)
    d=(ii .-lon0).^2 .+ (jj .-lat0).^2
    findall(d.==minimum(d))[1]
end

"""
    SST_coarse_grain.lowres_position(ii, jj, kdf)

Convert coarse-grid indices `ii`,`jj` (as found in `kdf`) to their cell-
center `(longitude, latitude)` coordinates, given the coarse resolution
implied by `kdf`'s index range (`dx = 360/maximum(i)`).

See also [`lowres_index`](@ref) (the inverse: coordinates → index).
"""
lowres_position(ii,jj,kdf) = begin
    (i,j)=([x.i for x in kdf],[x.j for x in kdf])
	dx=Int(360/maximum(i))
	(dx*ii.-dx/2,dx*jj.-dx/2 .-90)
end

"""
    SST_coarse_grain.lowres_cell(lon0, lat0, gdf)

Find the coarse-grid cell of the grouped table `gdf` (from
[`lowres_read`](@ref)) that **contains** `(lon0,lat0)`.

Returns the `GroupKey` `(i,j)` of that cell, or `nothing` if the point lies
over land / in a cell absent from the table (cells with no valid ocean data
in the first file are never stored), or if `lat0` is outside `[-90,90]`.

Differs from [`lowres_index`](@ref), which returns the *nearest stored cell
center* (so a land point silently maps to a distant ocean cell) and does
not wrap longitude. Here `lon0` is wrapped modulo 360, so `-155` and `205`
are equivalent. Cell width `dx` is inferred as in `lowres_index`
(`dx = 360/maximum(i)`), i.e. 10° for the standard table.

The returned key is opaque and specific to this table; pass it to
[`lowres_series`](@ref) or `gdf[key]`.
"""
function lowres_cell(lon0,lat0,gdf)
    (-90 <= lat0 <= 90) || return nothing
    dx=Int(360/maximum(k.i for k in keys(gdf)))
    i=floor(Int,mod(lon0,360)/dx)+1
    j=min(floor(Int,(lat0+90)/dx)+1, Int(180/dx))
    haskey(gdf,(i=i,j=j)) ? keys(gdf)[findfirst(k->k.i==i&&k.j==j,keys(gdf))] : nothing
end

"""
    SST_coarse_grain.lowres_series(gdf, cell)

Daily time series of the coarse-grained cell `cell` (a `GroupKey` from
[`lowres_cell`](@ref)). Returns `(time, sst)`: `time` is a `Vector{Date}`
decoded from the table's `YYYYMMDD` integer column `t`, and `sst` the
area-mean sea surface temperature (degC) over the cell's ocean points
(`SST_coarse_grain.areamean` of the raw OISST `sst` variable, i.e. **not**
an anomaly).
"""
function lowres_series(gdf,cell)
    g=gdf[cell]
    time=[Date(t÷10000,(t÷100)%100,t%100) for t in g.t]
    (time=time,sst=g.sst)
end

"""
    SST_coarse_grain.point_series(lon0, lat0, gdf)

Daily coarse-grained SST series at `(lon0,lat0)`; the SST counterpart of
`ECCO_procs.ECCO_map_point` (figure-free). Returns `nothing` over land or
outside the table (see [`lowres_cell`](@ref)), otherwise a `NamedTuple`:

- `time`, `sst`: see [`lowres_series`](@ref) (daily, gap-free in the standard table)
- `cell`: opaque `GroupKey` identifying the coarse cell
- `lon`, `lat`: center of that cell, from [`lowres_position`](@ref)
- `requested`: `(lon=lon0, lat=lat0)` as given
- `time_axis`: `:multi_year` (contrast with ECCO's `:seasonal`)
- `metadata`: `(variable, units, dlon, source)`

!!! note
    The values are the mean over a whole `dlon`×`dlon` box (10°×10°), not
    the value at the requested point.
"""
function point_series(lon0,lat0,gdf)
    cell=lowres_cell(lon0,lat0,gdf)
    isnothing(cell) && return nothing
    s=lowres_series(gdf,cell)
    (lon,lat)=lowres_position(cell.i,cell.j,keys(gdf))
    dlon=Float64(360/maximum(k.i for k in keys(gdf)))
    (time=s.time,sst=s.sst,cell=cell,lon=lon,lat=lat,
     requested=(lon=lon0,lat=lat0),time_axis=:multi_year,
     metadata=(variable="sst",units="degC",dlon=dlon,
       source="NOAA OISST v2.1 AVHRR daily, coarse-grained (area mean of 0.25 degree cells)"))
end

"""
    SST_coarse_grain.table_info(df; list=nothing, newest=nothing)

Summary and consistency checks for a coarse-grained table `df` (the first
output of [`lowres_read`](@ref); long format with columns `i`, `j`, `t`,
`sst`, where `t` is an integer `YYYYMMDD`). Returns a `NamedTuple`:

- `t0`, `t1`: first and last date in the table
- `ndays`: number of distinct days
- `gap_free`: `true` if `ndays` equals the number of calendar days from `t0` to `t1`
- `ncells`: number of coarse cells `(i,j)`
- `list_rows`, `list_ok`: if a file `list` is given (e.g. from `SST_FILES.read_files_list`),
  its number of rows and whether it covers the table (`rows >= ndays`; see
  `SST_timeseries.calc`, which aligns the table to the list by position); else `nothing`
- `lag_days`: if `newest` (a `Date`, e.g. the last daily file on disk) is given,
  the number of days by which the table ends before it (`0` if up to date); else `nothing`

```julia
d = SST_FILES.local_days(path)
info = SST_coarse_grain.table_info(df; list=list, newest=last(d.dates))
info.lag_days > 0 && @warn "table is stale; rerun SST_processing.coarse_grain"
```
"""
function table_info(df;list=nothing,newest=nothing)
    d=Date.(string.(sort!(unique(df.t))),dateformat"yyyymmdd")
    expected=Dates.value(d[end]-d[1])+1
    (t0=d[1],t1=d[end],ndays=length(d),gap_free=length(d)==expected,
     ncells=nrow(unique(df[!,[:i,:j]])),
     list_rows=isnothing(list) ? nothing : size(list,1),
     list_ok=isnothing(list) ? nothing : size(list,1)>=length(d),
     lag_days=isnothing(newest) ? nothing : Dates.value(newest-d[end]))
end

end

##

module SST_processing

using Distributed, Dataverse, DataFrames
import Dataverse.downloads: Downloads
import Climatology: SST_FILES, SST_coarse_grain, read_Dataset
import Climatology: SST_demo_path, to_monthly_file, write_SST_climatology

"""
    SST_processing.download_ersst_files(; path=SST_demo_path, dir=joinpath(path,"files_ersst"),
                                          last_month=nothing, ntasks=8, verbose=false)

Download the monthly ERSST v5 NetCDF files that are not yet in `dir`
(created if needed). The file lists are regenerated first via
[`SST_FILES.ersst_file_lists`](@ref), which also defines `last_month` (default:
previous calendar month; ERSST has no preliminary files, so a month appears
online only after it ends).

Existing files are never re-downloaded. Each file is downloaded to a
`.part` temporary name and renamed on success, so an interrupted download
cannot leave a truncated `.nc` that later looks complete. Downloads run as
`ntasks` concurrent tasks in the current process (the files are small; no
`Distributed` workers are needed). Months that are not online are skipped
(reported if `verbose`) and remain in `ersst_to_get_file_list.csv`.

Returns the local paths of all months now available.
"""
function download_ersst_files(;path=SST_demo_path,dir=joinpath(path,"files_ersst"),
        last_month=nothing,ntasks=8,verbose=false)
    mkpath(dir)
    SST_FILES.ersst_file_lists(path=path,dir=dir,last_month=last_month)
    list=SST_FILES.CSV.read(joinpath(path,"ersst_whole_file_list.csv"),DataFrame)

    todo=list[list.todo.==true,:]
    verbose && println("$(nrow(todo)) of $(nrow(list)) months to download")

    asyncmap(eachrow(todo);ntasks=ntasks) do r
        tmp=r.fil*".part"
        try
            Downloads.download(r.url,tmp)
            mv(tmp,r.fil;force=true)
        catch
            rm(tmp;force=true)
            verbose && println("not available online : "*basename(r.fil))
        end
    end

    # refresh the to-get list so it reflects what is still missing
    list.todo=[!isfile(f) for f in list.fil]
    SST_FILES.CSV.write(joinpath(path,"ersst_to_get_file_list.csv"),list[list.todo,:])

    list.fil[.!list.todo]
end

"""
    SST_processing.download_files(; path=SST_demo_path, short_demo=false, verbose=false)

Download the OISST daily NetCDF files listed by `SST_FILES.file_lists`/
`SST_FILES.read_files_list`, distributing work across available Julia
workers.

If `path` doesn't exist, it's created. The file list is regenerated via
`SST_FILES.file_lists(path=path)`; when `short_demo` is `true`, only the
most recent 30 files are downloaded (for quick testing). Work is split
evenly across `nworkers()` via `@distributed`; each missing file is
downloaded via `Downloads.download`, falling back to a
`"_preliminary.nc"`-suffixed URL/filename if the primary download fails
(OISST publishes near-real-time files under a `_preliminary` suffix
before the finalized file is available), and silently skipping (with an
optional `verbose` message) if neither is found.

Returns the list of successfully-available local file paths (preferring
the finalized file over the preliminary one where both exist), excluding
any still-missing entries.
"""
function download_files(;path=SST_demo_path,short_demo=false,verbose=false)
    !ispath(path) ? mkdir(path) : nothing
    fil,_=SST_FILES.file_lists(path=path)
    list=SST_FILES.read_files_list(path=path)
    list=(short_demo ? list[end-29:end,:] : list)
    n_per_workwer=Int(ceil(length(list.fil)/nworkers()))

    if !isempty(list.fil)

    @sync @distributed for m in 1:nworkers()
        n0=n_per_workwer*(m-1)+1
        n1=min(n_per_workwer*m,length(list.fil))
        verbose ? println("$(n0),$(n1)") : nothing
        for r in eachrow(list[n0:n1,:])
            !isdir(dirname(r.fil)) ? mkdir(dirname(r.fil)) : nothing
            if !isfile(r.fil)
            verbose ? println(r.fil) : nothing
            try
                Downloads.download(r.url,r.fil)
            catch
                try
                    Downloads.download(r.url[1:end-3]*"_preliminary.nc",r.fil[1:end-3]*"_preliminary.nc")
                catch
                    verbose ? println("file not found online : "*r.fil[1:end-3]) : nothing
                end
            end
            end
        end
    end

    else

        verbose ? println("no more files to process") : nothing

    end

    nl=length(list.fil)
    tst=fill("",nl)
    for ll in 1:nl
        if isfile(list.fil[ll])
            tst[ll]=list.fil[ll]
        elseif isfile(list.fil[ll][1:end-3]*"_preliminary.nc")
            list.fil[ll][1:end-3]*"_preliminary.nc"
            tst[ll]=list.fil[ll][1:end-3]*"_preliminary.nc"
        else
            tst[ll]=""
        end
    end

    tst[findall((!isempty).(tst))]

end

## 

"""
    SST_processing.coarse_grain(; datname="oisst", varname="sst", dlon=10.0,
                                  path=SST_demo_path, short_demo=false)

Coarse-grain the downloaded OISST daily files (see
[`download_files`](@ref)) onto a `dlon`-degree grid, writing one CSV per
input file plus a single merged CSV of all coarse-grained values.

Reads the file list `"\$(datname)_whole_file_list.csv"` from `path` (when
`short_demo` is `true`, only the most recent 10 files); determines the
sparse set of non-empty coarse cells once via
`SST_coarse_grain.indices(list)` (reused for every file, since land/ocean
geography doesn't change over time). Distributes files evenly across
`nworkers()`: for each file, loads `varname` from the NetCDF (falling
back to the `"_preliminary.nc"` filename if the primary is missing),
computes the coarse-cell area means via `SST_coarse_grain.areamean`, and
writes the result to its own CSV under
`"\$(varname)_lowres_files/\$(varname)_lowres_<date>.csv"` (any pre-existing
output directory for `varname` is moved aside via `mv` to a temp path
before starting, rather than merged into).

After all files are processed, calls `SST_coarse_grain.merge_files` to
concatenate the per-file CSVs into a single
`"lowres_oisst_\$(varname)_\$(dlon).csv"`.
"""
function coarse_grain(;datname="oisst",varname="sst",dlon=10.0,
       path=SST_demo_path,short_demo=false)

    ## setup
    list=SST_FILES.read_files_list(file="$(datname)_whole_file_list.csv",path=path,add_ymd=false)
    list=(short_demo ? list[end-9:end,:] : list)

    ind=SST_coarse_grain.indices(list)
    nt=length(list.fil)
    n_per_workwer=Int(ceil(nt/nworkers()))

    file_root=SST_coarse_grain.file_root(variable=varname,path=path)
    isdir(dirname(file_root)) ? mv(dirname(file_root),tempname()) : nothing
    mkdir(dirname(file_root))

    ## distributed computation
    @sync @distributed for m in 1:nworkers()
        n0=n_per_workwer*(m-1)+1
        n1=min(n_per_workwer*m,length(list.fil))
        dnl=Int(dlon/0.25)
        nnl=Int(720/dnl)
        println("$(n0),$(n1)")
        for n in n0:n1
            r=list[n,:]
            fil=(isfile(r.fil) ? r.fil : r.fil[1:end-3]*"_preliminary.nc")
            if isfile(fil)
                #calculate
                ds=read_Dataset(fil)
                tmp=ds[varname][:,:]
                sst=[SST_coarse_grain.areamean(tmp,ii,jj,dnl) for ii in 1:nnl*2, jj in 1:nnl]
                #save to csv
                df=SST_FILES.DataFrame(i=ind.i,j=ind.j,sst=Float32.(sst[ind.k]))
                tmp=split(basename(r.fil),".")[2]
                SST_FILES.CSV.write(file_root*tmp*".csv",df)
            end
        end
    end

    ## write to final file
    SST_coarse_grain.merge_files(variable=varname,path=path,dlon=dlon)
end

##

"""
    SST_processing.monthly_climatology(; datname="oisst", varname="sst", path=SST_demo_path)

Compute the 1992–2011 monthly climatology (mean SST and mean anomaly, per
calendar month) from the downloaded OISST daily files, and write it to a
single climatology file via `write_SST_climatology`.

Reads the file list `"\$(datname)_whole_file_list.csv"` from `path`,
selects the 1992–2011 subset, and groups it by calendar month. For each
of `"sst"` and `"anom"`, computes the across-years mean for each of the
12 calendar months via `SST_FILES.monthlymean` (distributed across
`nworkers()`), writing each month's mean field via `to_monthly_file` to a
temporary output directory.

Finally combines the 12 per-month files into the single climatology file
via `write_SST_climatology(output_path, year0, year1, lon, lat)`
(`year0=1992`, `year1=2011`), and returns that file's path.

!!! note
    The function's own `varname` keyword (default `"sst"`) is shadowed
    by an internal loop variable of the same name iterating over
    `("sst","anom")` — the keyword argument itself has no effect on the
    computation; both `sst` and `anom` climatologies are always computed
    regardless of what's passed in.
"""
function monthly_climatology(;datname="oisst",varname="sst",path=SST_demo_path)
    year0=1992; year1=2011
    list=SST_FILES.read_files_list(file="$(datname)_whole_file_list.csv",path=path,add_ymd=true)
    lon,lat=SST_FILES.read_lon_lat(list.fil[1])

    sel=findall([(f.year>=year0 && f.year<=year1) for f in eachrow(list)])
    suf="$(year0)_$(year1)_"
    gdf=groupby(list[sel,:],:month)

    output_path=tempname(); mkdir(output_path)
    println("output path="*output_path)
    
    n_per_workwer=Int(ceil(12/nworkers()))
    n_per_workwer*nworkers()!==12 ? println("need nworkers to divide 12") : nothing

    for varname in ("sst","anom")
        @sync @distributed for m in 1:nworkers()
            for mm in 1:n_per_workwer
                month=(m-1)*n_per_workwer+mm
                tmp=SST_FILES.monthlymean(gdf,month,varname=varname)
                to_monthly_file(tmp,month,varname=varname,output_path=output_path)
            end
        end
    end

    output_file=write_SST_climatology(output_path,year0,year1,lon,lat)
end

end

##

module SST_timeseries

using DataFrames, Statistics, Dates
import Climatology: SST_coarse_grain

"""
    SST_timeseries.calc(input, list; title="", gdf=nothing)

Compute a full SST time-series diagnostic `NamedTuple` — raw values,
day-of-year climatology, anomaly, and extreme-warm quantile bands —
suitable for the Makie extension's `by_time`/`by_year`/`MHW`/
`local_and_global` plots (via `X.options.timeseries`).

`input` is either a raw SST vector, or a `DataFrames.GroupKey` into `gdf`
(a grouped `DataFrame`, e.g. grouped by grid cell `(i,j)`) — in the
latter case `gdf[input].sst` supplies the series. `list` is the full
file/date list (as from `SST_FILES.read_files_list`), used to align
`year`/`month`/`day` with each point in the series and to compute the
climatology.

Internally:
1. `repeatclim` computes the 1992–2011 day-of-year climatology (via
   `clim`/`gdf_clim`) and repeats it across the full series length,
   giving `clim`.
2. `anom` computes `sst - climatology`, re-centered to the climatology's
   own median (so `anom`'s scale matches `sst`, not a zero-centered
   anomaly).
3. `calc_quantile` computes, for each day of year, the 10th/90th
   percentile of the 1992–2011 anomaly (a ±2-day window around each
   calendar day) — returned as `low`/`high`, used to flag extreme
   warm/cool periods (e.g. the Makie extension's `MHW` plot).

Returns `(sst, clim, anom, title, year, month, day, low, high)`. `title`
defaults to `"SST time series"` unless overridden.
"""
function calc(input,list; title="", gdf=nothing)
	if isa(input,DataFrames.GroupKey)
		sst1=gdf[input].sst[:]
	else
		sst1=input[:]
	end
    nt=size(sst1,1)
    sst2=repeatclim(sst1,list[1:nt,:])
	sst3=anom(sst1,list[1:nt,:])

	ttl="SST time series"
	#isa(input,DataFrames.GroupKey) ? ttl=ttl*"for i="*string(input.i)*", j="*string(input.j) : nothing 
	!isempty(title) ?  ttl=title : nothing

    ts=(sst=sst1,clim=sst2,anom=sst3,title=ttl,
    year=list.year[1:nt],month=list.month[1:nt],day=list.day[1:nt])

    tmp1=calc_quantile(ts)

    merge(ts,tmp1)
end

function gdf_clim(list)
	sel=findall([(f.year>=1992 && f.year<=2011) for f in eachrow(list)])
	groupby(list[sel,:],[:month,:day])
end

@inline clim(sst,list) = [mean(sst[a.t[:]]) for a in gdf_clim(list)]

@inline function anom(sst,list)
    c=clim(sst,list)
    a=0*sst
    for t in 1:length(list.t)
        (y,m,d)=(list.year[t],list.month[t],list.day[t])
        tt=min(1+(Date(y,m,d)-Date(y,1,1)).value,365)
        a[t]=sst[t]-c[tt]
    end
    a.+median(c)
end

@inline function repeatclim(sst,list)
    c=clim(sst,list)
    a=0*sst
    for t in 1:length(list.t)
        (y,m,d)=(list.year[t],list.month[t],list.day[t])
        tt=min(1+(Date(y,m,d)-Date(y,1,1)).value,365)
        a[t]=c[tt]
    end
    a
end

##

@inline function calc_quantile(x,msk,yearday,yd)
	d0=yearday[yd]
	d1=[sum(mod1.( d0 .+ (-2:2),365) .==dd)==1 for dd in yearday]
	
	sel=findall(msk .&& d1)
	quantile(x[sel], [0.1, 0.9])
end

@inline function calc_quantile(ts)
    x=ts.sst-ts.clim
	msk=(ts.year.>=1992 .&& ts.year.<=2011)
	
	yearday=Date.(ts.year,ts.month,ts.day)-Date.(ts.year,1,1)
	yearday=min.(1 .+ [yd.value for yd in yearday],365)

	ts_low=zeros(365)
	ts_high=zeros(365)
	for yd in 1:365
		ts_low[yd],ts_high[yd]=calc_quantile(x,msk,yearday,yd)
	end
	
	(low=ts_low[yearday],high=ts_high[yearday])
end

"""
    SST_timeseries.global_mean(df, G, list; dnl=40, title="Global Mean SST")

Global-mean SST time series from the coarse-grained table `df` (see
`SST_coarse_grain.lowres_read`), passed through [`calc`](@ref) so the result
has the same fields (`sst`, `clim`, `anom`, `year`, `month`, `day`, `low`,
`high`, `title`) as a local series. `G` is `SST_coarse_grain.grid(file)`,
`list` the file list (`SST_FILES.read_files_list`), and `dnl` the number of
native 0.25 degree cells per coarse cell (`dlon/0.25`).

!!! note
    The area weights are computed once, from the rows of the **last day**,
    and applied to every day. This assumes the table has the same cell order
    on every day, as in `sst_anomaly_notebook.jl`. That assumption has not
    been checked against the table.
"""
function global_mean(df,G::NamedTuple,list;dnl=40,title="Global Mean SST")
    gdf1=groupby(df,:t)
    last_day=gdf1[end]
    area_tmp=[SST_coarse_grain.areaintegral(G.msk,x.i,x.j,G,dnl) for x in eachrow(last_day)]
    glmsst=[sum(g.sst[:].*area_tmp)/sum(area_tmp) for g in gdf1]
    calc(glmsst,list,title=title)
end

"""
    SST_timeseries.zonal_anomaly(G, df, dnl, list)

Zonal-mean SST anomaly as an `(nlat, ntime)` array: `SST_coarse_grain.calc_zm(G,df,dnl)`
with each latitude band (row) replaced by `sst - clim` from [`calc`](@ref)
(the 1992-2011 day-of-year climatology). Rows without data stay `NaN`.
This is the `zonal_mean` input of `plot(SSTdiag(...,:TimeLat))`.

!!! note
    Inherits the `dnl` behaviour of [`SST_coarse_grain.calc_zm`](@ref), which
    passes `dnl` through even when it is `missing`.
"""
function zonal_anomaly(G::NamedTuple,df,dnl,list)
    z=SST_coarse_grain.calc_zm(G,df,dnl)
    for l in 1:size(z,1)
        if !isnan(z[l,1])
            tmp=calc(z[l,:],list,title="Zonal Mean SST")
            z[l,:].=tmp.sst-tmp.clim
        end
    end
    z
end

end

##

module SST_scenarios
	function read_temp(fil)
	
		log=readlines(fil)
	
		ii=findall([occursin("tas=",i) for i in log])
		nt=length(ii)
		tas=zeros(nt)
		year=zeros(nt)
	
		for i in 1:nt
			tmp=split(log[ii[i]],"=")[2]
			tas[i]=parse(Float64,split(tmp,"degC")[1])
			year[i]=parse(Float64,split(tmp,"in")[2])
		end
	
		year,tas
	end
	
	function calc_offset(year_sst,ny,scenario=245)
		year1=year_sst+ny
		hector_fil="hector_scenarios/temperature_ssp$(scenario).log"
		hector_year,hector_tas=read_temp(hector_fil)
		y0=findall(hector_year.==year_sst)[1]
		y1=findall(hector_year.==year1)[1]
		hector_tas[y1]-hector_tas[y0]
	end
end

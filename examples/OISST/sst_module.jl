module sst_files

using Printf, DataFrames, CSV, Dates, NCDatasets, Glob

read_files_list(;path=tempdir(),file="oisst_whole_file_list.csv",add_ymd=true) = begin
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
    url0="https://www.ncei.noaa.gov/thredds/fileServer/OisstBase/NetCDF/V2.1/AVHRR/"

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


function ersst_file_lists(;path=tempdir())
    url0="https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/netcdf/"

    nmonths=(2023-1854)*12+7
    file_list=DataFrame(fil=String[],url=String[],todo=Bool[])
    for t in 1:nmonths
        dd=Date(1854,1,1)+Dates.Month(t-1)
        y=year(dd)
        m=month(dd)
        d=day(dd)
        url=@sprintf "%s%s%04i%02i.nc" url0 "ersst.v5." y m
        fil=@sprintf "files_ersst/ersst.v5.%04i%02i.nc" y m
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
function test_files(list,ii=[])
    test=zeros(1,length(list.fil))
    isempty(ii) ? jj=collect(1:length(list.fil)) : jj=ii
    for f in jj
       try
        ds=Dataset(list.fil[f])
        close(ds)
       catch e
        println(basename(list.fil[f]))
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

function monthlymean(gdf,m;path0=pwd(),varname="sst")
    list=joinpath.(path0,gdf[m].fil)
    ds=Dataset(list[1])
    tmp=0*ds[varname][:,:,1,1]
    [tmp.+=Dataset(f)[varname][:,:,1,1] for f in list]
    tmp./length(list)
end

function to_monthly_file(arr,m; varname="sst",output_path=tempdir())
    fil=joinpath(output_path,"$(varname)_month$(m).nc")
    ds = Dataset(fil,"c")
    defDim(ds,"i",size(arr,1))
    defDim(ds,"j",size(arr,2))
    v = defVar(ds,varname,Float32,("i","j"))
    arr[ismissing.(arr)].=NaN
    v[:,:] = arr
    close(ds)
    return fil
end

###

read_lon_lat(fil) = begin
    lon=Dataset(fil)["lon"][:]
    lat=Dataset(fil)["lat"][:]
    lon,lat
end

"""
    write_climatology(output_path,year0,year1,lon,lat)

Consolidate monhtly fields into one file with 
- 12 months
- both sst and anom
- coordinate variables
- some metadata
"""
function write_climatology(output_path,year0,year1,lo,la)
    arr=zeros(1440,720,12,2)
    for m in 1:12
        arr[:,:,m,1].=Dataset(joinpath(output_path,"sst_month$(m).nc"))["sst"][:,:]
        arr[:,:,m,2].=Dataset(joinpath(output_path,"anom_month$(m).nc"))["anom"][:,:]
    end

    fi=joinpath(output_path,"OISST_mean_monthly_$(year0)_$(year1).nc")
    #
    ds = NCDataset(fi,"c")
    ds.attrib["title"] = "OISST climatology for $(year0) to $(year1)"
    ds.attrib["author"] = "Gael Forget"
    defDim(ds,"lon",1440); defDim(ds,"lat",720); defDim(ds,"month",12); 
    #
    lon = defVar(ds,"lon",Float32,("lon",))
    lat = defVar(ds,"lat",Float32,("lat",))
    mon = defVar(ds,"month",Float32,("month",))
    sst = defVar(ds,"sst",Float32,("lon","lat","month"))
    anom = defVar(ds,"anom",Float32,("lon","lat","month"))
    #
    lon[:] = lo[:]
    lat[:] = la[:]
    mon[:] = 1:12
    sst[:,:,:] = arr[:,:,:,1]
    anom[:,:,:] = arr[:,:,:,2]
    #
    close(ds)
    fi
end

end 

## 

module coarse_grain

using Statistics, DataFrames, CSV, NCDatasets, Glob

nl=720
#dnl=40 #for 10 degree squares
dnl=8 #for 2 degree squares
nnl=Int(nl/dnl)

@inline areamean(arr,ii,jj) = 
    mean(skipmissing(
        arr[(ii-1)*dnl.+collect(1:dnl),(jj-1)*dnl.+collect(1:dnl)]
        ))

function indices(list)
    arr=Dataset(list.fil[1])["sst"][:,:]

    ii=[ii for ii in 1:nnl*2, jj in 1:nnl]
    jj=[jj for ii in 1:nnl*2, jj in 1:nnl]    
    tmp=[areamean(arr,ii,jj) for ii in 1:nnl*2, jj in 1:nnl]
    kk=findall((!isnan).(tmp))
    (i=ii[kk],j=jj[kk],k=kk)
end

"""
    grid(fil)

Return `(lon=lon,lat=lat,msk=msk,area=area)` based on `fil`.
"""
function grid(fil)
    ds=NCDataset(fil,"r")
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

@inline areaintegral(arr,i::Int,j::Int,G::NamedTuple) = begin
    ii=(i-1)*dnl.+collect(1:dnl)
    jj=(j-1)*dnl.+collect(1:dnl)
    nansum(arr[ii,jj].*G.msk[ii,jj].*G.area[ii,jj])
end

function calc_zm(G::NamedTuple,df)
    gdf_tim=groupby(df, :t)
    arr=NaN*zeros(maximum(df.j),length(gdf_tim))
    for k in minimum(df.j):maximum(df.j)
        area_tmp=[areaintegral(G.msk,x.i,x.j,G) for x in eachrow(gdf_tim[1])]
        area_tmp[gdf_tim[1].j.!==k].=0
        tmp1=[sum(tmp1.sst[:].*area_tmp)/sum(area_tmp) for tmp1 in gdf_tim]
        arr[k,:].=tmp1
    end
    return arr
end

"""
    lowres_merge(;path=dirname(file_root()))

Merge all files found in chosen path.
"""
function merge_files(;path=tempdir(),outputfile="sst_lowres.csv", nam="sst_lowres")
    file_list=glob("$(nam)*csv",path)

    df=DataFrame(i=Int[],j=Int[],t=Int[],sst=Float32[])
    [lowres_append!(df,f) for f in file_list]

    CSV.write(joinpath(tempdir(),outputfile),df)
end

function lowres_append!(df,f)
    tmp=CSV.read(f,DataFrame)
    tmp.t.=parse(Int,split(basename(f),"_")[end][1:8])
    append!(df,tmp)
    return tmp
end

file_root(subdir="sst_lowres_files",filesuff="sst_lowres_") = joinpath(tempdir(),subdir,filesuff)

"""
    lowres_read(;path=tempdir())

Read `sst_lowres.csv`
"""
function lowres_read(;path=tempdir(),fil="lowres_oisst_sst.csv")
    fil=joinpath(path,fil)
    df=CSV.read(fil,DataFrame)
    gdf=groupby(df, [:i, :j])
    kdf=keys(gdf)
    return (df,gdf,kdf)
end

function lowres_index(lon0,lat0,kdf)
    (i,j)=([x.i for x in kdf],[x.j for x in kdf])
	dx=Int(360/maximum(i))
    (ii,jj)=(dx*i.-dx/2,dx*j.-dx/2 .-90)
    d=(ii .-lon0).^2 .+ (jj .-lat0).^2
    findall(d.==minimum(d))[1]
end

lowres_position(ii,jj,kdf) = begin
    (i,j)=([x.i for x in kdf],[x.j for x in kdf])
	dx=Int(360/maximum(i))
	(dx*ii.-dx/2,dx*jj.-dx/2 .-90)
end

end

##

module scenarios
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

##

module timeseries

using DataFrames, Statistics, Dates

function calc(input,list; title="", gdf=nothing)
	if isa(input,DataFrames.GroupKey)
		sst1=gdf[input].sst[:]
	else
		sst1=input[:]
	end
	sst2=repeatclim(sst1,list)
	sst3=anom(sst1,list)

	ttl="SST time series"
	#isa(input,DataFrames.GroupKey) ? ttl=ttl*"for i="*string(input.i)*", j="*string(input.j) : nothing 
	!isempty(title) ?  ttl=title : nothing

    ts=(sst=sst1,clim=sst2,anom=sst3,title=ttl,
    year=list.year,month=list.month,day=list.day)

    tmp1=timeseries.calc_quantile(ts)

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

##

end

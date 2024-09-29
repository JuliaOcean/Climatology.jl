
module SLA_MAIN

using Dataverse
import Climatology: SeaLevelAnomaly, read_Dataset, Dates, write_SLA_PODAAC, write_SLA_CMEMS
import Base: read

#fil=["sla_podaac.nc","sla_cmems.nc"]
function read(x::SeaLevelAnomaly)
    ID=x.name
    path=x.path
    
    fil=string(ID)*".nc"
    sla_file=joinpath(path,fil)
    !isdir(path) ? mkdir(path) : nothing
    if !isfile(sla_file)
        DOI="doi:10.7910/DVN/OYBLGK"
        lst=Dataverse.file_list(DOI)
        Dataverse.file_download(lst,fil,path)
    end

    ds=read_Dataset(sla_file)
    op=(dates=sla_dates(sla_file),)
    SeaLevelAnomaly(name=x.name,path=path,data=[ds],options=op)
end

podaac_date(n)=Dates.Date("1992-10-05")+Dates.Day(5*n)
podaac_sample_dates=podaac_date.(18:73:2190)
cmems_date(n)=Dates.Date("1993-01-01")+Dates.Day(1*n)
podaac_all_dates=podaac_date.(1:2190)
cmems_all_dates=cmems_date.(1:10632)

sla_dates(fil) = ( fil=="sla_podaac.nc" ? podaac_all_dates : cmems_all_dates)

end

##

module SLA_PODAAC

using Dates, DataStructures
import Climatology: Downloads, read_Dataset, write_SLA_PODAAC

#note : this need up-to-date credentials in ~/.netrc and ~/.ncrc
url0="https://opendap.earthdata.nasa.gov/collections/C2270392799-POCLOUD/granules/"
##url0="https://podaac-tools.jpl.nasa.gov/drive/files/allData/merged_alt/L4/cdr_grid/"
path0=joinpath(pwd(),"SEA_SURFACE_HEIGHT_ALT_GRIDS_L4_2SATS_5DAY_6THDEG_V_JPL2205")*"/"

#url1="https://opendap.earthdata.nasa.gov/collections/C2102959417-POCLOUD/granules/"
#url=url1*"oscar_currents_interim_20230101.nc"
#path1=joinpath(pwd(),"OSCAR_L4_OC_INTERIM_V2.0")*"/"


"""
    get_grid(;url=url0,file="",range_lon=360.0.+(-35.0,-22),range_lat=(34.0,45))

"""
function get_grid(;url=url0,file="",range_lon=360.0.+(-35.0,-22),range_lat=(34.0,45))

    if !isempty(file)
        fil=file
        ds=read_Dataset(fil)
        lon=ds["lon"][:]
        lat=ds["lat"][:]
    else
        url=url*"ssh_grids_v2205_1992101012.dap.nc"
    #    fil=joinpath(tempdir(),"ssh_grids_v2205_1992101012.dap.nc")
        fil=Downloads.download(url)
        ds=read_Dataset(fil)
        lon=Float64.(ds["Longitude"][:])
        lat=Float64.(ds["Latitude"][:])        
    end

    ii=findall( (lon.>range_lon[1]) .& (lon.<range_lon[2]) )
    jj=findall( (lat.>range_lat[1]) .& (lat.<range_lat[2]) )

    (lon=lon,lat=lat,ii=ii,jj=jj,nt=2190,file=fil)
end

function file_name(n)
    d0=Date("1992-10-05")
    d=d0+Dates.Day(n*5)
    dtxt=string(d)
    "ssh_grids_v2205_"*dtxt[1:4]*dtxt[6:7]*dtxt[9:10]*"12.nc" #".dap.nc"
end

function read_slice(url,gr)
    #fil=Downloads.download(url)
    #ds=read_Dataset(fil)
    ds=read_Dataset(url)
    SLA=ds["SLA"][gr.ii,gr.jj,1]
    SLA[ismissing.(SLA)].=NaN
    Float64.(SLA)
end

"""
    SLA_PODAAC.subset()

For download directions, see [this site](https://podaac.jpl.nasa.gov/dataset/SEA_SURFACE_HEIGHT_ALT_GRIDS_L4_2SATS_5DAY_6THDEG_V_JPL2205)

```
SLA_PODAAC.subset(; read_from_file=SLA.file)
```
"""
function subset(;
    path0="SEA_SURFACE_HEIGHT_ALT_GRIDS_L4_2SATS_5DAY_6THDEG_V_JPL2205/",
    username="unknown",
    password="unknown",
    range_lon=360.0.+(-35.0,-22),
    range_lat=(34.0,45),
    read_from_file="",
    save_to_file=false,
    )

    if !isempty(read_from_file)
        gr=SLA_PODAAC.get_grid(file=read_from_file)
        ds=read_Dataset(read_from_file)
        i0=1; i1=gr.nt
        data=ds["SLA"][:,:,:]
    else
        gr=get_grid(range_lon=range_lon,range_lat=range_lat)
        i0=1; i1=gr.nt    
        data=zeros(length(gr.ii),length(gr.jj),i1-i0+1)
        for n=i0:i1
            mod(n,100)==0 ? println(n) : nothing
            data[:,:,n-i0+1]=read_slice(path0*file_name(n),gr)
        end
    end

    #show(gr)
    save_to_file ? write_SLA_PODAAC(gr,data) : data
end

end #module SLA_PODAAC

module SLA_CMEMS

using URIs, DataStructures
import Climatology: Downloads, read_Dataset, write_SLA_CMEMS

"""
    SLA_CMEMS.subset()

For download directions, see [this site](https://marine.copernicus.eu)

For data documentation, see [this page](https://data.marine.copernicus.eu/product/SEALEVEL_GLO_PHY_L4_MY_008_047/description)

```
SLA_CMEMS.subset(username=username,password=password)
```
"""
function subset(;
    var="cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1D",
    username="unknown",
    password="unknown",
    range_lon=(-35.0,-22),
    range_lat=(34.0,45),
    read_from_file="",
    save_to_file=false,
    )

    if !isempty(read_from_file)
        ds=read_Dataset(read_from_file)
        SSH=ds["SLA"]
        lon=ds["lon"][:]
        lat=ds["lat"][:]
    else
        url="https://my.cmems-du.eu/thredds/dodsC/"*var
        url2 = string(URI(URI(url),userinfo = string(username,":",password)))
        ds = read_Dataset(url2)
        SSH=ds["sla"]
        lon=ds["longitude"][:]
        lat=ds["latitude"][:]
    end

    ii=findall( (lon.>range_lon[1]) .& (lon.<range_lon[2]) )
    jj=findall( (lat.>range_lat[1]) .& (lat.<range_lat[2]) )
    data = SSH[ii,jj,:]

    #show(gr)
    save_to_file ? write_SLA_CMEMS(lon[ii],lat[jj],data) : data
end

end #module SLA_CMEMS

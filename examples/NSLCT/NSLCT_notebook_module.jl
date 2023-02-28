
module NSLCT

using HTTP, JSON
using DataFrames, Dates, Statistics

slcp_base_url = "https://sealevel-dataservices.jpl.nasa.gov/nexus"
timeSeriesSpark_base_url="https://sealevel-dataservices.jpl.nasa.gov/nexus/timeSeriesSpark"

sea_base_url = "https://sealevel-dataservices.jpl.nasa.gov/edge/ws/search/tidegauge"
ipcc_base_url = "https://sealevel-dataservices.jpl.nasa.gov/edge/ws/search/projection"
taskforce_base_url = "https://sealevel-dataservices.jpl.nasa.gov/edge/ws/search/projection"

#tile_base_url = "https://sealevel-dataservices.jpl.nasa.gov/onearth/wmts"
#https://sealevel-dataservices.jpl.nasa.gov/onearth/wmts/epsg4326/all//SSH_ECCO_version4_release4/default/2016-01-15/4km/2/

"""
    get_list()

Get list of data sets from `slcp_base_url`, parse, and return as a `DataFrame`.
"""
function get_list()
    response = HTTP.get(slcp_base_url*"/list")
    tmp=JSON.parse(String(response.body))
    tmp[1]

    list=DataFrame([Symbol(s) => [x[s] for x in tmp] for s in keys(tmp[1])])
    list.t0=[Date(t.iso_start[1:10]) for t in eachrow(list)]
    list.t1=[Date(t.iso_end[1:10]) for t in eachrow(list)]
    list.dt=list.t1-list.t0

    sort!(list,:dt,rev=true)
    list[:,[:t0,:t1,:title]]
end

"""
    search()

Search list of data sets .

```
list=NSLCT.get_list()
NSLCT.search(list,"sea","surface","height")
NSLCT.search(list,"sea surface height")
```
"""
#search(list,word::String) = list[occursin.(Ref(lowercase(word)),lowercase.(list.title)),:]
search(list,word::AbstractString) = occursin(" ",word) ? search(list,split(word)) : list[occursin.(Ref(lowercase(word)),lowercase.(list.title)),:]
search(list,words...) = length(words)==1 ? search(list,words[1]) : search(search(list,words[1]),words[2:end])
search(list,words) = length(words)==1 ? search(list,words[1]) : search(search(list,words[1]),words[2:end])

#pre-processing
pp_ssh_ts(x)=x.-median(x)

#ecco_ssh_ts=get_ssh_ts("SSH_ECCO_version4_release4")
#measures_ssh_ds=get_ssh_ts("SEA_SURFACE_HEIGHT_ALT_GRIDS_L4_2SATS_5DAY_6THDEG_V_JPL2205_Monthly")
function get_ssh_ts(dataset_name="SSH_ECCO_version4_release4")
    # Temporal bounds
    start_time = DateTime(1992,10,1)
    start_time_txt = Dates.format(start_time,"yyyy-mm-ddTHH:MM:SSZ")
    end_time = DateTime(2018,1,1)
    end_time_txt = Dates.format(end_time,"yyyy-mm-ddTHH:MM:SSZ")

    # Spatio-Temporal bounding box
    limits = (  min_lon=-150, max_lon=-77, min_lat=-8, max_lat=8, 
                min_tim=start_time_txt, max_tim=end_time_txt )

    # Download data
    params = Dict("ds"=>dataset_name, 
    "minLon"=>limits.min_lon, "minLat"=>limits.min_lat, 
    "maxLon"=>limits.max_lon, "maxLat"=>limits.max_lat,
    "startTime"=>limits.min_tim, "endTime"=>limits.max_tim)        
    response=HTTP.request("GET", timeSeriesSpark_base_url; query=params)

    #url0=timeSeriesSpark_base_url*"?ds=$(ecco_ssh_ds)&"*
    #    "minLon=$(limits.min_lon)&maxLon=$(limits.max_lon)&"*
    #    "minLat=$(limits.min_lat)&maxLat=$(limits.max_lat)&"*
    #    "startTime=$(limits.min_tim)&endTime=$(limits.max_tim)"
    #response = HTTP.get(url0)

    #url00
    #response=HTTP.get(url00)

    # Reorganize data into DataFrame
    tmp=JSON.parse(String(response.body))
    ssh_ts=DataFrame([Symbol(s) => [x[1][s] for x in tmp["data"]] for s in keys(tmp["data"][1][1])])
    tmp=[Date(x[1:10]) for x in ssh_ts.iso_time]
    ssh_ts.tim=year.(tmp)+month.(tmp)./12 .-0.5/12
    ssh_ts.sla=pp_ssh_ts(ssh_ts.mean)

    ssh_ts[:,[:tim,:sla]]
end

function get_projection()
    ne_va_id=299
    params = (psmsl_id=ne_va_id,task_force=true)
    #params = Dict("psmsl_id"=>ne_va_id,"task_force"=>true)

    response=HTTP.request("GET", taskforce_base_url; query=params)
end

end #module NSLCT

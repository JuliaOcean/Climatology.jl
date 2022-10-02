# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.3
#   kernelspec:
#     display_name: Julia 1.8.1
#     language: julia
#     name: julia-1.8
# ---

using HTTP, JSON
using DataFrames, Dates
using GLMakie

slcp_base_url = "https://sealevel-dataservices.jpl.nasa.gov/nexus"
response = HTTP.get(slcp_base_url*"/list")
tmp=JSON.parse(String(response.body))
tmp[1]

# +
list=DataFrame([Symbol(s) => [x[s] for x in tmp] for s in keys(tmp[1])])
list.t0=[Date(t.iso_start[1:10]) for t in eachrow(list)]
list.t1=[Date(t.iso_end[1:10]) for t in eachrow(list)]
list.dt=list.t1-list.t0

sort!(list,:dt,rev=true)
list=list[:,[:t0,:t1,:title]]
summary(list)

# +
#search(list,word::String) = list[occursin.(Ref(lowercase(word)),lowercase.(list.title)),:]
search(list,word::AbstractString) = occursin(" ",word) ? search(list,split(word)) : list[occursin.(Ref(lowercase(word)),lowercase.(list.title)),:]
search(list,words...) = length(words)==1 ? search(list,words[1]) : search(search(list,words[1]),words[2:end])
search(list,words) = length(words)==1 ? search(list,words[1]) : search(search(list,words[1]),words[2:end])

search(list,"sea","surface","height")
search(list,"sea surface height")

# +
# Dataset names ingested into data analysis tool
ecco_ssh_ds = "SSH_ECCO_version4_release4"
measures_ssh_ds = "SEA_SURFACE_HEIGHT_ALT_GRIDS_L4_2SATS_5DAY_6THDEG_V_JPL2205_Monthly"

# Temporal bounds
start_time = DateTime(1992,10,1)
end_time = DateTime(2018,1,1)

# Spatial bounding box
bb = ( min_lon=-150, max_lon=-77, min_lat=-8, max_lat=8 )
# -

# HTTP.get("https://sealevel-dataservices.jpl.nasa.gov/nexus/timeSeriesSpark?ds=SSH_ECCO_version4_release4&minLon=-150&minLat=-8&maxLon=-77&maxLat=8&startTime=1992-10-01T00:00:00Z&endTime=2018-01-01T00:00:00Z")

# +
#1. get "spatial" mean from timeSeriesSpark
if true
    params = Dict("ds"=>ecco_ssh_ds, "minLon"=>bb.min_lon, "minLat"=>bb.min_lat, 
    "maxLon"=>bb.max_lon, "maxLat"=>bb.max_lat,
    "startTime"=>"1992-10-01T00:00:00Z","endTime"=>"2018-01-01T00:00:00Z")
    timeSeriesSpark_base_url="https://sealevel-dataservices.jpl.nasa.gov/nexus/timeSeriesSpark"
    response=HTTP.request("GET", timeSeriesSpark_base_url; query=params)
else
    url0="https://sealevel-dataservices.jpl.nasa.gov/nexus/timeSeriesSpark?"*
        "ds=$(ecco_ssh_ds)&minLon=$(bb.min_lon)&maxLon=$(bb.max_lon)&"*
        "minLat=$(bb.min_lat)&maxLat=$(bb.max_lat)&"*
        "startTime=$(start_time)Z&endTime=$(end_time)Z"
    response = HTTP.get(url0)
end
tmp=JSON.parse(String(response.body))

#2. reorganize data into DataFrame
ecco_ssh_ts=DataFrame([Symbol(s) => [x[1][s] for x in tmp["data"]] for s in keys(tmp["data"][1][1])])
tmp=[Date(x[1:10]) for x in ecco_ssh_ts.iso_time]
ecco_ssh_ts.t=year.(tmp)+month.(tmp)./12 .-0.5/12

#3. visualize data data
names(ecco_ssh_ts)
# -

lines(ecco_ssh_ts.t,ecco_ssh_ts.mean)

# +
ne_va_id=299
params = (psmsl_id=ne_va_id,task_force=true)
#params = Dict("psmsl_id"=>ne_va_id,"task_force"=>true)

taskforce_base_url = "https://sealevel-dataservices.jpl.nasa.gov/edge/ws/search/projection"
response=HTTP.request("GET", taskforce_base_url; query=params)

# +
slcp_base_url = "https://sealevel-dataservices.jpl.nasa.gov/nexus"
sea_base_url = "https://sealevel-dataservices.jpl.nasa.gov/edge/ws/search/tidegauge"
ipcc_base_url = "https://sealevel-dataservices.jpl.nasa.gov/edge/ws/search/projection"

#tile_base_url = "https://sealevel-dataservices.jpl.nasa.gov/onearth/wmts"
#https://sealevel-dataservices.jpl.nasa.gov/onearth/wmts/epsg4326/all//SSH_ECCO_version4_release4/default/2016-01-15/4km/2/
# -





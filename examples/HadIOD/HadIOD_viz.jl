
#data source : https://www.metoffice.gov.uk/hadobs/hadiod/
#reference : Atkinson, C. P., N. A. Rayner, J. J. Kennedy and S. A. Good, 2014. An Integrated Database of Ocean Temperature and Salinity Observations, Journal of Geophysical Research: Oceans, 119, 7139-7163, doi:10.1002/2014JC010053

#data download : https://www.metoffice.gov.uk/hadobs/hadiod/download-hadiod1-2-0-0.html#hadiod_historic

using NCDatasets, GLMakie, Random, Downloads

## Download file if needed

pth=joinpath(tempdir(),"hadiod_sample")
fil="hadiod1200_19841224.nc" #file to be read
prof=6872 #profile to be displayed

if !isfile(joinpath(pth,fil))
    url="https://www.metoffice.gov.uk/hadobs/hadiod/data/hadiod1200/main/HadIOD1200.data.1984.zip"
    zipfile=Downloads.download(url)
    !isdir(pth) ? mkdir(pth) : nothing
    run(`unzip $(zipfile) $(fil) -d $(pth)`)
end

#see also, for alternative corrections:
#fil1="1984/hadiod1200_btcorrs_19841224.nc"
#fil1="1984/hadiod1200_shipcorrs_19841224.nc"

ds=Dataset(joinpath(pth,fil))

## platform ID, position, time, number of samples

plat_id=[string(ds["plat_id"][:,i]...) for i in 1:6876]
lon=ds["lon"][:]
lat=ds["lat"][:]
time=ds["time"][:]
rowSize=ds["rowSize"][:]

## all SST values + corrections

tmp=loadragged(ds["temp"],:)
tmp_type_corr=loadragged(ds["temp_type_corr"],:)
tmp_plat_corr=loadragged(ds["temp_plat_corr"],:)
sst=[tmp[i][1] for i in 1:length(tmp)]
BT=[tmp_type_corr[i][1] for i in 1:length(tmp)]
BT[ismissing.(BT)].=0.0
BP=[tmp_plat_corr[i][1] for i in 1:length(tmp)]
BP[ismissing.(BP)].=0.0

#sst=Float64.(skipmissing(sst+BT+BP))
sst=Float64.(sst)
sst_corr=Float64.(BT+BP)

## one profile + corrections

BT=loadragged(ds["temp_type_corr"],prof)[1]
BT[ismissing.(BT)].=0.0
BP=loadragged(ds["temp_plat_corr"],prof)[1]
BP[ismissing.(BP)].=0.0

temp=Float64.(loadragged(ds["temp"],prof)[1])
potemp=Float64.(loadragged(ds["potemp"],prof)[1])

#sal=Float64.(loadragged(ds["sal"],prof)[1])
sal=Float64[]

depth=Float64.(loadragged(ds["depth"],prof)[1])
DC=Float64.(loadragged(ds["depth_corr"],prof)[1])
DC[ismissing.(DC)].=0.0
depth=depth+DC

close(ds)

## plotting

fig=Figure()
ax1=Axis(fig[1,1],xlabel="potemp",ylabel="depth",
    title="profile no "*string(prof))
lines!(ax1,potemp,-depth,label="uncorrected")
scatter!(ax1,potemp+BT+BP,-depth,marker=:x,markersize=10.0,label="corrected")
axislegend(ax1,position=:rb)

ax2=Axis(fig[1,2],xlabel="latitude",ylabel="temp[1] (i.e., sst)",
    title="sample from "*basename(fil))
nn=length(sst)
nn>200 ? ii=shuffle(1:nn)[1:200] : ii=1:nn
scatter!(ax2,lat[ii],sst[ii],label="uncorrected")
scatter!(ax2,lat[ii],sst[ii]+sst_corr[ii],marker=:x,markersize=20.0,label="corrected")
axislegend(ax2,position=:rt)

fig


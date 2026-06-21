
"""
    write_SST_climatology(output_path,year0,year1,lon,lat)

Consolidate monhtly fields into one file with 
- 12 months
- both sst and anom
- coordinate variables
- some metadata
"""
function write_SST_climatology(output_path,year0,year1,lo,la)
    arr=zeros(1440,720,12,2)
    for m in 1:12
        arr[:,:,m,1].=Dataset(joinpath(output_path,"sst_month$(m).nc"))["sst"][:,:]
        arr[:,:,m,2].=Dataset(joinpath(output_path,"anom_month$(m).nc"))["anom"][:,:]
    end

    fi=joinpath(output_path,"OISST_mean_monthly_$(year0)_$(year1).nc")
    #
    ds = Dataset(fi,"c")
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

function to_monthly_file(arr,m; varname="sst",output_path=SST_demo_path)
    fil=joinpath(output_path,"$(varname)_month$(m).nc")
    ds = read_Dataset(fil,"c")
    defDim(ds,"i",size(arr,1))
    defDim(ds,"j",size(arr,2))
    v = defVar(ds,varname,Float32,("i","j"))
    arr[ismissing.(arr)].=NaN
    v[:,:] = arr
    close(ds)
    return fil
end

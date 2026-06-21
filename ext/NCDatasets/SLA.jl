
function write_SLA_PODAAC(gr,data)
    fil=joinpath(tempdir(),"podaac_sla_dev.nc")
    Dataset(fil,"c",attrib = OrderedDict("title" => "Azores Regional Subset")) do ds
        defVar(ds,"SLA",data,("lon","lat","time"), attrib = OrderedDict(
            "units" => "m", "long_name" => "Sea Level Anomaly",
            "comments" => "source is https://sealevel.nasa.gov/data/dataset/?identifier=SLCP_SEA_SURFACE_HEIGHT_ALT_GRIDS_L4_2SATS_5DAY_6THDEG_V_JPL2205_2205")),
        defVar(ds,"lon",gr.lon[gr.ii],("lon",), attrib = OrderedDict(
            "units" => "degree", "long_name" => "Longitude"))
        defVar(ds,"lat",gr.lat[gr.jj],("lat",), attrib = OrderedDict(
            "units" => "degree", "long_name" => "Latitude"))
        end
    println("File name :")
    fil
end

function write_SLA_CMEMS(lon,lat,data)
    fil=joinpath(tempdir(),"cmems_sla_dev.nc")
    read_Dataset(fil,"c",attrib = OrderedDict("title" => "Azores Regional Subset")) do ds
        defVar(ds,"SLA",data,("lon","lat","time"), attrib = OrderedDict(
            "units" => "m", "long_name" => "Sea Level Anomaly",
            "comments" => "source is https://my.cmems-du.eu")),
        defVar(ds,"lon",lon,("lon",), attrib = OrderedDict(
            "units" => "degree", "long_name" => "Longitude"))
        defVar(ds,"lat",lat,("lat",), attrib = OrderedDict(
            "units" => "degree", "long_name" => "Latitude"))
        end
    println("File name :")
    fil
end


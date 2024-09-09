module ClimatologyNCDatasetsExt

    import Climatology: ECCO, load, read_Dataset, ECCOdiags_to_nc, write_SST_climatology, SST_demo_path, to_monthly_file
    import MeshArrays, Printf 
    import MeshArrays: GridSpec, Tiles, GridLoadVar, GRID_LLC90    
    import NCDatasets: Dataset, defDim, defVar

    read_Dataset(args...;kwargs...)=Dataset.(args...;kwargs...)

    function ECCOdiags_to_nc(;path_in=".",file_out=tempname()*".nc",
            year1=1960,nt=771,title="this is a test file")
        nsec=23
        nlatMT=179
        nlatZM=90

        ds = Dataset(file_out,"c")

        ## dimension definitions

        defDim(ds,"x",30); defDim(ds,"y",30); defDim(ds,"tile",117)
        #defDim(ds,"lon",360); defDim(ds,"lat",180)
        defDim(ds,"latMT",nlatMT); defDim(ds,"latZM",nlatZM)
        defDim(ds,"section",23); defDim(ds,"depth",50)
        defDim(ds,"time",nt); defDim(ds,"date",nt)
        defDim(ds,"month",12)
        defDim(ds,"time_clim",14)
        defDim(ds,"level_clim",6)

        ## dimension variables

        latMT=-89.0:89.0
        v = defVar(ds,"latMT",Float32,("latMT",)); v[:]=latMT
        dlat=2.0; latZM=vec(-90+dlat/2:dlat:90-dlat/2)
        v = defVar(ds,"latZM",Float32,("latZM",)); v[:]=latZM
        v = defVar(ds,"time",Float32,("time",)); v[:]=year1-0.5/12 .+ (1:nt)/12
        v = defVar(ds,"month",Float32,("month",)); v[:]=1:12

        time_clim=string.([:J,:F,:M,:A,:M,:J,:J,:A,:S,:O,:N,:D, :annual, :std])
        v = defVar(ds,"time_clim",String,("time_clim",));  v[:]=time_clim
        level_clim=[1 10 20 29 38 44]
        v = defVar(ds,"level_clim",Int,("level_clim",));  v[:]=level_clim

        ds.attrib["title"] = title

        ## simple array data

        list_in=ECCO.diagnostics_set1(path_in)
        for k in 1:size(list_in,1)
    #        println(list_in[k,"name"])
            tmp=load(list_in[k,"file"])["single_stored_object"]
            v = defVar(ds,list_in[k,"name"],Float32,list_in[k,"dims"])
            v.attrib["units"] = list_in[k,"units"]
            if isa(tmp[1],NamedTuple)
                [v[s,:,:] = tmp[s].val for s in 1:nsec]
                w = defVar(ds,"transport_name",String,("section",))
                [w[s] = split(tmp[s].nam,".")[1] for s in 1:nsec]
            else
                if (ndims(tmp)>1)&&(size(tmp,1)==nt)
                    v[:] = permutedims(tmp,ndims(tmp):-1:1)
                else
                    v[:] = tmp
                end
            end
        end

        ## MeshArrays data

        γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
        τ=Tiles(γ,30,30)

        list_dims=("x","y","tile")
        lon_clim = Tiles(τ,GridLoadVar("XC",γ))
        v = defVar(ds,"lon_clim",Float64,list_dims)
        [v[:,:,tile]=lon_clim[tile] for tile in 1:117]
        lat_clim = Tiles(τ,GridLoadVar("YC",γ))
        v = defVar(ds,"lat_clim",Float64,list_dims)
        [v[:,:,tile]=lat_clim[tile] for tile in 1:117]

        list_in=ECCO.diagnostics_set2(path_in)
        for k in 1:size(list_in,1)
    #        println(list_in[k,"name"])
            tmp=load(list_in[k,"file"])
            list_dims=("x","y","tile","time_clim")
            v = defVar(ds,list_in[k,"name"],Float64,list_dims)

            for time in 1:12
                tmp1=Tiles(τ,tmp["mon"][:,time]); 
                [v[:,:,tile,time]=tmp1[tile] for tile in 1:117]
            end
            tmp1=Tiles(τ,tmp["mean"]); [v[:,:,tile,13]=tmp1[tile] for tile in 1:117]
            tmp1=Tiles(τ,tmp["std"]); [v[:,:,tile,14]=tmp1[tile] for tile in 1:117]
        end

        list_in=ECCO.diagnostics_set3(path_in)
        for k in 1:size(list_in,1)
            #println(list_in[k,"name"])
            list_dims=("x","y","tile","level_clim","time_clim")
            v = defVar(ds,list_in[k,"name"],Float64,list_dims)
            for level in 1:length(level_clim)
                suff=Printf.@sprintf("%02d.jld2",level_clim[level])
                file=list_in[k,"file"][1:end-7]*suff
                tmp=load(file)
                #println("file = "*file)
                for time in 1:12
                    tmp1=Tiles(τ,tmp["mon"][:,time]); 
                    [v[:,:,tile,level,time]=tmp1[tile] for tile in 1:117]
                end
                tmp1=Tiles(τ,tmp["mean"])
                [v[:,:,tile,level,13]=tmp1[tile] for tile in 1:117]
                tmp1=Tiles(τ,tmp["std"]) 
                [v[:,:,tile,level,14]=tmp1[tile] for tile in 1:117]
            end
        end

        close(ds)

        file_out
    end

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
    
end



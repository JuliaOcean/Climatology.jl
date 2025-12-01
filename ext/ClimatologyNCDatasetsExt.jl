module ClimatologyNCDatasetsExt

    import Climatology: ECCO, load, read_Dataset, ECCOdiags_to_nc
    import Climatology: write_SST_climatology, SST_demo_path, to_monthly_file
    import Climatology: write_SLA_PODAAC, write_SLA_CMEMS
    import MeshArrays, Printf 
    import MeshArrays: GridSpec, Tiles, GridLoadVar    
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

        γ=GridSpec(ID=:LLC90)
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

    ##

    import Climatology: read_IAP, file_IAP, write_H_to_T

    import NCDatasets, MeshArrays
    using DataStructures: OrderedDict
    using MeshArrays: gridmask, Integration

    """
        file_IAP(path,y,m)
    """
    file_IAP(path,y,m)=begin
        mm=(m<10 ? "0$m" : "$m")        
        joinpath(path,"IAPv4_Temp_monthly_1_6000m_year_$(y)_month_$(mm).nc")
    end

    """
        read_IAP(F,var,tim,tmp=[])

    ```
    using Climatology, NCDatasets, MeshArrays

    p0="IAPv4_IAP_Temperature_gridded_1month_netcdf/monthly/"
    fil=Climatology.file_IAP(p0,"2023","12")
    depth=Dataset(fil)["depth_std"][:]
    temp=Climatology.read_IAP(fil,"temp",1,[])
    mask=1.0*(!ismissing).(temp)

    G=Gris_simple.GridLoad_lonlatdep(depth,mask)
    tmp=zeros(G.XC.grid)*ones(length(depth))
    Climatology.read_IAP(fil,"temp",1,tmp)
    ```
    """
    function read_IAP(F,var,tim,tmp=[])
        fil=F

        ds=Dataset(fil)
        temp=permutedims(ds[var][:,:,:],(2,3,1))
        close(ds)
        temp[findall(ismissing.(temp))].=0
        temp[findall(isnan.(temp))].=0
        if !isempty(tmp)
            tmp.=read(Float32.(temp),tmp)
            tmp
        else
            temp
        end
    end

    """
        write_H_to_T(file::String,M::gridmask,G::NamedTuple,H::Array)

    Write `H / Integration.volumes(M,G)` to file.

    ```
    using Climatology, NCDatasets, MeshArrays

    G=MeshArrays.Grids_simple.GridLoad_lonlatdep(depth,mask)
    M=Integration.define_sums(grid=G,regions=(10,5))
    H=ones(length(M.names),length(M.depths),3)
    V=Integration.volumes(M,G)
    Climatology.write_H_to_T(tempname()*".nc",M,G,H,V)
    ```
    """
    function write_H_to_T(file::String,M::gridmask,G::NamedTuple,H::Array,V::Array)
        nb,nz,nt=size(H)
        inv_vol=1.0./V
        #inv_vol[V.==0].=0
        pos=gridpos(M,(10,5))
        arr2d=zeros(36,32)
        arr3d=zeros(36,32,nz)
        arr4d=zeros(36,32,nz,nt)
        println(nz)

        ds = Dataset(file,"c")
        defDim(ds,"lon",36); defDim(ds,"lat",32); 
        defDim(ds,"dep",size(H,2)); defDim(ds,"tim",size(H,3));
        ds.attrib["title"] = "this is a test file"

        dlo=10; dla=5;
        lons=collect(-180:dlo:180); lons=0.5*(lons[1:end-1]+lons[2:end])
        lats=[-90 ; -75:dla:75 ; 90]; lats=0.5*(lats[1:end-1]+lats[2:end])
        vlo=defVar(ds,"lon",lons,("lon",), attrib = OrderedDict("units" => "degree", "long_name" => "Longitude"))
        vla=defVar(ds,"lat",lats,("lat",), attrib = OrderedDict("units" => "degree", "long_name" => "Latitude"))

        v1 = defVar(ds,"volume",Float32,("lon","lat","dep"), attrib = OrderedDict("units" => "m^3",))

        arr3d.=0 
        for ii in 1:nb
        i,j=pos[ii]
        [arr3d[i,j,k]=V[ii,k] for k in 1:nz]
        end

        v1[:,:,:] = arr3d

        v = defVar(ds,"temperature",Float32,("lon","lat","dep","tim"), attrib = OrderedDict("units" => "degree Celsius",))
        v.attrib["comments"] = "this is a string attribute with Unicode Ω ∈ ∑ ∫ f(x) dx"

        arr4d.=0
        for t in 1:nt
        for ii in 1:nb
        i,j=pos[ii]
        [arr4d[i,j,k,t]=H[ii,k,t]*inv_vol[ii,k] for k in 1:nz]
        end
        end

        v[:,:,:,:] = arr4d

        close(ds)
        file
    end

    """
        gridpos(M::gridmask,res::Tuple)

    ```
    gridpos(M,(10,5))
    ```
    """
    gridpos(M::gridmask,res::Tuple)=begin
        n=length(M.names)
        allpos=fill((0,0),n)
        for i in 1:n
            t1=split(M.names[i],"Nto")
            t2=split(t1[2],"N_")
            t3=split(t2[2],"Eto")
            t4=split(t3[2],"E")
            tt=[t1[1] t2[1] t3[1] t4[1]]
            tt=parse.(Ref(Int),tt)

            dlo=res[1]; dla=res[2]
            lons=collect(-180:dlo:180)
            lats=[-90 ; -75:dla:75 ; 90]
            thispos=(findall(lons.==tt[3])[1],findall(lats.==tt[1])[1])

            allpos[i]=thispos
        end
        allpos
    end

##

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


end

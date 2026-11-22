# SST (OISST / ERSST) tests. Included from runtests.jl; see TEST_GROUP and SKIP_DOWNLOADS there.
#
# With SKIP_DOWNLOADS=true no daily files are downloaded, so only these parts run:
#   A. always:                    file lists, read-only file discovery, table_info on synthetic data
#   B. downloaded daily files:    skipped
#   C. cached Zenodo table:       runs if the DataDep OISST-stats1 is already in the scratch space
#   D. needs both B and C:        skipped
# Without SKIP_DOWNLOADS everything runs, as before.

# true if the DataDep is already in the scratch space (checking never downloads)
_datadep_cached(name,file) = isfile(joinpath(Climatology.DataDeps.determine_save_path(name),file))

@testset "SST" begin

    input_path=tempname()

    list_downloaded=String[]
    if !_SKIP_DOWNLOADS
        list_downloaded=SST_processing.download_files(path=input_path,short_demo=true)
        @test ispath(input_path)

        output_path=SST_processing.coarse_grain(path=input_path,short_demo=true)
        @test isfile(output_path)
    end
    have_files=!isempty(list_downloaded)
    have_stats=!_SKIP_DOWNLOADS || _datadep_cached("OISST-stats1","lowres_oisst_sst_10.0.csv")

    #@everywhere using Climatology, NCDatasets
    #output_path=SST_processing.monthly_climatology(path=input_path)
    #mv(output_file,joinpath(input_path,basename(output_file)))

    ## A. file lists, read-only file discovery, table checks on synthetic data

    (fil1,fil2)=SST_FILES.file_lists(path=input_path)
    whole_list=SST_FILES.CSV.read(fil1,SST_FILES.DataFrame)

    list_pb=SST_FILES.test_files(whole_list)
    @test isa(list_pb,Vector)

    (fil1,fil2)=SST_FILES.ersst_file_lists(path=input_path)
    @test isfile(fil1)

    pld=mktempdir()
    for f in ("198201/oisst-avhrr-v02r01.19820101.nc",                 # final only
              "202611/oisst-avhrr-v02r01.20261114.nc",                 # both final and preliminary
              "202611/oisst-avhrr-v02r01.20261114_preliminary.nc",
              "202611/oisst-avhrr-v02r01.20261115_preliminary.nc",     # preliminary only
              "202611/notes.txt","202611/oisst-avhrr-v02r01.2026111.nc",   # stray names
              "other/oisst-avhrr-v02r01.20200101.nc")                  # not a YYYYMM folder
        mkpath(dirname(joinpath(pld,f))); touch(joinpath(pld,f))
    end
    ld=SST_FILES.local_days(pld)
    @test ld.dates==[Date(1982,1,1),Date(2026,11,14),Date(2026,11,15)]
    @test ld.preliminary==[false,false,true]
    @test ld.files[3]==joinpath(pld,"202611","oisst-avhrr-v02r01.20261115.nc") && !isfile(ld.files[3])
    @test ld.actual[2]==ld.files[2] && ld.actual[3]==ld.files[3][1:end-3]*"_preliminary.nc"
    @test all(isfile,ld.actual)
    @test isempty(SST_FILES.local_days(mktempdir()).dates)

    ## A2. ERSST map data from synthetic monthly files (no downloads)

    erd=mktempdir()
    for (fname,val) in (("ersst.v5.202001.nc",1.5f0),("ersst.v5.202002.nc",-2.0f0))
        ds=NCDatasets.Dataset(joinpath(erd,fname),"c")
        NCDatasets.defDim(ds,"lon",4); NCDatasets.defDim(ds,"lat",3)
        NCDatasets.defDim(ds,"lev",1); NCDatasets.defDim(ds,"time",1)
        NCDatasets.defVar(ds,"lon",Float32.(0:3),("lon",))
        NCDatasets.defVar(ds,"lat",Float32.(-1:1),("lat",))
        v=NCDatasets.defVar(ds,"ssta",Float32,("lon","lat","lev","time");fillvalue=Float32(-9999))
        arr=Array{Union{Missing,Float32}}(fill(val,4,3,1,1)); arr[1,1,1,1]=missing
        v[:,:,:,:]=arr
        close(ds)
    end

    md=SST_FILES.ersst_map_data(Date(2020,1,15);path=erd)
    @test md.metadata.fallback==false
    @test md.title=="ersst.v5.202001.nc"
    @test md.metadata.file==joinpath(erd,"ersst.v5.202001.nc")
    @test isnan(md.field[1,1]) && md.field[2,2]==1.5
    @test size(md.field)==(4,3)

    md2=SST_FILES.ersst_map_data(Date(2020,5,1);path=erd)
    @test md2.metadata.fallback==true
    @test md2.metadata.file==joinpath(erd,"ersst.v5.202002.nc")
    @test occursin("not on disk",md2.title)

    md3=SST_FILES.ersst_map_data(Date(2019,1,1);path=erd)
    @test md3.metadata.file==joinpath(erd,"ersst.v5.202001.nc")

    @test_throws ErrorException SST_FILES.ersst_map_data(Date(2020,1,1);path=mktempdir())

    if !_SKIP_DOWNLOADS
        @test isa(plot(SSTdiag("",md.title,:ersst_map;map_data=md)),Figure)
    end

    ## A3. map_point: land (no series) and decimal years (no downloads, no basemap)

    mpl=(field=Union{Missing,Float64}[1.0 missing; -1.0 0.5],lon=[205.0,206.0],lat=[44.0,45.0],
         colorrange=(-4.0,4.0),colormap=:curl,title="synthetic",metadata=(units="degC",))
    @test isa(plot(SSTdiag("","x",:map_point;map_data=mpl,point=nothing,timeseries=nothing,lon=15.0,lat=25.0)),Figure)
    @test isa(plot(SSTdiag("","x",:map_point;map_data=mpl,point=nothing,timeseries=nothing,lon=15.0,lat=25.0,show_anom=true)),Figure)

    SP=Base.get_extension(Climatology,:ClimatologyMakieExt).SST_plots
    @test SP.decimal_years([Date(2024,2,29)])[1]==2024+59/366
    @test SP.decimal_years([Date(2023,1,1),Date(2023,12,31)])==[2023.0,2023+364/365]

    ## B. downloaded daily files

    if have_files
        fil=list_downloaded[end]

        lon,lat=SST_FILES.read_lon_lat(fil)
        @test isa(lon,Vector)

        gr=SST_coarse_grain.grid(fil)
        @test isa(gr,NamedTuple)

        (df,gdf,kdf)=SST_coarse_grain.lowres_read(path=input_path)
        kdf0=kdf[SST_coarse_grain.lowres_index(205,25,kdf)]
        (lon1,lat1)=SST_coarse_grain.lowres_position(kdf0.i,kdf0.j,kdf)
        #ts=SST_timeseries.calc(kdf0,whole_list,gdf=gdf)

        @test isa(df,SST_FILES.DataFrame)
    end

    ## C. coarse-grained table from Zenodo (DataDep OISST-stats1)

    if have_stats
        path_OISST_stats=Climatology.downloads.OISST_stats_download()

        dlon=10.0
        dnl=Int(dlon/0.25)
        (df,gdf,kdf)=SST_coarse_grain.lowres_read(fil="lowres_oisst_sst_$(dlon).csv",path=path_OISST_stats)

        lon0=205; lat0=25
        list=SST_FILES.read_files_list(path=input_path)[1:length(unique(df.t)),:]

        kdf0=kdf[SST_coarse_grain.lowres_index(lon0,lat0,kdf)]
        (lon1,lat1)=SST_coarse_grain.lowres_position(kdf0.i,kdf0.j,kdf)
        ts=SST_timeseries.calc(kdf0,list,gdf=gdf)

        f1=plot(SSTdiag(input_path,"",:by_year; timeseries=ts))
        @test isa(f1,Figure)

        f2=plot(SSTdiag(input_path,"",:by_time; timeseries=ts,show_anom=false,show_clim=false))
        @test isa(f2,Figure)

        f3=plot(SSTdiag(input_path,"",:MHW; timeseries=ts))
        @test isa(f3,Figure)

        ## point lookup / series on the coarse-grained table (SST_coarse_grain.point_series)
        c=SST_coarse_grain.lowres_cell(205,45,gdf)
        @test (c.i,c.j)==(21,14)
        @test SST_coarse_grain.lowres_cell(-155,45,gdf)==c    # longitude wraps
        @test isnothing(SST_coarse_grain.lowres_cell(15,25,gdf))   # Sahara: no stored cell
        @test isnothing(SST_coarse_grain.lowres_cell(0,95,gdf))
        ps=SST_coarse_grain.point_series(205,45,gdf)
        @test ps.cell==c && (ps.lon,ps.lat)==(205.0,45.0) && ps.time_axis==:multi_year
        @test length(ps.time)==length(ps.sst) && all(diff(ps.time).==Day(1))
        @test isnothing(SST_coarse_grain.point_series(15,25,gdf))

        ## map_point: map plus 10° box series, ocean (both anomaly settings) and land
        mpm=(field=Union{Missing,Float64}[1.0 missing; -1.0 0.5],lon=[205.0,206.0],lat=[44.0,45.0],
             colorrange=(-4.0,4.0),colormap=:curl,title="synthetic",metadata=(units="degC",))
        tsp=SST_timeseries.calc(ps.cell,list,gdf=gdf)
        @test isa(plot(SSTdiag("",mpm.title,:map_point;map_data=mpm,point=ps,timeseries=tsp,lon=205.0,lat=45.0,show_anom=true)),Figure)
        @test isa(plot(SSTdiag("",mpm.title,:map_point;map_data=mpm,point=ps,timeseries=tsp,lon=205.0,lat=45.0,show_anom=false)),Figure)

        ## table checks (Zenodo table)
        info=SST_coarse_grain.table_info(df)
        @test (info.t0,info.t1,info.ndays)==(Date(1982,1,1),Date(2024,9,6),15590)
        @test info.gap_free && info.ncells==546==length(kdf)
        @test isnothing(info.list_rows) && isnothing(info.list_ok) && isnothing(info.lag_days)
        info=SST_coarse_grain.table_info(df;list=list,newest=Date(2024,9,10))
        @test info.list_rows==size(list,1) && info.list_ok==(size(list,1)>=info.ndays)
        @test info.lag_days==4
    end

    ## D. needs the downloaded files (grid, map file) and the Zenodo table

    if have_files && have_stats
        gdf1=SST_FILES.groupby(df, :t)
        tmp1=gdf1[end]
        area_tmp=[SST_coarse_grain.areaintegral(gr.msk,x.i,x.j,gr,dnl) for x in eachrow(tmp1)]
        glmsst=[sum(tmp1.sst[:].*area_tmp)/sum(area_tmp) for tmp1 in gdf1]
        ts_global=SST_timeseries.calc(glmsst,list,title="Global Mean SST")
        @test isequal(SST_timeseries.global_mean(df,gr,list;dnl=dnl),ts_global)

        x=SSTdiag(input_path,"local and global SST anomalies",:local_and_global;
                    timeseries=ts,timeseries_global=ts_global)
        f4=plot(x)
        @test isa(f4,Figure)

        ##

        file_climatology=joinpath(path_OISST_stats,"OISST_mean_monthly_1992_2011.nc")
        to_map=(field=SST_FILES.read_map(variable="anom",file=fil,file_climatology=file_climatology),
                colorrange=4 .*(-1.0,1.0),colormap=:thermal,
                lon=gr.lon,lat=gr.lat,lon1=lon1,lat1=lat1,showgrid=false)
        f7=plot(SSTdiag(input_path,"test",:map; map_data=to_map))
        @test isa(f7,Figure)

        ## figure-free map data (SST_FILES.map_data) agrees with read_map / grid
        for v in ("sst","anom","anom_recompute")
            md=SST_FILES.map_data(variable=v,file=fil,file_climatology=file_climatology)
            @test isequal(md.field,SST_FILES.read_map(variable=v,file=fil,file_climatology=file_climatology))
            @test md.lon==gr.lon && md.lat==gr.lat
            @test size(md.field)==(length(md.lon),length(md.lat))
            @test isa(plot(SSTdiag(input_path,md.title,:map; map_data=merge(md,(lon1=lon1,lat1=lat1,showgrid=false)))),Figure)
        end
        @test_throws ErrorException SST_FILES.map_data(variable="bogus",file=fil)
        @test_throws ErrorException SST_FILES.map_data(variable="anom_recompute",file=fil)

        ##

        zm=SST_coarse_grain.calc_zm(gr,df,dnl)

        za=SST_timeseries.zonal_anomaly(gr,df,dnl,list)
        @test size(za)==size(zm) && isequal(isnan.(za),isnan.(zm))
        l=findfirst(!isnan,zm[:,1])
        tmp=SST_timeseries.calc(zm[l,:],list,title="Zonal Mean SST")
        @test za[l,:]==tmp.sst-tmp.clim
        f5=plot(SSTdiag(input_path,"OISST anomaly",:TimeLat; timeseries=list,zonal_mean=zm))
        @test isa(f5,Figure)
    end

    ## table_info on synthetic data (no files needed)

    t=[20240228,20240229,20240301]                                  # leap day, month change
    dfs=SST_FILES.DataFrame(i=repeat([1,2],inner=3),j=1,t=repeat(t,2),sst=1.0)
    i=SST_coarse_grain.table_info(dfs;list=SST_FILES.DataFrame(x=1:2),newest=Date(2024,3,5))
    @test i.gap_free && i.ndays==3 && i.ncells==2
    @test i.list_rows==2 && i.list_ok==false && i.lag_days==4
    @test !SST_coarse_grain.table_info(dfs[dfs.t.!=20240229,:]).gap_free

end

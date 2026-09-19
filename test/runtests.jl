using Test, Climatology, Statistics, MITgcm, CairoMakie, Suppressor
import NCDatasets, NetCDF, MeshArrays
import Climatology.Dates: DateTime, Month

ENV["DATADEPS_ALWAYS_ACCEPT"]=true

#use this environment variable to bypass downloads / API calls that require server access
_SKIP_DOWNLOADS = parse(Bool,get(ENV, "SKIP_DOWNLOADS", "false"))

p=dirname(pathof(Climatology))

@testset "NCDatasetsExt" begin
    G=MeshArrays.Grids_simple.GridLoad_lonlatdep(collect(1000:1000:5000),ones(360,180,5))
    M=MeshArrays.Integration.define_sums(grid=G,regions=(10,5))
    H=ones(length(M.names),length(M.depths),3)
    V=MeshArrays.Integration.volumes(M,G)
    Climatology.write_H_to_T(tempname()*".nc",M,G,H,V)
end

@testset "Time Series Analysis" begin
    dates = collect((DateTime(2020,1,16):Month(1):DateTime(2024,12,31)));
    t_years = Climatology.datetime_to_years(dates);
    dates_reconstructed = Climatology.years_to_datetime(t_years)

    nt=length(dates); data = randn(nt) .+ sin.(2π .* (1:nt) ./ 12)
    data_smc = simple_monthly_climatology(dates, data)
    data_fit = fit_time_series(dates, data, order_season=3)

    @test isa(data_fit,Vector)
end

@testset "Climatology.jl" begin

    ## 1. SST

    input_path=tempname()

    list_downloaded=SST_processing.download_files(path=input_path,short_demo=true)
    @test ispath(input_path)

    output_path=SST_processing.coarse_grain(path=input_path,short_demo=true)
    @test isfile(output_path)

    #@everywhere using Climatology, NCDatasets
    #output_path=SST_processing.monthly_climatology(path=input_path)
    #mv(output_file,joinpath(input_path,basename(output_file)))

    (fil1,fil2)=SST_FILES.file_lists(path=input_path)
    whole_list=SST_FILES.CSV.read(fil1,SST_FILES.DataFrame)
    fil=list_downloaded[end]

    lon,lat=SST_FILES.read_lon_lat(fil)
    @test isa(lon,Vector)

    gr=SST_coarse_grain.grid(fil)
    @test isa(gr,NamedTuple)

    list_pb=SST_FILES.test_files(whole_list)
    @test isa(list_pb,Vector)

    (fil1,fil2)=SST_FILES.ersst_file_lists(path=input_path)
    @test isfile(fil1)

    (df,gdf,kdf)=SST_coarse_grain.lowres_read(path=input_path)
	kdf0=kdf[SST_coarse_grain.lowres_index(205,25,kdf)]
    (lon1,lat1)=SST_coarse_grain.lowres_position(kdf0.i,kdf0.j,kdf)
    #ts=SST_timeseries.calc(kdf0,whole_list,gdf=gdf)

    @test isa(df,SST_FILES.DataFrame)

    ###

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

	gdf1=SST_FILES.groupby(df, :t)
	tmp1=gdf1[end]
	area_tmp=[SST_coarse_grain.areaintegral(gr.msk,x.i,x.j,gr,dnl) for x in eachrow(tmp1)]
	glmsst=[sum(tmp1.sst[:].*area_tmp)/sum(area_tmp) for tmp1 in gdf1]
    ts_global=SST_timeseries.calc(glmsst,list,title="Global Mean SST")

    x=SSTdiag(input_path,"local and global SST anomalies",:local_and_global;
                timeseries=ts,timeseries_global=ts_global)
    f4=plot(x)
    @test isa(f4,Figure)

    ##

    path_OISST_stats=Climatology.downloads.OISST_stats_download()
    file_climatology=joinpath(path_OISST_stats,"OISST_mean_monthly_1992_2011.nc")
    to_map=(field=SST_FILES.read_map(variable="anom",file=fil,file_climatology=file_climatology),
            colorrange=4 .*(-1.0,1.0),colormap=:thermal,
            lon=gr.lon,lat=gr.lat,lon1=lon1,lat1=lat1,showgrid=false)
    f7=plot(SSTdiag(input_path,"test",:map; map_data=to_map))
    @test isa(f7,Figure)

    ##

    zm=SST_coarse_grain.calc_zm(gr,df,dnl)
    f5=plot(SSTdiag(input_path,"OISST anomaly",:TimeLat; timeseries=list,zonal_mean=zm))
    @test isa(f5,Figure)

    ## 2. ECCO

    γ=MeshArrays.GridSpec(ID=:LLC90)
    _SKIP_DOWNLOADS ? nothing : Climatology.get_ecco_files(γ,"oceQnet")
    tmp=read_nctiles(joinpath(ScratchSpaces.ECCO,"oceQnet/oceQnet"),"oceQnet",γ,I=(:,:,1))
    
    tmp=[mean(tmp[j][findall((!isnan).(tmp[j]))]) for j=1:5]
    ref=[19.88214831145215,47.63055475475805,-44.1122401210416,
         3.4402271721659816,30.14270126344508]
    @test tmp==ref

    _SKIP_DOWNLOADS ? nothing : get_occa_velocity_if_needed()
    _SKIP_DOWNLOADS ? nothing : get_occa_variable_if_needed("DDuvel")
    @test isfile(joinpath(ScratchSpaces.OCCA,"DDuvel.0406clim.nc"))

    _SKIP_DOWNLOADS ? nothing : get_ecco_velocity_if_needed()
    _SKIP_DOWNLOADS ? nothing : get_ecco_variable_if_needed("UVELMASS")
    @test isdir(joinpath(ScratchSpaces.ECCO,"UVELMASS"))

    ##

    if !_SKIP_DOWNLOADS
    Climatology.MITPROFclim_download()
#   Climatology.CBIOMESclim_download()
    Climatology.ECCOdiags_add("release2")
    @test true 
    end

    ##

    if true
        var_list3d=("THETA","SALT","UVELMASS","VVELMASS",
        "ADVx_TH","ADVy_TH","DFxE_TH","DFyE_TH")
        var_list2d=("MXLDEPTH","SIarea","sIceLoad","ETAN")
        _SKIP_DOWNLOADS ? nothing : [get_ecco_variable_if_needed(v) for v in var_list3d]
        _SKIP_DOWNLOADS ? nothing : [get_ecco_variable_if_needed(v) for v in var_list2d]
    else
        _SKIP_DOWNLOADS ? nothing : get_ecco_variable_if_needed("MXLDEPTH")
    end

    pth=ECCO.standard_analysis_setup(ScratchSpaces.ECCO)
    list0=ECCO_helpers.standard_list_toml("")
    P0=ECCO_helpers.parameters(pth,"r2",list0[4])

    !isdir(dirname(P0.pth_out)) ? mkdir(dirname(P0.pth_out)) : nothing
    pth_trsp=joinpath(pth,P0.sol,"ECCO_transport_lines")
    isdir(pth_trsp) ? mv(pth_trsp,tempname()) : nothing
    ECCO_helpers.transport_lines(P0.Γ,pth_trsp)
    
    for k in [collect(1:8)...,12,13,25,26,27,28]
        P=ECCO_helpers.parameters(P0,list0[k])
        !isdir(P.pth_out) ? mkdir(P.pth_out) : nothing
        ECCO_diagnostics.driver(P)
    end

    fil0=joinpath(P0.pth_out,"zonmean2d.jld2")
    @test isfile(fil0)

    ##

    [ECCO_procs.years_min_max(sol) for sol in ("ECCOv4r3","ECCOv4r4","ECCOv4r5","OCCA2HR1","OCCA2HR2")]

    sol="ECCO4R2"
    year0,year1=ECCO_procs.years_min_max(sol)
    pth_out=Climatology.downloads.ECCOdiags_add(sol)

    ECCOdiags_to_nc(path_in=datadep"ECCO4R2-stdiags",year1=1992,nt=240)

    using CairoMakie

    P=ECCO_procs.parameters()
    nammap=P.clim_longname[11]
    statmap="mean"
    timemap=1

    plot(ECCOdiag(pth_out, nammap, :ECCO_map; P=P, statistic=statmap, time=timemap))
    
    l0=1; l1=90

    plot(ECCOdiag(pth_out, "THETA_clim", :ECCO_TimeLat;
        P=P, period=(year0,year1), level=1, ylims=(l0,l1), colormap_factor=1.0))

    plot(ECCOdiag(pth_out, "THETA_clim", :ECCO_TimeLatAnom;
        P=P, period=(year0,year1), level=1, ylims=(l0,l1), colormap_factor=1.0))

    k0=1; k1=30
    plot(ECCOdiag(pth_out, "THETA_clim", :ECCO_DepthTime;
        P=P, period=(year0,year1), factor=1.0, level=28, klims=(k0,k1)))

    plot(ECCOdiag(pth_out, "THETA", :ECCO_GlobalMean; level=0, period=(year0,year1)))

    plot(ECCOdiag(pth_out, "OHT", :ECCO_OHT1))

    plot(ECCOdiag(pth_out, "overturn", :ECCO_Overturn2; grid=P.Γ))

    plot(ECCOdiag(pth_out, "overturn", :ECCO_Overturn1;
        level=29, low1="auto", period=(year0,year1)))

    ntr1=P.list_trsp[1]
    fig=plot(ECCOdiag(pth_out, "trsp", :ECCO_Transports;
        namtrs=[ntr1], ncols=1, list_trsp=P.list_trsp, period=(year0,year1)))

    @test isa(fig,Figure)
    @test ispath(pth_out)

    ## 3. SSH/SLA

    if !_SKIP_DOWNLOADS
    SLA=read(SeaLevelAnomaly(name="sla_podaac"))
    f3=plot(SLA)
    @test isa(f3,Figure)

    file=joinpath(SLA.path,SLA.name*".nc")
    gr=SLA_PODAAC.get_grid(file=file)
    data=SLA_PODAAC.read_slice(file,gr)
    sub=SLA_PODAAC.subset(; read_from_file=file,save_to_file=true)
    @test isa(sub,String)

    SLA=read(SeaLevelAnomaly(name="sla_cmems"))
    file=joinpath(SLA.path,SLA.name*".nc")
    sub=SLA_CMEMS.subset(; read_from_file=file,save_to_file=true)
    @test isa(sub,String)
    end

end
##

@testset "SurfaceFluxes" begin
    data=ERA5.read_sample()
    sst=ERA5.interpolate_sst(data.OISST.sst,data.tim)
    df=ERA5.surface_balance(data.ERA5,sst)
    @test isa(df,Climatology.DataFrame)

    da=Climatology.SurfaceFluxDiag((plot_type=:default,),df)
    plot(da)
    da=Climatology.SurfaceFluxDiag((plot_type=:surface_balance,),(df=df,tim=data.tim,sst=sst))
    plot(da)
    da=Climatology.SurfaceFluxDiag((plot_type=:Qnet_cumsum,),(df=df,tim=data.tim,sst=sst))
    fi=plot(da)
    @test isa(fi,Figure)
end
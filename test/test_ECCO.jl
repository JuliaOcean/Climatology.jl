# ECCO tests (downloads gated by _SKIP_DOWNLOADS where noted). Included from runtests.jl.
@testset "ECCO" begin

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

    # data-level accessor: no Figure, native field retained, wrapper unchanged
    dmap=ECCO_procs.ECCO_map_data(ECCOdiag(pth_out, nammap, :ECCO_map; P=P, statistic=statmap, time=timemap))
    wmap=ECCO_procs.ECCO_map(ECCOdiag(pth_out, nammap, :ECCO_map; P=P, statistic=statmap, time=timemap))
    @test keys(wmap)==(:λ,:field,:levels,:title)
    @test isequal(dmap.field,wmap.field) && dmap.levels==wmap.levels && dmap.title==wmap.title
    @test isa(dmap.native,MeshArrays.AbstractMeshArray) && dmap.location==:C
    @test size(dmap.field)==size(P.λ.lon)
    @test any(isnan,dmap.field) && any(!isnan,dmap.field)
    @test sum(count(isnan,dmap.native[f]) for f in 1:5)==sum(count(isnan,(P.μ*dmap.native)[f]) for f in 1:5) #accessor adds no masking to native
    dmon=ECCO_procs.ECCO_map_data(ECCOdiag(pth_out, nammap, :ECCO_map; P=P, statistic="mon", time=3))
    @test dmon.metadata.time==3 && dmap.metadata.time===nothing
    @test_throws ErrorException ECCO_procs.ECCO_map_data(ECCOdiag(pth_out, "not a product", :ECCO_map; P=P))

    # lon/lat -> native wet cell
    c=ECCO_procs.nearest_wet_cell(P,-30.0,40.0)
    @test !isnothing(c) && !isnan(P.μ[c.f][c.i,c.j])
    @test abs(c.lon+30.0)<3.0 && abs(c.lat-40.0)<3.0
    @test c==ECCO_procs.nearest_wet_cell(P,330.0,40.0) #longitude wrapping
    # cell -> monthly climatology at that cell
    Xmap=ECCOdiag(pth_out, P.clim_longname[11], :ECCO_map; P=P)
    ps=ECCO_procs.map_point_series(Xmap,c)
    @test length(ps.values)==12 && ps.month==1:12 && all(!isnan,ps.values)
    @test isapprox(sum(ps.values)/12,ps.mean) #stored mean = average of the 12 months
    @test ps.std>0 && (ps.cell.f,ps.cell.i,ps.cell.j)==(c.f,c.i,c.j)
    dnat=ECCO_procs.ECCO_map_data(Xmap)
    @test ps.mean==dnat.native[c.f][c.i,c.j] #consistent with the map accessor
    @test_throws ErrorException ECCO_procs.map_point_series(Xmap,(f=1,i=0,j=1))
    @test_throws ErrorException ECCO_procs.map_point_series(Xmap,(f=9,i=1,j=1))
    iland=findfirst(isnan,P.μ[1]) #first land cell of face 1
    @test_throws ErrorException ECCO_procs.map_point_series(Xmap,(f=1,i=iland[1],j=iland[2]))
    @test isnothing(ECCO_procs.nearest_wet_cell(P,10.0,25.0)) #Sahara: all four neighbors are land
    @test isnothing(ECCO_procs.nearest_wet_cell(P,0.0,-89.75)) #Antarctica

    # :ECCO_map_point plot type (data + figure)
    Xpt=ECCOdiag(pth_out, P.clim_longname[11], :ECCO_map_point; P=P, lon=-30.0, lat=40.0)
    rpt=ECCO_procs.ECCO_map_point(Xpt)
    @test rpt.cell==c && rpt.series.values==ps.values && rpt.requested==(lon=-30.0,lat=40.0)
    @test isa(plot(Xpt),Figure)
    Xland=ECCOdiag(pth_out, P.clim_longname[11], :ECCO_map_point; P=P, lon=10.0, lat=25.0)
    rland=ECCO_procs.ECCO_map_point(Xland)
    @test isnothing(rland.cell) && isnothing(rland.series) #land: no error
    @test isa(plot(Xland),Figure)

    # request resolution and diagnostic constructors (same objects as the standard notebook)
    @test ECCO_procs.years_min_max("ECCO4R3")==ECCO_procs.years_min_max("ECCOv4r3")==(1992,2015)
    @test ECCO_procs.years_min_max("ECCO4R5")==(1992,2019) && ECCO_procs.years_min_max("ECCO4R2")==(1992,2011)
    @test ECCO_procs.available_period("custom")==(1992,2027)
    @test_throws ErrorException ECCO_procs.available_period("nope")
    @test_throws ErrorException ECCO_procs.diagnostic_path("custom")
    @test_throws ErrorException ECCO_procs.diagnostic_path("custom"; custom_path=joinpath(tempdir(),"does_not_exist"))
    @test ECCO_procs.diagnostic_path("custom"; custom_path=tempdir())==tempdir()

    req=(product=:ECCO, solution="ECCO4R2", period=(1992,2011))
    d=ECCO_procs.resolve_request(req; P=P)
    @test d.pth_out==pth_out && d.period==(1992,2011) && d.grid===P.Γ
    @test d.metadata.source==:Climatology_DataDeps && d.metadata.available_period==(1992,2011)
    @test d.metadata.available_map_products==P.clim_longname
    @test_throws ErrorException ECCO_procs.resolve_request((product=:ECCO, solution="ECCO4R2", period=(1992,2020)); P=P)
    @test_throws ErrorException ECCO_procs.resolve_request((product=:ECCO, solution="ECCO4R2", period=(2000,1999)); P=P)
    @test ECCO_procs.resolve_request((product=:ECCO, solution="custom", period=(1992,2000)); P=P, custom_path=tempdir()).metadata.source==:local

    same(a,b)=(a.path==b.path && a.name==b.name && isequal(a.options,b.options)) #isequal: P.μ has NaN (land)
    @test same(ECCO_procs.map_diagnostic(d; name=P.clim_longname[11], statistic="mon", month=3),
        ECCOdiag(pth_out, P.clim_longname[11], :ECCO_map; P=P, statistic="mon", time=3))
    @test_throws ErrorException ECCO_procs.map_diagnostic(d; name="x")
    @test_throws ErrorException ECCO_procs.map_diagnostic(d; name=P.clim_longname[1], statistic="bad")
    @test_throws ErrorException ECCO_procs.map_diagnostic(d; name=P.clim_longname[1], month=13)
    @test same(ECCO_procs.timelat_diagnostic(d; variable="THETA", level=2, ylims=(-50,50), colormap_factor=0.5),
        ECCOdiag(pth_out, "THETA_clim", :ECCO_TimeLat; P=P, period=(1992,2011), level=2, ylims=(-50,50), colormap_factor=0.5))
    @test same(ECCO_procs.timelat_anomaly_diagnostic(d; variable="THETA"),
        ECCOdiag(pth_out, "THETA_clim", :ECCO_TimeLatAnom; P=P, period=(1992,2011), level=1, ylims=(-90.0,90.0), colormap_factor=1.0))
    @test same(ECCO_procs.depthtime_diagnostic(d; variable="THETA", latitude_index=28),
        ECCOdiag(pth_out, "THETA_clim", :ECCO_DepthTime; P=P, period=(1992,2011), factor=1.0, level=28, klims=(1,30)))
    @test same(ECCO_procs.global_mean_diagnostic(d; variable="THETA", level=0),
        ECCOdiag(pth_out, "THETA", :ECCO_GlobalMean; level=0, period=(1992,2011)))
    @test same(ECCO_procs.overturning_mean_diagnostic(d), ECCOdiag(pth_out, "overturn", :ECCO_Overturn2; grid=P.Γ))
    @test same(ECCO_procs.overturning_time_diagnostic(d; level=29),
        ECCOdiag(pth_out, "overturn", :ECCO_Overturn1; level=29, low1="auto", period=(1992,2011)))
    @test same(ECCO_procs.transports_diagnostic(d; names=P.list_trsp[1:2]),
        ECCOdiag(pth_out, "trsp", :ECCO_Transports; namtrs=P.list_trsp[1:2], ncols=1, list_trsp=P.list_trsp, period=(1992,2011)))
    @test_throws ErrorException ECCO_procs.transports_diagnostic(d; names=String[])
    @test_throws ErrorException ECCO_procs.transports_diagnostic(d; names=["not a section"])

    # parameters(path=...): same metadata for a standard release, empty transports for a bare directory
    Pp=ECCO_procs.parameters(path=pth_out)
    @test Pp.clim_longname==P.clim_longname && Pp.list_trsp==P.list_trsp
    bare=mktempdir(); mkpath(joinpath(bare,"THETA_clim")); touch(joinpath(bare,"THETA_clim","THETA_k01.jld2"))
    Pb=ECCO_procs.parameters(path=bare)
    @test isempty(Pb.list_trsp) && length(Pb.clim_name)==1

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

end

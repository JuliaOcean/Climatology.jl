using Test, Climatology, Statistics, MITgcm
import NCDatasets, NetCDF, MeshArrays

ENV["DATADEPS_ALWAYS_ACCEPT"]=true

p=dirname(pathof(Climatology))

@testset "Climatology.jl" begin
    γ=MeshArrays.GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
    Climatology.get_ecco_files(γ,"oceQnet")
    tmp=read_nctiles(joinpath(ScratchSpaces.ECCO,"oceQnet/oceQnet"),"oceQnet",γ,I=(:,:,1))
    
    tmp=[mean(tmp[j][findall((!isnan).(tmp[j]))]) for j=1:5]
    ref=[19.88214831145215,47.63055475475805,-44.1122401210416,
         3.4402271721659816,30.14270126344508]
    @test tmp==ref

    get_occa_velocity_if_needed()
    get_occa_variable_if_needed("DDuvel")
    @test isfile(joinpath(ScratchSpaces.OCCA,"DDuvel.0406clim.nc"))

    get_ecco_velocity_if_needed()
    get_ecco_variable_if_needed("UVELMASS")
    @test isdir(joinpath(ScratchSpaces.ECCO,"UVELMASS"))

    ##

    Climatology.MITPROFclim_download()
    Climatology.CBIOMESclim_download()
    Climatology.ECCOdiags_add("release2")
    @test true 

    ##

    if true
        var_list3d=("THETA","SALT","UVELMASS","VVELMASS",
        "ADVx_TH","ADVy_TH","DFxE_TH","DFyE_TH")
        var_list2d=("MXLDEPTH","SIarea","sIceLoad","ETAN")
        [get_ecco_variable_if_needed(v) for v in var_list3d]
        [get_ecco_variable_if_needed(v) for v in var_list2d]
    else
        get_ecco_variable_if_needed("MXLDEPTH") 
    end

    MeshArrays.GRID_LLC90_download()
    pth=ECCO.standard_analysis_setup(ScratchSpaces.ECCO)
    list0=ECCO_helpers.standard_list_toml("")
    P0=ECCO_helpers.parameters(pth,"r2",list0[4])

    !isdir(dirname(P0.pth_out)) ? mkdir(dirname(P0.pth_out)) : nothing
    pth_trsp=joinpath(pth,P0.sol,"ECCO_transport_lines")
    !isdir(pth_trsp) ? ECCO_helpers.transport_lines(P0.Γ,pth_trsp) : nothing
    
    for k in [collect(1:8)...,12,13,25,26,27,28]
        P=ECCO_helpers.parameters(P0,list0[k])
        !isdir(P.pth_out) ? mkdir(P.pth_out) : nothing
        ECCO_diagnostics.driver(P)
    end

    fil0=joinpath(P0.pth_out,"zonmean2d.jld2")
    @test isfile(fil0)

    ##

    sol="ECCO4R2"
    year0,year1=ECCO_procs.years_min_max(sol)
    pth_out=Climatology.downloads.ECCOdiags_add(sol)

    ECCOdiags_to_nc(path_in=datadep"ECCO4R2-stdiags",year1=1992,nt=240)

    using CairoMakie

    P=ECCO_procs.parameters()
    nammap=P.clim_longname[11]
    statmap="mean"
    timemap=1
    plot(ECCOdiag(path=pth_out,name="tbd",options=
        (plot_type=:ECCO_map,nammap=nammap,P=P,statmap=statmap,timemap=timemap)))
    
    plot(ECCOdiag(path=pth_out,name="THETA_clim",options=
        (plot_type=:ECCO_TimeLat,year0=year0,year1=year1,cmap_fac=1.0,
        k=1,P=P,years_to_display=[year0 year1+1])))

    l0=1; l1=90
    plot(ECCOdiag(path=pth_out,name="THETA_clim",options=
        (plot_type=:ECCO_TimeLatAnom,year0=year0,year1=year1,cmap_fac=1.0,
        k=1,l0=l0,l1=l1,P=P,years_to_display=[year0 year1+1])))

    k0=1; k1=30
    plot(ECCOdiag(path=pth_out,name="THETA_clim",options=
        (plot_type=:ECCO_DepthTime,facA=1.0,l=28,year0=year0,year1=year1,
        k0=k0,k1=k1,P=P,years_to_display=[year0 year1+1])))

    plot(ECCOdiag(path=pth_out,name="THETA",options=
        (plot_type=:ECCO_GlobalMean,k=0,year0=year0,year1=year1,
        years_to_display=[year0 year1+1])))

    plot(ECCOdiag(path=pth_out,name="OHT",options=(plot_type=:ECCO_OHT1,)))

    plot(ECCOdiag(path=pth_out,name="overturn",options=(plot_type=:ECCO_Overturn2,grid=P.Γ)))

    plot(ECCOdiag(path=pth_out,name="overturn",options=
    (plot_type=:ECCO_Overturn1,kk=29,low1="auto",year0=year0,year1=year1,
    years_to_display=[year0 year1+1])))

    ntr1=P.list_trsp[1]
    plot(ECCOdiag(path=pth_out,name="trsp",options=
    (plot_type=:ECCO_Transports,namtrs=[ntr1],ncols=1,list_trsp=P.list_trsp,
    year0=year0,year1=year1,years_to_display=[year0 year1+1])))

    @test ispath(pth_out)

end

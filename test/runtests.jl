using OceanStateEstimation, MeshArrays, Statistics
using Test

using Dataverse
pyDataverse.APIs(do_install=true)

p=dirname(pathof(OceanStateEstimation))

@testset "OceanStateEstimation.jl" begin
    γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
    tmp=OceanStateEstimation.get_ecco_files(γ,"oceQnet")
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

    OceanStateEstimation.MITPROFclim_download()
    OceanStateEstimation.CBIOMESclim_download()
    OceanStateEstimation.ECCOdiags_add("release2")
    OceanStateEstimation.ECCOdiags_add("interp_coeffs")
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

end

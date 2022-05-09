using OceanStateEstimation, MeshArrays, Statistics
using Test

p=dirname(pathof(OceanStateEstimation))

@testset "OceanStateEstimation.jl" begin
    γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
    tmp=OceanStateEstimation.get_ecco_files(γ,"oceQnet")
    tmp=[mean(tmp[j][findall((!isnan).(tmp[j]))]) for j=1:5]
    ref=[19.88214831145215,47.63055475475805,-44.1122401210416,
         3.4402271721659816,30.14270126344508]
    @test tmp==ref

    lst=joinpath(p,"../examples/nctiles_climatology.csv")
    lists=dataverse_lists(lst)

    get_occa_velocity_if_needed()
    get_occa_variable_if_needed("DDuvel")
    @test isfile(joinpath(OCCAclim_path,"DDuvel.0406clim.nc"))

    get_ecco_velocity_if_needed()
    get_ecco_variable_if_needed("UVELMASS")
    @test isdir(joinpath(ECCOclim_path,"UVELMASS"))

    isdir(MITPROFclim_path)
    isdir(CBIOMESclim_path)
    @test true 

    if false
        var_list3d=("THETA","SALT","UVELMASS","VVELMASS",
        "ADVx_TH","ADVy_TH","DFxE_TH","DFyE_TH")
        var_list2d=("MXLDEPTH","SIarea","sIceLoad","ETAN")
        [get_ecco_variable_if_needed(v) for v in var_list3d]
        [get_ecco_variable_if_needed(v) for v in var_list2d]
    else
        get_ecco_variable_if_needed("MXLDEPTH") 
    end

    MeshArrays.GRID_LLC90_download()
    pth=ECCO.standard_analysis_setup(ECCOclim_path)
    list0=ECCO_helpers.standard_list_toml("")
    P=ECCO_helpers.parameters(pth,"r2",list0,4)

    !isdir(dirname(P.pth_out)) ? mkdir(dirname(P.pth_out)) : nothing
    !isdir(P.pth_out) ? mkdir(P.pth_out) : nothing
    ECCO_diagnostics.driver(P)

#    include(joinpath("..","examples","ECCO","ECCO_standard_loop.jl"))
#    Pkg.activate(pth)

    fil0=joinpath(P.pth_out,"zonmean2d.jld2")
    @test isfile(fil0)

end

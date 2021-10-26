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

end

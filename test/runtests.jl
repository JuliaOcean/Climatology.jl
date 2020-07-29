using OceanStateEstimation, MeshArrays, Statistics
using Test

@testset "OceanStateEstimation.jl" begin
    γ=GridSpec("LatLonCap","./")
    tmp=get_ecco_files(γ,"oceQnet")
    tmp=[mean(tmp[j][findall((!isnan).(tmp[j]))]) for j=1:5]
    ref=[19.88214831145215,47.63055475475805,-44.1122401210416,
         3.4402271721659816,30.14270126344508]
    @test tmp==ref
end

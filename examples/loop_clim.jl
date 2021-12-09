
@everywhere sol0="r3"
@everywhere nam="THETA"
@everywhere calc="clim"
#@everywhere list_kk=[1 10 20 29 38 44]
#@everywhere list_kk=[1 10 20]
@everywhere list_kk=[29 38 44]

@everywhere using JLD2

for loc_kk in list_kk
 save("kk.jld2","kk",loc_kk)
 @everywhere kk=load("kk.jld2","kk")
 include("ECCO_standard_analysis.jl")
end


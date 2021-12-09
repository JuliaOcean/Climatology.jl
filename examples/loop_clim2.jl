
@everywhere sol0="r5"
#@everywhere list_nam=["SSH","MXLDEPTH","SIarea"]
@everywhere list_nam=["BSF"]
@everywhere calc="clim"
#@everywhere list_kk=[1 10 20 29 38 44]

@everywhere using JLD2

for loc_nam in list_nam
 save("nam.jld2","nam",loc_nam)
 @everywhere nam=load("nam.jld2","nam")
 include("ECCO_standard_analysis.jl")
end


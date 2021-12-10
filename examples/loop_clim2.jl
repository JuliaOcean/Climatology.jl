
@everywhere include("ECCO_pkg_grid_etc.jl")

@everywhere sol0="r2"
@everywhere list_nam=["BSF","SSH","MXLDEPTH","SIarea"]
@everywhere calc="clim"
@everywhere kk=1

for loc_nam in list_nam
 save("nam.jld2","nam",loc_nam)
 @everywhere nam=load("nam.jld2","nam")

 @everywhere pth_in,pth_out,pth_tmp,sol,nt,list_steps=ECCO_path_etc(sol0,calc,nam)
 !isdir(pth_out) ? mkdir(pth_out) : nothing
 !isdir(pth_tmp) ? mkdir(pth_tmp) : nothing

 include("ECCO_standard_analysis.jl")
end


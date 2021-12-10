
@everywhere include("ECCO_pkg_grid_etc.jl")

@everywhere sol0="r2"
@everywhere nam="THETA"
@everywhere calc="clim"
@everywhere list_kk=[1 10 20 29 38 44]

for loc_kk in list_kk
 save("kk.jld2","kk",loc_kk)
 @everywhere kk=load("kk.jld2","kk")

 @everywhere pth_in,pth_out,pth_tmp,sol,nt,list_steps=ECCO_path_etc(sol0,calc,nam)
 !isdir(pth_out) ? mkdir(pth_out) : nothing
 !isdir(pth_tmp) ? mkdir(pth_tmp) : nothing

 include("ECCO_standard_analysis.jl")
end


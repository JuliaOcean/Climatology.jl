
@everywhere include("ECCO_pkg_grid_etc.jl")
@everywhere include("ECCO_standard_analysis.jl")

@everywhere sol0="r2"
@everywhere list0=TOML.parsefile("ECCO_diags/ECCO_standard_list.toml")

for ff in 1:length(list0["kk"])
    save("ECCO_diags/taskID.jld2","ID",ff)
    
    @sync @everywhere gg=load("ECCO_diags/taskID.jld2","ID")

    @everywhere calc=list0["calc"][gg]
    @everywhere nam=list0["nam"][gg]
    @everywhere kk=list0["kk"][gg]
    
    @everywhere pth_in,pth_out,pth_tmp,sol,nt,list_steps=ECCO_path_etc(sol0,calc,nam)
    !isdir(pth_out) ? mkdir(pth_out) : nothing
    !isdir(pth_tmp) ? mkdir(pth_tmp) : nothing

    main_function(calc,sol,nam,kk)

end


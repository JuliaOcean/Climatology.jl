
#tmp=readlines("ECCO_standard_Project.toml")[2:end];
#tmp=[split(i," = ")[1] for i in tmp];
#[Pkg.add(i) for i in tmp]

#import MeshArrays
#MeshArrays.GRID_LLC90_download()

#mkdir("ECCO_diags")

@everywhere include("ECCO_pkg_grid_etc.jl")
@everywhere include("ECCO_standard_analysis.jl")

@everywhere sol0="r2"
@everywhere sol="ECCOv4"*sol0*"_analysis"

pth0=joinpath(pth,sol)
!isdir(pth0) ? mkdir(pth0) : nothing

@everywhere fil=joinpath(pth,sol,"ECCO_standard_list.toml")
!isfile(fil) ? ECCO_standard_list_toml(fil) : nothing
@everywhere list0=TOML.parsefile(joinpath(pth,sol,"ECCO_standard_list.toml"))

@everywhere pth0=joinpath(pth,sol,"ECCO_transport_lines")
!isdir(pth0) ? ECCO_transport_lines(pth0) : nothing
@everywhere list_trsp,msk_trsp,ntr=reload_transport_lines(pth0)

for ff in 1:length(list0["kk"])
    save(joinpath(pth,sol,"taskID.jld2"),"ID",ff)
    
    @sync @everywhere gg=load(joinpath(pth,sol,"taskID.jld2"),"ID")

    @everywhere calc=list0["calc"][gg]
    @everywhere nam=list0["nam"][gg]
    @everywhere kk=list0["kk"][gg]
    
    @everywhere pth_in,pth_out,pth_tmp,nt,list_steps=ECCO_path_etc(pth,sol0,calc,nam)
    !isdir(pth_tmp) ? mkdir(pth_tmp) : nothing

    main_function(calc,sol,nam,kk)

end


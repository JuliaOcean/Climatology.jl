
using Distributed, SharedArrays, JLD2

#@everywhere include("ECCO_pkg_grid_etc.jl")
#@everywhere include("ECCO_standard_analysis.jl")

@everywhere begin
    using OceanStateEstimation, Printf, TOML, MeshArrays

    γ,Γ,LC=ECCO_helper_functions.GridLoad_Main()
    nr=length(Γ.DRF)
    nl=length(LC)

    sol0="r2"
    sol="ECCOv4"*sol0*"_analysis" 
    
    fil=joinpath(pth,sol,"ECCO_standard_list.toml")
    pth_trsp=joinpath(pth,sol,"ECCO_transport_lines")
end

pth0=joinpath(pth,sol)
!isdir(pth0) ? mkdir(pth0) : nothing

!isfile(fil) ? ECCO_helper_functions.standard_list_toml(fil) : nothing
@everywhere list0=TOML.parsefile(joinpath(pth,sol,"ECCO_standard_list.toml"))

!isdir(pth_trsp) ? ECCO_helper_functions.transport_lines(pth_trsp) : nothing
@everywhere list_trsp,msk_trsp,ntr=ECCO_helper_functions.reload_transport_lines(pth_trsp)

list1=collect(1:length(list0["kk"]))
#zz=collect(1:6)
#zz=[7,8,12,13]
#zz=[25,26,27,28]

for ff in zz
    save(joinpath(pth,sol,"taskID.jld2"),"ID",ff)
    
    @sync @everywhere gg=load(joinpath(pth,sol,"taskID.jld2"),"ID")

    @everywhere calc=list0["calc"][gg]
    @everywhere nam=list0["nam"][gg]
    @everywhere kk=list0["kk"][gg]
    
    @everywhere pth_in,pth_out,pth_tmp,nt,list_steps=ECCO_helper_functions.path_etc(pth,sol0,calc,nam)
    !isdir(pth_tmp) ? mkdir(pth_tmp) : nothing

    ECCO_diagnostics.main_function(calc,sol,nam,kk; pth_in=pth_in, pth_out=pth_tmp)

end


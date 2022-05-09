
using Distributed, SharedArrays, JLD2

#@everywhere include("ECCO_pkg_grid_etc.jl")
#@everywhere include("ECCO_standard_analysis.jl")

@everywhere begin
    using OceanStateEstimation, Printf, TOML, MeshArrays

    #γ,Γ,LC=ECCO_helper_functions.GridLoad_Main()

    sol0="r2"
    sol="ECCOv4"*sol0*"_analysis" 
    
    fil=joinpath(pth,sol,"ECCO_standard_list.toml")
    pth_trsp=joinpath(pth,sol,"ECCO_transport_lines")
end

!isdir(joinpath(pth,sol)) ? mkdir(joinpath(pth,sol)) : nothing

!isfile(fil) ? ECCO_helper_functions.standard_list_toml(fil) : nothing
list0=TOML.parsefile(joinpath(pth,sol,"ECCO_standard_list.toml"))

!isdir(pth_trsp) ? ECCO_helper_functions.transport_lines(pth_trsp) : nothing
#@everywhere list_trsp,msk_trsp,ntr=ECCO_helper_functions.reload_transport_lines(pth_trsp)

list1=collect(1:length(list0["kk"]))
#list1=collect(1:6)
#list1=[7,8,12,13]
#list1=[25,26,27,28]

for ff in list1
    save(joinpath(pth,sol,"taskID.jld2"),"ID",ff)
    
    @sync @everywhere gg=load(joinpath(pth,sol,"taskID.jld2"),"ID")
    @everywhere P=ECCO_helper_functions.parameters(pth,sol0,list0,gg)            
    
    !isdir(P.pth_out) ? mkdir(P.pth_out) : nothing

    ECCO_diagnostics.main_function(P)

end


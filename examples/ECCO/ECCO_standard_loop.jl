using Distributed, OceanStateEstimation

@everywhere begin
    using Pkg; Pkg.activate("./")
    
    using OceanStateEstimation, Printf, TOML, MeshArrays
    using Distributed, SharedArrays, JLD2

    pth=joinpath(tempdir(),"ECCO_diags_dev")

    sol0="r2"
    sol="ECCOv4"*sol0*"_analysis" 
    
    fil=joinpath(pth,sol,"ECCO_standard_list.toml")
    pth_trsp=joinpath(pth,sol,"ECCO_transport_lines")
end

!isdir(joinpath(pth,sol)) ? mkdir(joinpath(pth,sol)) : nothing

!isfile(fil) ? ECCO_helper_functions.standard_list_toml(fil) : nothing
@everywhere list0=TOML.parsefile(joinpath(pth,sol,"ECCO_standard_list.toml"))

P0=ECCO_helper_functions.parameters(pth,sol0,list0,1)
!isdir(pth_trsp) ? ECCO_helper_functions.transport_lines(P0.Γ,pth_trsp) : nothing

list1=collect(1:length(list0["kk"]))
#list1=collect(1:6)
#list1=[7,8,12,13]
#list1=[25,26,27,28]

for ff in list1
#    save(joinpath(pth,sol,"taskID.jld2"),"ID",ff)
#    @sync @everywhere gg=load(joinpath(pth,sol,"taskID.jld2"),"ID")
    P=ECCO_helper_functions.parameters(pth,sol0,list0,ff)
    !isdir(P.pth_out) ? mkdir(P.pth_out) : nothing
    ECCO_diagnostics.main_function(P)
end


using Distributed, OceanStateEstimation

@everywhere begin
    using Pkg, OceanStateEstimation
    pth0=pwd()
    pth=ECCO.standard_analysis_setup(pth0)
    Pkg.activate(pth)

    sol0="r2"
    list0=ECCO_helpers.standard_list_toml("")
    P0=ECCO_helpers.parameters(pth,sol0,list0[1])
end

!isdir(joinpath(pth,P0.sol)) ? mkdir(joinpath(pth,P0.sol)) : nothing
pth_trsp=joinpath(pth,P0.sol,"ECCO_transport_lines")
!isdir(pth_trsp) ? ECCO_helpers.transport_lines(P0.Γ,pth_trsp) : nothing

#list1=collect(1:length(list0))
list1=collect(3:5)
#list1=[7,8,12,13]
#list1=[25,26,27,28]

for ff in list1
#    save(joinpath(pth,sol,"taskID.jld2"),"ID",ff)
#    @sync @everywhere gg=load(joinpath(pth,sol,"taskID.jld2"),"ID")
#
#    P=ECCO_helpers.parameters(pth,sol0,list0[ff])
    P=ECCO_helpers.parameters(P0,list0[ff])
    !isdir(P.pth_out) ? mkdir(P.pth_out) : nothing
    println("starting calc,sol,nam=$(P.calc),$(P.sol),$(P.nam) ...")

    ECCO_diagnostics.driver(P)
end


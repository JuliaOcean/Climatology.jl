using Distributed, Climatology
import NCDatasets, MITgcm, NetCDF

@everywhere begin
    using Pkg, Climatology
    import NCDatasets, MITgcm, NetCDF
    pth0=ScratchSpaces.ECCO
    #pth0=ENV["custom_path"]
    pth=ECCO.standard_analysis_setup(pth0)
    Pkg.activate(pth)

    sol0="r2"
    #sol0="custom_solution"
    list0=ECCO_helpers.standard_list_toml("")
    P0=ECCO_helpers.parameters(pth,sol0,list0[1])
end

!isdir(joinpath(pth,P0.sol)) ? mkdir(joinpath(pth,P0.sol)) : nothing
pth_trsp=joinpath(pth,P0.sol,"ECCO_transport_lines")
!isdir(pth_trsp) ? ECCO_helpers.transport_lines(P0.Γ,pth_trsp) : nothing

#list1=collect(1:length(list0))
list1=collect(1:5)

for ff in list1
    P=ECCO_helpers.parameters(P0,list0[ff])
    !isdir(P.pth_out) ? mkdir(P.pth_out) : nothing
    println("starting calc,sol,nam=$(P.calc),$(P.sol),$(P.nam) ...")

    ECCO_diagnostics.driver(P)
end


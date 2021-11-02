
#using OceanStateEstimation, MeshArrays, NCTiles
using JLD2
import CairoMakie as Mkie

pth_output=joinpath(tempdir(),"OptimalTransport_example")

function EMD_plot(fil)
    d=load(fil)["d"];
    d[findall(d.==0.0)].=NaN;
    fig = Mkie.Figure(resolution = (600,400), backgroundcolor = :grey95, fontsize=12)
    ax = Mkie.Axis(fig[1,1])
    hm=Mkie.heatmap!(d)
    Mkie.Colorbar(fig[1,2], hm, height = Mkie.Relative(0.65))
    fig
end

function EMD_plot_all(pth=pth_output)
    fil1=joinpath(pth,"ModToMod.jld2")
    fil2=joinpath(pth,"SatToSat.jld2")
    fil3=joinpath(pth,"ModToSat.jld2")

    d1=load(fil1)["d"]; d1[findall(d1.==0.0)].=NaN;
    d2=load(fil2)["d"]; d2[findall(d2.==0.0)].=NaN;
    d3=load(fil3)["d"]; d3[findall(d3.==0.0)].=NaN;
    #just to check the alignment of dimensions d3[1:end,1].=NaN

    #cr=(0.07, 0.15)
    cr=(0.0, 10.0)

    fig = Mkie.Figure(resolution = (600,400), backgroundcolor = :grey95, fontsize=12)
    ax = Mkie.Axis(fig[1,1])
    hm=Mkie.heatmap!(d1, colorrange = cr, colormap=:inferno)
    Mkie.ylims!(ax, (12.5, 0.5)); Mkie.xlims!(ax, (0.5,12.5))  
    ax = Mkie.Axis(fig[1,2])
    hm=Mkie.heatmap!(transpose(d3), colorrange = cr, colormap=:inferno)
    Mkie.ylims!(ax, (12.5, 0.5)); Mkie.xlims!(ax, (0.5,12.5))  
    ax = Mkie.Axis(fig[2,1])
    hm=Mkie.heatmap!(d3, colorrange = cr, colormap=:inferno)
    Mkie.ylims!(ax, (12.5, 0.5)); Mkie.xlims!(ax, (0.5,12.5))  
    ax = Mkie.Axis(fig[2,2])
    hm=Mkie.heatmap!(d2, colorrange = cr, colormap=:inferno)
    Mkie.ylims!(ax, (12.5, 0.5)); Mkie.xlims!(ax, (0.5,12.5))  
    Mkie.Colorbar(fig[1:2,3], hm, height = Mkie.Relative(0.65))
    fig
end


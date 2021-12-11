using Documenter, OceanStateEstimation, Pkg
import PlutoSliderServer, CairoMakie
Pkg.precompile()

makedocs(;
    modules=[OceanStateEstimation],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Functionalities" => "API.md",
    ],
    repo="https://github.com/gaelforget/OceanStateEstimation.jl/blob/{commit}{path}#L{line}",
    sitename="OceanStateEstimation.jl",
    authors="gaelforget <gforget@mit.edu>",
    assets=String[],
)

OceanStateEstimation.CBIOMESclim_download()
OceanStateEstimation.ECCOdiags_download()
#lst=("CBIOMES_climatogy_create.jl",)
lst=("CBIOMES_climatology_plot.jl","ECCO_standard_plots.jl")

pth_out=joinpath(@__DIR__,"build", "examples")
!isdir(pth_out) ? mkdir(pth_out) : nothing

if true #skip generating file automatically 
    for i in lst
    println("PlutoSliderServer($i) started")
    fil_in=joinpath(@__DIR__,"..", "examples",i)
    fil_out=joinpath(pth_out,i[1:end-2]*"html")
    PlutoSliderServer.export_notebook(fil_in)
    mv(fil_in[1:end-2]*"html",fil_out)
    cp(fil_in,fil_out[1:end-4]*"jl")
    println("PlutoSliderServer($i) completed")
end
end

deploydocs(;
    repo="github.com/gaelforget/OceanStateEstimation.jl",
)

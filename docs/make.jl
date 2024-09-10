using Documenter, Climatology
import MITgcm, MeshArrays, NCDatasets, Pkg
import PlutoSliderServer, CairoMakie
import Climatology: Downloads

Pkg.precompile()

makedocs(;
    modules=[Climatology],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
        "Data Sets" => "examples.md",
        "Data Files" => "API.md",
        "Internals (ECCO)" => "ECCO.md",
    ],
    repo="https://github.com/JuliaOcean/Climatology.jl/blob/{commit}{path}#L{line}",
    sitename="Climatology.jl",
    authors="JuliaOcean <gforget@mit.edu>",
    warnonly = [:cross_references,:missing_docs],
)

Climatology.CBIOMESclim_download()
Climatology.ECCOdiags_add("release2")
Climatology.ECCOdiags_add("OCCA2HR1")
withenv(MeshArrays.interpolation_setup,"DATADEPS_ALWAYS_ACCEPT"=>true)

#temporary fix, which should no longer be needed after 0.5.9 is released
earth_jpg=joinpath(tempdir(),"Blue_Marble.jpg")
url="https://upload.wikimedia.org/wikipedia/commons/5/56/Blue_Marble_Next_Generation_%2B_topography_%2B_bathymetry.jpg"
!isfile(earth_jpg) ? Downloads.download(url,earth_jpg) : nothing

lst=("SSH/SatelliteAltimetry.jl",
     "OISST/sst_anomaly_notebook.jl","CBIOMES/CBIOMES_climatology_plot.jl",
     "ECCO/ECCO_standard_plots.jl","NSLCT/NSLCT_notebook.jl",
     "OptimalTransport/OptimalTransport_demo.jl","HadIOD/HadIOD_viz.jl")

pth_out=joinpath(@__DIR__,"build", "examples")
!isdir(pth_out) ? mkdir(pth_out) : nothing

#lst=[]
for i in lst
    println("PlutoSliderServer($i) started")
    fil_in=joinpath(@__DIR__,"..", "examples",i)
    fil_out=joinpath(pth_out,basename(i)[1:end-2]*"html")
    PlutoSliderServer.export_notebook(fil_in)
    mv(fil_in[1:end-2]*"html",fil_out)
    cp(fil_in,fil_out[1:end-4]*"jl")
    println("PlutoSliderServer($i) completed")
end

deploydocs(;
    repo="github.com/JuliaOcean/Climatology.jl",
)

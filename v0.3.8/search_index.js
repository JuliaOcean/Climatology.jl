var documenterSearchIndex = {"docs":
[{"location":"ECCO/","page":"ECCO","title":"ECCO","text":"The workflow presented here is as follows.","category":"page"},{"location":"ECCO/","page":"ECCO","title":"ECCO","text":"set up recommended for running analyses of ECCO estimates.\nrun one computation loop on the ECCO monthly files.","category":"page"},{"location":"ECCO/","page":"ECCO","title":"ECCO","text":"ECCO.standard_analysis_setup","category":"page"},{"location":"ECCO/#OceanStateEstimation.ECCO.standard_analysis_setup","page":"ECCO","title":"OceanStateEstimation.ECCO.standard_analysis_setup","text":"ECCO.standard_analysis_setup(pth0::String)\n\nCreate temporary run folder pth where data folder pth0 will be linked. \n\nData folder pth0 should be the path to ECCO data.\n\nFor example:\n\nusing OceanStateEstimation, Pkg\npth=ECCO.standard_analysis_setup(ScratchSpaces.ECCO)\n\nThe Project.toml file found in pth provides an environment ready for ECCO analyses. \n\nThis environment can be activated and instantiated:\n\nPkg.activate(pth)\nPkg.instantiate()\n\n\n\n\n\n","category":"function"},{"location":"ECCO/","page":"ECCO","title":"ECCO","text":"Here is an example of parameters P to compute zonal mean temperatures at level 5.","category":"page"},{"location":"ECCO/","page":"ECCO","title":"ECCO","text":"ECCO_helpers.parameters","category":"page"},{"location":"ECCO/#OceanStateEstimation.ECCO_helpers.parameters","page":"ECCO","title":"OceanStateEstimation.ECCO_helpers.parameters","text":"parameters(P0,params)\n\nPrepare parameter NamedTuple for use in ECCO_diagnostics.driver.\n\nP1=parameters(P0,p) \n\nis faster than e.g. parameters(pth,\"r2\",p) as grid, etc get copied from P0 to P1.\n\n\n\n\n\nparameters(pth0::String,sol0::String,params)\n\nPrepare parameter NamedTuple for use in ECCO_diagnostics.driver.\n\nFor example, to compute zonal mean temperatures at level 5:\n\np=(calc = \"zonmean\", nam = \"THETA\", lev = 5)\npth=ECCO.standard_analysis_setup(ScratchSpaces.ECCO)\nP0=ECCO_helpers.parameters(pth,\"r2\",p)\n\nor, from a predefined list:\n\nlist0=ECCO_helpers.standard_list_toml(\"\")\npth=ECCO.standard_analysis_setup(ScratchSpaces.ECCO)\nP1=ECCO_helpers.parameters(pth,\"r2\",list0[1])\n\n\n\n\n\n","category":"function"},{"location":"ECCO/","page":"ECCO","title":"ECCO","text":"The computation loop, over all months, can then be carried out as follows.","category":"page"},{"location":"ECCO/","page":"ECCO","title":"ECCO","text":"ECCO_diagnostics.driver","category":"page"},{"location":"ECCO/#OceanStateEstimation.ECCO_diagnostics.driver","page":"ECCO","title":"OceanStateEstimation.ECCO_diagnostics.driver","text":"driver(P)\n\nCall main computation loop as specified by parameters P.\n\nThe main computation loop choice depends on the P parameter values. Methods include:\n\nmain_clim\nmain_glo\nmain_zonmean\nmain_overturn\nmain_MHT\nmain_trsp\n\n\n\n\n\n","category":"function"},{"location":"ECCO/","page":"ECCO","title":"ECCO","text":"Modules = [OceanStateEstimation.ECCO_io]","category":"page"},{"location":"ECCO/#OceanStateEstimation.ECCO_io.read_monthly-Tuple{Any, Any, Any}","page":"ECCO","title":"OceanStateEstimation.ECCO_io.read_monthly","text":"read_monthly(P,nam,t)\n\nRead record t for variable nam from file locations specified via parameters P.\n\nThe method used to read nam is selected based on nam's value. Methods include:\n\nread_monthly_default\nread_monthly_SSH\nread_monthly_MHT\nread_monthly_BSF\n\n\n\n\n\n","category":"method"},{"location":"examples/#Physical-Oceanography","page":"Examples","title":"Physical Oceanography","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"ECCO_standard_plots.jl (➭ code link) : explore the Ocean climatology, currents, and more.\nNSLCT_notebook.jl (➭ code link) : plot global mean sea level data from NASA and Dataverse.","category":"page"},{"location":"examples/#Detail","page":"Examples","title":"Detail","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"ECCO_standard_plots.jl lets you explore climate indices and climatologies derived (via this and that) from gridded ocean climatologies. The data is retrieved from dataverse.org, and intermediate results from zenodo.org.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"NSLCT_notebook.jl lets you access sea level data from NASA and Dataver portals (HTTP.jl, Dataverse.jl), organize it into tables (DataFrames.jl), and plot it (Makie.jl).","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"note: Note\nFor more on these estimates, and how to use them in Julia, please refer to the following documentation and links therein.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"MeshArrays.jl : gridded Earth variables, domain decomposition, C-grid support; Ocean Circulation, Geography tutorials.\nMITgcmTools.jl : framework to interact with MITgcm (setup, run, output, plot, etc) and ECCO output.\nIndividualDisplacements.jl : simulation and analysis of materials moving through oceanic and atmospheric flows.","category":"page"},{"location":"examples/#Marine-Ecosystems","page":"Examples","title":"Marine Ecosystems","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"CBIOMES_climatology_plot (➭ code link) : visualize ocean colour and biomass climatologies\nCBIOMES_climatology_create (➭ code link) : recreate the CBIOMES-global climatology files","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The CBIOMES-global climatology (alpha version) is a global ocean state estimate that covers the period from 1992 to 2011 (ECCO). It is based on Forget et al 2015 for ocean physics MIT general circulation model and on Dutkiewicz et al 2015 for marine biogeochemistry and ecosystems Darwin Project model.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"note: Note\nFor more notebooks involving CBIOMES and related efforts, take a look at the following pages.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Marine Ecosystem Notebooks : Darwin Model, Ocean Color data, Gradients field program, and more.\nJuliaCon2021 workshop : Modeling Marine Ecosystems At Multiple Scales Using Julia.\nPlanktonIndividuals.jl : simulate the life cycle of ocean phytoplankton cells and their environment.","category":"page"},{"location":"examples/#More-Notebooks","page":"Examples","title":"More Notebooks","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"OptimalTransport_demo.jl : using optimal transport for e.g. model-data comparison\nHadIOD_viz.jl : download, read, and plot a subset of the HadIOD T/S database","category":"page"},{"location":"examples/#User-Manual","page":"Examples","title":"User Manual","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"To run the notebook on a local computer or in the cloud, please refer to the Pluto docs. Directions are also provided in the following pages.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"ECCO_standard_plots.jl (see final section)\nJuliaClimate How-To \nECCO/Julia storymap\nvideo demonstration","category":"page"},{"location":"examples/#References","page":"Examples","title":"References","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"OCCA : Forget 2010\nECCO v4 : Forget et al 2015\nCBIOMES-global : Forget 2018","category":"page"},{"location":"#OceanStateEstimation.jl","page":"Home","title":"OceanStateEstimation.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package is currently focused on serving and deriving climatologies from ocean state estimates. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"See Physical Oceanography and Marine Ecosystems for examples.","category":"page"},{"location":"","page":"Home","title":"Home","text":"It is in early development stage; breaking changes remain likely.","category":"page"},{"location":"API/#Intro","page":"Files","title":"Intro","text":"","category":"section"},{"location":"API/","page":"Files","title":"Files","text":"Climatologies are readily accessed using the Scratch.jl artifact system via OceanStateEstimation.jl as explained below. ","category":"page"},{"location":"API/#Use-Examples","page":"Files","title":"Use Examples","text":"","category":"section"},{"location":"API/#ECCO","page":"Files","title":"ECCO","text":"","category":"section"},{"location":"API/","page":"Files","title":"Files","text":"ECCO climatology files can downloaded using get_ecco_files. These files are for version 4 release 2, on the native model grid.","category":"page"},{"location":"API/","page":"Files","title":"Files","text":"using OceanStateEstimation, MeshArrays\nγ=GridSpec(\"LatLonCap\",MeshArrays.GRID_LLC90)\ntmp=OceanStateEstimation.get_ecco_files(γ,\"ETAN\")","category":"page"},{"location":"API/","page":"Files","title":"Files","text":"Precomputed quantities shown in ECCO_standard_plots.jl can be downloaded separately.","category":"page"},{"location":"API/","page":"Files","title":"Files","text":"OceanStateEstimation.ECCOdiags_add(\"release2\")\nreaddir(ScratchSpaces.ECCO)","category":"page"},{"location":"API/#OCCA","page":"Files","title":"OCCA","text":"","category":"section"},{"location":"API/","page":"Files","title":"Files","text":"get_occa_variable_if_needed(\"SIarea\")\nreaddir(ScratchSpaces.OCCA)","category":"page"},{"location":"API/#CBIOMES","page":"Files","title":"CBIOMES","text":"","category":"section"},{"location":"API/","page":"Files","title":"Files","text":"To retrieve the CBIOMES climatology, in the julia REPL for example :","category":"page"},{"location":"API/","page":"Files","title":"Files","text":"OceanStateEstimation.CBIOMESclim_download()\nreaddir(ScratchSpaces.CBIOMES)","category":"page"},{"location":"API/","page":"Files","title":"Files","text":"And the files, now found in ScratchSpaces.CBIOMES, can then be read using other libraries.","category":"page"},{"location":"API/","page":"Files","title":"Files","text":"using NCDatasets\nfil_out=joinpath(ScratchSpaces.CBIOMES,\"CBIOMES-global-alpha-climatology.nc\")\nnc=NCDataset(fil_out,\"r\")\nkeys(nc)","category":"page"},{"location":"API/#Path-Names","page":"Files","title":"Path Names","text":"","category":"section"},{"location":"API/","page":"Files","title":"Files","text":"Gridded fields are mostly retrieved from Harvard Dataverse. These can be relatively large files, compared to the package codes, so they are handled lazily (only downloaded when needed). Precomputed diagnostics have also been archived on zenodo.org.","category":"page"},{"location":"API/","page":"Files","title":"Files","text":"Artifact Name File Type Download Method\nScratchSpaces.ECCO NetCDF lazy, by variable, dataverse\nScratchSpaces.ECCO JLD2 lazy, whole, zenodo\nScratchSpaces.MITprof binary lazy, whole, zenodo\nScratchSpaces.OCCA NetCDF lazy, by variable, dataverse\nScratchSpaces.CBIOMES NetCDF lazy, whole, zenodo","category":"page"},{"location":"API/#Functions-Reference","page":"Files","title":"Functions Reference","text":"","category":"section"},{"location":"API/","page":"Files","title":"Files","text":"Modules = [OceanStateEstimation.downloads]","category":"page"},{"location":"API/#OceanStateEstimation.downloads.CBIOMESclim_download-Tuple{}","page":"Files","title":"OceanStateEstimation.downloads.CBIOMESclim_download","text":"CBIOMESclim_download()\n\nDownload lazy artifact to scratch space.\n\n\n\n\n\n","category":"method"},{"location":"API/#OceanStateEstimation.downloads.ECCOdiags_add-Tuple{String}","page":"Files","title":"OceanStateEstimation.downloads.ECCOdiags_add","text":"ECCOdiags_add(nam::String)\n\nAdd data to the scratch space folder. Known options for nam include  \"release1\", \"release2\", \"release3\", \"release4\", and \"interp_coeffs\".\n\n\n\n\n\n","category":"method"},{"location":"API/#OceanStateEstimation.downloads.MITPROFclim_download-Tuple{}","page":"Files","title":"OceanStateEstimation.downloads.MITPROFclim_download","text":"MITPROFclim_download()\n\nDownload lazy artifact to scratch space.\n\n\n\n\n\n","category":"method"},{"location":"API/#OceanStateEstimation.downloads.get_ecco_files","page":"Files","title":"OceanStateEstimation.downloads.get_ecco_files","text":"get_ecco_files(γ::gcmgrid,v::String,t=1)\n\nusing MeshArrays, OceanStateEstimation, MITgcmTools\nγ=GridSpec(\"LatLonCap\",MeshArrays.GRID_LLC90)\ntmp=OceanStateEstimation.get_ecco_files(γ,\"oceQnet\")\n\n\n\n\n\n","category":"function"},{"location":"API/#OceanStateEstimation.downloads.get_ecco_variable_if_needed-Tuple{String}","page":"Files","title":"OceanStateEstimation.downloads.get_ecco_variable_if_needed","text":"get_ecco_variable_if_needed(v::String)\n\nDownload ECCO output for variable v to scratch space if needed\n\n\n\n\n\n","category":"method"},{"location":"API/#OceanStateEstimation.downloads.get_ecco_velocity_if_needed-Tuple{}","page":"Files","title":"OceanStateEstimation.downloads.get_ecco_velocity_if_needed","text":"get_ecco_velocity_if_needed()\n\nDownload ECCO output for u,v,w to scratch space if needed\n\n\n\n\n\n","category":"method"},{"location":"API/#OceanStateEstimation.downloads.get_occa_variable_if_needed-Tuple{String}","page":"Files","title":"OceanStateEstimation.downloads.get_occa_variable_if_needed","text":"get_occa_variable_if_needed(v::String)\n\nDownload OCCA output for variable v to scratch space if needed\n\n\n\n\n\n","category":"method"},{"location":"API/#OceanStateEstimation.downloads.get_occa_velocity_if_needed-Tuple{}","page":"Files","title":"OceanStateEstimation.downloads.get_occa_velocity_if_needed","text":"get_occa_velocity_if_needed()\n\nDownload OCCA output for u,v,w to scratch space if needed\n\n\n\n\n\n","category":"method"}]
}
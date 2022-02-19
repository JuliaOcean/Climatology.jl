var documenterSearchIndex = {"docs":
[{"location":"examples/#Physical-Oceanography","page":"Examples","title":"Physical Oceanography","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"ECCO_standard_plots.jl (➭ code link) depicts a selection of climate-relevant variables and indices as an example. These were derived from gridded estimates of the ocean state for physical variables like temperature, salinity, and currents (see ECCO_standard_analysis.jl) and retrieved from dataverse or zenodo.org (see OceanStateEstimation.ECCOdiags_download, and OceanStateEstimation.ECCOdiags_add).","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"To run the notebook proceed as follows. Once Pluto opens, in your web browser, paste the code link url (not the julia code itself) and click open. ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using Pluto\n\nimport OceanStateEstimation\nOceanStateEstimation.ECCOdiags_download()\nOceanStateEstimation.ECCOdiags_add(\"interp_coeffs\")\nOceanStateEstimation.ECCOdiags_add(\"release4\")\n\nPluto.run()","category":"page"},{"location":"examples/#Bio-Geo-Chemical-Climatology","page":"Examples","title":"Bio-Geo-Chemical Climatology","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"CBIOMES-global (alpha version) is a global ocean state estimate that covers the period from 1992 to 2011. It is based on Forget et al 2015 for ocean physics MIT general circulation model and on Dutkiewicz et al 2015 for marine biogeochemistry and ecosystems Darwin Project model.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"CBIOMES_climatology_plot (➭ code link) : visualize the climatology maps interactively using Pluto.jl\nCBIOMES_climatology_create (➭ code link) : recreate the climatology file. The original is archived here in zenodo.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"To retrieve this climatology, in the julia REPL for example :","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using OceanStateEstimation, NCTiles\nOceanStateEstimation.CBIOMESclim_download()\nfil_out=joinpath(CBIOMESclim_path,\"CBIOMES-global-alpha-climatology.nc\")\nnc=NCTiles.NCDataset(fil_out,\"r\")","category":"page"},{"location":"examples/#References","page":"Examples","title":"References","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"OCCA : Forget 2010\nECCO v4 : Forget et al 2015\nCBIOMES-global : Forget 2018","category":"page"},{"location":"#OceanStateEstimation.jl","page":"Home","title":"OceanStateEstimation.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package is currently focused on serving and deriving climatologies from ocean state estimates. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"See Physical Oceanography and Bio-Geo-Chemical Climatology for examples.","category":"page"},{"location":"","page":"Home","title":"Home","text":"It is in early development stage; breaking changes remain likely.","category":"page"},{"location":"API/#Data-Access","page":"Functionalities","title":"Data Access","text":"","category":"section"},{"location":"API/","page":"Functionalities","title":"Functionalities","text":"The gridded fields used in the ECCO examples can be retrieved from ecco-group.org and, for the ECCOv4r2 estimate, from Harvard Dataverse or zenodo.org. Two monthly climatologies (ECCOv4r2 and OCCA) are also readily available using the Julia artifact system as explained below. These can be relatively large files, compared to the package codes, so they are handled lazily (only downloaded when needed). ","category":"page"},{"location":"API/","page":"Functionalities","title":"Functionalities","text":"Artifact path File Type Download Method\nECCOclim_path NetCDF lazy, by variable, dataverse\nOCCAclim_path NetCDF lazy, by variable, dataverse\nMITPROFclim_path binary lazy, whole, zenodo\nECCOdiags_path JLD2 lazy, whole, zenodo","category":"page"},{"location":"API/#Basic-Usage","page":"Functionalities","title":"Basic Usage","text":"","category":"section"},{"location":"API/","page":"Functionalities","title":"Functionalities","text":"For ECCO :","category":"page"},{"location":"API/","page":"Functionalities","title":"Functionalities","text":"using OceanStateEstimation, MeshArrays\nγ=GridSpec(\"LatLonCap\",MeshArrays.GRID_LLC90)\ntmp=OceanStateEstimation.get_ecco_files(γ,\"ETAN\")","category":"page"},{"location":"API/","page":"Functionalities","title":"Functionalities","text":"For OCCA :","category":"page"},{"location":"API/","page":"Functionalities","title":"Functionalities","text":"using OceanStateEstimation\nget_occa_variable_if_needed(\"SIarea\")\nreaddir(OCCAclim_path)","category":"page"},{"location":"API/#Functions-Reference","page":"Functionalities","title":"Functions Reference","text":"","category":"section"},{"location":"API/","page":"Functionalities","title":"Functionalities","text":"Modules = [OceanStateEstimation]","category":"page"},{"location":"API/#OceanStateEstimation.CBIOMESclim_download-Tuple{}","page":"Functionalities","title":"OceanStateEstimation.CBIOMESclim_download","text":"CBIOMESclim_download()\n\nDownload lazy artifact to CBIOMESclim_path.\n\n\n\n\n\n","category":"method"},{"location":"API/#OceanStateEstimation.ECCOdiags_add-Tuple{String}","page":"Functionalities","title":"OceanStateEstimation.ECCOdiags_add","text":"ECCOdiags_add(nam::String)\n\nAdd data to the ECCOdiagspath folder. Known options for nam include  \"release1\", \"release3\", \"release4\", and \"interpcoeffs\". Note that  \"release2\" is the estimate that's readily donwloaded by ECCOdiags_download().\n\n\n\n\n\n","category":"method"},{"location":"API/#OceanStateEstimation.ECCOdiags_download-Tuple{}","page":"Functionalities","title":"OceanStateEstimation.ECCOdiags_download","text":"ECCOdiags_download()\n\nDownload lazy artifact to ECCOdiags_path.\n\n\n\n\n\n","category":"method"},{"location":"API/#OceanStateEstimation.MITPROFclim_download-Tuple{}","page":"Functionalities","title":"OceanStateEstimation.MITPROFclim_download","text":"MITPROFclim_download()\n\nDownload lazy artifact to MITPROFclim_path.\n\n\n\n\n\n","category":"method"},{"location":"API/#OceanStateEstimation.dataverse_lists-Tuple{String}","page":"Functionalities","title":"OceanStateEstimation.dataverse_lists","text":"dataverse_lists(lst::String)\n\nRead and derive lists (ID,name,URL) from csv file (ID,name) and return as tuple\n\nlists=dataverse_lists(lst)\n\n\n\n\n\n","category":"method"},{"location":"API/#OceanStateEstimation.get_ecco_files","page":"Functionalities","title":"OceanStateEstimation.get_ecco_files","text":"get_ecco_files(γ::gcmgrid,v::String,t=1)\n\nusing MeshArrays, OceanStateEstimation\nγ=GridSpec(\"LatLonCap\",MeshArrays.GRID_LLC90)\ntmp=get_ecco_files(γ,\"oceQnet\")\n\n\n\n\n\n","category":"function"},{"location":"API/#OceanStateEstimation.get_ecco_variable_if_needed-Tuple{String}","page":"Functionalities","title":"OceanStateEstimation.get_ecco_variable_if_needed","text":"get_ecco_variable_if_needed(v::String)\n\nDownload ECCO output for variable v to ECCOclim_path if needed\n\n\n\n\n\n","category":"method"},{"location":"API/#OceanStateEstimation.get_ecco_velocity_if_needed-Tuple{}","page":"Functionalities","title":"OceanStateEstimation.get_ecco_velocity_if_needed","text":"get_ecco_velocity_if_needed()\n\nDownload ECCO output for u,v,w to ECCOclim_path if needed\n\n\n\n\n\n","category":"method"},{"location":"API/#OceanStateEstimation.get_from_dataverse-Tuple{String, String, String}","page":"Functionalities","title":"OceanStateEstimation.get_from_dataverse","text":"get_from_dataverse(lst::String,nam::String,pth::String)\n\nusing OceanStateEstimation, CSV, DataFrames\npth=dirname(pathof(OceanStateEstimation))\nlst=joinpath(pth,\"../examples/OCCA_climatology.csv\")\nnams=CSV.read(lst,DataFrame)\n[get_from_dataverse(lst,nam,OCCAclim_path) for nam in nams.name[:]]\n\n\n\n\n\n","category":"method"},{"location":"API/#OceanStateEstimation.get_occa_variable_if_needed-Tuple{String}","page":"Functionalities","title":"OceanStateEstimation.get_occa_variable_if_needed","text":"get_occa_variable_if_needed(v::String)\n\nDownload OCCA output for variable v to OCCAclim_path if needed\n\n\n\n\n\n","category":"method"},{"location":"API/#OceanStateEstimation.get_occa_velocity_if_needed-Tuple{}","page":"Functionalities","title":"OceanStateEstimation.get_occa_velocity_if_needed","text":"get_occa_velocity_if_needed()\n\nDownload OCCA output for u,v,w to OCCAclim_path if needed\n\n\n\n\n\n","category":"method"}]
}


## Physical Oceanography

- [Sea Surface Temperature](sst_anomaly_notebook.html) (➭ [code link](https://raw.githubusercontent.com/JuliaOcean/Climatology.jl/master/examples/OISST/sst_anomaly_notebook.jl)) : plot global mean and regional sea surface temperature (NOAA's OISST).
- [Sea Level Estimates](NSLCT_notebook.html) (➭ [code link](https://raw.githubusercontent.com/JuliaOcean/Climatology.jl/master/examples/NSLCT/NSLCT_notebook.jl)) : plot global mean and regional sea level data from NASA (NASA's NSLCT and ECCO)).
- [Ocean State Estimates](ECCO_standard_plots.html) (➭ [code link](https://raw.githubusercontent.com/JuliaOcean/Climatology.jl/master/examples/ECCO/ECCO_standard_plots.jl)) : explore ocean transports, climate indices, siub-surface temperature, and many more variables from full ocean state estimates (ECCO, OCCA).

#### Detail

- [sst\_anomaly\_notebook.jl](sst_anomaly_notebook.html) plot maps as well as time series of SST and SST anomalies.
- [NSLCT\_notebook.jl](NSLCT_notebook.html) lets you access sea level data from NASA and Dataver portals (`HTTP.jl`, `Dataverse.jl`), organize it into tables (`DataFrames.jl`), and plot it (`Makie.jl`).
- [ECCO\_standard\_plots.jl](ECCO_standard_plots.html) lets you explore climate indices and climatologies derived (via [this](https://github.com/JuliaOcean/Climatology.jl/blob/master/examples/ECCO/ECCO_standard_calcs.jl) and [that](https://github.com/JuliaOcean/Climatology.jl/blob/master/examples/ECCO/ECCO_standard_loop.jl)) from gridded, time-variable ocean climatologies (ECCO4, OCCA2). The data is retrieved from [dataverse.org](https://dataverse.harvard.edu/dataverse/ECCO), and intermediate results from [zenodo.org](https://zenodo.org).


## Marine Ecosystems

- [Plankton, Chemistry, and Light](CBIOMES_climatology_plot.html) (➭ [code link](https://raw.githubusercontent.com/JuliaOcean/Climatology.jl/master/examples/CBIOMES/CBIOMES_climatology_plot.jl)) : visualize ocean colour and biomass climatologies estimated using the Darwin3 model.

#### Detail

- The [CBIOMES1](https://github.com/CBIOMES/global-ocean-model) climatology (alpha version) is a global ocean state estimate that covers the period from 1992 to 2011 (ECCO). It is based on Forget et al 2015 for ocean physics [MIT general circulation model](https://mitgcm.readthedocs.io/en/latest/#) and on Dutkiewicz et al 2015 for marine biogeochemistry and ecosystems [Darwin Project model](https://darwin3.readthedocs.io/en/latest/phys_pkgs/darwin.html).
- [CBIOMES\_climatology\_create](https://JuliaOcean.github.io/Climatology.jl/v0.1.13/examples/CBIOMES_model_climatogy.html) (➭ [code link](https://raw.githubusercontent.com/JuliaOcean/Climatology.jl/master/examples/CBIOMES/CBIOMES_climatology_create.jl)) : recreate the CBIOMES-global climatology files
- [OptimalTransport\_demo.jl](OptimalTransport_demo.html) : using optimal transport for model-data comparison (CBIOMES1 vs satellite data).

## Other Notebooks

- [HadIOD\_viz.jl](HadIOD_viz.html) : download, read, and plot a subset of the [HadIOD](https://www.metoffice.gov.uk/hadobs/hadiod/) T/S database
- [SatelliteAltimetry.html](https://juliaocean.github.io/OceanRobots.jl/dev/examples/SatelliteAltimetry.html) (➭ [code link](https://raw.githubusercontent.com/JuliaOcean/OceanRobots.jl/master/examples/SatelliteAltimetry.jl)) : gridded satellite data for sea level anomaly
- the suite of examples provided in [OceanRobots.jl](https://juliaocean.github.io/OceanRobots.jl/dev/examples/) that focus on observations collected at sea.

## References

- OCCA1 : [Forget 2010](https://doi.org/10.1175/2009JPO4043.1)
- ECCO4 : [Forget et al 2015](https://gmd.copernicus.org/articles/8/3071/2015/)
- CBIOMES1: [Forget 2018](https://zenodo.org/record/2653669#.YbwAUi1h0ow)
- OCCA2 : [Forget 2024](https://doi.org/10.21203/rs.3.rs-3979671/v1)
	
## Notes

!!! note
    For more on these estimates, and how to use them in Julia, please refer to the following documentation and links therein.

- [MeshArrays.jl](https://juliaclimate.github.io/MeshArrays.jl/dev/) : gridded Earth variables, domain decomposition, C-grid support; [Ocean Circulation](https://juliaclimate.github.io/MeshArrays.jl/dev/tutorials/vectors.html), [Geography](https://juliaclimate.github.io/MeshArrays.jl/dev/tutorials/geography.html) tutorials.
- [MITgcm.jl](https://gaelforget.github.io/MITgcm.jl/dev/) : framework to interact with MITgcm (setup, run, output, plot, etc), CBIOMES, and ECCO output.
- [IndividualDisplacements.jl](https://juliaclimate.github.io/IndividualDisplacements.jl/dev/) : simulation and analysis of materials moving through oceanic and atmospheric flows.

!!! note
    For more notebooks involving [CBIOMES](https://cbiomes.org) and related efforts, take a look at the following pages.

- [Marine Ecosystem Notebooks](https://github.com/JuliaOcean/MarineEcosystemNotebooks) : Darwin3 model, ocean color data, gradients field program, and more.
- [JuliaCon2021 workshop](https://github.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl) : _Modeling Marine Ecosystems At Multiple Scales Using Julia_.
- [PlanktonIndividuals.jl](https://juliaocean.github.io/PlanktonIndividuals.jl/dev/) : simulate the life cycle of ocean phytoplankton cells and their environment.

!!! note
    To run the notebook on a local computer or in the cloud, please refer to the [Pluto docs](https://github.com/fonsp/Pluto.jl/wiki). Directions are also provided in the following pages.

- [ECCO\_standard\_plots.jl](https://JuliaOcean.github.io/Climatology.jl/dev/examples/ECCO_standard_plots.html)
- [JuliaClimate How-To](https://juliaclimate.github.io/Notebooks/#directions) 
- [ECCO/Julia storymap](https://ecco-group.org/storymaps.htm?id=69)
- [video demonstration](https://www.youtube.com/watch?v=mZevMagHatc&list=PLXO7Tdh24uhPFZ5bph6Y_Q3-CRSfk5cDU)


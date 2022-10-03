
## Physical Oceanography

- [NSLCT\_notebook.jl](NSLCT_notebook.html) (➭ [code link](https://raw.githubusercontent.com/gaelforget/OceanStateEstimation.jl/master/examples/NSLCT/NSLCT_notebook.jl)) : download and plot sea level data from NASA portal
- [ECCO\_standard\_plots.jl](ECCO_standard_plots.html) (➭ [code link](https://raw.githubusercontent.com/gaelforget/OceanStateEstimation.jl/master/examples/ECCO/ECCO_standard_plots.jl)) : explore the Global Ocean (climatologies, currents, ...)

[NSLCT\_notebook.jl](NSLCT_notebook.html) lets you access sea level data from NASA portal (`HTTP.jl`), organize it (`DataFrames.jl`), and plot it (`Makie.jl`).

[ECCO\_standard\_plots.jl](ECCO_standard_plots.html) lets you explore climate indices and climatologies derived (via [this](https://github.com/gaelforget/OceanStateEstimation.jl/blob/master/examples/ECCO/ECCO_standard_calcs.jl) and [that](https://github.com/gaelforget/OceanStateEstimation.jl/blob/master/examples/ECCO/ECCO_standard_loop.jl)) from gridded ocean climatologies. The data is retrieved from [dataverse.org](https://dataverse.harvard.edu/dataverse/ECCO), and intermediate results from [zenodo.org](https://zenodo.org).

!!! note
    For more on these gridded estimates, and how to use them, see :

- [MeshArrays.jl](https://juliaclimate.github.io/MeshArrays.jl/dev/) : gridded Earth variables, domain decomposition, C-grid support; [Ocean Circulation](https://juliaclimate.github.io/MeshArrays.jl/dev/tutorials/vectors.html), [Geography](https://juliaclimate.github.io/MeshArrays.jl/dev/tutorials/geography.html) tutorials.
- [MITgcmTools.jl](https://juliaclimate.github.io/MiTgcmTools.jl/dev/) : framework to interact with MITgcm (setup, run, output, plot, etc) and ECCO output.
- [IndividualDisplacements.jl](https://juliaclimate.github.io/IndividualDisplacements.jl/dev/) : simulation and analysis of materials moving through oceanic and atmospheric flows.

## Marine Ecosystems

- [CBIOMES\_climatology\_plot](CBIOMES_climatology_plot.html) (➭ [code link](https://raw.githubusercontent.com/gaelforget/OceanStateEstimation.jl/master/examples/CBIOMES/CBIOMES_climatology_plot.jl)) : visualize ocean colour and biomass climatologies
- [CBIOMES\_climatology\_create](https://gaelforget.github.io/OceanStateEstimation.jl/v0.1.13/examples/CBIOMES_model_climatogy.html) (➭ [code link](https://raw.githubusercontent.com/gaelforget/OceanStateEstimation.jl/master/examples/CBIOMES/CBIOMES_climatology_create.jl)) : recreate the CBIOMES-global climatology files

The [CBIOMES-global](https://github.com/CBIOMES/global-ocean-model) climatology (alpha version) is a global ocean state estimate that covers the period from 1992 to 2011 (ECCO). It is based on Forget et al 2015 for ocean physics [MIT general circulation model](https://mitgcm.readthedocs.io/en/latest/#) and on Dutkiewicz et al 2015 for marine biogeochemistry and ecosystems [Darwin Project model](https://darwin3.readthedocs.io/en/latest/phys_pkgs/darwin.html).

!!! note
    For more notebooks involving [CBIOMES](https://cbiomes.org) and related efforts :

- [Marine Ecosystem Notebooks](https://github.com/JuliaOcean/MarineEcosystemNotebooks) : Darwin Model, Ocean Color data, Gradients field program, and more.
- [JuliaCon2021 workshop](https://github.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl) : _Modeling Marine Ecosystems At Multiple Scales Using Julia_.
- [PlanktonIndividuals.jl](https://juliaocean.github.io/PlanktonIndividuals.jl/dev/) : simulate the life cycle of ocean phytoplankton cells and their environment.

## More Notebooks

- [OptimalTransport\_demo.jl](https://gaelforget.github.io/OceanStateEstimation.jl/dev/examples/OptimalTransport_demo.html) : using optimal transport for e.g. model-data comparison
- [HadIOD\_viz.jl](https://github.com/gaelforget/OceanStateEstimation.jl/blob/master/examples/HadIOD/HadIOD_viz.jl) : download, read, and plot a subset of the [HadIOD](https://www.metoffice.gov.uk/hadobs/hadiod/) T/S database

## User Manual

To run the notebook on a local computer or in the cloud, please refer to directions provided in :

- [ECCO\_standard\_plots.jl](https://gaelforget.github.io/OceanStateEstimation.jl/dev/examples/ECCO_standard_plots.html) (see final section)
- [JuliaClimate How-To](https://juliaclimate.github.io/Notebooks/#directions) 
- [ECCO/Julia storymap](https://ecco-group.org/storymaps.htm?id=69)
- [video demonstration](https://www.youtube.com/watch?v=mZevMagHatc&list=PLXO7Tdh24uhPFZ5bph6Y_Q3-CRSfk5cDU)

## References

- OCCA : [Forget 2010]()
- ECCO v4 : [Forget et al 2015](https://gmd.copernicus.org/articles/8/3071/2015/)
- CBIOMES-global : [Forget 2018](https://zenodo.org/record/2653669#.YbwAUi1h0ow)
	


## Physical Oceanography

[ECCO\_standard\_plots.jl](ECCO_standard_plots.html) (➭ [code link](https://raw.githubusercontent.com/gaelforget/OceanStateEstimation.jl/master/examples/ECCO/ECCO_standard_plots.jl)) depicts a few climate-relevant variables and indices as an example. These quantities were derived (see [this notebook](https://github.com/gaelforget/OceanStateEstimation.jl/blob/master/examples/ECCO/ECCO_standard_calcs.jl) and [this script](https://github.com/gaelforget/OceanStateEstimation.jl/blob/master/examples/ECCO/ECCO_standard_loop.jl)) from estimates of the ocean state that provide gridded fields for physical variables like temperature, salinity, and currents. The initial fields can be retrieved from the [ECCO dataverse](https://dataverse.harvard.edu/dataverse/ECCO), and intermediate results from [zenodo.org](https://zenodo.org).

If you want to run the notebook on your local computer or in the cloud, please refer to the directions provided at the [notebook final section](https://gaelforget.github.io/OceanStateEstimation.jl/dev/examples/ECCO_standard_plots.html). Or more generally, to : 

- [JuliaClimate Notebooks How-To](https://juliaclimate.github.io/Notebooks/#directions) 
- [ecco-group.org story map](https://ecco-group.org/storymaps.htm?id=69)
- [video demonstration](https://www.youtube.com/watch?v=mZevMagHatc&list=PLXO7Tdh24uhPFZ5bph6Y_Q3-CRSfk5cDU)
	
For more notebooks on these gridded estimates, and how to analyze them, see :

- [MeshArrays.jl](https://juliaclimate.github.io/MeshArrays.jl/dev/) : gridded Earth variables, domain decomposition, C-grid support; [Ocean Circulation](https://juliaclimate.github.io/MeshArrays.jl/dev/tutorials/vectors.html), [Geography](https://juliaclimate.github.io/MeshArrays.jl/dev/tutorials/geography.html) tutorials.
- [IndividualDisplacements.jl](https://juliaclimate.github.io/IndividualDisplacements.jl/dev/) : simulation and analysis of materials moving through oceanic and atmospheric flows.
- [MITgcmTools.jl](https://juliaclimate.github.io/MiTgcmTools.jl/dev/) : framework to interact with MITgcm (setup, run, output, plot, etc) and ECCO.

## Marine Ecosystems

The [CBIOMES-global](https://github.com/CBIOMES/global-ocean-model) climatology (alpha version) is a global ocean state estimate that covers the period from 1992 to 2011 (ECCO). It is based on Forget et al 2015 for ocean physics [MIT general circulation model](https://mitgcm.readthedocs.io/en/latest/#) and on Dutkiewicz et al 2015 for marine biogeochemistry and ecosystems [Darwin Project model](https://darwin3.readthedocs.io/en/latest/phys_pkgs/darwin.html).

- [CBIOMES\_climatology\_plot](CBIOMES_climatology_plot.html) (➭ [code link](https://raw.githubusercontent.com/gaelforget/OceanStateEstimation.jl/master/examples/CBIOMES/CBIOMES_climatology_plot.jl)) : visualize the climatology maps interactively using [Pluto.jl](https://github.com/fonsp/Pluto.jl/wiki/🔎-Basic-Commands-in-Pluto)
- [CBIOMES\_climatology\_create](https://gaelforget.github.io/OceanStateEstimation.jl/v0.1.13/examples/CBIOMES_model_climatogy.html) (➭ [code link](https://raw.githubusercontent.com/gaelforget/OceanStateEstimation.jl/master/examples/CBIOMES/CBIOMES_climatology_create.jl)) : recreate the climatology file. The original is archived [here in zenodo](https://doi.org/10.5281/zenodo.5598417).

For more notebooks involving [CBIOMES](https://cbiomes.org) and related efforts :

- [Marine Ecosystem Notebooks](https://github.com/JuliaOcean/MarineEcosystemNotebooks) : Darwin Model, Ocean Color data, Gradients field program, and more.
- [PlanktonIndividuals.jl](https://juliaocean.github.io/PlanktonIndividuals.jl/dev/) : simulate the life cycle of ocean phytoplankton cells and their environment.
- [OptimalTransport\_demo.jl](https://gaelforget.github.io/OceanStateEstimation.jl/dev/examples/OptimalTransport_demo.html) : application of optimal transport to compare model and data.
- [JuliaCon2021 workshop](https://github.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl) : _Modeling Marine Ecosystems At Multiple Scales Using Julia_.

To retrieve this climatology, in the `julia REPL` for example :

```julia
using OceanStateEstimation, NCTiles
OceanStateEstimation.CBIOMESclim_download()
fil_out=joinpath(CBIOMESclim_path,"CBIOMES-global-alpha-climatology.nc")
nc=NCTiles.NCDataset(fil_out,"r")
```

## References

- OCCA : [Forget 2010]()
- ECCO v4 : [Forget et al 2015](https://gmd.copernicus.org/articles/8/3071/2015/)
- CBIOMES-global : [Forget 2018](https://zenodo.org/record/2653669#.YbwAUi1h0ow)
	

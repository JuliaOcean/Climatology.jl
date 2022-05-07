
## Physical Oceanography

[ECCO\_standard\_plots.jl](ECCO_standard_plots.html) (➭ [code link](https://raw.githubusercontent.com/gaelforget/OceanStateEstimation.jl/master/examples/ECCO/ECCO_standard_plots.jl)) depicts a selection of climate-relevant variables and indices as an example. These were derived from gridded estimates of the ocean state for physical variables like temperature, salinity, and currents (see `ECCO_standard_analysis.jl`) and retrieved from [dataverse](https://dataverse.harvard.edu/dataverse/ECCO) or [zenodo.org](https://zenodo.org) (see [`OceanStateEstimation.ECCOdiags_download`](@ref), and [`OceanStateEstimation.ECCOdiags_add`](@ref)).

To run the notebook proceed as follows. Once Pluto opens, in your web browser, paste the `code link` url (not the julia code itself) and click open. 

```
using Pluto

import OceanStateEstimation
OceanStateEstimation.ECCOdiags_download()
OceanStateEstimation.ECCOdiags_add("interp_coeffs")
OceanStateEstimation.ECCOdiags_add("release4")

Pluto.run()
```

## Bio-Geo-Chemical Climatology

CBIOMES-global (alpha version) is a global ocean state estimate that covers the period from 1992 to 2011. It is based on Forget et al 2015 for ocean physics MIT general circulation model and on Dutkiewicz et al 2015 for marine biogeochemistry and ecosystems Darwin Project model.

- [CBIOMES\_climatology\_plot](CBIOMES_climatology_plot.html) (➭ [code link](https://raw.githubusercontent.com/gaelforget/OceanStateEstimation.jl/master/examples/CBIOMES/CBIOMES_climatology_plot.jl)) : visualize the climatology maps interactively using [Pluto.jl](https://github.com/fonsp/Pluto.jl/wiki/🔎-Basic-Commands-in-Pluto)
- [CBIOMES\_climatology\_create](https://gaelforget.github.io/OceanStateEstimation.jl/v0.1.13/examples/CBIOMES_model_climatogy.html) (➭ [code link](https://raw.githubusercontent.com/gaelforget/OceanStateEstimation.jl/master/examples/CBIOMES/CBIOMES_climatology_create.jl)) : recreate the climatology file. The original is archived [here in zenodo](https://doi.org/10.5281/zenodo.5598417).
- [OptimalTransport\_demo](https://gaelforget.github.io/OceanStateEstimation.jl/dev/examples/OptimalTransport_demo.html) (➭ [code link](https://raw.githubusercontent.com/gaelforget/OceanStateEstimation.jl/master/examples/OptimalTransport/OptimalTransport_demo.jl)) :  tbd

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
	

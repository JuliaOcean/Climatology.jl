## Path Variables

The gridded fields used in the ECCO examples can be retrieved from [ecco-group.org](https://ecco-group.org/products.htm) and, for the `ECCOv4r2` estimate, from [Harvard Dataverse](https://dataverse.harvard.edu) or [zenodo.org](https://zenodo.org). Two monthly climatologies (`ECCOv4r2` and `OCCA`) are also readily available using the `Julia` artifact system as explained below. These can be relatively large files, compared to the package codes, so they are handled `lazily` (only downloaded when needed). More files are handled in similar fashion.

| Artifact path | File Type  | Download Method |
|:----------------|:----------------:|-----------------:|
| ECCOclim_path             | NetCDF              | lazy, by variable, [dataverse](https://dataverse.harvard.edu/dataverse/ECCO?q=&types=dataverses&sort=dateSort&order=desc&page=1) |
| OCCAclim_path             | NetCDF              |lazy, by variable, [dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/RNXA2A) |
| MITPROFclim_path             | binary    | lazy, whole, [zenodo](https://zenodo.org/record/5101243#.YXiEci1h1qs) |
| ECCOdiags_path             | JLD2    | lazy, whole, [zenodo](https://zenodo.org/record/5773401#.YbQmhS1h3Pg) |
| CBIOMESclim_path             | NetCDF    | lazy, whole, [zenodo](https://zenodo.org/record/5598417#.YoW46C-B3MU) |

## Basic Usage

#### ECCO

ECCO climatology files can downloaded using `get_ecco_files`. These files are for version 4 release 2, on the native model grid.

```julia
using OceanStateEstimation, MeshArrays
γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
tmp=OceanStateEstimation.get_ecco_files(γ,"ETAN")
```

Precomputed quantities shown in [ECCO\_standard\_plots.jl](examples/ECCO_standard_plots.html) can be downloaded separately.

```
OceanStateEstimation.ECCOdiags_download()
OceanStateEstimation.ECCOdiags_add("interp_coeffs")
```

### OCCA

```julia
using OceanStateEstimation
get_occa_variable_if_needed("SIarea")
readdir(OCCAclim_path)
```

### CBIOMES

To retrieve the CBIOMES climatology, in the `julia REPL` for example :

```julia
using OceanStateEstimation
OceanStateEstimation.CBIOMESclim_download()
```

And the files, now found in `CBIOMESclim_path`, can then be read using other libraries.

```
using NCDatasets
fil_out=joinpath(CBIOMESclim_path,"CBIOMES-global-alpha-climatology.nc")
nc=NCDataset(fil_out,"r")
```


## Functions Reference

```@autodocs
Modules = [OceanStateEstimation.downloads]
```

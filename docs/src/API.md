## Path Variables

The gridded fields used in the ECCO examples can be retrieved from [ecco-group.org](https://ecco-group.org/products.htm) and, for the `ECCOv4r2` estimate, from [Harvard Dataverse](https://dataverse.harvard.edu) or [zenodo.org](https://zenodo.org). Two monthly climatologies (`ECCOv4r2` and `OCCA`) are also readily available using the `Scratch.jl` artifact system as explained below. These can be relatively large files, compared to the package codes, so they are handled `lazily` (only downloaded when needed).

| Artifact Name | File Type  | Download Method |
|:----------------|:----------------:|-----------------:|
| ECCO clim             | NetCDF              | lazy, by variable, [dataverse](https://dataverse.harvard.edu/dataverse/ECCO?q=&types=dataverses&sort=dateSort&order=desc&page=1) |
| ECCO diags             | JLD2    | lazy, whole, [zenodo](https://zenodo.org/record/5773401#.YbQmhS1h3Pg) |
| MITprof clim             | binary    | lazy, whole, [zenodo](https://zenodo.org/record/5101243#.YXiEci1h1qs) |
| OCCA clim             | NetCDF              |lazy, by variable, [dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/RNXA2A) |
| CBIOMES clim             | NetCDF    | lazy, whole, [zenodo](https://zenodo.org/record/5598417#.YoW46C-B3MU) |

## Basic Usage

#### ECCO

ECCO climatology files can downloaded using `get_ecco_files`. These files are for version 4 release 2, on the native model grid.

```@example 1
using OceanStateEstimation, MeshArrays
γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
tmp=OceanStateEstimation.get_ecco_files(γ,"ETAN")
```

Precomputed quantities shown in [ECCO\_standard\_plots.jl](examples/ECCO_standard_plots.html) can be downloaded separately.

```@example 1
OceanStateEstimation.ECCOdiags_add("release2")
readdir(ScratchSpaces.ECCO)
```

### OCCA

```@example 1
get_occa_variable_if_needed("SIarea")
readdir(ScratchSpaces.OCCA)
```

### CBIOMES

To retrieve the CBIOMES climatology, in the `julia REPL` for example :

```@example 1
OceanStateEstimation.CBIOMESclim_download()
readdir(ScratchSpaces.CBIOMES)
```

And the files, now found in `ScratchSpaces.CBIOMES`, can then be read using other libraries.

```@example 1
using NCDatasets
fil_out=joinpath(ScratchSpaces.CBIOMES,"CBIOMES-global-alpha-climatology.nc")
nc=NCDataset(fil_out,"r")
keys(nc)
```


## Functions Reference

```@autodocs
Modules = [OceanStateEstimation.downloads]
```

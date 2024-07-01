## Intro

Climatologies are readily downloaded and accessed using the [Scratch.jl](https://github.com/JuliaPackaging/Scratch.jl#readme) artifact system as explained below. 

## Use Examples

### ECCO

ECCO climatology files can downloaded using `get_ecco_files`. These files are for version 4 release 2, on the native model grid.

```@example 1
using Climatology
get_ecco_variable_if_needed("ETAN")

using MeshArrays, MITgcm, NetCDF
path=joinpath(ScratchSpaces.ECCO,"ETAN/ETAN")
γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
tmp=read_nctiles(path,"ETAN",γ,I=(:,:,1))
```

Precomputed quantities shown in [ECCO\_standard\_plots.jl](examples/ECCO_standard_plots.html) can be downloaded separately.

```@example 1
Climatology.ECCOdiags_add("release2")
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
withenv("DATADEPS_ALWAYS_ACCEPT"=>true) do
	path_clim1=datadep"CBIOMES-clim1"
	readdir(path_clim1)
end
```

And the files, now found in `datadep"CBIOMES-clim1"`, can then be read using other libraries.

```@example 1
using NCDatasets
path_clim1=datadep"CBIOMES-clim1"
fil=joinpath(path_clim1,"CBIOMES-global-alpha-climatology.nc")
nc=NCDataset(fil,"r")
keys(nc)
```

### MITprof

To retrieve the MITprof climatologies :

```@example 1
withenv("DATADEPS_ALWAYS_ACCEPT"=>true) do
	readdir(datadep"MITprof-clim1")
end
```

## Path Names

Gridded fields are mostly retrieved from [Harvard Dataverse](https://dataverse.harvard.edu). These can be relatively large files, compared to the package codes, so they are handled `lazily` (only downloaded when needed). Precomputed diagnostics have also been archived on [zenodo.org](https://zenodo.org).

| Artifact Name | File Type  | Download Method |
|:----------------|:----------------:|-----------------:|
| `ScratchSpaces.ECCO`             | NetCDF              | lazy, by variable, [dataverse](https://dataverse.harvard.edu/dataverse/ECCO?q=&types=dataverses&sort=dateSort&order=desc&page=1) |
| `ScratchSpaces.ECCO`             | JLD2    | lazy, whole, [zenodo](https://zenodo.org/record/5773401#.YbQmhS1h3Pg) |
| `datadep"MITprof-clim1"`             | binary    | lazy, whole, [zenodo](https://zenodo.org/record/5101243#.YXiEci1h1qs) |
| `ScratchSpaces.OCCA`             | NetCDF              |lazy, by variable, [dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/RNXA2A) |
| `datadep"CBIOMES-clim1"`             | NetCDF    | lazy, whole, [zenodo](https://zenodo.org/record/5598417#.YoW46C-B3MU) |

## Functions Reference

```@autodocs
Modules = [Climatology.downloads]
```

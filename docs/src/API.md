## Data Access

The gridded fields used in the ECCO examples can be retrieved from [ecco-group.org](https://ecco-group.org/products.htm) and, for the `ECCOv4r2` estimate, from [Harvard Dataverse](https://dataverse.harvard.edu) or [zenodo.org](https://zenodo.org). Two monthly climatologies (`ECCOv4r2` and `OCCA`) are also readily available using the `Julia` artifact system as explained below. These can be relatively large files, compared to the package codes, so they are handled `lazily` (only downloaded when needed). 

| Artifact path | File Type  | Download Method |
|:----------------|:----------------:|-----------------:|
| ECCOclim_path             | NetCDF              | lazy, by variable, [dataverse](https://dataverse.harvard.edu/dataverse/ECCO?q=&types=dataverses&sort=dateSort&order=desc&page=1) |
| OCCAclim_path             | NetCDF              |lazy, by variable, [dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/RNXA2A) |
| MITPROFclim_path             | binary    | lazy, whole, [zenodo](https://zenodo.org/record/5101243#.YXiEci1h1qs) |
| ECCOdiags_path             | JLD2    | lazy, whole, [zenodo](https://zenodo.org/record/5773401#.YbQmhS1h3Pg) |

#### Basic Usage

For ECCO :

```julia
using OceanStateEstimation, MeshArrays
γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
tmp=OceanStateEstimation.get_ecco_files(γ,"ETAN")
```

For OCCA :

```julia
using OceanStateEstimation
get_occa_variable_if_needed("SIarea")
readdir(OCCAclim_path)
```

## Functions Reference

```@autodocs
Modules = [OceanStateEstimation,ECCO]
```

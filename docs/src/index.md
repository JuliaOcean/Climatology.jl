# OceanStateEstimation.jl

This package is still at a very early stage of development. Currently it just provides the following functions and path variables.

```@index
```

## Physical Oceanography

A selection of ocean climatologies / state estimates are provided via the artifact system; see:

- `ECCOclim_path` (netcdf; downloaded on demand)
- `OCCAclim_path` (netcdf; downloaded on demand)
- `MITPROFclim_path` (binaries; readily downloaded)

## CBIOMES-global climatology

[See this page](examples/CBIOMES_model_climatogy.html)

```@autodocs
Modules = [OceanStateEstimation]
```
